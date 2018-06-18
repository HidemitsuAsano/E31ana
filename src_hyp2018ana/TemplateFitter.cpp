#include "TemplateFitter.h"

using namespace std;

TemplateFitter::TemplateFitter()
  : fFitter(0), fChi2(-1.0), fNDF(-1), fData(0), fMC(0)
{ 
}

TemplateFitter::TemplateFitter(TH1 *data, TObjArray *mc)
  : fFitter(0), fChi2(-1.0), fNDF(-1), fData(data), fMC(mc)
{
  fName.resize(fMC->GetEntries(), "");
  fFraction.resize(fMC->GetEntries(), -1.0);
  fError.resize(fMC->GetEntries(), -1.0);
  fFix.resize(fMC->GetEntries(), false);

  int ndf=0;
  for( int binx=1; binx<=fData->GetNbinsX(); binx++ ){
    for( int biny=1; biny<=fData->GetNbinsY(); biny++ ){
      for( int binz=1; binz<=fData->GetNbinsZ(); binz++ ){
	bool flag=false;
	if( fData->GetBinContent(binx, biny, binz)>0 ) flag=true;

	for( int i=0; i<fMC->GetEntries(); i++ ){
	  TH1* h1=(TH1*)fMC->At(i);
	  if( h1->GetBinContent(binx, biny, binz)>0 ) flag=true;
	}
	if( flag ) ndf++;
      }
    }
  }
  fNDF=ndf;
}

bool TemplateFitter::getResult(int i, double &p, double &e)
{
  if( !fIsFitted ) return false;

  p=fFraction[i];
  e=fError[i];
  return true;
}

bool TemplateFitter::getScaleErr(int i, double &scale, double &err)
{
  if( !fIsFitted ) return false;

  double nData=getEntries(fData);
  double nHist=getEntries((TH1*)fMC->At(i));

  scale=nData*fFraction[i]/nHist;
  err=nData*fError[i]/nHist;
  return true;
}

void TemplateFitter::setScale(int index, double val)
{
  double nData=getEntries(fData);
  double nHist=getEntries((TH1*)fMC->At(index));

  fFraction[index]=nHist*val/nData;
}

bool TemplateFitter::fit(int dump_level)
{
  if( dump_level>0 ){
    std::cout<<"===== TemplateFitter::Fit START ====="<<std::endl;
    std::cout<<"      dump level : "<<dump_level<<endl;
  }
  int nFitHist=0;
  vector<int> index_map(fMC->GetEntries(), 0);
  for( int i=0; i<fMC->GetEntries(); i++ ){
    if( fFraction[i]==0 ) continue;

    // cout<<i<<" "<<fFraction[i]<<endl;
    // cout<<" Entry:"<<getEntries(dynamic_cast<TH1*>(fMC->At(i)))<<endl;
    if( getEntries(dynamic_cast<TH1*>(fMC->At(i)))>0 ){
      index_map[nFitHist]=i;
      nFitHist++;
    }    
  }
  if( nFitHist==0 ){
    if( dump_level>0 ) cout<<"  !!! non fit hist !!!"<<endl;
    return false;
  }

  TObjArray *mc=new TObjArray(nFitHist);
  for( int i=0; i<nFitHist; i++ ) mc-> Add(fMC->At(index_map[i]));

  if( fFitter ) delete fFitter;
  if( dump_level>1 ) fFitter = new TFractionFitter(fData, mc);
  else  fFitter = new TFractionFitter(fData, mc, "Q");
  int nFreeParam=0;
  for( int i=0; i<nFitHist; i++ ){
    if( fFraction[index_map[i]]>0 ){
      TString name=fName[index_map[i]];
      if( name=="" ) name=Form("param%d", i);
      fFitter->GetFitter()->SetParameter(i, name, fFraction[index_map[i]], 0.05, 0., 1.0);
      if( fFix[index_map[i]] ){
	if( fFraction[index_map[i]]<=0 ){
	  if( dump_level>0 ) cout<<"  !!! Fix parameter not set !!!"<<endl;
	  return false;
	}
	fFitter->GetFitter()->FixParameter(i);
      }
      else nFreeParam++;
    }
  }
  gMinuit-> Command("SET STRategy 2");

  if( nFreeParam==0 ){
    if( dump_level>0 ) cout<<"  !!! non free parameter !!!"<<endl;
    return false;
  }

  while( fFitter->Fit()!=0 ){
    if( dump_level>0 ){
      cout<<"> Refitting"<<endl;
    }

    double param[nFreeParam];
    double val=1;
    for( int i=0; i<nFreeParam-1; i++ ){
      val*=gRandom->Uniform();
      param[i]=val;
    }
    param[nFreeParam-1]=1.0-val;

    delete fFitter;
    if( dump_level>1 ) fFitter = new TFractionFitter(fData, mc);
    else  fFitter = new TFractionFitter(fData, mc, "Q");

    int index=0;
    for( int i=0; i<nFitHist; i++ ){
      TString name=fName[index_map[i]];
      if( name=="" ) name=Form("param%d", i);
      if( fFix[index_map[i]] ){
	fFitter->GetFitter()->SetParameter(i, name, fFraction[index_map[i]], 0.05, 0., 1.0);
	fFitter->GetFitter()->FixParameter(i);
	continue;
      }

      fFitter->GetFitter()->SetParameter(i, name, param[index], 0.05, 0., 1.0);
      index++;
    }
  }

  if( dump_level>0 ){
    std::cout<<"===== TemplateFitter::Fit scuccessfully FINISH ====="<<std::endl;
  }

  double frac[nFitHist];
  double frac_sum=0;
  for( int i=0; i<nFitHist; i++ ){
    double param, err;
    fFitter-> GetResult(i, param, err);
    frac_sum+=param;
    frac[i]=param;
    if( fFix[index_map[i]] ){
      if( err!=0 && dump_level>0 ){
	cout<<"    ! Fix parameter : "<<i<<" error is not 0"<<endl;
      }
      if( fFraction[index_map[i]]!=param && dump_level>0 ){
	cout<<"    ! Fix parameter : "<<i<<" not match "<<fFraction[index_map[i]]<<"  "<<param<<endl;
      }    
    }
    fFraction[index_map[i]]=param;
  }
  fChi2 = fFitter-> GetChisquare();

  fIsFitted=true;

  for( int i=0; i<nFreeParam; i++ ){
    double drive[nFitHist];
    for( int j=0; j<nFitHist; j++ ){
      if( i==j ) drive[j]=frac[i]/(frac_sum*frac_sum);
      else drive[j]=-(frac_sum-frac[i])/(frac_sum*frac_sum);
    }

    double err2=0;
    for( int j=0; j<nFreeParam; j++ ){
      for( int k=0; k<nFreeParam; k++ ){
	double mat_elem=fFitter->GetFitter()->GetCovarianceMatrixElement(j, k);
	err2+=drive[j]*drive[k]*mat_elem;
      }
    }
    double propagate_err=sqrt(err2);

    double param, err;
    fFitter-> GetResult(i, param, err);
    fError[index_map[i]]=err;
    // if( propagate_err>fFraction[index_map[i]] ){
    //   fFitter-> GetFitter()->PrintResults(1, 0.0);
    // }

    fError[index_map[i]]=propagate_err;
    if( fFix[index_map[i]] ) fError[index_map[i]]=0;
  }


  return true;
}

double TemplateFitter::getEntries(TH1 *h1)
{
  double sum=0;
  for( int binx=1; binx<=h1->GetNbinsX(); binx++ ){
    for( int biny=1; biny<=h1->GetNbinsY(); biny++ ){
      for( int binz=1; binz<=h1->GetNbinsZ(); binz++ ){
	sum+=h1->GetBinContent(binx, biny, binz);
      }
    }
  }
  return sum;
}
