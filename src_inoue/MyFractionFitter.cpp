#include "MyFractionFitter.h"

using namespace std;

const double MIN_THRE=1e-8;

MyFractionFitter *gFractionFitter=0;

static void fractionFCN( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  if( !gFractionFitter ){
    cout<<"!!!!! FractionFitter do not setted !!!!!"<<endl;
    exit(0);
  }

  f=gFractionFitter->fcn(npar, par);
}

int MyFractionFitter::Fit()
{
  cout<<"===== MyFractionFitter::Fit() ====="<<endl;
  fMinuit-> SetFCN(fractionFCN);
  double arglist[20];
  arglist[0]=3000; arglist[1]=1;
  int ierflag=0;
  fMinuit-> mnexcm("SIMPLEX", arglist, 1, ierflag);
  fMinuit-> mnexcm("MIGARD", arglist, 1, ierflag);
  fMinuit-> mnexcm("HESSE", arglist, 1, ierflag);

  return -1;
}

MyFractionFitter::MyFractionFitter(TH1 *h1, TObjArray *mc) : fMinuit(0), fData(0), fPlot(0)
{
  cout<<"===== MyFractionFitter Constructor START  ====="<<endl;
  fData=h1;
  fPlot=dynamic_cast<TH1*>(h1->Clone(Form("%s_true", h1->GetName())));
  fPlot->Reset();

  for( int i=0; i<mc->GetEntries(); i++ ){
    fMC.push_back(dynamic_cast<TH1*>(mc->At(i)->Clone(Form("%s_mc%d", mc->At(i)->GetName(), i))));
    fPred.push_back(dynamic_cast<TH1*>(mc->At(i)->Clone(Form("%s_mc%d_true", mc->At(i)->GetName(), i))));
    fPred[i]->Reset();
  }
  fFraction.resize(fMC.size());
  fError.resize(fMC.size());
  check();

  gFractionFitter=this;

  if( fMinuit ) delete fMinuit;

  fMinuit=new TMinuit((int)fMC.size());
  fMinuit-> SetPrintLevel(-1);

  double all_mc_ent=0;
  for( int i=0; i<(int)fMC.size(); i++ ) all_mc_ent+=fMC[i]->GetEntries();
  cout<<"> Data entries : "<<fData->GetEntries()<<"  MC sum : "<<all_mc_ent<<endl;

  for( int i=0; i<(int)fMC.size(); i++ ){
    double frac=fMC[i]->GetEntries()/all_mc_ent;
    fFraction.at(i)=frac;
    fError.at(i)=0.01*frac;
    int ierflag=0;
    fMinuit->mnparm(i, (TString)fMC[i]->GetName(), frac, 0.01*frac, 0.0, 2.0*frac, ierflag);
  }

  cout<<"===== MyFractionFitter Constructor FINISH ====="<<endl;
}

MyFractionFitter::~MyFractionFitter()
{
  if( fPlot ) delete fPlot;
  if( fMinuit ) delete fMinuit;
}

double MyFractionFitter::fcn(int npar, double *param)
{
  //  cout<<"    > MyFractionFitter::fcn()"<<endl;
  for( int i=0; i<npar; i++ ){
    fFraction[i]=param[i]*fData->GetEntries()/fMC[i]->GetEntries();
    if( fFraction[i]==0 ){
      cout<<"!!!!! MyFractionFitter  fraction["<<i<<"]==0 !!!!!"<<endl;
      exit(0);
    }
  }

  double result=0;
  for( int bin=1; bin<=fData->GetNbinsX(); bin++ ){
    if( fData->GetBinContent(bin)==0 ){
      double sum=0;
      for( int i=0; i<npar; i++ ){
	if( fMC[i]->GetBinContent(bin)>0 ){
	  fPred[i]->SetBinContent(bin, fMC[i]->GetBinContent(bin)/(1.0+fFraction[i]));
	}
	else{
	  fPred[i]->SetBinContent(bin, 0);
	}
	sum+=fPred[i]->GetBinContent(bin);
      }
      fPlot-> SetBinContent(bin, sum);
    }
    else{
      int n_max=0;
      double contents_max=0;
      double maxFrac=0;
      for( int i=0; i<npar; i++ ){
	if( fFraction[i]>maxFrac ){
	  n_max=1;
	  maxFrac=fFraction[i];
	  contents_max=fMC[i]->GetBinContent(bin);
	}
	else if( fFraction[i]==maxFrac ){
	  n_max++;
	  contents_max+=fMC[i]->GetBinContent(bin);
	}
      }
      double t_min=-1/maxFrac;

      if( contents_max==0 ){
	double A=fData->GetBinContent(bin)/(1.0+maxFrac);
	for( int i=0; i<npar; i++ ){
	  if( fFraction[i]==maxFrac ) continue;
	  A -= fMC[i]->GetBinContent(bin)*fFraction[i]/(maxFrac-fFraction[i]);
	}
	
	if( A>0 ){
	  double sum=0;
	  for( int i=0; i<npar; i++ ){
	    if( fFraction[i]==maxFrac ){
	      fPred[i]->SetBinContent(bin, A);
	    }
	    else{
	      if( fMC[i]->GetBinContent(bin)>0 ){
		fPred[i]->SetBinContent(bin, fMC[i]->GetBinContent(bin)/(1.0+fFraction[i]));
	      }
	      else{
		fPred[i]->SetBinContent(bin, 0);
	      }
	    }
	    sum+=fPred[i]->GetBinContent(bin);
	  }
	  fPlot-> SetBinContent(bin, sum);
	}
	else{
	  FindPrediction(bin, t_min);
	}
      }
      FindPrediction(bin, t_min);
    }

    result-=fPlot->GetBinContent(bin);
    if( fData-> GetBinContent(bin)>0 && fPlot->GetBinContent(bin)>0 ){
      result+=fData->GetBinContent(bin)*log(fPlot->GetBinContent(bin));
    }
    for( int i=0; i<npar; i++ ){
      result-=fMC[i]->GetBinContent(bin);
      if( fMC[i]->GetBinContent(bin) && fPred[i]->GetBinContent(bin) ){
	result+=fMC[i]->GetBinContent(bin)*log(fPred[i]->GetBinContent(bin));
      }
    }
  }

  return result;
}

double MyFractionFitter::func(int bin, double x)
{
  double val=fData->GetBinContent(bin)/(1-x);
  for( int i=0; i<(int)fMC.size(); i++ ){
    val-=fMC[i]->GetBinContent(bin)/(x+1.0/fFraction[i]);
  }
  return val;
}

void MyFractionFitter::FindPrediction(int bin, const double min)
{
  const double max=1;
  const double min_step=MIN_THRE*0.1*(max-min);
  double step=0.01*(max-min);

  double x=min;
  int sign=1;
  for( int ite=0; ite<=100000; ite++ ){
    double val1=func(bin, x);
    x+=sign*step;
    double val2=func(bin, x);

    if( val1*val2<0 ){
      //      cout<<"Sign change  ite="<<ite<<"  step="<<step<<" x="<<x<<"  val1="<<val1<<"  val2="<<val2<<endl;
      if( step<min_step || (fabs(val1)<MIN_THRE && fabs(val2)<MIN_THRE) ){
	// cout<<"      x="<<x<<"  step="<<step<<"  min="<<min<<"  max="<<max<<endl;
	// cout<<"            val1="<<val1<<"  val2="<<val2<<endl;
	break;
      }
      step*=0.1;
      sign*=-1;
    }

    if( ite==100000 ){
      cout<<"!!!!! MyFractionFitter::FindPrediction Iteration>100000 !!!!!"<<endl;
      cout<<"      x="<<x<<"  step="<<step<<"  min="<<min<<"  max="<<max<<endl;
      cout<<"            val1="<<val1<<"  val2="<<val2<<endl;
      exit(0);
    }
  }

  double sum=0;
  for( int i=0; i<(int)fMC.size(); i++ ){
    if( fMC[i]->GetBinContent(bin)>0 ){
      fPred[i]->SetBinContent(bin, fMC[i]->GetBinContent(bin)/(1.0+fFraction[i]));
    }
    else{
      fPred[i]->SetBinContent(bin, 0);
    }
    sum+=fPred[i]->GetBinContent(bin);
  }
  fPlot-> SetBinContent(bin, sum);
}


void MyFractionFitter::check()
{
  {
    double entries=0;
    for( int bin=1; bin<=(int)fData->GetNbinsX(); bin++ ){
      entries+=fData->GetBinContent(bin);
    }
    if( fData->GetEntries()!=entries ){
      cout<<"  !!! Data Entries not match !!!"<<endl;
      cout<<"      "<<fData->GetEntries()<<"  "<<entries<<endl;
      fData->SetEntries(entries);
    }

  }
  for( int i=0; i<(int)fMC.size(); i++ ){
    if( fData->GetNbinsX()!=fMC[i]->GetNbinsX() ){
      cout<<"!!!!! MyFractionFitter::check not match nbin !!!!!"<<endl;
      exit(-1);
    }

    double entries=0;
    for( int bin=1; bin<=(int)fMC[i]->GetNbinsX(); bin++ ) entries+=fMC[i]->GetBinContent(bin);
    if( fMC[i]->GetEntries()!=entries ){
      cout<<"  !!! MC["<<i<<"] Entries not match !!!"<<endl;
      cout<<"      "<<fMC[i]->GetEntries()<<"  "<<entries<<endl;
      fMC[i]->SetEntries(entries);
    }
  }
}
