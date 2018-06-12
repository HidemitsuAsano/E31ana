#include "FPAnaMan.h"

using namespace std;

FPAnaMan::FPAnaMan()
{
}

void FPAnaMan::makeCS()
{
  double x[fBinMin.size()];
  double xerr[fBinMin.size()];
  double ls_0[fBinMin.size()];
  double ls_1[fBinMin.size()];
  double ls_err_0[fBinMin.size()];
  double ls_err_1[fBinMin.size()];
  double ls_err_0_stat[fBinMin.size()];
  double ls_err_1_stat[fBinMin.size()];

  double tot_err=sqrt(pow(err_of_beam/num_of_beam, 2)
                      +pow(err_DAQ/eff_DAQ, 2)
                      +2*pow(err_CDC/eff_CDC, 2)
                      +pow(err_FDC1/eff_FDC1,2)+
                      +pow(err_trigC/eff_trigC, 2));
  //  double tot_err=0;

  double tot_factor=barn*micro/(num_of_beam*num_of_target*eff_DAQ*eff_CDC*eff_CDC*eff_FDC1);
  //  double tot_factor=barn*micro/(num_of_beam*num_of_target*eff_DAQ*eff_CDC*eff_CDC*eff_trigC*eff_FDC1);
  //  std::cout<<"Tot eff : "<<eff_DAQ*eff_CDC*eff_CDC*eff_trigC*eff_FDC1<<" +- "<<eff_DAQ*eff_CDC*eff_CDC*eff_trigC*eff_FDC1*tot_err<<std::endl;
  //  std::cout<<"err : "<<tot_err<<std::endl;
  fParamInfo->setVal(0, "CS_factor", tot_factor);
  fParamInfo->setVal(0, "CS_err", tot_err);
  for( int i=0; i<fBinMin.size(); i++ ){
    fHistInfo->setVal(i, "NumData0", fNumData[0][i]);
    fHistInfo->setVal(i, "NumData1", fNumData[1][i]);
  }

  for( int i=0; i<fBinMin.size(); i++ ){
    x[i]=0.5*(fBinMax[i]+fBinMin[i]);
    xerr[i]=0.5*(fBinMax[i]-fBinMin[i]);

    if( fEff[0][i]>0 && fSolidAng[0][i] && fNumData[0][i]>0 ){
      ls_0[i]=fNumData[0][i]*fTrig[0][i]/(fEff[0][i]*fSolidAng[0][i]);
      ls_0[i]=fNumData[0][i]*fTrig[0][i]/(fEff[0][i]*fSolidAng[0][i]);
      ls_0[i]*=1./(num_of_beam*num_of_target);
      ls_0[i]*=1./(eff_DAQ*eff_CDC*eff_CDC);
      ls_0[i]*=barn*micro;
      ls_0[i]*=1./eff_FDC1;
      //      ls_0[i]*=tot_factor;
      ls_0[i]*=1./(1000.*(fBinMax[i]-fBinMin[i]));
      ls_err_0[i]=ls_0[i]*sqrt(1./fNumData[0][i]+1./fEff[0][i]+tot_err*tot_err);
      ls_err_0_stat[i]=ls_0[i]*sqrt(1./fNumData[0][i]+1./fEff[0][i]);
    }
    else{
      ls_0[i]=0;
      ls_err_0[i]=0;
    }

    if( fEff[1][i]>0 && fSolidAng[1][i] && fNumData[1][i]>0 ){
      ls_1[i]=fNumData[1][i]*fTrig[1][i]/(fEff[1][i]*fSolidAng[1][i]);
      ls_1[i]*=1./(num_of_beam*num_of_target);
      ls_1[i]*=1./(eff_DAQ*eff_CDC*eff_CDC);
      ls_1[i]*=barn*micro;
      ls_1[i]*=1./eff_FDC1;
      //      ls_1[i]*=tot_factor;
      ls_1[i]*=1./(1000.*(fBinMax[i]-fBinMin[i]));
      ls_err_1[i]=ls_1[i]*sqrt(1./fNumData[1][i]+1./fEff[1][i]+tot_err*tot_err);
      ls_err_1_stat[i]=ls_1[i]*sqrt(1./fNumData[1][i]+1./fEff[1][i]);
    }
    else{
      ls_1[i]=0;
      ls_err_1[i]=0;
    }
  }

  TGraphErrors *gra0=new TGraphErrors(fBinMin.size(), x, ls_0, xerr, ls_err_0);
  gra0-> Write("CrossSection_sim0");
  TGraphErrors *gra0_stat=new TGraphErrors(fBinMin.size(), x, ls_0, xerr, ls_err_0_stat);
  gra0_stat-> Write("CrossSection_sim0_stat");
  TGraphErrors *gra1=new TGraphErrors(fBinMin.size(), x, ls_1, xerr, ls_err_1);
  gra1-> Write("CrossSection_sim1");
  TGraphErrors *gra1_stat=new TGraphErrors(fBinMin.size(), x, ls_1, xerr, ls_err_1_stat);
  gra1_stat-> Write("CrossSection_sim1_stat");

}

void FPAnaMan::makeLineShape()
{
  double x[fBinMin.size()];
  double xerr[fBinMin.size()];
  double ls_0[fBinMin.size()];
  double ls_1[fBinMin.size()];
  double ls_err_0[fBinMin.size()];
  double ls_err_1[fBinMin.size()];

  for( int i=0; i<fBinMin.size(); i++ ){
    x[i]=0.5*(fBinMax[i]+fBinMin[i]);
    xerr[i]=0.5*(fBinMax[i]-fBinMin[i]);

    if( fEff[0][i]>0 && fSolidAng[0][i] ){
      ls_0[i]=fNumData[0][i]*fTrig[0][i]/(fEff[0][i]*fSolidAng[0][i]);
      ls_err_0[i]=ls_0[i]*(1./sqrt(fNumData[0][i])+1./sqrt(fEff[0][i]));
    }
    else{
      ls_0[i]=0;
      ls_err_0[i]=0;
    }

    if( fEff[1][i]>0 && fSolidAng[1][i] ){
      ls_1[i]=fNumData[1][i]*fTrig[1][i]/(fEff[1][i]*fSolidAng[1][i]);
      ls_err_1[i]=ls_1[i]*(1./sqrt(fNumData[1][i])+1./sqrt(fEff[1][i]));
    }
    else{
      ls_1[i]=0;
      ls_err_1[i]=0;
    }
  }

  TGraphErrors *gra0=new TGraphErrors(fBinMin.size(), x, ls_0, xerr, ls_err_0);
  gra0-> Write("LineShape_sim0");
  TGraphErrors *gra1=new TGraphErrors(fBinMin.size(), x, ls_1, xerr, ls_err_1);
  gra1-> Write("LineShape_sim1");
}

void FPAnaMan::makeLineShape_2step()
{
  double x[fBinMin.size()];
  double xerr[fBinMin.size()];
  double num[fBinMin.size()];
  double err[fBinMin.size()];
  double ls[fBinMin.size()];
  double ls_err[fBinMin.size()];

  TH1F *h1 = (TH1F*)f_MC[0]->Get("KP_MM_pimS0");
  for( int bin=1; bin<=h1->GetNbinsX(); bin++ ){
    double mm=h1->GetBinCenter(bin);
    int mm_index=getMMIndex(mm);
    if( mm_index<0 ) continue;
    num[mm_index]+=h1->GetBinContent(bin);
  }

  for( int i=0; i<fBinMin.size(); i++ ){
    x[i]=0.5*(fBinMax[i]+fBinMin[i]);
    xerr[i]=0.5*(fBinMax[i]-fBinMin[i]);

    if( fEff[0][i]>0 && fSolidAng[0][i] ){
      ls[i]=num[i]*fTrig[0][i]/(fEff[0][i]*fSolidAng[0][i]);
      ls_err[i]=ls[i]*(1./sqrt(num[i])+1./sqrt(fEff[0][i]));
    }
  }
  TGraphErrors *gra0=new TGraphErrors(fBinMin.size(), x, ls, xerr, ls_err);
  gra0-> Write("LineShape_pimS0_offshell");


}

void FPAnaMan::readData_tmp()
{
  TH1F *h1 = (TH1F*)f_data->Get("KP_MM_pimL");
  for( int bin=1; bin<=h1->GetNbinsX(); bin++ ){
    double mm=h1->GetBinCenter(bin);
    int mm_index=getMMIndex(mm);
    if( mm_index<0 ) continue;
    fNumData[0][mm_index]+=h1->GetBinContent(bin);
  }

  h1 = (TH1F*)f_data->Get("KP_MM_pimS0");
  for( int bin=1; bin<=h1->GetNbinsX(); bin++ ){
    double mm=h1->GetBinCenter(bin);
    int mm_index=getMMIndex(mm);
    if( mm_index<0 ) continue;
    fNumData[1][mm_index]+=h1->GetBinContent(bin);
  }

  // for( int i=0; i<fBinMin.size(); i++ ){
  //   cout<<fNumData[0][i]<<"  "<<fNumData[1][i]<<endl;
  // }

}

void FPAnaMan::fillHist()
{
  TTree *tree =(TTree*)f_data->Get("ppimpim_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);

  cout<<"===== p pi- pi- event Read entry : "<<tree->GetEntries()<<" ====="<<endl;
  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    tree-> GetEntry(ev);
    TLorentzVector tgt_lmom;
    tgt_lmom.SetVectM(TVector3(0., 0., 0.), dMass);
    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();

    ForwardChargeInfo *fcInfo=anaInfo->forwardCharge(0);
    TLorentzVector fp_lmom=fcInfo->lmom();

    CDSInfo *pim0=anaInfo->CDS(CDS_PiMinus, 0);
    CDSInfo *pim1=anaInfo->CDS(CDS_PiMinus, 1);
    if( pim0->dca()>pim1->dca() ) swap(pim0, pim1);

    TLorentzVector pim_lmom0=pim0->lmom();
    TLorentzVector pim_lmom1=pim1->lmom();

    double kp_mm= (beam_lmom+tgt_lmom-fp_lmom).M();
    double kppimpim_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom0-pim_lmom1).M();
    double kppim0_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom0).M();
    double kppim1_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom1).M();

    bool mmL_flag=false, mmS0_flag=false, mmP_flag=false, mmPgamma_flag=false;
    if( FC_mmL_MIN<kppim0_mm && kppim0_mm<FC_mmL_MAX ) mmL_flag=true;
    if( FC_mmL_MIN<kppim1_mm && kppim1_mm<FC_mmL_MAX ) mmL_flag=true;
    if( FC_mmS0_MIN<kppim0_mm && kppim0_mm<FC_mmS0_MAX ) mmS0_flag=true;
    if( FC_mmS0_MIN<kppim1_mm && kppim1_mm<FC_mmS0_MAX ) mmS0_flag=true;
    if( FC_mmP_MIN<kppimpim_mm && kppimpim_mm<FC_mmP_MAX ) mmP_flag=true;
    if( 0.99<kppimpim_mm && kppimpim_mm<1.06 ) mmPgamma_flag=true;

    if( mmL_flag && mmP_flag ){
      int mm_index=getMMIndex(kp_mm);
      if( mm_index>=0 ) fNumData[0][mm_index]+=1.0;
    }
    if( mmS0_flag && mmPgamma_flag ){
      int mm_index=getMMIndex(kp_mm);
      if( mm_index>=0 ) fNumData[1][mm_index]+=1.0;
    }
  }
}

void FPAnaMan::readSim(int index)
{
  TFile *f=0;
  if( index==0 ) f=f_sim1;
  else if( index==1 ) f=f_sim2;
  else{
    cout<<"!!!!! Invaild index="<<index<<" !!!!!"<<endl;
    exit(0);
  }

  for( int i=0; i<fBinMin.size(); i++ ){
    fTrig[index][i]=0;
    fEff[index][i]=0;
  }

  TNtuple *tup=(TNtuple*)f->Get("tup_trig_mass_fp");
  Float_t mm;
  tup-> SetBranchAddress("m", &mm);
  for( int i=0; i<tup->GetEntries(); i++ ){
    tup-> GetEntry(i);
    //    cout<<mm<<endl;
    int mm_index = getMMIndex(mm);
    if( mm_index>=0 ) fTrig[index][mm_index]+=1.0;
  }

  tup=(TNtuple*)f->Get("tup_eff_mass_fp_pipi_mm2");
  tup-> SetBranchAddress("m", &mm);
  for( int i=0; i<tup->GetEntries(); i++ ){
    tup-> GetEntry(i);
    //    cout<<mm<<endl;
    int mm_index = getMMIndex(mm);
    if( mm_index>=0 ) fEff[index][mm_index]+=1.0;
  }

  // for( int i=0; i<fBinMin.size(); i++ ){
  //   cout<<"sim"<<index<<" : "<<fTrig[index][i]<<"  "<<fEff[index][i]<<endl;
  // }

  tup=(TNtuple*)f->Get("tup_acc_fp");
  Float_t sa_x, sa_y, flag;
  tup-> SetBranchAddress("m", &mm);
  tup-> SetBranchAddress("SA_x", &sa_x);
  tup-> SetBranchAddress("SA_y", &sa_y);
  tup-> SetBranchAddress("flag", &flag);

  for( int i=0; i<tup->GetEntries(); i++ ){
    tup-> GetEntry(i);
    //    cout<<mm<<"  sa : "<<sa_x<<", "<<sa_y<<" flag : "<<flag<<endl;
    int mm_index=getMMIndex(mm);
    if( mm_index<0 ) continue;

    TH2F *h2 = (TH2F*)f_out->Get(Form("Gen_%d_sim%d", mm_index, index));
    h2-> Fill(sa_x, sa_y);
    if( flag>0.5 ){
      TH2F *h2 = (TH2F*)f_out->Get(Form("Acc_%d_sim%d", mm_index, index));
      h2-> Fill(sa_x, sa_y);
    }
  }

  double x[fBinMin.size()];
  double xerr[fBinMin.size()];
  double y[fBinMin.size()];
  for( int mm_index=0; mm_index<fBinMin.size(); mm_index++ ){
    x[mm_index]=0.5*(fBinMax[mm_index]+fBinMin[mm_index]);
    xerr[mm_index]=0.5*(fBinMax[mm_index]-fBinMin[mm_index]);
    y[mm_index]=0;

    TH2F *h2_1 = (TH2F*)f_out->Get(Form("Gen_%d_sim%d", mm_index, index));
    TH2F *h2_2 = (TH2F*)f_out->Get(Form("Acc_%d_sim%d", mm_index, index));
    TH2F *h2_3 = (TH2F*)f_out->Get(Form("Acceptance_%d_sim%d", mm_index, index));
    
    for( int binx=1; binx<=h2_1->GetNbinsX(); binx++ ){
      for( int biny=1; biny<=h2_1->GetNbinsY(); biny++ ){
	double gen=h2_1->GetBinContent(binx, biny);
	double acc=h2_2->GetBinContent(binx, biny);

	if( gen>0 ){	
	  if( acc/gen>0.5 ){
	    y[mm_index]+=h2_1->GetXaxis()->GetBinWidth(binx)*h2_1->GetYaxis()->GetBinWidth(biny);
	    fSolidAng[index][mm_index]+=h2_1->GetXaxis()->GetBinWidth(binx)*h2_1->GetYaxis()->GetBinWidth(biny);
	  }
	  h2_3-> SetBinContent(binx, biny, acc/gen);
	}
      }
    }
    fHistInfo->setVal(mm_index, "SolidAng0", fSolidAng[0][mm_index]);
    fHistInfo->setVal(mm_index, "SolidAng1", fSolidAng[1][mm_index]);
  }
  TGraphErrors *gra=new TGraphErrors(fBinMin.size(), x, y, xerr, 0);
  gra-> Write(Form("gra_FC_sa_sim%d", index));
}

void FPAnaMan::initHist(int index)
{
  for( int i=0; i<fBinMin.size(); i++ ){
    new TH2F(Form("Gen_%d_sim%d", i, index), "", 200, -0.5, 0.5, 200, -0.5, 0.5); 
    new TH2F(Form("Acc_%d_sim%d", i, index), "", 200, -0.5, 0.5, 200, -0.5, 0.5); 
    new TH2F(Form("Acceptance_%d_sim%d", i, index), "", 200, -0.5, 0.5, 200, -0.5, 0.5); 
  }
}

void FPAnaMan::makeAcc()
{
  double x[fBinMin.size()];
  double xerr[fBinMin.size()];
  double y[fBinMin.size()];
  double yerr[fBinMin.size()];
  for( int index=0; index<2; index++ ){
    for( int i=0; i<fBinMin.size(); i++ ){
      x[i]=0.5*(fBinMin[i]+fBinMax[i]);
      xerr[i]=0.5*(fBinMax[i]-fBinMin[i]);

      if( fTrig[index][i]>0 ){      
	y[i]=fEff[index][i]/fTrig[index][i];
	yerr[i]=y[i]/sqrt(fEff[index][i]);
      }
      else{
	y[i]=0;
	yerr[i]=0;
      }
    }
    TGraphErrors *gra=new TGraphErrors(fBinMin.size(), x, y, xerr, yerr);
    gra-> Write(Form("gra_Acc%d", index));
  }
  for( int i=0; i<fBinMin.size(); i++ ){
    fHistInfo-> setVal(i, "Trig0", fTrig[0][i]);
    fHistInfo-> setVal(i, "Trig1", fTrig[1][i]);
    fHistInfo-> setVal(i, "Eff0", fEff[0][i]);
    fHistInfo-> setVal(i, "Eff1", fEff[1][i]);
  }
}

void FPAnaMan::init(TString conffile)
{
  ifstream ifs(conffile.Data());

  string str, str2;
  string outfilename;
  int num;

  while( getline(ifs, str) ){
    if( str[0]=='#' ) continue;

    if( str.find("DataFile:")==0 ){
      stringstream ss(str); ss>>str2; ss>>str2;
      f_data =new TFile(str2.c_str());
    }
    if( str.find("Sim1:")==0 ){
      stringstream ss(str); ss>>str2; ss>>str2;
      f_sim1 =new TFile(str2.c_str());
    }
    if( str.find("Sim2:")==0 ){
      stringstream ss(str); ss>>str2; ss>>str2;
      f_sim2 =new TFile(str2.c_str());
    }
    if( str.find("MC")==0 ){
      stringstream ss(str); ss>>str2; ss>>str2;
      f_MC.push_back(new TFile(str2.c_str()));
    }

    if( str.find("OutRoot:")==0 ){
      stringstream ss(str); ss>>str2; ss>>str2;
      outfilename=str2;
    }

    if( str.find("NBin:")==0 ){
      stringstream ss(str); ss>>str2; ss>>num;

      double min, max;
      for( int i=0; i<num; i++ ){
        getline(ifs, str);
        if( sscanf(str.c_str(), "%lf %lf", &min, &max)==2 ){
          fBinMin.push_back(min);
          fBinMax.push_back(max);

	  fNumData[0].push_back(0);
	  fNumData[1].push_back(0);

	  fTrig[0].push_back(0);
	  fTrig[1].push_back(0);
	  fEff[0].push_back(0);
	  fEff[1].push_back(0);
	  fSolidAng[0].push_back(0);
	  fSolidAng[1].push_back(0);
	}
      }
    }
  }

  f_out=new TFile(outfilename.c_str(), "recreate");
  std::vector<double> bin=fBinMin;
  bin.push_back(fBinMax[fBinMax.size()-1]);

  std::vector<TString> names;
  names.push_back("NumData0");
  names.push_back("NumData1");
  names.push_back("Trig0");
  names.push_back("Eff0");
  names.push_back("Trig1");
  names.push_back("Eff1");
  names.push_back("SolidAng0");
  names.push_back("SolidAng1");
  fHistInfo=new HistInfo(bin, names);

  std::vector<double> tmp; tmp.push_back(0); tmp.push_back(1);
  std::vector<TString> paramNames;
  paramNames.push_back("CS_factor");
  paramNames.push_back("CS_err");
  fParamInfo=new HistInfo(tmp, paramNames);

  initHist(0);
  initHist(1);

  readSim(0);
  readSim(1);

  //  fillHist();
  readData_tmp();
  makeAcc();
  makeLineShape();
  makeCS();
}

void FPAnaMan::finit()
{
  cout<<"===== FPAnaMan::finit start ====="<<endl;
  fHistInfo->dump();
  fHistInfo->Write("info_kp");
  fParamInfo->Write("info_param");
  f_out-> Write();
  f_out-> Close();
  cout<<"===== FPAnaMan::finit finish ====="<<endl;
}

int FPAnaMan::getMMIndex(double mm)
{
  for( int i=0; i<fBinMin.size(); i++ ){
    if( fBinMin[i]<=mm && mm<fBinMax[i] ) return i;
  }

  return -1;
}
