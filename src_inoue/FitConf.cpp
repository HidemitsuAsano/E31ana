#include "FitConf.h"

using namespace std;

void FitConf::fillHist_Sim(const int id)
{
  TFile *f=0;
  if( id==1 ) f=f_sim1;
  if( id==2 ) f=f_sim2;
  TString suffix2=Form("sim%d", id);
  cout<<"===== fill Hist for sim  TFile : "<<f->GetName()<<"  suffix : "<<suffix2<<" ====="<<endl;

  TTree *tree = (TTree*)f->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);
  cout<<"      Entries : "<<tree->GetEntries()<<endl;

  TH1F *h1;
  h1=(TH1F*)f_out->Get(Form("FitIM_%s", suffix2.Data())); h1-> Reset();
  for( int i=0; i<fBinMin.size(); i++ ){
    h1=(TH1F*)f_out->Get(Form("FitKNpi_%s_%d", suffix2.Data(), i)); h1->Reset();
  }

  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    tree-> GetEntry(ev);
    TLorentzVector target_lmom;
    target_lmom.SetVectM(TVector3(0, 0, 0), dMass);

    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TLorentzVector n_lmom=anaInfo->forwardNeutral(0)->lmom();
    TLorentzVector pim_lmom=anaInfo->CDS(CDS_PiMinus, 0)->lmom();
    TLorentzVector pip_lmom=anaInfo->CDS(CDS_PiPlus, 0)->lmom();
    TLorentzVector pipi_lmom=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->lmom();

    double kn_mm=(beam_lmom+target_lmom-n_lmom).M();
    double knpipi_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom-pip_lmom).M();
    double knpim_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom).M();
    double knpip_mm=(beam_lmom+target_lmom-n_lmom-pip_lmom).M();

    double npim_im=(n_lmom+pim_lmom).M();
    double npip_im=(n_lmom+pip_lmom).M();
    double pipi_im=pipi_lmom.M();

    h1=(TH1F*)f_out->Get(Form("FitIM_%s", suffix2.Data()));
    if( !h1 ){
      cout<<"  !!! "<<Form("FitIM_%s", suffix2.Data())<<" !!!"<<endl;
      exit(0);
    }
    double weight=getScale(id, kn_mm);
    h1-> Fill(weight*pipi_im);
    h1-> Fill(weight*npim_im);
    h1-> Fill(1.0+weight*npip_im);

    if( Npipi_N_MIN<knpipi_mm && knpipi_mm<Npipi_N_MAX ){
      bool K0_flag=false; bool Sm_flag=false; bool Sp_flag=false;
      if( Npipi_K0_MIN<pipi_im && pipi_im<Npipi_K0_MAX ) K0_flag=true;
      if( Npipi_Sm_MIN<npim_im && npim_im<Npipi_Sm_MAX ) Sm_flag=true;
      if( Npipi_Sp_MIN<npip_im && npip_im<Npipi_Sp_MAX ) Sp_flag=true;

      if( !K0_flag && !Sm_flag && !Sp_flag ){
	int bin_index=getMMIndex(kn_mm);
	if( bin_index>=0 ){
	  h1=(TH1F*)f_out->Get(Form("FitKNpi_%s_%d", suffix2.Data(), bin_index));
	  if( !h1 ){
	    cout<<"  !!! "<<Form("FitKNpi_%s_%d", suffix2.Data(), bin_index)<<" !!!"<<endl;
	  }
	  h1-> Fill(npim_im-1.0);
	  h1-> Fill(npip_im);
	}
      }
    }
  }
}


void FitConf::fillHist(TFile *f, const char *suffix)
{
  cout<<"===== fill Hist  TFile : "<<f->GetName()<<"  suffix : "<<suffix<<" ====="<<endl;
  TString suffix2=suffix;

  TTree *tree = (TTree*)f->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);
  cout<<"      Entries : "<<tree->GetEntries()<<endl;

  TH1F *h1;
  h1=(TH1F*)f_out->Get(Form("FitIM_%s", suffix2.Data())); h1-> Reset();
  for( int i=0; i<fBinMin.size(); i++ ){
    h1=(TH1F*)f_out->Get(Form("FitKNpi_%s_%d", suffix2.Data(), i)); h1->Reset();
  }

  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    tree-> GetEntry(ev);
    TLorentzVector target_lmom;
    target_lmom.SetVectM(TVector3(0, 0, 0), dMass);

    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TLorentzVector n_lmom=anaInfo->forwardNeutral(0)->lmom();
    TLorentzVector pim_lmom=anaInfo->CDS(CDS_PiMinus, 0)->lmom();
    TLorentzVector pip_lmom=anaInfo->CDS(CDS_PiPlus, 0)->lmom();
    TLorentzVector pipi_lmom=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->lmom();

    double kn_mm=(beam_lmom+target_lmom-n_lmom).M();
    double knpipi_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom-pip_lmom).M();
    double knpim_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom).M();
    double knpip_mm=(beam_lmom+target_lmom-n_lmom-pip_lmom).M();

    double npim_im=(n_lmom+pim_lmom).M();
    double npip_im=(n_lmom+pip_lmom).M();
    double pipi_im=pipi_lmom.M();

    h1=(TH1F*)f_out->Get(Form("FitIM_%s", suffix2.Data()));
    if( !h1 ){
      cout<<"  !!! "<<Form("FitIM_%s", suffix2.Data())<<" !!!"<<endl;
      exit(0);
    }
    h1-> Fill(pipi_im);
    h1-> Fill(npim_im);
    h1-> Fill(1.0+npip_im);

    if( Npipi_N_MIN<knpipi_mm && knpipi_mm<Npipi_N_MAX ){
      bool K0_flag=false; bool Sm_flag=false; bool Sp_flag=false;
      if( Npipi_K0_MIN<pipi_im && pipi_im<Npipi_K0_MAX ) K0_flag=true;
      if( Npipi_Sm_MIN<npim_im && npim_im<Npipi_Sm_MAX ) Sm_flag=true;
      if( Npipi_Sp_MIN<npip_im && npip_im<Npipi_Sp_MAX ) Sp_flag=true;

      if( !K0_flag && !Sm_flag && !Sp_flag ){
	int bin_index=getMMIndex(kn_mm);
	if( bin_index>=0 ){
	  h1=(TH1F*)f_out->Get(Form("FitKNpi_%s_%d", suffix2.Data(), bin_index));
	  if( !h1 ){
	    cout<<"  !!! "<<Form("FitKNpi_%s_%d", suffix2.Data(), bin_index)<<" !!!"<<endl;
	  }
	  h1-> Fill(npim_im-1.0);
	  h1-> Fill(npip_im);
	}
      }
    }
  }
}

void FitConf::readData()
{
  TTree *tree = (TTree*)f_data->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);

  for( int i=0; i<fBinMin.size(); i++ ) fNumData[i]=0;
  fNumDataK0=0;
  fNumDataSm=0;
  fNumDataSp=0;

  cout<<"Read Data Entries : "<<tree->GetEntries()<<endl;
  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    tree-> GetEntry(ev);
    TLorentzVector target_lmom;
    target_lmom.SetVectM(TVector3(0, 0, 0), dMass);
    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TLorentzVector n_lmom=anaInfo->forwardNeutral(0)->lmom();
    TLorentzVector pim_lmom=anaInfo->CDS(CDS_PiMinus, 0)->lmom();
    TLorentzVector pip_lmom=anaInfo->CDS(CDS_PiPlus, 0)->lmom();
    TLorentzVector pipi_lmom=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->lmom();

    double kn_mm=(beam_lmom+target_lmom-n_lmom).M();
    // cout<<"beam : "<<beam_lmom.M()<<endl;
    // cout<<"fn : "<<n_lmom.M()<<endl;
    // cout<<"KN : "<<kn_mm<<endl;
    double knpipi_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom-pip_lmom).M();
    double knpim_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom).M();
    double knpip_mm=(beam_lmom+target_lmom-n_lmom-pip_lmom).M();

    double npim_im=(n_lmom+pim_lmom).M();
    double npip_im=(n_lmom+pip_lmom).M();
    double pipi_im=pipi_lmom.M();

    if( Npipi_N_MIN<knpipi_mm && knpipi_mm<Npipi_N_MAX ){

      bool K0_flag=false; bool Sm_flag=false; bool Sp_flag=false;
      if( Npipi_K0_MIN<pipi_im && pipi_im<Npipi_K0_MAX ) K0_flag=true;
      if( Npipi_Sm_MIN<npim_im && npim_im<Npipi_Sm_MAX ) Sm_flag=true;
      if( Npipi_Sp_MIN<npip_im && npip_im<Npipi_Sp_MAX ) Sp_flag=true;

      if( K0_flag ) fNumDataK0+=1.0;
      if( Sm_flag ) fNumDataSm+=1.0;
      if( Sp_flag ) fNumDataSp+=1.0;

      if( !K0_flag && !Sm_flag && !Sp_flag ){
	int bin_index=getMMIndex(kn_mm);
	if( bin_index>=0 ) fNumData[bin_index] += 1.0;
      }
    }
  }

  for( int i=0; i<fBinMin.size(); i++ ){
    cout<<"d(K-, n) ["<<fBinMin[i]<<"  "<<fBinMax[i]<<"] = "<<fNumData[i]<<endl;
  }

  delete anaInfo;
}

void FitConf::readSim1()
{
  TTree *tree = (TTree*)f_sim1->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);

  for( int i=0; i<fBinMin.size(); i++ ) fNumSim1[i]=0;

  cout<<"Read Sim1 Entries : "<<tree->GetEntries()<<endl;
  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    tree-> GetEntry(ev);
    TLorentzVector target_lmom;
    target_lmom.SetVectM(TVector3(0, 0, 0), dMass);
    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TLorentzVector n_lmom=anaInfo->forwardNeutral(0)->lmom();
    TLorentzVector pim_lmom=anaInfo->CDS(CDS_PiMinus, 0)->lmom();
    TLorentzVector pip_lmom=anaInfo->CDS(CDS_PiPlus, 0)->lmom();
    TLorentzVector pipi_lmom=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->lmom();

    double kn_mm=(beam_lmom+target_lmom-n_lmom).M();
    double knpipi_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom-pip_lmom).M();
    double knpim_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom).M();
    double knpip_mm=(beam_lmom+target_lmom-n_lmom-pip_lmom).M();

    double npim_im=(n_lmom+pim_lmom).M();
    double npip_im=(n_lmom+pip_lmom).M();
    double pipi_im=pipi_lmom.M();

    if( Npipi_N_MIN<knpipi_mm && knpipi_mm<Npipi_N_MAX ){

      bool K0_flag=false; bool Sm_flag=false; bool Sp_flag=false;
      if( Npipi_K0_MIN<pipi_im && pipi_im<Npipi_K0_MAX ) K0_flag=true;
      if( Npipi_Sm_MIN<npim_im && npim_im<Npipi_Sm_MAX ) Sm_flag=true;
      if( Npipi_Sp_MIN<npip_im && npip_im<Npipi_Sp_MAX ) Sp_flag=true;

      if( !K0_flag && !Sm_flag && !Sp_flag ){
	int bin_index=getMMIndex(kn_mm);
	if( bin_index>=0 ) fNumSim1[bin_index] += 1.0;
      }
    }
  }

  for( int i=0; i<fBinMin.size(); i++ ){
    cout<<"d(K-, n) ["<<fBinMin[i]<<"  "<<fBinMax[i]<<"] = "<<fNumSim1[i]<<endl;
  }
  delete anaInfo;
}

void FitConf::readSim2()
{
  TTree *tree = (TTree*)f_sim2->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);

  for( int i=0; i<fBinMin.size(); i++ ) fNumSim2[i]=0;

  cout<<"Read Sim2 Entries : "<<tree->GetEntries()<<endl;
  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    tree-> GetEntry(ev);
    TLorentzVector target_lmom;
    target_lmom.SetVectM(TVector3(0, 0, 0), dMass);
    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TLorentzVector n_lmom=anaInfo->forwardNeutral(0)->lmom();
    TLorentzVector pim_lmom=anaInfo->CDS(CDS_PiMinus, 0)->lmom();
    TLorentzVector pip_lmom=anaInfo->CDS(CDS_PiPlus, 0)->lmom();
    TLorentzVector pipi_lmom=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->lmom();

    double kn_mm=(beam_lmom+target_lmom-n_lmom).M();
    double knpipi_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom-pip_lmom).M();
    double knpim_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom).M();
    double knpip_mm=(beam_lmom+target_lmom-n_lmom-pip_lmom).M();

    double npim_im=(n_lmom+pim_lmom).M();
    double npip_im=(n_lmom+pip_lmom).M();
    double pipi_im=pipi_lmom.M();

    if( Npipi_N_MIN<knpipi_mm && knpipi_mm<Npipi_N_MAX ){

      bool K0_flag=false; bool Sm_flag=false; bool Sp_flag=false;
      if( Npipi_K0_MIN<pipi_im && pipi_im<Npipi_K0_MAX ) K0_flag=true;
      if( Npipi_Sm_MIN<npim_im && npim_im<Npipi_Sm_MAX ) Sm_flag=true;
      if( Npipi_Sp_MIN<npip_im && npip_im<Npipi_Sp_MAX ) Sp_flag=true;

      if( !K0_flag && !Sm_flag && !Sp_flag ){
	int bin_index=getMMIndex(kn_mm);
	if( bin_index>=0 ) fNumSim2[bin_index] += 1.0;
      }
    }
  }

  for( int i=0; i<fBinMin.size(); i++ ){
    cout<<"d(K-, n) ["<<fBinMin[i]<<"  "<<fBinMax[i]<<"] = "<<fNumSim2[i]<<endl;
  }
  delete anaInfo;
}

void FitConf::initHist(char *suffix)
{
  for( int i=0;  i<fBinMin.size(); i++ ){
    int nbin=2000/fRebinF;
    new TH1F(Form("FitKNpi_%s_%d", suffix, i), Form("FitKNpi_%s_%d", suffix, i), nbin, 0, 2.0); // fro d(K-, n pi-) d(K-, n pi+)
  }
  new TH1F(Form("FitIM_%s", suffix), Form("FitIM_%s", suffix), 3000, 0, 3.0);
}

void FitConf::init(char *filename)
{
  ifstream ifs(filename);

  string outfilename;
  string str, str2;
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
    if( str.find("K0:")==0 ){
      stringstream ss(str); ss>>str2; ss>>num;
      for( int i=0; i<num; i++ ){
	getline(ifs, str);
	f_K0.push_back(new TFile(str.c_str()));
      }
    }
    if( str.find("Sm:")==0 ){
      stringstream ss(str); ss>>str2; ss>>num;
      for( int i=0; i<num; i++ ){
	getline(ifs, str);
	f_Sm.push_back(new TFile(str.c_str()));
      }
    }
    if( str.find("Sp:")==0 ){
      stringstream ss(str); ss>>str2; ss>>num;
      for( int i=0; i<num; i++ ){
	getline(ifs, str);
	f_Sp.push_back(new TFile(str.c_str()));
      }
    }
    if( str.find("BG:")==0 ){
      stringstream ss(str); ss>>str2; ss>>num;
      for( int i=0; i<num; i++ ){
	getline(ifs, str);
	f_BG.push_back(new TFile(str.c_str()));
      }
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
	  fNumData.push_back(0);
	  fNumSim1.push_back(0);
	  fNumSim2.push_back(0);

	  fEff1.push_back(0);
	  fEff2.push_back(0);
	  fTrig1.push_back(0);
	  fTrig2.push_back(0);
	  fFrac1.push_back(0);
	  fFrac2.push_back(0);
	  fErr1.push_back(0);	
	  fErr2.push_back(0);	
	  fErr12.push_back(0);
	  fChi2.push_back(0);	
	  fNDF.push_back(0);
	}
	else{ 
	  std::cout<<"  !!! read bin range error !!!"<<std::endl;
	  exit(0);
	}
      }
    }
    if( str.find("Rebin:")==0 ){
      stringstream ss(str); ss>>str2; ss>>num;
      fRebinF=num;
    }
  }
  fScaleK0.resize(f_K0.size());
  fScaleSm.resize(f_Sm.size());
  fScaleSp.resize(f_Sp.size());

  f_out=new TFile(outfilename.c_str(), "recreate");
  f_out-> cd();
  initHist("data");
  initHist("sim1");
  initHist("sim2");
  for( int i=0; i<f_K0.size(); i++ ) initHist(Form("K0_%d", i));
  for( int i=0; i<f_Sm.size(); i++ ) initHist(Form("Sp_%d", i));
  for( int i=0; i<f_Sp.size(); i++ ) initHist(Form("Sm_%d", i));
  for( int i=0; i<f_BG.size(); i++ ) initHist(Form("BG_%d", i));

  cout<<"===== init Hist Finish ====="<<endl;
  //  f_out->ls();

  // readData();
  // readSim1();
  // readSim2();

  fillHist(f_data, "data");
  fillHist_Sim(1);
  fillHist_Sim(2);
  for( int i=0; i<f_K0.size(); i++ ) fillHist(f_K0[i], Form("K0_%d", i));
  for( int i=0; i<f_Sm.size(); i++ ) fillHist(f_Sm[i], Form("Sp_%d", i));
  for( int i=0; i<f_Sp.size(); i++ ) fillHist(f_Sp[i], Form("Sm_%d", i));
  for( int i=0; i<f_BG.size(); i++ ) fillHist(f_BG[i], Form("BG_%d", i));

  TH1F *h1 = (TH1F*)f_out->Get("FitIM_K0_0");
  fScaleK0[0]=h1->GetEntries()/(3.0*fNumDataK0);
  fScaleSm[0]=h1->GetEntries()/(3.0*fNumDataSm);
  fScaleSp[0]=h1->GetEntries()/(3.0*fNumDataSp);
}

void FitConf::finit()
{
  cout<<"===== FitConf::finit ====="<<endl;
  if( f_out ){
    f_out-> Write();
    f_out-> Close();
    delete f_out;
  }

}

int FitConf::getMMIndex(double mm)
{
  for( int i=0; i<fBinMin.size(); i++ ){
    if( fBinMin[i]<=mm && mm<fBinMax[i] ) return i;
  }

  return -1;
}

double FitConf::getScale(const int id, double mm)
{
  if( id<1 || 2<id ){
    cout<<"  !!! FitConf::getScale error id="<<id<<"  !!!"<<endl;
  }

  for( int i=0; i<fBinMin.size(); i++ ){
    if( fBinMin[i]<=mm && mm<fBinMax[i] ){
      if( id==1 ) return fNumData[i]*fFrac1[i]/fNumSim1[i];
      if( id==2 ) return fNumData[i]*fFrac2[i]/fNumSim2[i];
    }
  }

  //  cout<<"  !!! FitConf::getScale"<<id<<"  out of range mm:"<<mm<<" !!!"<<endl;
  return 0.0;
}

