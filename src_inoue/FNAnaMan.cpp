#include "FNAnaMan.h"

using namespace std;

FNAnaMan::FNAnaMan() : fRebinF(1)
{
}

void FNAnaMan::postAna(){
  postAna(fFileData, "data");
  for( int i=0; i<fFileSim.size(); i++ ) postAna(fFileSim[i], Form("Sim%d", i), i);
  for( int i=0; i<fFileK0.size();  i++ ) postAna(fFileK0[i],  Form("K0_%d", i));
  for( int i=0; i<fFileSm.size();  i++ ) postAna(fFileSm[i],  Form("Sm_%d", i));
  for( int i=0; i<fFileSp.size();  i++ ) postAna(fFileSp[i],  Form("Sp_%d", i));
}

void FNAnaMan::fillHist(TDirectory *dir, const std::string &hist_name, double x, double y, int bin_index, int index){
  TH2F *h2=dynamic_cast<TH2F*>(dir->Get(hist_name.c_str()));
  if( index<0 ) h2-> Fill(x, y);
  else{
    if( bin_index<0 ) return;
    double weight=fInfoSim->val(bin_index, Form("Scale%d", index));
    h2-> Fill(x, y, weight);    
  }
}

void FNAnaMan::fillHist(TDirectory *dir, const std::string &hist_name, double x, double y){
  TH2F *h2=dynamic_cast<TH2F*>(dir->Get(hist_name.c_str()));
  h2-> Fill(x, y);    
}

void FNAnaMan::fillHist(TDirectory *dir, const std::string &hist_name, double x){
  TH1F *h1=dynamic_cast<TH1F*>(dir->Get(hist_name.c_str()));
  h1-> Fill(x);    
}

void FNAnaMan::fillHist(TDirectory *dir, const std::string &hist_name, double x, int bin_index, int index){
  TH1F *h1=dynamic_cast<TH1F*>(dir->Get(hist_name.c_str()));
  if( index<0 ) h1-> Fill(x);
  else{
    if( bin_index<0 ) return;
    double weight=fInfoSim->val(bin_index, Form("Scale%d", index));
    h1-> Fill(x, weight);    
  }
}

bool FNAnaMan::fitKNpi_all(int dump_level)
{
  int nerr=0;
  for( int i=0; i<fInfoSim->nBin(); i++ ){
    if( !fitKNpi(i, dump_level) ) nerr++;
  }
  cout<<"===== FNAnaMan::fitKNpi all  n err : "<<nerr<<endl;
  if( nerr==0 ) return true;
  else return false;
}

void FNAnaMan::fitAll(int dump_level)
{
  fitKNpi_all(dump_level);
  fitIM(dump_level);
}

bool FNAnaMan::fitIM(int dump_level)
{
  cout<<"===== FNAnaMan::fitIM Fit START ====="<<endl;
  cout<<"      dump level : "<<dump_level<<endl;
  TH1 *h1_data=(TH1*)gFile->Get("FitIM_data");
  vector<TH1*> vec_mc;
  vector<double> scale;
  vec_mc.push_back((TH1*)gFile->Get("FitIM_K0_0"));
  scale.push_back(fInfoK0->val(0,"Scale"));
  vec_mc.push_back((TH1*)gFile->Get("FitIM_Sm_0"));
  scale.push_back(fInfoSm->val(0,"Scale"));
  vec_mc.push_back((TH1*)gFile->Get("FitIM_Sp_0"));
  scale.push_back(fInfoSp->val(0,"Scale"));
  for( int i=0; i<fFileSim.size(); i++ ){
    vec_mc.push_back((TH1*)gFile->Get(Form("FitIM_Sim%d", i)));
    scale.push_back(1.0);
  }

  int nHist=vec_mc.size();
  TObjArray *mc=new TObjArray(nHist);
  for( int i=0; i<nHist; i++ ) mc->Add(vec_mc[i]);

  TemplateFitter *fitter = new TemplateFitter(h1_data, mc);
  for( int i=0; i<3; i++ ){
    if( scale[i]>0 ) fitter-> setScale(i, scale[i]);
    else fitter-> setScale(i, 1./3.);
  }
  for( int i=3; i<nHist; i++ ){
    fitter->setFixScale(i, scale[i]);
  }

  bool flag=fitter->fit(dump_level);
  if( !flag ){
    cout<<"!!!!! FNAnaMan::fitIM fault !!!!!"<<endl;
    return false;
  }

  double param, err;
  fitter->getScaleErr(0, param, err);
  fInfoK0->setVal(0, "Scale", param);
  fInfoK0->setVal(0, "Err", err);

  fitter->getScaleErr(1, param, err);
  fInfoSm->setVal(0, "Scale", param);
  fInfoSm->setVal(0, "Err", err);

  fitter->getScaleErr(2, param, err);
  fInfoSp->setVal(0, "Scale", param);
  fInfoSp->setVal(0, "Err", err);

  fLog-> AddLine(Form("FNAnaMan::fitIM  Chi2=%lf  NDF=%d", fitter->chi2(), fitter->ndf()));
  cout<<"===== FNAnaMan::fitIM Fit successfully FINISH ====="<<endl;
  return flag;
}

bool FNAnaMan::fitNpipiIM_K0(int dump_level)
{
  cout<<"===== FNAnaMan::fitNpipi_IM w/ K0 Fit START ====="<<endl;
  cout<<"      dump level : "<<dump_level<<endl;
  TH1 *h1_data=(TH1*)gFile->Get("FitIM_npipi_wN_wK0_data");
  vector<TH1*> vec_mc;
  vector<double> scale;
  for( int i=0; i<fFileK0.size(); i++ ){
    vec_mc.push_back((TH1*)gFile->Get(Form("FitIM_npipi_wN_wK0_K0_%d", i)));
    scale.push_back(fInfoK0->val(i,"Scale"));
  }
  for( int i=0; i<fFileSm.size(); i++ ){
    vec_mc.push_back((TH1*)gFile->Get(Form("FitIM_npipi_wN_wK0_Sm_%d", i)));
    scale.push_back(fInfoSm->val(i,"Scale"));
  }
  for( int i=0; i<fFileSp.size(); i++ ){
    vec_mc.push_back((TH1*)gFile->Get(Form("FitIM_npipi_wN_wK0_Sp_%d", i)));
    scale.push_back(fInfoSp->val(i,"Scale"));
  }
  for( int i=0; i<fFileSim.size(); i++ ){
    vec_mc.push_back((TH1*)gFile->Get(Form("FitIM_npipi_wN_wK0_Sim%d", i)));
    scale.push_back(1.0);
  }

  for( int i=0; i<vec_mc.size(); i++ ){
    cout<<vec_mc[i]->GetEntries()<<endl;
  }

  int nHist=vec_mc.size();
  TObjArray *mc=new TObjArray(nHist);
  for( int i=0; i<nHist; i++ ) mc->Add(vec_mc[i]);

  TemplateFitter *fitter = new TemplateFitter(h1_data, mc);
  for( int i=0; i<fFileK0.size(); i++ ){
    if( scale[i]>0 ) fitter-> setScale(i, scale[i]);
    else fitter-> setScale(i, scale[i]/fFileK0.size());
  }
  for( int i=fFileK0.size(); i<nHist; i++ ){
    fitter->setFixScale(i, scale[i]);
  }

  bool flag=fitter->fit(dump_level);
  if( !flag ){
    cout<<"!!!!! FNAnaMan::fitNpipiIM_K0 fault !!!!!"<<endl;
    return false;
  }

  for( int i=0; i<fFileK0.size(); i++ ){
    double param, err;
    fitter->getScaleErr(i, param, err);
    cout<<param<<"   "<<err<<endl;
    fInfoK0->setVal(i, "Scale", param);
    fInfoK0->setVal(i, "Err", err);
  }

  cout<<"===== FNAnaMan::fitIM Fit successfully FINISH ====="<<endl;
  return flag;
}

bool FNAnaMan::fitKNpi(int index, int dump_level)
{
  if( dump_level>0 ){
    cout<<"===== FNAnaMan::fitKNpi("<<index<<") Fit START ====="<<endl;
    cout<<"      dump level : "<<dump_level<<endl;
  }

  if( index<0 || fInfoSim->nBin()<index ) return false;
  TH1 *h1_data=(TH1*)gFile->Get(Form("FitKNpi_data_%d", index));
  vector<TH1*> vec_mc;
  vector<double> scale;
  vec_mc.push_back((TH1*)gFile->Get(Form("FitKNpi_Sim0_%d", index)));
  scale.push_back(fInfoSim->val(index, "Scale0"));
  vec_mc.push_back((TH1*)gFile->Get(Form("FitKNpi_Sim1_%d", index)));
  scale.push_back(fInfoSim->val(index, "Scale1"));
  for( int i=0; i<fFileK0.size(); i++ ){
    vec_mc.push_back((TH1*)gFile->Get(Form("FitKNpi_K0_%d_%d", i, index)));
    scale.push_back(fInfoK0->val(i, "Scale"));
  }
  for( int i=0; i<fFileSm.size(); i++ ){
    vec_mc.push_back((TH1*)gFile->Get(Form("FitKNpi_Sm_%d_%d", i, index)));
    scale.push_back(fInfoSm->val(i, "Scale"));
  }
  for( int i=0; i<fFileSp.size(); i++ ){
    vec_mc.push_back((TH1*)gFile->Get(Form("FitKNpi_Sp_%d_%d", i, index)));
    scale.push_back(fInfoSp->val(i, "Scale"));
  }

  int nHist=vec_mc.size();
  TObjArray *mc=new TObjArray(nHist);
  for( int i=0; i<nHist; i++ ) mc->Add(vec_mc[i]);

  TemplateFitter *fitter = new TemplateFitter(h1_data, mc);
  if( scale[0]>0 ) fitter-> setScale(0, scale[0]);
  else fitter-> setScale(0, 0.5);

  if( scale[1]>0 ) fitter-> setScale(1, scale[1]);
  else fitter-> setScale(1, 0.5);
  for( int i=2; i<nHist; i++ ){
    fitter->setFixScale(i, scale[i]);
  }
  bool flag=fitter->fit(dump_level);
  if( !flag ){
    cout<<"!!!!! FitKNpi("<<index<<") fit fault !!!!!"<<endl;
  }

  for( int i=0; i<2; i++ ){
    double scale, err;
    fitter->getScaleErr(i, scale, err);
    fInfoSim->setVal(index, Form("Scale%d", i), scale);
    fInfoSim->setVal(index, Form("Err%d", i), err);
  }
  fInfoSim->setVal(index, "Chi2", fitter->chi2());
  fInfoSim->setVal(index, "NDF", fitter->ndf());

  TFractionFitter *ffitter=fitter->fitter();
  fInfoSim->setVal(index, "Mat00", ffitter->GetFitter()->GetCovarianceMatrixElement(0,0));
  fInfoSim->setVal(index, "Mat01", ffitter->GetFitter()->GetCovarianceMatrixElement(0,1));
  fInfoSim->setVal(index, "Mat11", ffitter->GetFitter()->GetCovarianceMatrixElement(1,1));
  std::cout<<"(0,2)="<<ffitter->GetFitter()->GetCovarianceMatrixElement(0,2)<<std::endl;
  std::cout<<"(1,2)="<<ffitter->GetFitter()->GetCovarianceMatrixElement(1,2)<<std::endl;
  std::cout<<"(2,2)="<<ffitter->GetFitter()->GetCovarianceMatrixElement(2,2)<<std::endl;


  delete fitter;

  fLog-> AddLine(Form("FNAnaMan::fitKNpi(%d)  Chi2=%lf  NDF=%d", index, fitter->chi2(), fitter->ndf()));
  if( dump_level>0 ) cout<<"===== FNAnaMan::fitKNpi("<<index<<") scucessfully FINISH  ====="<<endl;
  return flag;
}

void FNAnaMan::fillHistMC()
{
  clearHist("MCsum");
  for( int i=0; i<fFileSim.size(); i++ ) clearHist(Form("Sim%d", i));
  for( int i=0; i<fFileK0.size(); i++ ) clearHist(Form("K0_%d", i));
  for( int i=0; i<fFileSm.size(); i++ ) clearHist(Form("Sm_%d", i));
  for( int i=0; i<fFileSp.size(); i++ ) clearHist(Form("Sp_%d", i));

  for( int i=0; i<fFileSim.size(); i++ ) fillHist(fFileSim[i], Form("Sim%d", i), i);
  for( int i=0; i<fFileK0.size(); i++ ) fillHist(fFileK0[i], Form("K0_%d", i));
  for( int i=0; i<fFileSm.size(); i++ ) fillHist(fFileSm[i], Form("Sm_%d", i));
  for( int i=0; i<fFileSp.size(); i++ ) fillHist(fFileSp[i], Form("Sp_%d", i));
}

void FNAnaMan::fillHist(TFile *f, string suffix, int index)
{
  TH1F* h1;
  TH2F* h2;
  TTree *tree = (TTree*)f->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);

  std::vector<int> num(fInfoSim->nBin(), 0.0);

  // cout<<"===== FillHist ======================================="<<endl;
  // cout<<"> File : "<<f->GetName()<<endl;
  // cout<<"> suffix : "<<suffix<<endl;
  // cout<<"> index : "<<index<<endl;
  // cout<<">   Entries : "<<tree->GetEntries()<<endl;  
  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    tree-> GetEntry(ev);
    TLorentzVector target_lmom;
    target_lmom.SetVectM(TVector3(0, 0, 0), dMass);
    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TVector3 labToCM=-(beam_lmom+target_lmom).BoostVector();
    TLorentzVector n_lmom=anaInfo->forwardNeutral(0)->lmom();
    TLorentzVector pim_lmom=anaInfo->CDS(CDS_PiMinus, 0)->lmom();
    TLorentzVector pip_lmom=anaInfo->CDS(CDS_PiPlus, 0)->lmom();
    TLorentzVector pipi_lmom=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->lmom();
    TLorentzVector pim_lmom_CM=pim_lmom, pip_lmom_CM=pip_lmom;    
    pim_lmom_CM.Boost(labToCM), pip_lmom_CM.Boost(labToCM);

    double kn_mm=(beam_lmom+target_lmom-n_lmom).M();
    double knpipi_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom-pip_lmom).M();
    double knpim_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom).M();
    double knpip_mm=(beam_lmom+target_lmom-n_lmom-pip_lmom).M();

    double npim_im=(n_lmom+pim_lmom).M();
    double npip_im=(n_lmom+pip_lmom).M();
    double pipi_im=pipi_lmom.M();
    double npipi_im=(n_lmom+pim_lmom+pip_lmom).M();

    if( Npipi_N_MIN<knpipi_mm && knpipi_mm<Npipi_N_MAX ){
      double weight=1.0;
      int bin_index=fInfoSim->index(kn_mm);
      if( index>=0 ){
	if( bin_index<0 ) weight=0.0;
	else weight=fInfoSim->val(bin_index, Form("Scale%d", index));
      }

      bool K0_flag=isK0(pipi_im); bool Sm_flag=isIM_Sm(npim_im); bool Sp_flag=isIM_Sp(npip_im);
      h1=(TH1F*)gFile->Get(Form("FitIM_npipi_wN_%s", suffix.c_str())), h1-> Fill(npipi_im);
      if( isK0(pipi_im) ) h1=(TH1F*)gFile->Get(Form("FitIM_npipi_wN_wK0_%s", suffix.c_str())), h1-> Fill(npipi_im, weight);
      if( isIM_Sm(npim_im) ) h1=(TH1F*)gFile->Get(Form("FitIM_npipi_wN_wSm_%s", suffix.c_str())), h1-> Fill(npipi_im, weight);
      if( isIM_Sp(npip_im) ) h1=(TH1F*)gFile->Get(Form("FitIM_npipi_wN_wSp_%s", suffix.c_str())), h1-> Fill(npipi_im, weight);

      h1=(TH1F*)gFile->Get(Form("FitIM_%s", suffix.c_str()));
      h1-> Fill(pipi_im, weight);
      h1-> Fill(npim_im, weight);
      h1-> Fill(1.0+npip_im, weight);
      h1=(TH1F*)gFile->Get(Form("IM_pipi_%s", suffix.c_str())), h1-> Fill(pipi_im, weight);
      h1=(TH1F*)gFile->Get(Form("IM_npim_%s", suffix.c_str())), h1-> Fill(npim_im, weight);
      h1=(TH1F*)gFile->Get(Form("IM_npip_%s", suffix.c_str())), h1-> Fill(npip_im, weight);

      if( isSignal(pipi_im, npim_im, npip_im) ){
	h1=(TH1F*)gFile->Get(Form("pim_ang_wN_woAll_CM_%s", suffix.c_str())), h1-> Fill(pim_lmom_CM.CosTheta(), weight);
	h1=(TH1F*)gFile->Get(Form("pip_ang_wN_woAll_CM_%s", suffix.c_str())), h1-> Fill(pip_lmom_CM.CosTheta(), weight);

	if( bin_index>=0 ){
	  h1=(TH1F*)gFile->Get(Form("FitKNpi_%s_%d", suffix.c_str(), bin_index));
	  h1-> Fill(knpim_mm-1.0);
	  h1-> Fill(knpip_mm);
	  
	  h1=(TH1F*)gFile->Get(Form("KNpim_MM_%s_%d", suffix.c_str(), bin_index)), h1-> Fill(knpim_mm);
	  h1=(TH1F*)gFile->Get(Form("KNpip_MM_%s_%d", suffix.c_str(), bin_index)), h1-> Fill(knpip_mm);
	  
	  h2 = (TH2F*)gFile->Get(Form("KNpim_KNpip_MM_%s_%d", suffix.c_str(), bin_index)), h2-> Fill(knpim_mm, knpip_mm);
	}
      }
    }
  }
}

void FNAnaMan::initHist(string suffix)
{
  int nbin=2000/fRebinF;
  for( int i=0;  i<fInfoSim->nBin(); i++ ){
    new TH2F(Form("KNpim_KNpip_MM_%s_%d", suffix.c_str(), i), "", nbin, 0, 2.0, nbin, 0, 2.0);
    new TH1F(Form("FitKNpi_%s_%d", suffix.c_str(), i), Form("FitKNpi_%s_%d", suffix.c_str(), i), nbin, 0, 2.0);

    new TH1F(Form("KNpim_MM_%s_%d", suffix.c_str(), i), Form("KNpim_MM_%s_%d", suffix.c_str(), i), nbin, 0, 2.0);
    new TH1F(Form("KNpip_MM_%s_%d", suffix.c_str(), i), Form("KNpip_MM_%s_%d", suffix.c_str(), i), nbin, 0, 2.0);
  }
  new TH1F(Form("FitIM_%s", suffix.c_str()), Form("FitIM_%s", suffix.c_str()), 600, 0, 3.0);

  new TH1F(Form("IM_pipi_%s", suffix.c_str()), "", 200, 0, 1.0);
  new TH1F(Form("IM_npim_%s", suffix.c_str()), "", 200, 1.0, 2.0);
  new TH1F(Form("IM_npip_%s", suffix.c_str()), "", 200, 1.0, 2.0);

  new TH1F(Form("FitIM_npipi_wN_%s",       suffix.c_str()), "", nbin, 0.0, 2.0);
  new TH1F(Form("FitIM_npipi_wN_wK0_%s",   suffix.c_str()), "", nbin, 0.0, 2.0);
  new TH1F(Form("FitIM_npipi_wN_wSm_%s",   suffix.c_str()), "", nbin, 0.0, 2.0);
  new TH1F(Form("FitIM_npipi_wN_wSp_%s",   suffix.c_str()), "", nbin, 0.0, 2.0);
  new TH1F(Form("FitIM_npipi_wN_woAll_%s", suffix.c_str()), "", nbin, 0.0, 2.0);

  new TH1F(Form("pim_ang_wN_woAll_CM_%s", suffix.c_str()), "", 1000, -1.0, 1.0);
  new TH1F(Form("pip_ang_wN_woAll_CM_%s", suffix.c_str()), "", 1000, -1.0, 1.0);
}

void FNAnaMan::clearHist(string suffix)
{
  TH1 *h1;
  for( int i=0;  i<fInfoSim->nBin(); i++ ){
    h1=(TH1*)gFile->Get(Form("KNpim_KNpip_MM_%s_%d", suffix.c_str(), i)); h1->Reset();
    h1=(TH1*)gFile->Get(Form("FitKNpi_%s_%d", suffix.c_str(), i)); h1->Reset();
    
    h1=(TH1*)gFile->Get(Form("KNpim_MM_%s_%d", suffix.c_str(), i)); h1->Reset();
    h1=(TH1*)gFile->Get(Form("KNpip_MM_%s_%d", suffix.c_str(), i)); h1->Reset();
  }
  h1=(TH1*)gFile->Get(Form("FitIM_%s", suffix.c_str())), h1->Reset();
  h1=(TH1*)gFile->Get(Form("IM_pipi_%s", suffix.c_str())), h1->Reset();
  h1=(TH1*)gFile->Get(Form("IM_npim_%s", suffix.c_str())), h1->Reset();
  h1=(TH1*)gFile->Get(Form("IM_npip_%s", suffix.c_str())), h1->Reset();

  h1=(TH1*)gFile->Get(Form("FitIM_npipi_wN_%s",       suffix.c_str())), h1->Reset();
  h1=(TH1*)gFile->Get(Form("FitIM_npipi_wN_wK0_%s",   suffix.c_str())), h1->Reset();
  h1=(TH1*)gFile->Get(Form("FitIM_npipi_wN_wSm_%s",   suffix.c_str())), h1->Reset();
  h1=(TH1*)gFile->Get(Form("FitIM_npipi_wN_wSp_%s",   suffix.c_str())), h1->Reset();
  h1=(TH1*)gFile->Get(Form("FitIM_npipi_wN_woAll_%s", suffix.c_str())), h1->Reset();

  h1=(TH1*)gFile->Get(Form("pim_ang_wN_woAll_CM_%s",   suffix.c_str())), h1->Reset();
  h1=(TH1*)gFile->Get(Form("pip_ang_wN_woAll_CM_%s", suffix.c_str())), h1->Reset();
}

void FNAnaMan::write(TString dirname)
{
  TDirectory *dir=0;
  if( dirname.Length()>0 ){
    dir=fFileOut->mkdir(dirname);
  }
  else{
    dir=fFileOut;
  }

  TIter nextKey(fFileDummy->GetListOfKeys());
  TKey *key;

  dir-> cd();
  while( (key=(TKey*)nextKey()) ){
    TObject *obj=(TObject*)fFileDummy->Get(key->GetName());
    obj->Write();
  }
  fInfoSim->Write("Info_npipin");
  fInfoK0->Write("Info_K0");
  fInfoSm->Write("Info_Sm");
  fInfoSp->Write("Info_Sp");
  TGraphErrors *gra=gra_chi2();

  gra-> Write("fitKNpi_chi2");

  fFileDummy->cd();
}

void FNAnaMan::readSim(int index)
{
  TTree *tree = (TTree*)fFileSim[index]->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);

  std::vector<int> num(fInfoSim->nBin(), 0.0);

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
      if( isSignal(pipi_im, npim_im, npip_im) ){
        int index=fInfoSim->index(kn_mm);
        if( index>=0 ) num[index]++;
      }
    }
  }

  for( int bin=0; bin<fInfoSim->nBin(); bin++ ){
    fInfoSim->setVal(bin, Form("NumSim%d",index), num[bin]);
  }
}

void FNAnaMan::readData()
{
  TTree *tree = (TTree*)fFileData->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);

  std::vector<int> num(fInfoSim->nBin(), 0.0);
  int n_All=0;
  int n_K0=0;
  int n_Sm=0;
  int n_Sp=0;

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
      n_All++;

      bool K0_flag=false; bool Sm_flag=false; bool Sp_flag=false;
      if( Npipi_K0_MIN<pipi_im && pipi_im<Npipi_K0_MAX ) K0_flag=true;
      if( Npipi_Sm_MIN<npim_im && npim_im<Npipi_Sm_MAX ) Sm_flag=true;
      if( Npipi_Sp_MIN<npip_im && npip_im<Npipi_Sp_MAX ) Sp_flag=true;

      if( isK0(pipi_im)    ) n_K0++;
      if( isIM_Sm(npim_im) ) n_Sm++;
      if( isIM_Sp(npip_im) ) n_Sp++;

      if( isSignal(pipi_im, npim_im, npip_im) ){
        int index=fInfoSim->index(kn_mm);
        if( index>=0 ) num[index]++;
      }
    }
  }
  
  TH1F *h1 = (TH1F*)gFile->Get("event_num_data");
  h1-> SetBinContent(0, n_All);
  h1-> SetBinContent(1, n_K0);
  h1-> SetBinContent(2, n_Sm);
  h1-> SetBinContent(3, n_Sp);

  for( int bin=0; bin<fInfoSim->nBin(); bin++ ){
    fInfoSim->setVal(bin, "NumData", num[bin]);
  }
}

void FNAnaMan::readAcc()
{
  for( int i=0; i<2; i++ ){
    vector<double> trig(fInfoSim->nBin(), 0.0);
    vector<double> eff(fInfoSim->nBin(), 0.0);

    Float_t mm, K0, Sm, Sp, knpim, knpip;
    TNtuple *tup=(TNtuple*)fFileSim[i]->Get("tup_trig_mass_fn");
    tup->SetBranchAddress("m", &mm);
    for( int ev=0; ev<tup->GetEntries(); ev++ ){
      tup->GetEntry(ev);
      int index=fInfoSim->index(mm);
      if( index>=0 ) trig.at(index)+=1.0;
    }

    tup=(TNtuple*)fFileSim[i]->Get("tup_eff_mass_fn_pipi_mmN_K0_Sm_Sp_knpi_mm");
    tup->SetBranchAddress("m",  &mm);
    tup->SetBranchAddress("K0", &K0);
    tup->SetBranchAddress("Sm", &Sm);
    tup->SetBranchAddress("Sp", &Sp);
    tup->SetBranchAddress("knpim_nn", &knpim);
    tup->SetBranchAddress("knpip_mm", &knpip);
    for( int ev=0; ev<tup->GetEntries(); ev++ ){
      tup->GetEntry(ev);
      int index=fInfoSim->index(mm);
      if( isSignal(K0, Sm, Sp) && index>=0 ) eff.at(index)+=1.0;
    }

    for( int bin=0; bin<fInfoSim->nBin(); bin++ ){
      fInfoSim->setVal(bin, Form("Eff%d", i), eff[bin]);
      fInfoSim->setVal(bin, Form("Trig%d", i), trig[bin]);
    }
  }
  //  fInfoSim->dump();
}

void FNAnaMan::init(TString filename)
{
  fLog=new TMacro("log", "d(K^{-} n #pi^{+} #pi^{-}) fitting log");
  fLog-> AddLine(Form("ConfFileName : %s", filename.Data()));
  printParam(fLog);
  ifstream ifs(filename);
  int num;
  string str;
  double val;
  char c_str[2000];
  string outfilename;
  vector<double> bin;
  while( getline(ifs, str) ){
    fLog-> AddLine(str.c_str());
    if( sscanf(str.c_str(), "DataFile: %s", c_str)==1 ){
      fFileData = new TFile(c_str);
    }
    else if( sscanf(str.c_str(), "OutFile: %s", c_str)==1 ){
      outfilename=c_str;
    }
    else if( sscanf(str.c_str(), "SimFile: %s", c_str)==1 ){
      fFileSim.push_back(new TFile(c_str));
    }
    else if( sscanf(str.c_str(), "K0File: %s", c_str)==1 ){
      fFileK0.push_back(new TFile(c_str));
    }
    else if( sscanf(str.c_str(), "SmFile: %s", c_str)==1 ){
      fFileSm.push_back(new TFile(c_str));
    }
    else if( sscanf(str.c_str(), "SpFile: %s", c_str)==1 ){
      fFileSp.push_back(new TFile(c_str));
    }
    else if( sscanf(str.c_str(), "Rebin: %d", &num)==1 ){
      fRebinF=num;
    }
    else if( sscanf(str.c_str(), "Bin: %lf", &val)==1 ){
      bin.push_back(val);
    }
  }
  std::cout<<"OutFile : "<<outfilename<<std::endl;
  fFileOut=new TFile(outfilename.c_str(), "recreate");
  fFileDummy=new TFile("tmp.root", "recreate");
  vector<TString> names;
  names.push_back("NumData");
  names.push_back("NumSim0");
  names.push_back("NumSim1");
  names.push_back("Trig0");
  names.push_back("Eff0");
  names.push_back("Trig1");
  names.push_back("Eff1");
  names.push_back("Scale0");
  names.push_back("Err0");
  names.push_back("Scale1");
  names.push_back("Err1");
  names.push_back("Chi2");
  names.push_back("NDF");
  names.push_back("Mat00");
  names.push_back("Mat01");
  names.push_back("Mat11");

  new TH1F("event_num_data", "event number", 10, 0, 10);

  fInfoSim=new HistInfo(bin, names);

  vector<TString> names2;
  names2.push_back("Num");
  names2.push_back("Scale");
  names2.push_back("Err");
  std::vector<double> bin2;
  bin2.push_back(0);
  for( int i=0; i<fFileK0.size(); i++ ){
    bin2.push_back(i+1);
  }
  fInfoK0=new HistInfo(bin2, names2);

  bin2.clear();
  bin2.push_back(0);
  for( int i=0; i<fFileSm.size(); i++ ){
    bin2.push_back(i+1);
  }
  fInfoSm=new HistInfo(bin2, names2);

  bin2.clear();
  bin2.push_back(0);
  for( int i=0; i<fFileSp.size(); i++ ){
    bin2.push_back(i+1);
  }
  fInfoSp=new HistInfo(bin2, names2);

  initHist("data");
  initHist("MCsum");
  for( int i=0; i<fFileSim.size(); i++ ) initHist(Form("Sim%d", i));
  for( int i=0; i<fFileK0.size(); i++ ) initHist(Form("K0_%d", i));
  for( int i=0; i<fFileSm.size(); i++ ) initHist(Form("Sm_%d", i));
  for( int i=0; i<fFileSp.size(); i++ ) initHist(Form("Sp_%d", i));

  fFileDummy->Write();

  readAcc();
  readData();
  for( int i=0; i<fFileSim.size(); i++ ){
    readSim(i);
  }

  fillHist(fFileData, "data");
  fillHistMC();

  TH1F *h1=(TH1F*)gFile->Get("event_num_data");
  double all_data=h1->GetBinContent(1);
  double K0_data=h1->GetBinContent(2);
  double Sm_data=h1->GetBinContent(3);
  double Sp_data=h1->GetBinContent(4);

  h1=(TH1F*)gFile->Get("IM_pipi_K0_0");
  double K0_sim=h1->GetEntries();

  h1=(TH1F*)gFile->Get("IM_npim_Sm_0");
  double Sm_sim=h1->GetEntries();

  h1=(TH1F*)gFile->Get("IM_npip_Sp_0");
  double Sp_sim=h1->GetEntries();

  fInfoK0->setVal(0, "Scale", K0_data/(all_data*K0_sim));
  fInfoSm->setVal(0, "Scale", Sm_data/(all_data*Sm_sim));
  fInfoSp->setVal(0, "Scale", Sp_data/(all_data*Sp_sim));
}

void FNAnaMan::finit()
{
  fFileOut-> cd();
  fLog-> Write();
  fFileOut->Close();
}

TGraphErrors *FNAnaMan::gra_chi2()
{
  double x[fInfoSim->nBin()];
  double y[fInfoSim->nBin()];
  double xerr[fInfoSim->nBin()];
  double yerr[fInfoSim->nBin()];

  // cout<<"index Chi2    /   NDF      =   Chi2/NDF"<<endl;
  // for( int i=0; i<fInfoSim->nBin(); i++ ){
  //   double chi2=fInfoSim->val(i, "Chi2");
  //   double ndf=fInfoSim->val(i, "NDF");
  //   double max=fInfoSim->bin_max(i);
  //   double min=fInfoSim->bin_max(i);

  //   x[i]=0.5*(max+min);
  //   xerr[i]=0.5*(max-min);
  //   y[i]=chi2/ndf;

  //   cout<<i<<"  "<<chi2<<"/"<<ndf<<"="<<chi2/ndf<<endl;
  // }
  TGraphErrors *gra =new TGraphErrors(fInfoSim->nBin(), x, y, xerr, y);
  return gra;
}
