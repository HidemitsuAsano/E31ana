void SubtractMix_H2(const int dEcut=2,const int sysud=0)
{
  //mode0 realdata
  //mode1 Sp
  //mode2 Sm
  //mode3 K0
  const int version = 14;

  TFile *fr;
  TFile *fmix;
  TFile *fsub;
  TString fnamer;
  TString fnamem;
  TString fnames;
   
  //real data
  fnamer = Form("evanaIMsigma_npi_h2_v%d_out_dE%d_iso_nostop.root",version,dEcut);
  fnamem = Form("evanaIMsigma_npi_h2_v%d_MIX_out_dE%d_iso_nostop_sys%d.root",version,dEcut,sysud);
  
  TIter *nexthistr;
  TIter *nexthistm;

  TKey *keyr;
  TKey *keym;

  TH1 *h1r = NULL;
  TH1 *h1m = NULL;
  TH1 *h1s = NULL;
  TObject *obj = NULL;
  fr = TFile::Open(fnamer);
  fmix = TFile::Open(fnamem);
  fr->Print();
  fmix->Print();
  fnames = fnamer;
  fnames.Replace(std::string(fnamer).size()-5,5,Form("_sub_sys%d.root",sysud));

  fsub = new TFile(fnames,"RECREATE");
  fsub->Print();

  nexthistr = new TIter(fr->GetListOfKeys());
  nexthistm = new TIter(fmix->GetListOfKeys());

  while(  (keyr= (TKey*)(*nexthistr)())
      && (keym= (TKey*)(*nexthistm)())
      ){
    TClass *clr = gROOT->GetClass(keyr->GetClassName());
    TClass *clm = gROOT->GetClass(keym->GetClassName());
    if(!clr->InheritsFrom("TH1")) continue;
    if(!clm->InheritsFrom("TH1")) continue;
    h1r = (TH1*)keyr->ReadObj();
    h1m = (TH1*)keym->ReadObj();
    h1r->Add(h1m,-1.0);
    if(h1r->InheritsFrom("TH2")) {
      //h1r->SetMinimum(0);
    }
    h1r->Write();
  }
  fsub->Close();
  fr->Close();
  fmix->Close();

}
