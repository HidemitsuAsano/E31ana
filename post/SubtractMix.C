void SubtractMix(const int mode=0)
{
  //mode0 realdata
  //mode1 Sp
  //mode2 Sm
  //mode3 K0

  TFile *fr[4]={NULL};
  TFile *fmix[4]={NULL};
  TFile *fsub[4]={NULL};
  TString fnamer[4];
  TString fnamem[4];
  TString fnames[4];
   
  //real data
  if(mode==0){
    const int version = 206;
    /*
    fnamer[0] = Form("evanaIMpisigma_npippim_v%d_out_iso.root",version);
    fnamem[0] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_iso.root",version);
    fnamer[1] = Form("evanaIMpisigma_npippim_v%d_out_iso_qlo.root",version);
    fnamem[1] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_iso_qlo.root",version);
    fnamer[2] = Form("evanaIMpisigma_npippim_v%d_out_iso_qhi.root",version);
    fnamem[2] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_iso_qhi.root",version);
    fnamer[3] = Form("evanaIMpisigma_npippim_v%d_out_iso_theta15.root",version);
    fnamem[3] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_iso_theta15.root",version);
    */
    fnamer[0] = Form("evanaIMpisigma_npippim_v%d_out_iso_nostop.root",version);
    fnamem[0] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_iso_nostop.root",version);
    fnamer[1] = Form("evanaIMpisigma_npippim_v%d_out_iso_qlo_nostop.root",version);
    fnamem[1] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_iso_qlo_nostop.root",version);
    fnamer[2] = Form("evanaIMpisigma_npippim_v%d_out_iso_qhi_nostop.root",version);
    fnamem[2] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_iso_qhi_nostop.root",version);
    fnamer[3] = Form("evanaIMpisigma_npippim_v%d_out_iso_theta15_nostop.root",version);
    fnamem[3] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_iso_theta15_nostop.root",version);
  }else if(mode==1){
    const int version = 132;
    fnamer[0] = Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_out_iso.root",version);
    fnamem[0] = Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_MIX_cut4_out_iso.root",version);
    fnamer[1] = Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_out_iso_qlo.root",version);
    fnamem[1] = Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_MIX_cut4_out_iso_qlo.root",version);
    fnamer[2] = Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_out_iso_qhi.root",version);
    fnamem[2] = Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_MIX_cut4_out_iso_qhi.root",version);
    fnamer[3] = Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_out_iso_theta15.root",version);
    fnamem[3] = Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_MIX_cut4_out_iso_theta15.root",version);
  }else if(mode==2){
    const int version = 132;
    fnamer[0] = Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_out_iso.root",version);
    fnamem[0] = Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_MIX_cut4_out_iso.root",version);
    fnamer[1] = Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_out_iso_qlo.root",version);
    fnamem[1] = Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_MIX_cut4_out_iso_qlo.root",version);
    fnamer[2] = Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_out_iso_qhi.root",version);
    fnamem[2] = Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_MIX_cut4_out_iso_qhi.root",version);
    fnamer[3] = Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_out_iso_theta15.root",version);
    fnamem[3] = Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_MIX_cut4_out_iso_theta15.root",version);
  }else if(mode==3){
    const int version = 11;
    fnamer[0] = Form("../simpost/simIMpisigma_K0nn_pippimn_v%d_out_iso.root",version);
    fnamem[0] = Form("../simpost/simIMpisigma_K0nn_pippimn_v%d_MIX_cut4_out_iso.root",version);
    fnamer[1] = Form("../simpost/simIMpisigma_K0nn_pippimn_v%d_out_iso_qlo.root",version);
    fnamem[1] = Form("../simpost/simIMpisigma_K0nn_pippimn_v%d_MIX_cut4_out_iso_qlo.root",version);
    fnamer[2] = Form("../simpost/simIMpisigma_K0nn_pippimn_v%d_out_iso_qhi.root",version);
    fnamem[2] = Form("../simpost/simIMpisigma_K0nn_pippimn_v%d_MIX_cut4_out_iso_qhi.root",version);
    fnamer[3] = Form("../simpost/simIMpisigma_K0nn_pippimn_v%d_out_iso_theta15.root",version);
    fnamem[3] = Form("../simpost/simIMpisigma_K0nn_pippimn_v%d_MIX_cut4_out_iso_theta15.root",version);
  }
  
  TIter *nexthistr[4];
  TIter *nexthistm[4];
  const int nqcut = 4;

  TKey *keyr[4];
  TKey *keym[4];

  TH1 *h1r = NULL;
  TH1 *h1m = NULL;
  TH1 *h1s = NULL;
  TObject *obj = NULL;
  for(int iq=0;iq<nqcut;iq++){
    fr[iq] = TFile::Open(fnamer[iq]);
    fmix[iq] = TFile::Open(fnamem[iq]);
    fr[iq]->Print();
    fmix[iq]->Print();
    fnames[iq] = fnamer[iq];
    fnames[iq].Replace(std::string(fnamer[iq]).size()-5,5,"_sub.root");

    fsub[iq] = new TFile(fnames[iq],"RECREATE");
    fsub[iq]->Print();

    nexthistr[iq] = new TIter(fr[iq]->GetListOfKeys());
    nexthistm[iq] = new TIter(fmix[iq]->GetListOfKeys());

    while(  (keyr[iq]= (TKey*)nexthistr[iq]())
          && (keym[iq]= (TKey*)nexthistm[iq]())
        ){
      TClass *clr = gROOT->GetClass(keyr[iq]->GetClassName());
      TClass *clm = gROOT->GetClass(keym[iq]->GetClassName());
      if(!clr->InheritsFrom("TH1")) continue;
      if(!clm->InheritsFrom("TH1")) continue;
      h1r = (TH1*)keyr[iq]->ReadObj();
      h1m = (TH1*)keym[iq]->ReadObj();
      h1r->Add(h1m,-1.0);
      if(h1r->InheritsFrom("TH2")) {
        h1r->SetMinimum(0);
      }
      h1r->Write();
    }
    fsub[iq]->Close();
    fr[iq]->Close();
    fmix[iq]->Close();
  }


}
