void SubtractMix(const int mode=0)
{
  //mode0 realdata
  //mode1 Sp
  //mode2 Sm
  //mode3 K0
  TH1::SetDefaultSumw2();
  const int version = 239;
  const int dEcut[3]={2,4,6};
  const int sysud[3]={-1,0,1};
  TFile *fr[4][3]={NULL};
  TFile *fmix[4][3][3]={NULL};
  TFile *fsub[4][3][3]={NULL};
  TString fnamer[4][3];
  TString fnamem[4][3][3];
  TString fnames[4][3][3];
   
  //real data
  for(int iEcut=0;iEcut<3;iEcut++){
    fnamer[0][iEcut] = Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_nostop.root",version,dEcut[iEcut]);
    fnamer[1][iEcut] = Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_qlo_nostop.root",version,dEcut[iEcut]);
    fnamer[2][iEcut] = Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_qhi_nostop.root",version,dEcut[iEcut]);
    fnamer[3][iEcut] = Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_theta15_nostop.root",version,dEcut[iEcut]);
    for(int isys=0;isys<3;isys++){
      fnamem[0][iEcut][isys] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_dE%d_iso_nostop_sys%d.root",version,dEcut[iEcut],sysud[isys]);
      fnamem[1][iEcut][isys] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_dE%d_iso_qlo_nostop_sys%d.root",version,dEcut[iEcut],sysud[isys]);
      fnamem[2][iEcut][isys] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_dE%d_iso_qhi_nostop_sys%d.root",version,dEcut[iEcut],sysud[isys]);
      fnamem[3][iEcut][isys] = Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_dE%d_iso_theta15_nostop_sys%d.root",version,dEcut[iEcut],sysud[isys]);
    }
  }

  TIter *nexthistr;
  TIter *nexthistm;
  const int nqcut = 4;

  TKey *keyr;
  TKey *keym;

  TH1 *h1r = NULL;
  TH1 *h1m = NULL;
  TH1 *h1s = NULL;
  TObject *obj = NULL;
  for(int iq=0;iq<nqcut;iq++){
    for(int iEcut=0;iEcut<3;iEcut++){
      fr[iq][iEcut] = TFile::Open(fnamer[iq][iEcut]);
      fr[iq][iEcut]->Print();
      for(int isys=0;isys<3;isys++){
        fmix[iq][iEcut][isys] = TFile::Open(fnamem[iq][iEcut][isys],"READ");
        if(!fmix[iq][iEcut][isys]) continue;
        fmix[iq][iEcut][isys]->Print();
        
        fnames[iq][iEcut][isys] = fnamer[iq][iEcut];
        fnames[iq][iEcut][isys].Replace(std::string(fnamer[iq][iEcut]).size()-5,10,Form("_sys%d_sub.root",sysud[isys]));

        fsub[iq][iEcut][isys] = new TFile(fnames[iq][iEcut][isys],"RECREATE");
        fsub[iq][iEcut][isys]->Print();

        nexthistr = new TIter(fr[iq][iEcut]->GetListOfKeys());
        nexthistm = new TIter(fmix[iq][iEcut][isys]->GetListOfKeys());

        while(  (keyr= (TKey*)nexthistr())
            && (keym= (TKey*)nexthistm())
            ){
          TClass *clr = gROOT->GetClass(keyr->GetClassName());
          TClass *clm = gROOT->GetClass(keym->GetClassName());
          if(!clr->InheritsFrom("TH1")) continue;
          if(!clm->InheritsFrom("TH1")) continue;
          h1r = (TH1*)keyr->ReadObj();
          h1m = (TH1*)keym->ReadObj();
          h1r->Add(h1m,-1.0);
          if(h1r->InheritsFrom("TH2")) {
            h1r->SetMinimum(0);
          }
          h1r->Write();
        }
        fsub[iq][iEcut][isys]->Close();
        fmix[iq][iEcut][isys]->Close();
      }
      fr[iq][iEcut]->Close();
    }
  }

}
