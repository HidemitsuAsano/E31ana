#include "HistManwMC.h"

void HistManwMC::initCalib()
{
  //*****************//
  //*** Beam Line ***//
  //*****************//
  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("T0NC_offset_T0%d", seg), Form("T0-NC offset T0 seg%d", seg), 1000, -5, 5);

    new TH2F(Form("T0_eu_ctm_%d_NC", seg), Form("T0 seg%d dE_{up} vs offset by NC", seg),   500, 0, 20, 500,  -5, 5);
    new TH2F(Form("T0_ed_ctm_%d_NC", seg), Form("T0 seg%d dE_{down} vs offset by NC", seg), 500, 0, 20, 500,  -5, 5);
  }

  new TH1F("nBPD", "BPD Multiplicity", 71, -0.5, 70.5);
  new TH1F("hitpatBPD", "BPD hit pattern", 70, 0.5, 70.5);
  new TH1F("hitpatBPD_calib", "BPD hit pattern", 70, 0.5, 70.5);
  new TH1F("hitpatBPD_calib_pi", "BPD hit pattern", 70, 0.5, 70.5);
  new TH1F("hitpatBPD_calib_k", "BPD hit pattern", 70, 0.5, 70.5);
  for( int seg=1; seg<=70; seg++ ){
    new TH1F(Form("BPD_UT_%d", seg), Form("BPD TDC_{up} seg%d", seg),   4000, 0, 4000);
    new TH1F(Form("BPD_DT_%d", seg), Form("BPD TDC_{down} seg%d", seg), 4000, 0, 4000);
    new TH1F(Form("BPD_UA_%d", seg), Form("BPD ADC_{up} seg%d", seg),   4000, 0, 4000);
    new TH1F(Form("BPD_DA_%d", seg), Form("BPD ADC_{down} seg%d", seg), 4000, 0, 4000);
    new TH1F(Form("BPD_UAwT_%d", seg), Form("BPD ADC_{up} seg%d", seg),   4000, 0, 4000);
    new TH1F(Form("BPD_DAwT_%d", seg), Form("BPD ADC_{down} seg%d", seg), 4000, 0, 4000);
    new TH1F(Form("BPD_UAwT2_%d", seg), Form("BPD ADC_{up} seg%d", seg),   4000, 0, 4000);
    new TH1F(Form("BPD_DAwT2_%d", seg), Form("BPD ADC_{down} seg%d", seg), 4000, 0, 4000);

    new TH1F(Form("T0BPD_offset_BPD%d", seg), Form("T0-BPD offset BPD seg%d", seg), 1000, -50, 50);

    new TH2F(Form("BPD_eu_offset_%d", seg), Form("BPD seg%d dE_{up} vs offset by NC", seg),   500, 0, 20, 500,  -10, 10);
    new TH2F(Form("BPD_ed_offset_%d", seg), Form("BPD seg%d dE_{down} vs offset by NC", seg), 500, 0, 20, 500,  -10, 10);
  }

  new TH1F("BHDT0_TOF_pi", "T0-BHD seg%d TOF #pi", 1000, 0.0, 5.0); 
  new TH1F("BHDT0_TOF_k",  "T0-BHD seg%d TOF K", 1000, 0.0, 5.0); 
  for( int seg=1; seg<=20; seg++ ){
    new TH1F(Form("BHDT0_offset_pi_BHD%d", seg), Form("T0-BHD seg%d TOF #pi", seg), 1000, -20.0, 20.0); 
    new TH1F(Form("BHDT0_offset_k_BHD%d", seg),  Form("T0-BHD seg%d TOF K", seg), 1000, -20.0, 20.0); 

    new TH2F(Form("BHD_eu_offset_%d", seg), Form("BHD seg%d dE_{up} vs offset", seg),   500, 0, 20, 500, -5, 5);
    new TH2F(Form("BHD_ed_offset_%d", seg), Form("BHD seg%d dE_{down} vs offset", seg), 500, 0, 20, 500, -5, 5);
  }

  for( int lay=1; lay<=8; lay++ ){
    new TH1F(Form("BLC1a_res_%d", lay), Form("BLC1a residual layer%d", lay), 1000, -1, 1);
    new TH1F(Form("BLC1b_res_%d", lay), Form("BLC1b residual layer%d", lay), 1000, -1, 1);
    new TH1F(Form("BLC2a_res_%d", lay), Form("BLC2a residual layer%d", lay), 1000, -1, 1);
    new TH1F(Form("BLC2b_res_%d", lay), Form("BLC2b residual layer%d", lay), 1000, -1, 1);
    new TH1F(Form("BPC_res_%d",   lay), Form("BPC residual layer%d",   lay), 1000, -1, 1);
    new TH1F(Form("FDC1_res_%d",  lay), Form("FDC1 residual layer%d",  lay), 1000, -1, 1);
  }

  //***************//
  //*** Forward ***//
  //***************//
  new TH2F("NC_eu_ctm", "NC seg%d dE_{up} vs offset",   500, 0, 100, 500, -5, 5);
  new TH2F("NC_ed_ctm", "NC seg%d dE_{down} vs offset", 500, 0, 100, 500, -5, 5);

  for( int seg=1; seg<=112; seg++ ){
    new TH1F(Form("T0NC_offset_NC%d_All", seg), Form("T0-NC offset NC seg%d", seg), 1000, -5, 5);
    new TH1F(Form("T0NC_offset_NC%d", seg), Form("T0-NC offset NC seg%d", seg), 1000, -5, 5);

    new TH2F(Form("NC_eu_ctm_%d", seg), Form("NC seg%d dE_{up} vs offset", seg),   500, 0, 100, 500, -5, 5);
    new TH2F(Form("NC_ed_ctm_%d", seg), Form("NC seg%d dE_{down} vs offset", seg), 500, 0, 100, 500, -5, 5);

    new TH1F(Form("NC_eu_%d", seg), Form("NC seg%d dE_{up}", seg),   1000, -10.0, 90);
    new TH1F(Form("NC_ed_%d", seg), Form("NC seg%d dE_{down}", seg), 1000, -10.0, 90);

    new TH1F(Form("NC_eu_wt_%d", seg), Form("NC seg%d dE_{up} w/ TDC", seg),   1000, -10.0, 90);
    new TH1F(Form("NC_ed_wt_%d", seg), Form("NC seg%d dE_{down} w/ TDC", seg), 1000, -10.0, 90);

    new TH1F(Form("NC_tsub_%d", seg), Form("NC tsub seg%d", seg), 1000, -50, 50);
    new TH1F(Form("NC_hitpos_%d", seg), Form("NC hitpos seg%d", seg), 1000, -150, 150);
  }

  //***********//
  //*** CDS ***//
  //***********//
  for( int seg=1; seg<=36; seg++ ){
    new TH1F(Form("CDH_offset_%d", seg), Form("CDH seg%d offset", seg), 1000, -10, 10);

    new TH1F(Form("CDH_eu_%d", seg), Form("CDH seg%d dE_{up}", seg),   1000, -10.0, 90);
    new TH1F(Form("CDH_ed_%d", seg), Form("CDH seg%d dE_{down}", seg), 1000, -10.0, 90);

    new TH1F(Form("CDH_eu_wt_%d", seg), Form("CDH seg%d dE_{up}", seg),   1000, -10.0, 90);
    new TH1F(Form("CDH_ed_wt_%d", seg), Form("CDH seg%d dE_{down}", seg), 1000, -10.0, 90);

    new TH1F(Form("CDH_eu_wpi_%d", seg), Form("CDH seg%d dE_{up}", seg),   1000, -10.0, 90);
    new TH1F(Form("CDH_ed_wpi_%d", seg), Form("CDH seg%d dE_{down}", seg), 1000, -10.0, 90);
  }
}

void HistManwMC::fillCalib(EventHeader* header)
{
  TH1F *h1;
  TH2F *h2;

  int nBPD=0;
  for( int i=0; i<blMan->nBPD(); i++ ){
    HodoscopeLikeHit *hit = blMan->BPD(i);
    int seg = hit-> seg();
    int adcu = hit->adcu(), adcd = hit->adcd();
    int tdcu = hit->tdcu(), tdcd = hit->tdcd();
    h1 = (TH1F*)rtFile-> Get(Form("BPD_UA_%d", seg)), h1-> Fill(adcu);
    h1 = (TH1F*)rtFile-> Get(Form("BPD_DA_%d", seg)), h1-> Fill(adcd);
    if( hit->CheckRange() ){
      nBPD++;
      h1 = (TH1F*)rtFile-> Get("hitpatBPD"), h1-> Fill(seg);
      h1 = (TH1F*)rtFile-> Get(Form("BPD_UAwT_%d", seg)), h1-> Fill(adcu);
      h1 = (TH1F*)rtFile-> Get(Form("BPD_DAwT_%d", seg)), h1-> Fill(adcd);

      h1 = (TH1F*)rtFile-> Get(Form("BPD_UT_%d", seg)), h1-> Fill(tdcu);
      h1 = (TH1F*)rtFile-> Get(Form("BPD_DT_%d", seg)), h1-> Fill(tdcd);
    }
  }
  h1 = (TH1F*)rtFile->Get("nBPD"), h1-> Fill(nBPD);

  if( fT0_hit.size()==1 && fBeamPID==Beam_Kaon ){
    for( int i=0; i<bltrackMan->ntrackBLC1a(); i++ ){
      LocalTrack *track = bltrackMan->trackBLC1a(i);
      for( int j=0; j<track->nhit(); j++ ){
	ChamberLikeHit *hit = track->hit(j);
	int lay = hit-> layer();
	double res = hit->resl();
	h1 = (TH1F*)rtFile->Get(Form("BLC1a_res_%d", lay)), h1-> Fill(res);
      }
    }
    for( int i=0; i<bltrackMan->ntrackBLC1b(); i++ ){
      LocalTrack *track = bltrackMan->trackBLC1b(i);
      for( int j=0; j<track->nhit(); j++ ){
	ChamberLikeHit *hit = track->hit(j);
	int lay = hit-> layer();
	double res = hit->resl();
	h1 = (TH1F*)rtFile->Get(Form("BLC1b_res_%d", lay)), h1-> Fill(res);
      }
    }
    for( int i=0; i<bltrackMan->ntrackBLC2a(); i++ ){
      LocalTrack *track = bltrackMan->trackBLC2a(i);
      for( int j=0; j<track->nhit(); j++ ){
	ChamberLikeHit *hit = track->hit(j);
	int lay = hit-> layer();
	double res = hit->resl();
	h1 = (TH1F*)rtFile->Get(Form("BLC2a_res_%d", lay)), h1-> Fill(res);
      }
    }
    for( int i=0; i<bltrackMan->ntrackBLC2b(); i++ ){
      LocalTrack *track = bltrackMan->trackBLC2b(i);
      for( int j=0; j<track->nhit(); j++ ){
	ChamberLikeHit *hit = track->hit(j);
	int lay = hit-> layer();
	double res = hit->resl();
	h1 = (TH1F*)rtFile->Get(Form("BLC2b_res_%d", lay)), h1-> Fill(res);
      }
    }
    for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
      LocalTrack *track = bltrackMan->trackBPC(i);
      for( int j=0; j<track->nhit(); j++ ){
	ChamberLikeHit *hit = track->hit(j);
	int lay = hit-> layer();
	double res = hit->resl();
	h1 = (TH1F*)rtFile->Get(Form("BPC_res_%d", lay)), h1-> Fill(res);
      }
    }
    for( int i=0; i<bltrackMan->ntrackFDC1(); i++ ){
      LocalTrack *track = bltrackMan->trackFDC1(i);
      for( int j=0; j<track->nhit(); j++ ){
	ChamberLikeHit *hit = track->hit(j);
	int lay = hit-> layer();
	double res = hit->resl();
	h1 = (TH1F*)rtFile->Get(Form("FDC1_res_%d", lay)), h1-> Fill(res);
      }
    }
  }

  if( fT0_hit.size()==1 && fBHD_hit.size()==1 ){
    double T0time = fT0_hit[0]-> ctmean();
    double BHDtime = fBHD_hit[0]-> ctmean();
    double tof = BHDtime-T0time;

    double calc_beta = 1./sqrt(1.+fD5mom*fD5mom/(parMass[fBeamPID]*parMass[fBeamPID]));
    double offset = tof-770./(100.*calc_beta*Const);

    int BHDseg = fBHD_hit[0]-> seg();
    if( fBeamPID==Beam_Pion ){
      h1 = (TH1F*)rtFile-> Get(Form("BHDT0_offset_pi_BHD%d", BHDseg)), h1-> Fill(offset);
    }
    if( fBeamPID==Beam_Kaon ){
      h1 = (TH1F*)rtFile-> Get(Form("BHDT0_offset_k_BHD%d", BHDseg)), h1-> Fill(offset);
      h2 = (TH2F*)rtFile-> Get(Form("BHD_eu_offset_%d", BHDseg)), h2-> Fill(offset, fBHD_hit[0]->eu());
      h2 = (TH2F*)rtFile-> Get(Form("BHD_ed_offset_%d", BHDseg)), h2-> Fill(offset, fBHD_hit[0]->ed());
    }
  }

  if( fT0_hit.size()==1 && fBPD_hit.size()==1 && bltrackMan->ntrackBLC2()==1 ){
    double T0time = fT0_hit[0]-> ctmean();
    double BPDtime = fBPD_hit[0]-> ctmean();
    double tof = BPDtime-T0time;
    int BPDseg = fBPD_hit[0]-> seg();
    TVector3 BPDpos;
    confMan-> GetGeomMapManager()-> GetGPos(CID_BPD, BPDseg, BPDpos);
    BPDpos = bltrackMan-> trackBLC2(0)-> GetPosatZ(BPDpos.Z()-0.25);
    double fl = (BPDpos-fT0pos).Mag();

    double calc_beta = 1./sqrt(1.+fD5mom*fD5mom/(parMass[fBeamPID]*parMass[fBeamPID]));
    double offset = tof - fl/(100.*Const*calc_beta);

    h1 = (TH1F*)rtFile-> Get("hitpatBPD_calib"), h1-> Fill(BPDseg);
    if( fBeamPID==Beam_Pion ) h1 = (TH1F*)rtFile-> Get("hitpatBPD_calib_pi"), h1-> Fill(BPDseg);
    if( fBeamPID==Beam_Kaon ){
      h1 = (TH1F*)rtFile-> Get("hitpatBPD_calib_k"), h1-> Fill(BPDseg);
      h1 = (TH1F*)rtFile-> Get(Form("T0BPD_offset_BPD%d", BPDseg)), h1-> Fill(offset);

      int adcu = fBPD_hit[0]-> adcu(), adcd = fBPD_hit[0]-> adcd(); 
      double eu = fBPD_hit[0]->eu(), ed = fBPD_hit[0]->ed(); 
      h1 = (TH1F*)rtFile-> Get(Form("BPD_UAwT2_%d", BPDseg)), h1-> Fill(adcu);
      h1 = (TH1F*)rtFile-> Get(Form("BPD_DAwT2_%d", BPDseg)), h1-> Fill(adcd);

      h2 = (TH2F*)rtFile-> Get(Form("BPD_eu_offset_%d", BPDseg)), h2-> Fill(eu, offset);
      h2 = (TH2F*)rtFile-> Get(Form("BPD_ed_offset_%d", BPDseg)), h2-> Fill(ed, offset);
    }
  }

  if( header ){
    if( header-> IsTrig(Trig_Neutral) ){
      for( int i=0; i<blMan->nNC(); i++ ){
	HodoscopeLikeHit *hit = blMan->NC(i);
	int seg = hit-> seg();
	double eu = hit->eu(), ed = hit->ed();
	h1 = (TH1F*)rtFile-> Get(Form("NC_eu_%d", seg)), h1-> Fill(eu);
	h1 = (TH1F*)rtFile-> Get(Form("NC_ed_%d", seg)), h1-> Fill(ed);
	if( hit-> CheckRange() ){
	  h1 = (TH1F*)rtFile-> Get(Form("NC_eu_wt_%d", seg)), h1-> Fill(eu);
	  h1 = (TH1F*)rtFile-> Get(Form("NC_ed_wt_%d", seg)), h1-> Fill(ed);
	}
      }
    }
    if( header->IsTrig(Trig_KCDH1) ){
      for( int i=0; i<cdsMan->nCDH(); i++ ){
	HodoscopeLikeHit *hit = cdsMan->CDH(i);
	int seg = hit->seg();
	double eu = hit->eu(), ed = hit->ed();
	h1 = (TH1F*)rtFile-> Get(Form("CDH_eu_%d", seg)), h1-> Fill(eu);
	h1 = (TH1F*)rtFile-> Get(Form("CDH_ed_%d", seg)), h1-> Fill(ed);
	if( hit->CheckRange() ){
	  h1 = (TH1F*)rtFile-> Get(Form("CDH_eu_wt_%d", seg)), h1-> Fill(eu);
	  h1 = (TH1F*)rtFile-> Get(Form("CDH_ed_wt_%d", seg)), h1-> Fill(ed);
	}
      }
    }
  }

  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *track = cdstrackMan->GoodTrack(i);
    for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
      for( int i=0; i<track->nTrackHit(lay); i++ ){
	CDCHit *hit = track->hit(cdsMan, lay, i);
	h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d", lay)), h1-> Fill(hit->resl());
	if( track->PID()==CDS_PiMinus ) h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d_pim", lay)), h1-> Fill(hit->resl());
	if( track->PID()==CDS_PiPlus  ) h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d_pip", lay)), h1-> Fill(hit->resl());
	if( track->PID()==CDS_Kaon    ) h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d_km", lay)), h1-> Fill(hit->resl());
	if( track->PID()==CDS_Proton  ) h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d_p", lay)), h1-> Fill(hit->resl());
      }

      if( track->nCDHHit()==1 ){
	if( track->PID()==CDS_PiMinus || track->PID()==CDS_PiPlus ){
	  if( track->Momentum()>0.3 ){
	    HodoscopeLikeHit *hit = track->CDHHit(cdsMan, 0);
	    int seg=hit->seg();
	    double eu=hit->eu(), ed=hit->ed();
	    h1 = (TH1F*)rtFile-> Get(Form("CDH_eu_wpi_%d", seg)), h1-> Fill(eu);
	    h1 = (TH1F*)rtFile-> Get(Form("CDH_ed_wpi_%d", seg)), h1-> Fill(ed);
	  }
	}
      }
    }
  }
}
