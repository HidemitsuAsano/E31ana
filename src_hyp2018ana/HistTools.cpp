#include "HistTools.h"

void fillCDC(CDSTrack *track, const double &beta, CDSHitMan* cdsMan, ConfMan *conf)
{
  TH1F *h1;
  TH2F *h2;

  bool SINGLE=true;
  for( int l=1; l<=NumOfCDCLayers; l++ ){
    if( track->nTrackHit(l)!=1 ) SINGLE=false;
    for( int i=0; i<track->nTrackHit(l); i++ ){
      CDCHit * cdc=track->TrackHit(cdsMan, l, i);
      int wire=cdc->wire();

      h1 = (TH1F*)gFile-> Get(Form("hitpatCDC_%d_track", l)), h1-> Fill(wire);
    }
  }
  double param[5];
  track->GetParameters(param);

  if( SINGLE && fabs(param[3])<12 && fabs(param[0])<6 ){
    for( int l=1; l<=NumOfCDCLayers; l++ ){
      for( int i=0; i<track->nTrackHit(l); i++ ){
	CDCHit * cdc=track->TrackHit(cdsMan, l, i);
	int wire=cdc->wire();
	double res = cdc->resl();
	double dt=cdc->dt();
	double dl=cdc->dl();
	double dlr=dl-res;
	double dxdt=conf-> GetXTMapManager()-> CalcDxDt(CID_CDC, l, 0, dt);

	// h1 = (TH1F*)gFile-> Get(Form("CDC_dt_l%d_w%d", l, wire)), h1-> Fill(dt);
	// h1 = (TH1F*)gFile-> Get(Form("CDC_res_l%d_w%d", l, wire)), h1-> Fill(res);
	// h2 = (TH2F*)gFile-> Get(Form("CDC_dt_dl_%d", l)), h2-> Fill(dt, dlr);
	// h2 = (TH2F*)gFile-> Get(Form("CDC_dt_res_%d", l)), h2-> Fill(dt, res);
	// h2 = (TH2F*)gFile-> Get("CDC_ob2_res"), h2-> Fill(1/(beta*beta), res/dxdt);
	// h2 = (TH2F*)gFile-> Get(Form("CDC_ob2_res_%d", l)), h2-> Fill(1/(beta*beta), res/dxdt);
      }
    }
  }
}

void fillBLC1a(LocalTrack *track)
{
  TH1F *h1;
  TH2F *h2;

  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int lay  = hit->layer();
    int wire = hit->wire();
    double dt = hit-> dt();
    double res = hit-> resl();
    double dl=hit-> dl();
    if( hit->leftright()==0 ) dl *= -1.0;

    //    h2 = (TH2F*)gFile-> Get(Form("BLC1a_l%d_dx_res", lay)), h2-> Fill(dl, res);
    h1 = (TH1F*)gFile-> Get(Form("BLC1a_l%d_res", lay)), h1-> Fill(res);
    h1 = (TH1F*)gFile-> Get(Form("BLC1a_l%d_dt", lay)), h1-> Fill(dt);
    h1 = (TH1F*)gFile-> Get(Form("BLC1a_l%d_w%d_dt", lay, wire)), h1-> Fill(dt);
  }
}

void fillBLC1b(LocalTrack *track)
{
  TH1F *h1;
  TH2F *h2;

  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int lay  = hit->layer();
    int wire = hit->wire();
    double dt = hit-> dt();
    double res = hit-> resl();
    double dl=hit-> dl();
    if( hit->leftright()==0 ) dl *= -1.0;

    //    h2 = (TH2F*)gFile-> Get(Form("BLC1b_l%d_dx_res", lay)), h2-> Fill(dl, res);
    h1 = (TH1F*)gFile-> Get(Form("BLC1b_l%d_res", lay)), h1-> Fill(res);
    h1 = (TH1F*)gFile-> Get(Form("BLC1b_l%d_dt", lay)), h1-> Fill(dt);
    h1 = (TH1F*)gFile-> Get(Form("BLC1b_l%d_w%d_dt", lay, wire)), h1-> Fill(dt);
  }
}

void fillBLC2a(LocalTrack *track)
{
  TH1F *h1;
  TH2F *h2;

  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int lay  = hit->layer();
    int wire = hit->wire();
    double dt = hit-> dt();
    double res = hit-> resl();
    double dl=hit-> dl();
    if( hit->leftright()==0 ) dl *= -1;

    //    h2 = (TH2F*)gFile-> Get(Form("BLC2a_l%d_dx_res", lay)), h2-> Fill(dl, res);
    h1 = (TH1F*)gFile-> Get(Form("BLC2a_l%d_res", lay)), h1-> Fill(res);
    h1 = (TH1F*)gFile-> Get(Form("BLC2a_l%d_dt", lay)), h1-> Fill(dt);
    h1 = (TH1F*)gFile-> Get(Form("BLC2a_l%d_w%d_dt", lay, wire)), h1-> Fill(dt);
  }
}

void fillBLC2b(LocalTrack *track)
{
  TH1F *h1;
  TH2F *h2;

  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int lay  = hit->layer();
    int wire = hit->wire();
    double dt = hit-> dt();
    double res = hit-> resl();
    double dl=hit-> dl();
    if( hit->leftright()==0 ) dl *= -1;

    //    h2 = (TH2F*)gFile-> Get(Form("BLC2b_l%d_dx_res", lay)), h2-> Fill(dl, res);
    h1 = (TH1F*)gFile-> Get(Form("BLC2b_l%d_res", lay)), h1-> Fill(res);
    h1 = (TH1F*)gFile-> Get(Form("BLC2b_l%d_dt", lay)), h1-> Fill(dt);
    h1 = (TH1F*)gFile-> Get(Form("BLC2b_l%d_w%d_dt", lay, wire)), h1-> Fill(dt);
  }
}



void fillBPC(LocalTrack *track)
{
  TH1F *h1;
  TH2F *h2;

  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int cid  = hit-> cid();
    int lay  = hit->layer();
    int wire = hit->wire();
    double dt = hit-> dt();
    double res = hit-> resl();
    double dl = hit-> dl();
    if( hit->leftright()==0 ) dl *= -1;

    //    h2 = (TH2F*)gFile-> Get(Form("BPC_l%d_dx_res", lay)), h2-> Fill(dl, res);
    h1 = (TH1F*)gFile-> Get(Form("BPC_l%d_dt", lay)), h1-> Fill(dt);
    h1 = (TH1F*)gFile-> Get(Form("BPC_l%d_res", lay)), h1-> Fill(res);
    h1 = (TH1F*)gFile-> Get(Form("BPC_l%d_w%d_dt", lay, wire)), h1-> Fill(dt);
  }
}

void fillFDC1(LocalTrack *track)
{
  TH1F *h1;
  TH2F *h2;

  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int cid  = hit-> cid();
    int lay  = hit->layer();
    int wire = hit->wire();
    double dt = hit-> dt();
    double res = hit-> resl();
    double dl=hit-> dl();
    if( hit->leftright()==0 ) dl *= -1;

    //    h2 = (TH2F*)gFile-> Get(Form("FDC1_l%d_dx_res", lay)), h2-> Fill(dl, res);
    h1 = (TH1F*)gFile-> Get(Form("FDC1_l%d_dt", lay)), h1-> Fill(dt);
    h1 = (TH1F*)gFile-> Get(Form("FDC1_l%d_res", lay)), h1-> Fill(res);
    h1 = (TH1F*)gFile-> Get(Form("FDC1_l%d_w%d_dt", lay, wire)), h1-> Fill(dt);
  }
}

void initHistCDH()
{
  new TH1F("CDH_trig", "CDH trigger",   36, 0.5, 36.5);
  new TH1F("CDH_eff",  "CDH effective", 36, 0.5, 36.5);

  new TH2F("CDH_offset_mom",     "CDH offset vs mom",         400, -5, 5, 400, 0, 1.0);
  new TH2F("CDH_offset_mom_pim", "CDH offset vs mom #pi^{-}", 400, -5, 5, 400, 0, 1.0);
  new TH2F("CDH_offset_mom_pip", "CDH offset vs mom #pi^{+}", 400, -5, 5, 400, 0, 1.0);
  new TH2F("CDH_offset_mom_km",  "CDH offset vs mom K^{-}",   400, -5, 5, 400, 0, 1.0);
  new TH2F("CDH_offset_mom_p",   "CDH offset vs mom p",       400, -5, 5, 400, 0, 1.0);

  for( int seg=1; seg<=36; seg++ ){
    new TH1F(Form("CDH_time_%d", seg),    Form("CDH seg%d time", seg), 1050, -5, 100);
    new TH1F(Form("CDH_dE_%d", seg),      Form("CDH seg%d dE", seg), 1050, -5, 100);
    new TH1F(Form("CDH_dE_woT_%d", seg),  Form("CDH seg%d dE", seg), 1050, -5, 100);
    new TH1F(Form("CDH_dE_wCDC_%d", seg), Form("CDH seg%d dE", seg), 1050, -5, 100);

    // new TNtuple(Form("CDH%d_slewing_info",seg), Form("CDH%d_slewing_infp", seg), "eu:ed:tu:td:ctu:ctd:ctm");
    new TH1F(Form("CDH_offset0_%d", seg),     Form("CDH seg%d offset", seg), 1000, -10, 10);

    new TH1F(Form("CDH_offset_%d", seg),     Form("CDH seg%d offset", seg), 1000, -10, 10);
    new TH1F(Form("CDH_offset_%d_pim", seg), Form("CDH seg%d offset", seg), 1000, -10, 10);
    new TH1F(Form("CDH_offset_%d_pip", seg), Form("CDH seg%d offset", seg), 1000, -10, 10);
    new TH1F(Form("CDH_offset_%d_km", seg),  Form("CDH seg%d offset", seg), 1000, -10, 10);
    new TH1F(Form("CDH_offset_%d_p", seg),   Form("CDH seg%d offset", seg), 1000, -10, 10);

    new TH1F( Form("CDH_deg_diff_%d", seg), Form("CDH seg%d angle diff", seg), 1000, -45, 45);

    new TH2F( Form("CDH_Z_tsub_%d", seg), Form("CDH z vs tsub seg%d", seg), 500, -200, 200, 500, -25, 25);

    for( int seg2=1; seg2<=5; seg2++ ){
      new TH1F(Form("CDH_offset_%d_T0%d", seg, seg2), Form("CDH offset seg%d T0 seg%d", seg, seg2), 1000, -10, 10); 
    }
  }
  for( int seg2=1; seg2<=5; seg2++ ){
    new TH1F(Form("CDH_offset_T0%d", seg2), Form("CDH offset T0 seg%d", seg2), 1000, -10, 10); 
  }
}

void initHistBLDC()
{
  new TH1F("BLC1a_eff", "BLC1a eff", 10, 0, 10);
  new TH1F("BLC1b_eff", "BLC1a eff", 10, 0, 10);
  new TH1F("BLC2a_eff", "BLC1a eff", 10, 0, 10);
  new TH1F("BLC2b_eff", "BLC1a eff", 10, 0, 10);
  new TH1F("BPC_eff_k",   "BPC eff",   10, 0, 10);
  new TH1F("BPC_eff_pi",   "BPC eff",   10, 0, 10);

  for( int lay=1; lay<=8; lay++ ){
    new TH1F(Form("BLC1a_l%d_res", lay), Form("BLC1a layer%d residual", lay), 1000, -1.0, 1.0);
    new TH1F(Form("BLC1b_l%d_res", lay), Form("BLC1b layer%d residual", lay), 1000, -1.0, 1.0);
    new TH1F(Form("BLC2a_l%d_res", lay), Form("BLC2a layer%d residual", lay), 1000, -1.0, 1.0);
    new TH1F(Form("BLC2b_l%d_res", lay), Form("BLC2b layer%d residual", lay), 1000, -1.0, 1.0);
    new TH1F(Form("BPC_l%d_res", lay),   Form("BPC layer%d residual", lay), 1000, -1.0, 1.0);
    new TH1F(Form("FDC1_l%d_res", lay),  Form("FDC1 layer%d residual", lay), 1000, -1.0, 1.0);

    new TH1F(Form("BLC1a_l%d_dt", lay), Form("BLC1a layer%d dt", lay), 1000, -100, 400);
    new TH1F(Form("BLC1b_l%d_dt", lay), Form("BLC1a layer%d dt", lay), 1000, -100, 400);
    new TH1F(Form("BLC2a_l%d_dt", lay), Form("BLC1a layer%d dt", lay), 1000, -100, 400);
    new TH1F(Form("BLC2b_l%d_dt", lay), Form("BLC1a layer%d dt", lay), 1000, -100, 400);
    new TH1F(Form("BPC_l%d_dt",   lay), Form("BPC layer%d dt",   lay), 1000, -100, 400);
    new TH1F(Form("FDC1_l%d_dt",  lay), Form("FDC1 layer%d dt",  lay), 1000, -100, 400);

    // new TH2F(Form("BLC1a_l%d_dx_res", lay), Form("BLC1a layer%d dx vs res", lay), 100, -0.5, 0.5, 100, -0.1, 0.1);
    // new TH2F(Form("BLC1b_l%d_dx_res", lay), Form("BLC2a layer%d dx vs res", lay), 100, -0.5, 0.5, 100, -0.1, 0.1);
    // new TH2F(Form("BLC2a_l%d_dx_res", lay), Form("BLC2a layer%d dx vs res", lay), 100, -0.5, 0.5, 100, -0.1, 0.1);
    // new TH2F(Form("BLC2b_l%d_dx_res", lay), Form("BLC2b layer%d dx vs res", lay), 100, -0.5, 0.5, 100, -0.1, 0.1);
    // new TH2F(Form("BPC_l%d_dx_res",   lay), Form("BPC layer%d dx vs res",   lay), 100, -0.5, 0.5, 100, -0.1, 0.1);
    // new TH2F(Form("FDC1_l%d_dx_res",  lay), Form("FDC1 layer%d dx vs res",  lay), 100, -0.5, 0.5, 100, -0.1, 0.1);

    for( int wire=1; wire<=32; wire++ ){
      new TH1F(Form("BLC1a_l%d_w%d_dt", lay, wire), Form("BLC1a l%d w%d dt", lay, wire), 1000, -100, 400);
      new TH1F(Form("BLC1b_l%d_w%d_dt", lay, wire), Form("BLC1b l%d w%d dt", lay, wire), 1000, -100, 400);
      new TH1F(Form("BLC2a_l%d_w%d_dt", lay, wire), Form("BLC2a l%d w%d dt", lay, wire), 1000, -100, 400);
      new TH1F(Form("BLC2b_l%d_w%d_dt", lay, wire), Form("BLC2b l%d w%d dt", lay, wire), 1000, -100, 400);
    }

    for( int wire=1; wire<=15; wire++ ){
      new TH1F(Form("BPC_l%d_w%d_dt", lay, wire), Form("BPC l%d w%d dt", lay, wire), 1000, -100, 400);
    }

    for( int wire=1; wire<=62; wire++ ){
      new TH1F(Form("FDC1_l%d_w%d_dt", lay, wire), Form("FDC1 l%d w%d dt", lay, wire), 1000, -100, 400);
    }
  }
}

void initHistT0NC()
{
  for( int seg=1; seg<=5; seg++ ){
    //    new TNtuple(Form("T0%d_slewing_info",seg), Form("T0%d_slewing_info",seg), "eu:ed:tu:td:ctu:ctd:ctm");
    new TH1F(Form("T0NC_offset_T0%d", seg), Form("T0-NC offset T0 seg%d", seg), 1000, -5, 5);
  }

  for( int seg=1; seg<=112; seg++ ){
    //    new TNtuple(Form("NC%d_slewing_info",seg), Form("NC%d_slewing_info",seg), "eu:ed:tu:td:ctu:ctd:ctm");
    new TH1F(Form("T0NC_offset_NC%d_All", seg), Form("T0-NC offset NC seg%d", seg), 1000, -5, 5);
    new TH1F(Form("T0NC_offset_NC%d", seg), Form("T0-NC offset NC seg%d", seg), 1000, -5, 5);

    for( int seg2=1; seg2<=5; seg2++ ){
      // new TNtuple(Form("NC%d_slewing_info_T0%d",seg, seg2), Form("NC%d_slewing_info_T0%d",seg, seg2), "eu:ed:tu:td:ctu:ctd:ctm");
      // new TNtuple(Form("T0%d_slewing_info_NC%d",seg2, seg), Form("T0%d_slewing_info_NC%d",seg2, seg), "eu:ed:tu:td:ctu:ctd:ctm");
      new TH1F(Form("T0NC_offset_T0%d_NC%d", seg2, seg), Form("T0-NC offset T0 seg%d NC seg%d", seg2, seg), 1000, -5, 5);
    }

    new TH1F(Form("NC_tsub_%d", seg), Form("NC tsub seg%d", seg), 1000, -50, 50);
    new TH1F(Form("NC_hitpos_%d", seg), Form("NC hitpos seg%d", seg), 1000, -150, 150);
  }
}

void initHistCDC()
{
  new TH2F("CDC_ob2_res", "CDC 1/#beta^{2} vx residual", 500,   0,  100, 500,  -25,  25);
  for( int l=1; l<=NumOfCDCLayers; l++ ){
    new TH1F(Form("hitpatCDC_%d", l), Form("CDC hit pattern l%d", l), NumOfCDCWiresInLayer[l-1], 0.5, NumOfCDCWiresInLayer[l-1]);
    new TH1F(Form("hitpatCDC_%d_track", l), Form("CDC hit pattern l%d", l), NumOfCDCWiresInLayer[l-1], 0.5, NumOfCDCWiresInLayer[l-1]);

    // new TH2F(Form("CDC_dt_dl_%d", l), Form("CDC dt vs dl l%d", l),                  400, -20, 380, 600, -0.1,  1.1);
    // new TH2F(Form("CDC_dt_res_%d", l), Form("CDC dt vs residual l%d", l),           400, -20, 380, 500, -0.25, 0.25);
    // new TH2F(Form("CDC_ob2_res_%d", l), Form("CDC 1/#beta^{2} vx residual l%d", l), 500,   0, 100, 500,  -25,  25);

    // for( int w=1; w<=NumOfCDCWiresInLayer[l-1]; w++ ){
    //   new TH1F(Form("CDC_dt_l%d_w%d", l ,w),  Form("CDC dt layer%d wire%d", l, w),       1000, -100, 400);
    //   new TH1F(Form("CDC_res_l%d_w%d", l ,w), Form("CDC residual layer%d wire%d", l, w), 1000, -0.1, 0.1);
    // }
  }
}

void initHistIH()
{
  new TH1F("IH_trig", "IH trigger",   24, 0.5, 24.5);
  new TH1F("IH_eff",  "IH effective", 24, 0.5, 24.5);

  for( int seg=1; seg<=24; seg++ ){
    new TH1F( Form("IH_ADC_%d", seg), Form("IH ADC seg%d", seg), 4000,  0, 4000);
    new TH1F( Form("IH_ADC_%d_woT", seg), Form("IH ADC seg%d", seg), 4000,  0, 4000);
    new TH1F( Form("IH_ADC_%d_cut", seg), Form("IH ADC seg%d", seg), 4000,  0, 4000);
    new TH1F( Form("IH_dE_%d", seg), Form("IH dE seg%d", seg),   1000, -0.5, 9.5);
    new TH1F( Form("IH_dE_%d_woT", seg), Form("IH dE seg%d", seg),   1000, -0.5, 9.5);
    new TH1F( Form("IH_dE_%d_cut", seg), Form("IH dE seg%d", seg),   1000, -0.5, 9.5);
    new TH1F( Form("IH_offset_%d", seg), Form("IH offset seg%d", seg), 5000, -1000, 1000);
    new TH1F( Form("IH_time_%d", seg), Form("IH time seg%d", seg), 1000, -50, 50);

    new TH2F( Form("IH_Z_dE_%d_cut", seg), Form("IH z pos vs dE seg%d", seg), 500, -100, 100, 500, -0.5, 9.5);

    new TH1F( Form("IH_deg_diff_%d", seg), Form("IH seg%d angle diff", seg), 1000, -45, 45);

    new TH2F( Form("IH_dE_offset_%d", seg), Form("IH dE vs offset seg%d", seg), 500, 0, 10, 500, -25, 25);
  }
}

void fillHistIH(CDSHitMan *cdsMan)
{
  TH1F *h1;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan->CDH(i);
    int seg = hit-> seg();
    double dE = hit-> emean(), time = hit->ctmean();
    h1 = (TH1F*)gFile-> Get(Form("CDH_dE_woT_%d", seg)), h1-> Fill(dE);
    if( hit-> CheckRange() ){
      h1 = (TH1F*)gFile-> Get(Form("CDH_dE_%d", seg)), h1-> Fill(dE);
      h1 = (TH1F*)gFile-> Get(Form("CDH_time_%d", seg)), h1-> Fill(time);
    }
  }

  for( int i=0; i<cdsMan->nIH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan->IH(i);
    int seg=hit-> seg(), adc=hit->adcu(), tdc = hit->tdcu();
    double dE = hit->eu(), time = hit->ctu();
    h1 = (TH1F*)gFile-> Get(Form("IH_ADC_%d_woT", seg)), h1-> Fill(adc);
    h1 = (TH1F*)gFile-> Get(Form("IH_dE_%d_woT", seg)), h1-> Fill(dE);
    if( 0<tdc && tdc<4000 ){
      h1 = (TH1F*)gFile-> Get(Form("IH_ADC_%d", seg)), h1-> Fill(adc);
      h1 = (TH1F*)gFile-> Get(Form("IH_dE_%d", seg)), h1-> Fill(dE);
      h1 = (TH1F*)gFile-> Get(Form("IH_time_%d", seg)), h1-> Fill(time);
    }
  }
}
