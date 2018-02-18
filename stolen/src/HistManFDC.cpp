#include "HistManFDC.h"

HistManFDC::HistManFDC(TFile *f) : fFile(f)
{
  fFile-> cd();
  for( int lay=1; lay<=6; lay++ ){
    new TH1F(Form("hitpat_FDC_l%d", lay), Form("FDC layer%d, Hit Pattern", lay), 64, 0.5, 64.5);
    new TH1F(Form("FDC_tdc_l%d",lay), Form("FDC layer%d TDC", lay), 4000, -0.5, 3999.5);
    new TH1F(Form("FDC_time_l%d",lay), Form("FDC layer%d time", lay), 4000, -500, 500);
    new TH1F(Form("FDC_time_l%d_1hit" ,lay), Form("FDC layer%d time 1hit", lay), 4000, -500, 500);
    new TH1F(Form("mul_FDC_l%d", lay), Form("FDC mutiplicity layer%d", lay), 65, -0.5, 64.5);
    new TH1F(Form("FDC_dl_l%d", lay), Form("FDC drift length layer%d", lay), 5000, -0.05, 0.45);
    for( int wire=1; wire<=64; wire++ ){
      new TH1F(Form("FDC_tdc_l%d_w%d", lay, wire), Form("FDC Layer%d, wire%d TDC", lay, wire), 4000, -0.5, 3999.5);
      new TH1F(Form("FDC_time_l%d_w%d", lay, wire), Form("FDC Layer%d, wire%d time", lay, wire), 4000, -500, 500);
      new TH1F(Form("FDC_time_l%d_w%d_1hit", lay, wire), Form("FDC Layer%d, wire%d time 1hit", lay, wire), 4000, -500, 500);
    }
  }
  new TH1F("FDC_eff_16", "FDC Efficiency 1&6 layer hit", 8, -0.5, 7.5);
  new TH1F("FDC_ohit_trig", "FDC Trigger",   8, 0.5, 8.5);
  new TH1F("FDC_ohit_eff",  "FDC Efficency", 8, 0.5, 8.5);

  new TH1F("FDC_trig_tracking", "FDC Trigger Event", 6, 0.5, 6.5);
  new TH1F("FDC_eff_tracking", "FDC Effective Event", 6, 0.5, 6.5);

  for( int lay=0; lay<=6; lay++ ){
    new TH2F(Form("FDC_prof_l%d", lay), Form("FDC profile l%d", lay), 500, -10, 10, 500, -10, 10);
  }

  new TH1F("FDC_status", "FDC Tracking Status", 10, -0.5, 9.5);

  //### for Trigger ###//
  new TH1F("TrigT_ADC",    "Trigger Top ADC",    4000, -0.5, 3999.5);
  new TH1F("TrigT_ADC_wT", "Trigger Top ADC",    4000, -0.5, 3999.5);
  new TH1F("TrigT_TDC",    "Trigger Top TDC",    4000, -0.5, 3999.5);
  new TH1F("TrigB_ADC",    "Trigger Buttom ADC", 4000, -0.5, 3999.5);
  new TH1F("TrigB_ADC_wT", "Trigger Buttom ADC", 4000, -0.5, 3999.5);
  new TH1F("TrigB_TDC",    "Trigger Buttom TDC", 4000, -0.5, 3999.5);
}

void HistManFDC::fill(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, DCEffMan *effMan, TKOHitCollection *tko)
{
  TH1* h1;
  TH2 *h2;
  int adc1=-999, adc2=-999, tdc1=-999, tdc2=-999;
  for( int i=0; i<tko->entries(); i++ ){
    TKOHit* hit = tko->hit(i);
    if( hit-> cr()==2 && hit->sl()==21 && hit->ch()==30 ) adc1=hit-> data();
    if( hit-> cr()==2 && hit->sl()==21 && hit->ch()==31 ) adc2=hit-> data();
    if( hit-> cr()==2 && hit->sl()==11 && hit->ch()==13 ) tdc1=hit-> data();
    if( hit-> cr()==2 && hit->sl()==11 && hit->ch()==14 ) tdc2=hit-> data();
  }
  h1 = (TH1F*)fFile-> Get("TrigT_ADC"), h1-> Fill(adc1);
  h1 = (TH1F*)fFile-> Get("TrigB_ADC"), h1-> Fill(adc2);
  h1 = (TH1F*)fFile-> Get("TrigT_TDC"), h1-> Fill(tdc1);
  h1 = (TH1F*)fFile-> Get("TrigB_TDC"), h1-> Fill(tdc2);
  if( tdc1>0 ) h1 = (TH1F*)fFile-> Get("TrigT_ADC_wT"), h1-> Fill(adc1);
  if( tdc2>0 ) h1 = (TH1F*)fFile-> Get("TrigB_ADC_wT"), h1-> Fill(adc2);


  //  if( tdc1>-1 && tdc2>-1 ){
  {
    h1 = getTH("FDC_status"), h1->Fill(bltrackMan->status(CID_FDC1));
    for( int lay=1; lay<=6; lay++ ){
      h1 = getTH(Form("mul_FDC_l%d", lay)), h1-> Fill(blMan->nFDC1(lay));
      for( int i=0; i<blMan->nFDC1(lay); i++ ){
	ChamberLikeHit *hit = blMan-> FDC1(lay, i);
	int wire = hit-> wire();
	int tdc = hit-> tdc();
	double dl = hit->dl();
	double time = hit-> dt();
	h1 = getTH(Form("FDC_tdc_l%d", lay)), h1-> Fill(tdc);
	h1 = getTH(Form("FDC_time_l%d", lay)), h1-> Fill(time);
	h1 = getTH(Form("FDC_tdc_l%d_w%d", lay, wire)), h1-> Fill(tdc);
	h1 = getTH(Form("FDC_time_l%d_w%d", lay, wire)), h1-> Fill(time);
	h1 = getTH(Form("hitpat_FDC_l%d", lay)), h1-> Fill(wire);
	h1 = getTH(Form("FDC_dl_l%d", lay)), h1-> Fill(dl);
	if( effMan-> all1hit() ){
	  h1 = getTH(Form("FDC_time_l%d_1hit", lay)), h1-> Fill(time);
	  h1 = getTH(Form("FDC_time_l%d_w%d_1hit", lay, wire)), h1-> Fill(time);
	}
      }
    }
    
    h1 = getTH("FDC_eff_16");
    h1-> Fill(0);
    if( effMan->trig18() ){
      h1 = getTH("FDC_eff_16");
      for( int lay=1; lay<=6; lay++ ){
	if( blMan->nFDC1(lay)>0 ) h1->Fill(lay);
      }
      if( bltrackMan->ntrackFDC1()>0 ) h1-> Fill(7);
    }
    
    for( int lay=1; lay<=6; lay++ ){
      if( effMan-> trig(lay) ){
	h1 = getTH("FDC_ohit_trig"), h1->Fill(lay);
	if( effMan-> eff(lay, 5) ) h1 = getTH("FDC_ohit_eff"), h1->Fill(lay);
      }
    }

    LocalTrack *track = effMan->getTrack();
    if( track!=0 ){
      for( int i=0; i<track->nhit(); i++ ){
	ChamberLikeHit *c_hit = track->hit(i);
	int lay = c_hit->layer();
	double wz = c_hit-> wz();
	TVector3 hitpos = track-> GetPosatZ(wz);

	h2 = getTH2(Form("FDC_prof_l%d", lay)), h2-> Fill(hitpos.X(), hitpos.Y());
      }
    }

    double default_res[8];
    for( int lay=1; lay<=6; lay++ ){
      double tresl, eresl;
      conf-> GetReslMapManager()-> GetParam(CID_FDC1, lay, 1, tresl, eresl);
      default_res[lay-1] = tresl;
    }

    for( int lay=1; lay<=6; lay++ ){
      bltrackMan-> Clear();
      for( int wire=1; wire<=64; wire++ ) conf->GetReslMapManager()-> SetParam(CID_FDC1, lay, wire, 1.0e10, 0.0);
      bltrackMan-> LocalTracking(blMan, conf, CID_FDC1, "notiming noslope");
      if( bltrackMan->ntrackFDC1()==1 ){
	LocalTrack *track = bltrackMan-> trackFDC1(0);
	int nw, xy;
	double z, xy0, dxy, wl, tilt, ra;
	conf-> GetBLDCWireMapManager()-> GetParam(CID_FDC1, lay, nw, z, xy, xy0, dxy, wl, tilt, ra);

	double local_x, local_y;
	track-> XYLocalPosatZ(z, local_x, local_y, true);


	bool eff_flag = false;
	for( int i=0; i<blMan->nFDC1(lay); i++ ){
	  ChamberLikeHit *hit = blMan-> FDC1(lay, i);
	  double wx = hit-> wx();
	  double dl = hit-> dl();

	  if( fabs(wx-local_x)<dxy ) eff_flag = true;
	}
	h1 = getTH("FDC_trig_tracking"), h1-> Fill(lay);
	if( eff_flag ) h1 = getTH("FDC_eff_tracking"), h1-> Fill(lay);
      }

      for( int wire=1; wire<=64; wire++ ) conf->GetReslMapManager()-> SetParam(CID_FDC1, lay, wire, default_res[lay-1], 0.0);
    }

  }
}

TH1* HistManFDC::getTH(const char *name)
{
  TH1* h1 = (TH1*)fFile->Get(name);
  if( h1==0 ){
    std::cout<<" !!! getTH21("<<name<<") !!!"<<std::endl;
    exit(-1);
  }
  std::string classname = h1-> ClassName();
  if( classname.find("TH")==std::string::npos ){
    std::cout<<"  !!! HistManFDC::getTH1 error name="<<name<<" !!!"<<std::endl;
    exit(-1);
  }
  else return h1;
}

TH2* HistManFDC::getTH2(const char *name)
{
  TH2* h2 = (TH2*)fFile->Get(name);
  if( h2==0 ){
    std::cout<<" !!! getTH2("<<name<<") !!!"<<std::endl;
    exit(-1);
  }
  std::string classname = h2-> ClassName();
  if( classname.find("TH2")==std::string::npos ){
    std::cout<<"  !!! HistManFDC::getTH1 error name="<<name<<" !!!"<<std::endl;
    exit(-1);
  }
  else return h2;
}
