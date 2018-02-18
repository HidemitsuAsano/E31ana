#include "HistManBPC.h"

HistManBPC::HistManBPC(TFile *f) : fFile(f)
{
  fFile-> cd();
  for( int lay=1; lay<=8; lay++ ){
    new TH1F(Form("hitpat_BPC_l%d", lay), Form("BPC layer%d, Hit Pattern", lay), 15, 0.5, 15.5);
    new TH1F(Form("BPC_tdc_l%d",lay), Form("BPC layer%d TDC", lay), 4000, -0.5, 3999.5);
    new TH1F(Form("BPC_time_l%d",lay), Form("BPC layer%d time", lay), 4000, -500, 500);
    new TH1F(Form("BPC_time_l%d_1hit" ,lay), Form("BPC layer%d time 1hit", lay), 4000, -500, 500);
    new TH1F(Form("mul_BPC_l%d", lay), Form("BPC mutiplicity layer%d", lay), 16, -0.5, 15.5);
    new TH1F(Form("BPC_dl_l%d", lay), Form("BPC drift length layer%d", lay), 5000, -0.05, 0.45);
    new TH1F(Form("BPC_res_l%d", lay), Form("BPC residual layer%d", lay), 1000, -1.0, 1.0);
    new TH2F(Form("BPC_xt_l%d", lay), Form("BPC dl vs time layer%d", lay), 1000, -1.0, 1.0, 1000, -10, 140);
    for( int wire=1; wire<16; wire++ ){
      new TH1F(Form("BPC_tdc_l%d_w%d", lay, wire), Form("BPC Layer%d, wire%d TDC", lay, wire), 4000, -0.5, 3999.5);
      new TH1F(Form("BPC_time_l%d_w%d", lay, wire), Form("BPC Layer%d, wire%d time", lay, wire), 4000, -500, 500);
      new TH1F(Form("BPC_time_l%d_w%d_1hit", lay, wire), Form("BPC Layer%d, wire%d time 1hit", lay, wire), 4000, -500, 500);
    }
  }
  new TH1F("BPC_eff_18", "BPC Efficiency 1&8 layer hit", 10, -0.5, 9.5);
  new TH1F("BPC_ohit_trig", "BPC Trigger",   8, 0.5, 8.5);
  new TH1F("BPC_ohit_eff",  "BPC Efficency", 8, 0.5, 8.5);

  for( int lay=0; lay<=8; lay++ ){
    new TH2F(Form("BPC_prof_l%d", lay), Form("BPC profile l%d", lay), 500, -10, 10, 500, -10, 10);
  }

  new TH1F("BPC_status", "BPC Tracking Status", 10, -0.5, 9.5);
  new TH1F("BPC_trig_tracking", "BPC Triger Event", 8, 0.5, 8.5);
  new TH1F("BPC_eff_tracking", "BPC Effective Event", 8, 0.5, 8.5);

  new TH1F("BPC_trig_tracking_cut", "BPC Triger Event", 8, 0.5, 8.5);
  new TH1F("BPC_eff_tracking_cut", "BPC Effective Event", 8, 0.5, 8.5);

  //### for Trigger ###//
  new TH1F("TrigT_ADC",    "Trigger Top ADC",    4000, -0.5, 3999.5);
  new TH1F("TrigT_ADC_wT", "Trigger Top ADC",    4000, -0.5, 3999.5);
  new TH1F("TrigT_TDC",    "Trigger Top TDC",    4000, -0.5, 3999.5);
  new TH1F("TrigB_ADC",    "Trigger Buttom ADC", 4000, -0.5, 3999.5);
  new TH1F("TrigB_ADC_wT", "Trigger Buttom ADC", 4000, -0.5, 3999.5);
  new TH1F("TrigB_TDC",    "Trigger Buttom TDC", 4000, -0.5, 3999.5);
}

void HistManBPC::fill(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, DCEffMan *effMan)
{
  TH1* h1;
  TH2 *h2;
  //### for Trigger ###//
  HodoscopeLikeHit *hit=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->seg()==1 ) hit=blMan->T0(i);
  }
  int tdc1=-999, tdc2=-999, adc1=-999, adc2=-999;
  if( hit ) tdc1 = hit-> tdcu(), tdc2 = hit-> tdcd(), adc1 = hit-> adcu(), adc2 = hit-> adcd();
  h1 = getTH("TrigT_ADC"), h1->Fill(adc1);
  if( tdc1>0 ) h1 = getTH("TrigT_ADC_wT"), h1->Fill(adc1);
  h1 = getTH("TrigB_ADC"), h1->Fill(adc2);
  if( tdc2>0 ) h1 = getTH("TrigB_ADC_wT"), h1->Fill(adc2);
  h1 = getTH("TrigT_TDC"), h1->Fill(tdc1);
  h1 = getTH("TrigB_TDC"), h1->Fill(tdc2);


  if( adc1>90 && adc2>140 ){
    h1 = getTH("BPC_status"), h1->Fill(bltrackMan->status(CID_BPC));
    for( int lay=1; lay<=8; lay++ ){
      h1 = getTH(Form("mul_BPC_l%d", lay)), h1-> Fill(blMan->nBPC(lay));
      for( int i=0; i<blMan->nBPC(lay); i++ ){
	ChamberLikeHit *hit = blMan-> BPC(lay, i);
	int wire = hit-> wire();
	int tdc = hit-> tdc();
	double dl = hit->dl();
	double time = hit-> dt();
	h1 = getTH(Form("BPC_tdc_l%d", lay)), h1-> Fill(tdc);
	h1 = getTH(Form("BPC_time_l%d", lay)), h1-> Fill(time);
	h1 = getTH(Form("BPC_tdc_l%d_w%d", lay, wire)), h1-> Fill(tdc);
	h1 = getTH(Form("BPC_time_l%d_w%d", lay, wire)), h1-> Fill(time);
	h1 = getTH(Form("hitpat_BPC_l%d", lay)), h1-> Fill(wire);
	h1 = getTH(Form("BPC_dl_l%d", lay)), h1-> Fill(dl);
	if( effMan-> all1hit() ){
	  h1 = getTH(Form("BPC_time_l%d_1hit", lay)), h1-> Fill(time);
	  h1 = getTH(Form("BPC_time_l%d_w%d_1hit", lay, wire)), h1-> Fill(time);
	}
      }
    }
    
    h1 = getTH("BPC_eff_18");
    h1-> Fill(0);
    if( effMan->trig18() ){
      h1 = getTH("BPC_eff_18");
      for( int lay=1; lay<=8; lay++ ){
	if( blMan->nBPC(lay)>0 ) h1->Fill(lay);
      }
      if( bltrackMan->ntrackBPC()>0 ) h1-> Fill(9);    
    }

    for( int lay=1; lay<=8; lay++ ){
      if( effMan-> trig(lay) ){
	h1 = getTH("BPC_ohit_trig"), h1->Fill(lay);
	if( effMan-> eff(lay, 5) ) h1 = getTH("BPC_ohit_eff"), h1->Fill(lay);
      }
    }

    LocalTrack *track = effMan->getTrack();
    if( track!=0 ){
      for( int i=0; i<track->nhit(); i++ ){
	ChamberLikeHit *c_hit = track->hit(i);
	double time = c_hit-> dt();
	int lay = c_hit->layer();
	double wz = c_hit-> wz();
	double dl = c_hit-> dl();
	double dt = c_hit-> dt();
	double resl = c_hit-> resl();
	double local_x, local_y;
	track-> XYLocalPosatZ(wz, local_x, local_y);
	double dltrack;
	if( c_hit-> xy()==0 ) dltrack = c_hit->x()-c_hit->wx();
	else if( c_hit-> xy()==1 ) dltrack = c_hit->y()-c_hit->wy();

	h2 = getTH2(Form("BPC_prof_l%d", lay)), h2-> Fill(local_x, local_y);
	h1 = getTH(Form("BPC_res_l%d",lay)), h1-> Fill(resl);
	h2 = getTH2(Form("BPC_xt_l%d",lay)), h2-> Fill(dltrack, dt);
      }
    }

    double default_res[8];
    for( int lay=1; lay<=8; lay++ ){
      double tresl, eresl;
      conf-> GetReslMapManager()-> GetParam(CID_BPC, lay, 1, tresl, eresl);
      default_res[lay-1] = tresl;
    }


    for( int lay=1; lay<=8; lay++ ){
      bltrackMan-> Clear();
      for( int wire=1; wire<=15; wire++ ){
	conf-> GetReslMapManager()-> SetParam(CID_BPC, lay, wire, 1.0e10, 0.0);
      }

      bool trig_flag=false;
      bltrackMan-> LocalTracking(blMan, conf, CID_BPC, "notiming noslope");
      if( bltrackMan-> ntrackBPC()==1 ){
	LocalTrack *track = bltrackMan->trackBPC(0);
	double local_x, local_y;
	int nw, xy;
	double z, xy0, dxy, wirel, tilt, ra;
	conf-> GetBLDCWireMapManager()-> GetParam(CID_BPC, lay, nw, z, xy, xy0, dxy, wirel, tilt, ra);
	track-> XYLocalPosatZ(z, local_x, local_y);
	double r = sqrt(local_x*local_x+local_y*local_y);
	if( r<4.5 ) trig_flag=true;
	else trig_flag = false;

	bool eff_flag = false;
	for( int i=0; i<blMan->nBPC(lay); i++ ){
	  ChamberLikeHit *hit = blMan->BPC(lay, i);
	  if( hit->xy()==0 ){
	    if( fabs(hit->wx()-local_x)<1.0*dxy ) eff_flag = true;
	  }
	  else{
	    if( fabs(hit->wy()-local_y)<1.0*dxy ) eff_flag = true;
	  }
	}
	h1 = getTH("BPC_trig_tracking"), h1-> Fill(lay);
	if( eff_flag) h1 = getTH("BPC_eff_tracking"), h1-> Fill(lay);
	if( trig_flag ){
	  h1 = getTH("BPC_trig_tracking_cut"), h1-> Fill(lay);
	  if( eff_flag) h1 = getTH("BPC_eff_tracking_cut"), h1-> Fill(lay);
	}
      }

      for( int wire=1; wire<=15; wire++ ){
	conf-> GetReslMapManager()-> SetParam(CID_BPC, lay, wire, default_res[lay-1], 0.0);
      }
    }
  }
}

TH1* HistManBPC::getTH(const char *name)
{
  TH1* h1 = (TH1*)fFile->Get(name);
  if( h1==0 ){
    std::cout<<" !!! getTH21("<<name<<") !!!"<<std::endl;
    exit(-1);
  }
  std::string classname = h1-> ClassName();
  if( classname.find("TH")==std::string::npos ){
    std::cout<<"  !!! HistManBPC::getTH1 error name="<<name<<" !!!"<<std::endl;
    exit(-1);
  }
  else return h1;
}

TH2* HistManBPC::getTH2(const char *name)
{
  TH2* h2 = (TH2*)fFile->Get(name);
  if( h2==0 ){
    std::cout<<" !!! getTH2("<<name<<") !!!"<<std::endl;
    exit(-1);
  }
  std::string classname = h2-> ClassName();
  if( classname.find("TH2")==std::string::npos ){
    std::cout<<"  !!! HistManBPC::getTH1 error name="<<name<<" !!!"<<std::endl;
    exit(-1);
  }
  else return h2;
}
