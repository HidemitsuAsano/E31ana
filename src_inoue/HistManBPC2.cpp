#include "HistManBPC2.h"

HistManBPC2::HistManBPC2(TFile *f) : fFile(f)
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


  //### for Trigger ###//
  new TH1F("TrigT_ADC",    "Trigger Top ADC",    4000, -0.5, 3999.5);
  new TH1F("TrigT_ADC_wT", "Trigger Top ADC",    4000, -0.5, 3999.5);
  new TH1F("TrigT_TDC",    "Trigger Top TDC",    4000, -0.5, 3999.5);
  new TH1F("TrigB_ADC",    "Trigger Buttom ADC", 4000, -0.5, 3999.5);
  new TH1F("TrigB_ADC_wT", "Trigger Buttom ADC", 4000, -0.5, 3999.5);
  new TH1F("TrigB_TDC",    "Trigger Buttom TDC", 4000, -0.5, 3999.5);
}

void HistManBPC2::fill(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, DCEffMan *effMan)
{
  TH1* h1;
  TH2 *h2;
  //### for Trigger ###//
  bool def_hit=false;
  for( int i=0; i<blMan->nDEF(); i++ ){
    HodoscopeLikeHit *hit = blMan->DEF(i);
    if( hit-> CheckRange() ) def_hit=true;
  }

  if( def_hit ){
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
	TVector3 hitpos = track-> GetPosatZ(wz);

	h2 = getTH2(Form("BPC_prof_l%d", lay)), h2-> Fill(hitpos.X(), hitpos.Y());
	if( c_hit-> xy()==0 ){
	  double wx = c_hit->wx();
	  double dl = c_hit->dl();
	  double res;
	  if( c_hit->leftright()==0 ) res = hitpos.x()-(wx-dl);
	  if( c_hit->leftright()==1 ) res = hitpos.x()-(wx+dl);
	  h1 = getTH(Form("BPC_res_l%d",lay)), h1-> Fill(res);
	  h2 = getTH2(Form("BPC_xt_l%d",lay)), h2-> Fill(hitpos.X()-wx, time);
	}
	else if( c_hit-> xy()==1 ){
	  double wy = c_hit->wy();
	  double dl = c_hit->dl();
	  double res;
	  if( c_hit->leftright()==0 ) res = hitpos.Y()-(wy-dl);
	  if( c_hit->leftright()==1 ) res = hitpos.Y()-(wy+dl);
	  h1 = getTH(Form("BPC_res_l%d",lay)), h1-> Fill(res);
	  h2 = getTH2(Form("BPC_xt_l%d",lay)), h2-> Fill(hitpos.Y()-wy, time);
	}
      }
    }
  }
}

TH1* HistManBPC2::getTH(const char *name)
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

TH2* HistManBPC2::getTH2(const char *name)
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
