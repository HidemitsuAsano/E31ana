#include "MyHistBHDT0.h"

using namespace std;

void initHistBHDT0()
{
  new TH1F("BHDT0",    "BHD-T0 TOF", 2500, 0.0, 50);
  new TH1F("BHDT0_K",  "BHD-T0 TOF trigK", 2500, 0.0, 50);
  new TH1F("BHDT0_pi", "BHD-T0 TOF trig#pi", 2500, 0.0, 50);
  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("BHDT0_T0%d",    seg), Form("BHD-T0 T0 seg%d",         seg), 2500, 0.0, 50);
    new TH1F(Form("BHDT0_K_T0%d",  seg), Form("BHD-T0 T0 trigK seg%d",   seg), 2500, 0.0, 50);
    new TH1F(Form("BHDT0_pi_T0%d", seg), Form("BHD-T0 T0 trig#pi seg%d", seg), 2500, 0.0, 50);

    for( int seg2=1; seg2<=20; seg2++ ){
      new TH1F(Form("BHDT0_BHD%d_T0%d",    seg2, seg), Form("BHD-T0 BHD seg%d T0 seg%d",         seg2, seg), 2500, 0.0, 50);
      new TH1F(Form("BHDT0_K_BHD%d_T0%d",  seg2, seg), Form("BHD-T0 trigK BHD seg%d T0 seg%d",   seg2, seg), 2500, 0.0, 50);
      new TH1F(Form("BHDT0_pi_BHD%d_T0%d", seg2, seg), Form("BHD-T0 trig#pi BHD seg%d T0 seg%d", seg2, seg), 2500, 0.0, 50);
    }
  }

  for( int seg=1; seg<=20; seg++ ){
    new TH1F(Form("BHDT0_BHD%d",    seg), Form("BHD-T0 BHD seg%d",         seg), 2500, 0.0, 50);
    new TH1F(Form("BHDT0_K_BHD%d",  seg), Form("BHD-T0 trigK BHD seg%d",   seg), 2500, 0.0, 50);
    new TH1F(Form("BHDT0_pi_BHD%d", seg), Form("BHD-T0 trig#pi BHD seg%d", seg), 2500, 0.0, 50);
  }
}

void fillBHDT0(EventHeader *header, BeamLineHitMan *blMan)
{
  vector<HodoscopeLikeHit*> T0hits = MyTools::getHodo(blMan, CID_T0);
  vector<HodoscopeLikeHit*> BHDhits = MyTools::getHodo(blMan, CID_BHD);

  if( T0hits.size()!=1 ) return;
  HodoscopeLikeHit *T0hit=T0hits[0];

  for( int i=0; i<BHDhits.size(); i++ ){
    HodoscopeLikeHit *BHDhit=BHDhits[i];
    double tof=T0hit->ctmean()-BHDhit->ctmean();
    int T0seg=T0hit->seg(), BHDseg=BHDhit->seg();

    MyHistTools::fillTH("BHDT0", tof);
    MyHistTools::fillTH(Form("BHDT0_T0%d", T0seg), tof);
    MyHistTools::fillTH(Form("BHDT0_BHD%d", BHDseg), tof);
    MyHistTools::fillTH(Form("BHDT0_BHD%d_T0%d", BHDseg, T0seg), tof);
    if( header->IsTrig(Trig_Kaon) ){
      MyHistTools::fillTH("BHDT0_K", tof);
      MyHistTools::fillTH(Form("BHDT0_K_T0%d", T0seg), tof);
      MyHistTools::fillTH(Form("BHDT0_K_BHD%d", BHDseg), tof);
      MyHistTools::fillTH(Form("BHDT0_K_BHD%d_T0%d", BHDseg, T0seg), tof);
    }
    if( header->IsTrig(Trig_Pion) ){
      MyHistTools::fillTH("BHDT0_pi", tof);
      MyHistTools::fillTH(Form("BHDT0_pi_T0%d", T0seg), tof);
      MyHistTools::fillTH(Form("BHDT0_pi_BHD%d", BHDseg), tof);
      MyHistTools::fillTH(Form("BHDT0_pi_BHD%d_T0%d", BHDseg, T0seg), tof);
    }
  }
}

