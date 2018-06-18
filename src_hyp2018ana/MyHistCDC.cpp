#include "MyHistCDC.h"

void initHistCDC(ConfMan *conf)
{
  new TH1F("CDC_chi2", "CDC chi-square", 1000, 0, 500);

  for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
    int nwire=conf->GetCDCWireMapManager()->nwires(lay);
    new TH1F(Form("mulCDC_%d", lay), Form("CDC layer%d multiplicity", lay), nwire+1, -0.5, nwire+0.5);
    new TH1F(Form("hitpatCDC_%d", lay), Form("CDC layer%d hit pattern", lay), nwire,  0.5, nwire+0.5);
    new TH1F(Form("CDC_res_%d", lay), Form("CDC layer%d residual", lay), 1000, -1.0, 1.0);
    new TH1F(Form("CDC_TDC_%d", lay), Form("CDC layer%d TDC", lay), 2000, 0, 2000);
    new TH1F(Form("CDC_dt_%d", lay), Form("CDC layer%d dt", lay), 500, -100, 400);

    new TH1F(Form("CDC_TDC_%d_1hit", lay), Form("CDC layer%d TDC", lay), 2000, 0, 2000);
    new TH1F(Form("CDC_dt_%d_1hit", lay), Form("CDC layer%d dt", lay), 500, -100, 400);

    new TH2F(Form("CDC_dt_dx_%d", lay), Form("CDC layer%d dt vs dx", lay), 100, -50, 450, 100, -0.3, 1.2);
    new TH2F(Form("CDC_dt_res_%d", lay), Form("CDC layer%d dt vs res", lay), 100, -50, 450, 100, -0.3, 0.3);

    for( int wire=1; wire<=nwire; wire++ ){
      new TH1F(Form("CDC_res_%d_%d", lay, wire), Form("CDC layer%d wire%d residual", lay, wire), 1000, -1.0, 1.0);
      new TH1F(Form("CDC_TDC_%d_%d", lay, wire), Form("CDC layer%d wire%d TDC", lay, wire),      2000, 0, 2000);
      new TH1F(Form("CDC_dt_%d_%d",  lay, wire), Form("CDC layer%d wire%d dt", lay, wire),       500, -100, 400);

      new TH1F(Form("CDC_TDC_%d_%d_1hit", lay, wire), Form("CDC layer%d wire%d TDC", lay, wire),      2000, 0, 2000);
      new TH1F(Form("CDC_dt_%d_%d_1hit",  lay, wire), Form("CDC layer%d wire%d dt", lay, wire),       500, -100, 400);

      new TH1F(Form("CDC_TDC_%d_%d_track", lay, wire), Form("CDC layer%d wire%d TDC", lay, wire),      2000, 0, 2000);
      new TH1F(Form("CDC_dt_%d_%d_track",  lay, wire), Form("CDC layer%d wire%d dt", lay, wire),       500, -100, 400);

      new TH2F(Form("CDC_dt_dx_%d_%d", lay, wire), Form("CDC layer%d wire%d dt vs dx", lay, wire),   100, -50, 450, 100, -0.3, 1.2);
      new TH2F(Form("CDC_dt_res_%d_%d", lay, wire), Form("CDC layer%d wire%d dt vs res", lay, wire), 100, -50, 450, 100, -0.3, 0.3);
    }
  }
}

void fillCDC(CDSTrackingMan *cdstrackMan, CDSHitMan *cdsMan)
{
  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *track=cdstrackMan->GoodTrack(i);
    MyHistTools::fillTH("CDC_chi2", track->Chi());

    if( track->Chi()<30 ){
      for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
	for( int j=0; j<track->nTrackHit(lay); j++ ){
	  CDCHit *hit = track->hit(cdsMan, lay, j);
	  MyHistTools::fillTH(Form("CDC_res_%d", lay), hit->resl());
	  MyHistTools::fillTH(Form("CDC_res_%d_%d", lay, hit->wire()), hit->resl());
	  
	  MyHistTools::fillTH(Form("CDC_TDC_%d_%d_track", lay, hit->wire()), hit->tdc());
	  MyHistTools::fillTH(Form("CDC_dt_%d_%d_track", lay, hit->wire()), hit->dt());

	  MyHistTools::fillTH(Form("CDC_dt_dx_%d", lay), hit->dt(), hit->dl()+hit->resl());
	  MyHistTools::fillTH(Form("CDC_dt_dx_%d_%d", lay, hit->wire()), hit->dt(), hit->dl()+hit->resl());

	  MyHistTools::fillTH(Form("CDC_dt_res_%d", lay), hit->dt(), hit->resl());
	  MyHistTools::fillTH(Form("CDC_dt_res_%d_%d", lay, hit->wire()), hit->dt(), hit->resl());
	}
      }
    }
  }
}

void fillCDC(CDSHitMan *cdsMan)
{
  bool all1hit=MyTools::CDCall1hit(cdsMan);

  for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
    MyHistTools::fillTH(Form("mulCDC_%d", lay), cdsMan->nCDC(lay));
    for( int i=0; i<cdsMan->nCDC(lay); i++ ){
      CDCHit *hit=cdsMan->CDC(lay, i);
      int wire=hit->wire(), tdc=hit->tdc();
      double dt=hit->dt();
      MyHistTools::fillTH(Form("hitpatCDC_%d", lay), wire);
      MyHistTools::fillTH(Form("CDC_TDC_%d", lay), tdc);
      MyHistTools::fillTH(Form("CDC_dt_%d", lay), dt);
      MyHistTools::fillTH(Form("CDC_TDC_%d_%d", lay, wire), tdc);
      MyHistTools::fillTH(Form("CDC_dt_%d_%d", lay, wire), dt);

      if( all1hit ){
	MyHistTools::fillTH(Form("CDC_TDC_%d_%d_1hit", lay, wire), tdc);
	MyHistTools::fillTH(Form("CDC_dt_%d_%d_1hit", lay, wire), dt);
      }
    }
  }
}
