#include "MyHistBLC2a.h"

void initHistBLC2a()
{
  const int nlay=8;
  const int nw=32;

  for( int lay=1; lay<=nlay; lay++ ){
    new TH1F(Form("hitpatBLC2a_%d", lay), Form("BLC2a hit pattern lay%d", lay), nw, 0.5, nw+0.5);
    new TH1F(Form("mulBLC2a_%d", lay), Form("BLC2a multiplicity lay%d", lay), nw+1, -0.5, nw+0.5);

    new TH1F(Form("BLC2a_TDC_%d", lay), Form("BLC2a TDC lay%d", lay), 2000, 0, 2000);
    new TH1F(Form("BLC2a_dt_%d", lay), Form("BLC2a dt lay%d", lay), 600, -200, 400);

    new TH1F(Form("BLC2a_TDC_%d_1hit", lay), Form("BLC2a TDC lay%d", lay), 2000, 0, 2000);
    new TH1F(Form("BLC2a_dt_%d_1hit", lay), Form("BLC2a dt lay%d", lay), 600, -200, 400);
    new TH1F(Form("BLC2a_dt_%d_track", lay), Form("BLC2a dt lay%d", lay), 600, -200, 400);

    new TH1F(Form("BLC2a_res_%d", lay), Form("BLC2a res lay%d", lay), 500, -0.1, 0.1);

    for( int w=1; w<=nw; w++ ){
      new TH1F(Form("BLC2a_TDC_%d_%d", lay, w), Form("BLC2a TDC lay%d wire%d", lay, w), 2000, 0, 2000);
      new TH1F(Form("BLC2a_dt_%d_%d", lay, w), Form("BLC2a dt lay%d wire%d", lay, w), 500, -100, 400);

      new TH1F(Form("BLC2a_TDC_%d_%d_1hit", lay, w), Form("BLC2a TDC lay%d wire%d", lay, w), 2000, 0, 2000);
      new TH1F(Form("BLC2a_dt_%d_%d_1hit",  lay, w), Form("BLC2a dt lay%d wire%d", lay, w), 500, -100, 400);
      new TH1F(Form("BLC2a_dt_%d_%d_track", lay, w), Form("BLC2a dt lay%d wire%d", lay, w), 500, -100, 400);

      new TH1F(Form("BLC2a_res_%d_%d", lay, w), Form("BLC2a res lay%d", lay, w), 500, -0.1, 0.1);
    }

    if( lay%2==1 ){
      new TH2F(Form("BLC2a_%d_tsub_tmean", lay), Form("BLC2a lay%d-%d time sub vs time mean", lay, lay+1), 500, -250, 250, 500, -100, 400);
      new TH2F(Form("BLC2a_%d_tsub_ctm", lay), Form("BLC2a lay%d-%d time sub vs corr. time mean", lay, lay+1), 500, -250, 250, 500, -100, 400);
    }

  }
}

void fillBLC2a(LocalTrack *track)
{
  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int lay=hit->layer(), wire=hit->wire(), tdc=hit->tdc();
    double dt=hit->dt();
    double res=hit->resl();

    MyHistTools::fillTH(Form("BLC2a_res_%d", lay), res);
    MyHistTools::fillTH(Form("BLC2a_res_%d_%d", lay, wire), res);
  }
}


void fillBLC2a(BeamLineHitMan *blMan)
{
  const int nlay=8;
  bool all_1hit=true;
  for( int lay=1; lay<=nlay; lay++ ){
    if( blMan->nBLC2a(lay)!=1 ) all_1hit=false;
    MyHistTools::fillTH(Form("mulBLC2a_%d", lay), blMan->nBLC2a(lay));

    for( int i=0; i<blMan->nBLC2a(lay); i++ ){
      ChamberLikeHit *hit=blMan->BLC2a(lay, i);
      int wire=hit->wire();
      int tdc=hit->tdc();
      double dt=hit->dt();

      MyHistTools::fillTH(Form("hitpatBLC2a_%d", lay), wire);
      MyHistTools::fillTH(Form("BLC2a_TDC_%d", lay), tdc);
      MyHistTools::fillTH(Form("BLC2a_dt_%d", lay), dt);
      MyHistTools::fillTH(Form("BLC2a_TDC_%d_%d", lay, wire), tdc);
      MyHistTools::fillTH(Form("BLC2a_dt_%d_%d", lay, wire), dt);
    }
  }

  if( all_1hit ){
    for( int lay=1; lay<=nlay; lay++ ){
      ChamberLikeHit *hit=blMan->BLC2a(lay, 0);
      int wire=hit->wire();
      int tdc=hit->tdc();
      double dt=hit->dt();

      MyHistTools::fillTH(Form("BLC2a_TDC_%d_1hit", lay), tdc);
      MyHistTools::fillTH(Form("BLC2a_dt_%d_1hit", lay), dt);
      MyHistTools::fillTH(Form("BLC2a_TDC_%d_%d_1hit", lay, wire), tdc);
      MyHistTools::fillTH(Form("BLC2a_dt_%d_%d_1hit", lay, wire), dt);
    }
  }
}
