#include "MyHistBPC.h"

using namespace std;

void initHistBPC()
{
  const int nlay=8;
  const int nw=32;

  for( int lay=1; lay<=nlay; lay++ ){
    new TH1F(Form("hitpatBPC_%d", lay), Form("BPC hit pattern lay%d", lay), nw, 0.5, nw+0.5);
    new TH1F(Form("mulBPC_%d", lay), Form("BPC multiplicity lay%d", lay), nw+1, -0.5, nw+0.5);

    new TH1F(Form("BPC_TDC_%d", lay), Form("BPC TDC lay%d", lay), 2000, 0, 2000);
    new TH1F(Form("BPC_dt_%d", lay), Form("BPC dt lay%d", lay), 500, -100, 400);

    new TH1F(Form("BPC_TDC_%d_1hit", lay), Form("BPC TDC lay%d", lay), 2000, 0, 2000);
    new TH1F(Form("BPC_dt_%d_1hit", lay), Form("BPC dt lay%d", lay), 500, -100, 400);
    new TH1F(Form("BPC_dt_%d_track", lay), Form("BPC dt lay%d", lay), 500, -100, 400);

    new TH1F(Form("BPC_res_%d", lay), Form("BPC res lay%d", lay), 500, -0.1, 0.1);

    if( lay%2==1 ){
      new TH2F(Form("BPC_%d_tsub_tmean", lay), Form("BPC lay%d-%d time sub vs time average", lay, lay+1), 500, -250, 250, 500, -100, 400);
      new TH2F(Form("BPC_%d_tsub_ctm", lay), Form("BPC lay%d-%d time sub vs corr. time average", lay, lay+1), 500, -250, 250, 500, -100, 400);
    }

    for( int w=1; w<=nw; w++ ){
      new TH1F(Form("BPC_TDC_%d_%d", lay, w), Form("BPC TDC lay%d wire%d", lay, w), 2000, 0, 2000);
      new TH1F(Form("BPC_dt_%d_%d", lay, w), Form("BPC dt lay%d wire%d", lay, w), 500, -100, 400);

      new TH1F(Form("BPC_TDC_%d_%d_1hit", lay, w), Form("BPC TDC lay%d wire%d", lay, w), 2000, 0, 2000);
      new TH1F(Form("BPC_dt_%d_%d_1hit", lay, w), Form("BPC dt lay%d wire%d", lay, w), 500, -100, 400);
      new TH1F(Form("BPC_dt_%d_%d_track", lay, w), Form("BPC dt lay%d wire%d", lay, w), 500, -100, 400);
      new TH1F(Form("BPC_res_%d_%d", lay, w), Form("BPC res lay%d wire%d", lay, w), 500, -0.1, 0.1);
    }
  }

  new TH1F("ntrackBPC_1", "BPC ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_2", "BPC ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBPC", "BPC ntrack", 10, -0.5, 9.5);
  new TH1F("BPC_time", "BPC time", 500, -200, 300);
  new TH1F("BPC_chi2", "BPC chi-square", 500, 0.0, 500);
}

void fillBPC(LocalTrack *track)
{
  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int lay=hit->layer(), wire=hit->wire(), tdc=hit->tdc();
    double dt=hit->dt();
    double res=hit->resl();

    MyHistTools::fillTH(Form("BPC_res_%d", lay), res);
    MyHistTools::fillTH(Form("BPC_res_%d_%d", lay, wire), res);
  }
}

void fillBPC(BeamLineHitMan *blMan)
{
  const int nlay=8;
  bool all_1hit=true;
  for( int lay=1; lay<=nlay; lay++ ){
    if( blMan->nBPC(lay)!=1 ) all_1hit=false;
    MyHistTools::fillTH(Form("mulBPC_%d", lay), blMan->nBPC(lay));

    for( int i=0; i<blMan->nBPC(lay); i++ ){
      ChamberLikeHit *hit=blMan->BPC(lay, i);
      int wire=hit->wire();
      int tdc=hit->tdc();
      double dt=hit->dt();

      MyHistTools::fillTH(Form("hitpatBPC_%d", lay), wire);
      MyHistTools::fillTH(Form("BPC_TDC_%d", lay), tdc);
      MyHistTools::fillTH(Form("BPC_dt_%d", lay), dt);
      MyHistTools::fillTH(Form("BPC_TDC_%d_%d", lay, wire), tdc);
      MyHistTools::fillTH(Form("BPC_dt_%d_%d", lay, wire), dt);
    }
  }

  if( all_1hit ){
    for( int lay=1; lay<=nlay; lay++ ){
      ChamberLikeHit *hit=blMan->BPC(lay, 0);
      int wire=hit->wire();
      int tdc=hit->tdc();
      double dt=hit->dt();

      MyHistTools::fillTH(Form("BPC_TDC_%d_1hit", lay), tdc);
      MyHistTools::fillTH(Form("BPC_dt_%d_1hit", lay), dt);
      MyHistTools::fillTH(Form("BPC_TDC_%d_%d_1hit", lay, wire), tdc);
      MyHistTools::fillTH(Form("BPC_dt_%d_%d_1hit", lay, wire), dt);
    }
  }
}

void fillBPC(BeamLineTrackMan *bltrackMan)
{
  int nBPC_0=0, nBPC_1=0, nBPC_2=0;
  LocalTrack *trackBPC=0;
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    double time=bltrackMan->trackBPC(i)->GetTrackTime();
    MyHistTools::fillTH("BPC_time", time);
    if( -10<time && time<10  ){
      nBPC_0++;
      trackBPC=bltrackMan->trackBPC(i);
    }
    if( -30<time && time<100 ) nBPC_1++;
    if( -50<time && time<200 ) nBPC_2++;
  }

  MyHistTools::fillTH("ntrackBPC", bltrackMan->ntrackBPC());
  MyHistTools::fillTH("ntrackBPC_1", nBPC_1);
  MyHistTools::fillTH("ntrackBPC_2", nBPC_2);

  if( nBPC_0==1 && nBPC_1==1 ){
    MyHistTools::fillTH("ntrackBPC", trackBPC->GetTrackTime());
    MyHistTools::fillTH("BPC_chi2", trackBPC->chi2all());
    if( trackBPC->chi2all()<10 ){
      for( int i=0; i<trackBPC->nhit(); i++ ){
	ChamberLikeHit *hit=trackBPC->hit(i);
	MyHistTools::fillTH(Form("BPC_dt_%d_%d_track",  hit->layer(),  hit->wire()), hit->dt());
	MyHistTools::fillTH(Form("BPC_dt_%d_track",  hit->layer()), hit->dt());
      }

      for( int xy=0; xy<2; xy++ ){
	for( int i=0; i<trackBPC->ncluster(xy); i++){
	  BLDCCluster *cluster=trackBPC->cluster(xy, i);

	  if( cluster->nhit()!=2 ) continue;
	  int lay=cluster->hit(0)->layer();
	  if( lay%2==0 ) lay=cluster->hit(1)->layer();
	  if( lay%2==0 ) continue;

	  MyHistTools::fillTH(Form("BPC_%d_tsub_tmean", lay), cluster->GetTimeSub(), cluster->GetTimeMean());
	  MyHistTools::fillTH(Form("BPC_%d_tsub_ctm", lay), cluster->GetTimeSub(), cluster->GetCTime());
	}
      }
    }
  }
}
