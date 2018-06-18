#include "MyHistBLC1.h"

void initHistBLC1()
{
  new TH1F("ntrackBLC1_1", "BLC1 ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_2", "BLC1 ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1", "BLC1 ntrack", 10, -0.5, 9.5);
  new TH1F("BLC1_time", "BLC1 time", 500, -200, 300);
  new TH1F("BLC1_chi2", "BLC1 chi-square", 500, 0.0, 500);
}

void fillBLC1(BeamLineTrackMan *bltrackMan)
{
  int nBLC1_0=0, nBLC1_1=0, nBLC1_2=0;
  LocalTrack *trackBLC1=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    double time=bltrackMan->trackBLC1(i)->GetTrackTime();
    MyHistTools::fillTH("BLC1_time", time);
    if( -10<time && time<10  ){
      nBLC1_0++;
      trackBLC1=bltrackMan->trackBLC1(i);
    }
    if( -30<time && time<100 ) nBLC1_1++;
    if( -50<time && time<200 ) nBLC1_2++;
  }

  MyHistTools::fillTH("ntrackBLC1", bltrackMan->ntrackBLC1());
  MyHistTools::fillTH("ntrackBLC1_1", nBLC1_1);
  MyHistTools::fillTH("ntrackBLC1_2", nBLC1_2);

  if( nBLC1_0==1 && nBLC1_1==1 ){
    MyHistTools::fillTH("BLC1_chi2", trackBLC1->chi2all());
    if( trackBLC1->chi2all()<10 ){
      for( int i=0; i<trackBLC1->nhit(); i++ ){
	ChamberLikeHit *hit=trackBLC1->hit(i);
	if( hit->cid()==CID_BLC1a ){
	  MyHistTools::fillTH(Form("BLC1a_dt_%d_%d_track",  hit->layer(),  hit->wire()), hit->dt());
	  MyHistTools::fillTH(Form("BLC1a_dt_%d_track",  hit->layer()), hit->dt());
	}
	else{
	  MyHistTools::fillTH(Form("BLC1b_dt_%d_%d_track",  hit->layer(),  hit->wire()), hit->dt());
	  MyHistTools::fillTH(Form("BLC1b_dt_%d_track",  hit->layer()), hit->dt());
	}
      }

      for( int xy=0; xy<2; xy++ ){
	for( int i=0; i<trackBLC1->ncluster(xy); i++ ){
	  BLDCCluster *cluster = trackBLC1->cluster(xy, i);
	  if( cluster->nhit()!=2 ) continue;
	  int lay=trackBLC1->hit(0)->layer();
	  if( lay%2==0 ) lay=trackBLC1->hit(1)->layer();
	  if( lay%2==0 ) continue;

	  if( cluster->hit(0)->cid()==CID_BLC1a ){
	    MyHistTools::fillTH(Form("BLC1a_%d_tsub_tmean", lay), cluster->GetTimeSub(), cluster->GetTimeMean());
	    MyHistTools::fillTH(Form("BLC1a_%d_tsub_ctm", lay), cluster->GetTimeSub(), cluster->GetCTime());
	  }
	  else if( cluster->hit(0)->cid()==CID_BLC1b ){
	    MyHistTools::fillTH(Form("BLC1b_%d_tsub_tmean", lay), cluster->GetTimeSub(), cluster->GetTimeMean());
	    MyHistTools::fillTH(Form("BLC1b_%d_tsub_ctm", lay), cluster->GetTimeSub(), cluster->GetCTime());
	  }
	}
      }
    }
  }
}
