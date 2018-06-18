#include "MyHistBLC2.h"

void initHistBLC2()
{
  new TH1F("ntrackBLC2_1", "BLC2 ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_2", "BLC2 ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2", "BLC2 ntrack", 10, -0.5, 9.5);
  new TH1F("BLC2_time", "BLC2 time", 500, -200, 300);
  new TH1F("BLC2_chi2", "BLC2 chi-square", 500, 0.0, 500);
}

void fillBLC2(BeamLineTrackMan *bltrackMan)
{
  int nBLC2_0=0, nBLC2_1=0, nBLC2_2=0;
  LocalTrack *trackBLC2=0;
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    double time=bltrackMan->trackBLC2(i)->GetTrackTime();
    MyHistTools::fillTH("BLC2_time", time);
    if( -10<time && time<10  ){
      nBLC2_0++;
      trackBLC2=bltrackMan->trackBLC2(i);
    }
    if( -30<time && time<100 ) nBLC2_1++;
    if( -50<time && time<200 ) nBLC2_2++;
  }

  MyHistTools::fillTH("ntrackBLC2", bltrackMan->ntrackBLC2());
  MyHistTools::fillTH("ntrackBLC2_1", nBLC2_1);
  MyHistTools::fillTH("ntrackBLC2_2", nBLC2_2);

  if( nBLC2_0==1 && nBLC2_1==1 ){
    MyHistTools::fillTH("BLC2_chi2", trackBLC2->chi2all());
    if(trackBLC2->chi2all()<10){
      for( int i=0; i<trackBLC2->nhit(); i++ ){
        ChamberLikeHit *hit=trackBLC2->hit(i);
        if( hit->cid()==CID_BLC2a ){
	  MyHistTools::fillTH(Form("BLC2a_dt_%d_%d_track", hit->layer(), hit->wire()), hit->dt());
	  MyHistTools::fillTH(Form("BLC2a_dt_%d_track", hit->layer()), hit->dt());
        }
        else{
	  MyHistTools::fillTH(Form("BLC2b_dt_%d_%d_track", hit->layer(), hit->wire()), hit->dt());
	  MyHistTools::fillTH(Form("BLC2b_dt_%d_track", hit->layer()), hit->dt());
        }
      }

      for( int xy=0; xy<2; xy++ ){
        for( int i=0; i<trackBLC2->ncluster(xy); i++ ){
          BLDCCluster *cluster = trackBLC2->cluster(xy, i);
          if( cluster->nhit()!=2 ) continue;
          int lay=trackBLC2->hit(0)->layer();
          if( lay%2==0 ) lay=trackBLC2->hit(1)->layer();
          if( lay%2==0 ) continue;

          if( cluster->hit(0)->cid()==CID_BLC2a ){
	    MyHistTools::fillTH(Form("BLC2a_%d_tsub_tmean", lay), cluster->GetTimeSub(), cluster->GetTimeMean());
	    MyHistTools::fillTH(Form("BLC2a_%d_tsub_ctm", lay), cluster->GetTimeSub(), cluster->GetCTime());
          }
          else if( cluster->hit(0)->cid()==CID_BLC2b ){
	    MyHistTools::fillTH(Form("BLC2b_%d_tsub_tmean", lay), cluster->GetTimeSub(), cluster->GetTimeMean());
	    MyHistTools::fillTH(Form("BLC2b_%d_tsub_ctm", lay), cluster->GetTimeSub(), cluster->GetCTime());
          }
        }
      }
    }
  }
}
