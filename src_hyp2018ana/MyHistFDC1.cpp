#include "MyHistFDC1.h"

using namespace std;

void initHistFDC1()
{
  const int nlay=6;
  const int nw=64;

  for( int lay=1; lay<=nlay; lay++ ){
    new TH1F(Form("hitpatFDC1_%d", lay), Form("FDC1 hit pattern lay%d", lay), nw, 0.5, nw+0.5);
    new TH1F(Form("mulFDC1_%d", lay), Form("FDC1 multiplicity lay%d", lay), nw+1, -0.5, nw+0.5);

    new TH1F(Form("FDC1_TDC_%d", lay), Form("FDC1 TDC lay%d", lay), 2000, 0, 2000);
    new TH1F(Form("FDC1_dt_%d", lay), Form("FDC1 dt lay%d", lay), 600, -200, 400);

    new TH1F(Form("FDC1_TDC_%d_1hit", lay), Form("FDC1 TDC lay%d", lay), 2000, 0, 2000);
    new TH1F(Form("FDC1_dt_%d_1hit", lay), Form("FDC1 dt lay%d", lay), 600, -200, 400);
    new TH1F(Form("FDC1_dt_%d_track", lay), Form("FDC1 dt lay%d", lay), 600, -200, 400);
    new TH1F(Form("FDC1_dt_%d_track_wBVC", lay), Form("FDC1 dt lay%d", lay), 500, -100, 400);

    new TH1F(Form("FDC1_dt_%d_select", lay), Form("FDC1 dt lay%d", lay), 500, -100, 400);

    new TH1F(Form("FDC1_res_%d", lay), Form("FDC1 res lay%d", lay), 500, -0.1, 0.1);
    for( int w=1; w<=nw; w++ ){
      new TH1F(Form("FDC1_TDC_%d_%d", lay, w), Form("FDC1 TDC lay%d wire%d", lay, w), 2000, 0, 2000);
      new TH1F(Form("FDC1_dt_%d_%d", lay, w), Form("FDC1 dt lay%d wire%d", lay, w), 500, -100, 400);

      new TH1F(Form("FDC1_TDC_%d_%d_1hit", lay, w), Form("FDC1 TDC lay%d wire%d", lay, w), 2000, 0, 2000);
      new TH1F(Form("FDC1_dt_%d_%d_1hit", lay, w), Form("FDC1 dt lay%d wire%d", lay, w), 500, -100, 400);
      new TH1F(Form("FDC1_dt_%d_%d_track", lay, w), Form("FDC1 dt lay%d wire%d", lay, w), 500, -100, 400);
      new TH1F(Form("FDC1_dt_%d_%d_track_wBVC", lay, w), Form("FDC1 dt lay%d wire%d", lay, w), 500, -100, 400);

      new TH1F(Form("FDC1_dt_%d_%d_select", lay, w), Form("FDC1 dt lay%d wire%d", lay, w), 500, -100, 400);

      new TH1F(Form("FDC1_res_%d_%d", lay, w), Form("FDC1 res lay%d wire%d", lay, w), 500, -0.1, 0.1);
    }
    if( lay%2==1 ){
      new TH2F(Form("FDC1_%d_tsub_tmean", lay), 
	       Form("FDC1 lay%d-%d time sub vs time mean", lay, lay+1), 500, 
	       -250, 250, 500, -100, 400);
      new TH2F(Form("FDC1_%d_tsub_ctm", lay), 
	       Form("FDC1 lay%d-%d time sub vs corr. time mean", lay, lay+1), 
	       500, -250, 250, 500, -100, 400);

      new TH2F(Form("FDC1_%d_tsub_tmean_select", lay), 
	       Form("FDC1 lay%d-%d time sub vs time mean", lay, lay+1), 
	       500, -250, 250, 500, -100, 400);
      new TH2F(Form("FDC1_%d_tsub_ctm_select", lay), 
	       Form("FDC1 lay%d-%d time sub vs corr. time mean", lay, lay+1), 
	       500, -250, 250, 500, -100, 400);
    }
  }


  new TH1F("ntrackFDC1_1", "FDC1 ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackFDC1_2", "FDC1 ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackFDC1", "FDC1 ntrack", 10, -0.5, 9.5);
  new TH1F("FDC1_time", "FDC1 time", 500, -200, 300);
  new TH1F("FDC1_chi2", "FDC1 chi-square", 500, 0.0, 500);
}

void fillFDC1(LocalTrack *track)
{
  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit=track->hit(i);
    int lay=hit->layer(), wire=hit->wire(), tdc=hit->tdc();
    double dt=hit->dt();
    double res=hit->resl();

    //    cout<<" FDC1 l"<<lay<<" w"<<wire<<"  resi="<<res<<endl;
    MyHistTools::fillTH(Form("FDC1_res_%d", lay), res);
    MyHistTools::fillTH(Form("FDC1_res_%d_%d", lay, wire), res);
  }
}


void fillFDC1(BeamLineHitMan *blMan)
{
  const int nlay=6;
  bool all_1hit=true;
  for( int lay=1; lay<=nlay; lay++ ){
    if( blMan->nFDC1(lay)!=1 ) all_1hit=false;
    MyHistTools::fillTH(Form("mulFDC1_%d", lay), blMan->nFDC1(lay));

    for( int i=0; i<blMan->nFDC1(lay); i++ ){
      ChamberLikeHit *hit=blMan->FDC1(lay, i);
      int wire=hit->wire();
      int tdc=hit->tdc();
      double dt=hit->dt();

      MyHistTools::fillTH(Form("hitpatFDC1_%d", lay), wire);
      MyHistTools::fillTH(Form("FDC1_TDC_%d", lay), tdc);
      MyHistTools::fillTH(Form("FDC1_dt_%d", lay), dt);
      MyHistTools::fillTH(Form("FDC1_TDC_%d_%d", lay, wire), tdc);
      MyHistTools::fillTH(Form("FDC1_dt_%d_%d", lay, wire), dt);
    }
  }

  if( all_1hit ){
    for( int lay=1; lay<=nlay; lay++ ){
      ChamberLikeHit *hit=blMan->FDC1(lay, 0);
      int wire=hit->wire();
      int tdc=hit->tdc();
      double dt=hit->dt();

      MyHistTools::fillTH(Form("FDC1_TDC_%d_1hit", lay), tdc);
      MyHistTools::fillTH(Form("FDC1_dt_%d_1hit", lay), dt);
      MyHistTools::fillTH(Form("FDC1_TDC_%d_%d_1hit", lay, wire), tdc);
      MyHistTools::fillTH(Form("FDC1_dt_%d_%d_1hit", lay, wire), dt);
    }
  }
}

void fillFDC1(EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan)
{
  int nFDC1_0=0, nFDC1_1=0, nFDC1_2=0;
  LocalTrack *trackFDC1=0;
  for( int i=0; i<bltrackMan->ntrackFDC1(); i++ ){
    double time=bltrackMan->trackFDC1(i)->GetTrackTime();
    MyHistTools::fillTH("FDC1_time", time);
    if( -10<time && time<10  ){
      nFDC1_0++;
      trackFDC1=bltrackMan->trackFDC1(i);
    }
    if( -30<time && time<100 ) nFDC1_1++;
    if( -50<time && time<200 ) nFDC1_2++;
  }

  MyHistTools::fillTH("ntrackFDC1", bltrackMan->ntrackFDC1());
  MyHistTools::fillTH("ntrackFDC1_1", nFDC1_1);
  MyHistTools::fillTH("ntrackFDC1_2", nFDC1_2);

  if( nFDC1_0==1 && nFDC1_1==1 ){
    MyHistTools::fillTH("FDC1_chi2", trackFDC1->chi2all());
  }
  if( bltrackMan->ntrackFDC1()==1 ){
    LocalTrack *track = bltrackMan->trackFDC1(0);

    if( 25<track->GetTrackTime() || track->GetTrackTime()<40 ){
      if( track->chi2all()<10 ){
	vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
	if( BVChits.size()==1 ){
	  for( int i=0; i<track->nhit(); i++ ){
	    ChamberLikeHit *hit=track->hit(i);
	    MyHistTools::fillTH(Form("FDC1_dt_%d_track_wBVC", hit->layer()), hit->dt());
	    MyHistTools::fillTH(Form("FDC1_dt_%d_%d_track_wBVC", hit->layer(), hit->wire()), hit->dt());
	  }
	  if( header->IsTrig(Trig_Kf) || header->IsTrig(Trig_Beam) ){
	    for( int i=0; i<track->nhit(); i++ ){
	      ChamberLikeHit *hit=track->hit(i);
	      MyHistTools::fillTH(Form("FDC1_dt_%d_select", hit->layer()), hit->dt());
	      MyHistTools::fillTH(Form("FDC1_dt_%d_%d_select", hit->layer(), hit->wire()), hit->dt());
	    }

	    for( int xy=0; xy<1; xy++ ){
	      for( int i=0; i<track->ncluster(xy); i++){
		BLDCCluster *cluster=track->cluster(xy, i);
		
		if( cluster->nhit()!=2 ) continue;
		int lay=cluster->hit(0)->layer();
		if( lay%2==0 ) lay=cluster->hit(1)->layer();
		if( lay%2==0 ) continue;
		
		MyHistTools::fillTH(Form("FDC1_%d_tsub_tmean_select", lay), cluster->GetTimeSub(), cluster->GetTimeMean());
		MyHistTools::fillTH(Form("FDC1_%d_tsub_ctm_select", lay), cluster->GetTimeSub(), cluster->GetCTime());
	      }
	    }
	  }
	}

	for( int i=0; i<track->nhit(); i++ ){
	  ChamberLikeHit *hit=track->hit(i);
	  MyHistTools::fillTH(Form("FDC1_dt_%d_track", hit->layer()), hit->dt());
	  MyHistTools::fillTH(Form("FDC1_dt_%d_%d_track", hit->layer(), hit->wire()), hit->dt());
	}
	
	for( int xy=0; xy<1; xy++ ){
	  for( int i=0; i<track->ncluster(xy); i++){
	    BLDCCluster *cluster=track->cluster(xy, i);
	    
	    if( cluster->nhit()!=2 ) continue;
	    int lay=cluster->hit(0)->layer();
	    if( lay%2==0 ) lay=cluster->hit(1)->layer();
	    if( lay%2==0 ) continue;
	    
	    MyHistTools::fillTH(Form("FDC1_%d_tsub_tmean", lay), cluster->GetTimeSub(), cluster->GetTimeMean());
	    MyHistTools::fillTH(Form("FDC1_%d_tsub_ctm", lay), cluster->GetTimeSub(), cluster->GetCTime());
	  }
	}
      }
    }
  }
}
