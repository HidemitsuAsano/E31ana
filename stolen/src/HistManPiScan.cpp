#include "HistManBeamAna.h"

static const double pi_v = 0.992046*Const*m; // 1.1GeV/c
//static const double pi_v = 0.9904*Const*m;   // 1.0GeV/c
//static const double pi_v = 0.988188*Const*m; // 0.9GeV/c
static const double USWK_Z = 250; //cm


void HistManBeamAna::initPiScan()
{
  fFile->cd();
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("T0CVC%d_offset", seg), Form("CVC seg%d offset by T0-CVC", seg), 20000, -100, 100);
  }
  for( int seg=1; seg<=27; seg++ ){
    new TH1F(Form("T0PC%d_offset", seg), Form("PC seg%d offset by T0-PC", seg), 20000, -100, 100);
  }
  for( int seg=1; seg<=112; seg++ ){
    new TH1F(Form("T0NC%d_offset", seg), Form("NC seg%d offset by T0-NC", seg), 20000, -100, 100);
  }
}

void HistManBeamAna::fillPiScan(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan)
{
  TH1 *h1;
  TH2 *h2;

  std::cout<<"start"<<std::endl;
  std::cout<<"BHD mul : "<<fBHD_hit.size()<<std::endl;
  std::cout<<"T0 mul : "<<fT0_hit.size()<<std::endl;
  std::cout<<"BLC2 track : "<<bltrackMan->ntrackBLC2()<<std::endl;
  std::cout<<"BPC track  : "<<bltrackMan->ntrackBPC()<<std::endl;
  std::cout<<"FDC1 track : "<<bltrackMan->ntrackFDC1()<<std::endl;

  if( fT0_hit.size()==1 && bltrackMan->ntrackFDC1()==1 && bltrackMan->ntrackBPC()==1 ){
    std::cout<<"aaaaa"<<std::endl;
    //    LocalTrack *trackBLC2 = bltrackMan->trackBLC2(0);
    LocalTrack *trackBPC  = bltrackMan->trackBPC(0);
    LocalTrack *trackFDC  = bltrackMan->trackFDC1(0);
    TVector3 gBPCpos, gFDCpos, gBPCdir, gFDCdir;
    conf-> GetBLDCWireMapManager()-> GetGParam(CID_FDC1, gFDCpos, gFDCdir);
    conf-> GetBLDCWireMapManager()-> GetGParam(CID_BPC, gBPCpos, gBPCdir);
    TVector3 pos;
    conf->GetGeomMapManager()-> GetPos(CID_T0, fT0_hit[0]->seg(), pos);
    //    TVector3 T0pos = trackBLC2-> GetPosatZ(pos.Z()-0.5);
    TVector3 T0pos = trackBPC-> GetPosatZ(pos.Z()-0.5);
    TVector3 FDCpos = trackFDC-> GetPosatZ(gFDCpos.Z());
    TVector3 BPCpos = trackBPC-> GetPosatZ(gBPCpos.Z());
    TVector3 dir = (FDCpos-BPCpos).Unit();
    TVector3 USWKpos = FDCpos+(250-FDCpos.Z())*dir;

    if( fCVC_hit.size()+fPC_hit.size()==1 ){
      HodoscopeLikeHit *charge_hit =0;
      TVector3 charge_pos;
      double tof;
      if( fCVC_hit.size()==1 ){
	tof = fCVC_hit[0]->ctmean()-fT0_hit[0]->ctmean();
	conf-> GetGeomMapManager()-> GetGPos(CID_CVC, fCVC_hit[0]->seg(), charge_pos);
      }
      else{
	tof = fPC_hit[0]->ctmean()-fT0_hit[0]->ctmean();
	conf-> GetGeomMapManager()-> GetGPos(CID_PC, fPC_hit[0]->seg(), charge_pos);
      }
      double bend_angle = (USWKpos-BPCpos).Angle(charge_pos-USWKpos);

      double leff = 100;
      double r = leff/sqrt(2*(1-cos(bend_angle)));
      double larc = r*bend_angle;
      double line = 2*leff/sqrt(2*(1+cos(bend_angle)));
      double diff = line-larc;
      double fl = (T0pos-BPCpos).Mag()+(BPCpos-USWKpos).Mag()+(USWKpos-charge_pos).Mag()-diff;
      fl -= 1.5;
      double offset = tof - fl/pi_v;

      if( fCVC_hit.size()==1 ) h1 = getTH(Form("T0CVC%d_offset", fCVC_hit[0]->seg())), h1->Fill(offset);
      else  h1 = getTH(Form("T0PC%d_offset", fPC_hit[0]->seg())), h1->Fill(offset);

      for( int i=0; i<blMan->nNC(); i++ ){
	if( blMan->NC(i)->CheckRange() ){
	  HodoscopeLikeHit *nc_hit = blMan->NC(i);
	  tof = nc_hit->ctmean()-fT0_hit[0]->ctmean();
	  TVector3 nc_pos;
	  conf-> GetGeomMapManager()-> GetGPos(CID_NC, nc_hit->seg(), nc_pos);
	  bend_angle = (USWKpos-BPCpos).Angle(nc_pos-USWKpos);
	  r = leff/sqrt(2*(1-cos(bend_angle)));
	  larc = r*bend_angle;
	  line = 2*leff/sqrt(2*(1+cos(bend_angle)));
	  diff = line-larc;
	  fl = (T0pos-BPCpos).Mag()+(BPCpos-USWKpos).Mag()+(USWKpos-nc_pos).Mag()-diff;
	  fl -= 2.5;
	  double offset = tof - fl/pi_v;	 
	  h1 = getTH(Form("T0NC%d_offset", nc_hit->seg())), h1-> Fill(offset);
	}
      }
    }
  }
}
