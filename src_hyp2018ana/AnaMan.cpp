#include "AnaMan.h"

static const double BLDC_TimeWindowMIN=-20;
static const double BLDC_TimeWindowMAX=20;
static const double BLDC_MAX_CHI=30;
static const double BEAM_MOM_MAX_CHI=30;
static const double Beam_K_MIN=27;
static const double Beam_K_MAX=31;
static const double Beam_PI_MIN=24;
static const double Beam_PI_MAX=27;
static const double NC_GAMMA_MIN=0.95;

AnaMan::AnaMan(ConfMan *conf)
{
  fBeamSpec = new BeamSpectrometer(conf);
}

void AnaMan::set(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, CDSHitMan *cdsMan)
{
  //*** Fill BL Counter ***//
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ) fBHD_hit.push_back(blMan->BHD(i));
  }

  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ) fT0_hit.push_back(blMan->T0(i));
  }

  for( int i=0; i<blMan->nBPD(); i++ ){
    if( blMan->BPD(i)->CheckRange() ) fBPD_hit.push_back(blMan->BPD(i));
  }

  for( int i=0; i<blMan->nDEF(); i++ ){
    if( blMan->DEF(i)->CheckRange() ) fDEF_hit.push_back(blMan->DEF(i));
  }

  for( int i=0; i<blMan->nBVC(); i++ ){
    if( blMan->BVC(i)->CheckRange() ) fBVC_hit.push_back(blMan->BVC(i));
  }

  for( int i=0; i<blMan->nCVC(); i++ ){
    if( blMan->CVC(i)->CheckRange() ) fCVC_hit.push_back(blMan->CVC(i));
  }

  for( int i=0; i<blMan->nPC(); i++ ){
    if( blMan->PC(i)->CheckRange() ) fPC_hit.push_back(blMan->PC(i));
  }

  for( int i=0; i<blMan->nBD(); i++ ){
    if( blMan->BD(i)->CheckRange() ) fBD_hit.push_back(blMan->BD(i));
  }

  for( int i=0; i<blMan->nLB(); i++ ){
    if( blMan->LB(i)->CheckRange() ) fLB_hit.push_back(blMan->LB(i));
  }
  for( int i=0; blMan->nNC(); i++ ){
    HodoscopeLikeHit *hit = blMan-> NC(i);
    int seg = hit-> seg();
    int lay = 1+(seg-1)/14;
    int seg2 = 1+(seg-1)%14;
    if( hit-> CheckRange() ){
      fNC_hit[lay-1].push_back(hit);
    }
  }
  for( int i=0; i<blMan->nWVC(); i++ ){
    if( 0<blMan->WVC(i)->tdcu() && blMan->WVC(i)->tdcu()<4000 ) fWVC_hit.push_back(blMan->WVC(i));
  }

  //*** Fill CDS Counter ***//
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    if( cdsMan->CDH(i)->CheckRange() ) fCDH_hit.push_back(cdsMan->CDH(i));
  }


  //*** Fill BLDC ***//
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    double track_time = bltrackMan->trackBLC1(i)->GetTrackTime();
    fTrackBLC1.push_back(bltrackMan->trackBLC1(i));
    if( -50<track_time && track_time<200 ) fTrackBLC1_1.push_back(bltrackMan->trackBLC1(i));
    if( -30<track_time && track_time<100 ) fTrackBLC1_2.push_back(bltrackMan->trackBLC1(i));
  }

  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    double track_time = bltrackMan->trackBLC2(i)->GetTrackTime();
    fTrackBLC2.push_back(bltrackMan->trackBLC2(i));
    if( -50<track_time && track_time<200 ) fTrackBLC2_1.push_back(bltrackMan->trackBLC2(i));
    if( -30<track_time && track_time<100 ) fTrackBLC2_2.push_back(bltrackMan->trackBLC2(i));
  }

  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    double track_time = bltrackMan->trackBPC(i)->GetTrackTime();
    fTrackBPC.push_back(bltrackMan->trackBPC(i));
    if( -50<track_time && track_time<200 ) fTrackBPC_1.push_back(bltrackMan->trackBPC(i));
    if( -30<track_time && track_time<100 ) fTrackBPC_2.push_back(bltrackMan->trackBPC(i));
  }

  //*** BHD-T0 TOF Analysis ***//
  if( fT0_hit.size()==1 ){
    if( fT0_hit[0]->seg()!=1 ){
      fT0time = fT0_hit[0]->ctmean();
      fStatus = 1;
      for( int i=0; i<fBHD_hit.size(); i++ ){
	if( 5<fBHD_hit[i]->seg() && fBHD_hit[i]->seg()<16 ){
	  double tof = fT0time-fBHD_hit[i]->ctmean();
	  if( fBeamPID!=Beam_Kaon ){
	    if( Beam_K_MIN<tof && tof<Beam_K_MAX ) fBeamPID==Beam_Kaon;
	  }
	  if( fBeamPID==Beam_Other ){
	    if( Beam_PI_MIN<tof && tof<Beam_PI_MAX ) fBeamPID==Beam_Pion;
	  }
	}
      }
      if( fBeamPID==Beam_Kaon ) fStatus=2;
    }
  }

  //*** Beam Momentum  Analysis ***//
  bool BLC1_flag= false;
  bool BLC2_flag=false;
  bool BPC_flag = false;
  if( fTrackBLC1_2.size()==1 ){
    double BLC1_tracktime = fTrackBLC1_2[0]-> GetTrackTime();
    double BLC1_chi2 = fTrackBLC1_2[0]-> chi2all();
    if( BLDC_TimeWindowMIN<BLC1_tracktime && BLC1_tracktime<BLDC_TimeWindowMAX && BLC1_chi2<BLDC_MAX_CHI ) BLC1_flag=true;
  }
  if( fTrackBLC2_2.size()==1 ){
    double BLC2_tracktime = fTrackBLC2_2[0]-> GetTrackTime();
    double BLC2_chi2 = fTrackBLC2_2[0]-> chi2all();
    if( BLDC_TimeWindowMIN<BLC2_tracktime && BLC2_tracktime<BLDC_TimeWindowMAX && BLC2_chi2<BLDC_MAX_CHI ) BLC2_flag=true;
  }
  if( fTrackBPC_2.size()==1 ){
    double BPC_tracktime = fTrackBPC_2[0]-> GetTrackTime();
    double BPC_chi2 = fTrackBPC_2[0]-> chi2all();
    if( BLDC_TimeWindowMIN<BPC_tracktime && BPC_tracktime<BLDC_TimeWindowMAX && BPC_chi2<BLDC_MAX_CHI ) BPC_flag=true;
  }
  if( fStatus==2 ){
    if( BLC1_flag ){
      fStatus=3;
      if( BLC2_flag ){
	fStatus=4;
      }
    }
  }
  if( BLC1_flag && BLC2_flag ){
    fBeamSpec-> TMinuitFit(fTrackBLC1_2[0], fTrackBLC2_2[0], conf);
    if( fBeamSpec-> chisquare()<BEAM_MOM_MAX_CHI ){
      fD5mom = fBeamSpec->mom();
      fStatus=5;
    }
  }
  if( fStatus=5 && BPC_flag ){
    fT0pos = fTrackBPC[0]-> GetPosatZ(-110.5);
    fStatus=6;
  }
}

void AnaMan::set(ConfMan *conf, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan)
{
  if( fStatus==6 && cdstrackMan->nGoodTrack()>0 ) fStatus=7;
  if( fBeamPID==Beam_Pion && fBeamPID==Beam_Kaon && fD5mom>0.1 && fTrackBPC_2.size()==1 ){
    for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
      CDS1Data cdsData(i);
      int status = cdsData.set(fBeamPID, fT0time, fD5mom, fTrackBPC_2[0], cdstrackMan->GoodTrack(i), cdsMan);
      if( status==0 ) fCDS1.push_back(cdsData);
    }

    for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
      for( int j=0; j<cdstrackMan->nGoodTrack(); j++ ){
	CDS2Data cdsData(i, j);
	int status = cdsData.set(fBeamPID, fD5mom, fTrackBPC_2[0], cdstrackMan->GoodTrack(i), cdstrackMan->GoodTrack(j));
	if( status==0 ) fCDS2.push_back(cdsData);
      }
    }

    double max_dis = DBL_MAX;
    for( int i=0; i<fCDS1.size(); i++ ){
      if( fCDS1[i].dis()<max_dis ){
	fVtxID = i;
	max_dis = fCDS1[i].dis();
      }
    }

    //*** search NC hit ***//
    if( fBVC_hit.size()==0 && fCVC_hit.size()==0 ){
      for( int lay=0; lay<8; lay++ ){
	if( fNC_hit[lay].size()>0 ){
	  double NCtime = DBL_MAX;
	  HodoscopeLikeHit *hit = 0;
	  for( int i=0; i<fNC_hit[lay].size(); i++ ){
	    if( fNC_hit[lay][i]->ctmean()<NCtime ){
	      NCtime = fNC_hit[lay][i]->ctmean();
	      hit = fNC_hit[lay][i];
	    }
	  }

	  fFdE = hit-> emean();
	  fFtime = hit->ctmean();
	  conf-> GetGeomMapManager()-> GetGPos(CID_NC, hit->seg(), fFHitPos);
	  if( fVtxID>=0 ){
	    TVector3 vtx = fCDS1[fVtxID].vtxBeam();
	    double tof = fFtime-fT0time;
	    double fl = (vtx-fFHitPos).Mag();
	    fFBeta = fl/(tof*100.*Const);
	    if( fFBeta>NC_GAMMA_MIN ) fFPID=F_Gamma;
	    else{
	      fFPID=F_Neutron;
	      double mom = nMass*fFBeta/sqrt(1-fFBeta*fFBeta);
	      fFMom = (fFHitPos-vtx).Unit();
	      fFMom.SetMag(mom);
	    }
	  }
	  break;
	}
      }
    }
  }
}

void AnaMan::clear()
{
  fBeamSpec-> Clear();

  fStatus=0;
  fT0time = DBL_MIN;
  fT0pos = DEFVECT;
  fBeamPID=Beam_Other;
  fD5mom = DBL_MIN;

  fBHD_hit.clear();
  fT0_hit.clear();
  fBPD_hit.clear();
  fDEF_hit.clear();
  fCDH_hit.clear();
  fBVC_hit.clear();
  fCVC_hit.clear();
  fPC_hit.clear();
  fBD_hit.clear();
  fLB_hit.clear();
  fWVC_hit.clear();
  for( int i=0; i<8; i++ ) fNC_hit[i].clear();

  fTrackBLC1.clear();
  fTrackBLC2.clear();
  fTrackBPC.clear();

  fTrackBLC1_1.clear();
  fTrackBLC2_1.clear();
  fTrackBPC_1.clear();

  fTrackBLC1_2.clear();
  fTrackBLC2_2.clear();
  fTrackBPC_2.clear();

  fVtxID=-1;
  fCDS1.clear();
  fCDS2.clear();

  fFPID = F_Other;
  fFMom = DEFVECT;
  fFHitPos = DEFVECT;
  fFdE = DBL_MIN;
  fFBeta = DBL_MIN;
  fFtime = DBL_MIN;
}
