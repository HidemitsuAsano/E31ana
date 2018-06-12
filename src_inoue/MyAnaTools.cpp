# include "MyAnaTools.h"

using namespace std;

double MyAnaTools::T0BVC_tof(AnaInfo *anaInfo, BeamLineHitMan *blMan)
{
  if( anaInfo->nBeam()!=1 ) return -9999;
  BeamInfo *beam=anaInfo->beam(0);
  vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
  if( BVChits.size()==1 ){
    HodoscopeLikeHit *BVChit=BVChits[0];
    TVector3 T0pos = beam->T0pos();
    return (BVChit->ctmean()-beam->T0time());
  }

  return -9999;
}

TLorentzVector MyAnaTools::target_lmom()
{
  TLorentzVector lmom;
  string mat=DetectorList::GetInstance()->GetMaterial(CID_Fiducial);
  if( mat=="LHydrogen" ) lmom.SetVectM(TVector3(0, 0, 0), pMass);
  else if( mat=="LDeuterium" ) lmom.SetVectM(TVector3(0, 0, 0), dMass);
  else if( mat=="LHelium-3" ) lmom.SetVectM(TVector3(0, 0, 0), ThreeHeMass);
  return lmom;
}

int MyAnaTools::targetNA()
{
  string mat=DetectorList::GetInstance()->GetMaterial(CID_Fiducial);
  if( mat=="LHydrogen" ) return 1;
  else if( mat=="LDeuterium" ) return 2;
  else if( mat=="LHelium-3" ) return 3;
  return 0;
}

bool MyAnaTools::goodBeam(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return false;
  if( anaInfo->beam(0)->T0seg()<0 ) return false;
  if( !isTOFKaon(anaInfo->beam(0)->tof()) ) return false;
  if( !anaBLC1(anaInfo) ) return false;
  if( !anaBLC2(anaInfo) ) return false;
  if( !anaD5(anaInfo) ) return false;
  if( !connectD5BHD(anaInfo) ) return false;
  if( !anaBPC(anaInfo) ) return false;
  if( !connectBLC2BPC(anaInfo) ) return false;

  anaInfo->beam(0)->SetFlag(true);
  return true;
}

bool MyAnaTools::trigBLC1(AnaInfo *anaInfo, ConfMan *conf)
{
  if( anaInfo->nBeam()!=1 ) return false;
  if( anaInfo->beam(0)->T0seg()<0 ) return false;
  if( !isTOFKaon(anaInfo->beam(0)->tof()) ) return false;
  if( !anaBLC2(anaInfo) ) return false;
  if( !anaBPC(anaInfo) ) return false;
  if( !connectBLC2BPC(anaInfo) ) return false;
  if( !beamFiducial(anaInfo, conf) ) return false;

  return true;
}

bool MyAnaTools::trigBLC2(AnaInfo *anaInfo, ConfMan *conf)
{
  if( anaInfo->nBeam()!=1 ) return false;
  if( anaInfo->beam(0)->T0seg()<0 ) return false;
  if( !isTOFKaon(anaInfo->beam(0)->tof()) ) return false;
  if( !anaBLC1(anaInfo) ) return false;
  if( !anaBPC(anaInfo) ) return false;
  if( !beamFiducial(anaInfo, conf) ) return false;

  return true;
}

bool MyAnaTools::trigBPC(AnaInfo *anaInfo, ConfMan *conf, BeamLineHitMan *blMan)
{
  //   cout<<"trigBPC check "<<endl;
  bool DEFflag=false;
  std::vector<HodoscopeLikeHit*> DEFhits=MyTools::getHodo(blMan, CID_DEF);
  for( int i=0; i<DEFhits.size(); i++ ){
    if( 1<DEFhits[i]->seg() && DEFhits[i]->seg()<8 ) DEFflag=true;
  }
  if( !DEFflag ) return false;

  BeamInfo *beam=anaInfo->beam(0);
  if( anaInfo->nBeam()!=1 ) return false;
  if( anaInfo->beam(0)->T0seg()<0 ) return false;
  if( !isTOFKaon(anaInfo->beam(0)->tof()) ) return false;
  if( !anaBLC1(anaInfo) ) return false;
  int nBLC2=0;
  BLDCTrackInfo infoBLC2;
  for( int i=0; i<beam->nBLC2(); i++ ){
    BLDCTrackInfo tmpBLC2=beam->BLC2(i);
    if( BLC2_TIME_WINDOW_MIN<tmpBLC2.time() && tmpBLC2.time()<BLC2_TIME_WINDOW_MAX ){
      nBLC2++;
      infoBLC2=tmpBLC2;
    }
  }
  if( nBLC2!=1 ) return false;

  if( infoBLC2.chi2()>BLC2_CHI2_MAX ) return false;
  TVector3 pos=infoBLC2.GetPosatZ(-20.3);
  TVector3 dir=infoBLC2.dir();
  // cout<<"dir x : "<<dir.X()<<endl;
  // cout<<"dir y : "<<dir.Y()<<endl;

  double r=sqrt(pow(pos.X()+0.04, 2)+pow(pos.Y()-0.05, 2));
  if( fabs(r)>1.0 ) return false;

  if( !anaD5(anaInfo) ) return false;
  if( !connectD5BHD(anaInfo) ) return false;
  if( !anaFDC1(anaInfo) ) return false;
  double T0BVCtof=MyAnaTools::T0BVC_tof(anaInfo, blMan);
  if( T0BVCtof<5 || 20<T0BVCtof ) return false; 

  return true;
}

bool MyAnaTools::trigFDC1(AnaInfo *anaInfo, ConfMan *conf, BeamLineHitMan *blMan)
{
  if( !goodBeam(anaInfo) ) return false;
  if( !beamFiducial(anaInfo, conf) ) return false;
  vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
  bool BVCflag=false;
  for( int i=0; i<BVChits.size(); i++ ){
    double T0BVCtof=BVChits[i]->ctmean()-anaInfo->beam(0)->T0time();
    if( 1<BVChits[i]->seg() && BVChits[i]->seg()<8 ){
      if( 5<T0BVCtof && T0BVCtof<20 ) BVCflag=true;
    }
  }
  if( !BVCflag ) return false;

  return true;
};

bool MyAnaTools::beamFiducial(AnaInfo *anaInfo, ConfMan *conf)
{
  if( anaInfo->nBeam()!=1 ) return false;
  BeamInfo *beam=anaInfo->beam(0);

  int nBPC=0;
  BLDCTrackInfo infoBPC;
  for( int i=0; i<beam->nBPC(); i++ ){
    BLDCTrackInfo tmpBPC=beam->BPC(i);
    if( BPC_TIME_WINDOW_MIN<tmpBPC.time() && tmpBPC.time()<BPC_TIME_WINDOW_MAX ){
      nBPC++;
      infoBPC=tmpBPC;
    }
  }
  if( nBPC!=1 ) return false;

  TVector3 pos;
  conf-> GetGeomMapManager()->GetPos(CID_TarCell, 0, pos);

  if( GeomTools::GetID(infoBPC.GetPosatZ(pos.Z()))==CID_Fiducial ) return true;
  return false;
}

bool MyAnaTools::connectBLC2BPC(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return false;

  BeamInfo *beam=anaInfo->beam(0);
  int nBLC2=0;
  BLDCTrackInfo infoBLC2;
  for( int i=0; i<beam->nBLC2(); i++ ){
    BLDCTrackInfo tmpBLC2=beam->BLC2(i);
    if( BLC2_TIME_WINDOW_MIN<tmpBLC2.time() && tmpBLC2.time()<BLC2_TIME_WINDOW_MAX ){
      nBLC2++;
      infoBLC2=tmpBLC2;
    }
  }
  if( nBLC2!=1 ) return false;

  int nBPC=0;
  BLDCTrackInfo infoBPC;
  for( int i=0; i<beam->nBPC(); i++ ){
    BLDCTrackInfo tmpBPC=beam->BPC(i);
    if( BPC_TIME_WINDOW_MIN<tmpBPC.time() && tmpBPC.time()<BPC_TIME_WINDOW_MAX ){
      nBPC++;
      infoBPC=tmpBPC;
    }
  }
  if( nBPC!=1 ) return false;

  double z = 0.5*(-130-20.3);
  TVector3 dir_diff=infoBLC2.dir()-infoBPC.dir();
  TVector3 pos_diff=infoBLC2.GetPosatZ(z)-infoBPC.GetPosatZ(z);

  if( pos_diff.X()<BLC2BPC_X_MIN || BLC2BPC_X_MAX<pos_diff.X() ){ return false; }
  if( pos_diff.Y()<BLC2BPC_Y_MIN || BLC2BPC_Y_MAX<pos_diff.Y() ){ return false; }
  if( dir_diff.X()<BLC2BPC_dX_MIN || BLC2BPC_dX_MAX<dir_diff.X() ){ return false; }
  if( dir_diff.Y()<BLC2BPC_dY_MIN || BLC2BPC_dY_MAX<dir_diff.Y() ){ return false; }
 
  return true;
}

bool MyAnaTools::isTOFKaon(double tof){  return (BEAM_TOF_K_MIN<tof && tof<BEAM_TOF_K_MAX); };

bool MyAnaTools::anaD5(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return false;
  if( anaInfo->beam(0)->D5chi2()<-1 || D5_CHI2_MAX<anaInfo->beam(0)->D5chi2() ) return false;

  return true;
}

bool MyAnaTools::connectD5BHD(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return false;
  if( anaInfo->beam(0)->D5chi2()<-1 || D5_CHI2_MAX<anaInfo->beam(0)->D5chi2() ) return false;
  BeamInfo *beam=anaInfo->beam(0);
  double D5mom=beam->D5mom();
  for( int i=0; i<beam->nBHD(); i++ ){
    int BHDseg;
    double BHDtime;
    beam->getBHD(i, BHDseg, BHDtime);
    if( BHDseg<6 || 15<BHDseg ) continue;
    if( BHD_MATCH_MIN[BHDseg-6]<D5mom && D5mom<BHD_MATCH_MAX[BHDseg-6] ) return true;
  }

  return false;
}

bool MyAnaTools::anaBLC1(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return false;
  BeamInfo *beam=anaInfo->beam(0);
  int nBLC1=0;
  BLDCTrackInfo infoBLC1;
  for( int i=0; i<beam->nBLC1(); i++ ){
    BLDCTrackInfo tmpBLC1=beam->BLC1(i);
    if( BLC1_TIME_WINDOW_MIN<tmpBLC1.time() && tmpBLC1.time()<BLC1_TIME_WINDOW_MAX ){
      nBLC1++;
      infoBLC1=tmpBLC1;
    }
  }
  if( nBLC1!=1 ) return false;

  if( BLC1_TIME_MIN<infoBLC1.time() && infoBLC1.time()<BLC1_TIME_MAX ){
    if( infoBLC1.chi2()<BLC1_CHI2_MAX ) return true;
  }
  return false;
}

bool MyAnaTools::anaFDC1(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return false;
  BeamInfo *beam=anaInfo->beam(0);
  int nFDC1=0;
  BLDCTrackInfo infoFDC1;
  for( int i=0; i<beam->nFDC1(); i++ ){
    BLDCTrackInfo tmpFDC1=beam->FDC1(i);
    if( FDC1_TIME_WINDOW_MIN<tmpFDC1.time() && tmpFDC1.time()<FDC1_TIME_WINDOW_MAX ){
      nFDC1++;
      infoFDC1=tmpFDC1;
    }
  }
  if( nFDC1!=1 ) return false;

  if( FDC1_TIME_MIN<infoFDC1.time() && infoFDC1.time()<FDC1_TIME_MAX ){
    if( infoFDC1.chi2()<FDC1_CHI2_MAX ) return true;
  }
  return false;
}

bool MyAnaTools::anaBLC2(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return false;
  BeamInfo *beam=anaInfo->beam(0);
  int nBLC2=0;
  BLDCTrackInfo infoBLC2;
  for( int i=0; i<beam->nBLC2(); i++ ){
    BLDCTrackInfo tmpBLC2=beam->BLC2(i);
    if( BLC2_TIME_WINDOW_MIN<tmpBLC2.time() && tmpBLC2.time()<BLC2_TIME_WINDOW_MAX ){
      nBLC2++;
      infoBLC2=tmpBLC2;
    }
  }
  if( nBLC2!=1 ) return false;

  if( BLC2_TIME_MIN<infoBLC2.time() && infoBLC2.time()<BLC2_TIME_MAX ){
    if( infoBLC2.chi2()<BLC2_CHI2_MAX ) return true;
  }
  return false;
}

bool MyAnaTools::anaBPC(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return false;
  BeamInfo *beam=anaInfo->beam(0);
  int nBPC=0;
  BLDCTrackInfo infoBPC;
  for( int i=0; i<beam->nBPC(); i++ ){
    BLDCTrackInfo tmpBPC=beam->BPC(i);
    if( BPC_TIME_WINDOW_MIN<tmpBPC.time() && tmpBPC.time()<BPC_TIME_WINDOW_MAX ){
      nBPC++;
      infoBPC=tmpBPC;
    }
  }
  if( nBPC!=1 ) return false;

  if( BPC_TIME_MIN<infoBPC.time() && infoBPC.time()<BPC_TIME_MAX ){
    if( infoBPC.chi2()<BPC_CHI2_MAX ) return true;
  }
  return false;
}
