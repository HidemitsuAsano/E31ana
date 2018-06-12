#include "AnaData.h"

static const double BEAM_PI_MIN = 24.5;
static const double BEAM_PI_MAX = 27.25;
static const double BEAM_K_MIN = 27.25;
static const double BEAM_K_MAX = 30.0;
static const double BEAM_P_MIN = 30.0;
static const double BEAM_P_MAX = 35.0;

AnaData::AnaData()
{
  clear();
}

bool AnaData::set(EventHeader *header, ConfMan *conf, const double &D5, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, CDSHitMan *cdsMan,  CDSTrackingMan *cdstrackMan)
{
  fEventNumber = header-> ev();
  int nT0=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
      fT0time = blMan->T0(i)->ctmean();
    }
  }
  if( nT0!=0 ) return false;

  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      double tof = fT0time-blMan->BHD(i)->ctmean();
      if( fBeamPID==Beam_Other ){
	if( BEAM_P_MIN<tof  && tof<BEAM_P_MAX  ) fBeamPID==Beam_Proton;
	if( BEAM_PI_MIN<tof && tof<BEAM_PI_MAX ) fBeamPID==Beam_Kaon;
      }

      if( fBeamPID!=Beam_Kaon ){
	if( BEAM_K_MIN<tof  && tof<BEAM_K_MAX  ) fBeamPID==Beam_Kaon;
      }
    }
  }
  if( fBeamPID==Beam_Other ) return false;

  if( bltrackMan->ntrackBLC1()!=1 && bltrackMan->ntrackBLC2()!=1 ) return false;

  fD5mom=D5;
  if( bltrackMan->ntrackBPC()!=1 ) return false;
  LocalTrack *trackBPC = bltrackMan->trackBPC(0);
  TVector3 T0pos = trackBPC-> GetPosatZ(-110.5);

  double dis=DBL_MAX;
  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *goodtrack = cdstrackMan-> GoodTrack(i);
    CDS1Data cdsData(i);
    if( cdsData.set(fBeamPID, fT0time, fD5mom, trackBPC, goodtrack, cdsMan) ){
      fCDS1Data.push_back(cdsData);
      if( cdsData.dis()<dis ){
	fVtxFlag = true;
	dis = cdsData.dis();
	fVtx = cdsData.vtxBeam();
      }
    }
  }
  if( !fVtxFlag ) return false;
  else{
    TVector3 beam_p = trackBPC-> GetMomDir();
    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(T0pos, fVtx, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
    beam_p.SetMag(beam_out);
    fBeamLmom.SetVectM(beam_p, parMass[fBeamPID]);
  }


  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    for( int j=i+1; cdstrackMan->nGoodTrack(); j++ ){
      CDSTrack *cds1 = cdstrackMan->GoodTrack(i);
      CDSTrack *cds2 = cdstrackMan->GoodTrack(j);

      CDS2Data cdsData(i, j);
      if( cdsData.set(fBeamPID, fD5mom, trackBPC, cds1, cds2) ){
	fCDS2Data.push_back(cdsData);
      }
    }
  }

  std::vector<HodoscopeLikeHit*> CVChits;
  std::vector<HodoscopeLikeHit*> PChits;
  bool NC_hit_flag = false;
  std::vector<HodoscopeLikeHit*> NChits[8];
  for( int i=0; i<blMan->nCVC(); i++ ){
    if( blMan->CVC(i)->CheckRange() ) CVChits.push_back(blMan->CVC(i));
  }
  for( int i=0; i<blMan->nPC(); i++ ){
    if( blMan->PC(i)->CheckRange() ) PChits.push_back(blMan->PC(i));
  }
  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->CheckRange() ){
      NC_hit_flag = true;
      int seg = blMan->NC(i)->seg();
      int lay = (seg-1)/14;
      NChits[lay].push_back(blMan->NC(i));
    }
  }

  if( CVChits.size()==0 && NC_hit_flag ){
    for( int lay=0; lay<8; lay++ ){
      if( NChits[lay].size()!=0 ){
	double NCtime = DBL_MAX;
	int NCseg = -1;
	for( int i=0; i<NChits[lay].size(); i++ ){
	  if( NChits[lay][i]->ctmean()<NCtime ){
	    NCtime = NChits[lay][i]->ctmean();
	    NCseg = NChits[lay][i]-> seg();
	    fFCID = CID_NC;
	    fFseg = NChits[lay][i]->seg();
	    fFdE = NChits[lay][i]->emean();
	  }
	}
	break;
      }
    }
  }


  return true;
}

CDS1Data* AnaData::CDS1(const int &pid, const int &n)
{
  int index=0;
  for( int i=0; i<fCDS1Data.size(); i++ ){
    if( pid==fCDS1Data[i].pid() ){
      if( n==index ) return &fCDS1Data[i];
      index++;
    }
  }

  std::cout<<"  !!! AnaData::CDS1("<<pid<<", "<<n<<") not found !!!"<<std::endl;
  return 0;
}

CDS2Data* AnaData::CDS2(const int &pid1, const int &pid2, const int &n)
{
  int index=0;
  for( int i=0; i<fCDS2Data.size(); i++ ){
    if( pid1==fCDS2Data[i].pid1() && pid2==fCDS2Data[i].pid2() ){
      if( n==index ) return &fCDS2Data[i];
      index++;
    }
    else if( pid2==fCDS2Data[i].pid1() && pid1==fCDS2Data[i].pid2() ){
      if( n==index ) return &fCDS2Data[i];
      index++;
    }
  }

  std::cout<<"  !!! AnaData::CDS2("<<pid1<<", "<<pid2<<", "<<n<<") not found !!!"<<std::endl;
  return 0;
}

void AnaData::clear()
{
  fEventNumber = -1;
  fT0time = DBL_MIN;
  fD5mom = 0.0;
  fBeamPID = Beam_Other;
  fBeamLmom.SetXYZT(0., 0., 0., 0.);
  fVtxFlag = false;
  fVtx.SetXYZ(DBL_MIN, DBL_MIN, DBL_MIN);

  fCDS1Data.clear();
  fCDS2Data.clear();

  fForwardPID = F_Other;
  fFCID=-1;
  fFseg=-1;
  fFdE=0.0;
  fFBeta=0.0;
  fForwardLmom.SetXYZT(0., 0., 0., 0.);
}
