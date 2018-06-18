#include "BeamInfo.h"

using namespace std;

BeamInfo::BeamInfo() : 
  fFlag(false), fPID(Beam_Other), fBHDseg(DEFAULTI), fBHDtime(DEFAULTD), fT0seg(DEFAULTI), fT0time(DEFAULTD), 
  fD5mom(DEFAULTD), fD5chisquare(DEFAULTD), fBLC1id(DEFAULTI), fBLC2id(DEFAULTI), fBPCid(DEFAULTI),
  fT0pos(DEFVECT), fVertex(DEFVECT), fVertexMom(DEFAULTD)
{
}

TLorentzVector BeamInfo::lmom() const
{
  TLorentzVector lmom;
  TVector3 mom=fBPCinfo[fBPCid].dir();
  mom.SetMag(fVertexMom);
  lmom.SetVectM(mom, parMass[fPID]);
  return lmom;
}

TLorentzVector BeamInfo::lmom(const TVector3 &pos) const
{
  TLorentzVector lmom;
  double beam_out, tmp_tof;
  ELossTools::CalcElossBeamTGeo(fT0pos, pos, fD5mom, parMass[fPID], beam_out, tmp_tof);
  TVector3 mom=fBPCinfo[fBPCid].dir();
  mom.SetMag(beam_out);
  lmom.SetVectM(mom, parMass[fPID]);
  return lmom;
}

HodoscopeLikeHit* BeamInfo::T0(BeamLineHitMan *blMan) const
{
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->seg()==fT0seg ) return blMan->T0(i);
  }
  cout<<"  !!! T0 seg"<<fT0seg<<" not found"<<endl;
  return 0;
}

HodoscopeLikeHit* BeamInfo::BHD(BeamLineHitMan *blMan) const
{
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->seg()==fBHDseg ) return blMan->BHD(i);
  }
  return 0;
}

void BeamInfo::SetBPCid(const int &id, ConfMan *conf)
{
  fBPCid=id;
  double T0z=-110.5;
  if( conf && fT0seg>0 ){
    TVector3 gpos;
    conf-> GetGeomMapManager()-> GetGPos(CID_T0, fT0seg, gpos);
    T0z=gpos.z()-0.5;
  }
  fT0pos=fBPCinfo[fBPCid].GetPosatZ(T0z);
}

void BeamInfo::SetBLDC(BeamLineTrackMan *bltrackMan)
{
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    fBLC1info.push_back(BLDCTrackInfo(i, bltrackMan->trackBLC1(i)));
  }
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    fBLC2info.push_back(BLDCTrackInfo(i, bltrackMan->trackBLC2(i)));
  }
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    fBPCinfo.push_back(BLDCTrackInfo(i, bltrackMan->trackBPC(i)));
  }
  for( int i=0; i<bltrackMan->ntrackFDC1(); i++ ){
    fFDC1info.push_back(BLDCTrackInfo(i, bltrackMan->trackFDC1(i)));
  }
}

bool BeamInfo::getBHD(const int &n, int &seg, double &time)
{
  if( n<0 || (int)fBHDtimes.size()<n ) return false;
  std::map<int, double>::iterator ite=fBHDtimes.begin();
  for( int i=0; i<n; i++ ) ++ite;
  seg=ite->first;
  time=ite->second;
  return true;
}

void BeamInfo::dump() const
{
  cout<<"T0 seg"<<fT0seg<<" T0 time : "<<fT0time<<endl;
  cout<<"BHD seg"<<fBHDseg<<"  BHD time : "<<fBHDtime<<endl;
  cout<<"BHD-T0 tof : "<<fT0time-fBHDtime<<" [ns]  trigger : ";
  if( fPID==Beam_Kaon        ) cout<<"Kaon"<<endl;
  else if( fPID==Beam_Pion   ) cout<<"Pion"<<endl;
  else if( fPID==Beam_Proton ) cout<<"Proton"<<endl;
  else if( fPID==Beam_Proton ) cout<<"Deuteron"<<endl;
  else if( fPID==Beam_Other  ) cout<<"Other"<<endl;

  cout<<"D5 mom : "<<fD5mom<<"[GeV/c]"<<endl;
}

double BeamInfo::calc_tof() const
{
  const double fl=770;
  const double velocity=100.*Const*fD5mom/(parMass[fPID]*sqrt(1+fD5mom*fD5mom/(parMass[fPID]*parMass[fPID])));
  return fl/velocity;
}

void BeamInfo::SetVertex(const TVector3 &pos)
{
  fVertex=pos;
  double beam_out, tmp_tof;
  ELossTools::CalcElossBeamTGeo(fT0pos, pos, fD5mom, parMass[fPID], beam_out, tmp_tof);
  fVertexMom=beam_out;
}

void BeamInfo::SetBLMan(BeamLineHitMan *blMan)
{
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      HodoscopeLikeHit *hit=blMan->BHD(i);
      fBHDtimes[hit->seg()]=hit->ctmean();
    }
  }
}

BLDCTrackInfo::BLDCTrackInfo()
{
}

BLDCTrackInfo::BLDCTrackInfo(int id, LocalTrack *track) :
  fID(id), fTime(track->GetTrackTime()), fChi2(track->chi2all()), 
  fPos(track->GetPosatZ(track->hit(0)->gz())),  fDir(track->GetMomDir())
{
}
