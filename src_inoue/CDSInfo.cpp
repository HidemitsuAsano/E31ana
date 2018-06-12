#include "CDSInfo.h"

using namespace std;

CDSInfo::CDSInfo() :
  fFlag(false), fTrackID(-1), fPID(CDS_DEFAULT), fCDHseg(-1), fIHseg(-1), 
  fMom(DEFAULTD), fMass2(DEFAULTD), fBeta(DEFAULTD), fOffset(DEFAULTD),
  fVertexCDS(DEFVECT), fVertexBeam(DEFVECT), fMomentum(DEFVECT)
{
}

HodoscopeLikeHit *CDSInfo::CDH(CDSHitMan *cdsMan) const
{
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    if( cdsMan->CDH(i)->seg()==fCDHseg ) return cdsMan->CDH(i);
  }
  return 0;
}

CDSTrack *CDSInfo::track(ConfMan *conf, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan) const
{
  CDSTrack *track=cdstrackMan->Track(fTrackID);
  bool calc_flag=false;
  for( int lay=1; lay<=15; lay++ ){
    if( track->nTrackHit(lay)==1 ){
      CDCHit *hit=track->hit(cdsMan, lay, 0);
      if( hit->resl()<-998 ) calc_flag=true;
      break;
    }
  }
  if( calc_flag ){
    if( track->FittingLevel()==5 )  track->Retiming(cdsMan, conf, beta(), true);
    else if( track->FittingLevel()==4 ){
      if( conf->GetSlewingMapManager()->isParam(CID_CDC, 1, 0) || conf->GetSlewingMapManager()->isParam(CID_CDC, 1, 1) ){
	track->Retiming(cdsMan, conf, beta(), true);
      }
      track->Retiming(cdsMan, conf, beta(), false);
    }
    else track->Retiming(cdsMan, conf, beta(), false);
    track-> SetHitPos(cdsMan);
  }

  return track;
}

HodoscopeLikeHit *CDSInfo::IH(CDSHitMan *cdsMan) const
{
  for( int i=0; i<cdsMan->nIH(); i++ ){
    if( cdsMan->IH(i)->seg()==fCDHseg ) return cdsMan->IH(i);
  }
  return 0;
}

TLorentzVector CDSInfo::lmom() const
{
  TLorentzVector lmom;
  lmom.SetVectM(fMomentum, pdgMass());
  return lmom;
}

void CDSInfo::dump() const
{
  cout<<"===== CDSInfo::dump ====="<<endl;
  cout<<"> Flag : "<<boolalpha<<fFlag<<"  TrackID : "<<fTrackID<<endl;
  if( fPID==CDS_PiPlus   ) cout<<"> Pion+"<<endl;
  else if( fPID==CDS_Proton   ) cout<<"> Proton"<<endl;
  else if( fPID==CDS_Deuteron ) cout<<"> Deutern"<<endl;
  else if( fPID==CDS_Helium3  ) cout<<"> He3"<<endl;
  else if( fPID==CDS_Triton   ) cout<<"> Triton"<<endl;
  else if( fPID==CDS_PiMinus  ) cout<<"> Pion-"<<endl;
  else if( fPID==CDS_Kaon     ) cout<<"> Kaon-"<<endl;
  else if( fPID==CDS_Other    ) cout<<"> Other"<<endl;
  else if( fPID==CDS_DEFAULT  ) cout<<"> DEFAULT"<<endl;
  else cout<<"> unknown"<<endl;

  cout<<"> mass2 : "<<fMass2<<"[(GeV/c2)2]  mom : "<<fMom<<"[GeV/c]  beta : "<<fBeta<<"  offset : "<<fOffset<<endl;
  cout<<"> Vtx Beam "<<Form("(%lf, %lf, %lf) [cm]", fVertexBeam.X(), fVertexBeam.Y(), fVertexBeam.Z())<<endl;
  cout<<"> Vtx CDS  "<<Form("(%lf, %lf, %lf) [cm]", fVertexCDS.X(), fVertexCDS.Y(), fVertexCDS.Z())<<endl;
  cout<<"> DCA : "<<dca()<<" [cm]"<<endl;
  cout<<"> Momentum "<<Form("(%lf, %lf, %lf) [cm]", fMomentum.X(), fMomentum.Y(), fMomentum.Z())<<endl;
  if( fCDHseg>0 ) cout<<"> CDH seg"<<fCDHseg<<"  "<<endl;
  else            cout<<"> ! CDH no hit  "<<endl;

  if( fIHseg>0 ) cout<<"> IH seg"<<fIHseg<<"  "<<endl;
  else           cout<<"> ! IH no hit  "<<endl;
}
