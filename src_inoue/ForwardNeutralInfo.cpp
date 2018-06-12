#include "ForwardNeutralInfo.h"

using namespace std;

ForwardNeutralInfo::ForwardNeutralInfo()
  : fPID(F_Other), fSeg(-1), fOffset(DEFAULTD), fTime(DEFAULTD), fdE(DEFAULTD), 
    fVertex(DEFVECT), fHitPos(DEFVECT), fBeta(DEFAULTD), fMomentum(DEFVECT)
{
}

void ForwardNeutralInfo::calc(const BeamInfo *beam, bool sim)
{
  double beam_out, tmp_tof;
  //                                                                         kpMass
  ELossTools::CalcElossBeamTGeo(beam->T0pos(), beam->vertex(), beam->D5mom(), beam->mass(), beam_out, tmp_tof);
  double tof=fTime-beam->T0time()-tmp_tof;
  double fl=(fHitPos-fVertex).Mag();
  if( sim ) fl-=2.5;
  fBeta=fl/(tof*100*Const);
  if( fPID==F_Other){
    if( fBeta<0.95 ){
      fPID=F_Neutron;
      double mom=nMass*fBeta/sqrt(1-fBeta*fBeta);
      fMomentum=fHitPos-fVertex;
      fMomentum.SetMag(mom);
    }
    else{
      // cout<<Form("T0pos (%lf, %lf, %lf)", beam->T0pos().X(), beam->T0pos().Y(), beam->T0pos().Z())<<endl;
      // cout<<Form("Vtx (%lf, %lf, %lf)", beam->vertex().X(), beam->vertex().Y(), beam->vertex().Z())<<endl;
      fPID=F_Gamma;
      fOffset=fl/(100.*Const)+tmp_tof;
    }
  }
  else if( fPID==F_Neutron ){
    double mom=nMass*fBeta/sqrt(1-fBeta*fBeta);
    fMomentum=fHitPos-fVertex;
    fMomentum.SetMag(mom);
  }
}

void ForwardNeutralInfo::SetHodo(const HodoscopeLikeHit *hit)
{
  fSeg=hit->seg();
  fTime=hit->ctmean();
  fdE=hit->emean();
  fHitPos=hit->pos();
  fHitPos.SetY(hit->hitpos());
}

HodoscopeLikeHit* ForwardNeutralInfo::NC(BeamLineHitMan *blMan) const
{
  for( int i=0; i<blMan->nNC(); i++ ){
    if( fSeg==blMan->NC(i)->seg() ) return blMan->NC(i);
  }
  return 0;
}

HodoscopeLikeHit* ForwardNeutralInfo::NC(const int &id, BeamLineHitMan *blMan) const
{
  if( id<0 || (int)fClusterSeg.size()<id ) return 0;
  for( int i=0; i<blMan->nNC(); i++ ){
    if( fClusterSeg[id]==blMan->NC(i)->seg() ) return blMan->NC(i);
  }
  return 0;
}

TLorentzVector ForwardNeutralInfo::lmom() const
{
  TLorentzVector lmom;
  lmom.SetVectM(fMomentum, particleMass[fPID]);
  return lmom;
}

void ForwardNeutralInfo::dump() const
{
  cout<<"===== ForwardNeutralInfo::dump() ====="<<endl;
  if( fPID==F_Kaon ) cout<<"> Kaon"<<endl;
  else if( fPID==F_Pion     ) cout<<"> Pion"<<endl;
  else if( fPID==F_Proton   ) cout<<"> Proton"<<endl;
  else if( fPID==F_Deuteron ) cout<<"> Deuteron"<<endl;
  else if( fPID==F_Neutron  ) cout<<"> Neutron"<<endl;
  else if( fPID==F_Gamma    ) cout<<"> Gamma"<<endl;
  else if( fPID==F_Other    ) cout<<"> Other"<<endl;
  else cout<<"> unknown"<<endl;
  cout<<"> Time : "<<fTime<<"[ns]  seg"<<fSeg<<" dE : "<<fdE<<"[MeVee]  beta : "<<fBeta<<endl;
  cout<<"> Vertex "<<Form("(%lf, %lf, %lf) [cm]", fVertex.X(), fVertex.Y(), fVertex.Z())<<endl;
  cout<<"> HitPos "<<Form("(%lf, %lf, %lf) [cm]", fHitPos.X(), fHitPos.Y(), fHitPos.Z())<<endl;
  cout<<"> Momentum "<<Form("(%lf, %lf, %lf) [GeV/c]", fMomentum.X(), fMomentum.Y(), fMomentum.Z())<<endl;

  string out[7];
  for( int lay=0; lay<7; lay++ ){
    out[lay]="| | | | | | | | | | | | | | | | |";
  }
  for( int i=0; i<(int)fClusterSeg.size(); i++ ){
    int lay=(fClusterSeg[i]-1)/16;
    int seg2=(fClusterSeg[i]-1)%16;
    if( fClusterSeg[i]==fSeg ){
      out[lay][1+2*seg2]='#';
    }
    else{
      out[lay][1+2*seg2]='*';
    }
  }
  for(  int lay=0; lay<7; lay++ ){
    cout<<out[lay]<<endl;
  }
}
