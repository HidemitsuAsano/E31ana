#include "CDS2Info.h"

using namespace std;

CDS2Info::CDS2Info() :
  fFlag(false), fTrackID1(-1), fTrackID2(-1), fPID1(CDS_DEFAULT), fPID2(CDS_DEFAULT), 
  fVertex1(DEFVECT), fVertex2(DEFVECT), fVertexBeam(DEFVECT), fMomentum1(DEFVECT), fMomentum2(DEFVECT)
{
}

TLorentzVector CDS2Info::lmom1() const
{
  TLorentzVector lmom;
  lmom.SetVectM(fMomentum1, pdgMass1());
  return lmom;
}

TLorentzVector CDS2Info::lmom2() const
{
  TLorentzVector lmom;
  lmom.SetVectM(fMomentum2, pdgMass2());
  return lmom;
}

CDS2Info CDS2Info::swap() const
{
  CDS2Info swap;
  swap.SetTrackID1(trackID2());
  swap.SetTrackID2(trackID1());
  swap.SetPID1(pid2());
  swap.SetPID2(pid1());
  swap.SetVertex1(vertex2());
  swap.SetVertex2(vertex1());
  swap.SetVertexBeam(vertexBeam());
  swap.SetMomentum1(momentum2());
  swap.SetMomentum2(momentum1());
  return swap;
}

void CDS2Info::dump() const
{
  cout<<"===== CDS2Info::dump ====="<<endl;
  cout<<"> Flag : "<<boolalpha<<fFlag<<"  TrackID1 : "<<fTrackID1<<" TrackID2 : "<<fTrackID2<<endl;
  if( fPID1==CDS_PiPlus   )      cout<<"> Pion+   ";
  else if( fPID1==CDS_Proton   ) cout<<"> Proton  ";
  else if( fPID1==CDS_Deuteron ) cout<<"> Deutern ";
  else if( fPID1==CDS_Helium3  ) cout<<"> He3     ";
  else if( fPID1==CDS_Triton   ) cout<<"> Triton  ";
  else if( fPID1==CDS_PiMinus  ) cout<<"> Pion-   ";
  else if( fPID1==CDS_Kaon     ) cout<<"> Kaon-   ";
  else if( fPID1==CDS_Other    ) cout<<"> Other   ";
  else if( fPID1==CDS_DEFAULT  ) cout<<"> DEFAULT ";
  else                           cout<<"> unknown ";

  if( fPID2==CDS_PiPlus   )      cout<<"> Pion+   "<<endl;
  else if( fPID2==CDS_Proton   ) cout<<"> Proton  "<<endl;
  else if( fPID2==CDS_Deuteron ) cout<<"> Deutern "<<endl;
  else if( fPID2==CDS_Helium3  ) cout<<"> He3     "<<endl;
  else if( fPID2==CDS_Triton   ) cout<<"> Triton  "<<endl;
  else if( fPID2==CDS_PiMinus  ) cout<<"> Pion-   "<<endl;
  else if( fPID2==CDS_Kaon     ) cout<<"> Kaon-   "<<endl;
  else if( fPID2==CDS_Other    ) cout<<"> Other   "<<endl;
  else if( fPID2==CDS_DEFAULT  ) cout<<"> DEFAULT "<<endl;
  else                           cout<<"> unknown "<<endl;

  cout<<"> Vtx1 "<<Form("(%lf, %lf, %lf) [cm]", fVertex1.X(), fVertex1.Y(), fVertex1.Z())<<endl;
  cout<<"> Vtx2 "<<Form("(%lf, %lf, %lf) [cm]", fVertex2.X(), fVertex2.Y(), fVertex2.Z())<<endl;
  cout<<"> VtxBeam "<<Form("(%lf, %lf, %lf) [cm]", fVertexBeam.X(), fVertexBeam.Y(), fVertexBeam.Z())<<endl;
  cout<<"> DCA : "<<dca()<<"[cm]  displaced : "<<displaced()<<"[cm]"<<endl;
  cout<<"> momentum1 "<<Form("(%lf, %lf, %lf) [GeV/c]", fMomentum1.X(), fMomentum1.Y(), fMomentum1.Z())<<endl;
  cout<<"> momentum2 "<<Form("(%lf, %lf, %lf) [GeV/c]", fMomentum2.X(), fMomentum2.Y(), fMomentum2.Z())<<endl;
}
