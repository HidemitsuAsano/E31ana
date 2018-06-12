#ifndef CDS2INFO_HH
#define CDS2INFO_HH 1

#include "CDSTrackingMan.h"
#include "TLorentzVector.h"

class CDS2Info : public TObject
{
 public:
  CDS2Info();
  virtual ~CDS2Info(){};

 private:
  bool fFlag;
  int fTrackID1;
  int fTrackID2;
  int fPID1;
  int fPID2;
  TVector3 fVertex1;
  TVector3 fVertex2;
  TVector3 fVertexBeam;
  TVector3 fMomentum1;
  TVector3 fMomentum2;

 public:
  void SetFlag(const bool &flag){ fFlag=flag; };
  void SetTrackID1(const int &id){ fTrackID1=id; };
  void SetTrackID2(const int &id){ fTrackID2=id; };
  void SetPID1(const int &pid){ fPID1=pid; };
  void SetPID2(const int &pid){ fPID2=pid; };
  void SetVertex1(const TVector3 &pos){ fVertex1=pos; };
  void SetVertex2(const TVector3 &pos){ fVertex2=pos; };
  void SetVertexBeam(const TVector3 &pos){ fVertexBeam=pos; };
  void SetMomentum1(const TVector3 &mom){ fMomentum1=mom; };
  void SetMomentum2(const TVector3 &mom){ fMomentum2=mom; };

  bool flag() const { return fFlag; };
  int trackID1() const { return fTrackID1; };
  int trackID2() const { return fTrackID2; };
  int pid1() const { return fPID1; };
  int pid2() const { return fPID2; };
  CDSTrack* track1(CDSTrackingMan *cdstrackMan){ return cdstrackMan->Track(fTrackID1); };
  CDSTrack* track2(CDSTrackingMan *cdstrackMan){ return cdstrackMan->Track(fTrackID2); };
  double pdgMass1() const { return (0<=fPID1 && fPID1<(int)(sizeof(cdsMass)/sizeof(cdsMass[0]))) ? cdsMass[fPID1] : 0.0; };
  double pdgMass2() const { return (0<=fPID2 && fPID2<(int)(sizeof(cdsMass)/sizeof(cdsMass[0]))) ? cdsMass[fPID2] : 0.0; };
  TVector3 vertex1() const { return fVertex1; };
  TVector3 vertex2() const { return fVertex2; };
  TVector3 vertexBeam() const { return fVertexBeam; };
  TVector3 momentum1() const { return fMomentum1; };
  TVector3 momentum2() const { return fMomentum2; };

  TVector3 vertexMean() const { return 0.5*(fVertex1+fVertex2); };
  double dca() const { return (fVertex1-fVertex2).Mag(); };
  double displaced() const { return (vertexMean()-fVertexBeam).Mag(); };

  double im() const { return (lmom1()+lmom2()).M(); };
  TLorentzVector lmom() const { return (lmom1()+lmom2()); }
  TLorentzVector lmom1() const;
  TLorentzVector lmom2() const;

  CDS2Info swap() const; // return swaped Object

  void dump() const;

  ClassDef(CDS2Info, 1);
};
#endif
