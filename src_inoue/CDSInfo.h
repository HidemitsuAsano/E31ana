#ifndef CDSINFO_HH
#define CDSINFO_HH 1

#include "CDSTrackingMan.h"
#include "BeamLineTrackMan.h"
#include "TLorentzVector.h"

class CDSInfo : public TObject 
{
 public:
  CDSInfo();
  virtual ~CDSInfo(){};

 private:
  bool fFlag;
  int fTrackID;
  int fPID;
  int fCDHseg;
  int fIHseg;
  double fMom;
  double fMass2;
  double fBeta;
  double fOffset;
  TVector3 fVertexCDS;
  TVector3 fVertexBeam;
  TVector3 fMomentum;

 public:
  void SetFlag(const bool &flag){ fFlag=flag; };
  void SetTrackID(const int &id){ fTrackID=id; };
  void SetPID(const int &pid){ fPID=pid; };
  void SetCDHseg(const int &seg){ fCDHseg=seg; };
  void SetIHseg(const int &seg){ fIHseg=seg; };
  void SetMom(const double &mom){ fMom=mom; };
  void SetMass2(const double &mass2){ fMass2=mass2; };
  void SetBeta(const double &beta){ fBeta=beta; };
  void SetOffset(const double &offset){ fOffset=offset; };
  void SetVertexCDS(const TVector3 &pos){ fVertexCDS=pos; };
  void SetVertexBeam(const TVector3 &pos){ fVertexBeam=pos; };
  void SetMomentum(const TVector3 &mom){ fMomentum=mom; };

  bool flag() const { return fFlag; };
  CDSTrack *track(CDSTrackingMan *cdstrackMan) const { return cdstrackMan->Track(fTrackID); };
  CDSTrack *track(ConfMan *conf, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan) const;
  int CDHseg() const { return fCDHseg; };
  int IHseg() const { return fIHseg; };
  HodoscopeLikeHit *CDH(CDSHitMan *cdsMan) const;
  HodoscopeLikeHit *IH(CDSHitMan *cdsMan) const;
  int trackID() const { return fTrackID; };
  int pid() const { return fPID; };
  double pdgMass() const { return (0<=fPID && fPID<(int)(sizeof(cdsMass)/sizeof(cdsMass[0]))) ? cdsMass[fPID] : 0.0; };
  double mom() const { return fMom; };
  double mass2() const { return fMass2; };
  double beta() const { return fBeta; };
  double offset() const { return fOffset; };
  TVector3 vertexCDS() const { return fVertexCDS; };
  TVector3 vertexBeam() const { return fVertexBeam; };
  TVector3 momentum() const { return fMomentum; };

  double dca() const { return (fVertexCDS-fVertexBeam).Mag(); };
  TLorentzVector lmom() const;

  void dump() const;

  ClassDef(CDSInfo, 1);
};
#endif
