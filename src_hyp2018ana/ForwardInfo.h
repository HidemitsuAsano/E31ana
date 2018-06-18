#ifndef FORWARDINFO_HH
#define FORWARDINFO_HH 1

#include "BeamLineTrackMan.h"
#include "TLorentzVector.h"

class ForwardInfo : public TObject
{
 public:
  ForwardInfo();
  virtual ~ForwardInfo(){};

 private:
  int fPID;
  int fCounterInfo; // CID : SEG 3digits
  double fBeta;
  double fTime;
  double fFlightLength;
  TVector3 fVertex;
  TVector3 fMomentum;

 public:
  void SetPID(const int &id){ fPID=id; };
  void SetCounter(const int &cid, const int &seg){ fCounterInfo=1000*cid+seg; };
  void SetBeta(const double &beta){ fBeta=beta; };
  void SetTime(const double &time){ fTime=time; };
  void SetFL(const double &fl){ fFlightLength=fl; };
  void SetVertex(const TVector3 &pos){ fVertex=pos; };
  void SetMomentum(const TVector3 &mom){ fMomentum=mom; };

  int pid() const { return fPID; };
  int cid() const { return fCounterInfo/1000; };
  int seg() const { return fCounterInfo%1000; };
  double beta() const { return fBeta; };
  double time() const { return fTime; };
  double fl() const { return fFlightLength; };
  double pdgMass() const { return (0<=fPID && fPID<sizeof(particleMass)/sizeof(particleMass[0])) ? particleMass[0] : 0.0; };
  TVector3 vertex() const { return fVertex; };
  TVector3 momentum() const { return fMomentum; };
  TLorentzVector lmom() const;

  ClassDef(ForwardInfo, 1);
};
#endif
