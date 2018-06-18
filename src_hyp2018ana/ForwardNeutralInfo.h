#ifndef FORWARDNEUTRALINFO_HH
#define FORWARDNEUTRALINFO_HH 1

#include "ConfMan.h"
#include "TLorentzVector.h"
#include "BeamLineHitMan.h"
#include "BeamInfo.h"

class ForwardNeutralInfo : public TObject
{
 public:
  ForwardNeutralInfo();
  virtual ~ForwardNeutralInfo(){};

 private:
  int fPID;
  int fSeg;  // hit diciding timing
  double fOffset;
  double fTime;
  double fdE;
  TVector3 fVertex;
  TVector3 fHitPos;
  double fBeta;
  TVector3 fMomentum;
  std::vector<int> fClusterSeg;

 public:
  void SetPID(const int &pid){ fPID=pid; }
  void SetHodo(const HodoscopeLikeHit *hit);
  void SetVertex(const TVector3 &vtx){ fVertex=vtx; }
  void SetHitPos(const TVector3 &pos){ fHitPos=pos; }
  void SetCluster(const int &seg){ fClusterSeg.push_back(seg); }

  int pid() const { return fPID; }
  int seg() const { return fSeg; }
  int nCluster() const { return fClusterSeg.size(); }
  double offset() const { return fOffset; }
  double beta() const { return fBeta; }
  double time() const { return fTime; }
  double dE() const { return fdE; }
  double fl() const { return (fVertex-fHitPos).Mag(); }
  TVector3 vertex() const { return fVertex; }
  TVector3 hitpos() const { return fHitPos; }
  TVector3 momentum() const { return fMomentum; }
  TLorentzVector lmom() const;
  HodoscopeLikeHit* NC(BeamLineHitMan *blMan) const;
  HodoscopeLikeHit* NC(const int &i, BeamLineHitMan *blMan) const;

  void calc(const BeamInfo *beam, bool sim=false);

  void dump() const;

  ClassDef(ForwardNeutralInfo, 1);
};
#endif
