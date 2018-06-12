#ifndef FORWARDCHARGINFO
#define FORWARDCHARGINFO_HH 1

#include "BeamLineTrackMan.h"
#include "BeamInfo.h"
#include "CDSInfo.h"
#include "ProtonArm.h"

class ForwardChargeInfo : public TObject
{
 public:
  ForwardChargeInfo();
  virtual ~ForwardChargeInfo(){};

 private:
  TVector3 fVertex;
  TVector3 fFitVtx;
  TVector3 fFitFDC1;
  double fDiffCounter;
  BLDCTrackInfo fFDC1info;

  int fSeg;
  double fTime;

  int fStep;
  bool fPCflag;
  int fPID;
  double fMass2ByAng;
  double fMass2ByRK;
  double fFlightLengthByArc;
  double fFlightLengthByRK;
  double fBeta;
  double fMomByTOF;
  double fMomByAng;
  double fMomByRK;
  double fCalcTOF;
  double fOffset;
  TVector3 fHitPos;
  TVector3 fMomentum;

 public:
  void SetStep(const int &step){ fStep=step; }
  void SetPID(const int &id){ fPID=id; }
  void SetVertex(const TVector3 &pos){ fVertex=pos; }
  void SetFDC1(const BLDCTrackInfo &track){ fFDC1info=track; };
  void SetHodo(HodoscopeLikeHit *hit);

  TVector3 vertex() const { return fVertex; }
  TVector3 fitVertex() const { return fFitVtx; }
  BLDCTrackInfo FDC1() const { return fFDC1info; };

  double time() const { return fTime; }
  int step() const { return fStep; }
  bool isPC() const { return fPCflag; }
  int seg() const { return fSeg; }
  int pid() const { return fPID; }
  double diffCounter() const { return fDiffCounter; }
  TVector3 fitFDC1() const { return fFitFDC1; }

  double mass2byAng() const { return fMass2ByAng; }
  double mass2byRK() const { return fMass2ByRK; }
  double fl() const { return fFlightLengthByRK>0 ? fFlightLengthByRK : fFlightLengthByArc; }
  double flByArc() const { return fFlightLengthByArc; }
  double flByRK() const { return fFlightLengthByRK; }
  double beta() const { return fBeta; }
  double angle(ConfMan *conf) const;
  double radius(ConfMan *conf) const { return (100./sqrt(2*(1-cos(angle(conf))))); }

  double momByTOF() const { return fMomByTOF; }
  double momByAng() const { return fMomByAng; }
  double momByRK() const { return fMomByRK; }
  double calc_tof() const { return fCalcTOF; }
  double offset() const { return fOffset; }
  TVector3 momentum() const { return fMomentum; }
  TVector3 hitpos() const { return fHitPos; }

  HodoscopeLikeHit *hodo(BeamLineHitMan *blMan) const;

  TLorentzVector lmom() const { TLorentzVector lmom; lmom.SetVectM(fMomentum, particleMass[fPID]); return lmom; };

  bool calc_forward(BeamLineHitMan *blMan, BeamInfo *beam, CDSInfo *cds, ConfMan *conf);
  bool fit_forward(BeamLineHitMan *blMan, BeamInfo *beam, CDSInfo *cds, ConfMan *conf, bool refit=false);
  bool calcMomByTOF(BeamInfo *beam, CDSInfo *cds);

  bool calc_simple_beam_through(const BeamInfo &beam, ConfMan *conf, BeamLineHitMan *blMan);
  bool calc_simple_p_through(const BeamInfo &beam, ConfMan *conf, BeamLineHitMan *blMan);
  std::vector<std::vector<double> > calc_beam_through_wUSWK(const BeamInfo &beam, ConfMan *conf, BeamLineHitMan *blMan);
  std::vector<std::vector<double> > fit_beam_through_wPoints(const BeamInfo &beam, ConfMan *conf);
  void dump() const;
  double findMass2(const double &fl, const double &mom, const double &tof);

  ClassDef(ForwardChargeInfo, 1);
};
#endif
