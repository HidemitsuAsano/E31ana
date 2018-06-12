#ifndef CDS1Data_h
#define CDS1Data_h

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "LocalTrack.h"

#include "float.h"

#include "GlobalVariables.h"
#include "TrackTools.h"

class CDS1Data : public TObject
{
 public:
  CDS1Data();
  CDS1Data(const int &id);
  virtual ~CDS1Data(){};

 private:
  int fID;
  int fPID;
  int fBeamPID;
  int fCDHseg;
  double fTOF;
  double fMass2;
  double fMom;
  double fBeta;
  double fdE;
  TVector3 fBeamMom;
  TVector3 fCDSMom;

  TVector3 fVtxBeam;
  TVector3 fVtxCDS;

 public:
  int id() const { return fID; };
  int pid() const { return fPID; };
  int CDHseg() const { return fCDHseg; };
  double tof() const { return fTOF; };
  double mass2() const { return fMass2; };
  double mom() const { return fMom; };
  double beta() const { return fBeta; };
  double de() const { return fdE; };
  TLorentzVector lmom() const { 
    TLorentzVector lmom;
    lmom.SetVectM(fCDSMom, cdsMass[fPID]);
    return lmom; 
  };
  TLorentzVector lmomBeam() const { 
    TLorentzVector lmom;
    lmom.SetVectM(fBeamMom, parMass[fBeamPID]);
    return lmom; 
  };

  TVector3 vtxBeam() const { return fVtxBeam; };
  TVector3 vtxCDS() const { return fVtxCDS; };
  double dis() const { return (fVtxBeam-fVtxCDS).Mag(); };

  double mm(const TLorentzVector &tgt) const {
    TLorentzVector beam_lmom, cds_lmom;
    beam_lmom.SetVectM(fBeamMom, parMass[fBeamPID]);
    cds_lmom.SetVectM(fCDSMom, cdsMass[fPID]);
    return (beam_lmom+tgt-cds_lmom).M();
  }

  // This method return status 0: succeed  1: no CDH hit  2: false FindMass2C  3: false CalcVertexTimeLength  4: false GetMomentum
  int set(const int &beam_pid, const double &t0time, const double &mom, LocalTrack *beam, CDSTrack *cds, CDSHitMan *cdsMan);

  ClassDef(CDS1Data, 1);
};
#endif
