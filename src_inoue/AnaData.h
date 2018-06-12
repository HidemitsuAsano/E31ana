#ifndef AnaData_h
#define AnaData_h

#include "CDS1Data.h"
#include "CDS2Data.h"
#include "EventHeader.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

class AnaData : public TObject
{
 public:
  AnaData();
  virtual ~AnaData(){};

 private:
  int fEventNumber;
  bool fVtxFlag;
  double fD5mom;
  double fT0time;
  int fBeamPID;
  TLorentzVector fBeamLmom;
  TVector3 fVtx;

  std::vector<CDS1Data> fCDS1Data;
  std::vector<CDS2Data> fCDS2Data;

  int fForwardPID;
  int fFCID;
  int fFseg;
  double fFdE;
  double fFBeta;
  TLorentzVector fForwardLmom;

 public:
  int ev() const { return fEventNumber; };
  int beam_pid() const { return fBeamPID; };
  TLorentzVector beam_lmom() const { return fBeamLmom; };
  bool isVtx() const { return fVtxFlag; };
  TVector3 vtx() const { return fVtx; };
  int forward_pid() const { return fForwardPID; };
  TLorentzVector forward_lmom() const { return fForwardLmom; };

  int nCDS1() const { return fCDS1Data.size(); };
  int nCDS2() const { return fCDS2Data.size(); };
  CDS1Data *CDS1(const int &i){ return &fCDS1Data[i]; };
  CDS2Data *CDS2(const int &i){ return &fCDS2Data[i]; };
  CDS1Data *CDS1(const int &pid, const int &n);
  CDS2Data *CDS2(const int &pid1, const int &pid2, const int &n);

  int fpid() const { return fForwardPID; };
  int fcid() const { return fFCID; };
  int fseg() const { return fFseg; };
  double fdE() const { return fFdE; };
  double fbeta() const { return fFBeta; };
  TLorentzVector flmom() const { return fForwardLmom; };

  bool set(EventHeader *header, ConfMan *conf, const double &D5, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan);

  void clear();

  ClassDef(AnaData, 1);
};
#endif
