#ifndef IMPISIGMAUTIL
#define IMPISIGMAUTIL 1

#include <iostream>
#include "ConfMan.h"
#include "HitMan.h"
#include "CDSHitMan.h"
#include "LocalTrack.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "CDSTrackingMan.h"
#include "IMPiSigmaAnaPar.h"
#include <TLorentzVector.h>

namespace Util
{ 

  int GetCDHMul(CDSHitMan *cdsman,const int ntrack=0, const bool MCFlag=false);
  bool IsForwardCharge(BeamLineHitMan *blman);
  int GetCDHNeighboringNHits(const std::vector <int> &seg, const std::vector <int> &allhit);
  int GetNHitsCDCOuter(const TVector3 PosCDH, CDSHitMan *cdsman);
  double AnaBeamSpec(ConfMan *confman,BeamLineTrackMan *bltrackman,const int blc1id, const int blc2id);
  int CDSChargedAna(const bool docdcretiming,
                    LocalTrack *bpctrack,
                    CDSHitMan *cdsman,
                    CDSTrackingMan *trackman,
                    ConfMan *confman,
                    const TLorentzVector beam,
                    const double ctmt0, 
                    std::vector <int> &cdhseg,
                    std::vector <int> &pimid,
                    std::vector <int> &pipid,
                    std::vector <int> &kmid,
                    std::vector <int> &protonid,
                    const bool MCFlag=false
                    );
  double AnalyzeT0(BeamLineHitMan *blman, ConfMan *confman);
  int BeamPID(EventHeader *header, const double ctmt0, BeamLineHitMan *blman); 
  int EveSelectBeamline(BeamLineTrackMan *bltrackman, 
                        CDSTrackingMan *trackman,
                        ConfMan *confman,
                        int &blc1id, 
                        int &blc2id, 
                        int &bpcid);

  void AnaCDHHitPos(const double meas_tof, const double beta_calc, 
                 LocalTrack *bpc,
                 const TLorentzVector LVec_beambf,
                 CDSTrack *track,
                 CDSHitMan *cdsman,
                 ConfMan *confman);

};

#endif
