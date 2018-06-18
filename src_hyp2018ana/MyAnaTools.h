#ifndef MYANATOOLS_HH
#define MYANATOOLS_HH 1

#include "MyTools.h"
#include "MyBeamAnaParam.h"
#include "ConfMan.h"
#include "TGeoPhysicalNode.h"

namespace MyAnaTools
{
  bool goodBeam(AnaInfo *anaInfo);
  bool trigBLC1(AnaInfo *anaInfo, ConfMan *conf);
  bool trigBLC2(AnaInfo *anaInfo, ConfMan *conf);
  bool trigBPC(AnaInfo *anaInfo, ConfMan *conf, BeamLineHitMan *blMan);
  bool trigFDC1(AnaInfo *anaInfo, ConfMan *conf, BeamLineHitMan *blMan);

  bool isTOFKaon(double tof);
  bool anaBLC1(AnaInfo *anaInfo);
  bool anaBLC2(AnaInfo *anaInfo);
  bool anaBPC(AnaInfo *anaInfo);
  bool anaD5(AnaInfo *anaInfo);
  bool connectBLC2BPC(AnaInfo *anaInfo);
  bool connectD5BHD(AnaInfo *anaInfo);
  bool beamFiducial(AnaInfo *anaInfo, ConfMan *conf);

  bool anaFDC1(AnaInfo *anaInfo);

  TLorentzVector target_lmom();
  int targetNA();

  double T0BVC_tof(AnaInfo *anaInfo, BeamLineHitMan *blMan);
};
#endif
