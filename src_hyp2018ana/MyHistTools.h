#ifndef MYHISTTOOLS
#define MYHISTTOOLS_H 1

#include "ConfMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "AnaInfo.h"
#include "TNtuple.h"
#include "MyTools.h"

namespace MyHistTools
{
  void fillTH(TString name, double val2);
  void fillTH(TString name, double val1, double val2);

  void initLpim();
  void fillLpim(EventHeader *header, AnaInfo *anaInfo, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan);

  void initPC();
  void fillPC(BeamLineHitMan *blMan);
  void fillT0PC(const BeamInfo &beam, const ForwardChargeInfo &forwoard, BeamLineHitMan *blMan);

  void initCVC();
  void fillCVC(BeamLineHitMan *blMan);
  void fillT0CVC(const BeamInfo &beam, const ForwardChargeInfo &forwoard, BeamLineHitMan *blMan);

  void initT0();
  void fillT0(BeamLineHitMan *blMan);

  void initBHD();
  void fillBHD(BeamLineHitMan *blMan);

  void initBHDT0();
  void fillBHDT0(double tof, EventHeader *header);
  void fillBHDT0(const BeamInfo &beam, BeamLineHitMan *blMan);

  void initFDC1();
  void fillFDC1(LocalTrack* track);
  void fillFDC1(BeamLineHitMan *blMan);

  void initNC();
  void fillNC(AnaInfo *info, BeamLineHitMan *blMan);

  void initFC();
  void fillFC(ConfMan *conf, EventHeader *header, AnaInfo *anaInfo, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan);

  void initCDS();
  void fillCDS(AnaInfo *info, CDSTrackingMan *cdstrackMan);
};
#endif
