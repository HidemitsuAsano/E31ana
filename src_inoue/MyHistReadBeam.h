#ifndef MYHISTREADBEAM_HH
#define MYHISTREADBEAM_HH

#include "AnaInfo.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "BeamLineHitMan.h"
#include "MyAnaTools.h"

void initHistReadBeam();
void fillReadBeam_T0(EventHeader *header, BeamLineHitMan *blMan);
void fillReadBeam_BHDT0(EventHeader *header, AnaInfo *anaInfo);
void fillReadBeam_BLC1(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo);
void fillReadBeam_BLC2(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo);
void fillReadBeam_BPC(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo, BeamLineHitMan *blMan);
void fillReadBeam_D5(EventHeader *header, AnaInfo *anaInfo);
void fillReadBeam_D5BHD(EventHeader *header, AnaInfo *anaInfo);
void fillReadBeam_BLC2BPC(EventHeader *header, AnaInfo *anaInfo);
void fillReadBeam_Profile(ConfMan *conf, EventHeader *header, AnaInfo *anaInfo);

void fillReadBeam_FDC1(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo, BeamLineHitMan *blMan);

#endif
