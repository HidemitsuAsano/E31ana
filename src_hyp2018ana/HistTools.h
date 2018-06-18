#ifndef HISTTOOLS_HH
#define HISTTOOLS_HH 1

#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"

#include "LocalTrack.h"
#include "ConfMan.h"
#include "CDSTrack.h"

void initHistBLDC();
void initHistCDH();
void initHistT0NC();
void initHistCDC();
void initHistIH();

void fillCDC(CDSTrack *track, const double &beta, CDSHitMan *cdsMan, ConfMan *conf);
void fillBLC1a(LocalTrack *track);
void fillBLC1b(LocalTrack *track);
void fillBLC2a(LocalTrack *track);
void fillBLC2b(LocalTrack *track);
void fillBPC(LocalTrack *track);
void fillFDC1(LocalTrack *track);
void fillHistIH(CDSHitMan *cdsMan);

#endif
