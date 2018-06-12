#ifndef MYHISTTRIGCHK_H
#define MYHISTTRIGCHK_H 1

#include "AnaInfo.h"
#include "EventHeader.h"
#include "ConfMan.h"
#include "MyHistTools.h"
#include "MyTools.h"
#include "MyAnaTools.h"

void initHistTrigChk();
void fillHistTrigChk(EventHeader *header, ConfMan *conf, BeamLineHitMan *blMan, CDSHitMan *cdsMan, AnaInfo *anaInfo);


#endif
