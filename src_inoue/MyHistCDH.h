#ifndef MYHISTCDH_HH
#define MYHISTCDH_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistCDH();
void fillCDH_ADC(CDSHitMan *cdsMan);
void fillCDH(BeamInfo *beam, CDSInfo *cdsInfo, CDSHitMan *cdsMan);

#endif
