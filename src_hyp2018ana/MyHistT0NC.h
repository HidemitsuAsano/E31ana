#ifndef MYHISTT0NC_HH
#define MYHISTT0NC_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistT0NC();
void fillT0NC(BeamLineHitMan *blMan, AnaInfo *anaInfo);
void fillSlewingNC(BeamLineHitMan *blMan, AnaInfo *anaInfo, ForwardNeutralInfo *fnInfo);


#endif
