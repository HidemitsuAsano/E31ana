#ifndef HYHISTPCCVC_HH
#define HYHISTPCCVC_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistT0PCCVC();
void fillT0PCCVC_BT(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamInfo *beam);

#endif
