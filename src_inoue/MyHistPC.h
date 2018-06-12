#ifndef MYHISTPC_HH
#define MYHISTPC_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistPC();
void fillPC_ADC(BeamLineHitMan *blMan);

#endif
