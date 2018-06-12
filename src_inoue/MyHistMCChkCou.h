#ifndef MYHISTMCCHKCOU_H
#define MYHISTMCCHKCOU_H 1

#include "KnuclRootData.h"
#include "MyHistTools.h"

void initHistMCChkCou();
void fillHistMCChkCou(EventHeaderMC *headerMC, DetectorData *detData, MCData *mcDatax);

#endif
