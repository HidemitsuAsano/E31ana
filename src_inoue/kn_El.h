#ifndef KN_EL_H
#define KN_EL_H 1

#include "SimTools.h"
#include "TFile.h"
#include "GlobalVariables.h"
#include "KnuclRootData.h"
#include "TH1.h"
#include "TH2.h"

void initHist();
int fillHist(DetectorData *detData, MCData* mcData, ReactionData *reacData);
void writeHist();

#endif
