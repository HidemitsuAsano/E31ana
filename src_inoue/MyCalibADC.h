#ifndef MYCALIBADC_H
#define MYCALIBADC_H 1

#include "GlobalVariables.h"
#include "ConfMan.h"
#include "MyCalibTools.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLine.h>

//return [0]:up offset [1]:up gain [2]:up threshold
//       [3]:up offset [4]:up gain [5]:up threshold
std::vector<double> CalibADC(ConfMan *conf, TFile *of, int cid, int seg, double dE);

#endif
