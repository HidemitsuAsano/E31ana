#ifndef MYTDCCALIB_H
#define MYTDCCALIB_H 1

#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TCanvas.h>

#include "ConfMan.h"

int TDCCalib(ConfMan *conf, int cr, int sl, int ch, double period, TFile *of);
// 0: successfully finish   1: some error happen parameter is updated   2~: error happen parameter is not updated

#endif
