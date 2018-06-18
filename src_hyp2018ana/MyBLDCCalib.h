#ifndef MYBLDCCALIB_H
#define MYBLDCCALIB_H 1

#include "GlobalVariables.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TLine.h>

#include "ConfMan.h"
#include "MyCalibTools.h"

void Chamber_offset0(ConfMan *conf, TFile *of, int cid, int lay, int wire=0);
void Chamber_dxdt(ConfMan *conf, std::ostream &ofs=std::cout);
void Chamber_dxdt(ConfMan *conf, int cid, std::ostream &ofs=std::cout);

#endif
