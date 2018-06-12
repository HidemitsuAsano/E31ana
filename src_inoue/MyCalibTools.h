#ifndef MYCALIBTOOLS_H
#define MYCALIBTOOLS_H 1

#include <iostream>
#include <TH1.h>
#include <TGraph.h>
#include <TF1.h>
#include <TROOT.h>

#include "ConfMan.h"

namespace MyCalibTools
{
  bool setTimeOffset(ConfMan *conf, int cid, int lay, int wire, double offset);
  bool setTimeOffset(ConfMan *conf, int cid, int seg, double offset);

  TF1* fitGaus_3sigma(TH1F *h1);
  TF1* fitGaus_3sigma(TH1F *h1, double mean, double sigma);

  TF1* fitGaus_25sigma(TH1F *h1);
  TF1* fitGaus_25sigma(TH1F *h1, double mean, double sigma);

  TF1* fitGaus_2sigma(TH1F *h1);
  TF1* fitGaus_2sigma(TH1F *h1, double mean, double sigma);

  TF1* fitLandau(TH1F *h1, double min);
};
#endif
