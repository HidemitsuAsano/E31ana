#ifndef FITCONF_H
#define FITCONF_H 1

#include <sstream>

#include "AnaInfo.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "MyParam.h"
#include "MyHistTools.h"

class FitConf
{
 public:

 private:
  TFile *f_data;
  TFile *f_sim1;
  TFile *f_sim2;
  std::vector<TFile*> f_K0;
  std::vector<TFile*> f_Sm;
  std::vector<TFile*> f_Sp;
  std::vector<TFile*> f_BG;
  std::vector<double> fScaleK0;
  std::vector<double> fScaleSm;
  std::vector<double> fScaleSp;
  std::vector<double> fScaleBG;

  TFile *f_out;
  int fRebinF;

  double fNumDataK0;
  double fNumDataSm;
  double fNumDataSp;

  std::vector<double> fBinMin;
  std::vector<double> fBinMax;
  std::vector<double> fNumData;
  std::vector<double> fNumSim1;
  std::vector<double> fNumSim2;

  std::vector<double> fTrig1;
  std::vector<double> fTrig2;
  std::vector<double> fEff1;
  std::vector<double> fEff2;

  std::vector<double> fFrac1;
  std::vector<double> fFrac2;
  std::vector<double> fErr1;
  std::vector<double> fErr2;
  std::vector<double> fErr12;

  std::vector<double> fChi2;
  std::vector<double> fNDF;

 public:
  void init(char *filename);
  void initHist(char *suffix);
  void readData();
  void readSim1();
  void readSim2();

  void fillHist(TFile *f, const char *suffix);
  void fillHist_Sim(const int id);

  void finit();

 public:
  int getMMIndex(double mm);
  double getScale(const int i, double mm);
};
#endif
