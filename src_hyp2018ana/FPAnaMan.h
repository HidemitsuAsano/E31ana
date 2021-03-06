#ifndef FPANAMAN_H
#define FPANAMAN_H 1

#include <fstream>
#include <string>
#include <vector>

#include "AnaInfo.h"
#include "TFile.h"
#include "TString.h"
#include "TNtuple.h"
#include "TGraphErrors.h"

#include "MyAnaParam.h"
#include "HistInfo.h"

#include "CS_Param.h"

class FPAnaMan
{
 public:
  FPAnaMan();
  ~FPAnaMan(){};

 public:
  void init(TString conffile);
  void readSim(int index);
  void finit();

 private:
  HistInfo *fParamInfo;
  HistInfo *fHistInfo;
  TFile *f_data;
  TFile *f_sim1;
  TFile *f_sim2;
  std::vector<TFile*> f_MC;
  TFile *f_out;

  std::vector<double> fBinMin;
  std::vector<double> fBinMax;
  std::vector<double> fNumData[2];

  std::vector<double> fTrig[2];
  std::vector<double> fEff[2];
  std::vector<double> fSolidAng[2];

 public:
  void readData_tmp();
  void initHist(int index);
  void makeAcc();
  void makeLineShape();
  void makeLineShape_2step();
  void makeCS();
  void fillHist();

  int getMMIndex(double mm);

};
#endif
