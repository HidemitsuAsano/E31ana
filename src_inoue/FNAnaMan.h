#ifndef FNANAMAN_H
#define FNANAMAN_H 1

#include "HistInfo.h"
#include "AnaInfo.h"
#include "MyParam.h"
//#include "MyAnaParam.h"
#include "TemplateFitter.h"

#include "KNpipi_cut.h"

#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TVirtualFitter.h"
#include "TMacro.h"
#include "TObjString.h"

#include "MyTools.h"

#include <climits>
#include <iostream>
#include <fstream>

class FNAnaMan : public TObject
{
 public:
  FNAnaMan();
  virtual ~FNAnaMan(){};

 private:
  int fRebinF;

  HistInfo *fInfoK0;
  HistInfo *fInfoSm;
  HistInfo *fInfoSp;
  std::vector<TFile*> fFileK0;
  std::vector<TFile*> fFileSm;
  std::vector<TFile*> fFileSp;

  HistInfo *fInfoSim;
  std::vector<TFile*> fFileSim;

  TFile *fFileData;
  TFile *fFileOut;
  TFile *fFileDummy;
  TMacro *fLog;
  void initHist(std::string str);
  void clearHist(std::string str);
  void readAcc();
  void readData();
  void readSim(int index);
  bool fitKNpi(int index, int dump_level=0);
  
 public:
  void init(TString filename);
  void finit();
  void postAna();
  void postAna(TFile *f, const std::string &dirname, int index=-1);
  void fillHist(TDirectory *dir, const std::string &hist_name, double x, double y, int bin_index, int index);
  void fillHist(TDirectory *dir, const std::string &hist_name, double data, int bin_index, int index);

  void fillHist(TDirectory *dir, const std::string &hist_name, double x, double y);
  void fillHist(TDirectory *dir, const std::string &hist_name, double data);
  
  void write(TString dirname);
  void fillHistMC();

  bool fitIM(int dump_level=0);
  bool fitNpipiIM_K0(int dump_level=0);

  void fillHist(TFile *f, std::string suffix, int index=-1);
  bool fitKNpi_all(int dump_level=0);
  void fitAll(int dump_level=0);

  TGraphErrors *gra_chi2();
};
#endif
