#ifndef MyFractionFitter_h
#define MyFractionFitter_h 1

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TH1.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"

class MyFractionFitter : public TObject
{
 public:
  MyFractionFitter(TH1* data, TObjArray *mc);
  ~MyFractionFitter();

 private:
  TVirtualFitter *fFitter;
  TH1* fData;
  TH1* fPlot;
  std::vector<double> fFraction;
  std::vector<double> fError;

  std::vector<TH1*> fMC;
  std::vector<TH1*> fPred;

 public:
  int Fit();
  double fcn(int npar, double *param);

 private:
  double func(int bin, double x);
  void FindPrediction(int bin, const double min);
  void check();

  ClassDef(MyFractionFitter, 0);
};
#endif
