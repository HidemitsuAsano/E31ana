#ifndef TemplateFitter_H
#define TemplateFitter_H 1

#include "TObject.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TRandom.h"
#include "TMinuit.h"
#include "TString.h"

#include <deque>
#include <iostream>

class TemplateFitter : public TObject
{
 public:
  TemplateFitter();
  TemplateFitter(TH1 *data, TObjArray *mc);
  virtual ~TemplateFitter(){};

 private:
  TFractionFitter *fFitter;
  bool fIsFitted;
  TH1 *fData;
  TObjArray *fMC;

  double fChi2;
  int fNDF;

  std::vector<TString> fName;
  std::vector<double> fFraction;
  std::vector<double> fError;
  std::deque<bool> fFix;

 public:
  bool isFitted() const { return fIsFitted; }
  TFractionFitter *fitter(){ return fFitter; }
  int n() const { return fMC->GetEntries(); }
  double frac(int i) const { return fFraction[i]; }
  double err(int i) const { return fError[i]; }
  bool isFix(int i) const { return fFix[i]; }
  bool getResult(int i, double &p, double &e);
  bool getScaleErr(int i, double &scale, double &err);
  double chi2() const { return fChi2; }
  int ndf() const { return fNDF; }

  void fix(int i){ fFix[i]=true; }
  void fix(int i, double val){ fFix[i]=true; fFraction[i]=val; }
  void setScale(int i, double val);
  void setFixScale(int i, double val){ fFix[i]=true; setScale(i, val); }
  void setFrac(int i, double val){ fFraction[i]=val; }

  bool fit(int dump_level=0);

  double getEntries(TH1* h1);

  ClassDef(TemplateFitter, 1);
};
#endif
