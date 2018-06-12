#ifndef MCINFO_HH
#define MCINFO_HH 1

#include "TObject.h"
#include "TLorentzVector.h"

class MCInfo : public TObject
{
 public:
  MCInfo();
  virtual ~MCInfo(){};

 private:
  std::vector<TLorentzVector> fMeasured;
  std::vector<TLorentzVector> fMissing;

  std::vector<TLorentzVector> fMCMeasured;
  std::vector<TLorentzVector> fMCMissing;

 public:
  int n_measured(){ (int)fMeasured.size(); };
  int n_missing(){ (int)fMissing.size(); };

  TLorentzVector measured_lmom(const int &i){ return fMeasured[i]; };
  TLorentzVector missing_lmom(const int &i=0){ return fMissing[i]; };

  TLorentzVector mc_measured_lmom(const int &i){ return fMCMeasured[i]; };
  TLorentzVector mc_missing_lmom(const int &i=0){ return fMCMissing[i]; };

  void set_measured(const TLorentzVector &data, const TLorentzVector &mc){ fMeasured.push_back(data); fMCMeasured.push_back(mc); };
  void set_missing(const TLorentzVector &data, const TLorentzVector &mc){ fMissing.push_back(data); fMCMissing.push_back(mc); };

  void Clear();

  ClassDef(MCInfo, 1);
};
#endif
