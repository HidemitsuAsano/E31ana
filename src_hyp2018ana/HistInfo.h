#ifndef HISTINFO_H
#define HISTINFO_h

#include "TObject.h"
#include "TString.h"

#include <vector>
#include <iostream>

class HistInfo : public TObject
{
 public:
  HistInfo();
  HistInfo(std::vector<double> bin, int n);
  HistInfo(std::vector<double> bin, std::vector<TString> names);
  virtual ~HistInfo(){};

 private:
  std::vector<double> fBin;
  std::vector<TString> fName;
  std::vector<double> fUnder;
  std::vector<double> fOver;
  std::vector<std::vector<double> > fVal;

 public:
  int nBin() const { return fBin.size()-1; };
  int nVal() const { return fName.size(); };
  int index(double val);
  TString name(int i) const { return fName.at(i); }
  double under(int i) const { return fUnder.at(i); }
  double bin_min(int i) const { return fBin.at(i); }
  double bin_max(int i) const { return fBin.at(i+1); }
  double over(int i) const { return fOver.at(i); }
  double val(int i, int j=0) const { return fVal.at(i).at(j); }
  double val(int i, TString name) const;

  void setName(int i, const TString &name){ fName.at(i)=name; }
  void setVal(int i, int j, double val){ fVal.at(i).at(j)=val; }
  bool setVal(int i, TString name, double val);

  void dump();

  ClassDef(HistInfo, 1);
};
#endif
