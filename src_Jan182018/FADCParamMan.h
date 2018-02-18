// FADCParamMan.h

#ifndef FADCParamMan_h
#define FADCParamMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"

class FADCParamMan : public TObject
{
 public:
  FADCParamMan();
  FADCParamMan( const std::string & filename );
  FADCParamMan( const FADCParamMan &right );
  ~FADCParamMan();

  void SetFileName( const std::string & filename );
  bool Initialize();

 private:

  std::string FileName;
  int BASE_L;
  int BASE_U;
  int PEAK_L;
  int PEAK_U;
  int POST_L;
  int POST_U;
  int NPoint;
  int Interval; //in ns
  std::string BaseFunc;
  std::string PeakFunc;
  std::string PostFunc;

 public:
  std::string GetFileName() { return FileName; }
  std::string GetBaseFunc() { return BaseFunc; }
  std::string GetPeakFunc() { return PeakFunc; }
  std::string GetPostFunc() { return PostFunc; }
  int GetUpperLimitBase() { return BASE_U; }
  int GetLowerLimitBase() { return BASE_L; }
  int GetUpperLimitPeak() { return PEAK_U; }
  int GetLowerLimitPeak() { return PEAK_L; }
  int GetUpperLimitPost() { return POST_U; }
  int GetLowerLimitPost() { return POST_L; }
  int GetNPoint() { return NPoint; }
  int GetInterval() { return Interval; }

  void Clear();
  
  ClassDef( FADCParamMan, 1 );
};

#endif
