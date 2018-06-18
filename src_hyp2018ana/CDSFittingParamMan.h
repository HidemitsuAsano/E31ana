// CDSFittingParamMan.h

#ifndef CDSFittingParamMan_h
#define CDSFittingParamMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"
#include "TF1.h"

class CDSFittingParamMan : public TObject
{
 public:
  CDSFittingParamMan();
  CDSFittingParamMan( const std::string & filename );
  CDSFittingParamMan( const CDSFittingParamMan &right );
  ~CDSFittingParamMan();

  void SetFileName( const std::string & filename );
  bool Initialize();

 private:

  std::string FileName;
  int MAXCDCHIT;
  double MAGFIELD;
  int MAXCHI;
  int CDCTDC_U;
  int CDCTDC_L;

  std::map<int, TF1*> PIDFuncMin;
  std::map<int, TF1*> PIDFuncMax;

 public:
  std::string GetFileName() { return FileName; }
  int GetMaxCDCHit() { return MAXCDCHIT; }
  double GetMagneticField() { return MAGFIELD; }
  int GetMaxChi() { return MAXCHI; }

  int GetUpperLimitCDCTDC() { return CDCTDC_U; }
  int GetLowerLimitCDCTDC() { return CDCTDC_L; }

  void SetMagneticField(double f){ MAGFIELD = f;}
  int PID(double mass2, double mom);
  bool PID(int pid, double mass2, double mom);

  TF1* PIDfunc_min(int pid){ return PIDFuncMin.at(pid); };
  TF1* PIDfunc_max(int pid){ return PIDFuncMax.at(pid); };

  void Clear();
  
  ClassDef( CDSFittingParamMan, 1 );
};

#endif
