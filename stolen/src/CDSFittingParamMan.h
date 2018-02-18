// CDSFittingParamMan.h

#ifndef CDSFittingParamMan_h
#define CDSFittingParamMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"

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

 public:
  std::string GetFileName() { return FileName; }
  int GetMaxCDCHit() { return MAXCDCHIT; }
  double GetMagneticField() { return MAGFIELD; }
  int GetMaxChi() { return MAXCHI; }

  int GetUpperLimitCDCTDC() { return CDCTDC_U; }
  int GetLowerLimitCDCTDC() { return CDCTDC_L; }

  void SetMagneticField(double f){ MAGFIELD = f;}

  void Clear();
  
  ClassDef( CDSFittingParamMan, 1 );
};

#endif
