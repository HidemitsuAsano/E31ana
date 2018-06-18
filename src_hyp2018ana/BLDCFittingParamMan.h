// BLDCFittingParamMan.h

#ifndef BLDCFittingParamMan_h
#define BLDCFittingParamMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"

class BLDCFittingParamMan : public TObject
{
 public:
  BLDCFittingParamMan();
  BLDCFittingParamMan( const std::string & filename );
  BLDCFittingParamMan( const BLDCFittingParamMan &right );
  ~BLDCFittingParamMan();

  void SetFileName( const std::string & filename );
  bool Initialize();

 private:

  std::string FileName;
  int MAXBLDCHIT;
  int MAXHITinLAYER;
  int MAXHITinTRACK;
  double MAGFIELD;
  int MAXCHI;
  int MAXCHIPRE;
  double MAXSLOPE;
  int BLDCTDC_U;
  int BLDCTDC_L;

  int MinHitXBLC1;
  int MinHitXBLC1a;
  int MinHitXBLC1b;
  int MinHitXBLC2;
  int MinHitXBLC2a;
  int MinHitXBLC2b;
  int MinHitXBPC;
  int MinHitXFDC1;

  int MinHitYBLC1;
  int MinHitYBLC1a;
  int MinHitYBLC1b;
  int MinHitYBLC2;
  int MinHitYBLC2a;
  int MinHitYBLC2b;
  int MinHitYBPC;
  int MinHitYFDC1;

  bool LayerDeadBLC1a[8];
  bool LayerDeadBLC1b[8];
  bool LayerDeadBLC2a[8];
  bool LayerDeadBLC2b[8];
  bool LayerDeadBPC[8];
  bool LayerDeadFDC1[8];
  bool LayerKilledBLC1a[8];
  bool LayerKilledBLC1b[8];
  bool LayerKilledBLC2a[8];
  bool LayerKilledBLC2b[8];
  bool LayerKilledBPC[8];
  bool LayerKilledFDC1[8];

 public:
  std::string GetFileName() { return FileName; }
  int GetMaxBLDCHit() { return MAXBLDCHIT; }
  int GetMaxHitInLayer() { return MAXHITinLAYER; }
  int GetMaxHitInTrack() { return MAXHITinTRACK; }
  double GetMagneticField() { return MAGFIELD; }
  int GetMaxChi() { return MAXCHI; }
  int GetMaxChiMWPCFit() { return MAXCHIPRE; }
  double GetMaxSlope() { return MAXSLOPE; }

  int GetMinHitXBLC1()  { return MinHitXBLC1; }
  int GetMinHitXBLC1a() { return MinHitXBLC1a; }
  int GetMinHitXBLC1b() { return MinHitXBLC1b; }
  int GetMinHitXBLC2()  { return MinHitXBLC2; }
  int GetMinHitXBLC2a() { return MinHitXBLC2a; }
  int GetMinHitXBLC2b() { return MinHitXBLC2b; }
  int GetMinHitXBPC() { return MinHitXBPC; }
  int GetMinHitXFDC1() { return MinHitXFDC1; }
  int GetMinHitX(const int &cid);

  int GetMinHit( const int &xy, const int &cid);

  int GetMinHitYBLC1()  { return MinHitYBLC1; }
  int GetMinHitYBLC1a() { return MinHitYBLC1a; }
  int GetMinHitYBLC1b() { return MinHitYBLC1b; }
  int GetMinHitYBLC2()  { return MinHitYBLC2; }
  int GetMinHitYBLC2a() { return MinHitYBLC2a; }
  int GetMinHitYBLC2b() { return MinHitYBLC2b; }
  int GetMinHitYBPC() { return MinHitYBPC; }
  int GetMinHitYFDC1() { return MinHitYFDC1; }
  int GetMinHitY(const int &cid);

  int GetUpperLimitBLDCTDC() { return BLDCTDC_U; }
  int GetLowerLimitBLDCTDC() { return BLDCTDC_L; }

  bool layerdead( const int &id, const int &i );
  bool layerdeadBLC1a( const int &i ) { return LayerDeadBLC1a[i-1]; }
  bool layerdeadBLC1b( const int &i ) { return LayerDeadBLC1b[i-1]; }
  bool layerdeadBLC2a( const int &i ) { return LayerDeadBLC2a[i-1]; }
  bool layerdeadBLC2b( const int &i ) { return LayerDeadBLC2b[i-1]; }
  bool layerdeadBPC( const int &i ) { return LayerDeadBPC[i-1]; }
  bool layerdeadFDC1( const int &i ) { return LayerDeadFDC1[i-1]; }

  bool layerkilled( const int &id, const int &i );
  bool layerkilledBLC1a( const int &i ) { return LayerKilledBLC1a[i-1]; }
  bool layerkilledBLC1b( const int &i ) { return LayerKilledBLC1b[i-1]; }
  bool layerkilledBLC2a( const int &i ) { return LayerKilledBLC2a[i-1]; }
  bool layerkilledBLC2b( const int &i ) { return LayerKilledBLC2b[i-1]; }
  bool layerkilledBPC( const int &i ) { return LayerKilledBPC[i-1]; }
  bool layerkilledFDC1( const int &i ) { return LayerKilledFDC1[i-1]; }

  void Clear();
  
  ClassDef( BLDCFittingParamMan, 1 );
};

#endif
