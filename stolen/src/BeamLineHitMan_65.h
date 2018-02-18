#ifndef BeamLineHitMan_h
#define BeamLineHitMan_h 1

#include <vector>
#include <map>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "HitMan.h"

class BeamLineHitMan : public HitMan
{
 public:
  BeamLineHitMan();
  virtual ~BeamLineHitMan() {};

  HodoscopeLikeContainer *HodoContainer(const int &cid);
  ChamberLikeContainer   *ChmContainer(const int &cid, const int &layer);
  ChamberLikeContainer   *BLDCContainer(const int &cid, const int &layer){ return ChmContainer(cid,layer); }
  CherenkovLikeContainer *ChereContainer(const int &cid);
  MTDCLikeContainer      *MTDCContainer(const int &cid);

 private:
  HodoscopeLikeContainer BHDContainer;  
  HodoscopeLikeContainer BHDpostContainer;  
  HodoscopeLikeContainer T0Container;   
  HodoscopeLikeContainer T0preContainer;   
  HodoscopeLikeContainer T0postContainer;   
  HodoscopeLikeContainer E0Container;   
  HodoscopeLikeContainer CVCContainer;  
  HodoscopeLikeContainer BPDContainer;  
  HodoscopeLikeContainer NCContainer;  
  HodoscopeLikeContainer PCContainer;  
  HodoscopeLikeContainer BVContainer;  
  HodoscopeLikeContainer BDContainer;  
  HodoscopeLikeContainer LBContainer;  
  HodoscopeLikeContainer WVContainer;  
  HodoscopeLikeContainer HVC1Container;  
  HodoscopeLikeContainer HVC2Container;  
  HodoscopeLikeContainer Temp1Container;  
  HodoscopeLikeContainer Temp2Container;  
  HodoscopeLikeContainer StartContainer;  

  MTDCLikeContainer T0MTDCContainer;
  MTDCLikeContainer BHDMTDCContainer;
  MTDCLikeContainer HVC1MTDCContainer;
  MTDCLikeContainer HVC2MTDCContainer;
  MTDCLikeContainer BVCMTDCContainer;
  MTDCLikeContainer REFMTDCContainer;

  CherenkovLikeContainer LC1Container;  
  CherenkovLikeContainer LC2Container;  
  CherenkovLikeContainer ACContainer;   
  CherenkovLikeContainer WCContainer;   
  CherenkovLikeContainer GCContainer;   

  static const int NumOfBLCLayers=8;
  ChamberLikeContainer   BLC1aContainer[NumOfBLCLayers];
  ChamberLikeContainer   BLC1bContainer[NumOfBLCLayers];
  ChamberLikeContainer   BLC2aContainer[NumOfBLCLayers];
  ChamberLikeContainer   BLC2bContainer[NumOfBLCLayers];
  ChamberLikeContainer   BPCContainer[NumOfBLCLayers];
  ChamberLikeContainer   FDC1Container[NumOfBLCLayers];

 public:
  // BHD
  int  nBHD() const { return BHDContainer.size(); }
  HodoscopeLikeHit *BHD( const int &i ) { return &BHDContainer[i]; }

  // BHD
  int  nBHDpost() const { return BHDpostContainer.size(); }
  HodoscopeLikeHit *BHDpost( const int &i ) { return &BHDpostContainer[i]; }

  // T0
  int  nT0() const { return T0Container.size(); }
  HodoscopeLikeHit *T0( const int &i ) { return &T0Container[i]; }

  int  nT0pre() const { return T0preContainer.size(); }
  HodoscopeLikeHit *T0pre( const int &i ) { return &T0preContainer[i]; }

  int  nT0post() const { return T0postContainer.size(); }
  HodoscopeLikeHit *T0post( const int &i ) { return &T0postContainer[i]; }

  // BPD
  int  nBPD() const { return BPDContainer.size(); }
  HodoscopeLikeHit *BPD( const int &i ) { return &BPDContainer[i]; }

  // E0
  int  nE0() const { return E0Container.size(); }
  HodoscopeLikeHit *E0( const int &i ) { return &E0Container[i]; }
  int  nDEF() const { return E0Container.size(); }
  HodoscopeLikeHit *DEF( const int &i ) { return &E0Container[i]; }

  // CVC
  int  nTOF() const { return CVCContainer.size(); }
  int  nCVC() const { return CVCContainer.size(); }
  HodoscopeLikeHit *TOF( const int &i ) { return &CVCContainer[i]; }
  HodoscopeLikeHit *CVC( const int &i ) { return &CVCContainer[i]; }

  // Neutron counter
  int  nNC() const { return NCContainer.size(); }
  HodoscopeLikeHit *NC( const int &i ) { return &NCContainer[i]; }

  // Proton counter
  int  nPC() const { return PCContainer.size(); }
  HodoscopeLikeHit *PC( const int &i ) { return &PCContainer[i]; }

  // BVC counter
  int  nBVC() const { return BVContainer.size(); }
  HodoscopeLikeHit *BVC( const int &i ) { return &BVContainer[i]; }

  // WVC counter
  int  nWVC() const { return WVContainer.size(); }
  HodoscopeLikeHit *WVC( const int &i ) { return &WVContainer[i]; }

  // HVC1 counter
  int  nHVC1() const { return HVC1Container.size(); }
  int  nHVC2() const { return HVC2Container.size(); }
  HodoscopeLikeHit *HVC1( const int &i ) { return &HVC1Container[i]; }
  HodoscopeLikeHit *HVC2( const int &i ) { return &HVC2Container[i]; }

  // longbar counter
  int  nLongbar() const { return LBContainer.size(); }
  int  nLB() const { return LBContainer.size(); }
  HodoscopeLikeHit *Longbar( const int &i ) { return &LBContainer[i]; }
  HodoscopeLikeHit *LB( const int &i ) { return &LBContainer[i]; }

  // Beamdump
  int  nBeamdump() const { return BDContainer.size(); }
  int  nBD() const { return BDContainer.size(); }
  HodoscopeLikeHit *Beamdump( const int &i ) { return &BDContainer[i]; }
  HodoscopeLikeHit *BD( const int &i ) { return &BDContainer[i]; }

  // Temp
  int  nTemp1() const { return Temp1Container.size(); }
  int  nTemp2() const { return Temp2Container.size(); }
  HodoscopeLikeHit *Temp1( const int &i ) { return &Temp1Container[i]; }
  HodoscopeLikeHit *Temp2( const int &i ) { return &Temp2Container[i]; }

  // Start
  int  nStart() const { return StartContainer.size(); }
  HodoscopeLikeHit *Start( const int &i ) { return &StartContainer[i]; }

  // LC1/2
  int  nLC1() const { return LC1Container.size(); }
  int  nLC2() const { return LC2Container.size(); }
  CherenkovLikeHit *LC1( const int &i ) { return &LC1Container[i]; }
  CherenkovLikeHit *LC2( const int &i ) { return &LC2Container[i]; }

  // AC
  int  nAC() const { return ACContainer.size(); }
  CherenkovLikeHit *AC( const int &i ) { return &ACContainer[i]; }

  // WC
  int  nWC() const { return WCContainer.size(); }
  CherenkovLikeHit *WC( const int &i ) { return &WCContainer[i]; }

  // GC
  int  nGC() const { return GCContainer.size(); }
  CherenkovLikeHit *GC( const int &i ) { return &GCContainer[i]; }

  // BLC1
  int nBLC1a( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers) ? BLC1aContainer[layer-1].size() : 0; }
  ChamberLikeHit *BLC1a( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BLC1aContainer[layer-1][i] : 0; }
  int nBLC1b( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers) ? BLC1bContainer[layer-1].size() : 0; }
  ChamberLikeHit *BLC1b( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BLC1bContainer[layer-1][i] : 0; }

  // BLC2
  int nBLC2a( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers) ? BLC2aContainer[layer-1].size() : 0; }
  ChamberLikeHit *BLC2a( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BLC2aContainer[layer-1][i] : 0; }
  int nBLC2b( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers) ? BLC2bContainer[layer-1].size() : 0; }
  ChamberLikeHit *BLC2b( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BLC2bContainer[layer-1][i] : 0; }

  // BPC
  int nBPC( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers ) ? BPCContainer[layer-1].size() : 0; }
  ChamberLikeHit *BPC( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BPCContainer[layer-1][i] : 0; }

  // FDC1
  int nFDC1( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers ) ? FDC1Container[layer-1].size() : 0; }
  ChamberLikeHit *FDC1( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &FDC1Container[layer-1][i] : 0; }

  int nBLDC( const int &cid, const int &layer){ return nChm(cid,layer); }
  ChamberLikeHit *BLDC( const int &id, const int &layer, const int &i){ return Chm(id,layer,i); }

 public:
  bool Calc( ConfMan *conf, double toffs=0.);
  void CalcBLDC( ConfMan *conf, const double &toffs=0. );

 private:
  double T0Time;
  int T0Seg;

 public:
  void CalcT0Time();

  int t0seg() const { return T0Seg; }
  double t0time() const { return T0Time; }
  void Clear();
  void CheckContainerSize();
  void RemoveNoTDCData();

  // simulation
  bool SetSimulatedResolution( ConfMan *conf );

  ClassDef( BeamLineHitMan, 1 );
};

#endif
