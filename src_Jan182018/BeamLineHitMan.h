#ifndef BeamLineHitMan_h
#define BeamLineHitMan_h 1

#include <vector>
#include <map>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TKO.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "GlobalVariables.h"
#include "ConfMan.h"

class BeamLineHitMan : public TObject
{
 public:
  BeamLineHitMan();
  virtual ~BeamLineHitMan() {};

 private:
  typedef std::vector <HodoscopeLikeHit> HodoscopeLikeContainer;
  HodoscopeLikeContainer BHDContainer;  
  HodoscopeLikeContainer BHDpostContainer;  
  HodoscopeLikeContainer PAContainer;   
  HodoscopeLikeContainer T0Container;   
  HodoscopeLikeContainer T0preContainer;   
  HodoscopeLikeContainer T0postContainer;   
  HodoscopeLikeContainer E0Container;   
  HodoscopeLikeContainer B1Container;   
  HodoscopeLikeContainer B2Container;   
  HodoscopeLikeContainer RangeContainer;
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
  HodoscopeLikeContainer HVC3Container;  
  HodoscopeLikeContainer Temp1Container;  
  HodoscopeLikeContainer Temp2Container;  
  HodoscopeLikeContainer Temp3Container;  

  typedef std::vector <CherenkovLikeHit> CherenkovLikeContainer;
  CherenkovLikeContainer LC1Container;  
  CherenkovLikeContainer LC2Container;  
  CherenkovLikeContainer ACContainer;   
  CherenkovLikeContainer WCContainer;   
  CherenkovLikeContainer GCContainer;   

  typedef std::vector <ChamberLikeHit> ChamberLikeContainer;
  ChamberLikeContainer   BLC1aContainer[NumOfBLCLayers];
  ChamberLikeContainer   BLC1bContainer[NumOfBLCLayers];
  ChamberLikeContainer   BLC2aContainer[NumOfBLCLayers];
  ChamberLikeContainer   BLC2bContainer[NumOfBLCLayers];
  ChamberLikeContainer   BPCContainer[NumOfBPCLayers];
  ChamberLikeContainer   FDC1Container[NumOfFDC1Layers];

  typedef std::map <int,int> SegmentContainer;
  SegmentContainer HodoscopeSegContainer;
  SegmentContainer CherenkovSegContainer;

 public:

  // BHD
  int  nBHD() const { return BHDContainer.size(); }
  HodoscopeLikeHit *BHD( const int &i ) { return &BHDContainer[i]; }

  // BHD
  int  nBHDpost() const { return BHDpostContainer.size(); }
  HodoscopeLikeHit *BHDpost( const int &i ) { return &BHDpostContainer[i]; }

  // PA
  int  nPA() const { return PAContainer.size(); }
  HodoscopeLikeHit *PA( const int &i ) { return &PAContainer[i]; }

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

  // B1
  int  nB1() const { return B1Container.size(); }
  HodoscopeLikeHit *B1( const int &i ) { return &B1Container[i]; }

  // B2
  int  nB2() const { return B2Container.size(); }
  HodoscopeLikeHit *B2( const int &i ) { return &B2Container[i]; }

  // Range
  int  nRange() const { return RangeContainer.size(); }
  HodoscopeLikeHit *Range( const int &i ) { return &RangeContainer[i]; }

  // TOF
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

  // CV counter
  //  int  nCV() const { return CVContainer.size(); }
  //  HodoscopeLikeHit *CV( const int &i ) { return &CVContainer[i]; }

  // BVC counter
  int  nBVC() const { return BVContainer.size(); }
  HodoscopeLikeHit *BVC( const int &i ) { return &BVContainer[i]; }

  // WVC counter
  int  nWVC() const { return WVContainer.size(); }
  HodoscopeLikeHit *WVC( const int &i ) { return &WVContainer[i]; }

  // HVC1 counter
  int  nHVC1() const { return HVC1Container.size(); }
  HodoscopeLikeHit *HVC1( const int &i ) { return &HVC1Container[i]; }
  // HVC2 counter
  int  nHVC2() const { return HVC2Container.size(); }
  HodoscopeLikeHit *HVC2( const int &i ) { return &HVC2Container[i]; }
  // HVC2 counter
  int  nHVC3() const { return HVC3Container.size(); }
  HodoscopeLikeHit *HVC3( const int &i ) { return &HVC3Container[i]; }

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
  HodoscopeLikeHit *Temp1( const int &i ) { return &Temp1Container[i]; }

  int  nTemp2() const { return Temp2Container.size(); }
  HodoscopeLikeHit *Temp2( const int &i ) { return &Temp2Container[i]; }

  int  nTemp3() const { return Temp3Container.size(); }
  HodoscopeLikeHit *Temp3( const int &i ) { return &Temp3Container[i]; }


  // LC1
  int  nLC1() const { return LC1Container.size(); }
  CherenkovLikeHit *LC1( const int &i ) { return &LC1Container[i]; }

  // LC2
  int  nLC2() const { return LC2Container.size(); }
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

  int nHodo(const int &cid);
  HodoscopeLikeHit *Hodoi( const int &cid, const unsigned int &i );

  int  nChere(const int &cid);
  CherenkovLikeHit *Cherei( const int &cid, const unsigned int &i );

  HodoscopeLikeHit *Hodo( const int &cid, const int &seg );
  CherenkovLikeHit *Chere( const int &cid, const int &seg );

  int nBLDC( const int &cid, const int &layer);
  ChamberLikeHit *BLDC( const int &id, const int &layer, const int &i);

  // BLC1
  int nBLC1a();
  int nBLC1a( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers) ? BLC1aContainer[layer-1].size() : 0; }
  ChamberLikeHit *BLC1a( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BLC1aContainer[layer-1][i] : 0; }
  int nBLC1b();
  int nBLC1b( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers) ? BLC1bContainer[layer-1].size() : 0; }
  ChamberLikeHit *BLC1b( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BLC1bContainer[layer-1][i] : 0; }

  // BLC2
  int nBLC2a();
  int nBLC2a( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers) ? BLC2aContainer[layer-1].size() : 0; }
  ChamberLikeHit *BLC2a( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BLC2aContainer[layer-1][i] : 0; }
  int nBLC2b();
  int nBLC2b( const int &layer ) { return (0<layer&&layer<=NumOfBLCLayers) ? BLC2bContainer[layer-1].size() : 0; }
  ChamberLikeHit *BLC2b( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBLCLayers) ? &BLC2bContainer[layer-1][i] : 0; }

  // BPC
  int nBPC();
  int nBPC( const int &layer ) { return (0<layer&&layer<=NumOfBPCLayers ) ? BPCContainer[layer-1].size() : 0; }
  ChamberLikeHit *BPC( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfBPCLayers) ? &BPCContainer[layer-1][i] : 0; }

  // FDC1
  int nFDC1();
  int nFDC1( const int &layer ) { return (0<layer&&layer<=NumOfFDC1Layers ) ? FDC1Container[layer-1].size() : 0; }
  ChamberLikeHit *FDC1( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfFDC1Layers) ? &FDC1Container[layer-1][i] : 0; }


 public:

  void AddHit( const HodoscopeLikeHit &hit );
  void AddHit( const CherenkovLikeHit &hit );
  void AddHit( const ChamberLikeHit &hit );

  bool Convert( TKOHitCollection *tko, ConfMan *conf );
  bool Convert( const int &c, const int &n, const int &a, const int &data, ConfMan *conf );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, HodoscopeLikeContainer *container );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, CherenkovLikeContainer *container );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &lay, const int &wire,
		const int &data, ChamberLikeContainer *container );

  bool Calc( ConfMan *conf );
  void CalcBLDC( ConfMan *conf, const double &toffs=0. );

  bool FindTrackPDC();

 private:
  double T0Time;
  int T0Seg;

 public:
  void CalcT0Time();

  int t0seg() const { return T0Seg; }
  void Clear();
  void CheckContainerSize();
  
  // simulation
  bool SetSimulatedResolution( ConfMan *conf );

  ClassDef( BeamLineHitMan, 1 );
};

#endif
