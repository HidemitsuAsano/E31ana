/* DataContainer.h */
#ifndef DataContainer_h
#define DataContainer_h 1

#include "GlobalVariables.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

class DataContainer : public TObject
{
 public:
  DataContainer();
  ~DataContainer();

  void SetHodoscope(BeamLineHitMan *blMan, CDSHitMan *cdsMan);
  bool DoTrackingBL(BeamLineHitMan *blMan, ConfMan *confMan);
  bool ExecuteCDS(CDSHitMan *cdsMan, ConfMan *confMan);
  void CalcVertex_beam(ConfMan *confMan);

 private:
  CDSTrackingMan   trackingMan;
  BeamLineTrackMan bltrackMan;

  typedef std::vector <HodoscopeLikeHit> HodoscopeLikeContainer;
  // Beam Line
  HodoscopeLikeContainer BHDContainer;
  HodoscopeLikeContainer BHDpostContainer;
  HodoscopeLikeContainer T0Container;
  HodoscopeLikeContainer T0preContainer;
  HodoscopeLikeContainer T0postContainer;
  HodoscopeLikeContainer E0Container;
  HodoscopeLikeContainer CVCContainer;
  HodoscopeLikeContainer BPDContainer;
  HodoscopeLikeContainer BVCContainer;
  HodoscopeLikeContainer PCContainer;
  HodoscopeLikeContainer BVContainer;
  HodoscopeLikeContainer BDContainer;
  HodoscopeLikeContainer LBContainer;
  HodoscopeLikeContainer WVCContainer;
  HodoscopeLikeContainer HVC1Container;
  HodoscopeLikeContainer HVC2Container;
  HodoscopeLikeContainer Temp1Container;
  HodoscopeLikeContainer Temp2Container;

  // Neutron Counter
  static const int NumOfNCLayer=7;
  static const int NumOfNCSegmentsInLayer=16;
  HodoscopeLikeContainer NCContainer[NumOfNCLayer];

  // CDS 
  HodoscopeLikeContainer CDHContainer;
  HodoscopeLikeContainer IHContainer;

  typedef std::vector <CherenkovLikeHit> CherenkovLikeContainer;
  CherenkovLikeContainer LC1Container;
  CherenkovLikeContainer LC2Container;
  CherenkovLikeContainer ACContainer;
  CherenkovLikeContainer WCContainer;
  CherenkovLikeContainer GCContainer;

 public:
  // for Tracking Manager
  CDSTrackingMan *GetCDSTrackingMan(){ return &trackingMan; };
  BeamLineTrackMan *GetBLTrackMan(){ return &bltrackMan; };

  // for BeamLine Hodoscopes
  int  nBHD() const { return BHDContainer.size(); }
  HodoscopeLikeHit *BHD( const int &i ) { return &BHDContainer[i]; }

  int  nBHDpost() const { return BHDpostContainer.size(); }
  HodoscopeLikeHit *BHDpost( const int &i ) { return &BHDpostContainer[i]; }

  int  nT0() const { return T0Container.size(); }
  HodoscopeLikeHit *T0( const int &i ) { return &T0Container[i]; }

  int  nT0pre() const { return T0preContainer.size(); }
  HodoscopeLikeHit *T0pre( const int &i ) { return &T0preContainer[i]; }

  int  nT0post() const { return T0postContainer.size(); }
  HodoscopeLikeHit *T0post( const int &i ) { return &T0postContainer[i]; }

  int  nBPD() const { return BPDContainer.size(); }
  HodoscopeLikeHit *BPD( const int &i ) { return &BPDContainer[i]; }

  int  nE0() const { return E0Container.size(); }
  HodoscopeLikeHit *E0( const int &i ) { return &E0Container[i]; }
  int  nDEF() const { return E0Container.size(); }
  HodoscopeLikeHit *DEF( const int &i ) { return &E0Container[i]; }

  int  nTOF() const { return CVCContainer.size(); }
  int  nCVC() const { return CVCContainer.size(); }
  HodoscopeLikeHit *TOF( const int &i ) { return &CVCContainer[i]; }
  HodoscopeLikeHit *CVC( const int &i ) { return &CVCContainer[i]; }

  int  nPC() const { return PCContainer.size(); }
  HodoscopeLikeHit *PC( const int &i ) { return &PCContainer[i]; }

  int  nBVC() const { return BVContainer.size(); }
  HodoscopeLikeHit *BVC( const int &i ) { return &BVContainer[i]; }

  int  nWVC() const { return WVCContainer.size(); }
  HodoscopeLikeHit *WVC( const int &i ) { return &WVCContainer[i]; }

  int  nHVC1() const { return HVC1Container.size(); }
  HodoscopeLikeHit *HVC1( const int &i ) { return &HVC1Container[i]; }

  int  nHVC2() const { return HVC2Container.size(); }
  HodoscopeLikeHit *HVC2( const int &i ) { return &HVC2Container[i]; }

  int  nLongbar() const { return LBContainer.size(); }
  int  nLB() const { return LBContainer.size(); }
  HodoscopeLikeHit *Longbar( const int &i ) { return &LBContainer[i]; }
  HodoscopeLikeHit *LB( const int &i ) { return &LBContainer[i]; }

  int  nBeamdump() const { return BDContainer.size(); }
  int  nBD() const { return BDContainer.size(); }
  HodoscopeLikeHit *Beamdump( const int &i ) { return &BDContainer[i]; }
  HodoscopeLikeHit *BD( const int &i ) { return &BDContainer[i]; }

  int  nTemp1() const { return Temp1Container.size(); }
  HodoscopeLikeHit *Temp1( const int &i ) { return &Temp1Container[i]; }

  int  nTemp2() const { return Temp2Container.size(); }
  HodoscopeLikeHit *Temp2( const int &i ) { return &Temp2Container[i]; }

  int nNC() const;
  int nNC(const int &layer) const { return (0<layer&&layer<=NumOfNCLayer) ? NCContainer[layer-1].size() : 0; }
  HodoscopeLikeHit *NC(const int &i);
  HodoscopeLikeHit *NC(const int &layer, const int &i){ return (0<layer&&layer<NumOfNCLayer) ? &NCContainer[layer-1][i] : 0; }

  // for CDS Hodoscopes
  int nCDH() const { return CDHContainer.size(); };
  HodoscopeLikeHit *CDH(const int &i){ return &CDHContainer[i]; };

  int nIH() const { return IHContainer.size(); };
  HodoscopeLikeHit *IH(const int &i){ return &IHContainer[i]; };

  // for Cherenlov Counter
  int  nLC1() const { return LC1Container.size(); }
  CherenkovLikeHit *LC1( const int &i ) { return &LC1Container[i]; }

  int  nLC2() const { return LC2Container.size(); }
  CherenkovLikeHit *LC2( const int &i ) { return &LC2Container[i]; }

  int  nAC() const { return ACContainer.size(); }
  CherenkovLikeHit *AC( const int &i ) { return &ACContainer[i]; }

  int  nWC() const { return WCContainer.size(); }
  CherenkovLikeHit *WC( const int &i ) { return &WCContainer[i]; }

  int  nGC() const { return GCContainer.size(); }
  CherenkovLikeHit *GC( const int &i ) { return &GCContainer[i]; }

  void Clear();

  ClassDef(DataContainer, 1)
};

#endif
