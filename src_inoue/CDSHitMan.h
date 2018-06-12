#ifndef CDSHitMan_h
#define CDSHitMan_h 1

#include <vector>
#include <map>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "HitMan.h"

class CDSHitMan : public HitMan
{
 public:
  CDSHitMan();
  virtual ~CDSHitMan() {};

  HodoscopeLikeContainer *HodoContainer(const int &cid);
  CDCLikeContainer   *CDCHitContainer(const int &cid, const int &layer);
  
 private:
  HodoscopeLikeContainer CDHContainer;
  HodoscopeLikeContainer IHContainer;

  CDCLikeContainer CDCContainer[NumOfCDCLayers];

 public:
  // CDH
  int nCDH() const { return CDHContainer.size(); }
  HodoscopeLikeHit *CDH( const int &i ) { return &CDHContainer[i]; }

  // IH
  int nIH() const { return IHContainer.size(); }
  HodoscopeLikeHit *IH( const int &i ) { return &IHContainer[i]; }

  // CDC
  int nCDC();
  int nCDC( const int &layer ) { return (0<layer&&layer<=NumOfCDCLayers) ? CDCContainer[layer-1].size() : 0; }
  CDCHit *CDC( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfCDCLayers) ? &CDCContainer[layer-1][i] : 0; }

 public:
  bool Calc( ConfMan *conf, double toffs=0 );
  void RemoveNoTDCData();

  void Clear();
  void CheckContainerSize();
  
  // simulation
  bool SetSimulatedResolution( ConfMan *conf );

  ClassDef( CDSHitMan, 1 );
};

#endif

