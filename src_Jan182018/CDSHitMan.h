#ifndef CDSHitMan_h
#define CDSHitMan_h 1

#include <vector>
#include <map>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TKO.h"
#include "HodoscopeLikeHit.h"
#include "CDCHit.h"
#include "ConfMan.h"

class CDSHitMan : public TObject
{
 public:
  CDSHitMan();
  virtual ~CDSHitMan() {};

 private:
  typedef std::vector <HodoscopeLikeHit> HodoscopeLikeContainer;
  HodoscopeLikeContainer CDHContainer;
  HodoscopeLikeContainer IHContainer;

  typedef std::vector <CDCHit> CDCHitContainer;
  CDCHitContainer CDCContainer[NumOfCDCLayers];

  typedef std::map <int,int> SegmentContainer;
  SegmentContainer HodoscopeSegContainer;

 public:
  // CDH
  int nCDH() const { return CDHContainer.size(); }
  HodoscopeLikeHit *CDH( const int &i ) { return &CDHContainer[i]; }

  // IH
  int nIH() const { return IHContainer.size(); }
  HodoscopeLikeHit *IH( const int &i ) { return &IHContainer[i]; }

  int nHodo(const int &cid);
  HodoscopeLikeHit *Hodo( const int &cid, const int &seg );
  HodoscopeLikeHit *Hodoi( const int &cid, const unsigned int &i );

  // CDC
  int nCDC();
  int nCDC( const int &layer ) { return (0<layer&&layer<=NumOfCDCLayers) ? CDCContainer[layer-1].size() : 0; }
  CDCHit *CDC( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfCDCLayers) ? &CDCContainer[layer-1][i] : 0; }

 public:
  void AddHit( const HodoscopeLikeHit &hit );
  void AddHit( const CDCHit &hit );
  
  bool Convert( TKOHitCollection *tko, ConfMan *conf );
  bool Convert( const int &c, const int &n, const int &a, const int &data, ConfMan *conf );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, HodoscopeLikeContainer *container );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &lay, const int &wire,
		const int &data, CDCHitContainer *container );

  bool Calc( ConfMan *conf );
  bool RemoveNoTDCData();

  void Clear();
  void CheckContainerSize();
  
  // simulation
  bool SetSimulatedResolution( ConfMan *conf );

  ClassDef( CDSHitMan, 1 );
};

#endif

