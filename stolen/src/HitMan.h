#ifndef HitMan_h
#define HitMan_h 1

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
#include "MTDCHit.h"
#include "CDCHit.h"
#include "GlobalVariables.h"
#include "ConfMan.h"

typedef std::vector <HodoscopeLikeHit> HodoscopeLikeContainer;
typedef std::vector <CherenkovLikeHit> CherenkovLikeContainer;
typedef std::vector <ChamberLikeHit>   ChamberLikeContainer;
typedef std::vector <CDCHit>  CDCLikeContainer;
typedef std::vector <MTDCHit> MTDCLikeContainer;  
typedef std::map <int,int> SegmentContainer;

class HitMan : public TObject
{
 public:
  HitMan();
  virtual ~HitMan() {}
  
 protected:
  SegmentContainer HodoscopeSegContainer;
  SegmentContainer CherenkovSegContainer;

 public:
  virtual HodoscopeLikeContainer *HodoContainer(const int &cid){ return 0; }
  virtual ChamberLikeContainer   *ChmContainer(const int &cid, const int &layer){ return 0; }
  virtual CDCLikeContainer       *CDCHitContainer(const int &cid, const int &layer){ return 0; }
  virtual CherenkovLikeContainer *ChereContainer(const int &cid){ return 0; }
  virtual MTDCLikeContainer      *MTDCContainer(const int &cid){ return 0; }

  
  int nHodo(const int &cid);
  HodoscopeLikeHit *Hodo( const int &cid, const int &seg );
  HodoscopeLikeHit *Hodoi( const int &cid, const unsigned int &i );
  
  //  int nChm( const int &cid);
  int nChm( const int &cid, const int &layer );
  ChamberLikeHit *Chm( const int &cid,const int &layer, const int &i );

  int  nChere(const int &cid);
  CherenkovLikeHit *Cherei( const int &cid, const unsigned int &i );
  CherenkovLikeHit *Chere( const int &cid, const int &seg );

  //  int nMTDC( const int &cid);
  //  MTDCHit *MTDC( const int &id, const int &i);
  
 public:
  void AddHit( const HodoscopeLikeHit &hit );
  void AddHit( const CDCHit &hit );
  void AddHit( const ChamberLikeHit &hit );
  void AddHit( const CherenkovLikeHit &hit );
  void AddHit( const MTDCHit &hit );
  
  bool Convert( TKOHitCollection *tko, ConfMan *conf, double toffs=0. );
  bool Convert( const int &c, const int &n, const int &a, const int &data, ConfMan *conf );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, HodoscopeLikeContainer *container );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, CherenkovLikeContainer *container );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, MTDCLikeContainer *container );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &lay, const int &wire,
		const int &data, ChamberLikeContainer *container );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &lay, const int &wire,
		const int &data, CDCLikeContainer *container );

  virtual bool Calc( ConfMan *conf, double toffs )=0;
  virtual void RemoveNoTDCData(){}
  virtual void Clear();
  virtual void CheckContainerSize(){}
  
  // simulation
  virtual bool SetSimulatedResolution( ConfMan *conf ){ return true; }

 private:
  static const unsigned int KEYMASK  = 0x000F;
  static const unsigned int SEGMASK    = 0x01FF;      /* SEG Mask 9 Bits (0-511) */
  static const unsigned int CIDMASK    = 0x007F;      /* CID Mask 7 Bits (0-127) */
  static const int          SEGSHIFT   = 8;
  static const int          CIDSHIFT   = 24;
  static const unsigned int KEYFLAG  = 0x0003;
  
 protected:
  const int KEY(const int &cid,const int &seg){ return ((((cid)&CIDMASK)<<CIDSHIFT) | (((seg)&SEGMASK)<<SEGSHIFT) | KEYFLAG ); }

  ClassDef( HitMan, 1 );
};

#endif

