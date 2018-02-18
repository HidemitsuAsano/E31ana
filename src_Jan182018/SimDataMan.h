#ifndef SimDataMan_h
#define SimDataMan_h 1

#include <vector>
#include <string>

#include "TObject.h"

#include "ConfMan.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "HodoscopeLikeHit.h"
#include "ChamberLikeHit.h"
#include "CDCHit.h"

class SimDataMan : public TObject
{
 public:
  SimDataMan();
  SimDataMan( ConfMan *conf );
  virtual ~SimDataMan() {delete CDSHit; delete BeamLineHit;};

 private:
  ConfMan *Conf;
  CDSHitMan *CDSHit;
  BeamLineHitMan *BeamLineHit;

 public:
  CDSHitMan *GetCDSHit() { return CDSHit; }
  BeamLineHitMan *GetBeamLineHit() { return BeamLineHit; }
  void SetChamberLikeHit( int cid,int layer, int wire, double dl, double toffset=0 );
  void SetHodoscopeLikeHit(int cid, int seg, double time, double dene, double hitpos );


  void SetCDCHit( int layer, int wire, double dl, double toffset=0 );
  void SetCDHHit( int seg, double time, double dene, double hitpos );

  void SetBLCHit( int tag, int layer, int wire, double dl, double toffset=0 );
                   // tag: BLC1=1, BLC2=0
  void SetT0Hit( int seg, double time, double dene, double hitpos );

  void Clear();

  ClassDef( SimDataMan, 1 );
};

#endif

