#ifndef SDDHitMan_h
#define SDDHitMan_h 1

#include <vector>
#include <map>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TKO.h"
#include "SDDHit.h"
#include "SDDHitCom.h"
#include "HodoscopeLikeHit.h"
#include "GlobalVariables.h"
#include "ConfMan.h"

class SDDHitMan : public TObject
{
 public:
  SDDHitMan();
  virtual ~SDDHitMan() {};

 private:
  typedef std::vector <SDDHit> SDDHitContainer;
  SDDHitContainer SDDContainer;
  typedef std::vector <SDDHitCom> SDDComContainer;
  SDDComContainer SDDCommon;
  //  typedef std::map <int,int> SegmentContainer;
  //  SegmentContainer SDDSegContainer;

  //  typedef std::vector <HodoscopeLikeHit> DefineContainer;
  //  DefineContainer Define;


 public:
  // SDD
  int  nSDD() const { return SDDContainer.size(); }
  SDDHit *SDD( const int &i ) { return &SDDContainer[i]; }
  SDDHitCom *SDDCom() { return &SDDCommon[0]; }
  //  HodoscopeLikeHit *DEF() { return &Define[0]; }

 private:
  int ADef,TDef;

 public:
  int adef() const { return ADef; }
  int tdef() const { return TDef; }

  void AddHit( const SDDHit &hit );
  void AddHit( const SDDHitCom &hit );
  //  void AddHit( const HodoscopeLikeHit &hit );

  int nHit(int threshold=700);
  int reset(int nth);
  int highth(int nth);

  bool Convert( TKOHitCollection *tko, ConfMan *conf );
  bool Convert( const int &c, const int &n, const int &a, const int &data, ConfMan *conf );
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, SDDHitContainer *container );
  /*
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, DefineContainer *container );
  */
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data, SDDComContainer *container );
  bool SetData( const int &seg, EventStruct &vme, SDDHitContainer *container );
  bool SetData( EventStruct &vme);
  //  bool SetData( const int &c, const int &n, const int &a,
  //		const int &cid, const int &seg, const int &at, const int &ud,
  //		const int &data);

  bool Calc( ConfMan *conf );

  void Clear();
  void CheckContainerSize();
  
  ClassDef( SDDHitMan, 1 );
};

#endif
