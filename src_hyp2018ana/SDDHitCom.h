#ifndef SDDHitCom_h
#define SDDHitCom_h 1

#include <string>
#include <iostream>
#include <cmath>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ConfMan.h"

class SDDHitCom : public TObject
{
 private:
  int CounterID;
  int HitID; 
  
 public:
  SDDHitCom();
  SDDHitCom( const SDDHitCom &hit );
  virtual ~SDDHitCom() {};
  int cid() const { return CounterID; }
  void SetCID( const int &cid ){ CounterID = cid; }

 private:
  int T0t,E0t1,E0t2,E0t3,First1,First2; 
  int CrateT0t, SlotT0t, ChannelT0t;
  int CrateE0t1, SlotE0t1, ChannelE0t1;
  int CrateE0t2, SlotE0t2, ChannelE0t2;
  int CrateE0t3, SlotE0t3, ChannelE0t3;
  int CrateFirst1, SlotFirst1, ChannelFirst1;
  int CrateFirst2, SlotFirst2, ChannelFirst2;

 public:
  int t0t() const { return T0t; }
  int e0t1() const { return E0t1; }
  int e0t2() const { return E0t2; }
  int e0t3() const { return E0t3; }
  int first1() const { return First1; }
  int first2() const { return First2; }

  int crt0t() const { return CrateT0t; }
  int slt0t() const { return SlotT0t; }
  int cht0t() const { return ChannelT0t; }
  int cre0t1() const { return CrateE0t1; }
  int sle0t1() const { return SlotE0t1; }
  int che0t1() const { return ChannelE0t1; }
  int cre0t2() const { return CrateE0t2; }
  int sle0t2() const { return SlotE0t2; }
  int che0t2() const { return ChannelE0t2; }
  int cre0t3() const { return CrateE0t3; }
  int sle0t3() const { return SlotE0t3; }
  int che0t3() const { return ChannelE0t3; }
  int crfirst1() const { return CrateFirst1; }
  int slfirst1() const { return SlotFirst1; }
  int chfirst1() const { return ChannelFirst1; }
  int crfirst2() const { return CrateFirst2; }
  int slfirst2() const { return SlotFirst2; }
  int chfirst2() const { return ChannelFirst2; }

  void SetT0t(       const int & v   ) { T0t = v; }
  void SetE0t1(      const int & v   ) { E0t1 = v; }
  void SetE0t2(      const int & v   ) { E0t2 = v; }
  void SetE0t3(      const int & v   ) { E0t3 = v; }
  void SetFirst1(      const int & v   ) { First1 = v; }
  void SetFirst2(      const int & v   ) { First2 = v; }

  void SetCrateT0t(   const int & cr  ) { CrateT0t = cr; }
  void SetSlotT0t(    const int & sl  ) { SlotT0t = sl; }
  void SetChannelT0t( const int & ch  ) { ChannelT0t = ch; }
  void SetCrateE0t1(   const int & cr  ) { CrateE0t1 = cr; }
  void SetSlotE0t1(    const int & sl  ) { SlotE0t1 = sl; }
  void SetChannelE0t1( const int & ch  ) { ChannelE0t1 = ch; }
  void SetCrateE0t2(   const int & cr  ) { CrateE0t2 = cr; }
  void SetSlotE0t2(    const int & sl  ) { SlotE0t2 = sl; }
  void SetChannelE0t2( const int & ch  ) { ChannelE0t2 = ch; }
  void SetCrateE0t3(   const int & cr  ) { CrateE0t3 = cr; }
  void SetSlotE0t3(    const int & sl  ) { SlotE0t3 = sl; }
  void SetChannelE0t3( const int & ch  ) { ChannelE0t3 = ch; }
  void SetCrateFirst1(   const int & cr  ) { CrateFirst1 = cr; }
  void SetSlotFirst1(    const int & sl  ) { SlotFirst1 = sl; }
  void SetChannelFirst1( const int & ch  ) { ChannelFirst1 = ch; }
  void SetCrateFirst2(   const int & cr  ) { CrateFirst2 = cr; }
  void SetSlotFirst2(    const int & sl  ) { SlotFirst2 = sl; }
  void SetChannelFirst2( const int & ch  ) { ChannelFirst2 = ch; }

  
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data );
  
  //  bool CheckRange();
  
 private:
  // conversion data
  double T0time,E0time1,E0time2,E0time3,Firstime1,Firstime2;
  // corrected time
  double CT0time,CE0time1,CE0time2,CE0time3,CFirstime1,CFirstime2; 

 public:
  double t0time()  const { return T0time; }
  double e0time1() const { return E0time1; }
  double e0time2() const { return E0time2; }
  double e0time3() const { return E0time3; }
  double firstime1() const { return Firstime1; }
  double firstime2() const { return Firstime2; }
  double ct0time()  const { return CT0time; }
  double ce0time1() const { return CE0time1; }
  double ce0time2() const { return CE0time2; }
  double ce0time3() const { return CE0time3; }
  double cfirstime1() const { return CFirstime1; }
  double cfirstime2() const { return CFirstime2; }


  void SetT0time(  const double &t ) { T0time = t; }
  void SetE0time1( const double &t ) { E0time1 = t; }
  void SetE0time2( const double &t ) { E0time2 = t; }
  void SetE0time3( const double &t ) { E0time3 = t; }
  void SetFirstime1( const double &t ) { Firstime1 = t; }
  void SetFirstime2( const double &t ) { Firstime2 = t; }
  void SetCT0time(  const double &t ) { CT0time = t; }
  void SetCE0time1( const double &t ) { CE0time1 = t; }
  void SetCE0time2( const double &t ) { CE0time2 = t; }
  void SetCE0time3( const double &t ) { CE0time3 = t; }
  void SetCFirstime1( const double &t ) { CFirstime1 = t; }
  void SetCFirstime2( const double &t ) { CFirstime2 = t; }

  bool Calc( ConfMan *conf );
  void Clear();

 public:
  //  void Reverse( ConfMan *conf );

  ClassDef(SDDHitCom, 1 );
};

#endif
