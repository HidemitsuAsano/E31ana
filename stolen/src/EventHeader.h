// EventHeader.h

#ifndef EVENTHEADER_H
#define EVENTHEADER_H

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include <fstream>
#include "TKO.h"
#include "ConfMan.h"
#include "TTree.h"
//#include "vme_common.h"
//#include "ana_common.h"

class EventHeader : public TObject
{
 public:
  EventHeader();
  virtual ~EventHeader() {};

 private:
  int RunNumber;
  int EventNumber;
  int VEventNumber;
  int BlockEventNumber;

  short TriggerPattern[20];
  double TriggerTime[20];

 public:
  int run()  const { return RunNumber; }
  int ev()   const { return EventNumber; }
  int vev()   const { return VEventNumber; }
  int blev() const { return BlockEventNumber; }

  int trigmode(    ConfMan* conf=0);
  bool trigmode( const int &trig,    ConfMan* conf=0);
  bool trigmode2( const int &trig,    ConfMan* conf=0);
  int trigparticle(ConfMan* conf=0);
  int trigncdh(    ConfMan* conf=0);

  int pattern( const int &i ) { return ((i<20) ? (int)TriggerPattern[i] : -1); }
  double time( const int &i ) { return ((i<20) ? TriggerTime[i] : -1); }
  bool IsTrig(int id, ConfMan* conf=0);
  bool beam(ConfMan* conf=0);
  bool kaon(ConfMan* conf=0);
  bool pion(ConfMan* conf=0);
  bool proton(ConfMan* conf=0);
  bool electron(ConfMan* conf=0);
  bool cds(ConfMan* conf=0);
  bool sim(ConfMan* conf=0);
  int simid() { return sim() ? pattern(Trig_SIM) : -1 ; }

  void SetRunNumber( const int &run ) { RunNumber = run; }
  void SetEventNumber( const int &evnum ) { EventNumber = evnum; }
  void SetVEventNumber( const int &evnum ) { VEventNumber = evnum; }
  void SetBlockEventNumber( const int &evnum ) { BlockEventNumber = evnum; }
  void SetTriggerPattern( const int &i, const int &tdc ) { TriggerPattern[i] = tdc; }
  void SetSimReacID( int i ) { TriggerPattern[Trig_SIM] = i; }
  
 private:
  int Utime;
  short Spillv;
  short Evidv;
  short Spillt[2];
  short Evidt[2];
  short SmpSw;
  short TkoSw;
  short AdcCnt;
  bool VMESYNC;

 public:
  int utime() const { return Utime; }
  int spillv() const { return Spillv; }
  int evidv() const { return Evidv; }
  int spillt(const int &i) const { return Spillt[i]; }
  int evidt(const int &i) const { return Evidt[i]; }
  bool sync() const { return VMESYNC; }
  void sync_set(bool flag){VMESYNC = flag;}

  void SetUtime( const int &time ){ Utime=time; }
  void Convert( TKOHitCollection *tko, ConfMan *conf );
  /* void SetVMEinfo( EventStruct &vevent ); */
  /* bool EventMatch(TTree *vme, EventStruct &veve, ConfMan *conf); */
  /* bool EventMatch(TTree *vme, EventStruct &veve, ConfMan *conf, std::ofstream &fout); */
  //  int EventMatch(TTree *vme, int &spill, int &evid);

  void Initialize();
  void Clear();

  ClassDef( EventHeader, 1 );
};

#endif

