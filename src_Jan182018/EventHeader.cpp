// EventHeader.cpp
#include <stdlib.h>
#include "EventHeader.h"
#include "GlobalVariables.h"

ClassImp(EventHeader);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
EventHeader::EventHeader() : TObject()
{
  Initialize();
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int EventHeader::trigmode(ConfMan *conf)
{
  if(IsTrig(Trig_Cosmic))  return Mode_Cosmic;
  if(IsTrig(Trig_Reject))  return Mode_Reject;
  if(IsTrig(Trig_Kf,conf)) return Mode_Kf;
  if(IsTrig(Trig_KCDH1f,conf)) return Mode_KCDH1f;
  if(IsTrig(Trig_KCDH2,conf))  return Mode_KCDH2;
  if(IsTrig(Trig_KCDH1,conf)&&IsTrig(Trig_Neutral,conf)) return Mode_KCDH1N;
  if(IsTrig(Trig_KCDH1,conf)&&IsTrig(Trig_Charged,conf)) return Mode_KCDH1C;
  if(IsTrig(Trig_PivBVC,conf)&&IsTrig(Trig_Neutral,conf)) return Mode_PiN;
  if(IsTrig(Trig_PivBVC,conf)&&IsTrig(Trig_Charged,conf)) return Mode_PiC;
  if(IsTrig(Trig_PiCDH1,conf)) return Mode_PiCDH1;
  return Mode_Unknown;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool EventHeader::trigmode(const int &mode,ConfMan *conf)
{
  if(mode==Mode_Cosmic)return IsTrig(Trig_Cosmic);
  if(mode==Mode_Reject)return IsTrig(Trig_Reject);
  if(mode==Mode_Kf)    return IsTrig(Trig_Kf,conf);
  if(mode==Mode_KCDH1f)return IsTrig(Trig_KCDH1f,conf);
  if(mode==Mode_KCDH2) return IsTrig(Trig_KCDH2,conf);
  if(mode==Mode_KCDH1N)return IsTrig(Trig_KCDH1,conf)&&IsTrig(Trig_Neutral,conf);
  if(mode==Mode_KCDH1C)return IsTrig(Trig_KCDH1,conf)&&IsTrig(Trig_Charged,conf);
  if(mode==Mode_PiN)   return IsTrig(Trig_PivBVC,conf)&&IsTrig(Trig_Neutral,conf);
  if(mode==Mode_PiC)   return IsTrig(Trig_PivBVC,conf)&&IsTrig(Trig_Charged,conf);
  if(mode==Mode_PiCDH1)return IsTrig(Trig_PiCDH1,conf);
  else return false;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int EventHeader::trigparticle(ConfMan *conf)
{
  if(IsTrig(Trig_Kaon,conf)) return Beam_Kaon;
  if(IsTrig(Trig_Pion,conf)) return Beam_Pion;
  if(IsTrig(Trig_Proton,conf)) return Beam_Proton;
  return Beam_Other;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int EventHeader::trigncdh(ConfMan *conf)
{
  if(IsTrig(Trig_KCDH1,conf))//&&IsTrig(Trig_PiCDH1,conf))
    return 1;
  if(IsTrig(Trig_KCDH2,conf))//&&IsTrig(Trig_PiCDH2,conf))
    return 2;
  else return 0;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool EventHeader::IsTrig(int id, ConfMan* conf)
{
  if( conf )
    if(conf->GetGateMapManager())
      return conf->GetGateMapManager()
	->CheckRange( CID_MISC, id, 0, 1, 0, (double)TriggerPattern[id] );
	
  if( 0<TriggerPattern[id] && TriggerPattern[id]<4095 )
    return true;
  else
    return false;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool EventHeader::beam(ConfMan *conf)
{
  return IsTrig(Trig_Beam,conf);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool EventHeader::kaon(ConfMan *conf)
{
  return IsTrig(Trig_Kaon,conf);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool EventHeader::pion(ConfMan *conf)
{
  return IsTrig(Trig_Pion,conf);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool EventHeader::proton(ConfMan *conf)
{
  return IsTrig(Trig_Proton,conf);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool EventHeader::electron(ConfMan *conf)
{
  return IsTrig(Trig_Electron,conf);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool EventHeader::cds(ConfMan *conf)
{
  return false;//IsTrig(Trig_CDS,conf);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void EventHeader::Convert( TKOHitCollection *tko, ConfMan *conf )
{
#if 0
  std::cout << " Enter EventHeader::Convert(TKO,conf) " << std::endl;
#endif
  for( int ih=0; ih<tko->entries(); ih++ ){
    int c, n, a, data;
    int tmp;
    int cid, seg, at, ud;
    c = tko->hit(ih)->cr(); n = tko->hit(ih)->sl(); a = tko->hit(ih)->ch(); data = tko->hit(ih)->data();
#if 0
    std::cout << " c:" << c << " n:" << n << " a:" << a << " data:" << data << std::endl;
#endif
    conf->GetCounterMapManager()->GetInfo( c, n, a, cid, tmp, seg, at, ud );
    //    if( cid != CID_MISC ) continue;
    switch( cid ){
    case CID_MISC:
      if(seg<20&&seg>=0){
	TriggerPattern[seg] = data;
	if( conf->GetGainMapManager() )
	  if( 0<data ) 
	    TriggerTime[seg] = conf->GetGainMapManager()->CalcCValue( c, n, a, 1, data );
      }
#if 0
      std::cout << " c:" << c << " n:" << n << " a:" << a << " data:" << data << " seg:"<<seg<<std::endl;
#endif
      /*
	if( seg == Trig_Beam )     TriggerPattern[seg] = data;
	if( seg == Trig_Kaon )     TriggerPattern[seg] = data;
	if( seg == Trig_Pion )     TriggerPattern[seg] = data;
	if( seg == Trig_Proton )   TriggerPattern[seg] = data;
	if( seg == Trig_Electron ) TriggerPattern[seg] = data;
	if( seg == Trig_CDS )      TriggerPattern[seg] = data;
      */
      break;
    case CID_GPIO:
      if (at == 0){
	Evidt[seg]  = 0xFFF&data;
	//	std::cout<<"Evidt: "<<Evidt<<std::endl;
      }
      if (at == 1){
	Spillt[seg] = 0xFFF&data;
	//	std::cout<<"Spillt: "<<Spillt<<std::endl;
      }
      break;
    default:
      break;
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
// void EventHeader::SetVMEinfo( EventStruct &vevent )
// {
//   Utime=vevent.utime;
//   Spillv=vevent.spill_v;
//   Evidv=vevent.evid_v;
//   SmpSw=vevent.smpSw;
//   TkoSw=vevent.tkoSw;
//   AdcCnt=vevent.adcCnt;
// }
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //

// int EventHeader::EventMatch(TTree *mtdcTree, int &spill_gpio, int &evid_gpio)
// {
//   mtdcTree->GetEntry(VEventNumber);
//   Evidv=evid_gpio;
//   Spillv=spill_gpio;
  
//   if(Spillv==Spillt && Evidv==Evidt){
//     VMESYNC=true;
//     return 1;
//   }else if(Evidv>Evidt && Evidv<Evidt+5){ //for MTDC lost events
//     VMESYNC=false;    
//     return 2;
//   }else if(Evidv>=Evidt+5){ //for MTDC lost events
//     VMESYNC=false;    
//     return 20;
//   }else if(Evidv<Evidt && Evidv+5>Evidt){ //for TKO lost events
//     VMESYNC=false;    
//     return 3;
//   }else if(Evidv+5<=Evidt){ //for TKO lost events
//     VMESYNC=false;    
//     return 3;
//   }  
// }
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
// bool EventHeader::EventMatch(TTree *vtree,EventStruct &veve, ConfMan *conf)
// {
//   //  for(int i=VEventNumber-1;i<(int)vtree->GetEntries();i++){

//   for(int i=VEventNumber;i<(int)vtree->GetEntries();i++){
//     vtree->GetEntry(i);
//     //  std::cout<<i<<std::endl;
//     //    std::cout<<veve.spill_v<<std::endl;
//     Evidv=veve.evid_v;
//     Spillv=veve.spill_v;
//     if(Spillv==Spillt){
//       if(Evidv==Evidt){
// #if 0
// 	std::cout<<"EventMatched.  "<<i<<std::endl;
// 	std::cout<<Form("Evidt=%d : Evidv=%d : Spillt=%d : Spillv=%d",Evidt,Evidv,Spillt,Spillv)<<std::endl;
// #endif
// 	VMESYNC=true;
// 	VEventNumber=i;
// 	return true;
//       }
//     }else if(!(Spillv<Spillt||Spillv-Spillt>220)){
//       VEventNumber -= 10;      
//       return false;
//     }
//   }
//   VEventNumber = -1;      
//   return false;
// }
// // + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
// bool EventHeader::EventMatch(TTree *vtree,EventStruct &veve, ConfMan *conf,std::ofstream &fout)
// {
//   //  for(int i=VEventNumber-1;i<(int)vtree->GetEntries();i++){
  
//   for(int i=VEventNumber;i<(int)vtree->GetEntries();i++){
//     vtree->GetEntry(i);
//     //  std::cout<<i<<std::endl;
//     //    std::cout<<veve.spill_v<<std::endl;
//     Evidv=veve.evid_v;
//     Spillv=veve.spill_v;
//     if(Spillv==Spillt&&Evidv==Evidt){
// #if 0
//       std::cout<<"EventMatched.  "<<i<<std::endl;
//       std::cout<<Form("Evidt=%d : Evidv=%d : Spillt=%d : Spillv=%d",Evidt,Evidv,Spillt,Spillv)<<std::endl;
// #endif
//       fout<<"1\t"<<Evidv<<"\t"<<Spillv<<"\t"<<std::endl;
//       VMESYNC=true;
//       VEventNumber=i;
//       return true;
//     }else if(!(Spillv<Spillt||Spillv-Spillt>220)){
// 	VEventNumber -= 1;      
// 	return false;
//     }
//     fout<<"0\t"<<Evidv<<"\t"<<Spillv<<"\t"<<std::endl;
//   }
//   VEventNumber -= 1;      
//   return false;
// }
void EventHeader::Initialize()
{
  VEventNumber = -1;
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void EventHeader::Clear()
{
  RunNumber = -1;
  EventNumber = -1;
  BlockEventNumber = -1;
  Utime = -1;
  Spillv = -1;
  Evidv = -1;
  for(int i=0;i<2;i++){
    Spillt[i] = -1;
    Evidt[i] = -1;
  }
  SmpSw = -1;
  TkoSw = -1;
  AdcCnt = -1;
  VMESYNC=false;
  for( int i=0; i<20; i++ ) TriggerPattern[i] = -1;
}
