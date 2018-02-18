//////////////////////////////////
//        vme_common.h 
//         2011/09/28  
//              T. Hashimoto 
//////////////////////////////////
#ifndef _VME_COMMON_H_
#define _VME_COMMON_H_ 1

#define NCH         8
#define MaxFADCPoint 32

typedef struct{
  int  utime;
  int  spill_v;
  int   evid_v;
  int  spill_t;
  int   evid_t;
  int  nsample[NCH];
  int  padcl[NCH];
  int  padch[NCH];
  int  fadc[NCH][MaxFADCPoint];
  int  smpSw;
  int  tkoSw;
  int  userB;      // User bit at sis                                          
  int  adcCnt;     // CAEN inner counter                                       
  int  tadc[NCH];
} EventStruct;

#endif
