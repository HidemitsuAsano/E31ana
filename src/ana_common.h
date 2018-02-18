//////////////////////////////////
//        ana_common.h 
//         2010/11/09  
//              H. Shi 
//////////////////////////////////
#include <stdint.h>

#ifndef _ANA_COMMON_H_
#define _ANA_COMMON_H_

#define  NSMPEVT           256
#define  NCH                8
#define  NGPIORM            1
#define  NSIS3301           1
#define  NV1785N            1
#define  NV1290             2
#define  MAXHITS            50
//number of events from one v1290 buffer
#define  NEvent             20 

#define  SIS_DATABIT_MASK    0x3fff
#define  PADC_DATABIT_MASK   0x0fff

// Structures for raw data format
typedef struct{
   unsigned int size;
   int tepoch;
   int len;    // header length -2
   int iev;    // 0xffffffff & iev: event number
   int iev32;  // 0xffffffff & (iev >> 32)
   int serial;   // 0x0123fedc
   int utepoch;  // originally blank1
   int blank2;
   int sis_sample;   // for user bit
   int blank3;
   int rpv100L;
   int gpioRML;
   int sisL;
   int v1290L;
   int blank4;
   int blank5;
} HeaderStruct;

typedef struct{
   int len1;
   int tag1;     // 0x48467935
   int num;      // number of rpv100 modules
   int blank1;
   int len2;     // 8 words for 8 channel counters
   int blank2;
   int i;        // rpv100 ID number
   int tag2;     // 0x39185326
   int count[8];
} Rpv100Struct;


typedef struct{
   int len1;
   int tag1;    //0x37641234
   int num;     // number of gpioRM modules
   int blank1;
   int len2;     // 4 words data
   int blank2;
   int i;        // gpioRM ID number
   int tag2;     // 0xabcdec12
   int evid;     // 12 bit
   int spill;    // 8 bit
   int serial;   // 31 bit (not recorded)
   int blank3;
} GpioStruct;

typedef struct{
   int len1;
   int tag1;         // 0x12901290
   int num;          // number of v1290 modules
   int blank1;
   int len2;         // ?
   int stored_event;   
   int i;            // v1290 ID number
   int tag2;         // 0x1290aaaa
   int data_buf[1024];
} V1290Struct;

// Structure for root file
typedef struct{
  int  utime;
  int  evid_daq;
  int  spill_gpio;
  int  evid_gpio;
  int  rpv100[8];     // rpv100 scaler
  int  v1290_stored_event[NV1290]; //use to count for event
  int  v1290_gheader[NV1290]; //use to count for event
  int  v1290_gtrailer[NV1290]; //use to count for event
  int  v1290_theader[NV1290]; //should be 4 times gheader
  int  v1290_ttrailer[NV1290]; //should be 4 times gtrailer
  int  mtdc[NV1290][32][MAXHITS]; //mtdc[nch][nhit], tdc value

} EventStruct;

#endif
