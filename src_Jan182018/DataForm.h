// DataForm.h

#ifndef DataForm_h
#define DataForm_h 1

/* unidaq */
const int UniRecTypeNormal = 0;
const int UniRecTypeBegin  = 1;
const int UniRecTypePause  = 2;
const int UniRecTypeResume = 3;
const int UniRecTypeEnd    = 4;

#define CHAR0(x) (((x)&0x000000FF))
#define CHAR1(x) (((x)&0x0000FF00)>>8)
#define CHAR2(x) (((x)&0x00FF0000)>>16)
#define CHAR3(x) (((x)&0xFF000000)>>24)

/* Data Header               */
/* structure:                */
/* 0xAAAA0000 | NSMPA&0xFFFF */
const int Data_Header = 0xAAAA0000;
const int Data_NumSMPMask = 0x0000FFFF;

/* SMP Header                    */
/* structure:                    */
/* 1. SMP_HEADER                 */
/* 2. SMP_ADDRESS | NDATA&0xFFFF */
const int SMP_Header = 0xABCDABCD;
const int SMP_Footer = 0xDCBADCBA;
const int SMP_AddressMask = 0xFFFF0000;
const int SMP_NumDataMask = 0x0000FFFF;

/* TKO Sparce Data */
const int SD_ModuleDownBit = 0xFF000000;
#define SD_ModuleDownBitMask(x) ( (x)&0xFF800000 )
#define SD_GetNumEv(x)    ( (x)     &0xFFFF)
#define SD_GetSerialNo(x) (((x)>>16)&0xFFFF)
#define SD_GetData(x)     ( (x)     &0xFFFF)
#define SD_GetSA(x)       (((x)>>16)&0xFF  )
#define SD_GetMA(x)       (((x)>>27)&0x1F  )

#define DrT_GetData(x)    ( (x)     &0x7FF )
#define DrT_GetCH(x)      (((x)>>11)&0x1F  )

/* Modules */
const int HRTDC_RANGE_MIN = 0;
const int HRTDC_RANGE_MAX = 4095;


#endif
