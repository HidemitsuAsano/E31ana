// GlobalVariables.h
#ifndef Globals_h
#define Globals_h 1

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <complex>

#include <TH1.h>
#include <TH2.h>
#include <TObject.h>


const double piMass = 0.13957;
const double pMass = 0.938272;
const double nMass = 0.939565;
const double dMass = 1.87561;
const double lMass = 1.115683;
const double kpMass = 0.4936;
const double ThreeHeMass = 2.80839;


const double Const = 0.299792458; // =m/10^9
//const double Const=29.97;// [cm/ns]

//const int Event_Max_Size = 8192;
const int Event_Max_Size = 0xFFFF;

const int NumOfCDCLayers        = 15;
const int NumOfCDHSegments      = 36;
const int NumOfIHSegments       = 24;
const int NumOfBHDSegments      = 20;
const int NumOfPASegments       = 8;
const int NumOfT0Segments       = 5;
const int NumOfE0Segments       = 3;
const int NumOfDEFSegments       = 8;
const int NumOfB1Segments       = 1;
const int NumOfB2Segments       = 1;
const int NumOfLC1Segments      = 1;
const int NumOfLC2Segments      = 2;
const int NumOfACSegments       = 1;
const int NumOfWCSegments       = 1;
const int NumOfGCSegments       = 1;
const int NumOfRangeSegments    = 10;
//const int NumOfTOFstopSegments  = 32;
const int NumOfTOFstopSegments  = 34;
const int NumOfCVCSegments  = 34;
const int NumOfBVCSegments  = 8;
const int NumOfBDSegments        = 5;
const int NumOfLBSegments        = 2;
const int NumOfWVCSegments        = 2;
const int NumOfHVC1Segments        = 4;
const int NumOfHVC2Segments        = 4;
const int NumOfHVC3Segments        = 2;

const int NumOfNCLayer          = 7;
const int NumOfNCSegmentsInLayer= 16;
const int NumOfNCSegments       = 112;
const int NumOfPCSegments       = 27;
const int NumOfBPDSegments      = 70;
const int NumOfBLCLayers        = 8;
const int NumOfBPCLayers        = 8;
const int NumOfFDC1Layers        = 8;
const int NumOfSDDs             = 8;

const int NumOfBLCWiresInLayer=32;
const int NumOfBPCWiresInLayer=15;
const int NumOfFDC1WiresInLayer=64;

const int NumOfCDCWiresInLayer[15]={81,81,81,99,99,108,108,126,126,153,153,162,162,180,180};
// Data type
enum gBeamParticle { Beam_Kaon     = 0,
		     Beam_Pion     = 1,
		     Beam_Proton   = 2,
		     Beam_Other    = 3
};
const double Mass_Beam[4]={kpMass,piMass,pMass,0.};
const double particleMass[6]={piMass,pMass,dMass,piMass,kpMass,0.};
enum gCDSParticle { CDS_PiPlus     = 0,
		    CDS_Proton     = 1,
		    CDS_Deuteron   = 2,
		    CDS_PiMinus    = 3,
		    CDS_Kaon       = 4,
		    CDS_Other      = 5

};

enum gDataType { Type_CDS1  = 0,
		 Type_CDS2  = 1, // reserved
		 Type_CDS3  = 2, // reserved
		 Type_CDS4  = 3, // reserved
		 Type_CDS5  = 4, // reserved
		 Type_BL1   = 5,
		 Type_BL2   = 6, // reserved
		 Type_BL3   = 7, // reserved
		 Type_BL4   = 8, // reserved
		 Type_BL5   = 9, // reserved
		 Type_E15_1 = 10,
		 Type_E15_2 = 11, // reserved
		 Type_E15_3 = 12, // reserved
		 Type_E15_4 = 13, // reserved
		 Type_E15_5 = 14, // reserved
		 Type_SDD1 = 15, // nov beam
		 Type_SDD2 = 16, // sdd daq
		 Type_SDD3 = 17, // reserved
		 Type_SDD4 = 18,  // reserved
		 Type_SDD5 = 19  // reserved
};

// Counter ID
enum gCounterID { CID_CDC     = 0,
		  CID_CDH     = 1,
		  CID_BHD     = 2,
		  CID_PA      = 3,
		  CID_T0      = 4,
		  CID_E0      = 5,
		  CID_DEF     = 5,
		  CID_B1      = 6,
		  CID_LC1     = 7,
		  CID_LC2     = 8,
		  CID_AC      = 9,
		  CID_WC      = 10,
		  CID_GC      = 11,
		  CID_Range   = 12,
		  CID_B2      = 13,
		  CID_TOFstop = 14,
		  CID_CVC     = 14,
		  CID_PDC1    = 15,
		  CID_BLC1a   = 15,
		  CID_PDC2    = 16,
		  CID_BLC1b   = 16,
		  CID_BLC2a   = 17,
		  CID_BLC2b   = 18,
		  CID_SDD     = 19,
		  CID_BLC1    = 21,
		  CID_BLC2    = 22,
		  CID_FDC1    = 23,
		  CID_FDC2    = 24,
		  CID_ZVC     = 30,
		  CID_KDV     = 31,
		  CID_NC      = 32,
		  CID_BVC     = 33,
		  CID_PC      = 35,
		  CID_Longbar = 36,
		  CID_LB      = 36,
		  CID_WVC     = 37,
		  CID_BPC     = 40,
		  CID_BPD     = 41,
		  CID_IH      = 42,
		  CID_T0pre   = 51,
		  CID_T0post  = 52,
		  CID_BHDpost = 56,
		  CID_HVC1    = 61,
		  CID_HVC2    = 62,
		  CID_HVC3    = 63,
		  CID_BD      = 90,
		  CID_BeamDump= 90,
		  CID_VCC     = 91,
		  CID_TEMP1   = 91,
		  CID_TEMP2   = 92,
		  CID_TEMP3   = 93,
		  CID_GPIO    = 97,
		  CID_MISC    = 98,
		  CID_TEMP    = 99
};

// TriggerPattern
enum gTriggerPattern { Trig_Beam     = 1, // <-- correct ??
		       Trig_Kaon     = 2,
		       Trig_Electron = 3,
		       Trig_KCDH1f   = 3,
		       Trig_Pion     = 4,
		       Trig_Proton   = 5,
		       Trig_KCDH1    = 6,
		       Trig_KCDH2    = 7,
		       Trig_PivBVC   = 8,
		       Trig_PiCDH1   = 9,
		       Trig_PiCDH2   = 10,
		       Trig_Kf       = 11,
		       Trig_1stMix   = 12,
		       Trig_Charged  = 13,
		       Trig_Neutral  = 14,
		       Trig_Cosmic   = 15,
		       Trig_Reject   = 16
};

enum gTriggerMode { Mode_Beam     = 1, // <-- correct ??
		    Mode_Kf       = 2,
		    Mode_KCDH1f   = 3,
		    Mode_PiN      = 4,
		    Mode_PiC      = 5,
		    Mode_KCDH1N   = 6,
		    Mode_KCDH1C   = 7,
		    Mode_KCDH2    = 8,
		    Mode_Cosmic   = 9,
		    Mode_Reject   = 10,
		    Mode_PiCDH1   = 11,
		    Mode_Unknown  = 15
};

// Crate Number
enum gCrateNumber { Crate_PDC     = 0,
		    Crate_BLC     = 1,
		    Crate_BLHodo  = 2,
		    Crate_CDC1    = 3,
		    Crate_CDC2    = 4,
		    Crate_CDC3    = 5,
		    Crate_CDSHodo = 6,
		    Crate_NC      = 7,
		    Crate_NCPC    = 8,
		    Crate_BPD     = 9,
		    Crate_FDC1    = 10
};

const std::string DefaultFileName = "None!!";

#include "TObject.h"

class GlobalVariables : public TObject
{
 public:
  GlobalVariables();
  ~GlobalVariables();

  ClassDef(GlobalVariables, 1);
};

#endif
