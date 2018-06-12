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
#include <TMath.h>
//ROOT default unit: MeV, ns, cm , degree
const double piMass = 0.13957;
const double pMass = 0.938272;
const double nMass = 0.939565;
const double dMass = 1.87561;
const double lMass = 1.115683;
const double s0Mass= 1.192642;
const double spMass= 1.18937;
const double smMass= 1.197449;
const double kpMass = 0.4936;
const double k0Mass = 0.497614;
const double ThreeHeMass = 2.80839;

const double m      = 100.0;
const double cm     = 1.0;
const double mm     = 0.1;
const double um     = 0.0001;
const double MeV     = 1./1000.;
const double GeV     = 1.0;
const double ns     = 1.0;
const double ps     = 1/.1000;
const double degree = 1.0;
const double deg    = 1.0;
const double TwoPi   = TMath::Pi()*2.0;
const double Deg2Rad = TMath::Pi()/180.;
//const double Rad2Deg = 180./TMath::Pi();
const double Const = 0.299792458; // =m/10^9
//const double Const=29.97;// [cm/ns]

//const int Event_Max_Size = 8192;
const int Event_Max_Size = 0xFFFF;

const int NumOfCDCLayers=15;
const int NumOfCDCWiresInLayer[15]={81,81,81,99,99,108,108,126,126,153,153,162,162,180,180};
enum enumArgv{ Argv_ProgramName=0,
	       Argv_ConfFileName=1,
	       Argv_RunNum=2,
	       Argv_OutFileName=3,
	       Argv_InFileName=4,
	       Argv_SummaryFileName=5,
	       Argv_TrackFileName=6,
	       Argv_MTDCFileName=7
};
enum enumSetupType{ Setup_E15=0, Setup_HEATES=1, Setup_PSI=2 };
// Data type
enum gBeamParticle { Beam_Kaon     = 0,
		     Beam_Pion     = 1,
		     Beam_Proton   = 2,
		     Beam_Deuteron = 3,
		     Beam_Other    = 4
};
enum gForwardParticle { F_Kaon     = 0,
			F_Pion     = 1,
			F_Proton   = 2,
			F_Deuteron = 3,
			F_Neutron  = 4,
			F_Gamma    = 5,
			F_Other    = 6
};
const double parMass[7]={kpMass,piMass,pMass,dMass,nMass,0.,0.};
const double particleMass[7]={kpMass,piMass,pMass,dMass,nMass,0.,0.};
const double cdsMass[9]={piMass,pMass,dMass,ThreeHeMass,ThreeHeMass,piMass,kpMass,0.0001,0.0001};
enum gCDSParticle { CDS_PiPlus     = 0,
		    CDS_Proton     = 1,
		    CDS_Deuteron   = 2,
		    CDS_Triton     = 3,
		    CDS_Helium3    = 4,
		    CDS_PiMinus    = 5,
		    CDS_Kaon       = 6,
		    CDS_Other      = 7,
		    CDS_DEFAULT    = 8
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
enum gCounterID { CID_CDC      = 0,
		  CID_CDH      = 1,
		  CID_BHD      = 2,
		  //		  CID_PA      = 3,
		  CID_T0       = 4,
		  CID_DEF      = 5,
		  CID_E0       = 6,
		  //		  CID_B1      = 6, //
		  CID_LC1      = 7,
		  CID_LC2      = 8,
		  CID_AC       = 9,
		  CID_WC       = 10,
		  CID_GC       = 11,
		  //		  CID_Range   = 12,
		  //		  CID_B2      = 13,
		  CID_TOFstop  = 14,
		  CID_CVC      = 14,
		  //		  CID_PDC1    = 15,
		  CID_BLC1a    = 15,
		  //	 	  CID_PDC2    = 16,
		  CID_BLC1b    = 16,
		  CID_BLC2a    = 17,
		  CID_BLC2b    = 18,
		  CID_SDD      = 19,
		  CID_TES      = 20,
		  CID_BLC1     = 21,
		  CID_BLC2     = 22,
		  CID_FDC1     = 23,
		  CID_FDC2     = 24,
		  CID_ZVC      = 30,
		  CID_KDV      = 31,
		  CID_NC       = 32,
		  CID_BVC      = 33,
		  CID_PC       = 35,
		  CID_Longbar  = 36,
		  CID_LB       = 36,
		  CID_WVC      = 37,
		  CID_BPC      = 40,
		  CID_BPD      = 41,
		  CID_IH       = 42,
		  CID_T0pre    = 51,
		  CID_T0post   = 52,
		  CID_BHDpost  = 56,
		  CID_HVC1     = 61,
		  CID_HVC2     = 62,
		  CID_BC1     = 71,
		  CID_BC2     = 72,
		  CID_BC3     = 73,
		  CID_BC4     = 74,
		  CID_BHDmul   = 81,
		  CID_T0mul    = 82,
		  CID_BVCmul   = 83,
		  CID_HVC1mul  = 84,
		  CID_HVC2mul  = 85,
		  CID_REFmul   = 86,
		  CID_BD       = 90,
		  //	 	  CID_BDC     = 90,
		  CID_TEMP1    = 91,
		  CID_TEMP2    = 92,
		  CID_GPIO     = 97,
		  CID_MISC     = 98,
		  CID_TRIG     = 98,
		  CID_TEMP     = 99,
		  CID_Hall     = 100,
		  CID_Floor    = 101,
		  CID_BeamDump = 110,
		  CID_SideDump = 111,
		  CID_NShield  = 112,
		  CID_SideCon  = 113,
		  CID_DoorCon  = 114,
		  CID_Doraemon = 120,
		  CID_USWK     = 121,
		  CID_CDCCFRP  = 130,
		  CID_CDCMylar = 131,
		  //		  CID_TarChm   = 140,
		  CID_TarSys    =140,
		  CID_RadS     = 141,
		  CID_TarCFRP  = 142,
		  CID_TarCap   = 143,
		  CID_TarRing  = 144,
		  CID_TarChm   = 145,
		  CID_TarCell  = 150,
		  CID_Target   = 151,
		  CID_CellBe   = 152,
		  CID_CellTube   = 152,
		  CID_CellAlBe = 153,
		  CID_CellFlange = 153,
		  CID_BShield  = 154,
		  CID_CellWindow = 154,
		  CID_BFrange  = 155,
		  CID_CellRing = 156,
		  CID_Fiducial = 160,
		  CID_DegC =  200, 
		  CID_Degrader1 =  200, 
		  CID_DegCu = 201,
		  CID_Degrader2 =  201, 
		  CID_Collimator = 202,
		  CID_Shield1 = 203,
		  CID_Shield2 = 204,
		  CID_Shield3 = 205,
		  CID_Shield4 = 206,
		  CID_TESSYS = 210,
		  CID_RS60K = 211,
		  CID_RS1K = 212,
		  CID_RS50mK = 213,
		  CID_TESWindow = 214,
		  CID_TESChamber = 215,
		  CID_RSwindow = 216,
		  CID_TESCol = 220,
		  CID_TESSi = 221,
		  CID_TESGrid = 222,
		  CID_TESHS = 223,  //heat sink ??
		  CID_TESAntico = 224,
		  CID_TESdummy = 225,
		  CID_BeamProf = 230
};
const int BHodoIDList[]={CID_BHD, CID_T0, CID_DEF, CID_CVC, CID_NC,
			 CID_BVC, CID_PC, CID_LB, CID_WVC, CID_BPD,
			 CID_T0pre,CID_T0post,CID_BHDpost,CID_HVC1,CID_HVC2,
			 CID_BD,CID_TEMP1,CID_TEMP2};
const int nBHodoIDList=sizeof(BHodoIDList)/sizeof(int);

const int ChereIDList[]={CID_LC1,CID_LC2,CID_AC,CID_WC,CID_GC};
const int nChereIDList=sizeof(ChereIDList)/sizeof(int);

const int BLDCIDList[]={CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_BPC,
			CID_FDC1,CID_BLC1,CID_BLC2};
const int nBLDCIDList=sizeof(BLDCIDList)/sizeof(int);

const int MTDCIDList[]={ CID_BHDmul,  CID_T0mul ,  CID_BVCmul,	  CID_HVC1mul,  CID_HVC2mul,
			 CID_REFmul };
const int nMTDCIDList=sizeof(MTDCIDList)/sizeof(int);

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
		       Trig_KCDH2f   = 9,
		       Trig_PiCDH2   = 10,
		       Trig_KCDH3    = 10,
		       Trig_Kf       = 11,
		       Trig_1stMix   = 12,
		       Trig_Charged  = 13,
		       Trig_Neutral  = 14,
		       Trig_Cosmic   = 15,
		       Trig_KvBVC    = 16,
		       Trig_Reject   = 16,
		       Trig_SIM      = 17
};

enum gTriggerMode { Mode_Beam     = 1, // <-- correct ??
		    Mode_Kf       = 2,
		    Mode_KCDH1f   = 3,
		    Mode_PiN      = 4,
		    Mode_PiC      = 5,
		    Mode_KCDH1N   = 6,
		    Mode_KCDH1C   = 7,
		    Mode_KCDH2    = 8,
		    Mode_KvBVCN   = 9,
		    Mode_KvBVCC   = 10,
		    Mode_PiCDH1N  = 11,
		    Mode_PiCDH1C  = 12,
		    Mode_KCDH3    = 13,
		    Mode_Cosmic   = 15,
		    Mode_Reject   = 16,
		    Mode_Unknown  = 17,
		    Mode_SIM      = 18
};

// Crate Number
enum gCrateNumber { Crate_BLC1    = 0,
		    Crate_BLC2    = 1,
		    Crate_BLHodo  = 2,
		    Crate_CDC1    = 3,
		    Crate_CDC2    = 4,
		    Crate_CDC3    = 5,
		    Crate_CDSHodo = 6,
		    Crate_NC      = 7,
		    Crate_NCPC    = 8,
		    Crate_FDCBPD  = 9
};

const std::string DefaultFileName = "None!!";

const int DEFAULTI=-1;
const double DEFAULTD=-999.;
const int ERRORI=-11;
const double ERRORD=-9999.;
const double DEFVECT[3]={-999.,-999.,-999.};
#include "TObject.h"

class GlobalVariables : public TObject
{
 public:
  GlobalVariables();
  ~GlobalVariables();

  ClassDef(GlobalVariables, 1);
};

#endif
