//===========================================================
//root macro to create geometry file for K1.8BR spectrometer
// y.ma@riken June 10, 2013
//units: degree, g/cm3, mm
//===========================================================
// #if defined(__CINT__) && !defined(__MAKECINT__)
// {
//    Info("k18br_geom.C",
//         "Has to be run in compiled mode, esp. if you want to pass parameters.");
//    gSystem->CompileMacro("k18br_geom.C");
//    k18br_geom();
// }
// #else
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include "Getline.h"
#include "TGeoMaterial.h"
#include "TGeoManager.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "TEveGeoShapeExtract.h"
#include "TEveGeoShape.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TEveSelection.h"
#include "TEvePointSet.h"
#include "TEveRGBAPalette.h"
#include "TColor.h"
#include "TRandom.h"
#include "TMath.h"
#include "TTree.h"
#include "TKey.h"
#include "TGLAnnotation.h"
#include "TGLViewer.h"
#include "TGLViewerBase.h"
#include "TEveStraightLineSet.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TEveViewer.h"
#include "TEveTrackPropagator.h"

//ROOT default unit: cm , degree
const Double_t m      = 100.0;
const Double_t cm     = 1.0;
const Double_t mm     = 0.1;
const Double_t degree = 1.0;
const Double_t deg    = 1.0;
const Double_t twopi   = TMath::Pi()*2.0;
const Double_t deg2rad = twopi/360.0;
const Double_t rad2deg = 360.0/twopi;
enum gCounterID { CID_CDC     = 0,
		  CID_CDH     = 1,
		  CID_BHD     = 2,
		  //		  CID_PA      = 3,
		  CID_T0      = 4,
		  CID_E0      = 5,
		  CID_DEF     = 5,
		  //		  CID_B1      = 6, //
		  CID_LC1     = 7,
		  CID_LC2     = 8,
		  CID_AC      = 9,
		  CID_WC      = 10,
		  CID_GC      = 11,
		  //		  CID_Range   = 12,
		  //		  CID_B2      = 13,
		  CID_TOFstop = 14,
		  CID_CVC     = 14,
		  //		  CID_PDC1    = 15,
		  CID_BLC1a   = 15,
		  //		  CID_PDC2    = 16,
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
		  CID_BHDmul  = 81,
		  CID_T0mul   = 82,
		  CID_BVCmul  = 82,
		  CID_HVC1mul = 83,
		  CID_HVC2mul = 84,
		  CID_REFmul  = 85,
		  CID_BD      = 90,
		  CID_BeamDump= 90,
		  CID_TEMP1   = 91,
		  CID_TEMP2   = 92,
		  CID_GPIO    = 97,
		  CID_MISC    = 98,
		  CID_TEMP    = 99
};
gCounterID CID;

//Double_t pos_global[CIDs]={x, y, z, dx, dy, dz};
Double_t pos_global[100][6];
//Double_t pos_seg[CIDs][id]={x, y, z, dx, dy, dz, l, w, t};
Double_t pos_seg[100][200][9];
TGeoMedium *Vacuum;
TGeoMedium *LHe3;
TGeoMedium *Beryllium;
TGeoMedium *Air;
TGeoMedium *CFRP;
TGeoMedium *Mylar;
TGeoMedium *Scintillator;
TGeoMedium *Aluminium;
TGeoMedium *Concrete;
TGeoMedium *Iron;
TGeoMedium *Tungsten;
TGeoVolume *hadron_hall;

void k18br_parser(const char* conf_file);
void constructMaterial();
void constructHall(TGeoManager * k18br_geom);
void constructNC(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructCVC(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructPC(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructUSWK(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructFDC(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructDORA(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructCDC(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructCDH(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructIH(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructBPC(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructBPD(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructDEF(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructTGT(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructT0(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);
void constructBLC2ab(TGeoManager * k18br_geom, TGeoVolume *hadron_hall);

void k18br_geom(const char* conf_file="conf/Run49c/analyzer.conf")
//void k18br_geom(const char* conf_file)
{
  //  gSystem->Load("libGeom");
  //  gSystem->Load("libGdml"); 
  //parse conf file and save geometry information to 
  //Double_t pos_global[CIDs][6];
  //Double_t pos_seg[CIDs][IDs][9];
  k18br_parser(conf_file);
  for(Int_t i=0; i<6; i++)
    cout<<"pos_global[16][i] = "<<pos_global[16][i]<<endl;
  for(Int_t i=0; i<8; i++)
    cout<<"pos_seg[16][2][i] = "<<pos_seg[16][2][i]<<endl;
  
  TGeoManager *k18br_geom = new TGeoManager("K1.8BR spectrometer", 
					    "geometry file for K1.8BR spectrometer");
  constructMaterial();
  constructHadronHall(k18br_geom);
  //optional based on configuration
  constructNC(k18br_geom, hadron_hall);
  constructCVC(k18br_geom, hadron_hall);
  //  constructPC(k18br_geom, hadron_hall);
  constructUSWK(k18br_geom, hadron_hall);
  constructFDC(k18br_geom, hadron_hall);
  constructDORA(k18br_geom, hadron_hall);
  constructCDC(k18br_geom, hadron_hall);
  constructCDH(k18br_geom, hadron_hall);
  constructIH(k18br_geom, hadron_hall);
  constructBPC(k18br_geom, hadron_hall);
  //  constructBPD(k18br_geom, hadron_hall);
  //  constructDEF(k18br_geom, hadron_hall);
  constructTGT(k18br_geom, hadron_hall);
  constructT0(k18br_geom, hadron_hall);
  constructBLC2ab(k18br_geom, hadron_hall);

  //--- close the geometry
  k18br_geom->CloseGeometry();
  k18br_geom->Export("k18br_geom.root");
  k18br_geom->SetVisLevel(4);
  k18br_geom->SetVisOption(0);
  k18br_geom->GetVolume("hadron_hall")->Draw("ogl");
  //  k18br_geom->GetVolume("ushiwaka_assembly")->Draw("ogl,same");
  k18br_geom->CheckOverlaps(0.1, "d");
  k18br_geom->PrintOverlaps();
  //  k18br_geom->GetVolume("hadron_hall")->Raytrace();
}

void k18br_parser(const char* conf_file)
{
  const Int_t MAX_CHARS_PER_LINE = 512;
  const Int_t MAX_TOKENS_PER_LINE = 20;
  char* const DELIMITER = " ";
  Int_t nglobal = 6;
  Int_t nseg = 9;
  ifstream fin_conf;
  char BLDCWireMap_param_name[50];
  char GeomBL_param_name[50];
  char GeomCDS_param_name[50];
  fin_conf.open(conf_file);
  if (!fin_conf.good()){ 
    cout<<"cant open file: "<<conf_file<<endl;
    getchar();
  }
  while (!fin_conf.eof()){
    char buf[MAX_CHARS_PER_LINE];
    fin_conf.getline(buf, MAX_CHARS_PER_LINE);
    int n = 0;
    const char* token[MAX_TOKENS_PER_LINE] = {};
    token[0] = strtok(buf, DELIMITER);
    if (token[0]) {
      for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
	token[n] = strtok(0, DELIMITER);
	if (!token[n]) break;
      }
    }
    for (int i = 0; i < n; i++){
      if( !strcmp(token[0], "BLDCWireMap:") ){
	strcpy(BLDCWireMap_param_name, token[1]);
      }else if( !strcmp(token[0], "GeomBL:") ){
	strcpy(GeomBL_param_name, token[1]);
      }else if( !strcmp(token[0], "GeomCDS:") ){
	strcpy(GeomCDS_param_name, token[1]);
      }
    }
  }
  cout<<"BLDCWireMap_param_name = "<<BLDCWireMap_param_name<<endl;
  cout<<"GeomBL_param_name = "<<GeomBL_param_name<<endl;
  cout<<"GeomCDS_param_name = "<<GeomCDS_param_name<<endl;
  //getchar();

  ifstream fin_BLDCWireMap_param;
  fin_BLDCWireMap_param.open(BLDCWireMap_param_name);
  if (!fin_BLDCWireMap_param.good()){ 
    cout<<"cant open file: "<<BLDCWireMap_param_name<<endl;
    getchar();
  }
  while (!fin_BLDCWireMap_param.eof()){
    char buf[MAX_CHARS_PER_LINE];
    fin_BLDCWireMap_param.getline(buf, MAX_CHARS_PER_LINE);
    int n = 0;
    const char* token[MAX_TOKENS_PER_LINE] = {};
    token[0] = strtok(buf, DELIMITER);
    if (token[0]) {
      for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
	token[n] = strtok(0, DELIMITER); 
	if (!token[n]) break; 
      }
    }
    for (int i = 0; i < n; i++){
      if( !strcmp(token[0], "GPOS:") ){
	//cout << "Token[" << i << "] = " << token[i] << endl;
  	if(i >= 2 && i < nglobal+2 )
  	  pos_global[atoi(token[1])][i-2]=atof(token[i]);
      }else if( strcmp(token[0], "#") ){
	//cout << "Token[" << i << "] = " << token[i] << endl;
  	if(i >= 2 && i < nseg + 2 )
  	  pos_seg[atoi(token[0])][atoi(token[1])][i-2]=atof(token[i]);
      }
    }
  }

  ifstream fin_GeomBL_param;
  fin_GeomBL_param.open(GeomBL_param_name);
  if (!fin_GeomBL_param.good()){ 
    cout<<"cant open file: "<<GeomBL_param_name<<endl;
    getchar();
  }
  while (!fin_GeomBL_param.eof()){
    char buf[MAX_CHARS_PER_LINE];
    fin_GeomBL_param.getline(buf, MAX_CHARS_PER_LINE);
    int n = 0;
    const char* token[MAX_TOKENS_PER_LINE] = {};
    token[0] = strtok(buf, DELIMITER);
    if (token[0]) {
      for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
	token[n] = strtok(0, DELIMITER); 
	if (!token[n]) break; 
      }
    }
    for (int i = 0; i < n; i++){
      if( !strcmp(token[0], "GPOS:") ){
	//  	cout << "Token[" << i << "] = " << token[i] << endl;
  	if(i >= 2 && i < nglobal+2 )
  	  pos_global[atoi(token[1])][i-2]=atof(token[i]);
      }else if( strcmp(token[0], "#") ){
	// 	cout << "Token[" << i << "] = " << token[i] << endl;
  	if(i >= 2 && i < nseg+2 )
  	  pos_seg[atoi(token[0])][atoi(token[1])][i-2]=atof(token[i]);
      }
    }
  }
  ifstream fin_GeomCDS_param;
  fin_GeomCDS_param.open(GeomCDS_param_name);
  if (!fin_GeomCDS_param.good()){ 
    cout<<"cant open file: "<<GeomBL_param_name<<endl;
    getchar();
  }
  while (!fin_GeomCDS_param.eof()){
    char buf[MAX_CHARS_PER_LINE];
    fin_GeomCDS_param.getline(buf, MAX_CHARS_PER_LINE);
    int n = 0;
    const char* token[MAX_TOKENS_PER_LINE] = {};
    token[0] = strtok(buf, DELIMITER);
    if (token[0]) {
      for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
	token[n] = strtok(0, DELIMITER); 
	if (!token[n]) break; 
      }
    }
    for (int i = 0; i < n; i++){
      if( !strcmp(token[0], "GPOS:") ){
	//	cout << "Token[" << i << "] = " << token[i] << endl;
  	if(i >= 2 && i < nglobal+2 )
  	  pos_global[atoi(token[1])][i-2]=atof(token[i]);
      }else if( strcmp(token[0], "#") ){
	//	cout << "Token[" << i << "] = " << token[i] << endl;
  	if(i >= 2 && i < nseg+2 )
  	  pos_seg[atoi(token[0])][atoi(token[1])][i-2]=atof(token[i]);
      }
    }
  }
}

void constructMaterial(){
  //====================================================
  //--- define materials
  //====================================================
  TGeoMaterial *matVacuum    = new TGeoMaterial("Vacuum",    1.008,  1,  1.e-25);
  TGeoMaterial *matHydrogen  = new TGeoMaterial("Hydrogen",  1.01,   1,  0.089/1000.);
  TGeoMaterial *matHelium    = new TGeoMaterial("Helium",    4.00,   2,  0.179/1000.);
  TGeoMaterial *matLHe3      = new TGeoMaterial("LHe3",      4.00,   2,  0.080);
  TGeoMaterial *matArgon     = new TGeoMaterial("Argon",     39.95,  18, 1.784/1000.);
  TGeoMaterial *matOxygen    = new TGeoMaterial("Oxygen",    15.99,  8,  1.429/1000.);
  TGeoMaterial *matNitrogen  = new TGeoMaterial("Nitrogen",  14.01,  7,  1.251/1000.);
  TGeoMaterial *matBeryllium = new TGeoMaterial("Beryllium", 9.01,   4,  1.85);
  TGeoMaterial *matCarbon    = new TGeoMaterial("Carbon",    12.01,  6,  2.267);
  TGeoMaterial *matAluminium = new TGeoMaterial("Aluminium", 26.98,  13, 2.70);
  TGeoMaterial *matSilicon   = new TGeoMaterial("Silicon",   28.09,  14, 2.32);
  TGeoMaterial *matChlorine  = new TGeoMaterial("Chlorine",  35.45,  17, 3.2/1000.);
  TGeoMaterial *matCalcium   = new TGeoMaterial("Calcium",   40.08,  20, 1.55);
  TGeoMaterial *matIron      = new TGeoMaterial("Iron",      55.85,  26, 7.87);
  TGeoMaterial *matTungsten  = new TGeoMaterial("Tungsten",  183.84, 74, 19.25);
  TGeoMixture  *matAir       = new TGeoMixture("Air", 3, 1.20/1000.);  
  matAir->AddElement(matNitrogen, 0.78);
  matAir->AddElement(matOxygen,   0.21);
  matAir->AddElement(matArgon,    0.01);
  //carbon fiber
  TGeoMixture  *matCFRP      = new TGeoMixture("CFRP", 1, 1.70);  
  matCFRP->AddElement(matCarbon,  1);
  //carbon dioxcide
  TGeoMixture  *matCO2      = new TGeoMixture("CO2", 2, 1.97/1000.);  
  matCO2->AddElement(matCarbon,  1.0/3.0);
  matCO2->AddElement(matOxygen,  2.0/3.0);
  //Methane
  TGeoMixture  *matCH4      = new TGeoMixture("CH4", 2, 0.656/1000.);  
  matCH4->AddElement(matCarbon,    1.0/5.0);
  matCH4->AddElement(matHydrogen,  4.0/5.0);
  //ethane
  TGeoMixture  *matC2H6      = new TGeoMixture("C2H6", 2, 1.356/1000.);  
  matC2H6->AddElement(matCarbon,    2.0/8.0);
  matC2H6->AddElement(matHydrogen,  6.0/8.0);
  //isobutane
  TGeoMixture  *matIsobutane = new TGeoMixture("Isobutane", 2, 2.51/1000.);  
  matIsobutane->AddElement(matCarbon,    4.0/14.0);
  matIsobutane->AddElement(matHydrogen,  10.0/14.0);
  //liquid scintillator
  TGeoMixture  *matNE213 = new TGeoMixture("NE213", 2, 0.935);  
  matNE213->AddElement(matCarbon,    8.0/16.0);
  matNE213->AddElement(matHydrogen,  8.0/16.0);
  //plastic scintillator
  TGeoMixture  *matScintillator = new TGeoMixture("Scintillator", 2, 1.032);  
  matScintillator->AddElement(matCarbon,    8.0/16.0);
  matScintillator->AddElement(matHydrogen,  8.0/16.0);
  //PET
  TGeoMixture  *matPET = new TGeoMixture("PET", 3, 1.35);  
  matPET->AddElement(matCarbon,    5.0/11.0);
  matPET->AddElement(matHydrogen,  4.0/11.0);
  matPET->AddElement(matOxygen,    2.0/11.0);
  //Mylar
  TGeoMixture  *matMylar = new TGeoMixture("Mylar", 3, 1.39);  
  matMylar->AddElement(matCarbon,    5.0/11.0);
  matMylar->AddElement(matHydrogen,  4.0/11.0);
  matMylar->AddElement(matOxygen,    2.0/11.0);
  //Kapton
  TGeoMixture  *matKapton = new TGeoMixture("Kapton", 3, 1.42);  
  matKapton->AddElement(matCarbon,    22.0/39.0);
  matKapton->AddElement(matHydrogen,  10.0/39.0);
  matKapton->AddElement(matOxygen,    5.0/39.0);
  matKapton->AddElement(matNitrogen,  2.0/39.0);
  //Polyethylene
  TGeoMixture  *matPolyethylene = new TGeoMixture("Polyethylene", 2, 0.92);  
  matPolyethylene->AddElement(matCarbon,    2.0/6.0);
  matPolyethylene->AddElement(matHydrogen,  4.0/6.0);
  //Concrete
  TGeoMixture *matConcrete = new TGeoMixture("Concrete", 6, 2.3);
  matConcrete->AddElement(matSilicon,     0.227915);
  matConcrete->AddElement(matOxygen,      0.605410);
  matConcrete->AddElement(matHydrogen,    0.099720);
  matConcrete->AddElement(matCalcium,     0.049860);
  matConcrete->AddElement(matAluminium,   0.014245);
  matConcrete->AddElement(matIron,        0.002850);
  //Aerogel
  TGeoMixture  *matAerogel = new TGeoMixture("Aerogel", 2, 0.2);  
  matAerogel->AddElement(matSilicon,   1.0/3.0);
  matAerogel->AddElement(matOxygen,    2.0/3.0);
  //Polyvinyl Chloride
  TGeoMixture  *matPolyvinyl = new TGeoMixture("Polyvinyl", 3, 1.4);  
  matPolyvinyl->AddElement(matCarbon,   2.0/6.0);
  matPolyvinyl->AddElement(matHydrogen, 3.0/6.0);
  matPolyvinyl->AddElement(matChlorine, 1.0/6.0);
  //====================================================
  //--- define tracking media
  //====================================================
  Vacuum       = new TGeoMedium("Vacuum",       0,  matVacuum);
  LHe3         = new TGeoMedium("LHe3",         1,  matLHe3);
  Beryllium    = new TGeoMedium("Beryllium",    2,  matLHe3);
  Air          = new TGeoMedium("Air",          3,  matAir);
  CFRP         = new TGeoMedium("CFRP",         4,  matCFRP);
  Mylar        = new TGeoMedium("Mylar",        5,  matMylar);
  Scintillator = new TGeoMedium("Scintillator", 6,  matScintillator);
  Aluminium    = new TGeoMedium("Aluminium",    7,  matAluminium);
  Concrete     = new TGeoMedium("Concrete",     8,  matConcrete);
  Iron         = new TGeoMedium("Iron",         9,  matIron);
  Tungsten     = new TGeoMedium("Tungsten",     10, matTungsten);
}

void constructHadronHall(TGeoManager *k18br_geom){
  Double_t beam_line_height  = 2.0*m;
  //====================================================
  //--- hadron hall infrastructure
  //====================================================
  //====================================================
  //--- define hadron hall
  //x,y,z = half of real size
  //====================================================
  Double_t hadron_hall_x = 20.0*m/2.0;
  Double_t hadron_hall_y = 5.0*m/2.0;
  Double_t hadron_hall_z = 50.0*m/2.0;
  Double_t hadron_hall_org[3] = {0., 0., 0.};
  hadron_hall = k18br_geom->MakeBox("hadron_hall", 
				    Air, 
				    hadron_hall_x,
				    hadron_hall_y,
				    hadron_hall_z);
  k18br_geom->SetTopVolume(hadron_hall);
  //====================================================
  //--- define floor
  //====================================================
  Double_t floor_x = 16.0*m/2.0;
  Double_t floor_y = 0.2*m/2.0;
  Double_t floor_z = 24.0*m/2.0;
  TGeoVolume *floor = k18br_geom->MakeBox("floor", 
					  Concrete, 
					  floor_x,
					  floor_y,
					  floor_z);
  Double_t floor_pos_x =  0.0*m;
  Double_t floor_pos_y = -beam_line_height - floor_y;
  Double_t floor_pos_z =  7.0*m;
  Double_t floor_angle =  0.0*degree; 
  TGeoRotation *floor_rot = new TGeoRotation();
  floor_rot->RotateY(floor_angle);
  TGeoCombiTrans *floor_trans = new TGeoCombiTrans(floor_pos_x,
						   floor_pos_y,
						   floor_pos_z,
						   floor_rot);
  floor->SetLineColor(kYellow);
  hadron_hall->AddNode(floor, 0, floor_trans);
  //====================================================
  //--- define operator
  //====================================================
  Double_t operator_assembly_x = 60*cm/2.0;
  Double_t operator_assembly_y = 180*cm/2.0;
  Double_t operator_assembly_z = 30.0*cm/2.0;
  TGeoVolume *operator_assembly = k18br_geom->MakeBox("operator_assembly", 
						      Air, 
						      operator_assembly_x,
						      operator_assembly_y,
						      operator_assembly_z);
  Double_t operator_assembly_pos_x =  3.0*m;
  Double_t operator_assembly_pos_y = -beam_line_height + operator_assembly_y;
  Double_t operator_assembly_pos_z =  7.0*m;
  Double_t operator_assembly_angle =  90.0*degree; 
  TGeoRotation *operator_assembly_rot = new TGeoRotation();
  operator_assembly_rot->RotateY(operator_assembly_angle);
  TGeoCombiTrans *operator_assembly_trans = new TGeoCombiTrans(operator_assembly_pos_x,
							       operator_assembly_pos_y,
							       operator_assembly_pos_z,
							       operator_assembly_rot);
  operator_assembly->SetLineColor(kBlue);
  hadron_hall->AddNode(operator_assembly, 0, operator_assembly_trans);
  Double_t operator_head_rmin     = 0.0*cm;
  Double_t operator_head_rmax     = 15.0*cm;
  Double_t operator_head_thetamin = 0.0*degree;
  Double_t operator_head_thetamax = 180.0*degree;
  Double_t operator_head_phimin   = 0.0*degree;
  Double_t operator_head_phimax   = 360.0*degree;
  TGeoVolume *operator_head = k18br_geom->MakeSphere("operator_head", 
						     Air,
						     operator_head_rmin,
						     operator_head_rmax,
						     operator_head_thetamin,
						     operator_head_thetamax,
						     operator_head_phimin,
						     operator_head_phimax);
  Double_t operator_head_pos_x =  0.0*m;
  Double_t operator_head_pos_y =  0.75*m;
  Double_t operator_head_pos_z =  0.0*m;
  Double_t operator_head_angle =  0.0*degree; 
  TGeoRotation *operator_head_rot = new TGeoRotation();
  operator_head_rot->RotateY(operator_head_angle);
  TGeoCombiTrans *operator_head_trans = new TGeoCombiTrans(operator_head_pos_x,
							   operator_head_pos_y,
							   operator_head_pos_z,
							   operator_head_rot);
  operator_head->SetLineColor(kBlue);
  operator_assembly->AddNode(operator_head, 0, operator_head_trans);
  Double_t operator_body_x = 40.0*cm/2.0;
  Double_t operator_body_y = 50.0*cm/2.0;
  Double_t operator_body_z = 10.0*cm/2.0;
  TGeoVolume *operator_body = k18br_geom->MakeBox("operator_body", 
						  Air,
						  operator_body_x,
						  operator_body_y,
						  operator_body_z);
  Double_t operator_body_pos_x =  0.0*m;
  Double_t operator_body_pos_y =  0.35*m;
  Double_t operator_body_pos_z =  0.0*m;
  Double_t operator_body_angle =  0.0*degree; 
  TGeoRotation *operator_body_rot = new TGeoRotation();
  operator_body_rot->RotateY(operator_body_angle);
  TGeoCombiTrans *operator_body_trans = new TGeoCombiTrans(operator_body_pos_x,
							   operator_body_pos_y,
							   operator_body_pos_z,
							   operator_body_rot);
  operator_body->SetLineColor(kBlue);
  operator_assembly->AddNode(operator_body, 0, operator_body_trans);
  Double_t operator_arm_rmin = 0.0*cm;
  Double_t operator_arm_rmax = 5.0*cm;
  Double_t operator_arm_z    = 70.0*cm/2.0;
  TGeoVolume *operator_arm = k18br_geom->MakeTube("operator_arm", 
						  Air,
						  operator_arm_rmin,
						  operator_arm_rmax,
						  operator_arm_z);
  for(Int_t i=0; i<2; i++){
    Double_t operator_arm_pos_x;
    if( i==0 )
      operator_arm_pos_x =  20.0*cm + operator_arm_rmax;
    if( i==1 )
      operator_arm_pos_x = -20.0*cm - operator_arm_rmax;
    Double_t operator_arm_pos_y =  0.25*m;
    Double_t operator_arm_pos_z =  0.0*m;
    Double_t operator_arm_angle =  90.0*degree; 
    TGeoRotation *operator_arm_rot = new TGeoRotation();
    operator_arm_rot->RotateX(operator_arm_angle);
    TGeoCombiTrans *operator_arm_trans = new TGeoCombiTrans(operator_arm_pos_x,
							    operator_arm_pos_y,
							    operator_arm_pos_z,
							    operator_arm_rot);
    operator_arm->SetLineColor(kBlue);
    operator_assembly->AddNode(operator_arm, i, operator_arm_trans);
  }
  Double_t operator_leg_rmin = 0.0*cm;
  Double_t operator_leg_rmax = 10.0*cm;
  Double_t operator_leg_z    = 100.0*cm/2.0;
  TGeoVolume *operator_leg = k18br_geom->MakeTube("operator_leg", 
						  Air,
						  operator_leg_rmin,
						  operator_leg_rmax,
						  operator_leg_z);
  for(Int_t i=0; i<2; i++){
    Double_t operator_leg_pos_x;
    if( i==0 )
      operator_leg_pos_x =  operator_leg_rmax;
    if( i==1 )
      operator_leg_pos_x = -operator_leg_rmax;
    Double_t operator_leg_pos_y =  -40*cm;
    Double_t operator_leg_pos_z =   0.0*m;
    Double_t operator_leg_angle =   90.0*degree; 
    TGeoRotation *operator_leg_rot = new TGeoRotation();
    operator_leg_rot->RotateX(operator_leg_angle);
    TGeoCombiTrans *operator_leg_trans = new TGeoCombiTrans(operator_leg_pos_x,
							    operator_leg_pos_y,
							    operator_leg_pos_z,
							    operator_leg_rot);
    operator_leg->SetLineColor(kBlue);
    operator_assembly->AddNode(operator_leg, i, operator_leg_trans);
  }
  //====================================================
  //--- define beam dump
  //====================================================
  Double_t beam_dump_x = 8.9492/2.*m;
  Double_t beam_dump_y = (beam_line_height*2.0)/2.0;
  Double_t beam_dump_z = 1.0/2.*m;
  TGeoVolume *beam_dump = k18br_geom->MakeBox("beam_dump", 
					      Iron, 
					      beam_dump_x,
					      beam_dump_y,
					      beam_dump_z);
  Double_t beam_dump_pos_x = 0.6*m;
  Double_t beam_dump_pos_y = 0.0*m;
  Double_t beam_dump_pos_z = 16.8*m;
  Double_t beam_dump_angle = 14*degree; 
  TGeoRotation *beam_dump_rot = new TGeoRotation();
  beam_dump_rot->RotateY(beam_dump_angle);
  TGeoCombiTrans *beam_dump_trans = new TGeoCombiTrans(beam_dump_pos_x,
						      beam_dump_pos_y,
						      beam_dump_pos_z,
						      beam_dump_rot);
  beam_dump->SetLineColor(kGreen);
  hadron_hall->AddNode(beam_dump, 0, beam_dump_trans);
  //====================================================
  //--- define side dump
  //====================================================
  Double_t side_dump_x = 1.0/2.*m;
  Double_t side_dump_y = 2.0/2.*m;
  Double_t side_dump_z = 4.0/2.*m;
  TGeoVolume *side_dump = k18br_geom->MakeBox("side_dump", 
					      Iron, 
					      side_dump_x,
					      side_dump_y,
					      side_dump_z);
  Double_t side_dump_pos_x = 3.9*m;
  Double_t side_dump_pos_y = 0.0*m;
  Double_t side_dump_pos_z = 13.4*m;
  Double_t side_dump_angle = 14*degree; 
  TGeoRotation *side_dump_rot = new TGeoRotation();
  side_dump_rot->RotateY(side_dump_angle);
  TGeoCombiTrans *side_dump_trans = new TGeoCombiTrans(side_dump_pos_x,
						      side_dump_pos_y,
						      side_dump_pos_z,
						      side_dump_rot);
  side_dump->SetLineColor(kGreen);
  hadron_hall->AddNode(side_dump, 0, side_dump_trans);
  //====================================================
  //--- define neutron shield
  //====================================================
  Double_t neutron_shield_x = 0.3/2.*m;
  Double_t neutron_shield_y = (beam_line_height*2.0)/2.0;
  Double_t neutron_shield_z = 4.0/2.*m;
  TGeoVolume *neutron_shield = k18br_geom->MakeBox("neutron_shield", 
						   Iron, 
						   neutron_shield_x,
						   neutron_shield_y,
						   neutron_shield_z);
  Double_t neutron_shield_pos_x = 2.4*m;
  Double_t neutron_shield_pos_y = 0.0*m;
  Double_t neutron_shield_pos_z = 13.774*m;
  Double_t neutron_shield_angle = 14*degree; 
  TGeoRotation *neutron_shield_rot = new TGeoRotation();
  neutron_shield_rot->RotateY(neutron_shield_angle);
  TGeoCombiTrans *neutron_shield_trans = new TGeoCombiTrans(neutron_shield_pos_x,
							    neutron_shield_pos_y,
							    neutron_shield_pos_z,
							    neutron_shield_rot);
  neutron_shield->SetLineColor(kGreen);
  hadron_hall->AddNode(neutron_shield, 0, neutron_shield_trans);
  //====================================================
  //--- define side concrete
  //====================================================
  Double_t side_concrete_x = 1.0/2.*m;
  Double_t side_concrete_y = (beam_line_height*2.0)/2.0;
  Double_t side_concrete_z = 9.0/2.*m;
  TGeoVolume *side_concrete = k18br_geom->MakeBox("side_concrete", 
						  Concrete, 
						  side_concrete_x,
						  side_concrete_y,
						  side_concrete_z);
  Double_t side_concrete_pos_x = -5.2*m;
  Double_t side_concrete_pos_y =  0.0*m;
  Double_t side_concrete_pos_z =  14.1*m;
  Double_t side_concrete_angle =  14*degree; 
  TGeoRotation *side_concrete_rot = new TGeoRotation();
  side_concrete_rot->RotateY(side_concrete_angle);
  TGeoCombiTrans *side_concrete_trans = new TGeoCombiTrans(side_concrete_pos_x,
							   side_concrete_pos_y,
							   side_concrete_pos_z,
							   side_concrete_rot);
  side_concrete->SetLineColor(kGreen);
  hadron_hall->AddNode(side_concrete, 0, side_concrete_trans);
  //====================================================
  //--- define door concrete
  //====================================================
  Double_t door_concrete_x = 3.0/2.*m;
  Double_t door_concrete_y = (beam_line_height*2.0)/2.0;
  Double_t door_concrete_z = 1.0/2.*m;
  TGeoVolume *door_concrete = k18br_geom->MakeBox("door_concrete", 
  						  Concrete, 
  						  door_concrete_x,
  						  door_concrete_y,
  						  door_concrete_z);
  Double_t door_concrete_pos_x = -4.2*m;
  Double_t door_concrete_pos_y =  0.0*m;
  Double_t door_concrete_pos_z =  9.8*m;
  Double_t door_concrete_angle =  14*degree; 
  TGeoRotation *door_concrete_rot = new TGeoRotation();
  door_concrete_rot->RotateY(door_concrete_angle);
  TGeoCombiTrans *door_concrete_trans = new TGeoCombiTrans(door_concrete_pos_x,
  							   door_concrete_pos_y,
  							   door_concrete_pos_z,
  							   door_concrete_rot);
  door_concrete->SetLineColor(kGreen);
  hadron_hall->AddNode(door_concrete, 0, door_concrete_trans);
}

void constructCVC(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define charge veto counter (CVC)
  //====================================================
  CID=CID_CVC; //14
  Double_t CVC_assembly_x = (abs(pos_seg[CID][1][0]-pos_seg[CID][34][0]) + pos_seg[CID][1][7])*cm/2.0;
  Double_t CVC_assembly_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t CVC_assembly_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *CVC_assembly = k18br_geom->MakeBox("CVC_assembly", 
  						 Vacuum, 
  						 CVC_assembly_x,
  						 CVC_assembly_y,
  						 CVC_assembly_z);
  Double_t CVC_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t CVC_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t CVC_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *CVC_assembly_rot   = new TGeoRotation();
  if( pos_global[CID][3]==0 && pos_global[CID][4]==0 )
    CVC_assembly_rot->RotateY(0);
  TGeoCombiTrans *CVC_assembly_trans = new TGeoCombiTrans(CVC_assembly_pos_x,
  							  CVC_assembly_pos_y,
  							  CVC_assembly_pos_z,
  							  CVC_assembly_rot);
  CVC_assembly->SetLineColor(kBlack);
  hadron_hall->AddNode(CVC_assembly, 0, CVC_assembly_trans);			
  Double_t CVC_segment_x = pos_seg[CID][1][7]*cm/2.0;
  Double_t CVC_segment_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t CVC_segment_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *CVC_segment = k18br_geom->MakeBox("CVC_segment", 
  						Scintillator, 
  						CVC_segment_x,
  						CVC_segment_y,
  						CVC_segment_z);
  for (Int_t i=1; i<=34; i++){
    Double_t CVC_segment_pos_x = pos_seg[CID][i][0]*cm;
    Double_t CVC_segment_pos_y = pos_seg[CID][i][1]*cm;
    Double_t CVC_segment_pos_z = pos_seg[CID][i][2]*cm;
    TGeoRotation *CVC_segment_rot   = new TGeoRotation();
    if( pos_seg[CID][i][3]==0 && pos_seg[CID][i][4]==0 )    
      CVC_segment_rot->RotateY(0);
    TGeoCombiTrans *CVC_segment_trans = new TGeoCombiTrans(CVC_segment_pos_x,
  							   CVC_segment_pos_y,
  							   CVC_segment_pos_z,
  							   CVC_segment_rot);
    CVC_segment->SetLineColor(kYellow);
    CVC_assembly->AddNode(CVC_segment, i, CVC_segment_trans);
  }
}

void constructNC(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define neutron counter (NC)
  //====================================================
  CID=CID_NC; //32
  Double_t NC_assembly_x = (abs(pos_seg[CID][1][0]-pos_seg[CID][16][0]) + pos_seg[CID][1][7])*cm/2.0;
  Double_t NC_assembly_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t NC_assembly_z = (abs(pos_seg[CID][1][2]-pos_seg[CID][112][2]) + pos_seg[CID][1][8])*cm/2.0;
  TGeoVolume *NC_assembly = k18br_geom->MakeBox("NC_assembly", 
  						Vacuum, 
  						NC_assembly_x,
  						NC_assembly_y,
  						NC_assembly_z);
  Double_t NC_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t NC_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t NC_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *NC_assembly_rot   = new TGeoRotation();
  if( pos_global[CID][3]==0 && pos_global[CID][4]==0 )
    NC_assembly_rot->RotateY(0);
  TGeoCombiTrans *NC_assembly_trans = new TGeoCombiTrans(NC_assembly_pos_x,
  							 NC_assembly_pos_y,
  							 NC_assembly_pos_z,
  							 NC_assembly_rot);
  NC_assembly->SetLineColor(kBlack);
  hadron_hall->AddNode(NC_assembly, 0, NC_assembly_trans);			
  Double_t NC_segment_x = pos_seg[CID][1][7]*cm/2.0;
  Double_t NC_segment_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t NC_segment_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *NC_segment = k18br_geom->MakeBox("NC_segment", 
  					       Scintillator, 
  					       NC_segment_x,
  					       NC_segment_y,
  					       NC_segment_z);
  for (Int_t i=1; i<=112; i++){
    Double_t NC_segment_pos_x = pos_seg[CID][i][0]*cm;
    Double_t NC_segment_pos_y = pos_seg[CID][i][1]*cm;
    Double_t NC_segment_pos_z = pos_seg[CID][i][2]*cm;
    TGeoRotation *NC_segment_rot   = new TGeoRotation();
    if( pos_seg[CID][i][3]==0 && pos_seg[CID][i][4]==0 )
      NC_assembly_rot->RotateY(0);
    TGeoCombiTrans *NC_segment_trans = new TGeoCombiTrans(NC_segment_pos_x,
  							  NC_segment_pos_y,
  							  NC_segment_pos_z,
  							  NC_segment_rot);
    NC_segment->SetLineColor(kBlue);
    NC_assembly->AddNode(NC_segment, i, NC_segment_trans);
  }
}

void constructPC(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define proton counter (PC)
  //====================================================
  CID=CID_PC; //35
  Double_t PC_assembly_x = (abs(pos_seg[CID][1][0]-pos_seg[CID][32][0]) + pos_seg[CID][1][7])*cm/2.0;
  Double_t PC_assembly_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t PC_assembly_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *PC_assembly = k18br_geom->MakeBox("PC_assembly", 
  						Vacuum, 
  						PC_assembly_x,
  						PC_assembly_y,
  						PC_assembly_z);
  Double_t PC_assembly_pos_x = -3.079*m;
  Double_t PC_assembly_pos_y =  0.0*m;
  Double_t PC_assembly_pos_z =  13.496*m;
  Double_t PC_assembly_angle = -15.843*degree;
  TGeoRotation *PC_assembly_rot   = new TGeoRotation();
  PC_assembly_rot->RotateY(PC_assembly_angle);
  TGeoCombiTrans *PC_assembly_trans = new TGeoCombiTrans(PC_assembly_pos_x,
  							 PC_assembly_pos_y,
  							 PC_assembly_pos_z,
  							 PC_assembly_rot);
  PC_assembly->SetLineColor(kBlack);
  hadron_hall->AddNode(PC_assembly, 0, PC_assembly_trans);	   
  Double_t PC_segment_x = 10.0/2.*cm;
  Double_t PC_segment_y = 1.5/2.*m;
  Double_t PC_segment_z = 3.0/2.*cm;
  TGeoVolume *PC_segment = k18br_geom->MakeBox("PC_segment", 
  					       Scintillator, 
  					       PC_segment_x,
  					       PC_segment_y,
  					       PC_segment_z);
  for (Int_t i=1; i<=27; i++){
    Double_t PC_segment_pos_x = (14-i)*PC_segment_x*2.0;
    Double_t PC_segment_pos_y = 0.0;
    Double_t PC_segment_pos_z = 0.0;
    Double_t PC_segment_angle = 0.0*degree;
    TGeoRotation *PC_segment_rot   = new TGeoRotation();
    PC_segment_rot->RotateY(PC_segment_angle);
    TGeoCombiTrans *PC_segment_trans = new TGeoCombiTrans(PC_segment_pos_x,
  							  PC_segment_pos_y,
  							  PC_segment_pos_z,
  							  PC_segment_rot);
    PC_segment->SetLineColor(kRed);
    PC_assembly->AddNode(PC_segment, i, PC_segment_trans);
  }
}

void constructUSWK(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define ushiwaka yoke
  //====================================================
  Double_t ushiwaka_assembly_x = 2.4/2.*m;
  Double_t ushiwaka_assembly_y = 1.4/2.*m;
  Double_t ushiwaka_assembly_z = 1.4/2.*m;
  TGeoVolume *ushiwaka_assembly = k18br_geom->MakeBox("ushiwaka_assembly", 
  						      Vacuum, 
  						      ushiwaka_assembly_x,
  						      ushiwaka_assembly_y,
  						      ushiwaka_assembly_z);
  Double_t ushiwaka_assembly_pos_x = 0.0*m;
  Double_t ushiwaka_assembly_pos_y = 0.0*m;
  Double_t ushiwaka_assembly_pos_z = 2.5*m;
  Double_t ushiwaka_assembly_angle = 0.0;
  TGeoRotation *ushiwaka_assembly_rot = new TGeoRotation();
  ushiwaka_assembly_rot->RotateY(ushiwaka_assembly_angle);
  TGeoCombiTrans *ushiwaka_assembly_trans = new TGeoCombiTrans(ushiwaka_assembly_pos_x,
  							       ushiwaka_assembly_pos_y,
  							       ushiwaka_assembly_pos_z,
  							       ushiwaka_assembly_rot);
  ushiwaka_assembly->SetLineColor(kYellow);
  hadron_hall->AddNode(ushiwaka_assembly, 0, ushiwaka_assembly_trans);
  Double_t ushiwaka_yoke_ud_x     = 2.4/2.*m;
  Double_t ushiwaka_yoke_ud_y     = 0.5/2.*m;
  Double_t ushiwaka_yoke_ud_z     = 1.4/2.*m;
  TGeoVolume *ushiwaka_yoke_ud = k18br_geom->MakeBox("ushiwaka_yoke_ud", 
  						     Iron, 
  						     ushiwaka_yoke_ud_x,
  						     ushiwaka_yoke_ud_y,
  						     ushiwaka_yoke_ud_z);
  for(Int_t i = 0; i<2; i++){
    Double_t ushiwaka_yoke_ud_pos_x = 0.0*m;
    Double_t ushiwaka_yoke_ud_pos_y;
    if(i == 0)
      ushiwaka_yoke_ud_pos_y = 0.45*m;
    if(i == 1)
      ushiwaka_yoke_ud_pos_y = -0.45*m;
    Double_t ushiwaka_yoke_ud_pos_z = 0.0*m;
    Double_t ushiwaka_yoke_ud_angle = 0.0;
    TGeoRotation *ushiwaka_yoke_ud_rot   = new TGeoRotation();
    ushiwaka_yoke_ud_rot->RotateY(ushiwaka_yoke_ud_angle);
    TGeoCombiTrans *ushiwaka_yoke_ud_trans = new TGeoCombiTrans(ushiwaka_yoke_ud_pos_x,
  							       ushiwaka_yoke_ud_pos_y,
  							       ushiwaka_yoke_ud_pos_z,
  							       ushiwaka_yoke_ud_rot);
    ushiwaka_yoke_ud->SetLineColor(kYellow);
    ushiwaka_assembly->AddNode(ushiwaka_yoke_ud, i, ushiwaka_yoke_ud_trans);
  }
  Double_t ushiwaka_yoke_lr_x     = 0.79/2.*m;
  Double_t ushiwaka_yoke_lr_y     = 0.4/2.*m;
  Double_t ushiwaka_yoke_lr_z     = 1.4/2.*m;
  TGeoVolume *ushiwaka_yoke_lr = k18br_geom->MakeBox("ushiwaka_yoke_lr", 
  						     Iron, 
  						     ushiwaka_yoke_lr_x,
  						     ushiwaka_yoke_lr_y,
  						     ushiwaka_yoke_lr_z);
  for(Int_t i = 0; i<2; i++){
    Double_t ushiwaka_yoke_lr_pos_x;
    if(i == 0)
      ushiwaka_yoke_lr_pos_x = 0.805*m;
    if(i == 1)
      ushiwaka_yoke_lr_pos_x = -0.805*m;
    Double_t ushiwaka_yoke_lr_pos_y = 0.0*m;
    Double_t ushiwaka_yoke_lr_pos_z = 0.0*m;
    Double_t ushiwaka_yoke_lr_angle = 0.0;
    TGeoRotation *ushiwaka_yoke_lr_rot   = new TGeoRotation();
    ushiwaka_yoke_lr_rot->RotateY(ushiwaka_yoke_lr_angle);
    TGeoCombiTrans *ushiwaka_yoke_lr_trans = new TGeoCombiTrans(ushiwaka_yoke_lr_pos_x,
  							       ushiwaka_yoke_lr_pos_y,
  							       ushiwaka_yoke_lr_pos_z,
  							       ushiwaka_yoke_lr_rot);
    ushiwaka_yoke_lr->SetLineColor(kYellow);
    ushiwaka_assembly->AddNode(ushiwaka_yoke_lr, i, ushiwaka_yoke_lr_trans);
  }
}

void constructFDC(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define FDC
  //====================================================
  CID=CID_FDC1; //23
  Double_t FDC_assembly_x = 60.0*cm/2.0;
  Double_t FDC_assembly_y = 28.0*cm/2.0;
  Double_t FDC_assembly_z = 12.0*cm/2.0;
  TGeoVolume *FDC_assembly = k18br_geom->MakeBox("FDC_assembly", 
						 Vacuum, 
						 FDC_assembly_x,
						 FDC_assembly_y,
						 FDC_assembly_z);
  Double_t FDC_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t FDC_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t FDC_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *FDC_assembly_rot   = new TGeoRotation();
  FDC_assembly_rot->RotateX( TMath::ATan(pos_global[CID][4])*rad2deg );
  FDC_assembly_rot->RotateY( TMath::ATan(pos_global[CID][5])*rad2deg );
  TGeoCombiTrans *FDC_assembly_trans = new TGeoCombiTrans(FDC_assembly_pos_x,
							  FDC_assembly_pos_y,
							  FDC_assembly_pos_z,
							  FDC_assembly_rot);
  FDC_assembly->SetLineColor(kWhite);
  hadron_hall->AddNode(FDC_assembly, 0, FDC_assembly_trans);
  Double_t FDC_sense_wire_inner_r = 0.0*cm;
  Double_t FDC_sense_wire_outer_r = 0.002*cm;
  //  Double_t FDC_sense_wire_outer_r = 1.0*mm;
  Double_t FDC_sense_wire_z = pos_seg[CID][1][5]*cm/2.0;
  TGeoVolume *FDC_sense_wire = k18br_geom->MakeTube("FDC_sense_wire",
						    Tungsten,
						    FDC_sense_wire_inner_r,
						    FDC_sense_wire_outer_r,
						    FDC_sense_wire_z);
  FDC_sense_wire->SetLineColor(kRed);
  Double_t FDC_field_wire_inner_r = 0.0*cm;
  Double_t FDC_field_wire_outer_r = 0.008*cm;
  //  Double_t FDC_field_wire_outer_r = 1.0*mm;
  Double_t FDC_field_wire_z = pos_seg[CID][1][5]*cm/2.0;
  TGeoVolume *FDC_field_wire = k18br_geom->MakeTube("FDC_field_wire",
						    Tungsten,
						    FDC_field_wire_inner_r,
						    FDC_field_wire_outer_r,
						    FDC_field_wire_z);
  FDC_field_wire->SetLineColor(kBlue);
  for(Int_t i=1; i<=6; i++){
    for(Int_t j=1; j<=pos_seg[CID][i][0]; j++){
      Double_t FDC_sense_wire_pos_x = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm;
      Double_t FDC_sense_wire_pos_y = 0*cm;
      Double_t FDC_sense_wire_pos_z = pos_seg[CID][i][1]*cm;
      TGeoRotation *FDC_sense_wire_rot = new TGeoRotation();
      FDC_sense_wire_rot->RotateX(90.0);
      FDC_sense_wire_rot->RotateZ(pos_seg[CID][i][7]*degree);
      TGeoCombiTrans *FDC_sense_wire_trans = new TGeoCombiTrans(FDC_sense_wire_pos_x,
								  FDC_sense_wire_pos_y,
								  FDC_sense_wire_pos_z,
								  FDC_sense_wire_rot);
      FDC_assembly->AddNode(FDC_sense_wire, i*100+j, FDC_sense_wire_trans);
    }
  }
}

void constructDORA(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define DORA magnet
  //====================================================
  Double_t DORA_yoke_inner_r = 590.0*mm;
  Double_t DORA_yoke_outer_r = 985.0*mm;
  Double_t DORA_yoke_z = (1170.0*mm)/2.0;
  TGeoVolume *DORA_yoke = k18br_geom->MakeTube("DORA_yoke", 
  					      Iron, 
  					      DORA_yoke_inner_r,
  					      DORA_yoke_outer_r,
  					      DORA_yoke_z);
  Double_t DORA_yoke_pos_x = 0.0*m;
  Double_t DORA_yoke_pos_y = 0.0*m;
  Double_t DORA_yoke_pos_z = 0.0*m;
  Double_t DORA_yoke_angle = 0.0*degree;
  TGeoRotation *DORA_yoke_rot = new TGeoRotation();
  DORA_yoke_rot->RotateY(DORA_yoke_angle);
  TGeoCombiTrans *DORA_yoke_trans = new TGeoCombiTrans(DORA_yoke_pos_x,
  						      DORA_yoke_pos_y,
  						      DORA_yoke_pos_z,
  						      DORA_yoke_rot);
  DORA_yoke->SetLineColor(kBlue);
  hadron_hall->AddNode(DORA_yoke, 0, DORA_yoke_trans);			  
  Double_t DORA_endcap_length = 155.0*mm;
  Double_t DORA_endcap_z = DORA_endcap_length/2.0;
  Double_t DORA_endcap_inner_r = 150.0*mm;
  Double_t DORA_endcap_outer_r = 985.0*mm;
  TGeoVolume *DORA_endcap = k18br_geom->MakeTube("DORA_endcap", 
						 Iron, 
						 DORA_endcap_inner_r,
						 DORA_endcap_outer_r,
						 DORA_endcap_z);
  for(Int_t i=0; i<2; i++){
    Double_t DORA_endcap_pos_x = 0.0*m;
    Double_t DORA_endcap_pos_y = 0.0*m;
    Double_t DORA_endcap_pos_z;
    if( i==0 )
      DORA_endcap_pos_z = -DORA_yoke_z - DORA_endcap_z;
    if( i==1 )
      DORA_endcap_pos_z =  DORA_yoke_z + DORA_endcap_z;
    Double_t DORA_endcap_angle = 0.0*degree;
    TGeoRotation *DORA_endcap_rot = new TGeoRotation();
    DORA_endcap_rot->RotateY(DORA_endcap_angle);
    TGeoCombiTrans *DORA_endcap_trans = new TGeoCombiTrans(DORA_endcap_pos_x,
							   DORA_endcap_pos_y,
							   DORA_endcap_pos_z,
							   DORA_endcap_rot);
    DORA_endcap->SetLineColor(kBlue);
    hadron_hall->AddNode(DORA_endcap, i, DORA_endcap_trans);			  
  }
}

void constructCDH(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  CID=CID_CDH; //1
  Double_t CDH_assembly_inner_r = 534.0*mm;
  Double_t CDH_assembly_outer_r = 584.0*mm;
  Double_t CDH_assembly_z       = pos_seg[CID][1][6]*cm/2.0;
  TGeoVolume *CDH_assembly = k18br_geom->MakeTube("CDH_assembly", 
  						  Vacuum, 
  						  CDH_assembly_inner_r,
  						  CDH_assembly_outer_r,
  						  CDH_assembly_z);
  Double_t CDH_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t CDH_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t CDH_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *CDH_assembly_rot = new TGeoRotation();
  if( pos_global[CID][3]==0 && pos_global[CID][4]==0 )  
    CDH_assembly_rot->RotateY(0);
  TGeoCombiTrans *CDH_assembly_trans = new TGeoCombiTrans(CDH_assembly_pos_x,
  							  CDH_assembly_pos_y,
  							  CDH_assembly_pos_z,
  							  CDH_assembly_rot);
  CDH_assembly->SetLineColor(kBlack);
  hadron_hall->AddNode(CDH_assembly, 0, CDH_assembly_trans);
  Double_t CDH_segment_x = pos_seg[CID][1][8]*cm/2.0;
  Double_t CDH_segment_y = pos_seg[CID][1][7]*cm/2.0;
  Double_t CDH_segment_z = pos_seg[CID][1][6]*cm/2.0;
  TGeoVolume *CDH_segment = k18br_geom->MakeBox("CDH_segment", 
  						Scintillator, 
  						CDH_segment_x,
  						CDH_segment_y,
  						CDH_segment_z);
  for (Int_t i=1; i<=36; i++){
    Double_t CDH_segment_pos_x = pos_seg[CID][i][0]*cm;
    Double_t CDH_segment_pos_y = pos_seg[CID][i][1]*cm;
    Double_t CDH_segment_pos_z = pos_seg[CID][i][2]*cm;
    TGeoRotation *CDH_segment_rot   = new TGeoRotation();
    if( i<=18 )
      CDH_segment_rot->RotateZ( TMath::ACos(pos_seg[CID][i][3])*rad2deg );
    else
      CDH_segment_rot->RotateZ( -TMath::ACos(pos_seg[CID][i][3])*rad2deg );
    TGeoCombiTrans *CDH_segment_trans = new TGeoCombiTrans(CDH_segment_pos_x,
							   CDH_segment_pos_y,
							   CDH_segment_pos_z,
							   CDH_segment_rot);
    CDH_segment->SetLineColor(kRed);
    CDH_assembly->AddNode(CDH_segment, i, CDH_segment_trans);
  }
}

Double_t  CDC_wire_pos_r[67];
void constructCDC(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define cylindrical drift chamber (CDC)
  //====================================================
  for (int i=0; i<67; i++){
    if(i==0){
      CDC_wire_pos_r[i] = 17.5*cm; 
    }else{
      CDC_wire_pos_r[i] = CDC_wire_pos_r[i-1] + 0.45*cm;
    }
    if(i==1 || i==12 || i==13 || i==21 || i==22 || i==30 || i==31 ||
       i==39 || i==40 || i==48 || i==49 || i==57 || i==58 || i==66){
      CDC_wire_pos_r[i] += 0.2*cm;
    }
  }
  Int_t get_CDC_numb_of_wires(Int_t i){
    if( i<12 ) {
      return 72;
    }
    else if( i>=12 && i<21 ) {
      return  90;
    }
    else if( i>=21 && i<30 ) { 
      return 100;
    }
    else if( i>=30 && i<39 ) { 
      return 120;
    }
    else if( i>=39 && i<48 ) { 
      return 150;
    }
    else if( i>=48 && i<57 ) { 
      return 160;
    }
    else if( i>=57 && i<=66 ) {
      return 180;
    }
  }
  Bool_t get_CDC_sense_wire_flag(Int_t i){
    if (i== 3 || i== 6 || i== 9 || i==15 || i==18 || i==24 || i==27 || i==33 || i==36 || i==42 || i==45 || i==51 || i==54 || i==60 || i==63 )
      return true;
    else 
      return false;
  }
  Bool_t get_idd_flag(Int_t i){
    if ( i== 0 || i== 1 || i== 3 || i== 5 || i== 7 || i== 9|| i==11 || i==12 || i==13 || i==15 || i==17 || i==19 || i==21 || i==22 || i==24 || i==26 || i==28|| i==30 || i==31 || i==33 || i==35 || i==37 || i==39 || i==40 || i==42 || i==44 || i==46 || i==48|| i==49 || i==51 || i==53 || i==55 || i==57 || i==58 || i==60 || i==62 || i==64 || i==66 ) 
      return false;
    else 
      return true;
  }
  Double_t get_CDC_wire_pos_r(Int_t i){
    return CDC_wire_pos_r[i];
  }
  Double_t get_CDC_wire_tilt(Int_t i, Double_t CDC_length){
    Double_t CDC_OFFSET = 3.0;
    if( i<12 || (i>=30 && i<39) || (i>=57 && i<=66))
      return 0.0;
    else 
      return (twopi/4-atan2(2*CDC_length, 2*get_CDC_wire_pos_r(i)*sin(twopi/get_CDC_numb_of_wires(i)*CDC_OFFSET/2.0)))*360.0/twopi;
  }
  Double_t get_CDC_wire_length(Int_t i, Double_t CDC_length){
    return (CDC_length/2.0)/cos(fabs(get_CDC_wire_tilt(i, CDC_length)));  
  }
  Double_t get_CDC_wire_pos_x(Int_t i, Int_t j){
    Double_t CDC_OFFSET = 3.0;
    Double_t unit_angle = twopi/get_CDC_numb_of_wires(i);
    Double_t Rwire      = get_CDC_wire_pos_r(i);
    Double_t x = 0.0;
    Double_t ang = 0.0;
    if( !get_idd_flag( i ) ) {
      ang = unit_angle*(Double_t)j;
      x   = Rwire*cos(ang);
    }else if( get_idd_flag( i ) ){
      ang = unit_angle*(Double_t)j + unit_angle/2.0;
      x   = Rwire*cos(ang);
    }
    return x;
  }
  Double_t get_CDC_wire_pos_y(Int_t i, Int_t j){
    Double_t CDC_OFFSET = 3.0;
    Double_t unit_angle = twopi/get_CDC_numb_of_wires(i);
    Double_t Rwire      = get_CDC_wire_pos_r(i);
    Double_t y = 0.0;
    Double_t ang = 0.0;
    if( !get_idd_flag( i ) ) {
      ang = unit_angle*(Double_t)j;
      y   = Rwire*sin(ang); 
    }else if( get_idd_flag( i ) ){
      ang = unit_angle*(Double_t)j + unit_angle/2.0;
      y   = Rwire*sin(ang); 
    }
    return y;
  }
  Double_t get_CDC_wire_pos_angle(Int_t i, Int_t j){
    Double_t CDC_OFFSET = 3.0;
    Double_t unit_angle = twopi/get_CDC_numb_of_wires(i);
    Double_t angle = 0.0;
    if( !get_idd_flag( i ) ) {
      angle = unit_angle*(Double_t)j;
    }else if( get_idd_flag( i ) ){
      angle = unit_angle*(Double_t)j + unit_angle/2.0;
    }
    return angle*360.0/twopi;
  }
  Double_t CDC_length           = 838.8*mm;
  Double_t CDC_assembly_inner_r = 150.0*mm;
  Double_t CDC_assembly_outer_r = 531.0*mm;
  Double_t CDC_assembly_z = CDC_length/2.0;
  TGeoVolume *CDC_assembly = k18br_geom->MakeTube("CDC_assembly", 
  						  Vacuum, 
  						  CDC_assembly_inner_r,
  						  CDC_assembly_outer_r,
  						  CDC_assembly_z);
  Double_t CDC_assembly_pos_x = 0.0*m;
  Double_t CDC_assembly_pos_y = 0.0*m;
  Double_t CDC_assembly_pos_z = 0.0*m;
  Double_t CDC_assembly_angle = 0.0*degree;
  TGeoRotation *CDC_assembly_rot = new TGeoRotation();
  CDC_assembly_rot->RotateY(CDC_assembly_angle);
  TGeoCombiTrans *CDC_assembly_trans = new TGeoCombiTrans(CDC_assembly_pos_x,
  							  CDC_assembly_pos_y,
  							  CDC_assembly_pos_z,
  							  CDC_assembly_rot);
  CDC_assembly->SetLineColor(kWhite);
  hadron_hall->AddNode(CDC_assembly, 0, CDC_assembly_trans);			  
  Double_t CDC_inner_window_inner_r = 150.0*mm;
  Double_t CDC_inner_window_outer_r = 151.0*mm;
  Double_t CDC_inner_window_z = CDC_length/2.0;
  TGeoVolume *CDC_inner_window = k18br_geom->MakeTube("CDC_inner_window", 
  						      CFRP, 
  						      CDC_inner_window_inner_r,
  						      CDC_inner_window_outer_r,
  						      CDC_inner_window_z);
  Double_t CDC_inner_window_pos_x = 0.0*m;
  Double_t CDC_inner_window_pos_y = 0.0*m;
  Double_t CDC_inner_window_pos_z = 0.0*m;
  Double_t CDC_inner_window_angle = 0.0*degree;
  TGeoRotation *CDC_inner_window_rot = new TGeoRotation();
  CDC_inner_window_rot->RotateY(CDC_inner_window_angle);
  TGeoCombiTrans *CDC_inner_window_trans = new TGeoCombiTrans(CDC_inner_window_pos_x,
  							      CDC_inner_window_pos_y,
  							      CDC_inner_window_pos_z,
  							      CDC_inner_window_rot);
  CDC_inner_window->SetLineColor(kWhite);
  CDC_assembly->AddNode(CDC_inner_window, 0, CDC_inner_window_trans);		
  Double_t CDC_outer_window_inner_r = 530.0*mm;
  Double_t CDC_outer_window_outer_r = 531.0*mm;
  Double_t CDC_outer_window_z = CDC_length/2.0;
  TGeoVolume *CDC_outer_window = k18br_geom->MakeTube("CDC_outer_window", 
  						      Mylar, 
  						      CDC_outer_window_inner_r,
  						      CDC_outer_window_outer_r,
  						      CDC_outer_window_z);
  Double_t CDC_outer_window_pos_x = 0.0*m;
  Double_t CDC_outer_window_pos_y = 0.0*m;
  Double_t CDC_outer_window_pos_z = 0.0*m;
  Double_t CDC_outer_window_angle = 0.0*degree;
  TGeoRotation *CDC_outer_window_rot = new TGeoRotation();
  CDC_outer_window_rot->RotateY(CDC_outer_window_angle);
  TGeoCombiTrans *CDC_outer_window_trans = new TGeoCombiTrans(CDC_outer_window_pos_x,
  							      CDC_outer_window_pos_y,
  							      CDC_outer_window_pos_z,
  							      CDC_outer_window_rot);
  CDC_outer_window->SetLineColor(kWhite);
  CDC_assembly->AddNode(CDC_outer_window, 0, CDC_outer_window_trans);		
  for(Int_t i=0; i<67; i++){
    for (Int_t j=0; j<get_CDC_numb_of_wires(i); j++){
      if( get_CDC_sense_wire_flag( i ) ){
  	Double_t CDC_sense_wire_inner_r = 0.0*mm;
  	Double_t CDC_sense_wire_outer_r = 0.03*mm;
  	//Double_t CDC_sense_wire_outer_r = 1.0*mm;
  	Double_t CDC_sense_wire_z = get_CDC_wire_length(i, CDC_length);
  	TGeoVolume *CDC_sense_wire = k18br_geom->MakeTube("CDC_sense_wire", 
  							  Tungsten,
  							  CDC_sense_wire_inner_r,
  							  CDC_sense_wire_outer_r,
  							  CDC_sense_wire_z);
  	Double_t CDC_sense_wire_pos_x = get_CDC_wire_pos_x(i, j);
  	Double_t CDC_sense_wire_pos_y = get_CDC_wire_pos_y(i, j);
  	Double_t CDC_sense_wire_pos_z = 0.0*mm;
  	Double_t CDC_sense_wire_angle = get_CDC_wire_tilt(i, CDC_length);
  	TGeoRotation *CDC_sense_wire_rot = new TGeoRotation();
  	CDC_sense_wire_rot->RotateX(-CDC_sense_wire_angle);
  	Double_t angle = get_CDC_wire_pos_angle(i, j);
  	CDC_sense_wire_rot->RotateZ(angle);
  	TGeoCombiTrans *CDC_sense_wire_trans = new TGeoCombiTrans(CDC_sense_wire_pos_x,
  								  CDC_sense_wire_pos_y,
  								  CDC_sense_wire_pos_z,
  								  CDC_sense_wire_rot);
  	CDC_sense_wire->SetLineColor(kRed);
  	CDC_assembly->AddNode(CDC_sense_wire, i*1000+j, CDC_sense_wire_trans);	   
      }else{
  	Double_t CDC_field_wire_inner_r = 0.0*mm;
  	Double_t CDC_field_wire_outer_r = 0.10*mm;
  	//Double_t CDC_field_wire_outer_r = 1.0*mm;
  	Double_t CDC_field_wire_z = get_CDC_wire_length(i, CDC_length);
  	TGeoVolume *CDC_field_wire = k18br_geom->MakeTube("CDC_field_wire", 
  							  Aluminium,
  							  CDC_field_wire_inner_r,
  							  CDC_field_wire_outer_r,
  							  CDC_field_wire_z);
  	Double_t CDC_field_wire_pos_x = get_CDC_wire_pos_x(i, j);
  	Double_t CDC_field_wire_pos_y = get_CDC_wire_pos_y(i, j);
  	Double_t CDC_field_wire_pos_z = 0.0*mm;
  	Double_t CDC_field_wire_angle = get_CDC_wire_tilt(i, CDC_length);
  	TGeoRotation *CDC_field_wire_rot = new TGeoRotation();
  	CDC_field_wire_rot->RotateX(-CDC_field_wire_angle);
  	Double_t angle = get_CDC_wire_pos_angle(i, j);
  	CDC_field_wire_rot->RotateZ(angle);
  	TGeoCombiTrans *CDC_field_wire_trans = new TGeoCombiTrans(CDC_field_wire_pos_x,
  								  CDC_field_wire_pos_y,
  								  CDC_field_wire_pos_z,
  								  CDC_field_wire_rot);
  	CDC_field_wire->SetLineColor(kGreen);
  	CDC_assembly->AddNode(CDC_field_wire, i*1000+j, CDC_field_wire_trans);	
      }
    }
  }
}

void constructIH(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  CID=CID_IH; //42
  //wait for proper param file
  Double_t IH_length         = 600.0*mm;
  Double_t IH_assembly_inner_r = 130.0*mm;
  Double_t IH_assembly_outer_r = 150.0*mm;
  Double_t IH_assembly_z       = IH_length/2.0;
  TGeoVolume *IH_assembly = k18br_geom->MakeTube("IH_assembly", 
  						  Vacuum, 
  						  IH_assembly_inner_r,
  						  IH_assembly_outer_r,
  						  IH_assembly_z);
  Double_t IH_assembly_pos_x = 0.0*m;
  Double_t IH_assembly_pos_y = 0.0*m;
  Double_t IH_assembly_pos_z = 0.0*m;
  Double_t IH_assembly_angle = 0.0*degree;
  TGeoRotation *IH_assembly_rot = new TGeoRotation();
  IH_assembly_rot->RotateY(IH_assembly_angle);
  TGeoCombiTrans *IH_assembly_trans = new TGeoCombiTrans(IH_assembly_pos_x,
  							 IH_assembly_pos_y,
  							 IH_assembly_pos_z,
  							 IH_assembly_rot);
  IH_assembly->SetLineColor(kBlack);
  hadron_hall->AddNode(IH_assembly, 0, IH_assembly_trans);			  
  Double_t IH_segment_x = 37.0*mm/2.0;
  Double_t IH_segment_y = 3.0*mm/2.0;
  Double_t IH_segment_z  = 600.0*mm/2.0;
  TGeoVolume *IH_segment = k18br_geom->MakeBox("IH_segment", 
  					       Scintillator, 
  					       IH_segment_x,
  					       IH_segment_y,
  					       IH_segment_z);
  for(Int_t i=0; i<24; i++){
    Double_t IH_segment_pos_x = 139.0*mm*cos(-5.0*deg2rad + i*15*deg2rad);
    Double_t IH_segment_pos_y = 139.0*mm*sin(-5.0*deg2rad + i*15*deg2rad);
    Double_t IH_segment_pos_z = 0.0*m;
    Double_t IH_segment_angle = i*15.0*deg + 90*deg;
    TGeoRotation *IH_segment_rot = new TGeoRotation();
    IH_segment_rot->RotateZ(IH_segment_angle);
    TGeoCombiTrans *IH_segment_trans = new TGeoCombiTrans(IH_segment_pos_x,
  							  IH_segment_pos_y,
  							  IH_segment_pos_z,
  							  IH_segment_rot);
    IH_segment->SetLineColor(kRed);
    IH_assembly->AddNode(IH_segment, i, IH_segment_trans);			  
  }
}

void constructBPC(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define BPC
  //====================================================
  CID=CID_BPC; //40
  // Double_t BPC_assembly_inner_r = 0.0*cm;
  // Double_t BPC_assembly_outer_r = 100.0*cm;
  // Double_t BPC_assembly_z = 5.0*cm/2.0;
  // TGeoVolume *BPC_assembly = k18br_geom->MakeTube("BPC_assembly", 
  // 						  Vacuum, 
  // 						  BPC_assembly_inner_r,
  // 						  BPC_assembly_outer_r,
  // 						  BPC_assembly_z);
  //temparory setting, waiting for proper param file
  Double_t BPC_assembly_x = 11.7*cm/2.0;
  Double_t BPC_assembly_y = 11.7*cm/2.0;
  Double_t BPC_assembly_z = 6.0*cm/2.0;
  TGeoVolume *BPC_assembly = k18br_geom->MakeBox("BPC_assembly", 
						 Vacuum, 
						 BPC_assembly_x,
						 BPC_assembly_y,
						 BPC_assembly_z);
  Double_t BPC_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t BPC_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t BPC_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *BPC_assembly_rot   = new TGeoRotation();
  BPC_assembly_rot->RotateX( TMath::ATan(pos_global[CID][4])*rad2deg );
  BPC_assembly_rot->RotateY( TMath::ATan(pos_global[CID][5])*rad2deg );
  TGeoCombiTrans *BPC_assembly_trans = new TGeoCombiTrans(BPC_assembly_pos_x,
							  BPC_assembly_pos_y,
							  BPC_assembly_pos_z,
							  BPC_assembly_rot);
  BPC_assembly->SetLineColor(kWhite);
  hadron_hall->AddNode(BPC_assembly, 0, BPC_assembly_trans);
  //temparory setting, waiting for proper param file
  Double_t BPC_mylar_x = 11.7*cm/2.0;
  Double_t BPC_mylar_y = 11.7*cm/2.0;
  Double_t BPC_mylar_z = 0.009*cm/2.0;
  TGeoVolume *BPC_mylar = k18br_geom->MakeBox("BPC_mylar", 
					      Mylar, 
					      BPC_mylar_x,
					      BPC_mylar_y,
					      BPC_mylar_z);
  BPC_mylar->SetLineColor(kWhite);
  Double_t BPC_sense_wire_inner_r = 0.0*cm;
  Double_t BPC_sense_wire_outer_r = 0.00125*cm;
  //  Double_t BPC_sense_wire_outer_r = 1.0*mm;
  Double_t BPC_sense_wire_z = pos_seg[CID][1][5]*cm/2.0;
  TGeoVolume *BPC_sense_wire = k18br_geom->MakeTube("BPC_sense_wire",
						    Tungsten,
						    BPC_sense_wire_inner_r,
						    BPC_sense_wire_outer_r,
						    BPC_sense_wire_z);
  BPC_sense_wire->SetLineColor(kRed);
  Double_t BPC_field_wire_inner_r = 0.0*cm;
  Double_t BPC_field_wire_outer_r = 0.005*cm;
  //  Double_t BPC_field_wire_outer_r = 1.0*mm;
  Double_t BPC_field_wire_z = pos_seg[CID][1][5]*cm/2.0;
  TGeoVolume *BPC_field_wire = k18br_geom->MakeTube("BPC_field_wire",
						    Tungsten,
						    BPC_field_wire_inner_r,
						    BPC_field_wire_outer_r,
						    BPC_field_wire_z);
  BPC_field_wire->SetLineColor(kBlue);
  for(Int_t i=1; i<=8; i++){
    for(Int_t j=1; j<=pos_seg[CID][i][0]; j++){
      Double_t BPC_sense_wire_pos_x;
      Double_t BPC_sense_wire_pos_y;
      Double_t BPC_sense_wire_pos_z;
      TGeoRotation *BPC_sense_wire_rot = new TGeoRotation();
      Double_t BPC_field_wirr_pos_x;
      Double_t BPC_field_wire_pos_y;
      Double_t BPC_field_wire_pos_z;
      TGeoRotation *BPC_field_wire_rot = new TGeoRotation();
      if( pos_seg[CID][i][2] == 0 ){
	BPC_sense_wire_pos_x = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm;
	BPC_sense_wire_pos_y = 0*cm;
	BPC_sense_wire_pos_z = pos_seg[CID][i][1]*cm;
	BPC_sense_wire_rot->RotateX(90.0);
	BPC_sense_wire_rot->RotateZ(0.0);
	BPC_field_wire_pos_x = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm + 0.36*cm;
	BPC_field_wire_pos_y = 0*cm;
	BPC_field_wire_pos_z = pos_seg[CID][i][1]*cm;
	BPC_field_wire_rot->RotateX(90.0);
	BPC_field_wire_rot->RotateZ(0.0);
      }else{
	BPC_sense_wire_pos_x = 0*cm;
	BPC_sense_wire_pos_y = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm;
	BPC_sense_wire_pos_z = pos_seg[CID][i][1]*cm;
	BPC_sense_wire_rot->RotateX(90.0);
	BPC_sense_wire_rot->RotateZ(90.0);
	BPC_field_wire_pos_x = 0*cm;
	BPC_field_wire_pos_y = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm + 0.36*cm;
	BPC_field_wire_pos_z = pos_seg[CID][i][1]*cm;
	BPC_field_wire_rot->RotateX(90.0);
	BPC_field_wire_rot->RotateZ(90.0);
      }
      TGeoCombiTrans *BPC_sense_wire_trans = new TGeoCombiTrans(BPC_sense_wire_pos_x,
								BPC_sense_wire_pos_y,
								BPC_sense_wire_pos_z,
								BPC_sense_wire_rot);
      BPC_assembly->AddNode(BPC_sense_wire, i*100+j, BPC_sense_wire_trans);
      TGeoCombiTrans *BPC_field_wire_trans = new TGeoCombiTrans(BPC_field_wire_pos_x,
								BPC_field_wire_pos_y,
								BPC_field_wire_pos_z,
								BPC_field_wire_rot);
      BPC_assembly->AddNode(BPC_field_wire, i*100+j, BPC_field_wire_trans);
    }
    Double_t BPC_mylar_pos_x = 0.0*cm;
    Double_t BPC_mylar_pos_y = 0.0*cm;
    Double_t BPC_mylar_pos_z = pos_seg[CID][i][1]*cm - 0.36*cm;
    TGeoRotation *BPC_mylar_rot = new TGeoRotation();
    BPC_mylar_rot->RotateX(0.0);
    TGeoCombiTrans *BPC_mylar_trans = new TGeoCombiTrans(BPC_mylar_pos_x,
    							 BPC_mylar_pos_y,
    							 BPC_mylar_pos_z,
    							 BPC_mylar_rot);
    BPC_assembly->AddNode(BPC_mylar, i, BPC_mylar_trans);
  }
}

void constructBPD(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define BPD
  //====================================================
  CID=CID_BPD; //41
  Double_t BPD_assembly_x = (abs(pos_seg[CID][1][0]-pos_seg[CID][32][0]) + pos_seg[CID][1][7])*cm/2.0;
  Double_t BPD_assembly_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t BPD_assembly_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *BPD_assembly = k18br_geom->MakeBox("BPD_assembly", 
  						 Vacuum, 
  						 BPD_assembly_x,
  						 BPD_assembly_y,
  						 BPD_assembly_z);
  Double_t BPD_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t BPD_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t BPD_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *BPD_assembly_rot   = new TGeoRotation();
  if( pos_global[CID][3]==0 && pos_global[CID][4]==0 )
    BPD_assembly_rot->RotateY(0);
  TGeoCombiTrans *BPD_assembly_trans = new TGeoCombiTrans(BPD_assembly_pos_x,
  							  BPD_assembly_pos_y,
  							  BPD_assembly_pos_z,
  							  BPD_assembly_rot);
  BPD_assembly->SetLineColor(kBlack);
  hadron_hall->AddNode(BPD_assembly, 0, BPD_assembly_trans);			
  Double_t BPD_segment_x = pos_seg[CID][1][7]*cm/2.0;
  Double_t BPD_segment_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t BPD_segment_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *BPD_segment = k18br_geom->MakeBox("BPD_segment", 
  						Scintillator, 
  						BPD_segment_x,
  						BPD_segment_y,
  						BPD_segment_z);
  for (Int_t i=1; i<=8; i++){
    Double_t BPD_segment_pos_x = pos_seg[CID][i][0]*cm;
    Double_t BPD_segment_pos_y = pos_seg[CID][i][1]*cm;
    Double_t BPD_segment_pos_z = pos_seg[CID][i][2]*cm;
    TGeoRotation *BPD_segment_rot   = new TGeoRotation();
    if( pos_seg[CID][i][3]==0 && pos_seg[CID][i][4]==0 )    
      BPD_segment_rot->RotateY(0);
    TGeoCombiTrans *BPD_segment_trans = new TGeoCombiTrans(BPD_segment_pos_x,
  							   BPD_segment_pos_y,
  							   BPD_segment_pos_z,
  							   BPD_segment_rot);
    BPD_segment->SetLineColor(kYellow);
    BPD_assembly->AddNode(BPD_segment, i, BPD_segment_trans);
  }
}

void constructDEF(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define DEF
  //====================================================
  CID=CID_DEF; //5
  Double_t DEF_assembly_x = (abs(pos_seg[CID][1][0]-pos_seg[CID][32][0]) + pos_seg[CID][1][7])*cm/2.0;
  Double_t DEF_assembly_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t DEF_assembly_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *DEF_assembly = k18br_geom->MakeBox("DEF_assembly", 
  						 Vacuum, 
  						 DEF_assembly_x,
  						 DEF_assembly_y,
  						 DEF_assembly_z);
  Double_t DEF_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t DEF_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t DEF_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *DEF_assembly_rot   = new TGeoRotation();
  if( pos_global[CID][3]==0 && pos_global[CID][4]==0 )
    DEF_assembly_rot->RotateY(0);
  TGeoCombiTrans *DEF_assembly_trans = new TGeoCombiTrans(DEF_assembly_pos_x,
  							  DEF_assembly_pos_y,
  							  DEF_assembly_pos_z,
  							  DEF_assembly_rot);
  DEF_assembly->SetLineColor(kBlack);
  hadron_hall->AddNode(DEF_assembly, 0, DEF_assembly_trans);			
  Double_t DEF_segment_x = pos_seg[CID][1][7]*cm/2.0;
  Double_t DEF_segment_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t DEF_segment_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *DEF_segment = k18br_geom->MakeBox("DEF_segment", 
  						Scintillator, 
  						DEF_segment_x,
  						DEF_segment_y,
  						DEF_segment_z);
  for (Int_t i=1; i<=8; i++){
    Double_t DEF_segment_pos_x = pos_seg[CID][i][0]*cm;
    Double_t DEF_segment_pos_y = pos_seg[CID][i][1]*cm;
    Double_t DEF_segment_pos_z = pos_seg[CID][i][2]*cm;
    TGeoRotation *DEF_segment_rot   = new TGeoRotation();
    if( pos_seg[CID][i][3]==0 && pos_seg[CID][i][4]==0 )    
      DEF_segment_rot->RotateY(0);
    TGeoCombiTrans *DEF_segment_trans = new TGeoCombiTrans(DEF_segment_pos_x,
  							   DEF_segment_pos_y,
  							   DEF_segment_pos_z,
  							   DEF_segment_rot);
    DEF_segment->SetLineColor(kYellow);
    DEF_assembly->AddNode(DEF_segment, i, DEF_segment_trans);
  }
}

void constructTGT(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define TGT
  //make new entry for TGT shift and rotation?
  //====================================================
  Double_t target_length     = 260.0*mm;
  Double_t target_assembly_inner_r = 0.0*mm;
  Double_t target_assembly_outer_r = 75.0*mm;
  Double_t target_assembly_z       = (target_length + 1*mm)/2.;
  TGeoVolume *target_assembly = k18br_geom->MakeTube("target_assembly", 
  						     Vacuum, 
  						     target_assembly_inner_r,
  						     target_assembly_outer_r,
  						     target_assembly_z);
  Double_t target_assembly_pos_x = 0.0*m;
  Double_t target_assembly_pos_y = 0.0*m;
  Double_t target_assembly_pos_z = 0.0*m;
  Double_t target_assembly_angle = 0.0*degree;
  TGeoRotation *target_assembly_rot = new TGeoRotation();
  target_assembly_rot->RotateY(target_assembly_angle);
  TGeoCombiTrans *target_assembly_trans = new TGeoCombiTrans(target_assembly_pos_x,
  							     target_assembly_pos_y,
  							     target_assembly_pos_z,
  							     target_assembly_rot);
  target_assembly->SetLineColor(kWhite);
  hadron_hall->AddNode(target_assembly, 0, target_assembly_trans);		    
  Double_t target_LHe3_inner_r = 0.0*mm;
  Double_t target_LHe3_outer_r = 34.0*mm;
  Double_t target_LHe3_z       = target_length/2.0;
  TGeoVolume *target_LHe3 = k18br_geom->MakeTube("target_LHe3", 
  						LHe3, 
  						target_LHe3_inner_r,
  						target_LHe3_outer_r,
  						target_LHe3_z);
  Double_t target_LHe3_pos_x = 0.0*mm;
  Double_t target_LHe3_pos_y = 0.0*mm;
  Double_t target_LHe3_pos_z = 0.0*mm;
  Double_t target_LHe3_angle = 0.0*deg;
  TGeoRotation *target_LHe3_rot = new TGeoRotation();
  target_LHe3_rot->RotateZ(target_LHe3_angle);
  TGeoCombiTrans *target_LHe3_trans = new TGeoCombiTrans(target_LHe3_pos_x,
  							target_LHe3_pos_y,
  							target_LHe3_pos_z,
  							target_LHe3_rot);
  target_LHe3->SetLineColor(kRed);
  target_assembly->AddNode(target_LHe3, 0, target_LHe3_trans);		  
  Double_t target_cell_inner_r = 34.0*mm;
  Double_t target_cell_outer_r = 34.3*mm;
  Double_t target_cell_z       = target_length/2.0;
  TGeoVolume *target_cell = k18br_geom->MakeTube("target_cell", 
  						 Beryllium, 
  						 target_cell_inner_r,
  						 target_cell_outer_r,
  						 target_cell_z);
  Double_t target_cell_pos_x = 0.0*mm;
  Double_t target_cell_pos_y = 0.0*mm;
  Double_t target_cell_pos_z = 0.0*mm;
  Double_t target_cell_angle = 0.0*deg;
  TGeoRotation *target_cell_rot = new TGeoRotation();
  target_cell_rot->RotateZ(target_cell_angle);
  TGeoCombiTrans *target_cell_trans = new TGeoCombiTrans(target_cell_pos_x,
  							 target_cell_pos_y,
  							 target_cell_pos_z,
  							 target_cell_rot);
  target_cell->SetLineColor(kBlack);
  target_assembly->AddNode(target_cell, 0, target_cell_trans);		  
  Double_t target_support_inner_r = 66.0*mm;
  Double_t target_support_outer_r = 66.3*mm;
  Double_t target_support_z       = 20.0*cm/2.0;
  TGeoVolume *target_support = k18br_geom->MakeTube("target_support", 
  						 Aluminium, 
  						 target_support_inner_r,
  						 target_support_outer_r,
  						 target_support_z);
  Double_t target_support_pos_x = 0.0*mm;
  Double_t target_support_pos_y = 0.0*mm;
  Double_t target_support_pos_z = 0.0*mm;
  Double_t target_support_angle = 0.0*deg;
  TGeoRotation *target_support_rot = new TGeoRotation();
  target_support_rot->RotateZ(target_support_angle);
  TGeoCombiTrans *target_support_trans = new TGeoCombiTrans(target_support_pos_x,
  							    target_support_pos_y,
  							    target_support_pos_z,
  							    target_support_rot);
  target_cell->SetLineColor(kBlack);
  target_assembly->AddNode(target_support, 0, target_support_trans);		  
  Double_t target_CFRP_inner_r = 74.0*mm;
  Double_t target_CFRP_outer_r = 75.0*mm;
  Double_t target_CFRP_z       = 26.0*cm/2.0;
  TGeoVolume *target_CFRP = k18br_geom->MakeTube("target_CFRP", 
  						 CFRP, 
  						 target_CFRP_inner_r,
  						 target_CFRP_outer_r,
  						 target_CFRP_z);
  Double_t target_CFRP_pos_x = 0.0*mm;
  Double_t target_CFRP_pos_y = 0.0*mm;
  Double_t target_CFRP_pos_z = 0.0*mm;
  Double_t target_CFRP_angle = 0.0*deg;
  TGeoRotation *target_CFRP_rot = new TGeoRotation();
  target_CFRP_rot->RotateZ(target_CFRP_angle);
  TGeoCombiTrans *target_CFRP_trans = new TGeoCombiTrans(target_CFRP_pos_x,
  							 target_CFRP_pos_y,
  							 target_CFRP_pos_z,
  							 target_CFRP_rot);
  target_cell->SetLineColor(kBlack);
  target_assembly->AddNode(target_CFRP, 0, target_CFRP_trans);		  
  Double_t target_endcap_inner_r = 0.0*mm;
  Double_t target_endcap_outer_r = 75.0*mm;
  Double_t target_endcap_z       = 0.5*mm/2.0;
  TGeoVolume *target_endcap = k18br_geom->MakeTube("target_endcap", 
  						   Aluminium, 
  						   target_endcap_inner_r,
  						   target_endcap_outer_r,
  						   target_endcap_z);
  for(Int_t i=0; i<2; i++){
    Double_t target_endcap_pos_x = 0.0*mm;
    Double_t target_endcap_pos_y = 0.0*mm;
    Double_t target_endcap_pos_z;
    if( i==0 )
      target_endcap_pos_z = -target_CFRP_z - target_endcap_z; 
    if( i==1 )
      target_endcap_pos_z =  target_CFRP_z + target_endcap_z; 
    Double_t target_endcap_angle = 0.0*deg;
    TGeoRotation *target_endcap_rot = new TGeoRotation();
    target_endcap_rot->RotateZ(target_endcap_angle);
    TGeoCombiTrans *target_endcap_trans = new TGeoCombiTrans(target_endcap_pos_x,
  							     target_endcap_pos_y,
  							     target_endcap_pos_z,
  							     target_endcap_rot);
    target_cell->SetLineColor(kBlack);
    target_assembly->AddNode(target_endcap, i, target_endcap_trans);
  }
}

void constructT0(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define T0
  //====================================================
  CID=CID_T0; //4
  Double_t T0_assembly_x = (fabs(pos_seg[CID][1][0]-pos_seg[CID][5][0]) + pos_seg[CID][1][7])*cm/2.0;
  Double_t T0_assembly_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t T0_assembly_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *T0_assembly = k18br_geom->MakeBox("T0_assembly", 
  						 Vacuum, 
  						 T0_assembly_x,
  						 T0_assembly_y,
  						 T0_assembly_z);
  Double_t T0_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t T0_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t T0_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *T0_assembly_rot   = new TGeoRotation();
  if( pos_global[CID][3]==0 && pos_global[CID][4]==0 )
    T0_assembly_rot->RotateY(0);
  TGeoCombiTrans *T0_assembly_trans = new TGeoCombiTrans(T0_assembly_pos_x,
  							 T0_assembly_pos_y,
  							 T0_assembly_pos_z,
  							 T0_assembly_rot);
  T0_assembly->SetLineColor(kBlack);
  hadron_hall->AddNode(T0_assembly, 0, T0_assembly_trans);			
  Double_t T0_segment_x = pos_seg[CID][1][7]*cm/2.0;
  Double_t T0_segment_y = pos_seg[CID][1][6]*cm/2.0;
  Double_t T0_segment_z = pos_seg[CID][1][8]*cm/2.0;
  TGeoVolume *T0_segment = k18br_geom->MakeBox("T0_segment", 
  						Scintillator, 
  						T0_segment_x,
  						T0_segment_y,
  						T0_segment_z);
  for (Int_t i=1; i<=5; i++){
    Double_t T0_segment_pos_x = pos_seg[CID][i][0]*cm;
    Double_t T0_segment_pos_y = pos_seg[CID][i][1]*cm;
    Double_t T0_segment_pos_z = pos_seg[CID][i][2]*cm;
    TGeoRotation *T0_segment_rot   = new TGeoRotation();
    if( pos_seg[CID][i][3]==0 && pos_seg[CID][i][4]==0 )    
      T0_segment_rot->RotateY(0);
    TGeoCombiTrans *T0_segment_trans = new TGeoCombiTrans(T0_segment_pos_x,
  							  T0_segment_pos_y,
  							  T0_segment_pos_z,
  							  T0_segment_rot);
    T0_segment->SetLineColor(kRed);
    T0_assembly->AddNode(T0_segment, i, T0_segment_trans);
  }
}

void constructBLC2ab(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
  //====================================================
  //--- define BLC2ab
  //====================================================
  CID=CID_BLC2a; //17
  Double_t BLC2ab_assembly_x = 34.0*cm/2.0;
  Double_t BLC2ab_assembly_y = 34.0*cm/2.0;
  Double_t BLC2ab_assembly_z = 32.0*cm/2.0;
  TGeoVolume *BLC2ab_assembly = k18br_geom->MakeBox("BLC2ab_assembly", 
						   Vacuum, 
						   BLC2ab_assembly_x,
						   BLC2ab_assembly_y,
						   BLC2ab_assembly_z);
  Double_t BLC2ab_assembly_pos_x = pos_global[CID][0]*cm;
  Double_t BLC2ab_assembly_pos_y = pos_global[CID][1]*cm;
  Double_t BLC2ab_assembly_pos_z = pos_global[CID][2]*cm;
  TGeoRotation *BLC2ab_assembly_rot   = new TGeoRotation();
  BLC2ab_assembly_rot->RotateX( TMath::ATan(pos_global[CID][4])*rad2deg );
  BLC2ab_assembly_rot->RotateY( TMath::ATan(pos_global[CID][5])*rad2deg );
  BLC2ab_assembly_rot->RotateZ( pos_seg[CID][1][7]*degree );
  TGeoCombiTrans *BLC2ab_assembly_trans = new TGeoCombiTrans(BLC2ab_assembly_pos_x,
							     BLC2ab_assembly_pos_y,
							    BLC2ab_assembly_pos_z,
							    BLC2ab_assembly_rot);
  BLC2ab_assembly->SetLineColor(kWhite);
  hadron_hall->AddNode(BLC2ab_assembly, 0, BLC2ab_assembly_trans);
  Double_t BLC2ab_mylar_x = 34.0*cm/2.0;
  Double_t BLC2ab_mylar_y = 34.0*cm/2.0;
  Double_t BLC2ab_mylar_z = 0.075*cm/2.0;
  TGeoVolume *BLC2ab_mylar = k18br_geom->MakeBox("BLC2ab_mylar", 
						 Mylar, 
						 BLC2ab_mylar_x,
						 BLC2ab_mylar_y,
						 BLC2ab_mylar_z);
  BLC2ab_mylar->SetLineColor(kWhite);

  //====================================================
  //--- define BLC2a
  //====================================================
  Double_t BLC2a_sense_wire_inner_r = 0.0*mm;
  Double_t BLC2a_sense_wire_outer_r = 0.012*mm;
  //  Double_t BLC2a_sense_wire_outer_r = 1.0*mm;
  Double_t BLC2a_sense_wire_z = pos_seg[CID][1][5]*cm;
  TGeoVolume *BLC2a_sense_wire = k18br_geom->MakeTube("BLC2a_sense_wire",
						      Tungsten,
						      BLC2a_sense_wire_inner_r,
						      BLC2a_sense_wire_outer_r,
						      BLC2a_sense_wire_z);
  BLC2a_sense_wire->SetLineColor(kRed);
  Double_t BLC2a_field_wire_inner_r = 0.0*mm;
  Double_t BLC2a_field_wire_outer_r = 0.075*mm;
  //  Double_t BLC2a_field_wire_outer_r = 1.0*mm;
  Double_t BLC2a_field_wire_z = pos_seg[CID][1][5]*cm;
  TGeoVolume *BLC2a_field_wire = k18br_geom->MakeTube("BLC2a_field_wire",
						      Tungsten,
						      BLC2a_field_wire_inner_r,
						      BLC2a_field_wire_outer_r,
						      BLC2a_field_wire_z);
  BLC2a_field_wire->SetLineColor(kBlue);
  for(Int_t i=1; i<=8; i++){
    for(Int_t j=1; j<=pos_seg[CID][i][0]; j++){
      Double_t BLC2a_sense_wire_pos_x;
      Double_t BLC2a_sense_wire_pos_y;
      Double_t BLC2a_sense_wire_pos_z;
      TGeoRotation *BLC2a_sense_wire_rot = new TGeoRotation();
      Double_t BLC2a_field_wirr_pos_x;
      Double_t BLC2a_field_wire_pos_y;
      Double_t BLC2a_field_wire_pos_z;
      TGeoRotation *BLC2a_field_wire_rot = new TGeoRotation();
      if( pos_seg[CID][i][2] == 0 ){
	BLC2a_sense_wire_pos_x = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm;
	BLC2a_sense_wire_pos_y = 0*cm;
	BLC2a_sense_wire_pos_z = pos_seg[CID][i][1]*cm;
	BLC2a_sense_wire_rot->RotateX(90.0);
	BLC2a_sense_wire_rot->RotateZ(0.0);
	BLC2a_field_wire_pos_x = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm + 0.5*cm;
	BLC2a_field_wire_pos_y = 0*cm;
	BLC2a_field_wire_pos_z = pos_seg[CID][i][1]*cm;
	BLC2a_field_wire_rot->RotateX(90.0);
	BLC2a_field_wire_rot->RotateZ(0.0);
      }else{
	BLC2a_sense_wire_pos_x = 0*cm;
	BLC2a_sense_wire_pos_y = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm;
	BLC2a_sense_wire_pos_z = pos_seg[CID][i][1]*cm;
	BLC2a_sense_wire_rot->RotateX(90.0);
	BLC2a_sense_wire_rot->RotateZ(90.0);
	BLC2a_field_wire_pos_x = 0*cm;
	BLC2a_field_wire_pos_y = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm + 0.5*cm;
	BLC2a_field_wire_pos_z = pos_seg[CID][i][1]*cm;
	BLC2a_field_wire_rot->RotateX(90.0);
	BLC2a_field_wire_rot->RotateZ(90.0);
      }
      TGeoCombiTrans *BLC2a_sense_wire_trans = new TGeoCombiTrans(BLC2a_sense_wire_pos_x,
								  BLC2a_sense_wire_pos_y,
								  BLC2a_sense_wire_pos_z,
								  BLC2a_sense_wire_rot);
      BLC2ab_assembly->AddNode(BLC2a_sense_wire, i*100+j, BLC2a_sense_wire_trans);
      TGeoCombiTrans *BLC2a_field_wire_trans = new TGeoCombiTrans(BLC2a_field_wire_pos_x,
								  BLC2a_field_wire_pos_y,
								  BLC2a_field_wire_pos_z,
								  BLC2a_field_wire_rot);
      BLC2ab_assembly->AddNode(BLC2a_field_wire, i*100+j, BLC2a_field_wire_trans);
    }
    Double_t BLC2a_mylar_pos_x = 0.0*cm;
    Double_t BLC2a_mylar_pos_y = 0.0*cm;
    Double_t BLC2a_mylar_pos_z = pos_seg[CID][i][1]*cm - 0.25*cm;
    TGeoRotation *BLC2a_mylar_rot = new TGeoRotation();
    BLC2a_mylar_rot->RotateX(0.0);
    TGeoCombiTrans *BLC2a_mylar_trans = new TGeoCombiTrans(BLC2a_mylar_pos_x,
							   BLC2a_mylar_pos_y,
							   BLC2a_mylar_pos_z,
							   BLC2a_mylar_rot);
    BLC2ab_assembly->AddNode(BLC2ab_mylar, i, BLC2a_mylar_trans);
  }

  //====================================================
  //--- define BLC2b
  //====================================================
  CID=CID_BLC2b; //17
  Double_t BLC2b_sense_wire_inner_r = 0.0*mm;
  Double_t BLC2b_sense_wire_outer_r = 0.012*mm;
  //  Double_t BLC2b_sense_wire_outer_r = 1.0*mm;
  Double_t BLC2b_sense_wire_z = pos_seg[CID][1][5];
  TGeoVolume *BLC2b_sense_wire = k18br_geom->MakeTube("BLC2b_sense_wire",
						      Tungsten,
						      BLC2b_sense_wire_inner_r,
						      BLC2b_sense_wire_outer_r,
						      BLC2b_sense_wire_z);
  BLC2b_sense_wire->SetLineColor(kRed);
  Double_t BLC2b_field_wire_inner_r = 0.0*mm;
  Double_t BLC2b_field_wire_outer_r = 0.075*mm;
  //  Double_t BLC2b_field_wire_outer_r = 1.0*mm;
  Double_t BLC2b_field_wire_z = pos_seg[CID][1][5];
  TGeoVolume *BLC2b_field_wire = k18br_geom->MakeTube("BLC2b_field_wire",
						      Tungsten,
						      BLC2b_field_wire_inner_r,
						      BLC2b_field_wire_outer_r,
						      BLC2b_field_wire_z);
  BLC2b_field_wire->SetLineColor(kBlue);
  for(Int_t i=1; i<=8; i++){
    for(Int_t j=1; j<=pos_seg[CID][i][0]; j++){
      Double_t BLC2b_sense_wire_pos_x;
      Double_t BLC2b_sense_wire_pos_y;
      Double_t BLC2b_sense_wire_pos_z;
      TGeoRotation *BLC2b_sense_wire_rot = new TGeoRotation();
      Double_t BLC2b_field_wirr_pos_x;
      Double_t BLC2b_field_wire_pos_y;
      Double_t BLC2b_field_wire_pos_z;
      TGeoRotation *BLC2b_field_wire_rot = new TGeoRotation();
      if( pos_seg[CID][i][2] == 0 ){
	BLC2b_sense_wire_pos_x = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm;
	BLC2b_sense_wire_pos_y = 0*cm;
	BLC2b_sense_wire_pos_z = pos_seg[CID][i][1]*cm;
	BLC2b_sense_wire_rot->RotateX(90.0);
	BLC2b_sense_wire_rot->RotateZ(0.0);
	BLC2b_field_wire_pos_x = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm + 0.5*cm;
	BLC2b_field_wire_pos_y = 0*cm;
	BLC2b_field_wire_pos_z = pos_seg[CID][i][1]*cm;
	BLC2b_field_wire_rot->RotateX(90.0);
	BLC2b_field_wire_rot->RotateZ(0.0);
      }else{
	BLC2b_sense_wire_pos_x = 0*cm;
	BLC2b_sense_wire_pos_y = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm;
	BLC2b_sense_wire_pos_z = pos_seg[CID][i][1]*cm;
	BLC2b_sense_wire_rot->RotateX(90.0);
	BLC2b_sense_wire_rot->RotateZ(90.0);
	BLC2b_field_wire_pos_x = 0*cm;
	BLC2b_field_wire_pos_y = (j-1)*pos_seg[CID][i][4]*cm + pos_seg[CID][i][3]*cm + 0.5*cm;
	BLC2b_field_wire_pos_z = pos_seg[CID][i][1]*cm;
	BLC2b_field_wire_rot->RotateX(90.0);
	BLC2b_field_wire_rot->RotateZ(90.0);
      }
      TGeoCombiTrans *BLC2b_sense_wire_trans = new TGeoCombiTrans(BLC2b_sense_wire_pos_x,
								  BLC2b_sense_wire_pos_y,
								  BLC2b_sense_wire_pos_z,
								  BLC2b_sense_wire_rot);
      BLC2ab_assembly->AddNode(BLC2b_sense_wire, i*100+j, BLC2b_sense_wire_trans);
      TGeoCombiTrans *BLC2b_field_wire_trans = new TGeoCombiTrans(BLC2b_field_wire_pos_x,
								  BLC2b_field_wire_pos_y,
								  BLC2b_field_wire_pos_z,
								  BLC2b_field_wire_rot);
      BLC2ab_assembly->AddNode(BLC2b_field_wire, i*100+j, BLC2b_field_wire_trans);
    }
    Double_t BLC2b_mylar_pos_x = 0.0*cm;
    Double_t BLC2b_mylar_pos_y = 0.0*cm;
    Double_t BLC2b_mylar_pos_z = pos_seg[CID][i][1]*cm - 0.25*cm;
    TGeoRotation *BLC2b_mylar_rot = new TGeoRotation();
    BLC2b_mylar_rot->RotateX(0.0);
    TGeoCombiTrans *BLC2b_mylar_trans = new TGeoCombiTrans(BLC2b_mylar_pos_x,
							   BLC2b_mylar_pos_y,
							   BLC2b_mylar_pos_z,
							   BLC2b_mylar_rot);
    BLC2ab_assembly->AddNode(BLC2ab_mylar, i, BLC2b_mylar_trans);
  }
}

//#endif

// void constructFDC(TGeoManager *k18br_geom, TGeoVolume *hadron_hall){
//   //====================================================
//   //--- define forward drift chamber (FDC)
//   //====================================================
//   Bool_t get_FDC_sense_wire_flag(Int_t i){
//     if( i==3 || i==6 || i==12 || i==15 || i==21 || i==24 )
//       return true;
//     else
//       return false;
//   }
//   Int_t get_FDC_numb_of_wires(Int_t i){
//     if( i==3 || i==6 || i==12 || i==15 || i==21 || i==24 || i==0 || i==9 || i==18 || i==27 )
//       return 64;
//     else
//       return 65;
//   }
//   Double_t get_FDC_wire_length(Int_t i){
//     if( i==12 || i==15 || i==0 || i==9 || i==18 || i==27 )
//       return 264.0*mm;
//     else
//       return 273.3*mm;
//   }
//   Double_t get_FDC_wire_tilt(Int_t i){
//     if( i==12 || i==15 || i==0 || i==9 || i==18 || i==27 )
//       return 0.0*degree;
//     else if ( i>=1 && i<=8 )
//       return 15.0*degree;
//     else
//       return -15.0*degree;
//   }
//   Double_t get_FDC_wire_pos_x(Int_t i, Int_t j){
//     if( i==0 || i==10 || i==12 || i==14 || i==16 || i==18 )
//       return -1.5*mm + (j-32)*6.0*mm;
//     if( i==9 || i==11 || i==13 || i==15 || i==17 || i==27 )
//       return +1.5*mm + (j-32)*6.0*mm;
//     if( i==1 || i==3 || i==5 || i==7 )
//       return +1.5*mm + (j-32)*6.0/cos(15.0*deg2rad)*mm;
//     if( i==2 || i==4 || i==6 || i==8 )
//       return -1.5*mm + (j-32)*6.0/cos(15.0*deg2rad)*mm;
//     if( i==19 || i==21 || i==23 || i==25 )
//       return +1.5*mm + (j-32)*6.0/cos(15.0*deg2rad)*mm;
//     if( i==20 || i==22 || i==24 || i==26 )
//       return -1.5*mm + (j-32)*6.0/cos(15.0*deg2rad)*mm;
//   }
//   Double_t get_FDC_wire_pos_y(Int_t i, Int_t j){
//     return 0.0*mm;
//   }
//   Double_t get_FDC_wire_pos_z(Int_t i, Int_t j){
//     if( i==0 )
//       return -30.0*mm;
//     else if( i==9 )
//       return -10.0*mm;
//     else if( i==18 )
//       return 10.0*mm;
//     else if( i==27 )
//       return 30.0*mm;
//     else if( 0<i&&i<9 ){
//       return -20.0*mm + ((i-4)-0.5)*sqrt(3.0)*mm;
//     }else if( 9<i&&i<18 ){
//       return 0.0*mm + ((i-4)-0.5)*sqrt(3.0)*mm;
//     }else if( 18<i&&i<27 ){
//       return 20.0*mm + ((i-4)-0.5)*sqrt(3.0)*mm;
//     }
//   }
//   Double_t FDC_assembly_x = 518.0*mm/2.0;
//   Double_t FDC_assembly_y = 280.0*mm/2.0;
//   Double_t FDC_assembly_z = 121.0*mm/2.0;
//   TGeoVolume *FDC_assembly = k18br_geom->MakeBox("FDC_assembly", 
//   						 Vacuum, 
//   						 FDC_assembly_x,
//   						 FDC_assembly_y,
//   						 FDC_assembly_z);
//   Double_t FDC_assembly_pos_x =  -18.9*cm;
//   Double_t FDC_assembly_pos_y =  0.0*m;
//   Double_t FDC_assembly_pos_z =  1696.0*mm;
//   Double_t FDC_assembly_angle =  0.0*degree;
//   TGeoRotation *FDC_assembly_rot = new TGeoRotation();
//   FDC_assembly_rot->RotateY(FDC_assembly_angle);
//   TGeoCombiTrans *FDC_assembly_trans = new TGeoCombiTrans(FDC_assembly_pos_x,
//   							  FDC_assembly_pos_y,
//   							  FDC_assembly_pos_z,
//   							  FDC_assembly_rot);
//   FDC_assembly->SetLineColor(kWhite);
//   hadron_hall->AddNode(FDC_assembly, 0, FDC_assembly_trans);			
//   for(Int_t i=0; i<28; i++){
//     for (Int_t j=0; j<get_FDC_numb_of_wires(i); j++){
//       if( get_FDC_sense_wire_flag( i ) ){
//   	Double_t FDC_sense_wire_inner_r = 0.0*mm;
//   	Double_t FDC_sense_wire_outer_r = 0.02*mm;
//   	//  	Double_t FDC_sense_wire_outer_r = 2.0*mm;
//   	Double_t FDC_sense_wire_z = get_FDC_wire_length(i)/2.0;
//   	TGeoVolume *FDC_sense_wire = k18br_geom->MakeTube("FDC_sense_wire", 
//   							  Tungsten,
//   							  FDC_sense_wire_inner_r,
//   							  FDC_sense_wire_outer_r,
//   							  FDC_sense_wire_z);
//   	Double_t FDC_sense_wire_pos_x = get_FDC_wire_pos_x(i, j);
//   	Double_t FDC_sense_wire_pos_y = get_FDC_wire_pos_y(i, j);
//   	Double_t FDC_sense_wire_pos_z = get_FDC_wire_pos_z(i, j);
//   	Double_t FDC_sense_wire_angle = get_FDC_wire_tilt(i);
//   	TGeoRotation *FDC_sense_wire_rot = new TGeoRotation();
//   	FDC_sense_wire_rot->RotateX(90.0*degree);
//   	FDC_sense_wire_rot->RotateZ(FDC_sense_wire_angle);
//   	TGeoCombiTrans *FDC_sense_wire_trans = new TGeoCombiTrans(FDC_sense_wire_pos_x,
//   								  FDC_sense_wire_pos_y,
//   								  FDC_sense_wire_pos_z,
//   								  FDC_sense_wire_rot);
//   	FDC_sense_wire->SetLineColor(kRed);
//   	FDC_assembly->AddNode(FDC_sense_wire, i*1000+j, FDC_sense_wire_trans);	   
//       }else{
//   	Double_t FDC_field_wire_inner_r = 0.0*mm;
//   	Double_t FDC_field_wire_outer_r = 0.08*mm;
//   	//  	Double_t FDC_field_wire_outer_r = 2.0*mm;
//   	Double_t FDC_field_wire_z = get_FDC_wire_length(i)/2.0;
//   	TGeoVolume *FDC_field_wire = k18br_geom->MakeTube("FDC_field_wire", 
//   							  Aluminium,
//   							  FDC_field_wire_inner_r,
//   							  FDC_field_wire_outer_r,
//   							  FDC_field_wire_z);
//   	Double_t FDC_field_wire_pos_x = get_FDC_wire_pos_x(i, j);
//   	Double_t FDC_field_wire_pos_y = get_FDC_wire_pos_y(i, j);
//   	Double_t FDC_field_wire_pos_z = get_FDC_wire_pos_z(i, j);
//   	Double_t FDC_field_wire_angle = get_FDC_wire_tilt(i);
//   	TGeoRotation *FDC_field_wire_rot = new TGeoRotation();
//   	FDC_field_wire_rot->RotateX(90.0*degree);
//   	FDC_field_wire_rot->RotateZ(FDC_field_wire_angle);
//   	TGeoCombiTrans *FDC_field_wire_trans = new TGeoCombiTrans(FDC_field_wire_pos_x,
//   								  FDC_field_wire_pos_y,
//   								  FDC_field_wire_pos_z,
//   								  FDC_field_wire_rot);
//   	FDC_field_wire->SetLineColor(kGreen);
//   	FDC_assembly->AddNode(FDC_field_wire, i*1000+j, FDC_field_wire_trans);	
//       }
//     }
//   }
// }

