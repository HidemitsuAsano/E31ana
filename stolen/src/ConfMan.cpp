// ConfMan.cpp

#include <string>
#include <iostream>
#include <cstdio>
#include <new>

#include "ConfMan.h"
#include "Tools.h"
#include "GeomTools.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <time.h>
#include <sys/time.h>

#ifndef NAN
#define NAN (0.0/0.0)
#endif

ClassImp(ConfMan)

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::ConfMan()
  : TObject()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::ConfMan( const std::string &filename )
  : TObject()
{
  Clear();
  ConfFileName = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::ConfMan( const std::string &filename,
                  int   runnum 
                )
  : TObject()
{
  Clear();
  ConfFileName = filename;
  RunNumber    = runnum;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::ConfMan( const std::string &filename,
                  int   runnum,
		  const std::string &outfilename) 
  : TObject()
{
  Clear();
  ConfFileName = filename;
  RunNumber    = runnum;
  OutFileName  = outfilename;
}
ConfMan::ConfMan( const std::string &filename,
                  int   runnum,
		  const std::string &outfilename, 
		  const std::string &vmefilename )
  : TObject()
{
  Clear();
  ConfFileName     = filename;
  RunNumber        = runnum;
  OutFileName      = outfilename;
  CDSTrackFileName = vmefilename;
}
ConfMan::ConfMan( const std::string &filename,
                  int   runnum,
		  const std::string &outfilename, 
		  const std::string &cdstrackfilename, 
		  const std::string &mtdcfilename )
  : TObject()
{
  Clear();
  ConfFileName     = filename;
  RunNumber        = runnum;
  OutFileName      = outfilename;
  CDSTrackFileName = cdstrackfilename;
  MTDCFileName     = mtdcfilename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::~ConfMan()
{
  if(CounterMapManager) delete CounterMapManager;
  if(BLDCWireMapManager) delete BLDCWireMapManager;
  if(CDCWireMapManager) delete CDCWireMapManager;
  if(GainMapManager) delete GainMapManager;
  if(GateMapManager) delete GateMapManager;
  if(XTMapManager) delete XTMapManager;
  if(GeomMapManager) delete GeomMapManager;
  if(SlewingMapManager) delete SlewingMapManager;
  if(ScalerMapManager) delete ScalerMapManager;
  if(ReslMapManager) delete ReslMapManager;
  if(CDSFittingParamManager) delete CDSFittingParamManager;
  if(BLDCFittingParamManager) delete BLDCFittingParamManager;
  if(FADCParamManager) delete FADCParamManager;
  if(DCTimeCorrManager) delete DCTimeCorrManager;
  if(TransferMatrixManager) delete TransferMatrixManager;
  //  End();
}

void ConfMan::Clear()
{
  ConfFileName = DefaultFileName;
  DataType = Type_CDS1;
  StartEvNum = StopEvNum = 0;
  StartBlockEvNum = StopBlockEvNum = 0;
  OutFileName = "tmp.root";
  MTDCFileName = DefaultFileName;
  CDSTrackFileName = DefaultFileName;
  CounterMapFileName = DefaultFileName;
  BLDCWireMapFileName = DefaultFileName;
  CDCGeometryFileName = DefaultFileName;
  CDCASDMapFileName = DefaultFileName;
  CDCChannelMapFileName = DefaultFileName;
  GainMapFileNameCDS = DefaultFileName;
  GainMapFileNameBL = DefaultFileName;
  GainMapFileNameSDD = DefaultFileName;
  GateMapFileNameCDS = DefaultFileName;
  GateMapFileNameBL = DefaultFileName;
  GateMapFileNameSDD = DefaultFileName;
  XTMapFileNameCDS = DefaultFileName;
  XTMapFileNameBL = DefaultFileName;
  XTMapROOTNameBL = DefaultFileName;
  GeometryMapFileNameCDS = DefaultFileName;
  GeometryMapFileNameBL = DefaultFileName;
  GeometryMapFileNameHall = DefaultFileName;
  SlewingMapFileNameCDC = DefaultFileName;
  SlewingMapFileNameCDS = DefaultFileName;
  SlewingMapFileNameBL = DefaultFileName;
  SlewingMapFileNameSDD = DefaultFileName;
  ReslMapFileName = DefaultFileName;
  CDSFittingParamFileName = DefaultFileName;
  BLDCFittingParamFileName = DefaultFileName;
  FADCParamFileName = DefaultFileName;
  TransferMatrixFileName = DefaultFileName;
  DCTimeCorrFileNameCDC = DefaultFileName;
  DCTimeCorrFileNameBLDC = DefaultFileName;
  DetectorListFileName = DefaultFileName;  
  MaterialListFileName = DefaultFileName;
  RunNumber         = (int)NULL;
  CurrentInSolenoid = -9999;
  CurrentInUshiwaka = -9999;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
static const int MAXCHAR = 144;
bool ConfMan::ReadConfFile(std::string conffilename)
{
  FILE *fp;
  char str[MAXCHAR], str1[MAXCHAR];
  int evnum;
  int typ;

  DB_Manager* database = new DB_Manager();
  int run_begin =   1;
  int run_end   =1000;

  bool status_CounterMap_List     = false; 
  bool status_GainMapBL_List      = false; 
  bool status_GateMapBL_List      = false; 
  bool status_GainMapCDS_List     = false; 
  bool status_DCTimeCorrBLDC_List = false; 
  bool status_DetectorList_List = false; 
  bool status_XTMapBLROOT_List    = false; 
  bool status_XTMapCDS_List       = false; 
  bool status_XTMapBL_List       = false; 
  bool status_SlewingMapBL_List   = false; 
  bool status_SlewingMapCDS_List  = false; 
  bool status_ReslMap_List        = false; 
  bool status_TransferMatrix_List        = false; 
  bool status_CDSFittingParam_List        = false; 
  bool status_CDCASDMap_List     = false; 

  if( (fp=fopen(conffilename.c_str(), "r"))==0 ){
    std::cerr << " file open fail : " << conffilename << std::endl;
    exit(-1);
  }
  std::cout << " file opend : " << conffilename << std::endl;

  while( fgets( str, MAXCHAR, fp )!=0 ){
    if( str[0]=='#' ) continue;
    
    if( sscanf(str,"include: %s", str1)==1 )
      ReadConfFile(str1);
    if( sscanf(str,"CounterMap: %s", str1)==1 )
      CounterMapFileName = str1;
    else if( sscanf(str,"BLDCWireMap: %s", str1)==1 )
      BLDCWireMapFileName = str1;
    else if( sscanf(str,"CDCGeom: %s", str1)==1 )
      CDCGeometryFileName = str1;
    else if( sscanf(str,"CDCASDMap: %s", str1)==1 )
      CDCASDMapFileName = str1;
    else if( sscanf(str,"CDCChannelMap: %s", str1)==1 )
      CDCChannelMapFileName = str1;
    else if( sscanf(str,"GainMapCDS: %s", str1)==1 )
      GainMapFileNameCDS = str1;
    else if( sscanf(str,"GainMapBL: %s", str1)==1 )
      GainMapFileNameBL  = str1;
    else if( sscanf(str,"GainMapBLprefix: %s", str1)==1 ){
      if(getenv("RUN")!=NULL)
	GainMapFileNameBL = Form("%s%s.param",str1,getenv("RUN"));
    }
    else if( sscanf(str,"GainMapSDD: %s", str1)==1 )
      GainMapFileNameSDD  = str1;
    else if( sscanf(str,"GateMapCDS: %s", str1)==1 )
      GateMapFileNameCDS = str1;
    else if( sscanf(str,"GateMapBL: %s", str1)==1 )
      GateMapFileNameBL  = str1;
    else if( sscanf(str,"GateMapBLprefix: %s", str1)==1 ){
      if(getenv("RUN")!=NULL)
	GateMapFileNameBL = Form("%s%s.param",str1,getenv("RUN"));
    }
    else if( sscanf(str,"GateMapSDD: %s", str1)==1 )
      GateMapFileNameSDD  = str1;
    else if( sscanf(str,"XTMapCDS: %s", str1)==1 )
      XTMapFileNameCDS = str1;
    else if( sscanf(str,"XTMapBL: %s", str1)==1 )
      XTMapFileNameBL = str1;
    else if( sscanf(str,"XTMapBLROOT: %s", str1)==1 )
      XTMapROOTNameBL = str1;
    else if( sscanf(str,"XTMapBLROOTprefix: %s", str1)==1 ){
      if(getenv("RUN")!=NULL)
	XTMapROOTNameBL = Form("%s%s.root",str1,getenv("RUN"));
    }
    else if( sscanf(str,"GeomCDS: %s", str1)==1 )
      GeometryMapFileNameCDS = str1;
    else if( sscanf(str,"GeomBL: %s", str1)==1 )
      GeometryMapFileNameBL = str1;
    else if( sscanf(str,"GeomHall: %s", str1)==1 )
      GeometryMapFileNameHall = str1;
    else if( sscanf(str,"SlewingMapCDC: %s", str1)==1 )
      SlewingMapFileNameCDC = str1;
    else if( sscanf(str,"SlewingMapCDS: %s", str1)==1 )
      SlewingMapFileNameCDS = str1;
    else if( sscanf(str,"SlewingMapBL: %s", str1)==1 )
      SlewingMapFileNameBL = str1;
    else if( sscanf(str,"SlewingMapBLprefix: %s", str1)==1 ){
      if(getenv("RUN")!=NULL)
	SlewingMapFileNameBL = Form("%s%s.param",str1,getenv("RUN"));
    }
    else if( sscanf(str,"SlewingMapSDD: %s", str1)==1 )
      SlewingMapFileNameSDD = str1;
    else if( sscanf(str,"ScalerMap: %s", str1)==1 )
      ScalerMapFileName = str1;
    else if( sscanf(str,"ReslMap: %s", str1)==1 )
      ReslMapFileName = str1;
    else if( sscanf(str,"CDSFittingParam: %s", str1)==1 )
      CDSFittingParamFileName = str1;
    else if( sscanf(str,"BLDCFittingParam: %s", str1)==1 )
      BLDCFittingParamFileName = str1;
    else if( sscanf(str,"FADCParam: %s", str1)==1 )
      FADCParamFileName = str1;
    else if( sscanf(str,"TransferMatrix: %s", str1)==1 )
      TransferMatrixFileName = str1;
    else if( sscanf(str,"DCTimeCorrBLDC: %s", str1)==1 )
      DCTimeCorrFileNameBLDC = str1;
    else if( sscanf(str,"DCTimeCorrCDC: %s", str1)==1 )
      DCTimeCorrFileNameCDC = str1;
    else if( sscanf(str,"DetectorList: %s", str1)==1 )
      DetectorListFileName = str1;
    else if( sscanf(str,"MaterialList: %s", str1)==1 )
      MaterialListFileName = str1;
    else if( sscanf(str,"DataType: %d", &typ)==1 )
      DataType = typ;
    else if( sscanf(str,"StartEvNum: %d", &evnum)==1 )
      StartEvNum = evnum;
    else if( sscanf(str,"StopEvNum: %d", &evnum)==1 )
      StopEvNum = evnum;
    else if( sscanf(str,"SkipEvNum: %d", &evnum)==1 )
      SkipEvNum.push_back(evnum);
    else if( sscanf(str,"StartBlockEvNum: %d", &evnum)==1 )
      StartBlockEvNum = evnum;
    else if( sscanf(str,"StopBlockEvNum: %d", &evnum)==1 )
      StopBlockEvNum = evnum;
    //===============================================================//
    // Store run dependent parameters
    //===============================================================//
    else if( sscanf(str,"RUN_SUMMARY: %s",str1)==1 ){
      std::string run_summary = str1;
      database->ConvertRunSummaryCVStoDB(run_summary.data()); 
    } 
    else if( sscanf(str,"CounterMap_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_CounterMap_List)     status_CounterMap_List     = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_CounterMap,     run_begin, run_end,calib_data.data());
    }
    else if( sscanf(str,"GainMapBL_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_GainMapBL_List)      status_GainMapBL_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_GainMapBL,      run_begin, run_end,calib_data.data());
    }
    else if( sscanf(str,"GateMapBL_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_GateMapBL_List)      status_GateMapBL_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_GateMapBL,      run_begin, run_end,calib_data.data());
    }
    else if( sscanf(str,"GainMapCDS_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_GainMapCDS_List)      status_GainMapCDS_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_GainMapCDS,     run_begin, run_end,calib_data.data());
    }
    else if( sscanf(str,"DCTimeCorrBLDC_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_DCTimeCorrBLDC_List)      status_DCTimeCorrBLDC_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_DCTimeCorrBLDC,  run_begin, run_end,calib_data.data());
    }
    else if( sscanf(str,"XTMapBLROOT_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_XTMapBLROOT_List)      status_XTMapBLROOT_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_XTMapBLROOT,    run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"XTMapCDS_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_XTMapCDS_List)      status_XTMapCDS_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_XTMapCDS,    run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"XTMapBL_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_XTMapBL_List)      status_XTMapBL_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_XTMapBL,    run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"SlewingMapBL_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_SlewingMapBL_List)      status_SlewingMapBL_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_SlewingMapBL,  run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"SlewingMapCDS_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_SlewingMapCDS_List)      status_SlewingMapCDS_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_SlewingMapCDS,  run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"ReslMap_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_ReslMap_List)      status_ReslMap_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_ReslMap,  run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"TransferMatrix_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_TransferMatrix_List)      status_TransferMatrix_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_TransferMatrix,  run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"CDSFittingParam_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_CDSFittingParam_List)      status_CDSFittingParam_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_CDSFittingParam,  run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"DetectorList_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_DetectorList_List)      status_DetectorList_List      = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_DetectorList,  run_begin, run_end,calib_data.data());
    } 
    else if( sscanf(str,"CDCASDMap_List: %s %d %d",str1, &run_begin, &run_end )==3 ){
      if (!status_CDCASDMap_List)     status_CDCASDMap_List     = true;  
      std::string calib_data = str1;
      database->ConvertFromFile(ID_CDCASDMap,     run_begin, run_end,calib_data.data());
    }
  }
 
  if (status_CounterMap_List) 
  {
    DB_Calibration* CounterMap     = database->FetchCalibData (ID_CounterMap,    RunNumber);
    if (CounterMap){
      CounterMapFileName.clear();    
      CounterMapFileName = CounterMap->GetContent();
    }
  }     
  if (status_GainMapBL_List)      
  {
    DB_Calibration* GainMapBL      = database->FetchCalibData (ID_GainMapBL,     RunNumber);
    if (GainMapBL) {
      GainMapFileNameBL.clear();     
      GainMapFileNameBL  = GainMapBL->GetContent();
    }
  }     
  if (status_GateMapBL_List)       
  {
    DB_Calibration* GateMapBL      = database->FetchCalibData (ID_GateMapBL,     RunNumber);
    if (GateMapBL) {
      GateMapFileNameBL.clear();     
      GateMapFileNameBL = GateMapBL->GetContent();
    }
  }     
  if (status_GainMapCDS_List)       
  {
    DB_Calibration* GainMapCDS     = database->FetchCalibData (ID_GainMapCDS,    RunNumber);
    if (GainMapCDS) {
      GainMapFileNameCDS.clear();     
      GainMapFileNameCDS = GainMapCDS->GetContent();
    }
  }     
  if (status_DCTimeCorrBLDC_List)   
  {
    DB_Calibration* DCTimeCorrBLDC = database->FetchCalibData (ID_DCTimeCorrBLDC,RunNumber);
    if (DCTimeCorrBLDC) {
      DCTimeCorrFileNameBLDC.clear(); 
      DCTimeCorrFileNameBLDC = DCTimeCorrBLDC->GetContent();
    }
  }     
  if (status_XTMapBLROOT_List)      
  {
    DB_Calibration* XTMapBLROOT    = database->FetchCalibData (ID_XTMapBLROOT,   RunNumber);
    if (XTMapBLROOT) {
      XTMapROOTNameBL.clear();       
      XTMapROOTNameBL = XTMapBLROOT->GetContent();
    }
  }     
  if (status_XTMapCDS_List)       
  {
    DB_Calibration* XTMapCDS       = database->FetchCalibData (ID_XTMapCDS,      RunNumber);
    if (XTMapCDS) {
      XTMapFileNameCDS.clear();      
      XTMapFileNameCDS = XTMapCDS->GetContent();
    }
  }     
  if (status_XTMapBL_List)       
  {
    DB_Calibration* XTMapBL       = database->FetchCalibData (ID_XTMapBL,      RunNumber);
    if (XTMapBL) {
      XTMapFileNameBL.clear();      
      XTMapFileNameBL = XTMapBL->GetContent();
    }
  }     
  if (status_SlewingMapBL_List)     
  {
    DB_Calibration* SlewingMapBL   = database->FetchCalibData (ID_SlewingMapBL,  RunNumber);
    if (SlewingMapBL) {
      SlewingMapFileNameBL.clear();  
      SlewingMapFileNameBL = SlewingMapBL->GetContent();
    }
  }     
  if (status_SlewingMapCDS_List)   
  {
    DB_Calibration* SlewingMapCDS  = database->FetchCalibData (ID_SlewingMapCDS, RunNumber);
    if (SlewingMapCDS) {
      SlewingMapFileNameCDS.clear(); 
      SlewingMapFileNameCDS = SlewingMapCDS->GetContent();
    }
  }     
  if (status_ReslMap_List)        
  {
    DB_Calibration* ReslMap        = database->FetchCalibData (ID_ReslMap,       RunNumber);
    if (ReslMap) {
      ReslMapFileName.clear();       
      ReslMapFileName = ReslMap->GetContent();
    }
  }     
  if (status_TransferMatrix_List)        
  {
    DB_Calibration* TransferMatrix        = database->FetchCalibData (ID_TransferMatrix,       RunNumber);
    if (TransferMatrix) {
      TransferMatrixFileName.clear();       
      TransferMatrixFileName = TransferMatrix->GetContent();
    }
  }     
  if (status_CDSFittingParam_List)        
  {
    DB_Calibration* CDSFittingParam        = database->FetchCalibData (ID_CDSFittingParam,       RunNumber);
    if (CDSFittingParam) {
      CDSFittingParamFileName.clear();       
      CDSFittingParamFileName = CDSFittingParam->GetContent();
    }
  }     
  if (status_DetectorList_List)        
  {
    DB_Calibration* DetectorList        = database->FetchCalibData (ID_DetectorList,       RunNumber);
    if (DetectorList) {
      DetectorListFileName.clear();       
      DetectorListFileName = DetectorList->GetContent();
    }
  }     
  if (status_CDCASDMap_List) 
  {
    DB_Calibration* CDCASDMap     = database->FetchCalibData (ID_CDCASDMap,    RunNumber);
    if (CDCASDMap){
      CDCASDMapFileName.clear();    
      CDCASDMapFileName = CDCASDMap->GetContent();
    }
  }     

  DB_RunSummary*  summary        = database->FetchRunSummaryData (RunNumber);
  if (summary) {
    summary->ShowData();
    CurrentInSolenoid = summary->GetSolenoidCurr(); 
    CurrentInUshiwaka = summary->GetUshiwakaCurr();
    BeamMomentum = summary->GetBeamMomentum();
  }
 
  delete database;
  fclose(fp);
  return true;
}
bool ConfMan::Initialize(bool PRINT, bool GEOM)
{
  try{
    static const std::string funcname = "ConfMan::Initialize";
    std::cout << "[" << funcname << "] Initialize. <--" << ConfFileName << std::endl;

    ReadConfFile(ConfFileName);
    
    if(PRINT){
      std::cout << "############################################" << std::endl;
      std::cout << "          DataType : " << DataType << std::endl;
      std::cout << "        CounterMap : " << CounterMapFileName << std::endl;
      std::cout << "       BLDCWireMap : " << BLDCWireMapFileName << std::endl;
      std::cout << "       CDCGeometry : " << CDCGeometryFileName << std::endl;
      std::cout << "         CDCASDMap : " << CDCASDMapFileName << std::endl;
      std::cout << "     CDCChannelMap : " << CDCChannelMapFileName << std::endl;
      std::cout << "        GainMapCDS : " << GainMapFileNameCDS << std::endl;
      std::cout << "         GainMapBL : " << GainMapFileNameBL << std::endl;
      //    std::cout << "        GainMapSDD : " << GainMapFileNameSDD << std::endl;
      std::cout << "        GateMapCDS : " << GateMapFileNameCDS << std::endl;
      std::cout << "         GateMapBL : " << GateMapFileNameBL << std::endl;
      //    std::cout << "        GateMapSDD : " << GateMapFileNameSDD << std::endl;
      std::cout << "          XTMapCDS : " << XTMapFileNameCDS << std::endl;
      std::cout << "           XTMapBL : " << XTMapFileNameBL << std::endl;
      std::cout << "       XTMapBLROOT : " << XTMapROOTNameBL << std::endl;
      std::cout << "       GeometryCDS : " << GeometryMapFileNameCDS << std::endl;
      std::cout << "        GeometryBL : " << GeometryMapFileNameBL << std::endl;
      std::cout << "      GeometryHall : " << GeometryMapFileNameHall << std::endl;
      std::cout << "     SlewingMapCDC : " << SlewingMapFileNameCDC << std::endl;
      std::cout << "     SlewingMapCDS : " << SlewingMapFileNameCDS << std::endl;
      std::cout << "      SlewingMapBL : " << SlewingMapFileNameBL << std::endl;
      //    std::cout << "     SlewingMapSDD : " << SlewingMapFileNameSDD << std::endl;
      std::cout << "         ScalerMap : " << ScalerMapFileName << std::endl;
      std::cout << "           ReslMap : " << ReslMapFileName << std::endl;
      std::cout << "   CDSFittingParam : " << CDSFittingParamFileName << std::endl;
      std::cout << "  BLDCFittingParam : " << BLDCFittingParamFileName << std::endl;
      //    std::cout << "         FADCParam : " << FADCParamFileName << std::endl;
      std::cout << "    TransferMatrix : " << TransferMatrixFileName << std::endl;
      std::cout << "    DCTimeCorrBLDC : " << DCTimeCorrFileNameBLDC << std::endl;
      std::cout << "     DCTimeCorrCDC : " << DCTimeCorrFileNameCDC << std::endl;
      std::cout << "      DetectorList : " << DetectorListFileName << std::endl;
      std::cout << "      MaterialList : " << MaterialListFileName << std::endl;
      std::cout << "     StartEventNum : " << StartEvNum << std::endl;
      std::cout << "      StopEventNum : " << StopEvNum << std::endl;
      for( int i=0; i<(int)SkipEvNum.size(); i++ )
	std::cout << "      SkipEventNum : " << SkipEvNum[i] << std::endl;
      std::cout << "StartBlockEventNum : " << StartBlockEvNum << std::endl;
      std::cout << " StopBlockEventNum : " << StopBlockEvNum << std::endl;
      std::cout << "############################################" << std::endl;
    }
    InitializeParameterFiles();
    std::cout<<"[ConfMan::Initialize] successfully finished !!!"<<std::endl;
  }
  catch(std::exception &e){
    std::cout<<"ConfMan::Initialize catch exception: "<<e.what()<<std::endl;
    std::cout<<"                               Type: "<<typeid(e).name()<<std::endl;
  }
  if(GEOM)
    GeomTools::MakeGeometry(this);

  struct timeval time;
  gettimeofday( &time, NULL );
  gRandom->SetSeed(time.tv_sec+time.tv_usec);
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ConfMan::InitializeParameterFiles()
{
  CounterMapManager=0;
  BLDCWireMapManager=0;
  CDCWireMapManager=0;
  GainMapManager=0;
  GateMapManager=0;
  XTMapManager=0;
  GeomMapManager=0;
  SlewingMapManager=0;
  ScalerMapManager=0;
  ReslMapManager=0;
  CDSFittingParamManager=0;
  BLDCFittingParamManager=0;
  FADCParamManager=0;
  TransferMatrixManager=0;
  DCTimeCorrManager=0;
  
  if( CounterMapFileName != DefaultFileName ){
    CounterMapManager = new CounterMapMan( CounterMapFileName );
  }
  if( CounterMapManager!=0 ) CounterMapManager->Initialize();

  if( BLDCWireMapFileName != DefaultFileName )
    BLDCWireMapManager = new BLDCWireMapMan();
  if( BLDCWireMapManager!=0 ){
    BLDCWireMapManager->SetFileName( BLDCWireMapFileName );
    BLDCWireMapManager->Initialize();
  }
  
  if( CDCGeometryFileName != DefaultFileName &&
      CDCASDMapFileName != DefaultFileName &&
      CDCChannelMapFileName != DefaultFileName )
    CDCWireMapManager = new CDCWireMapMan( CDCASDMapFileName, CDCChannelMapFileName, CDCGeometryFileName );
  if( CDCWireMapManager!=0 ) CDCWireMapManager->Initialize();
  

  if( GainMapFileNameCDS != DefaultFileName || GainMapFileNameBL  != DefaultFileName || GainMapFileNameSDD  != DefaultFileName )
    GainMapManager = new GainMapMan();
  if( GainMapManager!=0 ){
    GainMapManager->SetFileNameCDS( GainMapFileNameCDS );
    GainMapManager->SetFileNameBL(  GainMapFileNameBL );
    GainMapManager->SetFileNameSDD(  GainMapFileNameSDD );
    GainMapManager->Initialize();
  }

  GateMapManager = new GateMapMan();
  if( GateMapManager!=0 ){
    GateMapManager->SetFileNameCDS( GateMapFileNameCDS );
    GateMapManager->SetFileNameBL(  GateMapFileNameBL );
    GateMapManager->SetFileNameSDD(  GateMapFileNameSDD );
    GateMapManager->Initialize();
  }
  
  if( XTMapFileNameCDS != DefaultFileName || XTMapFileNameBL  != DefaultFileName || XTMapROOTNameBL  != DefaultFileName)
    XTMapManager = new XTMapMan();
  if( XTMapManager!=0 ){
    XTMapManager->SetFileNameCDS( XTMapFileNameCDS );
    XTMapManager->SetFileNameBL(  XTMapFileNameBL );
    XTMapManager->SetROOTNameBL(  XTMapROOTNameBL );
    XTMapManager->Initialize();
  }

  if( GeometryMapFileNameCDS != DefaultFileName || GeometryMapFileNameBL  != DefaultFileName )
    GeomMapManager = new GeomMapMan();
  if( GeomMapManager!=0 ){
    GeomMapManager->SetFileNameCDS( GeometryMapFileNameCDS );
    GeomMapManager->SetFileNameBL(  GeometryMapFileNameBL );
    GeomMapManager->SetFileNameHall(  GeometryMapFileNameHall );
    GeomMapManager->Initialize();
  }

  if( SlewingMapFileNameCDC != DefaultFileName ||  SlewingMapFileNameCDS != DefaultFileName || SlewingMapFileNameBL  != DefaultFileName || SlewingMapFileNameSDD  != DefaultFileName )
    SlewingMapManager = new SlewingMapMan();
  if( SlewingMapManager!=0 ){
    SlewingMapManager->SetFileNameCDC( SlewingMapFileNameCDC );
    SlewingMapManager->SetFileNameCDS( SlewingMapFileNameCDS );
    SlewingMapManager->SetFileNameBL(  SlewingMapFileNameBL );
    SlewingMapManager->SetFileNameSDD(  SlewingMapFileNameSDD );
    SlewingMapManager->Initialize();
  }
  
  if( ScalerMapFileName != DefaultFileName ){
    ScalerMapManager = new ScalerMapMan( ScalerMapFileName );
  }
  if( ScalerMapManager!=0 ) ScalerMapManager->Initialize();

  if( ReslMapFileName != DefaultFileName ){
    ReslMapManager = new ReslMapMan( ReslMapFileName );
  }
  if( ReslMapManager!=0 ) ReslMapManager->Initialize();

  if( CDSFittingParamFileName != DefaultFileName ){
    CDSFittingParamManager = new CDSFittingParamMan( CDSFittingParamFileName );
  }
  if( CDSFittingParamManager!=0 ) CDSFittingParamManager->Initialize();

  if( BLDCFittingParamFileName != DefaultFileName ){
    BLDCFittingParamManager = new BLDCFittingParamMan( BLDCFittingParamFileName );
  }
  if( BLDCFittingParamManager!=0 ) BLDCFittingParamManager->Initialize();

  if( FADCParamFileName != DefaultFileName ){
    FADCParamManager = new FADCParamMan( FADCParamFileName );
  }
  if( FADCParamManager!=0 ) FADCParamManager->Initialize();

  if( TransferMatrixFileName != DefaultFileName ){
    TransferMatrixManager = new TransferMatrixMan( TransferMatrixFileName );
  }
  if( TransferMatrixManager!=0 ) TransferMatrixManager->Initialize();

  if( DCTimeCorrFileNameCDC != DefaultFileName || DCTimeCorrFileNameBLDC  != DefaultFileName )
    DCTimeCorrManager = new DCTimeCorrMan();
  if( DCTimeCorrManager!=0 ){
    DCTimeCorrManager->SetFileNameCDC( DCTimeCorrFileNameCDC );
    DCTimeCorrManager->SetFileNameBLDC(  DCTimeCorrFileNameBLDC );
    DCTimeCorrManager->Initialize();
  }

  if( DetectorListFileName  != DefaultFileName ){
    DetectorList* dlist=DetectorList::GetInstance();
    dlist->Initialize(DetectorListFileName.data());
  }else{
    std::cout<< "Invalid DetectorListFileName "<<DetectorListFileName<<" !!!!!"<<std::endl;
    exit(-1);
  }
  
  if ( GetSolenoidCurrent()!=-9999 && GetUshiwakaCurrent()!=-9999 ) {
    std::cout << "Solenoid/Ushiwaka = " << GetSolenoidCurrent() << "[A], " << GetUshiwakaCurrent() << "[A] " << std::endl;
    std::cout << "magnetic field strength in Solenoid [Before] = " << CDSFittingParamManager->GetMagneticField() << std::endl;
    double field = fabs(CDSFittingParamManager->GetMagneticField());
    if ( GetSolenoidCurrent()>0 ) CDSFittingParamManager->SetMagneticField( 1.0*field);
    else if ( GetSolenoidCurrent()<0 ) CDSFittingParamManager->SetMagneticField(-1.0*field);
    else if ( GetSolenoidCurrent()==0 )CDSFittingParamManager->SetMagneticField(0);
    std::cout << "magnetic field strength in Solenoid [After]  = " << CDSFittingParamManager->GetMagneticField() << std::endl;

  } else {
    std::cout << "Solenoid/Ushiwaka current is fail to set. default value in CDSFittingParam.param will be used" << std::endl;
  }

  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int ConfMan::CheckEvNum( int evnum, int blocknum )
{
  int status = 0;
  if( 0<StartBlockEvNum && blocknum<StartBlockEvNum ) status = 1;
  if( 0<StopBlockEvNum  && StopBlockEvNum<blocknum ) status = 2;

  if( 0<StartEvNum && evnum<StartEvNum ) status = 1;
  if( 0<StopEvNum  && StopEvNum<evnum ) status = 2;

  return status;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ConfMan::End()
{
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ConfMan::SaveCode(){
  const char *dirname="src";
  char *slash = (char*)strrchr(dirname,'/');
  char *locdir;
  if (slash) locdir = slash+1;
  else       locdir = (char*)dirname;
  printf("processing dir %s\n",dirname);
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->mkdir(locdir);
  adir->cd();
  void *dirp = gSystem->OpenDirectory(dirname);
  if (!dirp) return;
  char *direntry;
  Long_t id, size,flags,modtime;
  //loop on all entries of this directory
  while ((direntry=(char*)gSystem->GetDirEntry(dirp))) {
    TString afile = Form("%s/%s",dirname,direntry);
    gSystem->GetPathInfo(afile,&id,&size,&flags,&modtime);
    if (direntry[0] == '.')             continue; //forget the "." and ".." special cases
    if (!strcmp(direntry,"CVS"))        continue; //forget some special directories
    if (!strcmp(direntry,"htmldoc"))    continue;
    if (strstr(dirname,"root/include")) continue;
    if (strstr(direntry,"G__"))         continue;
    if (strstr(direntry,".o"))          continue;
    if (strstr(direntry,".d"))          continue;
    if (strstr(direntry,"Dict.h"))      continue;
    if (strstr(direntry,"Dict.cpp"))    continue;
    if (strstr(direntry,"#"))           continue;
    if (strstr(direntry,"~"))           continue;
    if (strstr(direntry,"Dict.cpp"))    continue;
    if (strstr(direntry,".cpp")    ||
	strstr(direntry,".h")    ||
	strstr(direntry,".c")) {
      Tools::WriteFile(afile, direntry);
    } else {
      if (flags != 3)                     continue; //must be a directory
      //we have found a valid sub-directory. Process it
      //      importdir(afile);
    }
  }
  gSystem->FreeDirectory(dirp);
  savdir->cd();
  return;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ConfMan::SaveCDSParam(){
  const char* dirname="cds";
  char *slash = (char*)strrchr(dirname,'/');
  char *locdir;
  if (slash) locdir = slash+1;
  else       locdir = (char*)dirname;
  printf("processing dir %s\n",dirname);
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->mkdir(locdir);
  adir->cd();
  Tools::WriteFile(  ConfFileName        );
  Tools::WriteFile(  XTMapFileNameCDS    );
  Tools::WriteFile(  GainMapFileNameCDS  );
  Tools::WriteFile(  GateMapFileNameCDS  );
  Tools::WriteFile(  SlewingMapFileNameCDS  );
  Tools::WriteFile(  CDSFittingParamFileName);
  Tools::WriteFile(  DCTimeCorrFileNameCDC  );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ConfMan::SaveParams(){
  const char* dirname="param";
  char *slash = (char*)strrchr(dirname,'/');
  char *locdir;
  if (slash) locdir = slash+1;
  else       locdir = (char*)dirname;
  printf("processing dir %s\n",dirname);
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->mkdir(locdir);
  adir->cd();
  Tools::WriteFile(  ConfFileName  );
  Tools::WriteFile(  CounterMapFileName  );
  Tools::WriteFile(  BLDCWireMapFileName  );
  Tools::WriteFile(  CDCGeometryFileName  );
  Tools::WriteFile(  CDCASDMapFileName  );
  Tools::WriteFile(  CDCChannelMapFileName  );
  Tools::WriteFile(  GainMapFileNameCDS  );
  Tools::WriteFile(  GainMapFileNameBL  );
  Tools::WriteFile(  GainMapFileNameSDD  );
  Tools::WriteFile(  GateMapFileNameCDS  );
  Tools::WriteFile(  GateMapFileNameBL  );
  Tools::WriteFile(  GateMapFileNameSDD  );
  Tools::WriteFile(  XTMapFileNameCDS  );
  Tools::WriteFile(  XTMapFileNameBL  );
  Tools::WriteFile(  XTMapROOTNameBL  );
  Tools::WriteFile(  GeometryMapFileNameCDS  );
  Tools::WriteFile(  GeometryMapFileNameBL  );
  Tools::WriteFile(  SlewingMapFileNameCDS  );
  Tools::WriteFile(  SlewingMapFileNameBL  );
  Tools::WriteFile(  SlewingMapFileNameSDD  );
  Tools::WriteFile(  ReslMapFileName  );
  Tools::WriteFile(  CDSFittingParamFileName  );
  Tools::WriteFile(  BLDCFittingParamFileName  );
  Tools::WriteFile(  FADCParamFileName  );
  Tools::WriteFile(  TransferMatrixFileName  );
  Tools::WriteFile(  DCTimeCorrFileNameBLDC  );
  Tools::WriteFile(  DCTimeCorrFileNameCDC  );
  Tools::WriteFile(  DetectorListFileName  );
  Tools::WriteFile(  MaterialListFileName  );
  savdir->cd();
  return;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
