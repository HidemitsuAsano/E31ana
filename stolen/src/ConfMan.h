// ConfMan.h

#ifndef ConfMan_h
#define ConfMan_h 1

#include <string>
#include <vector>

#include <math.h>

#include "TObject.h"
#include "GlobalVariables.h"
#include "CounterMapMan.h"
#include "CDCWireMapMan.h"
#include "BLDCWireMapMan.h"
#include "GainMapMan.h"
#include "GateMapMan.h"
#include "XTMapMan.h"
#include "GeomMapMan.h"
#include "SlewingMapMan.h"
#include "ScalerMapMan.h"
#include "ReslMapMan.h"
#include "CDSFittingParamMan.h"
#include "BLDCFittingParamMan.h"
#include "FADCParamMan.h"
#include "TransferMatrixMan.h"
#include "DCTimeCorrMan.h"
#include "DetectorList.h"

#include "DB_Manager.h"

#include "TDirectory.h"
#include "TMacro.h"
#include "TSystem.h"

class ConfMan : public TObject
{
 private:
  std::string ProgramName;
  std::string ConfFileName;
  std::string OutFileName;
  std::string CDSTrackFileName;
  std::string MTDCFileName;

  int   RunNumber;
  float   BeamMomentum;
  int   CurrentInSolenoid;
  int   CurrentInUshiwaka;

 public:
  ConfMan();
  ConfMan( const std::string &filename );
  ConfMan( const std::string &filename, int runnum);
  ConfMan( const std::string &filename, int runnum, const std::string &outfilename);
  ConfMan( const std::string &filename, int runnum, const std::string &outfilename,
	   const std::string &cdstrackfilename );
  ConfMan( const std::string &filename, int runnum, const std::string &outfilename,
	   const std::string &cdstrackfilename, const std::string &mtdcfilename );
  ~ConfMan();

 private:
  ConfMan( const ConfMan &conf );

 private:
  std::string CounterMapFileName;
  CounterMapMan *CounterMapManager;

  std::string BLDCWireMapFileName;
  BLDCWireMapMan *BLDCWireMapManager;

  std::string CDCGeometryFileName;
  std::string CDCASDMapFileName;
  std::string CDCChannelMapFileName;
  CDCWireMapMan *CDCWireMapManager;

  std::string GainMapFileNameCDS;
  std::string GainMapFileNameBL;
  std::string GainMapFileNameSDD;
  GainMapMan *GainMapManager;

  std::string GateMapFileNameCDS;
  std::string GateMapFileNameBL;
  std::string GateMapFileNameSDD;
  GateMapMan *GateMapManager;

  std::string XTMapFileNameCDS;
  std::string XTMapFileNameBL;
  std::string XTMapROOTNameBL;
  XTMapMan *XTMapManager;

  std::string GeometryMapFileNameCDS;
  std::string GeometryMapFileNameBL;
  std::string GeometryMapFileNameHall;
  GeomMapMan *GeomMapManager;

  std::string SlewingMapFileNameCDS;
  std::string SlewingMapFileNameCDC;
  std::string SlewingMapFileNameBL;
  std::string SlewingMapFileNameSDD;
  SlewingMapMan *SlewingMapManager;

  std::string ScalerMapFileName;
  ScalerMapMan *ScalerMapManager;

  std::string ReslMapFileName;
  ReslMapMan *ReslMapManager;

  std::string CDSFittingParamFileName;
  CDSFittingParamMan *CDSFittingParamManager;

  std::string BLDCFittingParamFileName;
  BLDCFittingParamMan *BLDCFittingParamManager;

  std::string FADCParamFileName;
  FADCParamMan *FADCParamManager;

  std::string TransferMatrixFileName;
  TransferMatrixMan *TransferMatrixManager;

  std::string DCTimeCorrFileNameCDC;
  std::string DCTimeCorrFileNameBLDC;
  DCTimeCorrMan *DCTimeCorrManager;

  std::string DetectorListFileName;  
  std::string MaterialListFileName;  

  int DataType;
  int StartEvNum, StopEvNum;
  typedef std::vector < int > SkipEvNumContainer;
  SkipEvNumContainer SkipEvNum;

  int StartBlockEvNum, StopBlockEvNum;

 public:
  bool Initialize(bool PRINT=true,bool GEOM=true);
  bool InitializeParameterFiles();
  bool ReadConfFile(std::string name);
  bool End();
  
  void SetProgramName( const std::string &name ) { ProgramName = name; }
  void SetConfFileName( const std::string &name ) { ConfFileName = name; }
  void SetOutFileName( const std::string &name ) { OutFileName = name; }
  void SetCDSTrackFileName( const std::string &name ) { CDSTrackFileName = name; }
  void SetMTDCFileName( const std::string &name ) { MTDCFileName = name; }

  std::string GetProgramName(      ) const { return ProgramName; }
  //  char*       GetProgramName(      ) const { return ProgramName.data(); }
  std::string GetConfFileName(     ) const { return ConfFileName; }
  std::string GetOutFileName(      ) const { return OutFileName; }
  std::string GetCDSTrackFileName( ) const { return CDSTrackFileName; }
  std::string GetMTDCTrackFileName() const { return MTDCFileName; }
  
  void SetDataType( const int &type ) { DataType = type; }
  void SetStartStopEventNum( const int &start, const int &stop ) { StartEvNum = start; StopEvNum = stop; }
  void SetStartStopBlockNum( const int &start, const int &stop ) { StartBlockEvNum = start; StopBlockEvNum = stop; }
  void SetSkipEventNum( const std::vector<int> &skip ) { SkipEvNum = skip; }

  CounterMapMan  *GetCounterMapManager() {  return CounterMapManager; }
  BLDCWireMapMan *GetBLDCWireMapManager() { return BLDCWireMapManager; }
  CDCWireMapMan  *GetCDCWireMapManager() {  return CDCWireMapManager; }
  GainMapMan     *GetGainMapManager() {     return GainMapManager; }
  GateMapMan     *GetGateMapManager() {     return GateMapManager; }
  XTMapMan       *GetXTMapManager() {       return XTMapManager; }
  GeomMapMan     *GetGeomMapManager() {     return GeomMapManager; }
  SlewingMapMan  *GetSlewingMapManager() {  return SlewingMapManager; }
  ScalerMapMan   *GetScalerMapManager() {   return ScalerMapManager; }
  ReslMapMan     *GetReslMapManager() {     return ReslMapManager; }
  CDSFittingParamMan *GetCDSFittingParamManager() { return CDSFittingParamManager; }
  BLDCFittingParamMan *GetBLDCFittingParamManager() { return BLDCFittingParamManager; }
  FADCParamMan   *GetFADCParamManager() { return FADCParamManager; }
  TransferMatrixMan *GetTransferMatrixManager() { return TransferMatrixManager; }
  DCTimeCorrMan  *GetDCTimeCorrManager() { return DCTimeCorrManager; }
 
  void SetRunNumber(int irun){RunNumber = irun;}
  void SetSolenoidCurrent(int i){CurrentInSolenoid = i;}
  void SetUshiwakaCurrent(int i){CurrentInUshiwaka = i;}
  void SetBeamMomentum(float i){BeamMomentum = i;}

  int  GetRunNumber()       const {return RunNumber;}
  int  GetSolenoidCurrent() const {return CurrentInSolenoid;}
  int  GetUshiwakaCurrent() const {return CurrentInUshiwaka;}
  int  GetBeamMomentum()    const {return BeamMomentum;}
 
  int GetDataType() const { return DataType; }
  int GetStartEvNum() const { return StartEvNum; }
  int GetStopEvNum() const { return StopEvNum; }
  int GetNSkipEvNum() const { return SkipEvNum.size(); }
  int GetSkipEvNum( int i ) const { return SkipEvNum[i]; }
  int GetStartBlockEvNum() const { return StartBlockEvNum; }
  int GetStopBlockEvNum() const { return StopBlockEvNum; }

  int  CheckEvNum( int evnum, int blocknum );
  void SaveParams();
  void SaveCode();
  void SaveCDSParam();
  void Clear();

  ClassDef( ConfMan, 1 );
};

#endif
