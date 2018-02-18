// ConfMan.h

#ifndef ConfMan_h
#define ConfMan_h 1

#include <string>
#include <vector>

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

#include "TDirectory.h"
#include "TMacro.h"
#include "TSystem.h"

class ConfMan : public TObject
{
 private:
  std::string ConfFileName;
  std::string OutFileName;
  std::string VMEFileName;

 public:
  ConfMan();
  ConfMan( const std::string &filename );
  ConfMan( const std::string &filename, const std::string &outfilename);
  ConfMan( const std::string &filename, const std::string &outfilename,
	   const std::string &vmefilename );
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
  GeomMapMan *GeomMapManager;
  std::string SlewingMapFileNameCDS;
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
  
  int DataType;
  int StartEvNum, StopEvNum;
  typedef std::vector < int > SkipEvNumContainer;
  SkipEvNumContainer SkipEvNum;

  int StartBlockEvNum, StopBlockEvNum;

 public:
  bool Initialize(bool PRINT=true);
  bool InitializeParameterFiles();
  bool ReadConfFile(std::string name);
  bool End();

  void SetConfFileName( const std::string &name ) { ConfFileName = name; }
  void SetCounterMapFileName( const std::string &name ) { CounterMapFileName = name; }
  void SetBLDCWireMapFileNames( const std::string &name ){ BLDCWireMapFileName = name; }
  void SetCDCWireMapFileNames( const std::string &geoname, const std::string &asdname,
			       const std::string &chname )
  { CDCGeometryFileName = geoname; CDCASDMapFileName = asdname; CDCChannelMapFileName = chname; }
  void SetGainMapFileNameCDS( const std::string &name ) { GainMapFileNameCDS = name; }
  void SetGainMapFileNameBL(  const std::string &name ) { GainMapFileNameBL = name; }
  void SetGainMapFileNameSDD(  const std::string &name ) { GainMapFileNameSDD = name; }
  void SetGateMapFileNameCDS( const std::string &name ) { GateMapFileNameCDS = name; }
  void SetGateMapFileNameBL(  const std::string &name ) { GateMapFileNameBL = name; }
  void SetGateMapFileNameSDD(  const std::string &name ) { GateMapFileNameSDD = name; }
  void SetXTMapFileNameCDS( const std::string &name ) { XTMapFileNameCDS = name; }
  void SetXTMapFileNameBL(  const std::string &name ) { XTMapFileNameBL = name; }
  void SetXTMapROOTNameBL(  const std::string &name ) { XTMapROOTNameBL = name; }
  void SetGeomMapFileNameCDS( const std::string &name ) { GeometryMapFileNameCDS = name; }
  void SetGeomMapFileNameBL(  const std::string &name ) { GeometryMapFileNameBL = name; }
  void SetSlewingMapFileNameCDS( const std::string &name ) { SlewingMapFileNameCDS = name; }
  void SetSlewingMapFileNameBL(  const std::string &name ) { SlewingMapFileNameBL = name; }
  void SetSlewingMapFileNameSDD(  const std::string &name ) { SlewingMapFileNameSDD = name; }
  void SetScalerMapFileName( const std::string &name ) { ScalerMapFileName = name; }
  void SetReslMapFileName( const std::string &name ) { ReslMapFileName = name; }
  void SetCDSFittingParamFileName( const std::string &name ) { CDSFittingParamFileName = name; }
  void SetBLDCFittingParamFileName( const std::string &name ) { BLDCFittingParamFileName = name; }
  void SetFADCParamFileName( const std::string &name ) { FADCParamFileName = name; }
  void SetTransferMatrixFileName( const std::string &name ) { TransferMatrixFileName = name; }

  void SetCounterMapManager( CounterMapMan *map ) { CounterMapManager = map; }
  void SetBLDCWireMapManager( BLDCWireMapMan *map ) { BLDCWireMapManager = map; }
  void SetCDCWireMapManager( CDCWireMapMan *map ) { CDCWireMapManager = map; }
  void SetGainMapManager( GainMapMan *map ) { GainMapManager = map; }
  void SetGateMapManager( GateMapMan *map ) { GateMapManager = map; }
  void SetXTMapManager( XTMapMan *map ) { XTMapManager = map; }
  void SetGeomMapManager( GeomMapMan *map ) { GeomMapManager = map; }
  void SetSlewingMapManager( SlewingMapMan *map ) { SlewingMapManager = map; }
  void SetScalerMapManager( ScalerMapMan *map ) { ScalerMapManager = map; }
  void SetReslMapManager( ReslMapMan *map ) { ReslMapManager = map; }
  void SetCDSFittingParamManager( CDSFittingParamMan *map ) { CDSFittingParamManager = map; }
  
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
  FADCParamMan *GetFADCParamManager() { return FADCParamManager; }
  TransferMatrixMan *GetTransferMatrixManager() { return TransferMatrixManager; }
  
  std::string GetOutFileName() const { return OutFileName; }
  std::string GetVMEFileName() const { return VMEFileName; }
  std::string GetConfFileName() const { return ConfFileName; }
  std::string GetCounterMapFileName() const { return CounterMapFileName; }
  std::string GetBLDCWireMapFileName() const { return BLDCWireMapFileName; }
  std::string GetCDCGeometryFileName() const { return CDCGeometryFileName; }
  std::string GetCDCASDMapFileName() const { return CDCASDMapFileName; }
  std::string GetCDCChannelMapFileName() const { return CDCChannelMapFileName; }
  std::string GetGainMapFileNameCDS() const { return GainMapFileNameCDS; }
  std::string GetGainMapFileNameBL() const { return GainMapFileNameBL; }
  std::string GetGainMapFileNameSDD() const { return GainMapFileNameSDD; }
  std::string GetGateMapFileNameCDS() const { return GateMapFileNameCDS; }
  std::string GetGateMapFileNameBL() const { return GateMapFileNameBL; }
  std::string GetGateMapFileNameSDD() const { return GateMapFileNameSDD; }
  std::string GetXTMapFileNameCDS() const { return XTMapFileNameCDS; }
  std::string GetXTMapFileNameBL() const { return XTMapFileNameBL; }
  std::string GetXTMapROOTNameBL() const { return XTMapROOTNameBL; }
  std::string GetGeometryFileNameCDS() const { return GeometryMapFileNameCDS; }
  std::string GetGeometryFileNameBL() const { return GeometryMapFileNameBL; }
  std::string GetSlewingMapFileNameCDS() const { return SlewingMapFileNameCDS; }
  std::string GetSlewingMapFileNameBL() const { return SlewingMapFileNameBL; }
  std::string GetSlewingMapFileNameSDD() const { return SlewingMapFileNameSDD; }
  std::string GetScalerMapFileName() const { return ScalerMapFileName; }
  std::string GetReslMapFileName() const { return ReslMapFileName; }
  std::string GetCDSFittingParamFileName() const { return CDSFittingParamFileName; }
  std::string GetBLDCFittingParamFileName() const { return BLDCFittingParamFileName; }
  std::string GetFADCParamFileName() const { return FADCParamFileName; }
  std::string GetTransferMatrixFileName() const { return TransferMatrixFileName; }
  int GetDataType() const { return DataType; }
  int GetStartEvNum() const { return StartEvNum; }
  int GetStopEvNum() const { return StopEvNum; }
  int GetNSkipEvNum() const { return SkipEvNum.size(); }
  int GetSkipEvNum( int i ) const { return SkipEvNum[i]; }
  int GetStartBlockEvNum() const { return StartBlockEvNum; }
  int GetStopBlockEvNum() const { return StopBlockEvNum; }

  int CheckEvNum( int evnum, int blocknum );
  void SaveParams();
  void SaveCode();
  void WriteFile(TString a, TString b=0);

  ClassDef( ConfMan, 1 );
};

#endif
