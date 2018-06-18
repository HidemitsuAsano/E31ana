/* DB_Manager.h */
#ifndef DB_Manager_h
#define DB_Manager_h 1

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <time.h>

#include "TObject.h"
#include "TTimeStamp.h"
#include "TString.h"

typedef enum  { ID_CounterMap, ID_BLDCWireMap, ID_CDCGeom,        ID_CDCASDMap,   ID_CDCChannelMap, ID_GainMapCDS, //6
                ID_GainMapBL,  ID_GateMapBL,   ID_DCTimeCorrBLDC, ID_XTMapCDS,    ID_XTMapBLROOT, ID_XTMapBL, //12
                ID_GeomCDS,    ID_GeomBL,      ID_SlewingMapCDS,  ID_SlewingMapBL,  //16
                ID_ScalerMap,  ID_TransferMatrix, ID_ReslMap, ID_DetectorList, ID_CDSFittingParam
} Calib_type; 

typedef std::vector <TString> CalibDataCont;
class DB_Calibration : public TObject
{
 private:
  CalibDataCont CaibrationData;
  TString       CalibName;
  int           ValidFrom; 
  int           ValidTo; 
  TTimeStamp    inserted_time;
  
 public:
  DB_Calibration();
  DB_Calibration(TString name, int start, int stop, TTimeStamp t0, CalibDataCont data);
  ~DB_Calibration(){}
  int  GetBegin() {return ValidFrom;}
  int  GetEnd()   {return ValidTo;  }   
  void DumpCalibData(int run_id);
  const char* GetContent(){return CaibrationData[0].Data();}

  void       ShowData();  

  ClassDef(DB_Calibration, 1)

};

class DB_RunSummary : public TObject
{
 private:
  int        run_number;
  TTimeStamp start_time;
  TTimeStamp end_time;  
  float      MR_Power; 
  float      K18_Momentum;
  TString    BeamParticle;
  TString    BeamPolarization;
  TString    TargetMaterial;
  int        CurrentInSolenoid;
  int        CurrentInUshiwaka;
  TString    TrigCondition;
  TString    Comments;      
  TTimeStamp inserted_time;

 public:
  DB_RunSummary();
  DB_RunSummary(int   i0,  TTimeStamp t0,  TTimeStamp t1, TTimeStamp t2, float f0, float f1,
                char* st0, char* st1,  char* st2, char* st3, int i1, int i2, char* st4);
  ~DB_RunSummary();

  void        SetRunNumber(int i)       {run_number       =  i;}
  void        SetMR_Power(float f)      {MR_Power         =  f;} 
  void        SetBeamMomentum(float f)  {K18_Momentum     =  f;}
  void        SetBeamParticle(char* st) {BeamParticle     = st;}
  void        SetBeamPol(char* st)      {BeamPolarization = st;}
  void        SetTargetType(char* st)   {TargetMaterial   = st;}
  void        SetSolenoidCurr(int i)    {CurrentInSolenoid=  i;}   
  void        SetUshiwakaCurr(int i)    {CurrentInUshiwaka=  i;}   
  void        SetTrigCondition(char* st){TrigCondition    = st;}
  void        SetComments(char* st)     {Comments         = st;}

  int         GetRunNumber()    {return run_number;}
  float       GetMR_Power()     {return MR_Power;} 
  float       GetBeamMomentum() {return K18_Momentum;}
  const char* GetBeamParticle() {return BeamParticle.Data();}
  const char* GetBeamPol()      {return BeamPolarization.Data();}
  const char* GetTargetType()   {return TargetMaterial.Data();}
  const char* GetTrigCondition(){return TrigCondition.Data();}
  const char* GetComments()     {return Comments.Data();}
  int         GetSolenoidCurr() {return CurrentInSolenoid;}   
  int         GetUshiwakaCurr() {return CurrentInUshiwaka;}   

  void        ShowData();  
  void        CreateCDSFittingParameter();
  void        CreateBLDCFittingParam();
  void        CreateAnalysisConf();

  ClassDef(DB_RunSummary, 1)
};

typedef std::vector <DB_RunSummary>    RunSummaryContainer;
typedef std::vector <DB_Calibration>   CalibrationContainer;
class DB_Manager : public TObject
{
private:
  RunSummaryContainer RunSummaryCont;
  CalibrationContainer CalibData[30];

  int tranceport(char *a, char *b); 

 public:
  DB_Manager();
  virtual ~DB_Manager(); 
  

//
// RunSummary
//
  void ConvertRunSummaryCVStoDB(const char* st="./RunSummary_run49c.csv");
  int  GetEntriesRunSummary(){return RunSummaryCont.size();    }
  DB_RunSummary*   FetchRunSummaryData(int run_number);  

//
//  Calibration Data
//
  void ConvertFromFile(Calib_type id, int i0=0, int i1=0, const char* st="calibdata.dat");
  int  GetEntriesCalibData(Calib_type id) {return CalibData[id].size();}
  DB_Calibration*  FetchCalibData (Calib_type id, int run_number);

  const char* GetCalibName(Calib_type id);

  ClassDef(DB_Manager, 1)

};

#endif
