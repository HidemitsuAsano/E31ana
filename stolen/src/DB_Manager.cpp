#include "DB_Manager.h"

ClassImp(DB_Manager)

ClassImp(DB_RunSummary)

ClassImp(DB_Calibration)

DB_Calibration::DB_Calibration() : TObject()
{
  CaibrationData.clear();
  CalibName.Clear(); 
  ValidFrom     = -999;
  ValidTo       = -999;
  inserted_time.Set(2000,1,1,0,0,0,0,0,0);
}

DB_Calibration::DB_Calibration(TString name, int start, int stop, TTimeStamp t0, CalibDataCont data)
{
  CalibName      = name;
  ValidFrom      = start;
  ValidTo        = stop;
  inserted_time  = t0;
  CaibrationData = data;
}

void DB_Calibration::DumpCalibData(int run_number)
{
  // 
  // Creating file name
  //
  char number[4];
  sprintf(number,"%4.4d",run_number);

  TString dname = "./param_run";
  dname += number;

  mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
  if ( mkdir(dname.Data(), mode) != 0 )
   {
     //printf("directry %s exist!!!\n", dname.Data());
   }
  //
  //
  TString fname = dname;
  fname += "/E15_Calib_"; 
  fname += CalibName.Data();
  fname += ".param";

  int ith = 0;
  int i_entries = (int)CaibrationData.size();
  
  std::cout << "Calib Type = " << CalibName.Data() << ":            " << fname.Data() << std::endl;
  if (strcmp(CalibName.Data(),"XTMapBLROOT")!=0 && strcmp(CalibName.Data(),"DCTimeCorrBLDC")!=0  ) { 
    std::fstream out(fname.Data(),std::ios::out);
    if (run_number >= GetBegin() && 
        run_number <= GetEnd()      ) {
       for (int i=0; i<i_entries; i++){
         ith++;
         out << CaibrationData[i].Data() ;
         //std::cout << CaibrationData[i].Data() ;
       } 
       return;
    }
  } else {
    if (run_number >= GetBegin() && 
        run_number <= GetEnd()      ) {
       for (int i=0; i<i_entries; i++){
         ith++;
         //CaibrationData[i].Data() ;
         //std::cout << CaibrationData[i].Data() ;
         std::string command = "cp ";
         command += CaibrationData[i].Data();
         command += " "; 
         command += dname;
         command += "/E15_Calib_";
         command += CalibName.Data();
         command += ".root";
         //std::cout << command << std::endl;
         system( command.data() );
       } 
       return;
    }
  }
  if (ith>=i_entries)  {
    printf("GainMapBLData :: No Data in Database for RUN%4.4d\n",run_number);
  }
    
}

DB_RunSummary::DB_RunSummary(){
  run_number     =  -999;
  start_time.Set   (2000,1,1,0,0,0,0,0,0);
  end_time.Set     (2000,1,1,0,0,0,0,0,0);  
  inserted_time.Set(2000,1,1,0,0,0,0,0,0);  
  MR_Power         = 0.0;
  K18_Momentum     = 0.0;
  BeamParticle     = "Higgs";
  BeamPolarization = "neutral";
  TargetMaterial   = "Gold";
  TrigCondition    = "slepton";
  Comments         = "Default!!!!";
  CurrentInSolenoid=0;
  CurrentInUshiwaka=0;
}

DB_RunSummary::DB_RunSummary(int   i0,  TTimeStamp t0,  TTimeStamp t1, TTimeStamp t2, float f0, float f1,
                char* st0, char* st1,  char* st2, char* st3, int i1, int i2, char* st4)
{
  run_number       = i0;
  start_time       = t0;
  end_time         = t1;  
  inserted_time    = t2;  
  MR_Power         = f0;
  K18_Momentum     = f1;
  BeamParticle     = st0;
  BeamPolarization = st1;
  TargetMaterial   = st2;
  TrigCondition    = st3;
  CurrentInSolenoid= i1;
  CurrentInUshiwaka= i2;
  Comments         = st4;
}

DB_RunSummary::~DB_RunSummary()
{
}

void DB_RunSummary::ShowData(){
  printf("##################################################################\n");
  printf("Run Number       = %4d\n",run_number);
  printf("##################################################################\n");
  printf("Run Started      = %s \n",start_time.AsString("lc"));
  printf("Run Ended        = %s \n",end_time.AsString("lc"));
  printf("Trig Condition   = %s \n",TrigCondition.Data());  
  printf("Beam Particle    = %s \n",BeamParticle.Data());
  printf("Beam Polarization= %s \n",BeamPolarization.Data());
  printf("Target           = %s \n",TargetMaterial.Data());
  printf("Solenoid         = %6d[A]\n",CurrentInSolenoid);
  printf("Ushiwaka         = %6d[A]\n",CurrentInUshiwaka);
  printf("##################################################################\n");
  printf("Data inserted on = %s \n",inserted_time.AsString("lc"));
  printf("##################################################################\n");
}

DB_Manager::DB_Manager() : TObject () 
{
  RunSummaryCont.clear();
}

DB_Manager::~DB_Manager() 
{
}

void DB_Manager::ConvertRunSummaryCVStoDB(const char* st)
{
   FILE *fp;
   std::string fname(st);
   char  line_tmp[512];
   char  line    [512];

   fp = fopen( fname.data(), "r" );
   if( fp == NULL ){
     printf( "cannot open file = %s\n", fname.data() );
     exit(1);
   }

   int ncol = 55; // number of cols in RUN49c run summary table 

   int ith   = 0;
   int nelem = 0;
   while ( fgets(line_tmp, 512, fp)!=NULL) {
     char  run_id[8];     char  dummy[24];     char  start_time[36]; char  end_time [36];
     char  AccPow[8];     char  ExTarg[8];     char  HDTarg[8];      char  BeamPol[8];
     char  K18D1[8];      char  K18Mom[8];     char  K11D1[8];       char  ESS[8];
     char  CM[8];         char  Beam[8];       char  Trigger[128];   char  TrigReq[32];
     char  TrigAcc[32];   char  Solenoid[8];   char  Ushiwaka[8];    char  FT[12];
     char  SYIM[12];      char  TM[12];        char  nbeam[12];      char  kaon[12];
     char  pion[12];      char  MomSp[12];     char  MomSm[12];      char  MSp[12];
     char  MSm[12];       char  IFHp[12];      char  IFHm[12];       char  IFVp[12];
     char  IFVm[12];      char  DEFinTrig[12]; char  comment[128];

     std::string command("%[^,]");

     for (int i=0; i<(ncol-1); i++){ command.append(",%[^,]"); }

     std::string command_insert;

     nelem = tranceport(line_tmp, line);
     if ( nelem != ncol ) {
       std::cout << "Data structure for run summary is inconsistent "   << std::endl;
       std::cout << "with RUN49c run summary table structure (default)" << std::endl;
       std::cout << "Read data from CSV might be inncorrect. " << std::cout;
       std::cout << "Please check data structure " << std::endl; exit(1);
     }

     ith++;
     sscanf(line,command.c_str(),
         run_id,dummy,dummy,dummy,dummy,start_time,end_time,dummy,AccPow,ExTarg,HDTarg,BeamPol,K18D1,
         K18Mom,K11D1,ESS,CM, Beam,Trigger,TrigReq,TrigAcc,Solenoid,Ushiwaka,FT,SYIM,TM,nbeam,kaon,
         pion,MomSp,MomSm,MSp,MSm,IFHp,IFHm,IFVp,IFVm,DEFinTrig,comment,
         dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy );

      std::string ushiwaka_cur(Ushiwaka);

      std::string s1(start_time);
      std::string s2(end_time);
 
      std::vector<char> v1(s1.begin(), s1.end());   std::vector<char> v2(s2.begin(), s2.end());
      std::replace(v1.begin(), v1.end(), '/', ' '); std::replace(v2.begin(), v2.end(), '/', ' ');
      std::replace(v1.begin(), v1.end(), ':', ' '); std::replace(v2.begin(), v2.end(), ':', ' ');
      std::string tmpTime1(v1.begin(), v1.end());   std::string tmpTime2(v2.begin(), v2.end());

      char start_1[5],start_2[5],start_3[5],start_4[5],start_5[5];
      char end_t_1[5],end_t_2[5],end_t_3[5],end_t_4[5],end_t_5[5];
      sscanf(tmpTime1.data(),"%s %s %s %s %s ",start_1,start_2,start_3,start_4,start_5);
      sscanf(tmpTime2.data(),"%s %s %s %s %s ",end_t_1,end_t_2,end_t_3,end_t_4,end_t_5);

      TTimeStamp st_time(atoi(start_1),atoi(start_2),atoi(start_3),atoi(start_4),atoi(start_5), 0, kFALSE, 0);
      TTimeStamp sp_time(atoi(end_t_1),atoi(end_t_2),atoi(end_t_3),atoi(end_t_4),atoi(end_t_5), 0, kFALSE, 0);
      TTimeStamp in_time;

      int UshiwakaCur = 0;
      if ( strcmp(ushiwaka_cur.data()," scan") != 0 ) {
	UshiwakaCur = atoi(Ushiwaka);
      } else {
        UshiwakaCur = 0; 
      }

      DB_RunSummary run_summary(atoi(run_id),st_time,sp_time,in_time,atof(AccPow),atof(K18Mom),
                                Beam,BeamPol,ExTarg,Trigger,atoi(Solenoid),atoi(Ushiwaka),comment);
      RunSummaryCont.push_back(run_summary);
   }   
   
}

void DB_Manager::ConvertFromFile(Calib_type id, int i0, int i1, const char* st)
{
   std::string fname("");
   fname+=st;
#if 0
   std::cout << "Reading ! " << fname << std::endl;
#endif
   //FILE *fp;
   char  line_tmp[128];

   typedef std::vector <TString> CalibDataCont;
   CalibDataCont calData;

   TString name(GetCalibName(id));
   TTimeStamp in_time;

   sprintf(line_tmp,"%s",fname.data());
   calData.push_back(line_tmp);

   //   std::cout<<name<<"  "<<id<<"  "<<i0<<"  "<<i1<<"  "<<line_tmp<<std::endl;    
   DB_Calibration DB_Data(name ,i0, i1, in_time, calData);
   CalibData[id].push_back(DB_Data);

}

DB_RunSummary* DB_Manager::FetchRunSummaryData(int run_number){
  int ith = 0;
  int nRunsInDB = (int)RunSummaryCont.size(); 
  for (int i=0; i<nRunsInDB; i++){
     ith++;
     if (run_number == RunSummaryCont[i].GetRunNumber()) {
       return &RunSummaryCont[i];
     }
  }
  if (ith>=nRunsInDB)  {
    printf("RunSummary :: No Data in Database for RUN%4.4d\n",run_number); 
    printf("Default value will be used for run info.\n");
    return 0;
  }
  return 0;
}

DB_Calibration* DB_Manager::FetchCalibData(Calib_type id, int run_number){
  int ith = 0;
  int i_entries = (int)CalibData[id].size();
#if 0
  printf("DB_Manager::Data in DB %20s : entries = %3d\n",GetCalibName(id),i_entries);
#endif
  for (int i=0; i<i_entries; i++){
    ith++;
    if (run_number >= CalibData[id][i].GetBegin() && 
        run_number <= CalibData[id][i].GetEnd()      ) {
       return &CalibData[id][i];
    }
  }
  if (ith>=i_entries)  {
    printf("DB_Calib  :: No Data in Database for RUN%4.4d\n",run_number);
    printf("Default value will be used for calibration parameters\n");
    return 0;
  }
  return 0;
}

const char* DB_Manager::GetCalibName(Calib_type id){
  std::string name;
  if      (id == ID_CounterMap     ) { name = "CounterMap";    }
  else if (id == ID_BLDCWireMap    ) { name = "BLDCWireMap";   }
  else if (id == ID_CDCGeom        ) { name = "CDCGeom";       }
  else if (id == ID_CDCASDMap      ) { name = "CDCASDMap";     }
  else if (id == ID_CDCChannelMap  ) { name = "CDCChannelMap"; } 
  else if (id == ID_GainMapCDS     ) { name = "GainMapCDS";    } 
  else if (id == ID_GainMapBL      ) { name = "GainMapBL";     } 
  else if (id == ID_GateMapBL      ) { name = "GateMapBL";     } 
  else if (id == ID_DCTimeCorrBLDC ) { name = "DCTimeCorrBLDC";} 
  else if (id == ID_XTMapCDS       ) { name = "XTMapCDS";      } 
  else if (id == ID_XTMapBLROOT    ) { name = "XTMapBLROOT";   } 
  else if (id == ID_XTMapBL        ) { name = "XTMapBL";   } 
  else if (id == ID_GeomCDS        ) { name = "GeomCDS";       } 
  else if (id == ID_GeomBL         ) { name = "GeomBL";        } 
  else if (id == ID_SlewingMapCDS  ) { name = "SlewingMapCDS"; } 
  else if (id == ID_SlewingMapBL   ) { name = "SlewingMapBL";  } 
  else if (id == ID_ScalerMap      ) { name = "ScalerMap";     } 
  else if (id == ID_TransferMatrix ) { name = "TransferMatrix";} 
  else if (id == ID_ReslMap        ) { name = "ReslMap";       }
  else if (id == ID_DetectorList   ) { name = "DetectorList";       }
  else if (id == ID_CDSFittingParam) { name = "CDSFittingParam";  }
  else                               { std::cout << "No Calibration for " << id << " are available"  << std::endl;}               

  //std::cout <<  id << " " << name.Data() << std::endl;
  return name.data(); 
}


void DB_Calibration::ShowData(){
  
  printf("##################################################################\n");
  printf("DB_NAME          = %s \n",CalibName.Data());
  printf("##################################################################\n");
  int ncol = (int)CaibrationData.size();
  for (int i=0; i<ncol; i++){
    std::cout << CaibrationData[i].Data() << std::endl;
  }
  printf("##################################################################\n");
  printf("Data inserted on = %s \n",inserted_time.AsString("lc"));
  printf("##################################################################\n");

}

void DB_RunSummary::CreateBLDCFittingParam(){

  // 
  // Creating file name
  //
  char number[4];
  sprintf(number,"%4.4d",run_number);
  //
  TString dname = "./param_run";
  dname += number;
  //
  mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
  if ( mkdir(dname.Data(), mode) != 0 )
   {
     //printf("directry %s exist!!!\n", dname.Data());
    }
  TString fname = dname;
  fname += "/E15_Calib_BLDCFittingParam.param";

  std::fstream out(fname.Data(),std::ios::out);
 
  out << " #        \n";
  out << " # BLDCFittingParam        \n";
  out << " #        \n";
  out << " #        \n";
  out << " ########        \n";
  out << " # MaxNumber of BLDCHit for analysis        \n";
  out << " MaxBLDCHit: 48        \n";
  out << " MaxHitInLayer: 16        \n";
  out << " MaxHitInTrack: 16        \n";
  out << " ########        \n";
  out << " # Limit of BLDCTDC#        \n";
  out << " UpperLimitBLDCTDC: 850        \n";
  out << " LowerLimitBLDCTDC: 350        \n";
  out << " #        \n";
  out << " ########        \n";
  out << " # Limit of Chi#        \n";
  out << " MaxChi: 500        \n";
  out << " MaxChiPreTracking: 5        \n";
  out << " MaxSlope: 0.4        \n";
  out << " #        \n";
  out << " ########        \n";
  out << " ######        \n";
  out << " #Layer Condition: 1(alive),0(killed),-1(dead)        \n";
  out << " LayerBLC1a: 1 1 1 1 1 1 1 1        \n";
  out << " LayerBLC1b: 1 1 1 1 1 1 1 1        \n";
  out << " LayerBLC2a: 1 1 1 1 1 1 1 1        \n";
  out << " LayerBLC2b: 1 1 1 1 1 1 1 1        \n";
  out << " LayerBPC:   1 1 1 1 1 1 1 1        \n";
  out << " LayerFDC1:  1 1 1 1 1 1 1 1        \n";
  out << " ########        \n";
  out << " # MinNumber for tracking        \n";
  out << " MinHitXBLC1: 5        \n";
  out << " MinHitXBLC1a: 3        \n";
  out << " MinHitXBLC1b: 3        \n";
  out << " MinHitXBLC2: 5        \n";
  out << " MinHitXBLC2a: 3        \n";
  out << " MinHitXBLC2b: 3        \n";
  out << " MinHitXBPC:  3        \n";
  out << " MinHitYBLC1: 5        \n";
  out << " MinHitYBLC1a: 3        \n";
  out << " MinHitYBLC1b: 3        \n";
  out << " MinHitYBLC2: 5        \n";
  out << " MinHitYBLC2a: 3        \n";
  out << " MinHitYBLC2b: 3        \n";
  out << " MinHitYBPC:  3        \n";
  out << " MinHitXFDC1: 6        \n";

}

void DB_RunSummary::CreateCDSFittingParameter()
{
  float magnetic_field = (float)CurrentInSolenoid/983.*0.7; 

  // 
  // Creating file name
  //
  char number[4];
  sprintf(number,"%4.4d",run_number);
  //
  TString dname = "./param_run";
  dname += number;
  //
  mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
  if ( mkdir(dname.Data(), mode) != 0 )
   {
     //printf("directry %s exist!!!\n", dname.Data());
    }
  TString fname = dname;
  fname += "/E15_Calib_CDSFittingParam.param";

  std::fstream out(fname.Data(),std::ios::out);
 
  out << "########                                \n";
  out << "#                                       \n";
  out << "# CDSFittingParam                       \n";
  out << "#                                       \n";
  out << "########                                \n";
  out << "# MaxNumber of CDCHit for analysis      \n";
  out << "MaxCDCHit: 150                          \n";
  out << "########                                \n";
  out << "# Limit of CDCTDC#                      \n";
  out << "UpperLimitCDCTDC: 950                   \n";
  out << "LowerLimitCDCTDC: 350                   \n";
  out << "#                                       \n";
  out << "########                                \n";
  out << "# Limit of Chi#                         \n";
  out << "MaxChi: 100                             \n";
  out << "#                                       \n";
  out << "########                                \n";
  out << "# Magnetic Field of Solenoid#           \n";
  out << "MagneticField: " <<  magnetic_field << "\n";
  out << "#                                       \n";
  out << "######                                  \n";
}

void DB_RunSummary::CreateAnalysisConf()
{

  // 
  // Creating file name
  //
  char number[4];
  sprintf(number,"%4.4d",run_number);
  //
  TString dname = "./param_run";
  dname += number;
  //
  mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
  if ( mkdir(dname.Data(), mode) != 0 )
   {
     //printf("directry %s exist!!!\n", dname.Data());
    }
  TString fname = dname;
  fname += "/analyzer.conf";

  std::fstream out(fname.Data(),std::ios::out);

  out << "### Type                                           \n";
  out << "### CDS1 : 0, BL1 : 5                              \n";
  out << "DataType:      5                                   \n";
  out << "###                                                \n";
  out << "StartBlockEvNum: 0                                 \n";
  out << "StopBlockEvNum:  0                                 \n";
  out << "StartEvNum: 0                                      \n";
  out << "StopEvNum:  -5000                                  \n";
  out << "SkipEvNum:  0                                      \n";
  out << "SkipEvNum:  0                                      \n";
  out << "###                                                \n";
  out << "CounterMap:       E15_Calib_CounterMap.param       \n";
  out << "BLDCWireMap:      E15_Calib_BLDCWireMap.param      \n";
  out << "CDCGeom:          E15_Calib_CDCGeom.param          \n";
  out << "CDCASDMap:        E15_Calib_CDCASDMap.param        \n";
  out << "CDCChannelMap:    E15_Calib_CDCChannelMap.param    \n";
  out << "GainMapCDS:       E15_Calib_GainMapCDS.param       \n";
  out << "GainMapBL:        E15_Calib_GainMapBL.param        \n"; 
  out << "GateMapBL:        E15_Calib_GateMapBL.param        \n";
  out << "DCTimeCorrBLDC:   E15_Calib_DCTimeCorrBLDC.root    \n";
  out << "XTMapCDS:         E15_Calib_XTMapCDS.param         \n";
  out << "XTMapBLROOT:      E15_Calib_XTMapBLROOT.root       \n";
  out << "GeomCDS:          E15_Calib_GeomCDS.param          \n";
  out << "GeomBL:           E15_Calib_GeomBL.param           \n";
  out << "SlewingMapCDS:    E15_Calib_SlewingMapCDS.param    \n";
  out << "SlewingMapBL:     E15_Calib_SlewingMapBL.param     \n";
  out << "ScalerMap:        E15_Calib_ScalerMap.param        \n";
  out << "CDSFittingParam:  E15_Calib_CDSFittingParam.param  \n";
  out << "BLDCFittingParam: E15_Calib_BLDCFittingParam.param \n";
  out << "TransferMatrix:   E15_Calib_TransferMatrix.param   \n";
  out << "ReslMap:          E15_Calib_ReslMap.param          \n";

}


int DB_Manager::tranceport(char *a, char *b)
{
  int doublequote = 0;
  int ncommma     = 0;
  while((*b++ = *a)) {
     char tmp = *a++;
     if(tmp == '"') {
       doublequote++;
     }
     if      (tmp == ',' && (doublequote%2)==0) {*b++ = ' '; ncommma++;}
     else if (tmp == ',' && (doublequote%2)==1) {*b-- = ';';}
  }
  return ncommma;
}

