// ConfMan.cpp

#include <string>
#include <iostream>
#include <cstdio>
#include <new>

#include "ConfMan.h"

ClassImp(ConfMan)

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::ConfMan()
  : TObject()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::ConfMan( const std::string &filename )
  : TObject()
{
  ConfFileName = filename;
  DataType = Type_CDS1;
  StartEvNum = StopEvNum = 0;
  StartBlockEvNum = StopBlockEvNum = 0;
  OutFileName = "tmp.root";
  VMEFileName = DefaultFileName;
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
  SlewingMapFileNameCDS = DefaultFileName;
  SlewingMapFileNameBL = DefaultFileName;
  SlewingMapFileNameSDD = DefaultFileName;
  ReslMapFileName = DefaultFileName;
  CDSFittingParamFileName = DefaultFileName;
  BLDCFittingParamFileName = DefaultFileName;
  FADCParamFileName = DefaultFileName;
  TransferMatrixFileName = DefaultFileName;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::ConfMan( const std::string &filename,
		  const std::string &outfilename) 
  : TObject()
{
  ConfFileName = filename;
  DataType = Type_CDS1;
  StartEvNum = StopEvNum = 0;
  StartBlockEvNum = StopBlockEvNum = 0;
  OutFileName = outfilename;
  VMEFileName = DefaultFileName;
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
  SlewingMapFileNameCDS = DefaultFileName;
  SlewingMapFileNameBL = DefaultFileName;
  SlewingMapFileNameSDD = DefaultFileName;
  ReslMapFileName = DefaultFileName;
  CDSFittingParamFileName = DefaultFileName;
  BLDCFittingParamFileName = DefaultFileName;
  FADCParamFileName = DefaultFileName;
  TransferMatrixFileName = DefaultFileName;
}
ConfMan::ConfMan( const std::string &filename,
		  const std::string &outfilename, 
		  const std::string &vmefilename )
  : TObject()
{
  ConfFileName = filename;
  DataType = Type_CDS1;
  StartEvNum = StopEvNum = 0;
  StartBlockEvNum = StopBlockEvNum = 0;
  OutFileName = outfilename;
  VMEFileName = vmefilename;
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
  SlewingMapFileNameCDS = DefaultFileName;
  SlewingMapFileNameBL = DefaultFileName;
  SlewingMapFileNameSDD = DefaultFileName;
  ReslMapFileName = DefaultFileName;
  CDSFittingParamFileName = DefaultFileName;
  BLDCFittingParamFileName = DefaultFileName;
  FADCParamFileName = DefaultFileName;
  TransferMatrixFileName = DefaultFileName;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ConfMan::~ConfMan()
{
  End();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
static const int MAXCHAR = 144;
bool ConfMan::ReadConfFile(std::string conffilename)
{
  FILE *fp;
  char str[MAXCHAR], str1[MAXCHAR];
  int evnum;
  int typ;


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
  }
  
  fclose(fp);
  return true;
}
bool ConfMan::Initialize(bool PRINT)
{
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
    std::cout << "        GainMapSDD : " << GainMapFileNameSDD << std::endl;
    std::cout << "        GateMapCDS : " << GateMapFileNameCDS << std::endl;
    std::cout << "         GateMapBL : " << GateMapFileNameBL << std::endl;
    std::cout << "        GateMapSDD : " << GateMapFileNameSDD << std::endl;
    std::cout << "          XTMapCDS : " << XTMapFileNameCDS << std::endl;
    std::cout << "           XTMapBL : " << XTMapFileNameBL << std::endl;
    std::cout << "       XTMapBLROOT : " << XTMapROOTNameBL << std::endl;
    std::cout << "       GeometryCDS : " << GeometryMapFileNameCDS << std::endl;
    std::cout << "        GeometryBL : " << GeometryMapFileNameBL << std::endl;
    std::cout << "     SlewingMapCDS : " << SlewingMapFileNameCDS << std::endl;
    std::cout << "      SlewingMapBL : " << SlewingMapFileNameBL << std::endl;
    std::cout << "     SlewingMapSDD : " << SlewingMapFileNameSDD << std::endl;
    std::cout << "         ScalerMap : " << ScalerMapFileName << std::endl;
    std::cout << "           ReslMap : " << ReslMapFileName << std::endl;
    std::cout << "   CDSFittingParam : " << CDSFittingParamFileName << std::endl;
    std::cout << "  BLDCFittingParam : " << BLDCFittingParamFileName << std::endl;
    std::cout << "         FADCParam : " << FADCParamFileName << std::endl;
    std::cout << "    TransferMatrix : " << TransferMatrixFileName << std::endl;
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
    GeomMapManager->Initialize();
  }

  if( SlewingMapFileNameCDS != DefaultFileName || SlewingMapFileNameBL  != DefaultFileName || SlewingMapFileNameSDD  != DefaultFileName )
    SlewingMapManager = new SlewingMapMan();
  if( SlewingMapManager!=0 ){
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
//   delete &CounterMapManager;
//   delete &WireMapManager;
//   delete &GainMapManager;
//   delete &XTMapManager;
//   delete &CDHGeomMapManager;

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
    if (strstr(direntry,".o"))         continue;
    if (strstr(direntry,".d"))         continue;
    if (strstr(direntry,"Dict.h"))         continue;
    if (strstr(direntry,"Dict.cpp"))         continue;
    if (strstr(direntry,"#"))         continue;
    if (strstr(direntry,"~"))         continue;
    if (strstr(direntry,"Dict.cpp"))         continue;
    if (strstr(direntry,".cpp")    ||
	strstr(direntry,".h")    ||
	strstr(direntry,".c")) {
      WriteFile(afile, direntry);
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
  WriteFile(  ConfFileName  );
  WriteFile(  CounterMapFileName  );
  WriteFile(  BLDCWireMapFileName  );
  WriteFile(  CDCGeometryFileName  );
  WriteFile(  CDCASDMapFileName  );
  WriteFile(  CDCChannelMapFileName  );
  WriteFile(  GainMapFileNameCDS  );
  WriteFile(  GainMapFileNameBL  );
  WriteFile(  GainMapFileNameSDD  );
  WriteFile(  GateMapFileNameCDS  );
  WriteFile(  GateMapFileNameBL  );
  WriteFile(  GateMapFileNameSDD  );
  WriteFile(  XTMapFileNameCDS  );
  WriteFile(  XTMapFileNameBL  );
  WriteFile(  XTMapROOTNameBL  );
  WriteFile(  GeometryMapFileNameCDS  );
  WriteFile(  GeometryMapFileNameBL  );
  WriteFile(  SlewingMapFileNameCDS  );
  WriteFile(  SlewingMapFileNameBL  );
  WriteFile(  SlewingMapFileNameSDD  );
  WriteFile(  ReslMapFileName  );
  WriteFile(  CDSFittingParamFileName  );
  WriteFile(  BLDCFittingParamFileName  );
  WriteFile(  FADCParamFileName  );
  WriteFile(  TransferMatrixFileName  );
  savdir->cd();
  return;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ConfMan::WriteFile(TString afile,TString direntry){
  if(afile==DefaultFileName) return;
  //  std::cout<<afile<<std::endl;
  TMacro *m = new TMacro(afile);
  if(direntry)
    m->Write(direntry);
  else
    m->Write(afile);
  delete m;
  return;
}
