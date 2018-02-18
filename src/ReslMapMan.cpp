// ReslMapMan.cpp

#include "ReslMapMan.h"

ClassImp(ReslMap);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ReslMap::ReslMap()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ReslMap::~ReslMap()
{
}

ClassImp(ReslMapMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ReslMapMan::ReslMapMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ReslMapMan::ReslMapMan( const std::string & filename )
{
  FileName = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ReslMapMan::~ReslMapMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ReslMapMan::SetFileName( const std::string & filename )
{
  FileName = filename;
}

//const unsigned int KEYMASK  = 0x000F;
// |0111|1111|0001|1111|0000|1111|1111|0011|
const unsigned int WMASK    = 0x00FF;      /* Wire Mask 6 Bits (0-255) */
const unsigned int LMASK    = 0x001F;      /* Layer Mask 5 Bits (0-31) */
const unsigned int CMASK    = 0x007F;      /* CID Mask 7 Bits (0-31) */
const int          WSHIFT   =  4;
const int          LSHIFT   = 16;
const int          CSHIFT   = 24;
const unsigned int KEYFLAG  = 0x0003;

#define KEY(cid,layer,wire) \
((((cid)&CMASK)<<CSHIFT) | (((layer)&LMASK)<<LSHIFT) | (((wire)&WMASK)<<WSHIFT) | KEYFLAG )

const int MAXCHAR = 144;

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ReslMapMan::Initialize()
{
  static const std::string funcname = "ReslMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ... " << std::endl;

  int cid;
  int layer, wire;
  double tresl,eresl;
  int nd;
  unsigned int key;
  char str[MAXCHAR];
  FILE *fp;

  if( (fp=fopen(FileName.c_str(),"r"))==0 ){
    std::cerr << " File open fail. [" << FileName << "]" << std::endl;
    exit(-1);
  }
  
  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;

    if( (nd=sscanf(str,"%d %d %d %lf %lf", &cid, &layer, &wire, &tresl, &eresl ))==5 ){
#if 0
      std::cout << cid << "  " << layer << "  " << wire << "  " << tresl << "  " << eresl << std::endl;
#endif
      key = KEY(cid,layer,wire);
      ReslMap *rmap = new ReslMap();
      rmap->SetTResl(tresl); rmap->SetEResl(eresl);
      reslContainer[key] = *rmap;
      delete rmap;
    }
    else{
      std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
      std::cerr << std::string(str) << std::endl;
    }
  }

  fclose(fp);

  std::cout << "[" << funcname << "] Initialization finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ReslMapMan::SetParam( const int &cid, const int &layer, const int &wire, const double &tresl, const double &eresl )
{
  static const std::string funcname = "ReslMapMan::SetParam";
  unsigned int key;

  key = KEY(cid,layer,wire);
  ReslContainer::iterator ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){ // replace parameter
    (ir->second).SetTResl(tresl);
    (ir->second).SetEResl(eresl);
    //(ir->second) = resl;
    return true;
  }
  else{ // add parameter
      ReslMap *rmap = new ReslMap();
      rmap->SetTResl(tresl); rmap->SetEResl(eresl);
      reslContainer[key] = *rmap;
      delete rmap;
    //reslContainer[key] = resl;
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ReslMapMan::GetResolution( const int &cid, const int &layer, const int &wire, double &tresl, double &eresl )
{
  static const std::string funcname = "ReslMapMan::CalcResolution";
  unsigned int key;
  ReslContainer::iterator ir;

  // wire# != 0
  key = KEY(cid,layer,wire);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = gRandom->Gaus(0,(ir->second).tresl() );
    eresl = gRandom->Gaus(0,(ir->second).eresl() );
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return;
  }

  // wire# = 0
  key = KEY(cid,layer,0);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = gRandom->Gaus(0,(ir->second).tresl() );
    eresl = gRandom->Gaus(0,(ir->second).eresl() );
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return;
  }

  // layer# = 0
  key = KEY(cid,0,0);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = gRandom->Gaus(0,(ir->second).tresl() );
    eresl = gRandom->Gaus(0,(ir->second).eresl() );
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return;
  }

  // no param
  tresl = eresl = 0;
  return;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ReslMapMan::GetParam( const int &cid, const int &layer, const int &wire, double &tresl, double &eresl )
{
  static const std::string funcname = "ReslMapMan::GetParam";
  unsigned int key;
  ReslContainer::iterator ir;
  //  std::cout<<cid<<"  "<<layer<<"  "<<wire<<std::endl;
  // wire# != 0
  key = KEY(cid,layer,wire);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = (ir->second).tresl();
    eresl = (ir->second).eresl();
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return true;
  }

  // wire# = 0
  key = KEY(cid,layer,0);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = (ir->second).tresl();
    eresl = (ir->second).eresl();
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return true;
  }

  // layer# = 0
  key = KEY(cid,0,0);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = (ir->second).tresl();
    eresl = (ir->second).eresl();
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return true;
  }

  // no param
  tresl = eresl = 0;
  return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ReslMapMan::PrintMap( const int &Cid, std::ostream &p_out )
{
  static const std::string funcname = "ReslMapMan::PrintMap";
  unsigned int key;
  ReslMap rmap;

  std::cout << " ---- " << funcname << " ---- " << std::endl;

  for( ReslContainer::const_iterator ir=reslContainer.begin();
       ir!=reslContainer.end(); ir++ ){
    key = ir->first;
    rmap = ir->second;
    p_out << std::setw(5) << ((key<<CSHIFT)&CMASK)
	  << std::setw(8) << ((key<<LSHIFT)&LMASK)
	  << std::setw(8) << ((key<<WSHIFT)&WMASK)
	  << std::setw(10) << rmap.tresl()
	  << std::setw(10) << rmap.eresl()
	  << std::endl;
  }
}
