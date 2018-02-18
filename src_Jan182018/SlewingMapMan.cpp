// SlewingMapMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>
#include <cmath>

#include "SlewingMapMan.h"
#include "GlobalVariables.h"

ClassImp(SlewingMap);
ClassImp(SlewingMapMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SlewingMap::SlewingMap()
{
  //for( int i=0; i<MaxParam; i++ )
  //Par[i] = 0.0;
  //nPar = 0;
  Par.clear();
  Type = 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMap::SetParam( const int &npar, double *par )
{
  //if( MaxParam < npar ){
  //std::cout << " !!! #parameters is too large !!! " << std::endl;
  //return;
  //}

  //nPar = npar;
  Par.clear();
  for( int i=0; i<npar; i++ ){
    //Par[i] = par[i];
    Par.push_back(par[i]);
  }
  return;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMap::SetParam( const int &npar, const std::vector <double> &par )
{
  Par.clear();
  Par = par;
  return;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SlewingMapMan::SlewingMapMan()
{
  FileNameCDS = DefaultFileName;
  FileNameBL = DefaultFileName;
  FileNameSDD = DefaultFileName;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SlewingMapMan::SlewingMapMan( const SlewingMapMan &right )
{
  FileNameCDS = right.FileNameCDS;
  FileNameBL = right.FileNameBL;
  FileNameSDD = right.FileNameSDD;
  for( SlewingMapContainer::const_iterator i=right.slewingContainer.begin();
       i!=right.slewingContainer.end(); i++ ){
    slewingContainer[i->first] = i->second;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SlewingMapMan::~SlewingMapMan()
{
  slewingContainer.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMapMan::SetFileNameCDS( const std::string & filename )
{
  FileNameCDS = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMapMan::SetFileNameBL( const std::string & filename )
{
  FileNameBL = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMapMan::SetFileNameSDD( const std::string & filename )
{
  FileNameSDD = filename;
}

const unsigned int KEYMASK  = 0x000F;
// |1111|1110|1111|1011|1111|1110|1111|0011|
const unsigned int IMASK    = 0x000F;      /* I Mask 4 Bits (0-15) */
const unsigned int SMASK    = 0x01FF;      /* S Mask 9 Bits (0-511) */
const unsigned int UMASK    = 0x001F;      /* U Mask 5 Bits (0-31) */
const unsigned int CMASK    = 0x007F;      /* C Mask 7 Bits (0-127) */
const int          ISHIFT   =  4;
const int          SSHIFT   =  9;
const int          USHIFT   = 19;
const int          CSHIFT   = 25;
const unsigned int KEYFLAG  = 0x0003;

#define KEY(cid,seg,ud,i) \
((((cid)&CMASK)<<CSHIFT) | (((seg)&SMASK)<<SSHIFT) | (((ud)&UMASK)<<USHIFT) | (((i)&IMASK)<<ISHIFT) | KEYFLAG )

const int MAXCHAR = 144;

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SlewingMapMan::Initialize()
{
  static const std::string funcname = "SlewingMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";

  int cid,seg,ud,ith;
  //double par0,par1;
  int type;
  int npar;
  double par[MaxParam];
  int nd;
  unsigned int key;
  char str[MAXCHAR];
  FILE *fp;

  slewingContainer.clear();

  if( FileNameCDS!=DefaultFileName ){
    if( (fp=fopen(FileNameCDS.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameCDS << "]" << std::endl;
      exit(-1);
    }
    
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      
      //if( (nd=sscanf(str,"%d %d %d %d %lf %lf", &cid, &seg, &ud, &ith, &par0, &par1)) == 6 ) {
      if( (nd=sscanf(str,"%d %d %d %d %d %d %lf %lf %lf %lf", &cid, &seg, &ud, &type, &ith, &npar, &par[0], &par[1], &par[2], &par[3])) >= 7 ) {
#if 0
	std::cout << cid << "  " << seg << "  " << ud << "  " << type << "  " << ith << "  " << npar << "  ";
	for( int i=0; i<npar; i++ )
	  std::cout << par[i] << "  ";
	std::cout << std::endl;
#endif
	key = KEY(cid,seg,ud,ith);
	SlewingMap *aslewing = new SlewingMap();
	//aslewing->SetParam(par0,par1);
	aslewing->SetParam( npar, par );
	aslewing->SetType( type );
	slewingContainer[key] = *aslewing;
	delete aslewing;
      }
      else{
	std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	std::cerr << std::string(str) << std::endl;
      }
    }
    fclose(fp);
    //    std::cout << "[" << funcname << "] CDS is initialized." << std::endl;
    std::cout <<" CDS.";
  }

  if( FileNameBL!=DefaultFileName ){
    if( (fp=fopen(FileNameBL.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameBL << "]" << std::endl;
      exit(-1);
    }
    
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      
      //if( (nd=sscanf(str,"%d %d %d %d %lf %lf", &cid, &seg, &ud, &ith, &par0, &par1)) == 6 ) {
      if( (nd=sscanf(str,"%d %d %d %d %d %d %lf %lf %lf %lf", &cid, &seg, &ud, &type, &ith, &npar, &par[0], &par[1], &par[2], &par[3])) >= 7 ) {
#if 0
	std::cout << cid << "  " << seg << "  " << ud << "  " << type << "  " << ith << "  " << npar << "  ";
	for( int i=0; i<npar; i++ )
	  std::cout << par[i] << "  ";
	std::cout << std::endl;
#endif
	key = KEY(cid,seg,ud,ith);
	SlewingMap *aslewing = new SlewingMap();
	//aslewing->SetParam(par0,par1);
	aslewing->SetParam( npar, par );
	aslewing->SetType( type );
	slewingContainer[key] = *aslewing;
	delete aslewing;
      }
      else{
	std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	std::cerr << std::string(str) << std::endl;
      }
    }
    fclose(fp);
    //    std::cout << "[" << funcname << "] BeamLine is initialized." << std::endl;
    std::cout << " BeamLine.";
  }
  if( FileNameSDD!=DefaultFileName ){
    if( (fp=fopen(FileNameSDD.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameSDD << "]" << std::endl;
      exit(-1);
    }
    
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      
      //if( (nd=sscanf(str,"%d %d %d %d %lf %lf", &cid, &seg, &ud, &ith, &par0, &par1)) == 6 ) {
      if( (nd=sscanf(str,"%d %d %d %d %d %d %lf %lf %lf %lf", &cid, &seg, &ud, &type, &ith, &npar, &par[0], &par[1], &par[2], &par[3])) >= 7 ) {
#if 0
	std::cout << cid << "  " << seg << "  " << ud << "  " << type << "  " << ith << "  " << npar << "  ";
	for( int i=0; i<npar; i++ )
	  std::cout << par[i] << "  ";
	std::cout << std::endl;
#endif
	key = KEY(cid,seg,ud,ith);
	SlewingMap *aslewing = new SlewingMap();
	//aslewing->SetParam(par0,par1);
	aslewing->SetParam( npar, par );
	aslewing->SetType( type );
	slewingContainer[key] = *aslewing;
	delete aslewing;
      }
      else{
	std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	std::cerr << std::string(str) << std::endl;
      }
    }
    fclose(fp);
    std::cout << "[" << funcname << "] SDD is initialized." << std::endl;
  }
  //PrintMap();

  //  std::cout << "[" << funcname << "] Initialization finish." << std::endl;
  std::cout << " finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SlewingMapMan::SetParam( const int &cid, const int &seg, const int &ud, const int &ith, const int &type,
			      const int &npar, double *par )
{
  static const std::string funcname = "SlewingMapMan::SetParam";
  unsigned int key;

  key = KEY(cid,seg,ud,ith);
  SlewingMapContainer::iterator is = slewingContainer.find(key);
  if( is != slewingContainer.end() && npar<=MaxParam ){ // replace parameters
    (is->second).SetType(type);
    (is->second).SetParam(npar,par);
    return true;
  }
  else{
    key = KEY(cid,seg,ud,ith-1);
    is = slewingContainer.find(key);
    if( 0<ith && is!=slewingContainer.end() && npar<=MaxParam ){ // add parameters
      SlewingMap *aslewing = new SlewingMap();
      aslewing->SetParam( npar, par );
      aslewing->SetType( type );
      key = KEY(cid,seg,ud,ith);
      slewingContainer[key] = *aslewing;
      delete aslewing;
      return true;
    }
    else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid << " seg:" << seg << " ud:" << ud << " ith:" << ith << " type:" << type
	      << std::endl;
    return false;
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SlewingMapMan::SetParam( const int &cid, const int &seg, const int &ud, const int &ith, const int &type,
			      const int &npar, const std::vector <double> &par )
{
  static const std::string funcname = "SlewingMapMan::SetParam";
  unsigned int key;

  key = KEY(cid,seg,ud,ith);
  SlewingMapContainer::iterator is = slewingContainer.find(key);
  if( is != slewingContainer.end() && npar<=MaxParam ){ // replace parameters
    (is->second).SetType(type);
    (is->second).SetParam(npar,par);
    return true;
  }
  else{
    key = KEY(cid,seg,ud,ith-1);
    is = slewingContainer.find(key);
    if( 0<ith && is!=slewingContainer.end() && npar<=MaxParam ){ // add parameters
      SlewingMap *aslewing = new SlewingMap();
      aslewing->SetParam( npar, par );
      aslewing->SetType( type );
      key = KEY(cid,seg,ud,ith);
      slewingContainer[key] = *aslewing;
      delete aslewing;
      return true;
    }
    else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid << " seg:" << seg << " ud:" << ud << " ith:" << ith << " type:" << type
	      << std::endl;
    return false;
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SlewingMapMan::GetParam( const int &cid, const int &seg, const int &ud, const int &ith, int &type,
			      int &npar, std::vector <double> &par )
{
  static const std::string funcname = "SlewingMapMan::GetParam";
  unsigned int key;
  bool status;
  key = KEY(cid,seg,ud,ith);
  SlewingMapContainer::iterator is = slewingContainer.find(key);
  if( is != slewingContainer.end() ){
    par = (is->second).GetPar();
    npar = par.size();
    status = true;
  }
  else{
    status = false;
  }
  return status;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SlewingMapMan::GetParam( const int &cid, const int &seg, const int &ud, const int &ith, int &type,
			      int &npar, double *par )
{
  static const std::string funcname = "SlewingMapMan::GetParam";
  unsigned int key;
  bool status;
  key = KEY(cid,seg,ud,ith);
  SlewingMapContainer::iterator is = slewingContainer.find(key);
  if( is != slewingContainer.end() ){
    npar = (is->second).GetNPar();
    for(int i=0;i<npar;i++)
      par[i] = (is->second).GetPar(i);
    status = true;
  }
  else{
    status = false;
  }
  return status;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double SlewingMapMan::CalcCValue( const int &cid, const int &seg, const int &ud,
				  const double &time, const double &de )
{
  static const std::string funcname = "SlewingMapMan::CalcCValue";

  double ctime = time;

  if( fabs(de)< 0.00000001 ){
#if 0
    std::cout << " too small de[" << de << "]" << std::endl;
#endif
    return time;
  }

  unsigned int key;
  
  int type;
  int ith = 0;
  key = KEY(cid,seg,ud,ith);
  SlewingMapContainer::iterator is;
  double p0,p1,p2;
  while( (is=(slewingContainer.find(key))) != slewingContainer.end() ){
    type = (is->second).GetType();
    switch( type ){
    case 1:
      {
	p0 = (is->second).GetPar(0);
	p1 = (is->second).GetPar(1);
	ctime = ctime + p0/sqrt(de) + p1;
	ith++;
	key = KEY(cid,seg,ud,ith);
	break;
      }
    case 2:
      {
	p0 = (is->second).GetPar(0);
	p1 = (is->second).GetPar(1);
	p2 = (is->second).GetPar(2);
	ctime = ctime + p0/sqrt(de) + p1+p2*de;
	ith++;
	key = KEY(cid,seg,ud,ith);
	break;
      }
    case 3:
      {
	p0 = (is->second).GetPar(0);
	p1 = (is->second).GetPar(1);
	//	p2 = (is->second).GetPar(2);
	ctime = ctime + p0/sqrt(de-p1);
	ith++;
	key = KEY(cid,seg,ud,ith);
	break;
      }
    default:
      break;
    }
  }
  if( ith==0 ){
    std::cout << " cannot find parameters "
	      << " cid:" << cid << " seg:" << seg << " ud:" << ud << " time:" << time << " de:" << de
	      << std::endl;
  }
  
  return ctime;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double SlewingMapMan::CalcDATValue( const int &cid, const int &seg, const int &ud,
				    const double &ctime, const double &de )
{
  static const std::string funcname = "SlewingMapMan::CalcDATValue";
  double time = ctime;

  if( fabs(de)< 0.00000001 ){
#if 0
    std::cout << " too small de[" << de << "]" << std::endl;
#endif
    return time;
  }

  unsigned int key;
  
  int type;
  int ith_start = 0xf;
  int ith = ith_start;
  key = KEY(cid,seg,ud,ith);
  SlewingMapContainer::iterator is;
  while( 0<=ith ){
    if( (is=(slewingContainer.find(key))) != slewingContainer.end() ){
      type = (is->second).GetType();
      switch( type ){
      case 1:
	{
	  double p0 = (is->second).GetPar(0);
	  double p1 = (is->second).GetPar(1);
	  time = time - p0/sqrt(de) - p1;
	  break;
	}
      default:
	break;
      }
    }
    ith--;
    key = KEY(cid,seg,ud,ith);
  }
  if( ith==ith_start ){
    std::cout << " cannot find parameters "
	      << " cid:" << cid << " seg:" << seg << " ud:" << ud << " time:" << time << " de:" << de
	      << std::endl;
  }
  
  return time;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMapMan::PrintMapHeader( std::ostream &p_out )
{
  static const std::string funcname = "SlewingMapMan::PrintMapHeader";

  std::cout << " ---- " << funcname << " ---- " << std::endl;
  p_out<< "##" << std::endl
       << "#  Parameters for Slewing Correction" << std::endl
       << "#" << std::endl
       << "### type : 1" << std::endl
       << "#  ctime = time + par0/sqrt(dE) + par1" << std::endl
       << "#  par2 and par3 are not used" << std::endl
       << "###" << std::endl
       << "#" << std::endl
       << "###" << std::endl
       << "#  Different parameters sets can be set for same counter by numbering \"nth\" parameter in order to calculate iteratively." << std::endl
       << "#  The \"nth\" can be set up to 15" << std::endl
       << "###" << std::endl
       << "#" << std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMapMan::PrintMap( const int &Cid, std::ostream &p_out )
{
  static const std::string funcname = "SlewingMapMan::PrintMap";
  unsigned int key;
  SlewingMap aslewing;
  int cid, cid_old=-1;

  std::cout << " ---- " << funcname << " ---- " << std::endl;
  for( SlewingMapContainer::const_iterator i=slewingContainer.begin();
       i!=slewingContainer.end(); i++ ){
    key = i->first;
    aslewing = i->second;
    cid = ((key>>CSHIFT)&CMASK);
    if( !( cid==Cid || Cid==-1 ) ) continue;
    if( cid!=cid_old ){
      p_out << "#  CID  Seg  UD   correction-type  nth  nPar    Par0    Par1     Par2    Par3" << std::endl;
      cid_old = cid;
    }

    p_out << std::setw(5) << cid
	  << std::setw(6) << ((key>>SSHIFT)&SMASK)
	  << std::setw(4) << ((key>>USHIFT)&UMASK)
	  << std::setw(18) << aslewing.GetType()
	  << std::setw(5) << ((key>>ISHIFT)&IMASK)
	  << std::setw(6) << aslewing.GetNPar();
    p_out.setf(std::ios::showpoint);
    int i=0;
    for( i=0; i<3; i++ ){
      if(i<aslewing.GetNPar())
	p_out  << std::setprecision(5) << std::setw(12) << aslewing.GetPar(i);
      else
	p_out  << std::setprecision(5) << std::setw(12) << 0.;
    }
    p_out << std::endl;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMapMan::PrintMapBL( std::ostream &p_out )
{
  static const std::string funcname = "SlewingMapMan::PrintMapBL";
  PrintMapHeader(p_out);
  PrintMap(CID_BHD,p_out);
  PrintMap(CID_PA,p_out);
  PrintMap(CID_T0,p_out);
  PrintMap(CID_E0,p_out);
  PrintMap(CID_B1,p_out);
  PrintMap(CID_B2,p_out);
  PrintMap(CID_CVC,p_out);
  PrintMap(CID_NC,p_out);
  PrintMap(CID_CVC,p_out);
  PrintMap(CID_BVC,p_out);
  PrintMap(CID_LB,p_out);
  PrintMap(CID_PC,p_out);
  PrintMap(CID_BPD,p_out);
  PrintMap(CID_BD,p_out);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SlewingMapMan::PrintMapCDS( std::ostream &p_out )
{
  static const std::string funcname = "SlewingMapMan::PrintMapCDS";
  PrintMapHeader(p_out);
  PrintMap(CID_CDH,p_out);
  PrintMap(CID_IH,p_out);
}
