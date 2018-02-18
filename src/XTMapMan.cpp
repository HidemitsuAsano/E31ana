// XTMapMan.cpp

#include <string>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <new>

#include "XTMapMan.h"
#include "GlobalVariables.h"
#include "TFile.h"
#include "BLDCWireMapMan.h"

ClassImp(XTMap);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
XTMap::XTMap()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double XTMap::GetParam( int i )
{
  if( 0<=i && i<(int)Par.size() )
    return Par[i];
  else
    return -999.;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMap::Clear()
{
  Par.clear();
}

ClassImp(XTMapMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
XTMapMan::XTMapMan():verbosity(0)
{
  FileNameCDS = DefaultFileName;
  FileNameBL  = DefaultFileName;
  ROOTNameBL  = DefaultFileName;
  if(verbosity){
    std::cout << "verbosity level of " << __FILE__ << " : " << verbosity << std::endl;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
XTMapMan::XTMapMan( const XTMapMan &right ):verbosity(0)
{
  FileNameCDS = right.FileNameCDS;
  FileNameBL  = right.FileNameBL;
  ROOTNameBL  = right.ROOTNameBL;
  for( XTMapContainer::const_iterator i=right.xtContainer.begin();
       i!=right.xtContainer.end(); i++ ){
    xtContainer[i->first] = i->second;
  }
  for( XTMapContainer2::const_iterator i=right.xtContainer2.begin();
       i!=right.xtContainer2.end(); i++ ){
    xtContainer2[i->first] = i->second;
  }
  for( XTHistContainer::const_iterator i=right.xthistContainer.begin();
       i!=right.xthistContainer.end(); i++ ){
    xthistContainer[i->first] = i->second;
  }
  if(verbosity){
    std::cout << "verbosity level of " << __FILE__ << " : " << verbosity << std::endl;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
XTMapMan::~XTMapMan()
{
  xtContainer.clear();
  xtContainer2.clear();
  xthistContainer.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::SetFileNameCDS( const std::string & filename )
{
  FileNameCDS = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::SetFileNameBL( const std::string & filename )
{
  FileNameBL = filename;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::SetROOTNameBL( const std::string & filename )
{
  ROOTNameBL = filename;
}

const unsigned int KEYMASK  = 0x000F;
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

const int MAXCHAR = 512;

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool XTMapMan::Initialize(BLDCWireMapMan *wiremapman)
{
 
  static const std::string funcname = "XTMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ..."/* << std::endl*/;

  int cid;
  int layer,wire,npar;
  double par[10];
  double cellsize;
  int nd;
  unsigned int key;
  char str[MAXCHAR];
  FILE *fp;

  xtContainer.clear();
  xtContainer2.clear();

  if( FileNameCDS!=DefaultFileName ){
    if( (fp=fopen(FileNameCDS.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameCDS << "]" << std::endl;
      exit(-1);
    }
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      
      if( (nd=sscanf(str,"%d %d %d %d %lf %lf %lf %lf %lf %lf", 
		     &cid, &layer, &wire, &npar, 
		     &par[0], &par[1], &par[2], 
		     &par[3], &par[4], &par[5] )) >=6 ) {
#if 0
	std::cout << cid << layer << "  " << wire << "  " << npar;
	for( int i=0; i<npar; i++ ) std::cout << "  " << par[i];
	std::cout << std::endl;
#endif
	key = KEY(cid,layer,wire);
	XTMap *xtmap = new XTMap();
	for( int i=0; i<npar; i++ ) xtmap->SetParam(par[i]);
	xtContainer[key] = *xtmap;
	delete xtmap;
      }
      else{
	std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	std::cerr << std::string(str) << std::endl;
      }
    }
    fclose(fp);
    //    std::cout << "[" << funcname << "] CDS is initialized." << std::endl;
    std::cout <<" CDS." ;//<< std::endl;
  }

  std::vector<double> threshold;
  bool OLDPARAM=true;
  const int MAXCHAR = 512;
  const int MAXTOKEN = 30;
  const char* DELIMITER = " ";
  if( FileNameBL!=DefaultFileName ){
    if( (fp=fopen(FileNameBL.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameBL << "]" << std::endl;
      exit(-1);
    }
    
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      if( OLDPARAM&&(nd=sscanf(str,"%d %d %d %d %lf %lf %lf %lf %lf %lf", 
		     &cid, &layer, &wire, &npar, 
		     &par[0], &par[1], &par[2], 
		     &par[3], &par[4], &par[5] )) >=6 ) {
#if 0
	std::cout << cid << layer << "  " << wire << "  " << npar;
	for( int i=0; i<npar; i++ ) std::cout << "  " << par[i];
	std::cout << std::endl;
#endif
	key = KEY(cid,layer,wire);
	XTMap *xtmap = new XTMap();
	for( int i=0; i<npar; i++ ) xtmap->SetParam(par[i]);
	xtContainer[key] = *xtmap;
	delete xtmap;
      }else if( (nd=sscanf(str,"XTParam: %d %lf %d ",&cid,&cellsize,&npar))==3 ){
	//	std::cout<<cid<<"  "<<npar<<"  "<<cellsize<<std::endl;
	OLDPARAM=false;	
      }else if(!OLDPARAM){
	const char* token[MAXTOKEN] = {};
	int n=0;
	token[0] = strtok(str, DELIMITER);
	if (token[0] == NULL) continue;
	if( !strcmp(token[0], "#") ) continue;
	for (n = 1; n < MAXTOKEN; n++){
	  token[n] = strtok( NULL, DELIMITER);
	  if (token[n] == NULL ) break;
	}
	//	std::cout<<token[0]<<"  "<<n<<std::endl;
	if( !strcmp(token[0], "Threshold:") ){
	  threshold.clear();
	  for(int i=0;i<npar;i++){
	    double tmpth=atof(token[i+1])*cellsize/100.;
	    threshold.push_back(tmpth);
	  }
	}else if(n==npar+3){
	  if(atoi(token[0])!=cid){
	    std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	    std::cerr << std::string(str) << std::endl;
	    exit(-1);
	  }
	  layer=atoi(token[1]);
	  wire=atoi(token[2]);
	  key = KEY(cid,layer,wire);
	  TGraph* gr=new TGraph(npar);
	  for (int i = 0; i < npar; i++){
	    //	    std::cout<<cid<<"  "<<layer<<"  "<<wire<<"  "<<atof(token[2+i])<<"  "<<threshold[i]<<std::endl;
	    gr->SetPoint(i,atof(token[3+i]),threshold[i]);
	  }

	  xtContainer2[key] = *gr;
	  gr->SetName(Form("xt_%d_%d_%d",cid,layer,wire));
	  gROOT->Append(gr);
	}
      }
      else{
	std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	std::cerr << std::string(str) << std::endl;
      }
    }
    fclose(fp);
    std::cout <</* "[" << funcname << "]*/" BeamLine.";// << std::endl;
  }
  
  if( ROOTNameBL!=DefaultFileName ){
    TFile *f=new TFile(ROOTNameBL.c_str());
    if(!f->IsOpen()){
      std::cout<<"!!!!!\t"<<ROOTNameBL<<" does not exist !!!!!"<<std::endl; 
      ROOTNameBL=DefaultFileName;
      return true;
    }
    TH1F* h1;

    for(int lay=1;lay<=8;lay++){
      int nwireblc1a = wiremapman->GetNWire(CID_BLC1a,lay );
      if(verbosity>0){
        std::cout << __FILE__ << " L." << __LINE__ << std::endl;
        std::cout << "CID of BLC1a " << CID_BLC1a << std::endl;
        std::cout << "Number of wires of BLC1a: " << nwireblc1a << " in layer " << lay << std::endl;
      }
      for(int wire=1;wire<=nwireblc1a;wire++){
	h1=(TH1F*)f->Get(Form("xt_PDC1_%d_%d",lay,wire));
	if(!h1)  h1=(TH1F*)f->Get(Form("xt_BLC1a_%d_%d",lay,wire));
	key = KEY(CID_BLC1a,lay,wire);
	xthistContainer[key] = *h1;
	xthistContainer[key].SetDirectory(gROOT);
	//	delete h1
      }
      int nwireblc1b = wiremapman->GetNWire(CID_BLC1b,lay );
      if(verbosity>0){
        std::cout << __FILE__ << " L." << __LINE__ << std::endl;
        std::cout << "CID of BLC1b " << CID_BLC1b << std::endl;
        std::cout << "Number of wires of BLC1b: " << nwireblc1b << " in layer " << lay << std::endl;
      }
      for(int wire=1;wire<=nwireblc1b;wire++){
	h1=(TH1F*)f->Get(Form("xt_PDC2_%d_%d",lay,wire));
	if(!h1)  h1=(TH1F*)f->Get(Form("xt_BLC1b_%d_%d",lay,wire));
	key = KEY(CID_BLC1b,lay,wire);
	xthistContainer[key] = *h1;
	xthistContainer[key].SetDirectory(gROOT);
	//	delete h1;
      }
      int nwireblc2a = wiremapman->GetNWire(CID_BLC2a,lay );
      if(verbosity>0){
        std::cout << __FILE__ << " L." << __LINE__ << std::endl;
        std::cout << "CID of BLC2a " << CID_BLC2a << std::endl;
        std::cout << "Number of wires of BLC2a: " << nwireblc2a << " in layer " << lay << std::endl;
      }
      for(int wire=1;wire<=nwireblc2a;wire++){
	h1=(TH1F*)f->Get(Form("xt_BLC1_%d_%d",lay,wire));
	if(!h1)  h1=(TH1F*)f->Get(Form("xt_BLC2a_%d_%d",lay,wire));
	key = KEY(CID_BLC2a,lay,wire);
	xthistContainer[key] = *h1;
	xthistContainer[key].SetDirectory(gROOT);
	//	delete h1;
      }
      int nwireblc2b = wiremapman->GetNWire(CID_BLC2b,lay );
      if(verbosity>0){
        std::cout << __FILE__ << " L." << __LINE__ << std::endl;
        std::cout << "CID of BLC2b " << CID_BLC2b << std::endl;
        std::cout << "Number of wires of BLC2b: " << nwireblc2b << " in layer " << lay << std::endl;
      }
      for(int wire=1;wire<=nwireblc2b;wire++){
	h1=(TH1F*)f->Get(Form("xt_BLC2_%d_%d",lay,wire));
	if(!h1)  h1=(TH1F*)f->Get(Form("xt_BLC2b_%d_%d",lay,wire));
	key = KEY(CID_BLC2b,lay,wire);
	xthistContainer[key] = *h1;
	xthistContainer[key].SetDirectory(gROOT);
	//	delete h1;
      }
      //  for(int wire=1;wire<=15;wire++){
      int nwirebpc = wiremapman->GetNWire(CID_BPC,lay );
      if(verbosity>0){
        std::cout << __FILE__ << " L." << __LINE__ << std::endl;
        std::cout << "CID of BPC " << CID_BPC << std::endl;
        std::cout << "Number of wires of BPC: " << nwirebpc << " in layer " << lay << std::endl;
      }
      for(int wire=1;wire<=nwirebpc;wire++){
	h1=(TH1F*)f->Get(Form("xt_BPC_%d_%d",lay,wire));
	key = KEY(CID_BPC,lay,wire);
	xthistContainer[key] = *h1;
	xthistContainer[key].SetDirectory(gROOT);
	//	delete h1;
      }
      int nwirefdc1 = wiremapman->GetNWire(CID_FDC1,lay );
      if(verbosity>0){
        std::cout << __FILE__ << " L." << __LINE__ << std::endl;
        std::cout << "CID of FDC1 " << CID_FDC1 << std::endl;
        std::cout << "Number of wires of FDC1: " << nwirefdc1 << " in layer " << lay << std::endl;
      }
      for(int wire=1;wire<=nwirefdc1;wire++){
	h1=(TH1F*)f->Get(Form("xt_FDC1_%d_%d",lay,wire));
	if(h1){
	  key = KEY(CID_FDC1,lay,wire);
	  xthistContainer[key] = *h1;
	  xthistContainer[key].SetDirectory(gROOT);
	  //	delete h1;
	}
      }
    }
    std::cout << "BLhist.";
    f->Close();
    delete f;
  }
  std::cout <</* "[" << funcname << "] Initialization*/" finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool XTMapMan::SetParam( const int &cid, const int &layer, const int &wire, const int &npar, double *par )
{
  static const std::string funcname = "XTMapMan::SetParam";
  unsigned int key;

  key = KEY(cid,layer,wire);
  XTMapContainer::iterator ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){ // replace parameters
    (ix->second).Clear();
    for( int i=0; i<npar; i++ ) (ix->second).SetParam( par[i] );
    return true;
  }
  else{
    if( 0<wire ){ // add parameters
      XTMap *xtmap = new XTMap();
      for( int i=0; i<npar; i++ ) xtmap->SetParam(par[i]);
      xtContainer[key] = *xtmap;
      delete xtmap;
      return true;
    }
    else{
      std::cout << " Invalid value!!! [" << funcname << "]"
		<< " cid:" << cid << " layer:" << layer << " wire:" << wire
		<< std::endl;
      return false;
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double XTMapMan::CalcDriftLength( const int &cid, const int &layer, const int &wire, const double &dt )
{
  static const std::string funcname = "XTMapMan::CalcDriftLength";
  unsigned int key;
  double t = 1.;
  double ctime = 0.;
  XTMapContainer::iterator ix ;
  XTHistContainer::iterator ixh ;
  XTMapContainer2::iterator ixg ;
  
#if 0
  std::cout << cid << "  " << layer << "  " << wire << std::endl;
  key = KEY(cid,layer,wire);
  ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    for( int i=0; i<(ix->second).nParam(); i++ ){
      std::cout << (ix->second).GetParam(i) << "   ";
    }
    std::cout << std::endl;
  }
#endif
  
  key = KEY(cid,layer,wire);

  ixg = xtContainer2.find(key);
  if( ixg != xtContainer2.end() ){
    ctime = (ixg->second).Eval(dt);// in mm  
    //    std::cout<<dt<<"  "<<ctime<<std::endl;
    return ctime*0.1; //mm -> cm
  }
  
  ixh = xthistContainer.find(key);
  if( ixh != xthistContainer.end() ){
    // bin = (ixh->second).FindBin(dt);
    // ctime = (ixh->second).GetBinContent(bin);// in mm  
    ctime = (ixh->second).Interpolate(dt);// in mm   // 201505
    return ctime*0.1; //mm -> cm
  }
  
  // wire# != 0
  ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    for( int i=0; i<(ix->second).nParam(); i++ ){
      t = 1.;
      for( int j=0; j<i; j++ ){ t *= dt; }
      ctime += ((ix->second).GetParam(i))*t;
    }
    return ctime;
  }
  
  // wire# = 0
  key = KEY(cid,layer,0);
  ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    for( int i=0; i<(ix->second).nParam(); i++ ){
      t = 1.;
      for( int j=0; j<i; j++ ){ t *= dt; }
      ctime += ((ix->second).GetParam(i))*t;
    }
    return ctime;
  }
  
  // no param
  return dt;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double XTMapMan::CalcDxDt( const int &cid, const int &layer, const int &wire, const double &dt )
{
  static const std::string funcname = "XTMapMan::CalcDxDt";
  unsigned int key;
  double t = 1.;
  double ctime = 0.;
  int bin = 0;
  XTMapContainer::iterator ix ;
  XTHistContainer::iterator ixh ;
  
#if 0
  std::cout << cid << "  " << layer << "  " << wire << std::endl;
  key = KEY(cid,layer,wire);
  ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    for( int i=0; i<(ix->second).nParam(); i++ ){
      std::cout << (ix->second).GetParam(i) << "   ";
    }
    std::cout << std::endl;
  }
#endif
  key = KEY(cid,layer,wire);
  ixh = xthistContainer.find(key);
  if( ixh != xthistContainer.end() ){
    bin = (ixh->second).FindBin(dt);
    ctime = (ixh->second).GetBinContent(bin);// in mm  
    return ctime*0.1; //mm -> cm
  }
  
  // wire# != 0
  key = KEY(cid,layer,wire);
  ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    for( int i=0; i<(ix->second).nParam()-1; i++ ){
      t = 1.;
      for( int j=0; j<i; j++ ){ t *= dt; }
      ctime += (i+1)*((ix->second).GetParam(i+1))*t;
    }
    return ctime;
  }
  
  // wire# = 0
  key = KEY(cid,layer,0);
  ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    for( int i=0; i<(ix->second).nParam()-1; i++ ){
      t = 1.;
      for( int j=0; j<i; j++ ){ t *= dt; }
      ctime += (i+1)*((ix->second).GetParam(i+1))*t;
    }
    return ctime;
  }
  
  // no param
  return dt;
}

double XTMapMan::CalcDriftTime( const int &cid, const int &layer, const int &wire, const double &dl )
{
  // '10/5/20 k.t.
  // this method is quite stupid.
  // this method calculates dt only in the range in which the correration between dt and dl is the one valued function.
  static const std::string funcname = "XTMapMan::CalcDriftTime";

  if(cid!=CID_CDC) return 0;
  
  double lcon[10000],tcon[10000];
  double stepsize = 0.05; // nsec
  int num=0;
  double ll,llold=-999;
  for( int i=0; i<10000; i++ ){
    tcon[i] = stepsize*i;
    ll = CalcDriftLength( cid, layer, wire, tcon[i] );
    if( ll<llold ){
      break;
    }
    num++;
    lcon[i]=llold=ll;
  }

  int il=0;
  for( int i=0; i<num; i++ ){
    if( lcon[i]<dl ) il = i;
    else break;
  }
  if(9999<=il) il=9998;

  double dt = tcon[il] + (dl-lcon[il])*(tcon[il+1]-tcon[il])/(lcon[il+1]-lcon[il]);
#if 0  
  std::cout << " num:" << num << " il:" << il << std::endl;
  std::cout << " lcon[il]:" << lcon[il] << " dl:" << dl << " lcon[il+1]:" << lcon[il+1] << std::endl;
  std::cout << " dt:" << dt << std::endl;
#endif
  return dt;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int XTMapMan::nparam( const int &cid, const int &layer, const int &wire )
{
  static const std::string funcname = "XTMapMan::nparam";
  unsigned int key;

  // wire# != 0
  key = KEY(cid,layer,wire);
  XTMapContainer::iterator ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    return (ix->second).nParam();
  }

  // wire# = 0
  key = KEY(cid,layer,0);
  ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    return (ix->second).nParam();
  }
  return 0;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double XTMapMan::param( const int &cid, const int &layer, const int &wire, const int &i )
{
  static const std::string funcname = "XTMapMan::param";
  unsigned int key;

  // wire# != 0
  key = KEY(cid,layer,wire);
  XTMapContainer::iterator ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    return (ix->second).GetParam(i);
  }

  // wire# = 0
  key = KEY(cid,layer,0);
  ix = xtContainer.find(key);
  if( ix != xtContainer.end() ){
    return (ix->second).GetParam(i);
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::PrintMapHeader( std::ostream &p_out )
{
  static const std::string funcname = "XTMapMan::PrintMapHeader";

  std::cout << " ---- " << funcname << " ---- " << std::endl;
  p_out << "#" << std::endl
	<< "# XT curve for Drift chamers." << std::endl
	<< "#" << std::endl
	<< "# units : [nsec] --> [cm]" << std::endl
	<< "# At first, XT relations are represented by n-th polynomials." << std::endl
	<< "# dx = p0 + p1*dt + p2*dt**2 + p3*dt**3 + p4*dt**4 + p5*dt**5" << std::endl
	<< "#" << std::endl
	<< "# The wire#=0 means that the value is applied to all wires." << std::endl
	<< "# If wire#!=0 exists after wire#=0, the value is uploaded. " << std::endl
	<< "#" << std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::PrintMap( const int &Cid, std::ostream &p_out )
{
  static const std::string funcname = "XTMapMan::PrintMap";
  unsigned int key;
  XTMap xtmap;
  int cid, cid_old=-1;

  std::cout << " ---- " << funcname << " ---- " << std::endl;

  for( XTMapContainer::const_iterator i=xtContainer.begin();
       i!=xtContainer.end(); i++ ){
    key = i->first;
    xtmap = i->second;
    cid = ((key>>CSHIFT)&CMASK);
    if( !( cid==Cid || Cid==-1 ) ) continue;
    if( cid!=cid_old ){
      p_out	<< "# CID  Layer    Wire    npar            p0             p1             p2             p3             p4             p5" << std::endl;
      cid_old = cid;
    }
    p_out << std::setw(5) << cid
	  << std::setw(8) << ((key>>LSHIFT)&LMASK)
	  << std::setw(8) << ((key>>WSHIFT)&WMASK)
	  << std::setw(7) << xtmap.nParam();
    p_out.setf(std::ios::showpoint);
    p_out.setf(std::ios::scientific, std::ios::floatfield);
    for( int j=0; j<xtmap.nParam(); j++ ){
      p_out  << std::setprecision(5) << std::setw(15) << xtmap.GetParam(j);
    }
    p_out << std::endl;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::PrintMapBL( std::ostream &p_out )
{
  static const std::string funcname = "XTMapMan::PrintMapBL";
  PrintMapHeader(p_out);
  PrintMap(CID_BLC1a,p_out);
  PrintMap(CID_BLC1b,p_out);
  PrintMap(CID_BLC2a,p_out);
  PrintMap(CID_BLC2b,p_out);
  PrintMap(CID_BPC,p_out);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::PrintMapCDS( std::ostream &p_out )
{
  static const std::string funcname = "XTMapMan::PrintMapCDS";
  PrintMapHeader(p_out);
  PrintMap(CID_CDC,p_out);
}
