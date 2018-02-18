// DCTimeCorrMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>
#include <cmath>

#include "DCTimeCorrMan.h"
#include "GlobalVariables.h"

ClassImp(DCTimeCorr);
ClassImp(DCTimeCorrMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
DCTimeCorr::DCTimeCorr()
{
  Par.clear();
  Type = 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void DCTimeCorr::SetParam( const int &npar, double *par )
{
  Par.clear();
  for( int i=0; i<npar; i++ ){
    Par.push_back(par[i]);
  }
  return;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
DCTimeCorrMan::DCTimeCorrMan():FileNameCDC(DefaultFileName),FileNameBLDC(DefaultFileName)
{
  //FileNameCDC = DefaultFileName;
  //FileNameBLDC = DefaultFileName;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
DCTimeCorrMan::DCTimeCorrMan( const DCTimeCorrMan &right ):
FileNameCDC(right.FileNameCDC),
FileNameBLDC(right.FileNameBLDC)
{
  //FileNameCDC = right.FileNameCDC;
  //FileNameBLDC = right.FileNameBLDC;
  for( DCTimeCorrHistContainer::const_iterator i=right.dctimecorrContainer.begin();
       i!=right.dctimecorrContainer.end(); ++i ){
    dctimecorrContainer[i->first] = i->second;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
DCTimeCorrMan::~DCTimeCorrMan()
{
  dctimecorrContainer.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void DCTimeCorrMan::SetFileNameCDC( const std::string & filename )
{
  FileNameCDC = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void DCTimeCorrMan::SetFileNameBLDC( const std::string & filename )
{
  FileNameBLDC = filename;
}
//const unsigned int KEYMASK  = 0x000F;
// |0111|1111|0001|1111|0000|1111|1111|0011|
const unsigned int WMASK    = 0x00FF;      /* Wire Mask 8 Bits (0-255) */
const unsigned int LMASK    = 0x001F;      /* Layer Mask 5 Bits (0-31) */
const unsigned int CMASK    = 0x007F;      /* CID Mask 7 Bits (0-127) */
const int          WSHIFT   =  4;
const int          LSHIFT   = 16;
const int          CSHIFT   = 24;
const unsigned int KEYFLAG  = 0x0003;

#define KEY(cid,layer,wire) \
((((cid)&CMASK)<<CSHIFT) | (((layer)&LMASK)<<LSHIFT) | (((wire)&WMASK)<<WSHIFT) | KEYFLAG )
//const int MAXCHAR = 144;
//const int MaxParam = 4;
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool DCTimeCorrMan::Initialize()
{
  static const std::string funcname = "DCTimeCorrMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";

  unsigned int key;
  dctimecorrContainer.clear();
  
  TFile *f=new TFile(FileNameBLDC.c_str());
  if(!f->IsOpen()){
    std::cout<<"!!!!!\t"<<FileNameBLDC<<" does not exist !!!!!"<<std::endl; 
    FileNameBLDC=DefaultFileName;
    return true;
  }
  TGraph* gr;
  TString name[6]={"BLC1a","BLC1b","BLC2a","BLC2b","BPC","FDC1"};
  int cid[6]={CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_BPC,CID_FDC1};
  for(int ic=0;ic<6;ic++){
    for(int lay=1;lay<=8;lay++){      
      gr=(TGraph*)f->Get(Form("DCTimeCorr_%s_%d",name[ic].Data(),lay));
      if(!gr) continue;
      //      std::cout<<Form("DCTimeCorr_%s_%d",name[ic].Data(),lay)<<std::endl;
      key = KEY(cid[ic],lay,0);
      dctimecorrContainer[key] = *gr;
      //      dctimecorrContainer[key].SetDirectory(gROOT);
    }
  }
  f->Close();
  delete f;
  std::cout << " finish." << std::endl;
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double DCTimeCorrMan::CalcCValue( const int &cid, const int &layer, const int &wire,
				  const double &timemean, const double &timesub )
{
  static const std::string funcname = "DCTimeCorrMan::CalcCValue";
 
  double ctime = timemean;  

  unsigned int key;
  
  key = KEY(cid,layer,wire);
  DCTimeCorrHistContainer::iterator is;
  if( (is=(dctimecorrContainer.find(key))) != dctimecorrContainer.end() ){
    ctime = timemean - ((is->second).Eval(timesub));
    // if(abs(timesub)>100){
    //    std::cout<<cid<<"\t"<<(is->second).GetName()<<std::endl;
    //    std::cout<<"sub,mean,corr,ctime\t"<<timesub<<"\t"<<timemean<<"\t"<<(is->second).Eval(timesub)<<"\t"<<ctime<<std::endl;
    // }
  }else if(cid==CID_FDC1){
    ctime=timemean;
  }else{
    std::cout << " cannot find parameters "
	      << " cid:" << cid << " layer:" << layer << " wire:"<< wire <<" timemean:" << timemean << " timesub:" << timesub
	      << std::endl;
  }  
  return ctime;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double DCTimeCorrMan::CalcDATValue( const int &cid, const int &seg, const int &ud,
				    const double &ctime, const double &de )
{
  static const std::string funcname = "DCTimeCorrMan::CalcDATValue";
  double time = ctime;
  return time;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void DCTimeCorrMan::PrintMapHeader( std::ostream &p_out )
{
  static const std::string funcname = "DCTimeCorrMan::PrintMapHeader";

  std::cout << " ---- " << funcname << " ---- " << std::endl;
  p_out<< "##" << std::endl
       << "#  Parameters for DCTime Correction" << std::endl
       << "#" << std::endl
       << "### type : 1" << std::endl
       << "#  ctime = timemean - pow(timesub,2) * par0" << std::endl
       << "###" << std::endl
       << "#" << std::endl
       << "###" << std::endl
       << "###" << std::endl
       << "#" << std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void DCTimeCorrMan::PrintMap( const int &Cid, std::ostream &p_out )
{
  // static const std::string funcname = "DCTimeCorrMan::PrintMap";
  // unsigned int key;
  // DCTimeCorr amap;
  // int cid, cid_old=-1;

  // std::cout << " ---- " << funcname << " ---- " << std::endl;
  // for( DCTimeCorrContainer::const_iterator i=dctimecorrContainer.begin();
  //      i!=dctimecorrContainer.end(); i++ ){
  //   key = i->first;
  //   amap = i->second;
  //   cid = ((key>>CSHIFT)&CMASK);
  //   if( !( cid==Cid || Cid==-1 ) ) continue;
  //   if( cid!=cid_old ){
  //     p_out << "#  CID  Layer  Wire   Type   nPar    Par0    Par1" << std::endl;
  //     cid_old = cid;
  //   }

  //   p_out << std::setw(5) << cid
  // 	  << std::setw(6) << ((key>>LSHIFT)&LMASK)
  // 	  << std::setw(4) << ((key>>WSHIFT)&WMASK)
  // 	  << std::setw(18) << amap.GetType()
  // 	  << std::setw(6)  << amap.GetNPar();
  //   p_out.setf(std::ios::showpoint);
  //   int i=0;
  //   for( i=0; i<2; i++ ){
  //     if(i<amap.GetNPar())
  // 	p_out  << std::setprecision(5) << std::setw(12) << amap.GetPar(i);
  //     else
  // 	p_out  << std::setprecision(5) << std::setw(12) << 0.;
  //   }
  //   p_out << std::endl;
  // }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void DCTimeCorrMan::PrintMapBLDC( std::ostream &p_out )
{
  static const std::string funcname = "DCTimeCorrMan::PrintMapBLDC";
  PrintMapHeader(p_out);
  PrintMap(CID_BLC1a,p_out);
  PrintMap(CID_BLC1b,p_out);
  PrintMap(CID_BLC2a,p_out);
  PrintMap(CID_BLC2b,p_out);
  PrintMap(CID_BPC,p_out);
  PrintMap(CID_FDC1,p_out);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void DCTimeCorrMan::PrintMapCDC( std::ostream &p_out )
{
  static const std::string funcname = "DCTimeCorrMan::PrintMapCDC";
  PrintMapHeader(p_out);
  PrintMap(CID_CDC,p_out);
}
