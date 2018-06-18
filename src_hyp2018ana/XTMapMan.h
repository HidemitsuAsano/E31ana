// XTMapMan.h

#ifndef XTMapMan_h
#define XTMapMan_h 1

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "BLDCWireMapMan.h"
#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include "TH1.h"
#include "TGraph.h"
#include "TROOT.h"

class XTMap : public TObject
{
 public:
  XTMap();
  ~XTMap() {};
 private:
  std::vector <double> Par;
 public:
  void SetParam( const double &par ) { Par.push_back(par); }
  int nParam() { return Par.size(); }
  std::vector <double> *GetParam() { return &Par; }
  double GetParam( int i );
  
  void Clear();

  ClassDef(XTMap,1);
};

class XTMapMan : public TObject
{
 public:
  XTMapMan();
  ~XTMapMan();

  void SetFileNameCDS( const std::string & filename );
  void SetFileNameBL( const std::string & filename );
  void SetROOTNameBL( const std::string & filename );
  bool Initialize(BLDCWireMapMan *wiremapman);
  void Verbosity(int level){verbosity = level;};

  XTMapMan( const XTMapMan &right );
 private:

  int verbosity;
  std::string FileNameCDS;
  std::string FileNameBL;
  std::string ROOTNameBL;

  typedef std::map <unsigned int, TGraph > XTMapContainer2;
  XTMapContainer2 xtContainer2;
  typedef std::map <unsigned int, XTMap> XTMapContainer;
  XTMapContainer xtContainer;
  typedef std::map <unsigned int, TH1F> XTHistContainer;
  XTHistContainer xthistContainer;
  
  double driftlengthmax[6];//max. drift length in each BLDC
 public:
  void SetGainMapMan( const XTMapContainer container )  { xtContainer = container; }

  int nparam( const int &cid, const int &layer, const int &wire );
  int nparamBL( const int &cid, const int &layer, const int &wire );//H.Asano
  double param( const int &cid, const int &layer, const int &wire, const int &i );
  double paramBL( const int &cid, const int &layer, const int &wire, const int &i );//H.Asano

  std::string GetFileNameCDS() { return FileNameCDS; }
  std::string GetFileNameBL()  { return FileNameBL; }
  std::string GetROOTNameBL()  { return ROOTNameBL; }

  double CalcDriftLength( const int &cid, const int &layer, const int &wire, const double &dt );
  double CalcDxDt( const int &cid, const int &layer, const int &wire, const double &dt );
  bool SetParam( const int &cid, const int &layer, const int &wire, const int &npar, double *par );
  bool SetParamBL( const int &cid, const int &layer, const int &wire, const int &npar, double *par, double *threshold);

  void PrintMapHeader( std::ostream &p_out = std::cout );
  void PrintMapHeaderBL( std::ostream &p_out = std::cout );
  void PrintMap( const int &Cid = -1, std::ostream &p_out = std::cout );
  void PrintMapBL( std::ostream &p_out = std::cout );
  void PrintMapCDS( std::ostream &p_out = std::cout );

  // simulation
  double CalcDriftTime( const int &cid, const int &layer, const int &wire, const double &dl );

  ClassDef( XTMapMan, 1 );
};

#endif
