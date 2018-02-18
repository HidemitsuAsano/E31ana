// XTMapMan.h

#ifndef XTMapMan_h
#define XTMapMan_h 1

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include "TH1.h"
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
  bool Initialize();

  XTMapMan( const XTMapMan &right );
 private:

  std::string FileNameCDS;
  std::string FileNameBL;
  std::string ROOTNameBL;

  typedef std::map <unsigned int, XTMap> XTMapContainer;
  XTMapContainer xtContainer;
  typedef std::map <unsigned int, TH1F> XTHistContainer;
  XTHistContainer xthistContainer;

 public:
  void SetGainMapMan( const XTMapContainer container )  { xtContainer = container; }
  int nparam( const int &cid, const int &layer, const int &wire );
  double param( const int &cid, const int &layer, const int &wire, const int &i );
  std::string GetFileNameCDS() { return FileNameCDS; }
  std::string GetFileNameBL()  { return FileNameBL; }
  double CalcDriftLength( const int &cid, const int &layer, const int &wire, const double &dt );
  bool SetParam( const int &cid, const int &layer, const int &wire, const int &npar, double *par );

  void PrintMapHeader( std::ostream &p_out = std::cout );
  void PrintMap( const int &Cid = -1, std::ostream &p_out = std::cout );
  void PrintMapBL( std::ostream &p_out = std::cout );
  void PrintMapCDS( std::ostream &p_out = std::cout );

  // simulation
  double CalcDriftTime( const int &cid, const int &layer, const int &wire, const double &dl );

  ClassDef( XTMapMan, 1 );
};

#endif
