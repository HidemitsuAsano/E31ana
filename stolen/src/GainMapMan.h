// GainMapMan.h

#ifndef GainMapMan_h
#define GainMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class GainMap : public TObject
{
 public:
  GainMap();
  ~GainMap() {};
 private:
  double Param[3]; // 
 public:
  void SetGain( const double &gain ) { Param[1] = gain; }
  double GetGain() { return Param[1]; }
  void SetOffset( const double &offset ) { Param[0] = offset; }
  double GetOffset() { return Param[0]; }
  void SetParam( const int &i,const double &offset ) { Param[i] = offset; }
  double GetParam( const int &i ) { return Param[i]; }
  void SetParams( const double param[3] ) { for(int i=0;i<3;i++) Param[i]=param[i]; }
  void GetParams( double param[3] ) { for(int i=0;i<3;i++) param[i]=Param[i]; }
  ClassDef(GainMap,1);
};

class GainMapMan : public TObject
{
 public:
  GainMapMan();
  ~GainMapMan();

  void SetFileNameCDS( const std::string & filename );
  void SetFileNameBL(  const std::string & filename );
  void SetFileNameSDD(  const std::string & filename );
  bool Initialize();
  bool ReadFile(const std::string &filename);

  GainMapMan( const GainMapMan &right );
 private:

  std::string FileNameCDS;
  std::string FileNameBL;
  std::string FileNameSDD;

  typedef std::map <unsigned int, GainMap> GainMapContainer;
  GainMapContainer gainContainer;

 public:
  void SetGainMapMan( const GainMapContainer container )  { gainContainer = container; }

  std::string GetFileNameCDS() { return FileNameCDS; }
  std::string GetFileNameBL() { return FileNameBL; }
  std::string GetFileNameSDD() { return FileNameSDD; }
  bool GetParam( const int &c, const int &n, const int &a, double &gain, double &offset );
  bool SetParam( const int &c, const int &n, const int &a, const double &gain, const double &offset );
  bool GetParam( const int &c, const int &n, const int &a, double &pol0, double &pol1, double &pol2 );
  bool SetParam( const int &c, const int &n, const int &a, const double &pol0, const double &pol1, const double &pol2 );
  double CalcCValue( const int &c, const int &n, const int &a, const int &at, const double &val, const bool &smear=false );

  void PrintMapHeader( std::ostream &p_out = std::cout );
  void PrintMap( const int &crate=-1, std::ostream &p_out = std::cout );
  void PrintMapCDS( std::ostream &p_out = std::cout );
  void PrintMapBL( std::ostream &p_out = std::cout );
  void PrintMapSDD( std::ostream &p_out = std::cout );

  // simulation
  int CalcDATValue( const int &c, const int &n, const int &a, const int &at, const double &val );

  ClassDef( GainMapMan, 1 );
};

#endif
