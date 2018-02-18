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
  double Gain, Offset;
 public:
  void SetGain( const double &gain ) { Gain = gain; }
  double GetGain() { return Gain; }
  void SetOffset( const double &offset ) { Offset = offset; }
  double GetOffset() { return Offset; }

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
