// GateMapMan.h

#ifndef GateMapMan_h
#define GateMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class GateMap : public TObject
{
 public:
  GateMap();
  ~GateMap() {};
 private:
  double Low, High;
  int Type;
 public:
  void SetLow( const double &low ) { Low = low; }
  double GetLow() { return Low; }
  void SetHigh( const double &high ) { High = high; }
  double GetHigh() { return High; }
  void SetType( const int &type ) { Type=type; }
  int GetType() { return Type; }

  ClassDef(GateMap,1);
};

typedef std::map <unsigned int, GateMap> GateMapContainer;
class GateMapMan : public TObject
{
 public:
  GateMapMan();
  ~GateMapMan();

  void SetFileNameCDS( const std::string & filename );
  void SetFileNameBL(  const std::string & filename );
  void SetFileNameSDD(  const std::string & filename );
  bool Initialize();

  GateMapMan( const GateMapMan &right );
 private:
  std::string FileNameCDS;
  std::string FileNameBL;
  std::string FileNameSDD;

  GateMapContainer gateContainer;

 public:
  void SetGateMapMan( const GateMapContainer container )  { gateContainer = container; }  
  std::string GetFileNameCDS() { return FileNameCDS; }
  std::string GetFileNameBL() { return FileNameBL; }
  std::string GetFileNameSDD() { return FileNameSDD; }
  bool GetParam( const int &cid, const int &seg, const int &ud, const int &at, int &type, double &low, double &hifht );
  bool SetParam( const int &cid, const int &seg, const int &ud, const int &at, const int &type, const double &low, const double &high );
  bool CheckRange( const int &cid, const int &seg, const int &ud, const int &at, const int &type,const double &val );

  bool GetParam( const int &cid, const int &layer, const int &wire, int &type, double &low, double &hifht );
  bool SetParam( const int &cid, const int &layer, const int &wire, const int &type, const double &low, const double &high );
  bool CheckRange( const int &cid, const int &layer, const int &wire, const int &type, const double &val );

  void PrintMapHeader( std::ostream &p_out = std::cout );
  void PrintMap( const int &cid=-1, std::ostream &p_out = std::cout );
  void PrintMapCDS( std::ostream &p_out = std::cout );
  void PrintMapBL( std::ostream &p_out = std::cout );
  void PrintMapSDD( std::ostream &p_out = std::cout );

  ClassDef( GateMapMan, 1 );
};

#endif
