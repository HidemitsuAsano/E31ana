// CounterMapMan.h

#ifndef CounterMapMan_h
#define CounterMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"

#define MMAXSMP 20
#define MAXSLOT 23
#define NORMAL 1
#define DRT 0
class CounterMapMan : public TObject
{
 public:
  CounterMapMan();
  CounterMapMan( const std::string & filename );
  ~CounterMapMan();

  void SetFileName( const std::string & filename );
  bool Initialize();

  bool GetInfo( int c, int n, int a, int &cid, int &lay, int &seg, int &at, int &ud );
  int GetCID( int c, int n, int a );
  std::string GetName( const int &c, const int &n, const int &a );

  // for Hodoscope or Cherenkov
  bool GetCNA( int cid, int seg, int at, int ud, int &c, int &n, int &a );
  std::string GetName( const int &cid, const int &seg, const int &at, const int &ud);
  // for DC
  bool GetCNA( int cid, int layer, int wire, int at, int ud, int &c, int &n, int &a );
  std::string GetName( const int &cid, const int &lay, const int &wire, const int &at, const int &ud );

  int GetCrateNum( int address );
  int GetSMPAddress( int c );
  
  CounterMapMan( const CounterMapMan &right );
 private:

  std::string FileName;
  int NSMP;
  int NCH_SCA;
  int CRATE_TYPE[MMAXSMP][MAXSLOT];

  typedef std::map <unsigned int, unsigned int> fCounterMapContainer;
  typedef std::map <unsigned int, unsigned int> bCounterMapContainer;
  fCounterMapContainer fContainer;
  bCounterMapContainer bContainer;

  typedef std::map <unsigned int, unsigned int> fCrateDefContainer;
  typedef std::map <unsigned int, unsigned int> bCrateDefContainer;
  fCrateDefContainer fCrateDef;
  bCrateDefContainer bCrateDef;

  typedef std::map <unsigned int, std::string> nameCNAMapContainer;
  typedef std::map <unsigned int, std::string> nameCounterMapContainer;
  nameCNAMapContainer nameCNAContainer;
  nameCounterMapContainer nameCounterContainer;

 public:
  std::string GetFileName() { return FileName; }
  int GetNumSMP() { return NSMP; }
  int GetNumScaler() { return NCH_SCA; }
  int GetCrateType( const int &cr, const int &sl ) { return CRATE_TYPE[cr][sl-1]; } //sl 1 origin
  void PrintSimpleMap( std::ostream &p_out = std::cout );
  void PrintMap();

  void Clear();
  
  ClassDef( CounterMapMan, 1 );
};

#endif
