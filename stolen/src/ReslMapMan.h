// ReslMapMan.h

#ifndef RESLMAN_h
#define RESLMAN_h 1

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include <TRandom.h>

class ReslMap : public TObject
{
 public:
  ReslMap();
  ~ReslMap();

 private:
  double TResl, EResl;

 public:
  double tresl() { return TResl; }
  double eresl() { return EResl; }
  void SetTResl( const double &t ) { TResl = t; }
  void SetEResl( const double &e ) { EResl = e; }

  ClassDef(ReslMap,1);
};

class ReslMapMan : public TObject
{
 public:
  ReslMapMan();
  ReslMapMan(const std::string & filename );
  ~ReslMapMan();

  void SetFileName( const std::string & filename );
  bool Initialize();

  //ReslMapMan( const ReslMapMan &right );

 private:
  std::string FileName;
  //typedef std::map <unsigned int, double> ReslContainer;
  typedef std::map <unsigned int, ReslMap> ReslContainer;
  ReslContainer reslContainer;

 public:
  std::string GetFileName() { return FileName; }
  void GetResolution( const int &cid, const int &layer, const int &wire, double &tresl, double &eresl );
  bool SetParam( const int &cid, const int &layer, const int &wire, const double &tresl, const double &eresl );
  bool GetParam( const int &cid, const int &layer, const int &wire, double &tresl, double &eresl );

  void PrintMap( const int &Cid = -1, std::ostream &p_out = std::cout );

  ClassDef( ReslMapMan, 1 );
};

#endif
