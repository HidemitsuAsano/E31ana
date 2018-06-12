// SlewingMapMan.h

#ifndef SlewingMapMan_h
#define SlewingMapMan_h 1

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"

class SlewingMap : public TObject
{
 public:
  SlewingMap();
  ~SlewingMap() {};
 private:
  //int nPar;
  //double Par[MaxParam];
  std::vector <double> Par;
  int Type;
 public:
  //void SetPar( const int &i, const double &par ) { Par[i] = par; }
  void SetPar( const int &i, const double &par ) { if( i<(int)Par.size() ) Par[i] = par; }
  //int GetNPar() const { return nPar; }
  int GetNPar() const { return Par.size(); }
  //double GetPar( const int &i ) { return Par[i]; }
  double GetPar( const int &i ) { return i<(int)Par.size() ? Par[i] : 0; }
  std::vector <double> GetPar() { return Par; }
  void SetType( const int &i ) { Type = i; }
  int GetType() const { return Type; }
  void SetParam( const int &npar, double *par );
  void SetParam( const int &npar, const std::vector <double> &par );
  
  ClassDef(SlewingMap,1);
};

typedef std::map <unsigned int, SlewingMap> SlewingMapContainer;
class SlewingMapMan : public TObject
{
 public:
  SlewingMapMan();
  ~SlewingMapMan();

  void SetFileNameCDS( const std::string & filename );
  void SetFileNameBL( const std::string & filename );
  void SetFileNameSDD( const std::string & filename );
  void SetFileName( const std::string &filename){ FileNames.push_back(filename); }
  bool Initialize();
  bool ReadFile(const std::string &filename);

  SlewingMapMan( const SlewingMapMan &right );
 private:

  std::string FileNameCDS;
  std::string FileNameBL;
  std::string FileNameSDD;
  std::vector<std::string> FileNames;

  SlewingMapContainer slewingContainer;

 public:
  void SetSlewingMapMan( const SlewingMapContainer container )  { slewingContainer = container; }

  std::string GetFileNameCDS() { return FileNameCDS; }
  std::string GetFileNameBL() { return FileNameBL; }
  std::string GetFileNameSDD() { return FileNameSDD; }
  int nFile() const { return FileNames.size(); };
  std::string GetFileName(const int &i){ return FileNames[i]; }

  bool isParam( const int &cid, const int &seg, const int &ud, const int &ith=0);
  bool GetParam( const int &cid, const int &seg, const int &ud, const int &ith, int &type,
		 int &npar, std::vector <double> &par );
  bool GetParam( const int &cid, const int &seg, const int &ud, const int &ith, int &type,
		 int &npar, double *par );
  bool SetParam( const int &cid, const int &seg, const int &ud, const int &ith, const int &type,
		 const int &npar, double *par );
  bool SetParam( const int &cid, const int &seg, const int &ud, const int &ith, const int &type,
		 const int &npar, const std::vector <double> &par );
  double CalcCValue( const int &cid, const int &seg, const int &ud,
		     const double &time, const double &de );

  double CalcDATValue( const int &cid, const int &seg, const int &ud,
		       const double &ctime, const double &de );

  void PrintMapHeader( std::ostream &p_out = std::cout );
  void PrintMap( const int &Cid=-1, std::ostream &p_out = std::cout );
  void PrintMapBL( std::ostream &p_out = std::cout );
  void PrintMapCDS( std::ostream &p_out = std::cout );

  ClassDef( SlewingMapMan, 1 );
};

#endif
