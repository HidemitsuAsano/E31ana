// DCTimeCorrMan.h

#ifndef DCTimeCorrMan_h
#define DCTimeCorrMan_h 1

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"
#include "TGraph.h"
#include "TFile.h"
#include "TROOT.h"

class DCTimeCorr : public TObject
{
 public:
  DCTimeCorr();
  ~DCTimeCorr() {};
 private:
  std::vector <double> Par;
  int Type;
 public:
  void SetPar( const int &i, const double &par ) { if( i<(int)Par.size() ) Par[i] = par; }
  int GetNPar() const { return Par.size(); }
  double GetPar( const int &i ) { return i<(int)Par.size() ? Par[i] : 0; }

  void SetType( const int &i ) { Type = i; }
  int GetType() const { return Type; }

  void SetParam( const int &npar, double *par );
  
  ClassDef(DCTimeCorr,1);
};

class DCTimeCorrMan : public TObject
{
 public:
  DCTimeCorrMan();
  ~DCTimeCorrMan();

  void SetFileNameCDC( const std::string & filename );
  void SetFileNameBLDC( const std::string & filename );
  bool Initialize();
  
  DCTimeCorrMan( const DCTimeCorrMan &right );
 private:

  std::string FileNameCDC;
  std::string FileNameBLDC;

  typedef std::map <unsigned int, TGraph> DCTimeCorrHistContainer;
  DCTimeCorrHistContainer dctimecorrContainer;
  
 public:
  void SetDCTimeCorrMan( const DCTimeCorrHistContainer container )  { dctimecorrContainer = container; }

  std::string GetFileNameCDC() { return FileNameCDC; }
  std::string GetFileNameBLDC() { return FileNameBLDC; }

  double CalcCValue( const int &cid, const int &layer, const int &wire,
		     const double &timemean, const double &timesub );
  double CalcDATValue( const int &cid, const int &layer, const int &wire,
		       const double &timemean, const double &timesub );

  void PrintMapHeader( std::ostream &p_out = std::cout );
  void PrintMap( const int &Cid=-1, std::ostream &p_out = std::cout );
  void PrintMapBLDC( std::ostream &p_out = std::cout );
  void PrintMapCDC( std::ostream &p_out = std::cout );

  ClassDef( DCTimeCorrMan, 1 );
};

#endif
