// TransferMatrixMan.h

#ifndef TransferMatrixMan_h
#define TransferMatrixMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"
#include "TMatrixD.h"

class TransferMatrixMan : public TObject
{
 public:
  TransferMatrixMan();
  TransferMatrixMan( const std::string & filename );
  TransferMatrixMan( const TransferMatrixMan &right );
  ~TransferMatrixMan();

  void SetFileName( const std::string & filename );
  bool Initialize();

 private:

  std::string FileName;
  double CentralMomentum;
  double BLC1VIele[36];
  TMatrixD BLC1VIMatrix;
  TMatrixD BLC2VIMatrix;
  TMatrixD D5Matrix;
  double D5Matrix1st[36];
  double D5Matrix2nd[6][36];
  TMatrixD BLC2BLC1Matrix1st;
  TMatrixD BLC2BLC1Matrix2nd[6];

 public:
  std::string GetFileName() { return FileName; }
  double GetCentralMomentum() const { return CentralMomentum; } 
  const double* GetD5Matrix() const { return D5Matrix1st; }
  const double* GetD5Matrix2nd(const int &i) const { return (i<6) ? D5Matrix2nd[i] : 0 ; }
  void CalcBLC1toBLC2(double *parblc1, double *parblc2,const int &order=1);
  void CalcBLC2toBLC1(double *parblc2, double *parblc1, const int &order=1);
  void CalcBLC1toVI(double *parblc1, double *parvi,const int &order=1);
  void CalcBLC2toVI(double *parblc2, double *parvi,const int &order=1);
  void Clear();
  
  ClassDef( TransferMatrixMan, 1 );
};

#endif
