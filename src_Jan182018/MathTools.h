// MathTools.h

#ifndef MathTools_h
#define MathTools_h 1

#include <string>
#include <iostream>
#include <cmath>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TVector3.h"
#include "GlobalVariables.h"

class MathTools: public TObject
{
 public:
  MathTools();
  virtual ~MathTools(){};

 public:
  TVector3 CalcHelixPos(const double par[5], double z);
  double CalcHelixDCA( const double par1[5], const double par2[5], TVector3 &vtx1, TVector3 &vtx2, TVector3 &vtx );
  TVector3 CalcHelixMom(const double par[5], double z);
  double CalcDeg(const double x,const double y);

  ClassDef(MathTools,1);
};

#endif
