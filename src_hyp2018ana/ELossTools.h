// ELossTools.h
#ifndef ELossTools_h
#define ELossTools_h 1

#include <string>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TMath.h"
#include "GlobalVariables.h"
#include "GeomTools.h"
#include "MathTools.h"
#include "CDSTrack.h"

namespace ELossTools
{
  bool CalcHelixElossToVertex(const double param[5],const TVector3 &vertex, const double &momin, const double &mass, double &momout, double &tof);
  bool CalcHelixElossToVertexTGeo(const double param[5],const TVector3 &vertex, const double &momin, const double &mass, double &momout, double &tof,const double &rstart=(15.1+(53.-15.1)/2.)); 
  bool CalcHelixElossPointToPointTGeo(const double param[5],const TVector3 &posin, const TVector3 &posout,const double &momin, const double &mass, double &momout, double &tof,double offs=0.);  

  bool CalcHelixElossPointToPointTGeo(const TVector3 &posin, const TVector3 &momin, const double &mass, const double &z, const double &field, TVector3 &posout, TVector3 &momout, double &flightlength, double &tof);  
  
  bool CalcHelixElossToNextBoundary(const double param[5],TVector3 &pos1,const double &momin, const double &mass,
				    double &momout,  double &tof, double &length,int &id, int sign=1);
  bool CalcElossBeamTGeo(const TVector3 &in, const TVector3 &out, const double &momin, const double &mass, double &momout,double &tof);
  bool CalcElossForwardTGeo(const TVector3 &in, const TVector3 &fdcpos,const double &flength, const double &momin, const double &mass, double &momout,double &tof);
  bool CalcdE(const double param[5], const TVector3 &vertex, const double &rmax, const double &rmin,
	      const double &momentum, const double &mass , const TString &material, double &momout, double &tof);
  
  bool CalcdE(const double &momentum, const double &mass, const double &length, const TString &material, double &momout, const int &sign, double &tof, bool CORR=true);
  double CalcdEdX(const double &beta, const TString &matname, bool CORR=true);  
  double CalcdEdX(const double &beta,const double &rho, const double &I, const double &Z_A, const double &Z,
		  const double &C0, const double &a, const double &M, const double &X0, const double &X1, bool CORR=true);  
  inline double CalcdEdX(const double &beta,const double &rho, const double &I, const double &Z_A);
};

inline double ELossTools::CalcdEdX(const double &beta,const double &rho, const double &I, const double &Z_A){
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511; // MeV/c^2
  double gamma_2=1/(1-pow(beta,2.0));
  double W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  int z=1;
  double logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0);
  double val=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return val/1000.; // MeV->GeV
}
#endif
