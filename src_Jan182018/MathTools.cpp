// MathTools.cpp

#include "MathTools.h"

ClassImp(MathTools);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
MathTools::MathTools()
{
}

TVector3 MathTools::CalcHelixPos(const double par[5], double z)
{
  double phi2 = par[2]/par[4]*(par[3]-z);
  double x    = (par[0]+1./par[2] )*cos(par[1])-1./par[2]*cos(par[1]+phi2);
  double y    = (par[0]+1./par[2] )*sin(par[1])-1./par[2]*sin(par[1]+phi2);
  return TVector3(x,y,z);
}

double MathTools::CalcHelixDCA( const double par1[5], const double par2[5], TVector3 &vtx1, TVector3 &vtx2, TVector3 &vtx )
{
  const int NUM = 2;
  const int npoint = 100;
  const double initz = 0; //cm
  double region = 60.0; //+/-cm
  TVector3 pos1[npoint+1];
  TVector3 pos2[npoint+1];
  int num = 1;
  TVector3 now_vtx1;
  TVector3 now_vtx2;
  TVector3 now_vtx;
  double nowz = initz;
  double minl = 999.0;
  while(num<=NUM){
    //cerr<<"---"<<num<<endl;
    for(int i=0; i<npoint+1; i++){
      double z=nowz+(2*region/npoint)*(-npoint/2+i);
      //cerr<<z<<endl;
      pos1[i] = CalcHelixPos(par1,z);
      pos2[i] = CalcHelixPos(par2,z);
    }
    for(int i=0; i<npoint+1; i++){
      for(int j=0; j<npoint+1; j++){
	TVector3 diff = pos1[i]-pos2[j];
        double l = diff.Mag();
        if(l<minl){
          minl = l;
          now_vtx1 = pos1[i];
          now_vtx2 = pos2[j];
          now_vtx = now_vtx1+now_vtx2;
	  now_vtx *=0.5;
          nowz = now_vtx.z();
        }
      }
    }
    num++;
    region = 2*region/10;
  }
  vtx1 = now_vtx1;
  vtx2 = now_vtx2;
  vtx = now_vtx;

  // fine search by k.t.
  double step1 = 0.5;
  double step2 = 0.5;
  double dl = minl;
  double dldiff = 9999;
  TVector3 v1[4],v2[4];
  double dlp[4];
  int counter = 0;
  while( 0.0001<dldiff && counter<10000 ){
    double z1 = vtx1.Z();
    double z2 = vtx2.Z();
    v1[0] = CalcHelixPos( par1, z1+step1 );   v2[0] = CalcHelixPos( par2, z2+step2 );  dlp[0] = (v1[0]-v2[0]).Mag();
    v1[1] = CalcHelixPos( par1, z1+step1 );   v2[1] = CalcHelixPos( par2, z2-step2 );  dlp[1] = (v1[1]-v2[1]).Mag();
    v1[2] = CalcHelixPos( par1, z1-step1 );   v2[2] = CalcHelixPos( par2, z2+step2 );  dlp[2] = (v1[2]-v2[2]).Mag();
    v1[3] = CalcHelixPos( par1, z1-step1 );   v2[3] = CalcHelixPos( par2, z2-step2 );  dlp[3] = (v1[3]-v2[3]).Mag();
    bool valnewed=false;
    for( int i=0; i<4; i++ ){
      if( dlp[i]<dl ){
	vtx1 = v1[i]; vtx2 = v2[i];
	vtx = (vtx1+vtx2)*0.5;
	dldiff = fabs(dl-dlp[i]); dl = dlp[i]; 
	valnewed=true;
      }
    }
    if( !valnewed ){
      step1 *= 0.5; step2 *= 0.5;
    }
    counter++;
#if 0    
    std::cout << " minl:" << minl << " dl:" << dl << " step1:" << step1 << " step2:" << step2 
	      << " vz:" << vtx.Z() << " vz1:" << vtx1.Z() << " vz2:" << vtx2.Z()
	      << std::endl;
#endif
  }
#if 0    
  std::cout << " finally:" << dl << std::endl;
#endif

  if( counter==10000 ){
    std::cout << "!!! too many loops !!! " << std::endl;
  }
  return dl;
  //return minl;
}

TVector3 MathTools::CalcHelixMom(const double par[5], double z)
{
  //  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = -1*0.5; //T, "-1" is needed.
  double phi2 = (par[3]-z)*par[2]/par[4];
  double pt = fabs(1/par[2])*(Const*dMagneticField)/100.;
  double px = pt*(-1*sin(par[1]+phi2));
  double py = pt*(cos(par[1]+phi2));
  double pz = pt*(par[4]);
  return TVector3(px,py,pz);
}

double MathTools::CalcDeg(const double x,const  double y)
{
  double theta=atan( y/x )/3.141592*180.;
  if(x>=0&&y>=0) theta=theta;
  else if(x<0&&y>=0) theta=180+theta;
  else if(x<0&&y<0) theta=180+theta;
  else if(x>=0&&y<0) theta=360+theta;
  return theta;
}
