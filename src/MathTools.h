// MathTools.h
#ifndef MathTools_h
#define MathTools_h 1

#include <string>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TMath.h"
#include "GlobalVariables.h"
#include "GeomTools.h"

namespace MathTools
{
  inline TVector3 CalcHelixPosatZ(const double par[5],const double &z);
  inline TVector3 CalcHelixPosatR(const double par[5],const double &z);
  inline TVector3 CalcHelixGPosatR(const double par[5],const double gpar[4], const double &z);
  inline TVector3 CalcLinePosatR(const TVector3 &pos, const TVector3 &dir, const double &z);
  inline void CalcHelixXYtoZ( const double param[5], const double &x,const  double &y,double &z );

  inline double CalcHelixPhiatZ(const double par[5],const double &z);
  inline double CalcHelixPhiatR(const double par[5],const double &z);

  inline TVector3 CalcHelixMom(const double par[5], double z);  
  inline double   CalcHelixPhi( const double &x,const double &y,const double *par);
  inline double   CalcHelixPhi( const TVector3 &pos,const double *par);

  inline double CalcHelixArc(const double par[5], const double &r1,const double &r2);
  inline double CalcHelixArc(const double par[5], const TVector3 &pos1, const TVector3 &pos2);
  inline TVector3 CalcHelixStep(const double par[5], const TVector3 &pos, const double &step);  

  inline TVector3 GetPosition(const double &helixphi,const double *par);
  inline bool GetHelix(const TVector3 &pos,const TVector3 &mom, int &charge, const double &field, double par[5]);

  double CalcHelixDCA( const double par1[5], const double par2[5], TVector3 &vtx1, TVector3 &vtx2, TVector3 &vtx );
  
  inline void ConvertCDCPointLtoG(const TVector3 &in, TVector3 &out, const double *gpar);
  
  inline double CalcDeg(const double x,const double y);
  inline double CalcRad(const double x,const double y);
  
  bool ChangePivot(const TVector3 &oldpivot, const TVector3 &newpivot,
		   const double oldpar[5], double newpar[5],const int &charge);
  
  
  inline void PointToCircle(const double &x,const double &y,
			    const double &radius,
			    const double &x_cen,const double &y_cen,
			    double &dis,
			    double &xest,double &yest);
  inline void PointToCircle(const double *par,
			    const double &x,const double &y,
			    double &xest,double &yest);
  inline void PointToCircle(const double *par,
			    const double &x,const double &y,
			    double &xest,double &yest, double &dist);
  inline void PointToLine( const TVector3 &p,
			   const TVector3 &x, const TVector3 &a,
			   double &dist,TVector3 &xest );
  
  inline void LineToLine( const TVector3 &x1, const TVector3 &a1,
			  const TVector3 &x2, const TVector3 &a2,
			  const double &dl,
			  double &dist,
			  TVector3 &xest, TVector3 &next );    
  bool LineToCircle(const double &a,const double &b, const double &c, 
		    const double &rho ,const double &xc, const double &yc, 
		    double &x_p,double &y_p,double &x_n,double &y_n);
  inline bool PointToHelix(const TVector3 &hitpos, const double *par,
			   TVector3 &fitpos,double &dis); 
  inline bool PointToHelix(const TVector3 &hitpos, TVector3 &fitpos,const double *par);   
  bool LineToHelix(const TVector3 &a, const TVector3 &dline, 
		   const double *par, TVector3 &lnest,
		   TVector3 &hnest, double &dis);
  bool HelixToHelix( const double *par1,const double *par2, TVector3 &fitpos1,TVector3 &fitpos2, double &dis);
  
  bool HelixToHelixWresl( const double *par1,const double *par2, TVector3 &fitpos1,TVector3 &fitpos2, double &dis);


  inline bool CircleFit(const double *x, const double *y,const double *weight, const int &numofhit, double *par, double &chi2);
  
  inline double dfunc_PTH(const TVector3 &pos,const double &helixphi,const double *par);
  double dfunc_LTH(const TVector3 &lpos,const TVector3 &dline,const double &helixphi,const double *par);  
};

inline bool MathTools::GetHelix(const TVector3 &pos,const TVector3 &mom, int &charge, const double &field, double par[5])
{
  if(mom.Pt()<1e-9) return false;
  charge/=TMath::Abs(charge);
  par[2]=charge/mom.Pt()*Const*field/100.;
  par[4]=TMath::Tan(TMath::ACos(mom.Pt()/mom.Mag()));

  double cx=pos.X()+mom.Y()/mom.Pt()/par[2];
  double cy=pos.Y()-mom.X()/mom.Pt()/par[2];
  
  par[0]=sqrt(cx*cx+cy*cy)+TMath::Abs(1/par[2]);
  par[1]=TMath::ATan2(cy,cx);
  if(par[2]>0){
    par[0]*=-1;
    par[1]+=TMath::Pi();
  }
  if(par[1]<0) par[1]+=2*TMath::Pi();
  else if(par[1]>2*TMath::Pi()) par[1]-=2*TMath::Pi();
  double phi= MathTools::CalcHelixPhi(pos.X(),pos.Y(),par);

  par[3]=pos.Z()+1./par[2]*par[4]*phi;

  /* for(int i=0;i<5;i++) */
  /*   std::cout<<par[i]<<"  "; */
  /* std::cout<<std::endl; */
}

inline void MathTools::ConvertCDCPointLtoG(const TVector3 &in, TVector3 &out, const double *gpar){
  out=in+TVector3(gpar[0],gpar[1],gpar[2]);
  out.RotateZ(gpar[3]);
}

inline bool MathTools::CircleFit(const double *x,const double *y, const double *weight,const int &npoints, double *par, double &chi2)
{
  chi2=9999;
  int i;
  double xx, yy, xx2, yy2;
  double f, g, h, p, q, t, g0, g02, a, b, c, d;
  double xroot, ff, fp, xd, yd, g1;
  double dx, dy, dradius2, xnom;
  
  double xgravity = 0.0;
  double ygravity = 0.0;
  double x2 = 0.0;
  double y2 = 0.0;
  double xy = 0.0;
  double xx2y2 = 0.0;
  double yx2y2 = 0.0;
  double x2y22 = 0.0;
  double radius2 = 0.0;
  
  double mVariance = 0.0;
  
  //  const int MAX_HIT_IN_TRACK = 20;
  /* if(npoints >= MAX_HIT_IN_TRACK){ */
  /*   //fprintf(stderr, "circleFit: npoints %d >= MAX\n",npoints); */
  /*   return false; */
  /* } */
  
  if (npoints <= 3){
    //fprintf(stderr, "circleFit: npoints %d <= 3\n",npoints);
    return false;
  }
  // step 1. transfer the origin of the coordinate system to the center of gravity of the se Pos[n]    
  for (i=0; i<npoints; i++) {
    xgravity += x[i];
    ygravity += y[i];
  }
  xgravity /= npoints;
  ygravity /= npoints;
  
  
  // step 2. calcurate the gauss-bracket values for [x^2] etc.
  for (i=0; i<npoints; i++) {
    xx  = x[i]-xgravity;
    yy  = y[i]-ygravity;
    xx2 = xx*xx;
    yy2 = yy*yy;
    x2  += xx2;
    y2  += yy2;
    xy  += xx*yy;
    xx2y2 += xx*(xx2+yy2);
    yx2y2 += yy*(xx2+yy2);
    x2y22 += (xx2+yy2)*(xx2+yy2);
  }
  if (xy == 0.){
    //fprintf(stderr, "circleFit: xy = %f,    grav=%f, %f\n",xy,xgravity,ygravity);
    return false;
  }
  
  // step 3. calculate coefficients for eq (11)
  f = (3.*x2+y2)/npoints;
  g = (x2+3.*y2)/npoints;
  h = 2*xy/npoints;
  p = xx2y2/npoints;
  q = yx2y2/npoints;
  t = x2y22/npoints;
  g0 = (x2+y2)/npoints;
  g02 = g0*g0;
  a = -4.0;
  b = (f*g-t-h*h)/g02;
  c = (t*(f+g)-2.*(p*p+q*q))/(g02*g0);
  d = (t*(h*h-f*g)+2.*(p*p*g+q*q*f)-4.*p*q*h)/(g02*g02);
  
  //step 4. solve (11) by the Newton method
  xroot = 1.0;
  for (i=0; i<5; i++) {
    ff = (((xroot+a)*xroot+b)*xroot+c)*xroot+d;
    fp = ((4.*xroot+3.*a)*xroot+2.*b)*xroot+c;
    xroot -= ff/fp;
  }
  
  //step 5. calculate parameters
  g1 = xroot*g0;
  xnom = (g-g1)*(f-g1)-h*h;
  if (xnom == 0.){
    //fprintf(stderr, "circleFit: xnom1 = %f\n",xnom);
    return false;
  }
  
  yd = (q*(f-g1)-h*p)/xnom;
  xnom = f-g1;
  if (xnom == 0.){
    //fprintf(stderr, "circleFit: xnom2 = %f\n",xnom);
    return false;
  }
  
  xd = (p-h*yd )/xnom;
  
  radius2 = xd*xd+yd*yd+g1;

  par[0] = xd+xgravity;
  par[1] = yd+ygravity;
  par[2] = sqrt(radius2);  

  for (i=0; i<npoints; i++) {
    dx = x[i]-(par[0]);
    dy = y[i]-(par[1]);
    dradius2 = dx*dx+dy*dy;      
    mVariance += (dradius2+radius2-2.*sqrt(dradius2*radius2))/weight[i]/weight[i]; 
  }
  chi2 = mVariance / (npoints-3);
  return true;
}

inline TVector3 MathTools::GetPosition(const double &helixphi,const double *par)
{
  // drho=par[0], phi0=par[1], rho(kappa)=par[2], dz=par[3], tlam=par[4];
  // x= drho * cos (phi0) +  1/rho *( cos (phi0) - cos( phi0 +helixphi );
  // y= drho * sin (phi0) +  1/rho *( sin (phi0) - sin( phi0 +helixphi );
  // z= dz - 1/kappa * tlam * helixphi;
  // par[2]>0, magfield>0 -> positive charge, propagate clockwize
  TVector3 pos;
  if(par[2]!=0){
    pos.SetX( par[0]*cos(par[1])+1./par[2]*( cos(par[1])-cos(par[1]+helixphi) ) );
    pos.SetY( par[0]*sin(par[1])+1./par[2]*( sin(par[1])-sin(par[1]+helixphi) ) );
    pos.SetZ( par[3]-1./par[2]*par[4]*helixphi );
  }
  else if(par[2]==0)
    {
      pos.SetXYZ(par[0]*cos(par[1])-helixphi*(-sin(par[1]) ),
                 par[0]*sin(par[1])-helixphi*( cos(par[1]) ),
                 par[3]+par[4]*helixphi);
    }

  return pos;
}

inline void MathTools::PointToCircle(const double *par,const double &x,const double &y,
				     double &xest,double &yest)	      
{
  Double_t cx, cy;
  cx = (par[0]+1./par[2])*cos(par[1]);
  cy = (par[0]+1./par[2])*sin(par[1]);

  Double_t nx, ny, nn;
  nx = x - cx;
  ny = y - cy;
  nn = sqrt( nx*nx + ny*ny );
  nx /= nn;
  ny /= nn;

  xest = cx +fabs( 1./par[2])*nx;
  yest = cy +fabs( 1./par[2])*ny;
}

inline void MathTools::PointToCircle(const double *par,const double &x,const double &y,
				     double &xest,double &yest,double &dis) 
{
  Double_t cx, cy;
  double radius=1./par[2];
  cx = (par[0]+radius)*cos(par[1]);
  cy = (par[0]+radius)*sin(par[1]);

  Double_t nx, ny, nn;
  nx = x - cx;
  ny = y - cy;
  nn = sqrt( nx*nx + ny*ny );
  nx /= nn;
  ny /= nn;

  xest = cx +fabs(radius)*nx;
  yest = cy +fabs(radius)*ny;

  dis=fabs(nn-fabs(radius));
}

inline void MathTools::PointToCircle(const double &x,const double &y,
				     const double &radius,
				     const double &x_cen,const double &y_cen,
				     double &dis,
				     double &xest,double &yest)	      
{
  Double_t nx, ny;
  nx = x - x_cen;
  ny = y - y_cen;
  double rtmp=sqrt( nx*nx + ny*ny );
  nx /= rtmp;
  ny /= rtmp;
  dis=fabs(rtmp-radius);  
  xest=x_cen+radius*nx;
  yest=y_cen+radius*ny;
}

inline void MathTools::CalcHelixXYtoZ( const double param[5], const double &x,const  double &y,double &z )
{
  if(param[2]!=0)
    {
      double phi=MathTools::CalcHelixPhi(x,y,param);
      z=param[3]-1./param[2]*param[4]*phi;
    }
  else z=-999;
}

inline TVector3 MathTools::CalcHelixPosatZ(const double par[5],const double &z)
{
  double phi2 = par[2]/par[4]*(par[3]-z);
  double x    = (par[0]+1./par[2] )*cos(par[1])-1./par[2]*cos(par[1]+phi2);
  double y    = (par[0]+1./par[2] )*sin(par[1])-1./par[2]*sin(par[1]+phi2);
  return TVector3(x,y,z);
}

inline TVector3 MathTools::CalcLinePosatR(const TVector3 &pos, const TVector3 &dir,const double &r)
{
  double posx=pos.X();
  double posy=pos.Y();
  double dirx=dir.X();
  double diry=dir.Y();
  // (posx+s*dirx)^2+(posy+diry*s)^2=r^2
  double a=dirx*dirx+diry*diry;
  double b=2*(posx*dirx+posy*diry);
  double c=posx*posx+posy*posy-r*r;
  double s=(-b+sqrt(b*b-4*a*c))/2./a;

  TVector3 out=pos+s*dir;
#if 0
  std::cout<<"----------------"<<std::endl;
  std::cout<<a<<"  "<<b<<"  "<<c<<"  "<<s<<std::endl;
  pos.Print();
  dir.Print();
  out.Print();
#endif
  return out;
}

inline double MathTools::CalcHelixPhiatZ(const double par[5],const  double &z)
{
  double phi2 = par[2]/par[4]*(par[3]-z);
  return phi2;
}

inline TVector3 MathTools::CalcHelixStep(const double par[5],const TVector3 &lastpos, const double &step){
  double sign=1;
  //  double phi2=CalcHelixPhiatZ(par,lastpos.Z());
  double phi = CalcHelixPhi(lastpos,par);
  /* if(TMath::Abs(phi-tmp)>0.01) */
  double dphi= (step*sign*par[2])/sqrt(1+par[4]*par[4]);
  //  std::cout<<phi<<"   "<<phi2<<"  "<<dphi<<std::endl;
  return MathTools::GetPosition(phi+dphi,par);
}
inline TVector3 MathTools::CalcHelixPosatR(const double par[5], const double &r)
{
  double cx=(par[0]+1./par[2])*cos(par[1]);
  double cy=(par[0]+1./par[2])*sin(par[1]);
  double rho=fabs(1./par[2]);
  TVector3 center(cx,cy,0);
  double angle=TMath::ACos( (center.Mag2()+rho*rho-r*r)/2./center.Mag()/rho);
  TVector3 tmp=-center*(1./center.Mag())*rho;
  if(par[2]<0) tmp.RotateZ(angle);
  if(par[2]>0) tmp.RotateZ(-angle);
  double z;
  tmp+=center;
  MathTools::CalcHelixXYtoZ(par,tmp.X(),tmp.Y(),z);
  tmp.SetZ(z);
  return tmp;
}

inline double MathTools::CalcHelixPhiatR(const double par[5],const double &r)
{
  double cx=( par[0]+1./par[2] )*cos(par[1]);
  double cy=( par[0]+1./par[2] )*sin(par[1]);
  double rho=fabs(1./par[2]);
  TVector3 center(cx,cy,0);
  double angle=TMath::ACos( (center.Mag2()+rho*rho-r*r)/2/center.Mag()/rho );
  //  if(par[2]>0 ) angle=angle;
  //  else  angle=TMath::Pi()-angle;
  return angle;
}

inline double MathTools::CalcHelixArc(const double par[5],const double &r1,const double &r2)
{
  double cx=( par[0]+1./par[2] )*cos(par[1]);
  double cy=( par[0]+1./par[2] )*sin(par[1]);
  double rho=fabs(1./par[2]);
  TVector3 center(cx,cy,0);
  double angle1=TMath::ACos( (center.Mag2()+rho*rho-r1*r1)/2/center.Mag()/rho );
  double angle2=TMath::ACos( (center.Mag2()+rho*rho-r2*r2)/2/center.Mag()/rho );
  double larc=TMath::Abs(rho*(angle1-angle2))*sqrt(1+par[4]*par[4]);
  return larc;
}
inline double MathTools::CalcHelixArc(const double par[5], const TVector3 &pos1,const TVector3 &pos2)
{
  /* double phi1 = par[2]/par[4]*(par[3]-pos1.Z()); */
  /* double phi2 = par[2]/par[4]*(par[3]-pos2.Z()); */
  double phi1=CalcHelixPhi(pos1,par);
  double phi2=CalcHelixPhi(pos2,par);
  double rho= fabs(1./par[2]);
  double larc=TMath::Abs(rho*(phi1-phi2))*sqrt(1+par[4]*par[4]);
  return larc;
}

inline TVector3 MathTools::CalcHelixGPosatR(const double par[5],const double gpar[4],const double &r)
{
  //
  // under development
  //
  double cx=(par[0]+1./par[2])*cos(par[1]+gpar[3])+gpar[0];
  double cy=(par[0]+1./par[2])*sin(par[1]+gpar[3])+gpar[1];
  double rho=fabs(1./par[2]);
  TVector3 center(cx,cy,0);
  double angle=TMath::ACos(( center.Mag2()+rho*rho-r*r)/2/center.Mag()/rho);
  TVector3 tmp=-center*(1./center.Mag())*rho;
  if(par[2]<0) tmp.RotateZ(angle);
  if(par[2]>0) tmp.RotateZ(-angle);
  double z;
  MathTools::CalcHelixXYtoZ(par,tmp.X(),tmp.Y(),z);
  tmp.SetZ(z+gpar[2]);
  return center+tmp;
}

inline TVector3 MathTools::CalcHelixMom(const double par[5], double z)
{
  //  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = -1*0.7; //T, "-1" is needed.
  double phi2 = (par[3]-z)*par[2]/par[4];
  double pt = fabs(1/par[2])*(Const*dMagneticField)/100.;
  double px = pt*(-1*sin(par[1]+phi2));
  double py = pt*(cos(par[1]+phi2));
  double pz = pt*(par[4]);
  return TVector3(px,py,pz);
}

inline double MathTools::CalcDeg(const double x,const  double y)
{
  double theta=atan( y/x )*TMath::RadToDeg();
  if(x>=0&&y>=0) theta=theta;
  else if(x<0&&y>=0) theta=180+theta;
  else if(x<0&&y<0) theta=180+theta;
  else if(x>=0&&y<0) theta=360+theta;
  return theta;
}
inline double MathTools::CalcRad(const double x,const  double y)
{
  double theta=atan( y/x );
  //  if(x>=0&&y>=0) theta=theta;
  if(x<0&&y>=0) theta=TMath::Pi()+theta;
  else if(x<0&&y<0) theta=TMath::Pi()+theta;
  else if(x>=0&&y<0) theta=2*TMath::Pi()+theta;
  return theta;
}

inline void MathTools::PointToLine( const TVector3 &p,
				    const TVector3 &x, const TVector3 &a,
				    double &dist,TVector3 &xest )
{
  double k=(p*a-x*a)/a.Mag2();
  xest=x+k*a;
  dist=(xest-p).Mag();
}

inline double MathTools::CalcHelixPhi(const TVector3 &pos, const double *par)
{
  return MathTools::CalcHelixPhi(pos.X(),pos.Y(),par);
}
inline double MathTools::CalcHelixPhi(const double &x,const double &y,const double *par)
{
  //  return value between -pi and pi 
  double sin_c,cos_c;
  double cx=(par[0]+1./par[2])*cos(par[1]);
  double cy=(par[0]+1./par[2])*sin(par[1]);

  cos_c=x-cx;  sin_c=y-cy;
  double phi; 
  phi=MathTools::CalcRad(cos_c,sin_c);
  //  std::cout<<"x,y  "<<x<<"  "<<y<<"  "<<phi<<std::endl;
  if(phi>TMath::Pi() ) phi-=2*TMath::Pi();  
  if(par[2]<0 ) phi-=par[1];
  else  phi-=(par[1]+TMath::Pi() );
  while(fabs(phi)>TMath::Pi())
    {
      if(phi>0) phi-=2*TMath::Pi();
      else phi+=2*TMath::Pi();
    }

  return phi;
}

inline double MathTools::dfunc_PTH(const TVector3 &pos,const double &helixphi,const double *par)
{
  double inv=1./par[2];
  double x  = par[0]*cos(par[1]) + inv*(cos(par[1]) - cos(par[1]+helixphi));
  double y  = par[0]*sin(par[1]) + inv*(sin(par[1]) - sin(par[1]+helixphi));
  double z  = par[3] - inv*par[4]*helixphi;
  double dx = inv*sin(par[1]+helixphi);
  double dy = -inv*cos(par[1]+helixphi);
  double dz = -inv*par[4];

  double S = 2.*( (x - pos.x())*dx + (y - pos.y())*dy + (z - pos.z())*dz);
  return S;
}

inline bool MathTools::PointToHelix(const TVector3 &hitpos, const double *par,TVector3 &fitpos ,double &dis)
{
  if(!PointToHelix(hitpos,fitpos,par)) return false;
  TVector3 tmp;
  tmp=fitpos-hitpos;
  dis=tmp.Mag();
  return true;  
}

inline bool MathTools::PointToHelix(const TVector3 &hitpos,TVector3 &fitpos, const double *par)
{

  double phi=0,phi_b=0,phi_a=0;
  phi=MathTools::CalcHelixPhi(hitpos.x(),hitpos.y(),par);

/*   int max=0; */
/*   while(max<4) */
/*     { */
/*       double z=par[3]-1./par[2]*par[4]*phi; */
/*       double z_b=par[3]-1./par[2]*par[4]*(phi-2*TMath::Pi()); */
/*       double z_a=par[3]-1./par[2]*par[4]*(phi+2*TMath::Pi()); */
/*       if(fabs(z-hitpos.z())<fabs(z_b-hitpos.z()) && fabs(z-hitpos.z())<fabs(z_a-hitpos.z())) break; */
/*       else if(fabs(z-hitpos.z())>fabs(z_b-hitpos.z())) phi-=2*TMath::Pi(); */
/*       else if(fabs(z-hitpos.z())>fabs(z_a-hitpos.z())) phi+=2*TMath::Pi(); */
/*       max++; */
/*     } */

  double philen=1./par[2]*sqrt(1+par[4]*par[4]);
  double dist=1.;
  int trial=14; //sigma_position<1.2 micron 
  while(dist<128)
    {
      phi_b=phi-dist/philen;
      phi_a=phi+dist/philen;
      double dlen_b=dfunc_PTH(hitpos,phi_b,par);
      double dlen_a=dfunc_PTH(hitpos,phi_a,par);
      if(dlen_b*dlen_a<=0) break;
      else {dist*=2;trial++;}
    }
  
  if(dist>=128) 
    {
#if DEBUG
      std::cout<<"Can not find PTH inital param!!"<<std::endl;
#endif
      fitpos.SetXYZ(999,999,999);
      return false;
    }
  
 //Bisection Method
  for(int i=0;i<trial;i++)
    {
      phi_b=phi-dist/philen;
      phi_a=phi+dist/philen;
      double dlen_b=dfunc_PTH(hitpos,phi_b,par);
      //      double dlen_a=dfunc_PTH(hitpos,phi_a,par);
      double dlen=dfunc_PTH(hitpos,phi,par);
      if(dlen*dlen_b<=0 ) {phi=(phi_b+phi)/2.0;dist=dist/2.0;}
      else  {phi=(phi_a+phi)/2.0;dist=dist/2.0;}      
    }
  //
  
  fitpos=MathTools::GetPosition(phi,par);
  return true;  
}

inline void MathTools::LineToLine( const TVector3 &x1, const TVector3 &a1,
				   const TVector3 &x2, const TVector3 &a2,
				   const double &dl,
				   double &dist,
				   TVector3 &xest, TVector3 &next )
{
  // X1 = x1 + t*a1
  // X2 = x2 + s*a2
  // (X1-X2)*a1=0
  // (X1-X2)*a2=0
  //    ||
  //    \/
  // a*t + b*s = A1
  // c*t + d*s = A2
  // xest : on x1,a1
  TVector3 x = x2-x1;   
  double a =  a1.Dot(a1); 
  double b = -a1.Dot(a2);
  double c =  a2.Dot(a1);
  double d = -a2.Dot(a2); 
  double A1 = a1.Dot(x); 
  double A2 = a2.Dot(x);
  
  double D = a*d-b*c;

  TVector3 x2p;
  if( fabs(D)<0.00000000000001 ){
    dist = sqrt(x.Mag2()-A1*A1);
  }
  else{
    double s = (a*A2-c*A1)/D;
    double t = (d*A1-b*A2)/D;
    xest  = x1 + t*a1;
    x2p   = x2 + s*a2;
    dist = (xest-x2p).Mag();
    next = x2p+(xest-x2p)*(dl/dist);
  }
#if 0
  std::cout << " dl:" << dl << " dist:" << dist
            << " dl2:" << (next-x2p).Mag()  << std::endl;
#endif
}

#endif
