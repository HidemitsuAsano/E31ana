// HelixFit.cpp
#include <iostream>
#include "HelixFit.h"
#include "TMath.h"

#define DEBUG 0
ClassImp(HelixFit);


// parameters for TMinuit

/*
static const Double_t  FitStep[5] = { 0.005, 0.005, 1e-7,0.005,0.0001 };
static const Double_t LowLimit[5] = { -50, -2*TMath::Pi(), -0.5,-50,-10 };
static const Double_t  UpLimit[5] = { 50, 2*TMath::Pi(), 0.5,50,10 };
*/

static const Double_t  FitStep[5] = { 1e-20, 1e-21, 1e-22,1e-20,1e-20 };
static const Double_t LowLimit[5] = { -50, -2*TMath::Pi(), -0.5,-50,-10 };
static const Double_t  UpLimit[5] = { 50, 2*TMath::Pi(), 0.5,50,10 };

// static const Double_t  FitStep[3] = { 0.001, 0.0001, 0.001 };
// static const Double_t LowLimit[3] = { -4000, -TMath::Pi(), -2000 };
// static const Double_t  UpLimit[3] = { 4000, TMath::Pi(), 2000 };

// global variables for TMinuit
static Int_t gNumOfHits;
static TVector3 gHitPos[MAX_NUM_OF_HITS_H];
static Double_t gWeight[MAX_NUM_OF_HITS_H];


// --------------------------------------------------------------//
// functions for TMinuit
static void NearestPos( const Double_t *par, const Double_t &hitx, 
			const Double_t &hity,
			Double_t &fitx, Double_t &fity)
{
  Double_t cx, cy;
  cx = ( par[0] + 1./par[2] )*cos( par[1] );
  cy = ( par[0] + 1./par[2] )*sin( par[1] );

  Double_t nx, ny, nn;
  nx = hitx - cx;
  ny = hity - cy;
  nn = sqrt( nx*nx + ny*ny );
  nx /= nn;
  ny /= nn;
  
  fitx = cx + fabs( 1./par[2] )*nx;
  fity = cy + fabs( 1./par[2] )*ny;
}



static double CalcHelixPhi(const double &x,const double &y,const double *par)
{
  
  double sin_c=0,cos_c=0;
  double cx=(par[0]+1./par[2])*cos(par[1]);
  double cy=(par[0]+1./par[2])*sin(par[1]);

  cos_c=x-cx;  sin_c=y-cy;
  double phi; 
  MathTools *tool=new MathTools();
  phi=tool->CalcDeg(cos_c,sin_c)*TMath::DegToRad();
  delete tool;
  if(phi>TMath::Pi() ) phi-=2*TMath::Pi();
  
  if(par[2]<0 ) 
    {phi-=par[1];}
  else 
    {phi-=(par[1]+TMath::Pi() );}

   
  while(fabs(phi)>TMath::Pi())
    {
      if(phi>0) phi-=2*TMath::Pi();
      else phi+=2*TMath::Pi();
    }
  
  return phi;
}


static TVector3 GetPosition(const double &helixphi,const double *par)
{
  TVector3 pos;
  pos.SetXYZ(par[0]*cos(par[1])+1./par[2]*( cos(par[1])-cos(par[1]+helixphi) ),
	     par[0]*sin(par[1])+1./par[2]*( sin(par[1])-sin(par[1]+helixphi) ),
	     par[3]-1./par[2]*par[4]*helixphi);

  return pos;
}


static double dfunc_PTH(const TVector3 &pos,const double &helixphi,const double *par)
{
  double x  = par[0]*cos(par[1]) + (1./par[2])*(cos(par[1]) - cos(par[1]+helixphi));
  double y  = par[0]*sin(par[1]) + (1./par[2])*(sin(par[1]) - sin(par[1]+helixphi));
  double z  = par[3] - (1./par[2])*par[4]*helixphi;
  double dx = (1./par[2])*sin(par[1]+helixphi);
  double dy = (-1./par[2])*cos(par[1]+helixphi);
  double dz = (-1./par[2])*par[4];

  double S = 2.*( (x - pos.x())*dx + (y - pos.y())*dy + (z - pos.z())*dz);
  return S;
}

static bool NearestPosForStereo( const Double_t *par, const TVector3 &hitpos, 
		        TVector3 &fitpos )
{
  double dis=999;
  double phi=0,phi_b=0,phi_a=0;
  phi=CalcHelixPhi(hitpos.x(),hitpos.y(),par);

  /*
  int max=0; 
  while(max<4)
    {
      double z=par[3]-1./par[2]*par[4]*phi;
      double z_b=par[3]-1./par[2]*par[4]*(phi-2*TMath::Pi());
      double z_a=par[3]-1./par[2]*par[4]*(phi+2*TMath::Pi());
      if(fabs(z-hitpos.z())<fabs(z_b-hitpos.z()) && fabs(z-hitpos.z())<fabs(z_a-hitpos.z())) break;
      else if(fabs(z-hitpos.z())>fabs(z_b-hitpos.z())) phi-=2*TMath::Pi();
      else if(fabs(z-hitpos.z())>fabs(z_a-hitpos.z())) phi+=2*TMath::Pi();
      max++;
    }
  */
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
  
  if(dis>=128) 
    {
#if DEBUG
      std::cout<<"Can not find PTH in FCN2!!"<<std::endl;
#endif
    fitpos.SetXYZ(999,999,999);
    return false;}
 
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

  fitpos=GetPosition(phi,par);

  TVector3 tmp;
  tmp=fitpos-hitpos;
  dis=tmp.Mag();
  return true;
  
}



/*
static void fcn( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  Double_t chisq=0.;
  Int_t dof = 0;

  for( Int_t i=0; i<gNumOfHits; i++ ){
    Double_t hitx, hity, fitx, fity;
    hitx = gHitPos[i].x();
    hity = gHitPos[i].y();
    NearestPos(par, hitx, hity, fitx, fity);
    chisq +=( (hitx-fitx)*(hitx-fitx) + (hity-fity)*(hity-fity) )/gWeight[i]/gWeight[i];
    dof++;
  }

  f = chisq/(dof-3);
}
*/
static void fcn2( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  Double_t chisq=0.;
  Int_t dof = 0;

  for( Int_t i=0; i<gNumOfHits; i++ ){
    Double_t hitx, hity, fitx, fity;
    if(gHitPos[i].z()==-999)
      {
	hitx = gHitPos[i].x();
	hity = gHitPos[i].y();
	NearestPos(par, hitx, hity, fitx, fity);
	chisq +=( (hitx-fitx)*(hitx-fitx) + (hity-fity)*(hity-fity) )/gWeight[i]/gWeight[i];
      }
    else
      {
	TVector3 fittmp;
	if(!NearestPosForStereo(par, gHitPos[i],fittmp) ) 
	  {	
	    //    std::cout<<"Miss Newton Method in fcn!!"<<std::endl;
	    chisq+=9999;
	    continue;
	  }
	fittmp -=gHitPos[i];
	chisq +=fittmp.Mag2()/gWeight[i]/gWeight[i];

      }
    dof++;
  }

  f = chisq/(dof-5);
}

// --------------------------------------------------------------//


HelixFit::HelixFit()
{
  for( int i=0; i<5; i++ ){
    Par[i] = Err[i] = -999.;
  }

  NumOfHits = FitDof = FitStat = -999;
  FitChi2 = -999.;

  for( int i=0; i<MAX_NUM_OF_HITS_H; i++ ){
    HitPos[i].SetXYZ(-999,-999,-999);
    Weight[i] = -999.;
  }
  
  minuit = new TMinuit(5);
  TROOT minexam("HelixFit","helix fit using TMinuit");
}

HelixFit::HelixFit( const double *initPar, const double *hitx,
		    const double *hity,const double *hitz,
		    const double *weight, const int &numofhit)
{
  NumOfHits = (Int_t)numofhit;

  for( int i=0; i<5; i++ ) {
    Par[i] = (Double_t)initPar[i];
    Err[i] = -999.;
  }

  for( int i=0; i<MAX_NUM_OF_HITS_H; i++ ) {
    if( i<NumOfHits ){
      HitPos[i].SetXYZ( hitx[i],hity[i],hitz[i]);
      Weight[i] = (Double_t)weight[i];
      
      if(hitz[i]!=-999)
	{
	  TVector3 FitPos;
	  CalcNearestPosForStereo(HitPos[i],FitPos);
	  //	  std::cout<<"HitPos= "<<HitPos[i].x()<<" "<<HitPos[i].y()<<" "<<HitPos[i].z()<<std::endl
	  //  <<"FitPos= "<<FitPos.x()<<" "<<FitPos.y()<<" "<<FitPos.z()<<std::endl<<std::endl;  
	}
      
    }
    else {
      HitPos[i].SetXYZ(-999,-999,-999);
      Weight[i] = -999.;
    }
  }

  FitDof = FitStat = -999;
  FitChi2 = -999.;

  minuit = new TMinuit(5);
  TROOT minexam("HelixFit","Helix fit using TMinuit");

  fit();
}

HelixFit::~HelixFit()
{
  delete minuit;
}

void HelixFit::SetParameters( const double *param )
{
  for( int i=0; i<5; i++ ) Par[i] = (Double_t)param[i];
}

void HelixFit::SetHitPos( const TVector3 hitpos,  const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )
          HitPos[i]=hitpos;
}

void HelixFit::SetWeight( const double *weight, const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )
    Weight[i] = (Double_t)weight[i];
}

void HelixFit::GetParameters( double *param )
{
  for( int i=0; i<5; i++ ) param[i] = (double)Par[i];
}


void HelixFit::SetGlobalVariables()
{
  gNumOfHits = NumOfHits;
  for( int i=0; i<MAX_NUM_OF_HITS_H; i++ ){
    gHitPos[i] = HitPos[i];
    gWeight[i] = Weight[i];
  }
}

void HelixFit::fit()
{
#if 0
  std::cout << "!!! HelixFit::fit() !!!" << std::endl;
#endif

#if DEBUG
  Int_t plevel=1;
#else
  Int_t plevel=-1;
#endif

  SetGlobalVariables();

  minuit->SetPrintLevel( plevel );
  minuit->SetFCN( fcn2 );

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  //  minuit->mnexcm("SET ERR", arglist,1,ierflg);
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  // Set starting values and step sizes for parameters
  TString name[5] = {"d_rho", "phi_0", "rho","d_z","tanL"};
  for( Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i],Par[i],FitStep[i],LowLimit[i],UpLimit[i],ierflg);
  }
  minuit->Command("SET STRategy 0");
  // Now ready for minimization step
  arglist[0] = 1000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  //minuit->Command("TMProve 100");

  // Print results
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm,  errdef, nvpar, nparx, icstat);
#if DEBUG
  minuit->mnprin(5,amin);
#endif
  
  Int_t err;
  Double_t bnd1, bnd2;
  for( Int_t i=0; i<5; i++ ){
    minuit->mnpout(i, name[i], Par[i], Err[i], bnd1, bnd2, err);
  }
  
  FitStat = icstat;
  CalcChi2();
}

void HelixFit::CalcNearestPos( const Double_t &hitx, const Double_t &hity,
				  Double_t &fitx, Double_t &fity)
{
  Double_t cx, cy;
  cx = ( Par[0] + 1./Par[2] )*cos( Par[1] );
  cy = ( Par[0] + 1./Par[2] )*sin( Par[1] );

  Double_t nx, ny, nn;
  nx = hitx - cx;
  ny = hity - cy;
  nn = sqrt( nx*nx + ny*ny );
  nx /= nn;
  ny /= nn;
  
  fitx = cx + fabs( 1./Par[2] )*nx;
  fity = cy + fabs( 1./Par[2] )*ny;
}

void HelixFit::CalcNearestPosForStereo( const TVector3 &hitpos,  TVector3 &fitpos )
{

  Double_t cx, cy;
  cx = ( Par[0] + 1./Par[2] )*cos( Par[1] );
  cy = ( Par[0] + 1./Par[2] )*sin( Par[1] );

  Double_t nx, ny, nn;
  nx = hitpos.x() - cx;
  ny = hitpos.y() - cy;
  nn = sqrt( nx*nx + ny*ny );
  nx /= nn;
  ny /= nn;

  double phi=0,phi2=0;
  
  if( nx>=0 && ny>=0) phi=atan(ny/nx);
  else if( nx< 0 && ny >=0) phi=TMath::Pi()+atan(ny/nx);
  else if( nx<0 && ny<0) phi=TMath::Pi()+atan(ny/nx);
  else if( nx>=0 && ny<0) phi=2*TMath::Pi()+atan(ny/nx);

  if(1./Par[2]>0) phi=phi-Par[1]-TMath::Pi();  
  else phi=phi-Par[1];

  while(fabs(phi)>2*TMath::Pi()-1.0)
    {
      if(phi>0) phi-=2*TMath::Pi();
      else phi+=2*TMath::Pi();
    }

  fitpos.SetX( cx - 1./Par[2]*cos(Par[1]+phi) );
  fitpos.SetY( cy - 1./Par[2]*sin(Par[1]+phi) );
  fitpos.SetZ( Par[3] - 1./Par[2] *Par[4]*phi );
  phi2=phi;
  //  std::cout<<"phi "<<phi<<std::endl;

  double eps=1.0e-8;
  int numnewton=0;
  int numnewton2=0;
  bool flag=false;
  double aphi[2]={phi,phi};
  while(numnewton2<30)
    {
      aphi[0]=aphi[1];
      fitpos.SetX( cx - 1./Par[2]*cos(Par[1]+aphi[0]) );
      fitpos.SetY( cy - 1./Par[2]*sin(Par[1]+aphi[0]) );
      fitpos.SetZ( Par[3] - 1./Par[2] *Par[4]*aphi[0] );
      
      TVector3 dfitpos,ddfitpos;
      
      dfitpos.SetX(1./Par[2]*sin(Par[1]+aphi[0]));
      dfitpos.SetY(-1./Par[2]*cos(Par[1]+aphi[0]));
      dfitpos.SetZ(-1./Par[2]*Par[4] );
            
      ddfitpos.SetX(1./Par[2]*cos(Par[1]+aphi[0]));
      ddfitpos.SetY(1./Par[2]*sin(Par[1]+aphi[0]));
      ddfitpos.SetZ( 0 );

      double func=fitpos*dfitpos-hitpos*dfitpos;
      double dfunc=dfitpos.Mag2()+fitpos*ddfitpos-hitpos*ddfitpos;

      //      if(dfunc==0){ std::cout<<"Can't Calc Newton method!!"<<std::endl; return;}
      if(dfunc==0){ aphi[1]+=0.001;continue;}
      aphi[1]=aphi[0]-func/dfunc;
      if(fabs(aphi[0]-aphi[1])<eps) {flag=true; break;}
      numnewton++;
      if(numnewton>30) {numnewton=0;numnewton2++;aphi[1]=phi-.15+numnewton2*0.1;}
    }
  if(flag)
    {  
      fitpos.SetX( (Par[0]+1./Par[2] )*cos(Par[1])-1./Par[2]*cos(Par[1]+aphi[1]) );
      fitpos.SetY( (Par[0]+1./Par[2] )*sin(Par[1])-1./Par[2]*sin(Par[1]+aphi[1]));
      fitpos.SetZ( Par[3] - 1./Par[2] *Par[4]*aphi[1] );
      //      std::cout<<"PtH numnewton "<<numnewton<<std::endl;
    }
  else
    {
      fitpos.SetX( -999 );
      fitpos.SetY( -999 );
      fitpos.SetZ( -999 );
      /*
      //   std::cout<<"Too many run Newton method PtH"<<std::endl;
      const int Trial=1000;
      const double dphi=TMath::Pi()/6.;
      double dis;
      for(int i=0;i<Trial;i++)
	{
	  double phitmp=phi-dphi/2.+dphi/Trial*i;
	  TVector3 fittmp;
	  fittmp.SetX( (Par[0]+1./Par[2] )*cos(Par[1])-1./Par[2]*cos(Par[1]+phitmp) );
	  fittmp.SetY( (Par[0]+1./Par[2] )*sin(Par[1])-1./Par[2]*sin(Par[1]+phitmp));
	  fittmp.SetZ( Par[3]- 1./Par[2] *Par[4]*phitmp );

	  fittmp-=hitpos;
	  double distmp=fittmp.Mag();
	  if(i==0) {phi2=phitmp;dis=distmp;}
	  else if(distmp<dis) {dis=distmp; phi2=phitmp;}
	}
      
      fitpos.SetX( (Par[0]+1./Par[2] )*cos(Par[1])-1./Par[2]*cos(Par[1]+phi2) );
      fitpos.SetY( (Par[0]+1./Par[2] )*sin(Par[1])-1./Par[2]*sin(Par[1]+phi2));
      fitpos.SetZ( Par[3]- 1./Par[2] *Par[4]*phi2 );
      */
    }
  
  //  return true;
}



void HelixFit::CalcChi2()
{  
  Double_t chisq=0.;
  Int_t dof = 0;

  for( Int_t i=0; i<NumOfHits; i++ ){
    Double_t hitx, hity, fitx, fity;

    if(HitPos[i].z()==-999)
      {
	hitx = HitPos[i].x();
	hity = HitPos[i].y();
	CalcNearestPos( hitx, hity, fitx, fity);
	chisq +=( (hitx-fitx)*(hitx-fitx) + (hity-fity)*(hity-fity) )/Weight[i];
      }
    else if(fabs( HitPos[i].z() )<50.)
      {
	TVector3 fittmp;
	CalcNearestPosForStereo( HitPos[i],fittmp);
	fittmp -=HitPos[i];
	chisq +=fittmp.Mag2()/Weight[i]/Weight[i];
      }
    dof++;
  }

  FitDof = dof - 5;
  FitChi2 = chisq;
}

