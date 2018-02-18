// LinearTrack.cpp
#include <map>

#include "GlobalVariables.h"
#include "LinearTrack.h"
#include "TVector2.h"
#include "TMath.h"
#include "TMatrixD.h"

static const double SpatialResolutionOfBLDC=0.020; // [cm]

ClassImp(LinearTrack);

#define DEBUG 0
#define DEBUG2 0
// ----------------------------- //
// class LinearTrack              //
// ------------------------------//
LinearTrack::LinearTrack()
{
  A=B=C=D=-999.;
  GA=GB=GC=GD-999.;
  Dof=-999;
  Chi=-999.;
}

LinearTrack::LinearTrack( const LinearTrack &right )
{
  Clear();
  A = right.A;  B = right.B;  C = right.C;
  D = right.D;
  GA = right.GA;  GB = right.GB;  GC = right.GC;
  GD = right.GD; 
  Dof = right.Dof;
  Chi = right.Chi;
  
  for( vLinearTrackHitContainer::const_iterator it=right.LinearTrackHitContainer.begin(); it!=right.LinearTrackHitContainer.end(); it++ ){
    LinearTrackHitContainer.push_back( (*it) );
  }
}

LinearTrack::~LinearTrack()
{
  this->Clear();
}

void LinearTrack::DeleteHit( const int &i )
{
  vLinearTrackHitContainer::iterator it=LinearTrackHitContainer.begin();
  for( int j=0; j<i; j++ ) it++;
  LinearTrackHitContainer.erase(it);
}

void LinearTrack::Clear()
{
  LinearTrackHitContainer.clear();
}

bool LinearTrack::XYLocalPosatZ( const double &z, double &x, double &y, const bool &ROT )
{
  x = A + B*z ;
  y = C + D*z ;  
  if(ROT){
    double rot=0.;
    if( nhit()>0 )    rot = hit(0)->rot();    
    TVector2 org(x,y);
    TVector2 rotated=org.Rotate(-rot/180.*TMath::Pi());
    x=rotated.X();
    y=rotated.Y();
  }
  return true;
}

bool LinearTrack::XYPosatZ( const double &z, double &x, double &y )
{
  
  double gx=0,gy=0,gz=0;  
  if( nhit()>0 ){
    ChamberLikeHit *tmphit = hit(0);
    gx = tmphit->gx();  gy = tmphit->gy();  gz = tmphit->gz();
    //    std::cout<<"gx,gy,gz:\t"<<gx<<"\t"<<gy<<"\t"<<gz<<"\t"<<std::endl;
  }
  x = GA + GB*(z-gz) ;
  y = GC + GD*(z-gz) ;   
  
  //  x = GA + GB*z ;
  //  y = GC + GD*z ;   
  return true;
}

void LinearTrack::ConvLocalToGlobal()  
{
  double dgx=0,dgy=0,dgz=0;
  double gx=0,gy=0,gz=0;  
  if( nhit()>0 ){
    ChamberLikeHit *tmphit = hit(0);
    gx = tmphit->gx();  gy = tmphit->gy();  gz = tmphit->gz();
    dgx = tmphit->dgx();  dgy = tmphit->dgy();  dgz = tmphit->dgz();
    GA = A + gx;// + B*dgz;
    GB = B + dgx;
    GC = C + gy;// + D*dgz;
    GD = D + dgy;
  }
}

bool LinearTrack::Calc( ConfMan *conf )
{
  for( int i=0; i<nhit(); i++ )
    hit(i)->Calc(conf);    
  LinearFit(conf);
  CalcHitPosition();
  CalcResidual();
  ConvLocalToGlobal();
  return true;
}

static const unsigned int LRMASK=0x01;
static unsigned int LR( int hid, unsigned int key )
{
  key = key >> hid;
  key = key & LRMASK;
  return key;
}
//####################################################

bool LinearTrack::LinearFit( ConfMan *conf, const bool &CHECKLR )
{
#if 0
  std::cout << "!!! LinearTrack::LinearFit()" << std::endl;
#endif
  int MaxNumOfHitsInTrack=64;
  if(conf){
    BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
    MaxNumOfHitsInTrack =BLDCParam->GetMaxHitInTrack()*2;
  }
  if( nhit() > MaxNumOfHitsInTrack ){
    //    std::cout<<"!!! [LinearTrack::LinearFit()] too many hits !!!"<<std::endl;
    return false;
  }
  double x[MaxNumOfHitsInTrack], y[MaxNumOfHitsInTrack], theta[MaxNumOfHitsInTrack];
  int np = nhit();
  int np_org = np;
  int nlr=0;
  int ihitlr[MaxNumOfHitsInTrack];
  for( int i=0; i<np_org;i++){
    ihitlr[i]=-1;
    if( !layerstatus(hit(i)->cid(),hit(i)->layer()) ){
      np--; continue;
    }
    if(hit(i)->dl()<SpatialResolutionOfBLDC){
      ihitlr[i]=nlr; nlr++;
    }
  } 
  //  std::cout<<nlr <<" / "<<np<<std::endl;
  if( np < 2 ) return false;
  double minchi = 1.0e+9;
  double canda=0., candb=0., candc=0.,candd=0.;
  unsigned int candkey=0x0;
  unsigned int key = 0x0;
  if(!CHECKLR) nlr=0;
  int nconb = (int)pow(2,nlr);
  for( int ic=0; ic<nconb; ic++ ){      
    int j=0;
    for( int i=0; i<np_org;i++){
      if( !layerstatus(hit(i)->cid(),hit(i)->layer()) ) continue;
      double dl = hit(i)->dl();
      double temptheta;
      hit(i)->gwpos(y[j],temptheta,x[j],1);      
      //      std::cout<<"z,x,theta "<<x[j]<<"\t"<<y[j]<<"\t"<<temptheta<<std::endl; 
      theta[i]=temptheta/180.*TMath::Pi();      
      unsigned int lr =hit(i)->leftright();
      if(CHECKLR){
	if(ihitlr[i]>=0)
	  lr=LR(ihitlr[i],key);
      }
      if( lr==0 ) y[j] += dl;
      else y[j] -= dl;
      j++;
    }
    double right[4];
    double c2=0.,s2=0.,cs=0.,zc2=0.,zs2=0.,zcs=0.,z2c2=0.,z2s2=0.,z2cs=0;
    for(int i=0;i<4;i++){
      right[i]=0.;
    }
    for( int i=0; i<np; i++ ){
      c2 += pow(TMath::Cos(theta[i]),2);
      s2 += pow(TMath::Sin(theta[i]),2);
      cs += TMath::Cos(theta[i])*TMath::Sin(theta[i]);
      zc2 += x[i]*pow(TMath::Cos(theta[i]),2);
      zs2 += x[i]*pow(TMath::Sin(theta[i]),2);
      zcs += x[i]*TMath::Cos(theta[i])*TMath::Sin(theta[i]);
      z2c2 += x[i]*x[i]*pow(TMath::Cos(theta[i]),2);
      z2s2 += x[i]*x[i]*pow(TMath::Sin(theta[i]),2);
      z2cs += x[i]*x[i]*TMath::Cos(theta[i])*TMath::Sin(theta[i]);
      right[0]+=y[i]*TMath::Cos(theta[i]);
      right[1]+=y[i]*x[i]*TMath::Cos(theta[i]);
      right[2]+=y[i]*TMath::Sin(theta[i]);
      right[3]+=y[i]*x[i]*TMath::Sin(theta[i]);
    }
    TMatrixD mat(4,4);
    mat[0][0]= c2; mat[0][1]=zc2; mat[0][2]=cs; mat[0][3]=zcs; 
    mat[1][0]= zc2; mat[1][1]=z2c2; mat[1][2]=zcs; mat[1][3]=z2cs; 
    mat[2][0]= cs; mat[2][1]=zcs; mat[2][2]=s2; mat[2][3]=zs2; 
    mat[3][0]= zcs; mat[3][1]=z2cs; mat[3][2]=zs2; mat[3][3]=z2s2; 

    TMatrixD aa;
    aa.Use(4,1,right);
    mat.InvertFast();
    TMatrixD par(4,1);
    par.Mult(mat,aa);

    double chi = 0;
    double tempx,tempy,trackpos;
    for( int i=0; i<np; i++ ){
      tempx = par[0][0]+par[1][0]*x[i];
      tempy = par[2][0]+par[3][0]*x[i];
      trackpos = tempx*TMath::Cos(theta[i])+tempy*TMath::Sin(theta[i]);
      //      std::cout<<"i,trackpos,hitpos,chi2\t"<<i<<"\t"<<trackpos<<"\t"<<y[i]<<"\t"<<pow((y[i]-trackpos)/SpatialResolutionOfBLDC,2)<<std::endl;
      chi += pow((y[i]-trackpos)/SpatialResolutionOfBLDC,2);
    }
    if(np<5) chi=-1;
    else chi /= (double)(np-4);

    if( (chi < minchi) ){
      canda = par[0][0];
      candb = par[1][0];
      candc = par[2][0];
      candd = par[3][0];
      candkey=key;
      minchi=chi;	
    }
    key++;
  }
  SetABCD(canda,candb,candc,candd);
  SetDof(np-4);
  SetChisqr(minchi);
  int j=0;
  if(CHECKLR){
    for( int i=0; i<np_org; i++ ){
      if(ihitlr[i]>=0){
	unsigned int lr = LR(ihitlr[i],candkey);
	hit(i)->SetLeftRight((int)lr);
      }
      //    double dl = hit(i)->dl();
      //    double wpos = hit(i)->wpos(hit(i)->xy());
      //    if( lr==0 ) hit(i)->SetHitPos(hit(i)->xy(), wpos + dl );
      //    else hit(i)->SetHitPos(hit(i)->xy(), wpos - dl );
      j++;
    }  
  }
  return true;
}

void LinearTrack::SetGABCD(const double &a,const double &b, const double &c, const double &d )
{
  double dgx=0,dgy=0,dgz=0;
  double gx=0,gy=0,gz=0;  
  GA = a;
  GB = b;
  GC = c;
  GD = d;
  if( nhit()>0 ){
    ChamberLikeHit *tmphit = hit(0);
    gx = tmphit->gx();  gy = tmphit->gy();  gz = tmphit->gz();
    dgx = tmphit->dgx();  dgy = tmphit->dgy();  dgz = tmphit->dgz();
    //    std::cout<<"gx,gy,gz: "<<gx<<"\t"<<gy<<"\t"<<gz<<std::endl;
    A = a - gx ;//- b*dgz;
    B = b - dgx;
    C = c - gy ;//- d*dgz;
    D = d - dgy;
  }
}
void LinearTrack::CalcHitPosition( )
{
  for( int i=0; i<nhit(); i++ ){
    double x,y;//,z;
    double wx,theta,wz;
    hit(i)->gwpos(wx,theta,wz,1,0);
    XYLocalPosatZ(wz,x,y);
    hit(i)->SetHitPosition(x,y,wz);
  }
}

void LinearTrack::CalcResidual(bool PRINT )
{
  for( int i=0; i<nhit(); i++ ){
    ChamberLikeHit *tmphit = hit(i);
    double wx,theta,wz;
    tmphit->gwpos(wx,theta,wz,1,0);
    double x=tmphit->x();
    double y=tmphit->y();
    double trax=x*TMath::Cos(theta/180.*TMath::Pi())+y*TMath::Sin(theta/180.*TMath::Pi());
    double dltrack = trax-wx;
    double dl = tmphit->dl();
    int lr= tmphit->leftright();
    double resl;
    if(lr==0) resl = dltrack - dl;
    else resl = dltrack + dl;
    tmphit->SetResolution( resl );
    if(PRINT)
      std::cout<<"x,y,theta,trax,wx,dl,resl:   "<<x<<"\t"<<y<<"\t"<<theta<<"\t"<<trax<<"\t"<<wx<<"\t"<<dl<<"\t"<<resl<<std::endl;
  }
}


double LinearTrack::GetCalcChisquare( )
{
  double chisq=0;
  int ndf=0;
  for( int i=0; i<nhit(); i++ ){
    if(!layerstatus(hit(i)->cid(),hit(i)->layer()) ) continue;
    chisq += pow(hit(i)->resl()/SpatialResolutionOfBLDC,2);
    ndf++;
  }
  return chisq/(ndf-4);
}

