// BeamLineTrackMan.cpp
#include <map>

#include "GlobalVariables.h"
#include "BeamLineTrackMan.h"
#include "TVector2.h"
#include "TMath.h"

static const double SpatialResolutionOfBLDC=0.02; // [cm]

ClassImp(BeamLineTrackMan);
ClassImp(BLDCClusterMan);
ClassImp(BLDCCluster);
ClassImp(LocalTrack);

#define DEBUG 0
#define DEBUG2 0
#define CHECKRANGE 1
#define SLOPECUT 1
#define DEL 1
#define KILLLAYER 1
// ----------------------------- //
// class LocalTrack              //
// ------------------------------//
LocalTrack::LocalTrack()
{
  A=B=C=D=E=F=-999.;
  GA=GB=GC=GD=GE=GF=-999.;
  xzDof=yzDof=-999;
  xzChi=yzChi=-999.;
  Vtx1=Vty1=Vtz1=Vtx2=Vty2=Vtz2=-999.;
}

LocalTrack::LocalTrack( const LocalTrack &right )
{
  A = right.A;  B = right.B;  C = right.C;
  D = right.D;  E = right.E;  F = right.F;
  GA = right.GA;  GB = right.GB;  GC = right.GC;
  GD = right.GD;  GE = right.GE;  GF = right.GF;
  xzDof = right.xzDof;  yzDof = right.yzDof;
  xzChi = right.xzChi;  yzChi = right.yzChi;
  Vtx1 = right.Vtx1;  Vty1 = right.Vty1;  Vtz1 = right.Vtz1;
  Vtx2 = right.Vtx2;  Vty2 = right.Vty2;  Vtz2 = right.Vtz2;
  
  for( LocalTrackHitContainer::const_iterator it=right.xzLocalTrackHitContainer.begin(); it!=right.xzLocalTrackHitContainer.end(); it++ ){
    xzLocalTrackHitContainer.push_back( (*it) );
  }

  for( LocalTrackHitContainer::const_iterator it=right.yzLocalTrackHitContainer.begin(); it!=right.yzLocalTrackHitContainer.end(); it++ ){
    yzLocalTrackHitContainer.push_back( (*it) );
  }
}

LocalTrack::~LocalTrack()
{
  this->Clear();
}

void LocalTrack::DeleteHitXZ( const int &i )
{
  LocalTrackHitContainer::iterator it=xzLocalTrackHitContainer.begin();
  for( int j=0; j<i; j++ ) it++;
  xzLocalTrackHitContainer.erase(it);
}

void LocalTrack::DeleteHitYZ( const int &i )
{
  LocalTrackHitContainer::iterator it=yzLocalTrackHitContainer.begin();
  for( int j=0; j<i; j++ ) it++;
  yzLocalTrackHitContainer.erase(it);
}

ChamberLikeHit* LocalTrack::hit( const int &i )
{
  return ( i<(this->nhitxz()) ? (this->hitxz(i)) : (this->hityz( i - (this->nhitxz()) )) );
}

void LocalTrack::Clear()
{
  xzLocalTrackHitContainer.clear();
  yzLocalTrackHitContainer.clear();
}

void LocalTrack::DeleteHit(const int &xy, const int &i )
{
  if(xy==0) DeleteHitXZ(i);
  else if(xy==1) DeleteHitYZ(i);
}

bool LocalTrack::XYLocalPosatZ( const double &z, double &x, double &y )
{
  if( fabs(A) < 1.0e-9 ) return false;
  if( fabs(D) < 1.0e-9 ) return false;
  
  x = -1.*( B*z + C )/A;
  y = -1.*( E*z + F )/D;
  
  return true;
}

bool LocalTrack::XYSemiLocalPosatZ( const double &z, double &x, double &y )
{
  if( fabs(A) < 1.0e-9 ) return false;
  if( fabs(D) < 1.0e-9 ) return false;
  double rotx=0.,roty=0.;
  int cid=0;
  double tiltx=0.,tilty=0.;
  if( nhitxz()>0 ){
    ChamberLikeHit *hit = hitxz(0);
    rotx = hit->rot()/180.*TMath::Pi();
    cid= hit->cid();
    tiltx=hit->tilt()/180.*TMath::Pi();;
  }
  if( nhityz()>0 ){
    ChamberLikeHit *hit = hityz(0);
    roty = hit->rot()/180.*TMath::Pi();
    cid= hit->cid();
    tilty=hit->tilt()/180.*TMath::Pi();;
  }
  if(cid==CID_BLC2a||cid==CID_BLC2b){
    rotx-=135/180.*TMath::Pi();
    roty-=135/180.*TMath::Pi();
  }
  if(cid==CID_BLC1a||cid==CID_BLC1b){
    rotx-=45/180.*TMath::Pi();
    roty-=45/180.*TMath::Pi();
  }
  x = -1.*( (B-tiltx)*z + C )/A;
  y = -1.*( (E-tilty)*z + F )/D;

  TVector2 pos(x,y);
  TVector2 dir(B/A-tiltx,E/D-tilty);
  TVector2 pos1=pos.Rotate(rotx);
  TVector2 pos2=pos.Rotate(roty);
  x = pos1.X();//-gx;
  y = pos2.Y();//-gy;

  return true;  
}

void LocalTrack::semiabcdef(double &a,double &b,double &c,double &d,double &e,double &f)
{
  if( fabs(A) < 1.0e-9 ) return ;
  if( fabs(D) < 1.0e-9 ) return ;
  double rotx=0.,roty=0.;
  int cid=0;
  double tiltx=0.,tilty=0.;
  if( nhitxz()>0 ){
    ChamberLikeHit *hit = hitxz(0);
    rotx = hit->rot()/180.*TMath::Pi();
    cid= hit->cid();
    tiltx=hit->tilt()/180.*TMath::Pi();;
  }
  if( nhityz()>0 ){
    ChamberLikeHit *hit = hityz(0);
    roty = hit->rot()/180.*TMath::Pi();
    cid= hit->cid();
    tilty=hit->tilt()/180.*TMath::Pi();;
  }
  if(cid==CID_BLC2a||cid==CID_BLC2b){
    rotx-=135/180.*TMath::Pi();
    roty-=135/180.*TMath::Pi();
  }
  if(cid==CID_BLC1a||cid==CID_BLC1b){
    rotx-=45/180.*TMath::Pi();
    roty-=45/180.*TMath::Pi();
  }

  TVector2 pos(C/A,F/D);
  TVector2 dir(B/A-tiltx,E/D-tilty);
  TVector2 pos1=pos.Rotate(rotx);
  TVector2 pos2=pos.Rotate(roty);
  TVector2 dir1=dir.Rotate(rotx);
  TVector2 dir2=dir.Rotate(roty);
  a=1.;
  b= dir1.X();
  c= pos1.X();//-gx;
  d= 1.;
  e= dir2.Y();
  f = pos2.Y();//-gy;
}
bool LocalTrack::XYPosatZ( const double &z, double &x, double &y )
{
  if( fabs(GA) < 1.0e-9 ) return false;
  if( fabs(GD) < 1.0e-9 ) return false;

  double gx=0,gy=0,gz=0;  
  if( nhitxz()>0 ){
    ChamberLikeHit *hit = hitxz(0);
    gx = hit->gx();  gy = hit->gy();  gz = hit->gz();
    //std::cout<<"gx,gy,gz: "<<gx<<"\t"<<gy<<"\t"<<gz<<std::endl;
  }

  x = -1.*( GB*(z-gz) + GC )/GA + gx ;
  y = -1.*( GE*(z-gz) + GF )/GD + gy ;   
  //  x = -1.*( GB*z + GC )/GA;
  //  y = -1.*( GE*z + GF )/GD;
  
  return true;
}

void LocalTrack::ConvLocalToGlobal()  
{
  //  double gx, gy, gz;
  double rot;
  //  std::cout<<"convL2G nhit "<<nhitxz()<<std::endl;
  double dgx=0,dgy=0,dgz=0;
  double tilt=0; 
  if( nhitxz()>0 ){
    ChamberLikeHit *hit = hitxz(0);
    rot = hit->rot()/180.*TMath::Pi();
    dgx = hit->dgx();  dgy = hit->dgy();  dgz = hit->dgz();
    tilt= hit->tilt()/180.*TMath::Pi();;
    TVector2 pos(C/A,F/D);
    TVector2 dir(B/A-tilt,E/D);
    TVector2 pos2=pos.Rotate(rot);
    TVector2 dir2=dir.Rotate(rot);
    GA = 1.; 
    GB = dir2.X()-dgx;
    GC = pos2.X();//-gx;
  }
  if( nhityz()>0 ){
    ChamberLikeHit *hit = hityz(0);
    rot = hit->rot()/180.*TMath::Pi();
    dgx = hit->dgx();  dgy = hit->dgy();  dgz = hit->dgz();
    tilt= hit->tilt()/180.*TMath::Pi();;
    TVector2 pos(C/A,F/D);
    TVector2 dir(B/A,E/D-tilt);
    TVector2 pos2=pos.Rotate(rot);
    TVector2 dir2=dir.Rotate(rot);
    GD = 1.; 
    GE = dir2.Y()-dgy;
    GF = pos2.Y();//-gx;
  }
  //  std::cout<<"check convL2G abc gagbgc: "<<A<<" "<<B<<" "<<C<<" ,"<<
  //    GA<<" "<<GB<<" "<<GC<<std::endl;
}

bool LocalTrack::Calc( ConfMan *conf )
{
  for(int xy=0;xy<2;xy++){
    for( int j=0; j<nhit(xy); j++ )
      hit(xy,j)->Calc(conf);    
    LeastSquareFit(conf, xy);
  }
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
bool LocalTrack::PreTracking( ConfMan *conf , const int &xy )
{
#if 0
  std::cout << "!!! LocalTrack::PreTracking()" << std::endl;
#endif
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  int MaxNumOfHitsInTrack=BLDCParam->GetMaxHitInTrack();
  
  if( nhit() > MaxNumOfHitsInTrack*2 ){
    //    std::cout<<"!!! [LinearTrack::LinearFit()] too many hits !!!"<<std::endl;
    return false;
  }
  double x[MaxNumOfHitsInTrack*2], y[MaxNumOfHitsInTrack*2];//, theta[MaxNumOfHitsInTrack*2];
  int np = nhit();
  int np_org = np;
  //  int nlr=0;
  //  int ihitlr[MaxNumOfHitsInTrack*2];
  for( int i=0; i<np_org;i++){
    if(BLDCParam->layerkilled(hit(xy,i)->cid(),hit(xy,i)->layer())){
      np--; continue;    
    }
    //    if( !LAYER[8*((hit(xy,i)->cid()-1)%2)+hit(xy,i)->layer()-1] ) {
    if( !layerstatus(hit(xy,i)->cid(),hit(xy,i)->layer()) ) {
      np--; continue;    
    }
#if 0
    std::cout<<"cid,layer,wire\t"<<
      hit(xy,i)->cid()<<"\t"<<hit(xy,i)->layer()<<"\t"<<hit(xy,i)->wire()<<std::endl;
#endif     
  } 
  if( np < 2 ) return false;

  int j=0;
  for( int i=0; i<np_org;i++){
    if(BLDCParam->layerkilled(hit(xy,i)->cid(),
			      hit(xy,i)->layer())) continue;
    if( !layerstatus(hit(xy,i)->cid(),hit(xy,i)->layer()) ) continue;    
    x[j]=hit(xy,i)->wz();
    y[j]=hit(xy,i)->wpos(xy);      
    j++;
  }
  if(j!=np) return false;
  double Sx, Sy, Sxy, Sxx;
  Sx=Sy=Sxy=Sxx=0;
  for( int i=0; i<np; i++ ){
    Sx  += x[i];
    Sy  += y[i];
    Sxy += x[i]*y[i];
    Sxx += x[i]*x[i];
  }
  double D = (double)np*Sxx - Sx*Sx;
  if( D==0 ) return false;
  double aa = (Sxx*Sy - Sx*Sxy)/D;
  double bb = ((double)np*Sxy - Sx*Sy)/D;
  double chi = 0;
  double maxdl=hit(xy,0)->dxy()/2.;
  for( int i=0; i<np; i++ ) {
    double diff=pow(y[i]-(aa+bb*x[i]),2)/(maxdl*maxdl/3);
    //    std::cout<<"x,y, calcy, maxdl, diff  "<<x[i]<<"\t"<<y[i]<<"\t"<<aa+bb*x[i]<<"\t"<<maxdl<<"\t"<<diff<<std::endl;
    chi += diff;
  }
  if(np<3) chi=-1;
  else chi /= (double)(np-2);
  if( chi > BLDCParam->GetMaxChiPreTracking()  ) 
    return false;
  if( TMath::Abs(bb) > BLDCParam->GetMaxSlope()  ) 
    return false;
#if 0
  for( int i=0; i<np_org;i++){
    if( !layerstatus(hit(xy,i)->cid(),hit(xy,i)->layer()) ) continue;    
    std::cout<<"cid,layer,wire\t"<<
      hit(xy,i)->cid()<<"\t"<<hit(xy,i)->layer()<<"\t"<<hit(xy,i)->wire()<<std::endl;
  } 
  std::cout<<"chisquare=  "<<chi<<std::endl;
  std::cout<<"slope=  "<<TMath::Abs(bb)<<std::endl;
#endif     
  if(xy)
    SetDEF(1,bb,aa);
  else
    SetABC(1,bb,aa);
  SetDof(xy,np-2);
  SetChisqr(xy,chi);
  return true;
}

bool LocalTrack::LeastSquareFit( ConfMan *conf, const int &xy )
{
#if 0
  std::cout << "!!! LocalTrack::LeastSquareFit()" << std::endl;
#endif
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  int MaxNumOfHitsInTrack=BLDCParam->GetMaxHitInTrack();

  if( nhit(xy) > MaxNumOfHitsInTrack ) return false;
  double x[MaxNumOfHitsInTrack], y[MaxNumOfHitsInTrack];//, theta[MaxNumOfHitsInTrack];
  int np = nhit(xy);
  int np_org = np;
  for( int i=0; i<np_org;i++){
    if(BLDCParam->layerkilled(hit(xy,i)->cid(),
			      hit(xy,i)->layer())){
      np--; continue;
    }    
    if( !layerstatus(hit(xy,i)->cid(),hit(xy,i)->layer()) ){
      np--; continue;
    }
  }
  if( np < 2 ) return false;
  double minchi = 1.0e+9;
  double canda=0., candb=0., candc=0.;
  unsigned int candkey=0x0;
  unsigned int key = 0x0;
  int nconb = (int)pow(2,np);
  
  if(xy==0||xy==1){
    for( int ic=0; ic<nconb; ic++ ){      
      int j=0;
      for( int i=0; i<np_org;i++){
	if(BLDCParam->layerkilled(hit(xy,i)->cid(),
				  hit(xy,i)->layer())) continue;
	if( !layerstatus(hit(xy,i)->cid(),hit(xy,i)->layer()) ) continue;
	 double dl = hit(xy,i)->dl();
	 x[j] = hit(xy,i)->wz();
	 unsigned int lr = LR(j,key);
	 if( lr==0 ) y[j] = (hit(xy,i)->wpos(xy)) + dl;
	 else y[j] = (hit(xy,i)->wpos(xy)) - dl;
	 j++;
       }
       double Sx, Sy, Sxy, Sxx;
       Sx=Sy=Sxy=Sxx=0;
       for( int i=0; i<np; i++ ){
	 Sx  += x[i];
	 Sy  += y[i];
	 Sxy += x[i]*y[i];
	 Sxx += x[i]*x[i];
       }
       double D = (double)np*Sxx - Sx*Sx;
       if( D==0 ) return false;
       double aa = (Sxx*Sy - Sx*Sxy)/D;
       double bb = ((double)np*Sxy - Sx*Sy)/D;
       double chi = 0;
       for( int i=0; i<np; i++ ) chi += (y[i]-aa-bb*x[i])*(y[i]-aa-bb*x[i])/(SpatialResolutionOfBLDC*SpatialResolutionOfBLDC);

       if(np<3) chi=-1;
       else chi /= (double)(np-2);
#if SLOPECUT
       if( TMath::Abs(bb) < 0.1  ) 
#endif
	 if( (chi < minchi) || (chi==minchi && candb*candb<bb*bb) ){
	   canda = 1.;
	   candb = -1.*bb;
	   candc = -1.*aa;
	   candkey=key;
	   minchi=chi;	
	 }
       key++;
    }
    if(xy)
      SetDEF(canda,candb,candc);
    else
      SetABC(canda,candb,candc);
    SetDof(xy,np-2);
    SetChisqr(xy,minchi);
    int j=0;
    for( int i=0; i<np_org; i++ ){
      double dl = hit(xy,i)->dl();
      double wpos = hit(xy,i)->wpos(xy);
      if(BLDCParam->layerkilled(hit(xy,i)->cid(),
				hit(xy,i)->layer())||
	 !layerstatus(hit(xy,i)->cid(),hit(xy,i)->layer()) ){
	double x,y,z;
	z = hit(xy,i)->wz();
	XYLocalPosatZ(z,x,y);
	unsigned int lr=0;
	if(fabs(x-(wpos+dl))>fabs(x-(wpos-dl))) lr=1;
	if( lr==0 ) hit(xy,i)->SetHitPos( xy, wpos + dl );
	else hit(xy,i)->SetHitPos( xy, wpos - dl );
	hit(xy,i)->SetLeftRight((int)lr);	
	continue;
      }
      unsigned int lr = LR(j,candkey);
      if( lr==0 ) hit(xy,i)->SetHitPos( xy, wpos + dl );
      else hit(xy,i)->SetHitPos(xy, wpos - dl );
      hit(xy,i)->SetLeftRight((int)lr);
      j++;
    }
  }
  else{
    std::cerr << "usage: xy = 0 or 1 in LocalTrack::LeastSquareFit()" << std::endl;
    return false;
  }
  return true;
}

void LocalTrack::CalcHitPosition( const bool &TILT )
{
  for(int xy=0;xy<2;xy++)
    for( int i=0; i<nhit(xy); i++ ){
      double x,y,z;
      z = hit(xy,i)->wz();
      if(TILT)
	hit(xy,i)->gwpos(x,y,z,1,0);
      XYLocalPosatZ(z,x,y);
      hit(xy,i)->SetHitPosition(x,y,z);
    }
}

void LocalTrack::CalcResidual( const bool &ROT)
{
  for(int xy=0 ; xy<2; xy++)
    for( int i=0; i<nhit(xy); i++ ){
      ChamberLikeHit *tmphit = hit(xy,i);
      double wpos = tmphit->wpos(xy);
      double hpos = tmphit->hitpos(xy);
      if(ROT){
	double theta,wz;
	tmphit->gwpos(wpos,theta,wz,1,0);
	hpos=tmphit->x()*TMath::Cos(theta/180.*TMath::Pi())
	  +tmphit->y()*TMath::Sin(theta/180.*TMath::Pi());
      }
      double dltrack = hpos - wpos;
      double dl = tmphit->dl();
      int lr= tmphit->leftright();
      if(lr==0) tmphit->SetResolution( dltrack - dl );
      else tmphit->SetResolution( dltrack + dl );
    }
}

bool LocalTrack::LinearFit( ConfMan *conf, const int &xy )
{
#if 0
  std::cout << "!!! LocalTrack::LinearFit()" << std::endl;
#endif
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  int MaxNumOfHitsInTrack=BLDCParam->GetMaxHitInTrack()*2;

  if( nhit() > MaxNumOfHitsInTrack ){
    //    std::cout<<"!!! [LinearTrack::LinearFit()] too many hits !!!"<<std::endl;
    return false;
  }
  double x[MaxNumOfHitsInTrack], y[MaxNumOfHitsInTrack];
  int np = nhit(xy);
  int np_org = np;
  for( int i=0; i<np_org;i++){
    if(BLDCParam->layerkilled(hit(xy,i)->cid(),hit(xy,i)->layer())){
      np--;continue;
    }
    if( !layerstatus(hit(xy,i)->cid(),hit(xy,i)->layer()) ){
      np--; continue;
    }
  } 
  if( np < 2 ) return false;

  double minchi = 1.0e+9;
  double canda=0., candb=0., candc=0.;//,candd;
  unsigned int candkey=0x0;
  unsigned int key = 0x0;
  double slope=b();
  double intercept=c();
  if(xy){
    slope=e();
    intercept=f();
  }

  int nconb = (int)pow(2,np);
  int lr;
  for( int ic=0; ic<nconb; ic++ ){      
    int j=0;
    for( int i=0; i<np_org;i++){
      if(BLDCParam->layerkilled(hit(xy,i)->cid(),
				hit(xy,i)->layer())) continue;
      if( !layerstatus(hit(xy,i)->cid(),hit(xy,i)->layer()) ) continue;
      hit(xy,i)->gwpos2(y[j],x[j],slope,intercept,1,1);      
      //      std::cout<<"z,x,theta "<<x[j]<<"\t"<<y[j]<<"\t"<<temptheta<<std::endl; 
      lr=LR(j,key);
      if( lr==0 ) y[j] += hit(xy,i)->dl();
      else y[j] -= hit(xy,i)->dl();
      j++;
    }
    double Sx, Sy, Sxy, Sxx;
    Sx=Sy=Sxy=Sxx=0;
    for( int i=0; i<np; i++ ){
      Sx  += x[i];
      Sy  += y[i];
      Sxy += x[i]*y[i];
      Sxx += x[i]*x[i];
    }
    double D = (double)np*Sxx - Sx*Sx;
    if( D==0 ) return false;
    double aa = (Sxx*Sy - Sx*Sxy)/D;
    double bb = ((double)np*Sxy - Sx*Sy)/D;
    double chi = 0;
    for( int i=0; i<np; i++ ) chi += (y[i]-aa-bb*x[i])*(y[i]-aa-bb*x[i])/(SpatialResolutionOfBLDC*SpatialResolutionOfBLDC);
    
    if(np<3) chi=-1;
    else chi /= (double)(np-2);
#if SLOPECUT
    if( TMath::Abs(bb) < BLDCParam->GetMaxSlope()  ) 
#endif
      if( (chi < minchi) || (chi==minchi && candb*candb<bb*bb) ){
	canda = 1.;
	candb = -1.*bb;
	candc = -1.*aa;
	candkey=key;
	minchi=chi;	
      }
    key++;
  }
  if(xy)
    SetDEF(canda,candb,candc);
  else
    SetABC(canda,candb,candc);
  SetDof(xy,np-2);
  SetChisqr(xy,minchi);
  int j=0;
  for( int i=0; i<np_org; i++ ){
    double dl = hit(xy,i)->dl();
    double zpos; 
    double wpos;
    double tmpx,tmpy;
    hit(xy,i)->gwpos2(wpos,zpos,slope,intercept,1,1);      
    if(BLDCParam->layerkilled(hit(xy,i)->cid(),
			      hit(xy,i)->layer())||
       !layerstatus( hit(xy,i)->cid(), hit(xy,i)->layer()) ){
      XYLocalPosatZ(zpos,tmpx,tmpy);
      unsigned int lr=0;
      if(fabs(tmpx-(wpos+dl))>fabs(tmpx-(wpos-dl))) lr=1;
      if( lr==0 ) hit(xy,i)->SetHitPos( xy, wpos + dl );
      else hit(xy,i)->SetHitPos( xy, wpos - dl );
      hit(xy,i)->SetLeftRight((int)lr);	
      continue;
    }
    unsigned int lr = LR(j,candkey);
    if( lr==0 ) hit(xy,i)->SetHitPos( xy, wpos + dl );
    else hit(xy,i)->SetHitPos(xy, wpos - dl );
    hit(xy,i)->SetLeftRight((int)lr);
    j++;
  }
  return true;
}

// ----------------------------- //
// class BLDCCluster             //
// ------------------------------//

BLDCCluster::BLDCCluster( const BLDCCluster &right )
{
  for( BLDCHitCluster::const_iterator it=right.bldcHitCluster.begin(); it!=right.bldcHitCluster.end(); it++ ){
    bldcHitCluster.push_back( (*it) );
  }
}

BLDCCluster::~BLDCCluster()
{
  this->Clear();
}


// ----------------------------- //
// class BLDCClusterMan          //
// ------------------------------//

BLDCClusterMan::BLDCClusterMan( const BLDCClusterMan &right )
{
  for(int ud=0;ud<2;ud++)
    for(int xy=0;xy<2;xy++)
      for( BLDCClusterContainer::const_iterator it=right.bldcClusterContainer[ud+2*xy].begin(); it!=right.bldcClusterContainer[ud+2*xy].end(); it++ ){
	bldcClusterContainer[ud+2*xy].push_back( (*it) );
      }
}

BLDCClusterMan::~BLDCClusterMan()
{
  this->Clear();
}

void BLDCClusterMan::DeleteCluster( const int &ud, const int &xy, const int &i )
{
  BLDCClusterContainer::iterator it=bldcClusterContainer[ud+2*xy].begin();
  for( int j=0; j<i; j++ ) it++;
  bldcClusterContainer[ud+2*xy].erase( it );
}

void BLDCClusterMan::DeleteClusterUpXZ( const int &i )
{
  DeleteCluster(0,0,i);
}

void BLDCClusterMan::DeleteClusterDownXZ( const int &i )
{
  DeleteCluster(1,0,i);
}

void BLDCClusterMan::DeleteClusterUpYZ( const int &i )
{
  DeleteCluster(0,1,i);
}

void BLDCClusterMan::DeleteClusterDownYZ( const int &i )
{
  DeleteCluster(1,1,i);
}

void BLDCClusterMan::Clear()
{
  for(int ud=0;ud<2;ud++)
    for(int xy=0;xy<2;xy++)
      bldcClusterContainer[ud+2*xy].clear();
}

// ----------------------------- //
// class BeamLineTrackMan        //
// ------------------------------//

BeamLineTrackMan::BeamLineTrackMan()
{
  Clear();
}

BeamLineTrackMan::BeamLineTrackMan( const BeamLineTrackMan &right )
{
  for( BeamLineTrackContainer::const_iterator it=right.BLC1TrackContainer.begin(); it!=right.BLC1TrackContainer.end(); it++ ){
    BLC1TrackContainer.push_back( (*it) );
  }

  for( BeamLineTrackContainer::const_iterator it=right.BLC1aTrackContainer.begin(); it!=right.BLC1aTrackContainer.end(); it++ ){
    BLC1aTrackContainer.push_back( (*it) );
  }

  for( BeamLineTrackContainer::const_iterator it=right.BLC1bTrackContainer.begin(); it!=right.BLC1bTrackContainer.end(); it++ ){
    BLC1bTrackContainer.push_back( (*it) );
  }

  for( BeamLineTrackContainer::const_iterator it=right.BLC2TrackContainer.begin(); it!=right.BLC2TrackContainer.end(); it++ ){
    BLC2TrackContainer.push_back( (*it) );
  }

  for( BeamLineTrackContainer::const_iterator it=right.BLC2aTrackContainer.begin(); it!=right.BLC2aTrackContainer.end(); it++ ){
    BLC2aTrackContainer.push_back( (*it) );
  }

  for( BeamLineTrackContainer::const_iterator it=right.BLC2bTrackContainer.begin(); it!=right.BLC2bTrackContainer.end(); it++ ){
    BLC2bTrackContainer.push_back( (*it) );
  }

  for( BeamLineTrackContainer::const_iterator it=right.BPCTrackContainer.begin(); it!=right.BPCTrackContainer.end(); it++ ){
    BPCTrackContainer.push_back( (*it) );
  }
  for(int i=0;i<7;i++) STATUS[i]=right.STATUS[i];
}

BeamLineTrackMan::~BeamLineTrackMan()
{
  this->Clear();
}

void BeamLineTrackMan::Clear()
{
  for(int i=0;i<8;i++) STATUS[i]=-1;
  BLC1TrackContainer.clear();
  BLC1LinearTrackContainer.clear();
  BLC1aTrackContainer.clear();
  BLC1bTrackContainer.clear();
  BLC1aLinearTrackContainer.clear();
  BLC1bLinearTrackContainer.clear();

  BLC2TrackContainer.clear();
  BLC2LinearTrackContainer.clear();
  BLC2aTrackContainer.clear();
  BLC2bTrackContainer.clear();
  BLC2aLinearTrackContainer.clear();
  BLC2bLinearTrackContainer.clear();

  BPCTrackContainer.clear();
  BPCLinearTrackContainer.clear();
}

int BeamLineTrackMan::status(const int &cid){
  if(cid==CID_BLC1a) return STATUS[0];
  if(cid==CID_BLC1b) return STATUS[1];
  if(cid==CID_BLC2a) return STATUS[2];
  if(cid==CID_BLC2b) return STATUS[3];
  if(cid==CID_BPC) return STATUS[4];
  if(cid==CID_BLC1) return STATUS[5];
  if(cid==CID_BLC2) return STATUS[6];
  if(cid==CID_FDC1) return STATUS[7];
  return -1;
}

void BeamLineTrackMan::SetStatus(const int &cid, const int &sta){
  if(cid==CID_BLC1a)  STATUS[0]=sta;
  if(cid==CID_BLC1b)  STATUS[1]=sta;
  if(cid==CID_BLC2a)  STATUS[2]=sta;
  if(cid==CID_BLC2b)  STATUS[3]=sta;
  if(cid==CID_BPC)   STATUS[4]=sta;
  if(cid==CID_BLC1)  STATUS[5]=sta;
  if(cid==CID_BLC2)  STATUS[6]=sta;
  if(cid==CID_FDC1)   STATUS[7]=sta;
}

void BeamLineTrackMan::DeleteTrackBLC1( const int &i )
{    
  BeamLineTrackContainer::iterator it=BLC1TrackContainer.begin();
  for ( int j=0; j<i; j++ ) it++;
  BLC1TrackContainer.erase( it );
}

void BeamLineTrackMan::DeleteTrackBLC1a( const int &i )
{    
  BeamLineTrackContainer::iterator it=BLC1aTrackContainer.begin();
  for ( int j=0; j<i; j++ ) it++;
  BLC1aTrackContainer.erase( it );
}

void BeamLineTrackMan::DeleteTrackBLC1b( const int &i )
{    
  BeamLineTrackContainer::iterator it=BLC1bTrackContainer.begin();
  for ( int j=0; j<i; j++ ) it++;
  BLC1bTrackContainer.erase( it );
}

void BeamLineTrackMan::DeleteTrackBLC2( const int &i )
{    
  BeamLineTrackContainer::iterator it=BLC2TrackContainer.begin();
  for ( int j=0; j<i; j++ ) it++;
  BLC2TrackContainer.erase( it );
}

void BeamLineTrackMan::DeleteTrackBLC2a( const int &i )
{    
  BeamLineTrackContainer::iterator it=BLC2aTrackContainer.begin();
  for ( int j=0; j<i; j++ ) it++;
  BLC2aTrackContainer.erase( it );
}

void BeamLineTrackMan::DeleteTrackBLC2b( const int &i )
{    
  BeamLineTrackContainer::iterator it=BLC2bTrackContainer.begin();
  for ( int j=0; j<i; j++ ) it++;
  BLC2bTrackContainer.erase( it );
}

void BeamLineTrackMan::DeleteTrackBPC( const int &i )
{    
  BeamLineTrackContainer::iterator it=BPCTrackContainer.begin();
  for ( int j=0; j<i; j++ ) it++;
  BPCTrackContainer.erase( it );
}

// ------------------------------------------------------------------

bool BeamLineTrackMan::DoTracking( BeamLineHitMan *blMan, ConfMan *conf)
{
#if 0
  std::cout << "!!! BeamLineTrackMan::DoTracking()" << std::endl;
#endif 
  STATUS[4]=LocalTracking( blMan, conf, CID_BPC );
  STATUS[0]=LocalTracking( blMan, conf, CID_BLC1a );
  STATUS[1]=LocalTracking( blMan, conf, CID_BLC1b );
  STATUS[2]=LocalTracking( blMan, conf, CID_BLC2a );
  STATUS[3]=LocalTracking( blMan, conf, CID_BLC2b );
  //LinearTracking( conf, CID_BLC1 );
  //  LinearTracking( conf, CID_BLC2 );

  return 1;
}
//===================2011 10/18 sada============================================
//===================2012 6/29 hashimoto============================================
int BeamLineTrackMan::LocalTracking( BeamLineHitMan *blMan, ConfMan *conf, const int &id )
{
#if DEBUG
  std::cout << "!!! BeamLineTrackMan::LocalTracking()\tcid=" <<id<< std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  int MaxNumOfHitsInBLDC=BLDCParam->GetMaxBLDCHit();
  //  int MaxNumOfHitsInLayer=BLDCParam->GetMaxHitInLayer();
  for(int i=0;i<2;i++)
    for(int j=0;j<8;j++) LAYER[i][j]=true; 
  
  int ntrack=0;
  int nallhit = 0;
  int NumOfBLDCLayers=0;
  BeamLineTrackContainer TmpTrackContainer[2];
  if(id==CID_BLC1a || id==CID_BLC1b || id==CID_BLC2a || id==CID_BLC2b) NumOfBLDCLayers=NumOfBLCLayers;
  else if(id==CID_BPC ) NumOfBLDCLayers=NumOfBPCLayers;

  for( int layer=1; layer<=NumOfBLDCLayers; layer++ ){
    if( BLDCParam->layerdead(id,layer)) { continue; }
#if KILLLAYER
    int tmp_nhit=0;
#if CHECKRANGE
    for( int i=0; i<blMan->nBLDC(id,layer); i++)
      if(blMan->BLDC(id,layer,i)->CheckRange()) tmp_nhit++;
#else 
    tmp_nhit = blMan->nBLDC(id,layer);
#endif
    if( BLDCParam->layerkilled(id,layer)||
	tmp_nhit>BLDCParam->GetMaxHitInLayer()){
      LAYER[0][layer-1] = false;
      continue;       
    }
    nallhit += tmp_nhit;
#else
    nallhit += blMan->nBLDC(id,layer);
#endif
  }
  if( nallhit > MaxNumOfHitsInBLDC ){
    //    std::cout<<"too many hits in BLDC !!!"<<std::endl;
    //return ntrack;
    return 2;
  }

  BLDCClusterMan *clMan = new BLDCClusterMan();
  if( !Clustering( blMan, clMan, conf, id ) ){
    delete clMan; 
    std::cout<<"clustring failed !!!"<<std::endl;
    //    return ntrack; 
    return 3;
  }
#if 0
  for(int xy=0;xy<2;xy++){
    std::cout<<"id,xy:  "<<id<<"\t"<<xy<<std::endl;
    std::cout<<"ncu:  "<<clMan->ncluster(0,xy)<<std::endl;
    for( int icu=0; icu<clMan->ncluster(0,xy); icu++ ){
      BLDCCluster *cl=clMan->cluster(0,xy,icu);
      for( int i=0; i<cl->nhit(); i++ ) {
	ChamberLikeHit *hit = cl->hit(i);
	std::cout<<"ihit,layer,wire:   "<<i<<"\t"<<hit->layer()<<"\t"<<hit->wire()<<std::endl;
      }
    }
    std::cout<<"ncd:  "<<clMan->ncluster(1,xy)<<std::endl;
    for( int icd=0; icd<clMan->ncluster(1,xy); icd++ ){
      BLDCCluster *cl=clMan->cluster(1,xy,icd);
      for( int i=0; i<cl->nhit(); i++ ) {
	ChamberLikeHit *hit = cl->hit(i);
	std::cout<<"ihit,layer,wire:   "<<i<<"\t"<<hit->layer()<<"\t"<<hit->wire()<<std::endl;
      }
    }
  }
#endif
  for(int xy=0;xy<2;xy++){
    int tmp_ntrack=0;
    bool ALLCL=false;
    bool MINHIT=false;
    bool SLOPE=false;
    bool FIT=false;
    bool CHI2=false;
#if DEL
    bool WHILE=true;
    while(clMan->ncluster(0,xy)>0&&clMan->ncluster(1,xy)>0&&WHILE){
      double minchi2=9999;
      LocalTrack *tmptra=0;
      int iclu,icld;
#endif
      for( int icu=0; icu<clMan->ncluster(0,xy); icu++ ){
	for( int icd=0; icd<clMan->ncluster(1,xy); icd++ ){
#if DEBUG
	  std::cout<<"start icu,icd:  "<<icu<<"\t"<<icd<<std::endl; 
#endif	
	  ALLCL=true;
	  BLDCCluster *cl[2];
	  cl[0] = clMan->cluster( 0, xy, icu );
	  cl[1] = clMan->cluster( 1, xy, icd );
	  LocalTrack *track = new LocalTrack();
	  track->SetCID(id,id);
	  track->SetLayerStatus(id,LAYER[0]);
	  int tmp_nhit=0;
	  for(int ud=0;ud<2;ud++){
	    for( int i=0; i<cl[ud]->nhit(); i++ ) {
	      ChamberLikeHit *hit = cl[ud]->hit(i);
	      track->SetHit( xy, (*hit) );
	      if( BLDCParam->layerkilled(hit->cid(),hit->layer())==0 && LAYER[0][hit->layer()-1])
		tmp_nhit++;
	    }
	  }
	  if( tmp_nhit<BLDCParam->GetMinHit(xy,id) ) {
	    delete track;	  continue;
	  }	
	  MINHIT=true;
	  if( !LeastSquareFit( track, conf, xy ) ) {
	    delete track;	  continue;
	  }	
	  FIT=true;
	  
#if SLOPECUT
	  double slope=track->b();
	  if(xy) slope=track->e();
	  if( TMath::Abs(slope) > 0.1  ) {
	    delete track;   continue;
	  }	  
	  SLOPE=true;
#endif
	  if( track->chi2(xy)>BLDCParam->GetMaxChi() ) {
	    delete track;   continue;
	  }	
	  CHI2=true;  
#if DEL
	  if( track->chi2(xy)>minchi2 ) {
	    delete track;   continue;
	  }else{
	    if(tmptra) delete tmptra;
	    tmptra=track;
	    minchi2=track->chi2(xy);
	    iclu=icu;
	    icld=icd;
	  }	
#else
	  TmpTrackContainer[xy].push_back(*track);
	  delete track;
	  tmp_ntrack++;
#endif
	}
      }    
#if DEL
      if(tmptra){
	TmpTrackContainer[xy].push_back(*tmptra);
	delete tmptra;
	clMan->DeleteCluster( 0, xy, iclu );
	clMan->DeleteCluster( 1, xy, icld );      
	tmp_ntrack++;
      }else{
	WHILE=false;
      }
    }// while
#endif

    if( tmp_ntrack==0 ) {
      TmpTrackContainer[0].clear();
      TmpTrackContainer[1].clear();
      delete clMan; 
      if(!ALLCL) return 4;
      if(!MINHIT) return 5;
      if(!FIT) return 6;
#if SLOPECUT
      if(!SLOPE) return 8;
#endif
      if(!CHI2) return 7;
      return 10;
    }
#if DEBUG
    std::cout<<"xy, ntrack   "<<xy<<tmp_ntrack<<std::endl;
#endif
  } //xy


  delete clMan; 
  // Set a track candidate
  for(int tr1=0;tr1<(int)TmpTrackContainer[0].size();tr1++)
    {
      for(int tr2=0;tr2<(int)TmpTrackContainer[1].size();tr2++)
	{	
	  LocalTrack *trcand = new LocalTrack();
	  LocalTrack *tmptr[2];
	  tmptr[0] = &(TmpTrackContainer[0][tr1]);
	  tmptr[1] = &(TmpTrackContainer[1][tr2]);
	  for(int xy=0;xy<2;xy++){
	    for( int i=0; i<tmptr[xy]->nhit(xy); i++ ){
	      ChamberLikeHit *hit = tmptr[xy]->hit(xy,i);
	      trcand->SetHit( xy, (*hit) );
	    }
	    int dof; double chi;
	    chi = tmptr[xy]->chi2(xy);
	    dof = tmptr[xy]->dof(xy);
	    trcand->SetChisqr( xy, chi );
	    trcand->SetDof( xy, dof );
	  }
	  double a,b,c,d,e,f;
	  tmptr[0]->abc(a,b,c);
	  tmptr[1]->def(d,e,f);
	  trcand->SetABC(a,b,c);
	  trcand->SetDEF(d,e,f);
	  trcand->CalcHitPosition();
	  trcand->CalcResidual();	  
	  if(id==CID_BLC2a)   SetTrackBLC2a( (*trcand) );
	  else if(id==CID_BLC2b) SetTrackBLC2b( (*trcand) );
	  else if(id==CID_BLC1a) SetTrackBLC1a( (*trcand) );
	  else if(id==CID_BLC1b) SetTrackBLC1b( (*trcand) );
	  else if(id==CID_BPC) SetTrackBPC( (*trcand) );
	  ntrack++;
	  delete trcand;
	}
    }
#if DEBUG
  std::cout << "!!! BeamLineTrackMan::Local Tracking finished" << std::endl;
#endif
  ConvLocalToGlobal(id);
  TmpTrackContainer[0].clear();
  TmpTrackContainer[1].clear();

  return 1;
  //  return ntrack;
}

int BeamLineTrackMan::SemiLocalTracking( BeamLineHitMan *blMan, ConfMan *conf, const int &id )
{
#if DEBUG2
  std::cout << "!!! BeamLineTrackMan::SemiLocalTracking()" << std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  int MaxNumOfHitsInBLDC=BLDCParam->GetMaxBLDCHit();
  //  int MaxNumOfHitsInLayer=BLDCParam->GetMaxHitInLayer();
  
  int ntrack=0;
  int NumOfBLDCLayers=0;
  int id2[2];
  BeamLineTrackContainer TmpTrackContainer[2];
  if(id==CID_BLC1 || id==CID_BLC2a || id ==CID_BLC2b ){
    id2[0]=CID_BLC2a; 
    id2[1]=CID_BLC2b; 
    NumOfBLDCLayers=NumOfBLCLayers;
  }
  else if(id==CID_BLC1 || id==CID_BLC1a || id == CID_BLC1b ){
    id2[0]=CID_BLC1a; 
    id2[1]=CID_BLC1b; 
    NumOfBLDCLayers=NumOfBLCLayers;
  }else return 0;

  BLDCClusterMan *clMan[2];
  for(int ii=0;ii<2;ii++){
    int nallhit = 0;
    for( int layer=1; layer<=NumOfBLDCLayers; layer++ ){
      if( BLDCParam->layerdead(id2[ii],layer)) continue;   
      nallhit += blMan->nBLDC(id2[ii],layer);
    }
    if( nallhit > MaxNumOfHitsInBLDC ){
      //    std::cout<<"too many hits in BLDC !!!"<<std::endl;
      return ntrack;
    }
    clMan[ii] = new BLDCClusterMan();
    if( !Clustering( blMan, clMan[ii], conf, id2[ii] ) ){
      delete clMan[ii]; 
      std::cout<<"clustring failed !!!"<<std::endl;
      return ntrack; 
    }
  }

  for(int xy=0;xy<2;xy++){
    int tmp_ntrack=0;
    for( int icu1=0; icu1<clMan[0]->ncluster(0,xy); icu1++ ){
      for( int icd1=0; icd1<clMan[0]->ncluster(1,xy); icd1++ ){
	for( int icu2=0; icu2<clMan[1]->ncluster(0,xy); icu2++ ){
	  for( int icd2=0; icd2<clMan[1]->ncluster(1,xy); icd2++ ){
	    BLDCCluster *cl[4];
	    cl[0] = clMan[0]->cluster( 0, xy, icu1 );
	    cl[1] = clMan[0]->cluster( 1, xy, icd1 );
	    cl[2] = clMan[1]->cluster( 0, xy, icu2 );
	    cl[3] = clMan[1]->cluster( 1, xy, icd2 );
	    LocalTrack *track = new LocalTrack();
	    int tmp_nhit=0;
	    for(int ud=0;ud<4;ud++){
	      for( int i=0; i<cl[ud]->nhit(); i++ ) {
		ChamberLikeHit *hit = cl[ud]->hit(i);
		track->SetHit( xy, (*hit) );
		if( BLDCParam->layerkilled(hit->cid(),hit->layer())==0 )
		  tmp_nhit++;
	      }
	    }
#if DEBUG2
	    std::cout<<"xy, icu1, icd1, icu2, icd2, nhit   "<<xy<<"\t"<<icu1<<icd1<<icu2<<icd2<<"\t"
		     <<tmp_nhit<<std::endl;
#endif
	    if( tmp_nhit<BLDCParam->GetMinHit(xy,id2[0]) ) {
	      delete track;	  continue;
	    }	
	    if( !LeastSquareFit( track, conf, xy ) ) {
	      delete track;	  continue;
	    }	
	    if( track->chi2(xy)>BLDCParam->GetMaxChi() ) {
	      delete track;   continue;
	    }	  
	    TmpTrackContainer[xy].push_back(*track);
	    delete track;
	    tmp_ntrack++;
	  }
	}    
      }
    }
    if( tmp_ntrack==0 ) {
      TmpTrackContainer[0].clear();
      TmpTrackContainer[1].clear();
      delete clMan[0]; 
      delete clMan[1]; 
      return ntrack;
    }
#if DEBUG2
    std::cout<<"xy, ntrack   "<<xy<<tmp_ntrack<<std::endl;
#endif
  } //xy
  delete clMan[0]; 
  delete clMan[1]; 
  // Set a track candidate
  for(int tr1=0;tr1<(int)TmpTrackContainer[0].size();tr1++)
    {
      for(int tr2=0;tr2<(int)TmpTrackContainer[1].size();tr2++)
	{	
	  LocalTrack *trcand = new LocalTrack();
	  LocalTrack *tmptr[2];
	  tmptr[0] = &(TmpTrackContainer[0][tr1]);
	  tmptr[1] = &(TmpTrackContainer[1][tr2]);
	  for(int xy=0;xy<2;xy++){
	    for( int i=0; i<tmptr[xy]->nhit(xy); i++ ){
	      ChamberLikeHit *hit = tmptr[xy]->hit(xy,i);
	      trcand->SetHit( xy, (*hit) );
	    }
	    int dof; double chi;
	    chi = tmptr[xy]->chi2(xy);
	    dof = tmptr[xy]->dof(xy);
	    trcand->SetChisqr( xy, chi );
	    trcand->SetDof( xy, dof );
	  }
	  double a,b,c,d,e,f;
	  tmptr[0]->abc(a,b,c);
	  tmptr[1]->def(d,e,f);
	  trcand->SetABC(a,b,c);
	  trcand->SetDEF(d,e,f);
	  trcand->CalcHitPosition();
	  trcand->CalcResidual();	  
	  if(id==CID_BLC2a||id==CID_BLC2b||id ==CID_BLC2)      SetTrackBLC2( (*trcand) );
	  else if(id==CID_BLC1a||id==CID_BLC1b|| id==CID_BLC1)     SetTrackBLC1( (*trcand) );
	  ntrack++;
	  delete trcand;
	}
    }
#if DEBUG
  std::cout << "!!! BeamLineTrackMan::Semi Local Tracking finished" << std::endl;
#endif
  ConvLocalToGlobal(id);
  return ntrack;
}

int BeamLineTrackMan::LinearTracking( BeamLineHitMan *blMan, ConfMan *conf, const int &id )
{
#if DEBUG
  std::cout << "!!! BeamLineTrackMan::LinearTracking()\tcid=" <<id<< std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  int MaxNumOfHitsInBLDC=BLDCParam->GetMaxBLDCHit();
  //  int MaxNumOfHitsInLayer=BLDCParam->GetMaxHitInLayer();
 
  for(int i=0;i<2;i++)
    for(int j=0;j<8;j++) LAYER[i][j]=true; 
  
  int nallhit = 0;
  int NumOfBLDCLayers=0;
  BeamLineTrackContainer TmpTrackContainer[2];
  int id2[3];
  if(id==CID_BLC2a || id==CID_BLC2b || id ==CID_BLC2 ){
    id2[0]=CID_BLC2a; 
    id2[1]=CID_BLC2b; 
    id2[2]=CID_BLC2;
    NumOfBLDCLayers=NumOfBLCLayers;
  }
  else if(id==CID_BLC1 || id==CID_BLC1a || id == CID_BLC1b ){
    id2[0]=CID_BLC1a; 
    id2[1]=CID_BLC1b;
    id2[2]=CID_BLC1; 
    NumOfBLDCLayers=NumOfBLCLayers;
  }else return 0;
  
  //  int tmp_ntrack=0;
  //----------------------------
  // initialization done
  //----------------------------
  //check total number and busy layers
  //----------------------------
  for(int ic=0;ic<2;ic++){
    for( int layer=1; layer<=NumOfBLDCLayers; layer++ ){
      if( BLDCParam->layerdead(id2[ic],layer)) continue; 
#if KILLLAYER
      int tmp_nhit=0;
#if CHECKRANGE
      for( int i=0; i<blMan->nBLDC(id2[ic],layer); i++)
	if(blMan->BLDC(id2[ic],layer,i)->CheckRange()) tmp_nhit++;
#else 
      tmp_nhit = blMan->nBLDC(id2[ic],layer);
#endif
      if( BLDCParam->layerkilled(id2[ic],layer)||
	  tmp_nhit>BLDCParam->GetMaxHitInLayer()){
	LAYER[ic][layer-1] = false;
	continue;       
      }
      nallhit += tmp_nhit;
#else
      nallhit += blMan->nBLDC(id2[ic],layer);
#endif
    }
  }
  if( nallhit > MaxNumOfHitsInBLDC * 2 ){
    //    std::cout<<"too many hits in BLDC !!!"<<std::endl;
    //return ntrack;
    return 2;
  }
  //----------------------------
  //clustering
  //----------------------------
  BLDCClusterMan *clMan[2];
  clMan[0] = new BLDCClusterMan();
  clMan[1] = new BLDCClusterMan();
  for(int ic=0;ic<2;ic++){
    if( !Clustering( blMan, clMan[ic], conf, id2[ic], ic ) ){
      delete clMan[0]; 
      delete clMan[1]; 
      std::cout<<"clustring failed !!!\tid"<<id2[ic]<<std::endl;
      //    return ntrack; 
      return 3;
    }
  }
#if 0
  for(int ic=0;ic<2;ic++){
    for(int xy=0;xy<2;xy++){
      std::cout<<"id,xy:  "<<id2[ic]<<"\t"<<xy<<std::endl;
      std::cout<<"ncu:  "<<clMan[ic]->ncluster(0,xy)<<std::endl;
      for( int icu=0; icu<clMan[ic]->ncluster(0,xy); icu++ ){
	BLDCCluster *cl=clMan[ic]->cluster(0,xy,icu);
	for( int i=0; i<cl->nhit(); i++ ) {
	  ChamberLikeHit *hit = cl->hit(i);
	  std::cout<<"ihit,layer,wire:   "<<i<<"\t"<<hit->layer()<<"\t"<<hit->wire()<<std::endl;
	}
      }
      std::cout<<"ncd:  "<<clMan[ic]->ncluster(1,xy)<<std::endl;
      for( int icd=0; icd<clMan[ic]->ncluster(1,xy); icd++ ){
	BLDCCluster *cl=clMan[ic]->cluster(1,xy,icd);
	for( int i=0; i<cl->nhit(); i++ ) {
	  ChamberLikeHit *hit = cl->hit(i);
	  std::cout<<"ihit,layer,wire:   "<<i<<"\t"<<hit->layer()<<"\t"<<hit->wire()<<std::endl;
	}
      }
    }
  }
#endif
  //----------------------------
  // require more than 3 slayers which have at least one cluster
  //----------------------------
  int ncomb[2]={1,1};
  bool CLUSTER[2]={true,true};
  int ncluster[2][4]={};
  
  for(int xy=0;xy<2;xy++){
    int nslay=0;
    for(int ic=0;ic<2;ic++) {
      for(int iud=0;iud<2;iud++){
	if(clMan[ic]->ncluster(iud,xy)>0){
	  ncomb[xy] *= clMan[ic]->ncluster(iud,xy);
	  ncluster[xy][2*ic+iud]=clMan[ic]->ncluster(iud,xy);
	  nslay++;
	}
      }
    }
    if(nslay<3) CLUSTER[xy]=false;
  }

  //----------------------------
  // search track candidate in X/Y layers
  //----------------------------
  int icl[4]={0};
  if(CLUSTER[0]&&CLUSTER[1]){
    for(int xy=0;xy<2;xy++){
      int tmp_ntrack=0;
      bool ALLCL=false;
      bool MINHIT=false;
      bool SLOPE=false;
      bool FIT=false;
      bool CHI2=false;

      for(int icomb=0;icomb<ncomb[xy];icomb++){
	BLDCCluster *cl[4];
	int tmp=icomb;
	int ncl=0;
	for(int ic=0;ic<2;ic++) {
	  for(int iud=0;iud<2;iud++){
	    int islay=2*ic+iud;
	    if(ncluster[xy][islay]>0){
	      icl[islay]=tmp%ncluster[xy][islay];
	      tmp=tmp/ncluster[xy][islay];
	      cl[ncl] = clMan[ic]->cluster( iud, xy, icl[islay] );
	      ncl++;
	    }
	  }
	}
	LocalTrack *track = new LocalTrack();
	track->SetCID(id2[0],id2[1]);
	track->SetLayerStatus(id2[0],LAYER[0]);
	track->SetLayerStatus(id2[1],LAYER[1]);
	int tmp_nhit=0;
	for(int icl=0;icl<ncl;icl++){		  
	  for( int i=0; i<cl[icl]->nhit(); i++ ) {
	    ChamberLikeHit *hit = cl[icl]->hit(i);
	    track->SetHit( xy, (*hit) );
	    if( BLDCParam->layerkilled(hit->cid(),hit->layer())==0 &&
		track->layerstatus(hit->cid(),hit->layer()) ){	      
	      tmp_nhit++;
	    }
	  }
	}
	//	if( tmp_nhit<BLDCParam->GetMinHit(xy,id) ) {
	if( tmp_nhit<5 ) {
	  delete track;	  continue;
	}	
	MINHIT=true;
	if( !track->PreTracking( conf, xy ) ) {
	  delete track;	  continue;
	}	
	FIT=true;	  
	//	if( track->chi2(xy)>BLDCParam->GetMaxChi() ) {
	//	  delete track;   continue;
	//	}	
	//	CHI2=true;  

	TmpTrackContainer[xy].push_back(*track);
	delete track;
	tmp_ntrack++;	 
      }//icomb
      if( tmp_ntrack==0 ) {
	TmpTrackContainer[0].clear();
	TmpTrackContainer[1].clear();
	delete clMan[0]; 
	delete clMan[1]; 
	if(!ALLCL) return 4;
	if(!MINHIT) return 5;
	if(!FIT) return 6;
#if SLOPECUT
	if(!SLOPE) return 8;
#endif
	if(!CHI2) return 7;
	return 10;
      }
#if DEBUG
      std::cout<<"xy, ntrack   "<<xy<<tmp_ntrack<<std::endl;
#endif
    }//xy
  }// 3 slayers
  

  delete clMan[0]; 
  delete clMan[1]; 
  // Set a track candidate


  int ntrack=0;
  for(int tr1=0;tr1<(int)TmpTrackContainer[0].size();tr1++)
    {
      for(int tr2=0;tr2<(int)TmpTrackContainer[1].size();tr2++)
	{	
	  LocalTrack *trcand = new LocalTrack();
	  LocalTrack *tmptr[2];
	  tmptr[0] = &(TmpTrackContainer[0][tr1]);
	  tmptr[1] = &(TmpTrackContainer[1][tr2]);
	  for(int xy=0;xy<2;xy++){
	    for( int i=0; i<tmptr[xy]->nhit(xy); i++ ){
	      ChamberLikeHit *hit = tmptr[xy]->hit(xy,i);
	      trcand->SetHit( xy, (*hit) );
	    }
	  }
	  double a,b,c,d,e,f;
	  tmptr[0]->abc(a,b,c);
	  tmptr[1]->def(d,e,f);
	  trcand->SetABC(a,b,c);
	  trcand->SetDEF(d,e,f);

	  //	  for(int xy=0;xy<2;xy++){
	  if( !trcand->LinearFit( conf, 0 ) ) {
	    std::cout<<"LocalTrack::LinearFit() failed for x"<<std::endl;
	    delete trcand;	  continue;
	  }	
	  if( !trcand->LinearFit( conf, 1 ) ) {
	    std::cout<<"LocalTrack::LinearFit() failed for y"<<std::endl;
	    delete trcand;	  continue;
	  }	
	  if( trcand->chi2all() > BLDCParam->GetMaxChi() ) {
	    delete trcand;   continue;
	  }	
	  //CHI2=true;
	  trcand->CalcHitPosition();
	  trcand->CalcResidual();	  	  
	  if(id==CID_BLC2)   SetTrackBLC2( (*trcand) );
	  else if(id==CID_BLC1) SetTrackBLC1( (*trcand) );
	  ntrack++;
	  delete trcand;	  
	} //tr2
    }//tr1
  TmpTrackContainer[0].clear();
  TmpTrackContainer[1].clear();
  ConvLocalToGlobal(id);
#if DEBUG
  std::cout << "!!! BeamLineTrackMan::Linear Tracking finished" << std::endl;
#endif
  return 1;  
}

bool BeamLineTrackMan::ConvertLocalToLinear( const int &cid ){
  if(cid==CID_BLC1)    BLC1LinearTrackContainer.clear();
  else if(cid==CID_BLC2)    BLC2LinearTrackContainer.clear();
  else return false;
  for(int i=0; i<ntrackBLDC(cid);i++){
    LinearTrack* ltrack=new LinearTrack();
    LocalTrack *track=trackBLDC(cid,i);
    ltrack->SetCID(track->GetCID(0),track->GetCID(1));
    bool tmp[8];
    track->GetLayerStatus(0,tmp);
    ltrack->SetLayerStatus(track->GetCID(0),tmp);
    track->GetLayerStatus(1,tmp);
    ltrack->SetLayerStatus(track->GetCID(1),tmp);
    for(int xy=0; xy<2; xy++){
      for( int i=0; i<track->nhit(xy); i++ ) {
	ChamberLikeHit *hit = track->hit(xy,i);
	ltrack->SetHit( (*hit) );
      }
    }
    if(ltrack->LinearFit(0,false)){
      ltrack->CalcHitPosition();
      ltrack->CalcResidual();
      ltrack->ConvLocalToGlobal();
      SetLinearTrack(cid, (*ltrack) );
      delete ltrack;
    }	      
  }
  if(cid==CID_BLC1)    BLC1TrackContainer.clear();
  else if(cid==CID_BLC2)    BLC2TrackContainer.clear();
  return true;
}
int BeamLineTrackMan::LinearTracking( ConfMan *conf, const int &id )
{
#if DEBUG2
  std::cout << "!!! BeamLineTrackMan::LinearTracking()" << std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  //  int MaxNumOfHitsInBLDC=BLDCParam->GetMaxBLDCHit();
  //  int MaxNumOfHitsInLayer=BLDCParam->GetMaxHitInLayer();
  
  int ntrack=0;
  int NumOfBLDCLayers=0;
  int id2[3];
  if(id==CID_BLC2a || id==CID_BLC2b || id ==CID_BLC2 ){
    id2[0]=CID_BLC2a; 
    id2[1]=CID_BLC2b; 
    id2[2]=CID_BLC2;
    NumOfBLDCLayers=NumOfBLCLayers;
  }
  else if(id==CID_BLC1a || id==CID_BLC1b || id == CID_BLC1 ){
    id2[0]=CID_BLC1a; 
    id2[1]=CID_BLC1b;
    id2[2]=CID_BLC1; 
    NumOfBLDCLayers=NumOfBLCLayers;
  }else return 0;
  
  //  int tmp_ntrack=0;

  //  std::cout<<"itr1,itr2: "<<ntrackBLDC(id2[0])<<"\t"<<ntrackBLDC(id2[1])<<std::endl;
  if(ntrackBLDC(id2[0])!=1||ntrackBLDC(id2[1])!=1){
    //    return ntrack;
    SetStatus(id2[2],2);
    return 2;
  }
  bool MINHIT=false;
  bool FIT=false;
  bool CHI2=false;

  for( int itr1=0; itr1<ntrackBLDC(id2[0]) ; itr1++ ){
    LocalTrack *tr[2];
    tr[0] = trackBLDC( id2[0], itr1 );
    for( int itr2=0; itr2<ntrackBLDC(id2[1]) ; itr2++ ){
      tr[1] = trackBLDC( id2[1], itr2 );
      LinearTrack *track = new LinearTrack();
      int tmp_nhit=0;
      for(int ic=0;ic<2;ic++){
	for(int xy=0; xy<2; xy++){
	  for( int i=0; i<tr[ic]->nhit(xy); i++ ) {
	    ChamberLikeHit *hit = tr[ic]->hit(xy,i);
	    track->SetHit( (*hit) );
	    if( BLDCParam->layerkilled(hit->cid(),hit->layer())==0 )
	      tmp_nhit++;
	  }
	}
      }
      if( tmp_nhit<BLDCParam->GetMinHit(0,id2[0]) * 2 ) {
	delete track;	  continue;
      }	
      MINHIT=true;
      if( !track->LinearFit( conf ) ) {
	delete track;	  continue;
      }
      FIT=true;
      //      std::cout<<"chi2all "<<track->chi2all()<<std::endl;	
      if( track->chi2all() >BLDCParam->GetMaxChi() ) {
	delete track;   continue;
      }	  
      CHI2=true;
      track->CalcHitPosition();
      track->CalcResidual();
      track->ConvLocalToGlobal();	  
      if(id==CID_BLC2a||id==CID_BLC2b||id ==CID_BLC2) 
	SetLinearTrackBLC2( (*track) );
      else if(id==CID_BLC1a||id==CID_BLC1b|| id==CID_BLC1) 
	SetLinearTrackBLC1( (*track) );
      ntrack++;
      delete track;
    }   //itr2       
  } //itr1
#if DEBUG2
  std::cout<<"ntrack   "<<ntrack<<std::endl;
#endif
#if DEBUG
  std::cout << "!!! BeamLineTrackMan::Linear Tracking finished" << std::endl;
#endif
  if(!MINHIT){
    SetStatus(id2[2],3);
    return 3;
  }
  if(!FIT){
    SetStatus(id2[2],4);
    return 4;
  }
  if(!CHI2){
    SetStatus(id2[2],5);
    return 5;
  }
  if(ntrack==1){
    SetStatus(id2[2],1);  
    return ntrack;
  }else{
    SetStatus(id2[2],10);
    return 10;
  }
 
}

int BeamLineTrackMan::LocalLinearTracking( ConfMan *conf, const int &id )
{
#if DEBUG2
  std::cout << "!!! BeamLineTrackMan::LocalLinearTracking()" << std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  //  int MaxNumOfHitsInBLDC=BLDCParam->GetMaxBLDCHit();
  //  int MaxNumOfHitsInLayer=BLDCParam->GetMaxHitInLayer();
  
  int ntrack=0;
  int NumOfBLDCLayers=0;
  if(id==CID_BLC2a || id==CID_BLC2b || id ==CID_BLC2 ){
    NumOfBLDCLayers=NumOfBLCLayers;
  }
  else if(id==CID_BLC1a || id==CID_BLC1b || id == CID_BLC1 ){
    NumOfBLDCLayers=NumOfBLCLayers;
  }else return 0;
  
  //  int tmp_ntrack=0;
  for( int itr1=0; itr1<ntrackBLDC(id) ; itr1++ ){
    LocalTrack *tr;
    tr = trackBLDC( id, itr1 );
    if(tr->chi2all()>20) continue;
    LinearTrack *track = new LinearTrack();
    int tmp_nhit=0;
    for(int xy=0; xy<2; xy++){
      for( int i=0; i<tr->nhit(xy); i++ ) {
	ChamberLikeHit *hit = tr->hit(xy,i);
	track->SetHit( (*hit) );
	if( BLDCParam->layerkilled(hit->cid(),hit->layer())==0 )
	  tmp_nhit++;
      }
    }//xy      
    if( tmp_nhit<BLDCParam->GetMinHit(0,id) ) {
      delete track;	  continue;
    }	
    if( !track->LinearFit( conf ) ) {
      delete track;	  continue;
    }
    //      std::cout<<"chi2all "<<track->chi2all()<<std::endl;	
    if( track->chi2all() >BLDCParam->GetMaxChi() ) {
      delete track;   continue;
    }	  
    track->CalcHitPosition();
    track->CalcResidual();
    track->ConvLocalToGlobal();	  
    SetLinearTrack(id, (*track) );
    ntrack++;
    delete track;
  } //itr1
#if DEBUG2
  std::cout<<"ntrack   "<<ntrack<<std::endl;
#endif
#if DEBUG
  std::cout << "!!! BeamLineTrackMan::LocalLinearTracking finished" << std::endl;
#endif
  return ntrack;
}

int BeamLineTrackMan::ntrackBLDC(const int &cid)
{
  if(cid==CID_BPC)
    return ntrackBPC();
  if(cid==CID_BLC1)
    return ntrackBLC1();
  if(cid==CID_BLC1a)
    return ntrackBLC1a();
  if(cid==CID_BLC1b)
    return ntrackBLC1b();
  if(cid==CID_BLC2)
    return ntrackBLC2();
  if(cid==CID_BLC2a)
    return ntrackBLC2a();
  if(cid==CID_BLC2b)
    return ntrackBLC2b();
  return 0;
}

int BeamLineTrackMan::nltrackBLDC(const int &cid)
{
  if(cid==CID_BPC)
    return nltrackBPC();
  if(cid==CID_BLC1)
    return nltrackBLC1();
  if(cid==CID_BLC1a)
    return nltrackBLC1a();
  if(cid==CID_BLC1b)
    return nltrackBLC1b();
  if(cid==CID_BLC2)
    return nltrackBLC2();
  if(cid==CID_BLC2a)
    return nltrackBLC2a();
  if(cid==CID_BLC2b)
    return nltrackBLC2b();
  return 0;
}

void BeamLineTrackMan::SetLinearTrack(const int &cid, LinearTrack track)
{
  if(cid==CID_BPC)
    BPCLinearTrackContainer.push_back(track);
  if(cid==CID_BLC1)
    BLC1LinearTrackContainer.push_back(track);
  if(cid==CID_BLC1a)
    BLC1aLinearTrackContainer.push_back(track);
  if(cid==CID_BLC1b)
    BLC1bLinearTrackContainer.push_back(track);
  if(cid==CID_BLC2)
    BLC2LinearTrackContainer.push_back(track);
  if(cid==CID_BLC2a)
    BLC2aLinearTrackContainer.push_back(track);
  if(cid==CID_BLC2b)
    BLC2bLinearTrackContainer.push_back(track);
}

LocalTrack* BeamLineTrackMan::trackBLDC(const int &cid, const unsigned int &i)
{
  switch( cid ){
  case CID_BPC:
    return (0<=i&&i<BPCTrackContainer.size()) ? &BPCTrackContainer[i] : 0;
  case CID_BLC1:
    return (0<=i&&i<BLC1TrackContainer.size()) ? &BLC1TrackContainer[i] : 0;
  case CID_BLC1a:
    return (0<=i&&i<BLC1aTrackContainer.size()) ? &BLC1aTrackContainer[i] : 0;
  case CID_BLC1b:
    return (0<=i&&i<BLC1bTrackContainer.size())  ? &BLC1bTrackContainer[i] : 0;
  case CID_BLC2:
    return (0<=i&&i<BLC2TrackContainer.size()) ? &BLC2TrackContainer[i] : 0;
  case CID_BLC2a:
    return (0<=i&&i<BLC2aTrackContainer.size()) ? &BLC2aTrackContainer[i] : 0;
  case CID_BLC2b:
    return (0<=i&&i<BLC2bTrackContainer.size())  ? &BLC2bTrackContainer[i] : 0;
  }
  return 0;
}
LinearTrack* BeamLineTrackMan::ltrackBLDC(const int &cid, const unsigned int &i)
{
  switch( cid ){
  case CID_BPC:
    return (0<=i&&i<BPCLinearTrackContainer.size()) ? &BPCLinearTrackContainer[i] : 0;
  case CID_BLC1:
    return (0<=i&&i<BLC1LinearTrackContainer.size()) ? &BLC1LinearTrackContainer[i] : 0;
  case CID_BLC1a:
    return (0<=i&&i<BLC1aLinearTrackContainer.size()) ? &BLC1aLinearTrackContainer[i] : 0;
  case CID_BLC1b:
    return (0<=i&&i<BLC1bLinearTrackContainer.size())  ? &BLC1bLinearTrackContainer[i] : 0;
  case CID_BLC2:
    return (0<=i&&i<BLC2LinearTrackContainer.size()) ? &BLC2LinearTrackContainer[i] : 0;
  case CID_BLC2a:
    return (0<=i&&i<BLC2aLinearTrackContainer.size()) ? &BLC2aLinearTrackContainer[i] : 0;
  case CID_BLC2b:
    return (0<=i&&i<BLC2bLinearTrackContainer.size())  ? &BLC2bLinearTrackContainer[i] : 0;
  }
  return 0;
}


//-----------------------Old Clustring-------------------------------------------
bool BeamLineTrackMan::Clustering( BeamLineHitMan *blMan, BLDCClusterMan *clMan, ConfMan *conf, const int &id, const int &id2 )
{
#if 0
  std::cout << "!!! BeamLineTrackMan::Clustering()" << std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  
  bool uxflag, dxflag, uyflag, dyflag;
  uxflag=dxflag=uyflag=dyflag=false;
  
  int maxwire=0;
  if(id==CID_BPC ) maxwire=NumOfBPCWiresInLayer;
  else if(id==CID_BLC1a || id==CID_BLC1b ) maxwire=NumOfBLCWiresInLayer;
  else if(id==CID_BLC2a ||id==CID_BLC2b ) maxwire=NumOfBLCWiresInLayer;

  int maxlayer=0;
  if(id==CID_BPC ) maxlayer=NumOfBPCLayers;
  else if(id==CID_BLC1a || id==CID_BLC1b ) maxlayer=NumOfBLCLayers;
  else if(id==CID_BLC2a ||id==CID_BLC2b ) maxlayer=NumOfBLCLayers;

  int hp[32];
  std::map<int,int> cid;
  int length;
  
  // Search clusters in layer#1 and #2


  for(int icl=0;icl<4;icl++){
    int xy=-1;
    for( int i=0; i<maxwire; i++ ) hp[i]=0;
    for( int j=0; j<2; j++ ){
      int layer=icl*2+j+1;
      if( BLDCParam->layerdead(id,layer) ) continue;
      if( !LAYER[id2][layer-1] ) continue;
      for( int i=0; i<blMan->nBLDC(id,layer); i++ ){
	ChamberLikeHit *hit = blMan->BLDC(id,layer, i);
#if CHECKRANGE
	if(!hit->CheckRange()) continue;
#endif
	int wire = hit->wire();
	//	std::cout<<"cid layer wire \t"<<id<<"\t"<<layer<<"\t"<<wire<<std::endl;
	(hp[wire-1])++;
      }
    }
    
    cid.clear();
    length=0;
    for( int i=0; i<maxwire; i++ ){
      if( hp[i] > 0 ){
	length++;
      }
      else{
	if( length > 0 ){
	  cid.insert( std::map<int,int>::value_type(i-length+1,length) );
	  length=0;
	}
      }
    }
    
    for( std::map<int,int>::iterator it=cid.begin(); it!=cid.end(); it++ ){
      int first_wire = (*it).first;
      int clength = (*it).second;
      int last_wire = first_wire + clength - 1;
      int nhit=0;
      BLDCCluster *cluster = new BLDCCluster();
      for( int j=0; j<2; j++ ){
	int layer=icl*2+j+1;
	if( BLDCParam->layerdead(id,layer) ) continue;
	if( !LAYER[id2][layer-1] ) continue;
	for( int i=0; i<blMan->nBLDC(id,layer); i++ ){
	  ChamberLikeHit *hit = blMan->BLDC(id,layer, i);
#if CHECKRANGE
	  if(!hit->CheckRange()) continue;
#endif
	  if(xy!=-1&&xy!=hit->xy()){
	    std::cout<<"xy layer strange"<<std::endl;
	    return false;
	  }
	  xy=hit->xy();
	  int wire = hit->wire();
	  if( first_wire <= wire && wire <= last_wire ){
	    cluster->SetHit( (*hit) );
	    nhit++;
	  }
	}
      }
      if( nhit < 1 || clength>2 ){
	delete cluster;
	continue;
      }
      if(icl==0||icl==1){
	clMan->SetCluster( 0, xy, (*cluster) );
      }else if(icl==2||icl==3){
	clMan->SetCluster( 1, xy, (*cluster) );
      }
      delete cluster;
    }
  }
  //  return (uxflag && uyflag && dxflag && dyflag);
  return true;
}

// ----------------------New Clustering 2011 10/17 sada--------------------------------------------
bool BeamLineTrackMan::Clustering2( BeamLineHitMan *blMan, BLDCClusterMan *clMan, ConfMan *conf, const int &id )
{
#if 0
  std::cout << "!!! BeamLineTrackMan::Clustering()" << std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();

    
  // Search clusters in layer#1 and #2
  int maxlayer=0;
  if(id==CID_BPC ) maxlayer=NumOfBPCLayers;
  else if(id==CID_BLC1a || id==CID_BLC1b ) maxlayer=NumOfBLCLayers;
  else if(id==CID_BLC2a ||id==CID_BLC2b ) maxlayer=NumOfBLCLayers;

  else{
    std::cerr << "ERROR: ID is not defined in BeamLineTrackMan::Clustering()." << std::endl;
    return false;
  }

  for( int layer=1; layer<=maxlayer; layer++ ){   
    if( BLDCParam->layerkilled(id,layer) ) continue;
    int nBLDCHit=0;
    nBLDCHit=blMan->nBLDC(id,layer);
    if(nBLDCHit> BLDCParam->GetMaxHitInLayer() ) continue;

    for( int i=0; i<nBLDCHit; i++ ){
      ChamberLikeHit *hit=0;

      hit=blMan->BLDC(id,layer,i);
      
      //      int wire = hit->wire();
      int xy=hit->xy();
      double pos=-999;
      if(xy==0) pos=hit->wx(); 
      if(xy==1) pos=hit->wy();

      if(layer%2==1 && BLDCParam->layerkilled(id, layer+1)) 
	{
	  BLDCCluster *cluster = new BLDCCluster();
	  cluster->SetHit( (*hit) );
	  if(xy==0 && layer==1 )clMan->SetClusterUpXZ( (*cluster) );
	  else if(xy==1 && layer==3 )clMan->SetClusterUpYZ( (*cluster) );
	  else if(xy==0 && layer==5 )clMan->SetClusterDownXZ( (*cluster) );
	  else if(xy==1 && layer==7 )clMan->SetClusterDownYZ( (*cluster) );
	  delete cluster;
	  continue;
	}
      else if(layer%2==0 && BLDCParam->layerkilled(id, layer-1)) 
	{
	  BLDCCluster *cluster = new BLDCCluster();
	  cluster->SetHit( (*hit) );
	  if(xy==0 && layer==2 )clMan->SetClusterUpXZ( (*cluster) );
	  else if(xy==1 && layer==4 )clMan->SetClusterUpYZ( (*cluster) );
	  else if(xy==0 && layer==6 )clMan->SetClusterDownXZ( (*cluster) );
	  else if(xy==1 && layer==8 )clMan->SetClusterDownYZ( (*cluster) );
	  delete cluster;
	  continue;
	}

	//###layer2#####
      if(layer%2==0) continue; 
      int layer2=layer+1;
      if( BLDCParam->layerkilled(id, layer2) )  continue;
    
      int nBLDCHit2=0;

      nBLDCHit2=blMan->nBLDC(id,layer2);


      if(nBLDCHit2> BLDCParam->GetMaxHitInLayer() ) continue;

      for( int ii=0; ii<nBLDCHit2; ii++ ){
	BLDCCluster *cluster = new BLDCCluster();
	ChamberLikeHit *hit2=0;

	hit2=blMan->BLDC(id,layer2,ii);

	//	int wire2 = hit2->wire();
	int xy2=hit2->xy();
	double pos2=-1;
	if(xy==0 ||xy2==0) pos2=hit2->wx(); 
	else if(xy==1 ||xy2==1) pos2=hit2->wy();
	else continue;	
	double wiredis=0;
	if(id==CID_BLC1a ||id==CID_BLC1b ) wiredis=0.8;
	else if(id==CID_BLC2a ||id==CID_BLC2b ) wiredis=0.5;
	else if(id==CID_BPC ) wiredis=0.72;

	if( fabs(pos-pos2)<wiredis ) 
	  {
	    cluster->SetHit( (*hit) );
	    cluster->SetHit( (*hit2) );
	    if(xy==0 && layer==1 )clMan->SetClusterUpXZ( (*cluster) );
	    else if(xy==1 && layer==3 )clMan->SetClusterUpYZ( (*cluster) );
	    else if(xy==0 && layer==5 )clMan->SetClusterDownXZ( (*cluster) );
	    else if(xy==1 && layer==7 )clMan->SetClusterDownYZ( (*cluster) );
	    delete cluster;
	  }
	else 
	  {
	    delete cluster;
	  }
      }//Hit layer2
    }//Hit layer1
  }//layer

  return true;
}

//####################################################

bool BeamLineTrackMan::LeastSquareFit( LocalTrack *track, ConfMan *conf,
				       const int &xy )
{
  return track->LeastSquareFit(conf,xy);
}

void BeamLineTrackMan::ConvLocalToGlobal(const int &cid)
{  
  for( int i=0; i<ntrackBLDC(cid); i++ ){
    LocalTrack *track = trackBLDC(cid,i);
    track->ConvLocalToGlobal();
  }  
}
