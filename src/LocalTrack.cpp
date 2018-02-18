// BeamLineTrackMan.cpp
#include <map>

#include "GlobalVariables.h"
#include "BeamLineTrackMan.h"
#include "TVector2.h"
#include "TMath.h"

static const double SpatialResolutionOfBLDC=0.02; // [cm]

ClassImp(LocalTrack);

#define DEBUG 0
#define DEBUG2 0
#define CHECKRANGE 0
#define DEL 1
// ----------------------------- //
// class LocalTrack              //
// ------------------------------//
LocalTrack::LocalTrack()
{
  A=B=C=D=E=F=-999.;
  GA=GB=GC=GD=GE=GF=-999.;
  xzDof=yzDof=0;
  xzChi=yzChi=-999.;

  TrackTime=TrackTimeRMS=TrackTimeX=TrackTimeY=-999.;
}

LocalTrack::~LocalTrack()
{
}

int LocalTrack::nhit(){
  return nhit(0)+nhit(1);
}

int LocalTrack::nhit(const int &xy){
  int tmpnhit=0;
  for(int i=0;i<(int)clusterContainer[xy].size();i++)
    tmpnhit+=this->cluster(xy,i)->nhit();
  return tmpnhit;
}

ChamberLikeHit* LocalTrack::hit( const int &i )
{
  return ( i<(this->nhit(0)) ? (this->hit(0,i)) : (this->hit(1, i - (this->nhit(0)) )) );
}

ChamberLikeHit* LocalTrack::hit( const int &xy, const int &i )
{
  int icl=0;
  int sum=0;
  while(icl<this->ncluster(xy)){
    if((i-sum)<this->cluster(xy,icl)->nhit())
      return this->cluster(xy,icl)->hit(i-sum);
    sum+=this->cluster(xy,icl)->nhit();
    icl++;
  }
  std::cout<<"!!! error in ChamberLikeHit::hit(xy,i) !!!"<<std::endl;
  return 0;
}

void LocalTrack::Clear()
{
  clusterContainer[0].clear();
  clusterContainer[1].clear();
}

int LocalTrack::SetCluster(const int &xy,BLDCCluster *clu){
  clusterContainer[xy].push_back(*clu);
  int tmpnhit=0;
  for(int i=0;i<clu->nhit();i++){
    ChamberLikeHit* hit=clu->hit(i);
    if( hit->status()==1 ) 
      tmpnhit++;
  }
  return tmpnhit;
}

bool LocalTrack::CompareCluster( LocalTrack *tr, const bool &HIT )
{
  bool check=false;
  for(int xy=0;xy<2;xy++){
    for(int i1=0;i1<this->ncluster(xy);i1++){
      int id1=this->cluster(xy,i1)->GetClusterID();
      int cid1=this->cluster(xy,i1)->hit(0)->cid();      
      for(int i2=0;i2<tr->ncluster(xy);i2++){
	int id2=this->cluster(xy,i2)->GetClusterID();
	int cid2=tr->cluster(xy,i2)->hit(0)->cid();
	if(id1==id2&&cid1==cid2){
	  check= true;
	  return true;
	}
	if(cid1==cid2&&HIT){
	  for(int ih1=0;ih1<this->ncluster(xy);ih1++){
	    int lay1=this->cluster(xy,i1)->hit(ih1)->layer();
	    int wire1=this->cluster(xy,i1)->hit(ih1)->wire();
	    for(int ih2=0;ih2<tr->ncluster(xy);ih2++){	      
	      int lay2=tr->cluster(xy,i2)->hit(ih2)->layer();
	      int wire2=tr->cluster(xy,i2)->hit(ih2)->wire();	      
	      if(cid1==cid2&&lay1==lay2&&wire1==wire2){
		check= true;
		return true;
	      }		
	    }
	  }
	}
      }
    }
  }
  return check;
}

bool LocalTrack::CompareTrackHit( LocalTrack *tr )
{
  bool check=false;
  for(int xy=0;xy<2;xy++){
    for(int i1=0;i1<this->nhit(xy);i1++){
      for(int i2=0;i2<tr->nhit(xy);i2++){
	int cid1=this->hit(xy,i1)->cid();      
	int lay1=this->hit(xy,i1)->layer();
	int wire1=this->hit(xy,i1)->wire();
	int cid2=tr->hit(xy,i2)->cid();  
	int lay2=tr->hit(xy,i2)->layer();
	int wire2=tr->hit(xy,i2)->wire();	    
	if(cid1==cid2&&lay1==lay2&&wire1==wire2){
	  check= true;
	  return true;
	}		
      }
    }
  }
  return check;
}

void LocalTrack::CalcTrackTime(){
  int n[2]={0,0};
  double time[20];
  double timexy[2][10];
  for(int xy=0;xy<2;xy++){
    double sumtime=0.;
    for(int i=0;i<this->ncluster(xy);i++){
      double tmptime=this->cluster(xy,i)->GetCTime();
      if(tmptime<-500) continue;
      sumtime+=tmptime;
      time[n[0]+n[1]]=tmptime;
      timexy[xy][n[xy]]=tmptime;
      n[xy]++;      
    }
  }
  if(n[0]>0)
    TrackTimeX=TMath::Mean(n[0],timexy[0]);
  if(n[1]>0)
    TrackTimeY=TMath::Mean(n[1],timexy[1]);
  if((n[0]+n[1])>0){
    TrackTime=TMath::Mean(n[0]+n[1],time);
    TrackTimeRMS=TMath::RMS(n[0]+n[1],time);
  }
  //  std::cout<<this->ncluster(0)+this->ncluster(1)<<"\t"<<TrackTime<<"\t"<<TrackTimeRMS<<std::endl;
}

void LocalTrack::DeleteOffTimingCluster(){
  int n=0;
  double time[20];
  double sumtime=0.;
  int del=-1;
  double tmprms=999.;
  int ncluster=this->ncluster(0)+this->ncluster(1);
  for(int idel=0;idel<ncluster;idel++){
    int icl=0;
    for(int xy=0;xy<2;xy++){
      for(int i=0;i<this->ncluster(xy);i++){       
	icl++;
	if(icl-1==idel) continue;
	double tmptime=this->cluster(xy,i)->GetCTime();
	if(tmptime<-500) continue;
	sumtime+=tmptime;
	time[n]=tmptime;
	n++;      
      }      
    }
    if(n>0){
      double rms=TMath::RMS(n,time);
      if(rms<tmprms) del=idel;
    }
  }
  if(del!=-1){
    if(del<this->ncluster(0)){
      BLDCClusterContainer::iterator it = clusterContainer[0].begin();
      for( int j=0; j<=del; j++ ) ++it;
      clusterContainer[0].erase( it ) ;
    }else{
      BLDCClusterContainer::iterator it = clusterContainer[1].begin();
      for( int j=0; j<=(del-this->ncluster(0)); j++ ) ++it;
      clusterContainer[1].erase( it ) ;
    }
  }  
  CalcTrackTime();
  //  std::cout<<TrackTime<<std::endl;
}

void LocalTrack::DeleteHit(const int &xy, const int &i )
{
  //  if(xy==0) DeleteHitXZ(i);
  //  else if(xy==1) DeleteHitYZ(i);
}

bool LocalTrack::XYLocalPosatZ( const double &z, double &x, double &y, const bool &TILT, const bool &ROT )
{
  if( fabs(A) < 1.0e-9 ) return false;
  if( fabs(D) < 1.0e-9 ) return false;

  double tmpB=B,tmpE=E;
  if(TILT){
    if( this->nhit(0)>0 && this->nhit(1)>0 ){
      ChamberLikeHit *hitx = this->hit(0,0);
      ChamberLikeHit *hity = this->hit(1,0);
      double tiltx= hitx->tilt()*Deg2Rad;
      double tilty= hity->tilt()*Deg2Rad;
      tmpB=-(-B/A+tiltx)*A;
      tmpE=-(-E/D+tilty)*D;
    }
  }
  
  x = -1.*( tmpB*z + C )/A;
  y = -1.*( tmpE*z + F )/D;
  
  if(ROT){
    double rot=0.;
    if( nhit()>0 )    rot = hit(0)->dgz();    
    TVector2 org(x,y);
    TVector2 rotated=org.Rotate(-rot/180.*TMath::Pi());
    x=rotated.X();
    y=rotated.Y();
  }
  return true;
}

bool LocalTrack::XYPosatZ( const double &z, double &x, double &y )
{
  if( fabs(GA) < 1.0e-9 ) return false;
  if( fabs(GD) < 1.0e-9 ) return false;

  double gx=0,gy=0,gz=0;  
  if( this->nhit(0)>0 ){
    ChamberLikeHit *hit = this->hit(0,0);
    gx = hit->gx();  gy = hit->gy();  gz = hit->gz();
    //    std::cout<<"cid "<<hit->cid()<<"  gx,gy,gz: "<<gx<<"\t"<<gy<<"\t"<<gz<<std::endl;
  }

  x = -1.*( GB*(z-gz) + GC )/GA;
  y = -1.*( GE*(z-gz) + GF )/GD;   
  //  x = -1.*( GB*z + GC )/GA;
  //  y = -1.*( GE*z + GF )/GD;
  
  return true;
}

TVector3 LocalTrack::GetPosatZ(const double &z)
{
  double tmpx,tmpy;
  XYPosatZ(z,tmpx,tmpy);
  return TVector3(tmpx,tmpy,z);
}
TVector3 LocalTrack::GetMomDir()
{
  TVector3 tmpdir(-GB/GA,-GE/GD,1);
  return tmpdir.Unit();
}

void LocalTrack::ConvLocalToGlobal()  
{
  //  std::cout<<"convL2G nhit "<<nhit(0)<<std::endl;
  if( this->nhit(0)>0 && this->nhit(1)>0 ){
    ChamberLikeHit *hitx = this->hit(0,0);
    ChamberLikeHit *hity = this->hit(1,0);
    double tiltx= hitx->tilt()*Deg2Rad;
    double tilty= hity->tilt()*Deg2Rad;
    TVector3 pos(-C/A,-F/D,0);
    pos.RotateZ(hitx->dgz()*Deg2Rad);
    pos.RotateX(hitx->dgx()*Deg2Rad);
    pos.RotateY(hitx->dgy()*Deg2Rad);
    TVector3 dir(-B/A+tiltx,-E/D+tilty,1);
    dir.RotateZ(hitx->dgz()*Deg2Rad);
    dir.RotateX(hitx->dgx()*Deg2Rad);
    dir.RotateY(hitx->dgy()*Deg2Rad);
    GA = 1.;
    GB = -dir.X()/dir.Z();
    GC = -(pos.X()-dir.X()/dir.Z()*pos.Z()+hitx->gx());
    GD = 1.; 
    GE = -dir.Y()/dir.Z();
    GF = -(pos.Y()-dir.Y()/dir.Z()*pos.Z()+hitx->gy());
  }
}

void LocalTrack::ConvLocalToGlobal2()  // LinearTrack::ConvLocalToGlobal()
{
  double dgx=0,dgy=0,dgz=0;
  double gx=0,gy=0,gz=0;  
  if( nhit()>0 ){
    ChamberLikeHit *tmphit = hit(0);
    gx = tmphit->gx();  gy = tmphit->gy();  gz = tmphit->gz();
    dgx = tmphit->dgx();  dgy = tmphit->dgy();  dgz = tmphit->dgz();
    TVector3 pos(-C/A,-F/D,0);
    pos.RotateX(dgx*Deg2Rad);
    pos.RotateY(dgy*Deg2Rad);
    TVector3 dir(-B/A,-E/D,1);
    dir.RotateX(dgx*Deg2Rad);
    dir.RotateY(dgy*Deg2Rad);
    GA = 1.; 
    GB = -dir.X()/dir.Z();
    GC = -(pos.X()-dir.X()/dir.Z()*pos.Z()+gx);
    GD = 1.; 
    GE = -dir.Y()/dir.Z();
    GF = -(pos.Y()-dir.Y()/dir.Z()*pos.Z()+gy);
  }
}

bool LocalTrack::Calc( ConfMan *conf, TString option )
{
  for(int xy=0;xy<2;xy++){
    for( int j=0; j<nhit(xy); j++ )
      hit(xy,j)->Calc(conf);    
    LeastSquareFit(conf, xy, option);
  }
  
  if(option.Contains("linear") ){
    LinearFit(conf);
    CalcHitPosition(true);
    CalcResidual(true);
    ConvLocalToGlobal2();
  }
  else{
    CalcHitPosition();
    CalcResidual();
    ConvLocalToGlobal();  
  }
  return true;
}

static const unsigned int LRMASK=0x01;
inline static unsigned int LR( int hid, unsigned int key )
{
  key = key >> hid;
  key = key & LRMASK;
  return key;
}
//####################################################
bool LocalTrack::LeastSquareFit( ConfMan *conf, const int &xy, TString option )
{
#if DEBUG
  std::cout << "!!! LocalTrack::LeastSquareFit()  "<<option << std::endl;
#endif
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  int MaxNumOfHitsInTrack=BLDCParam->GetMaxHitInTrack();
  if( nhit(xy) > MaxNumOfHitsInTrack )
    return false;
  
  double x[MaxNumOfHitsInTrack], y[MaxNumOfHitsInTrack],w[MaxNumOfHitsInTrack];//, theta[MaxNumOfHitsInTrack];
  int np = nhit(xy);
  int np_org = np;
  for( int i=0; i<np_org;i++){
    if(hit(xy,i)->status()!=1){
      np--; continue;
    }    
#if DEBUG
    if(hit(xy,i)->cid()==CID_BPC)
    std::cout<<"cid,layer,wire,dl,wz,wpos\t"<<
      hit(xy,i)->cid()<<"\t"<<hit(xy,i)->layer()<<"\t"<<hit(xy,i)->wire()<<"\t"<<hit(xy,i)->dl()<<"\t"<<hit(xy,i)->wz()<<"\t"<< hit(xy,i)->wpos(xy)<<std::endl;
#endif     

  }
  //  std::cout<<"np,np_org: "<<np<<"\t"<<np_org<<std::endl;
  if( np < 2 ){
#if DEBUG 
    std::cout<<"LocalTrack::LeastSquareFit():  too few hits "<<std::endl;
#endif 
    return false;
  }
  double sigma=SpatialResolutionOfBLDC;
  if(option.Contains("mwpc")) sigma=hit(xy,0)->dxy()/sqrt(12);
  
  double minchi = 1.0e+9;
  double canda=0., candb=0., candc=0.;
  unsigned int candkey=0x0;
  unsigned int key = 0x0;
  int nconb = (int)pow(2,np);
  if(option.Contains("mwpc")) nconb=1;  
  double tmp,tmp2;
  for( int ic=0; ic<nconb; ic++ ){      
    int j=0;
    for( int i=0; i<np_org;i++){
      if( hit(xy,i)->status()!=1 ) continue;
      x[j] = hit(xy,i)->wz();
      y[j] = hit(xy,i)->wpos(xy);
      if(!option.Contains("mwpc")){
	conf->GetReslMapManager()->GetParam(hit(xy,i)->cid(),hit(xy,i)->layer(),hit(xy,i)->wire(),tmp,tmp2);
	w[j]=tmp;
	double dl = hit(xy,i)->dl();
	unsigned int lr = LR(j,key);
	if( lr==0 ) y[j] += dl;
	else y[j] -= dl;
      }else{
	w[j]=sigma;
      }
      //      std::cout<<w[j]<<std::endl;
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
    if( D==0 ){
#if DEBUG 
    std::cout<<"LocalTrack::LeastSquareFit():  invalid determinal "<<std::endl;
#endif 
    return false;
    }
    double aa = (Sxx*Sy - Sx*Sxy)/D;
    double bb = ((double)np*Sxy - Sx*Sy)/D;
    double chi = 0;
    for( int i=0; i<np; i++ ) chi += (y[i]-aa-bb*x[i])*(y[i]-aa-bb*x[i])/(w[i]*w[i]);    
    if(np<3) chi=-1;
    else chi /= (double)(np-2);

    if( (chi < minchi) || (chi==minchi && candb*candb<bb*bb) ){
      canda = 1.;
      candb = -1.*bb;
      candc = -1.*aa;
      candkey=key;
      minchi=chi;	
    }
    key++;
  }
  if(minchi==1.0e+9){
#if DEBUG 
    std::cout<<"LocalTrack::LeastSquareFit():  no track candidate "<<std::endl;
#endif 
    return false;
  }
  if(xy)
    SetDEF(canda,candb,candc);
  else
    SetABC(canda,candb,candc);
  SetDof(xy,np-2);
  SetChisqr(xy,minchi);
  if(option.Contains("mwpc")) return true;
  int j=0;
  for( int i=0; i<np_org; i++ ){
    double dl = hit(xy,i)->dl();
    double zpos= hit(xy,i)->wz(); 
    double wpos = hit(xy,i)->wpos(xy);
    double tmpx,tmpy;
    if( hit(xy,i)->status()!=1 ){
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


bool LocalTrack::LinearFit( ConfMan *conf, const bool &CHECKLR )
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
  double x[MaxNumOfHitsInTrack], y[MaxNumOfHitsInTrack], theta[MaxNumOfHitsInTrack], w[MaxNumOfHitsInTrack];
  int np = nhit();
  int np_org = np;
  int nlr=0;
  int ihitlr[MaxNumOfHitsInTrack];
  for( int i=0; i<np_org;i++){
    ihitlr[i]=-1;
    if( hit(i)->status()!=1 ){
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
  nlr=np;
  if(!CHECKLR) nlr=0;
  int nconb = (int)pow(2,nlr);
  for( int ic=0; ic<nconb; ic++ ){      
    int j=0;
    for( int i=0; i<np_org;i++){
      if( hit(i)->status()!=1 ) continue;
      double dl = hit(i)->dl();
      double temptheta;
      hit(i)->gwpos(y[j],temptheta,x[j],true);      
      //      std::cout<<"z,x,theta "<<x[j]<<"\t"<<y[j]<<"\t"<<temptheta<<std::endl; 
      theta[i]=temptheta*TMath::DegToRad();
      double tmp,tmp2;
      conf->GetReslMapManager()->GetParam(hit(i)->cid(),hit(i)->layer(),hit(i)->wire(),tmp,tmp2);
      w[j]=tmp;
      unsigned int lr =hit(i)->leftright();
      if(CHECKLR){
	lr=LR(j,key);
	//	if(ihitlr[i]>=0)
	//	  lr=LR(ihitlr[i],key);
      }
      if( lr==0 ) y[j] += dl;
      else y[j] -= dl;
      //      std::cout<<"j/np, dl, y = "<<j<<"/"<<np_org<<"\t"<<dl<<"\t"<<y[j]<<std::endl;
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
      //      std::cout<<"trackpos = "<<trackpos<<std::endl;
      //      std::cout<<"i,trackpos,hitpos,chi2\t"<<i<<"\t"<<trackpos<<"\t"<<y[i]<<"\t"<<pow((y[i]-trackpos)/SpatialResolutionOfBLDC,2)<<std::endl;
      chi += pow((y[i]-trackpos)/w[i],2);
    }
    if(np<5) chi=-1;
    else chi /= (double)(np-4);
    //    std::cout<<"ic/ncomb,key,chi="<<ic<<"/"<<nconb<<"\t"<<key<<"\t"<<chi<<std::endl;
    if( (chi < minchi) ){
      canda = par[0][0];
      candb = par[1][0];
      candc = par[2][0];
      candd = par[3][0];
      candkey = key;
      //      std::cout<<"ic,candkey,key="<<ic<<"\t"<<candkey<<"\t"<<key<<std::endl;
      minchi=chi;	
    }
    key++;
  }
  //  std::cout<<"LinearFit::Chi2 = "<<minchi<<" / "<<np-4<<" = "<<minchi/(np-4)<<std::endl;
  SetABC(1.,-candb,-canda);
  SetDEF(1.,-candd,-candc);
  SetDof(0,np-4);
  SetChisqr(0,minchi);
  int j=0;
  if(CHECKLR){
    for( int i=0; i<np_org; i++ ){
      if( hit(i)->status()!=1 ) continue;
      //      if(ihitlr[i]>=0){
      //      unsigned int lr = LR(ihitlr[i],candkey);
      unsigned int lr = LR(j,candkey);
      //      std::cout<<"i,key,j,lr\t"<<i<<"\t"<<candkey<<"\t"<<j<<"\t"<<lr<<std::endl;
      hit(i)->SetLeftRight((int)lr);      
      j++;
    }  
  }
  return true;
}

void LocalTrack::SetGParam(const double &b,const double &c, const double &e, const double &f )
{
  // b: X position at Z=0
  // c: dX/dZ
  // e: Y position at Z=0
  // f: dY/dZ
  double dgx=0,dgy=0,dgz=0;
  double gx=0,gy=0,gz=0;  
  GA = 1;
  GB = -c;
  GC = -b;
  GD = 1;
  GE = -f;
  GF = -e;
  if( nhit()>0 ){
    ChamberLikeHit *tmphit = hit(0);
    gx = tmphit->gx();  gy = tmphit->gy();  gz = tmphit->gz();
    dgx = tmphit->dgx();  dgy = tmphit->dgy();  dgz = tmphit->dgz();
    //    std::cout<<"gx,gy,gz: "<<gx<<"\t"<<gy<<"\t"<<gz<<std::endl;
    TVector3 pos(-GC-gx,-GF-gy,0);
    pos.RotateY(-dgy*Deg2Rad);
    pos.RotateX(-dgx*Deg2Rad);

    
    TVector3 dir(-GB,-GE,1.);
    dir.RotateY(-dgy*Deg2Rad);
    dir.RotateX(-dgx*Deg2Rad);
    A=1;
    B=-(dir.X()/dir.Z());              //dx
    C=-(pos.X()-dir.X()/dir.Z()*pos.Z()); //x
    D=1;
    E=-(dir.Y()/dir.Z());              //dy
    F=-(pos.Y()-dir.Y()/dir.Z()*pos.Z()); //y
  }
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
#if DEBUG2
      std::cout<<"xy,i,z,x,y  "<<xy<<"  "<<i<<"  "<<z<<"  "<<x<<"  "<<y<<std::endl;
#endif
    }
}

void LocalTrack::CalcResidual( const bool &ROT, const TString &option)
{
  //  std::cout<<"----------------------------"<<std::endl;
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
      if(option.Contains("mwpc")){
	tmphit->SetResolution( dltrack );
      }else{
	double dl = tmphit->dl();
	int lr= tmphit->leftright();
	if(lr==0) tmphit->SetResolution( dltrack - dl );
	else tmphit->SetResolution( dltrack + dl );
      }
#if DEBUG2
      std::cout<<"layer,wire,dltrack,dl,resid"
	       <<"\t"<<tmphit->layer()
	       <<"\t"<<tmphit->wire()
	       <<"\t"<<dltrack
	       <<"\t"<<dl
	       <<"\t"<<tmphit->resl()
	       <<std::endl;
#endif
    }
  //  GetCalcChisquare( );
}

double LocalTrack::GetCalcChisquare( )
{
  double chisq=0;
  int ndf=0;
  for( int i=0; i<nhit(); i++ ){      
    if( hit(i)->status()!=1 ) continue;
    chisq += pow(hit(i)->resl()/SpatialResolutionOfBLDC,2);
    ndf++;
  }
  //  std::cout<<"chi2 = "<<chi2all()<<"\t->\t"<<chisq<<" / "<<ndf-4<<" = "<<chisq/(ndf-4)<<std::endl;
  return chisq/(ndf-4);
}

void LocalTrack::Print(){
  std::cout<< "DC_ID:" << this->hit(0)->cid() <<std::endl;
  std::cout<< "x: Chi2,Dof,a,b,c :" << this->chi2xz() << "\t" << this->dofxz() 
	   << "\t" << this->a() << "\t" << this->b() << "\t" << this->c() << std::endl;
  std::cout<< "           ga,b,c :" << this->chi2xz() << "\t" << this->dofxz() 
	   << "\t" << this->ga() << "\t" << this->gb() << "\t" << this->gc() << std::endl;
  std::cout<< "y: Chi2,Dof,a,b,c :" << this->chi2yz() << "\t" << this->dofyz() 
	   << "\t" << this->d() << "\t" << this->e() << "\t" << this->f() << std::endl;
  std::cout<< "           ga,b,c :" << this->chi2yz() << "\t" << this->dofyz() 
	   << "\t" << this->gd() << "\t" << this->ge() << "\t" << this->gf() << std::endl;
  std::cout<<"------------------------------------------------------"<<std::endl;
}
