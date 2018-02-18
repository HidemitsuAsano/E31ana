#include "CDSTrack.h"

ClassImp(CDSTrack);
#define DEBUG 0
const double ResolutionOfCDC=0.02;

CDSTrack::CDSTrack() : TObject()
{
  Clear();
}

void CDSTrack::Clear()
{
  for(int i=0;i<15;i++) {trackContainer[i].clear();}
  goodflag=false;
  CDHflag=false;
  fittinglevel=0;
  for(int i=0;i<5;i++) param[i]=0;
  Line_o.SetXYZ(0,0,0);  Line_d.SetXYZ(0,0,0);
}


void CDSTrack::GetParameters(double *aparam)
{
  for(int i=0;i<5;i++) aparam[i]=param[i];
}

void CDSTrack::GetGParameters(double *aparam)
{
  for(int i=0;i<5;i++) aparam[i]=param[i];
}


void CDSTrack::GetTmpParameters(const int &num,double *aparam)
{
  if(0<=num&&num<5)
    {
      for(int i=0;i<5;i++) aparam[i]=paramtmp[num][i];
    }
}


void CDSTrack::SetParameters(const double *aparam)
{
  for(int i=0;i<5;i++) param[i]=aparam[i];
}


void CDSTrack::XYatZ( double &x, double &y, const double &z )
{

  if(param[2]!=0)
    {
      double phi=param[2]/param[4]*(param[3]-z);
      x = (param[0]+1./param[2])*cos(param[1])-1./param[2]*cos(param[1]+phi);
      y = (param[0]+1./param[2])*sin(param[1])-1./param[2]*sin(param[1]+phi);
    }
  else if(param[2]==0)
    {
      double time =-(param[3]-z)/param[4];
      x=param[0]*cos(param[1])-time*(-sin(param[1]));
      y=param[0]*sin(param[1])-time*(cos(param[1]));
    }
}

void CDSTrack::XYtoZ( const double &x,const  double &y,double &z )
{
  double phi=CalcHelixPhi(x,y);
  z=param[3]-1./param[2]*param[4]*phi;
}


void CDSTrack::AddHit(const CDCHit &hit)
{
  int cid = hit.cid();
  int layer = hit.layer();
#if 0
  std::cout << " in Add Hit, layer:" << hit.layer() << " wire:" << hit.wire() <\< " CID:" << cid << std::endl;
#endif

  if(cid ==CID_CDC)
    {
      trackContainer[layer-1].push_back(hit);
      trackContainer[layer-1][ trackContainer[layer-1].size()-1 ].SetHitID( trackContainer[layer-1].size()-1 );
#if 0
      std::cout << "?? wire:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].wire()
	<< " hid:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].hid()
	<< " slayer:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].slayer()
	  << std::endl;
#endif
    }
}


void CDSTrack::SetCDHHit(const HodoscopeLikeHit &cdhhit)
{
  int cid = cdhhit.cid();

  if(cid ==CID_CDH)
    {
      CDHhit=cdhhit;
      CDHflag=true;
    }

}



bool CDSTrack::RemoveAllHitInLayer( const int &layer ) 
{
  if(!(1<=layer && layer<=15)) return false;


  for(int i=0;i < (int)trackContainer[layer-1].size();i++) trackContainer[layer-1].pop_back();
  return true;

}

bool CDSTrack::DeleteHit( const int &layer ,const int &i) 
{
  if(!(1<=layer && layer<=15)) return false;

  TrackHitContainer::iterator it = trackContainer[layer-1].begin();
  for( int j=0; j<i; j++ ) ++it;
  trackContainer[layer-1].erase( it ) ;
  return true;
}


bool CDSTrack::FirstCircleFitting()
{
    
    double x[40],y[40],weight[40];
    int numallhit=0;

    for(int numlayer=0;numlayer<7;numlayer++)
      {
	int layer=0;
	if(numlayer<3) layer=numlayer+1;
	else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
	else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
	for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	  {
	    double dl=trackContainer[layer-1][n].dl();
	    if(dl<ResolutionOfCDC) dl=ResolutionOfCDC;

	    x[numallhit]=trackContainer[layer-1][n].wx();
	    y[numallhit]=trackContainer[layer-1][n].wy();
	    weight[numallhit]=dl;
	    numallhit++;
	  }
      }
          
    CircleFit *cirfit=new CircleFit(param,x,y,weight,numallhit);
    cirfit->GetParameters(param);
    cirfit->GetParameters(paramtmp[1]);

    cirfit->CalcChi2();

    ChiSqr= cirfit->chisquare();
    dof=cirfit->dof();
    delete cirfit;
    ChiSqr=(double)ChiSqr/dof;
    ChiSqrtmp[0]=ChiSqr;
    if(!cirfit->stat()) return false;
    if( ChiSqr >99999 ) 
      {
#if DEBUG
	std::cout<<"Too big chi2/dof f"<<std::endl;
#endif
	return false;
      }
    if(param[1]>2*TMath::Pi() ) param[1]-=2*TMath::Pi();
    if(param[2]==0) return false;
    fittinglevel=1;
    CirRho=fabs(1./param[2]);
    CirCenterX=(param[0]+1./param[2])*cos(param[1]);
    CirCenterY=(param[0]+1./param[2])*sin(param[1]);
    return true;
}



bool CDSTrack::FirstHelixFitting()
{

  double disZ[40];
  int numstlayer=0;

  for(int stlayer=0 ;stlayer<8;stlayer++)
    {
      int layer=0;
      if(0<=stlayer && stlayer<=3) layer=stlayer+4;
      else if(4<=stlayer && stlayer<=7) layer=stlayer+6;

      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	{
	  CDCHit *cdc=&trackContainer[layer-1][n];
	  TVector3 wpos,wposp;
	  
	  wpos.SetXYZ(trackContainer[layer-1][n].wx(),
		      trackContainer[layer-1][n].wy(),
		      trackContainer[layer-1][n].wz() );
	  
	  wposp.SetXYZ(trackContainer[layer-1][n].wxp(),
		       trackContainer[layer-1][n].wyp(),
		       trackContainer[layer-1][n].wzp() );

	  double aw,bw,cw;
	  
	  if((wposp.x()-wpos.x() )!=0)
	    {
	      aw=(wposp.y()-wpos.y() )/(wposp.x()-wpos.x() );
	      bw=-1;
	      cw= -(aw*wpos.x()+bw*wpos.y() );
	    }
	    else
	      {
		bw=(wposp.x()-wpos.x() )/(wposp.y()-wpos.y() );
		aw=-1;
		cw= -(aw*wpos.x()+bw*wpos.y() );
	      }

	    double x_p,x_n,y_p,y_n;
	    TVector3 WPos;
	    if(!LineToCircle(aw,bw,cw,CirRho,CirCenterX,CirCenterY,x_p,y_p,x_n,y_n) )
	      continue;
 	    if(sqrt( (x_p-wpos.x())*(x_p-wpos.x())
		     +(y_p-wpos.y())*(y_p-wpos.y()))<
	       sqrt( (x_n-wpos.x())*(x_n-wpos.x())
		     +(y_n-wpos.y())*(y_n-wpos.y())) ){ WPos.SetX(x_p);WPos.SetY(y_p);}
	    else{ WPos.SetX(x_n);WPos.SetY(y_n);}

	    double dis=sqrt( ( wpos.x()-WPos.x() )*(wpos.x()-WPos.x())
			     +(wpos.y()-WPos.y())*(wpos.y()-WPos.y() ) );
	    double wdis=sqrt( (wpos.x()-wposp.x() )*(wpos.x()-wposp.x() )+(wpos.y()-wposp.y() )*(wpos.y()-wposp.y() ) );
	    WPos.SetZ( wpos.z()+(wposp.z()-wpos.z() )*dis/wdis );
	    trackContainer[layer-1][n].SetHitPosition(WPos.x(),WPos.y(),WPos.z());    

	    //########dis Z ###########//

	    double dl=cdc->dl();
	    double tilt=cdc->tilt()*TMath::DegToRad();
	    disZ[numstlayer]=fabs(dl/sin(tilt) );
	    numstlayer++;
	}

    }
  
  int checki=-1;	  
  for(int i =0;i<pow(2,numstlayer);i++)
    {
      double hitx[256],hity[256],hitz[256];
      double phi[256];
      //      double weight[256];

      int allhit=0;
      int numst=0;
      double s_x=0,s_y=0,s_xx=0,s_xy=0;
      
      for(int layer=1;layer<=15;layer++)
	{
	  for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	    {
	      CDCHit *cdc=&trackContainer[layer-1][n];
	      
	      if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
		{
		
		  hitx[allhit]= cdc->x();
		  hity[allhit]= cdc->y();
		  hitz[allhit]= cdc->z();
		  phi[allhit]=CalcHelixPhi(hitx[allhit],hity[allhit]);

		  phi[allhit]=phi[allhit]/param[2];
		  int check=i;
		  
		  if(check<pow(2,numst) ) hitz[allhit]+=disZ[numst];
		  else 
		    {
		      int cflag;
		      for(int n=0;n<numst;n++)
			{
			  cflag=check%2;
			  check=(check-cflag)/2;
			}
		      cflag=check%2;
		      if(cflag==0)  hitz[allhit]+=disZ[numst];
		      else if(cflag==1) hitz[allhit]-=disZ[numst];
		    }

		  s_x+=phi[allhit];
		  s_y+=hitz[allhit];
		  s_xx+=phi[allhit]*phi[allhit];
		  s_xy+=phi[allhit]*hitz[allhit];

		  numst++;
		  allhit++;
		  
		}
	      
	    }
	}
      
      
      double aa,bb,cc;
      aa=(allhit*s_xy-s_x*s_y)/(allhit*s_xx-s_x*s_x);
      cc=(s_xx*s_y-s_xy*s_x)/(allhit*s_xx-s_x*s_x);
      bb=-1;
      
      double Chitmp=0;
      double dis=0;
      for(int n=0;n<allhit;n++)
	{
	  dis=fabs(aa*phi[n]+bb*hitz[n]+cc)/sqrt(aa*aa+bb*bb);
	  Chitmp+=dis*dis;
	}
      Chitmp=(double)Chitmp/allhit;      

      if(i==0)
	{
	  param[3]=cc;param[4]=-aa;
	  ChiSqr= Chitmp;
	  checki=i;	  
	}	
      else if(Chitmp<ChiSqr)
	{
	  param[3]=cc;param[4]=-aa;
	  ChiSqr= Chitmp;
	  checki=i;	  
	}
      
    }
  
  int numst=0;
  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	{
	  if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
	    {
		  double hitx= trackContainer[layer-1][n].x();
		  double hity= trackContainer[layer-1][n].y();
		  double hitz= trackContainer[layer-1][n].z();
		  int check=checki;
		  if(check<pow(2,numst) )
		    hitz+=disZ[numst];
		  else 
		    {
		      int cflag;
		      for(int n=0;n<numst;n++)
			{
			  cflag=check%2;
			  check=(check-cflag)/2;
			}
		      cflag=check%2;
		      if(cflag==0)  hitz+=disZ[numst];
		      else if(cflag==1) hitz-=disZ[numst];
		    }
		  trackContainer[layer-1][n].SetHitPosition(hitx,hity,hitz);   
		  numst++;
	    }

	}
    }
  ChiSqrtmp[2]=ChiSqr;
  for(int n=0;n<5;n++) paramtmp[3][n]=param[n];
    if( ChiSqr >999999 ) 
      {
#if DEBUG
	std::cout<<"Too big chi2/dof fh"<<std::endl;
#endif
	return false;
      }
    fittinglevel=3;
    SetHitPos();      
    return true;
}



bool CDSTrack::CircleFitting()
{
    
    double x[40],y[40],weight[40];
    int numallhit=0;

    for(int numlayer=0;numlayer<7;numlayer++)
      {
	int layer=0;
	if(numlayer<3) layer=numlayer+1;
	else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
	else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
	for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	  {
	    double dl=trackContainer[layer-1][n].dl();
	    if(dl<0) dl=0.0001;

	    x[numallhit]=trackContainer[layer-1][n].wx();
	    y[numallhit]=trackContainer[layer-1][n].wy();

	    double rtmp=sqrt( (x[numallhit]-CirCenterX)* (x[numallhit]-CirCenterX)+(y[numallhit]-CirCenterY)*(y[numallhit]-CirCenterY) );
	    
	    double sintmp= (y[numallhit]-CirCenterY)/rtmp;
	    double costmp= (x[numallhit]-CirCenterX)/rtmp;

	    double sign;
	    if(rtmp<fabs(1./param[2])) sign=1;
	    else  sign=-1;
	    
	    x[numallhit]+= sign*dl*costmp;
	    y[numallhit]+= sign*dl*sintmp;
	    trackContainer[layer-1][n].SetHitPosX(x[numallhit]);
	    trackContainer[layer-1][n].SetHitPosY(y[numallhit]);
	    
	    weight[numallhit]=ResolutionOfCDC;
	    numallhit++;
	  }
      }
          
    CircleFit *cirfit=new CircleFit(param,x,y,weight,numallhit);
    cirfit->GetParameters(param);

    if(param[1]<0) param[1]+=2*TMath::Pi();
    else if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi();


    cirfit->CalcChi2();
    ChiSqr= cirfit->chisquare();
    dof=cirfit->dof();
    if(!cirfit->stat()) return false;
    delete cirfit;
    ChiSqr=(double)ChiSqr/dof;
     ChiSqrtmp[1]=ChiSqr;
      for(int n=0;n<5;n++) paramtmp[2][n]=param[n];
    if( ChiSqr >99999 ) 
      {
#if DEBUG
	std::cout<<"Too big chi2/dof t"<<std::endl;
#endif
	return false;
      }

    CirCenterX=(param[0]+1./param[2])*cos(param[1]);
    CirCenterY=(param[0]+1./param[2])*sin(param[1]);
    CirRho=fabs(1./param[2]);
     //################################//
      //###########Check Charge ##########//	
      //################################//
    MathTools *tool=new MathTools();
    double phi_c=tool->CalcDeg(CirCenterX,CirCenterY);
    delete tool;
    double phi_mid=theta;
    if(phi_mid>270 &&phi_c<90 ) phi_mid-=360;
      else if(phi_c>270 &&phi_mid<90 ) phi_c-=360;
    if( (phi_mid-phi_c)<0 && param[2]>0) 
      {param[0]=-param[0];param[2]=-param[2];param[1]+=TMath::Pi();}
    else if( phi_mid-phi_c>0 && param[2]<0) 
      {param[0]=-param[0];param[2]=-param[2];param[1]-=TMath::Pi();}
    if(param[1]<0) param[1]+=2*TMath::Pi();
    else if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi();
 

    for(int n=0;n<5;n++) paramtmp[2][n]=param[n];
    fittinglevel=2;
    SetHitPos();
    return true;
}

bool CDSTrack::HelixFitting()
{
   
  int allhit=0;
  
  double hitx[50],hity[50],hitz[50];
  double weight[50];

  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	  {
	    CDCHit *cdc=&trackContainer[layer-1][n];

	    hitx[allhit]= cdc->x();
	    hity[allhit]= cdc->y();
	    if( (1<=layer &&layer<=3) || (8<=layer && layer<=9) || (14<=layer && layer<=15) )
	      {
		hitz[allhit]= -999;
		weight[allhit]=ResolutionOfCDC;
	      }
	    else 
	      {
		hitz[allhit]= cdc->z();
		weight[allhit]=ResolutionOfCDC;

	      }
	    allhit++;
	  }
      
    }

  HelixFit *helixfit=new HelixFit(param,hitx,hity,hitz,weight,allhit);
  helixfit->GetParameters(param); 
  
  CirCenterX=(param[0]+1./param[2])*cos(param[1]);
  CirCenterY=(param[0]+1./param[2])*sin(param[1]);
  CirRho=1./param[2];   
  
  ChiSqr= helixfit->chisquare();
  dof= helixfit->dof();
  ChiSqr=(double)ChiSqr/dof;
  for(int n=0;n<5;n++) paramtmp[4][n]=param[n];
  ChiSqrtmp[3]=(double)ChiSqr/dof;
  // if(!helixfit->stat()) return false;
  delete helixfit;
    if( ChiSqr >9999 ) 
      {
#if DEBUG
	std::cout<<"Too big chi2/dof h"<<std::endl;
#endif
	return false;
      }

    fittinglevel=4;
    SetHitPos();
  
    return true;

}


void CDSTrack::SetHitPos()
{
    
    for(int layer=1;layer<=15;layer++)
      {
	for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	  {
	    double dl=trackContainer[layer-1][n].dl();
	    if(dl<0) dl=0.002;

	    TVector3 whpos,hpos,wpos,wposp,dline;
	    wpos.SetX(trackContainer[layer-1][n].wx() );
	    wpos.SetY(trackContainer[layer-1][n].wy() );
	    wpos.SetZ(trackContainer[layer-1][n].wz() );
	    wposp.SetX(trackContainer[layer-1][n].wxp() );
	    wposp.SetY(trackContainer[layer-1][n].wyp() );
	    wposp.SetZ(trackContainer[layer-1][n].wzp() );
	    hpos.SetX(trackContainer[layer-1][n].x() );
	    hpos.SetY(trackContainer[layer-1][n].y() );
	    hpos.SetZ(trackContainer[layer-1][n].z() );
	    
	    dline=wposp-wpos;
	    dline=dline.Unit();

	    double k=hpos*dline-wpos*dline;
	    whpos=wpos+k*dline;
	    
	    TVector3 lnest,hnest;
	    TVector3 hitpos;
	    double dis,resi;
	    if(param[2]==0)
	      {		
		LineToLine(Line_o,Line_d,wpos,dline,dl,dis,lnest,hitpos);
	      }
	    else 
	      {
		if(!LineToHelix(whpos,dline,param,lnest,hnest,dis) )
		  {
#if DEBUG
		    std::cout<<"Miss LineToHelix!!"<<std::endl;
#endif
		    continue;
		  }
		TVector3 tmp;
		tmp=hnest-lnest;
		tmp=tmp.Unit();
		hitpos=lnest+dl*tmp;
	      }
	    resi=(dl-dis);
	    trackContainer[layer-1][n].SetResolution(resi);
	    trackContainer[layer-1][n].SetHitPosX(hitpos.x());
	    trackContainer[layer-1][n].SetHitPosY(hitpos.y());
	    trackContainer[layer-1][n].SetHitPosZ(hitpos.z());
	    
	  }
      }

}

bool CDSTrack::LineToCircle(const double &a,const double &b, const double &c, const double &rho ,const double &xc, const double &yc, double &x_p,double &y_p,double &x_n,double &y_n)
{
  double dis=fabs(a*xc+b*yc+c)/sqrt(a*a+b*b);
  if(dis>rho){
#if DEBUG
    std::cout<<"no crossing point! rho dis "<<rho<<" "<<dis<<" "<<c<<std::endl;
#endif 
    return false;}

  double p[3];

  if(b!=0)
    {
      p[0]=a*a/(b*b)+1;
      p[1]=2*a*c/(b*b)-2*xc+2*yc*a/b;
      p[2]=xc*xc+c*c/(b*b)+2*yc*c/b+yc*yc-rho*rho;
      
      x_p=(-p[1]+sqrt(p[1]*p[1]-4*p[0]*p[2]) )/(2*p[0]);
      x_n=(-p[1]-sqrt(p[1]*p[1]-4*p[0]*p[2]) )/(2*p[0]);
      y_p= -(a*x_p+c)/b;
      y_n= -(a*x_n+c)/b;
    }
  else if(a!=0)
    {
      p[0]=b*b/(a*a)+1;
      p[1]=2*b*c/(a*a)-2*yc+2*xc*b/a;
      p[2]=yc*yc+c*c/(a*a)+2*xc*c/a+xc*xc-rho*rho;
      
      y_p=(-p[1]+sqrt(p[1]*p[1]-4*p[0]*p[2]) )/(2*p[0]);
      y_n=(-p[1]-sqrt(p[1]*p[1]-4*p[0]*p[2]) )/(2*p[0]);
      x_p= -(b*y_p+c)/a;
      x_n= -(b*y_n+c)/a;

    }
  double dis2=sqrt( (x_p-xc)*(x_p-xc)+(y_p-yc)*(y_p-yc) );
  if(fabs(dis2-rho)>0.01) {std::cout<<"not macth rho! missing calc!"<<std::endl; return false;}
  return true;
}



double CDSTrack::dfunc_LTH(const TVector3 &lpos,const TVector3 &dline,const double &helixphi,const double *par)
{

  TVector3 hpos,dhpos;
  hpos=GetPosition(helixphi,par);
  dhpos.SetXYZ( (1./par[2])*sin(par[1]+helixphi)  ,
		(-1./par[2])*cos(par[1]+helixphi) ,
		(-1./par[2])*par[4]                );

  double k=hpos*dline-lpos*dline;
  double dk=dhpos*dline;

  double S = 2.*(lpos+k*dline-hpos)*(dk*dline-dhpos);
  return S;
}


bool CDSTrack::LineToHelix(const TVector3 &a, const  TVector3 &dline ,
			   const double *par, TVector3 &lnest,TVector3 &hnest,
			   double &dis)
{

  TVector3 dlineu=dline.Unit();
  double phi=0,phi_b=0,phi_a=0;
  phi=CalcHelixPhi(a.x(),a.y(),par);

  int max=0;
  
  while(max<4)
    {
      double z=par[3]-1./par[2]*par[4]*phi;
      double z_b=par[3]-1./par[2]*par[4]*(phi-2*TMath::Pi());
      double z_a=par[3]-1./par[2]*par[4]*(phi+2*TMath::Pi());
      if(fabs(z-a.z())<fabs(z_b-a.z()) && fabs(z-a.z())<fabs(z_a-a.z())) break;
      else if(fabs(z-a.z())>fabs(z_b-a.z())) phi-=2*TMath::Pi();
      else if(fabs(z-a.z())>fabs(z_a-a.z())) phi+=2*TMath::Pi();
      max++;
    }
  
  double philen=1./par[2]*sqrt(1+par[4]*par[4]);
  double dist=1.;
  int trial=14; //sigma_position<1.2 micron 
  while(dist<128)
    {
      phi_b=phi-dist/philen;
      phi_a=phi+dist/philen;
      double dlen_b=dfunc_LTH(a,dlineu,phi_b,par);
      double dlen_a=dfunc_LTH(a,dlineu,phi_a,par);
      if(dlen_b*dlen_a<=0) break;
      else {dist*=2;trial++;}
    }
  
  if(dis>=128) 
    {
#if DEBUG
      std::cout<<"Can not find LTH inital param!!"<<std::endl;
#endif
    return false;
    }
 
 //Bisection Method
  for(int i=0;i<trial;i++)
    {
      phi_b=phi-dist/philen;
      phi_a=phi+dist/philen;
      double dlen_b=dfunc_LTH(a,dlineu,phi_b,par);
      //      double dlen_a=dfunc_LTH(a,dlineu,phi_a,par);
      double dlen=dfunc_LTH(a,dlineu,phi,par);
      if(dlen*dlen_b<=0 ) {phi=(phi_b+phi)/2.0;dist=dist/2.0;}
      else  {phi=(phi_a+phi)/2.0;dist=dist/2.0;}      
    }
  //
  
  hnest=GetPosition(phi,par);
  double k=hnest*dlineu-a*dlineu;
  lnest=a+k*dlineu;
  TVector3 tmp;
  tmp=hnest-lnest;
  dis=tmp.Mag();
  return true;

}


bool CDSTrack::GetMomentum(const TVector3 &pos,TVector3 &p)
{
  double phi=CalcHelixPhi(pos.x(),pos.y());
  TVector3 hpos=GetPosition(phi);
  if(sqrt( (pos.x()-hpos.x())*(pos.x()-hpos.x())+(pos.y()-hpos.y())*(pos.y()-hpos.y()) )>5 ) 
  {
    std::cout<<"$E GetMomentum() check input pos!"<<std::endl;
    p.SetXYZ(-999,-999,-999);
    return false;
  }
  int trial=0;
  while(fabs(pos.z()-hpos.z() )>3.0)
    {
      trial++;
      double tmpphi=phi-2*trial*TMath::Pi();
      hpos=GetPosition(tmpphi);
      if( fabs(pos.z()-hpos.z() )>3.0 )
	{
	  tmpphi=phi+2*trial*TMath::Pi();
	  hpos=GetPosition(tmpphi);
	}
      if(trial>5 && fabs(pos.z()-hpos.z() )>3.0)
	{
	  std::cout<<"$E3 GetMomentum() pos "<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<" "<<" hpos "<<hpos.x()<<" "<<hpos.y()<<" "<<hpos.z()<<std::endl;
	  p.SetXYZ(-999,-999,-999);
	  return false;
	}
    }

  //  const double dMagneticField = -1*0.5; //T, "-1" is needed.
  double pt = fabs(1/param[2]*(Const*magneticfield)/100.);
  double px = pt*(-1*sin(param[1]+phi));
  double py = pt*(cos(param[1]+phi));
  double pz = pt*(param[4]);
  p.SetXYZ(px,py,pz);
  return true;

}


TVector3 CDSTrack::GetPosition(const double &helixphi)
{
  TVector3 pos;
  pos.SetXYZ(param[0]*cos(param[1])+1./param[2]*( cos(param[1])-cos(param[1]+helixphi) ),
	     param[0]*sin(param[1])+1./param[2]*( sin(param[1])-sin(param[1]+helixphi) ),
	     param[3]-1./param[2]*param[4]*helixphi);

  return pos;
}

TVector3 CDSTrack::GetPosition(const double &helixphi,const double *par)
{
  TVector3 pos;
  pos.SetXYZ(par[0]*cos(par[1])+1./par[2]*( cos(par[1])-cos(par[1]+helixphi) ),
	     par[0]*sin(par[1])+1./par[2]*( sin(par[1])-sin(par[1]+helixphi) ),
	     par[3]-1./par[2]*par[4]*helixphi);

  return pos;
}


double CDSTrack::dfunc_PTH(const TVector3 &pos,const double &helixphi,const double *par)
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

bool CDSTrack::PointToHelix(const TVector3 &hitpos, const double *par,TVector3 &fitpos ,double &dis)
{

  double phi=0,phi_b=0,phi_a=0;
  phi=CalcHelixPhi(hitpos.x(),hitpos.y(),par);

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
      std::cout<<"Can not find PTH inital param!!"<<std::endl;
#endif
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

  fitpos=GetPosition(phi,par);
  TVector3 tmp;
  tmp=fitpos-hitpos;
  dis=tmp.Mag();
  return true;
  
}

void CDSTrack::PointToCircle(const double &x,const double &y,const double &radius,const double &x_cen,const double &y_cen,double &dis,double &xest,double &yest)
{
  double rtmp=sqrt( (x-x_cen)*(x-x_cen)+(y-y_cen)*(y-y_cen) );
  double costmp=(x-x_cen)/rtmp;
  double sintmp=(y-y_cen)/rtmp;
  dis=fabs(rtmp-radius);

  xest=x_cen+radius*costmp;
  yest=y_cen+radius*sintmp;
}

void CDSTrack::ReconstructHit(ConfMan *Conf)
{
  int max_w[15]={72,72,72,
		 90,90,100,100,
		 120,120,150,150,
		 160,160,180,180 };

  for(int layer=1;layer<=15;layer++)
    {
      CDCHit hit_r;
      double dis=-1;
      CDCHit hit;
      for(int wire=1;wire<=max_w[layer-1];wire++)
	{
	  hit.Clear();
	  hit.SetCounterID( CID_CDC );
	  hit.SetLayer( layer );
	  hit.SetWire( wire );
	  hit.SetDriftLength( 0.01 );
	  hit.SetTimeOffsetTrue( 0 );
	  
	  if( Conf!=0 ){
	    if( Conf->GetCDCWireMapManager() ){
	      int c,n,a;
	      int slayer, asdnum, ttype, asdch;
	      double rad, phi, tilt;
	      Conf->GetCDCWireMapManager()->GetWire( layer, wire, slayer, asdnum, ttype, asdch,
						     rad, phi, tilt, c, n, a );
	      hit.SetCrate(c); hit.SetSlot(n); hit.SetChannel(a);
	      hit.SetSuperLayer(slayer); hit.SetASDNum(asdnum); hit.SetTransType(ttype); hit.SetASDChannel(asdch);
	      hit.SetRadius(rad); hit.SetPhi(phi); hit.SetTiltAngle(tilt);
	
	    }
	    else{
	      std::cout << " cdc wire map was not found " << std::endl;
	    }
	  }
	  hit.Calc(Conf);
	   
	    //##########
	  double phidir=0;
	  if(param[2]<0) phidir=param[1]+3.151592/2.;
	  else if(param[2]>0) phidir=param[1]-3.151592/2.;

	  if(phidir<0) phidir+=2*3.141592;
	  else if(phidir>2*TMath::Pi()) phidir-=2*3.141592;
	  double def_deg=fabs( hit.phi()*3.141592/180. - phidir );
	  if(def_deg>2*TMath::Pi()-1.0) {def_deg-=2*TMath::Pi();def_deg=fabs(def_deg);}
	  if(def_deg>3.141592/2.)  continue; 
	    

	  if( (1<=layer &&layer<=3) || (8<=layer &&layer<=9) || (14<=layer &&layer<=15) )	 
	    {
	      TVector3 wpos;
	      wpos.SetXYZ(hit.wx(),hit.wy(),hit.wz() );
	      double distmp,xest,yest;
	      PointToCircle(wpos.x(),wpos.y(),CirRho,CirCenterX,CirCenterY,distmp,xest,yest);
	      /*
	      if(dis==-1) { dis=distmp;hit.SetDriftLength( dis ); 
	      hit.Reverse( Conf ); hit_r=hit;}
	      else if(dis>distmp) { dis=distmp;hit.SetDriftLength( dis ); 
	      hit.Reverse( Conf ); hit_r=hit;}*/
	      if(distmp<0.95){ dis=distmp;hit.SetDriftLength( dis ); 
	      hit.Reverse( Conf ); hit_r=hit;AddHit(hit_r); }
	    }


	  else if( (4<=layer && layer<=7) || (10<=layer&&layer<=13) )
	    {
	      TVector3 whpos,wpos,wposp,dline;
	      wpos.SetXYZ(hit.wx(),hit.wy(),hit.wz() );
	      wposp.SetXYZ(hit.wxp(),hit.wyp(),hit.wzp() );

	      double distmp,xest,yest;
	      PointToCircle(wpos.x(),wpos.y(),CirRho,CirCenterX,CirCenterY,distmp,xest,yest);
	      if(distmp>10) continue;
	      dline=wposp-wpos;
	      dline=dline.Unit();

	      TVector3 lnest,hnest;
	      //	      double distmp;
	      if(!LineToHelix(wpos,dline,param,lnest,hnest,distmp) )
		{std::cout<<"Miss LineToHelix!!"<<std::endl;continue;}
	      /*
	      if(dis==-1) { dis=distmp;hit.SetDriftLength( dis ); 
	      hit.Reverse( Conf ); hit_r=hit; }
	      else if(dis>distmp) { dis=distmp;hit.SetDriftLength( dis ); 
	      hit.Reverse( Conf ); hit_r=hit;}
	      */
	      if(distmp<0.95){ dis=distmp;hit.SetDriftLength( dis ); 
	      hit.Reverse( Conf ); hit_r=hit;AddHit(hit_r); }
	    }

	  //  delete hit;
	}
     
      //      AddHit(hit_r);
      //      delete hit_r;
    }

  //############ Set Some Parameter ############//
  int in_outnum=0;
  int midnum=0;
  TVector3 in_outpos[10];
  TVector3 midpos[10];
  for(int layer=1;layer<=15;layer++)
    {

      if(trackContainer[layer-1].size()==0) continue;
      if(trackContainer[layer-1][0].slayer()==1 ||trackContainer[layer-1][0].slayer()==7)
	{
	  in_outpos[in_outnum].SetXYZ(trackContainer[layer-1][0].wx(),
				      trackContainer[layer-1][0].wy(),
				      trackContainer[layer-1][0].wz());

	    in_outnum++;
	}

      if(trackContainer[layer-1][0].slayer()==4 )
	{

	  midpos[midnum].SetXYZ(trackContainer[layer-1][0].wx(),
				      trackContainer[layer-1][0].wy(),
				      trackContainer[layer-1][0].wz());

	    midnum++;
	}
           
    }
  double s_xy=0,s_x=0,s_xx=0,s_y=0;
  for(int n=0;n<in_outnum;n++)
    {
      s_x+=in_outpos[n].x();
      s_y+=in_outpos[n].y();
      s_xx+=in_outpos[n].x()*in_outpos[n].x();
      s_xy+=in_outpos[n].x()*in_outpos[n].y();
    }
  Al=(in_outnum*s_xy-s_x*s_y)/(in_outnum*s_xx-s_x*s_x);
  Cl=(s_xx*s_y-s_xy*s_x)/(in_outnum*s_xx-s_x*s_x);
  Bl=-1;
  s_x=0;
  s_y=0;
  for(int n=0;n<midnum;n++)
    {
      s_x+=midpos[n].x();
      s_y+=midpos[n].y();
    }

  Clusterx[3]=s_x/midnum;
  Clustery[3]=s_y/midnum;
  theta=atan(Clustery[3]/Clusterx[3]);
  jitter=fabs(Al*Clusterx[3]+Bl*Clustery[3]+Cl)/(Al*Al+Bl*Bl);

}

bool CDSTrack::CalcInitialParameters()
{
 
  if(fabs(jitter)<1e-100 ) return false;
  double rho=fabs( (15*15+jitter*jitter)/(2*jitter) );
  
  double ap=Bl/Al;
  double bp=-1;
  double cp=Clustery[3]-ap*Clusterx[3];
  double MidXp,MidYp;
  MidXp=-(Cl/Bl-cp/bp)/(Al/Bl-ap/bp);
  MidYp=ap*MidXp+cp;


  TVector3 dirtmp;
  dirtmp.SetXYZ(MidXp-Clusterx[3],MidYp-Clustery[3],0);
  dirtmp=dirtmp.Unit();
  double x_c=Clusterx[3]+rho*dirtmp.x();
  double y_c=Clustery[3]+rho*dirtmp.y();
  
  double d_c=sqrt(x_c*x_c+y_c*y_c);
#if 0
  double d_ctmp=sqrt( (x_c-Clusterx[3])*(x_c-Clusterx[3])+(y_c-Clustery[3])*(y_c-Clustery[3]));
  std::cout<<"Clusterx[3],y xp,yp, x_c,y_c, d_c "<<Clusterx[3]<<" "
	   <<Clustery[3]<<" "<<MidXp<<" "<<MidYp<<" "<<x_c<<" "
	   <<y_c<<" "<<d_c<<" "<<d_ctmp<<" "<<rho<<" "<<std::endl;
#endif
  double cos_c=-x_c/d_c;
  double sin_c=-y_c/d_c;
	
  //  double x_o=x_c+rho*cos_c;
  //  double y_o=y_c+rho*sin_c;
  //  param[0]=sqrt(x_o*x_o+y_o*y_o );
  param[0]=d_c-rho;
  /*
  if(rho<d_c) param[0]=param[0];
  else if(rho>d_c) param[0]=-param[0];
  */
  MathTools *tool=new MathTools();
  param[1]=tool->CalcDeg(cos_c,sin_c);
  delete tool;
  param[1]=param[1]*TMath::DegToRad();
  
  if(jitter<0 ) 
    {param[2]=-1./rho;param[0]=-1*param[0];param[1]=param[1];}
  else 
    {param[2]=1./rho;param[0]=param[0];param[1]-=TMath::Pi();}
  CirCenterX=(param[0]+1./param[2])*cos(param[1]);
  CirCenterY=(param[0]+1./param[2])*sin(param[1]);
  CirRho=fabs(1./param[2]);
  param[3]=0;
  param[4]=0;
  for(int n=0;n<5;n++) paramtmp[0][n]=param[n];

  return true;


}

double CDSTrack::CalcHelixPhi(const double &x,const double &y)
{
  
  double sin_c=0,cos_c=0;
  cos_c=x-CirCenterX;  sin_c=y-CirCenterY;
  double phi; 
  MathTools *tool=new MathTools();
  phi=tool->CalcDeg(cos_c,sin_c)*TMath::DegToRad();
  delete tool;
  if(phi>TMath::Pi() ) phi-=2*TMath::Pi();

  if(param[2]<0 ) 
    {phi-=param[1];}
  else 
    {phi-=(param[1]+TMath::Pi() );}

     
  while(fabs(phi)>TMath::Pi())
    {
      if(phi>0) phi-=2*TMath::Pi();
      else phi+=2*TMath::Pi();
    }

  return phi;
}


double CDSTrack::CalcHelixPhi(const double &x,const double &y,const double *par)
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


void CDSTrack::Refit(ConfMan *conf)
{
  //Re-Calc CDCHit
  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	{
#if DEBUG
	  std::cout<<"hit "
		   <<"lay "<<trackContainer[layer-1][n].layer()
		   <<" wire "<<trackContainer[layer-1][n].wire()
		   <<" dt "<<trackContainer[layer-1][n].dt()
		   <<" dl "<<trackContainer[layer-1][n].dl()
		   <<" hx "<<trackContainer[layer-1][n].x()
		   <<" hy "<<trackContainer[layer-1][n].y()<<std::endl;
#endif
	  trackContainer[layer-1][n].Calc(conf);
#if DEBUG
	  std::cout<<"af "
		   <<"lay "<<trackContainer[layer-1][n].layer()
		   <<" wire "<<trackContainer[layer-1][n].wire()
		   <<" dt "<<trackContainer[layer-1][n].dt()
		   <<" dl "<<trackContainer[layer-1][n].dl()
		   <<" hx "<<trackContainer[layer-1][n].x()
		   <<" hy "<<trackContainer[layer-1][n].y()<<std::endl;
#endif
	}
    }

  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  magneticfield=CDSParam->GetMagneticField();
  int selected_chi=CDSParam->GetMaxChi();

  if(magneticfield==0)
    {
      
      if(!FirstLineFitting() ) 
	{
#if DEBUG
	  std::cout<<"first line track!!"<<std::endl;
#endif
	  SetGoodFlag(false);  return;
	}
      
#if DEBUG
      std::cout<<" firstline "<<Chi();
#endif
      
      bool flag=true;
      for(int trial=0;trial<10;trial++) { if(!LineFitting() ) flag=false;}
      if(!flag)  
      	{
#if DEBUG
	  std::cout<<"line track!!"<<std::endl;
#endif
	  SetGoodFlag(false); return;
	}     
#if DEBUG
      std::cout<<" line "<<Chi();
#endif


      if(!FirstStereoLineFitting() ) 
	{
#if DEBUG
	  std::cout<<"first sline track!!"<<std::endl; 
#endif
	  SetGoodFlag(false);  return;
	}
      
#if DEBUG
      std::cout<<" firstsl "<<Chi();
#endif

    SetHitPos();       
    int test=0;
    flag=true;
    for(int trial=0;trial<5;trial++) { if(!StereoLineFitting() ) {flag=false;test=trial;break;}}
    if(!flag) 
      {
#if DEBUG
	std::cout<<"sl track ! "<<test<<std::endl; 
#endif
	SetGoodFlag(false); 
	return;
      }
    
    int numt=0;
    double Chi=999,Chitmp=999;
    bool hflag1=true;	
    while(  numt<20 )
      {
	Chi=Chitmp;
	for(int trial=0;trial<10;trial++)
	  { if(!StereoLineFitting() ){hflag1=false; break;}  }	
	Chitmp=CDSTrack::Chi();
	if( (Chi-Chitmp)<0.1 ) break;
	numt++;
      }
    if(!hflag1) { return;}

#if DEBUG    
    std::cout<<" stline "<<Chi()<<std::endl;
#endif


    }
  else
    {

      if(!CalcInitialParameters())
	{
#if DEBUG
	  std::cout<<"ini par!"<<std::endl;
#endif
	  SetGoodFlag(false);return;
	}
      
      if(!FirstCircleFitting() ) 
	{
#if DEBUG
	  std::cout<<"first cir track!!"<<std::endl;
#endif
	  SetGoodFlag(false);  return;
	}
      
#if DEBUG
      std::cout<<" firstcir "<<track->Chi();
#endif
      
      bool flag=true;
      for(int trial=0;trial<10;trial++) { if(!CircleFitting() ) flag=false;}
      if(!flag)  
      	{
#if DEBUG
	  std::cout<<"cir track!!"<<std::endl;
#endif
	  SetGoodFlag(false); return;
	}     
#if DEBUG
      std::cout<<" cir "<<Chi();
#endif


      if(!FirstHelixFitting() ) 
	{
#if DEBUG
	  std::cout<<"first hel track!!"<<std::endl; 
#endif
	  SetGoodFlag(false);  return;
	}
      
#if DEBUG
      std::cout<<" firsthel "<<track->Chi();
#endif
	if( Chi()>2*selected_chi ) 
	  {
#if DEBUG
	    std::cout<<"firsthel chi!"<<std::endl; 
#endif 
	    SetGoodFlag(false);  return;
	  }
		   
    SetHitPos();       
    int test=0;
    flag=true;
    for(int trial=0;trial<5;trial++) { if(!HelixFitting() ) {flag=false;test=trial;break;}}
    if(!flag) 
      {
#if DEBUG
	std::cout<<"hel track ! "<<test<<std::endl;
#endif
	SetGoodFlag(false); 
	return;
      }
    
    int numt=0;
    double Chi=999,Chitmp=999;
    bool hflag1=true;	
    while(  numt<20 )
      {
	Chi=Chitmp;
	for(int trial=0;trial<10;trial++)
	  { if(!HelixFitting() ){hflag1=false; break;}  }	
	Chitmp=CDSTrack::Chi();
	if( (Chi-Chitmp)<0.1 ) break;
	numt++;
      }
    if(!hflag1) { return;}

#if DEBUG    
    std::cout<<" helix "<<Chi()<<std::endl;
#endif

    }
  Calc(conf);
  return;
}



void CDSTrack::Calc(ConfMan *conf)
{
  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  magneticfield=CDSParam->GetMagneticField();
  if(magneticfield==0)
    {
      CirCenterX=0; CirCenterY=0; CirRho=0;      
      pt=0;
      mom=0;
      param[2]=0;
    }
  else
    {
      CirCenterX=(param[0]+1./param[2])*cos(param[1]);
      CirCenterY=(param[0]+1./param[2])*sin(param[1]);
      CirRho=fabs(1./param[2]);
      pt=magneticfield*Const/100./param[2];
      mom=pt*sqrt(1+param[4]*param[4]);
    }

  return;
}


bool CDSTrack::DoLineFit( const int &n, const double *x, const double *y,
                              const double *w,
                              double &a, double &b, double &c, double &chi )
{
  double Sxw=0.,Syw=0.,Sww=0.,Sxy=0.,Sxx=0.;
  double w2;
  for( int i=0; i<n; i++ )
    {
      w2=w[i]*w[i];
      Sxw+=x[i]/w2; Syw+=y[i]/w2; Sww+=1./w2; Sxy+=(x[i]*y[i])/w2; 
      Sxx+=(x[i]*x[i])/w2;
    }
  // ax+by+c=0
  // y = -a/b*x -c/b
  double D = Sww*Sxx-Sxw*Sxw;
  //  std::cout<<"D Swx Swy "<<D<<" "<<Sxw<<" "<<Syw<<std::endl;
  a = -(Sww*Sxy-Sxw*Syw)/D;
  b = 1.;
  c = -(Syw*Sxx-Sxw*Sxy)/D;
  chi = 0.;
  for( int i=0; i<n; i++ ){
    chi += ((-a/b*x[i]-c/b)-y[i])*((-a/b*x[i]-c/b)-y[i])/(w[i]*w[i]);
  }
  chi /= (n-2);
  return true;
}
   
bool CDSTrack::FirstLineFitting()
{
  
  double x[40],y[40],weight[40];
  int numallhit=0;  
  for(int numlayer=0;numlayer<7;numlayer++)
    {
      int layer=0;
      if(numlayer<3) layer=numlayer+1;
      else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
      else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	{
	  double dl=trackContainer[layer-1][n].dl();
	  if(dl<ResolutionOfCDC) dl=ResolutionOfCDC;
	  x[numallhit]=trackContainer[layer-1][n].wx();
	  y[numallhit]=trackContainer[layer-1][n].wy();
	  weight[numallhit]=dl;
	  numallhit++;
	}
    }
  double  meanphi=theta*TMath::DegToRad();  
  for( int i=0; i<numallhit; i++ ){
    double tmpx = x[i], tmpy = y[i];
    x[i] = tmpx*cos(meanphi) + tmpy*sin(meanphi);
    y[i] = -tmpx*sin(meanphi) + tmpy*cos(meanphi);
  }
  double a,b,c;//,chi;
  DoLineFit(numallhit,x,y,weight,a,b,c,ChiSqr);     
  double a2 = a, b2 = b;
  a = a2*cos(meanphi)-b2*sin(meanphi);
  b = a2*sin(meanphi)+b2*cos(meanphi);
  Al = a; Bl = b; Cl = c;

  dof=numallhit-2; 
  if( ChiSqr >99999 ) 
    {
#if DEBUG
	std::cout<<"Too big chi2/dof fline"<<std::endl;
#endif
	return false;
      }
  fittinglevel=1;

  //### Set some param#######
  double x0,y0,tx0=0,ty0=-c/b2;
  x0=tx0*cos(meanphi)-ty0*sin(meanphi);
  y0=tx0*sin(meanphi)+ty0*cos(meanphi);
  TVector3 tmp_o,zero;
  double d_x,td_x=30;  double d_y,td_y=(-a2/b2*30);
  d_x=td_x*cos(meanphi)-td_y*sin(meanphi);
  d_y=td_x*sin(meanphi)+td_y*cos(meanphi);
  Line_d.SetXYZ(d_x,d_y,0);Line_d.Unit();
  tmp_o.SetXYZ(x0,y0,0);  zero.SetXYZ(0,0,0);
  double dis;
  PointToLine(zero,tmp_o,Line_d,dis,Line_o);
  param[2]=0;
  param[0]=dis;
  MathTools *tool=new MathTools();
  double tmpdeg=tool->CalcDeg(Line_o.x(),Line_o.y());
  param[1]=tmpdeg*TMath::DegToRad();
  tmpdeg+=90;
  if(tmpdeg>360) tmpdeg-=360;
  tmpdeg-=theta;tmpdeg=fabs(tmpdeg);if(tmpdeg>180) tmpdeg-=360;
  if(fabs(tmpdeg)>90 ) {param[0]*=-1;param[1]+=TMath::Pi();}
  if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi(); 
  delete tool;
  SetHitPos();
  return true;  
}

   
bool CDSTrack::LineFitting()
{
  
  double x[40],y[40],weight[40];
  int numallhit=0;  
  for(int numlayer=0;numlayer<7;numlayer++)
    {
      int layer=0;
      if(numlayer<3) layer=numlayer+1;
      else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
      else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	{
	  double dl=trackContainer[layer-1][n].dl();
	  if(dl<0) dl=0.0001;
	  x[numallhit]=trackContainer[layer-1][n].x();
	  y[numallhit]=trackContainer[layer-1][n].y();
	  weight[numallhit]=ResolutionOfCDC;
	  numallhit++;
	}
    }
  double meanphi=theta*TMath::DegToRad();  
  for( int i=0; i<numallhit; i++ ){
    double tmpx = x[i], tmpy = y[i];
    x[i] = tmpx*cos(meanphi) + tmpy*sin(meanphi);
    y[i] = -tmpx*sin(meanphi) + tmpy*cos(meanphi);
  }
  double a,b,c;
  DoLineFit(numallhit,x,y,weight,a,b,c,ChiSqr);     
  double a2 = a, b2 = b;
  a = a2*cos(meanphi)-b2*sin(meanphi);
  b = a2*sin(meanphi)+b2*cos(meanphi);
  Al = a; Bl = b; Cl = c;
#if DEBUG
  std::cout<<"linechi "<<ChiSqr<<std::endl;
#endif
  dof=numallhit-2; 
  if( ChiSqr >99999 ) 
    {
#if DEBUG
	std::cout<<"Too big chi2/dof line"<<std::endl;
#endif
	return false;
      }
  fittinglevel=2;

  //### Set some param#######
  double x0,y0,tx0=0,ty0=-c/b2;
  x0=tx0*cos(meanphi)-ty0*sin(meanphi);
  y0=tx0*sin(meanphi)+ty0*cos(meanphi);
  TVector3 tmp_o,zero;
  double d_x,td_x=30;  double d_y,td_y=(-a2/b2*30);
  d_x=td_x*cos(meanphi)-td_y*sin(meanphi);
  d_y=td_x*sin(meanphi)+td_y*cos(meanphi);
  Line_d.SetXYZ(d_x,d_y,0);Line_d.Unit();
  tmp_o.SetXYZ(x0,y0,0);  zero.SetXYZ(0,0,0);
  double dis;
  PointToLine(zero,tmp_o,Line_d,dis,Line_o);
  param[2]=0;
  param[0]=dis;
  MathTools *tool=new MathTools();
  double tmpdeg=tool->CalcDeg(Line_o.x(),Line_o.y());
  param[1]=tmpdeg*TMath::DegToRad();
  tmpdeg+=90;
  if(tmpdeg>360) tmpdeg-=360;
  tmpdeg-=theta;tmpdeg=fabs(tmpdeg);if(tmpdeg>180) tmpdeg-=360;
  if(fabs(tmpdeg)>90 ) {param[0]*=-1;param[1]+=TMath::Pi();}
  if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi(); 
  delete tool;
  SetHitPos();

  return true;  
}


bool CDSTrack::TestLineFitting()
{
  
  double x[40],y[40],weight[40];
  int numallhit=0;  
  for(int numlayer=0;numlayer<7;numlayer++)
    {
      int layer=0;
      if(numlayer<3) layer=numlayer+1;
      else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
      else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	{
	  double dl=trackContainer[layer-1][n].dl();
	  if(dl<0) dl=0.0001;
	  x[numallhit]=trackContainer[layer-1][n].x();
	  y[numallhit]=trackContainer[layer-1][n].y();
	  weight[numallhit]=ResolutionOfCDC;
	  numallhit++;
	}
    }
  double meanphi=2.00*TMath::DegToRad();  //theta*TMath::DegToRad();  
  for( int i=0; i<numallhit; i++ ){
    double tmpx = x[i], tmpy = y[i];
    x[i] = tmpx*cos(meanphi) + tmpy*sin(meanphi);
    y[i] = -tmpx*sin(meanphi) + tmpy*cos(meanphi);
  }
  double a,b,c;
  DoLineFit(numallhit,x,y,weight,a,b,c,ChiSqr);     
  double a2 = a, b2 = b;
  a = a2*cos(meanphi)-b2*sin(meanphi);
  b = a2*sin(meanphi)+b2*cos(meanphi);
  Al = a; Bl = b; Cl = c;
#if DEBUG
  std::cout<<"linechi "<<ChiSqr<<std::endl;
#endif
  dof=numallhit-2; 
  if( ChiSqr >99999 ) 
    {
#if DEBUG
	std::cout<<"Too big chi2/dof line"<<std::endl;
#endif
	return false;
      }
  fittinglevel=2;

  //### Set some param#######
  double x0,y0,tx0=0,ty0=-c/b2;
  x0=tx0*cos(meanphi)-ty0*sin(meanphi);
  y0=tx0*sin(meanphi)+ty0*cos(meanphi);
  TVector3 tmp_o,zero;
  double d_x,td_x=30;  double d_y,td_y=(-a2/b2*30);
  d_x=td_x*cos(meanphi)-td_y*sin(meanphi);
  d_y=td_x*sin(meanphi)+td_y*cos(meanphi);
  Line_d.SetXYZ(d_x,d_y,0);Line_d.Unit();
  tmp_o.SetXYZ(x0,y0,0);  zero.SetXYZ(0,0,0);
  double dis;
  PointToLine(zero,tmp_o,Line_d,dis,Line_o);
  param[2]=0;
  param[0]=dis;
  MathTools *tool=new MathTools();
  double tmpdeg=tool->CalcDeg(Line_o.x(),Line_o.y());
  param[1]=tmpdeg*TMath::DegToRad();
  tmpdeg+=90;
  if(tmpdeg>360) tmpdeg-=360;
  tmpdeg-=theta;tmpdeg=fabs(tmpdeg);if(tmpdeg>180) tmpdeg-=360;
  if(fabs(tmpdeg)>90 ) {param[0]*=-1;param[1]+=TMath::Pi();}
  if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi(); 
  delete tool;
  //  SetHitPos();

  return true;  
}


bool CDSTrack::FirstStereoLineFitting()
{

  double disZ[40];
  int numstlayer=0;

  for(int stlayer=0 ;stlayer<8;stlayer++)
    {
      int layer=0;
      if(0<=stlayer && stlayer<=3) layer=stlayer+4;
      else if(4<=stlayer && stlayer<=7) layer=stlayer+6;

      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	{
	  CDCHit *cdc=&trackContainer[layer-1][n];
	  TVector3 wpos,wposp;
	  
	  wpos.SetXYZ(trackContainer[layer-1][n].wx(),
		      trackContainer[layer-1][n].wy(),
		      trackContainer[layer-1][n].wz() );
	  
	  wposp.SetXYZ(trackContainer[layer-1][n].wxp(),
		       trackContainer[layer-1][n].wyp(),
		       trackContainer[layer-1][n].wzp() );
	  double aw,bw,cw;  
	  if((wposp.x()-wpos.x() )!=0)
	    {
	      aw=(wposp.y()-wpos.y() )/(wposp.x()-wpos.x() );
	      bw=-1;
	      cw= -(aw*wpos.x()+bw*wpos.y() );
	    }
	    else
	      {
		bw=(wposp.x()-wpos.x() )/(wposp.y()-wpos.y() );
		aw=-1;
		cw= -(aw*wpos.x()+bw*wpos.y() );
	      }
	  if((Al*bw-aw*Bl)==0) continue;
	  double x=(Bl*cw-bw*Cl)/(Al*bw-aw*Bl);
	  double y=-(Al*cw-aw*Cl)/(Al*bw-aw*Bl);
	  TVector3 WPos;WPos.SetXYZ(x,y,0);

	  double dis=sqrt( ( wpos.x()-WPos.x() )*(wpos.x()-WPos.x())
			   +(wpos.y()-WPos.y())*(wpos.y()-WPos.y() ) );
	  double wdis=sqrt( (wpos.x()-wposp.x() )*(wpos.x()-wposp.x() )+(wpos.y()-wposp.y() )*(wpos.y()-wposp.y() ) );

	  WPos.SetZ( wpos.z()+(wposp.z()-wpos.z() )*dis/wdis );
	  trackContainer[layer-1][n].SetHitPosition(WPos.x(),WPos.y(),WPos.z());    

	  //########dis Z ###########//	  
	  double dl=cdc->dl();
	  double tilt=cdc->tilt()*TMath::DegToRad();
	  disZ[numstlayer]=fabs(dl/sin(tilt) );
	  //	  std::cout<<"z[ "<<numstlayer<<"] "<<WPos.z()<<" "
	  //   <<disZ[numstlayer]<<" "<<dis<<" "<<wdis<<std::endl;      
	  numstlayer++;

	}

    }
  
  int checki=-1;	  
  for(int i =0;i<pow(2,numstlayer);i++)
    {
      double hitx[64],hity[64],hitz[64];
      double hitr[64];
      double weight[64];

      int allhit=0;
      int numst=0;
      double meanphi=0;      
      for(int layer=1;layer<=15;layer++)
	{
	  for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	    {
	      CDCHit *cdc=&trackContainer[layer-1][n];
	      
	      if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
		{
		
		  hitx[allhit]= cdc->x();
		  hity[allhit]= cdc->y();
		  hitz[allhit]= cdc->z();
		  weight[allhit]=ResolutionOfCDC*10;
		  meanphi=theta*TMath::DegToRad();
		  hitr[allhit]=hitx[allhit]*cos(meanphi)+hity[allhit]*sin(meanphi);
		  int check=i;		  
		  if(check<pow(2,numst) ) hitz[allhit]+=disZ[numst];
		  else 
		    {
		      int cflag;
		      for(int n=0;n<numst;n++)
			{
			  cflag=check%2;
			  check=(check-cflag)/2;
			}
		      cflag=check%2;
		      if(cflag==0)  hitz[allhit]+=disZ[numst];
		      else if(cflag==1) hitz[allhit]-=disZ[numst];
		    }
		  numst++;
		  allhit++;
		  
		}//if stereo	      
	    }//cdchit
	}//layer
      
      double aa,bb,cc;      double chitmp=0;
      DoLineFit(allhit,hitr,hitz,weight,aa,bb,cc,chitmp);      
      //std::cout<<"allhit hitr hitz w "<<allhit<<" "<<hitr[0]<<" "<<hitz[0]<<" "<<weight[0]<<std::endl;
      if(i==0)
	{
	  //Ap*x+Bp*y+Cp*z+Dp=0;
	  Ap=aa*cos(meanphi); Bp=aa*sin(meanphi);
	  Cp=bb; Dp=cc;
	  ChiSqr= chitmp;
	  checki=i;	  
	}	
      else if(chitmp<ChiSqr)
	{
	  Ap=aa*cos(meanphi); Bp=aa*sin(meanphi);
	  Cp=bb; Dp=cc;
	  ChiSqr= chitmp;
	  checki=i;	  
	}
      
    }
  
  int numst=0;
  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	{
	  if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
	    {
	      double hitx= trackContainer[layer-1][n].x();
	      double hity= trackContainer[layer-1][n].y();
	      double hitz= trackContainer[layer-1][n].z();
	      int check=checki;
	      if(check<pow(2,numst) )
		hitz+=disZ[numst];
	      else 
		{
		  int cflag;
		  for(int n=0;n<numst;n++)
		    {
		      cflag=check%2;
			  check=(check-cflag)/2;
		    }
		  cflag=check%2;
		  if(cflag==0)  hitz+=disZ[numst];
		  else if(cflag==1) hitz-=disZ[numst];
		}
	      trackContainer[layer-1][n].SetHitPosition(hitx,hity,hitz);   
	      numst++;
	    }
	  
	}
    }
  ChiSqrtmp[2]=ChiSqr;

  if( ChiSqr >999999 ) 
    {
#if DEBUG
      std::cout<<"Too big chi2/dof fsl"<<std::endl;
#endif
      return false;
      }
  fittinglevel=3;      
  double z_o=-(Ap*Line_o.x()+Bp*Line_o.y()+Dp)/Cp;
  Line_o.SetZ(z_o);param[3]=z_o;
  TVector3 tmp;tmp=Line_o+30*Line_d;
  double z_d=-(Ap*tmp.x()+Bp*tmp.y()+Dp)/Cp;
  tmp.SetZ(z_d);
  Line_d=(tmp-Line_o).Unit();
  tmp.SetXYZ(Line_d.x(),Line_d.y(),0);
  param[4]=-Line_d.z()/tmp.Mag();
  for(int n=0;n<5;n++) paramtmp[3][n]=param[n];
  SetHitPos();
  return true;
}


   
bool CDSTrack::StereoLineFitting()
{
  double x[40],y[40],z[40],weight[40];
  int numallhit=0;  
  for(int layer=1;layer<=15;layer++)
    {
	      
      if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
	{ 
	  for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	    {
	      double dl=trackContainer[layer-1][n].dl();
	      if(dl<0) dl=0.0001;
	      x[numallhit]=trackContainer[layer-1][n].x();
	      y[numallhit]=trackContainer[layer-1][n].y();
	      z[numallhit]=trackContainer[layer-1][n].z();
	      double tilt=trackContainer[layer-1][n].tilt()*TMath::DegToRad();
	      weight[numallhit]=ResolutionOfCDC*sin(tilt);
	      numallhit++;
	    }
	}
    }

  double meanphi=theta*TMath::DegToRad();  
  for( int i=0; i<numallhit; i++ ){
    double tmpx = x[i], tmpy = y[i];
    x[i] = tmpx*cos(meanphi) + tmpy*sin(meanphi);
    y[i] = -tmpx*sin(meanphi) + tmpy*cos(meanphi);
  }
  double a,b,c,chi;
  DoLineFit(numallhit,x,z,weight,a,b,c,chi);     
  //  double a2 = a, b2 = b;

  //Ap*x+Bp*y+Cp*z+Dp=0;
  Ap=a*cos(meanphi); Bp=a*sin(meanphi);
  Cp=b; Dp=c;

  fittinglevel=4;
  //### Set some param#######
  double z_o=-(Ap*Line_o.x()+Bp*Line_o.y()+Dp)/Cp;
  Line_o.SetZ(z_o);param[3]=z_o;
  TVector3 tmp;tmp=Line_o+30*Line_d;
  double z_d=-(Ap*tmp.x()+Bp*tmp.y()+Dp)/Cp;
  tmp.SetZ(z_d);
  Line_d=(tmp-Line_o).Unit();
  tmp.SetXYZ(Line_d.x(),Line_d.y(),0);
  param[4]=-Line_d.z()/tmp.Mag();
  for(int n=0;n<5;n++) paramtmp[4][n]=param[n];
  CalcLineChiSqr();
  if( ChiSqr >99999 ) 
    {
#if DEBUG
	std::cout<<"Too big chi2/dof sl"<<std::endl;
#endif
	return false;
    }
  SetHitPos();
  return true;  
}

bool CDSTrack::TestStereoLineFitting()
{
  double x[40],y[40],z[40],weight[40];
  int numallhit=0;  
  for(int layer=1;layer<=15;layer++)
    {
	      
      if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
	{ 
	  for(int n=0;n<(int)trackContainer[layer-1].size();n++)
	    {
	      double dl=trackContainer[layer-1][n].dl();
	      if(dl<0) dl=0.0001;
	      x[numallhit]=trackContainer[layer-1][n].x();
	      y[numallhit]=trackContainer[layer-1][n].y();
	      z[numallhit]=trackContainer[layer-1][n].z();
	      double tilt=trackContainer[layer-1][n].tilt()*TMath::DegToRad();
	      weight[numallhit]=ResolutionOfCDC*sin(tilt);
	      numallhit++;
	    }
	}
    }

  double meanphi=2.00*TMath::DegToRad();//theta*TMath::DegToRad();  
  for( int i=0; i<numallhit; i++ ){
    double tmpx = x[i], tmpy = y[i];
    x[i] = tmpx*cos(meanphi) + tmpy*sin(meanphi);
    y[i] = -tmpx*sin(meanphi) + tmpy*cos(meanphi);
  }
  double a,b,c,chi;
  DoLineFit(numallhit,x,z,weight,a,b,c,chi);     
  //  double a2 = a, b2 = b;

  //Ap*x+Bp*y+Cp*z+Dp=0;
  Ap=a*cos(meanphi); Bp=a*sin(meanphi);
  Cp=b; Dp=c;

  fittinglevel=4;
  //### Set some param#######
  double z_o=-(Ap*Line_o.x()+Bp*Line_o.y()+Dp)/Cp;
  Line_o.SetZ(z_o);param[3]=z_o;
  TVector3 tmp;tmp=Line_o+30*Line_d;
  double z_d=-(Ap*tmp.x()+Bp*tmp.y()+Dp)/Cp;
  tmp.SetZ(z_d);
  Line_d=(tmp-Line_o).Unit();
  tmp.SetXYZ(Line_d.x(),Line_d.y(),0);
  param[4]=-Line_d.z()/tmp.Mag();
  for(int n=0;n<5;n++) paramtmp[4][n]=param[n];
  CalcLineChiSqr();
  if( ChiSqr >99999 ) 
    {
#if DEBUG
	std::cout<<"Too big chi2/dof sl"<<std::endl;
#endif
	return false;
    }
  SetHitPos();
  return true;  
}

void CDSTrack::CalcLineChiSqr()
{
  TVector3 P( 0., -Cl/Bl, (Bp*Cl/Bl-Dp)/Cp );
  TVector3 Q( Bl, -Al, Bl*(Bp*Al/Bl-Ap)/Cp );
  TVector3 WirePosition,WirePositionp, WireDirection;
  TVector3 Xest, Next;
  CDCHit *hit = 0;
  double dist;
  double chi = 0;
  int nallhit=0;
  for(int layer=1;layer<=15;layer++)
    {
      for( int i=0; i<(int)trackContainer[layer-1].size(); i++ ){
	hit = &trackContainer[layer-1][i];
 	double dl = fabs(hit->dl());
	WirePosition.SetXYZ( hit->wx(), hit->wy(), hit->wz() );
	WirePositionp.SetXYZ( hit->wxp(), hit->wyp(), hit->wzp() );
	WireDirection = (WirePositionp-WirePosition).Unit();
	LineToLine( P, Q, WirePosition, WireDirection, dl, dist, Xest, Next );
	//LineToLine( Line_o, Line_d, WirePosition, WireDirection, dl, dist, Xest, Next );
	chi += (dl-dist)*(dl-dist)/(ResolutionOfCDC*ResolutionOfCDC);
	//	std::cout<<"dl dis "<<dl<<" "<< dist<<std::endl;
	
	nallhit++;
      }
    }
  chi /= (double)(nallhit-4);
  ChiSqr=chi;
  dof=nallhit-4;
}


void CDSTrack::PointToLine( const TVector3 &p,
			   const TVector3 &x, const TVector3 &a,
			   double &dist,TVector3 &xest )
{
  double k=(p*a-x*a)/a.Mag2();
  xest=x+k*a;
  dist=(xest-p).Mag();
}

void CDSTrack::LineToLine( const TVector3 &x1, const TVector3 &a1,
			   const TVector3 &x2, const TVector3 &a2,
			   const double &dl,
			   double &dist,
			   TVector3 &xest, TVector3 &next )
{
  TVector3 x = x2-x1;   // x = x1 + t*a1
  double a =  a1.Dot(a1);  // x = x2 + s*a2
  double b = -a1.Dot(a2);  //    ||
  double c =  a2.Dot(a1);  //    \/
  double d = -a2.Dot(a2);  // a*t + b*s = A1
  double A1 = a1.Dot(x);   // c*t + d*s = A2
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


