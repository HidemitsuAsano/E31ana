#include "CDSTrackingMan.h"
#define DEBUG 0
ClassImp(CDSTrackingMan);
ClassImp(TrackVertex);

const double Inclradius=20.40,Outclradius=47.78,cellsize=0.845*2.;
//const double ResolutionOfCDC=0.02;


CDSTrackingMan::CDSTrackingMan(): TObject()
{
  Clear();
  KilledLayer=0;
  selected_chi=50;
}

void CDSTrackingMan::Clear()
{
  trackContainer.clear();
  Vertex.clear();
  Vertex_beam.clear();
  GoodTrack.clear();
  for(int i=0;i<7;i++) ClusterContainer[i].clear();
}

void CDSTrackingMan::AddTrack(const CDSTrack &atrack )
{}

bool CDSTrackingMan::SearchCluster(CDSHitMan *cdsMan,ConfMan *conf)
{
  double fac_cell=1.2;

  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  int cdctdc_u=CDSParam->GetUpperLimitCDCTDC();
  int cdctdc_l=CDSParam->GetLowerLimitCDCTDC();

  for(int num_slayer=1;num_slayer<=7;num_slayer++)
    {
      int inilayer=0;
      int size_slayer=0;
      if(num_slayer==1){inilayer=1;size_slayer=3;}
      else if(2<=num_slayer && num_slayer<=7){inilayer=num_slayer*2;size_slayer=2;}
      int layer2=inilayer+1; 
      int layer3=inilayer+2;
      int layer2hit[256]={0};

      for(int i=0;i<cdsMan->nCDC(inilayer); i++)
	{
	  if(inilayer==KilledLayer) break;
	  CDCHit *cdc1=cdsMan->CDC(inilayer,i);
	  int tdc = cdc1->tdc();
	  if( !(cdctdc_l<tdc && tdc<cdctdc_u) ) continue;	  
	  double x1=cdc1->wx();    double y1=cdc1->wy();
	  bool hit_layer2flag=false;
	  for(int ii=0;ii<cdsMan->nCDC(layer2);ii++)
	    {
	      if(layer2==KilledLayer) break;
	      CDCHit *cdc2=cdsMan->CDC(layer2,ii);
	      int wire2=cdc2->wire();
	      int tdc2 = cdc2->tdc();
	      if( !(cdctdc_l<tdc2 && tdc2<cdctdc_u) ) continue;	  	      
	      double x2=cdc2->wx();  double y2=cdc2->wy();
	      double dist=sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2) );
	      if( dist>1e-5 && dist < cellsize*fac_cell )
		{
		  layer2hit[wire2]++;
		  hit_layer2flag=true;
		  if(num_slayer!=1)
		    {
		      HitCluster TmpCluster;
		      TmpCluster.Clear();
		      TmpCluster.SetHit( *cdc1 ); TmpCluster.SetHit( *cdc2 );
		      ClusterContainer[num_slayer-1].push_back(TmpCluster);
		    }
		  if(num_slayer==1)		  
		    {

		      bool hit_layer3flag=false;
		      for(int iii=0;iii<cdsMan->nCDC(layer3);iii++)
			{
			  if(layer3==KilledLayer) break;     		 
			  CDCHit *cdc3=cdsMan->CDC(layer3,iii);
			  int tdc3 = cdc3->tdc();
			  if( !(cdctdc_l<tdc3 && tdc3<cdctdc_u) ) continue;  
			  double x3=cdc3->wx(); double y3=cdc3->wy();
			  double dist2=sqrt( (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3) );
			  if( dist2>1e-5 && dist2 < cellsize*fac_cell )
			    {
			      HitCluster TmpCluster;
			      TmpCluster.Clear();
			      TmpCluster.SetHit( *cdc1 );
			      TmpCluster.SetHit( *cdc2 );
			      TmpCluster.SetHit( *cdc3 );
			      ClusterContainer[num_slayer-1].push_back(TmpCluster);
			      hit_layer3flag=true;
			    }
			}//layer3
		      
		      if(!hit_layer3flag )
			{
			  HitCluster TmpCluster;
			  TmpCluster.Clear();
			  TmpCluster.SetHit( *cdc1 );
			  TmpCluster.SetHit( *cdc2 );
			  ClusterContainer[num_slayer-1].push_back(TmpCluster);
			}
		      
		    }//slayer==1
		}//near hit
	    }//layer2
	      
	  if(!hit_layer2flag)
	    {
	      if(num_slayer!=1)
		{
		  HitCluster TmpCluster;
		  TmpCluster.SetHit( *cdc1 );
		  ClusterContainer[num_slayer-1].push_back(TmpCluster);
		}
	      if(num_slayer==1)
		{
		  for(int iii=0;iii<cdsMan->nCDC(layer3);iii++)
		    {
		      if(layer3==KilledLayer) break;
		      CDCHit *cdc3=cdsMan->CDC(layer3,iii);
		      int tdc3 = cdc3->tdc();
		      if( !(cdctdc_l<tdc3 && tdc3<cdctdc_u) ) continue;	  

		      double x3=cdc3->wx(); double y3=cdc3->wy();
		      double dist3=sqrt( (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3) );
		      if( dist3>1e-5 && dist3 < 2*cellsize*fac_cell )
			{
			  HitCluster TmpCluster;
			  TmpCluster.Clear();
			  TmpCluster.SetHit( *cdc1 ); TmpCluster.SetHit( *cdc3 );
			  ClusterContainer[num_slayer-1].push_back(TmpCluster);
			}
		    }//layer3
		}//slayer==1
	    }//nohit_layer2
	}//layer1

      for(int ii=0;ii<cdsMan->nCDC(layer2);ii++)
	{
	  if(layer2==KilledLayer) break;
      	  CDCHit *cdc2=cdsMan->CDC(layer2,ii);
	  int tdc2 = cdc2->tdc();
	  if( !(cdctdc_l<tdc2 && tdc2<cdctdc_u) ) continue;	  

	  int wire=cdc2->wire();
	  double x2=cdc2->wx(); double y2=cdc2->wy();
	  bool hit_flag=false;
	  if(layer2hit[wire]>0) hit_flag=true;
	  if(!hit_flag)
	    {
	      if(num_slayer!=1)
		{
		  HitCluster TmpCluster;
		  TmpCluster.Clear();
		  TmpCluster.SetHit( *cdc2 );
		  ClusterContainer[num_slayer-1].push_back(TmpCluster);
		}
	      if(num_slayer==1)	      
		{
		  for(int iii=0;iii<cdsMan->nCDC(layer3);iii++)
		    {
		      if(layer3==KilledLayer) break;
		
		      CDCHit *cdc3=cdsMan->CDC(layer3,iii);
		      int tdc3 = cdc3->tdc();
		      if( !(cdctdc_l<tdc3 && tdc3<cdctdc_u) ) continue;	  
		      double x3=cdc3->wx(); double y3=cdc3->wy();
		      double dist=sqrt( (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3) );
		      if( dist>1e-5 && dist < cellsize*fac_cell )
			{
			  HitCluster TmpCluster;
			  TmpCluster.Clear();
			  TmpCluster.SetHit(*cdc2);TmpCluster.SetHit( *cdc3 );
			  ClusterContainer[num_slayer-1].push_back(TmpCluster);
			}
		    }//layer3
		}//slayer==1
	    }
	}//layer2

      for(int n=0;n<(int)ClusterContainer[num_slayer-1].size();n++)
	{
	  ClusterContainer[num_slayer-1][n].Calc();
	}

    }//slayer

  return true;
}
    
bool CDSTrackingMan::FindingTrackCandidate(CDSHitMan *cdsMan,ConfMan *conf)
{    
   
    CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
    double magneticfield=CDSParam->GetMagneticField();

    double JITTER;
    if(magneticfield==0) JITTER=3.0;
    else JITTER=15.0;

    int numincl=ClusterContainer[0].size();
    for(int incl=0;incl<numincl;incl++)
      {

	//	bool trackfin=false;
	int outclnum=0;
	CDSTrack TmpTrack;
	for(int outcl=0;outcl<(int)ClusterContainer[6].size(); outcl++)
	  {
	    TmpTrack.Clear();
	    double def_deg=fabs(ClusterContainer[0][incl].Meanphi()
				-ClusterContainer[6][outcl].Meanphi());
	    if(def_deg>180) def_deg=360-def_deg;
	    if(def_deg>90 ) continue;
	    
	    outclnum++;
	    double x_in=ClusterContainer[0][incl].Meanx();
	    double y_in=ClusterContainer[0][incl].Meany();
	    double x_out=ClusterContainer[6][outcl].Meanx();
	    double y_out=ClusterContainer[6][outcl].Meany();
	
	    double a=1,b,c; //a*x+b*y+c=0
	    b=-(x_in-x_out)/(y_in-y_out);
	    c=-(x_in+b*y_in);
	    TmpTrack.SetABC(a,b,c);
	    TmpTrack.SetClusterX(1,x_in);TmpTrack.SetClusterY(1,y_in);	
	    TmpTrack.SetClusterY(7,x_out);TmpTrack.SetClusterY(7,y_out);

//#############mid axial cluster################
	    int num_midcl=0;
	    if((int)ClusterContainer[3].size()==0) continue;
	    for(int midcl=0; midcl<(int)ClusterContainer[3].size(); midcl++)
	      {		
	
		double phi_mid=ClusterContainer[3][midcl].Meanphi();
		double x_mid=ClusterContainer[3][midcl].Meanx();
		double y_mid=ClusterContainer[3][midcl].Meany();
		double dis_mid =fabs(a*x_mid+b*y_mid+c)/sqrt(a*a+b*b);
		double def_deg2=fabs(ClusterContainer[6][outcl].Meanphi()
				     -phi_mid);
		if(def_deg2>180) def_deg2=360-def_deg2;
		if(dis_mid<JITTER && def_deg2<90)
		  {
		    MathTools *tool=new MathTools();
		    TmpTrack.SetClusterX(4,x_mid);TmpTrack.SetClusterY(4,y_mid);	
		    double theta_o=0,theta_m=0;
		    theta_o=tool->CalcDeg(x_out-x_in,y_out-y_in);
		    theta_m=tool->CalcDeg(x_mid-x_in,y_mid-y_in);
		    delete tool;
		    if(fabs(theta_o-theta_m)>330) {theta_o-=360;theta_m-=360;}
		    if(fabs(theta_o-theta_m)<0) dis_mid=dis_mid;
		    else dis_mid=-1*dis_mid;
		    TmpTrack.SetJitter(dis_mid);
		    double theta=tool->CalcDeg(x_mid,y_mid);
		    TmpTrack.SetTheta(theta);

		    for(int i=0;i< ClusterContainer[0][incl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[0][incl].CDC(i) ) );

		      }
		    for(int i=0;i<(int)ClusterContainer[3][midcl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[3][midcl].CDC(i) ) );
		      }
		    for(int i=0;i< ClusterContainer[6][outcl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[6][outcl].CDC(i) ) );
		      }

		    num_midcl++;		    
		    trackContainer.push_back(TmpTrack);
		    TmpTrack.Clear();       		       
		  }

	      }//midcl 
	  }//outcl
      }//incl
	
    return true;
}

void  CDSTrackingMan::SelectSharedHit()
{
  
  for(int itr=0;itr<(int)trackContainer.size()-1;itr++)
    {
      CDSTrack *track1=&trackContainer[itr];
      double chi1=track1->Chi();
      if( !(track1->FittingLevel()==2 || track1->FittingLevel()==4) ) continue;
      for(int layer=1;layer<=NumOfCDCLayers;layer++)
	{ for(int ih1=0;ih1<track1->nTrackHit(layer);ih1++)
	  {
	    bool flaghit1=true;
	    CDCHit *hit1=track1->TrackHit(layer,ih1);
	    int wire1=hit1->wire();
	    for(int itr2=itr+1;itr2<(int)trackContainer.size();itr2++)
	      { 
		if(!flaghit1) break;
		CDSTrack *track2=&trackContainer[itr2];
		double chi2=track2->Chi();		
		if( !(track2->FittingLevel()==2 || track2->FittingLevel()==4) ) continue;

		for(int ih2=0;ih2<track2->nTrackHit(layer);ih2++)
		  {
		    CDCHit *hit2=track2->TrackHit(layer,ih2);
		    int wire2=hit2->wire();
		    if(wire1==wire2)
		      {
			if(chi1<chi2) track2->DeleteHit(layer,ih2);
			else flaghit1=false;
		      }
		  }
	      }//track2
	    if(!flaghit1) track1->DeleteHit(layer,ih1);
	  }//Hit1
	}//layer
   }//track1

  //###Delete track with too less Hit
  TrackContainer::iterator ittr;
  for(ittr=trackContainer.end()-1;ittr!=trackContainer.begin()-1;--ittr)
    {
      int hitlayer=0;
      CDSTrack *track=&*ittr;
      for(int layer=1;layer<=NumOfCDCLayers;layer++)
	{ 
	  if(track->nTrackHit(layer)>0) hitlayer++;
	}
      //      std::cout<<"trackLevel "<<track->FittingLevel()<<" chi "<<track->Chi()<<std::endl;
      if( !(track->FittingLevel()==2 || track->FittingLevel()==4) )   trackContainer.erase(ittr);//delete track
      if(track->FittingLevel()==2 &&( hitlayer<= 4 || track->Chi()>3*selected_chi )) trackContainer.erase(ittr);//delete track
      else if(track->FittingLevel()==4 && hitlayer<= 10) trackContainer.erase(ittr);//delete track


    }
}

bool CDSTrackingMan::FindingStereoHit(CDSHitMan *cdsMan,ConfMan *conf)
{

   
  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  double magfield=CDSParam->GetMagneticField();
    
  //########Hit streo wire###########

  TrackContainer tmptrackContainer;
  for(int itr=0;itr<(int)trackContainer.size();itr++)
    {
      std::vector<int> tmpvec[4];
      CDSTrack *track=&trackContainer[itr];
      if(track->FittingLevel()!=2) continue;
      double rho=track->Rho();
      double xcen=track->CenterX();
      double ycen=track->CenterY();
      TVector3 Line_o,Line_d;
      if(magfield==0)
	{
	  Line_o=track->LineOrigin(); Line_d=track->LineDirection();
	}
      for( int streonum=0; streonum<4; streonum++ )
	{
	  int streoslayer=0;
	  if(0<=streonum && streonum<=1) streoslayer=streonum+1; 
	  else if(2<=streonum && streonum<=3) streoslayer=streonum+2; 

	  for( int k=0; k<(int)ClusterContainer[streoslayer].size(); k++ )
	    {
	      // int tdc = cdsMan->CDC(streolayer,k)->tdc();
	      // if( !(400<tdc && tdc<900)) continue;
	      //double radius=cdsMan->CDC(streolayer,k)->radius();
	      HitCluster *cluster=&ClusterContainer[streoslayer][k];
	      double phi=cluster->Meanphi();
	      double def_deg_st=fabs(phi-track->Theta());
	      if(def_deg_st>180 ) def_deg_st=360-def_deg_st;
	      if(def_deg_st>90) continue;
	    
	      double x=cluster->Meanx();
	      double y=cluster->Meany();
	      double xp=cluster->Meanxp();
	      double yp=cluster->Meanyp();
	      
	      double wiredis = sqrt( (x-xp)*(x-xp)+(y-yp)*(y-yp) );
	      double dis,dis2,xest,yest;
	      TVector3 pos,posp,lnest; pos.SetXYZ(x,y,0);posp.SetXYZ(xp,yp,0);
	      if(magfield==0){
		track->PointToLine(pos,Line_o,Line_d,dis,lnest); 
		track->PointToLine(posp,Line_o,Line_d,dis2,lnest);   	  
		}
	      else{
		track->PointToCircle(x,y,rho,xcen,ycen,dis,xest,yest);
		track->PointToCircle(xp,yp,rho,xcen,ycen,dis2,xest,yest);
	      }
	      if(dis<wiredis*1.2 && dis2<wiredis*1.2 )
		{tmpvec[streonum].push_back(k);}
	      
	    }      
	}//stereocluster

      if(  tmpvec[0].size()==0 ||tmpvec[1].size()==0
	   ||tmpvec[2].size()==0||tmpvec[3].size()==0) continue;
      if( tmpvec[0].size()*tmpvec[1].size()*tmpvec[2].size()*tmpvec[3].size()>50) continue;

      int is[4]={0};
      for(is[0]=0;is[0]<(int)tmpvec[0].size();is[0]++)
	{ for(is[1]=0;is[1]<(int)tmpvec[1].size();is[1]++)
	  { for( is[2]=0;is[2]<(int)tmpvec[2].size();is[2]++)
	    { for( is[3]=0;is[3]<(int)tmpvec[3].size();is[3]++)
	      {

		int allsthit=0;
		CDSTrack tmptrack=*track;
		for(int n=0;n<4;n++)
		  {
		    int stereo=0;
		    if(0<=n && n<=1)stereo=n+1;else if(2<=n&&n<=3)stereo=n+2;
		    int maxhit=(int)ClusterContainer[stereo]
 [ tmpvec[n][ is[n] ] ].nHit();
		    allsthit+=maxhit; 
		    for(int i=0;i< maxhit; i++)
		      {
			tmptrack.AddHit(*(ClusterContainer[stereo][tmpvec[n][is[n]]].CDC(i) ) );
		      }	    
		  }
		if(allsthit<=6) continue;
		tmptrackContainer.push_back(tmptrack);

	      }
	    }
	  }
	}
      

    }//track
  trackContainer=tmptrackContainer;
  tmptrackContainer.clear();
  return true;
  
}


void CDSTrackingMan::PointToCircle(const double &x,const double &y,const double &radius,const double &x_cen,const double &y_cen,double &dis,double &xest,double &yest)
{
  double rtmp=sqrt( (x-x_cen)*(x-x_cen)+(y-y_cen)*(y-y_cen) );
  double costmp=(x-x_cen)/rtmp;
  double sintmp=(y-y_cen)/rtmp;
  dis=fabs(rtmp-radius);

  xest=x_cen+radius*costmp;
  yest=y_cen+radius*sintmp;
}

void CDSTrackingMan::CalcBetaMass(TVector3 vertex,LocalTrack *beam,CDSTrack* cdctrack,
				  ConfMan *conf, const int &beam_pid, const double &tof,
				  double &beta, double &mass)
{
  double param[5];
  cdctrack->GetParameters(param);
  //	if(param[2]==0 ||param[4]==0   ) continue;
  //  double drho=param[0], phi0=param[1];
  double rho=cdctrack->Rho(), tlam=param[4];
  double mom = cdctrack->Momentum();

  double x1,y1,x2,y2;
  double gx,gy,gz,dx,dy,dz;
  conf->GetGeomMapManager()->GetGParam(CID_T0,gx,gy,gz,dx,dy,dz);
  double z_t0=gz,z_vtxb=vertex.Z();
  beam->XYPosatZ(z_t0,x1,y1);
  beam->XYPosatZ(z_vtxb,x2,y2);
  double beam_dis=sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z_t0-z_vtxb)*(z_t0-z_vtxb));
  //  std::cout<<"zt0,zvtxb="<<z_t0<<"\t"<<z_vtxb<<std::endl;  
  double beta_b=0.;
  if(beam_pid==0)  beta_b=1.0/sqrt(1.0*1.0+kpMass*kpMass);
  else if(beam_pid==1)  beta_b=1.0/sqrt(1.0*1.0+piMass*piMass);
  
  double time_beam=beam_dis/beta_b/(Const*100);
  
  //#####CDC dis#######		  	
  double mass1=0;	    
  TVector3 cdhvtx=cdctrack->CDHVertex();
  double phi1,phi2;
  phi1=cdctrack->CalcHelixPhi(vertex.x(),vertex.y());
  phi2=cdctrack->CalcHelixPhi(cdhvtx.x(),cdhvtx.y());

  double cdc_dis=rho*fabs(phi2-phi1)*sqrt(1+tlam*tlam);
  double beta_calc=cdc_dis/(tof-time_beam)/(Const*100);
  double calc_mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
  double calc_mass1=sqrt(calc_mass2);
  
  if(calc_mass1>0.07 && calc_mass1<0.22 ) mass1=piMass;
  else  if(calc_mass1>0.6 && calc_mass1<1.2 ) mass1=pMass;
  
  double beta1=fabs(mom)/sqrt(mass1*mass1+mom*mom);
  // double time_cdc=cdc_dis/beta1/(Const*100);
  /*  
  std::cout<<"beam_dis,cdc_dis,tof,time_beam,beta,mom,mass="
	   <<beam_dis<<"\t"<<cdc_dis<<"\t"<<tof<<"\t"<<time_beam
	   <<"\t"<<beta<<"\t"<<mom<<"\t"<<mass<<std::endl;
  */
  //####Calc time#####
  //  double calc_time=time_cdc+time_beam;
  //  double calc_time2=cdc_dis/(Const*100)+time_beam;
  
  //  double ctsub=cdh1->ctsub();
  //  double calc_time_u=time_cdc+time_beam+ctsub;
  //  double calc_time_d=time_cdc+time_beam-ctsub;
  //	    std::cout<<"CDH seg= "<<cdh1->seg()<<std::endl;	    
  beta = beta_calc;
  mass = calc_mass1;
}

void CDSTrackingMan::CalcVertex_beam(const int &trackid,BeamLineTrackMan *bltrackman,ConfMan *conf)
{
  for(int bi=0;bi<bltrackman->ntrackBPC();bi++)
    CalcVertex_beam(trackid,bltrackman->trackBPC(bi),conf,bi);
}

TrackVertex CDSTrackingMan::CalcVertex_beam(const int &trackid,LocalTrack *bc,ConfMan *conf, const int &bi)
{
  TrackVertex *Vertextmp=new TrackVertex();
     
  CDSTrack *cdstrack=&trackContainer[trackid];
  Vertextmp->SetTrackID1(trackid);
  Vertextmp->SetTrackID2(bi);
  
  TVector3 pos,fitpos1,fitpos2,tmp;
  TVector3 mom1,mom2;
  double param[5];
  cdstrack->GetParameters(param);

  if(param[2]==0 )
    {	  
      TVector3 x1,a1,x2,a2;     
      x1=cdstrack->LineOrigin(); a1=cdstrack->LineDirection();
      double xp1,yp1,xp2,yp2;
      bc->XYPosatZ(0,xp1,yp1);x2.SetXYZ(xp1,yp1,0);
      bc->XYPosatZ(20,xp2,yp2);a2.SetXYZ(xp2-xp1,yp2-yp1,20);
      a2=a2.Unit();
      x2=x2; a2=a2;
      double dl=0.2,dis;     
      cdstrack->LineToLine(x1,a1,x2,a2,dl,dis,fitpos1,tmp);
      TVector3 tmpd;
      tmpd=(tmp-fitpos1).Unit();
      fitpos2=fitpos1+dis*tmpd;
      
      Vertextmp->SetVertexPos1(fitpos1);
      Vertextmp->SetVertexPos2(fitpos2);
      if(!cdstrack->GetMomentum(fitpos1,mom1)){
#if 0
	std::cout << " err in calc vertex beam " <<std::endl;
	std::cout<<"x1"<<std::endl;
	x1.Print();
	std::cout<<"a1"<<std::endl;
	a1.Print();
	std::cout<<"x2"<<std::endl;
	x2.Print();
	std::cout<<"a2"<<std::endl;
	a2.Print();
#endif
      }
      mom2.SetXYZ(0,0,0.9);
      Vertextmp->SetMomentum1(mom1);
      Vertextmp->SetMomentum2(mom2);
      Vertextmp->SetVertexPos_mean(0.5*(fitpos1+fitpos2));	  
    }  
  else
    {	  
      TVector3 x2,a2;     
      double xp1,yp1,xp2,yp2;
      bc->XYPosatZ(param[3],xp1,yp1);x2.SetXYZ(xp1,yp1,param[3]);
      bc->XYPosatZ(param[3]+20,xp2,yp2);a2.SetXYZ((xp2-xp1),(yp2-yp1),param[3]+20);
      a2=a2.Unit();
      x2=x2; a2=a2;
      double dis;     
      if(!cdstrack->LineToHelix(x2,a2,param,fitpos2,fitpos1,dis))
      	std::cout << " line to helix failed " <<std::endl;
      Vertextmp->SetVertexPos1(fitpos1);
      Vertextmp->SetVertexPos2(fitpos2);
      if(!cdstrack->GetMomentum(fitpos1,mom1)){
#if 0
	std::cout << " err in calc vertex beam " <<std::endl;
	std::cout<<"fitpos1"<<std::endl;
	fitpos1.Print();
	std::cout<<"param: ";
	for(int i=0;i<5;i++) std::cout<<param[i]<<"\t";
	std::cout<<std::endl;
	std::cout<<"x2"<<std::endl;
	x2.Print();
	std::cout<<"a2"<<std::endl;
	a2.Print();
#endif
      }
      mom2.SetXYZ(0,0,0.9);
      Vertextmp->SetMomentum1(mom1);
      Vertextmp->SetMomentum2(mom2);
      Vertextmp->SetVertexPos_mean(0.5*(fitpos1+fitpos2));      
    }  
  Vertex_beam.push_back(*Vertextmp);
  TrackVertex vertex=*Vertextmp;
  delete Vertextmp;
  return vertex;
}


void CDSTrackingMan::CalcVertex(const int &id1,const  int &id2)
{
  TrackVertex *Vertextmp=new TrackVertex();
  int tmpid[2]={id1,id2};
  if(id1>id2) {tmpid[0]=id2;tmpid[1]=id1;}
   
  Vertextmp->SetTrackID1(tmpid[0]);
  Vertextmp->SetTrackID2(tmpid[1]);
  
  CDSTrack *track1=&trackContainer[tmpid[0]];
  CDSTrack *track2=&trackContainer[tmpid[1]];
 
  TVector3 pos,fitpos1,fitpos2,tmp;
  TVector3 mom1,mom2;
  double param[2][5];
  track1->GetParameters(param[0]);track2->GetParameters(param[1]);     
  
  if(param[0][2]==0 && param[1][2]==0)
    {
     
      TVector3 x1,a1,x2,a2;     
      x1=track1->LineOrigin(); a1=track1->LineDirection();
      x2=track2->LineOrigin();a2=track2->LineDirection();
      
      double dl=0.2,dis;     
      track1->LineToLine(x1,a1,x2,a2,dl,dis,fitpos1,tmp);
      TVector3 tmpd;
      tmpd=(tmp-fitpos1).Unit();
      fitpos2=fitpos1+dis*tmpd;
      
      Vertextmp->SetVertexPos1(fitpos1);
      Vertextmp->SetVertexPos2(fitpos2);

      //      track1->GetMomentum(fitpos1,mom1);
      //      track2->GetMomentum(fitpos2,mom2);
      mom1.SetXYZ(0,0,0);      mom2.SetXYZ(0,0,0);
      Vertextmp->SetMomentum1(mom1);
      Vertextmp->SetMomentum2(mom2);
      Vertextmp->SetVertexPos_mean(0.5*(fitpos1+fitpos2));
      
    }
  else
    {
      double dis=0,distmp=0;
      double rho[2],x_c[2],y_c[2],chi[2];  
      double posphi=0;
      chi[0]=track1->Chi();chi[1]=track2->Chi();
      rho[0]=track1->Rho();  rho[1]=track2->Rho();     
      x_c[0]=track1->CenterX(); x_c[1]=track2->CenterX();
      y_c[0]=track1->CenterY(); y_c[1]=track2->CenterY();
      for(int scale=0;scale<3;scale++)
	{
	  for(int n=0;n<20;n++)
	    {
	      double phi=posphi+(n-10)*3*pow(10,-scale)*param[0][2];
	      pos=track1->GetPosition(phi);
	      if(!track2->PointToHelix(pos,param[1],tmp,distmp) ) 
		{distmp=999;}
	      if(n==0) { dis=distmp; posphi=phi; fitpos2=tmp;}
	      else if(distmp<dis) { dis=distmp; posphi=phi; fitpos2=tmp;}
	    }
	}
      fitpos1=track1->GetPosition(posphi);
      Vertextmp->SetVertexPos1(fitpos1);
      Vertextmp->SetVertexPos2(fitpos2);


      track1->GetMomentum(fitpos1,mom1);
      track2->GetMomentum(fitpos2,mom2);
      Vertextmp->SetMomentum1(mom1);
      Vertextmp->SetMomentum2(mom2);
      Vertextmp->SetVertexPos_mean(0.5*(fitpos1+fitpos2));
      //  std::cout<<"vertex "<<posphi<<" "<<Vertex.x()<<" "<<Vertex.y()<<" "<<Vertex.z()<<" "<<std::endl;

    }
  
  Vertex.push_back(*Vertextmp);
  delete Vertextmp;  
}




bool CDSTrackingMan::SearchCDHHit(CDSHitMan *cdsMan)
{    
  //################################//
  //###########CDH Search ##########//	
  //################################//
  
  for( int itr=0; itr<(int)trackContainer.size(); itr++ )
    {
      double dis=999;  
      int CDHnum=-1;
      TVector3 CDHvertex;
      CDSTrack *track=&trackContainer[itr];
      if(!track->GoodFlag()) continue;
      if(cdsMan->nCDH()==0) return false;
      for( int ii=0; ii<cdsMan->nCDH(); ii++ )
	{
	  double param[5];
	  track->GetParameters(param);
	  double arho, ax_c, ay_c;
	  TVector3 line_o,line_d;
	  if(param[2]==0)
	    {
	      line_o=track->LineOrigin(); line_d=track->LineDirection();
	      line_o.SetXYZ(line_o.x(),line_o.y(),0);
	      line_d.SetXYZ(line_d.x(),line_d.y(),0);
	    }
	  else 
	    {
	      arho=track->Rho();
	      ax_c=track->CenterX();
	      ay_c=track->CenterY();
	    }
	  double cdh_x=cdsMan->CDH(ii)->x();
	  double cdh_y=cdsMan->CDH(ii)->y();
	  MathTools *tools=new MathTools();
	  double phi_cdh=tools->CalcDeg(cdh_x,cdh_y);
	  delete tools;
	  double phi_track=track->Theta();
	  double diff_phi=fabs(phi_track-phi_cdh);
	  if(diff_phi>180) diff_phi=360-diff_phi;
	  
	  if(diff_phi>90 ) continue;
	  double distmp,xest,yest;
	  TVector3 cdhp,nest;cdhp.SetXYZ(cdh_x,cdh_y,0);
	  if(param[2]==0)
	    {
	      track->PointToLine(cdhp,line_o,line_d,distmp,nest);
	      xest=nest.x(); yest=nest.y();
	    }
	  else
	    { 
	      track->PointToCircle(cdh_x,cdh_y,arho,ax_c,ay_c,distmp,xest,yest);
	    }
	  if(dis==999){CDHnum=ii;dis=distmp; CDHvertex.SetXYZ(xest,yest,0);}
	  else if(distmp<dis){CDHnum=ii;dis=distmp;CDHvertex.SetXYZ(xest,yest,0); }
	  
	}

      if( cdsMan->CDH(CDHnum)->CheckRange() )
	{
	  
	  track->SetCDHHit(*cdsMan->CDH(CDHnum) );
	  double z;
	  track->XYtoZ(CDHvertex.x(),CDHvertex.y(),z);
	  CDHvertex.SetZ(z);
	  track->SetCDHVertex(CDHvertex);
	}
    }
  
  return true;;
}


bool CDSTrackingMan::Execute(CDSHitMan *cdsMan,ConfMan *conf)
{

  /***********************/
  /**** Tracking **********/
  /***********************/
  Clear();
  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  selected_chi=CDSParam->GetMaxChi();
  int maxcdchit=CDSParam->GetMaxCDCHit();
  //  int cdctdc_u=CDSParam->GetUpperLimitCDCTDC();
  //  int cdctdc_l=CDSParam->GetLowerLimitCDCTDC();
  double magneticfield=CDSParam->GetMagneticField();
#if DEBUG
  double t0=clock();
#endif
  //###### event cut######
  int numCDC=0; 
  int numCDH=0;
  for( int i=0; i<cdsMan->nCDH(); i++ )
    {
      if( cdsMan->CDH(i)->CheckRange() )
	numCDH++;
    }

  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    for( int i=0; i<cdsMan->nCDC(layer); i++ ){
	numCDC++;
    }      
  }
  
  if( numCDC>maxcdchit ){
    Clear();
    return true;
  }

  SearchCluster(cdsMan,conf);
  FindingTrackCandidate(cdsMan,conf);

  //####cut too many container####//
  if(trackContainer.size()>100)
    {
      Clear();
      return true;
    }
#if DEBUG
  double t1=clock();
#endif
  //###############################
  //####### Axial fit#############
  //###############################
  if(magneticfield==0) FullAxalLineFit();
  else  FullCircleFit();

  SelectSharedHit();   
  FindingStereoHit(cdsMan,conf);
#if DEBUG
  double t2=clock();
#endif
  //########################
  //######Stereo Fittig######
  //########################
  if(magneticfield==0) FullStereoLineFit();
  else FullHelixFit();
  
  SelectSharedHit();
  Calc(conf);
  SearchCDHHit(cdsMan);
#if DEBUG
  double t3=clock();
#endif
    /***********************/
    /**********************/
#if DEBUG
  std::cout<<"##tracking time ## "<<std::endl;
  std::cout<<"search time : "<<t1-t0<<std::endl;
  std::cout<<"circle time : "<<t2-t1<<std::endl;
  std::cout<<"helix  time : "<<t3-t2<<std::endl;
#endif
  return true;

}


void CDSTrackingMan::FullCircleFit()
{
  //###############################
  //####### circle fit#############
  //###############################
  for(int i=0;i<(int)trackContainer.size();i++)
    {
      CDSTrack *track=&trackContainer[i];
      if(!track->CalcInitialParameters())
	{
#if DEBUG
	  std::cout<<"ini par!"<<std::endl;
#endif
	  track->SetGoodFlag(false);continue;
	}
      
      if(!track->FirstCircleFitting() ) 
	{
#if DEBUG
	  std::cout<<"first cir track!!"<<std::endl;
#endif
	  track->SetGoodFlag(false);  continue;
	}
      
#if DEBUG
      std::cout<<" firstcir "<<track->Chi();
#endif
      
      bool flag=true;
      for(int trial=0;trial<10;trial++) { if(!track->CircleFitting() ) flag=false;}
      if(!flag)  
      	{
#if DEBUG
	  std::cout<<"cir track!!"<<std::endl;
#endif
	  track->SetGoodFlag(false); continue;
	}     
#if DEBUG
      std::cout<<" cir "<<track->Chi();
#endif
      
      
    }

}


void CDSTrackingMan::FullAxalLineFit()
{
  //###############################
  //####### circle fit#############
  //###############################
  for(int i=0;i<(int)trackContainer.size();i++)
    {
      CDSTrack *track=&trackContainer[i];
      
      if(!track->FirstLineFitting() ) 
	{
#if DEBUG
	  std::cout<<"first line track!!"<<std::endl;
#endif
	  track->SetGoodFlag(false);  continue;
	}
      
#if DEBUG
      std::cout<<" firstline "<<track->Chi();
#endif
      
      bool flag=true;
      for(int trial=0;trial<10;trial++) { if(!track->LineFitting() ) flag=false;}
      if(!flag)  
      	{
#if DEBUG
	  std::cout<<"line track!!"<<std::endl;
#endif
	  track->SetGoodFlag(false); continue;
	}     
#if DEBUG
      std::cout<<" line "<<track->Chi();
#endif

    }

}


void CDSTrackingMan::FullStereoLineFit()
{
  //#############################
  //######StereoLine Fittig######
  //#############################
  for(int i=0;i<(int)trackContainer.size();i++)
    {

      CDSTrack *track=&trackContainer[i];            
      if(!track->FirstStereoLineFitting() ) 
	{
#if DEBUG
	  std::cout<<"first sline track!!"<<std::endl; 
#endif
	  track->SetGoodFlag(false);  continue;
	}
      
#if DEBUG
      std::cout<<" firstsl "<<track->Chi();
#endif
      /*
	if( track->Chi()>2*selected_chi ) 
	  {
#if DEBUG
	    std::cout<<"first sl chi!"<<std::endl; 
#endif 
	    track->SetGoodFlag(false);  continue;
	  }
      */		   
    track->SetHitPos();       
    int test=0;
    bool flag=true;
    for(int trial=0;trial<5;trial++) { if(!track->StereoLineFitting() ) {flag=false;test=trial;break;}}
    if(!flag) 
      {
#if DEBUG
	std::cout<<"sl track ! "<<test<<std::endl; track->SetGoodFlag(false); 
#endif
	continue;
      }
    
    int numt=0;
    double Chi=999,Chitmp=999;
    bool hflag1=true;	
    while(  numt<20 )
      {
	Chi=Chitmp;
	for(int trial=0;trial<10;trial++)
	  { if(!track->StereoLineFitting() ){hflag1=false; break;}  }	
	Chitmp=track->Chi();
	if( (Chi-Chitmp)<0.1 ) break;
	numt++;
      }
    if(!hflag1) { continue;}

#if DEBUG    
    std::cout<<" stline "<<track->Chi()<<std::endl;
#endif
	
    
    }
	
}


void CDSTrackingMan::FullHelixFit()
{
  //########################
  //######Helix Fittig######
  //########################

  for(int i=0;i<(int)trackContainer.size();i++)
    {
#if DEBUG
      double ht0=clock();
#endif
      CDSTrack *track=&trackContainer[i];            
      if(!track->FirstHelixFitting() ) 
	{
#if DEBUG
	  std::cout<<"first hel track!!"<<std::endl; 
#endif
	  track->SetGoodFlag(false);  continue;
	}
      
#if DEBUG
      std::cout<<" firsthel "<<track->Chi();
#endif
	if( track->Chi()>2*selected_chi ) 
	  {
#if DEBUG
	    std::cout<<"firsthel chi!"<<std::endl; 
#endif 
	    track->SetGoodFlag(false);  continue;
	  }
#if DEBUG
	double ht1=clock();
#endif
    track->SetHitPos();       
    int test=0;
    bool flag=true;
    for(int trial=0;trial<5;trial++) { if(!track->HelixFitting() ) {flag=false;test=trial;break;}}
    if(!flag) 
      {
#if DEBUG
	std::cout<<"hel track ! "<<test<<std::endl; track->SetGoodFlag(false); 
#endif
	continue;
      }
#if DEBUG
      double ht2=clock();
#endif
    int numt=0;
    double Chi=999,Chitmp=999;
    bool hflag1=true;	
    while(  numt<20 )
      {
	Chi=Chitmp;
	for(int trial=0;trial<10;trial++)
	  { if(!track->HelixFitting() ){hflag1=false; break;}  }	
	Chitmp=track->Chi();
	if( (Chi-Chitmp)<0.1 ) break;
	numt++;
      }
    if(!hflag1) { continue;}
#if DEBUG    
    double ht3=clock();
    std::cout<<" helix "<<track->Chi()<<std::endl;

    std::cout<<"     helix track time : "<<i<<std::endl;
    std::cout<<"        first helix fit: "<<ht1-ht0<<std::endl;
    std::cout<<"        helix fit: "<<ht2-ht1<<std::endl;
    std::cout<<"        while helix fit: "<<ht3-ht2<<std::endl;
#endif
    }

}

void CDSTrackingMan::Calc(ConfMan *conf)
{
  Vertex_beam.clear();
  Vertex.clear();
  GoodTrack.clear();

  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  selected_chi=CDSParam->GetMaxChi();
  //  int maxcdchit=CDSParam->GetMaxCDCHit();
  //  int cdctdc_u=CDSParam->GetUpperLimitCDCTDC();
  //  int cdctdc_l=CDSParam->GetLowerLimitCDCTDC();
  //  double magneticfield=CDSParam->GetMagneticField();


  for(int i=0;i<(int)trackContainer.size();i++)
    {
      CDSTrack *track=&trackContainer[i];
      track->SetGoodFlag(false);
      if(track->FittingLevel()==4&&track->Chi()<selected_chi) 
	{
	  track->SetGoodFlag(true);
	  GoodTrack.push_back(i);
	  track->Calc(conf);
	}

    }

  for(int i=0;i<(int)GoodTrack.size();i++)
    {
      //      CDSTrack *track=&trackContainer[GoodTrack[i]];
      for(int ii=i+1;ii<(int)GoodTrack.size();ii++)
	{
	  //	  CDSTrack *track2=&trackContainer[GoodTrack[ii]];
	  CalcVertex(GoodTrack[i],GoodTrack[ii]);
	}
    }

  return;
}


bool CDSTrackingMan::GetVertex(const int &trackid1,const int &trackid2,TrackVertex &vertex)
{
  int gid[2]={-1};

  for(int i=0;i<(int)GoodTrack.size();i++)
    {
      if(trackid1==GoodTrack[i]) gid[0]=i;
      if(trackid2==GoodTrack[i]) gid[1]=i;
    }
  if(gid[0]==-1 || gid[1]==-1 || gid[0]==gid[1]) return false;
  if(gid[0]>gid[1]) 
    {
      int tmp=gid[0];
      gid[0]=gid[1];
      gid[1]=tmp;
    }

  int nGoodTrack=(int)GoodTrack.size();
  // int gidtmp[2]={0};
  int vertexidtmp=0;
  int vertexid=0;
  for(int id1=0;id1<nGoodTrack;id1++)
    { for(int id2=id1+1;id2<nGoodTrack;id2++)
      {
	if(id1==gid[0]&&id2==gid[1]) vertexid=vertexidtmp; 
	vertexidtmp++;
      } 
    }
  vertex=Vertex[vertexid];
  return true;
  
}


bool CDSTrackingMan::GetVertex_beam(const int &trackid,const int &beamid, TrackVertex &vertex)
{
  if(Vertex_beam.size()==0) return false;

  int nbeam=(double)Vertex_beam.size()/(double)GoodTrack.size();
  int gid[2]={-1,-1};
  if(nbeam<beamid) return false;
  int vertexidtmp=0;
  int vertexid=0;
  for(int bi=0;bi<nbeam;bi++)
    {
      for(int i=0;i<(int)GoodTrack.size();i++)
	{
	  if(trackid==GoodTrack[i] && beamid==bi) {gid[0]=i;gid[1]=bi;vertexid=vertexidtmp;}
	  vertexidtmp++;
	}
    }
  if(gid[0]==-1 ||gid[1]==-1 ) return false;

  vertex=Vertex_beam[vertexid];
  return true;
  
}

//###########################
//#######TrackVertex#########
//###########################


TrackVertex::TrackVertex(): TObject()
{
  trackid[0]=-999;
  trackid[1]=-999;
  Vertex_mean.SetXYZ(-999,-999,-999);
  Vertex1.SetXYZ(-999,-999,-999);
  Vertex2.SetXYZ(-999,-999,-999);
}


void TrackVertex::Clear()
{
  trackid[0]=-999;
  trackid[1]=-999;
  Vertex_mean.SetXYZ(-999,-999,-999);
  Vertex1.SetXYZ(-999,-999,-999);
  Vertex2.SetXYZ(-999,-999,-999);
}


TVector3 TrackVertex::GetVertexPos(const int &aid)
{
  if(aid==trackid[0]) return Vertex1;
  else if(aid==trackid[1]) return Vertex2;
  else 
    {
      TVector3 tmp;tmp.SetXYZ(-999,-999,-999);
      std::cout<<"No ID in This TrackVertex! "<<std::endl;
      return tmp;
    }
}


TVector3 TrackVertex::GetMomentum(const int &aid)
{
  if(aid==trackid[0]) return Momentum1;
  else if(aid==trackid[1]) return Momentum2;
  else 
    {
      TVector3 tmp;tmp.SetXYZ(-999,-999,-999);
      std::cout<<"No ID in This TrackVertex! "<<std::endl;
      return tmp;
    }
}
