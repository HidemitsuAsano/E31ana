#include "CDSTrackingMan.h"
#define DEBUG 0
#define DEBUG2 0
#define NOMINUIT 0
#define TESTSHORT 0
#define REMOVECLUSTER 1
ClassImp(CDSTrackingMan);
ClassImp(TrackVertex);

const double Inclradius=20.40,Outclradius=47.78,cellsize=0.845*2.;
//const double ResolutionOfCDC=0.02;

CDSTrackingMan::CDSTrackingMan(): TObject()
{
  //  Clear();
  KilledLayer=0;
  selected_chi=50;
}

void CDSTrackingMan::Clear()
{
  trackContainer.clear();
  Vertex.clear();
  Vertex_beam.clear();
  goodtrackIDContainer.clear();
  shorttrackIDContainer.clear();
  for(int i=0;i<7;i++) ClusterContainer[i].clear();
}

void CDSTrackingMan::AddTrack(const CDSTrack &atrack )
{
  trackContainer.push_back(atrack);
}
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
	  double x1=cdc1->wx();
	  double y1=cdc1->wy();
	  bool hit_layer2flag=false;
	  for(int ii=0;ii<cdsMan->nCDC(layer2);ii++)
	    {
	      if(layer2==KilledLayer) break;
	      CDCHit *cdc2=cdsMan->CDC(layer2,ii);
	      int wire2= cdc2->wire();
	      int tdc2 = cdc2->tdc();
	      if( !(cdctdc_l<tdc2 && tdc2<cdctdc_u) ) continue;	  	      
	      double x2=cdc2->wx(); 
	      double y2=cdc2->wy();
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
	  ClusterContainer[num_slayer-1][n].Calc(cdsMan);
	}
    }//slayer  
  return true;
}

//search cluseter

void CDSTrackingMan::DeleteUsedCluster(CDSHitMan *cdsMan)
{
  for(int num_slayer=1;num_slayer<=7;num_slayer++)
    {
      for(int icl=0;icl<(int)ClusterContainer[num_slayer-1].size();icl++)
	{
	  HitCluster *cl;
	  cl=&ClusterContainer[num_slayer-1][icl];
	  for(int ihit=0;ihit<cl->nHit();ihit++)
	    {
	      
	      CDCHit *hit=cl->CDC(cdsMan,ihit);
	      int layer=hit->layer();
	      int hid=hit->hid();
	      int wire=hit->wire();
	      bool flagshared=false;
	      for(int itr=0;itr<(int)trackContainer.size();itr++)
		{
		  CDSTrack *track1=&trackContainer[itr];
		  //		  if( !track1->GoodFlag() ) continue;
		  for(int ih1=0;ih1<track1->nTrackHit(layer);ih1++)
		    {
		      CDCHit *hit1=track1->hit(cdsMan,layer,ih1);
		      int wire1=hit1->wire();
		      int hid1=hit1->hid();
		      if(wire1==wire && hid==hid1) 
			{flagshared=true;break;}
		    }
		  if(flagshared) break;		  
		}// track	      
	      if(flagshared) cl->DeleteHit(ihit);	      
	    }//hit in cluster
	}//end of icl      
      //###Delete cluster with too less Hit
      HitClusterContainer::iterator it;
      for(it=ClusterContainer[num_slayer-1].end()-1;it!=ClusterContainer[num_slayer-1].begin()-1;--it)
	{
	  HitCluster *cl=&*it;
	  int nhit=cl->nHit();
	  if(num_slayer==1 && nhit<3)  ClusterContainer[num_slayer-1].erase(it);
	  if(num_slayer>1 && nhit<2)   ClusterContainer[num_slayer-1].erase(it);
	}
    }//end of num slayer

  return;
}    
//DeleteSharedCluster


bool CDSTrackingMan::FindingShortTrackCandidate(CDSHitMan *cdsMan,ConfMan *conf)
{       
    CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
    double magneticfield=CDSParam->GetMagneticField();

    int numincl=ClusterContainer[0].size();
    for(int incl=0;incl<numincl;incl++)
      {
	CDSTrack TmpTrack;
	TVector3 inpos;
	inpos.SetXYZ(ClusterContainer[0][incl].Meanx(),
		     ClusterContainer[0][incl].Meany(),
		     0);

	for(int midcl=0;midcl<(int)ClusterContainer[3].size(); midcl++)
	  {

	    double def_deg=fabs(ClusterContainer[0][incl].Meanphi()
				-ClusterContainer[3][midcl].Meanphi());
	    if(def_deg>180) def_deg=360-def_deg;
	    if(def_deg>60 ) continue;
	    TVector3 midpos;
	    midpos.SetXYZ(ClusterContainer[3][midcl].Meanx(),
			 ClusterContainer[3][midcl].Meany(),
			 0);

	    TVector3 cenpos;
	    cenpos=0.5*(inpos+midpos);
	    double rlimit=0;
	    if(magneticfield==0) rlimit=1.7*4;//4cell
	    else    rlimit=(cenpos-inpos).Mag();

	    double x_in=ClusterContainer[0][incl].Meanx();
	    double y_in=ClusterContainer[0][incl].Meany();
	    double x_mid=ClusterContainer[3][midcl].Meanx();
	    double y_mid=ClusterContainer[3][midcl].Meany();
	
	    double a=1,b,c; //a*x+b*y+c=0
	    b=-(x_in-x_mid)/(y_in-y_mid);
	    c=-(x_in+b*y_in);
	    TmpTrack.SetABC(a,b,c);
	    TmpTrack.SetClusterX(1,x_in);TmpTrack.SetClusterY(1,y_in);	
	    TmpTrack.SetClusterX(4,x_mid);TmpTrack.SetClusterY(4,y_mid);

	    //######stereo cluster############
	    for(int st1cl=0;st1cl<(int)ClusterContainer[1].size(); st1cl++)
	      {
		TVector3 st1pos;
		st1pos.SetXYZ(ClusterContainer[1][st1cl].Meanx(),
			      ClusterContainer[1][st1cl].Meany(),
			      0);
		if( (st1pos-cenpos).Mag()>rlimit) continue;
		TmpTrack.SetClusterX(2,ClusterContainer[1][st1cl].Meanx() );
		TmpTrack.SetClusterY(2,ClusterContainer[1][st1cl].Meany() );
		for(int st2cl=0;st2cl<(int)ClusterContainer[2].size(); st2cl++)
		  {
		    TVector3 st2pos;
		    st2pos.SetXYZ(ClusterContainer[2][st2cl].Meanx(),
				  ClusterContainer[2][st2cl].Meany(),
				  0);

		    if( (st2pos-cenpos).Mag()>rlimit) continue;
		    if( (st2pos-st1pos).Mag()>10*1.7) continue;//10cell size

		    TVector3 stmean; stmean=0.5*(st1pos+st2pos);
		    TVector3 tmp[2];tmp[0]=(midpos-inpos).Unit();
		    tmp[1]=inpos+(tmp[0]*stmean-tmp[0]*inpos)*tmp[0];
		    
		    double jitter=(tmp[1]-inpos).Mag();
		    double theta_o=0,theta_m=0;
		    theta_o=MathTools::CalcDeg(midpos.x()-inpos.x(),midpos.y()-inpos.y());
		    theta_m=MathTools::CalcDeg(stmean.x()-inpos.x(),stmean.y()-inpos.y());
		    if(fabs(theta_o-theta_m)>330) {theta_o-=360;theta_m-=360;}
		    if(fabs(theta_o-theta_m)<0) jitter=jitter;
		    else jitter=-1*jitter;      
		    double theta=MathTools::CalcDeg(cenpos.x(),cenpos.y());
		    TmpTrack.SetJitter(jitter);
		    TmpTrack.SetTheta(theta);
		    TmpTrack.SetClusterX(3,ClusterContainer[2][st2cl].Meanx());
		    TmpTrack.SetClusterY(3,ClusterContainer[2][st2cl].Meany());
		    for(int i=0;i< ClusterContainer[0][incl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[0][incl].CDC(cdsMan,i) ) );
		      }
		    for(int i=0;i<(int)ClusterContainer[1][st1cl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[1][st1cl].CDC(cdsMan,i) ) );
		      }
		    for(int i=0;i<(int)ClusterContainer[2][st2cl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[2][st2cl].CDC(cdsMan,i) ) );
		      }

		    for(int i=0;i<(int)ClusterContainer[3][midcl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[3][midcl].CDC(cdsMan,i) ) );
		      }
		    TmpTrack.SetShort(true);
		    trackContainer.push_back(TmpTrack);
		    TmpTrack.Clear();       		       	      
		    
		  }//st2		
	      }//st1
	  }//midcl
      }//incl

    shorttrackIDContainer.clear();
    for(int i=0;i<(int)trackContainer.size();i++)
      {
	CDSTrack *track=&trackContainer[i];
	if(track->IsShort())  shorttrackIDContainer.push_back(i);
      }
    return true;
}	
//Findingshort track    


bool CDSTrackingMan::FindingTrackCandidate(CDSHitMan *cdsMan,ConfMan *conf)
{       
    CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
    double magneticfield=CDSParam->GetMagneticField();

    double JITTER;
    if(magneticfield==0) JITTER=3.0;
    //    else JITTER=15.0;
    else JITTER=25.0;

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
	    TmpTrack.SetClusterX(7,x_out);TmpTrack.SetClusterY(7,y_out);

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
		    TmpTrack.SetClusterX(4,x_mid);TmpTrack.SetClusterY(4,y_mid);	
		    double theta_o=0,theta_m=0;
		    theta_o=MathTools::CalcDeg(x_out-x_in,y_out-y_in);
		    theta_m=MathTools::CalcDeg(x_mid-x_in,y_mid-y_in);
		    if(fabs(theta_o-theta_m)>330) {theta_o-=360;theta_m-=360;}
		    if(fabs(theta_o-theta_m)<0) dis_mid=dis_mid;
		    else dis_mid=-1*dis_mid;
		    TmpTrack.SetJitter(dis_mid);
		    double theta=MathTools::CalcDeg(x_mid,y_mid);
		    TmpTrack.SetTheta(theta);		    
		    for(int i=0;i< ClusterContainer[0][incl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[0][incl].CDC(cdsMan,i) ) );
		      }
		    for(int i=0;i<(int)ClusterContainer[3][midcl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[3][midcl].CDC(cdsMan,i) ) );
		      }
		    for(int i=0;i< ClusterContainer[6][outcl].nHit(); i++)
		      {
			TmpTrack.AddHit(*(ClusterContainer[6][outcl].CDC(cdsMan,i) ) );
		      }
		    num_midcl++;
		    TmpTrack.SetShort(false);
		    trackContainer.push_back(TmpTrack);
		    TmpTrack.Clear();       		       
		  }
	      }//midcl 
	  }//outcl
      }//incl	
    return true;
}

// void CDSTrackingMan::SelectSharedHit(CDSHitMan *cdsMan)
// { 
//   double chimin=0;
//   for(int itr=0;itr<(int)trackContainer.size()-1;itr++)
//     {
//       CDSTrack *track1=&trackContainer[itr];
//       double chi1=track1->Chi();
//       //      if( !(track1->FittingLevel()==2 || track1->FittingLevel()==4) ) continue;
//       bool flaghit1=true;
//       for(int layer=1;layer<=NumOfCDCLayers;layer++)
// 	{ 
// 	  for(int ih1=0;ih1<track1->nTrackHit(layer);ih1++)
// 	    {
// 	      CDCHit *hit1=track1->hit(cdsMan,layer,ih1);
// 	      int wire1=hit1->wire();
// 	      for(int itr2=itr+1;itr2<(int)trackContainer.size();itr2++)
// 		{ 
// 		  if(!flaghit1) break;
// 		  CDSTrack *track2=&trackContainer[itr2];
// 		  double chi2=track2->Chi();		
// 		  //	  if( !(track2->FittingLevel()==2 || track2->FittingLevel()==4) ) continue;
		  
// 		  for(int ih2=0;ih2<track2->nTrackHit(layer);ih2++)
// 		    {
// 		      CDCHit *hit2=track2->TrackHit(cdsMan,layer,ih2);
// 		      int wire2=hit2->wire();
// 		      if(wire1==wire2)
// 			{
// 			  if(chi1<chi2 && chi2 > chimin) track2->DeleteHit(layer,ih2);
// 			  else flaghit1=false;
// 			}
// 		    }
// 		}//track2
// 	      if(!flaghit1 && chi1 > chimin) track1->DeleteHit(layer,ih1);
// 	    }//Hit1
// 	}//layer
//     }//track1
  
//   //###Delete track with too less Hit
//   TrackContainer::iterator ittr;
//   for(ittr=trackContainer.end()-1;ittr!=trackContainer.begin()-1;--ittr)
//     {
//       int hitlayer=0;
//       CDSTrack *track=&*ittr;
//       for(int layer=1;layer<=NumOfCDCLayers;layer++)
// 	{ 
// 	  if(track->nTrackHit(layer)>0) hitlayer++;
// 	}
//       //      std::cout<<"trackLevel "<<track->FittingLevel()<<" chi "<<track->Chi()<<" hitlayer "<<hitlayer<<std::endl;
//       //      if( !(track->FittingLevel()==2 || track->FittingLevel()==4) )   trackContainer.erase(ittr);//delete track
//       if(track->FittingLevel()==2 &&( hitlayer<= 4 || track->Chi()>3*selected_chi )) trackContainer.erase(ittr);//delete track
//       else if(track->FittingLevel()==3 && hitlayer<= 10) trackContainer.erase(ittr);//delete track
//       else if(track->FittingLevel()==4 && hitlayer<= 10) trackContainer.erase(ittr);//delete track
//     }
// }

void CDSTrackingMan::SelectSharedHit(CDSHitMan *cdsMan)
// new Aug 2013
{ 
  //  double chimin=0;
  //select best chi2 track
  TrackContainer tmptrackContainer;
  while(trackContainer.size()>0){
#if DEBUG2
    std::cout<<"container size = "<<trackContainer.size()<<std::endl;
    for(int itr=0;itr<(int)trackContainer.size();itr++)
      {
	CDSTrack *track=&trackContainer[itr];
	std::cout<<track->Chi()<<"  ";
	track->Print();
      }    
#endif
    double tmpchi=99999;
    int tmptr=0;
    for(int itr=0;itr<(int)trackContainer.size();itr++)
      {
	CDSTrack *track1=&trackContainer[itr];
	double chi1=track1->Chi();
	if(chi1<tmpchi){
	  tmpchi=chi1;
	  tmptr=itr;
	}
      }
    CDSTrack *track1=&trackContainer[tmptr];  
    tmptrackContainer.push_back(*track1);
    //      if( !(track1->FittingLevel()==2 || track1->FittingLevel()==4) ) continue;
    for(int layer=1;layer<=NumOfCDCLayers;layer++)
      { 
	for(int ih1=0;ih1<track1->nTrackHit(layer);ih1++)
	  {
	    CDCHit *hit1=track1->hit(cdsMan,layer,ih1);
	    int wire1=hit1->wire();
	    for(int itr2=0;itr2<(int)trackContainer.size();itr2++)
	      { 
		CDSTrack *track2=&trackContainer[itr2];
		for(int ih2=0;ih2<track2->nTrackHit(layer);ih2++)
		  {
		    CDCHit *hit2=track2->TrackHit(cdsMan,layer,ih2);
		    int wire2=hit2->wire();
		    if(wire1==wire2)
		      track2->DeleteHit(layer,ih2);
		  }
	      }//track2
	  }//Hit1
      }//layer
    
    //###Delete track with too less Hit
    TrackContainer::iterator ittr;
    for(ittr=trackContainer.end()-1;ittr!=trackContainer.begin()-1;--ittr)
      {
	int hitlayer=0;
	CDSTrack *track=&*ittr;
	int slay[7];
	for(int i=0;i<7;i++) slay[i]=0;
	for(int layer=1;layer<=NumOfCDCLayers;layer++)
	  { 
	    if(track->nTrackHit(layer)>0){
	      hitlayer++;
	      if(layer>=1 &&layer<=3 ) slay[0]++;
	      else if(layer>=4 &&layer<=5 ) slay[1]++;
	      else if(layer>=6 &&layer<=7 ) slay[2]++;
	      else if(layer>=8 &&layer<=9 ) slay[3]++;
	      else if(layer>=10&&layer<=11) slay[4]++;
	      else if(layer>=12&&layer<=13) slay[5]++;
	      else if(layer>=14&&layer<=15) slay[6]++;
	    }
	  }
	//      std::cout<<"trackLevel "<<track->FittingLevel()<<" chi "<<track->Chi()<<" hitlayer "<<hitlayer<<std::endl;
	//      if( !(track->FittingLevel()==2 || track->FittingLevel()==4) )   trackContainer.erase(ittr);//delete track
	if(track->FittingLevel()==2 &&( hitlayer<= 4 || track->Chi()>3*selected_chi )) trackContainer.erase(ittr);//delete track
	//	if(track->FittingLevel()==2 &&( hitlayer<= 4 )) trackContainer.erase(ittr);//delete track
	else if(track->FittingLevel()==3 && hitlayer<= 10) trackContainer.erase(ittr);//delete track
	else if(track->FittingLevel()==4 && hitlayer<= 10) trackContainer.erase(ittr);//delete track
	else if(slay[0]==0||slay[3]==0||slay[6]==0) trackContainer.erase(ittr);//delete track
	else if(track->FittingLevel()>2 && slay[1]+slay[2]+slay[4]+slay[5]<6) trackContainer.erase(ittr);//delete track
      }
  }
  for(int itr=0;itr<(int)tmptrackContainer.size();itr++)
    {
      CDSTrack *track1=&tmptrackContainer[itr];
      trackContainer.push_back(*track1);
    }
}
//SelectSharedHit

void  CDSTrackingMan::SelectSharedHitforShort(CDSHitMan *cdsMan)
{
  for(int itr=0;itr<(int)trackContainer.size()-1;itr++)
    {
      CDSTrack *track1=&trackContainer[itr];
      if(!track1->IsShort()) continue;
      double chi1=track1->Chi();
      if( !(track1->FittingLevel()==2 || track1->FittingLevel()==4) ) continue;
      for(int layer=1;layer<=NumOfCDCLayers;layer++)
	{ for(int ih1=0;ih1<track1->nTrackHit(layer);ih1++)
	  {
	    bool flaghit1=true;
	    CDCHit *hit1=track1->hit(cdsMan,layer,ih1);
	    int wire1=hit1->wire();
	    for(int itr2=itr+1;itr2<(int)trackContainer.size();itr2++)
	      { 
		if(!flaghit1) break;
		CDSTrack *track2=&trackContainer[itr2];
		if(!track2->IsShort()) continue;
		double chi2=track2->Chi();		
		if( !(track2->FittingLevel()==2 || track2->FittingLevel()==4) ) continue;

		for(int ih2=0;ih2<track2->nTrackHit(layer);ih2++)
		  {
		    CDCHit *hit2=track2->TrackHit(cdsMan,layer,ih2);
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
      if(!track->IsShort()) continue;
      for(int layer=1;layer<=NumOfCDCLayers;layer++)
	{ 
	  if(track->nTrackHit(layer)>0) hitlayer++;
	}
      //      std::cout<<"trackLevel "<<track->FittingLevel()<<" chi "<<track->Chi()<<" hitlayer "<<hitlayer<<std::endl;
      if( !(track->FittingLevel()==2 || track->FittingLevel()==4) )
	{ 
#if DEBUG2
	  std::cout<<"flag false "<<track->FittingLevel()<<std::endl;
#endif
	  trackContainer.erase(ittr);//delete track
	}
      if(track->FittingLevel()==2 &&(hitlayer<8 || track->Chi()>10*selected_chi )) 
	{
#if DEBUG2
	  std::cout<<"chi or hit false "<<hitlayer<<" "<<track->Chi()<<std::endl;
#endif
	  trackContainer.erase(ittr);//delete track
	}
      else if(track->FittingLevel()==4 && hitlayer<8 )
	{
#if DEBUG2
	  std::cout<<"f  hit false "<<hitlayer<<" "<<track->Chi()<<std::endl;
#endif
	  
	  trackContainer.erase(ittr);//delete track
	}
    }
  shorttrackIDContainer.clear();
  for(int i=0;i<(int)trackContainer.size();i++)
    {
      CDSTrack *track=&trackContainer[i];
      if(track->IsShort())  shorttrackIDContainer.push_back(i);
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
#if DEBUG2
      std::cout<<"---finding stereo hit -----  " <<itr<<std::endl;
#endif
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
		MathTools::PointToLine(pos,Line_o,Line_d,dis,lnest); 
		MathTools::PointToLine(posp,Line_o,Line_d,dis2,lnest);   	  
	      }
	      else{
		MathTools::PointToCircle(x,y,rho,xcen,ycen,dis,xest,yest);
		MathTools::PointToCircle(xp,yp,rho,xcen,ycen,dis2,xest,yest);
	      }
	      if(dis<wiredis*1.2 && dis2<wiredis*1.2 )
		{tmpvec[streonum].push_back(k);}
	      
	    }      
	}//stereocluster
      int tmpn=0;
      for(int ivec=0;ivec<4;ivec++)
	if(tmpvec[ivec].size()>0) tmpn++;
#if DEBUG2
      for(int ivec=0;ivec<4;ivec++)
	std::cout<<ivec<<" "<<tmpvec[ivec].size()<<std::endl;
      std::cout<<"tmpn "<<tmpn<<std::endl;
#endif
      if(tmpn!=4) continue;
      if( tmpvec[0].size()*tmpvec[1].size()*tmpvec[2].size()*tmpvec[3].size()>100) continue;

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
			tmptrack.AddHit(*(ClusterContainer[stereo][tmpvec[n][is[n]]].CDC(cdsMan,i) ) );
		      }	    
		  }
		if(allsthit<6) continue;//track eff 20130218
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
//finding stereo hit

bool CDSTrackingMan::FindingAdditionalHit(CDSHitMan *cdsMan,ConfMan *conf)
{   
  // CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  // double magfield=CDSParam->GetMagneticField();
    
  // //########Hit streo wire###########
  // TrackContainer tmptrackContainer;
  // for(int itr=0;itr<(int)shorttrackContainer.size();itr++)
  //   {
  //     std::vector<int> tmpvec[2];
  //     CDSTrack *track=&shorttrackContainer[itr];
  //     if(track->FittingLevel()!=2) continue;
  //     double rho=track->Rho();
  //     double xcen=track->CenterX();
  //     double ycen=track->CenterY();
  //     TVector3 Line_o,Line_d;
  //     if(magfield==0)
  // 	{
  // 	  Line_o=track->LineOrigin(); Line_d=track->LineDirection();
  // 	}

  //     for( int streonum=0; streonum<2; streonum++ )
  // 	{
  // 	  int streoslayer=0;
  // 	  if(0<=streonum && streonum<=1) streoslayer=streonum+4; 

  // 	  for( int k=0; k<(int)ClusterContainer[streoslayer].size(); k++ )
  // 	    {
  // 	      // int tdc = cdsMan->CDC(streolayer,k)->tdc();
  // 	      // if( !(400<tdc && tdc<900)) continue;
  // 	      //double radius=cdsMan->CDC(streolayer,k)->radius();
  // 	      HitCluster *cluster=&ClusterContainer[streoslayer][k];
  // 	      double phi=cluster->Meanphi();
  // 	      double def_deg_st=fabs(phi-track->Theta());
  // 	      if(def_deg_st>180 ) def_deg_st=360-def_deg_st;
  // 	      if(def_deg_st>90) continue;
	    
  // 	      double x=cluster->Meanx();
  // 	      double y=cluster->Meany();
  // 	      double xp=cluster->Meanxp();
  // 	      double yp=cluster->Meanyp();
	      
  // 	      double wiredis = sqrt( (x-xp)*(x-xp)+(y-yp)*(y-yp) );
  // 	      double dis,dis2,xest,yest;
  // 	      TVector3 pos,posp,lnest; pos.SetXYZ(x,y,0);posp.SetXYZ(xp,yp,0);
  // 	      if(magfield==0){
  // 		MathTools::PointToLine(pos,Line_o,Line_d,dis,lnest); 
  // 		MathTools::PointToLine(posp,Line_o,Line_d,dis2,lnest); 
  // 		}
  // 	      else{
  // 		MathTools::PointToCircle(x,y,rho,xcen,ycen,dis,xest,yest);
  // 		MathTools::PointToCircle(xp,yp,rho,xcen,ycen,dis2,xest,yest);
  // 	      }
  // 	      if(dis<wiredis*1.2 && dis2<wiredis*1.2 )
  // 		{tmpvec[streonum].push_back(k);}
	      
  // 	    }      
  // 	}//stereocluster

  //     if( tmpvec[0].size()==0 )
  // 	{
  // 	  tmptrackContainer.push_back(*track);		
  // 	  continue;
  // 	}
  //     int is[2];
  //     for(is[0]=0;is[0]<(int)tmpvec[0].size();is[0]++)
  // 	{
  // 	  int allsthit=0;
  // 	  CDSTrack tmptrack=*track;
  // 	  int stereo=0;
  // 	  stereo=4;
  // 	  int maxhit=(int)ClusterContainer[stereo][ tmpvec[0][ is[0] ] ].nHit();
  // 	  allsthit+=maxhit; 
  // 	  if(tmpvec[1].size()==0)
  // 	    {
  // 	      stereo=4;
  // 	      maxhit=(int)ClusterContainer[stereo][ tmpvec[0][ is[0] ] ].nHit();
  // 	      for(int i=0;i< maxhit; i++)
  // 		{
  // 		  tmptrack.AddHit(*(ClusterContainer[stereo][tmpvec[0][is[0]]].CDC(cdsMan,i) ) );
  // 		}
  // 	      tmptrackContainer.push_back(tmptrack);		
	     	      
  // 	    }
  // 	  else
  // 	    {
  // 	      for(is[1]=0;is[1]<(int)tmpvec[1].size();is[1]++)
  // 		{
  // 		  for(int n=0;n<2;n++)
  // 		    {		  
  // 		      stereo=4+n;
  // 		      maxhit=(int)ClusterContainer[stereo][ tmpvec[n][ is[n] ] ].nHit();
  // 		      for(int i=0;i< maxhit; i++)
  // 			{
  // 			  tmptrack.AddHit(*(ClusterContainer[stereo][tmpvec[n][is[n]]].CDC(cdsMan,i) ) );
  // 			}
  // 		    }	    	      	      
  // 		  tmptrackContainer.push_back(tmptrack);		
		  
  // 		}
  // 	    }
  // 	}
         
  //   }//track
  // shorttrackContainer=tmptrackContainer;
  // tmptrackContainer.clear();
  return true;  
}
//finding additional stereo for short

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
      MathTools::LineToLine(x1,a1,x2,a2,dl,dis,fitpos1,tmp);
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
      if(!MathTools::LineToHelix(x2,a2,param,fitpos2,fitpos1,dis))
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


bool CDSTrackingMan::CalcVertex(const int &id1,const  int &id2,TrackVertex &vertex)
{
  TrackVertex* Vertextmp=new TrackVertex();

  CDSTrack *track1;
  CDSTrack *track2;
  int tmpid[2]={id1,id2};
  if(id1>id2) {tmpid[0]=id2;tmpid[1]=id1;}
  
  Vertextmp->SetTrackID1(tmpid[0]);
  Vertextmp->SetTrackID2(tmpid[1]);
  
  track1=&trackContainer[tmpid[0]];
  track2=&trackContainer[tmpid[1]];

  TVector3 pos,fitpos1,fitpos2,tmp;
  TVector3 mom1,mom2;
  double param[2][5];
  track1->GetParameters(param[0]);
  track2->GetParameters(param[1]);     

  if(param[0][2]==0 && param[1][2]==0)
    {
     
      TVector3 x1,a1,x2,a2;     
      x1=track1->LineOrigin(); a1=track1->LineDirection();
      x2=track2->LineOrigin(); a2=track2->LineDirection();
      
      double dl=0.2,dis;     
      MathTools::LineToLine(x1,a1,x2,a2,dl,dis,fitpos1,tmp);
      TVector3 tmpd;
      tmpd=(tmp-fitpos1).Unit();
      fitpos2=fitpos1+dis*tmpd;
      
      Vertextmp->SetVertexPos1(fitpos1);
      Vertextmp->SetVertexPos2(fitpos2);

      //      track1->GetMomentum(fitpos1,mom1);
      //      track2->GetMomentum(fitpos2,mom2);
      mom1.SetXYZ(0,0,0);      mom2.SetXYZ(0,0,0);
      Vertextmp->SetVertexPos_mean(0.5*(fitpos1+fitpos2));      
      if(dis<10 ) Vertextmp->SetFlag(true);
      else  Vertextmp->SetFlag(false);
     
    }
  else
    {
      double dis;
      if(MathTools::HelixToHelix(param[0],param[1],fitpos1,fitpos2,dis))
	{
	  Vertextmp->SetVertexPos1(fitpos1);
	  Vertextmp->SetVertexPos2(fitpos2);	  
	  track1->GetMomentum(fitpos1,mom1);
	  track2->GetMomentum(fitpos2,mom2);
	  Vertextmp->SetMomentum1(mom1);
	  Vertextmp->SetMomentum2(mom2);
	  Vertextmp->SetVertexPos_mean(0.5*(fitpos1+fitpos2));	  
	  Vertextmp->SetFlag(true);
	}
      else 
	{  Vertextmp->SetFlag(false);    }
    }  
  Vertex.push_back(*Vertextmp);
  vertex=*Vertextmp;
  delete Vertextmp;
  return true;
}

bool CDSTrackingMan::Execute(CDSHitMan *cdsMan,ConfMan *conf,bool shortflag)
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
  //#if 
#if DEBUG2
  for(int ilay=0;ilay<7;ilay++){
    std::cout<<"######### islay "<<ilay<<std::endl;
    for(int itr=0;itr<ClusterContainer[ilay].size();itr++){
      std::cout<<"------------icluster "<<itr<<" / "<<ClusterContainer[ilay].size()<<std::endl;
      HitCluster *cl=&ClusterContainer[ilay][itr];     
      for(int ihit=0;ihit<cl->nHit();ihit++)
	{	  
	  CDCHit *hit=cl->CDC(cdsMan,ihit);
	  int layer=hit->layer();
	  int hid=hit->hid();
	  int wire=hit->wire();
	  std::cout<<"ihit,layer,hid,wire "<<ihit<<" "<<layer<<" "<<hid<<"  "<<wire<<std::endl;
	}
    }
  }	    
#endif

  FindingTrackCandidate(cdsMan,conf);
#if DEBUG2
  std::cout<<"finding :"<<trackContainer.size()<<std::endl;
#endif
#if 0
  //#if DEBUG2
  for(int itr=0;itr<trackContainer.size();itr++){
    CDSTrack *track1=&trackContainer[itr];
    std::cout<<"######### itrack/ntrack "<<itr<<"/"<<trackContainer.size()<<std::endl;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      std::cout<<"------------layer,nhit "<<layer<<" "<<track1->nTrackHit(layer)<<std::endl;
      for(int ih1=0;ih1<track1->nTrackHit(layer);ih1++){
	CDCHit *hit1=track1->TrackHit(cdsMan,layer,ih1);
	int wire1=hit1->wire();
	std::cout<<"ihit,hid,wire "<<ih1<<" "<<hit1->hid()<<" "<<wire1<<std::endl;
      }
    }
  }	    
#endif

  //####cut too many container####//
  if(trackContainer.size()>600)
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
  if(magneticfield==0) FullAxalLineFit(cdsMan);
  else  FullCircleFit2(cdsMan);
#if DEBUG2
  std::cout<<"circle :"<<trackContainer.size()<<std::endl;
  for(int itr=0;itr<(int)trackContainer.size();itr++)
    {
      CDSTrack *track=&trackContainer[itr];
      std::cout<<track->Chi()<<"  ";
      track->Print();
    }
#endif

  SelectSharedHit(cdsMan);   
#if DEBUG2
  std::cout<<"shared circle :"<<trackContainer.size()<<std::endl;
  for(int itr=0;itr<(int)trackContainer.size();itr++)
    {
      CDSTrack *track=&trackContainer[itr];
      std::cout<<track->Chi()<<"  ";
      track->Print();
    }
#endif

  FindingStereoHit(cdsMan,conf);
#if DEBUG2
  std::cout<<"findingst :"<<trackContainer.size()<<std::endl;
  for(int itr=0;itr<(int)trackContainer.size();itr++)
    {
      CDSTrack *track=&trackContainer[itr];
      std::cout<<track->Chi()<<"  ";
      track->Print();
    }
#endif


#if DEBUG
  double t2=clock();
#endif
  //########################
  //######Stereo Fittig######
  //########################
  if(magneticfield==0) FullStereoLineFit(cdsMan);
  else FullHelixFit(cdsMan);
  
#if DEBUG2
  std::cout<<"helfit :"<<trackContainer.size()<<std::endl;
  for(int itr=0;itr<(int)trackContainer.size();itr++)
    {
      CDSTrack *track=&trackContainer[itr];
      std::cout<<track->Chi()<<"  ";
      track->Print();
    }
#endif


  SelectSharedHit(cdsMan);
#if DEBUG2
  std::cout<<"shared hel :"<<trackContainer.size()<<std::endl;
  for(int itr=0;itr<(int)trackContainer.size();itr++)
    {
      CDSTrack *track=&trackContainer[itr];
      std::cout<<track->Chi()<<"  ";
      track->Print();
    }
#endif


  Calc(cdsMan,conf,false);
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

  if(shortflag)
    {
      DeleteUsedCluster(cdsMan);
      DoShortTracking(cdsMan,conf);
    }
  Calc(cdsMan,conf,false);
#if REMOVECLUSTER
   for(int i=0;i<7;i++)
     ClusterContainer[i].clear();
#endif

  return true;
}
//Execute

bool CDSTrackingMan::DoShortTracking(CDSHitMan *cdsMan,ConfMan *conf)
{
#if TESTSHORT
  //Initialize
  Clear();
  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  selected_chi=CDSParam->GetMaxChi();
  int maxcdchit=CDSParam->GetMaxCDCHit();
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

  //  SearchCluster(cdsMan,conf);
#endif

  FindingShortTrackCandidate(cdsMan,conf);
#if DEBUG2
  std::cout<<"sfinding :"<<shorttrackIDContainer.size()<<std::endl;
#endif
#if 0
  //#if DEBUG2
  for(int itr=0;itr<(int)trackContainer.size();itr++){
    CDSTrack *track1=&trackContainer[itr];
    std::cout<<"######### itrack/ntrack "<<itr<<"/"<<trackContainer.size()<<std::endl;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      std::cout<<"------------layer,nhit "<<layer<<" "<<track1->nTrackHit(layer)<<std::endl;
      for(int ih1=0;ih1<track1->nTrackHit(layer);ih1++){
	CDCHit *hit1=track1->TrackHit(cdsMan,layer,ih1);
	int wire1=hit1->wire();
	std::cout<<"ihit,hid,wire "<<ih1<<" "<<hit1->hid()<<" "<<wire1<<std::endl;
      }
    }
  }	    
#endif

  //####cut too many container####//
  if((int)shorttrackIDContainer.size()>200)
    {
      Clear();
      return true;
    }
  //#####fin initialize###

  for(int itr=0;itr<(int)shorttrackIDContainer.size();itr++){
    CDSTrack *strack=&trackContainer[shorttrackIDContainer[itr]];
    strack->CalcShortInitialParameters();
    if(!strack->FirstCircleFitting(cdsMan)) continue;
    if(!strack->CircleFitting2(cdsMan)) continue;
    strack->Calc(conf);
  }          
#if DEBUG2
  std::cout<<"scfinding :"<<shorttrackIDContainer.size()<<std::endl;
#endif


  SelectSharedHitforShort(cdsMan);  
  //  FindingAdditionalHit(cdsMan,conf);
#if DEBUG2
  std::cout<<"scfinding :"<<shorttrackCIDontainer.size()<<std::endl;
#endif

  for(int itr=0;itr<(int)shorttrackIDContainer.size();itr++){
    CDSTrack *strack=&trackContainer[shorttrackIDContainer[itr]];
    if(!strack->FirstHelixFitting(cdsMan)) continue;
    if( !strack->HelixFitting(cdsMan)) continue;
    strack->Calc(conf);
    strack->SearchHodoHit(cdsMan,conf,CID_IH);
  }
  SelectSharedHitforShort(cdsMan);  
#if DEBUG2
  std::cout<<"sffinding :"<<shorttrackIDContainer.size()<<std::endl;
#endif
#if REMOVECLUSTER
   for(int i=0;i<7;i++)
     ClusterContainer[i].clear();
#endif

  return true;
}

void CDSTrackingMan::FullCircleFit(CDSHitMan *cdsMan)
{
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
      if(!track->FirstCircleFitting(cdsMan) ) 
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
      //      std::cout<<"------------------"<<std::endl;
      for(int trial=0;trial<5;trial++) { 
	if(!track->CircleFitting(cdsMan) ){
	  flag=false;
	  break;
	}
	//	std::cout<<track->Chi()<<std::endl;
      }
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

void CDSTrackingMan::FullCircleFit2(CDSHitMan *cdsMan)
{
  for(int i=0;i<(int)trackContainer.size();i++)
    {
      CDSTrack *track=&trackContainer[i];
      if(!track->CircleFitting2(cdsMan) ) 
	{
#if DEBUG
	  std::cout<<"cir2 track!!"<<std::endl;
#endif
	  track->SetGoodFlag(false);  continue;
	}
#if DEBUG
      std::cout<<" cir "<<track->Chi();
#endif           
    }
}

void CDSTrackingMan::FullAxalLineFit(CDSHitMan *cdsMan)
{
  for(int i=0;i<(int)trackContainer.size();i++)
    {
      CDSTrack *track=&trackContainer[i];
      
      if(!track->FirstLineFitting(cdsMan) ) 
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
      for(int trial=0;trial<10;trial++) { if(!track->LineFitting(cdsMan) ) flag=false;}
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


void CDSTrackingMan::FullStereoLineFit(CDSHitMan *cdsMan)
{
  //#############################
  //######StereoLine Fittig######
  //#############################
  for(int i=0;i<(int)trackContainer.size();i++)
    {

      CDSTrack *track=&trackContainer[i];            
      if(!track->FirstStereoLineFitting(cdsMan) ) 
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
      track->SetHitPos(cdsMan);       
#if NOMINUIT
      int test=0;
      bool flag=true;
      for(int trial=0;trial<1;trial++) { if(!track->StereoLineFitting(cdsMan) ) {flag=false;test=trial;break;}}
      if(!flag) 
	{
#if DEBUG
	  std::cout<<"sl track ! "<<test<<std::endl; track->SetGoodFlag(false); 
#endif
	  continue;
	}
    
      int numt=0;
      double Chi=9999,Chitmp=9999;
      bool hflag1=true;	
      Chitmp=track->Chi();
      while(  numt<20 )
	{
	  Chi=Chitmp;
	  for(int trial=0;trial<10;trial++)
	    { if(!track->StereoLineFitting(cdsMan) ){hflag1=false; break;}  }	
	  Chitmp=track->Chi();
	  if( (Chi-Chitmp)<0.1 ) break;
	  numt++;
	}
      if(!hflag1) { continue;}
#if DEBUG    
      std::cout<<" stline "<<track->Chi()<<std::endl;
#endif           
#else  // NOMINUIT
      track->StereoLineFitting2(cdsMan);	
#endif
    }	

}

void CDSTrackingMan::FullHelixFit(CDSHitMan *cdsMan)
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
      if(!track->FirstHelixFitting(cdsMan) ) 
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
    }
#if DEBUG2
  std::cout<<"helfitfirst :"<<trackContainer.size()<<std::endl;
  for(int itr=0;itr<(int)trackContainer.size();itr++)
    {
      CDSTrack *track=&trackContainer[itr];
      std::cout<<track->Chi()<<"  ";
      track->Print();
    }
#endif

  //  SelectSharedHit(cdsMan);
#if DEBUG2
  std::cout<<"shared helfirst :"<<trackContainer.size()<<std::endl;
  for(int itr=0;itr<(int)trackContainer.size();itr++)
    {
      CDSTrack *track=&trackContainer[itr];
      std::cout<<track->Chi()<<"  ";
      track->Print();
    }
#endif

  for(int i=0;i<(int)trackContainer.size();i++)
    {
      CDSTrack *track=&trackContainer[i];            
      if(!track->HelixFitting(cdsMan) ) continue;	
      
#if DEBUG    
      double ht3=clock();
      std::cout<<" helix "<<track->Chi()<<std::endl;      
      std::cout<<"     helix track time : "<<i<<std::endl;
      std::cout<<"        first helix fit: "<<ht1-ht0<<std::endl;
      std::cout<<"        while helix fit: "<<ht3-ht1<<std::endl;
#endif
    }
  
}

void CDSTrackingMan::Calc(CDSHitMan *cdsMan,ConfMan *conf,bool VERTEX)
{
  Vertex_beam.clear();
  Vertex.clear();
  goodtrackIDContainer.clear();
  shorttrackIDContainer.clear();

  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  selected_chi=CDSParam->GetMaxChi();
  //  int maxcdchit=CDSParam->GetMaxCDCHit();
  //  int cdctdc_u=CDSParam->GetUpperLimitCDCTDC();
  //  int cdctdc_l=CDSParam->GetLowerLimitCDCTDC();
  //  double magneticfield=CDSParam->GetMagneticField();

  //  std::cout<<"selected_chi = "<<selected_chi<<std::endl;
  for(int i=0;i<(int)trackContainer.size();i++)
    {
      CDSTrack *track=&trackContainer[i];
      track->SetHitPos(cdsMan);       
      if(track->IsShort()){
	track->SearchHodoHit(cdsMan,conf);
	track->Calc(conf);
	shorttrackIDContainer.push_back(i);
	continue;
      }
      track->SetGoodFlag(false);
      //      std::cout<<i<<" / "<<trackContainer.size()<<std::endl;
	
      if(track->FittingLevel()==4&&track->Chi()<selected_chi) 
	{
	  track->SetGoodFlag(true);
	  goodtrackIDContainer.push_back(i);
	  track->Calc(conf);
	  track->SearchHodoHit(cdsMan,conf);
	  //	  track->Print();
	}
    }
  //  std::cout<<goodtrackIDContainer.size()<<std::endl;
  if(VERTEX)
    {
      for(int i=0;i<(int)goodtrackIDContainer.size();i++)
	{
	  for(int ii=i+1;ii<(int)goodtrackIDContainer.size();ii++)
	    {
	      TrackVertex vertextmp;
	      CalcVertex(goodtrackIDContainer[i],goodtrackIDContainer[ii],vertextmp);
	    }
	}
    }
  return;
}


bool CDSTrackingMan::GetVertex(const int &trackid1,const int &trackid2,TrackVertex &vertex)
{
  int gid[2]={-1};
  int vertexidtmp=0;
  int vertexid=0;
  
  for(int i=0;i<nGoodTrack();i++)
    {
      if(trackid1==goodtrackIDContainer[i]) gid[0]=i;
      if(trackid2==goodtrackIDContainer[i]) gid[1]=i;
    }
  if(gid[0]==-1 || gid[1]==-1 || gid[0]==gid[1]) return false;
  if(gid[0]>gid[1]) 
    {
      int tmp=gid[0];
      gid[0]=gid[1];
      gid[1]=tmp;
    }

  for(int id1=0;id1<nGoodTrack();id1++)
    { 
      for(int id2=id1+1;id2<nGoodTrack();id2++)
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

  int nbeam=(int)Vertex_beam.size()/(int)goodtrackIDContainer.size();
  int gid[2]={-1,-1};
  if(nbeam<beamid) return false;
  int vertexidtmp=0;
  int vertexid=0;
  for(int bi=0;bi<nbeam;bi++)
    {
      for(int i=0;i<(int)goodtrackIDContainer.size();i++)
	{
	  if(trackid==goodtrackIDContainer[i] && beamid==bi) {gid[0]=i;gid[1]=bi;vertexid=vertexidtmp;}
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
  Clear();
}

void TrackVertex::Clear()
{
  trackid[0]=-999;
  trackid[1]=-999;
  Vertex_mean.SetXYZ(-999,-999,-999);
  Vertex1.SetXYZ(-999,-999,-999);
  Vertex2.SetXYZ(-999,-999,-999);
  Momentum1.SetXYZ(-999,-999,-999);
  Momentum2.SetXYZ(-999,-999,-999);
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
