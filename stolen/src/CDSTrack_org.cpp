#include "CDSTrack.h"

ClassImp(CDSTrack);
#define DEBUG 0
//const double ResolutionOfCDC=0.02;
const double ResolutionOfCDC=0.018;

CDSTrack::CDSTrack() : TObject()
{
  Clear();
}

void CDSTrack::Clear()
{
  for(int i=0;i<15;i++) { hitidContainer[i].clear(); }
  CDHhitID.clear();
  IHhitID.clear();
  goodflag=false;
  CDHflag=false;
  IHflag=false;
  fittinglevel=0;
  for(int i=0;i<5;i++) param[i]=0;
  Line_o.SetXYZ(0,0,0);  Line_d.SetXYZ(0,0,0);
  pid=CDS_DEFAULT;
}

void CDSTrack::GetParameters(double *aparam)
{
  for(int i=0;i<5;i++) aparam[i]=param[i];
}
bool CDSTrack::GetGParameters(double *aparam)
{
  TVector3 tmp;
  return GetParameters(CID_CDC,aparam,tmp);
}

bool CDSTrack::GetParameters(const int &id, double *aparam, TVector3 &vtx)
{
  vtx.SetXYZ(999,999,999);
  if(parCont.find(id)!=parCont.end()){
    for(int i=0;i<5;i++) aparam[i]=parCont[id][i];
    if(vtxContainer.find(id)!=vtxContainer.end()){
      vtx=vtxContainer[id];
      return true;
    }
  }
  return false;
}

#if TMPPARAM
void CDSTrack::GetTmpParameters(const int &num,double *aparam)
{
  if(0<=num&&num<5)
    {
      for(int i=0;i<5;i++) aparam[i]=paramtmp[num][i];
    }
}
#endif
void CDSTrack::SetParameters(const double *aparam)
{
  for(int i=0;i<5;i++) param[i]=aparam[i];
}

void CDSTrack::AddParameters(const int &id,const double *aparam, const TVector3 &vtx)
{
  parContainer tmp;
  for(int i=0;i<5;i++) tmp.push_back(aparam[i]);
  parCont[id]=tmp;
  vtxContainer[id]=vtx;
#if DEBUG
  //#if 1
  std::cout<<"--------- ID = "<<id<<"  "<<
    DetectorList::GetInstance()->GetName(id)<<std::endl;
  for(int i=0;i<5;i++)
    std::cout<<parCont[id][i]<<",";
  std::cout<<std::endl;
  std::cout<<"TVector3("<<vtxContainer[id].X()<<","
	   <<vtxContainer[id].Y()<<","
	   <<vtxContainer[id].Z()<<")"<<std::endl;
#endif
}

double CDSTrack::Momentum(const int &id){
  double tmp[5];
  TVector3 tmpvtx;
  if(GetParameters(id,tmp,tmpvtx)){
    double tmppt=magneticfield*Const/100./tmp[2];
    return tmppt*sqrt(1+tmp[4]*tmp[4]);
  }
  return -9999.;
}

bool CDSTrack::GetNthParameters(int n, int &id, double *aparam, TVector3 &vtx)
{
  std::map<int,parContainer>::iterator it = parCont.begin();
  for( int j=0; j<n; j++ ) ++it;
  id=it->first;
  for(int i=0;i<5;i++) aparam[i]=(it->second)[i];
  vtx=vtxContainer[id];

  return true;
}

bool CDSTrack::GetVertex(const TVector3 &pos, const TVector3 &dir, TVector3 &lpos, TVector3 &hpos)
{
  TVector3 tmpvtx;
  TVector3 tmp1,tmp2;
  double tmp[5];
  int id;
  double dis, tmpdis=999;
  //  std::cout<<"---------------"<<std::endl;
  for(int i=0;i<nParamSets();i++){
    GetNthParameters(i,id,tmp,tmpvtx);
    if(!MathTools::LineToHelix(pos,dir,tmp,tmp1,tmp2,dis) ){
#if DEBUG
      std::cout<<"MathTools::LineToHelix failed in CDSTrack::GetVertex"<<std::endl;
#endif
      continue;
    }
#if DEBUG
    std::cout<<"=== id,getid,dis "<<id<<"  "<<GeomTools::GetID(tmp2)<<"  "<<dis<<std::endl;
    tmpvtx.Print();
    tmp1.Print();
    tmp2.Print();
#endif
    if(dis<tmpdis||GeomTools::GetID(tmp2)==id){
      lpos=tmp1;
      hpos=tmp2;
      tmpdis=dis;
      if(GeomTools::GetID(tmp2)==id) return true;
    }
  }
  if(tmpdis<10) return true;

  //  if(lpos.Mag()<0.01)  Print();
  return false;
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
  if(param[2]!=0)
    {
      double phi=MathTools::CalcHelixPhi(x,y,param);
      z=param[3]-1./param[2]*param[4]*phi;
    }
  else
    {
      double tmp=(x-Line_o.X())/Line_d.X() ;
      z=tmp*Line_d.Z()+Line_o.Z();
    }
}


void CDSTrack::AddHit(const CDCHit &hit)
{
  int cid = hit.cid();
  int layer = hit.layer();
  short hid = hit.hid();
#if 0
  std::cout << " in Add Hit, layer:" << hit.layer() << " wire:" << hit.wire() <\< " CID:" << cid << std::endl;
#endif

  if(cid ==CID_CDC)
    {
      hitidContainer[layer-1].push_back(hid);
#if 0
      std::cout << "?? wire:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].wire()
		<< " hid:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].hid()
		<< " slayer:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].slayer()
		<< std::endl;
#endif
    }
}

void CDSTrack::SetHodoHit(const HodoscopeLikeHit &hit, const TVector3 &vertex)
{
  int cid = hit.cid();
  int hid = hit.hid();
  //  std::cout<<cid<<" "<<hid<<"  "<<hit.seg()<<std::endl;
  if(cid ==CID_IH)
    {
      IHhitID.push_back(hid);
      if(hit.CheckRange())      IHflag=true;
      IHvertex=vertex;
    }
  else if(cid ==CID_CDH)
    {
      CDHhitID.push_back(hid);
      if(hit.CheckRange()) CDHflag=true;
      CDHvertex=vertex;
    }
}

bool CDSTrack::SearchHodoHit(CDSHitMan *cdsMan, ConfMan *conf,const double &cdhmaxdist, const double &ihmaxdist)
{
  // global position correction temporary implemented for only helix 
  // will be modified also for straight track
  // global on/off will be selectable
  // anyway CDH hit search algorithm should be totally modified.
  // -- allow upto neibouring 2 hit 
  // -- hit should be judged by the distance between track at CDH R and CDH hit position considering ctsub
  // -- essentially same algorithm should be also applied to IH hit ( without Z information )

  // clear CDH info before hit search 20130619 hashi
  //  std::cout<<"CDSTrack::SearchCDHHit() start !!!"<<std::endl;

  if(!GoodFlag())    return false;
  SearchHodoHit(cdsMan,conf,CID_IH,ihmaxdist);
  if(!IsShort())
    SearchHodoHit(cdsMan,conf,CID_CDH,cdhmaxdist);
  return true;
}
bool CDSTrack::SearchHodoHit(CDSHitMan *cdsMan, ConfMan *conf,const gCounterID &cid,double maxdphi, bool HITONLY)
{
  double rin=0,rout=0;
  if(cid==CID_CDH){
    CDHvertex.SetXYZ(-999,-999,-999);
    CDHflag=false;
    CDHhitID.clear();
    rin=54.4;
    rout=57.4;
    if(maxdphi==0) maxdphi=7;
  }else if(cid==CID_IH){
    IHvertex.SetXYZ(-999,-999,-999);
    IHflag=false;
    IHhitID.clear();
    rin=13.7;
    rout=14.1;
    if(maxdphi==0) maxdphi=10;
  }else return false;
  double tmpparam[5];
  TVector3 tmp;
  GetParameters(CID_CDC,tmpparam,tmp);   
  TVector3 inpos =MathTools::CalcHelixPosatR(tmpparam,rin);
  TVector3 outpos=MathTools::CalcHelixPosatR(tmpparam,rout);
  if(cid==CID_CDH)  CDHvertex=inpos;

  if(cdsMan->nHodo(cid)==0) return false;
  for( int ii=0; ii<cdsMan->nHodo(cid); ii++ )
    {
      HodoscopeLikeHit *hit=cdsMan->Hodoi(cid,ii);
      if( HITONLY && !hit->CheckRange() ) continue;
      int seg=cdsMan->Hodoi(cid,ii)->seg();
      TVector3 cdhp;
      if(conf->GetGeomMapManager())
	conf->GetGeomMapManager()->GetGPos(cid,seg,cdhp);
      else continue;
      double indphi=inpos.Phi()-cdhp.Phi();
      if(indphi>TMath::Pi()) indphi=TMath::Pi()*2-indphi;
      indphi*=TMath::RadToDeg();
      
      double outdphi=outpos.Phi()-cdhp.Phi();
      if(outdphi>TMath::Pi()) outdphi=TMath::Pi()*2-outdphi;
      outdphi*=TMath::RadToDeg();
      if(TMath::Abs(indphi)>maxdphi&&TMath::Abs(outdphi)>maxdphi) continue;
      SetHodoHit(*hit,inpos);
    }
  return true;
}

bool CDSTrack::SearchHodoHit2(CDSHitMan *cdsMan, ConfMan *conf)
{
  // global position correction is not installed
  // able to calc both strignt and helix  
  CDHvertex.SetXYZ(-999,-999,-999);
  CDHflag=false;
  CDHhitID.clear();
  IHvertex.SetXYZ(-999,-999,-999);
  IHflag=false;
  IHhitID.clear();

  if(!GoodFlag())    return false;

  int cid[2]={CID_CDH,CID_IH};
  //  double rin[2]={54.5,13.7};
  //  double rout[2]={57.5,14.1};
  double maxdist[2]={6.,3.};
  // 20130619 hashi
  double gparam[4];
  conf->GetCDCWireMapManager()->gparam(gparam);

  for(int ic=0; ic<2;ic++){
    if(cdsMan->nHodo(cid[ic])==0) continue;
    TVector3 inpos,outpos;
    //    if(!MathTools::RadiusToHelix(rin[ic],param,inpos) ) continue;
    //    if(!MathTools::RadiusToHelix(rout[ic],param,outpos) ) continue;
    for( int ii=0; ii<cdsMan->nHodo(cid[ic]); ii++ )
      {
	if(cid[ic]==CID_CDH && !cdsMan->Hodoi(cid[ic],ii)->CheckRange() ) continue;
	if(cid[ic]==CID_IH && !(cdsMan->Hodoi(cid[ic],ii)->tdcu()>0) ) continue;
	int seg=cdsMan->Hodoi(cid[ic],ii)->seg();
	TVector3 cdhp;
	if(conf->GetGeomMapManager())
	  conf->GetGeomMapManager()
	    ->GetPos(cid[ic],seg,cdhp);
	double indist2=pow(inpos.X()-cdhp.X(),2)+pow(inpos.Y()-cdhp.Y(),2);
	double outdist2=pow(outpos.X()-cdhp.X(),2)+pow(outpos.Y()-cdhp.Y(),2);
	if(indist2>maxdist[ic]*maxdist[ic]&&outdist2>maxdist[ic]*maxdist[ic]) continue;
	if(cid[ic]==CID_CDH){
	  CDHhitID.push_back(ii);
	  CDHflag=true;
	  CDHvertex=inpos;
	}
	if(cid[ic]==CID_IH){
	  IHhitID.push_back(ii);
	  IHflag=true;
	  IHvertex=inpos;
	}
      }
  }
  return true;
}


bool CDSTrack::RemoveAllHitInLayer( const int &layer ) 
{
  if(!(1<=layer && layer<=15)) return false;
  for(int i=0;i < (int)hitidContainer[layer-1].size();i++) hitidContainer[layer-1].pop_back();
  return true;
}

bool CDSTrack::DeleteHit( const int &layer ,const int &i) 
{
  if(!(1<=layer && layer<=15)) return false;

  TrackHitIDContainer::iterator it = hitidContainer[layer-1].begin();
  for( int j=0; j<i; j++ ) ++it;
  hitidContainer[layer-1].erase( it ) ;
  return true;
}


bool CDSTrack::FirstCircleFitting(CDSHitMan *cdsMan)
{    
    double x[40],y[40],weight[40];
    int numallhit=0;
    for(int numlayer=0;numlayer<7;numlayer++)
      {
	int layer=0;
	if(numlayer<3) layer=numlayer+1;
	else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
	else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
	for(int n=0;n<(int)nTrackHit(layer);n++)
	  {
	    double dl=hit(cdsMan,layer,n)->dl();
	    if(dl<ResolutionOfCDC) dl=ResolutionOfCDC;

	    x[numallhit]=hit(cdsMan,layer,n)->wx();
	    y[numallhit]=hit(cdsMan,layer,n)->wy();
	    weight[numallhit]=dl;
	    numallhit++;
	  }
      }
          
    CircleFit *cirfit=new CircleFit(param,x,y,weight,numallhit);
    cirfit->GetParameters(param);
    cirfit->CalcChi2();
    ChiSqr= cirfit->chisquare();
    dof=cirfit->dof();
    ChiSqr=(double)ChiSqr/dof;
#if TMPPARAM
    cirfit->GetParameters(paramtmp[1]);
    ChiSqrtmp[0]=ChiSqr;
#endif
    delete cirfit;
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
    CirCenterX=(param[0]+1./param[2])*cos(param[1]);
    CirCenterY=(param[0]+1./param[2])*sin(param[1]);
    CirRho=fabs(1./param[2]);
    fittinglevel=1;
    return true;
}



bool CDSTrack::FirstHelixFitting(CDSHitMan *cdsMan)
{
  double disZ[40];
  int numstlayer=0;
  // set nearest point on wire to the Circle for each stereo layer
  for(int stlayer=0 ;stlayer<8;stlayer++)
    {
      int layer=0;
      if(0<=stlayer && stlayer<=3) layer=stlayer+4;
      else if(4<=stlayer && stlayer<=7) layer=stlayer+6;

      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  CDCHit *cdc=hit(cdsMan,layer,n);
	  TVector3 wpos,wposp;	  
	  wpos=hit(cdsMan,layer,n)->wpos();
	  wposp=hit(cdsMan,layer,n)->wposp();
	  double aw,bw,cw;	 
	    
	  if(TMath::Abs( wposp.x()-wpos.x() )>0.00001)
	    {
	      if(TMath::Abs( wposp.x()-wpos.x() )<0.0001){
		wpos.Print();
		wposp.Print();
	      }
	      aw=(wposp.y()-wpos.y() )/(wposp.x()-wpos.x() );
	      bw=-1;
	      cw= -(aw*wpos.x()+bw*wpos.y() );
	    }
	  else
	    {
	      // check !!
	      bw=(wposp.x()-wpos.x() )/(wposp.y()-wpos.y() );
	      aw= -1;
	      cw= -(aw*wpos.x()+bw*wpos.y() ); 
	    }
	  
	  double x_p,x_n,y_p,y_n;
	  TVector3 WPos;

	  if(!(MathTools::LineToCircle(aw,bw,cw,CirRho,CirCenterX,CirCenterY,x_p,y_p,x_n,y_n)) ){
	    // std::cout<<"!!!  "<<std::endl;
	    // wposp.Print();
	    // wpos.Print();
	    continue;
	  }
	  TVector3 Pos_p(x_p,y_p,0);
	  TVector3 Pos_n(x_n,y_n,0);
	  if( (Pos_p-wpos).Perp()< (Pos_n-wpos).Perp() ){ WPos.SetX(x_p);WPos.SetY(y_p);}
	  else{ WPos.SetX(x_n);WPos.SetY(y_n);}
	  
	  double dis=(wpos-WPos).Perp();
	  double wdis=(wpos-wposp).Perp();
	  WPos.SetZ( wpos.z()+(wposp.z()-wpos.z() )*dis/wdis );
	  hit(cdsMan,layer,n)->SetHitPosition(WPos);    
	  
	  //########dis Z ###########//
	  double dl=cdc->dl();
	  double tilt=cdc->tilt()*TMath::DegToRad();
	  disZ[numstlayer]=fabs(dl/sin(tilt) );
	  numstlayer++;
	}
    }
  
  // search minimum chi2 combination of the sign of disZ 
  // in phi-Z plane tracking
  // z = - 1/kappa * dipangle * phi + dz
  // kappa = param[2] in CircleFit
  // dipangle = param[4]
  // dz = param[3]
#define HELLR 0
  int checki=-1;	  
#if HELLR
  for(int i =0;i<pow(2,numstlayer);i++)
    {
#else
      int i=0;
#endif
      double hitx[256],hity[256],hitz[256];
      double phi[256];
      //      double weight[256];
      int allhit=0;
      int numst=0;
      double s_x=0,s_y=0,s_xx=0,s_xy=0;
      
      for(int layer=1;layer<=15;layer++)
	{
	  for(int n=0;n<(int)nTrackHit(layer);n++)
	    {
	      CDCHit *cdc=hit(cdsMan,layer,n);	      
	      if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
		{
		  hitx[allhit]= cdc->x();
		  hity[allhit]= cdc->y();
		  hitz[allhit]= cdc->z();
		  phi[allhit]=MathTools::CalcHelixPhi(hitx[allhit],hity[allhit],param);
		  phi[allhit]=phi[allhit]/param[2];
#if HELLR
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
#endif
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
      // aa*phi+bb*z+cc=0
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
#if HELLR      
    }

  // best combination of disZ signs for phi-Z track determined
  
  // set hit Z position considering the sign of drift length.
  int numst=0;
  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
	    {
	      CDCHit *cdc=hit(cdsMan,layer,n);
	      double hitx= cdc->x();
	      double hity= cdc->y();
	      double hitz= cdc->z();
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
	      cdc->SetHitPosition(hitx,hity,hitz);   
	      numst++;
	    } // stereo layer
	}//ntrackhit
    }//layer
#endif
#if TMPPARAM
  ChiSqrtmp[2]=ChiSqr;
  for(int n=0;n<5;n++) paramtmp[3][n]=param[n];
#endif
  if( ChiSqr >999999 ) 
    {
#if DEBUG
      std::cout<<"Too big chi2/dof fh"<<std::endl;
#endif
      return false;
    }
  fittinglevel=3;
  SetHitPos(cdsMan);       
  return true;
}

bool CDSTrack::CircleFitting(CDSHitMan *cdsMan)
{
    
    double x[40],y[40],weight[40];
    int numallhit=0;
    
    for(int numlayer=0;numlayer<7;numlayer++)
      {
	int layer=0;
	if(numlayer<3) layer=numlayer+1;
	else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
	else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
	for(int n=0;n<(int)nTrackHit(layer);n++)
	  {
	    CDCHit *cdc=hit(cdsMan,layer,n);
	    double dl=cdc->dl();
	    if(dl<0) dl=0.0001;

	    x[numallhit]=cdc->wx();
	    y[numallhit]=cdc->wy();

	    double rtmp=sqrt( (x[numallhit]-CirCenterX)* (x[numallhit]-CirCenterX)+(y[numallhit]-CirCenterY)*(y[numallhit]-CirCenterY) );
	    double sintmp= (y[numallhit]-CirCenterY)/rtmp;
	    double costmp= (x[numallhit]-CirCenterX)/rtmp;
	    double sign;
	    if(rtmp<fabs(1./param[2])) sign=1;
	    else  sign=-1;
	    
	    x[numallhit]+= sign*dl*costmp;
	    y[numallhit]+= sign*dl*sintmp;
	    cdc->SetHitPosX(x[numallhit]);
	    cdc->SetHitPosY(y[numallhit]);
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

    if(!cirfit->stat()) {
      delete cirfit;
      return false;
    }
    delete cirfit;
    ChiSqr=(double)ChiSqr/dof;
#if TMPPARAM
    ChiSqrtmp[1]=ChiSqr;
     for(int n=0;n<5;n++) paramtmp[2][n]=param[n];
#endif
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

    CheckCharge();

#if TMPPARAM
    for(int n=0;n<5;n++) paramtmp[2][n]=param[n];
#endif
    fittinglevel=2;
    SetHitPos(cdsMan);
    return true;
}
//******circleFitting*********//

void CDSTrack::CheckCharge()
{
  // new!! on 20130620
  //################################//
  //###########Check Charge ##########//	
  //################################//
  double phi_in=MathTools::CalcDeg(Clusterx[0]-CirCenterX,Clustery[0]-CirCenterY);
  double phi_mid=MathTools::CalcDeg(Clusterx[3]-CirCenterX,Clustery[3]-CirCenterY);
  double phi_out=MathTools::CalcDeg(Clusterx[6]-CirCenterX,Clustery[6]-CirCenterY);
#if 0
  for(int i=0;i<7;i++)
    std::cout<<i<<"  "<<Clusterx[i]<<"  "<<Clustery[i]<<std::endl;
  std::cout<<"------"<<phi_in<<"  "<<phi_mid<<"  "<<phi_out<<std::endl;
#endif
  phi_mid=phi_out;

  if(phi_mid>270 && phi_in<90 ) phi_mid-=360;
  else if(phi_in>270 &&phi_mid<90 ) phi_in-=360;
  
  if( (phi_mid-phi_in)<0 && param[2]<0) {
    //    std::cout<<"+"<<std::endl;
    param[0]=-param[0];
    param[2]=-param[2];
    param[4]=-param[4];
    param[1]+=TMath::Pi();
  }
  else if( phi_mid-phi_in>0 && param[2]>0){
#if 0
    std::cout<<"-"<<std::endl;
    for(int i=0;i<5;i++)
      std::cout<<param[i]<<"  ";
    std::cout<<" before"<<std::endl;
#endif
    param[0]=-param[0];
    param[2]=-param[2];
    param[4]=-param[4];
    param[1]-=TMath::Pi();
#if 0
    for(int i=0;i<5;i++)
      std::cout<<param[i]<<"  ";
    std::cout<<" after"<<std::endl;
#endif
  }
  //  if(param[1]<0||param[1]>TwoPi) std::cout<<"========  abnormal param[1] !!!  "<<param[1]<<std::endl;
  if(param[1]<0) param[1]+=2*TMath::Pi();
  else if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi();
//   if(param[1]<-TMath::Pi()){
//     param[1]+=2*TMath::Pi();
//   }
//   else if(param[1]>TMath::Pi()){
//     param[1]-=2*TMath::Pi();
//   }
}

void CDSTrack::CheckCharge2()
{
  //################################//
  //###########Check Charge ##########//	
  //################################//
  double phi_c=MathTools::CalcDeg(CirCenterX,CirCenterY);
  double phi_mid=theta;
  double phi_diff=phi_mid-phi_c;
  if( phi_diff>180 ) phi_diff-=360;
  if( phi_diff<-180 ) phi_diff+=360;  

  //  if(phi_mid>270 &&phi_c<90 ) phi_mid-=360;
  //  else if(phi_c>270 &&phi_mid<90 ) phi_c-=360;
  
  if( phi_diff<0 && param[2]>0) {
    param[0]=-param[0];
    param[2]=-param[2];
    param[4]=-param[4];
    param[1]+=TMath::Pi();
  }
  else if( phi_diff>0 && param[2]<0){
    param[0]=-param[0];
    param[2]=-param[2];
    param[4]=-param[4];
    param[1]-=TMath::Pi();
  }
  if(param[1]<0) param[1]+=2*TMath::Pi();
  else if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi();
  //  std::cout<<" phi_diff "<<phi_diff<<" "<<phi_mid<<" "<<param[2]<<std::endl;
}

bool CDSTrack::CircleFitting2(CDSHitMan *cdsMan)
{
    
  int allhit=0;

  for(int numlayer=0;numlayer<7;numlayer++)
    {
      int layer=0;
      if(numlayer<3) layer=numlayer+1;
      else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
      else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  allhit++;
	}
    }
  int lrnum=0;
  int firstchi=0;
  double chilr=0;double chilrtmp=0;
  double paramlr[5];
  double x[40],y[40],weight[40];
  int numallhit=0;
  // FirstFitting with MWPC mode
  for(int numlayer=0;numlayer<7;numlayer++)
    {
      int layer=0;
      if(numlayer<3) layer=numlayer+1;
      else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
      else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  CDCHit* cdc=hit(cdsMan,layer,n);
	  x[numallhit]=cdc->wx();
	  y[numallhit]=cdc->wy();
	  weight[numallhit]=0.3;
	  numallhit++;
	}
    }
  if(!MathTools::CircleFit(x,y,weight,numallhit,paramlr,chilrtmp)){
    std::cout<<"circlefitting2: mwpc fitting failed!!!"<<std::endl;
    return false;
  }

  CirCenterX=paramlr[0];
  CirCenterY=paramlr[1];
  CirRho=paramlr[2];
  
  
  //fitting all LR combinations
  for(int ilr=0;ilr<pow(2,allhit);ilr++)
    {
      numallhit=0;	
      for(int numlayer=0;numlayer<7;numlayer++)
	{
	  int layer=0;
	  if(numlayer<3) layer=numlayer+1;
	  else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
	    else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
	  for(int n=0;n<(int)nTrackHit(layer);n++)
	    {
	      CDCHit* cdc=hit(cdsMan,layer,n);
	      double dl=cdc->dl();
	      if(dl<0) dl=0.0001;
	      
	      x[numallhit]=cdc->wx();
	      y[numallhit]=cdc->wy();
	      
	      double rtmp=sqrt( (x[numallhit]-CirCenterX)* (x[numallhit]-CirCenterX)+(y[numallhit]-CirCenterY)*(y[numallhit]-CirCenterY) );
	      double sintmp= (y[numallhit]-CirCenterY)/rtmp;
	      double costmp= (x[numallhit]-CirCenterX)/rtmp;
	      double sign(0.0);
	      int check=ilr;
	      
	      if(check<pow(2,numallhit) ) sign=1;
	      else 
		{
		  int cflag;
		  for(int n=0;n<numallhit;n++)
		    {
		      cflag=check%2;
		      check=(check-cflag)/2;
		    }
		  cflag=check%2;
		  if(cflag==0)  sign=1;
		  else if(cflag==1) sign=-1;
		}
	      x[numallhit]+= sign*dl*costmp;
	      y[numallhit]+= sign*dl*sintmp;
	      weight[numallhit]=ResolutionOfCDC;
	      numallhit++;
	    }
	}
      if(!MathTools::CircleFit(x,y,weight,numallhit,param,chilrtmp))
	{
	  if(ilr==firstchi) firstchi++;
	  continue;
	}
      if(chilrtmp<chilr || ilr==firstchi)
	{
	  lrnum=ilr;
	  chilr=chilrtmp;
	  for(int n=0;n<5;n++) paramlr[n]=param[n];
	}
    }//LR
  //  std::cout<<"-------------------------------------------"<<std::endl;
  //  std::cout<<"CirParam1: "<<paramlr[0]<<" , "<<paramlr[1]<<" , "<<paramlr[2]<<std::endl;
  param[0] = sqrt(paramlr[0]*paramlr[0]+paramlr[1]*paramlr[1]) - paramlr[2] ;
  param[1] = atan2(paramlr[1],paramlr[0]);
  param[2] = 1./paramlr[2];
  
  //tumenaosi    
  numallhit=0;    
  for(int numlayer=0;numlayer<7;numlayer++)
    {
      int layer=0;
      if(numlayer<3) layer=numlayer+1;
  	else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
  	else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
  	for(int n=0;n<(int)nTrackHit(layer);n++)
  	  {
  	    CDCHit *cdc=hit(cdsMan,layer,n);
  	    double dl=cdc->dl();
  	    if(dl<0) dl=0.0001;
  	    x[numallhit]=cdc->wx();
  	    y[numallhit]=cdc->wy();
  
  	    double rtmp=sqrt( (x[numallhit]-CirCenterX)* (x[numallhit]-CirCenterX)+(y[numallhit]-CirCenterY)*(y[numallhit]-CirCenterY) );
  	    double sintmp= (y[numallhit]-CirCenterY)/rtmp;
  	    double costmp= (x[numallhit]-CirCenterX)/rtmp;
  	    double sign(0.0);
  	    int check=lrnum;
  
  	    if(check<pow(2,numallhit) ) sign=1;
  	    else 
  	      {
  		int cflag;
  		for(int n=0;n<numallhit;n++)
  		  {
  		    cflag=check%2;
  		    check=(check-cflag)/2;
  		  }
  		cflag=check%2;
  		if(cflag==0)  sign=1;
  		else if(cflag==1) sign=-1;
  	      }
  	    x[numallhit]+= sign*dl*costmp;
  	    y[numallhit]+= sign*dl*sintmp;
  	    cdc->SetHitPosX(x[numallhit]);
  	    cdc->SetHitPosY(y[numallhit]);
  	    weight[numallhit]=ResolutionOfCDC;
  	    numallhit++;
  	  }
    }
  //
  if(param[1]<0) param[1]+=2*TMath::Pi();
  else if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi();
  ChiSqr=(double)chilr;
#if TMPPARAM
  ChiSqrtmp[1]=ChiSqr;
  for(int n=0;n<5;n++) paramtmp[2][n]=param[n];
#endif
  if( ChiSqr >99999 ) 
    {
#if DEBUG
      std::cout<<"Too big chi2/dof t"<<std::endl;
#endif
      return false;
    }

  CheckCharge();  
  
#if TMPPARAM
  for(int n=0;n<5;n++) paramtmp[2][n]=param[n];
#endif
  fittinglevel=2;
  SetHitPos(cdsMan);

  return true;
}


//*********Circle Fitng2************//
bool CDSTrack::HelixFitting(CDSHitMan *cdsMan)
{
  int allhit=0;
  TVector3 wirepos[50],wiredir[50];
  double weight[50],drift[50];

  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  CDCHit *cdc=hit(cdsMan,layer,n);
	  wirepos[allhit]=(cdc->wpos()+cdc->wposp())*0.5;
	  if( (1<=layer &&layer<=3) || (8<=layer && layer<=9) || (14<=layer && layer<=15) )
            {
              wirepos[allhit].SetZ(-999);
            }
	  wiredir[allhit]=(cdc->wpos()-cdc->wposp());
	  weight[allhit]=ResolutionOfCDC;
	  drift[allhit]=cdc->dl();
	  allhit++;
	}      
    }
  
  
  HelixFit *helixfit=new HelixFit(param,wirepos,wiredir,weight,drift,allhit);
  helixfit->GetParameters(param); 
  CirCenterX=(param[0]+1./param[2])*cos(param[1]);
  CirCenterY=(param[0]+1./param[2])*sin(param[1]);
  CirRho=fabs(1./param[2]);
  ChiSqr= helixfit->chisquare();
  dof= helixfit->dof();
  ChiSqr=(double)ChiSqr/dof;
#if 0
  std::cout<<"helix param : ";
  for(int n=0;n<5;n++) std::cout<<param[n]<<"  ";
  std::cout<<std::endl;
  std::cout<<"        x,y,r  "<<CirCenterX<<"  "<<CirCenterY<<"  "<<CirRho<<std::endl;
#endif
#if TMPPARAM
  for(int n=0;n<5;n++) paramtmp[4][n]=param[n];
  ChiSqrtmp[3]=(double)ChiSqr/dof;
#endif
  // if(!helixfit->stat()) return false;
  delete helixfit;
  if( ChiSqr >9999 ) 
    {
#if DEBUG
      std::cout<<"Too big chi2/dof h"<<std::endl;
#endif
      return false;
    }
  //  SetHitPos(cdsMan,true);
  SetHitPos(cdsMan);
  fittinglevel=4;
  return true;
}

void CDSTrack::SetHitPos(CDSHitMan *cdsMan, bool PRINT)
{    
#if 0
  for(int i=0;i<5;i++)
    std::cout<<param[i]<<"  ";
  std::cout<<std::endl;
#endif

  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  CDCHit *cdc=hit(cdsMan,layer,n);
	  double dl=cdc->dl();
	  if(dl<0) dl=0.002;
	  
	  TVector3 whpos,hpos,wpos,wposp,dline;
	  wpos=cdc->wpos();
	  wposp=cdc->wposp();
	  hpos=cdc->hitpos();
	  dline=wposp-wpos;
	  dline=dline.Unit();
	  
	  double k=hpos*dline-wpos*dline; 
	  whpos=wpos+k*dline; //nearest point to hpos on the line
	  if(hpos.Mag()>999)
	    whpos=(wpos+wposp)*0.5;
	  TVector3 lnest,hnest;
	  TVector3 hitpos;
	  double dis,resi;
	  if(param[2]==0)
	    {		
	      MathTools::LineToLine(Line_o,Line_d,wpos,dline,dl,dis,lnest,hitpos);
	    }
	  else 
	    {
	      if(!MathTools::LineToHelix(whpos,dline,param,lnest,hnest,dis) )
		{ 		  
		  cdc->SetHitPosition(-999,-999,-999);	    
		  cdc->SetResolution(999.);
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
	  // if(PRINT){
	  //   std::cout<<layer<<"  dl = "<<dis<<std::endl;
	  // }

	  resi=(dl-dis);

	  double phih=hitpos.Phi();
	  double phiw=whpos.Phi();
	  double dphi=phih-phiw;
	  if(dphi<TMath::Pi()) dphi+=2*TMath::Pi();
	  if(dphi>TMath::Pi()) dphi-=2*TMath::Pi();
	  if(dphi<0) cdc->SetLeftRight(-1);
	  else cdc->SetLeftRight(1);

	  if(PRINT){
	    std::cout<<"layer,i,dl,dist,resi  "<<layer<<"  "<<n<<"  "<<dl<<"  "<<dis<<"  "<<resi<<"  ";
	    hitpos.Print();
	  }
	  cdc->SetResolution(resi);
	  cdc->SetHitPosition(hitpos);	    
	}
    }
}

bool CDSTrack::GetMomentum(const TVector3 &orgpos,TVector3 &p, bool GLOBAL, bool ELOSS)//, ConfMan *conf)
{
  TVector3 pos=orgpos;
  TVector3 initpos;
  double tmppar[5];
  int volid;
  if(ELOSS){
    volid=GeomTools::GetID(pos);
    if(!GetParameters(volid,tmppar,initpos)){
#if 1
      std::cout<<"GetMomentum::volid  "<<volid<<"  "; pos.Print();
#endif
      return false;
    }
      //      if(!GetParameters(CID_CDC,tmppar,initpos)) return false;
  }
  else if(GLOBAL){
    if(!GetParameters(CID_CDC,tmppar,initpos)) return false;
  }
  else GetParameters(tmppar);
 
#if DEBUG
  std::cout<<"++++++ CDSTrack::GetMomentum ++++++"<<std::endl;
  std::cout<<"vertex pos "; orgpos.Print();
  std::cout<<"init id,pos  "<< volid <<"  "; initpos.Print(); 
#endif
  // if(sqrt( (pos.x()-hpos.x())*(pos.x()-hpos.x())+(pos.y()-hpos.y())*(pos.y()-hpos.y()) )>5 ) 
  //   {
  //     std::cout<<"$E GetMomentum() check input pos "<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<" "
  // 	       <<" hpos "<<hpos.x()<<" "<<hpos.y()<<" "<<hpos.z()<<" phi "<<phi<<std::endl;
  //     p.SetXYZ(-999,-999,-999);
  //     return false;
  //   }
  // int trial=0;
  // while(fabs(pos.z()-hpos.z() )>3.0)
  //   {
  //     trial++;
  //     double tmpphi=phi-2*trial*TMath::Pi();
  //     hpos=MathTools::GetPosition(tmpphi,tmppar);
  //     //      std::cout<<trial<<"\t"<<tmpphi<<"\t";
  //     if( fabs(pos.z()-hpos.z() )>3.0 )
  // 	{
  // 	  tmpphi=phi+2*trial*TMath::Pi();
  // 	  hpos=MathTools::GetPosition(tmpphi,tmppar);
  // 	}
  //     if(trial>5 && fabs(pos.z()-hpos.z() )>3.0)
  // 	{
  // 	  std::cout<<"$E3 GetMomentum() pos "<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<" "
  // 		   <<" hpos "<<hpos.x()<<" "<<hpos.y()<<" "<<hpos.z()<<" phi "<<phi<<std::endl;
  // 	  for(int i=0;i<5;i++)
  // 	    std::cout<<tmppar[i]<<"\t";
  // 	  std::cout<<std::endl;
  // 	  p.SetXYZ(-999,-999,-999);
  // 	  return false;
  // 	}
  //   }
  double momout;
  if(ELOSS){
    double momtmp = fabs(1/tmppar[2]*(Const*magneticfield)/100.)*sqrt(1+tmppar[4]*tmppar[4]);
    //    std::cout<<"eloss calc tgeo"<<std::endl;
    double tmptof;
    if(!ELossTools::CalcHelixElossPointToPointTGeo(tmppar,initpos,pos,momtmp,Mass(),momout,tmptof,0.)){
#if DEBUG
      std::cout<<" error in ELossTools::CalcHelixElossPointToPointTGeo()"<<std::endl;
#endif
      return false;
    }
    //    std::cout<<"eloss calc tgeo OK"<<std::endl;
    double newpar[5];
    MathTools::ChangePivot(TVector3(0,0,0),pos,tmppar,newpar,charge());
    newpar[2]=magneticfield*Const/100./(momout/sqrt(1+newpar[4]*newpar[4]));
    MathTools::ChangePivot(pos,TVector3(0,0,0),newpar,tmppar,charge());
  }
  double phi=MathTools::CalcHelixPhi(pos.x(),pos.y(),tmppar);
  TVector3 hpos=MathTools::GetPosition(phi,tmppar);

  double pttmp = fabs(1/tmppar[2]*(Const*magneticfield)/100.);
  double px = pttmp*(-1*sin(tmppar[1]+phi));
  double py = pttmp*(cos(tmppar[1]+phi));
  double pz = pttmp*(tmppar[4]);
  p.SetXYZ(px,py,pz);
  return true;

}

void CDSTrack::ReconstructHit(CDSHitMan *cdsMan,ConfMan *Conf)
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
	  if(param[2]<0) phidir=param[1]+TwoPi/4.;
	  else if(param[2]>0) phidir=param[1]-TwoPi/4.;

	  if(phidir<0) phidir+=TwoPi;
	  else if(phidir>TwoPi) phidir-=TwoPi;
	  double def_deg=fabs( hit.phi()*Deg2Rad - phidir );
	  if(def_deg>TwoPi-1.0) {def_deg-=TwoPi;def_deg=fabs(def_deg);}
	  if(def_deg>TwoPi/4.)  continue; 
	    

	  if( (1<=layer &&layer<=3) || (8<=layer &&layer<=9) || (14<=layer &&layer<=15) )	 
	    {
	      TVector3 wpos;
	      wpos.SetXYZ(hit.wx(),hit.wy(),hit.wz() );
	      double distmp,xest,yest;
	      MathTools::PointToCircle(wpos.x(),wpos.y(),CirRho,CirCenterX,CirCenterY,distmp,xest,yest);
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
	      MathTools::PointToCircle(wpos.x(),wpos.y(),CirRho,CirCenterX,CirCenterY,distmp,xest,yest);
	      if(distmp>10) continue;
	      dline=wposp-wpos;
	      dline=dline.Unit();

	      TVector3 lnest,hnest;
	      //	      double distmp;
	      if(!MathTools::LineToHelix(wpos,dline,param,lnest,hnest,distmp) )
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
      if(nTrackHit(layer)==0) continue;
      if(hit(cdsMan,layer,0)->slayer()==1 ||hit(cdsMan,layer,0)->slayer()==7)
	{
	  in_outpos[in_outnum]=hit(cdsMan,layer,0)->wpos();
	  in_outnum++;
	}
      
      if(hit(cdsMan,layer,0)->slayer()==4 )
	{	  
	  midpos[midnum]= hit(cdsMan,layer,0)->wpos();
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
  param[1]=MathTools::CalcDeg(cos_c,sin_c);
  param[1]=param[1]*TMath::DegToRad();
  
  if(jitter<0 ) 
    {param[2]=-1./rho;param[0]=-1*param[0];param[1]=param[1];}
  else 
    {param[2]=1./rho;param[0]=param[0];param[1]-=TMath::Pi();}
  param[3]=0;
  param[4]=0;
  CirCenterX=(param[0]+1./param[2])*cos(param[1]);
  CirCenterY=(param[0]+1./param[2])*sin(param[1]);
  CirRho=fabs(1./param[2]);
  
#if TMPPARAM
  for(int n=0;n<5;n++) paramtmp[0][n]=param[n];
#endif
  return true;
}
//CalcInitial

bool CDSTrack::CalcShortInitialParameters()
{
  
  TVector3 dir[2];
  TVector3 pos[2];
  TVector3 spos[3];

  spos[0].SetXYZ( Clusterx[0],
		  Clustery[0],
		  0);
  spos[2].SetXYZ( Clusterx[3],
		  Clustery[3],
		  0);
  spos[1].SetXYZ( (Clusterx[1]+Clusterx[2])/2.0,
		 (Clustery[1]+Clustery[2])/2.0,
		 0);

  pos[0]=(spos[0]+spos[1])*0.5;
  pos[1]=(spos[1]+spos[2])*0.5;
  dir[0]=(spos[0]-spos[1]).Orthogonal();
  dir[1]=(spos[1]-spos[2]).Orthogonal();
  dir[0]=dir[0].Unit();
  dir[1]=dir[1].Unit();
  
  double  k=(pos[1]-pos[0]).Cross(dir[1]).Mag()/dir[0].Cross(dir[1]).Mag();
  TVector3 cen;cen=pos[0]+k*dir[0];
  double rho=(cen-spos[0]).Mag();
  double d_c=sqrt(cen.x()*cen.x()+cen.y()*cen.y());
  double cos_c=-cen.x()/d_c;
  double sin_c=-cen.y()/d_c;

  param[0]=d_c-rho;
  param[1]=MathTools::CalcDeg(cos_c,sin_c);
  param[1]=param[1]*TMath::DegToRad();
  if(jitter<0 ) 
    { param[2]=-1./rho; param[0]=-1* param[0]; param[1]= param[1];}
  else 
    { param[2]=1./rho; param[0]= param[0]; param[1]-=TMath::Pi();}
  param[3]=0;
  param[4]=0;
  CirCenterX=(param[0]+1./param[2])*cos(param[1]);
  CirCenterY=(param[0]+1./param[2])*sin(param[1]);
  CirRho=fabs(1./param[2]);
  
#if TMPPARAM
  for(int n=0;n<5;n++) paramtmp[0][n]=param[n];
#endif
  return true;
  

}
//CalcShortInitial
bool CDSTrack::Retiming(CDSHitMan *cdsMan, ConfMan *conf,double beta,bool SLEW)
{
  double cdhtime;
  int seg;
  if(!GetCDHHit(cdsMan,seg,cdhtime)) return false;
  double obeta2=1./pow(beta,2);
  double corr=-3.27644-0.124793*obeta2+0.00119356*obeta2*obeta2+6.32108*exp(-0.295883*obeta2); 
  corr+=0.18691-0.0923064*obeta2+0.00113628*obeta2*obeta2+0.430231*exp(-0.753885*obeta2); 
  //  std::cout<<obeta2<<"  "<<corr<<std::endl;
  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  CDCHit *cdc= hit(cdsMan,layer,n);
	  TVector3 hitpos=cdc->hitpos();
	  double dis=MathTools::CalcHelixArc(param,CDHvertex,hitpos);
	  double timeminus=dis/beta/Const/100;
	  double ctime= cdhtime-timeminus-5; //(cdhtime-timeminus-t0time) + t0time;
	  //	  std::cout<<"chttime "<<ctime<<" "<<cdhtime<<" "<<timeminus<<" "<<CDHvertex.x()<<std::endl;
	  if(SLEW)    ctime+=corr;	
	  cdc->SetTimeOffset(conf,ctime);
	}
    }  
  return true;
}
void CDSTrack::Refit(CDSHitMan *cdsMan,ConfMan *conf)
{
  //Re-Calc CDCHit
  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  CDCHit *cdc= hit(cdsMan,layer,n);
#if DEBUG
	  std::cout<<"hit "
		   <<"lay "<<cdc->layer()
		   <<" wire "<<cdc->wire()
		   <<" dt "<<cdc->dt()
		   <<" dl "<<cdc->dl()
		   <<" hx "<<cdc->x()
		   <<" hy "<<cdc->y()<<std::endl;
#endif
	  cdc->Calc(conf);
#if DEBUG
	  std::cout<<"af "
		   <<"lay "<<cdc->layer()
		   <<" wire "<<cdc->wire()
		   <<" dt "<<cdc->dt()
		   <<" dl "<<cdc->dl()
		   <<" hx "<<cdc->x()
		   <<" hy "<<cdc->y()<<std::endl;
#endif
	}
    }

  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  magneticfield=CDSParam->GetMagneticField();
  int selected_chi=CDSParam->GetMaxChi();

  if(magneticfield==0)
    {
      
      if(!FirstLineFitting(cdsMan) ) 
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
      for(int trial=0;trial<10;trial++) { if(!LineFitting(cdsMan) ) flag=false;}
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


      if(!FirstStereoLineFitting(cdsMan) ) 
	{
#if DEBUG
	  std::cout<<"first sline track!!"<<std::endl; 
#endif
	  SetGoodFlag(false);  return;
	}
      
#if DEBUG
      std::cout<<" firstsl "<<Chi();
#endif

    SetHitPos(cdsMan);       
    int test=0;
    flag=true;
    for(int trial=0;trial<5;trial++) { if(!StereoLineFitting(cdsMan) ) {flag=false;test=trial;break;}}
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
	  { if(!StereoLineFitting(cdsMan) ){hflag1=false; break;}  }	
	Chitmp=CDSTrack::Chi();
	if( (Chi-Chitmp)<0.1 ) break;
	numt++;
      }
    if(!hflag1) { return;}

#if DEBUG    
    std::cout<<" stline "<<Chi<<std::endl;
#endif


    }
  else
    {
      
      if(!CircleFitting2(cdsMan) ) 
	{
#if DEBUG
	  std::cout<<"first cir track!!"<<std::endl;
#endif
	  SetGoodFlag(false);  return;
	}
      
#if DEBUG
      //      std::cout<<" firstcir "<<Chi;
#endif
      

      if(!FirstHelixFitting(cdsMan) ) 
	{
#if DEBUG
	  std::cout<<"first hel track!!"<<std::endl; 
#endif
	  SetGoodFlag(false);  return;
	}
      
#if DEBUG
      //      std::cout<<" firsthel "<<track->Chi();
#endif
	if( Chi()>2*selected_chi ) 
	  {
#if DEBUG
	    std::cout<<"firsthel chi!"<<std::endl; 
#endif 
	    SetGoodFlag(false);  return;
	  }
		   
    SetHitPos(cdsMan);       

    if(!HelixFitting(cdsMan) ) return;	

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
      param[2]=0;
    }
  else
    {
      CheckCharge();
      CirCenterX=(param[0]+1./param[2])*cos(param[1]);
      CirCenterY=(param[0]+1./param[2])*sin(param[1]);
      CirRho=fabs(1./param[2]);
      pt=magneticfield*Const/100./param[2];
      double global[5];
      conf->GetCDCWireMapManager()->HelixLocalToGlobal(param,global,charge());
      double tmpr=15.3+(53.-15.3)/2.;
      TVector3 tmppos=MathTools::CalcHelixPosatR(global,tmpr);
      parCont.clear();
      vtxContainer.clear();
      AddParameters(CID_CDC,global,tmppos);
#if 0
      std::cout<<"++++++++++++++"<<std::endl;
      for(int i=0;i<5;i++) std::cout<<", "<<param[i]; std::cout<<std::endl;
      for(int i=0;i<5;i++) std::cout<<", "<<global[i]; std::cout<<std::endl;
      std::cout<<"++++++++++++++"<<std::endl;
#endif
      //      CalcELoss();
    }

  return;
}


bool CDSTrack::DoLineFit( const int &n, const double *x, const double *y,
                              const double *w,
                              double &a, double &b, double &c, float &chi )
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
   
bool CDSTrack::FirstLineFitting(CDSHitMan *cdsMan)
{
  
  double x[40],y[40],weight[40];
  int numallhit=0;  
  for(int numlayer=0;numlayer<7;numlayer++)
    {
      int layer=0;
      if(numlayer<3) layer=numlayer+1;
      else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
      else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  double dl=hit(cdsMan,layer,n)->dl();
	  if(dl<ResolutionOfCDC) dl=ResolutionOfCDC;
	  x[numallhit]=hit(cdsMan,layer,n)->wx();
	  y[numallhit]=hit(cdsMan,layer,n)->wy();
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
  phi = meanphi;
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
  MathTools::PointToLine(zero,tmp_o,Line_d,dis,Line_o);
  param[2]=0;
  param[0]=dis;
  double tmpdeg=MathTools::CalcDeg(Line_o.x(),Line_o.y());
  param[1]=tmpdeg*TMath::DegToRad();
  tmpdeg+=90;
  if(tmpdeg>360) tmpdeg-=360;
  tmpdeg-=theta;tmpdeg=fabs(tmpdeg);if(tmpdeg>180) tmpdeg-=360;
  if(fabs(tmpdeg)>90 ) {param[0]*=-1;param[1]+=TMath::Pi();}
  if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi(); 
  SetHitPos(cdsMan);
  return true;  
}

   
bool CDSTrack::LineFitting(CDSHitMan* cdsMan)
{
  
  double x[40],y[40],weight[40];
  int numallhit=0;  
  for(int numlayer=0;numlayer<7;numlayer++)
    {
      int layer=0;
      if(numlayer<3) layer=numlayer+1;
      else if(3<=numlayer && numlayer<=4) layer=numlayer-3+8;
      else if(5<=numlayer && numlayer<=6) layer=numlayer-5+14;
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  double dl=hit(cdsMan,layer,n)->dl();
	  if(dl<0) dl=0.0001;
	  x[numallhit]=hit(cdsMan,layer,n)->x();
	  y[numallhit]=hit(cdsMan,layer,n)->y();
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
  MathTools::PointToLine(zero,tmp_o,Line_d,dis,Line_o);
  param[2]=0;
  param[0]=dis;

  double tmpdeg=MathTools::CalcDeg(Line_o.x(),Line_o.y());
  param[1]=tmpdeg*TMath::DegToRad();
  tmpdeg+=90;
  if(tmpdeg>360) tmpdeg-=360;
  tmpdeg-=theta;tmpdeg=fabs(tmpdeg);if(tmpdeg>180) tmpdeg-=360;
  if(fabs(tmpdeg)>90 ) {param[0]*=-1;param[1]+=TMath::Pi();}
  if(param[1]>2*TMath::Pi()) param[1]-=2*TMath::Pi(); 
  SetHitPos(cdsMan);

  return true;  
}


bool CDSTrack::FirstStereoLineFitting(CDSHitMan *cdsMan)
{

  double disZ[40];
  int numstlayer=0;

  for(int stlayer=0 ;stlayer<8;stlayer++)
    {
      int layer=0;
      if(0<=stlayer && stlayer<=3) layer=stlayer+4;
      else if(4<=stlayer && stlayer<=7) layer=stlayer+6;

      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  CDCHit *cdc= hit(cdsMan,layer,n);
	  TVector3 wpos,wposp;	  
	  wpos=cdc->wpos();	  
	  wposp=cdc->wposp();

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
	  hit(cdsMan,layer,n)->SetHitPosition(WPos.x(),WPos.y(),WPos.z());    

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
	  for(int n=0;n<(int)nTrackHit(layer);n++)
	    {
	      CDCHit *cdc=hit(cdsMan,layer,n);
	      
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
      
      double aa,bb,cc;      float chitmp=0;
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
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
	    {
	      double hitx= hit(cdsMan,layer,n)->x();
	      double hity= hit(cdsMan,layer,n)->y();
	      double hitz= hit(cdsMan,layer,n)->z();
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
	      hit(cdsMan,layer,n)->SetHitPosition(hitx,hity,hitz);   
	      numst++;
	    }
	  
	}
    }
#if TMPPARAM
  ChiSqrtmp[2]=ChiSqr;
#endif
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
#if TMPPARAM
  for(int n=0;n<5;n++) paramtmp[3][n]=param[n];
#endif
  SetHitPos(cdsMan);
  return true;
}


   
bool CDSTrack::StereoLineFitting(CDSHitMan *cdsMan)
{
  double x[40],y[40],z[40],weight[40];
  int numallhit=0;  
  for(int layer=1;layer<=15;layer++)
    {
	      
      if( (4<=layer &&layer<=7) || (10<=layer && layer<=13) )
	{ 
	  for(int n=0;n<(int)nTrackHit(layer);n++)
	    {
	      double dl=hit(cdsMan,layer,n)->dl();
	      if(dl<0) dl=0.0001;
	      x[numallhit]=hit(cdsMan,layer,n)->x();
	      y[numallhit]=hit(cdsMan,layer,n)->y();
	      z[numallhit]=hit(cdsMan,layer,n)->z();
	      double tilt=hit(cdsMan,layer,n)->tilt()*TMath::DegToRad();
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
  double a,b,c;
  float chi;
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
#if TMPPARAM
  for(int n=0;n<5;n++) paramtmp[4][n]=param[n];
#endif
  CalcLineChiSqr(cdsMan);
  if( ChiSqr >99999 ) 
    {
#if DEBUG
	std::cout<<"Too big chi2/dof sl"<<std::endl;
#endif
	return false;
    }
  SetHitPos(cdsMan);
  return true;  
}

bool CDSTrack::StereoLineFitting2(CDSHitMan *cdsMan)
{
  int allhit=0;
  TVector3 wirepos[50],wiredir[50];
  double weight[50],drift[50];
  
  for(int layer=1;layer<=15;layer++)
    {
      for(int n=0;n<(int)nTrackHit(layer);n++)
	{
	  CDCHit *cdc=hit(cdsMan,layer,n);
	  wirepos[allhit]=(cdc->wpos()+cdc->wposp())*0.5;
	  wiredir[allhit]=(cdc->wpos()-cdc->wposp());
	  weight[allhit]=ResolutionOfCDC;
	  drift[allhit]=cdc->dl();
	  allhit++;
	}      
    }
  double tmppar[6]={Line_o.X(),Line_o.Y(),Line_o.Z(),
		    Line_d.X(),Line_d.Y(),Line_d.Z()};
  LineFit *fit=new LineFit(tmppar,wirepos,wiredir,weight,drift,allhit);
  fit->GetParameters(tmppar); 
  delete fit;
  
  Line_o.SetXYZ(tmppar[0],tmppar[1],tmppar[2]);
  Line_d.SetXYZ(tmppar[3],tmppar[4],tmppar[5]);

  fittinglevel=4;
  param[2]=0.;
#if TMPPARAM
  for(int n=0;n<5;n++) paramtmp[4][n]=param[n];
#endif
  CalcLineChiSqr(cdsMan);
  if( ChiSqr >99999 ) 
    {
#if DEBUG
	std::cout<<"Too big chi2/dof sl"<<std::endl;
#endif
	return false;
    }
  SetHitPos(cdsMan);
  return true;  
}

void CDSTrack::CalcLineChiSqr(CDSHitMan* cdsMan)
{
//   TVector3 P( 0., -Cl/Bl, (Bp*Cl/Bl-Dp)/Cp );
//   TVector3 Q( Bl, -Al, Bl*(Bp*Al/Bl-Ap)/Cp );
  TVector3 WirePosition,WirePositionp, WireDirection;
  TVector3 Xest, Next;
  CDCHit *cdc = 0;
  double dist;
  double chi = 0;
  int nallhit=0;
  for(int layer=1;layer<=15;layer++)
    {
      for( int i=0; i<(int)nTrackHit(layer); i++ ){
	cdc = hit(cdsMan,layer,i);
 	double dl = fabs(cdc->dl());
	WirePosition=cdc->wpos();
	WirePositionp=cdc->wposp();
	WireDirection = (WirePositionp-WirePosition).Unit();
	//	MathTools::LineToLine( P, Q, WirePosition, WireDirection, dl, dist, Xest, Next );
	MathTools::LineToLine( Line_o, Line_d, WirePosition, WireDirection, dl, dist, Xest, Next );
	chi += (dl-dist)*(dl-dist)/(ResolutionOfCDC*ResolutionOfCDC);
	//	std::cout<<"dl dis "<<dl<<" "<< dist<<std::endl;
	
	nallhit++;
      }
    }
  chi /= (double)(nallhit-4);
  ChiSqr=chi;
  dof=nallhit-4;
}

void CDSTrack::CalcELoss(double mass){
#if DEBUG
  std::cout<<"==============================start CalcEloss()"<<std::endl;
#endif
  if(mass<0)  mass=Mass();
  if(mass<=0){
#if DEBUG
    std::cout<<"\t !!! mass= "<<mass<<std::endl;
#endif
    return ;
  }
  int sign=1; // add energy loss
  double momout=Momentum();
  double tmpmom=momout;
  if(TMath::Abs(tmpmom)>10){
#if DEBUG
    std::cout<<"\t !!! momentum= "<<momout<<std::endl;
#endif
    return;
  }
  //  TVector tmppos;
  double tmppar[5];
  //  for(int i=0;i<5;i++) tmppar[i]=param[i];
  TVector3 pos1;
  GetParameters(CID_CDC,tmppar,pos1);
  double helixphi=MathTools::CalcHelixPhi(pos1,tmppar);

  double newpar[5];
  double length=0;
  int id=CID_CDC;

  int count=0;
  double tmptof,tmpl;
  double aaaa[5];
  TVector3 bbbb;

  //  std::cout<<"=============================="<<std::endl;
  while(1){
    // std::cout<<"                       ";
    // for(int i=0;i<5;i++)
    //   std::cout<<tmppar[i]<<"  ";
    // std::cout<<std::endl;
    count++;
    if(count>50) return;
    if( !ELossTools::CalcHelixElossToNextBoundary(param,pos1,tmpmom,mass,momout,tmptof,tmpl,id,sign) ) return;
    tmpmom=momout;
    if( pos1.Perp()>20 || TMath::Abs(pos1.Z())>100 ) return;
    helixphi=MathTools::CalcHelixPhiatZ(tmppar,pos1.Z());
    //    std::cout<<std::setw(10)<<helixphi<<"  "<<std::setw(5)<<id<<"      "<<pos1.Perp()<<std::endl;
    if(GetParameters(id,aaaa,bbbb))
      if(TMath::Abs(helixphi)>TMath::Pi() || tmppar[2]/TMath::Abs(tmppar[2])*helixphi >0 )
	return;    
    MathTools::ChangePivot(TVector3(0,0,0),pos1,tmppar,newpar,charge());
    newpar[2]=magneticfield*Const/100./(momout/sqrt(1+newpar[4]*newpar[4]));
    MathTools::ChangePivot(pos1,TVector3(0,0,0),newpar,tmppar,charge());
    AddParameters((int)id,tmppar,pos1);
  }
}

bool CDSTrack::CalcVertexTimeLength(const TVector3 &pos,const TVector3 &dir,const double &mass,TVector3 &lpos, TVector3 &hpos,double &tof, double &length, bool ADDPAR)
{
  TString mat,newmat;
  //  if(mass<0)  mass=Mass();
  if(mass<=0) return false;
  int sign=1; // add energy loss
  double momout=Momentum();
  double tmpmom=momout;
  if(TMath::Abs(tmpmom)>10)  return false;
  double tmppar[5];
  int id=CID_CDC;
  TVector3 pos1,pos2;
  GetParameters(id,tmppar,pos1);
  GeomTools::GetIDMat(pos1,mat);
  double newpar[5];
  length=0; tof=0;

  int count=0;
  double tmptof,tmpl,dis;
  double aaaa[5];
  TVector3 bbbb;
  while(1){
    count++;
    if(count>30) return false;    
    if(!MathTools::LineToHelix(pos,dir,tmppar,lpos,hpos,dis) ) return false;
#if 0
    if(ADDPAR){
      std::cout<<"vertex r,z"<<hpos.Perp()<<"   "<<hpos.Z()<<std::endl;

    }
#endif
    // return new Volume ID
    if(!GeomTools::HelixStepToNextVolume(tmppar,pos1,pos2,tmpl,newmat,id)) return false;
    double tmpl2=MathTools::CalcHelixArc(tmppar,pos1,hpos);
    if(tmpl>tmpl2){
#if DEBUG
      std::cout<<tmpl<<" >  "<<tmpl2<<std::endl;
#endif
      break;
    }
    pos1=pos2;
   //    if(GeomTools::IsSameVolumeHelix(tmppar,pos1,hpos)) break;
    if( !ELossTools::CalcdE(tmpmom,mass,tmpl,mat,momout,sign,tmptof ) ) return false;
    //    std::cout<<mat<<"  "<<mass<<"  "<<tmpl<<"  "<<tmpmom<<"  -> "<<momout<<std::endl;
    //    if(!ELossTools::CalcHelixElossToNextBoundary(tmppar,pos1,tmpmom,mass,momout,tmptof,tmpl,id,sign) ) return false;
    tmpmom=momout;
    mat=newmat;
    //    std::cout<<"id,tmpl,tmptof, tottof: "<<id<<"  "<<tmpl<<"  "<<tmptof<<"  "<<tof<<std::endl;
    length+=tmpl;
    tof+=tmptof;
    if(id==CID_CDCCFRP||id==CID_CDC){
      length=0;
      tof=0;
    }
    MathTools::ChangePivot(TVector3(0,0,0),pos1,tmppar,newpar,charge());
    newpar[2]=magneticfield*Const/100./(momout/sqrt(1+newpar[4]*newpar[4]));
    MathTools::ChangePivot(pos1,TVector3(0,0,0),newpar,tmppar,charge());
    if(ADDPAR)
      AddParameters((int)id,tmppar,pos1);
  }
  mat=newmat;
  tmpl=MathTools::CalcHelixArc(tmppar,pos1,hpos);
  ELossTools::CalcdE(tmpmom,mass,tmpl,mat,momout,sign,tmptof);
  if(id!=CID_CDCCFRP&&id!=CID_CDC){
    tof+=tmptof;
    length+=tmpl;
  }
  return true;
}

bool CDSTrack::GetCDHHit(CDSHitMan *cdsMan, int &seg, double &time){
  seg=0;
  time=999;
  for(int i=0;i<nCDHHit();i++){
    HodoscopeLikeHit *hit=CDHHit(cdsMan,i);
    if(hit&&time>hit->ctmean()){
      time=hit->ctmean();
      seg=hit->seg();
    }
  }
  if(seg==0) return false;
  else return true;
}
bool CDSTrack::GetIHHit(CDSHitMan *cdsMan, int &seg, double &time){
  seg=0;
  time=999;
  for(int i=0;i<nIHHit();i++){
    HodoscopeLikeHit *hit=IHHit(cdsMan,i);
    if(hit&&time>hit->ctmean()){
      time=hit->ctmean();
      seg=hit->seg();
    }
  }
  if(seg==0) return false;
  else return true;
}
void CDSTrack::Print()
{
  for(int i=0;i<5;i++) std::cout<<param[i]<<"  ";
  std::cout<<std::endl;
  for(int i=0;i<nParamSets();i++){
    int id=-1;
    double tmp[5]={0,0,0,0,0};
    TVector3 pos;
    GetNthParameters(i,id,tmp,pos);
    std::cout<<"id  "<<id<<"  ";
    for(int i=0;i<5;i++) std::cout<<tmp[i]<<"  ";
    std::cout<<pos.Perp()<<"  ";
    pos.Print();   
  }
}
