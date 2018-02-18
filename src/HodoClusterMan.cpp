// HodoClusterMan.cpp
#include <map>
#include <algorithm>

#include "GlobalVariables.h"
#include "HodoClusterMan.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"

ClassImp(HodoHit);
ClassImp(HodoClusterMan);
ClassImp(HodoCluster);
//bool SortHodoHit(const HodoHit& lhs, const HodoHit& rhs) { return lhs.seg()<rhs.seg(); }
#define DEBUG 0
// ----------------------------- //
// class HodoCluster             //
// ------------------------------//
HodoHit::HodoHit( int seg, double time, double edep, const TVector3 &tmppos):
  Seg(seg),Time(time),Edep(edep),Pos(tmppos)
{
}
HodoCluster::HodoCluster(const int &cid, const int &cluid):
  CounterID(cid),ClusterID(cluid),CTimeMean(DEFAULTD),EdepTotal(DEFAULTD),aTrackID(DEFAULTI),FLAG(false)
{
}

void HodoCluster::SetHit( const HodoscopeLikeHit &hit )
{ 
  hitContainer.push_back(HodoHit(hit.seg(),hit.ctmean(),hit.emean(),hit.pos()));
  return;
}

void HodoCluster::Calc( )
{
  EdepTotal=0.;
  double tmptime=0.;
  for(int i=0;i<nhit();i++){
    EdepTotal+=hit(i)->edep();
    tmptime+=hit(i)->time();
  }
  CTimeMean=tmptime/nhit();
}

int HodoCluster::nhit( const double &threshold){
  int tmpn=0;
  for(int i=0;i<nhit();i++){
    HodoHit *tmphit=hit(i);
    if(tmphit&&tmphit->edep()>threshold) tmpn++;
  }
  return tmpn; 
}

bool HodoCluster::first(  double threshold , int &seg, double &time, TVector3 &tmppos)
{
  std::sort(hitContainer.begin(),hitContainer.end());//,SortHodoHit);
  for(int i=0;i<nhit();i++){
    HodoHit *tmphit=hit(i);
    if(tmphit&&tmphit->edep()>threshold){
      seg=tmphit->seg();
      time=tmphit->time();
      tmppos=tmphit->pos();
      return true;
    }
  }
  return 0;
}

bool HodoCluster::firsttimeinlayer(  double threshold , int &seg, double &time, TVector3 &tmppos)
{
  int minlayer=9;
  double mintime=999;
  std::sort(hitContainer.begin(),hitContainer.end());//,SortHodoHit);
  for(int i=0;i<nhit();i++){
    HodoHit *tmphit=hit(i);
    int tmpseg=tmphit->seg();
    int layer=(tmpseg-1)/16;
    //    std::cout<<tmpseg<<"  "<<layer<<std::endl;
    double tmptime=tmphit->time();
    if(tmphit->edep()>threshold&&(layer<minlayer||(layer==minlayer&&tmptime<mintime))){
      seg=tmpseg;
      time=tmptime;
      tmppos=tmphit->pos();
      mintime=tmptime;
      minlayer=layer;
      //      return true;
    }
  }
  if(mintime>300)   return 0;
  else return true;
}

bool HodoCluster::firsttime(  double threshold , int &seg, double &time, TVector3 &tmppos)
{
  double mintime=999;
  std::sort(hitContainer.begin(),hitContainer.end());//,SortHodoHit);
  for(int i=0;i<nhit();i++){
    HodoHit *tmphit=hit(i);
    int tmpseg=tmphit->seg();
    double tmptime=tmphit->time();
    double edep=tmphit->edep();
    //    std::cout<<"seg,time,edep  "<<tmpseg<<"  "<<tmptime<<"  "<<edep<<"  ";
    if(edep>threshold&&tmptime<mintime){
      seg=tmpseg;
      time=tmptime;
      tmppos=tmphit->pos();
      mintime=tmptime;
      //      return true;
    }
  }
  if(mintime>300)   return 0;
  else return true;
}

bool HodoCluster::firsttime(  double threshold , int &seg, double &time, TVector3 &tmppos, double &tmpedep)
{
  double mintime=999;
  std::sort(hitContainer.begin(),hitContainer.end());//,SortHodoHit);
  for(int i=0;i<nhit();i++){
    HodoHit *tmphit=hit(i);
    int tmpseg=tmphit->seg();
    double tmptime=tmphit->time();
    double edep=tmphit->edep();
    //    std::cout<<"seg,time,edep  "<<tmpseg<<"  "<<tmptime<<"  "<<edep<<"  ";
    if(edep>threshold&&tmptime<mintime){
      seg=tmpseg;
      time=tmptime;
      tmppos=tmphit->pos();
      mintime=tmptime;
      tmpedep=edep;
      //      return true;
    }
  }
  if(mintime>300)   return 0;
  else return true;
}

bool HodoCluster::second(  double threshold , int &seg, double &time, TVector3 &tmppos)
{
  std::sort(hitContainer.begin(),hitContainer.end());//,SortHodoHit);
  HodoHit *tmphit=hit(0);
  if(!tmphit) return 0;
  int tmpseg=tmphit->seg();
  int layer=tmpseg/16+1;
  //  std::cout<<tmpseg<<"  "<<layer<<std::endl;
  seg=-1;
  double tmpedep=0;
  for(int i=1;i<nhit();i++){
    HodoHit *tmphit=hit(i);
    if(tmphit){
      //       std::cout<<tmphit->seg()<<"   "<<layer<<std::endl;
      if((tmphit->seg()>layer*16) && (tmphit->seg()<=(layer+1)*16) && tmphit->edep()>threshold){
	if(tmpedep<tmphit->edep()){
	  seg=tmphit->seg();
	  time=tmphit->time();
	  tmppos=tmphit->pos();
	  tmpedep=tmphit->edep();
	}
      }
    }
  }
  if(seg>0) return true;
  return false;
}

void HodoCluster::Clear() { 
  hitContainer.clear();
  ClusterID=DEFAULTI;
  CTimeMean=DEFAULTD;  
  EdepTotal=DEFAULTD;  
  aTrackID=DEFAULTI;
  FLAG=false;
}

// ----------------------------- //
// class HodoClusterMan          //
// ------------------------------//
HodoClusterMan::HodoClusterMan( BeamLineHitMan *blMan, CDSHitMan* cdsMan, ConfMan *conf )
{
  Clustering(blMan,conf);
  Clustering(cdsMan,conf);
}

void HodoClusterMan::Calc(  )
{
  HodoClusterContainer::iterator it = Container.begin();
  while( it != Container.end() )
    {
      for( unsigned int i=0;i<it->second.size(); i++ ){
	  it->second[i].Calc();
      }
    }
}

void HodoClusterMan::DeleteCluster( const int &cid,const int &i )
{
  HodoClusterContainer::iterator it=Container.find(cid);
  if( it != Container.end() ){
    ClusterContainer::iterator it2=it->second.begin();
    for(int j=0;j<i;j++ ) it2++;
    it->second.erase(it2);
  }
}

int HodoClusterMan::ncluster( int cid, HitMan *hitMan, double threshold)
{
  int tmpn=0;
  HodoClusterContainer::iterator it=Container.find(cid);
  if( it != Container.end() )
    for( unsigned int i=0;i<it->second.size(); i++ )
      if( it->second[i].nhit(threshold)  > 0) tmpn++;
  return tmpn;
}

void HodoClusterMan::Clear()
{
  Container.clear();
  nccvcID.clear();
}

bool HodoClusterMan::Clustering( HitMan *hitMan , ConfMan *conf, const int &cid )
{
  TVector3 pos1,pos2;
  conf->GetGeomMapManager()->GetPos(cid,1,pos1);
  if(cid==CID_NC)
    conf->GetGeomMapManager()->GetPos(cid,18,pos2);
  else
    conf->GetGeomMapManager()->GetPos(cid,2,pos2);
  double dist=(pos2-pos1).Mag()+0.01;
  double wtime=3;
  if(cid==CID_DEF||cid==CID_IH) wtime=10;
  //  std::cout<<hitMan->nHodo(cid)<<std::endl;
  for(int i=0;i<hitMan->nHodo(cid);i++){
    HodoscopeLikeHit *hit=hitMan->Hodoi(cid,i);
    if(!hit || !hit->CheckRange()) continue;
#if DEBUG
    std::cout<<cid<<"  "<<hit->seg()<<"  "<<hit->ctmean()<<std::endl;
#endif
    hit->ctmean();
    TVector3 pos;//=hit->pos();
    conf->GetGeomMapManager()->GetPos(cid,hit->seg(),pos);
    int nclu=ncluster(cid);
    bool ADD=false;
    for(int iclu=0;iclu<nclu;iclu++){
      //      std::cout<<"clustertime: "<<cluster(cid,iclu)->GetCTimeMean()<<std::endl;
      if(TMath::Abs( cluster(cid,iclu)->GetCTimeMean() - hit->ctmean()) < 5. ){
	for(int ihit=0;ihit<cluster(cid,iclu)->nhit();ihit++){
	  //	  TVector3 tmppos=cluster(cid,iclu)->hit(blMan,ihit)->pos();
	  TVector3 tmppos;
	  conf->GetGeomMapManager()->GetPos(cid,cluster(cid,iclu)->hit(hitMan,ihit)->seg(),tmppos);
	  if((tmppos-pos).Mag()<dist){
	    cluster(cid,iclu)->SetHit(*hit);
	    cluster(cid,iclu)->Calc();
	    ADD=true;
	    break;
	  }
	}
      }
    }
    if(!ADD){
      HodoCluster* tmpclu=new HodoCluster(cid,nclu);
      tmpclu->SetHit(*hit);
      tmpclu->Calc();
      SetCluster(cid,*tmpclu);
      delete tmpclu;
    }
  }    				     
#if DEBUG
  for(int i=0;i<ncluster(cid);i++){
    std::cout<<"----cluster "<<i<<"  dist "<<dist<<std::endl;
    for(int ihit=0;ihit<cluster(cid,i)->nhit();ihit++){
      HodoscopeLikeHit *hit=cluster(cid,i)->hit(hitMan,ihit);
      std::cout<<cid<<"  "<<hit->seg()<<"  "<<hit->ctmean()<<"  "<<hit->emean()<<"  ";
      hit->pos().Print();
    }
  }
  std::cout<<"=========================="<<std::endl;
#endif
  return true;
}

bool HodoClusterMan::Clustering( BeamLineHitMan *blMan , ConfMan *conf )
{
  const int DetList[]={/*CID_BPD,CID_DEF,CID_BVC,*/CID_CVC,CID_NC};
  const int nDet=sizeof(DetList)/sizeof(int);
  for(int i=0;i<nDet;i++){
    Clustering(blMan,conf,DetList[i]);
  }
  return true;
}

bool HodoClusterMan::Clustering( CDSHitMan *cdsMan , ConfMan *conf )
{
  const int nDet=2;
  const int DetList[nDet]={CID_CDH,CID_IH};
  for(int i=0;i<nDet;i++){
    Clustering(cdsMan,conf,DetList[i]);
  }
  return true;
}

HodoCluster *HodoClusterMan::ncfirst( double threshold , int &seg, double &time, TVector3 &tmppos)
{
  double tmptime=999;
  int tmpnum=-1;
  for(int i=0;i<ncluster(CID_NC);i++){
    if(cluster(CID_NC,i)->first(threshold,seg,time,tmppos))
      if(time<tmptime){
	//	tmptime=tmptime; //bug fixed 20130907 
	tmptime=time;
	tmpnum=i;
      }
  }
  if(tmpnum==-1) return 0;
  else{
    cluster(CID_NC,tmpnum)->first(threshold,seg,time,tmppos);
    return cluster(CID_NC,tmpnum);
  }
}

HodoCluster *HodoClusterMan::ncfirsttime( double threshold , int &seg, double &time, TVector3 &tmppos)
{
  double tmptime=999;
  int tmpnum=-1;
  for(int i=0;i<ncluster(CID_NC);i++){
    if(cluster(CID_NC,i)->firsttime(threshold,seg,time,tmppos))
      if(time<tmptime){
	//	tmptime=tmptime; //bug fixed 20130907 
	tmptime=time;
	tmpnum=i;
      }
  }
  if(tmpnum==-1) return 0;
  else{
    cluster(CID_NC,tmpnum)->firsttime(threshold,seg,time,tmppos);
    return cluster(CID_NC,tmpnum);
  }
}

HodoCluster *HodoClusterMan::ncfirsttime( double threshold , int &seg, double &time, TVector3 &tmppos, double &edep)
{
  double tmptime=999;
  int tmpnum=-1;
  for(int i=0;i<ncluster(CID_NC);i++){
    if(cluster(CID_NC,i)->firsttime(threshold,seg,time,tmppos))
      if(time<tmptime){
	//	tmptime=tmptime; //bug fixed 20130907 
	tmptime=time;
	tmpnum=i;
      }
  }
  if(tmpnum==-1) return 0;
  else{
    cluster(CID_NC,tmpnum)->firsttime(threshold,seg,time,tmppos,edep);
    return cluster(CID_NC,tmpnum);
  }
}

