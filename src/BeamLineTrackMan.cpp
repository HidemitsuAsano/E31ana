// BeamLineTrackMan.cpp
#include <map>

#include "GlobalVariables.h"
#include "BeamLineTrackMan.h"
#include "TVector2.h"
#include "TMath.h"

static const double SpatialResolutionOfBLDC=0.02; // [cm]

ClassImp(BeamLineTrackMan);

#define DEBUG 0
#define DEBUG2 0
#define CHECKRANGE 0
#define SLOPECUT 1
#define DEL 1
#define KILLLAYER 0

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
  for( BeamLineTrackContainer::const_iterator it=right.FDC1TrackContainer.begin(); it!=right.FDC1TrackContainer.end(); it++ ){
    FDC1TrackContainer.push_back( (*it) );
  }
  for(int i=0;i<10;i++) STATUS[i]=right.STATUS[i];
}

BeamLineTrackMan::~BeamLineTrackMan()
{
  //  this->Clear();
}

void BeamLineTrackMan::Clear()
{
  for(int i=0;i<10;i++) STATUS[i]=-1;
  BLC1TrackContainer.clear();
  BLC1aTrackContainer.clear();
  BLC1bTrackContainer.clear();

  BLC2TrackContainer.clear();
  BLC2aTrackContainer.clear();
  BLC2bTrackContainer.clear();

  BPCTrackContainer.clear();
  FDC1TrackContainer.clear();
}

BeamLineTrackContainer* BeamLineTrackMan::BLDCTrackContainer(const int &cid){
  switch(cid){
  case CID_BLC1a: return &BLC1aTrackContainer;
  case CID_BLC1b: return &BLC1bTrackContainer;
  case CID_BLC2a: return &BLC2aTrackContainer;
  case CID_BLC2b: return &BLC2bTrackContainer;
  case CID_BPC:   return &BPCTrackContainer;
  case CID_BLC1:  return &BLC1TrackContainer;
  case CID_BLC2:  return &BLC2TrackContainer;
  case CID_FDC1:  return &FDC1TrackContainer;
  default:        return 0;
  }
}  
int BeamLineTrackMan::status(const int &cid){
  switch(cid){
  case CID_BLC1a: return STATUS[0];
  case CID_BLC1b: return STATUS[1];
  case CID_BLC2a: return STATUS[2];
  case CID_BLC2b: return STATUS[3];
  case CID_BPC:   return STATUS[4];
  case CID_BLC1:  return STATUS[5];
  case CID_BLC2:  return STATUS[6];
  case CID_FDC1:  return STATUS[7];
  default:        return -1;
  }
}

int BeamLineTrackMan::SetStatus(const int &cid, const int &sta){
  switch(cid){
  case CID_BLC1a:  STATUS[0]=sta;   break;
  case CID_BLC1b:  STATUS[1]=sta;   break;
  case CID_BLC2a:  STATUS[2]=sta;   break;
  case CID_BLC2b:  STATUS[3]=sta;   break;
  case CID_BPC:    STATUS[4]=sta;   break;
  case CID_BLC1:   STATUS[5]=sta;   break;
  case CID_BLC2:   STATUS[6]=sta;   break;
  case CID_FDC1:   STATUS[7]=sta;   break;
  }
  return sta;
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

bool BeamLineTrackMan::DoTracking( BeamLineHitMan *blMan, ConfMan *conf, const bool &SEMILOCAL, const bool &TIMING)
{
#if 0
  std::cout << "!!! BeamLineTrackMan::DoTracking()" << std::endl;
#endif 
  TString option="";
  if( !TIMING ) option += "notiming";
  LocalTracking( blMan, conf, CID_BPC  , option );
  LocalTracking( blMan, conf, CID_BLC1a, option );
  LocalTracking( blMan, conf, CID_BLC1b, option );
  LocalTracking( blMan, conf, CID_BLC2a, option );
  LocalTracking( blMan, conf, CID_BLC2b, option );
  LocalTracking( blMan, conf, CID_FDC1 , option + "noslope" + "notiming" );
  if(SEMILOCAL){
    option+="linear";
    LocalTracking( blMan, conf, CID_BLC1, option );
    LocalTracking( blMan, conf, CID_BLC2, option );
  }
  return 1;
}
// ------------------------------------------------------------------
int BeamLineTrackMan::LocalTracking( BeamLineHitMan *blMan, ConfMan *conf, const int &id, TString option )
{
#if 0
  std::cout << "!!! BeamLineTrackMan::LocalTracking()\tcid = " <<id<< std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
  int MaxNumOfHitsInBLDC=BLDCParam->GetMaxBLDCHit();
  //  int MaxNumOfHitsInLayer=BLDCParam->GetMaxHitInLayer();
 
  int nallhit = 0;
  int NumOfBLDCLayers=DetectorList::GetInstance()->GetNlayers(id);
  BeamLineTrackContainer TmpTrackContainer[2];
  BeamLineTrackContainer TmpTrackContainer2[2];
  BeamLineTrackContainer TmpTrackContainer3;
  int id2[2];

  int nbldc=1;
  int nxy=2;
  int nud=2;

  id2[0]=id;
  id2[1]=id;
  if( id ==CID_BLC2 ){
    nbldc=2;
    id2[0]=CID_BLC2a; 
    id2[1]=CID_BLC2b; 
  }
  else if(id==CID_BLC1){
    id2[0]=CID_BLC1a; 
    id2[1]=CID_BLC1b;
    nbldc=2;
  }
  else if(id==CID_FDC1){
    nxy=1;
    nud=3;
  }
  
  //----------------------------
  //check total number and busy layers
  //----------------------------
  for(int ic=0;ic<nbldc;ic++){
    for( int layer=1; layer<=NumOfBLDCLayers; layer++ ){
      int tmp_nhit=0;
      tmp_nhit = blMan->nBLDC(id2[ic],layer);
      if( BLDCParam->layerdead(id2[ic],layer)){
	for(int i=0;i<tmp_nhit;i++)
	  blMan->BLDC(id2[ic],layer,i)->SetStatus(-1);	  
	continue; 
      }
#if KILLLAYER
      if( BLDCParam->layerkilled(id2[ic],layer)||
	  tmp_nhit>BLDCParam->GetMaxHitInLayer()){
	for(int i=0;i<tmp_nhit;i++){
	  blman->BLDC(id2[ic],layer,i)->SetStatus(0);
	  continue;       
	}
      }
#endif
      nallhit += tmp_nhit;
#if DEBUG
      std::cout<<"layer,nhit,nallhit"
	       <<"\t"<<layer
	       <<"\t"<<tmp_nhit
	       <<"\t"<<nallhit
	       <<std::endl;
#endif
    }
  }
  if( nallhit > MaxNumOfHitsInBLDC * nbldc ){
#if DEBUG
    std::cout<<"too many hits in BLDC !!!"<<std::endl;
    std::cout<<nallhit<<"\t>\t"<<MaxNumOfHitsInBLDC*nbldc<<std::endl;
#endif
    return SetStatus(id,2);
  }
  //----------------------------
  //clustering
  //----------------------------
  BLDCClusterMan *clMan[2];
  clMan[0] = new BLDCClusterMan();
  clMan[1] = new BLDCClusterMan();

  for(int ic=0;ic<nbldc;ic++){
    if( !Clustering( blMan, clMan[ic], conf, id2[ic], ic ) ){
      delete clMan[0]; 
      delete clMan[1]; 
      std::cout<<"clustring failed !!!\tid"<<id2[ic]<<std::endl;
      //    return ntrack; 
    }
  }
#if DEBUG  // for debug
  //  if(id==CID_FDC1)
  {
    std::cout<<"clustering done"<<std::endl;
    for(int ic=0;ic<nbldc;ic++){
      for(int xy=0;xy<nxy;xy++){
	std::cout<<"id,xy:  "<<id2[ic]<<"\t"<<xy<<std::endl;
	for(int islay=0;islay<nud;islay++){
	  std::cout<<"islay,ncluster:  "<<islay<<"\t"<<clMan[ic]->ncluster(islay,xy)<<std::endl;
	  for( int icu=0; icu<clMan[ic]->ncluster(islay,xy); icu++ ){
	    BLDCCluster *cl=clMan[ic]->cluster(islay,xy,icu);
	    std::cout<<"icluster,time: "<<icu<<"\t"<<cl->GetTimeMean()<<std::endl;
	    for( int i=0; i<cl->nhit(); i++ ) {
	      ChamberLikeHit *hit = cl->hit(i);
	      std::cout<<"\ticl,ihit,cid,layer,wire,dt,dl:   "<<icu<<"\t"<<i<<"\t"<<hit->cid()<<"\t"<<hit->layer()<<"\t"<<hit->wire()<<"\t"<<hit->dt()<<"\t"<<hit->dl()<<std::endl;
	    }
	  }
	}
      }
    }
  }
#endif
  //----------------------------
  // require more than 3 slayers which have at least one cluster
  //----------------------------
  int ncomb[2]={1,1};
  bool CLUSTER=true;
  int ncluster[2][4]={};
  
  int minslay=2;
  if(nbldc==2||nxy==1) minslay=3;

  for(int xy=0;xy<nxy;xy++){
    int nslay=0;
    for(int ic=0;ic<nbldc;ic++) {
      for(int iud=0;iud<nud;iud++){
	if(clMan[ic]->ncluster(iud,xy)>0&&clMan[ic]->ncluster(iud,xy)<10){
	  ncomb[xy] *= clMan[ic]->ncluster(iud,xy);
	  ncluster[xy][2*ic+iud]=clMan[ic]->ncluster(iud,xy);
	  nslay++;
	}
      }
    }
    if(nslay<minslay) CLUSTER=false;
  }

  //----------------------------
  // search track candidate in X/Y layers w/ MWPC mode
  //----------------------------
  if(CLUSTER){
    for(int xy=0;xy<nxy;xy++){
      int clnum[4]={-1,-1,-1,-1};
      int tmp_ntrack=0;
      bool MINHIT=false;
      bool SLOPE=false;
      bool TIMING=false;
      bool FIT=false;
      bool CHI2=false;
      
      for(int icomb=0;icomb<ncomb[xy];icomb++){
	BLDCCluster *cl[4];
	int tmp=icomb;
	int ncl=0;
	LocalTrack *track = new LocalTrack();
	
	for(int ic=0;ic<nbldc;ic++){
	  for(int iud=0;iud<nud;iud++){
	    int islay=2*ic+iud;
	    if(ncluster[xy][islay]>0){
	      clnum[islay]=tmp%ncluster[xy][islay];
	      tmp=tmp/ncluster[xy][islay];
	      cl[ncl] = clMan[ic]->cluster( iud, xy, clnum[islay] );
	      ncl++;
	    }
	  }
	}
	int tmp_nhit=0;	
	for(int icl=0;icl<ncl;icl++)	
	  tmp_nhit+=track->SetCluster(xy,cl[icl]);
	//	for(int i=0;i<track->nhit();i++)
	//	  std::cout<<track->hit(i)->wire()<<"\t";
	if( tmp_nhit<BLDCParam->GetMinHit(xy,id) ) {
	  //	  if(id==CID_FDC1)
	  //	  std::cout<<tmp_nhit<<"\t"<<BLDCParam->GetMinHit(xy,id)<<std::endl;
	  //	if( tmp_nhit<5 ) {
	  //	  std::cout<<"too few hit number "<<tmp_nhit<<"  "<<BLDCParam->GetMinHit(xy,id)<<std::endl;
	  delete track;	  continue;
	}	
	MINHIT=true;
	track->CalcTrackTime();
#if DEBUG
	std::cout<<"\t";
	for(int i=0;i<track->nhit();i++)
	  std::cout<<track->hit(i)->wire()<<"\t";
	for(int i=0;i<track->ncluster(xy);i++)
	  std::cout<<track->cluster(xy,i)->GetTimeMean()<<"\t";
	std::cout<<std::endl;
#endif	
#if 0
	while( !option.Contains("notiming")
	       && track->GetTrackTimeRMS() > 20 ) {
	  track->DeleteOffTimingCluster();
	  std::cout<<"!!! TrackTime "<<track->GetTrackTimeRMS()<<std::endl;
	  if( (track->ncluster(0)+track->ncluster(1)) < minslay ) break;
	}	
	if( (track->ncluster(0)+track->ncluster(1)) < minslay ){
	  delete track; continue;
	}
#else
	if( !option.Contains("notiming")
	    && track->GetTrackTimeRMS() > 20 ) {
	  delete track; continue;
	}
#endif
	TIMING=true;
	
	if( !track->LeastSquareFit( conf, xy, "mwpc" ) ) {
	  std::cout<<"MWPC fit failed"<<std::endl;
	  delete track;	  continue;
	}
	FIT=true;
	
	double slope=track->b();
	if(xy) slope=track->e();
	
	//	std::cout<<"slope: "<<slope<<std::endl;
#if SLOPECUT
	if( !option.Contains("noslope") && TMath::Abs(slope) > 0.2  ) {
	  //	  std::cout<<"slope !!!"<<std::endl;
	  delete track;   continue;
	}	  
	SLOPE=true;
#endif
	if( track->chi2(xy)>BLDCParam->GetMaxChiMWPCFit() && id!=CID_FDC1) { 
	  //	  std::cout<<"mwpc chi !!!"<<std::endl;
	  delete track;   continue;
	}	
	CHI2=true;  
	
	TmpTrackContainer[xy].push_back(*track);
	//	std::cout<<"OK !!!"<<std::endl;
	delete track;
	tmp_ntrack++;	 
      }//icomb
      if( tmp_ntrack==0 ) {
	//sstd::cout<< "xy: "<<xy<<" no track !!!"<<std::endl;
	TmpTrackContainer[0].clear();
	TmpTrackContainer[1].clear();
	delete clMan[0]; 
	delete clMan[1]; 
	if(!MINHIT) return SetStatus(id,5);
	if(!TIMING) return SetStatus(id,6);
	if(!FIT) return SetStatus(id,7);
#if SLOPECUT
	if(!SLOPE)     return SetStatus(id,9);
#endif
	if(!CHI2)     return SetStatus(id,8);
	return SetStatus(id,10);
      }
#if 0
      if(id==CID_FDC1)
      std::cout<<"xy, ntrack   "<<xy<<"\t"<<tmp_ntrack<<std::endl;
#endif
    }//xy
  }else{
    TmpTrackContainer[0].clear();
    TmpTrackContainer[1].clear();
    delete clMan[0]; 
    delete clMan[1]; 
    return SetStatus(id,4);
  }// cluster OK
  delete clMan[0]; 
  delete clMan[1]; 
  // Search Minimum Chisquare Track in X/Y layers
  
  for(int xy=0;xy<nxy;xy++){
    bool FIRST=true;
    while(TmpTrackContainer[xy].size()>0){
      double minchi=9999;
      double mintrrms=9999;
      int trnum=-1;
      LocalTrack *tmptr;
      for(int tr=0;tr<(int)TmpTrackContainer[xy].size();tr++){
	tmptr = &(TmpTrackContainer[xy][tr]);
	if(FIRST&&!option.Contains("mwpc")) 
	  if(!tmptr->LeastSquareFit(conf,xy)) continue;	
	double tmpchi=tmptr->chi2(xy);
	double tmptrrms=tmptr->GetTrackTimeRMS();
#define SELECTCHECK 0
#if SELECTCHECK	
	double tmptime=tmptr->GetTrackTime();
	for(int i=0;i<tmptr->nhit();i++)
	  std::cout<<tmptr->hit(i)->wire()<<"\t";
	std::cout<<tmpchi<<"\t";
	std::cout<<tmptime<<"\t";
	std::cout<<tmptrrms<<"\t";
#endif
	if( (tmpchi < minchi &&  (tmptrrms < mintrrms || tmptrrms <10 || mintrrms > 10 || option.Contains("notiming") ) )
	    || ( tmptrrms < mintrrms && (tmpchi < minchi + 5 ) ) ){
#if SELECTCHECK
	  std::cout<<"!!!";
#endif
	  minchi =  tmptr->chi2(xy);
	  mintrrms = tmptr->GetTrackTimeRMS();
	  trnum=tr;
	}
#if SELECTCHECK
	std::cout<<std::endl;
#endif
      }//tr
      if(FIRST)	FIRST=false;
      if(trnum==-1){
	TmpTrackContainer[xy].clear();
	continue;
      }
      tmptr = &(TmpTrackContainer[xy][trnum]);
      //      std::cout<<"xy:"<<xy<<std::endl;
      TmpTrackContainer2[xy].push_back(*tmptr);
      tmptr = &(TmpTrackContainer2[xy].back());
#if DEBUG
      std::cout<<"\tbest track\t";
      for(int i=0;i<tmptr->nhit();i++)
	std::cout<<tmptr->hit(i)->wire()<<"\t";
      std::cout<<minchi<<"\t";
      std::cout<<tmptr->GetTrackTime()<<"\t";
      std::cout<<mintrrms<<"\t";
      std::cout<<std::endl;
#endif 
      BeamLineTrackContainer::iterator it=TmpTrackContainer[xy].begin();
      while(it!=TmpTrackContainer[xy].end()){
	LocalTrack* tmptr2 = &(*it);
#if ERASE
	for(int i=0;i<tmptr2->nhit();i++)
	  std::cout<<tmptr2->hit(i)->wire()<<"\t";
#endif
	if( tmptr->CompareTrackHit(tmptr2) ){
	  TmpTrackContainer[xy].erase(it);
#if ERASE
	  std::cout<<"!!!";
#endif
	}
	else
	  it++;
#if ERASE
	std::cout<<std::endl;
#endif
      }
    } 
  }
  // Set a track candidate
#if DEBUG
  if(id==CID_FDC1)
  for(int xy=0;xy<nxy;xy++)
    for(int tr1=0;tr1<(int)TmpTrackContainer2[xy].size();tr1++){
      LocalTrack *tr = &(TmpTrackContainer2[xy][tr1]);
      for(int i=0;i<tr->nhit();i++){
	ChamberLikeHit *hit=tr->hit(i);
	std::cout<<"cid,layer,wire,dt,lr= "
		 << hit->cid() << "\t"
		 << hit->layer() << "\t"
		 << hit->wire() << "\t"
		 << hit->dt() << "\t"
		 << hit->leftright() << "\t"
		 <<std::endl;
      }
    }
#endif

  int ntrack=0;

  if(nxy==2){
    for(int tr1=0;tr1<(int)TmpTrackContainer2[0].size();tr1++)
      {
	for(int tr2=0;tr2<(int)TmpTrackContainer2[1].size();tr2++)
	  {	
	    LocalTrack *trcand = new LocalTrack();
	    LocalTrack *tmptr[2];
	    tmptr[0] = &(TmpTrackContainer2[0][tr1]);
	    tmptr[1] = &(TmpTrackContainer2[1][tr2]);
	    for(int xy=0;xy<2;xy++){
	      for( int i=0; i<tmptr[xy]->ncluster(xy); i++ ){
		BLDCCluster *clu = tmptr[xy]->cluster(xy,i);
		trcand->SetCluster( xy, clu );
	      }
	      int dof; double chi;
	      chi = tmptr[xy]->chi2(xy);
	      dof = tmptr[xy]->dof(xy);
	      trcand->SetChisqr( xy, chi );
	      trcand->SetDof( xy, dof );
	    }
	    trcand->CalcTrackTime();

	    if( !option.Contains("notiming") &&
		trcand->GetTrackTimeRMS() > 15 ) {
	      delete trcand;	  continue;
	    }	
	    double a,b,c,d,e,f;
	    tmptr[0]->abc(a,b,c);
	    tmptr[1]->def(d,e,f);
	    trcand->SetABC(a,b,c);
	    trcand->SetDEF(d,e,f);
	    
	    TmpTrackContainer3.push_back(*trcand);
	    delete trcand;
	  } //tr2
      }//tr 1
    while(TmpTrackContainer3.size()>0){
      double minchi=9999;
      double mintrrms=9999;
      int besttr=-1;
      for(int itr=0;itr<(int)TmpTrackContainer3.size();itr++){
	LocalTrack* tmptr=&(TmpTrackContainer3[itr]);
	double tmpchi=tmptr->chi2all();
	double tmptrrms=tmptr->GetTrackTimeRMS();
#if 0
	for(int i=0;i<tmptr->nhit();i++)
	  std::cout<<tmptr->hit(i)->wire()<<"\t";
	std::cout<<std::endl<<"\t";
	std::cout<<tmpchi<<"\t";
	std::cout<<tmptrrms<<"\t";
	std::cout<<tmptr->GetTrackTime()<<"\t";
	std::cout<<minchi<<"\t";
	std::cout<<std::endl;
#endif 	
	if(tmpchi< minchi && 
	   (option.Contains("notiming")||tmptrrms<mintrrms) ){
	  minchi=tmpchi;
	  mintrrms=tmptrrms;
	  besttr=itr;
	}
      }
      if(besttr==-1){
	TmpTrackContainer3.clear();
	//	std::cout<<"no track candidate!!!"<<std::endl;
      }
      else{
	LocalTrack* trcand=&(TmpTrackContainer3[besttr]);
	// if(option.Contains("linear")){
	//   for(int xy=0;xy<2;xy++){
	//     if( !trcand->LeastSquareFit( conf, xy, "linear" ) ) {
	//       std::cout<<"LocalTrack::LinearFit() failed for xy="<<xy<<std::endl;
	//     }	
	//   }
	// }
	//	    if( trcand->chi2all() > BLDCParam->GetMaxChi() ) {
	//	      delete trcand;   continue;
	//	    }	
	//CHI2=true;
	trcand->CalcHitPosition();
	trcand->CalcResidual(false,option);	  	  
	SetTrack( id, (*trcand) );
	ntrack++;
	LocalTrack *tmptr = new LocalTrack(*trcand);
#if 0
	std::cout<<"-----------"<<std::endl;
	for(int i=0;i<tmptr->nhit();i++)
	  std::cout<<tmptr->hit(i)->wire()<<"\t";
	std::cout<<std::endl<<"\t";
	std::cout<<tmptr->chi2all()<<"\t";
	std::cout<<tmptr->GetTrackTime()<<"\t";
	std::cout<<tmptr->GetTrackTimeRMS()<<"\t";
	std::cout<<std::endl;
	std::cout<<"-----------"<<std::endl;
#endif 	

	BeamLineTrackContainer::iterator it=TmpTrackContainer3.begin();
	while(it!=TmpTrackContainer3.end()){
#if 0
	  std::cout<<"+++++"<<std::endl;
	  for(int i=0;i<tmptr->nhit();i++)
	    std::cout<<tmptr->hit(i)->wire()<<"\t";
	  std::cout<<std::endl<<"\t";
	  std::cout<<tmptr->chi2all()<<"\t";
	  std::cout<<tmptr->GetTrackTime()<<"\t";
	  std::cout<<tmptr->GetTrackTimeRMS()<<"\t";
	  std::cout<<std::endl;
	  std::cout<<"-----------"<<std::endl;
#endif 	
	  LocalTrack* tmptr2 = &(*it);
	  if( tmptr->CompareTrackHit(tmptr2) ){
	    TmpTrackContainer3.erase(it);
	  }
	  else
	    it++;
	}
	delete tmptr;
      }
    } // container size	  
  }else if(nxy==1){
    for(int tr1=0;tr1<(int)TmpTrackContainer2[0].size();tr1++)
      {	
	LocalTrack *tmptr = &(TmpTrackContainer2[0][tr1]);
	if(tmptr->LinearFit(conf,true)){
	  tmptr->CalcHitPosition(true);
	  tmptr->CalcResidual(true);
	  tmptr->ConvLocalToGlobal2();
	  SetTrack( id, (*tmptr) );
	}
      }
  }
  TmpTrackContainer[0].clear();
  TmpTrackContainer[1].clear();
  TmpTrackContainer2[0].clear();
  TmpTrackContainer2[1].clear();
  TmpTrackContainer3.clear();
  ConvLocalToGlobal(id);
  if(option.Contains("linear"))
    ConvertLocalToLinear(conf,id);
#if DEBUG
  std::cout << "!!! BeamLineTrackMan::Linear Tracking finished" << std::endl;
#endif 
  return SetStatus(id,1);
}

bool BeamLineTrackMan::ConvertLocalToLinear( ConfMan *conf, const int &cid ){
  for(int i=0; i<ntrackBLDC(cid);i++){
    LocalTrack *track=trackBLDC(cid,i);
    if(track->LinearFit(conf,false)){
      track->CalcHitPosition(true);
      track->CalcResidual(true);
      track->ConvLocalToGlobal2();
    }	      
  }
  return true;
}

void BeamLineTrackMan::SetTrack( const int &cid, LocalTrack track)
{
  BeamLineTrackContainer* container=BLDCTrackContainer(cid);
  if(container)
    container->push_back(track);
}

int BeamLineTrackMan::ntrackBLDC(const int &cid,const double &ll,const double &ul)
{
  int n=0;
  BeamLineTrackContainer* container=BLDCTrackContainer(cid);
  if(container)
    for(int i=0;i<(int)container->size();i++)
      if((*container)[i].CheckRange(ll,ul)) n++;
  return n;
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
  case CID_FDC1:
    return (0<=i&&i<FDC1TrackContainer.size())  ? &FDC1TrackContainer[i] : 0;
  default: return 0;
  }
}

bool BeamLineTrackMan::Clustering( BeamLineHitMan *blMan, BLDCClusterMan *clMan, ConfMan *conf, const int &id, const int &id2 )
{
#if 0
  std::cout << "!!! BeamLineTrackMan::Clustering()" << std::endl;
#endif 
  BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();    
  int ncl=4;
  double maxsub=9999.;

  if(id==CID_BPC )    maxsub=120.;
  else if(id==CID_BLC1a || id==CID_BLC1b )    maxsub=160.;
  else if(id==CID_BLC2a || id==CID_BLC2b )    maxsub=80.;
  else if(id==CID_FDC1 ){
    ncl=3;
    maxsub=120;
  }
  int clusterID=0;
  int ixy[2]={0,0};
  for(int icl=0;icl<ncl;icl++){
    int xy=-1;
    int layer1=icl*2+1;
    int layer2=icl*2+2;
    if( BLDCParam->layerdead(id,layer1) && BLDCParam->layerdead(id,layer2))
      continue;
    else if( BLDCParam->layerdead(id,layer1) || blMan->nBLDC(id,layer1)==0 ){
      for( int i=0; i<blMan->nBLDC(id,layer2); i++ ){
	BLDCCluster *cluster = new BLDCCluster();
	ChamberLikeHit *hit = blMan->BLDC(id,layer2, i);      
	xy=hit->xy();
	cluster->SetHit( (*hit) );
	clusterID++;
	cluster->SetClusterID(clusterID);
	clMan->SetCluster( ixy[xy], xy, (*cluster) );
	delete cluster;
      }
      continue;
    }else{
      for( int i1=0; i1<blMan->nBLDC(id,layer1); i1++ ){
	bool PAIRING=false;
	ChamberLikeHit *hit1 = blMan->BLDC(id,layer1, i1); 	
#if 0
	std::cout<<"cid,layer,wire,xy\t"
		 <<hit1->cid()<<"\t"
		 <<hit1->layer()<<"\t"
		 <<hit1->wire()<<"\t"
		 <<hit1->xy()<<"\t"
		 <<std::endl;
#endif
	xy=hit1->xy();
	double dxy=hit1->dxy();
	if( !BLDCParam->layerdead(id,layer2) ){
	  for( int i2=0; i2<blMan->nBLDC(id,layer2); i2++ ){
	    ChamberLikeHit *hit2 = blMan->BLDC(id,layer2, i2); 	
	    double pos1=-999,pos2=999;
	    if(xy==0){
	      pos1=hit1->wx();
	      pos2=hit2->wx();
	    }else if(xy==1){
	      pos1=hit1->wy();
	      pos2=hit2->wy();
	    }
	    double tsub=TMath::Abs((hit1->dt()-hit2->dt())/2.);
	    if(TMath::Abs(pos1-pos2)<dxy/1.5&&tsub<maxsub){
	      BLDCCluster *cluster = new BLDCCluster();
	      cluster->SetHit( (*hit1) );
	      cluster->SetHit( (*hit2) );
	      clusterID++;
	      cluster->SetClusterID(clusterID);
	      clMan->SetCluster( ixy[xy], xy, (*cluster) );
	      delete cluster;
	      PAIRING=true;
	    }
	  }
	}
	if(!PAIRING){
	  BLDCCluster *cluster = new BLDCCluster();
	  cluster->SetHit( (*hit1) );
	  clusterID++;
	  cluster->SetClusterID(clusterID);
	  clMan->SetCluster( ixy[xy],xy, (*cluster) );
	  delete cluster;	  
	}//paring
      }//i1
      for( int i2=0; i2<blMan->nBLDC(id,layer2); i2++ ){
	bool PAIRING=false;
	ChamberLikeHit *hit2 = blMan->BLDC(id,layer2, i2); 	
	xy=hit2->xy();
	double dxy=hit2->dxy();
	for( int i1=0; i1<blMan->nBLDC(id,layer1); i1++ ){
	  ChamberLikeHit *hit1 = blMan->BLDC(id,layer1, i1); 	
	  double pos1=-999,pos2=999;
	  if(xy==0){
	    pos1=hit1->wx();
	    pos2=hit2->wx();
	  }else if(xy==1){
	    pos1=hit1->wy();
	    pos2=hit2->wy();
	  }
	  double tsub=TMath::Abs((hit1->dt()-hit2->dt())/2.);
	  if(TMath::Abs(pos1-pos2)<dxy/1.5&&tsub<maxsub){
	    PAIRING=true;
	  } 
	}//i1
	if(!PAIRING){
	  BLDCCluster *cluster = new BLDCCluster();
	  cluster->SetHit( (*hit2) );
	  clusterID++;
	  cluster->SetClusterID(clusterID);
	  clMan->SetCluster( ixy[xy], xy, (*cluster) );
	  delete cluster;
	}//pair
      }//i2
    }//2hit	   
    if(xy!=-1) ixy[xy]++;
  }//icl
  clMan->Calc(conf);

  return true;
}
//####################################################
void BeamLineTrackMan::ConvLocalToGlobal(const int &cid)
{  
  for( int i=0; i<ntrackBLDC(cid); i++ ){
    LocalTrack *track = trackBLDC(cid,i);
    track->ConvLocalToGlobal();
  }  
}
