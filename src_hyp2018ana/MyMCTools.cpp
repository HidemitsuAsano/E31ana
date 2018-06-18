#include "MyMCTools.h"

using namespace std;

bool MyMCTools::effNC(DetectorData *detData)
{
  vector<DetectorHit*> NChit=NChits(detData);
  bool adc_flag=false;
  bool n_hit=false;
  for( int i=0; i<(int)NChit.size(); i++ ){
    if( NChit[i]->adc()>NC_THRE ) adc_flag=true;
    if( NChit[i]->pdg()==2112 && NChit[i]->parentID()==0 ) n_hit=true;
  }
  if( adc_flag && n_hit ) return true;
  else return false;
}

vector<DetectorHit*> MyMCTools::NChits(DetectorData* detData)
{
  vector<DetectorHit*> hits;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(i);
    if( hit->detectorID()==CID_NC ) hits.push_back(hit);
  }
  return hits;
}

bool MyMCTools::goodBeam(MCData *mcData, DetectorData *detData)
{
  Track *beam=initTrack(mcData, 321);
  bool T0_flag=false;
  bool BPC_flag[8], BLC2a_flag[8], BLC2b_flag[8];
  for( int i=0; i<8; i++ ) BPC_flag[i]=false, BLC2a_flag[i]=false, BLC2b_flag[i]=false;

  for( int i=0; i<beam->detectorHitLinkSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(beam->detectorHitLink(i));
    if( hit->detectorID()==CID_T0 ) T0_flag=true;
    if( hit->detectorID()==CID_BPC ) BPC_flag[hit->layerID()]=true;
    if( hit->detectorID()==CID_BLC2a ) BLC2a_flag[hit->layerID()]=true;
    if( hit->detectorID()==CID_BLC2b ) BLC2b_flag[hit->layerID()]=true;
  }

  int nBPC=0, nBLC2a=0, nBLC2b=0;
  for( int i=0; i<8; i++ ){
    if( BPC_flag[i] ) nBPC++;
    if( BLC2a_flag[i] ) nBLC2a++;
    if( BLC2b_flag[i] ) nBLC2b++;
  }

  if( T0_flag && nBPC==8 && nBLC2a==8 && nBLC2b==8 ) return true;

  // cout<<"T0 flag : "<<boolalpha<<T0_flag<<endl;
  // cout<<"nBPC    : "<<nBPC<<endl;
  // cout<<"nBLC2a  : "<<nBLC2a<<endl;
  // cout<<"nBLC2b  : "<<nBLC2b<<endl;
  return false;
}

bool MyMCTools::goodPi_beam(MCData *mcData, DetectorData *detData)
{
  Track *beam=initTrack(mcData, -211);
  if( !beam ) beam=initTrack(mcData, 211);
  bool T0_flag=false;
  bool BPC_flag[8], BLC2a_flag[8], BLC2b_flag[8];
  for( int i=0; i<8; i++ ) BPC_flag[i]=false, BLC2a_flag[i]=false, BLC2b_flag[i]=false;

  for( int i=0; i<beam->detectorHitLinkSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(beam->detectorHitLink(i));
    if( hit->detectorID()==CID_T0 ) T0_flag=true;
    if( hit->detectorID()==CID_BPC ) BPC_flag[hit->layerID()]=true;
    if( hit->detectorID()==CID_BLC2a ) BLC2a_flag[hit->layerID()]=true;
    if( hit->detectorID()==CID_BLC2b ) BLC2b_flag[hit->layerID()]=true;
  }

  int nBPC=0, nBLC2a=0, nBLC2b=0;
  for( int i=0; i<8; i++ ){
    if( BPC_flag[i] ) nBPC++;
    if( BLC2a_flag[i] ) nBLC2a++;
    if( BLC2b_flag[i] ) nBLC2b++;
  }

  if( T0_flag && nBPC==8 && nBLC2a==8 && nBLC2b==8 ) return true;

  // cout<<"T0 flag : "<<boolalpha<<T0_flag<<endl;
  // cout<<"nBPC    : "<<nBPC<<endl;
  // cout<<"nBLC2a  : "<<nBLC2a<<endl;
  // cout<<"nBLC2b  : "<<nBLC2b<<endl;
  return false;
}

bool MyMCTools::goodP_beam(MCData *mcData, DetectorData *detData)
{
  Track *beam=initTrack(mcData, 2212);
  bool T0_flag=false;
  bool BPC_flag[8], BLC2a_flag[8], BLC2b_flag[8];
  for( int i=0; i<8; i++ ) BPC_flag[i]=false, BLC2a_flag[i]=false, BLC2b_flag[i]=false;

  for( int i=0; i<beam->detectorHitLinkSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(beam->detectorHitLink(i));
    if( hit->detectorID()==CID_T0 ) T0_flag=true;
    if( hit->detectorID()==CID_BPC ) BPC_flag[hit->layerID()]=true;
    if( hit->detectorID()==CID_BLC2a ) BLC2a_flag[hit->layerID()]=true;
    if( hit->detectorID()==CID_BLC2b ) BLC2b_flag[hit->layerID()]=true;
  }

  int nBPC=0, nBLC2a=0, nBLC2b=0;
  for( int i=0; i<8; i++ ){
    if( BPC_flag[i] ) nBPC++;
    if( BLC2a_flag[i] ) nBLC2a++;
    if( BLC2b_flag[i] ) nBLC2b++;
  }

  if( T0_flag && nBPC==8 && nBLC2a==8 && nBLC2b==8 ) return true;

  // cout<<"T0 flag : "<<boolalpha<<T0_flag<<endl;
  // cout<<"nBPC    : "<<nBPC<<endl;
  // cout<<"nBLC2a  : "<<nBLC2a<<endl;
  // cout<<"nBLC2b  : "<<nBLC2b<<endl;
  return false;
}

Track* MyMCTools::initTrack(MCData *mcData, int pdg)
{
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track=mcData->track(i);
    if( track->parentTrackID()==0 && track->pdgID()==pdg ) return track;
  }
  cout<<"!!! trackSize="<<mcData->trackSize()<<endl;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track=mcData->track(i);
    cout<<"    PDG : "<<track->pdgID()<<endl;
  }

  return 0;
}

Track* MyMCTools::track(MCData *mcData, int id)
{
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track=mcData->track(i);
    if( track->trackID()==id ) return track;
  }
  cout<<" !!! trackSize="<<mcData->trackSize()<<"  track id="<<id<<" !!!"<<std::endl;
 
  return 0;
}

Track* MyMCTools::trackByPDG(MCData *mcData, int id)
{
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track=mcData->track(i);
    if( track->pdgID()==id ) return track;
  }
  // cout<<" !!! trackSize="<<mcData->trackSize()<<"  track id="<<id<<" !!!"<<std::endl;
  // for( int i=0; i<mcData->trackSize(); i++ ){
  //   Track *track=mcData->track(i);
  //   cout<<"    PDG : "<<track->pdgID()<<endl;
  // }

  return 0;
}

double MyMCTools::lmom(int pdg, ReactionData *reacData)
{
  for( int i=0; i<reacData->ParticleSize(); i++ ){
    if( reacData->PDG(i)==pdg ) return 0.001*reacData->GetParticle(i).M();
  }
  cout<<"  !!! MyMCTools::lmom("<<pdg<<") not found !!!"<<endl;
  return 0;
}

int MyMCTools::nHit(int cid, DetectorData *detData)
{
  int nhit=0;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *ddata=detData->detectorHit(i);
    if( ddata->detectorID()==cid ) nhit++;
  }
  return nhit;
}

bool MyMCTools::isHit(int cid, Track *track, DetectorData *detData)
{
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *ddata=detData->detectorHit(i);
    if( ddata->detectorID()==cid && ddata->trackID()==track->trackID() ) return true;
  }
  return false;
}

vector<Track*> MyMCTools::getInitTracks(MCData *mcData)
{
  vector<Track*> init;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track = mcData->track(i);
    if( track->parentTrackID()==0 ){
      init.push_back(track);
    }
  }
  return init;
}

vector<Track*> MyMCTools::getDaughter(Track *track, MCData* mcData)
{
  vector<Track*> daughter;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *tmp = mcData->track(i);
    if( tmp->parentTrackID()==track->trackID() ) daughter.push_back(tmp);
  }
  return daughter;
}

CDSInfo* MyMCTools::getCDSInfo(DetectorHit *hit, AnaInfo *anaInfo, CDSHitMan *cdsMan)
{
  if( hit->detectorID()!=CID_CDH ) return 0;

  for( int i=0; i<anaInfo->nCDS(); i++ ){
    CDSInfo *cds=anaInfo->CDS(i);
    HodoscopeLikeHit *CDHhit=cds->CDH(cdsMan);
    //    cout<<" seg : "<<CDHhit->seg()<<endl;
    TVector3 pos=CDHhit->pos();
    //    cout<<Form("CDH pos(%lf, %lf, %lf)", pos.X(), pos.Y(), pos.Z())<<endl;
    if( CDHhit->seg()==hit->channelID()+1 ) return cds;
  }

  cout<<anaInfo->nCDS()<<endl;
  for( int i=0; i<anaInfo->nCDS(); i++ ){
    CDSInfo *cds=anaInfo->CDS(i);
    HodoscopeLikeHit *CDHhit=cds->CDH(cdsMan);
    TVector3 pos=CDHhit->pos();
    cout<<" seg : "<<CDHhit->seg()<<endl;
    cout<<Form("CDH pos(%lf, %lf, %lf)", pos.X(), pos.Y(), pos.Z())<<endl;
    if( CDHhit->seg()==hit->channelID()+1 ) return cds;
  }

  return 0;
}

DetectorHit* MyMCTools::getHit(DetectorData *detData, HodoscopeLikeHit *hit)
{
  int cid=hit->cid();
  int seg=hit->seg();
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *dhit=detData->detectorHit(i);
    if( dhit->detectorID()==cid && dhit->channelID()+1==seg ) return dhit;
  }
  return 0;
}

bool MyMCTools::effFC(DetectorData *detData, MCData *mcData, BeamLineHitMan *blMan)
{
  int nPC=0;
  int nCVC=0;
  int nFDC1=0;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(i);
    if( hit->detectorID()==CID_CVC ) nCVC++;
    if( hit->detectorID()==CID_PC  ) nPC++;
    if( hit->detectorID()==CID_FDC1 ) nFDC1++;
  }

  bool fc_flag=false;
  if( nPC==1 && nCVC==0 ) fc_flag=true;
  if( nCVC==1 && nPC==0 ) fc_flag=true;

  if( nFDC1<6 ) fc_flag=false;
  if( !fc_flag ) return false;

  Track *ptrack=0;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *tmp=mcData->track(i);
    if( tmp->parentTrackID()==0 && tmp->pdgID()==2212 ){
      ptrack=tmp;
    }
  }

  bool hit_flag=false;
  int nFDC1mc=0;
  for( int i=0; i<ptrack->detectorHitLinkSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(ptrack->detectorHitLink(i));
    if( hit->detectorID()==CID_CVC || hit->detectorID()==CID_PC ) hit_flag=true;
    if( hit->detectorID()==CID_FDC1 ) nFDC1mc++;
  }

  if( fc_flag && hit_flag && nFDC1mc>=6 ) return true;
  return false;
}

bool MyMCTools::pimS0all_hit(MCData *mcData, DetectorData *detData)
{
  Track *ptrack=0;
  Track* Strack=0;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *tmp=mcData->track(i);
    if( tmp->parentTrackID()==0 && tmp->pdgID()==2212 ) ptrack=tmp;
    if( tmp->parentTrackID()==0 && tmp->pdgID()==3114 ) Strack=tmp;
  }
  Track *Strack2=0;
  Track *pim_track=0;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *tmp=mcData->track(i);
    if( tmp->parentTrackID()==Strack->trackID() && tmp->pdgID()==3212 ) Strack2=tmp;
    if( tmp->parentTrackID()==Strack->trackID() && tmp->pdgID()==-211 ) pim_track=tmp;
  }
  Track *Ltrack=0;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *tmp=mcData->track(i);
    if( tmp->parentTrackID()==Strack2->trackID() && tmp->pdgID()==3122 ) Ltrack=tmp;
  }
  Track *son_pi=0;
  Track *son_N=0;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *tmp=mcData->track(i);
    if( tmp->parentTrackID()==Ltrack->trackID() && tmp->pdgID()<1000 ) son_pi=tmp;
    if( tmp->parentTrackID()==Ltrack->trackID() && tmp->pdgID()<1000 ) son_N=tmp;
  }

  bool fc_flag=false;
  bool hit_flag=false;
  int nFDChit=0;
  for( int i=0; i<ptrack->detectorHitLinkSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(ptrack->detectorHitLink(i));
    if( hit->detectorID()==CID_CVC || hit->detectorID()==CID_PC ) hit_flag=true;
    if( hit->detectorID()==CID_FDC1 ) nFDChit++;
  }
  if( hit_flag && nFDChit==6 ) fc_flag=true;
  if( fc_flag ){
    if( son_pi && son_pi->pdgID()==-211 ){
      bool CDHhit1=false;
      bool CDHhit2=false;
      int nCDC1=0;
      int nCDC2=0;
      for( int i=0; i<pim_track->detectorHitLinkSize(); i++ ){
	DetectorHit *hit=detData->detectorHit(pim_track->detectorHitLink(i));
	if( hit->detectorID()==CID_CDH ) CDHhit1=true;
	if( hit->detectorID()==CID_CDC ) nCDC1++;
      }

      for( int i=0; i<son_pi->detectorHitLinkSize(); i++ ){
	DetectorHit *hit=detData->detectorHit(son_pi->detectorHitLink(i));
	if( hit->detectorID()==CID_CDH ) CDHhit2=true;
	if( hit->detectorID()==CID_CDC ) nCDC2++;
      }
      if( CDHhit1 && CDHhit2 ) return true;
    }
  }
  return false;
}

double MyMCTools::Ystar_mass(ReactionData *reacData)
{
  TLorentzVector lmom;
  for( int i=0; i<reacData->ParticleSize(); i++ ){
    if( reacData->PDG(i)==13122 || reacData->PDG(i)==3224 || reacData->PDG(i)==3214 || reacData->PDG(i)==3114 ){
      lmom=reacData->GetParticle(i);
    }
  }
  //  cout<<"LMOM : "<<lmom.M()<<endl;
  return 0.001*lmom.M();
}

double MyMCTools::S1385mass(ReactionData *reacData)
{
  TLorentzVector lmom;
  for( int i=0; i<reacData->ParticleSize(); i++ ){
    if( reacData->PDG(i)==3224 || reacData->PDG(i)==3214 || reacData->PDG(i)==3114 ){
      lmom=reacData->GetParticle(i);
    }
  }
  //  cout<<"LMOM : "<<lmom.M()<<endl;
  return 0.001*lmom.M();
}
