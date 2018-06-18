#include <TLorentzVector.h>


#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
#include "BeamSpectrometer.h"
//#include "Particle.h"
#include "GlobalVariables.h"
#include "TransferMatrixMan.h"
#include "HodoClusterMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "Tools.h"
#include "MathTools.h"
#include "TrackTools.h"

class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TFile *cdcFile;
  TTree *cdcTree;
  TTree *scaTree;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackMan;
  BeamLineHitMan *blMan;
  BeamLineTrackMan *bltrackMan;
  EventHeader *header;
  EventHeader *header2;
  int t0,t1;
  
  int AllGoodTrack;
  int nTrack;
  int CDC_Event_Number;
  
public:
  void Initialize( ConfMan *conf );
  int BeginOfAnEvent();
  bool EndOfAnEvent();
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysis::EventAnalysis()
  : EventTemp()
{
}

EventAnalysis::~EventAnalysis()
{
}

void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  TString cdcfname=confMan->GetCDSTrackFileName();
  cdcFile = new TFile(cdcfname);
  if(!cdcFile->IsOpen()){
    std::cout<<" failed to open " <<cdcfname<< "  !!!"<<std::endl;
    exit(false);
  }
  if(cdcFile->IsZombie()){
    std::cout<<cdcfname<< " is Zombie !!!"<<std::endl;
    exit(false);
  }

  trackMan=new CDSTrackingMan();
  header2 =new EventHeader();
  cdcTree=(TTree*)cdcFile->Get("EventTree");
  cdcTree->SetBranchAddress( "CDSTrackingMan", &trackMan );
  cdcTree->SetBranchAddress( "EventHeader" ,&header2);
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  rtFile->cd();

  header     = new EventHeader();
  cdsMan     = new CDSHitMan();
  blMan      = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();

  t0=time(0);
  AllGoodTrack=0;
  nTrack=0;
  CDC_Event_Number=0;
}

void EventAnalysis::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysis::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

int EventAnalysis::BeginOfAnEvent()
{
  header->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();

  int status=0;
  Event_Number++;
  { 
    status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if(status) return status;
  }
  if( Event_Number%5000==0)
    {
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0) << std::endl;
    }
  return status;
}

bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysis::UAna " << std::endl;
#endif

  int status= BeginOfAnEvent();
  if( status==1 ) return true;
  if( status==2 ) return false;

  // check event matching
  if(CDC_Event_Number>=cdcTree->GetEntries()) return false;  
  cdcTree->GetEntry(CDC_Event_Number);
  if(header2->ev()!=Event_Number){
    int tmpnum=CDC_Event_Number;
    while(header2->ev()<Event_Number&&tmpnum<=cdcTree->GetEntries()){
      tmpnum++;    
      cdcTree->GetEntry(tmpnum);
    }
    if(header2->ev()!=Event_Number)  return true;
    else CDC_Event_Number=tmpnum;
  }
  CDC_Event_Number++;

  // check trigger flag
  if( header2->IsTrig(Trig_Cosmic) ) return EndOfAnEvent();
  if(!header2->kaon()) return EndOfAnEvent();

  // check cdc track number
  int nGoodTrack=trackMan->nGoodTrack();
  int nallTrack=  trackMan->nTrack();
  AllGoodTrack+=nGoodTrack;
  nTrack+=nallTrack;
  if( nGoodTrack<1 ) return EndOfAnEvent();
  // convert raw data
  header->Convert( tko, confMan );
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);  
  cdsMan->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  trackMan->Calc( cdsMan, confMan );  

  // check T0 hit
  int nT0=0;
  int segT0=-1;    
  double ctmT0=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
      ctmT0 = blMan->T0(i)->ctmean();
      segT0=blMan->T0(i)->seg();
    }
  }
  if( nT0!=1 )   return EndOfAnEvent();

  //check BHD
  int nBHD=0;
  double ctmBHD=0;
  int segBHD=-1;
  int npid=3;
  int pid_beam=2;//0:pi 1:K 2:else
  double tofBHDT0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      nBHD++;
      ctmBHD = blMan->BHD(i)->ctmean();
      tofBHDT0 = ctmT0-ctmBHD;    
      if(header->kaon()&&tofBHDT0>27.5&&tofBHDT0<31.) pid_beam=Beam_Kaon;
      else if(header->pion()&&tofBHDT0<27.5&&tofBHDT0>24.5) pid_beam=Beam_Pion;
    }
  }
  //  std::cout<<"TOFBHDT0  "<<tofBHDT0<<std::endl;
  Tools::H1("tofBHDT0", tofBHDT0,200,20,40 );  
  if(pid_beam!=Beam_Kaon)  return EndOfAnEvent();

  //require 1 track for BPC
  bltrackMan->Clear();
  bltrackMan->LocalTracking(blMan,confMan,CID_BPC);
  int nbpc=0;
  int bpcid;
  for(int i=0;i<bltrackMan->ntrackBPC();i++)
    if(bltrackMan->trackBPC(i)->CheckRange(-30,100)){
      nbpc++;
      bpcid=i;
    }
  if(nbpc!=1) return EndOfAnEvent();
  LocalTrack *bpctrack=bltrackMan->trackBPC(bpcid);    

  //################################
  //####CDS 1 track analysis #######
  //################################
  
  Tools::H1("nGoodTrackCDC",trackMan->nGoodTrack(),10,-0.5,9.5);
  
  for( int it=0; it<trackMan->nGoodTrack(); it++ ){
    CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
    double param[5];
    track->GetParameters(param);
    double mom = track->Momentum();
    TVector3 vtxb1, vtxb2;
    double tmpdis;
    TrackTools::CalcLineHelixVertex(bpctrack,track,vtxb1,vtxb2,tmpdis);
    track->SetPID(-1);
    if(!track->CDHFlag() ) continue;

    double beta=-1.;
    double tof=999.;
    double mass2=-999.;

    for(int icdh=0;icdh<track->nCDHHit();icdh++){
      HodoscopeLikeHit *cdhhit=track->CDHHit(cdsMan,icdh);
      double tmptof    = cdhhit->ctmean()-ctmT0;      
      if(tmptof<tof||tof==999.){
	tof=tmptof;
      }
    }
    TrackTools::CalcBetaMass(vtxb1,bpctrack,track,confMan,pid_beam,tof,beta,mass2);
    double mass=sqrt(fabs(mass2));
    Tools::H2("mom_overbeta",1./beta,mom,200,0,10,240,-1.2,1.2);
    Tools::H2("mom_tof", tof, mom,200,-5,15,240,-1.2,1.2);
    int pid=TrackTools::PID(mom,mass2);      
    track->SetPID(pid);

    TVector3 cvtxb1,cvtxb2,cPd;
    double tofvtxcdc,tmpl;
    track->GetVertex(bpctrack->GetPosatZ(0),bpctrack->GetMomDir(),cvtxb2,cvtxb1);
    track->CalcVertexTimeLength(bpctrack->GetPosatZ(-110.5),
				bpctrack->GetMomDir(),
				cdsMass[pid],
				cvtxb2,cvtxb1,tofvtxcdc,tmpl,true);
    Tools::H1("VertexZ",cvtxb2.Z(),1000,-50,50);
    Tools::H2("VertexZX",cvtxb2.Z(),cvtxb2.X(),200,-25,25,200,-10,10);
    Tools::H2("VertexZY",cvtxb2.Z(),cvtxb2.Y(),200,-25,25,200,-10,10);
    Tools::H2("VertexXY",cvtxb2.X(),cvtxb2.Y(),200,-10,10,200,-10,10);
  }// cdc itrack


  //################################
  //####CDS 2 track analysis #######
  //################################
  for( int it1=0; it1<trackMan->nGoodTrack(); it1++ ){
    for( int it2=it1+1; it2<trackMan->nGoodTrack(); it2++ ){
      //      std::cout<<"---------------"<<std::endl;
      CDSTrack *track1 = trackMan->Track( trackMan->GoodTrackID(it1) );
      CDSTrack *track2 = trackMan->Track( trackMan->GoodTrackID(it2) );
      int pid1=track1->PID();
      int pid2=track2->PID();
      if(pid1<0||pid2<0) continue;
      double mass1=cdsMass[pid1];
      double mass2=cdsMass[pid2];
      TVector3 vtx1,vtx2;      
      TVector3 Pp1,Pp2;
      if( !TrackTools::Calc2HelixVertex2(track1,track2,vtx1,vtx2) ) continue;
      TVector3 vtx=(vtx1+vtx2)*0.5;
      if( GeomTools::GetID(vtx)!=CID_Fiducial ) continue;
      if( !track1->GetMomentum(vtx1,Pp1) )	continue;
      if( !track2->GetMomentum(vtx2,Pp2) )	continue;
      TVector3 Pp= Pp1+Pp2;
      TLorentzVector L1; L1.SetVectM( Pp1, mass1 );
      TLorentzVector L2; L2.SetVectM( Pp2, mass2 );
      double im = (L1+L2).M();	
      
      TVector3 cvtx1,cvtx2;
      if( !TrackTools::Calc2HelixVertex(track1,track2,cvtx1,cvtx2) )continue;
      TVector3 cPp1,cPp2;
      if( !track1->GetMomentum(cvtx1,cPp1,true,true) ) continue;
      if( !track2->GetMomentum(cvtx2,cPp2,true,true) ) continue;

      TVector3 cPp= cPp1+cPp2;
      TLorentzVector cL1; cL1.SetVectM( cPp1, mass1 );
      TLorentzVector cL2; cL2.SetVectM( cPp2, mass2 );
      double cim = (cL1+cL2).M();	

      int combid=pow(2,pid1)+pow(2,pid2);
      if(combid==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
	Tools::H1(Form("InvPiPi"), im, 1000,0,1 );
	Tools::H1(Form("cInvPiPi"), cim ,1000,0,1);
      }
      if(combid==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
	Tools::H1(Form("InvPiP"), im,1000,1,2 );
	Tools::H1(Form("cInvPiP"), cim,1000,1,2 );
      }
    }
  }	

  return EndOfAnEvent();
}

bool EventAnalysis::EndOfAnEvent(){
  header->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  return true;
}
void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;

  rtFile->cd();
  confMan->SaveCode();
  confMan->SaveParams();
  gFile->Write();
  gFile->Close();
  cdcFile->Close();

  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
