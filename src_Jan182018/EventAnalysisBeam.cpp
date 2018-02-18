#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#include "time.h"

#define TREE 0

class EventAnalysisBeam: public EventTemp
{
public:
  EventAnalysisBeam();
  ~EventAnalysisBeam();

private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  BeamLineHitMan *blMan;
  CDSHitMan *cdsMan;
  BeamLineTrackMan *trackMan;
  BeamSpectrometer *beam;
  EventHeader *header;
  ScalerMan *scaMan;
  int t0,t1;  

  void HistSemiLocal();
  void FillSemiLocal( const int &ic,LocalTrack *track);

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  void UTime( int time ){};
  bool UAna( TKOHitCollection *tko );
  void Finalize();

  void InitializeHistogram();
};

EventAnalysisBeam::EventAnalysisBeam()
  : EventTemp()
{
}

EventAnalysisBeam::~EventAnalysisBeam()
{
}

//const int MaxTreeSize = 1900000000000;
void EventAnalysisBeam::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisBeam::Initialize " << std::endl;
#endif

  confMan = conf;
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();
  evTree=new TTree( "EventTree", "EventTree" );  
  beam = new BeamSpectrometer(conf);
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "BeamLineHitMan", &blMan );
  trackMan = new BeamLineTrackMan();
  if( trackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "BeamLineTrackMan", &trackMan );
  cdsMan = new CDSHitMan();
  if( cdsMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSHitMan", &cdsMan );
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  t0=clock();
}

void EventAnalysisBeam::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisHodo::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool EventAnalysisBeam::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisBeam::UAna " << std::endl;
#endif
   
  Event_Number++;
  if(1){ 
    int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; 
  }
  t1=clock();
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s): " << (t1-t0)*1.0e-6 << std::endl;
    
  blMan->Convert( tko, confMan );
  header->Convert( tko, confMan );
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  rtFile->cd();  
  TH1F *h1;
  //  TH2F* h2;

  TString tmp[2]={"BLC1","BLC2"};
  int tmpid[2]={CID_BLC1,CID_BLC2};
  trackMan->DoTracking( blMan, confMan );
  trackMan->LinearTracking( blMan, confMan, CID_BLC1 );
  trackMan->LinearTracking( blMan, confMan, CID_BLC2 );
  trackMan->ConvertLocalToLinear(CID_BLC1);
  trackMan->ConvertLocalToLinear(CID_BLC2);
  if(trackMan->nltrackBLDC(tmpid[0])==1&&trackMan->nltrackBLDC(tmpid[1])==1){
    beam->TMinuitFit(trackMan->ltrackBLC1(0),trackMan->ltrackBLC2(0),confMan);
    h1=(TH1F*)gFile->Get("Chi2D5");
    h1->Fill(beam->chisquare());    
    h1=(TH1F*)gFile->Get("Mom");
    h1->Fill(beam->mom());
  }

  int ntr;
  for(int i=0;i<2;i++){
    ntr=trackMan->ntrackBLDC(tmpid[i]);
    h1=(TH1F*)gFile->Get(Form("nTrack%s",tmp[i].Data())); h1->Fill(ntr);  
    for(int itr=0;itr<ntr;itr++){
      LocalTrack *track=trackMan->trackBLDC(tmpid[i],itr);
      double chi2=track->chi2all();
      h1=(TH1F*)gFile->Get(Form("Chi2%s",tmp[i].Data())); h1->Fill(chi2);
      if(chi2>20) continue;
      FillSemiLocal(i,track);
    }
  }

#if TREE
  evTree->Fill();
#endif
  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  trackMan->Clear();
  beam->Clear();
  return true;
}

void EventAnalysisBeam::Finalize()
{
  std::cout << " Enter EventAnalysisBeam::Finalize " << std::endl;
  rtFile->cd();
  confMan->SaveCode();
  confMan->SaveParams();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete header;
  delete cdsMan;
  delete trackMan;
  delete beam;
}


void EventAnalysisBeam::FillSemiLocal(const int &c,LocalTrack* track)
{
  TH1F* h1;
  TH2F* h2;
  double a1,b1,c1,d1,e1,f1;
  double cid,dt,dltrack,dl;
  double resi=-99.;
  int wire,layer,lr;
  TString tmp[2]={"BLC1","BLC2"};
  int fac[2]={1,-1};
  TString strlr[2]={"L","R"};    
  double chi2=track->chi2all();  
  track->gabc(a1,b1,c1);
  track->gdef(d1,e1,f1);
  h2=(TH2F*)gFile->Get(Form("AB%s",tmp[c].Data())); h2->Fill(b1,e1);
  double tempx,tempy;	
  if(c==0)
    track->XYPosatZ(0.,tempx,tempy);//cm
  else
    track->XYPosatZ(-130.,tempx,tempy);//cm
  h2=(TH2F*)gFile->Get(Form("XY%s",tmp[c].Data())); h2->Fill(tempx,tempy);

  for( int ih=0;ih<track->nhit();ih++ ){
    ChamberLikeHit *hit=track->hit(ih);
    cid=hit->cid();
    layer=hit->layer();
    if(cid==16||cid==18)
      layer+=8;
    wire=hit->wire();
    dl=hit->dl();
    dt=hit->dt();
    if(hit->xy()) dltrack=hit->y()-hit->wy();
    else          dltrack=hit->x()-hit->wx();
    resi=hit->resl();
    lr=hit->leftright();

    h2=(TH2F*)gFile->Get(Form("h_dt_dltrack_%s_%d",tmp[c].Data(),layer)); 
    h2->Fill(dt,dltrack);
    h2=(TH2F*)gFile->Get(Form("h_dl_resi_%s_%d",tmp[c].Data(),layer)); 
    h2->Fill(dl*fac[lr],resi);
    h1=(TH1F*)gFile->Get(Form("h_resi_%s_%d",tmp[c].Data(),layer));
    h1->Fill(resi);
    h2=(TH2F*)gFile->Get(Form("h_dt_resi_%s_%d",tmp[c].Data(),layer));
    h2->Fill(dt*fac[lr],resi);
  }
}

void EventAnalysisBeam::HistSemiLocal(){
  TString tmp[2]={"BLC1","BLC2"};
  double dlmax[2]={0.4,0.25};
  for(int iname=0;iname<2;iname++){
    new TH1F( Form("Chi2%s",tmp[iname].Data()),
		Form("Chi square/Dof %s",tmp[iname].Data()),
		1000, 0,100 );
    new TH1F( Form("nTrack%s",tmp[iname].Data()),
		Form("ntrack %s",tmp[iname].Data()),
		20, 0,20 );
    new TH2F( Form("AB%s",tmp[iname].Data()),
		Form("B%s",tmp[iname].Data()),
		100,-0.1,0.1,100,-0.1,0.1 );
    new TH2F( Form("XY%s",tmp[iname].Data()),
		Form("XY%s",tmp[iname].Data()),
		100,-15,15,100,-15,15 );
    for(int i=0;i<16;i++){
      new TH2F(Form("h_dt_dltrack_%s_%d",tmp[iname].Data(),i+1),
	       Form("dt vs dltrack in %s layer%d",tmp[iname].Data(),i+1) ,
	       200,-50,200,200,-dlmax[iname]*1.2,dlmax[iname]*1.2);
      new TH2F(Form("h_dl_resi_%s_%d",tmp[iname].Data(),i+1),
	       Form("dl vs resi in %s layer%d",tmp[iname].Data(),i+1) ,
	       200,-dlmax[iname]*1.2,dlmax[iname]*1.2,200,-0.1,0.1);
      new TH2F(Form("h_dt_resi_%s_%d",tmp[iname].Data(),i+1),
	       Form("dt vs residu[cm] in %s layer%d",tmp[iname].Data(),i+1) ,
	       200,-200,200,200,-0.1,0.1);
      new TH1F(Form("h_resi_%s_%d",tmp[iname].Data(),i+1),
	       Form("residu[cm] in %s layer%d",tmp[iname].Data(),i+1) ,
	       300,-0.1,0.1);
    }
  }
}
void EventAnalysisBeam::InitializeHistogram()
{
  rtFile->cd();  
  HistSemiLocal();
  new TH1F("Chi2D5",  "Chi2D5", 500,0,100);
  new TH1F("Mom",  "Mom", 1000,900,1100);
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisBeam *event = new EventAnalysisBeam();
  return (EventTemp*)event;
}
