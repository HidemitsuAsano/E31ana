#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamSpectrometer.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#define TKOHIS 0
#define TREE 0
#define SUB 0
#define DEBUG 0
#define TRACK 1
class EventAnalysisMomentum: public EventTemp
{
public:
  EventAnalysisMomentum();
  ~EventAnalysisMomentum();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  BeamLineHitMan *blMan;
  BeamLineTrackMan *trackMan;
  BeamSpectrometer *beam;
  CDSHitMan *cdsMan;
  EventHeader *header;
  ScalerMan *scaMan;

  TString name[5];
  int cid[5];
  int fix[5];
  int nHodo[5];
  int Hodoseg[5];
  double ctmHodo[5];
  double ctsubHodo[5];
  double ctof[5][5];
  bool ELECTRON,PROTON;
  void HistSemiLocal();
  void FillLinear( const int &ic,LinearTrack *track);

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  double t0,t1;
  void InitializeHistogram();
};

EventAnalysisMomentum::EventAnalysisMomentum()
  : EventTemp()
{
}

EventAnalysisMomentum::~EventAnalysisMomentum()
{
}

//const int MaxTreeSize = 190000000000;
void EventAnalysisMomentum::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMomentum::Initialize " << std::endl;
#endif
  confMan = conf;  

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();
  rtFile ->cd();
  t0=clock();
  evTree=new TTree("EventTree","EventTree");
  scaTree=new TTree("ScalerTree","ScalerTree");

  beam = new BeamSpectrometer(conf);
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("EventHeader",&header);
#if CDS
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
#endif
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("BeamLineHitMan",&blMan);
  trackMan = new BeamLineTrackMan();
  if( trackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("BeamLineTrackMan",&trackMan);
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch("ScalerMan",&scaMan);
}

void EventAnalysisMomentum::USca( int nsca, unsigned int *sca )
{
}

bool EventAnalysisMomentum::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisMomentum::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  int t1=clock();
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number<< " Time (s): " << (t1-t0)*1.0e-6 << std::endl;
  
  if(getenv("RUN")!=NULL)  header->SetRunNumber(atoi(getenv("RUN")));
  else header->SetRunNumber(-1);
  header->SetEventNumber(Event_Number);

  TH1F *h1;
  TH2F *h2;
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
#if CDS
  cdsMan->Convert( tko, confMan );
#endif
#if TRACK
  trackMan->DoTracking( blMan, confMan );
  trackMan->LinearTracking( blMan, confMan, CID_BLC1 );
  trackMan->LinearTracking( blMan, confMan, CID_BLC2 );
  trackMan->ConvertLocalToLinear(  CID_BLC1  );
  trackMan->ConvertLocalToLinear(  CID_BLC2  );
#endif
  rtFile->cd();
  if(trackMan->nltrackBLC1()==1&&trackMan->nltrackBLC2()==1){
    beam->TMinuitFit(trackMan->ltrackBLC1(0),trackMan->ltrackBLC2(0),confMan);
    h1=(TH1F*)gFile->Get("Chi2D5");
    h1->Fill(beam->chisquare());    
    h1=(TH1F*)gFile->Get("Mom");
    h1->Fill(beam->mom());
    FillLinear(0,beam->blc1track());
    FillLinear(1,beam->blc2track());
  }
  int tmpid[2]={CID_BLC1,CID_BLC2};
  TString tmp[2]={"BLC1","BLC2"};
  for(int i=0;i<2;i++){
    int ntr=trackMan->nltrackBLDC(tmpid[i]);
    h1=(TH1F*)gFile->Get(Form("nTrack%s",tmp[i].Data())); h1->Fill(ntr);
    if(ntr!=1) continue;      
    for(int itr=0;itr<ntr;itr++){
      LinearTrack *track=trackMan->ltrackBLDC(tmpid[i],itr);
      double chi2=track->chi2all();
      h1=(TH1F*)gFile->Get(Form("Chi2%s",tmp[i].Data())); h1->Fill(chi2);
      if(chi2>20) continue;
      //  FillLinear(i,track);
    }
  }
  
#if TREE
  evTree->Fill();
#endif
  header->Clear();
  blMan->Clear();
  trackMan->Clear();
#if CDS
  cdsMan->Clear();
#endif
  return true;
}

void EventAnalysisMomentum::Finalize()
{
  std::cout << " Enter EventAnalysisMomentum::Finalize " << std::endl;
  
  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
#if CDS
  delete cdsMan;
#endif
  delete header;
  delete trackMan;
}

void EventAnalysisMomentum::FillLinear(const int &c,LinearTrack* track)
{
  TH1F* h1;
  TH2F* h2;
  double a1,b1,c1,d1;
  double cid,dt,dltrack,dl;
  double resi=-99.;
  int wire,layer,lr;
  TString tmp[2]={"BLC1","BLC2"};
  int fac[2]={1,-1};
  TString strlr[2]={"L","R"};    
  double chi2=track->chi2all();  
  track->gabcd(a1,b1,c1,d1);
  double theta=TMath::ATan(b1)*1000;
  double phi=TMath::ATan(d1)*1000;
  h2=(TH2F*)gFile->Get(Form("AB%s",tmp[c].Data())); h2->Fill(theta,phi);
  double tempx,tempy;	
  if(c==0)
    track->XYPosatZ(0.,tempx,tempy);//cm
  else
    track->XYPosatZ(-130.,tempx,tempy);//cm
  h2=(TH2F*)gFile->Get(Form("XY%s",tmp[c].Data())); h2->Fill(tempx,tempy);
  h2=(TH2F*)gFile->Get(Form("XA%s",tmp[c].Data())); h2->Fill(tempx,theta);
  h2=(TH2F*)gFile->Get(Form("YB%s",tmp[c].Data())); h2->Fill(tempy,phi);
  /**
  if(c==1){
    if(nT0==1){
      track->XYPosatZ(-100.,tempx,tempy);//cm
      //      std::cout<<"T0seg, x, y by BLC  "<<t0seg<<"\t"<<tempx<<"\t"<<tempy<<std::endl;
      h2=(TH2F*)gFile->Get(Form("%s_T0%d",tmp[c].Data(),t0seg));
      h2->Fill(tempx,tempy);
    }
    if(nCVC==1){
      track->XYPosatZ(1600.,tempx,tempy);//cm
      h1=(TH1F*)gFile->Get("HitPatTOF_ifBLC");
      h1->Fill(tofsseg);
      h2=(TH2F*)gFile->Get(Form("%s_CVC%d",tmp[c].Data(),tofsseg));
      h2->Fill(tempx,tempy);
    }
  }
  **/
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
#if DEBUG
    std::cout <<name<<"\tlayer" << layer
	      <<"\twire"<<wire << std::endl;    
#endif
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

void EventAnalysisMomentum::HistSemiLocal(){
  TString tmp[2]={"BLC1","BLC2"};
  double dlmax[2]={0.8,0.5};
  for(int iname=0;iname<2;iname++){
    new TH1F( Form("Chi2%s",tmp[iname].Data()),
		Form("Chi square/Dof %s",tmp[iname].Data()),
		1000, 0,100 );
    new TH1F( Form("nTrack%s",tmp[iname].Data()),
		Form("ntrack %s",tmp[iname].Data()),
		20, 0,20 );
    new TH1F( Form("nTrack%sifK",tmp[iname].Data()),
		Form("ntrack %sifK",tmp[iname].Data()),
		20, 0,20 );
    new TH2F( Form("AB%s",tmp[iname].Data()),
		Form("AB%s[mrad]",tmp[iname].Data()),
		100,-50.,50.,100,-50.,50. );
    new TH2F( Form("AB%s2",tmp[iname].Data()),
		Form("AB%s2[mrad]",tmp[iname].Data()),
		100,-50.,50.,100,-50.,50. );
    new TH2F( Form("XY%s",tmp[iname].Data()),
		Form("XY%s[cm]",tmp[iname].Data()),
		100,-15.,15.,100,-15.,15. );
    new TH2F( Form("XA%s",tmp[iname].Data()),
		Form("X[cm]A[mrad]%s",tmp[iname].Data()),
		100,-15.,15.,100,-50.,50. );
    new TH2F( Form("YB%s",tmp[iname].Data()),
		Form("Y[cm]B[mrad]%s",tmp[iname].Data()),
		100,-15.,15.,100,-50.,50. );
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
  for(int i=0;i<5;i++){
    new TH2F(Form("%s_T0%d",tmp[1].Data(),i+1),
	     Form("%s track T0 seg%d image",tmp[1].Data(),i+1) ,
	     100,-10,10,100,-10.,10.0);
  }    
  for(int i=0;i<34;i++){
    new TH2F(Form("%s_CVC%d",tmp[1].Data(),i+1),
	     Form("%s track CVC seg%d image",tmp[1].Data(),i+1) ,
	     800,-400,400,200,-100.,100.0);
  }    
  new TH1F("HitPatCVC_ifBLC","HitPatCVC_ifBLC",
	   34,0.5,34.5 );
  Tools::newTH2F("CorrX_BLC12",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrY_BLC12",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrA_BLC12",	   200,-50.,50.,200,-50.,50. );      
  Tools::newTH2F("CorrB_BLC12",	   200,-50.,50.,200,-50.,50. );      
  Tools::newTH2F("CorrX_BLC12_2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrY_BLC12_2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrXY_BLC12_2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrAB_BLC12_2",	   200,-10.,10.,200,-10.,10. );      
}

void EventAnalysisMomentum::InitializeHistogram()
{
  rtFile->cd();
  HistSemiLocal();  
  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );
  new TH1F("Mom",  "Mom", 1000,900,1100);
  new TH1F("Chi2D5",  "Chi2D5", 500,0,100);

}
EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMomentum *event = new EventAnalysisMomentum();
  return (EventTemp*)event;
}
