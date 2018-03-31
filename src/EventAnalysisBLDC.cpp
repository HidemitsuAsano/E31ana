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

#include "Tools.h"
#include "time.h"

#define TREE 0
#define TRACK 1
#define NCDH 0
#define WIRE 0
#define SEMILOCAL 0
#define LINEAR 0
#define BLC2DEAD 0
#define DEBUG 0

const int MAXEVENT=10000;
class EventAnalysisBLDC: public EventTemp
{
public:
  EventAnalysisBLDC();
  ~EventAnalysisBLDC();

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

  TString name[5];
  int nwire[5];
  int cid[5];
  double dlmax[5];
  bool FLAG[5];
  double tof_bhd;
  int nT0;
  int t0seg;
  int nCVC;
  int tofsseg;
  int nCDH;
  int trigregi[20];

  void HistWireCorr();
  void FillWireCorr();
  void HistBLDCCorr();
  void FillBLDCCorr();
  void HistBLCBPCCorr();
  void FillBLCBPCCorr();
  void HistTrig();
  void FillTrig();
  void HistRaw();
  void FillRaw(const int &c);
  void HistTrack();

  void HistTrackT0FF();
  void FillTrackT0FF(const int &ic,LocalTrack *track);
  void FillTrack( const int &ic,LocalTrack *track);
  void HistSemiLocal();
  void FillSemiLocal( const int &ic,LocalTrack *track);
  //  void HistLinear();
  void FillLinear( const int &ic,LinearTrack *track);

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  void UTime( int time ){};
  bool UAna( TKOHitCollection *tko );
  void Finalize();

  void InitializeHistogram();
};

EventAnalysisBLDC::EventAnalysisBLDC()
  : EventTemp()
{
}

EventAnalysisBLDC::~EventAnalysisBLDC()
{
}

//const int MaxTreeSize = 1900000000000;
void EventAnalysisBLDC::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisBLDC::Initialize " << std::endl;
#endif
  //TODO : avoid hard code (H.Asano)
  name[0]="BLC1a";  name[1]="BLC1b";  name[2]= "BLC2a";  name[3]="BLC2b";  name[4]="BPC";
  nwire[0]=32;  nwire[1]=32;  nwire[2]=32;  nwire[3]=32;  nwire[4]=32;
  cid[0]=15;  cid[1]=16;  cid[2]=17;  cid[3]=18;  cid[4]=40;
  dlmax[0]=0.4;  dlmax[1]=0.4;  dlmax[2]=0.25;  dlmax[3]=0.25; /* dlmax[4]=0.36;*/ dlmax[4]=0.30;

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

void EventAnalysisBLDC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisHodo::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool EventAnalysisBLDC::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisBLDC::UAna " << std::endl;
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
#if NCDH
  cdsMan->Convert( tko, confMan );  
#endif
  TH1F *h1;
#if TRACK
  trackMan->DoTracking( blMan, confMan );
  for(int i=0;i<5;i++){
    h1=(TH1F*)gFile->Get(Form("%sstatus",name[i].Data()));
    h1->Fill(0);
    h1->Fill(trackMan->status(cid[i]));
    if(trackMan->status(CID_BPC)==1){
      h1=(TH1F*)gFile->Get(Form("%sstatusifBPC",name[i].Data()));
      h1->Fill(0);
      h1->Fill(trackMan->status(cid[i]));
    }      
  }
  TH2F* h2;
  TString tmp[2]={"BLC1","BLC2"};
  int tmpid[2]={CID_BLC1,CID_BLC2};
#if LINEAR
  trackMan->LinearTracking( blMan, confMan, CID_BLC1 );
  trackMan->LinearTracking( blMan, confMan, CID_BLC2 );
  trackMan->ConvertLocalToLinear(  CID_BLC1  );
  trackMan->ConvertLocalToLinear(  CID_BLC2  );
  h1=(TH1F*)gFile->Get(Form("BLC1status"));
  h1->Fill(0);
  h1->Fill(trackMan->status(CID_BLC1));
  h1=(TH1F*)gFile->Get(Form("BLC2status"));
  h1->Fill(0);
  h1->Fill(trackMan->status(CID_BLC2));
#endif
#if 1
  //  if(header->proton()&&
  if(trackMan->nltrackBLDC(tmpid[0])==1&&trackMan->nltrackBLDC(tmpid[1])==1){
     //     &&trackMan->ltrackBLDC(tmpid[0],0)->chi2all()<10
     //     &&trackMan->ltrackBLDC(tmpid[1],0)->chi2all()<10){
    //    beam->TMinuitFit(trackMan->ltrackBLC1(0),trackMan->ltrackBLC2(0),confMan);
    //    std::cout<<"Momentum,Chi2: "<<beam->mom()<<"\t"<<beam->chisquare()<<std::endl;
    double blc1par[6];
    double parblc2[6];
    trackMan->ltrackBLDC(tmpid[0],0)->gabcd(blc1par[0],blc1par[1],blc1par[2],blc1par[3]);
    blc1par[1]=TMath::ATan(blc1par[1])*1000;
    blc1par[3]=TMath::ATan(blc1par[3])*1000;
    blc1par[4]=0.;
    blc1par[5]=0.;
    beam->CalcParBLC1toBLC2(blc1par,parblc2);
    //    std::cout<<parblc2[0]<<std::endl;
    double tempx,tempy;
    trackMan->ltrackBLDC(tmpid[1],0)->XYPosatZ(-130,tempx,tempy);
    double a,b,c,d;
    trackMan->ltrackBLDC(tmpid[1],0)->gabcd(a,b,c,d);
    h2=(TH2F*)gFile->Get("CorrX_BLC1BLC2_2");
    h2->Fill(tempx,parblc2[0]);
    h2=(TH2F*)gFile->Get("CorrY_BLC1BLC2_2");
    h2->Fill(tempy,parblc2[2]);
    h2=(TH2F*)gFile->Get("CorrXY_BLC1BLC2_2");
    h2->Fill(tempx-parblc2[0],tempy-parblc2[2]);
    h2=(TH2F*)gFile->Get("CorrAB_BLC1BLC2_2");
    h2->Fill(TMath::ATan(b)*1000-parblc2[1],TMath::ATan(d)*1000-parblc2[3]);
  }
#endif
#endif   
  rtFile->cd();

  int ntr;

  for(int i=0;i<5;i++)
    FLAG[i]=true;
  nT0=0; t0seg=0; nCVC=0; tofsseg=0; nCDH=0; tof_bhd=-1.;

  FillTrig();

  //  TH2F *h2;
  for(int c=0;c<5;c++){
#if DEBUG
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s): " << (t1-t0)*1.0e-6 << std::endl;
#endif
    for( int layer=1; layer<=NumOfBLCLayers; layer++ ){    
      h1 = (TH1F*)gFile->Get( Form("Mul%s_%d",name[c].Data(),layer) );
      h1->Fill( blMan->nBLDC(cid[c],layer) );
      if(blMan->nBLDC(cid[c],layer)!=1||!(blMan->BLDC(cid[c],layer,0))) FLAG[c]=false;
      //      if(blMan->nBLDC(cid[c],layer)>10){
#if DEBUG
      std::cout<<blMan->nBLDC(cid[c],layer)<<std::endl;
      for(int i=0;i<blMan->nBLDC(cid[c],layer);i++)
	std::cout<<name[c]<<"\t"<<layer<<"\t"<<blMan->BLDC(cid[c],layer,i)->wire()
		 <<"\t"<<blMan->BLDC(cid[c],layer,i)->dt()
		 <<"\t"<<blMan->BLDC(cid[c],layer,i)->dl()
		 <<"\t"<<blMan->BLDC(cid[c],layer,i)->CheckRange()
		 <<std::endl;
#endif
    }
    if(blMan->nBLDC(cid[c],1)==1&&blMan->nBLDC(cid[c],8)==1){
      h1=(TH1F*)gFile->Get( Form("Eff18%s",name[c].Data()) );
      for( int layer=1; layer<=NumOfBLCLayers; layer++ )
	if(blMan->nBLDC(cid[c],layer)>=1)
	  h1->Fill(layer);
      if(blMan->BLDC(cid[c],1,0)->wire()==(int)nwire[cid[c]]/2.&&
	 blMan->BLDC(cid[c],8,0)->wire()==(int)nwire[cid[c]]/2.){
	h1=(TH1F*)gFile->Get( Form("Eff18%scenter",name[c].Data()) );
	for( int layer=1; layer<=NumOfBLCLayers; layer++ )
	  if(blMan->nBLDC(cid[c],layer)==1)
	    h1->Fill(layer);
      }
    }
#if TRACK
    ntr=trackMan->ntrackBLDC(cid[c]);
    h1=(TH1F*)gFile->Get(Form("nTrack%s",name[c].Data())); h1->Fill(ntr);
    if(trigregi[11]>10&&trigregi[11]<2000&&tof_bhd>28){
      h1=(TH1F*)gFile->Get(Form("nTrack%sifK",name[c].Data())); h1->Fill(ntr);
    }
    if(trigregi[7]>0&&nCDH<2){
      h1=(TH1F*)gFile->Get(Form("nTrack%s_2",name[c].Data())); h1->Fill(ntr);
    }
    if(ntr!=1) continue;      
    for(int itr=0;itr<ntr;itr++){
      LocalTrack *track=trackMan->trackBLDC(cid[c],itr);
      double chi2=track->chi2all();
      h1=(TH1F*)gFile->Get(Form("Chi2%s",name[c].Data())); h1->Fill(chi2);
      if(chi2>20) continue;
      FillTrack(c,track);
      if(nT0!=1||ntr!=1) continue;
      FillTrackT0FF(c,track);
    }
#else
    if(FLAG[c]) FillRaw(c);
#endif
  }

#if TRACK
#if SEMILOCAL
  int ntrack[2];
  for(int i=0;i<2;i++){
    ntr=trackMan->ntrackBLDC(tmpid[i]);
    ntrack[i]=ntr;
    h1=(TH1F*)gFile->Get(Form("nTrack%s",tmp[i].Data())); h1->Fill(ntr);
    if(trigregi[11]>10&&trigregi[11]<1700&&tof_bhd>28){
      h1=(TH1F*)gFile->Get(Form("nTrack%sifK",tmp[i].Data())); h1->Fill(ntr);
    }
    if(trackMan->status(CID_BPC)==1){
      h1=(TH1F*)gFile->Get(Form("nTrack%sifBPC",tmp[i].Data())); h1->Fill(ntr);
    }
    if(ntr!=1) continue;      
    for(int itr=0;itr<ntr;itr++){
      LocalTrack *track=trackMan->trackBLDC(tmpid[i],itr);
      double chi2=track->chi2all();
      h1=(TH1F*)gFile->Get(Form("Chi2%s",tmp[i].Data())); h1->Fill(chi2);
      if(chi2>20) continue;
      FillSemiLocal(i,track);
    }
  }
  h2=(TH2F*)gFile->Get(Form("nTrackBLC1BLC2")); h2->Fill(ntrack[0],ntrack[1]);
  if(trackMan->status(CID_BPC)==1){
    h2=(TH2F*)gFile->Get(Form("nTrackBLC1BLC2ifBPC")); h2->Fill(ntrack[0],ntrack[1]);
  }
  if(trackMan->ntrackBLDC(tmpid[0])==1&&trackMan->ntrackBLDC(tmpid[1])==1){
    double blc1x,blc1y,blcx,blcy;
    trackMan->trackBLDC(CID_BLC1,0)->XYPosatZ(0.,blc1x,blc1y);
    trackMan->trackBLDC(CID_BLC2,0)->XYPosatZ(-130,blcx,blcy);
    h2=(TH2F*)gFile->Get("CorrX_BLC1BLC2");
    h2->Fill(blcx,blc1x);
    h2=(TH2F*)gFile->Get("CorrY_BLC1BLC2");
    h2->Fill(blcy,blc1y);
  }
#elif LINEAR
  for(int i=0;i<2;i++){
    ntr=trackMan->nltrackBLDC(tmpid[i]);
    h1=(TH1F*)gFile->Get(Form("nTrack%s",tmp[i].Data())); h1->Fill(ntr);
    if(trigregi[11]>10&&trigregi[11]<1700&&tof_bhd>28){
      h1=(TH1F*)gFile->Get(Form("nTrack%sifK",tmp[i].Data())); h1->Fill(ntr);
    }
    if(ntr!=1) continue;      
    for(int itr=0;itr<ntr;itr++){
      LinearTrack *track=trackMan->ltrackBLDC(tmpid[i],itr);
      double chi2=track->chi2all();
      h1=(TH1F*)gFile->Get(Form("Chi2%s",tmp[i].Data())); h1->Fill(chi2);
      if(chi2>20) continue;
      FillLinear(i,track);
    }
  }
  if(trackMan->nltrackBLDC(tmpid[0])==1&&trackMan->nltrackBLDC(tmpid[1])==1){
    double blc1x,blc1y,blcx,blcy;
    trackMan->ltrackBLDC(CID_BLC1,0)->XYPosatZ(0.,blc1x,blc1y);
    trackMan->ltrackBLDC(CID_BLC2,0)->XYPosatZ(-130,blcx,blcy);
    h2=(TH2F*)gFile->Get("CorrX_BLC1BLC2");
    h2->Fill(blcx,blc1x);
    h2=(TH2F*)gFile->Get("CorrY_BLC1BLC2");
    h2->Fill(blcy,blc1y);
  }

  if(trackMan->nltrackBLDC(CID_BLC2)==1&&trackMan->ntrackBLDC(CID_BPC)==1){
    double bpcx,bpcy,blcx,blcy;
    trackMan->ltrackBLDC(CID_BLC2,0)->XYPosatZ(-75,blcx,blcy);
    trackMan->trackBLDC(CID_BPC,0)->XYPosatZ(-75,bpcx,bpcy);
    //    std::cout<<"blcx,y bpcx,y\t"<<blcx<<"\t"<<blcy<<"\t"<<bpcx<<"\t"<<bpcy<<std::endl;
    h2=(TH2F*)gFile->Get("CorrX_BLC2BPC");
    h2->Fill(blcx,bpcx);
    h2=(TH2F*)gFile->Get("CorrY_BLC2BPC");
    h2->Fill(blcy,bpcy);
    FillBLCBPCCorr();
  }
#endif
#endif

  if(nT0==1) FillBLDCCorr();
  FillWireCorr();
#if TREE
  evTree->Fill();
#endif
  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  trackMan->Clear();
  return true;
}

void EventAnalysisBLDC::Finalize()
{
  std::cout << " Enter EventAnalysisBLDC::Finalize " << std::endl;
  rtFile->cd();
  confMan->SaveCode();
  confMan->SaveParams();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete header;
  delete cdsMan;
  delete trackMan;
}

void EventAnalysisBLDC::FillTrig(){
  TH1F *h1;
  //  TH2F *h2;
  for(int i=0;i<20;i++)
    trigregi[i]=-1;

  double t0ctm=0.;
  for( int i=0; i<blMan->nT0(); i++ ){
    HodoscopeLikeHit *hit = blMan->T0(i);
    int seg = hit->seg();
    if(hit->CheckRange()){
      t0seg=seg;
      nT0++;
      t0ctm=hit->ctmean();
    }
  }

  for( int i=0; i<blMan->nCVC(); i++ ){
    HodoscopeLikeHit *hit = blMan->CVC(i);
    int seg = hit->seg();
    if(hit->CheckRange()){
      tofsseg=seg;
      nCVC++;
    }
  }

  if(nT0==1){
    for( int i=0; i<blMan->nBHD(); i++ ){
      HodoscopeLikeHit *hit = blMan->BHD(i);
      if(hit->CheckRange()){
	tof_bhd=t0ctm-hit->ctmean();
      }
    }
  }
  h1=(TH1F*)gFile->Get("BHD_T0TOF");
  h1->Fill(tof_bhd);
  
  if(cdsMan!=0){
    for(int i=0; i<cdsMan->nCDH(); i++ ){
      HodoscopeLikeHit *hit=cdsMan->CDH(i);
      if(hit->CheckRange()) nCDH++;     
    }
  }
  for(int i=1;i<20;i++){
    trigregi[i]=header->pattern(i);
    if(trigregi[i]>0) {
      h1=(TH1F*)gFile->Get(Form("TrigRegi%d",i));
      h1->Fill(trigregi[i]);
      h1=(TH1F*)gFile->Get("HitPatTrigRegi");
      h1->Fill(i);
    }
  }
  h1=(TH1F*)gFile->Get(Form("TrigCheck"));
  h1->Fill(0);
  if(trigregi[6]>0||trigregi[8]>0||trigregi[9]>0||trigregi[11]>0)
    h1->Fill(1);
  if(trigregi[13]>0||trigregi[14]>0||trigregi[16]>0)
    h1->Fill(2);
  int trigpat=0;
  bool no2nd=true;
  if(trigregi[2]>0) trigpat=2;
  if(trigregi[4]>0) trigpat=8;
  if(trigregi[8]>0){
    if(trigregi[14]>0) 
      trigpat=3;
    else if(trigregi[13]>0) 
      trigpat=9;
    else if(trigregi[11]>0) 
      trigpat=10;
    //    if(header->pattern(11)>0) 
    //      trigpat=9;
    else if(trigregi[9]<0){
      trigpat=11;
      for(int i=1;i<20;i++){
	if(trigregi[i]>0) {
	  if(i!=8) no2nd=false;
	  h1=(TH1F*)gFile->Get("HitPatTrigRegiifpi");
	  h1->Fill(i);
	}
      }
      if(no2nd){
	h1=(TH1F*)gFile->Get("HitPatTrigRegiifpi");
	h1->Fill(0);
      }
    }
  }
  if(trigregi[9]>0) trigpat=1; // cdh1*pi
  if(trigregi[7]>0) trigpat=6; // cdh2*K
  if(trigregi[6]>0){ //cdh1*K
    if(trigregi[13]>0)     trigpat=4; // *charged
    if(trigregi[14]>0)     trigpat=5; // *neutral
  }

  if(trigregi[14]>0)
    for(int i=1;i<20;i++)
      if(trigregi[i]>0) {
	h1=(TH1F*)gFile->Get("HitPatTrigRegiifNtrig");
	h1->Fill(i);
      }
   
  h1=(TH1F*)gFile->Get("TrigPat");
  h1->Fill(trigpat);
  if(tof_bhd>28&&trackMan->ntrackBLDC(CID_BPC)>0){
    h1=(TH1F*)gFile->Get("TrigPatifKwBPC");
    h1->Fill(trigpat);
  }
}
void EventAnalysisBLDC::FillWireCorr(){
  //  TH1F *h1;
  TH2F *h2;
  for(int i=0;i<5;i++){
    if(FLAG[i]){
      for(int ilay1=1;ilay1<=8;ilay1++)
	for(int ilay2=ilay1+1;ilay2<=8;ilay2++){
	  h2=(TH2F*)gFile->Get(Form("Corr%s_%d_%s_%d",name[i].Data(),ilay1,name[i].Data(),ilay2));
	  h2->Fill(blMan->BLDC(cid[i],ilay1,0)->wire(),blMan->BLDC(cid[i],ilay2,0)->wire());
	}
    }  
  }
}
void EventAnalysisBLDC::FillRaw(const int &c){
  TH1F *h1;
  //  TH2F *h2;
  for( int layer=1; layer<=NumOfBLCLayers; layer++ ){    
    for( int i=0; i<blMan->nBLDC(cid[c],layer); i++ ){      
      ChamberLikeHit *hit = blMan->BLDC(cid[c],layer,i);
      int wire = hit->wire();
      double dt = hit->dt();
      int tdc = hit->tdc();
      //	    std::cout<<"wire: "<<wire<<"\ttdc: "<<tdc<<"\tdt: "<<dt<<std::endl;
      h1 = (TH1F*)gFile->Get( Form("HitPat%s_%d",name[c].Data(),layer) ); h1->Fill( wire );
      h1 = (TH1F*)gFile->Get( Form("h_dt_%s_%d",name[c].Data(),layer) ); h1->Fill( dt );
#if WIRE
      h1 = (TH1F*)gFile->Get( Form("h_dt_%s_%d_%d",name[c].Data(),layer,wire) ); h1->Fill( dt );
      h1 = (TH1F*)gFile->Get( Form("h_tdc_%s_%d_%d",name[c].Data(),layer,wire) ); h1->Fill( tdc );
#endif
    }
  }
}

void EventAnalysisBLDC::FillTrackT0FF(const int &c,LocalTrack* track){
  //  TH1F *h1;
  TH2F *h2;
  double tempx,tempy;
  track->XYLocalPosatZ(20.,tempx,tempy);//cm
  h2=(TH2F*)gFile->Get(Form("%s_T0%d",name[c].Data(),t0seg));
  h2->Fill(tempx,tempy);
  h2=(TH2F*)gFile->Get(Form("%s_T0",name[c].Data()));
  h2->Fill(tempx,tempy);
  if(trigregi[13]<0&&trigregi[14]<0&&trigregi[16]<0){
    track->XYLocalPosatZ(20.,tempx,tempy);//cm
    h2=(TH2F*)gFile->Get(Form("%s_GT0%d",name[c].Data(),t0seg));
    h2->Fill(tempx,tempy);
    h2=(TH2F*)gFile->Get(Form("%s_GT0",name[c].Data()));
    h2->Fill(tempx,tempy);
  }
  if(trigregi[7]>0&&nCDH>1){
    h2=(TH2F*)gFile->Get(Form("%s_GT0%d_2",name[c].Data(),t0seg));
    h2->Fill(tempx,tempy);
    h2=(TH2F*)gFile->Get(Form("%s_GT0_2",name[c].Data()));
    h2->Fill(tempx,tempy);
  }

  if(header->kaon()&&tof_bhd>28){
    track->XYPosatZ(0.,tempx,tempy);//cm
	h2=(TH2F*)gFile->Get(Form("%s_FF_kaon",name[c].Data()));
	h2->Fill(tempx,tempy);
  }
  if(header->pattern(11)>10&&header->pattern(11)<1700&&tof_bhd>28){
    track->XYPosatZ(0.,tempx,tempy);//cm
    h2=(TH2F*)gFile->Get(Form("%s_FF_kaon2",name[c].Data()));
    h2->Fill(tempx,tempy);    
  }
  if((header->beam()&&header->kaon())||header->pattern(11)>10){
    track->XYPosatZ(0.,tempx,tempy);//cm
    h2=(TH2F*)gFile->Get(Form("%s_FF_ktrig",name[c].Data()));
    h2->Fill(tempx,tempy);
  }
  if(header->pattern(5)>10&&tof_bhd>28){
    track->XYPosatZ(0.,tempx,tempy);//cm
    h2=(TH2F*)gFile->Get(Form("%s_FF_k_cdh",name[c].Data()));
    h2->Fill(tempx,tempy);
  }

  track->XYPosatZ(-100.,tempx,tempy);//cm
  h2=(TH2F*)gFile->Get(Form("%s_GT0%d_3",name[c].Data(),t0seg));
  h2->Fill(tempx,tempy);
  h2=(TH2F*)gFile->Get(Form("%s_GT0_3",name[c].Data()));
  h2->Fill(tempx,tempy);
}
void EventAnalysisBLDC::FillTrack(const int &c,LocalTrack* track)
{
  TH1F* h1;
  TH2F* h2;
  double a1,b1,c1,d1,e1,f1;
  double dt,dltrack,dl;
  double resi=-99.;
  int wire,layer,lr;
  int fac[2]={1,-1};
  TString strlr[2]={"L","R"};    
  double chi2=track->chi2all();  
  track->abc(a1,b1,c1);
  track->def(d1,e1,f1);
  h2=(TH2F*)gFile->Get(Form("AB%s",name[c].Data())); h2->Fill(b1,e1);
  double tempx,tempy;	
  track->XYLocalPosatZ(0.,tempx,tempy);//cm
  h2=(TH2F*)gFile->Get(Form("XY%s",name[c].Data())); h2->Fill(tempx,tempy);
  
  for( int ih=0;ih<track->nhit();ih++ ){
    ChamberLikeHit *hit=track->hit(ih);
    layer=hit->layer();
    wire=hit->wire();
    dl=hit->dl();
#if DEBUG
    std::cout <<name<<"\tlayer" << layer
	      <<"\twire"<<wire << std::endl;    
#endif
    dt=hit->dt();
    if(hit->xy()) dltrack=hit->y()-hit->wy();
    else          dltrack=hit->x()-hit->wx();
    resi=hit->resl();
    lr=hit->leftright();
    h1 = (TH1F*)gFile->Get( Form("HitPat%s_%d",name[c].Data(),layer) );
    h1->Fill( wire );
    h2=(TH2F*)gFile->Get(Form("h_dt_dltrack_%s_%d",name[c].Data(),layer)); 
    h2->Fill(dt,dltrack);
    h2=(TH2F*)gFile->Get(Form("h_dl_resi_%s_%d",name[c].Data(),layer)); 
    h2->Fill(dl*fac[lr],resi);
    h1=(TH1F*)gFile->Get(Form("h_resi_%s_%d",name[c].Data(),layer));
    h1->Fill(resi);
    h1=(TH1F*)gFile->Get(Form("h_resi_%s_%dA",name[c].Data(),layer)); 
    h1->Fill(resi*fac[lr]);
    h1=(TH1F*)gFile->Get(Form("h_resi_%s_%d%s",name[c].Data(),layer,strlr[lr].Data()));
    h1->Fill(resi);
    h2=(TH2F*)gFile->Get(Form("h_dt_resi_%s_%d",name[c].Data(),layer));
    h2->Fill(dt*fac[lr],resi);
#if WIRE
    h1=(TH1F*)gFile->Get(Form("h_resi_%s_%d_%dA",name[c].Data(),layer,wire)); 
    h1->Fill(resi*fac[lr]);
    h2=(TH2F*)gFile->Get(Form("h_dt_resi_%s_%d_%d",name[c].Data(),layer,wire));
    h2->Fill(dt*fac[lr],resi);
#endif
    if(chi2<5){
      h1=(TH1F*)gFile->Get(Form("h_dt_%s_%d",name[c].Data(),layer));
      h1->Fill(dt);
#if WIRE
      h1=(TH1F*)gFile->Get(Form("h_dt_%s_%d_%d",name[c].Data(),layer,wire));
      h1->Fill(dt);
      h1 = (TH1F*)gFile->Get( Form("h_tdc_%s_%d_%d",name[c].Data(),layer,wire) ); 
      h1->Fill( hit->tdc() );
#endif
    }
  }
}

void EventAnalysisBLDC::HistWireCorr(){
  for(int iname=0;iname<5;iname++){
    for(int ilay1=1;ilay1<=8;ilay1++)
      for(int ilay2=ilay1+1;ilay2<=8;ilay2++)
	new TH2F(Form("Corr%s_%d_%s_%d",name[iname].Data(),ilay1,name[iname].Data(),ilay2),
		 Form("Corr%s_%d_%s_%d",name[iname].Data(),ilay1,name[iname].Data(),ilay2),
		 nwire[iname]+1,0,nwire[iname]+1,nwire[iname]+1,0,nwire[iname]+1);
  }
}

void EventAnalysisBLDC::HistTrig(){
  for(int i=0;i<20;i++)
    new TH1F( Form("TrigRegi%d",i), Form("Trigger register sa%d",i),
	      4096,-0.5,4095.5);
  new TH1F( "HitPatTrigRegi","TrigPat",20,-0.5,19.5);
  new TH1F( "HitPatTrigRegiifpi","TrigPat",20,-0.5,19.5);
  new TH1F( "HitPatTrigRegiifNtrig","TrigPat",20,-0.5,19.5);
  new TH1F( "TrigPat","TrigPat",20,-0.5,19.5);
  new TH1F( "TrigPatifKwBPC","TrigPatifKwBPC",20,-0.5,19.5);
  new TH1F( "TrigPatifNtrig","TrigPatifNtrig",20,-0.5,19.5);
  new TH1F( "TrigCheck","TrigCheck",5,-0.5,4.5);

  new TH1F("BHD_T0TOF","BHD_T0 Time Of Flight",4000,0,100);
}

void EventAnalysisBLDC::HistRaw(){
  new TH1F( Form("BLC1status"),
	    Form("tracking status of BLC1"),
	    11, 0,11 );
  new TH1F( Form("BLC2status"),
	    Form("tracking status of BLC2"),
	    11, 0,11 );

  for(int iname=0;iname<5;iname++){
    new TH1F( Form("%sstatus",name[iname].Data()),
	      Form("tracking status of %s",name[iname].Data()),
	      11, 0,11 );
    new TH1F( Form("%sstatusifBPC",name[iname].Data()),
	      Form("tracking status of %s",name[iname].Data()),
	      11, 0,11 );
    new TH1F( Form("Eff18%s",name[iname].Data()),
	      Form("Efficiency if layer 18 hit %s",name[iname].Data()),
	      10, 0,10 );
    new TH1F( Form("Eff18%scenter",name[iname].Data()),
	      Form("Efficiency if layer 18 hit %s",name[iname].Data()),
	      10, 0,10 );
    for(int i=0;i<8;i++){
      new TH1F( Form("Mul%s_%d",name[iname].Data(),i+1),
		Form("Multiplicity %s Layer%d",name[iname].Data(),i+1),
		nwire[iname]+1, -0.5, nwire[iname]+0.5 );
      new TH1F( Form("HitPat%s_%d",name[iname].Data(),i+1),
		Form("Hit Pattern %s Layer%d",name[iname].Data(),i+1),
		nwire[iname], 0.5, nwire[iname]+0.5 );
      new TH1F(Form("h_dt_%s_%d",name[iname].Data(),i+1),
	       Form("dt in %s layer%d",name[iname].Data(),i+1) ,
	       3000,-200,400);
#if WIRE
      for(int j=1;j<=nwire[iname];j++){
	new TH1F(Form("h_dt_%s_%d_%d",name[iname].Data(),i+1,j),
		 Form("dt in %s layer%d wire%d",name[iname].Data(),i+1,j) ,
		 3000,-200,400);
	new TH1F(Form("h_tdc_%s_%d_%d",name[iname].Data(),i+1,j),
		 Form("tdc in %s layer%d wire%d",name[iname].Data(),i+1,j) ,
		 2048,-0.5,2047.5);	
      }
#endif
    }
  }
}

void EventAnalysisBLDC::HistTrackT0FF(){
  for(int iname=0;iname<5;iname++){
    new TH2F(Form("%s_GT0",name[iname].Data()),
	     Form("%s track GT0 image",name[iname].Data()) ,
	     100,-10,10,100,-10.,10.0);    
    new TH2F(Form("%s_GT0_2",name[iname].Data()),
	     Form("%s track GT0 image",name[iname].Data()) ,
	     100,-10,10,100,-10.,10.0);    
    new TH2F(Form("%s_GT0_3",name[iname].Data()),
	     Form("%s track GT0 image",name[iname].Data()) ,
	     100,-10,10,100,-10.,10.0);        
    for(int i=0;i<5;i++){
      new TH2F(Form("%s_T0%d",name[iname].Data(),i+1),
	       Form("%s track T0 seg%d image",name[iname].Data(),i+1) ,
	       100,-10,10,100,-10.,10.0);    
      new TH2F(Form("%s_GT0%d",name[iname].Data(),i+1),
	       Form("%s track GT0 seg%d image",name[iname].Data(),i+1) ,
	       100,-10,10,100,-10.,10.0);    
      new TH2F(Form("%s_GT0%d_2",name[iname].Data(),i+1),
	       Form("%s track GT0 seg%d image 2",name[iname].Data(),i+1) ,
	       100,-10,10,100,-10.,10.0);    
      new TH2F(Form("%s_GT0%d_3",name[iname].Data(),i+1),
	       Form("%s track GT0 seg%d image 3",name[iname].Data(),i+1) ,
	       100,-10,10,100,-10.,10.0);    
    }
    new TH2F(Form("%s_T0",name[iname].Data()),
	     Form("%s track T0 image",name[iname].Data()) ,
	     100,-10,10,100,-10.,10.0);    
  }
  for(int iname=0;iname<5;iname++){
    new TH2F(Form("%s_FF_kaon",name[iname].Data()),
	     Form("%s track FF kaon image ",name[iname].Data()) ,
	     300,-15,15,300,-15.,15.0);    
    new TH2F(Form("%s_FF_kaon2",name[iname].Data()),
	     Form("%s track FF kaon image ",name[iname].Data()) ,
	     300,-15,15,300,-15.,15.0);    
    new TH2F(Form("%s_FF_ktrig",name[iname].Data()),
	     Form("%s track FF kaon tirg image ",name[iname].Data()) ,
	     300,-15,15,300,-15.,15.0);    
    new TH2F(Form("%s_FF_k_cdh",name[iname].Data()),
	     Form("%s track FF kaon image with CDH hit",name[iname].Data()) ,
	     300,-15,15,300,-15.,15.0);    
  }

}
void EventAnalysisBLDC::HistTrack(){
  for(int iname=0;iname<5;iname++){
    new TH1F( Form("Chi2%s",name[iname].Data()),
		Form("Chi square/Dof %s",name[iname].Data()),
		1000, 0,100 );
    new TH1F( Form("nTrack%s",name[iname].Data()),
		Form("ntrack %s",name[iname].Data()),
		20, 0,20 );
    new TH1F( Form("nTrack%sifK",name[iname].Data()),
		Form("ntrack %sifK",name[iname].Data()),
		20, 0,20 );
    new TH1F( Form("nTrack%s_2",name[iname].Data()),
		Form("ntrack %sifK",name[iname].Data()),
		20, 0,20 );
    new TH2F( Form("AB%s",name[iname].Data()),
		Form("B%s",name[iname].Data()),
		100,-0.1,0.1,100,-0.1,0.1 );
    new TH2F( Form("XY%s",name[iname].Data()),
		Form("XY%s",name[iname].Data()),
		100,-15,15,100,-15,15 );
    for(int i=0;i<8;i++){
      new TH2F(Form("h_dt_dltrack_%s_%d",name[iname].Data(),i+1),
	       Form("dt vs dltrack in %s layer%d",name[iname].Data(),i+1) ,
	       200,-50,200,200,-dlmax[iname]*1.2,dlmax[iname]*1.2);
      new TH2F(Form("h_dl_resi_%s_%d",name[iname].Data(),i+1),
	       Form("dl vs resi in %s layer%d",name[iname].Data(),i+1) ,
	       200,-dlmax[iname]*1.2,dlmax[iname]*1.2,200,-0.1,0.1);
      new TH2F(Form("h_dt_resi_%s_%d",name[iname].Data(),i+1),
	       Form("dt vs residu[cm] in %s layer%d",name[iname].Data(),i+1) ,
	       200,-200,200,200,-0.1,0.1);
      new TH1F(Form("h_resi_%s_%d",name[iname].Data(),i+1),
	       Form("residu[cm] in %s layer%d",name[iname].Data(),i+1) ,
	       300,-0.1,0.1);
      new TH1F(Form("h_resi_%s_%dA",name[iname].Data(),i+1),
	       Form("residu[cm] in %s layer%dL",name[iname].Data(),i+1) ,
	       300,-0.1,0.1);
      new TH1F(Form("h_resi_%s_%dL",name[iname].Data(),i+1),
	       Form("residu[cm] in %s layer%dL",name[iname].Data(),i+1) ,
	       300,-0.1,0.1);
      new TH1F(Form("h_resi_%s_%dR",name[iname].Data(),i+1),
	       Form("residu[cm] in %s layer%dR",name[iname].Data(),i+1) ,
	       300,-0.1,0.1);
#if WIRE      
      for(int j=1;j<=nwire[iname];j++){
	new TH2F(Form("h_dt_resi_%s_%d_%d",name[iname].Data(),i+1,j),
		 Form("dt vs residu[cm] in %s layer%d wire %d",name[iname].Data(),i+1,j) ,
		 200,-200,200,200,-0.1,0.1);
	new TH1F(Form("h_resi_%s_%d_%dA",name[iname].Data(),i+1,j),
		 Form("residu[cm] in %s layer%d wire%d",name[iname].Data(),i+1,j) ,
		 300,-0.1,0.1);
      }
#endif
    }
  }
}
void EventAnalysisBLDC::HistBLDCCorr(){
  Tools::newTH2F( "BLC1a_bX", 100,-15,15,100,-1,1 );
  Tools::newTH2F( "BLC1a_bY", 100,-15,15,100,-1,1 );
  Tools::newTH2F( "BLC2a_bX", 100,-15,15,100,-1,1 );
  Tools::newTH2F( "BLC2a_bY", 100,-15,15,100,-1,1 );
  Tools::newTH2F( "BLC1a_bdxdz", 200,-0.1,0.1,200,-0.05,0.05 );
  Tools::newTH2F( "BLC1a_bdydz", 200,-0.1,0.1,200,-0.05,0.05 );
  Tools::newTH2F( "BLC2a_bdxdz", 200,-0.1,0.1,200,-0.05,0.05 );
  Tools::newTH2F( "BLC2a_bdydz", 200,-0.1,0.1,200,-0.05,0.05 );
  Tools::newTH2F( "BLC1Xa_bdy", 100,-15.,15.,100,-2.,2. );
  Tools::newTH2F( "BLC1Ya_bdx", 100,-15.,15.,100,-2.,2. );
  Tools::newTH2F( "BLC2Xa_bdy", 100,-15.,15.,100,-2.,2. );
  Tools::newTH2F( "BLC2Ya_bdx", 100,-15.,15.,100,-2.,2. );
  HistBLCBPCCorr();
}

void EventAnalysisBLDC::FillBLDCCorr(){
  double x1,y1,x2,y2;
  double a1,b1,c1,d1,e1,f1;
  double a2,b2,c2,d2,e2,f2;
  TH2F* h2;
  if(trackMan->ntrackBLC1a()==1&&trackMan->trackBLC1a(0)->chi2all()<5
     &&trackMan->ntrackBLC1b()==1&&trackMan->trackBLC1b(0)->chi2all()<5){    
    LocalTrack* track=trackMan->trackBLC1a(0);
    track->XYSemiLocalPosatZ(0.,x1,y1);//cm
    track->semiabcdef(a1,b1,c1,d1,e1,f1);
    track=trackMan->trackBLC1b(0);
    track->XYSemiLocalPosatZ(0.,x2,y2);//cm
    track->semiabcdef(a2,b2,c2,d2,e2,f2);
    h2=(TH2F*)gFile->Get("BLC1a_bX");
    h2->Fill(x1,x1-x2);
    h2=(TH2F*)gFile->Get("BLC1a_bY");
    h2->Fill(y1,y1-y2);
    h2=(TH2F*)gFile->Get("BLC1Xa_bdy");
    h2->Fill(x1,y1-y2);
    h2=(TH2F*)gFile->Get("BLC1Ya_bdx");
    h2->Fill(y1,x1-x2);
    h2=(TH2F*)gFile->Get("BLC1a_bdxdz");
    h2->Fill(b1,b1-b2);
    h2=(TH2F*)gFile->Get("BLC1a_bdydz");
    h2->Fill(e1,e1-e2);
  }
  if(trackMan->ntrackBLC2a()==1&&trackMan->trackBLC2a(0)->chi2all()<5
     &&trackMan->ntrackBLC2b()==1&&trackMan->trackBLC2b(0)->chi2all()<5){    
    LocalTrack* track=trackMan->trackBLC2a(0);
    track->XYSemiLocalPosatZ(0.,x1,y1);//cm
    track->semiabcdef(a1,b1,c1,d1,e1,f1);
    track=trackMan->trackBLC2b(0);
    track->XYSemiLocalPosatZ(0.,x2,y2);//cm
    track->semiabcdef(a2,b2,c2,d2,e2,f2);
    h2=(TH2F*)gFile->Get("BLC2a_bX");
    h2->Fill(x1,x1-x2);
    h2=(TH2F*)gFile->Get("BLC2a_bY");
    h2->Fill(y1,y1-y2);
    h2=(TH2F*)gFile->Get("BLC2Xa_bdy");
    h2->Fill(x1,y1-y2);
    h2=(TH2F*)gFile->Get("BLC2Ya_bdx");
    h2->Fill(y1,x1-x2);
    h2=(TH2F*)gFile->Get("BLC2a_bdxdz");
    h2->Fill(b1,b1-b2);
    h2=(TH2F*)gFile->Get("BLC2a_bdydz");
    h2->Fill(e1,e1-e2);
  }
}

void EventAnalysisBLDC::FillBLCBPCCorr(){
  double x1,y1,x2,y2;
  double a1,b1,c1,d1,e1,f1;
  double a2,b2,c2,d2,e2,f2;
  TH2F* h2;
  if(trackMan->nltrackBLC2()==1&&trackMan->ltrackBLC2(0)->chi2all()<5
     &&trackMan->ntrackBPC()==1&&trackMan->trackBPC(0)->chi2all()<5){    
    LinearTrack* ltrack=trackMan->ltrackBLC2(0);
    ltrack->XYPosatZ(-75.,x1,y1);//cm
    ltrack->gabcd(a1,b1,c1,d1);
    LocalTrack *track=trackMan->trackBPC(0);
    track->XYPosatZ(-75.,x2,y2);//cm
    track->gabc(a2,b2,c2);
    track->gdef(d2,e2,f2);
    h2=(TH2F*)gFile->Get("BLC2_BPCX");
    h2->Fill(x1,x1-x2);
    h2=(TH2F*)gFile->Get("BLC2_BPCY");
    h2->Fill(y1,y1-y2);
    h2=(TH2F*)gFile->Get("BLC2X_BPCdy");
    h2->Fill(x1,y1-y2);
    h2=(TH2F*)gFile->Get("BLC2Y_BPCdx");
    h2->Fill(y1,x1-x2);
    h2=(TH2F*)gFile->Get("BLC2_BPCdxdz");
    h2->Fill(b1,b1+b2);
    h2=(TH2F*)gFile->Get("BLC2_BPCdydz");
    h2->Fill(d1,d1+e2);
  }
}

void EventAnalysisBLDC::FillSemiLocal(const int &c,LocalTrack* track)
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

  track->XYPosatZ(0.,tempx,tempy);//cm
  h2=(TH2F*)gFile->Get(Form("XYatFF%s",tmp[c].Data())); h2->Fill(tempx,tempy);
  if(trackMan->ntrackBPC()>0){
    h2=(TH2F*)gFile->Get(Form("XYatFF%sifBPC",tmp[c].Data())); h2->Fill(tempx,tempy);
  }
  if(c==1){
    if(nT0==1){
      track->XYLocalPosatZ(20.,tempx,tempy);//cm
      //      std::cout<<"T0seg, x, y by BLC  "<<t0seg<<"\t"<<tempx<<"\t"<<tempy<<std::endl;
      h2=(TH2F*)gFile->Get(Form("%s_T0%d",tmp[c].Data(),t0seg));
      h2->Fill(tempx,tempy);
    }
    if(nCVC==1){
      track->XYPosatZ(1600.,tempx,tempy);//cm
      h1=(TH1F*)gFile->Get("HitPatCVC_ifBLC");
      h1->Fill(tofsseg);
      h2=(TH2F*)gFile->Get(Form("%s_CVC%d",tmp[c].Data(),tofsseg));
      h2->Fill(tempx,tempy);
    }
  }

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
void EventAnalysisBLDC::FillLinear(const int &c,LinearTrack* track)
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
  track->gabcd(a1,b1,c1,d1);
  h2=(TH2F*)gFile->Get(Form("AB%s",tmp[c].Data())); h2->Fill(b1,d1);
  double tempx,tempy;	
  if(c==0)
    track->XYPosatZ(0.,tempx,tempy);//cm
  else
    track->XYPosatZ(-130.,tempx,tempy);//cm
  h2=(TH2F*)gFile->Get(Form("XY%s",tmp[c].Data())); h2->Fill(tempx,tempy);
  track->XYPosatZ(0.,tempx,tempy);//cm
  h2=(TH2F*)gFile->Get(Form("XYatFF%s",tmp[c].Data())); h2->Fill(tempx,tempy);
  if(trackMan->ntrackBPC()>0){
    h2=(TH2F*)gFile->Get(Form("XYatFF%sifBPC",tmp[c].Data())); h2->Fill(tempx,tempy);
  }
  if(c==1){
    if(nT0==1){
      track->XYLocalPosatZ(20.,tempx,tempy);//cm
      //      std::cout<<"T0seg, x, y by BLC  "<<t0seg<<"\t"<<tempx<<"\t"<<tempy<<std::endl;
      h2=(TH2F*)gFile->Get(Form("%s_T0%d",tmp[c].Data(),t0seg));
      h2->Fill(tempx,tempy);
    }
    if(nCVC==1){
      track->XYPosatZ(1600.,tempx,tempy);//cm
      h1=(TH1F*)gFile->Get("HitPatCVC_ifBLC");
      h1->Fill(tofsseg);
      h2=(TH2F*)gFile->Get(Form("%s_CVC%d",tmp[c].Data(),tofsseg));
      h2->Fill(tempx,tempy);
    }
  }
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

void EventAnalysisBLDC::HistSemiLocal(){
  TString tmp[2]={"BLC1","BLC2"};
  Tools::newTH2F( "nTrackBLC1BLC2",  20, 0,20,20,0,20 );
  Tools::newTH2F( "nTrackBLC1BLC2ifBPC",  20, 0,20,20,0,20 );
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
    new TH1F( Form("nTrack%sifBPC",tmp[iname].Data()),
		Form("ntrack %sifBPC",tmp[iname].Data()),
		20, 0,20 );
    new TH2F( Form("AB%s",tmp[iname].Data()),
		Form("B%s",tmp[iname].Data()),
		100,-0.1,0.1,100,-0.1,0.1 );
    new TH2F( Form("XY%s",tmp[iname].Data()),
		Form("XY%s",tmp[iname].Data()),
		100,-15,15,100,-15,15 );
    new TH2F( Form("XYatFF%sifBPC",tmp[iname].Data()),
		Form("XYatFF%sifBPC",tmp[iname].Data()),
		100,-15,15,100,-15,15 );
    new TH2F( Form("XYatFF%s",tmp[iname].Data()),
		Form("XYatFF%s",tmp[iname].Data()),
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
  Tools::newTH1F("HitPatCVC_ifBLC",	   34,0.5,34.5 );
  Tools::newTH2F("CorrX_BLC1BLC2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrY_BLC1BLC2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrA_BLC1BLC2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrB_BLC1BLC2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrX_BLC1BLC2_2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrY_BLC1BLC2_2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrAB_BLC1BLC2_2",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrXY_BLC1BLC2_2",	   200,-10.,10.,200,-10.,10. );      
}

void EventAnalysisBLDC::HistBLCBPCCorr()
{
  Tools::newTH2F("CorrX_BLC2BPC",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrY_BLC2BPC",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrA_BLC2BPC",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F("CorrB_BLC2BPC",	   200,-10.,10.,200,-10.,10. );      
  Tools::newTH2F( "BLC2_BPCX", 100,-15,15,100,-5,5 );
  Tools::newTH2F( "BLC2_BPCY", 100,-15,15,100,-5,5 );
  Tools::newTH2F( "BLC2_BPCdxdz", 200,-0.1,0.1,200,-0.1,0.1 );
  Tools::newTH2F( "BLC2_BPCdydz", 200,-0.1,0.1,200,-0.1,0.1 );
  Tools::newTH2F( "BLC2X_BPCdy", 100,-15.,15.,100,-5.,5. );
  Tools::newTH2F( "BLC2Y_BPCdx", 100,-15.,15.,100,-5.,5. );
}
void EventAnalysisBLDC::InitializeHistogram()
{
  rtFile->cd();  
  HistTrig();
  HistRaw();
  HistWireCorr();
  HistTrack();
  HistTrackT0FF();
  HistBLDCCorr();
  HistSemiLocal();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisBLDC *event = new EventAnalysisBLDC();
  return (EventTemp*)event;
}
