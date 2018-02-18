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
  trackMan->DoTracking( blMan, confMan, true, true );
  if(trackMan->ntrackBLDC(tmpid[0])==1&&trackMan->ntrackBLDC(tmpid[1])==1){
    beam->TMinuitFit(trackMan->trackBLC1(0),trackMan->trackBLC2(0),confMan);
    h1=(TH1F*)gFile->Get("Chi2D5");
    h1->Fill(beam->chisquare());    
    h1=(TH1F*)gFile->Get("MomD5");
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

  int nT0=0,nPC=0,pcseg=-1,t0seg=-1;
  for(int i=0;i<blMan->nT0();i++)
    if(blMan->T0(i)->CheckRange()){
      nT0++;
      t0seg=blMan->T0(i)->seg();
    }
  for(int i=0;i<blMan->nPC();i++)
    if(blMan->PC(i)->CheckRange()){
      nPC++;
      pcseg=blMan->PC(i)->seg();
    }
  // ######################
  // PC analysis
  // #####################
  if(nT0==1&&nPC==1&&trackMan->ntrackBPC()==1){
    if(trackMan->ntrackFDC1()==1){
      TVector3 t0pos=  trackMan->trackBPC(0)->GetPosatZ(-110-0.5);
      TVector3 vertex=  trackMan->trackBPC(0)->GetPosatZ(0);
      vertex=  trackMan->trackBPC(0)->GetPosatZ(0);

      TVector3 pcpos,gfdcpos,gfdcdir,fdc1pos;
      confMan->GetGeomMapManager()->GetGPos(CID_PC,pcseg,pcpos);     
      confMan->GetBLDCWireMapManager()->GetGParam(CID_FDC1,gfdcpos,gfdcdir);
      fdc1pos=trackMan->trackFDC1(0)->GetPosatZ(gfdcpos.Z());

      double t0time=blMan->Hodo(CID_T0,t0seg)->ctmean();
      double pctime=blMan->Hodo(CID_PC,pcseg)->ctmean();
      double mass=pMass;  
      double betab=beam->mom()/sqrt(mass*mass+beam->mom()*beam->mom()); 
      double beamtime=((vertex-t0pos).Mag())/(Const*betab*100);
      
      //  TVector3 fdc1  
      TVector3 pdir=(fdc1pos-vertex).Unit();
      TVector3 uswkpos=fdc1pos+(250-fdc1pos.Z())*pdir;
      double tmpangle=(uswkpos-vertex).Angle(pcpos-uswkpos);
      
      double tofvpc=pctime-t0time-beamtime;
      
      //consider the difference between the lenght of 2 lines and arc trajectory in USWK
      double leff=100; //cm
      double r=leff/sqrt(2*(1-TMath::Cos(tmpangle)));
      double larc=r*tmpangle;      
      double line=2*leff/sqrt(2*(1+TMath::Cos(tmpangle)));
      double diff=line-larc;      
      double fl=(uswkpos-vertex).Mag()+(pcpos-uswkpos).Mag()-diff;
      fl-=1.5*cm;  //PC thickness
      double betap=fl/tofvpc/(Const*100);
      
      // ana mom from tof
      double tmpmom1=mass/sqrt(1/betap/betap-1);
      double tmpmom2;
      double tmpmom3;
      double tmptof1;
      double tmptof2;
      double mom1;
      double mom2;      
      double tmpdiff;
      double ctof=pctime-t0time;

      ELossTools::CalcElossBeamTGeo(t0pos,vertex,tmpmom1, pMass, tmpmom2,tmptof1);
      ELossTools::CalcElossForwardTGeo(vertex, fdc1pos, fl,tmpmom2, pMass, tmpmom3,tmptof2);
      
      //      mom1=2*tmpmom1-tmpmom3;
      mom1=tmpmom1;
      tmpdiff=tmptof1+tmptof2 - ctof;
      mom2=mom1;
      int count=0;
      do{
	mom2+=0.01*tmpdiff;
	if(!ELossTools::CalcElossBeamTGeo(t0pos,vertex,mom2, pMass, tmpmom2,tmptof1)) break;
	if(!ELossTools::CalcElossForwardTGeo(vertex, fdc1pos, fl,tmpmom2, pMass, tmpmom3,tmptof2)) break;
	diff=tmptof1+tmptof2 - ctof;
	count++;
	if(count>50) break; 
      } while (tmpdiff*diff>0);
      if(mom1>mom2){
	double tmp=mom1;
	mom1=mom2;
	mom2=tmp;
      }
      do{
	//	  std::cout<<mom1<<"  "<<mom2<<std::endl;
	double mommid=(mom1+mom2)*0.5;
	if(!ELossTools::CalcElossBeamTGeo(t0pos,vertex,mommid, pMass, tmpmom2,tmptof1))break;
	if(!ELossTools::CalcElossForwardTGeo(vertex, fdc1pos, fl,tmpmom2, pMass, tmpmom3,tmptof2))break;
	diff=tmptof1+tmptof2 - ctof;
	if(diff>0) mom1=mommid;
	else mom2=mommid;
	count++;
	if(count>50) break; 
      } while (TMath::Abs(mom2-mom1)>0.0001);
      
      h1=(TH1F*)gFile->Get("MomPC");
      h1->Fill(beam->mom());
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
  new TH1F("MomD5",  "MomD5", 1000,0.9,1.1);
  new TH1F("MomPC",  "MomPC", 1000,0.9,1.1);
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisBeam *event = new EventAnalysisBeam();
  return (EventTemp*)event;
}
