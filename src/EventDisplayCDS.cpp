#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <TRint.h>

#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "CDSTrackingMan.h"
#include "TrackTools.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#define TRACK 1
#define BLC 0
#define CDCTRACK 1
#define BEAM 1

class EventDisplay: public EventTemp
{
public:
  EventDisplay();
  ~EventDisplay();
private:
  Display *disp;
  TCanvas *canvas1;

  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  BeamLineTrackMan *blTrackMan;
  EventHeader *header;
  EventHeader *header2;
  ScalerMan *scaMan;
  CDSTrackingMan *trackMan;
  TTree *cdcTree;
  int CDC_Event_Number;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventDisplay::EventDisplay()
  : EventTemp()
{
}

EventDisplay::~EventDisplay()
{
}

const int MaxTreeSize = 1000000000;
void EventDisplay::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventDisplay::Initialize " << std::endl;
#endif

  int argc2=1;
  char *temp=(char *)malloc(100);
  strcpy(temp,conf->GetProgramName().data());
  char **argv2=&temp;
  new TRint( "theApp", &argc2, argv2 );
  
  confMan = conf;
  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  blTrackMan = new BeamLineTrackMan();
  //scaMan = new ScalerMan();

  trackMan=new CDSTrackingMan();
  header2 =new EventHeader();
  TString cdcfname=confMan->GetCDSTrackFileName();
  TFile *cdcFile = new TFile(cdcfname);
  if(!cdcFile->IsOpen()){
    std::cout<<" failed to open " <<cdcfname<< "  !!!"<<std::endl;
    exit(false);
  }
  if(cdcFile->IsZombie()){
    std::cout<<cdcfname<< " is Zombie !!!"<<std::endl;
    exit(false);
  }

  cdcTree=(TTree*)cdcFile->Get("EventTree");
  cdcTree->SetBranchAddress( "CDSTrackingMan", &trackMan );
  cdcTree->SetBranchAddress( "EventHeader" ,&header2);
  CDC_Event_Number=0;
  disp = new Display();
  canvas1 = new TCanvas( "canvas1", "canvas1", 1400, 1400 );
  std::cout << " disp:" << disp << " canvas1:" << canvas1 << std::endl;
  canvas1->Divide(2,2);
  disp->SetCDSFrameXY(-70,70,-70,70);
  disp->SetCDSFrameYZ(-70,70,-70,70);
  disp->SetBLDCFrameXZ(-25,25,-15,15);
  disp->SetBLDCFrameYZ(-25,25,-15,15);
}

void EventDisplay::USca( int nsca, unsigned int *sca )
{
}

bool EventDisplay::UAna( TKOHitCollection *tko )
{
#if 1
  std::cout << " Enter EventDisplay::UAna " << std::endl;
#endif

  Event_Number++;

  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%1==0 )
    std::cout << " Event# : " << Event_Number << std::endl;

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  blTrackMan->Clear();
  //  trackMan->Clear();
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  std::cout<< "convert" <<std::endl;  
  header->Convert( tko, confMan );
#if BEAM
  if(header->IsTrig(Trig_Cosmic)) return true;
#else
  if(!header->IsTrig(Trig_Cosmic)) return true;
#endif
  // check event matching
  //  bool CDCTRACK=false;
  cdcTree->GetEntry(CDC_Event_Number);
  if(header2->ev()==Event_Number){
    //    CDCTRACK=true;
    CDC_Event_Number++;
  }else{
    trackMan->Clear();
  }

  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

  std::vector<int> cdhseg;
  std::vector<int> ihseg;
  for(int i=0;i<cdsMan->nCDH();i++)
    if(cdsMan->CDH(i)->CheckRange()){
      cdhseg.push_back(cdsMan->CDH(i)->seg());
    }
  for(int i=0;i<cdsMan->nIH();i++)
    if(cdsMan->IH(i)->CheckRange()){
      ihseg.push_back(cdsMan->IH(i)->seg());
    }
  if(cdhseg.size()>=trackMan->nTrack()){
    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
  }
  blTrackMan->DoTracking(blMan,confMan,true,true);
  if(blTrackMan->ntrackBLDC(CID_BPC)!=1) return true;
  if(blTrackMan->ntrackBLDC(CID_BLC1)!=1) return true;
  if(blTrackMan->ntrackBLDC(CID_BLC2)!=1) return true;
  LocalTrack *bpctrack=0;
  if(blTrackMan->ntrackBLDC(CID_BPC)==1)
    bpctrack=blTrackMan->trackBPC(0);    
  int nT0=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
    }
  }
  int nBHD=0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      nBHD++;
    }
  }
  if(nT0!=1) return true;
  if(nBHD!=1) return true;

  std::cout<<"nCDH: "<<cdhseg.size()<<" ,nTrack: "<<trackMan->nTrack()<<std::endl;
  TVirtualPad *pad;
  pad = canvas1->GetPad(1);
  disp->DrawCDSFrameXY( pad );
  disp->DrawSegmentsXY( pad, confMan, CID_CDH );
  disp->DrawSegmentsXY( pad, confMan, CID_IH );
  disp->DrawCDCLayersXY( pad, confMan );
  disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_CDH );
  disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_IH );
  disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_CDC );

  pad = canvas1->GetPad(2);
  disp->DrawCDSFrameYZ( pad );
  disp->DrawCDCLayersYZ( pad, confMan );
  disp->DrawSegmentsYZ( pad, confMan, CID_CDH );
  disp->DrawCDSHitYZ( pad, confMan, cdsMan, CID_CDH );
  disp->DrawCDSHitYZ( pad, confMan, cdsMan, CID_CDC );

#if CDCTRACK
  std::cout<<" ntrack= "<<trackMan->nTrack()<<" goodtrack= "<<trackMan->nGoodTrack()<<std::endl; 
  double arho,ax_c,ay_c,aPt;
  double aparam[5];
  double gparam[5];
  TF1 *func_y;	
  TVector3 vertex(-999,-999,-999);
  double disvertex=999;  
  int nCDC=0;
  for(int i=0;i<trackMan->nTrack();i++)
    {
      CDSTrack *track=trackMan->Track(i);
      if(track->Chi()>30) continue;
      TVector3 cdhpos=track->CDHVertex();
      if(TMath::Abs(cdhpos.Z())<30) nCDC++;
      if(bpctrack){
	TVector3 vtxb1, vtxb2, vtxb;
	double tmpdis;      
	TrackTools::CalcLineHelixVertex(bpctrack,track,vtxb1,vtxb2,tmpdis);
	if(vtxb1.Mag()<1e-10){ std::cout<<"calc vertex failed"<<std::endl; continue; }
	if(tmpdis<disvertex){
	  disvertex=tmpdis;
	  vertex=vtxb1;
	}    
      }
    }
  if(nCDC<=cdhseg.size()) return true;
  if(GeomTools::GetID(vertex)!=CID_Fiducial) return true;

  for(int i=0;i<trackMan->nTrack();i++)
    {
      CDSTrack *track=trackMan->Track(i);
      //      if(track->Chi()>30) continue;
      track->SearchHodoHit(cdsMan,confMan,CID_CDH,7,false);
      track->SetHitPos(cdsMan);      
      bool goodflag=track->GoodFlag();
      std::cout<<"------------  "<<i<<"  "<<track->Chi()<<std::endl;
      //      if(!goodflag) continue;
      track->GetParameters(aparam);
      std::cout<<"cdhflag:  "<<track->CDHFlag()<<std::endl;      
      int nhit[15]={0};
      double jitter,theta,y_mid,chi2,seg;
      bool resiflag=false;

      for(int numlayer=1;numlayer<=15;numlayer++)
	{
	  for(int n=0;n<track->nTrackHit(numlayer);n++)
	    {
	      TVector3 hpos,wpos,wposp,cdhpos;
	      double dl=track->TrackHit(cdsMan,numlayer,n)->dl();	       
	      double resi=track->TrackHit(cdsMan,numlayer,n)->resl();		
	      //	      std::cout<<"layer= "<<numlayer<<" resid= "<<resi<<std::endl;
	      if(fabs(resi)>0.05) resiflag=true;
	      jitter=track->Jitter();		
	      //y_mid=track->MidY();		
	      theta=track->Theta();		
	      chi2=track->Chi();		
	      	      
	      hpos.SetXYZ(track->TrackHit(cdsMan,numlayer,n)->x(),
			  track->TrackHit(cdsMan,numlayer,n)->y(),
			  track->TrackHit(cdsMan,numlayer,n)->z() );
	      TMarker mhpos;
	      mhpos.SetMarkerStyle(i+20);
	      mhpos.SetMarkerColor(2+i);
	      mhpos.SetMarkerSize(0.5);
	      canvas1->cd(2);
	      mhpos.DrawMarker(hpos.z(),hpos.y());	      
	      canvas1->cd(1);
	      mhpos.DrawMarker(hpos.x(),hpos.y());	      
	      nhit[numlayer-1]++;
	    }
	  
	}
      
      //	if(!resiflag) continue;
      //      if(goodflag)
	{
	  
	  if(aparam[2]!=0)
	    {
	      arho=fabs(1./aparam[2]);
	      ax_c=(aparam[0]+1./aparam[2])*cos(aparam[1]);
	      ay_c=(aparam[0]+1./aparam[2])*sin(aparam[1]);
	      aPt=track->Pt();
	      std::cout<<"Pt"<<"= "<<aPt<<"GeV/c"<<std::endl;
	    }
	  std::cout<<" jitter= "<<jitter<<" theta= "<<theta<<" rchi2 "<<chi2<<std::endl;
	  
	}
      std::cout<<"param= "
	       <<aparam[0]<<" "
	       <<aparam[1]<<" "
	       <<aparam[2]<<" "
	       <<aparam[3]<<" "
	       <<aparam[4]<<" "
	       <<std::endl;
      std::cout<<"       x,y,r "<<ax_c<<" "<<ay_c<<" "<<arho<<std::endl;
      if(aparam[2]==0)
	{
	  pad = canvas1->GetPad(1);
	  pad->cd();
	  TLine lxy,lzy;
	  lxy.SetLineColor(2+i);
	  double a,b,c;
	  a=track->A(); b=track->B(); c=track->C();
	  
	  double x1,y1,z1,x2,y2,z2;
	  // x1=0;x2=50;
	  //y1=-c/b;y2=-(a*50+c)/b;
	  
	  x1=aparam[0]*cos(aparam[1])-(0)*sin(aparam[1]);
	  y1=aparam[0]*sin(aparam[1])+(0)*cos(aparam[1]);
	  z1=aparam[3]-(0)*aparam[4];
	  x2=aparam[0]*cos(aparam[1])-60*sin(aparam[1]);
	  y2=aparam[0]*sin(aparam[1])+60*cos(aparam[1]);
	  z2=aparam[3]-60*aparam[4];
	  
	  lxy.DrawLine(x1,y1,x2,y2);
	  
	  lzy.SetLineColor(2+i);	    
	  pad = canvas1->GetPad(2);
	  pad->cd();	
	  lzy.DrawLine(z1,y1,z2,y2);	  
	}
      else
	{    
	  pad = canvas1->GetPad(1);
	  pad->cd();
	  TArc arc;
	  arc.SetFillStyle(0);
	  arc.SetLineColor(2+i);
	  arc.DrawArc(ax_c,ay_c,arho,0,360);
	  //	  TVector3 cdhpos=MathTools::CalcHelixGPosatR(aparam,gparam,55.);
	  TVector3 cdhpos=track->CDHVertex();
	  TVector3 ihpos=track->IHVertex();
	  cdhpos.Print();
	  TMarker mark;
	  mark.SetMarkerColor(2+i);
	  mark.SetMarkerStyle(20);
	  mark.DrawMarker(cdhpos.X(),cdhpos.Y());
	  mark.DrawMarker(ihpos.X(),ihpos.Y());

	  pad = canvas1->GetPad(2);
	  pad->cd();	
	  mark.DrawMarker(cdhpos.Z(),cdhpos.Y());
	  
	  //	  double phi=MathTools::CalcDeg(cdhpos.X(),cdhpos.Y());
	  
	  TF1 *func_y=new TF1("func_y","([0])*sin([1])+1./[2]*(sin([1])-sin( [1]+([3]-x)/([4]*1./[2]) ) )",aparam[3]-1./aparam[2]*aparam[4]*(-TMath::Pi()/2. ),aparam[3]-1./aparam[2]*aparam[4]*(TMath::Pi()/2. ) );
	  func_y->SetParameters(aparam[0],aparam[1],aparam[2],aparam[3],aparam[4]);
	  func_y->SetLineColor(2+i);
	  func_y->SetLineWidth(1);
	  
	  func_y->Draw("LPSAME");
	}
      for(int i=0;i<track->nCDHHit();i++){
	HodoscopeLikeHit *hit = track->CDHHit(cdsMan,i);
	std::cout<<"CDH seg: "<<hit->seg()<<"  , eu: "<<hit->eu()<<"  , ed: "<<hit->ed()<<" , hit: "<<hit->CheckRange()<<std::endl;
      }
      //   ntrack++;
      
    }
  
#endif    
  
  
  canvas1->Update();
  canvas1->SaveAs("tmpdis.ps");
  bool status = disp->Wait();
  if( !status ) return status;

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
#if TRACK
  trackMan->Clear();
#endif
  return true;
}

void EventDisplay::Finalize()
{
  std::cout << " Enter EventDisplay::Finalize " << std::endl;

  delete blMan;
  delete cdsMan;
  delete header;

  gSystem->Exit(1);
  gROOT->GetApplication()->Run();

}

EventTemp *EventAlloc::EventAllocator()
{
  EventDisplay *event = new EventDisplay();
  return (EventTemp*)event;
}
