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
#include "EventHeader.h"
#include "ScalerMan.h"

#include "CDSTrackingMan.h"


#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#define TRACK 1

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
  EventHeader *header;
  ScalerMan *scaMan;
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

  int argc2;
  char **argv2;
  
  TRint *theApp = new TRint( "theApp", &argc2, argv2 );
  
  confMan = conf;
  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  //scaMan = new ScalerMan();

  disp = new Display();
  canvas1 = new TCanvas( "canvas1", "canvas1", 700, 700 );
  std::cout << " disp:" << disp << " canvas1:" << canvas1 << std::endl;
  canvas1->Divide(2,2);
  disp->SetCDSFrameXY(-70,70,-70,70);
  disp->SetCDSFrameYZ(-70,70,-70,70);
  disp->SetBLDCFrameXZ(-25,25,-15,15);
  disp->SetBLDCFrameYZ(-25,25,-15,15);
}

void EventDisplay::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventDisplay::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  //scaMan->SetBlockEventNumber( Block_Event_Number );
  //for( int i=0; i<nsca; i++ ){
  //  scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  //}
#if 0
  //std::cout << nsca<< std::endl;
  //for( int i=0; i<nsca; i++ ){
  //  std::cout << "  " << sca[i];
  //}
  //std::cout << std::endl;
#endif
  //scaMan->Clear();
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

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

#if TRACK  
  CDSTrackingMan *trackMan = new CDSTrackingMan();
  BeamLineTrackMan *blTrackMan=new BeamLineTrackMan();

  /***********************/
  /**** Tracking **********/
  /***********************/
  //  trackMan->SetKilledLayer(1);				       
  trackMan->Execute(cdsMan,confMan);
  blTrackMan->DoTracking(blMan,confMan);
  /***********************/
  /**********************/
    
  int nGoodTrack=trackMan->nGoodTrack();

  /*
  if( nGoodTrack<1 ){
    header->Clear();
    blMan->Clear();
    blTrackMan->Clear();
    cdsMan->Clear();
    trackMan->Clear();
    return true;
  }
  */
  for(int n=0;n<nGoodTrack;n++)
    {
      int GoodTrack=trackMan->GoodTrackID(n);
      trackMan->CalcVertex_beam(GoodTrack,blTrackMan,confMan);
    }
  std::cout<<"nTrack BPC "<<" "<<blTrackMan->ntrackBPC()<<std::endl;

#endif
  TVirtualPad *pad;
  pad = canvas1->GetPad(1);
  disp->DrawCDSFrameXY( pad );
  disp->DrawSegmentsXY( pad, confMan, CID_CDH );
  disp->DrawCDCLayersXY( pad, confMan );
  disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_CDH );
  disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_CDC );

  pad = canvas1->GetPad(2);
  disp->DrawCDSFrameYZ( pad );
  disp->DrawCDCLayersYZ( pad, confMan );
  disp->DrawCDSHitYZ( pad, confMan, cdsMan, CID_CDH );
  disp->DrawCDSHitYZ( pad, confMan, cdsMan, CID_CDC );

  pad = canvas1->GetPad(3);
  disp->DrawBLDCFrameXZ( pad );
  disp->DrawBLDCLayersXZ( pad, confMan, CID_BLC1 );
  disp->DrawBLDCLayersXZ( pad, confMan, CID_BLC2 );
  disp->DrawBLDCHit( pad, confMan,blMan, CID_BLC1,0 );
  disp->DrawBLDCHit( pad, confMan,blMan, CID_BLC2,0 );

  pad = canvas1->GetPad(4);
  disp->DrawBLDCFrameYZ( pad );
  disp->DrawBLDCLayersYZ( pad, confMan, CID_BLC1 );
  disp->DrawBLDCLayersYZ( pad, confMan, CID_BLC2 );
  disp->DrawBLDCHit( pad, confMan,blMan, CID_BLC1,1 );
  disp->DrawBLDCHit( pad, confMan,blMan, CID_BLC2,1 );

#if TRACK           
  
  canvas1->cd(3);
  std::cout<<"nTrack BLC,1,2 : "<<blTrackMan->ntrackBLC()<<" "<<blTrackMan->ntrackBLC1()<<" "<<blTrackMan->ntrackBLC2()<<std::endl;



   for(int itr=0;itr<blTrackMan->ntrackBLC1();itr++)
     {
       LocalTrack *blc1=blTrackMan->trackBLC1(itr);
       BLDCWireMapMan *BLwireman=confMan->GetBLDCWireMapManager();
       double gx,gy,gz,dx,dy,dz;
       BLwireman->GetGParam(CID_BLC1,gx,gy,gz,dx,dy,dz);
       TVector3 Gpos,dpos;
       Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);
       double z1,z2,x1,x2,y1,y2;
       z1=-25;z2=25;
       blc1->XYPosatZ(z1+Gpos.z(),x1,y1);
       blc1->XYPosatZ(z2+Gpos.z(),x2,y2);
       std::cout<<"BLC1 chi zx : zy : all "<<blc1->chi2xz()<<" : "<<blc1->chi2yz()<<" : "<<blc1->chi2all()<<std::endl;
       TLine line;
       line.SetLineStyle(2);
       line.SetLineColor(itr+2);
       line.DrawLine(z1,x1,z2,x2);
       double a,b,c;
       blc1->abc(a,b,c);
       std::cout<<"BLC1 a : b : c "<<a<<" : "<<b<<" : "<<c<<std::endl;
       for(int ih=0;ih<blc1->nhit();ih++)
 	{
 	  ChamberLikeHit *hit=blc1->hit(ih);
 	  TMarker blt_m;	
 	  double hpos,wpos,dltrack;
 	  if(hit->xy()==0) 
 	    {
 	      hpos=hit->x(); wpos=hit->wx();dltrack=fabs(hpos-wpos);
 	    }  
 	  else if(hit->xy()==1) 
 	    {
 	      hpos=hit->y(); wpos=hit->wy();dltrack=fabs(hpos-wpos);
 	    }  
 	  if(hit->xy()==0)
 	    {
 	      blt_m.SetMarkerStyle(20);
 	      blt_m.SetMarkerSize(0.2);
 	      blt_m.SetMarkerColor(3  );
 	      blt_m.DrawMarker(hit->wz(),hit->wx() );
 	    }
 	}
     }

   for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
     {
       LocalTrack *blc2=blTrackMan->trackBLC2(itr);
       BLDCWireMapMan *BLwireman=confMan->GetBLDCWireMapManager();
       double gx,gy,gz,dx,dy,dz;
       BLwireman->GetGParam(18,gx,gy,gz,dx,dy,dz);
       TVector3 Gpos,dpos;
       Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);
       double z1,z2,x1,x2,y1,y2;
       z1=-25;z2=25;
       blc2->XYPosatZ(z1+Gpos.z(),x1,y1);
       blc2->XYPosatZ(z2+Gpos.z(),x2,y2);
       std::cout<<"BLC2 chi zx : zy : all "<<blc2->chi2xz()<<" : "<<blc2->chi2yz()<<" : "<<blc2->chi2all()<<std::endl;
       TLine line;
       line.SetLineColor(itr+2);
       line.DrawLine(z1,x1,z2,x2);
       double a,b,c;
       blc2->abc(a,b,c);
       std::cout<<"BLC2 a : b : c "<<a<<" : "<<b<<" : "<<c<<std::endl;
       for(int ih=0;ih<blc2->nhit();ih++)
 	{
 	  ChamberLikeHit *hit=blc2->hit(ih);
 	  TMarker blt_m;	
 	  double hpos,wpos,dltrack;
 	  if(hit->xy()==0) 
 	    {
 	      hpos=hit->x(); wpos=hit->wx();dltrack=fabs(hpos-wpos);
 	    }  
 	  else if(hit->xy()==1) 
 	    {
 	      hpos=hit->y(); wpos=hit->wy();dltrack=fabs(hpos-wpos);
 	    }  
 	  if(hit->xy()==0)
 	    {
 	      blt_m.SetMarkerStyle(20);
 	      blt_m.SetMarkerSize(0.2);
 	      blt_m.SetMarkerColor(3  );
 	      blt_m.DrawMarker(hit->wz(),hit->wx() );
 	    }
 	}
     }
/*  
  for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
    {
      LocalTrack *blc=blTrackMan->trackBLC(itr);
      BLDCWireMapMan *BLwireman=confMan->GetBLDCWireMapManager();
      double gx,gy,gz,dx,dy,dz;
      BLwireman->GetGParam(CID_BLC1,gx,gy,gz,dx,dy,dz);
      TVector3 Gpos,dpos;
      Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);
      double z1,z2,x1,x2,y1,y2;
      z1=-25;z2=25;
      blc->XYPosatZ(z1+Gpos.z(),x1,y1);
      blc->XYPosatZ(z2+Gpos.z(),x2,y2);
      std::cout<<"BLC chi zx : zy : all "<<blc->chi2xz()<<" : "<<blc->chi2yz()<<" : "<<blc->chi2all()<<std::endl;
      TLine line;
      line.SetLineWidth(3);
      line.SetLineStyle(4);
      line.SetLineColor(itr+2);
      line.DrawLine(z1,x1,z2,x2);
      double a,b,c;
      blc->abc(a,b,c);
      //std::cout<<"BLC a : b : c "<<a<<" : "<<b<<" : "<<c<<std::endl;
      for(int ih=0;ih<blc->nhit();ih++)
	{
	  ChamberLikeHit *hit=blc->hit(ih);
	  TMarker blt_m;	
	  double hpos,wpos,dltrack;
	  if(hit->xy()==0) 
	    {
	      hpos=hit->x(); wpos=hit->wx();dltrack=fabs(hpos-wpos);
	    }  
	  else if(hit->xy()==1) 
	    {
	      hpos=hit->y(); wpos=hit->wy();dltrack=fabs(hpos-wpos);
	    }  
	  if(hit->xy()==0)
	    {
	      blt_m.SetMarkerStyle(20);
	      blt_m.SetMarkerSize(0.2);
	      blt_m.SetMarkerColor(1  );
	      blt_m.DrawMarker(hit->wz(),hit->wx() );
	    }
	}
    }
*/
  canvas1->cd(4);  

   for(int itr=0;itr<blTrackMan->ntrackBLC1();itr++)
     {
       LocalTrack *blc1=blTrackMan->trackBLC1(itr);
       BLDCWireMapMan *BLwireman=confMan->GetBLDCWireMapManager();
       double gx,gy,gz,dx,dy,dz;
       BLwireman->GetGParam(18,gx,gy,gz,dx,dy,dz);
       TVector3 Gpos,dpos;
       Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);
       double z1,z2,x1,x2,y1,y2;
       z1=-20;z2=20;
       blc1->XYPosatZ(z1+Gpos.z(),x1,y1);
       blc1->XYPosatZ(z2+Gpos.z(),x2,y2);
       TLine line;
       line.SetLineStyle(2);
       line.SetLineColor(itr+2);
       line.DrawLine(z1,y1,z2,y2);
       for(int ih=0;ih<blc1->nhit();ih++)
 	{
 	  ChamberLikeHit *hit=blc1->hit(ih);
 	  TMarker blt_m;	  
 	  if(hit->xy()==1)
 	    {
 	      blt_m.SetMarkerStyle(20);
 	      blt_m.SetMarkerSize(0.2);
 	      blt_m.SetMarkerColor(3  );
 	      blt_m.DrawMarker(hit->wz(),hit->wy() );
 	    }
 	}

     }

   for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
     {
       LocalTrack *blc2=blTrackMan->trackBLC2(itr);
       BLDCWireMapMan *BLwireman=confMan->GetBLDCWireMapManager();
       double gx,gy,gz,dx,dy,dz;
       BLwireman->GetGParam(18,gx,gy,gz,dx,dy,dz);
       TVector3 Gpos,dpos;
       Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);
       double z1,z2,x1,x2,y1,y2;
       z1=-20;z2=20;
       blc2->XYPosatZ(z1+Gpos.z(),x1,y1);
       blc2->XYPosatZ(z2+Gpos.z(),x2,y2);
       TLine line;
       line.SetLineColor(itr+2);
       line.DrawLine(z1,y1,z2,y2);
       for(int ih=0;ih<blc2->nhit();ih++)
 	{
 	  ChamberLikeHit *hit=blc2->hit(ih);
 	  TMarker blt_m;	  
 	  if(hit->xy()==1)
 	    {
 	      blt_m.SetMarkerStyle(20);
 	      blt_m.SetMarkerSize(0.2);
 	      blt_m.SetMarkerColor(3  );
 	      blt_m.DrawMarker(hit->wz(),hit->wy() );
 	    }
 	}
     }
   /*
  for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
    {
      LocalTrack *blc=blTrackMan->trackBLC(itr);
      BLDCWireMapMan *BLwireman=confMan->GetBLDCWireMapManager();
      double gx,gy,gz,dx,dy,dz;
      BLwireman->GetGParam(18,gx,gy,gz,dx,dy,dz);
      TVector3 Gpos,dpos;
      Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);
      double z1,z2,x1,x2,y1,y2;
      z1=-20;z2=20;
      blc->XYPosatZ(z1+Gpos.z(),x1,y1);
      blc->XYPosatZ(z2+Gpos.z(),x2,y2);
      TLine line;
      line.SetLineStyle(4);
      line.SetLineStyle(3);
      line.SetLineColor(itr+2);
      line.DrawLine(z1,y1,z2,y2);
      for(int ih=0;ih<blc->nhit();ih++)
	{
	  ChamberLikeHit *hit=blc->hit(ih);
	  TMarker blt_m;	  
	  if(hit->xy()==1)
	    {
	      blt_m.SetMarkerStyle(20);
	      blt_m.SetMarkerSize(0.2);
	      blt_m.SetMarkerColor(1  );
	      blt_m.DrawMarker(hit->wz(),hit->wy() );
	    }
	}

    }
   */
  

  canvas1->Update();


  
//   for(int i=0;i<blTrackMan->ntrackBLC2();i++)
//     { 
//       double gx,gy,gz,dx,dy,dz;
//       BLDCWireMapMan *BLwireman=confMan->GetBLDCWireMapManager();
//       BLwireman->GetGParam(18,gx,gy,gz,dx,dy,dz);
//       TVector3 Gpos,dpos;
//       Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);     
//       LocalTrack *blc2=blTrackMan->trackBLC2(i);
//       TVector3 x1,x2;
//       double xp1,xp2,yp1,yp2;
//       //      blc2->XYPosatZ(-60-Gpos.z(),xp1,yp1);x1.SetXYZ(xp1,yp1,-60-Gpos.z());
//       //blc2->XYPosatZ(60-Gpos.z(),xp2,yp2);x2.SetXYZ(xp2,yp2,60-Gpos.z());
//       //      x1=x1+Gpos; x2=x2+Gpos;
//       blc2->XYPosatZ(-60,xp1,yp1);x1.SetXYZ(xp1,yp1,-60);
//       blc2->XYPosatZ(60,xp2,yp2);x2.SetXYZ(xp2,yp2,60);
//       x1=x1; x2=x2;
//       TLine line;
//       TMarker bl_mk;
//       bl_mk.SetMarkerStyle(20);
//       bl_mk.SetMarkerSize(0.8);
//       canvas1->cd(1);	
//       line.DrawLine(x1.x(),x1.y(),x2.x(),x2.y());
//       bl_mk.DrawMarker(x1.x(),x1.y());
//       canvas1->cd(2);	
//       line.DrawLine(x1.z(),x1.y(),x2.z(),x2.y());
//       bl_mk.DrawMarker(x1.z(),x1.y());
//     }

  
//     for(int numslayer=1;numslayer<=7;numslayer++)
//       {
// 	for(int i=0;i<trackMan->nCluster(numslayer);i++)
// 	  {
// 	    HitCluster *cluster=trackMan->Cluster(numslayer,i);
// 	    double Meanphi=cluster->Meanphi();		
// 	    double Meanradius=cluster->Meanradius();		
	    
// 	    TVector3 cls_pos;
// 	    cls_pos.SetXYZ( Meanradius*cos(Meanphi/180*TMath::Pi()),
// 			    Meanradius*sin(Meanphi/180*TMath::Pi()),
// 			    0);
// 	    //		seg=track->CDHHit()->seg();		
// 	    TMarker cls_mk;
// 	    cls_mk.SetMarkerStyle(24);
// 	    cls_mk.SetMarkerColor(1);
// 	    cls_mk.SetMarkerSize(1.2);
// 	    canvas1->cd(1);
// 	    cls_mk.DrawMarker(cls_pos.x(),cls_pos.y());
// 	    //	    std::cout<<" ncl= "<<cluster->nHit()<<" phi "<<cluster->Meanphi()<<std::endl; 
// 	  }
//       }
/****  
       int ver_id=0;    
    for(int i=0;i<trackMan->nGoodTrack();i++)
      {
	
	for(int ii=i+1;ii<trackMan->nGoodTrack();ii++)
	  {
	    TrackVertex vertex;
	    //    if(!trackMan->GetVertex(trackMan->GoodTrackID(i),trackMan->GoodTrackID(ii),vertex) ) continue;
	    if(!trackMan->GetVertex(trackMan->GoodTrackID(i),trackMan->GoodTrackID(ii),vertex) ) continue;
	    TVector3 ver_mean=vertex.GetVertexPos_mean();
	    TMarker ver_mk;
	    ver_mk.SetMarkerStyle(20);
	    ver_mk.SetMarkerColor(2+ver_id);
	    ver_mk.SetMarkerSize(1.2);
	    canvas1->cd(1);
	    ver_mk.DrawMarker(ver_mean.x(),ver_mean.y());
	    canvas1->cd(2);
	    ver_mk.DrawMarker(ver_mean.z(),ver_mean.y());
	    ver_id++;
	    std::cout<<" vrx_cdc2= "<<"x y z "<<ver_mean.x()<<" "<<ver_mean.y()<<" "<<ver_mean.z()<<std::endl; 
	  }

	TrackVertex vertex_beam;
	if(!trackMan->GetVertex_beam(trackMan->GoodTrackID(i),0,vertex_beam)) continue;
	TVector3 verb_mean=vertex_beam.GetVertexPos_mean();
	TVector3 verb1=vertex_beam.GetVertexPos1();
	TVector3 verb2=vertex_beam.GetVertexPos2();
	TMarker verb_mk;
	verb_mk.SetMarkerStyle(22);
	verb_mk.SetMarkerColor(2+i);
	verb_mk.SetMarkerSize(1.1);
	canvas1->cd(1);
	verb_mk.DrawMarker(verb_mean.x(),verb_mean.y());
	canvas1->cd(2);
	verb_mk.DrawMarker(verb_mean.z(),verb_mean.y());
 
	std::cout<<"cdc vertex = "<<verb1.x()<<" "<<verb1.y()<<" "<<
	  verb1.z()<<std::endl;
	std::cout<<"blc vertex = "<<verb2.x()<<" "<<verb2.y()<<" "<<
	  verb2.z()<<std::endl;
****/
// 	for(int i=0;i<blTrackMan->ntrackBLC2();i++)
// 	  { 
// 	    LocalTrack *blc2=blTrackMan->trackBLC2(i);
// 	    TVector3 x1,x2;
// 	    double xp1,xp2,yp1,yp2,tan;
// 	    blc2->XYPosatZ(-74,xp1,yp1);x1.SetXYZ(xp1,yp1,-74);
// 	    blc2->XYPosatZ(74,xp2,yp2);x2.SetXYZ(xp2,yp2,74);
// 	    tan=148.0/sqrt( (xp1-xp2)*(xp1-xp2)+(yp1-yp2)*(yp1-yp2) );
// 	    double rho_b=-900.0/tan/0.5/2.99;
// 	    double phi=-(verb2.z()+74)/(rho_b*tan);
// 	    //  MathTools *tool=new MathTools();
// 	    //double phi0=tool->CalcDeg(xp1,yp1)*TMath::DegToRad();
// 	    // delete tool;
// 	    double d_r=sqrt(xp1*xp1+yp1*yp1);
// 	    double vx,vy;
// 	    //std::cout<<" d_r "<<d_r<<" phi0 "<<phi0<<" rho_b "<<rho_b<<" tan "<<tan<<"  phi_b "<<phi<<std::endl;
// 	    //vx=d_r*cos(phi0)+rho_b*(cos(phi0)-cos(phi0+phi));
// 	    //vy=d_r*sin(phi0)+rho_b*(sin(phi0)-sin(phi0+phi));
// 	    //	    std::cout<<"blc2 vertex = "<<vx<<" "<<vy<<" "<<
// 	    //     verb2.z()<<std::endl;
// 	  }
  /***
 	std::cout<<"vertex dis = "<<(verb1-verb2).Mag()<<std::endl;
      }
	
    std::cout<<" ntrack= "<<trackMan->nTrack()<<" goodtrack= "<<trackMan->nGoodTrack()<<std::endl; 
    double arho,ax_c,ay_c,aPt;
    double aparam[5];
    TF1 *func_y;	
    for(int i=0;i<trackMan->nTrack();i++)
      {
	CDSTrack *track=trackMan->Track(i);
	bool goodflag=track->GoodFlag();
	if(!goodflag) continue;
	track->GetParameters(aparam);

	int nhit[15]={0};
	double jitter,theta,y_mid,chi2,seg;
	bool resiflag=false;
	for(int numlayer=1;numlayer<=15;numlayer++)
	  {
	    for(int n=0;n<track->nTrackHit(numlayer);n++)
	      {
		TVector3 hpos,wpos,wposp,cdhpos;
		double dl=track->TrackHit(numlayer,n)->dl();		
		double resi=track->TrackHit(numlayer,n)->resl();		
		//		std::cout<<"layer= "<<numlayer<<" resid= "<<resi<<std::endl;
		if(fabs(resi)>0.05) resiflag=true;
		jitter=track->Jitter();		
		//y_mid=track->MidY();		
		theta=track->Theta();		
		chi2=track->Chi();		
	
		cdhpos.SetXYZ(track->CDHHit()->x(),
			    track->CDHHit()->y(),
			    track->CDHHit()->z() );
		//		seg=track->CDHHit()->seg();		
		TMarker mcdhpos;
		mcdhpos.SetMarkerStyle(20);
		if(!goodflag) mcdhpos.SetMarkerStyle(5);
	        mcdhpos.SetMarkerColor(2+i);
		mcdhpos.SetMarkerSize(0.7);
		canvas1->cd(1);
		mcdhpos.DrawMarker(cdhpos.x(),cdhpos.y());

		
		hpos.SetXYZ(track->TrackHit(numlayer,n)->x(),
			    track->TrackHit(numlayer,n)->y(),
			    track->TrackHit(numlayer,n)->z() );
		//std::cout<<"hit "<<hpos.x()<<" "<<hpos.y()<<" "<<hpos.z()<<std::endl;
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
	if(goodflag)
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
	    
	    TF1 *func_y=new TF1("func_y","([0])*sin([1])+1./[2]*(sin([1])-sin( [1]+([3]-x)/([4]*1./[2]) ) )",aparam[3]-1./aparam[2]*aparam[4]*(-TMath::Pi()/2. ),aparam[3]-1./aparam[2]*aparam[4]*(TMath::Pi()/2. ) );
	    func_y->SetParameters(aparam[0],aparam[1],aparam[2],aparam[3],aparam[4]);
	    func_y->SetLineColor(2+i);
	    func_y->SetLineWidth(1);
	    
	    pad = canvas1->GetPad(2);
	    pad->cd();	
	    
	    func_y->Draw("LPSAME");
	  }

	//   ntrack++;
	
      }
  ***/
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
