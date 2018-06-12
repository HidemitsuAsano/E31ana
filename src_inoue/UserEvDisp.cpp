
#include <iostream>
#include <string>
#include <fstream>
#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
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
#include "Display.h"
#include "Display3D.h"

/* ####2011/07/21########
const double MassThPPi = 0.8; // GeV/c2
const double TOFOffset = 3; // nsec
*/
/* #####2011/07/25 slew2######
const double MassThPPi = 0.7; // GeV/c2
const double TOFOffset = -2.5; // nsec
*/

const double MassThPPi = 0.7; // GeV/c2
const double TOFOffset = -3.5; // nsec

const double piMass = 0.138;
const double pMass = 0.934;
const double kpMass = 0.4936;
const double Const=29.97;// [cm/ns]
#define SIM 0


int main( int argc, char **argv )
{
 
  std::cout << " argc:" << argc << std::endl;
  if(argc != 4 )
    std::cout << "Plese set Conffile Outputfile Inputfile  "<< std::endl;

  for(int i=0; i<argc; i++ ){
    std::cout << "  " << argv[i] << std::endl;
  }

  std::string confFile,outFile,inFile;
  if( argv[1][0] != '!' ) confFile = argv[1];
  if( argv[2][0] != '!' ) outFile = argv[2];
  if( argv[3][0] != '!' ) inFile = argv[3];


 
  TRint *theApp = new TRint( "theApp", &argc, argv );
  //TROOT root( "GUI", "GUI" );
  //TApplication theApp( "App", &argc, argv );
  //gROOT->SetStyle( "Plain" );
  gROOT->cd();

  
  gSystem->Load("libPhysics.so");
  gSystem->Load( "./lib/libAll.so" );
  

  ConfMan *conf = new ConfMan(confFile);
  conf->Initialize();


  TFile *f;
  TTree *evtree;


  f = new TFile( inFile.c_str());
  evtree = (TTree*)f->Get( "EventTree" );
  
  CDSHitMan *cdsMan = 0;
  BeamLineHitMan *blMan = 0;
  BeamLineTrackMan *blTrackMan = 0;
  EventHeader *head = 0;
  CDSTrackingMan *trackMan = 0;

  
  TFile *fout;
  fout = new TFile( outFile.c_str(), "recreate" );

  bool AllData=false;
  TObjArray *objarr=evtree->GetListOfBranches();
  if(objarr->FindObject("CDSHitMan")!=0 && objarr->FindObject("CDSTrackingMan")!=0)
    {
      evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
      evtree->SetBranchAddress( "CDSTrackingMan", &trackMan );
      if(objarr->FindObject("BeamLineHitMan")!=0 && 
	 objarr->FindObject("BeamLineTrackMan")!=0 )
      {
	evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
	evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
	evtree->SetBranchAddress( "BeamLineTrackMan", &blTrackMan );
	evtree->SetBranchAddress( "CDSTrackingMan", &trackMan );

      }
      if(objarr->FindObject("EventHeader")!=0)
	{
	  evtree->SetBranchAddress( "EventHeader", &head );
	  AllData=true;
	}
    }
  else 
    {
      std::cout<<"error of EventTree branch"<<std::endl;
      return 0;
    }

  Display *disp = new Display();
  TCanvas *canvas1 = new TCanvas( "canvas1", "canvas1", 700, 700 );
  std::cout << " disp:" << disp << " canvas1:" << canvas1 << std::endl;
  canvas1->Divide(2,2);
  disp->SetCDSFrameXY(-70,70,-70,70);
  disp->SetCDSFrameYZ(-70,70,-70,70);
  disp->SetBLDCFrameXZ(-25,25,-15,15);
  disp->SetBLDCFrameYZ(-25,25,-15,15);

  int nev = evtree->GetEntries();
  std::cout <<" entry:" << nev << std::endl;

  for( int iev=0; iev<nev; iev++ ){
    //    if(iev<12000) continue;

    if( iev%2000 == 0 )
      std::cout << " Event : " << iev << std::endl;
    evtree->GetEvent(iev);

    /***********************/
    /**** Tracking **********/
    /***********************/
    //  trackMan->SetKilledLayer(15);				       
    //  trackMan->Execute(cdsMan,conf);
    blTrackMan->DoTracking(blMan,conf);
//     for(int n=0;n<nGoodTrack;n++)
//       {
// 	int GoodTrack=trackMan->GoodTrackID(n);
// 	trackMan->CalcVertex_beam(GoodTrack,blTrackMan,conf);
//       }
  
    /***********************/
    /**********************/

    
    int nGoodTrack=trackMan->nGoodTrack();

  //   if( nGoodTrack<1 ){
//       header->Clear();
//       blMan->Clear();
//       blTrackMan->Clear();
//       cdsMan->Clear();
//       trackMan->Clear();
//       return true;
//     }
    
  
    TVirtualPad *pad;
    pad = canvas1->GetPad(1);
    disp->DrawCDSFrameXY( pad );
    disp->DrawSegmentsXY( pad, conf, CID_CDH );
    disp->DrawCDCLayersXY( pad, conf );
    disp->DrawCDSHitXY( pad, conf, cdsMan, CID_CDH );
    disp->DrawCDSHitXY( pad, conf, cdsMan, CID_CDC );

    pad = canvas1->GetPad(2);
    disp->DrawCDSFrameYZ( pad );
    disp->DrawCDCLayersYZ( pad, conf );
    disp->DrawCDSHitYZ( pad, conf, cdsMan, CID_CDH );
    disp->DrawCDSHitYZ( pad, conf, cdsMan, CID_CDC );

    pad = canvas1->GetPad(3);
    disp->DrawBLDCFrameXZ( pad );
    disp->DrawBLDCLayersXZ( pad, conf, CID_BLC1 );
    disp->DrawBLDCLayersXZ( pad, conf, CID_BLC2 );
    disp->DrawBLDCHit( pad, conf,blMan, CID_BLC1,0 );
    disp->DrawBLDCHit( pad, conf,blMan, CID_BLC2,0 );

    pad = canvas1->GetPad(4);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawBLDCLayersYZ( pad, conf, CID_BLC1 );
    disp->DrawBLDCLayersYZ( pad, conf, CID_BLC2 );
    disp->DrawBLDCHit( pad, conf,blMan, CID_BLC1,1 );
    disp->DrawBLDCHit( pad, conf,blMan, CID_BLC2,1 );

  
    canvas1->cd(3);
    std::cout<<"nTrack BLC,1,2 : "<<blTrackMan->ntrackBLC()<<" "<<blTrackMan->ntrackBLC1()<<" "<<blTrackMan->ntrackBLC2()<<std::endl;




  for(int itr=0;itr<blTrackMan->ntrackBLC1();itr++)
    {
      LocalTrack *blc1=blTrackMan->trackBLC1(itr);
      BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
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
      //      std::cout<<"BLC1 a : b : c "<<a<<" : "<<b<<" : "<<c<<std::endl;
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
// 	      blt_m.SetMarkerStyle(20);
// 	      blt_m.SetMarkerSize(0.2);
// 	      blt_m.SetMarkerColor(3  );
// 	      blt_m.DrawMarker(hit->wz(),hit->wx() );
	    }
	}
    }

  for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
    {
      LocalTrack *blc2=blTrackMan->trackBLC2(itr);
      BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
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
      //      std::cout<<"BLC2 a : b : c "<<a<<" : "<<b<<" : "<<c<<std::endl;
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
// 	      blt_m.SetMarkerStyle(20);
// 	      blt_m.SetMarkerSize(0.2);
// 	      blt_m.SetMarkerColor(3  );
// 	      blt_m.DrawMarker(hit->wz(),hit->wx() );
	    }
	}
    }
  
  for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
    {
      LocalTrack *blc=blTrackMan->trackBLC(itr);
      BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
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

  canvas1->cd(4);  

//   for(int itr=0;itr<blTrackMan->ntrackBLC1();itr++)
//     {
//       LocalTrack *blc1=blTrackMan->trackBLC1(itr);
//       BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
//       double gx,gy,gz,dx,dy,dz;
//       BLwireman->GetGParam(18,gx,gy,gz,dx,dy,dz);
//       TVector3 Gpos,dpos;
//       Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);
//       double z1,z2,x1,x2,y1,y2;
//       z1=-20;z2=20;
//       blc1->XYPosatZ(z1+Gpos.z(),x1,y1);
//       blc1->XYPosatZ(z2+Gpos.z(),x2,y2);
//       TLine line;
//       line.SetLineStyle(2);
//       line.SetLineColor(itr+2);
//       line.DrawLine(z1,y1,z2,y2);
//       for(int ih=0;ih<blc1->nhit();ih++)
// 	{
// 	  ChamberLikeHit *hit=blc1->hit(ih);
// 	  TMarker blt_m;	  
// 	  if(hit->xy()==1)
// 	    {
// 	      blt_m.SetMarkerStyle(20);
// 	      blt_m.SetMarkerSize(0.2);
// 	      blt_m.SetMarkerColor(3  );
// 	      blt_m.DrawMarker(hit->wz(),hit->wy() );
// 	    }
// 	}

//     }

//   for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
//     {
//       LocalTrack *blc2=blTrackMan->trackBLC2(itr);
//       BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
//       double gx,gy,gz,dx,dy,dz;
//       BLwireman->GetGParam(18,gx,gy,gz,dx,dy,dz);
//       TVector3 Gpos,dpos;
//       Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);
//       double z1,z2,x1,x2,y1,y2;
//       z1=-20;z2=20;
//       blc2->XYPosatZ(z1+Gpos.z(),x1,y1);
//       blc2->XYPosatZ(z2+Gpos.z(),x2,y2);
//       TLine line;
//       line.SetLineColor(itr+2);
//       line.DrawLine(z1,y1,z2,y2);
//       for(int ih=0;ih<blc2->nhit();ih++)
// 	{
// 	  ChamberLikeHit *hit=blc2->hit(ih);
// 	  TMarker blt_m;	  
// 	  if(hit->xy()==1)
// 	    {
// 	      blt_m.SetMarkerStyle(20);
// 	      blt_m.SetMarkerSize(0.2);
// 	      blt_m.SetMarkerColor(3  );
// 	      blt_m.DrawMarker(hit->wz(),hit->wy() );
// 	    }
// 	}

//     }

  for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
    {
      LocalTrack *blc=blTrackMan->trackBLC(itr);
      BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
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

  

  canvas1->Update();


  
  for(int i=0;i<blTrackMan->ntrackBLC();i++)
    { 
      double gx,gy,gz,dx,dy,dz;
      BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
      BLwireman->GetGParam(18,gx,gy,gz,dx,dy,dz);
      TVector3 Gpos,dpos;
      Gpos.SetXYZ(gx,gy,gz);  dpos.SetXYZ(dx,dy,dz);     
      LocalTrack *blc=blTrackMan->trackBLC(i);
      TVector3 x1,x2;
      double xp1,xp2,yp1,yp2;
      //      blc2->XYPosatZ(-60-Gpos.z(),xp1,yp1);x1.SetXYZ(xp1,yp1,-60-Gpos.z());
      //blc2->XYPosatZ(60-Gpos.z(),xp2,yp2);x2.SetXYZ(xp2,yp2,60-Gpos.z());
      //      x1=x1+Gpos; x2=x2+Gpos;
      blc->XYPosatZ(-60,xp1,yp1);x1.SetXYZ(xp1,yp1,-60);
      blc->XYPosatZ(60,xp2,yp2);x2.SetXYZ(xp2,yp2,60);
      x1=x1; x2=x2;
      TLine line;
      TMarker bl_mk;
      bl_mk.SetMarkerStyle(20);
      bl_mk.SetMarkerSize(0.8);
      canvas1->cd(1);	
      line.DrawLine(x1.x(),x1.y(),x2.x(),x2.y());
      bl_mk.DrawMarker(x1.x(),x1.y());
      canvas1->cd(2);	
      line.DrawLine(x1.z(),x1.y(),x2.z(),x2.y());
      bl_mk.DrawMarker(x1.z(),x1.y());
    }


  for(int layer=1;layer<=15;layer++)
    {
      for(int nchit=0;nchit<cdsMan->nCDC(layer);nchit++)
	{
	  double dl=cdsMan->CDC(layer,nchit)->dl();		
	  double resi=cdsMan->CDC(layer,nchit)->resl();		
	  int wire=cdsMan->CDC(layer,nchit)->wire();		
	  std::cout<<"layer= "<<layer<<" wire= "<<wire<<" resid= "<<resi<<" dl= "<<dl<<std::endl;	  
	}
    }

  
    for(int numslayer=1;numslayer<=7;numslayer++)
      {
	for(int i=0;i<trackMan->nCluster(numslayer);i++)
	  {
	    HitCluster *cluster=trackMan->Cluster(numslayer,i);
	    double Meanphi=cluster->Meanphi();		
	    double Meanradius=cluster->Meanradius();		
	    
	    TVector3 cls_pos;
	    cls_pos.SetXYZ( Meanradius*cos(Meanphi/180*TMath::Pi()),
			    Meanradius*sin(Meanphi/180*TMath::Pi()),
			    0);
	    //		seg=track->CDHHit()->seg();		
	    TMarker cls_mk;
	    cls_mk.SetMarkerStyle(24);
	    cls_mk.SetMarkerColor(1);
	    cls_mk.SetMarkerSize(1.2);
	    canvas1->cd(1);
	    cls_mk.DrawMarker(cls_pos.x(),cls_pos.y());
	    //	    std::cout<<" ncl= "<<cluster->nHit()<<" phi "<<cluster->Meanphi()<<std::endl; 
	  }
      }
      
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
	for(int i=0;i<blTrackMan->ntrackBLC2();i++)
	  { 
	    LocalTrack *blc2=blTrackMan->trackBLC2(i);
	    TVector3 x1,x2;
	    double xp1,xp2,yp1,yp2,tan;
	    blc2->XYPosatZ(-74,xp1,yp1);x1.SetXYZ(xp1,yp1,-74);
	    blc2->XYPosatZ(74,xp2,yp2);x2.SetXYZ(xp2,yp2,74);
	    tan=148.0/sqrt( (xp1-xp2)*(xp1-xp2)+(yp1-yp2)*(yp1-yp2) );
	    double rho_b=-900.0/tan/0.5/2.99;
	    double phi=-(verb2.z()+74)/(rho_b*tan);
	    //  MathTools *tool=new MathTools();
	    //double phi0=tool->CalcDeg(xp1,yp1)*TMath::DegToRad();
	    // delete tool;
	    double d_r=sqrt(xp1*xp1+yp1*yp1);
	    double vx,vy;
	    //std::cout<<" d_r "<<d_r<<" phi0 "<<phi0<<" rho_b "<<rho_b<<" tan "<<tan<<"  phi_b "<<phi<<std::endl;
	    //vx=d_r*cos(phi0)+rho_b*(cos(phi0)-cos(phi0+phi));
	    //vy=d_r*sin(phi0)+rho_b*(sin(phi0)-sin(phi0+phi));
	    //	    std::cout<<"blc2 vertex = "<<vx<<" "<<vy<<" "<<
	    //     verb2.z()<<std::endl;
	  }
	
	std::cout<<"vertex dis = "<<(verb1-verb2).Mag()<<std::endl;
      }
	
    std::cout<<" ntrack= "<<trackMan->nTrack()<<" goodtrack= "<<trackMan->nGoodTrack()<<std::endl; 
    double arho,ax_c,ay_c,aPt;
    double aparam[5];
    double paramtmp3[5];
    double paramtmp4[5];
    TF1 *func_y;	
    for(int i=0;i<trackMan->nTrack();i++)
      {
	CDSTrack *track=trackMan->Track(i);
	bool goodflag=track->GoodFlag();
	if(!goodflag) continue;
	track->GetParameters(aparam);
	track->GetTmpParameters(3,paramtmp3);
	track->GetTmpParameters(4,paramtmp4);

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
		std::cout<<"layer= "<<numlayer<<" resid= "<<resi<<" dl= "<<dl<<std::endl;
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

	std::cout<<"paramtmp[3]= "
		 <<paramtmp3[0]<<" "
		 <<paramtmp3[1]<<" "
		 <<paramtmp3[2]<<" "
		 <<paramtmp3[3]<<" "
		 <<paramtmp3[4]<<" "
		 <<std::endl;

	std::cout<<"paramtmp[4]= "
		 <<paramtmp4[0]<<" "
		 <<paramtmp4[1]<<" "
		 <<paramtmp4[2]<<" "
		 <<paramtmp4[3]<<" "
		 <<paramtmp4[4]<<" "
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



  canvas1->Update();
  canvas1->SaveAs("tmpdis.ps");
  bool status = disp->Wait();
  if( !status ) return status;

  }

  //  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  trackMan->Clear();
  gFile->Write();
  gFile->Close();
  
  gSystem->Exit(1);
  gROOT->GetApplication()->Run();
  //  theApp.Run();

  return 0;

}
