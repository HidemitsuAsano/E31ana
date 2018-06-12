#include <iostream>
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
#include  <TGraphErrors.h> 

#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
#include "Display3D.h"

//const double MassThPPi = 0.8; // GeV/c2
//const double TOFOffset = 3; // nsec

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

int PID( double tof, double mom ) // 0:unknown, 1:pi+, 2;pi-, 3:proton
{
  //TF1 *fun1 = new TF1("fun1",Form("%lf/sqrt( pow((x+%lf),2)-1)",MassThPPi,TOFOffset),-2,14);
  int val=0;
  if( tof<TOFOffset ){
  }
  if( 0<mom ){
    if( mom < MassThPPi/sqrt(pow(tof+TOFOffset,2)-1) ){
      val = 1;
    }
    else{
      val = 3;
    }
  }
  else{
    if( fabs(mom) < MassThPPi/sqrt(pow(tof+TOFOffset,2)-1) ){
      val = 2;
    }
    else{
      val = 0;
    }
  }

  //std::cout << " tof:" << tof << " mom:" << mom << " type:" << val << std::endl;
  return val;
}

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

  int runnum = 3036;
  
  gSystem->Load("libPhysics.so");
  gSystem->Load( "./lib/libAll.so" );
  

  ConfMan *conf = new ConfMan( confFile );
    

  //  MathTools *mtool = new MathTools();

  TFile *f;
  TTree *evtree;
  f = new TFile( inFile.c_str() );

  evtree = (TTree*)f->Get( "EventTree" );
  

  CDSHitMan *cdsMan = 0;
  BeamLineHitMan *blMan = 0;
  BeamLineTrackMan *blTrackMan = 0;
  EventHeader *head = 0;
  CDSTrackingMan *trackMan = 0;
  
  TFile *fout;
  fout = new TFile( outFile.c_str(),"recreate" );
  //##### Creating Hist ###########
   
  for(int layer=1;layer<=NumOfCDCLayers;layer++)
    {
      new TH1F( Form("CDCresid%d",layer),Form("CDCresid%d",layer) ,200,-0.1,0.1);
      new TH1F( Form("CDCresid_pi%d",layer),Form("CDCresid_pi%d",layer) ,200,-0.1,0.1);
      new TH1F( Form("CDCresid_proton%d",layer),Form("CDCresid_proton%d",layer) ,200,-0.1,0.1);
      new TH2F( Form("CDCdt_resid%d",layer),Form("CDCdt_resid%d",layer) ,
		300,-20,280,200,-0.1,0.1);
      new TH1F( Form("CDCcellangle%d",layer),Form("CDCcellangle%d",layer) ,
		720,0,360);
      new TH1F( Form("CDCangle%d",layer),Form("CDCangle%d",layer) ,
		720,0,360);
	  for(int n=1;n<=6;n++)
	    {
	  new TH2F( Form("CDCdt_resid_angle%d_%d",n,layer),
		    Form("CDCdt_resid_angle%d_%d",n,layer) ,
		    300,-20,280,200,-0.1,0.1);
	  new TH2F( Form("CDCdt_resid_cell%d_%d",n,layer),
		    Form("CDCdt_resid_cell%d_%d",n,layer) ,
		    300,-20,280,200,-0.1,0.1);
	    }
      for(int wire=1;wire<=NumOfCDCWiresInLayer[layer-1];wire++)
	{
	  new TH1F( Form("CDCdt%d_%d",layer,wire),Form("CDCdt%d_wire",layer,wire) ,300,-60,240);
	  new TH1F( Form("CDCresid%d_%d",layer,wire),Form("CDCresid%d_wire",layer,wire) ,200,-0.2,0.2);
	  new TH2F( Form("CDCdt_resid%d_%d",layer,wire),Form("CDCdt_resid%d_%d",layer,wire) ,
		300,-20,280,200,-0.1,0.1);
	  new TH2F( Form("CDCdt_dl%d_%d",layer,wire),Form("CDCdt_dl%d_%d",layer,wire) ,
		300,-20,280,1000,-0.1,0.9);

	}
    }

   
  for(int layer=1;layer<=NumOfBLCLayers;layer++)
    {
      new TH1F( Form("BLC2resid%d",layer),Form("BLC2resid%d",layer) ,200,-0.1,0.1);
      new TH2F( Form("BLC2dt_resid%d",layer),Form("BLC2dt_resid%d",layer) ,
		200,-50,150,200,-0.1,0.1);
      for(int wire=1;wire<=NumOfBLCWiresInLayer;wire++)
	{
	  new TH1F( Form("BLC2dt%d_%d",layer,wire),Form("BLC2dt%d_%d",layer,wire) ,200,-50,150);
	  new TH1F( Form("BLC2resid%d_%d",layer,wire),Form("BLC2resid%d_wire",layer,wire) ,400,-0.2,0.2);
	  new TH2F( Form("BLC2dt_resid%d_%d",layer,wire),Form("BLC2dt_resid%d_%d",layer,wire) ,
		200,-50,150,400,-0.2,0.2);
	  new TH2F( Form("BLC2dt_dl%d_%d",layer,wire),Form("BLC2dt_dl%d_%d",layer,wire) ,
		200,-50,150,1000,-0.1,0.9);
	  new TH2F( Form("BLC2dt_dl2_%d_%d",layer,wire),
		    Form("BLC2dt_dl2_%d_%d",layer,wire) ,
		    200,-50,150,1000,-0.5,0.5);
	}
    }

    
  for(int layer=1;layer<=NumOfBLCLayers;layer++)
    {
      new TH1F( Form("BLC1resid%d",layer),Form("BLC1resid%d",layer) ,200,-0.1,0.1);
      new TH2F( Form("BLC1dt_resid%d",layer),Form("BLC1dt_resid%d",layer) ,
		200,-50,150,400,-0.2,0.2);
      for(int wire=1;wire<=NumOfBLCWiresInLayer;wire++)
	{
	  new TH1F( Form("BLC1dt%d_%d",layer,wire),Form("BLC1dt%d_%d",layer,wire) ,200,-50,150);
	  new TH1F( Form("BLC1resid%d_%d",layer,wire),Form("BLC1resid%d_wire",layer,wire) ,400,-0.2,0.2);
	  new TH2F( Form("BLC1dt_resid%d_%d",layer,wire),Form("BLC1dt_resid%d_%d",layer,wire) ,
		200,-50,150,200,-0.1,0.1);
	  new TH2F( Form("BLC1dt_dl%d_%d",layer,wire),
		    Form("BLC1dt_dl%d_%d",layer,wire) ,
		200,-50,150,1000,-0.1,0.9);
	  new TH2F( Form("BLC1dt_dl2_%d_%d",layer,wire),
		    Form("BLC1dt_dl2_%d_%d",layer,wire) ,
		    200,-50,150,1000,-0.5,0.5);
	}
    }
 

  //####  set Branch ###############

  TH1F *h1;
  TH2F *h2;

  bool AllData=false;
  TObjArray *objarr=evtree->GetListOfBranches();
  if(objarr->FindObject("CDSHitMan")!=0 && objarr->FindObject("CDSTrackingMan")!=0)
    {
      evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
      evtree->SetBranchAddress( "CDSTrackingMan", &trackMan );
      if(objarr->FindObject("BeamLineHitMan")!=0 && 
	 objarr->FindObject("BeamLineTrackMan")!=0 && 
	 objarr->FindObject("EventHeader")!=0 )
      {
	evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
	evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
	evtree->SetBranchAddress( "BeamLineTrackMan", &blTrackMan );
	evtree->SetBranchAddress( "EventHeader", &head );
	evtree->SetBranchAddress( "CDSTrackingMan", &trackMan );
	AllData=true;
      }
    }
  else 
    {
      std::cout<<"error of EventTree branch"<<std::endl;
      return 0;
    }



  int nev = evtree->GetEntries();
  std::cout << " AllEntry:" << nev << std::endl;
  
  for( int iev=0; iev<nev; iev++ ){
    if( iev%1000 == 0 )
      std::cout << " Event : " << iev << std::endl;
    evtree->GetEvent(iev);
    //cdsMan->Calc(conf);
    //cdsMan->CheckContainerSize();
      //blMan->CheckContainerSize();

    if(AllData)
      {
	int nT0=0;
	for( int i=0; i<blMan->nT0(); i++ ){
	  if( blMan->T0(i)->CheckRange() ) nT0++;
	}
	if( nT0!=1 ) continue;
	double ctmT0;
	for( int i=0; i<blMan->nT0(); i++ ){
	  if( blMan->T0(i)->CheckRange() ){
	    ctmT0 = blMan->T0(i)->ctmean();
	  }
	}
	
	int nBHD=0;
	for( int i=0; i<blMan->nBHD(); i++ ){
	  if( blMan->BHD(i)->CheckRange() ) nBHD++;
	}
	//       if( nBHD!=1 ) continue;
	double ctmBHD;
	for( int i=0; i<blMan->nBHD(); i++ ){
	  if( blMan->BHD(i)->CheckRange() ){
	    ctmBHD = blMan->BHD(i)->ctmean();
	  }
	}
	double tofBHDT0 = ctmT0-ctmBHD;
	
	int pid_beam;//0:pi 1:K 3:else
	
	if(-2<tofBHDT0 && tofBHDT0<2) pid_beam=0;
	else if(2<tofBHDT0 && tofBHDT0<5) pid_beam=1;
	else pid_beam=3;
	//if(pid_beam!=0) continue;
	//     std::cout << " nTrack:" << trackMan->nTrack() << std::endl;
      }

    int GoodTrack=trackMan->nGoodTrack();
    //    std::cout << " nGoodTrack:" << GoodTrack << std::endl;    

    for( int it=0; it<trackMan->nGoodTrack(); it++ ){
      CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );

      double chi = track->Chi();
      double param[5];
      track->GetParameters(param);
      //      if(param[2]==0  ) continue;
      double drho=param[0], phi0=param[1], rho=1./param[2], dz=param[3], tlam=param[4];
      double mom = track->Momentum();
      if(fabs(param[3])>15) continue;

      //##PID#############
      /*
      if(!track->CDHFlag()) continue;
      HodoscopeLikeHit *cdh = track->CDHHit();
      double tof = cdh->ctmean()-ctmT0;
      int ptype = PID( tof, mom );
      */
      //#####################
      for(int layer=1;layer<=NumOfCDCLayers;layer++)
	{
	  for(int nhit=0;nhit<track->nTrackHit(layer);nhit++)
	    {

	      CDCHit *cdc=track->TrackHit(layer,nhit);
	      double resi=cdc->resl();
	      double dt=cdc->dt();
	      double dl=cdc->dl();
	      double dlr=dl-resi;
	      int wire=cdc->wire();
	      double angle=cdc->phi();
	      TVector3 wpos,wposp,hpos,tmp,west,dir;
	      double dist;
	      wpos.SetXYZ(cdc->wx(),cdc->wy(),cdc->wz());
	      wposp.SetXYZ(cdc->wxp(),cdc->wyp(),cdc->wzp());
	      dir=(wposp-wpos).Unit();
	      hpos.SetXYZ(cdc->x(),cdc->y(),cdc->z());
	      if( (layer>=4 &&layer<=7)  || (layer>=10 &&layer<=13) )
		{
		  track->PointToLine( hpos,wpos,dir,dist,west );
		}
	      else west=wpos;
	      tmp=hpos-west;
	      MathTools *tool;
	      double cangle=tool->CalcDeg(tmp.x(),tmp.y());
	      cangle-=(angle-90);
	      if(cangle<0) cangle+=360;
	      if(cangle>360) cangle-=360;

	      int canglenum,anglenum;
	      if(angle<60) anglenum=1;
	      else  if(angle<120) anglenum=2;
	      else  if(angle<180) anglenum=3;
	      else  if(angle<240) anglenum=4;
	      else  if(angle<300) anglenum=5;
	      else  if(angle<360) anglenum=6;

	      if(cangle<30) canglenum=1;
	      else  if(cangle<90) canglenum=2;
	      else  if(cangle<150) canglenum=3;
	      else  if(cangle<210) canglenum=4;
	      else  if(cangle<270) canglenum=5;
	      else  if(cangle<330) canglenum=6;
	      else  if(cangle<360) canglenum=1;

	      h1 = (TH1F*)gFile->Get(Form("CDCcellangle%d",layer)); h1->Fill( cangle );		
	      h1 = (TH1F*)gFile->Get(Form("CDCangle%d",layer)); h1->Fill( angle );		
	      h1 = (TH1F*)gFile->Get(Form("CDCdt%d_%d",layer,wire)); h1->Fill( dt );		
	      h1 = (TH1F*)gFile->Get(Form("CDCresid%d",layer)); h1->Fill( resi );		

	      /*
	      if(ptype==1||ptype==2)
		{
		  h1 = (TH1F*)gFile->Get(Form("CDCresid_pi%d",layer)); h1->Fill( resi );		
		}
	      else if(ptype==3)
		{
		  h1 = (TH1F*)gFile->Get(Form("CDCresid_proton%d",layer)); h1->Fill( resi );		
		}
	      */
	      if(fabs(resi)<0.05) 
		{	
		  h2 = (TH2F*)gFile->Get(Form("CDCdt_resid%d",layer)); h2->Fill( dt,resi ); 
		  h2 = (TH2F*)gFile->Get(Form("CDCdt_resid%d_%d",layer,wire)); h2->Fill( dt,resi );		
		}
	      h1 = (TH1F*)gFile->Get(Form("CDCresid%d_%d",layer,wire)); h1->Fill( resi );		

	      h2 = (TH2F*)gFile->Get(Form("CDCdt_dl%d_%d",layer,wire)); h2->Fill( dt,dlr );		
	      
	      h2 = (TH2F*)gFile->Get(Form("CDCdt_resid_angle%d_%d",anglenum,layer,wire)); h2->Fill( dt,resi );		
	      h2 = (TH2F*)gFile->Get(Form("CDCdt_resid_cell%d_%d",canglenum,layer,wire)); h2->Fill( dt,resi );		
	    }
	}
    }
    
    if(AllData)
      {
	//#########BLHitMan #################
	for(int layer=0;layer<=NumOfBLCLayers;layer++)
	  {
	    for(int ih=0;ih<blMan->nBLC1(layer);ih++)
	      {
		ChamberLikeHit *hit=blMan->BLC1(layer,ih);
		double dt=hit->dt();
		int wire=hit->wire();
		h1 = (TH1F*)gFile->Get(Form("BLC1dt%d_%d",layer,wire)); h1->Fill( dt );		
	      }
	    for(int ih=0;ih<blMan->nBLC2(layer);ih++)
	      {
		ChamberLikeHit *hit=blMan->BLC2(layer,ih);
		double dt=hit->dt();
		int wire=hit->wire();
		h1 = (TH1F*)gFile->Get(Form("BLC2dt%d_%d",layer,wire)); h1->Fill( dt );		
	      }
	  }
  //######## BLCTrack  #################
	for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
	  {
	    LocalTrack *blc2=blTrackMan->trackBLC2(itr);
	    if(blc2->chi2all()>30) continue;
	    for(int ih=0;ih<blc2->nhit();ih++)
	      {
		ChamberLikeHit *hit=blc2->hit(ih);
		
		double resi= hit->resl();
		double dt=hit->dt();
		double dl=hit->dl();
		//	  double dlr=dl-resi;
		int cid=hit->cid();
		int wire=hit->wire();
		int layer=hit->layer();
		int blcnum=0;
		if(cid==17/*BLC1*/) blcnum=1;
		else 	if(cid==18/*BLC1*/) blcnum=2;
		
		
		double hpos,wpos,dlr;
		if(hit->xy()==0) 
		  {
		    hpos=hit->x(); wpos=hit->wx();dlr=(hpos-wpos);
		  }  
		else if(hit->xy()==1) 
		  {
		    hpos=hit->y(); wpos=hit->wy();dlr=(hpos-wpos);
		  }  
		
		//  if(dlr<0) continue;
		double fdlr=fabs(dlr);
		if(dlr<0 && fdlr<0.28) 
		  {
		    //	  h1 = (TH1F*)gFile->Get(Form("BLC%ddt%d_%d",blcnum,layer,wire)); h1->Fill( dt );		
		    h1 = (TH1F*)gFile->Get(Form("BLC%dresid%d",blcnum,layer)); h1->Fill( resi );		
		    h2 = (TH2F*)gFile->Get(Form("BLC%ddt_resid%d",blcnum,layer));
		    if(dt>120) h2->Fill( dt,dl-0.255 );
		    else if(fdlr<0.25) h2->Fill( dt,resi );
		    else  h2->Fill( dt,dl-0.255 );
		    h1 = (TH1F*)gFile->Get(Form("BLC%dresid%d_%d",blcnum,layer,wire)); h1->Fill( resi );		
		    h2 = (TH2F*)gFile->Get(Form("BLC%ddt_resid%d_%d",blcnum,layer,wire)); 
		    if(dt>120) h2->Fill( dt,dl-0.255 );
		    else if(fdlr<0.25) h2->Fill( dt,resi );
		    else  h2->Fill( dt,dl-0.255 );
		    
		    h2 = (TH2F*)gFile->Get(Form("BLC%ddt_dl%d_%d",blcnum,layer,wire)); h2->Fill( dt,fdlr );		
		  }
		h2 = (TH2F*)gFile->Get(Form("BLC%ddt_dl2_%d_%d",blcnum,layer,wire)); h2->Fill( dt,dlr );		
		
		
	      }
	  }
      }
    
    
  }
  
  //########  Calc Resid  ###############//
  
  TF1 *dgaus=new TF1("dgaus","gaus(0)+gaus(3)");
  TF1 *gaus1=new TF1("gaus1","gaus(0)");
  gaus1->SetParameters(100,0,0.02);
  dgaus->SetLineColor(2);
  double nlayer[15]={0};
  double resid[15]={0};
  double ex[15]={0};
  double ey[15]={0};
  for(int layer=1;layer<=15;layer++)
    {
      //      cout<<"Layer "<<layer<<endl;
      TH1F *h1=(TH1F*)gFile->Get(Form("CDCresid%d",layer) );
      h1->Fit("gaus1");
      double par[3],apar[6]; 
      gaus1->GetParameters(par);

      if(AllData)
	{
      dgaus->SetParameters(par[0],par[1],par[2]/2.0,par[0]/10.0,par[1],par[2]*2);
      //dgaus->SetParameters(500,0,0.02,100,0,0.05);
      h1->Fit("dgaus");
      dgaus->GetParameters(apar);
      nlayer[layer-1]=layer;
      resid[layer-1]=apar[2];
      ey[layer-1]=dgaus->GetParError(2);
	}
      else
	{
	  nlayer[layer-1]=layer;
	  resid[layer-1]=par[2];
	  ey[layer-1]=gaus1->GetParError(2);
	}

    }
  TGraphErrors *g1=new TGraphErrors(15,nlayer,resid,ex,ey);
  g1->SetLineColor(2);
  //  g1->Draw("ALP");
  g1->SetName("CDCResid");
  g1->Write();

  //##################
  
  gFile->Write();
  gFile->Close();
  
  gSystem->Exit(1);
  gROOT->GetApplication()->Run();
  //  theApp.Run();

  return 0;
}

