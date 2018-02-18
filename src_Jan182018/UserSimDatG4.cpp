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
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "SimDataMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
#include "Display3D.h"

  //#######Need knucl lib file!! Set your dir!! ###
#include "/home/sakuma/work/ana/geant/knucl3/include/KnuclRootData.h"
  //#######################################

#define DEBUG 0
#define TRACKING 1

const double ADC_TH_CVC = 0.5; // MeV
const double ADC_TH_TOFstop = 1.0; // MeV
const double ADC_TH_NC = 2.0; // MeV

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



  //  TROOT root( "GUI", "GUI" );
  // TApplication theApp( "App", &argc, argv );
  // gROOT->SetStyle( "Plain" );
  // gROOT->cd();

  TRint *theApp = new TRint( "theApp", &argc, argv );
  //TROOT root( "GUI", "GUI" );
  //TApplication theApp( "App", &argc, argv );
  //gROOT->SetStyle( "Plain" );
  gROOT->cd();

  const int residnum=250; 

  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("./lib/libAll.so");
  //#######Need knucl lib file!! Set your dir!! ###
  gSystem->Load("~/work/ana/geant/knucl3/libKnuclRootData.so");
  //#######################################

  /*** assign input file & call tree ***/

  TFile *f = new TFile( inFile.c_str() );
  TTree *evtree = (TTree*)f->Get("tree");
  
  /*** conf file for new parameter ***/
  
  ConfMan *conf = new ConfMan(confFile);
  conf->Initialize();

  Header* header = 0;
  DetectorData *detectorData = 0;
  MCData* mcData = 0;
  
  
  /*** declaration of classes ***/
  
  evtree->SetBranchAddress( "Header", &header );
  evtree->SetBranchAddress( "DetectorData", &detectorData );
  evtree->SetBranchAddress( "MCData", &mcData );

  
  /*** make output file & make tree ***/
 
  TFile *of = new TFile(outFile.c_str(),"recreate");
  
  TTree *otree = new TTree( "EventTree", "EventTree");
  //--- copy input to output ---//
  
  //otree->Branch( "header", &header );
  //otree->Branch( "detectorData", &detectorData );
  otree->Branch( "MCData", &mcData );
  //--- attach new branch on output ---//

  CDSHitMan* cdsMan = 0;
  CDSTrackingMan* cdsTrackMan = new CDSTrackingMan();
  otree->Branch( "CDSHitMan", &cdsMan );
  otree->Branch( "CDSTrackingMan", &cdsTrackMan );
  BeamLineHitMan* blMan = 0;
  otree->Branch( "BeamLineHitMan", &blMan );
  BeamLineTrackMan* blTrackMan = new BeamLineTrackMan();
  otree->Branch( "BeamLineTrackMan", &blTrackMan );


  int nev = evtree->GetEntries();
  int fev = 0;
  int tev = 0;
  //nev = 1000;
  std::cerr<<"# of events = "<<nev<<std::endl;
  for(int iev = 0; iev<nev; iev++ ){
    evtree->GetEvent(iev);
    SimDataMan *simMan = new SimDataMan(conf);
   
    
    if( iev%1000==0 )
      std::cout << " Event : " << iev << std::endl;

    //###########################//
    //### check detector hits ###//
    //###########################//
    // check hit size
    int nT0 = 0;
    int nCDH = 0;
    int nBLC1 = 0;
    int nBLC2 = 0;
    int nCDC = 0;
    int nNC = 0;
    int nTOFstop = 0;
    int nPC = 0;
    int nCVC = 0;
    double time0 = 0;
    for( int i=0; i<detectorData->detectorHitSize(); i++ ){
      DetectorHit *detectorHit = detectorData->detectorHit(i);
      int detectorID = detectorHit->detectorID();
      if( detectorID==CID_T0 ){
	nT0++;
	time0 = detectorHit->tdc();
      }
      else if( detectorID==CID_CDH ) nCDH++;
      else if( detectorID==CID_BLC1 ) nBLC1++;
      else if( detectorID==CID_BLC2 ) nBLC2++;
      else if( detectorID==CID_CDC ) nCDC++;
      else if( detectorID==CID_NC && detectorHit->adc()>ADC_TH_NC ) nNC++;
      else if( detectorID==CID_TOFstop && detectorHit->adc()>ADC_TH_TOFstop ) nTOFstop++;
      else if( detectorID==CID_PC ) nPC++;
      else if( detectorID==CID_CVC && detectorHit->adc()>ADC_TH_CVC ) nCVC++;
    }
#if DEBUG
    std::cerr<<"###### event number  = "<<iev<<std::endl;
    std::cerr<<"# of T0 hits   = "<<nT0<<std::endl;
    std::cerr<<"# of CDH hits  = "<<nCDH<<std::endl;
    std::cerr<<"# of BLC1 hits = "<<nBLC1<<std::endl;
    std::cerr<<"# of BLC2 hits = "<<nBLC2<<std::endl;
    std::cerr<<"# of CDC hits  = "<<nCDC<<std::endl;
    std::cerr<<"# of NC hits  = "<<nNC<<std::endl;
    std::cerr<<"# of TOFstop hits  = "<<nTOFstop<<std::endl;
    std::cerr<<"# of PC hits  = "<<nPC<<std::endl;
    std::cerr<<"# of CVC hits  = "<<nCVC<<std::endl;
#endif

    //if(!nT0==1) continue;
    //if(nCDC==0) continue;
    //if(nCDH==0) continue;
    //if(Ncdh<2) continue;
#if 1
    if ( !nT0 ) continue;
#else
    if(nT0==0)
      {
	simMan->SetT0Hit( 3, 0, 0, 0 );
	time0=0;
      }
#endif

    //##########################//
    //### fill detector hits ###//
    //##########################//
    DetectorHit *detectorHit = 0;
    double val0=0;
    double val1=0;
    double val2=0;
    for(int i=0; i<detectorData->detectorHitSize(); i++ ){
      detectorHit = detectorData->detectorHit(i);
      //--- T0 ---//
      if(detectorHit->detectorID()==CID_T0){
	Int_t seg = detectorHit->channelID()+1;
	val0=detectorHit->tdc()-time0;
	val1=detectorHit->adc();
	val2=detectorHit->pos().y();
	if(val1==0) val1=3.;//MeV
	simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
      }
      //--- CDH ---//
      if(detectorHit->detectorID()==CID_CDH){
	Int_t seg = detectorHit->channelID()+1;
	val0=detectorHit->tdc()+time0;
	val1=detectorHit->adc();
	val2=detectorHit->pos().z();
	if(val1==0) val1=6.;//MeV
	simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
      }
      //--- BLC ---//
      if(detectorHit->detectorID()==CID_BLC1 ||
	 detectorHit->detectorID()==CID_BLC2 ){
	Int_t layer = detectorHit->layerID()+1;
	Double_t wire = detectorHit->channelID()+1;
	int tag;
	if(detectorHit->detectorID()==CID_BLC1) tag=0;
	if(detectorHit->detectorID()==CID_BLC2) tag=1;
	val0=fabs(detectorHit->dx()/10.);
	val1=0;//time offset
	simMan->SetChamberLikeHit( detectorHit->detectorID(), layer, wire, val0, val1 );
      }
      //--- CDC ---//
      if(detectorHit->detectorID()==CID_CDC){
	Int_t layer = detectorHit->layerID()+1;
	Double_t wire = detectorHit->channelID()+1;
	val0=fabs(detectorHit->dx()/10.);
	val1=0;//time offset
	simMan->SetChamberLikeHit( detectorHit->detectorID(), layer, wire, val0, val1 );
      }
      //--- NC ---//
      if(detectorHit->detectorID()==CID_NC && detectorHit->adc()>ADC_TH_NC){
	Int_t layer = detectorHit->layerID()+1;
	Int_t channel = detectorHit->channelID()+1;
	Int_t seg = 16*(layer-1)+channel;
	val0=detectorHit->tdc()+time0;
	val1=detectorHit->adc();
	val2=detectorHit->pos().z();
	if(val1==0) val1=6.;//MeV
	simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
      }
      //--- TOFstop ---//
      if(detectorHit->detectorID()==CID_TOFstop && detectorHit->adc()>ADC_TH_TOFstop){
	Int_t seg = detectorHit->channelID()+1;
	val0=detectorHit->tdc()+time0;
	val1=detectorHit->adc();
	val2=detectorHit->pos().z();
	if(val1==0) val1=6.;//MeV
	simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
      }
      //--- PC ---//
      if(detectorHit->detectorID()==CID_PC){
	Int_t seg = detectorHit->channelID()+1;
	val0=detectorHit->tdc()+time0;
	val1=detectorHit->adc();
	val2=detectorHit->pos().z();
	if(val1==0) val1=6.;//MeV
	simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
      }
      //--- CVC ---//
      if(detectorHit->detectorID()==CID_CVC && detectorHit->adc()>ADC_TH_CVC){
	Int_t seg = detectorHit->channelID()+1;
	val0=detectorHit->tdc()+time0;
	val1=detectorHit->adc();
	val2=detectorHit->pos().z();
	if(val1==0) val1=6.;//MeV
	simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
      }
    } // for(int i=0; i<detectorData->detectorHitSize(); i++ ){

    cdsMan = simMan->GetCDSHit();
    blMan = simMan->GetBeamLineHit();

#if DEBUG
    std::cout<<"######Filled data#############"<<std::endl;
    for(int i=0; i<blMan->nT0(); i++){
      HodoscopeLikeHit *hod = blMan->T0(i);
      std::cerr<<"T0 : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
    }
    for(int i=0; i<cdsMan->nCDH(); i++){
      HodoscopeLikeHit *hod = cdsMan->CDH(i);
      std::cerr<<"CDH : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
    }
    for(int i=1; i<=NumOfBLCLayers; i++){
      for(int j=0; j<blMan->nBLC1(i); j++){
	ChamberLikeHit *chm = blMan->BLC1(i, j);
	std::cerr<<"BLC1 : "<<chm->hid()<<" "<<chm->layer()<<" "<<chm->wire()<<" "<<chm->dl()<<std::endl;
      }
    }
    for(int i=1; i<=NumOfBLCLayers; i++){
      for(int j=0; j<blMan->nBLC2(i); j++){
	ChamberLikeHit *chm = blMan->BLC2(i, j);
	std::cerr<<"BLC2 : "<<chm->hid()<<" "<<chm->layer()<<" "<<chm->wire()<<" "<<chm->dl()<<std::endl;
      }
    }
    for(int i=1; i<=NumOfCDCLayers; i++){
      for(int j=0; j<cdsMan->nCDC(i); j++){
	CDCHit *cdc = cdsMan->CDC(i, j);
	std::cerr<<"CDC : "<<cdc->hid()<<" "<<cdc->layer()<<" "<<cdc->wire()<<" "<<cdc->dl()<<std::endl;
      }
    }
    for(int i=0; i<blMan->nNC(); i++){
      HodoscopeLikeHit *hod = blMan->NC(i);
      std::cerr<<"NC : "<<hod->hid()<<" "<<hod->seg()<<" ("<<hod->seg()/16+1<<","<<hod->seg()%16<<") "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
    }
    for(int i=0; i<blMan->nTOF(); i++){
      HodoscopeLikeHit *hod = blMan->TOF(i);
      std::cerr<<"TOFstop : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
    }
    for(int i=0; i<blMan->nPC(); i++){
      HodoscopeLikeHit *hod = blMan->PC(i);
      std::cerr<<"PC : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
    }
    for(int i=0; i<blMan->nCV(); i++){
      HodoscopeLikeHit *hod = blMan->CV(i);
      std::cerr<<"CVC : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
    }
#endif

#if TRACKING
    //##########Tracking #############//
    cdsTrackMan->Execute(cdsMan,conf);
    blTrackMan->DoTracking(blMan,conf);
#if DEBUG
    std::cerr<<"ev-num = "<<iev<<" track-num = "<<cdsTrackMan->nGoodTrack()<<std::endl;
#endif
    if ( cdsTrackMan->nGoodTrack()>1 ) tev++;
    for(int n=0;n<cdsTrackMan->nGoodTrack();n++)
      {
	int GoodTrack=cdsTrackMan->GoodTrackID(n);
	cdsTrackMan->CalcVertex_beam(GoodTrack,blTrackMan,conf);
	CDSTrack *track=cdsTrackMan->Track(GoodTrack);
	for(int layer=1;layer<=15;layer++)
	  {
	    for(int m=0;m<track->nTrackHit(layer);m++)
	      {
		CDCHit *hit=track->TrackHit(layer,m);
		//		std::cout<<" Trackhit l : "<<layer<<" "<<"resi : "<<hit->resl()<<" dl : "<<hit->dl()<<
		//		  " "<<"hitpos : "<<hit->x()<<" "<<hit->y()<<" "<<hit->z()<<std::endl;
	      }
	  }
	double param[5];
	track->GetParameters(param);
	//	std::cout<<"track posz : "<<param[3]<<" chi "<<track->Chi()<<std::endl;
      }
    //###############################//
#endif

    otree->Fill();
    fev++;

    //######clear###########
    cdsMan->Clear();
    cdsTrackMan->Clear(); 
    blMan->Clear(); 
    blTrackMan->Clear();
    simMan->Clear();
  
  } // for(int iev = 0; iev<nev; iev++ ){

  std::cerr<<"filled events = "<<fev<<std::endl;
  std::cerr<<"2 or more tracks events = "<<tev<<std::endl;

  of->Write();
  of->Close();

  gSystem->Exit(1);
  gROOT->GetApplication()->Run();
  //  theApp.Run();

  return 0;
}
