#include<fstream.h>
#include<iostream.h>

// class knucl

class Header;
class DetectorData;
class MCData;

class DetectorHit;
class Track;

// class k18ana
class ConfMan;
class CDSHitMan;
class CDSHit;
class SinDataMan;


#define CDH_ID 6
#define CDC_ID 5
#define DEBUG 0


void SimDatManG4()
{

  const int residnum=250; 

  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("./lib/libAll.so");
  gSystem->Load("~/work/geant/knucl3/libKnuclRootData.so");


  /*** assign input file & call tree ***/

  //  TFile *f = new TFile("~/work/geant/knucl3/tmp.root");
  TFile *f = new TFile("~/work/geant/knucl3/SinglePi/pi150.root");
  TTree *evtree = (TTree*)f->Get("tree");

  /*** conf file for new parameter ***/

  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer3057-3060.conf");
  conf->Initialize();

  Header* header = 0;
  DetectorData *detectorData = 0;
  MCData* mcData = 0;


  /*** declaration of classes ***/

  evtree->SetBranchAddress( "Header", &header );
  evtree->SetBranchAddress( "DetectorData", &detectorData );
  evtree->SetBranchAddress( "MCData", &mcData );


  /*** make output file & make tree ***/

  //  TFile *of = new TFile("./simout.root","recreate");
  TFile *of = new TFile("./simdata/sim_pi150.root","recreate");
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

  
  int nev = evtree->GetEntries();
  int fev = 0;
  //nev = 1000;
  cerr<<"# of events = "<<nev<<endl;
  for(int iev = 0; iev<nev; iev++ ){
    // cdsMan->Clear();
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
    }
    cerr<<"# of T0 hits   = "<<nT0<<endl;
    cerr<<"# of CDH hits  = "<<nCDH<<endl;
    cerr<<"# of BLC1 hits = "<<nBLC1<<endl;
    cerr<<"# of BLC2 hits = "<<nBLC2<<endl;
    cerr<<"# of CDC hits  = "<<nCDC<<endl;
    if(!nT0==1) continue;
    if(nCDC==0) continue;
    if(nCDH==0) continue;
    //if(nCDH<2) continue;


    //##########################//
    //### fill detector hits ###//
    //##########################//
    DetectorHit *detectorHit = 0;
    Double_t val0 = 0;
    for(int i=0; i<detectorData->detectorHitSize(); i++ ){
      detectorHit = detectorData->detectorHit(i);
      //--- T0 ---//
      if(detectorHit->detectorID()==CID_T0){
	Int_t seg = detectorHit->channelID()+1;
	simMan->SetT0Hit( seg, val0, val0, val0 );
      }
      //--- CDH ---//
      if(detectorHit->detectorID()==CID_CDH){
	Int_t seg = detectorHit->channelID()+1;
	simMan->SetCDHHit( seg, val0, val0, val0 );
      }
      //--- BLC ---//
      if(detectorHit->detectorID()==CID_BLC1 ||
	 detectorHit->detectorID()==CID_BLC2 ){
	Int_t layer = detectorHit->layerID()+1;
	Double_t wire = detectorHit->channelID()+1;
	int tag = (detectorHit->detectorID()==CID_BLC1) ? 1 : 0;
	simMan->SetBLCHit( tag, layer, wire, val0, val0 );
      }
      //--- CDC ---//
      if(detectorHit->detectorID()==CID_CDC){
	Int_t layer = detectorHit->layerID()+1;
	Double_t wire = detectorHit->channelID()+1;
	simMan->SetCDCHit( layer, wire, val0, val0 );
      }
    } // for(int i=0; i<detectorData->detectorHitSize(); i++ ){

    cdsMan = simMan->GetCDSHit();
    cdsMan->Calc(conf);



    blMan = simMan->GetBeamLineHit();
    blMan->Calc(conf);

    //-------------------//
    //--- refill data ---//
    //-------------------//
    int mT0 = 0;
    int mCDH = 0;
    int mBLC1[NumOfBLCLayers]; for(int i=0; i<NumOfBLCLayers; i++) mBLC1[i] = 0;
    int mBLC2[NumOfBLCLayers]; for(int i=0; i<NumOfBLCLayers; i++) mBLC2[i] = 0;
    int mCDC[NumOfCDCLayers]; for(int i=0; i<NumOfCDCLayers; i++) mCDC[i] = 0;
    for(int i=0; i<detectorData->detectorHitSize(); i++ ){
      detectorHit = detectorData->detectorHit(i);
      //--- T0 ---//
      if(detectorHit->detectorID()==CID_T0){
	blMan->T0(mT0)->SetHitID(i);
	blMan->T0(mT0)->SetCTMean(detectorHit->tdc()-time0);
	blMan->T0(mT0)->SetEMean(detectorHit->adc());
	blMan->T0(mT0)->SetHitPosition(detectorHit->pos().y());
	blMan->T0(mT0)->SetADCu(100);
	blMan->T0(mT0)->SetADCd(100);
	blMan->T0(mT0)->SetTDCu(100);
	blMan->T0(mT0)->SetTDCd(100);
	mT0++;
      }
      //--- CDH ---//
      if(detectorHit->detectorID()==CID_CDH){
	cdsMan->CDH(mCDH)->SetHitID(i);
	cdsMan->CDH(mCDH)->SetCTMean(detectorHit->tdc()+time0);
	cdsMan->CDH(mCDH)->SetEMean(detectorHit->adc());
	cdsMan->CDH(mCDH)->SetHitPosition(detectorHit->pos().z());
	cdsMan->CDH(mCDH)->SetADCu(100);
	cdsMan->CDH(mCDH)->SetADCd(100);
	cdsMan->CDH(mCDH)->SetTDCu(100);
	cdsMan->CDH(mCDH)->SetTDCd(100);
	mCDH++;
      }
      //--- BLC1 ---//
      if(detectorHit->detectorID()==CID_BLC1){
	Int_t layer = detectorHit->layerID()+1;
	Double_t dl = detectorHit->dx();
	dl = dl/10.; // cm
	dl = fabs(dl);
	blMan->BLC1(layer,mBLC1[layer-1])->SetHitID(i);
	blMan->BLC1(layer,mBLC1[layer-1])->SetDriftLength(dl);
	blMan->BLC1(layer,mBLC1[layer-1])->SetTDC(500);
	blMan->BLC1(layer,mBLC1[layer-1])->SetDriftTime(1);
	mBLC1[layer-1]++;
      }
      //--- BLC2 ---//
      if(detectorHit->detectorID()==CID_BLC2){
	Int_t layer = detectorHit->layerID()+1;
	Double_t dl = detectorHit->dx();
	dl = dl/10.; // cm
	dl = fabs(dl);
	blMan->BLC2(layer,mBLC2[layer-1])->SetHitID(i);
	blMan->BLC2(layer,mBLC2[layer-1])->SetDriftLength(dl);
	blMan->BLC2(layer,mBLC2[layer-1])->SetTDC(500);
	blMan->BLC2(layer,mBLC2[layer-1])->SetDriftTime(1);
	mBLC2[layer-1]++;
      }
      //--- CDC ---//
      if(detectorHit->detectorID()==CID_CDC){
	Int_t layer = detectorHit->layerID()+1;
	Double_t dl = detectorHit->dx();

	dl = dl/10.; // cm
	dl = fabs(dl);

	cdsMan->CDC(layer,mCDC[layer-1])->SetHitID(i);
	cdsMan->CDC(layer,mCDC[layer-1])->SetDriftLength(dl);
	cdsMan->CDC(layer,mCDC[layer-1])->SetTDC(500);
	cdsMan->CDC(layer,mCDC[layer-1])->SetDriftTime(1);
	mCDC[layer-1]++;

      }
    }


#if DEBUG
    for(int i=0; i<blMan->nT0(); i++){
      HodoscopeLikeHit *hod = blMan->T0(i);
      cerr<<"T0 : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<endl;
    }
    for(int i=0; i<cdsMan->nCDH(); i++){
      HodoscopeLikeHit *hod = cdsMan->CDH(i);
      cerr<<"CDH : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<endl;
    }
    for(int i=1; i<=NumOfBLCLayers; i++){
      for(int j=0; j<blMan->nBLC1(i); j++){
	ChamberLikeHit *chm = blMan->BLC1(i, j);
	cerr<<"BLC1 : "<<chm->hid()<<" "<<chm->layer()<<" "<<chm->wire()<<" "<<chm->dl()<<endl;
      }
    }
    for(int i=1; i<=NumOfBLCLayers; i++){
      for(int j=0; j<blMan->nBLC2(i); j++){
	ChamberLikeHit *chm = blMan->BLC2(i, j);
	cerr<<"BLC2 : "<<chm->hid()<<" "<<chm->layer()<<" "<<chm->wire()<<" "<<chm->dl()<<endl;

      }
    }

    for(int i=1; i<=NumOfCDCLayers; i++){
      for(int j=0; j<cdsMan->nCDC(i); j++){
	CDCHit *cdc = cdsMan->CDC(i, j);
	cerr<<"CDC : "<<cdc->hid()<<" "<<cdc->layer()<<" "<<cdc->wire()<<" "<<cdc->dl()<<endl;

      }
    }
#endif

    cdsTrackMan->Execute(cdsMan,conf);
    otree->Fill();
    fev++;
  
  } // for(int iev = 0; iev<nev; iev++ ){

  cerr<<"filled events = "<<fev<<endl;

  of->Write();
  of->Close();
  return;  
}
