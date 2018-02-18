//----------------------------------------------------------------//
// ===== UserSimDatG4.cpp =====
// CDC tracking program for MC data
//----------------------------------------------------------------//
//  exe-file: sim [./sim $(ConfFile) $(OutFile) $(InFile)]
//  input:  conf-file & G4(knucl4.10)-file
//  output: CDC-tracking-file
//           <- including class EventHeader and CDSTrackingMan
//  *NOTE*: the same conf-file used to generate G4-file is needed
//----------------------------------------------------------------//
//  updated by F.S, 2016 12/28
//----------------------------------------------------------------//
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
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h> 
#include <TRandom.h>

#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "SimDataMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
#include "MathTools.h"
#include "Particle.h"
#include "Tools.h"
#include "TrackTools.h"
#include "KnuclRootData.h" // copy from knucl4.10

using namespace std;

int main( int argc, char **argv )
{
  if( argc!=4 ){
    cout<<argv[0]<<" $(ConfFile) $(OutFile) $(InFile)"<<endl;
    return 0;
  }

  ConfMan *conf = new ConfMan(argv[1]);
  conf->Initialize();

  TFile *infile = new TFile(argv[3]);
  TTree *tree2 = (TTree*)infile->Get("tree2");
  RunHeaderMC *runHeader = new RunHeaderMC();
  tree2-> SetBranchAddress("RunHeaderMC", &runHeader);
  if( tree2->GetEntries()==1 ) cout<<"  !!! tree2 entries==1 !!!"<<endl;
  else tree2->GetEntry(0);

  TTree *tree = (TTree*)infile->Get("tree");
  EventHeaderMC *evHeaderMC = new EventHeaderMC();
  DetectorData  *detData  = new DetectorData();
  ReactionData  *reacData = new ReactionData();
  MCData        *mcData   = new MCData();
  tree->SetBranchAddress("EventHeaderMC", &evHeaderMC);
  tree->SetBranchAddress("DetectorData", &detData);
  tree->SetBranchAddress("ReactionData", &reacData);
  tree->SetBranchAddress("MCData", &mcData);

  TFile *outfile = new TFile(argv[2], "recreate");
  TTree *evTree  = new TTree("EventTree", "EventTree");
  EventHeader      *header = new EventHeader();
  BeamLineHitMan   *blMan  = new BeamLineHitMan();
  CDSHitMan        *cdsMan = new CDSHitMan();
  BeamLineTrackMan *bltrackMan  = new BeamLineTrackMan();
  CDSTrackingMan   *cdstrackMan = new CDSTrackingMan();
  evTree->Branch("EventHeader", &header);
  evTree->Branch("CDSTrackingMan", &cdstrackMan);

  SimDataMan *simMan = new SimDataMan();

  int eventn = tree->GetEntries();
  int stopn  = conf->GetStopEvNum();
  int exen   = ( 0<stopn && stopn<eventn ) ? stopn : eventn;
  
  cout<<"===== CDC tracking in MC START ====="<<endl;
  cout<<"     # of All  Event in G4-file:       "<<eventn<<endl;
  cout<<"     # of Stop Event in analyzer.conf: "<<stopn<<endl;
  cout<<"     # of Exe  Event in this program:  "<<exen<<endl;
  cout<<"========================================================="<<endl;

  for( int ev=0; ev<exen; ev++ ){
    if( ev%100==0 ) cout<<">     Event : "<<ev<<" finished"<<endl;
   
    tree->GetEntry(ev); // Get MC data

    header->SetEventNumber(ev);

    simMan->Convert(detData, conf, blMan, cdsMan); // convert KnuclRootData->BeamLine/CDSHitMan

    cdstrackMan->Execute(cdsMan, conf);

#if 0
    //-------------------- event reductoin --------------------//
    // NOTE: 
    //  at least one chaged track in CDC is required.
    //  ==>> inclusive analysis is ignore
    if( !cdstrackMan->nGoodTrack() ){
      header->Clear();
      blMan->Clear();
      cdsMan->Clear();
      bltrackMan->Clear();
      cdstrackMan->Clear();
    }
    //-------------------- event reductoin --------------------//
#endif

    evTree->Fill();

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
  } // for( int ev=0; ev<exen; ev++ ){

  cout<<"===== CDC tracking in MC END ====="<<endl;

  outfile->Write();
  outfile->Close();

  infile->Close();

  delete conf;

  delete runHeader;
  delete evHeaderMC;
  delete detData;
  delete reacData;
  delete mcData;

  delete header;
  delete blMan;
  delete cdsMan;
  delete bltrackMan;
  delete cdstrackMan;

  delete simMan;

  delete infile;
  delete outfile;

  return 0;
}



