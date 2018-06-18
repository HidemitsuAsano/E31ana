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

  ConfMan *conf=new ConfMan(argv[1]);
  conf-> Initialize();

  TFile *infile = new TFile(argv[3]);
  TTree *tree2=(TTree*)infile->Get("tree2");
  RunHeaderMC *runHeader=0;
  tree2-> SetBranchAddress("RunHeaderMC", &runHeader);
  if( tree2->GetEntries()==1 ) cout<<"  !!! tree2 entries==0 !!!"<<endl;
  else tree2->GetEntry(0);

  TTree *tree=(TTree*)infile->Get("tree");
  EventHeaderMC *evHeaderMC=0;
  DetectorData *detData=0;
  ReactionData *reacData=0;
  MCData *mcData=0;
  tree-> SetBranchAddress("EventHeaderMC", &evHeaderMC);
  tree-> SetBranchAddress("DetectorData", &detData);
  tree-> SetBranchAddress("ReactionData", &reacData);
  tree-> SetBranchAddress("MCData", &mcData);

  TFile *outfile=new TFile(argv[2], "recreate");
  BeamLineHitMan *blMan=new BeamLineHitMan();
  CDSHitMan *cdsMan=new CDSHitMan();
  BeamLineTrackMan *bltrackMan=new BeamLineTrackMan();
  CDSTrackingMan *cdstrackMan=new CDSTrackingMan();
  TTree *evTree=new TTree("EventTree", "EventTree");
  evTree-> Branch("CDSTrackingMan", &cdstrackMan);

  SimDataMan *simMan = new SimDataMan();
  cout<<"===== knucl->k18ana Convert START ====="<<endl;
  cout<<"      All Event : "<<tree->GetEntries()<<endl;
  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    if( ev%100==0 ) cout<<">     Event : "<<ev<<" finished"<<endl;
    cout<<">     Event : "<<ev<<" finished"<<endl;
    tree-> GetEntry(ev); // Get MC data
    simMan-> Convert(detData, conf, blMan, cdsMan); // convert KnuclRootData->BeamLine/CDSHitMan

    bltrackMan-> DoTracking(blMan, conf, true, true);
    cdstrackMan-> Execute(cdsMan, conf);

    evTree->Fill();

    blMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
    cdsMan->Clear();
  }

  outfile->Write();
  outfile->Close();

  return 0;
}



