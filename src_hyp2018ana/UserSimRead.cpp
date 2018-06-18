#include "KnuclRootData.h"
#include "ConfMan.h"
#include "SimDataMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

using namespace std;

int main( int argc, char** argv )
{
  if( argc!=5 ){
    cout<<argv[0]<<" $(ConfFile) $(OutFile) $(KnuclFile) $(TrackingFile)"<<endl;
    return 0;
  }

  ConfMan *conf=new ConfMan(argv[1]);
  conf-> Initialize();

  TFile *simfile = new TFile(argv[3]);
  TTree *tree2=(TTree*)simfile->Get("tree2");
  RunHeaderMC *runHeader=0;
  tree2-> SetBranchAddress("RunHeaderMC", &runHeader);
  if( tree2->GetEntries()!=1 ){
    cout<<"  !!! tree2 entries!=1 !!!"<<endl;
    return 0;
  }
  else tree2->GetEntry(0);

  TTree *tree=(TTree*)simfile->Get("tree");
  EventHeaderMC *evHeaderMC=0;
  DetectorData *detData=0;
  ReactionData *reacData=0;
  MCData *mcData=0;
  tree-> SetBranchAddress("EventHeaderMC", &evHeaderMC);
  tree-> SetBranchAddress("DetectorData", &detData);
  tree-> SetBranchAddress("ReactionData", &reacData);
  tree-> SetBranchAddress("MCData", &mcData);

  TFile *anafile=new TFile(argv[4]);
  CDSHitMan *cdsMan=new CDSHitMan();
  BeamLineHitMan *blMan=new BeamLineHitMan();
  CDSTrackingMan *cdstrackMan = new CDSTrackingMan();
  BeamLineTrackMan *bltrackMan = new BeamLineTrackMan();

  TTree *evTree = (TTree*)anafile->Get("EventTree");
  evTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);

  SimDataMan *simMan = new SimDataMan();

  TFile *outfile=new TFile(argv[2], "recreate");
  outfile-> cd();
  if( evTree->GetEntries()!=tree->GetEntries() ){
    std::cout<<"  !!! TTree entries don't match !!!"<<std::endl;
    return 0;
  }

  std::cout<<"> TTree Entries : "<<evTree-> GetEntries()<<std::endl;
  for( int ev=0; ev<evTree-> GetEntries(); ev++ ){
    if( ev%1==0 ) std::cout<<"> Event Number "<<ev<<std::endl;
    tree-> GetEntry(ev);
    simMan-> Convert(detData, conf, blMan, cdsMan);
    bltrackMan-> DoTracking(blMan, conf, true, true);

    cdsMan-> Clear();
    blMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
  }
  return 0;
}

