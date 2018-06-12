#include "KnuclRootData.h"
#include "ConfMan.h"
#include "SimDataMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"
#include "AnaInfo.h"

#include "MyTools.h"
#include "AnaInfo.h"
#include "MyHistReadCalibCDC.h"

using namespace std;

int main( int argc, char** argv )
{
  if( argc!=5 ){
    cout<<argv[0]<<" $(ConfFile) $(OutFile) $(KnuclFile) $(AnaFile)"<<endl;
    return 0;
  }

  ConfMan *conf=new ConfMan(argv[1]);
  conf-> Initialize();

  SimDataMan *simMan = new SimDataMan();
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

  AnaInfo *anaInfo = new AnaInfo();
  TTree *evTree = (TTree*)anafile->Get("EventTree");
  evTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);
  evTree-> SetBranchAddress("AnaInfo", &anaInfo);

  if( evTree->GetEntries()!=tree->GetEntries() ){
    std::cout<<"  !!! TTree entries don't match !!!"<<std::endl;
    return 0;
  }

  TFile *outfile=new TFile(argv[2], "recreate");
  outfile-> cd();

  initHistReadCalibCDC(conf);

  std::cout<<"> TTree Entries : "<<evTree-> GetEntries()<<std::endl;
  for( int ev=0; ev<evTree-> GetEntries(); ev++ ){
    if( ev%1000==0 ) std::cout<<"> Event Number "<<ev<<std::endl;
    tree-> GetEntry(ev);
    evTree-> GetEntry(ev);
    simMan-> Convert(detData, conf, blMan, cdsMan);
    bltrackMan-> DoTracking(blMan, conf, true, true);

    fillHistReadCalibCDC(conf, blMan, cdsMan, cdstrackMan, anaInfo);

    cdsMan-> Clear();
    blMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
    anaInfo-> Clear();
  }
  cout<<"finish"<<endl;

  outfile-> Write();
  outfile-> Close();

  return 0;
}

