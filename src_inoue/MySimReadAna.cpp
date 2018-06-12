#include "KnuclRootData.h"
#include "ConfMan.h"
#include "SimDataMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"
#include "AnaInfo.h"

#include "MyTools.h"
#include "AnaInfo.h"
#include "MyHistReadCDS.h"
#include "MyHistReadNC.h"
#include "MyHistReadFC.h"
#include "MyHistReadKNpim.h"
#include "MyHistReadKNpip.h"
#include "MyHistReadKNkm.h"
#include "MyHistReadKNpipi_D2.h"
#include "MyHistReadKPpimpim_D2.h"
#include "MyHistReadCDS_Lpim_D2.h"

#include "MyHistMCAcc.h"
#include "MyHistMCChkCou.h"
#include "MyHistMCFC.h"
#include "MyHistMC_pLpim.h"

using namespace std;

int main( int argc, char** argv )
{
  if( argc!=5 ){
    cout<<argv[0]<<" $(ConfFile) $(OutFile) $(KnuclFile) $(AnaFile)"<<endl;
    return 0;
  }

  ConfMan *conf=new ConfMan(argv[1]);
  conf-> Initialize();
  ProtonArm::Initialize(conf, 1.0);

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

  TTree *npipi_tree=new TTree("npipi_event", "npipi_event");
  npipi_tree-> Branch("AnaInfo", &anaInfo);

  initHistReadCDS();
  initHistReadNC();
  initHistReadFC();
  initHistReadKNpim();
  initHistReadKNpip();
  initHistReadKNkm();

  initHistMCAcc();
  initHistMCChkCou();
  initHistMCFC();
  initHistMC_pLpim();

  if( MyAnaTools::targetNA() ){
    initHistReadKNpipi_D2(0, anaInfo);
    initHistReadKPpimpim_D2();
    initHistReadCDS_Lpim_D2();
  }

  std::cout<<"> TTree Entries : "<<evTree-> GetEntries()<<std::endl;
  int t0=time(0);

  for( int ev=0; ev<evTree-> GetEntries(); ev++ ){
    if( ev%1000==0 ){
      int t1=time(0);
      std::cout<<"> Event Number "<<ev<<"  Time : "<<t1-t0<<std::endl;
    }
    tree-> GetEntry(ev);
    evTree-> GetEntry(ev);
    simMan-> Convert(detData, conf, blMan, cdsMan);
    bltrackMan-> DoTracking(blMan, conf, true, true);

    anaInfo-> ClearFC();
    {
      ForwardChargeInfo fcInfo=MyTools::makeFC(conf, blMan, anaInfo);
      if( fcInfo.step()>=0 ){
	anaInfo-> AddCharge(fcInfo);
      }
    }
    // for( int i=0; i<anaInfo->nFCharge(); i++ ){
    //   if( anaInfo->forwardCharge(i)->pid()==F_Proton ){
    //     anaInfo->forwardCharge(i)->fit_forward(blMan, anaInfo->beam(0), anaInfo->minDCA(), conf);
    //   }
    // }

    fillHistReadCDS(anaInfo, cdsMan, cdstrackMan);
    fillHistReadNC(0, blMan, anaInfo);
    fillHistReadFC(0, conf, blMan, anaInfo);
    fillHistReadKNpim(0, blMan, anaInfo);
    fillHistReadKNpip(0, blMan, anaInfo);
    fillHistReadKNkm(0, blMan, anaInfo);

    fillHistMCAcc_fn(detData, mcData, reacData, anaInfo, blMan, cdstrackMan);
    fillHistMCAcc_fp(detData, mcData, reacData, anaInfo, blMan, cdstrackMan);
    fillHistMCChkCou(evHeaderMC, detData, mcData);
    fillHistMCFC(conf, anaInfo, cdsMan, blMan, cdstrackMan, bltrackMan, reacData, mcData, detData);
    fillHistMC_pLpim(detData, mcData, reacData, anaInfo, blMan, cdstrackMan);

    if( MyAnaTools::targetNA() ){
      fillHistReadKNpipi_D2(0, blMan, cdsMan, cdstrackMan, anaInfo);
      fillHistReadKPpimpim_D2(0, conf, blMan, cdsMan, cdstrackMan, anaInfo);
      fillHistReadCDS_Lpim_D2(0, blMan, anaInfo);
    }

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

