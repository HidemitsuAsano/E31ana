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
#include "CDHCluster.cpp"
#include "MyAnalysisBL.h"
#include "MyAnalysisCDS.h"
#include "MyAnalysisFWDNeutral.h"
#include "MyAnalysisdKnpipi.h"
#include "MyAnalysisdKnppim.h"
#include "MyAnalysisCheck.h"
#include "MyAnalysisHeKpipp.h"
#include "MyAnalysisCheckChi2.h"
#include "MyAnalysisHeKFWD.h"
#include "MyAnalysisHeKpippOneDummyLpair.h"
#include "MyAnalysisHeKpippMCReaction.h"
#include "MyAnalysisHeKpippMCKinFit.h"
#include "MyAnalysisHeKpippCheckResl.h"

//#define dKnpipi
//#define dKnppim
#define HeKpipp
//#define Check
//#define HeKFWD
//#define CheckChi2
//#define HeKpippOneDummyLpair
//#define HeKpippMCReaction
//#define HeKpippMCKinFit
//#define HeKpippCheckResl

#include "KnuclRootData.h"
//#include "/home/yamaga/geant/knucl4/include/KnuclRootData.h"
//#include "/home/yamaga/geant/knucl4/include/KnuclRootData.h"

int main( int argc, char **argv )
{
  //gSystem->Load(".lib/KnuclRootData_cc.so");
  gSystem->Load("/home/yamaga/geant/knucl4/build/KnuclRootData_cc.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");

  if(argc != 5 ){
    std::cout << "Plese set Conffile Outputfile Inputfile(MC data) Inputfile(SIM data)  "<< std::endl;
    return 0;
  }

  std::string confFile,outFile,inFile,simFile;
  if( argv[1][0] != '!' ){
    confFile = argv[1];
    if( confFile.find_last_of(".conf")==std::string::npos){
      std::cout<<"ConfFile should be *.conf"<<std::endl;
      return 0;
    }
  }
  if( argv[2][0] != '!' ){
    outFile = argv[2];
    if( outFile.find_last_of(".root")==std::string::npos ){
      std::cout<<"Oututfile should be *.root"<<std::endl;
      return 0;
    }
  }
  if( argv[3][0] != '!' ){
    inFile = argv[3];
    if( inFile.find_last_of(".root")==std::string::npos ){
      std::cout<<"Inputfile should be *.root"<<std::endl;
      return 0;
    }
  }
  if( argv[4][0] != '!' ){
    simFile = argv[4];
    if( simFile.find_last_of(".root")==std::string::npos ){
      std::cout<<"Inputfile should be *.root"<<std::endl;
      return 0;
    }
  }

  TFile *of = new TFile( outFile.c_str(), "recreate" ); 
  TFile *f = new TFile( inFile.c_str() );
  TFile *simf = new TFile( simFile.c_str() );

  ConfMan *conf = new ConfMan(confFile);
  conf->Initialize();
  conf->SetOutFileName(outFile);
  EventHeader* header = new EventHeader();
  EventHeaderMC* mcheader = new EventHeaderMC();
  SimDataMan *simMan = new SimDataMan(f, conf);
  CDSHitMan *cdsMan = new CDSHitMan();
  BeamLineHitMan *blMan = new BeamLineHitMan();
  BeamLineTrackMan *bltrackMan = new BeamLineTrackMan();
  CDSTrackingMan *cdstrackMan = new CDSTrackingMan();
  Particle* particle = new Particle();

  TTree* simtree = (TTree*)simf->Get("EventTree");
  simtree->SetBranchAddress("CDSHitMan", &cdsMan );
  simtree->SetBranchAddress("BeamLineHitMan", &blMan );
  simtree->SetBranchAddress("BeamLineTrackMan", &bltrackMan );
  simtree->SetBranchAddress("CDSTrackingMan", &cdstrackMan );
  simtree->SetBranchAddress("EventHeader",&header);
  simtree->SetBranchAddress("EventHeaderMC",&mcheader);
  simtree->SetBranchAddress("Particle",&particle);

#ifdef dKnpipi
  MyAnalysisdKnpipi* dKnpipiAna;
  dKnpipiAna  = new MyAnalysisdKnpipi(of, conf);
#endif
#ifdef dKnppim
  MyAnalysisdKnppim* dKnppimAna;
  dKnppimAna  = new MyAnalysisdKnppim(of, conf);
#endif
#ifdef HeKpipp
  MyAnalysisHeKpipp* HeKpippAna;
  HeKpippAna  = new MyAnalysisHeKpipp(of, conf);
#endif
#ifdef Check
  MyAnalysisCheck* CheckAna;
  CheckAna  = new MyAnalysisCheck(of, conf);
#endif
#ifdef HeKFWD
  MyAnalysisHeKFWD* HeKFWDAna;
  HeKFWDAna  = new MyAnalysisHeKFWD(of, conf);
#endif
#ifdef CheckChi2
  MyAnalysisCheckChi2* CheckChi2Ana;
  CheckChi2Ana  = new MyAnalysisCheckChi2(of, conf);
#endif
#ifdef HeKpippOneDummyLpair
  MyAnalysisHeKpippOneDummyLpair* HeKpippOneDummyLpairAna;
  HeKpippOneDummyLpairAna  = new MyAnalysisHeKpippOneDummyLpair(of, conf);
#endif
#ifdef HeKpippMCReaction
  MyAnalysisHeKpippMCReaction* HeKpippMCReactionAna;
  HeKpippMCReactionAna  = new MyAnalysisHeKpippMCReaction(of, conf);
#endif
#ifdef HeKpippMCKinFit
  MyAnalysisHeKpippMCKinFit* HeKpippMCKinFitAna;
  HeKpippMCKinFitAna  = new MyAnalysisHeKpippMCKinFit(of, conf);
#endif
#ifdef HeKpippCheckResl
  MyAnalysisHeKpippCheckResl* HeKpippCheckReslAna;
  HeKpippCheckReslAna  = new MyAnalysisHeKpippCheckResl(of, conf);
#endif

  int evn = simMan-> nEvent();
  if( conf->GetStopEvNum()>0 ){
    if( evn>conf->GetStopEvNum() ) evn=conf->GetStopEvNum();
  }

  std::cout<<" Read "<<simMan-> nEvent()<<" Event read"<<std::endl;
  bool flag = true;
  int SIM_Event_Number = 0;
  simtree -> GetEntry(SIM_Event_Number);
  for( int iev=0; iev<evn; iev++ ){
    simMan->Get(iev);
    if(iev%10000==0){
      std::cout<<" Event Number : "<<iev<<std::endl;
      std::cout<<" EventID : "<<simMan->GetEventHeader()->eventID()<<std::endl;
    }
    ReactionData* reaction = simMan->GetReactionData();
    DetectorData* detector = simMan->GetDetectorData();
    MCData* mcdata = simMan->GetMCData();
    EventHeaderMC* mchead = simMan->GetEventHeader();

    while( mcheader->eventID()<mchead->eventID() ){
      SIM_Event_Number++;
      if(SIM_Event_Number>simtree->GetEntries()){
        flag = false;
        break;
      }
      simtree -> GetEntry(SIM_Event_Number);
    }
    if(mcheader->eventID()>mchead->eventID()){
      continue;
    }
    if(!flag){ continue; }
    if(iev%10000==0){
      std::cout<<" SIMEventID : "<<mcheader->eventID()<<std::endl;
    }


    //int pronum = 3;
    //for(int itra=0; itra<mcdata->trackSize(); itra++){
    //  Track* tra = mcdata->track(itra);
    //  if(tra->pdgID()==13122){ pronum = tra->trackID(); }
    //  if(tra->parentTrackID()==0 || tra->parentTrackID()==pronum){
    //    simMan->PrintTrack(mcdata->track(itra));
    //  }
    //}

#ifdef dKnpipi
    dKnpipiAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle);
#endif
#ifdef dKnppim
    dKnppimAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle);
#endif
#ifdef HeKpipp
    HeKpippAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle);
#endif
#ifdef Check
    CheckAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle);
#endif
#ifdef HeKFWD
    HeKFWDAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle);
#endif
#ifdef CheckChi2
    CheckChi2Ana->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle);
#endif
#ifdef HeKpippOneDummyLpair
    HeKpippOneDummyLpairAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle);
#endif
#ifdef HeKpippMCReaction
    HeKpippMCReactionAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle, mcdata, reaction, detector);
#endif
#ifdef HeKpippMCKinFit
    HeKpippMCKinFitAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle, mcdata, reaction, detector);
#endif
#ifdef HeKpippCheckResl
    HeKpippCheckReslAna->DoAnalysis(conf, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle, mcdata, reaction, detector);
#endif

    header->Clear();
    mcheader->Clear();
    cdsMan-> Clear();
    blMan-> Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
    particle->Clear();
    simMan-> Clear();
#ifdef dKnpipi
    dKnpipiAna->Clear();
#endif
#ifdef dKnppim
    dKnppimAna->Clear();
#endif
#ifdef HeKpipp
    HeKpippAna->Clear();
#endif
#ifdef Check
    CheckAna->Clear();
#endif
#ifdef HeKFWD
    HeKFWDAna->Clear();
#endif
#ifdef CheckChi2
    CheckChi2Ana->Clear();
#endif
#ifdef HeKpippOneDummyLpair
    HeKpippOneDummyLpairAna->Clear();
#endif
#ifdef HeKpippMCReaction
    HeKpippMCReactionAna->Clear();
#endif
#ifdef HeKpippMCKinFit
    HeKpippMCKinFitAna->Clear();
#endif
#ifdef HeKpippCheckResl
    HeKpippCheckReslAna->Clear();
#endif
  }

#ifdef dKnpipi
  delete dKnpipiAna;
#endif
#ifdef dKnppim
  delete dKnppimAna;
#endif
#ifdef HeKpipp
  delete HeKpippAna;
#endif
#ifdef Check
  delete CheckAna;
#endif
#ifdef HeKFWD
  delete HeKFWDAna;
#endif
#ifdef CheckChi2
  delete CheckChi2Ana;
#endif
#ifdef HeKpippOneDummyLpair
  delete HeKpippOneDummyLpairAna;
#endif
#ifdef HeKpippMCReaction
  delete HeKpippMCReactionAna;
#endif
#ifdef HeKpippMCKinFit
  delete HeKpippMCKinFitAna;
#endif
#ifdef HeKpippCheckResl
  delete HeKpippCheckReslAna;
#endif

  std::cout << "Finish!!!" << std::endl;

  return 0;
}
