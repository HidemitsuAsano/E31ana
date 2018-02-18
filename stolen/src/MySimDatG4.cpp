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
#include "MyAnalysisBL.h"
#include "MyAnalysisCDS.h"
#include "MyAnalysisFWDNeutral.h"
#include "MyAnalysisFWDCharged.h"

#include "KnuclRootData.h"
//#include "/home/yamaga/geant/knucl4/include/KnuclRootData.h"

int main( int argc, char **argv )
{
  //gSystem->Load(".lib/KnuclRootData_cc.so");
  gSystem->Load("/home/yamaga/geant/knucl4/build/KnuclRootData_cc.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");

  if(argc != 4 ){
    std::cout << "Plese set Conffile Outputfile Inputfile  "<< std::endl;
    return 0;
  }

  std::string confFile,outFile,inFile;
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

  TFile *f = new TFile( inFile.c_str() );
  ConfMan *conf = new ConfMan(confFile);
  conf->Initialize();
  EventHeader* header = new EventHeader();
  EventHeaderMC* mcheader = new EventHeaderMC();
  RunHeaderMC* runheader = new RunHeaderMC();
  SimDataMan *simMan = new SimDataMan(f, conf);
  CDSHitMan *cdsMan = new CDSHitMan();
  BeamLineHitMan *blMan = new BeamLineHitMan();
  BeamLineTrackMan *bltrackMan = new BeamLineTrackMan();
  CDSTrackingMan *cdstrackMan = new CDSTrackingMan();
  Particle* particle = new Particle();

  TFile *of = new TFile( outFile.c_str(), "recreate" ); 
  TTree *otree = new TTree("EventTree","EventTree");
  otree->Branch( "CDSHitMan", &cdsMan );
  otree->Branch( "BeamLineHitMan", &blMan );
  otree->Branch( "BeamLineTrackMan", &bltrackMan );
  otree->Branch( "CDSTrackingMan", &cdstrackMan );
  otree -> Branch("EventHeader",&header);
  otree -> Branch("EventHeaderMC",&mcheader);
  otree -> Branch("RunHeaderMC",&runheader);
  otree -> Branch("Particle",&particle);

	//int seed = runheader->seed();
	//std::cout << "seed : " << seed << std::endl;
	//return 0;

  outFile.insert(outFile.find(".root"),"_anaAll");
  TFile *of2 = new TFile( outFile.c_str(), "recreate" ); 
  MyAnalysisBL* blAna = new MyAnalysisBL(of2, conf, true);
  MyAnalysisCDS* cdsAna  = new MyAnalysisCDS(of2, conf, true);
  MyAnalysisFWDNeutral* fwdnAna = new MyAnalysisFWDNeutral(of2, conf,true);
  MyAnalysisFWDCharged* fwdpAna = new MyAnalysisFWDCharged(of2, conf,true);

  int evn = simMan-> nEvent();
  if( conf->GetStopEvNum()>0 ){
    if( evn>conf->GetStopEvNum() ) evn=conf->GetStopEvNum();
  }

  std::cout<<" Read "<<simMan-> nEvent()<<" Event read"<<std::endl;
  for( int iev=0; iev<evn; iev++ ){
    ReactionData* reaction = simMan->GetReactionData();
    MCData* mcdata = simMan->GetMCData();
    simMan-> Get(iev, cdsMan, blMan);
    mcheader = simMan->GetEventHeader();
    cdstrackMan->Execute(cdsMan,conf);

		of2->cd();
		bool flag_bl=false;
		flag_bl = blAna->DoAnalysisSim(conf,header,blMan,bltrackMan,particle,simMan->beam_mom());
		bool flag_cds=false;
		flag_cds = cdsAna->DoAnalysisSim(conf,header,blMan,bltrackMan,cdsMan,cdstrackMan,particle);
		fwdnAna->DoAnalysis(conf,header,blMan,bltrackMan,cdsMan,cdstrackMan,particle);
		fwdpAna->DoAnalysis(conf,header,blMan,bltrackMan,cdsMan,cdstrackMan,particle);

		//int pronum = 3;
		//for(int itra=0; itra<mcdata->trackSize(); itra++){
		//  Track* tra = mcdata->track(itra);
		//  if(tra->pdgID()==13122){ pronum = tra->trackID(); }
		//  if(tra->parentTrackID()==0 || tra->parentTrackID()==pronum){
		//    simMan->PrintTrack(mcdata->track(itra));
		//  }
		//}
		if(flag_bl){
			if(flag_cds){
				of->cd();
				otree->Fill();
			}
			else{
				cdstrackMan->Clear();
				blMan-> Clear();
				cdsMan-> Clear();
				pBeam* tmp_beam = particle->beam(0);
				particle->Clear(); 
				particle->AddBeam(*tmp_beam);
				of->cd();
				otree->Fill();
			}
		}

		if(iev%1000==0){
			std::cout<<" Event Number : "<<iev<<std::endl;
			std::cout<<" EventID : "<<simMan->GetEventHeader()->eventID()<<std::endl;
			//simMan-> PrintEvent();
			//simMan-> PrintHit();
			//std::cout<<" Please input any char"<<std::endl;
			//std::cout<<" q: quit"<<std::endl;
			//char in;
			//std::cin>>in;
			//if( in=='q' ) break;
		}

		header->Clear();
		mcheader->Clear();
		cdsMan-> Clear();
		blMan-> Clear();
		bltrackMan->Clear();
		cdstrackMan->Clear();
		blAna->Clear();
		cdsAna->Clear();
		fwdnAna->Clear();
		fwdpAna->Clear();
		particle->Clear();
		simMan-> Clear();
	}

	of->Write();
	of->Close();

	of2->Write();
	of2->Close();

	std::cout << "Finish!!!" << std::endl;

	return 0;
}
