
#include <iostream>
#include <string>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>

#include "ConfMan.h"
#include "KnuclRootData.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "SimDataMan.h"
#include "AnaInfo.h"

#include "MyHistFC.h"
#include "MyMCTools.h"
#include "MyHistCDS.h"

using namespace std;

int main(int argc, char** argv)
{
  TString conffile;
  TString mcfile;
  TString datafile;
  TString outfile;

  if( argc==3 ){
    std::cout<<"  !!! "<<argv[0]<<" ${input_dir} ${outputfile}"<<std::endl;
    conffile = Form("%s/analyzer.conf", argv[1]);
    mcfile = Form("%s/test_mc.root", argv[1]);
    datafile = Form("%s/sim.root", argv[1]);
    outfile = argv[2];
  }
  else if( argc==4 ){
    conffile = Form("%s/analyzer.conf", argv[1]);
    mcfile = Form("%s/test_mc.root", argv[1]);
    datafile = Form("%s/%s", argv[1], argv[2]);
    outfile = argv[3];
  }
  else if( argc==5 ){
    conffile = argv[1];
    mcfile   = argv[3];
    datafile = argv[4];
    outfile  = argv[2];
  }
  else{
    std::cout<<"input "<<argv[0]<<" ${input_dir} ${outputfile}"<<std::endl;
    std::cout<<"input "<<argv[0]<<" ${input_dir} ${readfile} ${outputfile}"<<std::endl;
    std::cout<<"input "<<argv[0]<<" ${conffile} ${input1} ${input2} ${outputfile}"<<std::endl;
    return 0;
  }

  ConfMan *confMan = new ConfMan((std::string)conffile);
  confMan-> Initialize();
  std::cout<<"===== Geant Data Analysis Start ====="<<std::endl;
  std::cout<<"> Con fFile   : "<<conffile<<std::endl;
  std::cout<<"> MC File     : "<<mcfile<<std::endl;
  std::cout<<"> Data File   : "<<datafile<<std::endl;
  std::cout<<"> Output File : "<<outfile<<std::endl;

  BeamLineHitMan *blMan = new BeamLineHitMan();
  CDSHitMan *cdsMan = new CDSHitMan();
  CDSTrackingMan *cdstrackMan = new CDSTrackingMan();
  BeamLineTrackMan *bltrackMan = new BeamLineTrackMan();
  AnaInfo *anaInfo=new AnaInfo();

  TFile *f_mc = new TFile(mcfile);
  TTree *mcTree = (TTree*)f_mc-> Get("tree");
  TTree *mcTree2 = (TTree*)f_mc-> Get("tree2");
  MCData *mcData = new MCData();
  ReactionData *reacData = new ReactionData();
  DetectorData *detData = new DetectorData();
  mcTree-> SetBranchAddress("MCData", &mcData);
  mcTree-> SetBranchAddress("ReactionData", &reacData);
  mcTree-> SetBranchAddress("DetectorData", &detData);

  TFile *f_data = new TFile(datafile);
  TTree *evTree = (TTree*)f_data->Get("EventTree");
  evTree-> SetBranchAddress("AnaInfo", &anaInfo);
  evTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);

  SimDataMan *simData = new SimDataMan();
  TFile *of=new TFile(outfile, "recreate");
  of-> cd();
  initHistFC();
  initHistCDS();
  new TH1F("trig_mass", "trigger mass",   200, 1.0, 2.0);
  new TH1F("eff_p2pim_hit",  "p 2#pi^{-} effective", 200, 1.0, 2.0);
  new TH1F("eff_mass",  "effective mass", 200, 1.0, 2.0);

  if( mcTree->GetEntries()!=evTree->GetEntries() ){
    cout<<"!!!!! mcTree evTree Entry dont match !!!!!"<<endl;
    //    return 0;
  }

  for( int ev=0; ev<evTree->GetEntries(); ev++ ){
    if( ev%1000==0 ) cout<<"> Event Number : "<<ev<<endl;
    mcTree->GetEntry(ev);
    evTree-> GetEntry(ev);

    of-> cd();
    simData->Convert(detData, confMan, blMan, cdsMan);
    if( anaInfo->nFCharge()==1 && anaInfo->forwardCharge(0)->step()>1 ){
      fillFC(confMan, 0, blMan, cdsMan, cdstrackMan, anaInfo, anaInfo->forwardCharge(0));
    }

    if( anaInfo->nBeam()==1 && anaInfo->beam(0)->flag() ){
      for( int i=0; i<anaInfo->nCDS(); i++ ){
	if( anaInfo->CDS(i)->flag() ){
	  fillCDS(anaInfo->CDS(i), cdsMan, cdstrackMan);
	}
      }
    }

    if( MyMCTools::effFC(detData, mcData, blMan) && anaInfo->nBeam()==1 && anaInfo->beam(0)->flag() ){
      //      if( kd_pLpim(anaInfo) ) cout<<"> d(K-, p)\"L pi-\" event"<<endl;
      double Ymass = MyMCTools::S1385mass(reacData);
      MyHistTools::fillTH("trig_mass", Ymass);
      {

	Track *ptrack=0;
	Track* Strack=0;
	for( int i=0; i<mcData->trackSize(); i++ ){
	  Track *tmp=mcData->track(i);
	  if( tmp->parentTrackID()==0 && tmp->pdgID()==2212 ) ptrack=tmp;
	  if( tmp->parentTrackID()==0 && tmp->pdgID()==3114 ) Strack=tmp;
	}
	Track *Strack2=0;
	Track *pim_track=0;
	for( int i=0; i<mcData->trackSize(); i++ ){
	  Track *tmp=mcData->track(i);
	  if( tmp->parentTrackID()==Strack->trackID() && tmp->pdgID()==3212 ) Strack2=tmp;
	  if( tmp->parentTrackID()==Strack->trackID() && tmp->pdgID()==-211 ) pim_track=tmp;
	}
	Track *Ltrack=0;
	for( int i=0; i<mcData->trackSize(); i++ ){
	  Track *tmp=mcData->track(i);
	  if( tmp->parentTrackID()==Strack2->trackID() && tmp->pdgID()==3122 ) Ltrack=tmp;
	}
	Track *son_pi=0;
	Track *son_N=0;
	for( int i=0; i<mcData->trackSize(); i++ ){
	  Track *tmp=mcData->track(i);
	  if( tmp->parentTrackID()==Ltrack->trackID() && tmp->pdgID()<1000 ) son_pi=tmp;
	  if( tmp->parentTrackID()==Ltrack->trackID() && tmp->pdgID()<1000 ) son_N=tmp;
	}

	if( son_pi && son_pi->pdgID()==-211 ){
	  DetectorHit *CDHhit1=0;
	  DetectorHit *CDHhit2=0;
	  int nCDC1=0;
	  int nCDC2=0;
	  for( int i=0; i<pim_track->detectorHitLinkSize(); i++ ){
	    DetectorHit *hit=detData->detectorHit(pim_track->detectorHitLink(i));
	    if( hit->detectorID()==CID_CDH ) CDHhit1=hit;
	    if( hit->detectorID()==CID_CDC ) nCDC1++;
	  }

	  for( int i=0; i<son_pi->detectorHitLinkSize(); i++ ){
	    DetectorHit *hit=detData->detectorHit(son_pi->detectorHitLink(i))\
	      ;
	    if( hit->detectorID()==CID_CDH ) CDHhit2=hit;
	    if( hit->detectorID()==CID_CDC ) nCDC2++;
	  }
	  if( CDHhit1 && CDHhit2 ){
	    // MyHistTools::fillTH("eff_p2pim_hit", Ymass);
	    // cout<<"CDH hit 1 : "<<CDHhit1->channelID()<<"  CDH hit 2 : "<<CDHhit2->channelID()<<endl;
	    // cout<<"    n CDC : "<<nCDC1<<"   n CDC : "<<nCDC2<<endl;
	    // CDSInfo *cds1=MyMCTools::getCDSInfo(CDHhit1, anaInfo, cdsMan);
	    // CDSInfo *cds2=MyMCTools::getCDSInfo(CDHhit1, anaInfo, cdsMan);
	    // cout<<"  info1   : "<<cds1<<"     info2 : "<<cds2<<endl;
	    // if( cds1 && cds2 ){
	    //   cout<<"         "<<cds1->pid()<<"     "<<cds2->pid()<<endl;
	    // }
	  }
	}
      }

      if( kd_pS0pim(anaInfo) ){
	MyHistTools::fillTH("eff_mass", Ymass);
	//	cout<<"> d(K-, p)\"S0 pi-\" event"<<endl;
      }
    }
    //    if( ev>0 && ev%1000==0 ) break;
    blMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
    cdsMan->Clear();
    anaInfo->Clear();
  }
  of-> cd();
  of-> Write();
  of-> Close();
  std::cout<<"===== Geant Data Analysis Finish ====="<<std::endl;
}
