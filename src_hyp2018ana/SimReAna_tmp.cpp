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
  ProtonArm::Initialize(confMan);

  std::cout<<"===== Geant Data Analysis Start ====="<<std::endl;
  std::cout<<"> Con fFile   : "<<conffile<<std::endl;
  std::cout<<"> MC File     : "<<mcfile<<std::endl;
  std::cout<<"> Data File   : "<<datafile<<std::endl;
  std::cout<<"> Output File : "<<outfile<<std::endl;

  BeamLineHitMan *blMan = new BeamLineHitMan();
  CDSHitMan *cdsMan = new CDSHitMan();
  CDSTrackingMan *cdstrackMan = new CDSTrackingMan();
  BeamLineTrackMan *bltrackMan = new BeamLineTrackMan();

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
  evTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);

  SimDataMan *simData = new SimDataMan();
  TFile *of=new TFile(outfile, "recreate");
  TTree *outTree=new TTree("EventTree", "EventTree");
  AnaInfo *anaInfo=new AnaInfo();
  outTree-> Branch("CDSTrackingMan", &cdstrackMan);
  outTree-> Branch("AnaInfo", &anaInfo);

  of-> cd();
  initHistFC();
  initHistCDS();

  if( mcTree->GetEntries()!=evTree->GetEntries() ){
    cout<<"!!!!! mcTree evTree Entry dont match !!!!!"<<endl;
    //    return 0;
  }

  for( int ev=0; ev<evTree->GetEntries(); ev++ ){
    if( ev%1000==0 ) cout<<"> Event Number : "<<ev<<endl;
    mcTree->GetEntry(ev);
    evTree-> GetEntry(ev);

    simData->Convert(detData, confMan, blMan, cdsMan);
    bltrackMan-> DoTracking(blMan, confMan, true, true);

    BeamInfo beam;
    beam.SetPID(Beam_Kaon);
    beam.SetBLMan(blMan);
    beam.SetBLDC(bltrackMan);
    std::vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);
    if( T0hits.size()==1 ) beam.SetT0(T0hits[0]);
    if( bltrackMan->ntrackBPC()==1 ) beam.SetBPCid(0);
    if( bltrackMan->ntrackBLC2()==1 ) beam.SetBLC2id(0);
    TVector3 init_mom, init_pos;
    for( int i=0; i<mcData->trackSize(); i++ ){
      if( mcData->track(i)->parentTrackID()==0 && mcData->track(i)->pdgID()==321 ){
        init_pos=0.1*mcData->track(i)->vertex();
        init_mom=0.001*mcData->track(i)->momentum();
      }
    }
    beam.SetVtxMom(init_mom.Mag());
    double D5mom=init_mom.Mag();
    while( true ){
      double beam_out, beam_tof;
      ELossTools::CalcElossBeamTGeo(beam.T0pos(), init_pos, D5mom, parMass[Beam_Kaon], beam_out, beam_tof);
      if( fabs(init_mom.Mag()-beam_out)<1.0e-8 ){
        break;
      }
      D5mom += init_mom.Mag()-beam_out;
    }
    beam.SetD5mom(D5mom);

    if( T0hits.size()==1 && bltrackMan->ntrackBPC()==1 && bltrackMan->ntrackBLC2()==1 ) beam.SetFlag(true);
    anaInfo-> AddBeam(beam);

    if( anaInfo->beam(0)->flag() ){
      for( int i=0; i<cdstrackMan->nTrack(); i++ ){
        CDSTrack *track=cdstrackMan->Track(i);
        CDSInfo cds=MyTools::makeCDSInfo(i, cdsMan, cdstrackMan, anaInfo->beam(0), bltrackMan, confMan);
	anaInfo-> AddCDS(cds);
      }
      if( anaInfo->minDCA() ){ anaInfo->beam(0)-> SetVertex(anaInfo->minDCA()->vertexBeam()); }

      for( int i=0; i<cdstrackMan->nTrack(); i++ ){
        for( int j=i+1; j<cdstrackMan->nTrack(); j++ ){
          CDS2Info cds2=MyTools::makeCDS2Info(cdstrackMan, i, j, anaInfo);
          anaInfo-> AddCDS2(cds2);
        }
      }

      for( int i=0; i<anaInfo->nCDS(); i++ ){
	if( anaInfo->CDS(i)->flag() ){
	  fillCDS(anaInfo->CDS(i), cdsMan, cdstrackMan);
	}
      }

      ForwardNeutralInfo fnInfo=MyTools::makeFN(blMan, anaInfo, NC_THRE);
      if( fnInfo.pid()==F_Gamma || fnInfo.pid()==F_Neutron ) anaInfo->AddNeutral(fnInfo);

      ForwardChargeInfo fcInfo=MyTools::makeFC(confMan, blMan, anaInfo);
      if( fcInfo.step()>1 ){
        anaInfo->AddCharge(fcInfo);
	fillFC(confMan, 0, blMan, cdsMan, cdstrackMan, anaInfo, anaInfo->forwardCharge(0));
      }
    }
    outTree->Fill();

    blMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
    cdsMan->Clear();
    anaInfo->Clear();

    if( ev>25000 ) break;
  }
  of-> Write();
  of-> Close();
  std::cout<<"===== Geant Data Analysis Finish ====="<<std::endl;
}
