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
#include "KnuclRootData.h"
#include "AnaInfo.h"
#include "MyTools.h"
#include "MyMCTools.h"

using namespace std;

int main( int argc, char **argv )
{
  if( argc!=4 ){
    cout<<argv[0]<<" $(ConfFile) $(OutFile) $(InFile)"<<endl;
    return 0;
  }

  ConfMan *conf=new ConfMan(argv[1]);
  conf-> Initialize();
  ProtonArm::Initialize(conf);

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
  AnaInfo *anaInfo = new AnaInfo();
  TTree *evTree=new TTree("EventTree", "EventTree");
  evTree-> Branch("CDSTrackingMan", &cdstrackMan);
  evTree-> Branch("AnaInfo", &anaInfo);

  SimDataMan *simMan = new SimDataMan();
  cout<<"===== knucl->k18ana Convert START ====="<<endl;
  cout<<"      All Event : "<<tree->GetEntries()<<endl;
  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    if( ev%100==0 ) cout<<">     Event : "<<ev<<" finished"<<endl;
    cout<<">     Event : "<<ev<<" finished"<<endl;
    tree-> GetEntry(ev);
    simMan-> Convert(detData, conf, blMan, cdsMan);

    bltrackMan-> DoTracking(blMan, conf, true, true);
    cdstrackMan-> Execute(cdsMan, conf);

    BeamInfo beam;
    beam.SetPID(Beam_Kaon);
    beam.SetBLMan(blMan);
    beam.SetBLDC(bltrackMan);
    std::vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);

    if( T0hits.size()==1 && bltrackMan->ntrackBPC()==1 && bltrackMan->ntrackBLC2()==1 ){
      beam.SetT0(T0hits[0]); beam.SetBLC2id(0); beam.SetBPCid(0);

      DetectorHit *mcT0 = MyMCTools::getHit(detData, T0hits[0]);
      Track *beam_track=MyMCTools::initTrack(mcData, 321); // 321 is PDG ID of K+
      TVector3 in_pos=0.1*beam_track->vertex();
      TVector3 in_mom=0.001*beam_track->momentum();
      cout<<Form("init hit pos (%lf, %lf, %lf)", in_pos.X(), in_pos.Y(), in_pos.Z())<<endl;
      cout<<Form("init hit mom (%lf, %lf, %lf)", in_mom.X(), in_mom.Y(), in_mom.Z())<<endl;

      //      TVector3 pos=0.1*mcT0->pos();
      TVector3 pos=beam.T0pos();
      TVector3 mom=0.001*mcT0->momentum();
      cout<<Form("T0 hit pos (%lf, %lf, %lf)", pos.X(), pos.Y(), pos.Z())<<endl;
      cout<<Form("T0 hit mom (%lf, %lf, %lf)", mom.X(), mom.Y(), mom.Z())<<endl;
      cout<<"T0 hit time : "<<mcT0->time()<<endl;
      beam.SetFlag(true);
      double calc_tof, calc_mom;
      ELossTools::CalcElossBeamTGeo(in_pos, pos, in_mom.Mag(), kpMass, calc_mom, calc_tof);
      cout<<"mc tof : "<<mcT0->time()<<"  calc tof : "<<calc_tof<<"  dt : "<<mcT0->time()-calc_tof<<endl;
      cout<<"mc mom : "<<in_mom.Mag()<<"  calc mom : "<<calc_mom<<"  dp : "<<mom.Mag()-calc_mom<<endl;

      double D5mom=in_mom.Mag(); double calc_tof2; double out_mom;
      while( true ){

	ELossTools::CalcElossBeamTGeo(pos, in_pos, D5mom, kpMass, out_mom, calc_tof2);
	D5mom+=in_mom.Mag()-out_mom;
	if( fabs(in_mom.Mag()-out_mom)<10e-5 ) break;
      }
      ELossTools::CalcElossBeamTGeo(pos, in_pos, D5mom, kpMass, out_mom, calc_tof2);

      beam.SetT0time(T0hits[0]->ctmean()+calc_tof-calc_tof2); // Set tof difference
      beam.SetD5mom(D5mom); // Set D5 momentum
      cout<<" calc_tof : "<<calc_tof<<"  calc_tof2 : "<<calc_tof2<<endl;
      cout<<" T0 time (k18ana) : "<<T0hits[0]->ctmean()<<endl;
      cout<<" dt : "<<calc_tof-calc_tof2<<"  dp : "<<in_mom.Mag()-out_mom<<endl;
    }
    anaInfo-> AddBeam(beam);

    if( anaInfo->beam(0)->flag() ){
      for( int i=0; i<cdstrackMan->nTrack(); i++ ){
	CDSTrack *track=cdstrackMan->Track(i);
	CDSInfo cds=MyTools::makeCDSInfo(i, cdsMan, cdstrackMan, anaInfo->beam(0), bltrackMan, conf, true);
	anaInfo-> AddCDS(cds);
      }
      if( anaInfo->minDCA() ){ anaInfo->beam(0)-> SetVertex(anaInfo->minDCA()->vertexBeam()); }

      for( int i=0; i<cdstrackMan->nTrack(); i++ ){
	for( int j=i+1; j<cdstrackMan->nTrack(); j++ ){
	  CDS2Info cds2=MyTools::makeCDS2Info(cdstrackMan, i, j, anaInfo);
	  anaInfo-> AddCDS2(cds2);
	}
      }

      ForwardNeutralInfo fnInfo=MyTools::makeFN(blMan, anaInfo, NC_THRE);
      if( fnInfo.pid()==F_Gamma || fnInfo.pid()==F_Neutron ) anaInfo->AddNeutral(fnInfo);

      ForwardChargeInfo fcInfo=MyTools::makeFC(conf, blMan, anaInfo);
      if( fcInfo.step()>1 ){
	anaInfo->AddCharge(fcInfo);
      }
    }
    evTree->Fill();

    blMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
    cdsMan->Clear();
    anaInfo->Clear();
  }

  outfile->Write();
  outfile->Close();

  return 0;
}



