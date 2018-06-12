#include "ConfMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"
#include "SimDataMan.h"
#include "EventHeader.h"
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
  ProtonArm::Initialize(conf, 1.0);

  TFile *infile = new TFile(argv[3]);
  TTree *tree2=(TTree*)infile->Get("tree2");
  RunHeaderMC *runHeader=0;
  tree2-> SetBranchAddress("RunHeaderMC", &runHeader);
  if( tree2->GetEntries()==1 ) cout<<"  !!! tree2 entries!=1 !!!"<<endl;
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
  // TTree *evTree=new TTree("EventTree", "EventTree");
  // evTree-> Branch("CDSTrackingMan", &cdstrackMan);
  // evTree-> Branch("AnaInfo", &anaInfo);

  SimDataMan *simMan = new SimDataMan();
  cout<<"===== knucl->k18ana Convert START ====="<<endl;
  cout<<"      All Event : "<<tree->GetEntries()<<endl;
  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    if( ev%100==0 ) cout<<">     Event : "<<ev<<" finished"<<endl;

    tree-> GetEntry(ev);
    simMan-> Convert(detData, conf, blMan, cdsMan);

    bltrackMan-> DoTracking(blMan, conf, true, true);
    cdstrackMan-> Execute(cdsMan, conf);

    BeamInfo beam;
    beam.SetPID(Beam_Proton);
    beam.SetBLMan(blMan);
    beam.SetBLDC(bltrackMan);
    std::vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);

    if( T0hits.size()==1 && bltrackMan->ntrackBPC()==1 && bltrackMan->ntrackBLC2()==1 ){
      beam.SetT0(T0hits[0]); beam.SetBLC2id(0); beam.SetBPCid(0); beam.SetPID(Beam_Proton);
      {
	TVector3 pos=bltrackMan->trackBPC(0)->GetPosatZ(-20);
	std::cout<<Form("BPC pos (%lf, %lf, %lf)", pos.X(), pos.Y(), pos.Z())<<std::endl;
      }


      Track *beam_track=MyMCTools::initTrack(mcData, 2212); // 321 is PDG ID of K+
      TVector3 in_pos=0.1*beam_track->vertex();
      TVector3 in_mom=0.001*beam_track->momentum();
      DetectorHit *T0hit=MyMCTools::getHit(detData, T0hits[0]);

      TVector3 pos=beam.T0pos();
      beam.SetFlag(true);
      double calc_tof, calc_mom;
      double D5mom=gRandom->Gaus(in_mom.Mag(), 0.003*in_mom.Mag());
      beam.SetD5mom(D5mom);

      //      cout<<T0hits[0]->ctmean()+calc_tof-calc_tof2<<endl;
      //      cout<<D5mom<<endl;

      if( MyMCTools::goodP_beam(mcData, detData) ){
	beam.SetFlag(true);
	std::cout<<"Good Beam"<<std::endl;

	std::vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
	std::vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
	std::vector<HodoscopeLikeHit*> PChits=MyTools::getHodo(blMan, CID_PC);

	// std::cout<<"nFDC : "<<beam.nFDC1()<<std::endl;
	// std::cout<<"nBVC : "<<BVChits.size()<<std::endl;
	// std::cout<<"nCVC : "<<CVChits.size()<<std::endl;
	// std::cout<<"nPC : "<<PChits.size()<<std::endl;
	if( beam.nFDC1()==1 && BVChits.size()==1 ){
	  HodoscopeLikeHit *FChit=0;
	  if( CVChits.size()==0 && PChits.size()==1 ) FChit=PChits[0];
	  if( CVChits.size()==1 && PChits.size()==0 ) FChit=CVChits[0];
	  if( FChit ){
	    DetectorHit* FChit_mc=MyMCTools::getHit(detData, FChit);
	    //	    std::cout<<"Good BT"<<std::endl;
	    ForwardChargeInfo info;
	    info.SetPID(F_Proton);
	    info.SetVertex(beam.T0pos());
	    info.SetFDC1(beam.FDC1(0));
	    info.SetHodo(FChit);

	    //	    if( info.calc_beam_through_wUSWK(beam, conf).size()>0 ){
	    info.calc_simple_p_through(beam, conf, blMan);
	    //	    info.calc_beam_through_wUSWK(beam, conf, blMan);
	      std::cout<<"fl by data : "<<info.fl()<<std::endl;
	      TVector3 init_pos=beam_track->vertex();
	      TVector3 T0pos_mc=T0hit->pos();
	      double off_mc=0.1*(init_pos-T0pos_mc).Mag();
	      std::cout<<"  offset mc : "<<off_mc<<std::endl;
	      std::cout<<"  fl by mc : "<<0.1*beam_track->FlightLength()-off_mc<<std::endl;

	      std::cout<<"  TOF by data : "<<info.calc_tof()<<std::endl;
	      std::cout<<"  TOF by MC   : "<<FChit_mc->time()-T0hit->time()<<std::endl;

	      std::cout<<"  mom by TOF : "<<info.momByTOF()<<std::endl;
	      std::cout<<"  mom by MC  : "<<beam_track->momentum().Mag()<<std::endl;

	    // }
	    // else{
	    //   std::cout<<"!!! ForwardChargeInfo::calc_beam_through !!!"<<std::endl;
	    // }
	  }
	}
      }
    }


    blMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
    cdsMan->Clear();
    anaInfo->Clear();
  }

  outfile->Write();
  outfile->Close();

  std::cout<<"conv sim finish"<<std::endl;

  return 0;
}



