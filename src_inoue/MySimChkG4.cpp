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
    cout<<argv[0]<<" $(ConfFile) $(OutFile(dummy)) $(InFile)"<<endl;
    return 0;
  }

  ConfMan *conf=new ConfMan(argv[1]);
  conf-> Initialize();
  ProtonArm::Initialize(conf);

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

  BeamLineHitMan *blMan=new BeamLineHitMan();
  CDSHitMan *cdsMan=new CDSHitMan();
  BeamLineTrackMan *bltrackMan=new BeamLineTrackMan();
  CDSTrackingMan *cdstrackMan=new CDSTrackingMan();
  AnaInfo *anaInfo = new AnaInfo();

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
    beam.SetPID(Beam_Kaon);
    beam.SetBLMan(blMan);
    beam.SetBLDC(bltrackMan);
    std::vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);

    if( T0hits.size()==1 && bltrackMan->ntrackBPC()==1 && bltrackMan->ntrackBLC2()==1 ){
      beam.SetT0(T0hits[0]); beam.SetBLC2id(0); beam.SetBPCid(0); beam.SetPID(Beam_Kaon);

      Track *beam_track=MyMCTools::initTrack(mcData, 321); // 321 is PDG ID of K+
      TVector3 in_pos=0.1*beam_track->vertex();
      TVector3 in_mom=0.001*beam_track->momentum();

      TVector3 pos=beam.T0pos();
      beam.SetFlag(true);
      double calc_tof, calc_mom;
      ELossTools::CalcElossBeamTGeo(in_pos, pos, in_mom.Mag(), kpMass, calc_mom, calc_tof);

      double D5mom=in_mom.Mag(); double calc_tof2; double out_mom;
      while( true ){
	ELossTools::CalcElossBeamTGeo(pos, in_pos, D5mom, kpMass, out_mom, calc_tof2);
	D5mom+=in_mom.Mag()-out_mom;
	if( fabs(in_mom.Mag()-out_mom)<10e-5 ) break;
      }
      ELossTools::CalcElossBeamTGeo(pos, in_pos, D5mom, kpMass, out_mom, calc_tof2);

      //      cout<<T0hits[0]->ctmean()+calc_tof-calc_tof2<<endl;
      beam.SetT0time(T0hits[0]->ctmean()+calc_tof-calc_tof2); // Set tof difference
      //      cout<<D5mom<<endl;
      beam.SetD5mom(gRandom->Gaus(D5mom, 0.003*D5mom)); // Set D5 momentum w/ mom resolution w/ 0.3%

      if( MyMCTools::goodBeam(mcData, detData) ){
	beam.SetFlag(true);
      }
    }
    anaInfo-> AddBeam(beam);

    if( anaInfo->beam(0)->flag() ){
      for( int i=0; i<cdstrackMan->nTrack(); i++ ){
	CDSInfo cds=MyTools::makeCDSInfo(i, cdsMan, cdstrackMan, anaInfo->beam(0), bltrackMan, conf, true);
	anaInfo-> AddCDS(cds);
      }
      if( anaInfo->minDCA() ){ 
	anaInfo->beam(0)-> SetVertex(anaInfo->minDCA()->vertexBeam()); 
      }

      for( int i=0; i<cdstrackMan->nTrack(); i++ ){
	for( int j=i+1; j<cdstrackMan->nTrack(); j++ ){
	  CDS2Info cds2=MyTools::makeCDS2Info(cdstrackMan, i, j, anaInfo);
	  anaInfo-> AddCDS2(cds2);
	}
      }

      ForwardNeutralInfo fnInfo=MyTools::makeFN(blMan, anaInfo, NC_THRE, true);
      if( fnInfo.pid()==F_Gamma || fnInfo.pid()==F_Neutron ) anaInfo->AddNeutral(fnInfo);

      ForwardChargeInfo fcInfo=MyTools::makeFC(conf, blMan, anaInfo);
      if( fcInfo.step()>1 ){
	anaInfo->AddCharge(fcInfo);
      }
    }
    //    anaInfo-> dump();
    // if( anaInfo-> nFNeutral()==1 ){
    //   anaInfo->forwardNeutral(0)->dump();
    //   DetectorHit *dhit=MyMCTools::getHit(detData, anaInfo->forwardNeutral(0)->NC(blMan));
    //   cout<<" hit pdgID="<<dhit->pdg()<<"  parentID="<<dhit->parentID()<<endl;
    //   cout<<"     time="<<dhit->time()<<"  dE="<<dhit->adc()<<endl;
    // }

    blMan->Clear();
    bltrackMan->Clear();
    cdstrackMan->Clear();
    cdsMan->Clear();
    anaInfo->Clear();
  }

  return 0;
}



