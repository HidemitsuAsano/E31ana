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

// These Linked by Makefile.user
#include "HistManwMC.h"

static const double D5MomRes=0.001;

int main(int argc, char** argv)
{
  //  gSystem->Load("libMinuit.so");
  //  gSystem->Load("libPhysics.so");

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
    mcfile   = argv[2];
    datafile = argv[3];
    outfile  = argv[4];
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
  std::cout<<"===== Geant Data Analysis Finish ====="<<std::endl;

  BeamLineHitMan *blMan = new BeamLineHitMan();
  CDSHitMan *cdsMan = new CDSHitMan();
  CDSTrackingMan *cdstrackMan = new CDSTrackingMan();
  BeamLineTrackMan *bltrackMan = new BeamLineTrackMan();
  SimDataReader *simReader = new SimDataReader();

  TFile *f_mc = new TFile(mcfile);
  TTree *mcTree = (TTree*)f_mc-> Get("tree");
  TTree *mcTree2 = (TTree*)f_mc-> Get("tree2");

  simReader-> setTree(mcTree);
  simReader-> setTree2(mcTree2);

  TFile *f_data = new TFile(datafile);
  TTree *evTree = (TTree*)f_data->Get("EventTree");
  evTree-> SetBranchAddress("BeamLineHitMan", &blMan);
  evTree-> SetBranchAddress("CDSHitMan", &cdsMan);
  evTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);

  TFile *rtFile =new TFile(outfile, "recreate");
  HistManwMC *histMan = new HistManwMC(rtFile, confMan);
  histMan-> setSimReader(simReader);
  histMan-> setK18ana(blMan, cdsMan, bltrackMan, cdstrackMan);

  if( mcTree2 ){
    if( mcTree2->GetEntries()!=1 ){
      std::cout<<"  !!! tree2 entries not 1 !!!"<<std::endl;
      return 0;
    }
    mcTree2-> GetEntry(0);
  }
  std::cout<<"===== Initialize Histgram for Check START ====="<<std::endl;
  rtFile-> cd();
  new TH1F("EventReduction", "Event Reduction", 20, -0.5, 19.5);
  new TH1F("nT0", "T0 multiplicity", 6, -0.5, 5.5);
  new TH1F("ntrackBLC2", "BLC2 ntrack",         5, -0.5, 4.5);
  new TH1F("BLC2_chi2",  "BLC2 chi2",        1000,  0.0, 100);
  new TH1F("BLC2_time",  "BLC2 track time",  1000,  -20, 80);

  new TH1F("ntrackBPC", "BPC ntrack",         5, -0.5, 4.5);
  new TH1F("BPC_chi2",  "BPC chi2",        1000,  0.0, 100);
  new TH1F("BPC_time",  "BPC track time",  1000,  -20, 80);

  new TH2F("BLC2BPC_diff", "BLC2 BPC pos diff", 1000, -10, 10, 1000, -10, 10);
  new TH2F("BLC2BPC_dir_diff", "BLC2 BPC dir diff", 1000, -0.1, 0.1, 1000, -0.1, 0.1);
  histMan-> initHist();
  std::cout<<"===== Initialize Histgram for Check FINISH ====="<<std::endl;

  if( evTree->GetEntries()!=mcTree->GetEntries() ){
    std::cout<<"  !!! TTree entries don't match !!!"<<std::endl;
    return 0;
  }

  std::cout<<"> TTree Entries : "<<evTree-> GetEntries()<<std::endl;
  for( int ev=0; ev<evTree-> GetEntries(); ev++ ){
    bltrackMan-> Clear();
    if( ev%5000==0 )
      std::cout<<"> Event Number "<<ev<<std::endl;
      
    TH1F *h1;
    TH2F *h2;

    evTree-> GetEntry(ev);
    mcTree-> GetEntry(ev);

    bltrackMan-> DoTracking(blMan, confMan, true, true);
    //    bltrackMan-> DoTracking2(blMan, confMan);
    TH1F *h1_ER = (TH1F*)rtFile-> Get("EventReduction");
    h1_ER-> Fill(0);
    
    int nT0=0;
    int T0seg=-1;
    double T0time=DBL_MIN;
    for( int i=0; i<blMan->nT0(); i++ ){
      if( blMan->T0(i)-> CheckRange() ){
	T0time = blMan->T0(i)-> ctmean();
	T0seg = blMan->T0(i)-> seg();
	nT0++;
      }
    }
    h1 = (TH1F*)rtFile-> Get("nT0"), h1->Fill(nT0);
    if( nT0!=1 ){
      //      std::cout<<"  !!! nT0="<<nT0<<" !!!"<<std::endl;
      continue;
    }
    h1_ER-> Fill(1);
      
    MCData *mcData = simReader-> getMCData();
    DetectorData *detData = simReader->getDetectorData();

    int nKbeam=0;
    TVector3 vtx, mom;

    bool T0_flag=false;
    std::string hit_par="no hit";
    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit=detData->detectorHit(i);
      if( hit->detectorID()==CID_T0 && hit->channelID()+1==T0seg ){
	if( hit->pdg()==321 ) T0_flag=true;
	hit_par = simReader->parName(hit->pdg());
      }
    }
    if( !T0_flag ){
      //      std::cout<<" ev : "<<ev<<"  T0 don't hit beam  particle : "<<hit_par<<std::endl;
      continue;
    }

    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData-> track(i);
      if( track->pdgID()==321 && track->parentTrackID()==0 ){
	nKbeam++;
	vtx = 0.1*track->vertex();     // mm->cm
	mom = 0.001*track->momentum(); // MeV->GeV
      }
    }
    if( nKbeam!=1 ){
      std::cout<<"  !!! K beam find "<<nKbeam<<" !!!"<<std::endl;
      continue;
    }
    h1_ER-> Fill(2);
      
    TVector3 T0pos = vtx+(-110.5-vtx.Z())*mom.Unit()*(1./mom.CosTheta());
    double beam_mom_mc=mom.Mag();
    double D5mom = beam_mom_mc;
    double beam_out, beam_tof;
    while( true ){
      ELossTools::CalcElossBeamTGeo(T0pos, vtx, D5mom, parMass[Beam_Kaon], beam_out, beam_tof);
      if( fabs(beam_mom_mc-beam_out)<0.00001 ) break;
      D5mom += beam_mom_mc-beam_out;
    }
    //    std::cout<<"  D5 mom : "<<D5mom<<"  true mom : "<<beam_mom_mc<<"  calc mom : "<<beam_out<<std::endl;

    //*** Add resolution Beam Momentum
    D5mom *= gRandom-> Gaus(1.0, D5MomRes);
      
    h1 = (TH1F*)rtFile-> Get("ntrackBLC2"), h1-> Fill(bltrackMan->ntrackBLC2());
    for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
      LocalTrack *track = bltrackMan->trackBLC2(i);
      double tracktime = track-> GetTrackTime();
      double chi2 = track-> chi2all();
      h1 = (TH1F*)rtFile-> Get("BLC2_chi2"), h1-> Fill(chi2);
      h1 = (TH1F*)rtFile-> Get("BLC2_time"), h1-> Fill(tracktime);
    }
      
    h1 = (TH1F*)rtFile-> Get("ntrackBPC"), h1-> Fill(bltrackMan->ntrackBPC());
    for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
      LocalTrack *track = bltrackMan->trackBPC(i);
      double tracktime = track-> GetTrackTime();
      double chi2 = track-> chi2all();
      h1 = (TH1F*)rtFile-> Get("BPC_chi2"), h1-> Fill(chi2);
      h1 = (TH1F*)rtFile-> Get("BPC_time"), h1-> Fill(tracktime);
    }
      
    if( bltrackMan->ntrackBPC()==1 && bltrackMan->ntrackBLC2()==1 ){
      LocalTrack *trackBLC2 = bltrackMan-> trackBLC2(0);
      LocalTrack *trackBPC  = bltrackMan-> trackBPC(0);
      double center_z = 0.5*(-130-20.3);
      TVector3 BLC2_track_pos = trackBLC2-> GetPosatZ(center_z);
      TVector3 BPC_track_pos  = trackBPC-> GetPosatZ(center_z);
      TVector3 BLC2_mom_dir = trackBLC2-> GetMomDir();
      TVector3 BPC_mom_dir = trackBPC-> GetMomDir();
      h2 = (TH2F*)rtFile-> Get("BLC2BPC_diff"), h2-> Fill((BLC2_track_pos-BPC_track_pos).X(), (BLC2_track_pos-BPC_track_pos).Y());
      h2 = (TH2F*)rtFile-> Get("BLC2BPC_dir_diff"), h2-> Fill((BLC2_mom_dir-BPC_mom_dir).X(), (BLC2_mom_dir-BPC_mom_dir).Y());
    }

    if( bltrackMan-> ntrackBLC2()!=1 || bltrackMan-> ntrackBPC()!=1 ){
      //      std::cout<<"  !!! BPC Track : "<<bltrackMan->ntrackBPC()<<std::endl;
      continue;
    }

    LocalTrack *trackBLC2 = bltrackMan-> trackBLC2(0);
    double center_z = 0.5*(-130-20.3);
    TVector3 BLC2_track_pos = trackBLC2-> GetPosatZ(center_z);
    TVector3 BLC2_mom_dir = trackBLC2-> GetMomDir();
      
    DetectorHit *bp_hit = 0;
    int nBPD=0;
    int BPD_p_hit = 0;
    int BPD_pim_hit = 0;
    int BPD_pip_hit = 0;
    int BPD_beam_hit = 0;
    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit = detData->detectorHit(i);
      if( hit->detectorID()==CID_BPD && hit->adc()>0.1 ){
	nBPD++;
	if( hit->pdg()==2212 ) BPD_p_hit++, bp_hit=hit;
	else if( hit->pdg()==-211 ) BPD_pim_hit++;
	else if( hit->pdg()== 211 ) BPD_pip_hit++;
	else if( hit->pdg()== 321 && hit->parentID()==0 ) BPD_beam_hit++;
	//	else std::cout<<" BPD hit  seg : "<<hit->channelID()+1<<" par : "<<simReader->parName(hit->pdg())<<" dE : "<<hit->adc()<<std::endl;

	  
      }
    }

    //    std::cout<<"> n BPC track : "<<bltrackMan->ntrackBPC()<<std::endl; 
    int nBPC_beam[8];
    int BPC_beam_wire[8];
    int nBPC_p[8];
    int BPC_p_wire[8];
    for( int i=0; i<8; i++ ) nBPC_beam[i]=0, nBPC_p[i]=0;
    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit = detData->detectorHit(i);
      if( hit->detectorID()==CID_BPC ){
	if( hit->pdg()==321 && hit->parentID()==0 ){
	  nBPC_beam[hit->layerID()]++;
	  BPC_beam_wire[hit->layerID()]=hit->channelID()+1;
	}
	if( hit->pdg()==2212 ){
	  nBPC_p[hit->layerID()]++;
	  BPC_p_wire[hit->layerID()]=hit->channelID()+1;
	}
      }
    }
    for( int lay=1; lay<=8; lay++ ){
      //      std::cout<<"> BPC lay"<<lay<<" nHit : "<<blMan->nBPC(lay)<<"  nHit beam : "<<nBPC_beam[lay-1]<<"  nHit p : "<<nBPC_p[lay-1]<<std::endl;
    }
      
#if 0
    if( BPD_p_hit>0 ){
      Track *parent_track = simReader->getMCTrack(bp_hit->parentID());
      if( parent_track ){
	std::cout<<" PID : "<<simReader->parName(parent_track->pdgID())<<"  Reaction : "<<simReader->getReactionContainer()->getReaction(parent_track->trackID())<<std::endl;
	if( parent_track->pdgID()==3122 ){
	  std::cout<<"> n BPC track : "<<bltrackMan->ntrackBPC()<<std::endl; 
	  int nBPC_beam[8];
	  int BPC_beam_wire[8];
	  int nBPC_p[8];
	  int BPC_p_wire[8];
	  DetectorData *detData = simReader->getDetectorData();
	  for( int i=0; i<8; i++ ) nBPC_beam[i]=0, nBPC_p[i]=0;
	  for( int i=0; i<detData->detectorHitSize(); i++ ){
	    DetectorHit *hit = detData->detectorHit(i);
	    if( hit->detectorID()==CID_BPC ){
	      if( hit->pdg()==321 && hit->parentID()==0 ){
		nBPC_beam[hit->layerID()]++;
		BPC_beam_wire[hit->layerID()]=hit->channelID()+1;
	      }
	      if( hit->pdg()==2212 ){
		nBPC_p[hit->layerID()]++;
		BPC_p_wire[hit->layerID()]=hit->channelID()+1;
	      }
	    }
	  }
	  
	  for( int lay=1; lay<=8; lay++ ){
	    std::cout<<"> BPC lay"<<lay<<" nHit : "<<blMan->nBPC(lay)<<"  nHit beam : "<<nBPC_beam[lay-1]<<"  nHit p : "<<nBPC_p[lay-1]<<std::endl;
	  }
	}
      }
    }
#endif	
      
#if 0
    std::string input;
    std::cin>>input;
    if( input=="q" ){
      histMan-> Clear();
      break;
    }
#endif
 
#if 0
    if( bltrackMan-> ntrackBPC()!=1 ){
      //      std::cout<<"  !!! BPC Track : "<<bltrackMan->ntrackBPC()<<std::endl;
      continue;
    }
    h1_ER-> Fill(6);
#endif
    histMan-> setBeamPID(Beam_Kaon);
    histMan-> setBeamKaon(true);
    histMan-> setT0time(T0time);
    histMan-> setD5mom(D5mom);
    histMan-> setTrackBPC(bltrackMan->trackBPC(0));

    histMan-> ana(0);
    histMan-> fill(0);
   
    histMan->Clear();
  }

  TH1F *h1_ER = (TH1F*)rtFile-> Get("EventReduction");
  RunHeaderMC *runHeader = simReader->getRunHeader();
  h1_ER-> SetBinContent(1, runHeader->numGenerated());

  std::cout<<"> Event Analysis Finish"<<std::endl;

  delete confMan;
  delete blMan;
  delete cdsMan;
  delete cdstrackMan;
  delete bltrackMan;
  delete simReader;

  delete f_mc;
  delete f_data;

  histMan-> finit();

  std::cout<<"> Fill Closed"<<std::endl;
  return 0;
}
