//H. Asano
//This code is originated from: EventAnalysis_pipipnn_sakuma.cpp
//----------------------------------------------------------------//
//----------------------------------------------------------------//
//  input : raw data, conf-file, & CDC-tracking-file
//  output: when $(OutFile) is "tmp.root", the following 3 files are generated.
//     "tmp.root":      histogram file
//     "tmp_CDC.root": condensed event file from the CDC-tracking-file with p/p/pi selection
//           <- including class EventHeader and CDSTrackingMan
//     "tmp_npippim.root": basic information of pi+pi-n event is listed up in TTree

#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include <TLorentzVector.h>
#include "TrackTools.h"
#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "Tools.h"
#include "ELossTools.h"

#include <TDatabasePDG.h>
#include <Math/ProbFuncMathCore.h>

#define VTX_SIGMA 1 // 1: Sigma reconstruction,  vertex = K- & p
                    // 0: Lambda reconstruction, vertex = K- & pi+

#define DEBUG 0

#define KFDEBUG 0 // verbose level of the KinFitter
// 0: quiet, 1: print result, 2: print iterations, 3: print also matrices

//-- set run# --//

const int MaxTreeSize = 1000000000;

bool DoCDCRetiming = false;
int Verbosity = 0;

namespace blcuts{
  
  //Kaon selection from TOF T0-BHD
  const double beam_tof_k_min=27.8588;
  const double beam_tof_k_max=29.5663;
  const double beam_tof_pi_min=25.0;//rough 
  const double beam_tof_pi_max=27.0;//rough
  
  //FDC1 cuts
  const double fdc1_time_window_min=-30;
  const double fdc1_time_window_max=100;
  const double fdc1_time_min=-10;
  const double fdc1_time_max= 10;
  const double fdc1_chi2_max= 10;
  
  //BLC1 cuts
  const double blc1_time_window_min=-30;
  const double blc1_time_window_max=100;
  const double blc1_time_min=-10;
  const double blc1_time_max= 10;
  const double blc1_chi2_max= 10;
  
  //BLC2 cuts
  const double blc2_time_window_min=-30;
  const double blc2_time_window_max=100;
  const double blc2_time_min=-10;
  const double blc2_time_max= 10;
  const double blc2_chi2_max= 10;
  
  //BPC cuts
  const double bpc_time_window_min=-30;
  const double bpc_time_window_max=100;
  const double bpc_time_min=-10;
  const double bpc_time_max= 10;
  const double bpc_chi2_max= 10;
  
  //D5 cuts
  const double d5_chi2_max=30;
  
  //BLC2-BPC matching cuts
  const double blc2bpc_x_min=-1.1015;//Run78  rough
  const double blc2bpc_x_max=1.21206;//Run78  rough
  const double blc2bpc_y_min=-1.1015;//Run78  rough
  const double blc2bpc_y_max=1.21206;//Run78  rough 
  const double blc2bpc_dx_min=-0.0253846;//Run78 rough
  const double blc2bpc_dx_max=0.0242834;//Run78 rough
  const double blc2bpc_dy_min=-0.0246937;//Run78 rough
  const double blc2bpc_dy_max=0.02502;//Run78 rough
}

namespace cdscuts{
  const int cdhmulti = 3;
  const double tdc_cdh_max = 25; // ns
  const double cds_chi2_max=30;

}


class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;// histograms
  TFile *rtFile2;// condensed file
  TFile *rtFile3;// pi+,pi-,n event
  TFile *cdcFile;
  TTree *cdcTree;
  TTree *evTree;
  TTree *npippimTree;

  const EventHeader *header_CDC; // original in CDC-tracking-file
  CDSTrackingMan *trackMan_CDC; // original in CDC-tracking-file

  EventHeader *header;
  ScalerMan *scaMan;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  CDSTrackingMan *trackMan;
  BeamLineTrackMan *bltrackMan;
  int t0,t1;

  int ncheckev;
  int checkrun[26];
  int checkev[26];

  double scainit[40];
  double scaend[40];
  bool INIT;
  bool SCAOVERFLOW[40];


  int AllGoodTrack;
  int nTrack;
  int CDC_Event_Number;
  std::ofstream ofs;
  int PIDBeam;
  int blc1GoodTrackID;
  int blc2GoodTrackID;
  int bpcGoodTrackID;

  //** counters for filing **//
  int nFill_pippim;
  int nFill_npippim;
  //** counters for event abort **//
  int nAbort_nGoodTrack;
  int nAbort_nCDH;
  int nAbort_nT0;
  int nAbort_pid_beam;
  int nAbort_nbpc;
  int nAbort_bpctrack;
  int nAbort_singleBLtrack;
  int nAbort_fblc2bpc;
  int nAbort_flagbmom;
  int nAbort_ftarget;
  int nAbort_CDHiso;
  int nAbort_ppipi;
  int nAbort_end;

  //= = = = npippim final-sample tree = = = =//
  TLorentzVector mom_beam;   // 4-momentum(beam)
  TLorentzVector mom_target; // 4-momentum(target)
  TLorentzVector mom_pip;    // 4-momentum(pi+)
  TLorentzVector mom_pim;    // 4-momentum(pi-)
  TLorentzVector mom_p;      // 4-momentum(proton)
  TLorentzVector mom_n;      // 4-momentum(neutron)
  double NeutralBetaCDH; // veracity of neutral particle on CDH
  double dE;   // energy deposit on CDH
  TVector3 vtx_reaction; // vertex(reaction)
  int run_num;   // run number
  int event_num; // event number
  int block_num; // block number
  TLorentzVector kf1mom_beam;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector kf1mom_pip;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector kf1mom_pim;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector kf1mom_p;      // 4-momentum(proton) after kinematical refit for pi- Sigma+
  TLorentzVector kf1mom_n;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  double kf1_chi2;   // chi2 of kinematical refit
  double kf1_NDF;    // NDF of kinematical refit
  double kf1_status; // status of kinematical refit -> details can be found in this code
  double kf1_pvalue; // p-value of kinematical refit
  TLorentzVector kf2mom_beam;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
  TLorentzVector kf2mom_pip;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
  TLorentzVector kf2mom_pim;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
  TLorentzVector kf2mom_p;      // 4-momentum(proton) after kinematical refit for pi+ Sigma-
  TLorentzVector kf2mom_n;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
  double kf2_chi2;   // chi2 of kinematical refit
  double kf2_NDF;    // NDF of kinematical refit
  double kf2_status; // status of kinematical refit -> details can be found in this code
  double kf2_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc
  //= = = = npippim final-sample tree = = = =//

  
public:
  void Initialize( ConfMan *conf );
  void InitializeHistogram();
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  void Clear( int &nAbort);
};


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
EventAnalysis::EventAnalysis()
  : EventTemp()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
EventAnalysis::~EventAnalysis()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
}


void EventAnalysis::Initialize( ConfMan *conf )
{
  std::cout << " *** Enter EventAnalysis::Initialize " << std::endl;
  std::cout << " Verbosity Level " << Verbosity << std::endl;
  std::cout << " CDC Retiming ? " ;
  
  if(DoCDCRetiming) std::cout << " Yes" << endl;
  else              std::cout << "  No" << endl;
  std::cout << " CDH multiplicity cuts: " << cdscuts::cdhmulti << std::endl;
  

  INIT = true;
  //  spillinit = spillfini = -1;
  for( int i=0; i<40; i++ ){
    SCAOVERFLOW[i] = false;
    scaend[i] = 0;
  }

  //** CDS file open **//
  confMan = conf;
  cdcFile = new TFile( confMan->GetCDSTrackFileName().c_str() );
  if( !cdcFile->IsOpen() ){
    std::cerr<<" !!! failed to open " <<confMan->GetCDSTrackFileName().c_str()<< "  !!!"<<std::endl;
    exit( false );
  }
  
  gSystem->Load( "libPhysics.so" );

  //** Getting CDSTracking info. from CDCfile **//
  cdcTree = (TTree*)cdcFile->Get( "EventTree" );
  header_CDC = 0;
  trackMan_CDC = 0;
  cdcTree->SetBranchAddress( "CDSTrackingMan", &trackMan_CDC );
  cdcTree->SetBranchAddress( "EventHeader" ,&header_CDC );  

  //** output file 1 : histograms **//
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  rtFile->cd();
  InitializeHistogram();

  //** output file 2 : condensed event with pi+/pi-/n selection **//
  std::string outfile2 = confMan->GetOutFileName();
  outfile2.insert( outfile2.size()-5, "_CDC" );
  std::cout<<"CDC file "<<outfile2<<std::endl;
  rtFile2 = new TFile( outfile2.c_str(), "recreate" );
  rtFile2->cd();
  evTree = new TTree( "EventTree", "EventTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  evTree->Branch( "CDSTrackingMan", &trackMan_CDC ); //** = fill original CDSTrackingMan (w/o dE correction) **//

  //** output file 3 : npippim final-sample tree **//
  std::string outfile3 = confMan->GetOutFileName();
  outfile3.insert( outfile3.size()-5, "_npippim" );
  std::cout<<"npippim file "<<outfile3<<std::endl;
  rtFile3 = new TFile( outfile3.c_str(), "recreate" );
  rtFile3->cd();
  npippimTree = new TTree( "EventTree", "EventTree" );
  npippimTree->Branch( "mom_beam",   &mom_beam );
  npippimTree->Branch( "mom_target", &mom_target );
  npippimTree->Branch( "mom_pip", &mom_pip );
  npippimTree->Branch( "mom_pim", &mom_pim );
  npippimTree->Branch( "mom_p", &mom_p );
  npippimTree->Branch( "mom_n", &mom_n );
  npippimTree->Branch( "NeutralBetaCDH", &NeutralBetaCDH );
  npippimTree->Branch( "dE", &dE );
  npippimTree->Branch( "vtx_reaction", &vtx_reaction );
  npippimTree->Branch( "run_num", &run_num );
  npippimTree->Branch( "event_num", &event_num );
  //npippimTree->Branch( "block_num", &block_num );
  //npippimTree->Branch( "kf1mom_beam",   &kf1mom_beam );
  //npippimTree->Branch( "kf1mom_pip", &kf1mom_pip );
  //npippimTree->Branch( "kf1mom_pim", &kf1mom_pim );
  //npippimTree->Branch( "kf1mom_p", &kf1mom_p );
  //npippimTree->Branch( "kf1mom_n", &kf1mom_n );
  //npippimTree->Branch( "kf1_chi2", &kf1_chi2 );
  //npippimTree->Branch( "kf1_NDF", &kf1_NDF );
  //npippimTree->Branch( "kf1_status", &kf1_status );
  //npippimTree->Branch( "kf1_pvalue", &kf1_pvalue );
  //npippimTree->Branch( "kf2mom_beam",   &kf2mom_beam );
  //npippimTree->Branch( "kf2mom_pip", &kf2mom_pip );
  //npippimTree->Branch( "kf2mom_pim", &kf2mom_pim );
  //npippimTree->Branch( "kf2mom_p", &kf2mom_p );
  //npippimTree->Branch( "kf2mom_n", &kf2mom_n );
  //npippimTree->Branch( "kf2_chi2", &kf2_chi2 );
  //npippimTree->Branch( "kf2_NDF", &kf2_NDF );
  //npippimTree->Branch( "kf2_status", &kf2_status );
  //npippimTree->Branch( "kf2_pvalue", &kf2_pvalue );
  //npippimTree->Branch( "kf_flag", &kf_flag );

  trackMan = new CDSTrackingMan(); //** = dE correction is performed in this code **//
  if( trackMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }


  rtFile->cd();

  t0 = clock();
  ResetCounters();
}


void EventAnalysis::ResetCounters()
{
  AllGoodTrack = 0;
  nTrack = 0;
  CDC_Event_Number = 0;

  nFill_pippim = 0;
  nFill_npippim  = 0;

  nAbort_nGoodTrack = 0;
  nAbort_nCDH = 0;
  nAbort_nT0 = 0;
  nAbort_pid_beam = 0;
  nAbort_nbpc = 0;
  nAbort_bpctrack = 0;
  nAbort_singleBLtrack = 0;
  nAbort_fblc2bpc = 0;
  nAbort_flagbmom = 0;
  nAbort_ftarget = 0;
  nAbort_end = 0;
  nAbort_ppipi = 0;
  nAbort_CDHiso = 0;


  return;
}



//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
void EventAnalysis::USca( int nsca, unsigned int *sca )
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
#if 0
  std::cout << " *** Enter EventAnalysis::USca " << std::endl;
#endif

  rtFile->cd();
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

#if 0
  for( int i=0; i<11; i++ )
    std::cout<<sca[i]<<"\t";
  std::cout<<std::endl;
#endif

  if( INIT && sca[0]<5 ){
    for( int i=0; i<scaMan->nsca(); i++ ){
      scainit[i] = scaMan->sca(i)->val();
      std::cout<<i<<" "<<scainit[i]<<std::endl;
    }
    INIT = false;
  }

  if( INIT ){
    scaMan->Clear();
    return;
  }

  for( int i=0; i<scaMan->nsca(); i++ ){
    TString name = scaMan->sca(i)->name();
    if( scaend[i]>9.9e+07 && scaend[i]>scaMan->sca(i)->val() ) SCAOVERFLOW[i] = true;
    scaend[i] = scaMan->sca(i)->val();
    if( SCAOVERFLOW[i] ) scaend[i] += 1.0e+08;
  }
  
  t1 = clock();
  TH1F* h1 = (TH1F*)gFile->Get( "Time" );
  h1->Fill( Block_Event_Number, (t1-t0)/CLOCKS_PER_SEC );
  scaMan->Clear();
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
bool EventAnalysis::UAna( TKOHitCollection *tko )
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  // - - - - - - - - - - - - - - - - - - - - - - - - - - //
  // Event_Number: event number in "raw data"
  //  = header_CDC->ev() in cdc-tracking data
  // - - - - - - - - - - - - - - - - - - - - - - - - - - //
  // Block_Event_Number: block number in "raw data"
  //  = header_CDC->blev() in cdc-tracking data
  // - - - - - - - - - - - - - - - - - - - - - - - - - - //
  // CDC_Event_Number: event number tagged in this code
  // - - - - - - - - - - - - - - - - - - - - - - - - - - //

#if 0
  std::cout << " *** Enter EventAnalysis::UAna " << std::endl;
#endif

  //** fill event count **//
  Tools::Fill1D( Form("EventCheck"), 1 );
  Event_Number++;

  //** control of start and stop events **//
  {
    int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) {
      if( CDC_Event_Number>=cdcTree->GetEntries() ) return false;
      cdcTree->GetEntry( CDC_Event_Number );

      if( header_CDC->ev()<Event_Number ) return true;
      else if( header_CDC->ev()==Event_Number ){
	CDC_Event_Number++;    	  
	return true;
      }
    }
    if( status==2 ) return false; 
  }
  
#if DEBUG
  if( Event_Number>10000) return true; //** for debug **//
#endif
  
  if( Event_Number%1000==1 ){
    t1=clock();
    std::cout << "Run " << confMan->GetRunNumber()
	      << " Event# : " << Event_Number 
	      << " CDC_Event# : " << CDC_Event_Number
	      << " BlockEvent# : " << Block_Event_Number
	      << " AllTrack# : " << nTrack
	      << " GoodTrack# : " << AllGoodTrack
	      << " Time (s): " << (t1-t0)/CLOCKS_PER_SEC << std::endl;
  }


  //** Checking event number consistentcy with CDC tracking file **//
  if( CDC_Event_Number>=cdcTree->GetEntries() ) return false;
  cdcTree->GetEntry( CDC_Event_Number );
  if( header_CDC->ev()!=Event_Number ) return true;

  CDC_Event_Number++;
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  //** copy class CDSTrackingMan objects **//
  *trackMan = *trackMan_CDC; 

  //** Converting data=>tree and re-calc. of trackMan of CDS **//
  header->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  trackMan->Calc( cdsMan, confMan, true);  

  int nGoodTrack = trackMan->nGoodTrack();
  int nallTrack = trackMan->nTrack();
  AllGoodTrack += nGoodTrack;
  nTrack += nallTrack;
  Tools::Fill1D( Form("nGoodTrack"), nGoodTrack );

  //** # of good CDS tracks cut **//
  //allowing =>3 tracks ? OK? 
  if( nGoodTrack<2 ){//require pi+,pi-
    Clear( nAbort_nGoodTrack );
    return true;
  }

  //CDH-hits cut
  if(!EveSelectCDHMul()) return true;

  //beam line event selection
  if(!EveSelectBeamline()) return true;
  
  //PID beam
  if(PIDBeam!=Beam_Kaon) return true;

  //BLC1-D5-BLC2 analysis and chi2 selection
  double beammom = AnaBeamSpec(ConfMan *confMan);
  if(beammom <-100) return true;//chi2 cut

  //** beam momentum calculation **//
  TVector3 Pp_target; Pp_target.SetXYZ( 0, 0, 0 ); 
  TLorentzVector LVec_beambf;  // 4-Momentum(beam) in LAB
  TLorentzVector LVec_beam;    // 4-Momentum(beam) in LAB with dE correcion
  TLorentzVector LVec_target;  // 4-Momentum(deuteron-target) in LAB
  LVec_target.SetVectM( Pp_target, dMass );//deuteron target
  TLorentzVector LVec_targetP; // 4-Momentum(p-target) in LAB
  LVec_targetP.SetVectM( Pp_target, pMass );//H target
  TLorentzVector LVec_beambfCM;  // 4-Momentum(beam) in CM
  TLorentzVector LVec_beamCM;    // 4-Momentum(beam) in CM with dE correcion
  TLorentzVector LVec_targetCM;  // 4-Momentum(deuteron-target) in CM
  TLorentzVector LVec_targetPCM; // 4-Momentum(p-target) in CM
  
  double x1, y1, x2, y2;
  double z1 = 0, z2 = 20;
  LocalTrack *bpctrack = bltrackMan->trackBPC(bpcGoodTrackID);
  bpctrack->XYPosatZ( z1, x1, y1 );		  
  bpctrack->XYPosatZ( z2, x2, y2 );
  TVector3 lp; lp.SetXYZ( x1, y1, z1 );
  TVector3 ls; ls.SetXYZ( x2-x1, y2-y1, z2-z1); ls = ls.Unit();	       
  TVector3 Pp_beam = beammom*ls; 

  LVec_beambf.SetVectM( Pp_beam , kpMass );
  LVec_beam = LVec_beambf;
  TVector3 boost = (LVec_target+LVec_beam).BoostVector();
  LVec_beambfCM = LVec_beam;
  LVec_targetCM = LVec_target;
  LVec_targetPCM = LVec_targetP;
  //boost to CM frame
  LVec_beambfCM.Boost( -1*boost );
  LVec_targetCM.Boost( -1*boost );
  LVec_targetPCM.Boost( -1*boost );


  //** + + + + + + + + + + + + **//
  //**  PID in CDS             **//
  //** + + + + + + + + + + + + **//
  int CDHseg=-1;
  //** vectors for PID container **//
  std::vector <int> pip_ID;
  std::vector <int> pim_ID;
  std::vector <int> km_ID;
  std::vector <int> p_ID;
  std::vector <int> d_ID;

  std::vector <int> spip_ID;
  std::vector <int> spim_ID;
  std::vector <int> skm_ID;
  std::vector <int> sp_ID;

  std::vector <int> vCDHseg;


  //** PID of CDS tracks **//
  for( int it=0; it<trackMan->nGoodTrack(); it++ ){
    CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );

    //** chi2 cut can be applied in CDSTrackingMan with MaxChi in CDSFittingParam_posi.param **//
    Tools::Fill1D( Form("trackchi2_CDC"), track->Chi() );
    if( track->Chi()>cdscuts::cds_chi2_max ) continue; 
    
    //asano memo
    //checking if CDH hits at the projected position of the track ?
    if( !track->CDHFlag() ) continue;

    double mom = track->Momentum();
    TVector3 vtxbline, vtxhelix, vtxb;
    track->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtxbline, vtxbhelix );
    track->SetPID(-1);
    vtxb = (vtxbline+vtxbhelix)*0.5;//not used, so far

    double tof = 999.;
    double mass2 = -999.;
    //asano memo
    //perhaps, needs clustering ?
    int nCDHass = track->nCDHHit();
    Tools::Fill(Form("mul_CDH_assoc",nCDHass);
    for( int icdh=0; icdh<track->nCDHHit(); icdh++ ){
      HodoscopeLikeHit *cdhhit = track->CDHHit( cdsMan, icdh );
      double tmptof = cdhhit->ctmean()-ctmT0;
      if( tmptof<tof || tof==999. ){
        tof = tmptof;
        CDHseg = cdhhit->seg();
      }
    }
    
    //asano memo
    //AYASHII ?
    //check if the segment has been already used in the analysis.
    bool CDHflag = true;
    for( int icdhseg=0; icdhseg<(int)vCDHseg.size(); icdhseg++ ){
      if( CDHseg==vCDHseg[icdhseg] ) CDHflag = false;
    }
    if( !CDHflag ) continue;//go to next CDStrack
    vCDHseg.push_back( CDHseg );

    //** calculation of beta and squared-mass **//
    double tmptof, beta_calc;
    if( !TrackTools::FindMass2( track, bpctrack, tof, LVec_beam.Vect().Mag(),
				PIDBeam, beta_calc, mass2, tmptof ) ){
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      continue;
    }

    //** Retiming of CDC track by CDH info. **//
    if(DoCDCRetiming){
      track->Retiming( cdsMan, confMan, beta_calc, true );
      //asano memo
      //why need this ?
      for( int m=0; m<5; m++ ){
        track->HelixFitting( cdsMan ); //** 5 times iteration **//
      }
      track->Calc( confMan );
    }


    //** finalize PID **//
    if( !TrackTools::FindMass2( track, bpctrack, tof, LVec_beam.Vect().Mag(),
				PIDBeam, beta_calc, mass2, tmptof ) ){ //** not FindMass2C() [20170622] **//
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      continue;
    }
    int pid = -1;

    //if( GRUN=="49c" ){
    //  pid = TrackTools::PIDcorr( mom, mass2 ); //!! sada-D p.92 !!//
    //}else if ( GRUN=="65" ){
    //  pid = TrackTools::PIDcorr3( mom, mass2 ); //** from Yamaga library **// //!! Yamaga-D v03 !!//
    //}

    //RUN68 inoue tune
    pid = TrackTools::PIDcorr_wide(mom,mass2);

    track->SetPID( pid );
    Tools::Fill2D( "PID_CDS_beta", 1/beta_calc, mom );
    Tools::Fill2D( "PID_CDS", mass2, mom );

    //** energy loss calculation **//
    double tmpl;
    TVector3 vtx_beam, vtx_cds;
    if( !track->CalcVertexTimeLength( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), track->Mass(),
				      vtx_beam, vtx_cds, tmptof, tmpl, true ) ){
      std::cerr<<" !!! failure in energy loss calculation [CalcVertexTimeLength()] !!! "<<std::endl;
      continue;
    }

    if( pid==CDS_PiMinus ){
      pim_ID.push_back( trackMan->GoodTrackID(it) );
    }else if( pid==CDS_PiPlus ){
      pip_ID.push_back( trackMan->GoodTrackID(it) );
    }else if( pid==CDS_Proton ){
      p_ID.push_back( trackMan->GoodTrackID(it) );
    }else if( pid==CDS_Deuteron ){
      d_ID.push_back( trackMan->GoodTrackID(it) );
    }else if( pid==CDS_Kaon ){
      km_ID.push_back( trackMan->GoodTrackID(it) );
    }
  } // for( int it=0; it<trackMan->nGoodTrack(); it++ ){
  //** end of PID (except for neutron) **//

  Tools::Fill1D( Form("ntrack_CDS"), pip_ID.size()+p_ID.size()+d_ID.size()+pim_ID.size()+km_ID.size() );
  Tools::Fill1D( Form("ntrack_pi_plus"),  pip_ID.size() );
  Tools::Fill1D( Form("ntrack_proton"),   p_ID.size() );
  Tools::Fill1D( Form("ntrack_deuteron"), d_ID.size() );
  Tools::Fill1D( Form("ntrack_pi_minus"), pim_ID.size() );
  Tools::Fill1D( Form("ntrack_K_minus"),  km_ID.size() );
         
  //** charge veto with BVC, CVC (TOF=CVC), & PC **//
  int nBVC = 0;
  int nCVC = 0;
  int nPC  = 0;
  for( int i=0; i<blMan->nBVC(); i++ ){
    if( blMan->BVC(i)->CheckRange() ) nBVC++;
  }
  for( int i=0; i<blMan->nTOF(); i++ ){
    if( blMan->TOF(i)->CheckRange() ) nCVC++;
  }
  for( int i=0; i<blMan->nPC(); i++ ){
    if( blMan->PC(i)->CheckRange() ) nPC++;
  }
  Tools::Fill1D( Form("mul_BVC"), nBVC );
  Tools::Fill1D( Form("mul_CVC"), nCVC );
  Tools::Fill1D( Form("mul_PC"),  nPC );
  bool chargedhit = false;
  if( nBVC || nCVC || nPC ) chargedhit = true; 
  //if( GRUN=="49c" ){
  //  if( nBVC ) chargedhit = true;
  //}else if( GRUN=="65" ){
  //  if( nBVC || nCVC || nPC ) chargedhit = true;
  //}



  //TLorentzVector Lpipdef, Lpimdef;
  TVector3 vtx_react;
  //  pi+ pi- X event  
  //  with CDH multiplicity selection
  if( pim_ID.size()==1 && pip_ID.size()==1 &&
      trackMan->nGoodTrack()>1 && !chargedhit ){

    //=== condense pi+ pi- X candidates ===//
    rtFile2->cd();
    evTree->Fill();
    rtFile->cd();
    std::cout<<"### filled: Event_Number, Block_Event_Number, CDC_Event_Number = "
	     <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
    nFill_pippim++;
    //=== condense pi+ pi- X event ===//
    
    //** find CDH hit from neutral particles **//
    std::vector <int> NeutralCDHseg;//CDHhits - CDHhits used for charged particle tracking
    std::vector <int> CDHhit_list;
    for( int icdhhit=0; icdhhit<cdsMan->nCDH(); icdhhit++ ){
      if( cdsMan->CDH(icdhhit)->CheckRange() && 
          cdsMan->CDH(icdhhit)->ctmean()< cdscuts::tdc_cdh_max){
        CDHhit_list.push_back( cdsMan->CDH(icdhhit)->seg() );
      }
    }
    std::sort(vCDHseg.begin(), vCDHseg.end());
    std::sort(CDHhit_list.begin(), CDHhit_list.end());
    std::set_difference( CDHhit_list.begin(), CDHhit_list.end(),
			 vCDHseg.begin(), vCDHseg.end(),
			 std::back_inserter(NeutralCDHseg) );

    if(Verbosity){
      if( NeutralCDHseg.size()!=1 ){
        std::cerr<<" CDH neutral hit is not 1 :: "<<NeutralCDHseg.size()<<std::endl;
      }
      std::cerr<<"# of diff = "<<NeutralCDHseg.size()<<std::endl;
      std::cerr<<"CDH hits =   ";
      for( int n=0; n<(int)CDHhit_list.size(); n++ ){
        std::cerr<<CDHhit_list[n]<<" ";
      } std::cerr<<std::endl;
      std::cerr<<"track hits = ";
      for( int n=0; n<(int)vCDHseg.size(); n++ ){
        std::cerr<<vCDHseg[n]<<" ";
      } std::cerr<<std::endl;
      std::cerr<<"diff hits =  ";
      for( int n=0; n<(int)NeutralCDHseg.size(); n++ ){
        std::cerr<<NeutralCDHseg[n]<<" ";
      } std::cerr<<std::endl;
    }

    //** isolation cut **//
    int flag_isolation = 0;
    for( int ineuseg=0; ineuseg<(int)NeutralCDHseg.size(); ineuseg++ ){
      for( int ihit=0; ihit<(int)CDHhit_list.size(); ihit++ ){
	if( NeutralCDHseg[ineuseg]-CDHhit_list[ihit] ){
    Tools::Fill1D( Form("diff_CDH"), NeutralCDHseg[ineuseg]-CDHhit_list[ihit] );
  }
  //CDH has 36 segments. Requiring there is no hits on neighboring segments.
	if( abs(NeutralCDHseg[ineuseg]-CDHhit_list[m])==1 || abs(NeutralCDHseg[ineuseg]-CDHhit_list[m])==35 )
	  flag_isolation++;
      }
    }
    if( flag_isolation ){
      if(Verbosity) std::cerr<<"CDH hit candidate is NOT isolated !!!"<<std::endl;
      Clear( nAbort_CDHiso );
      return true;
    }

    //** copy neutral CDH hit candidate **//
    int icdh = -1;
    for( int ihit=0; ihit<cdsMan->nCDH(); ihit++ ){
      if( cdsMan->CDH(ihit)->seg()==NeutralCDHseg[0] ) icdh = ihit;
    }
    HodoscopeLikeHit *ncdhhit = cdsMan->CDH(icdh);

    //** charge veto using CDC **//
    TVector3 Pos_CDH;
    confMan->GetGeomMapManager()->GetPos( CID_CDH, ncdhhit->seg(), Pos_CDH );
    if(Verbosity){
      std::cerr<<"CDH candidate = "<<ncdhhit->seg()<<" -> "<<Pos_CDH.Phi()/TwoPi*360<<" deg"<<std::endl;
    }
    const double PhiMin = -15.0/360.*TwoPi; // rad
    const double PhiMax =  15.0/360.*TwoPi; // rad
    if(Verbosity){
      std::cerr<<"Min/Max = "<<PhiMin/TwoPi*360<<"/"<<PhiMax/TwoPi*360<<" deg"<<std::endl;
    }
    int nCDC = 0;
    for( int ilr=14; ilr<16; ilr++ ){ // charge veto using layer 14, 15
      for( int icdchit=0; icdchit<cdsMan->nCDC(ilr); icdchit++ ){
        CDCHit *cdc=cdsMan->CDC(ilr,icdchit);
        TVector3 Pos_CDC = cdc->wpos();
        Pos_CDC.SetZ(0); // only xy pos is used
        double angle = Pos_CDC.Angle(Pos_CDH); // rad
        if(Verbosity){std::cerr<<"CDC "<<ilr<<" "<<icdchit<<" "<<cdc->wire()<<" -> "<<Pos_CDC.Phi()/TwoPi*360.
          <<" deg :: diff = "<<angle/TwoPi*360<<" deg"<<std::endl;
        }
        Tools::Fill1D( Form("diff_CDH_CDC"), angle/TwoPi*360 );
        if( PhiMin<angle && angle<PhiMax ) nCDC++;
      }//icdchit
    }//ilr
    if(Verbosity){
      std::cerr<<"# of CDC hits for nCDH candidate = "<<nCDC<<std::endl;
    }
    Pos_CDH.SetZ(-1*ncdhhit->hitpos()); // (-1*) is correct in data analysis [20170926]
    //Pos_CDH.SetZ(ncdhhit->hitpos());


    //** neutral particle in CDH **//
    if( !nCDC ){
      CDSTrack *track_pip = trackMan->Track( pip_ID[0] ); // only 1 track
      CDSTrack *track_pim = trackMan->Track( pim_ID[0] ); // only 1 track

      //He3 analysis vertex calculation
      //TVector3 vtx_b; // Vertex(baem-particle)_on_beam
      //TVector3 vtx_p; // Vertex(baem-particle)_on_particle
      //#if VTX_SIGMA
      //      track_p->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      //#else
      //      track_pip->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      //#endif
      //track_p->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      //vtx_react = 0.5*(vtx_b+vtx_p); // reaction vertex
       
      //deuteron target
      TVector3 vtx_beam_wpip;//vertex(beam-pip) on beam
      TVector3 vtx_pip;//vertex(beam-pip) on beam
      track_pip->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wpip, vtx_pip ); 
      TVector3 vtx_beam_wpim;//vertex(beam-pip) on beam
      TVector3 vtx_pim;//vertex(beam-pip) on beam
      track_pim->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wpim, vtx_pim ); 
      
      double dcapipvtx =  (vtx_pip-vtx_beam_wpip).Mag();
      double dcapimvtx =  (vtx_pim-vtx_beam_wpim).Mag();
      

      //reaction vertex is determined from beam and nearest vtx 
      if(dcapipvtx <= dcapimvtx) vtx_react = 0.5*(vtx_pip+vtx_beam_wpip);
      else if (dcapipvtx > dcapimvtx) vtx_react = 0.5*(vtx_pim+vtx_beam_wpim);
      //vertex position from pi+/pi-
      //bool vtx_flag=TrackTools::Calc2HelixVertex(track_pip, track_pim, vtx1, vtx2);

      //why need this ?
      //double tof = 999.;
      //for( int icdh=0; icdh<track_p->nCDHHit(); icdh++ ){
      //  HodoscopeLikeHit *cdhhit = track_p->CDHHit( cdsMan, icdh );
      //  double tmptof = cdhhit->ctmean()-ctmT0;
      //  if( tmptof<tof || tof==999. ){
      //    tof = tmptof;
      //    CDHseg = cdhhit->seg();
      //  }
      //}
      //why need this ?

      //** beam kaon tof **//
      TVector3 Pos_T0;
      confMan->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
      double beamtof, momout;
      double z_pos = Pos_T0.Z();;
      //dE correction of beam  
      ELossTools::CalcElossBeamTGeo( bpctrack->GetPosatZ(z_pos), vtx_react,
				     LVec_beambf.Vect().Mag(), kpMass, momout, beamtof );
      LVec_beam.SetVectM( momout*LVec_beambf.Vect().Unit(), kpMass );
      double ntof = ncdhhit->ctmean()-ctmT0-beamtof;
      double nlen = (Pos_CDH-vtx_react).Mag();
      NeutralBetaCDH = nlen/ntof/(Const*100.);
      double tmp_mom = NeutralBetaCDH<1. ? nMass*NeutralBetaCDH/sqrt(1.-NeutralBetaCDH*NeutralBetaCDH) : 0;
      if(Verbosity) std::cerr<<"$$$ NeutralBetaCDH = "<<NeutralBetaCDH<<" mom_n = "<<tmp_mom<<std::endl; //" "<<1/sqrt(1+nMass*nMass)<<std::endl;

      //** reconstructoin of missing neutorn **//
      //TVector3 P_p;   // Momentum(p)
      TVector3 P_pim; // Momentum(pi-)
      TVector3 P_pip; // Momentum(pi+)
      TVector3 P_n;   // Momentum(n) (CDS)
      
      //TLorentzVector L_p;   // 4-Momentum(p)
      TLorentzVector L_pim; // 4-Momentum(pi-)
      TLorentzVector L_pip; // 4-Momentum(pi+)
      TLorentzVector L_n;   // 4-Momentum(n) (CDS)
      TLorentzVector L_nmiss; // 4-Momentum(n_miss)

      //track_p->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      //double dca_p  = (vtx_b-vtx_p).Mag(); // DCA(beam-p)
      //if( !track_p->GetMomentum( vtx_p, P_p, true, true ) ){
      //	std::cerr<<" !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      // }
      //track_pip->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      //double dca_pip  = (vtx_b-vtx_p).Mag(); // DCA(beam-pip)
      if( !track_pip->GetMomentum( vtx_pip, P_pip, true, true ) ){
        std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }
      //track_pim->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      //double dca_pim  = (vtx_b-vtx_p).Mag(); // DCA(beam-pim)
      if( !track_pim->GetMomentum( vtx_pim, P_pim, true, true ) ){
        std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }
      P_n = tmp_mom*(Pos_CDH-vtx_react).Unit();
      
      L_p.SetVectM(   P_p,   pMass );
      L_pim.SetVectM( P_pim, piMass );
      L_pip.SetVectM( P_pip, piMass );
      L_n.SetVectM(   P_n,   nMass );
      
      double mm_mass   = (LVec_target+LVec_beam-L_p-L_pim-L_pip-L_n).M();
      TVector3 P_missn = (LVec_target+LVec_beam-L_p-L_pim-L_pip-L_n).Vect();
      L_nmiss.SetVectM( P_missn, nMass );
      std::cerr<<"  missing mass = "<<mm_mass<<std::endl;

      TVector3 boost = (LVec_target+LVec_beam).BoostVector();
      TLorentzVector L_nmiss_CM = L_nmiss;
      TLorentzVector LVec_beam_CM = LVec_beam;
      L_nmiss_CM.Boost(-boost);
      LVec_beam_CM.Boost(-boost);
      double cos_n = L_nmiss_CM.Vect().Dot(LVec_beam_CM.Vect())/(L_nmiss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
      std::cerr<<"  missing mom | cos_CM = "<<cos_n<<std::endl;


      //** + + + + + + + + + + + + + **//
      //**  fill histograms & tree   **//
      //** + + + + + + + + + + + + + **//

      const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
      const double dE_MIN = 5.0; // 8.0MeVee * 3cm / 5cm;

      const double pipi_MIN = 0.485;
      const double pipi_MAX = 0.510;
      const double ppi_MIN = 1.1075;
      const double ppi_MAX = 1.1225;

      const double neutron_MIN = 0.85;
      const double neutron_MAX = 1.03;

      const double Sigmap_MIN = 1.18;
      const double Sigmap_MAX = 1.20;
      const double Sigmam_MIN = 1.19;
      const double Sigmam_MAX = 1.21;

      kf_flag = -1;

      Tools::Fill2D( Form("dE_betainv"), 1/NeutralBetaCDH, ncdhhit->emean() );
      Tools::Fill2D( Form("MMom_MMass"), mm_mass, P_missn.Mag() );
      
      if( GeomTools::GetID(vtx_react)==CID_Fiducial ){
	Tools::Fill2D( Form("dE_betainv_fiducial"), 1/NeutralBetaCDH, ncdhhit->emean() );
	Tools::Fill2D( Form("MMom_MMass_fiducial"), mm_mass, P_missn.Mag() );

	if(  NeutralBetaCDH<beta_MAX ){
	  Tools::Fill2D( Form("dE_betainv_fiducial_beta"), 1/NeutralBetaCDH, ncdhhit->emean() );
	  Tools::Fill2D( Form("MMom_MMass_fiducial_beta"), mm_mass, P_missn.Mag() );

	  if( dE_MIN<ncdhhit->emean() ){
	    Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE"), 1/NeutralBetaCDH, ncdhhit->emean() );
	    Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE"), mm_mass, P_missn.Mag() );

	    Tools::Fill1D( Form("IMpipi"), (L_pim+L_pip).M() );
	    Tools::Fill1D( Form("IMppi"), (L_p+L_pim).M() );

	    if( ((L_pim+L_pip).M()<pipi_MIN || pipi_MAX<(L_pim+L_pip).M()) &&
		((L_p+L_pim).M()<ppi_MIN || ppi_MAX<(L_p+L_pim).M()) ){ // K0 & Lambda subtraction
	      Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE_res"), 1/NeutralBetaCDH, ncdhhit->emean() );
	      Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE_res"), mm_mass, P_missn.Mag() );

	      if( neutron_MIN<mm_mass && mm_mass<neutron_MAX ){
		Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE_res_n"), 1/NeutralBetaCDH, ncdhhit->emean() );
		Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE_res_n"), mm_mass, P_missn.Mag() );

		Tools::Fill2D( Form("MMom_NMom"), P_n.Mag(), P_missn.Mag() );
		Tools::Fill2D( Form("IMnpim_IMnpip"), (L_n+L_pip).M(), (L_n+L_pim).M() );
		
		if( (Sigmap_MIN<(L_n+L_pip).M() && (L_n+L_pip).M()<Sigmap_MAX) ||
		    (Sigmam_MIN<(L_n+L_pim).M() && (L_n+L_pim).M()<Sigmam_MAX) ){
		  Tools::Fill2D( Form("IMmnpim_IMmnpip"), (L_nmiss+L_pip).M(), (L_nmiss+L_pim).M() );
		  Tools::Fill2D( Form("MMnppip_MMnppim"), (LVec_target+LVec_beam-L_p-L_pim-L_n).M(),
				 (LVec_target+LVec_beam-L_p-L_pip-L_n).M() );
		  
		  Tools::Fill2D( Form("Cosn_IMnppipi"), (L_n+L_p+L_pim+L_pip).M(), cos_n );
		  Tools::Fill2D( Form("Cosn_IMnpipi"), (L_n+L_pim+L_pip).M(), cos_n );
		  Tools::Fill2D( Form("IMnpipi_IMnppipi"), (L_n+L_p+L_pim+L_pip).M(), (L_n+L_pim+L_pip).M() );
				  
		  Tools::Fill1D( Form("DCA_p"), dca_p );
		  Tools::Fill1D( Form("DCA_pip"), dca_pip );
		  Tools::Fill1D( Form("DCA_pim"), dca_pim );
		}
	      } // if( neutron_MIN<mm_mass && mm_mass<neutron_MAX ){
	      
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
		// %%% Kinematical Fit using KinFitter %%% //
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
		//--- set TLorentzVector ---//
		// beam_K(K+), pi-/+, Sigma+/-, p, n, n from S, pi+/- from S 
		//  = 1) TLorentzVector LVec_beam, L_pim, (L_n+L_pip), L_p, L_nmiss, L_n, L_pip = for pi- Sigma+
		TLorentzVector TL_meas1[7]; // measured
		TL_meas1[0] = LVec_beam;
		TL_meas1[1] = L_pim;
		TL_meas1[2] = (L_n+L_pip);
		TL_meas1[3] = L_p;
		TL_meas1[4] = L_nmiss;
		TL_meas1[5] = L_n;
		TL_meas1[6] = L_pip;
		//  = 2) TLorentzVector LVec_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
		TLorentzVector TL_meas2[7]; // measured
		TL_meas2[0] = LVec_beam;
		TL_meas2[1] = L_pip;
		TL_meas2[2] = (L_n+L_pim);
		TL_meas2[3] = L_p;
		TL_meas2[4] = L_nmiss;
		TL_meas2[5] = L_n;
		TL_meas2[6] = L_pim;
		TLorentzVector TL_kfit1[7]; // kinematical fitted
		TLorentzVector TL_kfit2[7]; // kinematical fitted
		// LVec_target is defined as (0, 0, 0, M_3He)
		TVector3 TV_target = LVec_target.Vect();
		TVector3 TV_meas1[7];
		TVector3 TV_meas2[7];
		for( int i=0; i<7; i++ ){
		  TV_meas1[i] = TL_meas1[i].Vect();
		  TV_meas2[i] = TL_meas2[i].Vect();
		}

		TDatabasePDG *pdg = new TDatabasePDG();
		pdg->ReadPDGTable("pdg_table.txt");
		int PDG1[7] = {321, -211, 3222, 2212, 2112, 2112,  211}; // pi-Sigma+
		int PDG2[7] = {321,  211, 3112, 2212, 2112, 2112, -211}; // pi+Sigma-
		
		//--- KinFitter :: initialization ---//
		//  = 1) TLorentzVector LVec_beam, L_pim, (L_n+L_pip), L_p, L_nmiss, L_n, L_pip = for pi- Sigma+
		//  = 2) TLorentzVector LVec_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
		//*** definition of fit particles in cartesian coordinates ***//
		TString str_particle1[7] = {"L_beam", "L_pim", "L_Sp", "L_p", "L_mn", "L_n", "L_pip"};
		TString str_particle2[7] = {"L_beam", "L_pip", "L_Sm", "L_p", "L_mn", "L_n", "L_pim"};
		TFitParticlePxPyPz ParticleTgt = TFitParticlePxPyPz("target", "target", &TV_target,
								    pdg->GetParticle("He3")->Mass(), covZero);
		TFitParticlePxPyPz Particle1[7];
		TFitParticlePxPyPz Particle2[7];
		for( int i=0; i<7; i++ ){
		  Particle1[i] = TFitParticlePxPyPz(str_particle1[i], str_particle1[i], &TV_meas1[i],
						    pdg->GetParticle(PDG1[i])->Mass(), covParticle1[i]);
		  Particle2[i] = TFitParticlePxPyPz(str_particle2[i], str_particle2[i], &TV_meas2[i],
						    pdg->GetParticle(PDG2[i])->Mass(), covParticle2[i]);
		}
		//*** definition of constraints ***//
		// constraint :: mass of Sigma
		TFitConstraintM ConstMS1 = TFitConstraintM("M_Sp", "M_Sp", 0, 0, pdg->GetParticle(PDG1[2])->Mass());
		TFitConstraintM ConstMS2 = TFitConstraintM("M_Sm", "M_Sm", 0, 0, pdg->GetParticle(PDG2[2])->Mass());
		ConstMS1.addParticles1(&Particle1[5], &Particle1[6]);
		ConstMS2.addParticles1(&Particle2[5], &Particle2[6]);
		// constraint :: 4-momentum conservation
		TFitConstraintEp ConstEp1[4];
		TFitConstraintEp ConstEp2[4];
		TString str_constEp1[4]  = {"Px", "Py", "Pz", "E"};
		TString str_constEp2[4]  = {"Px", "Py", "Pz", "E"};
		for( int i=0; i<4; i++ ){
		  ConstEp1[i] = TFitConstraintEp(str_constEp1[i], str_constEp1[i], 0, TFitConstraintEp::component(i), 0);
		  ConstEp2[i] = TFitConstraintEp(str_constEp2[i], str_constEp2[i], 0, TFitConstraintEp::component(i), 0);
		  ConstEp1[i].addParticles1(&ParticleTgt, &Particle1[0]);
		  ConstEp2[i].addParticles1(&ParticleTgt, &Particle2[0]);
		  ConstEp1[i].addParticles2(&Particle1[1], &Particle1[3], &Particle1[4], &Particle1[5], &Particle1[6]);
		  ConstEp2[i].addParticles2(&Particle2[1], &Particle2[3], &Particle2[4], &Particle2[5], &Particle2[6]);
		}

		//--- KinFitter :: execution ---//
	    else if( ((L_pim+L_pip).M()<pipi_MIN || pipi_MAX<(L_pim+L_pip).M()) &&
		     ((ppi_MIN<(L_p+L_pim).M() && (L_p+L_pim).M()<ppi_MAX)) ){ // K0 subtraction & Lambda selection
	      //** reaction vertex would be obtained from K- & pi+, thus this analysis is temporary **//
	      Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE_Lambda"), 1/NeutralBetaCDH, ncdhhit->emean() );
	      Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE_Lambda"), mm_mass, P_missn.Mag() );

	      if( neutron_MIN<mm_mass && mm_mass<neutron_MAX ){
		Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE_Lambda_n"), 1/NeutralBetaCDH, ncdhhit->emean() );
		Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE_Lambda_n"), mm_mass, P_missn.Mag() );
		
		Tools::Fill2D( Form("MMom_NMom_Lambda"), P_n.Mag(), P_missn.Mag() );
		Tools::Fill2D( Form("Cosn_IMnppipi_Lambda"), (L_n+L_p+L_pim+L_pip).M(), cos_n );
		Tools::Fill2D( Form("Cosn_IMppipi_Lambda"), (L_p+L_pim+L_pip).M(), cos_n );
		Tools::Fill2D( Form("IMppipi_IMnppipi_Lambda"), (L_n+L_p+L_pim+L_pip).M(), (L_p+L_pim+L_pip).M() );

		Tools::Fill1D( Form("DCA_p_Lambda"), dca_p );
		Tools::Fill1D( Form("DCA_pip_Lambda"), dca_pip );
		Tools::Fill1D( Form("DCA_pim_Lambda"), dca_pim );
	      }
	    } // else if( ppi_MIN<(L_p+L_pim).M() && (L_p+L_pim).M()<ppi_MAX ){

	  } // if( dE_MIN<ncdhhit->emean() ){
	} // if(  beta<beta_MAX ){

	//** fill tree **//
	mom_beam   = LVec_beam;   // 4-momentum(beam)
	mom_target = LVec_target; // 4-momentum(target)
	mom_pip = L_pip;        // 4-momentum(pi+)
	mom_pim = L_pim;        // 4-momentum(pi-)
	mom_p = L_p;            // 4-momentum(proton)
	mom_n = L_n;            // 4-momentum(neutron)
	dE = ncdhhit->emean();
	// beta is already filled
	vtx_reaction = vtx_react; // vertex(reaction)
	run_num   = confMan->GetRunNumber(); // run number
	event_num = Event_Number;            // event number
	block_num = Block_Event_Number;      // block number

	std::cout<<"%%% npippim event: Event_Number, Block_Event_Number, CDC_Event_Number = "
		 <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
	rtFile3->cd();
	npippimTree->Fill();
	rtFile->cd();
	nFill_npippim++;
	//** fill tree **//

      } // if( GeomTools::GetID(vtx_react)==CID_Fiducial ){
    } // if( !nCDC ){
  }
  else{
    Clear( nAbort_ppipi );
    return true;
  }

  Clear( nAbort_end );
  return true;
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
void EventAnalysis::Finalize()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  std::cout << " *** Enter EventAnalysis::Finalize " << std::endl;

  rtFile2->Write();
  rtFile2->Close();

  rtFile3->Write();
  rtFile3->Close();

  rtFile->cd();
  for(int i=0;i<40;i++){
    TH1F* h1 = (TH1F*)gFile->Get("Scaler"); 
    std::cout<<i<<"  "<<scaend[i]<<std::endl;
    h1->Fill(i,scaend[i]-scainit[i]);
  }

  std::cout<<"====== Abort counter ========="<<std::endl;
  std::cout<<" nAbort_nGoodTrack    = "<<nAbort_nGoodTrack<<std::endl;
  std::cout<<" nAbort_nCDH          = "<<nAbort_nCDH<<std::endl;
  std::cout<<" nAbort_nT0           = "<<nAbort_nT0<<std::endl;
  std::cout<<" nAbort_pid_beam      = "<<nAbort_pid_beam<<std::endl;
  std::cout<<" nAbort_nbpc          = "<<nAbort_nbpc<<std::endl;
  std::cout<<" nAbort_bpctrack      = "<<nAbort_bpctrack<<std::endl;
  std::cout<<" nAbort_singleBLtrack = "<<nAbort_singleBLtrack<<std::endl;
  std::cout<<" nAbort_fblc2bpc      = "<<nAbort_fblc2bpc<<std::endl;
  std::cout<<" nAbort_flagbmom      = "<<nAbort_flagbmom<<std::endl;
  std::cout<<" nAbort_ftarget       = "<<nAbort_ftarget<<std::endl;
  std::cout<<" nAbort_CDHiso        = "<<nAbort_CDHiso<<std::endl;
  std::cout<<" nAbort_nAbort_ppipi  = "<<nAbort_ppipi<<std::endl;
  std::cout<<" nAbort_end           = "<<nAbort_end<<std::endl;
  std::cout<<"========= Abort counter ========="<<std::endl;
  std::cout<<"*** # of pi+ pi- p n n events = "<<nFill_npippim<<" ***"<<std::endl;

  //  confMan->SaveCDSParam();
  gFile->Write();
  gFile->Close();
  cdcFile->Close();

  delete cdsMan;
  delete header;
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
void EventAnalysis::InitializeHistogram()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  //** geneneral informantion **//
  Tools::newTH1F( Form("Time"), 3000, -0.5, 2999.5 );
  Tools::newTH1F( Form("EventCheck"), 20, 0, 20 );
  Tools::newTH1F( Form("Scaler"), 41, -0.5, 40.5 );

  //** CDC and CDH information from CDC-trackig file **//
  Tools::newTH1F( Form("nGoodTrack"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_CDH"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_CDH_assoc"), 11, -0.5, 10.5 );
  Tools::newTH2F( Form("CDHtime"), 36, -0.5, 35.5,1600,0,40);
 

  //** beam line **//
  Tools::newTH1F( Form("mul_BHD"), 12, -0.5, 11.5 );
  Tools::newTH1F( Form("mul_T0"),   6, -0.5, 5.5 );
  Tools::newTH1F( Form("tof_T0BHD"), 2000, 20, 40 );
  Tools::newTH1F( Form("tracktime_BPC"),  1200, -200, 400 );
  Tools::newTH1F( Form("trackchi2_BPC"),  200, 0, 20 );
  Tools::newTH1F( Form("ntrack_BPC"),  6, -0.5, 5.5 );
  Tools::newTH1F( Form("tracktime_BLC1"), 1200, -200, 400 );
  Tools::newTH1F( Form("tracktime_BLC2"), 1200, -200, 400 );
  Tools::newTH1F( Form("trackchi2_BLC1"), 200, 0, 20 );
  Tools::newTH1F( Form("trackchi2_BLC2"), 200, 0, 20 );
  Tools::newTH1F( Form("ntrack_BLC1"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_BLC2"), 6, -0.5, 5.5 );
  Tools::newTH2F( Form("dydx_BLC2BPC"),     130, -1.3, 1.3, 130, -1.3, 1.3 );
  Tools::newTH2F( Form("dydzdxdz_BLC2BPC"), 175, -0.035, 0.035, 175, -0.035, 0.035 );
  Tools::newTH1F( Form("trackchi2_beam"), 400, 0, 40 );
  Tools::newTH1F( Form("momentum_beam"), 180, 0.92, 1.10 );
  Tools::newTH1F( Form("PID_beam"),5,-1.5,3.5);

  //** CDS **//
  Tools::newTH1F( Form("trackchi2_CDC"), 1000, 0, 50 );
  Tools::newTH2F( Form("PID_CDS_beta"), 1000, 0, 5, 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS"), 1000, -0.6, 5, 1000, -1.2, 1.2 );
  Tools::newTH1F( Form("ntrack_CDS"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_pi_plus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_proton"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_deuteron"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_pi_minus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_K_minus"), 6, -0.5, 5.5 );

  //** forward counters **//
  Tools::newTH1F( Form("mul_BVC"), 9, -0.5, 8.5 );
  Tools::newTH1F( Form("mul_CVC"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_PC"), 11, -0.5, 10.5 );

  //** p pi+ pi- X event **//
  Tools::newTH1F( Form("diff_CDH"), 73, -36.5, 36.5 );
  Tools::newTH1F( Form("diff_CDH_CDC"), 180, 0, 180 );
  Tools::newTH2F( Form("dE_betainv"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fiducial"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fiducial_beta"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fiducial_beta_dE"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fiducial_beta_dE_res"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fiducial_beta_dE_res_n"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fiducial_beta_dE_Lambda"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fiducial_beta_dE_Lambda_n"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("MMom_MMass"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fiducial"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fiducial_beta"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fiducial_beta_dE"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fiducial_beta_dE_res"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fiducial_beta_dE_res_n"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fiducial_beta_dE_Lambda"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fiducial_beta_dE_Lambda_n"), 140, 0.4, 1.8, 100, 0, 1.5 );

  Tools::newTH1F( Form("IMpipi"), 40, 0.4, 0.6 );
  Tools::newTH1F( Form("IMppi"),  56, 1.06, 1.2 );

  Tools::newTH2F( Form("MMom_NMom"), 100, 0, 1.5, 100, 0, 1.5 );
  Tools::newTH2F( Form("IMnpim_IMnpip"), 70, 1, 1.7, 70, 1, 1.7 );
  Tools::newTH2F( Form("IMmnpim_IMmnpip"), 70, 1, 1.7, 70, 1, 1.7 );
  Tools::newTH2F( Form("MMnppip_MMnppim"), 70, 1, 1.7, 70, 1, 1.7 );
  Tools::newTH2F( Form("Cosn_IMnppipi"), 50, 2, 3, 50, -1, 1 );
  Tools::newTH2F( Form("Cosn_IMnpipi"), 100, 1, 2, 50, -1, 1 );
  Tools::newTH2F( Form("IMnpipi_IMnppipi"), 50, 2, 3, 100, 1, 2 );
  Tools::newTH1F( Form("DCA_p"), 200, 0, 2 );
  Tools::newTH1F( Form("DCA_pip"), 200, 0, 2 );
  Tools::newTH1F( Form("DCA_pim"), 200, 0, 2 );

  Tools::newTH2F( Form("MMom_NMom_Lambda"), 100, 0, 1.5, 100, 0, 1.5 );
  Tools::newTH2F( Form("Cosn_IMnppipi_Lambda"), 50, 2, 3, 50, -1, 1 );
  Tools::newTH2F( Form("Cosn_IMppipi_Lambda"), 100, 1, 2, 50, -1, 1 );
  Tools::newTH2F( Form("IMppipi_IMnppipi_Lambda"), 50, 2, 3, 100, 1, 2 );
  Tools::newTH1F( Form("DCA_p_Lambda"), 200, 0, 2 );
  Tools::newTH1F( Form("DCA_pip_Lambda"), 200, 0, 2 );
  Tools::newTH1F( Form("DCA_pim_Lambda"), 200, 0, 2 );

}

bool EvetnAnalysis::EveSelectCDHMul()
{
  //** # of CDH-hits cut **//
  int nCDH = 0;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    Tools::Fill2D(Form("CDHtime"),cdsMan->CDH(i)->seg(),cdsMan->CDH(i)->ctmean());
    //if( cdsMan->CDH(i)->CheckRange() ) nCDH++; //** only requirement of TDC **//
    if( cdsMan->CDH(i)->CheckRange() && cdsMan->CDH(i)->ctmean()<cdscuts::tdc_cdh_max ){
      nCDH++;
    }
  }
  Tools::Fill1D( Form("mul_CDH"), nCDH );
  if( nCDH != cdscuts::cdhmulti  ){ 
    Clear( nAbort_nCDH );
    return false;
  }

  return true;

}

//Beam PID + event selection of BLC1,2, BPC
bool EventAnalysis::EveSelectBeamline()
{
  //  beamline analysis & event selection      

  //** BLDC tracking **//
  bltrackMan->DoTracking(blMan,confMan,true,true);

  //** BHD & T0 **//
  int nBHD = 0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ) nBHD++;
  }
  int nT0 = 0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ) nT0++;
  }
  Tools::Fill1D( Form("mul_BHD"), nBHD );
  Tools::Fill1D( Form("mul_T0"),  nT0 );

  //** T0 = 1hit selection **//
  if( nT0!=1 ){  //!! sada-D p.72 !!//
    Clear( nAbort_nT0 );
    return false;      
  }

  //** Beam PID by T0-BHD TOF **//
  TVector3 vtxT0;
  double ctmT0 = 0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      ctmT0 = blMan->T0(i)->ctmean();
      confMan->GetGeomMapManager()->GetGPos(CID_T0, blMan->T0(i)->seg(), vtxT0);
    }
  }
  double ctmBHD;
  PIDBeam = -1; // 0:pi 1:K 3:else
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      ctmBHD = blMan->BHD(i)->ctmean();
      double tofBHDT0 = ctmT0-ctmBHD;
      Tools::Fill1D( Form("tof_T0BHD"), tofBHDT0 );      
      if( header->kaon() && blcuts::beam_tof_k_min <tofBHDT0 && tofBHDT0<blcuts::beam_tof_k_max  )
	PIDBeam = Beam_Kaon; 
      else if( header->pion() && blcuts::beam_tof_pi_min<tofBHDT0 && tofBHDT0<blcuts::beam_tof_pi_max )
	PIDBeam = Beam_Pion;
    }
  }
  Fill1D(Form("PID_beam"), PIDBeam);
  if( PIDBeam== -1  ){ //** unidentified particle is discarded (other than pi/K) **//
    Clear( nAbort_pid_beam );
    return false;
  }

  unsigned int nblc1 = 0;
  unsigned int nblc2 = 0;
  blc1GoodTrackID = -1;
  blc2GoodTrackID = -1;
  //** timing selection of BLC1/BLC2 **//
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    LocalTrack *blc1 = bltrackMan->trackBLC1(i);
    Tools::Fill1D( Form("tracktime_BLC1"), blc1->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BLC1"), blc1->chi2all() );
    if( blc1->CheckRange(blcuts::blc1_time_window_min, blcuts::blc1_time_window_max)){
      nblc1++;
      if( blc1->CheckRange(blcuts::blc1_time_min, blcuts::blc1_time_max) &&
	  bltrackMan->trackBLC1(i)->chi2all()<blcuts::blc1_chi2_max ) blc1GoodTrackID = i;
    }	
  }
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    LocalTrack *blc2 = bltrackMan->trackBLC2(i);
    Tools::Fill1D( Form("tracktime_BLC2"), blc2->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BLC2"), blc2->chi2all() );
    if( blc2->CheckRange(blcuts::blc2_time_window_min, blcuts::blc2_time_window_max) ){
      nblc2++;
      if( blc2->CheckRange(blcuts::blc2_time_min, blcuts::blc2_time_max) &&
	  bltrackMan->trackBLC2(i)->chi2all()<blc2_chi2_max ) blc2GoodTrackID = i;
    }	
  }
  Tools::Fill1D( Form("ntrack_BLC1"), nblc1 );
  Tools::Fill1D( Form("ntrack_BLC2"), nblc2 );

  //** single track selection in each BLC **//
  if( !(nblc1==1 && blc1GoodTrackID!=-1 && nblc2==1 && blc2GoodTrackID!=-1) ){ //** multi-good-beams event is discarded  **//
    Clear( nAbort_singleBLtrack );
    return false;
  }
  
  //** BPC track selection **//
  int nbpc = 0;
  bpcGoodTrackID = -1;
  double chibpc = 999;
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    Tools::Fill1D( Form("tracktime_BPC"), bltrackMan->trackBPC(i)->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BPC"), bltrackMan->trackBPC(i)->chi2all() );
    if( bltrackMan->trackBPC(i)->CheckRange(-30,100) ){
      nbpc++;
      bpcGoodTrackID = i;
      chibpc = bltrackMan->trackBPC(i)->chi2all();
    }
  }

  Tools::Fill1D( Form("ntrack_BPC"), nbpc );
  
  if( nbpc!=1 ){
    Clear( nAbort_nbpc );
    return false;
  }

  LocalTrack *bpctrack = bltrackMan->trackBPC(bpcGoodTrackID);
  if( !(bpctrack->CheckRange(-10,10)) || bpctrack->chi2all()>10 ){ //!! sada-D p.74 !!// //!! Yamaga-D v03 !!//
    Clear( nAbort_bpctrack );
    return false;
  }

  //** vertex calculation **//
  for( int it1=0; it1<trackMan->nGoodTrack(); it1++ ){
    trackMan->CalcVertex_beam( trackMan->GoodTrackID(it1), bltrackMan, confMan );
  }

  //** BLC2-BPC track matching **//
  bool fblc2bpc = false;
  for( int ii=0; ii<bltrackMan->ntrackBLC2(); ii++ ){
    if( ii!=blc2GoodTrackID ) continue;
    LocalTrack *blc2 = bltrackMan->trackBLC2(ii);
    double xblc2bpc[2]={0,0};
    double yblc2bpc[2]={0,0};	 
    double xpos[2]={0,0};
    double ypos[2]={0,0};

    TVector3 Pos_BPC, Pos_BLC2, rot;
    confMan->GetBLDCWireMapManager()->GetGParam( CID_BPC, Pos_BPC, rot );
    confMan->GetBLDCWireMapManager()->GetGParam( CID_BLC2a, Pos_BLC2, rot );
    double zPos_BPC = Pos_BPC.Z();
    double zPos_BLC2 = Pos_BLC2.Z();
    double zPos_BPC_BLC2 = (Pos_BPC.Z()+Pos_BLC2.Z())/2;

    bpctrack->XYPosatZ( zPos_BPC_BLC2, xblc2bpc[0], yblc2bpc[0] );
    bpctrack->XYPosatZ( zPos_BPC, xpos[0], ypos[0] );
    blc2->XYPosatZ( zPos_BPC_BLC2, xblc2bpc[1], yblc2bpc[1]);
    blc2->XYPosatZ( zPos_BLC2, xpos[1], ypos[1]);
    double dxdz[2], dydz[2];
    dxdz[0] = (xpos[0]-xblc2bpc[0]) / (zPos_BPC-zPos_BPC_BLC2);
    dxdz[1] = (xpos[1]-xblc2bpc[1]) / (zPos_BLC2-zPos_BPC_BLC2);
    dydz[0] = (ypos[0]-yblc2bpc[0]) / (zPos_BPC-zPos_BPC_BLC2);
    dydz[1] = (ypos[1]-yblc2bpc[1]) / (zPos_BLC2-zPos_BPC_BLC2);

    if( (xblc2bpc[1]-xblc2bpc[0])< blcuts::blc2bpc_x_min ||
	     (xblc2bpc[1]-xblc2bpc[0])> blcuts::blc2bpc_x_max ) fblc2bpc = false;
    else if( (yblc2bpc[1]-yblc2bpc[0])<blcuts::blc2bpc_y_min ||
	     (yblc2bpc[1]-yblc2bpc[0])>blcuts::blc2bpc_y_max ) fblc2bpc = false;
    else if( (dxdz[1]-dxdz[0])<blcuts::blc2bpc_dx_min ||
	     (dxdz[1]-dxdz[0])>blcuts::blc2bpc_dx_max ) fblc2bpc = false;
    else if( (dydz[1]-dydz[0])<blcuts::blc2bpc_dy_min ||
	     (dydz[1]-dydz[0])>blcuts::blc2bpc_dy_max ) fblc2bpc = false;
    else fblc2bpc = true;
    
    Tools::Fill2D( Form("dydx_BLC2BPC"), xblc2bpc[1]-xblc2bpc[0], yblc2bpc[1]-yblc2bpc[0] );
    Tools::Fill2D( Form("dydzdxdz_BLC2BPC"), dxdz[1]-dxdz[0], dydz[1]-dydz[0] );
  }
  if( !fblc2bpc ){
    Clear( nAbort_fblc2bpc );
    return false;
  }

  return true;
}


//BLC1-D5-BLC2 momentum analysis
double EventAnalysis::AnaBeamSpec(ConfMan *confMan){

  BeamSpectrometer *beamsp = new BeamSpectrometer( confMan );
  LocalTrack *blc1 = bltrackMan->trackBLC1(blc1GoodTrackID);
  LocalTrack *blc2 = bltrackMan->trackBLC2(blc2GoodTrackID);
  beamsp->TMinuitFit( blc1, blc2, confMan );
  double beammom = beamsp->mom();
  double bchi = beamsp->chisquare(); 
  delete beamsp;

  Tools::Fill1D( Form("trackchi2_beam"), bchi );
  if( bchi>blcuts::d5_chi2_max ){
    Clear( nAbort_flagbmom );
    return -9999.;
  }
  Tools::Fill1D( Form("momentum_beam"), beammom );

  return beammom;
}



//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
void EventAnalysis::Clear( int &nAbort)
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  header->Clear();
  cdsMan->Clear();
  blMan->Clear();
  trackMan->Clear();
  bltrackMan->Clear();
  nAbort++;
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
EventTemp *EventAlloc::EventAllocator()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
