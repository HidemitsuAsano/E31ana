//H. Asano
//This code is originated from: EventAnalysis_pipipnn_sakuma.cpp
//and modified to analyze k-d->Lambda pi- proton(missing) -> proton pi- pi- proton(missing)
//----------------------------------------------------------------//
//----------------------------------------------------------------//
//  input : raw data, conf-file, & CDC-tracking-file
//  output: when $(OutFile) is "tmp.root", the following 3 files are generated.
//     "tmp.root":      histogram file
//     "tmp_CDC.root": condensed event file from the CDC-tracking-file with pip/pim selection
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
#include "Tools.h"
#include "ELossTools.h"

#include <TDatabasePDG.h>
#include <Math/ProbFuncMathCore.h>

//kinfitter lib
#include <KinFitter/TKinFitter.h>
#include <KinFitter/TFitParticlePxPyPz.h>
#include <KinFitter/TFitConstraintM.h>
#include <KinFitter/TFitConstraintEp.h>


#include "IMPiSigmaAnaPar.h"
#include "IMPiSigmaHist.h"
#include "IMPiSigmaUtil.h"

#define KFDEBUG 0 // verbose level of the KinFitter
// 0: quiet, 1: print result, 2: print iterations, 3: print also matrices

//-- set run# --//

const unsigned int Verbosity = 0;
const bool DoCDCRetiming = false;
const bool DoKinFit = false;
const bool IsVtxDoubleCheck = false;
//-----------------------------------------//
//--- covariance matrices for KinFitter ---//
//-----------------------------------------//
// ### obtained from (p_meas[j]-p_gene[j])*(p_meas[k]-p_gene[k])
// ###  using G4-data with TH1F(Form("cov_%d_%d_%d", i, j, k), 100, -cov_MAX, cov_MAX);
//   evaluated using "Air" Dora MC
// 1) TLorentzVector LVec_beam, LVec_pim, (LVec_n+LVec_pip), LVec_nmiss, LVec_n, LVec_pip = for pi- Sigma+


class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  void ResetCounters();
  void InitKinFitMatrix();
  TFile *rtFile;// histograms
  TFile *rtFile2;// condensed file
  TFile *rtFile3;// pi+,pi-,n event
  TFile *cdcFile;
  TTree *cdcTree;
  TTree *evTree;
  TTree *ppimpimTree;

  const EventHeader *header_CDC; // original in CDC-tracking-file
  CDSTrackingMan *trackMan_CDC; // original in CDC-tracking-file

  EventHeader *header;
  ScalerMan *scaMan;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  CDSTrackingMan *trackMan;
  BeamLineTrackMan *bltrackMan;
  int t0,t1;

  double scainit[40];
  double scaend[40];
  bool INIT;
  bool SCAOVERFLOW[40];


  int AllGoodTrack;//
  int nTrack;//
  int CDC_Event_Number;//
  int blc1GoodTrackID;//event by event
  int blc2GoodTrackID;//event by event
  int bpcGoodTrackID;// event by event

  //** counters for filling **//
  int nFill_ppimpim;
  //** counters for event abort **//
  int nAbort_KCDH3trg;
  int nAbort_nGoodTrack;
  int nAbort_CDSPID;
  int nAbort_nCDH;
  int nAbort_nT0;
  int nAbort_pid_beam;
  int nAbort_nbpc;
  int nAbort_bpctrack;
  int nAbort_singleBLtrack;
  int nAbort_fblc2bpc;
  int nAbort_flagbmom;
  int nAbort_ftarget;
  int nAbort_pipip;
  int nAbort_end;

  // npippim final-sample tree Branch //

  // 4-momentum(beam) (reaction vtx determined by DCA)
  TLorentzVector mom_beam;
  TLorentzVector mom_target; // 4-momentum(target)
  //ordering of two pi- tracks filled in these vector is determined by chi2 of cds tracking
  TLorentzVector mom_pim1;    // 4-momentum(pi-)
  TLorentzVector mom_pim2;    // 4-momentum(pi-)
  TLorentzVector mom_p;      // 4-momentum(proton)
  TLorentzVector mom_p2;      // 4-momentum(proton)
  TVector3 vtx_reaction; //
  TVector3 vtx_displaced; //
  TVector3 vtx_pim1_beam; //
  TVector3 vtx_pim2_beam; //
  TVector3 vtx_p_beam; //
  TVector3 vtx_pim1_cdc;//
  TVector3 vtx_pim2_cdc;//
  TVector3 vtx_p_cdc;//
  TVector3 CA_pim1_pim1p;//closest approach of pim1-p pim1 side
  TVector3 CA_p_pim1p;//closest approach of pim1-p p side
  TVector3 CA_pim1_pim1p2;//closest approach of pim1-p pim1 side
  TVector3 CA_p2_pim1p2;//closest approach of pim1-p p side
  TVector3 CA_pim2_pim2p;
  TVector3 CA_p_pim2p;
  TVector3 CA_pim2_pim2p2;
  TVector3 CA_p2_pim2p2;
  TVector3 CA_pim1_pim1pim2;
  TVector3 CA_pim2_pim1pim2;
  TVector3 vtx_Lcan_p_pim1;
  TVector3 vtx_Lcan_p_pim2;
  TVector3 vtx_Lcan_p2_pim1;
  TVector3 vtx_Lcan_p2_pim2;
  int ForwardCharge;
  int run_num;   // run number
  int event_num; // event number
  int block_num; // block number
  double ctmT0;

  TMatrixD *covZero;
  TMatrixD *covParticle_Spmode[kin::npart];
  TMatrixD *covParticle_Smmode[kin::npart];
  TLorentzVector kfSpmode_mom_beam;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector kfSpmode_mom_pip;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector kfSpmode_mom_pim;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector kfSpmode_mom_n;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  double kfSpmode_chi2;   // chi2 of kinematical refit
  double kfSpmode_NDF;    // NDF of kinematical refit
  double kfSpmode_status; // status of kinematical refit, 0 :converged 1: not converged
  double kfSpmode_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc
  //= = = = npippim final-sample tree = = = =//
  TDatabasePDG *pdg;

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
  : EventTemp(),
    rtFile(nullptr),
    rtFile2(nullptr),
    rtFile3(nullptr),
    cdcFile(nullptr),
    cdcTree(nullptr),
    evTree(nullptr),
    ppimpimTree(nullptr)
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

  if(DoCDCRetiming) std::cout << " Yes" << std::endl;
  else              std::cout << "  No" << std::endl;

  std::cout << " Kinematic fit ? " ;
  if(DoKinFit) std::cout << " Yes" << std::endl;
  else         std::cout << " No"  << std::endl;

  std::cout << "Double Check VTX fid cut ? " ;
  if(IsVtxDoubleCheck) std::cout << " Yes" << std::endl;
  else         std::cout << " No"  << std::endl;

  std::cout << " CDH TDC cuts " << cdscuts_lpim::tdc_cdh_max << std::endl;
  std::cout << " CDH multiplicity cut: " << cdscuts_lpim::cdhmulti << std::endl;
  std::cout << " CDS # of good tracks cut: " << cdscuts_lpim::cds_ngoodtrack << std::endl;
  std::cout << " use closest pion for vertex " << cdscuts_lpim::useclosestpi << std::endl;
  std::cout << std::endl;
  std::cout << "##################################" << std::endl;
  std::cout << "CDS Neutron ID: beta_MAX " << anacuts::beta_MAX << std::endl;
  std::cout << "CDS Neutron ID: dE_MIN " << anacuts::dE_MIN << std::endl;

  INIT = true;
  //  spillinit = spillfini = -1;
  for( int i=0; i<40; i++ ) {
    SCAOVERFLOW[i] = false;
    scaend[i] = 0;
  }

  //** CDS file open **//
  confMan = conf;
  cdcFile = new TFile( confMan->GetCDSTrackFileName().c_str() );
  if( !cdcFile->IsOpen() ) {
    std::cerr<<" !!! failed to open " <<confMan->GetCDSTrackFileName().c_str()<< "  !!!"<<std::endl;
    exit( false );
  }

  gSystem->Load( "libPhysics.so" );

  //** Getting CDSTracking info. from CDCfile **//
  cdcTree = (TTree*)cdcFile->Get( "EventTree" );
  header_CDC = nullptr;
  trackMan_CDC = nullptr;
  cdcTree->SetBranchAddress( "CDSTrackingMan", &trackMan_CDC );
  cdcTree->SetBranchAddress( "EventHeader",&header_CDC );

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
  if( header==NULL ) {
    std::cerr << "!!!! EventHeader not found"   << std::endl;
    return;
  }
  evTree->Branch( "EventHeader", &header );
  evTree->Branch( "CDSTrackingMan", &trackMan_CDC ); //** = fill original CDSTrackingMan (w/o dE correction) **//

  //** output file 3 : ppimpim final-sample tree **//
  std::string outfile3 = confMan->GetOutFileName();
  outfile3.insert( outfile3.size()-5, "_ppimpim" );
  std::cout<<"npippim file "<<outfile3<<std::endl;
  rtFile3 = new TFile( outfile3.c_str(), "recreate" );
  rtFile3->cd();
  ppimpimTree = new TTree( "EventTree", "EventTree" );
  ppimpimTree->Branch( "mom_beam",   &mom_beam );
  ppimpimTree->Branch( "mom_target", &mom_target );
  ppimpimTree->Branch( "mom_pim1", &mom_pim1 );
  ppimpimTree->Branch( "mom_pim2", &mom_pim2 );
  ppimpimTree->Branch( "mom_p", &mom_p );
  ppimpimTree->Branch( "mom_p2", &mom_p2 );
  ppimpimTree->Branch( "vtx_reaction", &vtx_reaction );
  ppimpimTree->Branch( "vtx_displaced", &vtx_displaced );
  ppimpimTree->Branch( "vtx_pim1_beam", &vtx_pim1_beam );
  ppimpimTree->Branch( "vtx_pim2_beam", &vtx_pim2_beam );
  ppimpimTree->Branch( "vtx_p_beam", &vtx_p_beam );
  ppimpimTree->Branch( "vtx_pim1_cdc", &vtx_pim1_cdc );
  ppimpimTree->Branch( "vtx_pim2_cdc", &vtx_pim2_cdc );
  ppimpimTree->Branch( "vtx_p_cdc", &vtx_p_cdc );
  ppimpimTree->Branch( "CA_pim1_pim1p",&CA_pim1_pim1p);
  ppimpimTree->Branch( "CA_p_pim1p",&CA_p_pim1p);
  ppimpimTree->Branch( "CA_pim1_pim1p2",&CA_pim1_pim1p2);
  ppimpimTree->Branch( "CA_p2_pim1p2",&CA_p2_pim1p2);
  ppimpimTree->Branch( "CA_pim2_pim2p",&CA_pim2_pim2p);
  ppimpimTree->Branch( "CA_p_pim2p",&CA_p_pim2p);
  ppimpimTree->Branch( "CA_pim2_pim2p2",&CA_pim2_pim2p2);
  ppimpimTree->Branch( "CA_p2_pim2p2",&CA_p2_pim2p2);
  ppimpimTree->Branch( "CA_pim1_pim1pim2",&CA_pim1_pim1pim2);
  ppimpimTree->Branch( "CA_pim2_pim1pim2",&CA_pim2_pim1pim2);
  ppimpimTree->Branch( "vtx_Lcan_p_pim1",&vtx_Lcan_p_pim1);
  ppimpimTree->Branch( "vtx_Lcan_p_pim2",&vtx_Lcan_p_pim2);
  ppimpimTree->Branch( "vtx_Lcan_p2_pim1",&vtx_Lcan_p2_pim1);
  ppimpimTree->Branch( "vtx_Lcan_p2_pim2",&vtx_Lcan_p2_pim2);
  ppimpimTree->Branch( "ForwardCharge", &ForwardCharge);  
  //ppimpimTree->Branch( "run_num", &run_num );
  //ppimpimTree->Branch( "event_num", &event_num );
  //ppimpimTree->Branch( "block_num", &block_num );
  ppimpimTree->Branch( "kf_flag", &kf_flag );

  trackMan = new CDSTrackingMan(); //** = dE correction is performed in this code **//
  if( trackMan==NULL ) {
    std::cerr << "!!!! Can not find CDSTrackingMan" << std::endl;
    return;
  }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ) {
    std::cerr << "!!!! Can not find CDSHitMan" << std::endl;
    return;
  }
  blMan = new BeamLineHitMan();
  if( blMan==NULL ) {
    std::cerr << "!!!! Can not find BeamLineHitMan" << std::endl;
    return;
  }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan==NULL ) {
    std::cerr << "!!!! Can not find BeamLineTrackMan" << std::endl;
    return;
  }
  scaMan = new ScalerMan();
  if( scaMan==NULL ) {
    std::cerr << "!!!! Can not find ScalerMan" << std::endl;
    return;
  }


  rtFile->cd();

  t0 = clock();
  ResetCounters();
  InitKinFitMatrix();

  pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");
  //pdg->Print();
}


void EventAnalysis::ResetCounters()
{
  AllGoodTrack = 0;
  nTrack = 0;
  CDC_Event_Number = 0;

  nFill_ppimpim = 0;

  nAbort_nGoodTrack = 0;
  nAbort_CDSPID = 0;
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
  nAbort_pipip = 0;

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
  for( int i=0; i<nsca; i++ ) {
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

#if 0
  for( int i=0; i<11; i++ )
    std::cout<<sca[i]<<"\t";
  std::cout<<std::endl;
#endif

  if( INIT && sca[0]<5 ) {
    for( int i=0; i<scaMan->nsca(); i++ ) {
      scainit[i] = scaMan->sca(i)->val();
      std::cout<<i<<" "<<scainit[i]<<std::endl;
    }
    INIT = false;
  }

  if( INIT ) {
    scaMan->Clear();
    return;
  }

  for( int i=0; i<scaMan->nsca(); i++ ) {
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
  Event_Number++;

  //** control of start and stop events **//
  {
    int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) {
      if( CDC_Event_Number>=cdcTree->GetEntries() ) return false;
      cdcTree->GetEntry( CDC_Event_Number );

      if( header_CDC->ev()<Event_Number ) return true;
      else if( header_CDC->ev()==Event_Number ) {
        CDC_Event_Number++;
        return true;
      }
    }
    if( status==2 ) return false;
  }


  if( Event_Number%5000==1 ) {
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

  //trigger check
  //
  //unbiased kaon trigger,prescaled
  if( header->IsTrig(Trig_Kf)) Tools::Fill1D(Form("Trigger"),0);
  if( header->IsTrig(Trig_KCDH2f)) Tools::Fill1D(Form("Trigger"),2);
  //K x CDH3 trigger
  bool IsTrigKCDH3 = header->IsTrig(Trig_KCDH3);
  if(IsTrigKCDH3) {
    Tools::Fill1D(Form("Trigger"),1);
  } else {
    Clear(nAbort_KCDH3trg);
    return true;
    //return true;
  }

  const int nGoodTrack = trackMan->nGoodTrack();
  const int nallTrack = trackMan->nTrack();
  AllGoodTrack += nGoodTrack;
  nTrack += nallTrack;
  Tools::Fill1D( Form("nTrack"),nallTrack);
  Tools::Fill1D( Form("nGoodTrack"), nGoodTrack );

  Tools::Fill1D( Form("EventCheck"), 1 );

  //CDH-hits cut
  //if( Util::GetCDHMul(cdsMan,nGoodTrack,IsTrigKCDH3,false)!=cdscuts_lpim::cdhmulti){
  if( Util::GetCDHMul(cdsMan,nGoodTrack,IsTrigKCDH3,false)<cdscuts_lpim::cdhmulti ) {
    Clear( nAbort_nCDH );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 2 );

  //** # of good CDS tracks cut **//
  //if( nGoodTrack!=cdscuts_lpim::cds_ngoodtrack ) { //require pi+,pi-,proton
  if( nGoodTrack<cdscuts_lpim::cds_ngoodtrack ) { // dedicated for pi+ pi- event
    Clear( nAbort_nGoodTrack );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 3 );


  //beam line analysis and event selection

  //** BLDC tracking **//
  bltrackMan->DoTracking(blMan,confMan,true,true);
  //Get T0
  int t0seg=-1;
  ctmT0 = Util::AnalyzeT0(blMan,confMan,t0seg);
  if(ctmT0<-9000) {
    Clear( nAbort_nT0 );
    Tools::Fill1D( Form("EventCheck"), 15 );
    return true;
  }

  //PID beam
  const int beamPID = Util::BeamPID(header, ctmT0, blMan);
  if(beamPID!=Beam_Kaon) {
    Clear( nAbort_pid_beam );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 4 );

  const int blstatus = Util::EveSelectBeamline(bltrackMan,trackMan,confMan,blc1GoodTrackID,blc2GoodTrackID,bpcGoodTrackID);

  if(blstatus == -16) {
    Clear( nAbort_singleBLtrack );
    Tools::Fill1D( Form("EventCheck"), 16 );
    return true;
  }

  if(blstatus == -17) {
    Clear( nAbort_nbpc );
    Tools::Fill1D( Form("EventCheck"), 17 );
    return true;
  }

  if(blstatus == -18) {
    Clear( nAbort_bpctrack );
    Tools::Fill1D( Form("EventCheck"), 18 );
    return true;
  }

  if(blstatus == -19) {
    Clear( nAbort_fblc2bpc );
    Tools::Fill1D( Form("EventCheck"), 19 );
    return true;
  }

  Tools::Fill1D( Form("EventCheck"), 5 );


  //BLC1-D5-BLC2 analysis and chi2 selection
  const double beammom = Util::AnaBeamSpec(confMan,bltrackMan,blc1GoodTrackID,blc2GoodTrackID);
  if(beammom <-100) {
    Clear( nAbort_flagbmom );
    return true;//chi2 cut
  }
  Tools::Fill1D( Form("EventCheck"), 6 );

  //** beam momentum calculation **//
  TVector3 Pp_target(0,0,0);
  //Pp_target.SetXYZ( 0, 0, 0 );
  TLorentzVector LVec_beambf;  // 4-Momentum(beam) in LAB
  TLorentzVector LVec_beam;    // 4-Momentum(beam) in LAB with dE correcion
  TLorentzVector LVec_beam_vtx[2];    // 4-Momentum(beam) in LAB with dE correcion
  TLorentzVector LVec_target;  // 4-Momentum(deuteron-target) in LAB
  LVec_target.SetVectM( Pp_target, dMass );//deuteron target
  TLorentzVector LVec_targetP; // 4-Momentum(p-target) in LAB
  //LVec_targetP.SetVectM( Pp_target, pMass );//H target
  TLorentzVector LVec_beambfCM;  // 4-Momentum(beam) in CM
  TLorentzVector LVec_beamCM;    // 4-Momentum(beam) in CM with dE correcion
  TLorentzVector LVec_targetCM;  // 4-Momentum(deuteron-target) in CM
  //TLorentzVector LVec_targetPCM; // 4-Momentum(p-target) in CM

  double x1, y1, x2, y2;
  const double z1 = 0, z2 = 20;
  LocalTrack *bpctrack = bltrackMan->trackBPC(bpcGoodTrackID);
  bpctrack->XYPosatZ( z1, x1, y1 );
  bpctrack->XYPosatZ( z2, x2, y2 );
  TVector3 lp;
  lp.SetXYZ( x1, y1, z1 );
  TVector3 ls;
  ls.SetXYZ( x2-x1, y2-y1, z2-z1);
  ls = ls.Unit();
  const TVector3 Pp_beam = beammom*ls;

  LVec_beambf.SetVectM( Pp_beam, kpMass );
  LVec_beam = LVec_beambf;
  LVec_beam_vtx[0] = LVec_beambf;
  LVec_beam_vtx[1] = LVec_beambf;
  const TVector3 boost = (LVec_target+LVec_beam).BoostVector();
  LVec_beambfCM = LVec_beam;
  LVec_targetCM = LVec_target;
  //LVec_targetPCM = LVec_targetP;
  //boost to CM frame
  LVec_beambfCM.Boost( -1.*boost );
  LVec_targetCM.Boost( -1.*boost );


  //** + + + + + + + + + + + + **//
  //**  PID in CDS             **//
  //** + + + + + + + + + + + + **//
  //** vectors for PID container **//
  std::vector <int> pim_ID;
  std::vector <int> pip_ID;
  std::vector <int> km_ID;
  std::vector <int> p_ID;

  std::vector <int> vCDHseg;
  TVector3 pim_cdhprojected;
  TVector3 pip_cdhprojected;
  // PID of CDS tracks //
  const int nIDedTrack = Util::CDSChargedAna(
                           DoCDCRetiming,
                           bpctrack, cdsMan, trackMan, confMan,blMan,
                           LVec_beam, ctmT0,vCDHseg,pim_ID,pip_ID,km_ID,p_ID,pim_cdhprojected,pip_cdhprojected);
  if(nIDedTrack==-7) Tools::Fill1D( Form("EventCheck"), 7 );//chi2
  if(nIDedTrack==-8) Tools::Fill1D( Form("EventCheck"), 8 );//CDH hit
  if(nIDedTrack==-9) Tools::Fill1D( Form("EventCheck"), 9 );//no CDH sharing
  if(nIDedTrack==-10) Tools::Fill1D( Form("EventCheck"), 10 );//FindMass2_1
  if(nIDedTrack==-11) Tools::Fill1D( Form("EventCheck"), 11 );//FindMass2_2
  if(nIDedTrack==-12) Tools::Fill1D( Form("EventCheck"), 12 );//Energy loss
  if(nIDedTrack<0) {
    Clear(nAbort_CDSPID);
    return true;
  }
  Tools::Fill1D( Form("EventCheck"),13);
  Tools::Fill1D( Form("ntrack_CDS"), nIDedTrack );
  Tools::Fill1D( Form("ntrack_pi_plus"),  pip_ID.size() );
  Tools::Fill1D( Form("ntrack_proton"),   p_ID.size() );
  Tools::Fill1D( Form("ntrack_pi_minus"), pim_ID.size() );
  Tools::Fill1D( Form("ntrack_K_minus"),  km_ID.size() );

  //check forward charge
  //not used for event selection
  //bool forwardcharge = Util::IsForwardCharge(blMan);
  bool forwardcharge=false;

  int nPC =0 ;
  for( int i=0; i<blMan->nPC(); i++ ) {
    if( blMan->PC(i)->CheckRange() ) nPC++;
  }
  if(nPC) forwardcharge=true;

  //  pi+ pi- X event
  //  with CDH multiplicity selection
  //bool ppimpimFlag = false;
  //if( pim_ID.size()==2 && p_ID.size()==1) ppimpimFlag = true;
  if((pim_ID.size()==2) &&
      (
        ((p_ID.size()==1) && (3<=trackMan->nGoodTrack() && trackMan->nGoodTrack()<=4))
        ||((p_ID.size()==2) && (trackMan->nGoodTrack()==4))
      )
    ) {
    //=== pi+ pi- X candidates ===//
    rtFile2->cd();
    evTree->Fill();
    rtFile->cd();
    if(Verbosity) std::cout<<"### filled: Event_Number, Block_Event_Number, CDC_Event_Number = "
                             <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
    nFill_ppimpim++;

    CDSTrack *track_pim1 = trackMan->Track( pim_ID.at(0) ); //
    CDSTrack *track_pim2 = trackMan->Track( pim_ID.at(1) ); //
    CDSTrack *track_p    = trackMan->Track( p_ID.at(0) ); //
    CDSTrack *track_p2;
    if(p_ID.size()==2) {
      track_p2 = trackMan->Track( p_ID.at(1) ); //
    }

    if(Verbosity) {
      std::cout << "pim1 chi2  " << track_pim1->Chi() << std::endl;
      std::cout << "pim1 mom  " << track_pim1->Momentum() << std::endl;
      std::cout << "pim1 angle " << track_pim1->Theta() << std::endl;
      std::cout << "pim2 chi2  " << track_pim2->Chi() << std::endl;
      std::cout << "pim2 mom  " << track_pim2->Momentum() << std::endl;
      std::cout << "pim2 angle " << track_pim2->Theta() << std::endl;
    }

    //deuteron target
    TVector3 vtx_react;//reaction vertex
    TVector3 vtx_dis;//displaced vertex

    TVector3 vtx_beam_wpim1;//vertex(beam-pim1) on beam
    TVector3 vtx_pim1;//vertex(beam-pim) on cdc track
    track_pim1->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wpim1, vtx_pim1 );

    TVector3 vtx_beam_wpim2;//vertex(beam-pim) on beam
    TVector3 vtx_pim2;//vertex(beam-pim) on cdc track
    track_pim2->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wpim2, vtx_pim2 );

    TVector3 vtx_beam_wp;//vertex(beam-proton) on beam
    TVector3 vtx_p;//vertex(beam-proton) on cdc track
    track_p->GetVertex(bpctrack->GetPosatZ(0),bpctrack->GetMomDir(),vtx_beam_wp,vtx_p);

    TVector3 vtx_beam_wp2;//vertex(beam-proton) on beam
    TVector3 vtx_p2;//vertex(beam-proton) on cdc track
    if(p_ID.size()==2) {
      track_p2->GetVertex(bpctrack->GetPosatZ(0),bpctrack->GetMomDir(),vtx_beam_wp2,vtx_p2);
    }

    const double dcapim1vtx =  (vtx_pim1-vtx_beam_wpim1).Mag();
    const double dcapim2vtx =  (vtx_pim2-vtx_beam_wpim2).Mag();
    const TVector3 vtxpim1_mean = 0.5*(vtx_pim1+vtx_beam_wpim1);
    const TVector3 vtxpim2_mean = 0.5*(vtx_pim2+vtx_beam_wpim2);
    Tools::Fill1D( Form("DCA_pim1"), dcapim1vtx );
    Tools::Fill1D( Form("DCA_pim2"), dcapim2vtx );

    TVector3 ca_pim1_pim1p,ca_pim2_pim2p,ca_p_pim1p,ca_p_pim2p,ca_pim1_pim1pim2,ca_pim2_pim1pim2;

    bool vtx_flag1=TrackTools::Calc2HelixVertex(track_pim1, track_p, ca_pim1_pim1p, ca_p_pim1p);
    double dcapim1p=9999.;
    if(vtx_flag1) dcapim1p = (ca_pim1_pim1p-ca_p_pim1p).Mag();
    Tools::Fill1D( Form("DCA_pim1p"), dcapim1p);

    bool vtx_flag2=TrackTools::Calc2HelixVertex(track_pim2, track_p, ca_pim2_pim2p, ca_p_pim2p);
    double dcapim2p=9999.;
    if(vtx_flag2) dcapim2p = (ca_pim2_pim2p-ca_p_pim2p).Mag();
    Tools::Fill1D( Form("DCA_pim2p"), dcapim2p);

    bool vtx_flag3=TrackTools::Calc2HelixVertex(track_pim1, track_pim2, ca_pim1_pim1pim2, ca_pim2_pim1pim2);
    double dcapim1pim2=9999.;
    if(vtx_flag3) dcapim1pim2 = (ca_pim1_pim1pim2-ca_pim2_pim1pim2).Mag();
    Tools::Fill1D( Form("DCA_pim1pim2"), dcapim1pim2);

    TVector3 ca_pim1_pim1p2,ca_p2_pim1p2;
    if(p_ID.size()==2) {
      TrackTools::Calc2HelixVertex(track_pim1, track_p2, ca_pim1_pim1p2, ca_p2_pim1p2);
    }

    TVector3 ca_pim2_pim2p2,ca_p2_pim2p2;
    if(p_ID.size()==2) {
      TrackTools::Calc2HelixVertex(track_pim2, track_p2, ca_pim2_pim2p2, ca_p2_pim2p2);
    }

    //reaction vertex is determined from beam and nearest vtx of pi-
    TVector3 vtx_beam;
    //determine by pim-proton DCA
    if(dcapim1p < dcapim2p) {
      vtx_react = 0.5*(vtx_pim2+vtx_beam_wpim2);
      vtx_dis = ca_p_pim1p;
      vtx_beam = vtx_beam_wpim1;
    } else {
      vtx_react = 0.5*(vtx_pim1+vtx_beam_wpim1);
      vtx_dis = ca_p_pim2p;
      vtx_beam = vtx_beam_wpim2;
    }


    //** beam kaon tof **//
    TVector3 Pos_T0;
    confMan->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
    double beamtof=0;
    double momout=0;
    const double z_pos = Pos_T0.Z();;
    //const double zPos_T0 = Pos_T0.Z();
    //std::cout << "test" << std::endl;
    //std::cout << "bpctrack->GetPosatZ(0)" << std::endl;
    //std::cout << bpctrack->GetPosatZ(0).X() << "  " <<  bpctrack->GetPosatZ(0).Y()  << "  " << bpctrack->GetPosatZ(0).Z()  << std::endl;
    //std::cout << zPos_T0 << std::endl;
    //std::cout << bpctrack->GetPosatZ(zPos_T0 ).X() << "  " <<  bpctrack->GetPosatZ(zPos_T0).Y()  << "  " << bpctrack->GetPosatZ(zPos_T0).Z()  << std::endl;
    //Energy loss correction of beam
    ELossTools::CalcElossBeamTGeo( bpctrack->GetPosatZ(z_pos), vtx_react,
                                   LVec_beambf.Vect().Mag(), kpMass, momout, beamtof );
    LVec_beam.SetVectM( momout*LVec_beambf.Vect().Unit(), kpMass );


    //** reconstructoin of missing proton **//
    TVector3 P_pim1; // Momentum(pi-)
    TVector3 P_pim2; // Momentum(pi-)
    TVector3 P_p; // Momentum(proton)
    TVector3 P_p2; // Momentum(proton)

    TLorentzVector LVec_pim1; // 4-Momentum(pi-)
    TLorentzVector LVec_pim2; // 4-Momentum(pi-)
    TLorentzVector LVec_p;   // 4-Momentum(p) (CDS)
    TLorentzVector LVec_p2;   // 4-Momentum(p) (CDS)
    TLorentzVector LVec_pmiss; // 4-Momentum(n_miss)

    //energy loss correction and momentum correction using vertex info
    if( !track_pim1->GetMomentum( vtx_pim1, P_pim1, true, true ) ) {
      std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
    }

    if( !track_pim2->GetMomentum( vtx_pim2, P_pim2, true, true ) ) {
      std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
    }

    if( !track_p->GetMomentum( vtx_p, P_p, true, true ) ) {
      std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
    }

    if(p_ID.size()==2) {
      if( !track_p2->GetMomentum( vtx_p2, P_p2, true, true ) ) {
        std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }
    }

    LVec_pim1.SetVectM( P_pim1, piMass );
    LVec_pim2.SetVectM( P_pim2, piMass );
    LVec_p.SetVectM(   P_p,   pMass );//
    LVec_p2.SetVectM( P_p2, pMass);

    TVector3 P_ppim1 = P_p + P_pim1;
    TVector3 P_ppim2 = P_p + P_pim2;
    TVector3 P2_p2pim1 = P_p2 + P_pim1;
    TVector3 P2_p2pim2 = P_p2 + P_pim2;

    double dl=0;
    double dist=0;
    TVector3 nest;
    nest.SetXYZ(0,0,0);
    TVector3 CA_ppim1_center = 0.5*(ca_p_pim1p+ca_pim1_pim1p);
    MathTools::LineToLine(CA_ppim1_center,P_ppim1.Unit(),bpctrack->GetPosatZ(0), bpctrack->GetMomDir(),dl,dist,vtx_Lcan_p_pim1,nest);

    dl=0;
    dist=0;
    nest.SetXYZ(0,0,0);
    TVector3 CA_ppim2_center = 0.5*(ca_p_pim2p+ca_pim2_pim2p);
    MathTools::LineToLine(CA_ppim2_center,P_ppim2.Unit(),bpctrack->GetPosatZ(0), bpctrack->GetMomDir(),dl,dist,vtx_Lcan_p_pim2,nest);

    if(p_ID.size()==2) {
      dl=0;
      dist=0;
      nest.SetXYZ(0,0,0);
      TVector3 CA_p2pim1_center = 0.5*(ca_p2_pim1p2+ca_pim1_pim1p2);
      MathTools::LineToLine(CA_p2pim1_center,P2_p2pim1.Unit(),bpctrack->GetPosatZ(0), bpctrack->GetMomDir(),dl,dist,vtx_Lcan_p2_pim1,nest);

      dl=0;
      dist=0;
      nest.SetXYZ(0,0,0);
      TVector3 CA_p2pim2_center = 0.5*(ca_p2_pim2p2+ca_pim2_pim2p2);
      MathTools::LineToLine(CA_p2pim2_center,P2_p2pim2.Unit(),bpctrack->GetPosatZ(0), bpctrack->GetMomDir(),dl,dist,vtx_Lcan_p2_pim2,nest);
    }







    //const double pimphi = P_pim.Phi();
    //const double pipphi = P_pip.Phi();
    //Tools::Fill1D(Form("npimangle"),pimphi-cdhphi);
    //Tools::Fill1D(Form("npipangle"),pipphi-cdhphi);

    const double mm_mass   = (LVec_target+LVec_beam-LVec_pim1-LVec_pim2-LVec_p).M();
    //double mm_mass_vtx[2];
    //for(int ivtx=0;ivtx<2;ivtx++){
    //mm_mass_vtx[ivtx]= (LVec_target+LVec_beam_vtx[ivtx]-LVec_pim-LVec_pip-LVec_n_vtx[ivtx]).M();
    //}

    const TVector3 P_missp = (LVec_target+LVec_beam-LVec_pim1-LVec_pim2-LVec_p).Vect();
    //TVector3 P_missn_vtx[2];
    //for(int ivtx=0;ivtx<2;ivtx++){
    //  P_missn_vtx[ivtx] = (LVec_target+LVec_beam_vtx[ivtx]-LVec_pim-LVec_pip-LVec_n_vtx[ivtx]).Vect();
    //}

    LVec_pmiss.SetVectM( P_missp, pMass );
    //for(int ivtx=0;ivtx<2;ivtx++){
    //  LVec_nmiss_vtx[ivtx].SetVectM( P_missn_vtx[ivtx], nMass );
    //}

    if(Verbosity>10)std::cerr<<"  missing mass = "<<mm_mass<<std::endl;

    const TVector3 boost = (LVec_target+LVec_beam).BoostVector();
    TLorentzVector LVec_pmiss_CM = LVec_pmiss;
    TLorentzVector LVec_beam_CM = LVec_beam;
    LVec_pmiss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    //cos in CM frame
    const double cos_p = LVec_pmiss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pmiss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    if(Verbosity>10)std::cerr<<"  missing mom | cos_CM = "<<cos_p<<std::endl;


    //** + + + + + + + + + + + + + **//
    //**  fill histograms & tree   **//
    //** + + + + + + + + + + + + + **//
    kf_flag = -1;
    bool MissPFlag=false;
    bool LambdaFlag=false;

    Tools::Fill2D(Form("MMom_MMass"), mm_mass, P_missp.Mag() );
    Tools::Fill2D(Form("Vtx_ZX_nofid"),vtxpim1_mean.Z(),vtxpim1_mean.X());
    Tools::Fill2D(Form("Vtx_ZY_nofid"),vtxpim1_mean.Z(),vtxpim1_mean.Y());
    Tools::Fill2D(Form("Vtx_XY_nofid"),vtxpim1_mean.X(),vtxpim1_mean.Y());
    Tools::Fill2D(Form("Vtx_ZX_nofid"),vtxpim2_mean.Z(),vtxpim2_mean.X());
    Tools::Fill2D(Form("Vtx_ZY_nofid"),vtxpim2_mean.Z(),vtxpim2_mean.Y());
    Tools::Fill2D(Form("Vtx_XY_nofid"),vtxpim2_mean.X(),vtxpim2_mean.Y());
    //Fiducial cuts OK
    if( (!IsVtxDoubleCheck && (GeomTools::GetID(vtx_react)==CID_Fiducial)) ||
        ( IsVtxDoubleCheck &&
          (GeomTools::GetID(vtxpim1_mean)==CID_Fiducial) &&
          (GeomTools::GetID(vtxpim2_mean)==CID_Fiducial)))  {

      Tools::Fill2D(Form("Vtx_ZX_primfid"),vtx_react.Z(),vtx_react.X());
      Tools::Fill2D(Form("Vtx_ZY_primfid"),vtx_react.Z(),vtx_react.Y());
      Tools::Fill2D(Form("Vtx_XY_primfid"),vtx_react.X(),vtx_react.Y());
      Tools::Fill2D(Form("Vtx_ZX_fid"),vtxpim1_mean.Z(),vtxpim1_mean.X());
      Tools::Fill2D(Form("Vtx_ZY_fid"),vtxpim1_mean.Z(),vtxpim1_mean.Y());
      Tools::Fill2D(Form("Vtx_XY_fid"),vtxpim1_mean.X(),vtxpim1_mean.Y());
      Tools::Fill2D(Form("Vtx_ZX_fid"),vtxpim2_mean.Z(),vtxpim2_mean.X());
      Tools::Fill2D(Form("Vtx_ZY_fid"),vtxpim2_mean.Z(),vtxpim2_mean.Y());
      Tools::Fill2D(Form("Vtx_XY_fid"),vtxpim2_mean.X(),vtxpim2_mean.Y());


      Tools::Fill2D( Form("MMom_MMass_fid"), mm_mass, P_missp.Mag() );

      //missing mass proton ID
      if( anacuts_lpim::proton_MIN<mm_mass && mm_mass<anacuts_lpim::proton_MAX ) MissPFlag=true;

      //Lambda production in CDS
      if( (anacuts_lpim::ppi_MIN<(LVec_p+LVec_pim1).M() && (LVec_p+LVec_pim1).M()<anacuts_lpim::ppi_MAX)) LambdaFlag=true;
      if( (anacuts_lpim::ppi_MIN<(LVec_p+LVec_pim2).M() && (LVec_p+LVec_pim2).M()<anacuts_lpim::ppi_MAX)) LambdaFlag=true;

      if(MissPFlag) {
        Tools::Fill2D( Form("MMom_MMass_fid_p"), mm_mass, P_missp.Mag() );
        Tools::Fill2D( Form("MMom_PMom_fid_p"), P_p.Mag(), P_missp.Mag() );
        Tools::Fill2D( Form("IMppim1_IMppim2_p"), (LVec_p+LVec_pim1).M(), (LVec_p+LVec_pim2).M() );
        Tools::Fill2D( Form("q_IMppipi_p"), (LVec_p+LVec_pim1+LVec_pim2).M(), (LVec_beam.Vect()-LVec_pmiss.Vect()).Mag());
      }

      if(LambdaFlag) {
        Tools::Fill2D( Form("MMom_MMass_fid_wL"), mm_mass, P_missp.Mag() );
      }

      if(MissPFlag && LambdaFlag) {
        Tools::Fill2D( Form("q_IMppipi_wL_p"), (LVec_p+LVec_pim1+LVec_pim2).M(), (LVec_beam.Vect()-LVec_pmiss.Vect()).Mag());
      }

      /*
      if(DoKinFit) {
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
        // %%% Kinematical Fit using KinFitter %%% //
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
        //--- set TLorentzVector ---//
        // beam_K(K+), pi-/+, Sigma+/-,  missn, n from S, pi+/- from S
        //  = 1) TLorentzVector LVec_beam, LVec_pim, (LVec_n+LVec_pip), LVec_nmiss, LVec_n, LVec_pip = for pi- Sigma+
        TLorentzVector TL_meas_Spmode[kin::npart]; // measured
        //TL_meas_Spmode[kin::kmbeam] = LVec_beam;
        TL_meas_Spmode[kin::kmbeam] = LVec_beam_vtx[0];
        TL_meas_Spmode[kin::pim_g1] = LVec_pim;
        //TL_meas_Spmode[kin::Sp] = (LVec_n+LVec_pip);
        TL_meas_Spmode[kin::Sp] = (LVec_n_vtx[0]+LVec_pip);
        //TL_meas_Spmode[kin::nmiss] = LVec_nmiss;
        TL_meas_Spmode[kin::nmiss] = LVec_nmiss_vtx[0];
        //TL_meas_Spmode[kin::ncds] = LVec_n;
        TL_meas_Spmode[kin::ncds] = LVec_n_vtx[0];
        TL_meas_Spmode[kin::pip_g2] = LVec_pip;
        //  = 2) TLorentzVector LVec_beam, LVec_pip, (LVec_n+LVec_pim), LVec_nmiss, LVec_n, LVec_pim = for pi+ Sigma-
        TLorentzVector TL_meas_Smmode[kin::npart]; // measured
        //TL_meas_Smmode[kin::kmbeam] = LVec_beam;
        TL_meas_Smmode[kin::kmbeam] = LVec_beam_vtx[1];
        TL_meas_Smmode[kin::pip_g1] = LVec_pip;
        //TL_meas_Smmode[kin::Sm] = (LVec_n+LVec_pim);
        TL_meas_Smmode[kin::Sm] = (LVec_n_vtx[1]+LVec_pim);
        //TL_meas_Smmode[kin::nmiss] = LVec_nmiss;
        TL_meas_Smmode[kin::nmiss] = LVec_nmiss_vtx[1];
        //TL_meas_Smmode[kin::ncds] = LVec_n;
        TL_meas_Smmode[kin::ncds] = LVec_n_vtx[1];
        TL_meas_Smmode[kin::pim_g2] = LVec_pim;
        TLorentzVector TL_kfit_Spmode[kin::npart]; // kinematical fitted
        TLorentzVector TL_kfit_Smmode[kin::npart]; // kinematical fitted
        // LVec_target is defined as (0, 0, 0, M_d2)
        TVector3 TV_target = LVec_target.Vect();
        TVector3 TV_meas_Spmode[kin::npart];
        TVector3 TV_meas_Smmode[kin::npart];
        for( int i=0; i<kin::npart; i++ ) {
          TV_meas_Spmode[i] = TL_meas_Spmode[i].Vect();
          TV_meas_Smmode[i] = TL_meas_Smmode[i].Vect();
        }

        //These partcile IDs are defined in pythia6
        //see http://home.fnal.gov/~mrenna/lutp0613man2/node44.html
        //                    K-    pi-  S+     n     n     pi+
        const int PDG_Spmode[kin::npart] = {-321, -211, 3222, 2112, 2112,  211}; // pi-Sigma+
        //                    K-    pi+  S-     n     n     pi-
        const int PDG_Smmode[kin::npart] = {-321,  211, 3112, 2112, 2112, -211}; // pi+Sigma-


        //--- KinFitter :: initialization ---//
        //  = 1) TLorentzVector LVec_beam, LVec_pim, (LVec_n+LVec_pip), LVec_nmiss, LVec_n, LVec_pip = for pi- Sigma+
        //  = 2) TLorentzVector LVec_beam, LVec_pip, (LVec_n+LVec_pim), LVec_nmiss, LVec_n, LVec_pim = for pi+ Sigma-
        //definition of fit particles in cartesian coordinates //
        const TString str_particle_Spmode[kin::npart] = {"LVec_beam", "LVec_pim", "LVec_Sp", "LVec_mn", "LVec_n", "LVec_pip"};
        const TString str_particle_Smmode[kin::npart] = {"LVec_beam", "LVec_pip", "LVec_Sm", "LVec_mn", "LVec_n", "LVec_pim"};
        //asano memo
        //TFitParticlePxPyPz this is KinFitter class
        TFitParticlePxPyPz ParticleTgt = TFitParticlePxPyPz("target", "target", &TV_target,
                                         pdg->GetParticle("deuteron")->Mass(), covZero);
        TFitParticlePxPyPz Particle_Spmode[kin::npart];//-321,-211,3222,2112,2112, 211
        TFitParticlePxPyPz Particle_Smmode[kin::npart];//-321, 211,3112,2112,2112,-211
        for( int i=0; i<kin::npart; i++ ) {
          Particle_Spmode[i] = TFitParticlePxPyPz(str_particle_Spmode[i], str_particle_Spmode[i], &TV_meas_Spmode[i],
                                                  pdg->GetParticle(PDG_Spmode[i])->Mass(), covParticle_Spmode[i]);
          Particle_Smmode[i] = TFitParticlePxPyPz(str_particle_Smmode[i], str_particle_Smmode[i], &TV_meas_Smmode[i],
                                                  pdg->GetParticle(PDG_Smmode[i])->Mass(), covParticle_Smmode[i]);
        }//for i
        // definition of constraints //
        // constraint :: mass of Sigma
        TFitConstraintM ConstMS_Spmode = TFitConstraintM("M_Sp", "M_Sp", 0, 0, pdg->GetParticle(PDG_Spmode[kin::Sp])->Mass());
        TFitConstraintM ConstMS_Smmode = TFitConstraintM("M_Sm", "M_Sm", 0, 0, pdg->GetParticle(PDG_Smmode[kin::Sm])->Mass());
        ConstMS_Spmode.addParticles1(&Particle_Spmode[kin::ncds], &Particle_Spmode[kin::pip_g2]);
        ConstMS_Smmode.addParticles1(&Particle_Smmode[kin::ncds], &Particle_Smmode[kin::pim_g2]);
        // constraint :: 4-momentum conservation
        TFitConstraintEp ConstEp_Spmode[4];
        TFitConstraintEp ConstEp_Smmode[4];
        const TString str_constEp_Spmode[4]  = {"Px", "Py", "Pz", "E"};
        const TString str_constEp_Smmode[4]  = {"Px", "Py", "Pz", "E"};
        for( int i=0; i<4; i++ ) {
          ConstEp_Spmode[i] = TFitConstraintEp(str_constEp_Spmode[i], str_constEp_Spmode[i], 0, TFitConstraintEp::component(i), 0);
          ConstEp_Smmode[i] = TFitConstraintEp(str_constEp_Smmode[i], str_constEp_Smmode[i], 0, TFitConstraintEp::component(i), 0);
          ConstEp_Spmode[i].addParticles1(&ParticleTgt, &Particle_Spmode[kin::kmbeam]);
          ConstEp_Smmode[i].addParticles1(&ParticleTgt, &Particle_Smmode[kin::kmbeam]);
          ConstEp_Spmode[i].addParticles2(&Particle_Spmode[kin::pim_g1], &Particle_Spmode[kin::nmiss], &Particle_Spmode[kin::ncds], &Particle_Spmode[kin::pip_g2]);// pim, miss_n,cds_n, pip

          ConstEp_Smmode[i].addParticles2(&Particle_Smmode[kin::pip_g1], &Particle_Smmode[kin::nmiss], &Particle_Smmode[kin::ncds], &Particle_Smmode[kin::pim_g2]);// pip, miss_n,cds_n, pim
        }//for

        //--- KinFitter :: execution ---//
        // definition of the fitter //
        TKinFitter kinfitter_Spmode;
        TKinFitter kinfitter_Smmode;
        // add measured particles
        kinfitter_Spmode.addMeasParticles(&Particle_Spmode[kin::kmbeam], &Particle_Spmode[kin::pim_g1], &Particle_Spmode[kin::ncds], &Particle_Spmode[kin::pim_g2]); // K, pi-, n, pi+
        kinfitter_Smmode.addMeasParticles(&Particle_Smmode[kin::kmbeam], &Particle_Smmode[kin::pip_g1], &Particle_Smmode[kin::ncds], &Particle_Smmode[kin::pip_g2]); // K, pi+, n, pi-
        kinfitter_Spmode.addUnmeasParticles(&Particle_Spmode[kin::nmiss]); // missing-n
        kinfitter_Smmode.addUnmeasParticles(&Particle_Smmode[kin::nmiss]); // missing-n
        // add constraints
        kinfitter_Spmode.addConstraint(&ConstMS_Spmode); // mass of Sigma+
        kinfitter_Smmode.addConstraint(&ConstMS_Smmode); // mass of Sigma-
        for( int i=0; i<4; i++ ) {
          kinfitter_Spmode.addConstraint(&ConstEp_Spmode[i]); // 4-momentum conservation
          kinfitter_Smmode.addConstraint(&ConstEp_Smmode[i]); // 4-momentum conservation
        }

        // perform the fit //
        kinfitter_Spmode.setMaxNbIter(kin::maxitr);       // max number of iterations
        kinfitter_Smmode.setMaxNbIter(kin::maxitr);       // max number of iterations
        kinfitter_Spmode.setMaxDeltaS(kin::maxdchi2);     // max delta chi2
        kinfitter_Smmode.setMaxDeltaS(kin::maxdchi2);     // max delta chi2
        kinfitter_Spmode.setMaxF(kin::maxsumconst);          // max sum of constraints
        kinfitter_Smmode.setMaxF(kin::maxsumconst);          // max sum of constraints
        kinfitter_Spmode.setVerbosity(KFDEBUG);  // verbosity level
        kinfitter_Smmode.setVerbosity(KFDEBUG);  // verbosity level
        kinfitter_Spmode.fit();
        kinfitter_Smmode.fit();
        // copy fit results //
        for( int i=0; i<kin::npart; i++ ) {
          TL_kfit_Spmode[i] = (*Particle_Spmode[i].getCurr4Vec());
          TL_kfit_Smmode[i] = (*Particle_Smmode[i].getCurr4Vec());
        }
        TL_kfit_Spmode[kin::Sp] = TL_kfit_Spmode[kin::ncds]+TL_kfit_Spmode[kin::pip_g2];
        TL_kfit_Smmode[kin::Sm] = TL_kfit_Smmode[kin::ncds]+TL_kfit_Smmode[kin::pim_g2];


        Tools::Fill2D( Form("KFchi2_vs"), kinfitter_Spmode.getS()/kinfitter_Spmode.getNDF(),
                       kinfitter_Smmode.getS()/kinfitter_Smmode.getNDF() );

        if(Verbosity) {
          std::cerr<<"pi- S+ : status = "<<kinfitter_Spmode.getStatus()<<", chi2/NDF = "<<kinfitter_Spmode.getS()<<"/"<<kinfitter_Spmode.getNDF()<<std::endl;
          std::cerr<<"pi+ S- : status = "<<kinfitter_Smmode.getStatus()<<", chi2/NDF = "<<kinfitter_Smmode.getS()<<"/"<<kinfitter_Smmode.getNDF()<<std::endl;
        }

        // fill tree //
        kfSpmode_mom_beam   = TL_kfit_Spmode[kin::kmbeam];
        kfSpmode_mom_pip    = TL_kfit_Spmode[kin::pip_g2];
        kfSpmode_mom_pim    = TL_kfit_Spmode[kin::pim_g1];
        kfSpmode_mom_n      = TL_kfit_Spmode[kin::ncds];
        kfSpmode_chi2      = kinfitter_Spmode.getS();
        kfSpmode_NDF       = kinfitter_Spmode.getNDF();
        kfSpmode_status    = kinfitter_Spmode.getStatus();
        kfSpmode_pvalue    = ROOT::Math::chisquared_cdf_c(kinfitter_Spmode.getS(), kinfitter_Spmode.getNDF());
        kfSmmode_mom_beam   = TL_kfit_Smmode[kin::kmbeam];
        kfSmmode_mom_pip    = TL_kfit_Smmode[kin::pip_g1];
        kfSmmode_mom_pim    = TL_kfit_Smmode[kin::pim_g2];
        kfSmmode_mom_n      = TL_kfit_Smmode[kin::ncds];
        kfSmmode_chi2      = kinfitter_Smmode.getS();
        kfSmmode_NDF       = kinfitter_Smmode.getNDF();
        kfSmmode_status    = kinfitter_Smmode.getStatus();
        kfSmmode_pvalue    = ROOT::Math::chisquared_cdf_c(kinfitter_Smmode.getS(), kinfitter_Smmode.getNDF());
        kf_flag       = 1;

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
        // %%% Kinematical Fit using KinFitter %%% //
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
      }//if DoKinFit
      */

      mom_beam   = LVec_beam;   // 4-momentum(beam)
      mom_target = LVec_target; // 4-momentum(target)
      mom_pim1 = LVec_pim1;        // 4-momentum(pi+)
      mom_pim2 = LVec_pim2;        // 4-momentum(pi-)
      mom_p = LVec_p;            // 4-momentum(neutron)
      mom_p2 = LVec_p2;            // 4-momentum(neutron)
      vtx_reaction = vtx_react; // vertex(reaction)
      vtx_displaced = vtx_dis; // vertex(reaction)
      vtx_pim1_beam = vtx_beam_wpim1;
      vtx_pim2_beam = vtx_beam_wpim2;
      vtx_p_beam = vtx_beam_wp;
      vtx_pim1_cdc = vtx_pim1;
      vtx_pim2_cdc = vtx_pim2;
      vtx_p_cdc = vtx_p;
      CA_pim1_pim1p = ca_pim1_pim1p;
      CA_p_pim1p = ca_p_pim1p;
      CA_pim1_pim1p2 = ca_pim1_pim1p2;
      CA_p2_pim1p2 = ca_p2_pim1p2;
      CA_pim2_pim2p = ca_pim2_pim2p;
      CA_p_pim2p = ca_p_pim2p;
      CA_pim2_pim2p2 = ca_pim2_pim2p2;
      CA_p2_pim2p2 = ca_p2_pim2p2;
      CA_pim1_pim1pim2 = ca_pim1_pim1pim2;
      CA_pim2_pim1pim2 = ca_pim2_pim1pim2;
      ForwardCharge = forwardcharge;
      //run_num   = confMan->GetRunNumber(); // run number
      //event_num = Event_Number;            // event number
      //block_num = Block_Event_Number;      // block number

      if(Verbosity>10)std::cout<<"End loop :: %%% npippim event: Event_Number, Block_Event_Number, CDC_Event_Number = "
                                 <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
      rtFile3->cd();
      ppimpimTree->Fill();
      rtFile->cd();
      //** fill tree **//

    } // if( GeomTools::GetID(vtx_react)==CID_Fiducial )
  }//pi-,pi- p event
  else {
    Clear( nAbort_pipip );
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
  for(int i=0; i<40; i++) {
    TH1F* h1 = (TH1F*)gFile->Get("Scaler");
    std::cout<<i<<"  "<<scaend[i]<<std::endl;
    h1->Fill(i,scaend[i]-scainit[i]);
  }

  std::cout<<"====== Abort counter ========="<<std::endl;
  std::cout<<" nAbort_nGoodTrack    = "<<nAbort_nGoodTrack<<std::endl;
  std::cout<<" nAbort_CDSPID        = "<<nAbort_CDSPID<<std::endl;
  std::cout<<" nAbort_nCDH          = "<<nAbort_nCDH<<std::endl;
  std::cout<<" nAbort_nT0           = "<<nAbort_nT0<<std::endl;
  std::cout<<" nAbort_pid_beam      = "<<nAbort_pid_beam<<std::endl;
  std::cout<<" nAbort_nbpc          = "<<nAbort_nbpc<<std::endl;
  std::cout<<" nAbort_bpctrack      = "<<nAbort_bpctrack<<std::endl;
  std::cout<<" nAbort_singleBLtrack = "<<nAbort_singleBLtrack<<std::endl;
  std::cout<<" nAbort_fblc2bpc      = "<<nAbort_fblc2bpc<<std::endl;
  std::cout<<" nAbort_flagbmom      = "<<nAbort_flagbmom<<std::endl;
  std::cout<<" nAbort_ftarget       = "<<nAbort_ftarget<<std::endl;
  std::cout<<" nAbort_nAbort_pipip   = "<<nAbort_pipip<<std::endl;
  std::cout<<" nAbort_end           = "<<nAbort_end<<std::endl;
  std::cout<<"========= Abort counter ========="<<std::endl;
  std::cout<<"*** # of pi- pi- p events = "<<nFill_ppimpim<<" ***"<<std::endl;

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
  Tools::newTH1F( Form("EventCheck"), 20, 0, 20 );
  // Event Reduction Check memo
  // 1.  # of all events w/o cosmic trigger
  // 2.  after CDH cuts
  // 3.  after the # of CDS good tracks cut
  // 4.  after beam PID (kaon)
  // 5.  after beam line analysis
  // 6.  after BLC1-D5-BLC2 analysis
  // 7.  filled if CDS chi2 is bad
  // 8.  filled if a CDC track have no CDH hits
  // 9.  filled if CDC tracks share a CDH segments
  // 10. filled if FindMass2 is failed
  // 11. filled if FindMass2 is failed after Retiming
  // 12. filled if Energy loss calculation failed
  // 13. filled if CDSChargedAna is OK.
  // 14. filled if CDH isolation is OK.
  // 15. T0 multiplicity is NOT 1
  // 16. # of BLC1-BLC2 tracks is NOT 1
  // 17. # of BPC tracks is NOT 1
  // 18. BPC timing and track chi2 is NOT good
  // 19. BLC2-BPC track is not matched
  InitBasicHist();
  InitIMLambdaPimHist();


  return;
}




void EventAnalysis::InitKinFitMatrix()
{
  //-----------------------------------------//
  //--- covariance matrices for KinFitter ---//
  //-----------------------------------------//
  covZero = new TMatrixD(4, 4);
  covZero->Zero();
  covZero->ResizeTo(3, 3); // resize from 4x4 to 3x3
  for( int i=0; i<kin::npart; i++ ) {
    covParticle_Spmode[i] = new TMatrixD(4, 4);
    covParticle_Smmode[i] = new TMatrixD(4, 4);
    int n = 0;
    for( int j=0; j<4; j++ ) {
      for( int k=0; k<4; k++ ) {
        if( j==k ) {
          (*covParticle_Spmode[i])[j][k] = kin::covValSpmode[i][n]; // only diagonal elements
          (*covParticle_Smmode[i])[j][k] = kin::covValSmmode[i][n]; // only diagonal elements
        } else {
          (*covParticle_Spmode[i])[j][k] = 0;
          (*covParticle_Smmode[i])[j][k] = 0;
        }
        n++;
      }
    }
    covParticle_Spmode[i]->ResizeTo(3, 3); // resize from 4x4 to 3x3
    covParticle_Smmode[i]->ResizeTo(3, 3); // resize from 4x4 to 3x3
    covParticle_Spmode[i]->Print(); // Print all
    covParticle_Smmode[i]->Print(); // Print all
  }
  //-----------------------------------------//
  //--- covariance matrices for KinFitter ---//
  //-----------------------------------------//

  return;
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
  //beam line
  blc1GoodTrackID=-1;
  blc2GoodTrackID=-1;
  bpcGoodTrackID=-1;
  ctmT0 = 0;
  ForwardCharge = 0;

  //CDS
  vtx_reaction.SetXYZ(-9999.,-9999.,-9999.);
  vtx_displaced.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pim1_beam.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pim2_beam.SetXYZ(-9999.,-9999.,-9999.);
  vtx_p_beam.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pim1_cdc.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pim2_cdc.SetXYZ(-9999.,-9999.,-9999.);
  vtx_p_cdc.SetXYZ(-9999.,-9999.,-9999.);
  CA_pim1_pim1p.SetXYZ(-9999.,-9999.,-9999.);
  CA_p_pim1p.SetXYZ(-9999.,-9999.,-9999.);
  CA_pim1_pim1p2.SetXYZ(-9999.,-9999.,-9999.);
  CA_p2_pim1p2.SetXYZ(-9999.,-9999.,-9999.);
  CA_pim2_pim2p.SetXYZ(-9999.,-9999.,-9999.);
  CA_p_pim2p.SetXYZ(-9999.,-9999.,-9999.);
  CA_pim2_pim2p2.SetXYZ(-9999.,-9999.,-9999.);
  CA_p2_pim2p2.SetXYZ(-9999.,-9999.,-9999.);
  CA_pim1_pim1pim2.SetXYZ(-9999.,-9999.,-9999.);
  CA_pim2_pim1pim2.SetXYZ(-9999.,-9999.,-9999.);
  vtx_Lcan_p_pim1.SetXYZ(-9999.,-9999.,-9999.);
  vtx_Lcan_p_pim2.SetXYZ(-9999.,-9999.,-9999.);
  vtx_Lcan_p2_pim1.SetXYZ(-9999.,-9999.,-9999.);
  vtx_Lcan_p2_pim2.SetXYZ(-9999.,-9999.,-9999.);
  
  return;
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
EventTemp *EventAlloc::EventAllocator()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
