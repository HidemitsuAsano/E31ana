//H. Asano
//This code is originated from: EventAnalysis_pipipnn_sakuma.cpp
//and modified to analyze k-d->npi+Sigma-, npi-,Sigma+ events
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

const int MaxTreeSize = 1000000000;

const unsigned int Verbosity = 0;
const bool DoCDCRetiming = false;
const bool DoKinFit = true;
const bool IsVtxDoubleCheck = false;
const bool UseDecayVtx = true;
const unsigned int IsolationCutFlag = 1;
//-----------------------------------------//
//--- covariance matrices for KinFitter ---//
//-----------------------------------------//
// ### obtained from (p_meas[j]-p_gene[j])*(p_meas[k]-p_gene[k])
// ###  using G4-data with TH1F(Form("cov_%d_%d_%d", i, j, k), 100, -cov_MAX, cov_MAX);
//   evaluated using "Air" Dora MC
// 1) TLorentzVector LVec_beam, LVec_pim, (LVec_n+LVec_pip), LVec_pmiss, LVec_n, LVec_pip = for pi- Sigma+


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
  int nFill_pippim;
  int nFill_npippim;
  //** counters for event abort **//
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
  int nAbort_CDHiso;
  int nAbort_pipi;
  int nAbort_end;

  // npippim final-sample tree Branch //
  
  // 4-momentum(beam) (reaction vtx determined by DCA)
  TLorentzVector mom_beam;   
  // 4-momentum(beam) (reaction vtx Sp mode assumption)
  TLorentzVector mom_beam_Sp;  
  // 4-momentum(beam) (reaction vtx Sm mode assumption)
  TLorentzVector mom_beam_Sm;
  TLorentzVector mom_target; // 4-momentum(target)
  TLorentzVector mom_pip;    // 4-momentum(pi+)
  TLorentzVector mom_pim;    // 4-momentum(pi-)
  TLorentzVector mom_n;      // 4-momentum(neutron)
  TLorentzVector mom_n_Sp;  // 4-momentum(neutron),Sp mode assumption
  TLorentzVector mom_n_Sm;  // 4-momentum(neutron),Sm mode assumption
  double NeutralBetaCDH; // velocity of neutral particle on CDH from decay vtx 
  double NeutralBetaCDH_vtx[2];//1:pip_vtx,2:pim_vtx  
  double dE;   // energy deposit on CDH [MeVee]
  TVector3 vtx_reaction; // 
  TVector3 vtx_pip_beam; // 
  TVector3 vtx_pim_beam; // 
  TVector3 vtx_pip_cdc;//
  TVector3 vtx_pim_cdc;//
  TVector3 CA_pip;
  TVector3 CA_pim;
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
  TLorentzVector kfSmmode_mom_beam;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
  TLorentzVector kfSmmode_mom_pip;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
  TLorentzVector kfSmmode_mom_pim;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
  TLorentzVector kfSmmode_mom_n;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
  double kfSmmode_chi2;   // chi2 of kinematical refit
  double kfSmmode_NDF;    // NDF of kinematical refit
  double kfSmmode_status; // status of kinematical refit, 0 : convergetd ,1:not converged
  double kfSmmode_pvalue; // p-value of kinematical refit
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
  npippimTree(nullptr)
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
  
  std::cout << "Use Decay VTX for neutron ? " ;
  if(UseDecayVtx) std::cout << " Yes" << std::endl;
  else         std::cout << " No"  << std::endl;
  
  std::cout << "Isolation cut range ? " ;
  std::cout << IsolationCutFlag << "  segments" << std::endl;


  std::cout << " CDH TDC cuts " << cdscuts::tdc_cdh_max << std::endl;
  std::cout << " CDH multiplicity cut: " << cdscuts::cdhmulti << std::endl;
  std::cout << " CDS # of good tracks cut: " << cdscuts::cds_ngoodtrack << std::endl;
  std::cout << " use closest pion for vertex " << cdscuts::useclosestpi << std::endl;
  std::cout << std::endl;
  std::cout << "##################################" << std::endl;
  std::cout << "CDS Neutron ID: beta_MAX " << anacuts::beta_MAX << std::endl;
  std::cout << "CDS Neutron ID: dE_MIN " << anacuts::dE_MIN << std::endl;
  std::cout << "K0 rejection window "
            <<  anacuts::pipi_MIN << " - " << anacuts::pipi_MAX << std::endl;
  std::cout << "missing neutron window "
            << anacuts::neutron_MIN << " - " << anacuts::neutron_MAX << std::endl;
  std::cout << "Sigma Plus window "
            << anacuts::Sigmap_MIN << " - " << anacuts::Sigmap_MAX << std::endl;
  std::cout << "Sigma Minus window "
            << anacuts::Sigmam_MIN << " - " << anacuts::Sigmam_MAX << std::endl;
  std::cout << "##################################" << std::endl;

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

  //** output file 3 : npippim final-sample tree **//
  std::string outfile3 = confMan->GetOutFileName();
  outfile3.insert( outfile3.size()-5, "_npippim" );
  std::cout<<"npippim file "<<outfile3<<std::endl;
  rtFile3 = new TFile( outfile3.c_str(), "recreate" );
  rtFile3->cd();
  npippimTree = new TTree( "EventTree", "EventTree" );
  npippimTree->Branch( "mom_beam",   &mom_beam );
  npippimTree->Branch( "mom_beam_Sp",   &mom_beam_Sp  );
  npippimTree->Branch( "mom_beam_Sm",   &mom_beam_Sm  );
  npippimTree->Branch( "mom_target", &mom_target );
  npippimTree->Branch( "mom_pip", &mom_pip );
  npippimTree->Branch( "mom_pim", &mom_pim );
  npippimTree->Branch( "mom_n", &mom_n );
  npippimTree->Branch( "mom_n_Sp", &mom_n_Sp );
  npippimTree->Branch( "mom_n_Sm", &mom_n_Sm );
  npippimTree->Branch( "NeutralBetaCDH", &NeutralBetaCDH );
  npippimTree->Branch( "NeutralBetaCDH_vtx[2]", NeutralBetaCDH_vtx );
  npippimTree->Branch( "dE", &dE );
  npippimTree->Branch( "vtx_reaction", &vtx_reaction );
  npippimTree->Branch( "vtx_pip_beam", &vtx_pip_beam );
  npippimTree->Branch( "vtx_pim_beam", &vtx_pim_beam );
  npippimTree->Branch( "vtx_pip_cdc", &vtx_pip_cdc );
  npippimTree->Branch( "vtx_pim_cdc", &vtx_pim_cdc );
  npippimTree->Branch( "CA_pip",&CA_pip);
  npippimTree->Branch( "CA_pim",&CA_pim);
  //npippimTree->Branch( "run_num", &run_num );
  //npippimTree->Branch( "event_num", &event_num );
  //npippimTree->Branch( "block_num", &block_num );
  npippimTree->Branch( "kfSpmode_mom_beam",   &kfSpmode_mom_beam );
  npippimTree->Branch( "kfSpmode_mom_pip", &kfSpmode_mom_pip );
  npippimTree->Branch( "kfSpmode_mom_pim", &kfSpmode_mom_pim );
  npippimTree->Branch( "kfSpmode_mom_n", &kfSpmode_mom_n );
  npippimTree->Branch( "kfSpmode_chi2", &kfSpmode_chi2 );
  npippimTree->Branch( "kfSpmode_NDF", &kfSpmode_NDF );
  npippimTree->Branch( "kfSpmode_status", &kfSpmode_status );
  npippimTree->Branch( "kfSpmode_pvalue", &kfSpmode_pvalue );
  npippimTree->Branch( "kfSmmode_mom_beam",   &kfSmmode_mom_beam );
  npippimTree->Branch( "kfSmmode_mom_pip", &kfSmmode_mom_pip );
  npippimTree->Branch( "kfSmmode_mom_pim", &kfSmmode_mom_pim );
  npippimTree->Branch( "kfSmmode_mom_n", &kfSmmode_mom_n );
  npippimTree->Branch( "kfSmmode_chi2", &kfSmmode_chi2 );
  npippimTree->Branch( "kfSmmode_NDF", &kfSmmode_NDF );
  npippimTree->Branch( "kfSmmode_status", &kfSmmode_status );
  npippimTree->Branch( "kfSmmode_pvalue", &kfSmmode_pvalue );
  npippimTree->Branch( "kf_flag", &kf_flag );

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

  nFill_pippim = 0;
  nFill_npippim  = 0;

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
  nAbort_pipi = 0;
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

  const int nGoodTrack = trackMan->nGoodTrack();
  const int nallTrack = trackMan->nTrack();
  AllGoodTrack += nGoodTrack;
  nTrack += nallTrack;
  Tools::Fill1D( Form("nGoodTrack"), nGoodTrack );

  Tools::Fill1D( Form("EventCheck"), 1 );
  

  //CDH-hits cut
  if( Util::GetCDHMul(cdsMan,nGoodTrack)!=2){
    Clear( nAbort_nCDH );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 2 );

  //** # of good CDS tracks cut **//
  if( nGoodTrack!=2 ) { //require K
    Clear( nAbort_nGoodTrack );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 3 );


  //beam line analysis and event selection
  
  //** BLDC tracking **//
  bltrackMan->DoTracking(blMan,confMan,true,true);
  //Get T0
  ctmT0 = Util::AnalyzeT0(blMan,confMan);
  if(ctmT0<-9000){
    Clear( nAbort_nT0 );
    Tools::Fill1D( Form("EventCheck"), 15 );
    return true;
  }
  
  //PID beam
  const int beamPID = Util::BeamPID(header, ctmT0, blMan);
  if(beamPID!=Beam_Kaon){
    Clear( nAbort_pid_beam );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 4 );

  const int blstatus = Util::EveSelectBeamline(bltrackMan,trackMan,confMan,blc1GoodTrackID,blc2GoodTrackID,bpcGoodTrackID);
  
  if(blstatus == -16){
    Clear( nAbort_singleBLtrack );
    Tools::Fill1D( Form("EventCheck"), 16 );
    return true;
  }
  
  if(blstatus == -17){
    Clear( nAbort_nbpc );
    Tools::Fill1D( Form("EventCheck"), 17 );
    return true;
  }
  
  if(blstatus == -18){
    Clear( nAbort_bpctrack );
    Tools::Fill1D( Form("EventCheck"), 18 );
    return true;
  }
  
  if(blstatus == -19){
    Clear( nAbort_fblc2bpc );
    Tools::Fill1D( Form("EventCheck"), 19 );
    return true;
  }
  
  Tools::Fill1D( Form("EventCheck"), 5 );


  //BLC1-D5-BLC2 analysis and chi2 selection
  const double beammom = Util::AnaBeamSpec(confMan,bltrackMan,blc1GoodTrackID,blc2GoodTrackID);
  if(beammom <-100){
    Clear( nAbort_flagbmom );
    return true;//chi2 cut
  }
  Tools::Fill1D( Form("EventCheck"), 6 );

  //** beam momentum calculation **//
  TVector3 Pp_target(0,0,0);
  //Pp_target.SetXYZ( 0, 0, 0 );
  TLorentzVector LVec_beambf;  // 4-Momentum(beam) in LAB
  TLorentzVector LVec_beam;    // 4-Momentum(beam) in LAB with dE correcion
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
  const TVector3 boost = (LVec_target+LVec_beam).BoostVector();
  LVec_beambfCM = LVec_beam;
  LVec_targetCM = LVec_target;
  //LVec_targetPCM = LVec_targetP;
  //boost to CM frame
  LVec_beambfCM.Boost( -1*boost );
  LVec_targetCM.Boost( -1*boost );


  //** + + + + + + + + + + + + **//
  //**  PID in CDS             **//
  //** + + + + + + + + + + + + **//
  //** vectors for PID container **//
  std::vector <int> pim_ID;
  std::vector <int> pip_ID;
  std::vector <int> km_ID;
  std::vector <int> p_ID;

  std::vector <int> vCDHseg;
  // PID of CDS tracks //
  const int nIDedTrack = Util::CDSChargedAna(
    DoCDCRetiming,
    bpctrack, cdsMan, trackMan, confMan, 
    LVec_beam, ctmT0,vCDHseg,pim_ID,pip_ID,km_ID,p_ID);
  if(nIDedTrack==-7) Tools::Fill1D( Form("EventCheck"), 7 );
  if(nIDedTrack==-8) Tools::Fill1D( Form("EventCheck"), 8 );
  if(nIDedTrack==-9) Tools::Fill1D( Form("EventCheck"), 9 );
  if(nIDedTrack==-10) Tools::Fill1D( Form("EventCheck"), 10 );
  if(nIDedTrack==-11) Tools::Fill1D( Form("EventCheck"), 11 );
  if(nIDedTrack==-12) Tools::Fill1D( Form("EventCheck"), 12 );
  if(nIDedTrack<0){
    Clear(nAbort_CDSPID);
    return true;
  }
  Tools::Fill1D( Form("EventCheck"),13);
  Tools::Fill1D( Form("ntrack_CDS"), nIDedTrack );
  Tools::Fill1D( Form("ntrack_pi_plus"),  pip_ID.size() );
  Tools::Fill1D( Form("ntrack_proton"),   p_ID.size() );
  Tools::Fill1D( Form("ntrack_pi_minus"), pim_ID.size() );
  Tools::Fill1D( Form("ntrack_K_minus"),  km_ID.size() );

  //  pi+ pi- X event
  //  with CDH multiplicity selection
  bool kmFlag = false;
  if( km_ID.size()==1 && p_ID.size()==1 ) kmFlag = true;
  if( kmFlag &&
      (trackMan->nGoodTrack()==2) && !Util::IsForwardCharge(blMan)){
    //=== condense pi+ pi- X candidates ===//
    rtFile2->cd();
    evTree->Fill();
    rtFile->cd();
    if(Verbosity) std::cout<<"### filled: Event_Number, Block_Event_Number, CDC_Event_Number = "
                             <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
    nFill_pippim++;
    //=== condense pi+ pi- X event ===//

    //** find CDH hit from neutral particles **//
    std::vector <int> NeutralCDHseg;//CDHhits - CDHhits used for charged particle tracking
    std::vector <int> CDHhit_list;
    for( int icdhhit=0; icdhhit<cdsMan->nCDH(); icdhhit++ ) {
      if( cdsMan->CDH(icdhhit)->CheckRange() &&
          cdsMan->CDH(icdhhit)->ctmean()< cdscuts::tdc_cdh_max) {
        CDHhit_list.push_back( cdsMan->CDH(icdhhit)->seg() );
      }
    }
    std::sort(vCDHseg.begin(), vCDHseg.end());
    std::sort(CDHhit_list.begin(), CDHhit_list.end());
    std::set_difference( CDHhit_list.begin(), CDHhit_list.end(),
                         vCDHseg.begin(), vCDHseg.end(),
                         std::back_inserter(NeutralCDHseg) );

    if( NeutralCDHseg.size()!=1 ) {
      //if(Verbosity) {
        std::cerr<<" CDH neutral hit is not 1 :: "<<NeutralCDHseg.size()<<std::endl;
      //}
      //if(Verbosity){
        std::cerr<<"# of diff = "<<NeutralCDHseg.size()<<std::endl;
        std::cerr<<"CDH hits =   ";
        for( int n=0; n<(int)CDHhit_list.size(); n++ ) {
          std::cerr<<CDHhit_list[n]<<" ";
        }
        std::cerr<<std::endl;
        std::cerr<<"track hits = ";
        for( int n=0; n<(int)vCDHseg.size(); n++ ) {
          std::cerr<<vCDHseg[n]<<" ";
        }
        std::cerr<<std::endl;
        std::cerr<<"diff hits =  ";
        for( int n=0; n<(int)NeutralCDHseg.size(); n++ ) {
          std::cerr<<NeutralCDHseg[n]<<" ";
        }
        std::cerr<<std::endl;
      //}//if Verbosity
    }
    
    //** isolation cut **//
    int flag_isolation = 0;
    if(IsolationCutFlag==2){
      flag_isolation = Util::GetCDHNeighboringNHits(NeutralCDHseg,CDHhit_list,vCDHseg);
      flag_isolation+= Util::GetCDHTwoSegAwayNHits(NeutralCDHseg,CDHhit_list);
    }else if(IsolationCutFlag==1){
      flag_isolation = Util::GetCDHNeighboringNHits(NeutralCDHseg,CDHhit_list,vCDHseg);
    }else{
      //check cdh hit position anyway, but don't apply isolation cuts 
      flag_isolation = Util::GetCDHNeighboringNHits(NeutralCDHseg,CDHhit_list,vCDHseg);
      flag_isolation = 0;
    }

    if( flag_isolation ) {
      //if(Verbosity) std::cerr<<"CDH neutral hit candidate is NOT isolated !!!"<<std::endl;
      Clear( nAbort_CDHiso );
      return true;
    } else {
      if(Verbosity) std::cerr<<"CDH isolation cuts : OK " << std::endl;
      Tools::Fill1D( Form("EventCheck"), 14 );
    }

    // copy neutral CDH hit candidate
    int cdhcan = -1;
    for( int ihit=0; ihit<cdsMan->nCDH(); ihit++ ) {
      if( cdsMan->CDH(ihit)->seg()==NeutralCDHseg[0] ) cdhcan = ihit;
    }
    HodoscopeLikeHit *ncdhhit = cdsMan->CDH(cdhcan);

    TVector3 Pos_CDH;
    confMan->GetGeomMapManager()->GetPos( CID_CDH, ncdhhit->seg(), Pos_CDH );
    if(Verbosity) {
      std::cerr<<"CDH candidate = "<<ncdhhit->seg()<<" -> "<<Pos_CDH.Phi()/TwoPi*360.<<" deg"<<std::endl;
    }
    
    // charge veto using CDC
    const int nCDCforVeto = Util::GetNHitsCDCOuter(Pos_CDH,cdsMan,cdscuts::chargevetoangle);
    //Util::AnaPipPimCDCCDH(Pos_CDH,NeutralCDHseg,pip_ID[0],pim_ID[0],cdsMan,trackMan);
    
    Pos_CDH.SetZ(-1*ncdhhit->hitpos()); // (-1*) is correct in data analysis [20170926]
    //Pos_CDH.SetZ(ncdhhit->hitpos()); // (-1*) is correct in data analysis [20170926]


    //** neutral particle in CDH **//
    if( !nCDCforVeto ) {
      if(NeutralCDHseg.size()!=1) {
        std::cout << "L." << __LINE__ << " # of seg for neutral hits " << NeutralCDHseg.size() << std::endl;
      } else {
        Tools::Fill1D(Form("CDHNeutralSeg"),NeutralCDHseg.at(0));
      }

      //CDSTrack *track_pip = trackMan->Track( pip_ID.at(0) ); // should be only 1 track
      //CDSTrack *track_pim = trackMan->Track( pim_ID.at(0) ); // should be only 1 track
      CDSTrack *track_km = trackMan->Track( km_ID.at(0) ); // should be only 1 track

      //deuteron target
      TVector3 vtx_react;//reaction vertex
      TVector3 vtx_dis;//displaced vertex
      
      
      TVector3 vtx_beam_wkm;//vertex(beam-pip) on beam
      TVector3 vtx_km;//vertex(beam-pip) on beam
      track_km->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wkm, vtx_km );
      
      //reaction vertex is determined from beam and nearest vtx
      TVector3 vtx_beam;
      vtx_react = 0.5*(vtx_km+vtx_beam_wkm);
      vtx_beam = vtx_beam_wkm;


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
      const double ntof = ncdhhit->ctmean()-ctmT0-beamtof;
      double nlen;
      nlen = (Pos_CDH-vtx_react).Mag();

      //std::cout << "nlen "  << nlen << std::endl;
      //std::cout << "Pos_CDH x" << Pos_CDH.x() << std::endl;
      //std::cout << "Pos_CDH y" << Pos_CDH.y() << std::endl;
      //std::cout << "Pos_CDH z" << Pos_CDH.z() << std::endl;

      NeutralBetaCDH = nlen/ntof/(Const*100.);
      double tmp_mom = NeutralBetaCDH<1. ? nMass*NeutralBetaCDH/sqrt(1.-NeutralBetaCDH*NeutralBetaCDH) : 0;
      if(Verbosity) {
        std::cerr<<"L. " << __LINE__ ;
        std::cerr<<" NeutralBetaCDH = "<<NeutralBetaCDH<<" mom_n = "<<tmp_mom<<std::endl; //" "<<1/sqrt(1+nMass*nMass)<<std::endl;
      }

      //** reconstructoin of missing neutorn **//
      TVector3 P_km; // Momentum(pi-)

      TLorentzVector LVec_km; // 4-Momentum(pi-)
      TLorentzVector LVec_n;   // 4-Momentum(n) (CDS)
      TLorentzVector LVec_pmiss; // 4-Momentum(n_miss)
      
      //energy loss correction and momentum correction using vertex info
      if( !track_km->GetMomentum( vtx_km, P_km, true, true ) ) {
        std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }


      //Momentum (n CDS)
      TVector3 P_n;
      P_n = tmp_mom*(Pos_CDH-vtx_react).Unit();

      LVec_km.SetVectM( P_km, kpMass );
      LVec_n.SetVectM(   P_n,   nMass );//CDS n
      
      const double mm_mass   = (LVec_target+LVec_beam-LVec_km-LVec_n).M();
      const TVector3 P_missp = (LVec_target+LVec_beam-LVec_km-LVec_n).Vect();
      LVec_pmiss.SetVectM( P_missp, pMass );

      //** + + + + + + + + + + + + + **//
      //**  fill histograms & tree   **//
      //** + + + + + + + + + + + + + **//
      kf_flag = -1;
      bool MissPFlag=false;
      bool NBetaOK=false;
      bool NdEOK=false;

      Tools::Fill2D( Form("dE_betainv"), 1./NeutralBetaCDH, ncdhhit->emean() );
      Tools::Fill2D( Form("MMom_MMass"), mm_mass, P_missp.Mag() );
      
      //Fiducial cuts OK
      if( GeomTools::GetID(vtx_react)==CID_Fiducial){
        for( int i=0; i<cdsMan->nCDH(); i++ ) {
          Tools::Fill2D(Form("dE_CDHtime_pippimn"), cdsMan->CDH(i)->ctmean(), cdsMan->CDH(i)->emean());
        }

        Tools::Fill2D(Form("Vtx_ZX_fid"),vtx_react.Z(),vtx_react.X());
        Tools::Fill2D(Form("Vtx_ZY_fid"),vtx_react.Z(),vtx_react.Y());
        Tools::Fill2D(Form("Vtx_XY_fid"),vtx_react.X(),vtx_react.Y());


        Tools::Fill2D(Form("NeutraltimeEnergy"),ncdhhit->ctmean()-ctmT0-beamtof,ncdhhit->emean());
        Tools::Fill2D(Form("CDHzNeutraltime"),Pos_CDH.z(),ncdhhit->ctmean()-ctmT0-beamtof);
        Tools::Fill2D( Form("dE_betainv_fid"), 1./NeutralBetaCDH, ncdhhit->emean() );
        Tools::Fill2D( Form("MMom_MMass_fid"), mm_mass, P_missp.Mag() );

        if(NeutralBetaCDH<anacuts::beta_MAX) NBetaOK=true;

        if(NBetaOK) {
          Tools::Fill2D( Form("dE_betainv_fid_beta"), 1./NeutralBetaCDH, ncdhhit->emean() );
          Tools::Fill2D( Form("MMom_MMass_fid_beta"), mm_mass, P_missp.Mag() );
        }
        if(anacuts::dE_MIN<ncdhhit->emean()) NdEOK=true;

        if( NBetaOK && NdEOK ) {
          Tools::Fill2D( Form("dE_betainv_fid_beta_dE"), 1./NeutralBetaCDH, ncdhhit->emean() );
          Tools::Fill2D( Form("MMom_MMass_fid_beta_dE"), mm_mass, P_missp.Mag() );
        }

        //missing mass neutron ID
        if( anacuts::neutron_MIN<mm_mass && mm_mass<anacuts::neutron_MAX ) MissPFlag=true;

        if( NBetaOK && NdEOK ) {
          Tools::Fill2D( Form("dE_betainv_fid_beta_dE_woK0"), 1./NeutralBetaCDH, ncdhhit->emean() );
          Tools::Fill2D( Form("MMom_MMass_fid_beta_dE_woK0"), mm_mass, P_missp.Mag() );
          Tools::Fill2D(Form("CDHseg_MMass_fid_beta_dE_woK0"),ncdhhit->seg(),mm_mass);
          Tools::Fill2D(Form("CDHz_MMass_fid_beta_dE_woK0"),-1*ncdhhit->hitpos(),mm_mass);
          //Tools::Fill2D(Form("zVTX_MMass_fid_beta_dE_woK0"),vtx_react.z(),mm_mass);
          Tools::Fill2D( Form("dE_MMom_fid_beta_woK0"), P_missp.Mag(), ncdhhit->emean() );
          Tools::Fill2D( Form("dE_MMass_fid_beta_woK0"), mm_mass, ncdhhit->emean() );

          if(MissPFlag) {
            Tools::Fill2D( Form("dE_betainv_fid_beta_dE_woK0_n"), 1./NeutralBetaCDH, ncdhhit->emean() );
            Tools::Fill2D( Form("MMom_MMass_fid_beta_dE_woK0_n"), mm_mass, P_missp.Mag() );
            Tools::Fill2D( Form("NMom_NMom_fid_beta_dE_woK0_n"), P_n.Mag(), P_missp.Mag() );
          }

        }//MissPFlag && K0rejectFlag && (SigmaPFlag || SigmaMFlag)
          

        if(Verbosity>10)std::cout<<"%%% npippim event: Event_Number, Block_Event_Number, CDC_Event_Number = "
                                   <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
        rtFile3->cd();
        npippimTree->Fill();
        rtFile->cd();
        nFill_npippim++;
        //** fill tree **//
      } // if( GeomTools::GetID(vtx_react)==CID_Fiducial )
    } // if( !nCDCforVeto )
  }//pi+,pi-X event
  else {
    Clear( nAbort_pipi );
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
  std::cout<<" nAbort_CDHiso        = "<<nAbort_CDHiso<<std::endl;
  std::cout<<" nAbort_nAbort_pipi   = "<<nAbort_pipi<<std::endl;
  std::cout<<" nAbort_end           = "<<nAbort_end<<std::endl;
  std::cout<<"========= Abort counter ========="<<std::endl;
  std::cout<<"*** # of pi+ pi- n events = "<<nFill_npippim<<" ***"<<std::endl;

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
  InitIMPiSigmaHist();


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

  //CDS
  NeutralBetaCDH=-9999.;
  NeutralBetaCDH_vtx[0]=-9999.;
  NeutralBetaCDH_vtx[1]=-9999.;

  dE=-9999.;
  vtx_reaction.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pip_beam.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pim_beam.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pip_cdc.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pim_cdc.SetXYZ(-9999.,-9999.,-9999.);
  
  return;
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
EventTemp *EventAlloc::EventAllocator()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
