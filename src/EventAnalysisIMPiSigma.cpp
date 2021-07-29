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
const unsigned int IsolationCutFlag = 0;
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
  int nAbort_KCDH3trg;
  int nAbort_Kf;
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
  int nAbort_CDCInner3Lay;
  int nAbort_pipi;
  int nAbort_end;

  // npippim final-sample tree Branch //
  
  // 4-momentum(beam) (reaction vtx determined by DCA)
  TLorentzVector mom_beam;   
  // 4-momentum(beam) (reaction vtx Sp mode assumption)
  TLorentzVector mom_beam_Sp;  
  // 4-momentum(beam) (reaction vtx Sm mode assumption)
  TLorentzVector mom_beam_Sm;
  // 4-momentum(beam) (reaction vtx K0 mode assumption)
  TLorentzVector mom_beam_K0;
  TLorentzVector mom_target; // 4-momentum(target)
  TLorentzVector mom_pip;    // 4-momentum(pi+)
  TLorentzVector mom_pim;    // 4-momentum(pi-)
  TLorentzVector mom_n;      // 4-momentum(neutron)
  TLorentzVector mom_n_beam;      // 4-momentum(neutron)
  TLorentzVector mom_n_Sp;  // 4-momentum(neutron),Sp mode assumption
  TLorentzVector mom_n_Sm;  // 4-momentum(neutron),Sm mode assumption
  TLorentzVector mom_n_K0;  // 4-momentum(neutron),K0 mode assumption
  double NeutralBetaCDH; // velocity of neutral particle on CDH from decay vtx 
  double NeutralBetaCDH_beam; // velocity of neutral particle on CDH from decay vtx 
  double NeutralBetaCDH_vtx[3];//1:pip_vtx,2:pim_vtx  
  double tofpim;
  double tofpip;
  double tofn;
  double dE;   // energy deposit on CDH [MeVee] (neutral candidate)
  int neutralseg;   //neutral candidate segment
  int nhitOutCDC;
  int ForwardCharge;
  TVector3 vtx_reaction; // 
  TVector3 vtx_displaced; // 
  TVector3 vtx_pip_beam; // 
  TVector3 vtx_pim_beam; // 
  TVector3 vtx_pip_cdc;//
  TVector3 vtx_pim_cdc;//
  TVector3 CA_pip;//Closest Approach Point of CDS pip-pim tracks
  TVector3 CA_pim;//Closest Approach Point of CDS pip-pim tracks
  TVector3 CDH_Pos;//neutron candidate
  TVector3 CDH_Pos_pim;//pim CDC track projected position at CDH
  TVector3 CDH_Pos_pip;//pip CDC track projected position at CDH
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
  npippimTree->Branch( "mom_beam_K0",   &mom_beam_K0  );
  npippimTree->Branch( "mom_target", &mom_target );
  npippimTree->Branch( "mom_pip", &mom_pip );
  npippimTree->Branch( "mom_pim", &mom_pim );
  npippimTree->Branch( "mom_n", &mom_n );
  npippimTree->Branch( "mom_n_beam", &mom_n_beam );
  npippimTree->Branch( "mom_n_Sp", &mom_n_Sp );//decay point is pip-CA
  npippimTree->Branch( "mom_n_Sm", &mom_n_Sm );//decay point is pim-CA
  npippimTree->Branch( "mom_n_K0", &mom_n_K0 );//decay point = reaction point
  npippimTree->Branch( "NeutralBetaCDH", &NeutralBetaCDH );
  npippimTree->Branch( "NeutralBetaCDH_beam", &NeutralBetaCDH_beam );
  npippimTree->Branch( "NeutralBetaCDH_vtx[3]", NeutralBetaCDH_vtx );
  npippimTree->Branch( "tofpim",&tofpim);
  npippimTree->Branch( "tofpip",&tofpip);
  npippimTree->Branch( "tofn",&tofn);
  npippimTree->Branch( "dE", &dE );
  npippimTree->Branch( "neutralseg", &neutralseg );
  npippimTree->Branch( "nhitOutCDC", &nhitOutCDC );
  npippimTree->Branch( "ForwardCharge", &ForwardCharge);
  npippimTree->Branch( "vtx_reaction", &vtx_reaction );
  npippimTree->Branch( "vtx_displaced", &vtx_displaced );
  npippimTree->Branch( "vtx_pip_beam", &vtx_pip_beam );
  npippimTree->Branch( "vtx_pim_beam", &vtx_pim_beam );
  npippimTree->Branch( "vtx_pip_cdc", &vtx_pip_cdc );
  npippimTree->Branch( "vtx_pim_cdc", &vtx_pim_cdc );
  npippimTree->Branch( "CA_pip",&CA_pip);
  npippimTree->Branch( "CA_pim",&CA_pim);
  npippimTree->Branch( "CDH_Pos",&CDH_Pos);
  npippimTree->Branch( "CDH_Pos_pip",&CDH_Pos_pip);
  npippimTree->Branch( "CDH_Pos_pim",&CDH_Pos_pim);
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
  
  nAbort_KCDH3trg = 0;
  nAbort_Kf = 0;
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
  nAbort_CDCInner3Lay = 0;
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

  static bool FirstScalerRead=true;
  for( int i=0; i<scaMan->nsca(); i++ ) {
    TString name = scaMan->sca(i)->name();
    //test 
    int val =scaMan->sca(i)->val()-scaend[i];
    if(FirstScalerRead) {
      val -= scainit[i];
      std::cout << i << " init  " << scainit[i] << std::endl;
    }
    TH1F* h1 = (TH1F*)gFile->Get(Form("SCA%d",i));
    h1->Fill(Event_Number,val);
    
    if( scaend[i]>9.9e+07 && scaend[i]>scaMan->sca(i)->val() ) SCAOVERFLOW[i] = true;
    scaend[i] = scaMan->sca(i)->val();
    if( SCAOVERFLOW[i] ) scaend[i] += 1.0e+08;
  
  }
  FirstScalerRead =false;

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


  // fill event count //
  Event_Number++;

  // control of start and stop events //
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


  // Checking event number consistentcy with CDC tracking file //
  if( CDC_Event_Number>=cdcTree->GetEntries() ) return false;
  cdcTree->GetEntry( CDC_Event_Number );
  if( header_CDC->ev()!=Event_Number ) return true;

  CDC_Event_Number++;
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  // copy class CDSTrackingMan objects //
  *trackMan = *trackMan_CDC;

  // Converting data=>tree and re-calc. of trackMan of CDS //
  header->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  trackMan->Calc( cdsMan, confMan, true);
  
  //trigger check
  //
  //unbiased kaon trigger,prescaled
  //
  //if( header->IsTrig(Trig_Kf)){
  //  Tools::Fill1D(Form("Trigger"),0);
  //  Clear(nAbort_Kf);
  //  return true;
  //}
  //std::cout << header->trigmode() << std::endl;
  //if(header->IsTrig(Trig_Kf)){
  Tools::Fill1D(Form("Trigmode"),header->trigmode());
  if(header->trigmode(Mode_Kf)){
    Tools::Fill1D(Form("Trigger"),0);
    Tools::Fill1D(Form("Trigmode_Kf"),header->trigmode());
  }else{ 
  //std::cout << "abort" << std::endl;
  //  Clear(nAbort_Kf);
  //  return true;
  }
  Tools::Fill1D( Form("EventCheck"), 1 );
  if( header->IsTrig(Trig_KCDH2f)) Tools::Fill1D(Form("Trigger"),2);   
  //K x CDH3 trigger
  bool IsTrigKCDH3 = header->IsTrig(Trig_KCDH3);
  if(IsTrigKCDH3){
    Tools::Fill1D(Form("Trigger"),1);
  }else{
    Clear(nAbort_KCDH3trg);
    return true;
  }


  const int nGoodTrack = trackMan->nGoodTrack();
  const int nallTrack = trackMan->nTrack();
  AllGoodTrack += nGoodTrack;
  nTrack += nallTrack;
  Tools::Fill1D( Form("nTrack"),nallTrack);
  Tools::Fill1D( Form("nGoodTrack"), nGoodTrack);
  if(nGoodTrack==2){
    Tools::Fill1D( Form("nTrack_If2GoodTracks"),nallTrack);
  }
  Tools::Fill1D( Form("EventCheck"), 1 );

  //CDH emean recalibrator
  static bool isState=false;
  if(!isState){ 
    std::cout << "***********************************************" << std::endl;   
    std::cout << "L." << __LINE__ << " CDH emean is recalibrated " << std::endl;
    std::cout << "correction factor " << cdscuts::CDHemeanCal << std::endl;
    std::cout << "***********************************************" << std::endl;   
    isState=true;
  }
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    double emean = cdsMan->CDH(i)->emean();
    cdsMan->CDH(i)->SetEMean(emean*cdscuts::CDHemeanCal);
  }
  
  //CDH-hits cut
  if( Util::GetCDHMul(cdsMan,nGoodTrack,IsTrigKCDH3,false)!=cdscuts::cdhmulti){
    Clear( nAbort_nCDH );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 2 );

  // # of good CDS tracks cut //
  if( nGoodTrack!=cdscuts::cds_ngoodtrack  ) { //require pi+,pi-
  //if( nGoodTrack!=cdscuts::cds_ngoodtrack && nallTrack!=cdscuts::cds_ngoodtrack ) { //require pi+,pi-
    Clear( nAbort_nGoodTrack );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 3 );
  

  //beam line analysis and event selection
  
  //** BLDC tracking **//
  bltrackMan->DoTracking(blMan,confMan,true,true);
  //Get T0
  int t0seg=-1;
  //return -9999 if nhit T0 =>2 
  ctmT0 = Util::AnalyzeT0(blMan,confMan,t0seg);
  if(ctmT0<-9000){
    Clear( nAbort_nT0 );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"),15);
  
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
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 16 );

  if(blstatus == -17){
    Clear( nAbort_nbpc );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 17 );
  
  if(blstatus == -18){
    Clear( nAbort_bpctrack );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 18 );
  
  if(blstatus == -19){
    Clear( nAbort_fblc2bpc );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 19 );
  

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
  TLorentzVector LVec_beam_vtx[3];    // 4-Momentum(beam) in LAB with dE correcion
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
  
  Tools::H2(Form("bpcVtx_nofid"),x1,y1,500,-12.5,12.5,500,-12.5,12.5);
  if(GeomTools::GetID(lp)==CID_Fiducial){
    Tools::H2(Form("bpcVtx_fid"),x1,y1,500,-12.5,12.5,500,-12.5,12.5);
  }

  LVec_beambf.SetVectM( Pp_beam, kpMass );
  LVec_beam = LVec_beambf;
  LVec_beam_vtx[0] = LVec_beambf;
  LVec_beam_vtx[1] = LVec_beambf;
  LVec_beam_vtx[2] = LVec_beambf;
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

  std::vector <int> vCDHseg;//vector of CDH seg. 
  // PID of CDS tracks //
  TVector3 pim_cdhprojected;
  TVector3 pip_cdhprojected;
  const int nIDedTrack = Util::CDSChargedAna(
    DoCDCRetiming,
    bpctrack, cdsMan, trackMan, confMan,blMan, 
    LVec_beam, ctmT0,vCDHseg,pim_ID,pip_ID,km_ID,p_ID,pim_cdhprojected,pip_cdhprojected);
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
  bool pimpipFlag = false;
  if( pim_ID.size()==1 && pip_ID.size()==1) pimpipFlag = true;
  if( pimpipFlag &&
      (trackMan->nGoodTrack()==cdscuts::cds_ngoodtrack)) { //&& !Util::IsForwardCharge(blMan)){
    //=== pi+ pi- X candidates ===//
    rtFile2->cd();
    evTree->Fill();
    rtFile->cd();
    
    for( int i=0; i<cdsMan->nCDH(); i++ ){
      if((cdsMan->CDH(i)->CheckRange()) && (cdsMan->CDH(i)->ctmean()<cdscuts::tdc_cdh_max)){
        Tools::Fill2D( Form("CDHdE_pippim"),cdsMan->CDH(i)->seg(),cdsMan->CDH(i)->emean());
      }
    }

    //added Jul.28th,2019
    //purpose
    double tofpim_b = 0;
    double tofpip_b = 0;
    for( int it=0; it<trackMan->nGoodTrack(); it++ ) {
      CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );

      double mom = track->Momentum();//charge X momentum
      double tof = 999.;//TOF of CDH-T0 (slewing corrected)
      double mass2 = -999.;
      double correctedtof=0;//CDH-T0 (corrected by energy loss)
      double beta_calc=0;
      for( int icdh=0; icdh<track->nCDHHit(); icdh++ ) {
        HodoscopeLikeHit *cdhhit = track->CDHHit( cdsMan, icdh );
        double tmptof = cdhhit->ctmean()-ctmT0;
        if( tmptof<tof || tof==999. ) {
          tof = tmptof;
        }
      }
      if( !TrackTools::FindMass2( track, bpctrack, tof, LVec_beam.Vect().Mag(),
            Beam_Kaon, beta_calc, mass2, correctedtof ) ) {
      }
      //const int pid = TrackTools::PIDcorr_wide(mom,mass2);
      const int pid = TrackTools::PIDcorr(mom,mass2);
      track->SetPID( pid );
      Tools::Fill2D( "PID_CDS_beta_select", 1./beta_calc, mom );
      Tools::Fill2D( "PID_CDS_select", mass2, mom );
      if(pid == CDS_PiMinus){
        Tools::Fill2D("PID_CDS_PIM_beta_select",1./beta_calc,mom);
        Tools::Fill2D("PID_CDS_PIM_select",mass2,mom);
        tofpim = tof;//tree val.
      }else if(pid == CDS_PiPlus){
        Tools::Fill2D("PID_CDS_PIP_beta_select",1./beta_calc,mom);
        Tools::Fill2D("PID_CDS_PIP_select",mass2,mom);
        tofpip = tof;//tree val.
      }
      else if(pid == CDS_Proton) Tools::Fill2D("PID_CDS_Proton_select",mass2,mom);
      else if(pid == CDS_Kaon) Tools::Fill2D("PID_CDS_Kaon_select",mass2,mom);
    }//nGoodTrack
    if(Verbosity) std::cout<<"### filled: Event_Number, Block_Event_Number, CDC_Event_Number = "
                             <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
    nFill_pippim++;

    //** find CDH hit from neutral particles **//
    std::vector <int> NeutralCDHseg;//CDHhits - CDHhits used for charged particle tracking
    std::vector <int> CDHhitsegall;
    for( int icdhhit=0; icdhhit<cdsMan->nCDH(); icdhhit++ ) {
      if( cdsMan->CDH(icdhhit)->CheckRange() &&
          cdsMan->CDH(icdhhit)->ctmean()< cdscuts::tdc_cdh_max) {
        CDHhitsegall.push_back( cdsMan->CDH(icdhhit)->seg() );
      }
    }
    std::sort(vCDHseg.begin(), vCDHseg.end());
    std::sort(CDHhitsegall.begin(), CDHhitsegall.end());
    std::set_difference( CDHhitsegall.begin(), CDHhitsegall.end(),
                         vCDHseg.begin(), vCDHseg.end(),
                         std::back_inserter(NeutralCDHseg) );

    if( NeutralCDHseg.size()!=1 ) {
      //if(Verbosity) {
        std::cerr<<" CDH neutral hit is not 1 :: "<<NeutralCDHseg.size()<<std::endl;
      //}
      if(Verbosity){
        std::cerr<<"# of diff = "<<NeutralCDHseg.size()<<std::endl;
        std::cerr<<"CDH hits =   ";
        for( int n=0; n<(int)CDHhitsegall.size(); n++ ) {
          std::cerr<<CDHhitsegall[n]<<" ";
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
      }//if Verbosity
    }//NeutralCDHseg check
    
    //** isolation cut **//
    int flag_isolation = 0;
    if(IsolationCutFlag==2){
      flag_isolation = Util::GetCDHNeighboringNHits(NeutralCDHseg,CDHhitsegall,vCDHseg,cdsMan);
      flag_isolation+= Util::GetCDHTwoSegAwayNHits(NeutralCDHseg,CDHhitsegall);
    }else if(IsolationCutFlag==1){
      flag_isolation = Util::GetCDHNeighboringNHits(NeutralCDHseg,CDHhitsegall,vCDHseg,cdsMan);
    }else{
      //check cdh hit position anyway, but don't apply isolation cuts 
      flag_isolation = Util::GetCDHNeighboringNHits(NeutralCDHseg,CDHhitsegall,vCDHseg,cdsMan);
      flag_isolation = 0;
    }


    // copy neutral CDH hit candidate
    int cdhcan = -1;
    for( int ihit=0; ihit<cdsMan->nCDH(); ihit++ ) {
      if( cdsMan->CDH(ihit)->seg()==NeutralCDHseg[0] ) cdhcan = ihit;
    }
    HodoscopeLikeHit *ncdhhit = cdsMan->CDH(cdhcan);

    TVector3 Pos_CDH;
    const int CDHSeg = ncdhhit->seg();
    confMan->GetGeomMapManager()->GetPos( CID_CDH, CDHSeg, Pos_CDH );


    if(Verbosity) {
      std::cerr<<"CDH candidate = "<<CDHSeg<<" -> "<<Pos_CDH.Phi()/TwoPi*360.<<" deg"<<std::endl;
    }
    
    // charge veto using CDC
    Util::AnaPipPimCDCCDH(Pos_CDH,NeutralCDHseg,pip_ID[0],pim_ID[0],cdsMan,trackMan);
    
    //std::cout << __LINE__ << "  "  << -1.*ncdhhit->hitpos() << std::endl;
    if( flag_isolation ) {
      //if(Verbosity) std::cerr<<"CDH neutral hit candidate is NOT isolated !!!"<<std::endl;
      Clear( nAbort_CDHiso );
      return true;
    } else {
      if(Verbosity) std::cerr<<"CDH isolation cuts : OK " << std::endl;
      Tools::Fill1D( Form("EventCheck"), 14 );
    }
    const int nCDCInner3Lay = Util::GetNHitsCDCInner3Lay(cdsMan);
    Tools::Fill1D(Form("CDCInner3Mul"),nCDCInner3Lay);
     
    //if(nCDCInner3Lay>6){
    //  Clear(nAbort_CDCInner3Lay);
    //  return true;
    //}

    const int nCDCforVeto = Util::GetNHitsCDCOuterNoAss(Pos_CDH,cdsMan,trackMan,cdscuts::chargevetoangle);
    Tools::Fill1D(Form("NCDCOutHit"),nCDCforVeto);
    Pos_CDH.SetZ(-1.*ncdhhit->hitpos()); // (-1*) is correct in data analysis [20170926]
    //** neutral particle in CDH **//
    if(NeutralCDHseg.size()!=1) {
        std::cout << "L." << __LINE__ << " # of seg for neutral hits " << NeutralCDHseg.size() << std::endl;
    } else {
      Tools::Fill1D(Form("CDHNeutralSeg"),NeutralCDHseg.at(0));
    }

    CDSTrack *track_pip = trackMan->Track( pip_ID.at(0) ); // should be only 1 track
    CDSTrack *track_pim = trackMan->Track( pim_ID.at(0) ); // should be only 1 track

    //vertex calculation
    TVector3 vtx_react;//reaction vertex
    TVector3 vtx_dis;//displaced vertex

    TVector3 vtx_beam_wpip;//Closest approach(beam-pip) on beam
    TVector3 vtx_pip;//Closest approach(beam-pip) on CDC track
    track_pip->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wpip, vtx_pip );

    TVector3 vtx_beam_wpim;//Closest approach (beam-pim) on beam
    TVector3 vtx_pim;//Closest approach (beam-pim) on CDC track
    track_pim->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wpim, vtx_pim );

    const double dcapipvtx =  (vtx_pip-vtx_beam_wpip).Mag();
    const double dcapimvtx =  (vtx_pim-vtx_beam_wpim).Mag();
    const TVector3 vtxpip_mean = 0.5*(vtx_pip+vtx_beam_wpip);
    const TVector3 vtxpim_mean = 0.5*(vtx_pim+vtx_beam_wpim);
    Tools::Fill1D( Form("DCA_pip"), dcapipvtx );
    Tools::Fill1D( Form("DCA_pim"), dcapimvtx );

    TVector3 CA_pip_pippim,CA_pim_pippim;
    bool vtx_flag=TrackTools::Calc2HelixVertex(track_pip, track_pim, CA_pip_pippim, CA_pim_pippim);
    double dcapippim=-9999.;
    if(vtx_flag) dcapippim = (CA_pim_pippim-CA_pip_pippim).Mag();
    Tools::Fill1D( Form("DCA_pippim"), dcapippim);


    //reaction vertex is determined from beam and nearest vtx
    TVector3 vtx_beam;
    if(dcapipvtx < dcapimvtx) {
      //follows sakuma/sada's way , avg. of scattered particle ana beam particle [20180829]
      vtx_react = 0.5*(vtx_pip+vtx_beam_wpip);
      //if(cdscuts::useclosestpi) vtx_dis  = vtx_pip;
      //else              vtx_dis  = vtx_pim;
      //vtx_dis = CA_pim_pippim;
      vtx_dis = vtx_pim;
      vtx_beam = vtx_beam_wpip;
    } else {
      vtx_react = 0.5*(vtx_pim+vtx_beam_wpim);
      //if(cdscuts::useclosestpi) vtx_dis = vtx_pim;
      //else             vtx_dis = vtx_pip;
      //vtx_dis = CA_pip_pippim;
      vtx_dis = vtx_pip;
      vtx_beam = vtx_beam_wpim;
    }


    //** beam kaon tof **//
    TVector3 Pos_T0;
    confMan->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
    double beamtof=0;
    double momout=0;
    const double z_pos = Pos_T0.Z();;
    
    //std::cout << "test" << std::endl;
    //std::cout << "bpctrack->GetPosatZ(0)" << std::endl;
    //std::cout << bpctrack->GetPosatZ(0).X() << "  " <<  bpctrack->GetPosatZ(0).Y()  << "  " << bpctrack->GetPosatZ(0).Z()  << std::endl;
    //std::cout << zPos_T0 << std::endl;
    //std::cout << bpctrack->GetPosatZ(zPos_T0 ).X() << "  " <<  bpctrack->GetPosatZ(zPos_T0).Y()  << "  " << bpctrack->GetPosatZ(zPos_T0).Z()  << std::endl;
    
    //Energy loss correction of beam
    ELossTools::CalcElossBeamTGeo( bpctrack->GetPosatZ(z_pos), vtx_react,
        LVec_beambf.Vect().Mag(), kpMass, momout, beamtof );
    LVec_beam.SetVectM( momout*LVec_beambf.Vect().Unit(), kpMass );
    double beamtof_vtx[3];
    double momout_vtx[3];
    //Sp mode assumption
    ELossTools::CalcElossBeamTGeo( bpctrack->GetPosatZ(z_pos), 0.5*(vtx_pim+vtx_beam_wpim),
        LVec_beambf.Vect().Mag(), kpMass, momout_vtx[0], beamtof_vtx[0] );
    LVec_beam_vtx[0].SetVectM( momout_vtx[0]*LVec_beambf.Vect().Unit(), kpMass );
    //Sm mode assumption
    ELossTools::CalcElossBeamTGeo( bpctrack->GetPosatZ(z_pos), 0.5*(vtx_pip+vtx_beam_wpip),
        LVec_beambf.Vect().Mag(), kpMass, momout_vtx[1], beamtof_vtx[1] );
    LVec_beam_vtx[1].SetVectM( momout_vtx[1]*LVec_beambf.Vect().Unit(), kpMass );
    //K0 mode assumption
    ELossTools::CalcElossBeamTGeo( bpctrack->GetPosatZ(z_pos), 0.5*(vtx_beam_wpim+vtx_beam_wpip),
        LVec_beambf.Vect().Mag(), kpMass, momout_vtx[2], beamtof_vtx[2] );
    LVec_beam_vtx[2].SetVectM( momout_vtx[2]*LVec_beambf.Vect().Unit(), kpMass );


    //here flight time of Sigma is ignored.
    const double ntof = ncdhhit->ctmean()-ctmT0-beamtof;
    const double ntof_vtx[3] = {
      ncdhhit->ctmean()-ctmT0-beamtof_vtx[0], 
      ncdhhit->ctmean()-ctmT0-beamtof_vtx[1],
      ncdhhit->ctmean()-ctmT0-beamtof_vtx[2]
    };
    Tools::Fill1D(Form("CDH%d_T0%d_TOF_Neutral",CDHSeg,t0seg),ntof);
    double nlen;
    if(UseDecayVtx) nlen = (Pos_CDH-vtx_dis).Mag();  
    else nlen = (Pos_CDH-vtx_react).Mag();

    Tools::Fill2D(Form("ntof_nlen"),ntof,nlen);

    //subtract T0-target beam tof
    tofpim -=beamtof;
    tofpip -=beamtof;
    tofn = ntof;//tree val.
    //nlen_vtx[0] = (Pos_CDH-vtx_pip).Mag();
    //nlen_vtx[1] = (Pos_CDH-vtx_pim).Mag();
    double nlen_beam = (Pos_CDH-vtx_beam).Mag();

    double nlen_vtx[3];
    //nlen_vtx[0] = (Pos_CDH-CA_pip_pippim).Mag();
    //nlen_vtx[1] = (Pos_CDH-CA_pim_pippim).Mag();
    nlen_vtx[0] = (Pos_CDH-vtx_pip).Mag();
    nlen_vtx[1] = (Pos_CDH-vtx_pim).Mag();
    nlen_vtx[2] = (Pos_CDH-vtx_react).Mag();
    if(Verbosity>10) std::cout << "L." << __LINE__ << " flight length " << nlen << std::endl;
    NeutralBetaCDH = nlen/ntof/(Const*100.);
    NeutralBetaCDH_beam = nlen_beam/ntof/(Const*100.);
    for(int ivtx=0;ivtx<3;ivtx++){
      NeutralBetaCDH_vtx[ivtx] = nlen_vtx[ivtx]/ntof_vtx[ivtx]/(Const*100.);
    }
    double tmp_mom = NeutralBetaCDH<1. ? nMass*NeutralBetaCDH/sqrt(1.-NeutralBetaCDH*NeutralBetaCDH) : 0;
    double tmp_mom_beam = NeutralBetaCDH_beam<1. ? nMass*NeutralBetaCDH_beam/sqrt(1.-NeutralBetaCDH_beam*NeutralBetaCDH_beam) : 0;
    double tmp_mom_vtx[3];
    for(int ivtx=0;ivtx<3;ivtx++){
      tmp_mom_vtx[ivtx] = NeutralBetaCDH_vtx[ivtx]<1. ? nMass*NeutralBetaCDH_vtx[ivtx]/sqrt(1.-NeutralBetaCDH_vtx[ivtx]*NeutralBetaCDH_vtx[ivtx]) : 0;
    }
    if(Verbosity) {
      std::cerr<<"L. " << __LINE__ ;
      std::cerr<<" NeutralBetaCDH = "<<NeutralBetaCDH<<" mom_n = "<<tmp_mom<<std::endl; //" "<<1/sqrt(1+nMass*nMass)<<std::endl;
    }

    //** reconstructoin of missing neutorn **//
    TVector3 P_pim; // Momentum(pi-)
    TVector3 P_pip; // Momentum(pi+)

    TLorentzVector LVec_pim; // 4-Momentum(pi-)
    TLorentzVector LVec_pip; // 4-Momentum(pi+)
    TLorentzVector LVec_n;   // 4-Momentum(n) (CDS)
    TLorentzVector LVec_n_beam;   // 4-Momentum(n) (CDS)
    TLorentzVector LVec_n_vtx[3];   // 4-Momentum(n) (CDS) //decay vertex 0: pip (Spmode), 1:pim (Smmode), 2:reaction vertex (K0 mode)
    TLorentzVector LVec_nmiss; // 4-Momentum(n_miss)
    TLorentzVector LVec_nmiss_vtx[3]; // 4-Momentum(n_miss)

    //energy loss correction and momentum correction using vertex info
    if( !track_pip->GetMomentum( vtx_pip, P_pip, true, true ) ) {
      std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
    }

    if( !track_pim->GetMomentum( vtx_pim, P_pim, true, true ) ) {
      std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
    }

    //Momentum (n CDS)
    TVector3 P_n;
    if(UseDecayVtx){
      P_n = tmp_mom*((Pos_CDH-vtx_dis).Unit());
    }else{
      P_n = tmp_mom*((Pos_CDH-vtx_react).Unit());
    }
    TVector3 P_n_beam = tmp_mom_beam*((Pos_CDH-vtx_beam).Unit());
    TVector3 P_n_vtx[3];
    P_n_vtx[0] = tmp_mom_vtx[0]*((Pos_CDH-vtx_pip).Unit());
    P_n_vtx[1] = tmp_mom_vtx[1]*((Pos_CDH-vtx_pim).Unit());
    //P_n_vtx[0] = tmp_mom_vtx[0]*((Pos_CDH-CA_pip_pippim).Unit());
    //P_n_vtx[1] = tmp_mom_vtx[1]*((Pos_CDH-CA_pim_pippim).Unit());
    P_n_vtx[2] = tmp_mom_vtx[2]*((Pos_CDH-vtx_react).Unit());

    LVec_pim.SetVectM( P_pim, piMass );
    LVec_pip.SetVectM( P_pip, piMass );
    LVec_n.SetVectM( P_n, nMass );//CDS n
    LVec_n_beam.SetVectM( P_n_beam, nMass );//CDS n
    for(int ivtx=0;ivtx<3;ivtx++){
      LVec_n_vtx[ivtx].SetVectM( P_n_vtx[ivtx], nMass );//CDS n
    }
      
    /*
       std::cout << "NeutralBetaCDH " << NeutralBetaCDH << std::endl;
       std::cout << "NeutralBetaCDH_vtx[0] " << NeutralBetaCDH_vtx[0] << std::endl;
       std::cout << "NeutralBetaCDH_vtx[1] " << NeutralBetaCDH_vtx[1] << std::endl;
       std::cout << "mom:" << tmp_mom << std::endl;
       std::cout << "pip-n mom:" << tmp_mom_vtx[0] << std::endl;
       std::cout << "pim-n mom:" << tmp_mom_vtx[1] << std::endl;
       std::cout << "Lvec0 " << LVec_n_vtx[0].P() << std::endl;
       std::cout << "Lvec1 " << LVec_n_vtx[1].P() << std::endl;
    */

    const double cdhphi = Pos_CDH.Phi();
    const double pimphi = P_pim.Phi();
    const double pipphi = P_pip.Phi();
    Tools::Fill1D(Form("npimangle"),pimphi-cdhphi);
    Tools::Fill1D(Form("npipangle"),pipphi-cdhphi);

    const double mm_mass   = (LVec_target+LVec_beam-LVec_pim-LVec_pip-LVec_n).M();
    double mm_mass_vtx[3];
    for(int ivtx=0;ivtx<3;ivtx++){
      mm_mass_vtx[ivtx]= (LVec_target+LVec_beam_vtx[ivtx]-LVec_pim-LVec_pip-LVec_n_vtx[ivtx]).M();       
    }

    const TVector3 P_missn = (LVec_target+LVec_beam-LVec_pim-LVec_pip-LVec_n).Vect();
    TVector3 P_missn_vtx[3];
    for(int ivtx=0;ivtx<3;ivtx++){
      P_missn_vtx[ivtx] = (LVec_target+LVec_beam_vtx[ivtx]-LVec_pim-LVec_pip-LVec_n_vtx[ivtx]).Vect();  
    }

    LVec_nmiss.SetVectM( P_missn, nMass );
    for(int ivtx=0;ivtx<3;ivtx++){
      LVec_nmiss_vtx[ivtx].SetVectM( P_missn_vtx[ivtx], nMass );
    }

    if(Verbosity>10)std::cerr<<"  missing mass = "<<mm_mass<<std::endl;

    const TVector3 boost = (LVec_target+LVec_beam).BoostVector();
    TLorentzVector LVec_nmiss_CM = LVec_nmiss;
    TLorentzVector LVec_beam_CM = LVec_beam;
    LVec_nmiss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    //cos in CM frame
    const double cos_n = LVec_nmiss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_nmiss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    if(Verbosity>10)std::cerr<<"  missing mom | cos_CM = "<<cos_n<<std::endl;


    // + + + + + + + + + + + + + //
    //  fill histograms & tree   //
    // + + + + + + + + + + + + + //
    kf_flag = -1;
    bool K0rejectFlag=false;
    bool MissNFlag=false;
    bool SigmaPFlag=false;
    bool SigmaMFlag=false;
    bool NBetaOK=false;
    bool NdEOK=false;

    Tools::Fill2D( Form("dE_betainv"), 1./NeutralBetaCDH, ncdhhit->emean() );
    Tools::Fill2D( Form("MMom_MMass"), mm_mass, P_missn.Mag() );


    Tools::Fill2D(Form("Vtx_ZX_nofid"),vtxpip_mean.Z(),vtxpip_mean.X());
    Tools::Fill2D(Form("Vtx_ZY_nofid"),vtxpip_mean.Z(),vtxpip_mean.Y());
    Tools::Fill2D(Form("Vtx_XY_nofid"),vtxpip_mean.X(),vtxpip_mean.Y());
    Tools::Fill2D(Form("Vtx_ZX_nofid"),vtxpim_mean.Z(),vtxpim_mean.X());
    Tools::Fill2D(Form("Vtx_ZY_nofid"),vtxpim_mean.Z(),vtxpim_mean.Y());
    Tools::Fill2D(Form("Vtx_XY_nofid"),vtxpim_mean.X(),vtxpim_mean.Y());
    //Fiducial cuts OK
    if( (!IsVtxDoubleCheck && (GeomTools::GetID(vtx_react)==CID_Fiducial)) || 
        ( IsVtxDoubleCheck && 
          (GeomTools::GetID(vtxpim_mean)==CID_Fiducial) &&
          (GeomTools::GetID(vtxpip_mean)==CID_Fiducial)))  {


      for( int i=0; i<cdsMan->nCDH(); i++ ) {
        if((cdsMan->CDH(i)->CheckRange()) && (cdsMan->CDH(i)->ctmean()<cdscuts::tdc_cdh_max)){
          Tools::Fill2D(Form("dE_CDHtime_pippimn"), cdsMan->CDH(i)->ctmean(), cdsMan->CDH(i)->emean());
        }
      }

      Tools::Fill2D(Form("Vtx_ZX_primfid"),vtx_react.Z(),vtx_react.X());
      Tools::Fill2D(Form("Vtx_ZY_primfid"),vtx_react.Z(),vtx_react.Y());
      Tools::Fill2D(Form("Vtx_XY_primfid"),vtx_react.X(),vtx_react.Y());
      Tools::Fill2D(Form("Vtx_ZX_fid"),vtxpip_mean.Z(),vtxpip_mean.X());
      Tools::Fill2D(Form("Vtx_ZY_fid"),vtxpip_mean.Z(),vtxpip_mean.Y());
      Tools::Fill2D(Form("Vtx_XY_fid"),vtxpip_mean.X(),vtxpip_mean.Y());
      Tools::Fill2D(Form("Vtx_ZX_fid"),vtxpim_mean.Z(),vtxpim_mean.X());
      Tools::Fill2D(Form("Vtx_ZY_fid"),vtxpim_mean.Z(),vtxpim_mean.Y());
      Tools::Fill2D(Form("Vtx_XY_fid"),vtxpim_mean.X(),vtxpim_mean.Y());


      Tools::Fill2D( Form("NMomCDHtime"),ncdhhit->ctmean()-ctmT0-beamtof,P_n.Mag()); 
      Tools::Fill2D( Form("NMomCDHtime%d",CDHSeg),ncdhhit->ctmean()-ctmT0-beamtof,P_n.Mag()); 
      Tools::Fill2D( Form("NeutraltimeEnergy"),ncdhhit->ctmean()-ctmT0-beamtof,ncdhhit->emean());
      Tools::Fill2D( Form("CDHzNeutraltime"),Pos_CDH.z(),ncdhhit->ctmean()-ctmT0-beamtof);
      Tools::Fill2D( Form("CDH%dzNeutraltime",CDHSeg),Pos_CDH.z(),ncdhhit->ctmean()-ctmT0-beamtof);
      Tools::Fill2D( Form("dE_betainv_fid"), 1./NeutralBetaCDH, ncdhhit->emean() );
      Tools::Fill2D( Form("MMom_MMass_fid"), mm_mass, P_missn.Mag() );

      if(NeutralBetaCDH<anacuts::beta_MAX) NBetaOK=true;

      if(NBetaOK) {
        Tools::Fill2D( Form("dE_betainv_fid_beta"), 1./NeutralBetaCDH, ncdhhit->emean() );
        Tools::Fill2D( Form("MMom_MMass_fid_beta"), mm_mass, P_missn.Mag() );
      }
      if(anacuts::dE_MIN<ncdhhit->emean()) NdEOK=true;

      if( NBetaOK && NdEOK ) {
        Tools::Fill2D( Form("dE_betainv_fid_beta_dE"), 1./NeutralBetaCDH, ncdhhit->emean() );
        Tools::Fill2D( Form("MMom_MMass_fid_beta_dE"), mm_mass, P_missn.Mag() );
        Tools::Fill1D( Form("IMpipi_dE"), (LVec_pim+LVec_pip).M() );
        Tools::Fill2D( Form("IMpipi_NMom_dE"),P_n.Mag(), (LVec_pim+LVec_pip).M());


        //added Jul.28th,2019
        //purpose 
        for( int it=0; it<trackMan->nGoodTrack(); it++ ) {
          CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );

          double mom = track->Momentum();//charge X momentum
          double tof = 999.;//TOF of CDH-T0 (slewing corrected)
          double mass2 = -999.;
          double correctedtof=0;//CDH-T0 (corrected by energy loss)
          double beta_calc=0;
          for( int icdh=0; icdh<track->nCDHHit(); icdh++ ) {
            HodoscopeLikeHit *cdhhit = track->CDHHit( cdsMan, icdh );
            double tmptof = cdhhit->ctmean()-ctmT0;
            if( tmptof<tof || tof==999. ) {
              tof = tmptof;
            }
          }
          if( !TrackTools::FindMass2( track, bpctrack, tof, LVec_beam.Vect().Mag(),
                Beam_Kaon, beta_calc, mass2, correctedtof ) ) {

          }
          //const int pid = TrackTools::PIDcorr_wide(mom,mass2);
          const int pid = TrackTools::PIDcorr(mom,mass2);
          track->SetPID( pid );
          Tools::Fill2D( "PID_CDS_beta_select2", 1./beta_calc, mom );
          Tools::Fill2D( "PID_CDS_select2", mass2, mom );
          if(pid == CDS_PiMinus){
            Tools::Fill2D("PID_CDS_PIM_beta_select2",1./beta_calc,mom);
            Tools::Fill2D("PID_CDS_PIM_select2",mass2,mom);
          }else if(pid == CDS_PiPlus){
            Tools::Fill2D("PID_CDS_PIP_beta_select2",1./beta_calc,mom);
            Tools::Fill2D("PID_CDS_PIP_select2",mass2,mom);
          }
          else if(pid == CDS_Proton) Tools::Fill2D("PID_CDS_Proton_select2",mass2,mom);
          else if(pid == CDS_Kaon) Tools::Fill2D("PID_CDS_Kaon_select2",mass2,mom);
        }//for it
      }

      if( ((LVec_pim+LVec_pip).M()<anacuts::pipi_MIN || anacuts::pipi_MAX<(LVec_pim+LVec_pip).M())) K0rejectFlag=true;

      //missing mass neutron ID
      if( anacuts::neutron_MIN<mm_mass && mm_mass<anacuts::neutron_MAX ) MissNFlag=true;

      //Sigma+ production in CDS
      if( (anacuts::Sigmap_MIN<(LVec_n+LVec_pip).M() && (LVec_n+LVec_pip).M()<anacuts::Sigmap_MAX)) SigmaPFlag=true;

      //Sigma- production in CDS
      if( (anacuts::Sigmam_MIN<(LVec_n+LVec_pim).M() && (LVec_n+LVec_pim).M()<anacuts::Sigmam_MAX)) SigmaMFlag=true;

      if( NBetaOK && NdEOK ) {
        //K0rejection
        if(K0rejectFlag) {
          // K0 rejection
          Tools::Fill2D( Form("dE_betainv_fid_beta_dE_woK0"), 1./NeutralBetaCDH, ncdhhit->emean() );
          Tools::Fill2D( Form("MMom_MMass_fid_beta_dE_woK0"), mm_mass, P_missn.Mag() );
          Tools::Fill2D( Form("CDHseg_MMass_fid_beta_dE_woK0"),ncdhhit->seg(),mm_mass);
          Tools::Fill2D( Form("CDHz_MMass_fid_beta_dE_woK0"),Pos_CDH.z(),mm_mass);
          //Tools::Fill2D(Form("zVTX_MMass_fid_beta_dE_woK0"),vtx_react.z(),mm_mass);

          Tools::Fill2D( Form("IMnpim_IMnpip_dE_woK0"), (LVec_n+LVec_pip).M(), (LVec_n+LVec_pim).M() );
          Tools::Fill2D( Form("dE_MMom_fid_beta_woK0"), P_missn.Mag(), ncdhhit->emean() );
          Tools::Fill2D( Form("dE_MMass_fid_beta_woK0"), mm_mass, ncdhhit->emean() );
        } else {
          // K0 selection
          Tools::Fill2D( Form("dE_betainv_fid_beta_dE_wK0"), 1./NeutralBetaCDH, ncdhhit->emean() );
          Tools::Fill2D( Form("MMom_MMass_fid_beta_dE_wK0"), mm_mass, P_missn.Mag() );
        }//K0 rejection

        if(K0rejectFlag && MissNFlag) {
          Tools::Fill2D( Form("dE_betainv_fid_beta_dE_woK0_n"), 1./NeutralBetaCDH, ncdhhit->emean() );
          Tools::Fill2D( Form("MMom_MMass_fid_beta_dE_woK0_n"), mm_mass, P_missn.Mag() );
          Tools::Fill2D( Form("NMom_NMom_fid_beta_dE_woK0_n"), P_n.Mag(), P_missn.Mag() );
          Tools::Fill2D( Form("IMnpim_IMnpip_dE_woK0_n"), (LVec_n+LVec_pip).M(), (LVec_n+LVec_pim).M() );
          Tools::Fill2D( Form("CDHz_IMnpip_fid_beta_dE_woK0_n"),Pos_CDH.z(),(LVec_n+LVec_pip).M());
          Tools::Fill2D( Form("CDHz_IMnpim_fid_beta_dE_woK0_n"),Pos_CDH.z(),(LVec_n+LVec_pim).M());
        }

        if(K0rejectFlag && (SigmaPFlag || SigmaMFlag)) {
          Tools::Fill2D(Form("MMom_MMass_fid_beta_dE_woK0_wSid"),mm_mass, P_missn.Mag());
          for( int i=0; i<cdsMan->nCDH(); i++ ){
            if((cdsMan->CDH(i)->CheckRange()) && (cdsMan->CDH(i)->ctmean()<cdscuts::tdc_cdh_max)){
              Tools::Fill2D( Form("CDHdE_woK0_wSid"),cdsMan->CDH(i)->seg(),cdsMan->CDH(i)->emean());
            }
          }
        }

        if(K0rejectFlag && SigmaMFlag) {
          for( int i=0; i<cdsMan->nCDH(); i++ ){
            if((cdsMan->CDH(i)->CheckRange()) && (cdsMan->CDH(i)->ctmean()<cdscuts::tdc_cdh_max)){
              Tools::Fill2D( Form("CDHdE_woK0_wSmid"),cdsMan->CDH(i)->seg(),cdsMan->CDH(i)->emean());
            }
          }
        }


        if( MissNFlag ) {
          Tools::Fill1D( Form("IMnpipi_n"), (LVec_n+LVec_pim+LVec_pip).M() );
          Tools::Fill2D( Form("MMnmiss_IMnpipi_n"),(LVec_n+LVec_pim+LVec_pip).M(), LVec_nmiss.M());
        }

        if( MissNFlag && (SigmaPFlag || SigmaMFlag)) {
          Tools::Fill1D( Form("IMnpipi_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M() );
          Tools::Fill2D( Form("MMnmiss_IMnpipi_wSid_n"),(LVec_n+LVec_pim+LVec_pip).M(), LVec_nmiss.M());
        }

        if( MissNFlag && K0rejectFlag && (SigmaPFlag || SigmaMFlag)) {
          Tools::Fill2D( Form("IMmnpim_IMmnpip_woK0_wSid_n"), (LVec_nmiss+LVec_pip).M(), (LVec_nmiss+LVec_pim).M() );
          Tools::Fill2D( Form("MMnpip_MMnpim_woK0_wSid_n"), (LVec_target+LVec_beam-LVec_pim-LVec_n).M(),
              (LVec_target+LVec_beam-LVec_pip-LVec_n).M() );

          Tools::Fill1D( Form("IMnpipi_woK0_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M() );

          Tools::Fill2D( Form("dE_IMnpipi_woK0_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M(), ncdhhit->emean());
          //cos theta
          Tools::Fill2D( Form("Cosn_IMnpipi_woK0_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M(), cos_n );
          //
          Tools::Fill2D( Form("MMnmiss_IMnpipi_woK0_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M(), LVec_nmiss.M());
          Tools::Fill2D( Form("nmom_IMnpipi_woK0_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M(), LVec_n.P());

          //momentum transfer
          Tools::Fill2D( Form("q_IMnpipi_woK0_wSid_n"),(LVec_n+LVec_pim+LVec_pip).M(), (LVec_beam.Vect()-LVec_nmiss.Vect()).Mag());

          Tools::Fill1D( Form("DCA_pip_SigmaPM"),dcapipvtx);
          Tools::Fill1D( Form("DCA_pim_SigmaPM"),dcapimvtx);
          Tools::Fill1D( Form("DCA_pippim_SigmaPM"),dcapippim);
          //vertex position from pi+/pi-
            
            for( int i=0; i<cdsMan->nCDH(); i++ ){
              if((cdsMan->CDH(i)->CheckRange()) && (cdsMan->CDH(i)->ctmean()<cdscuts::tdc_cdh_max)){
                Tools::Fill2D( Form("CDHdE_woK0_wSid_n"),cdsMan->CDH(i)->seg(),cdsMan->CDH(i)->emean());
              }
            }
          }//MissNFlag && K0rejectFlag && (SigmaPFlag || SigmaMFlag)
          
          if(K0rejectFlag && MissNFlag && SigmaPFlag) {
            Tools::Fill1D( Form("DCA_pip_SigmaP"),dcapipvtx);
            Tools::Fill1D( Form("DCA_pim_SigmaP"),dcapimvtx);
            Tools::Fill1D( Form("DCA_pippim_SigmaP"),dcapippim);
          }

          if(K0rejectFlag && MissNFlag && SigmaMFlag) {
            Tools::Fill1D( Form("DCA_pip_SigmaM"),dcapipvtx);
            Tools::Fill1D( Form("DCA_pim_SigmaM"),dcapimvtx);
            Tools::Fill1D( Form("DCA_pippim_SigmaM"),dcapippim);
          }

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
            //*** definition of fit particles in Cartesian coordinates ***//
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
            //*** definition of constraints ***//
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
            //*** definition of the fitter ***//
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

            //*** perform the fit ***//
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
            //*** copy fit results ***//
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

            //** fill tree **//
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
        } // if( dE_MIN<ncdhhit->emean() )

        mom_beam   = LVec_beam;   // 4-momentum(beam)
        mom_beam_Sp = LVec_beam_vtx[0];
        mom_beam_Sm = LVec_beam_vtx[1];
        mom_beam_K0 = LVec_beam_vtx[2];
        mom_target = LVec_target; // 4-momentum(target)
        mom_pip = LVec_pip;        // 4-momentum(pi+)
        mom_pim = LVec_pim;        // 4-momentum(pi-)
        mom_n = LVec_n;            // 4-momentum(neutron)
        mom_n_beam = LVec_n_beam;            // 4-momentum(neutron)
        mom_n_Sp = LVec_n_vtx[0];
        mom_n_Sm = LVec_n_vtx[1];
        mom_n_K0 = LVec_n_vtx[2];
        dE = ncdhhit->emean();
        neutralseg = (ncdhhit->seg());
        nhitOutCDC = nCDCforVeto;
        ForwardCharge = Util::IsForwardCharge(blMan);
        // beta is already filled
        vtx_reaction = vtx_react; // vertex(reaction)
        vtx_displaced = vtx_dis; // vertex(reaction)
        vtx_pip_beam = vtx_beam_wpip;
        vtx_pim_beam = vtx_beam_wpim;
        vtx_pip_cdc = vtx_pip;
        vtx_pim_cdc = vtx_pim;
        CA_pip = CA_pip_pippim;
        CA_pim = CA_pim_pippim;
        CDH_Pos = Pos_CDH;//neutron hit pos
        CDH_Pos_pim = pim_cdhprojected;
        CDH_Pos_pip = pip_cdhprojected;
        run_num   = confMan->GetRunNumber(); // run number
        event_num = Event_Number;            // event number
        block_num = Block_Event_Number;      // block number

        if(Verbosity>10)std::cout<<"%%% npippim event: Event_Number, Block_Event_Number, CDC_Event_Number = "
                                   <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
        if(header->IsTrig(Trig_KCDH3)){
          Tools::H1(Form("Trig_npippim"),1,10,0,10);
          rtFile3->cd();
          npippimTree->Fill();
          nFill_npippim++;
          //** fill tree **//
        }else if(header->IsTrig(Trig_KCDH2f)){
          Tools::H1(Form("Trig_npippim"),2,10,0,10);
        }else if(header->IsTrig(Trig_Kf)){
          Tools::H1(Form("Trig_npippim"),3,10,0,10);
        }else{
          Tools::H1(Form("Trig_npippim"),4,10,0,10);
        }
        rtFile->cd();
      } // if( GeomTools::GetID(vtx_react)==CID_Fiducial )
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
  std::cout<<" nAbort_KCDH3trg      = "<<nAbort_KCDH3trg<<std::endl;
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
  std::cout<<" nAbort_CDCInner3Lay  = "<<nAbort_CDCInner3Lay<<std::endl;
  std::cout<<" nAbort_pipi          = "<<nAbort_pipi<<std::endl;
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
   for( int i=0; i<40; i++ ){
     Tools::newTH1F(Form("SCA%d",i), 2000000, -0.5, 1999999.5);
   }


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
  NeutralBetaCDH_beam=-9999.;
  NeutralBetaCDH_vtx[0]=-9999.;
  NeutralBetaCDH_vtx[1]=-9999.;
  NeutralBetaCDH_vtx[2]=-9999.;
  tofpim=-9999.;
  tofpip=-9999.;
  tofn=-9999.;

  dE=-9999.;
  neutralseg=-1;
  nhitOutCDC=-1;
  ForwardCharge=-1;
  vtx_reaction.SetXYZ(-9999.,-9999.,-9999.);
  vtx_displaced.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pip_beam.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pim_beam.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pip_cdc.SetXYZ(-9999.,-9999.,-9999.);
  vtx_pim_cdc.SetXYZ(-9999.,-9999.,-9999.);
  CA_pip.SetXYZ(-9999.,-9999.,-9999.);
  CA_pim.SetXYZ(-9999.,-9999.,-9999.);
  CDH_Pos.SetXYZ(-9999.,-9999.,-9999.);
  CDH_Pos_pim.SetXYZ(-9999.,-9999.,-9999.);
  CDH_Pos_pip.SetXYZ(-9999.,-9999.,-9999.);

  return;
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
EventTemp *EventAlloc::EventAllocator()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
