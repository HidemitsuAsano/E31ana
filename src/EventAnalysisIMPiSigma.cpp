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


#define KFDEBUG 0 // verbose level of the KinFitter
// 0: quiet, 1: print result, 2: print iterations, 3: print also matrices

//-- set run# --//

const int MaxTreeSize = 1000000000;

const bool DoCDCRetiming = false;
const int Verbosity = 0;
const bool DoKinFit = true;
const bool AddQAplots = true;

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
  bool EveSelectCDHMul();
  bool EveSelectBeamline();
  double AnaBeamSpec(ConfMan *confMan);
  bool IsForwardCharge();
  bool IsForwardNeutron();
  void InitKinFitMatrix();
  int CDSChargedAna(LocalTrack *bpc,
                    const TLorentzVector beam,
                    std::vector <int> &cdhseg,
                    std::vector <int> &pimid,
                    std::vector <int> &pipid,
                    std::vector <int> &kmid,
                    std::vector <int> &protonid);

  TFile *rtFile;// histograms
  TFile *rtFile2;// condensed file
  TFile *rtFile3;// pi+,pi-,n event
  TFile *cdcFile;
  TTree *cdcTree;
  TTree *evTree;
  TTree *npippimTree;

  Float_t tlogprob1;
  Float_t tlogprob2;

  Float_t twoimprob1;
  Float_t twoimprob2;

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

  //** counters for filling **//
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
  int nAbort_pipi;
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
  double kfSpmode_status; // status of kinematical refit -> details can be found in this code
  double kfSpmode_pvalue; // p-value of kinematical refit
  TLorentzVector kfSmmode_mom_beam;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
  TLorentzVector kfSmmode_mom_pip;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
  TLorentzVector kfSmmode_mom_pim;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
  TLorentzVector kfSmmode_mom_n;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
  double kfSmmode_chi2;   // chi2 of kinematical refit
  double kfSmmode_NDF;    // NDF of kinematical refit
  double kfSmmode_status; // status of kinematical refit -> details can be found in this code
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

  if(DoCDCRetiming) std::cout << " Yes" << std::endl;
  else              std::cout << "  No" << std::endl;

  std::cout << " Kinematic fit ? " ;
  if(DoKinFit) std::cout << " Yes" << std::endl;
  else         std::cout << " No"  << std::endl;


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
  header_CDC = 0;
  trackMan_CDC = 0;
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
    std::cerr << "!!!!" << std::endl;
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
    std::cerr << "!!!!" << std::endl;
    return;
  }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ) {
    std::cerr << "!!!!" << std::endl;
    return;
  }
  blMan = new BeamLineHitMan();
  if( blMan==NULL ) {
    std::cerr << "!!!!" << std::endl;
    return;
  }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan==NULL ) {
    std::cerr << "!!!!" << std::endl;
    return;
  }
  scaMan = new ScalerMan();
  if( scaMan==NULL ) {
    std::cerr << "!!!!" << std::endl;
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


  if( Event_Number%1000==1 ) {
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

  Tools::Fill1D( Form("EventCheck"), 1 );

  //CDH-hits cut
  if(!EveSelectCDHMul()) return true;
  Tools::Fill1D( Form("EventCheck"), 2 );

  //** # of good CDS tracks cut **//
  if( nGoodTrack!=cdscuts::cds_ngoodtrack ) { //require pi+,pi-
    Clear( nAbort_nGoodTrack );
    return true;
  }
  Tools::Fill1D( Form("EventCheck"), 3 );

  //beam line analysis and event selection
  if(!EveSelectBeamline()) return true;
  Tools::Fill1D( Form("EventCheck"), 4 );

  //PID beam
  if(PIDBeam!=Beam_Kaon) return true;
  Tools::Fill1D( Form("EventCheck"), 5 );

  //BLC1-D5-BLC2 analysis and chi2 selection
  double beammom = AnaBeamSpec(confMan);
  if(beammom <-100) return true;//chi2 cut
  Tools::Fill1D( Form("EventCheck"), 6 );

  //** beam momentum calculation **//
  TVector3 Pp_target;
  Pp_target.SetXYZ( 0, 0, 0 );
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
  double z1 = 0, z2 = 20;
  LocalTrack *bpctrack = bltrackMan->trackBPC(bpcGoodTrackID);
  bpctrack->XYPosatZ( z1, x1, y1 );
  bpctrack->XYPosatZ( z2, x2, y2 );
  TVector3 lp;
  lp.SetXYZ( x1, y1, z1 );
  TVector3 ls;
  ls.SetXYZ( x2-x1, y2-y1, z2-z1);
  ls = ls.Unit();
  TVector3 Pp_beam = beammom*ls;

  LVec_beambf.SetVectM( Pp_beam, kpMass );
  LVec_beam = LVec_beambf;
  TVector3 boost = (LVec_target+LVec_beam).BoostVector();
  LVec_beambfCM = LVec_beam;
  LVec_targetCM = LVec_target;
  //LVec_targetPCM = LVec_targetP;
  //boost to CM frame
  LVec_beambfCM.Boost( -1*boost );
  LVec_targetCM.Boost( -1*boost );
  //LVec_targetPCM.Boost( -1*boost );


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
  int nIDedTrack = CDSChargedAna(bpctrack,LVec_beam,vCDHseg,pim_ID,pip_ID,km_ID,p_ID);
  if(nIDedTrack<0) return true;
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
      trackMan->nGoodTrack()==cdscuts::cds_ngoodtrack && !IsForwardCharge() ) {
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

    if(Verbosity) {
      if( NeutralCDHseg.size()!=1 ) {
        std::cerr<<" CDH neutral hit is not 1 :: "<<NeutralCDHseg.size()<<std::endl;
      } else {
        Tools::Fill1D( Form("CDHNeutralseg"),NeutralCDHseg.at(0));
      }
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
    }

    //** isolation cut **//
    int flag_isolation = 0;
    for( int ineuseg=0; ineuseg<(int)NeutralCDHseg.size(); ineuseg++ ) {
      for( int ihit=0; ihit<(int)CDHhit_list.size(); ihit++ ) {
        if( NeutralCDHseg[ineuseg]-CDHhit_list[ihit] ) {
          Tools::Fill1D( Form("diff_CDH"), NeutralCDHseg[ineuseg]-CDHhit_list[ihit] );
        }
        //CDH has 36 segments. Requiring there is no hits on neighboring segments.
        if( abs(NeutralCDHseg[ineuseg]-CDHhit_list[ihit])==1 || abs(NeutralCDHseg[ineuseg]-CDHhit_list[ihit])==35 )
          flag_isolation++;
      }
    }

    if( flag_isolation ) {
      if(Verbosity) std::cerr<<"CDH hit candidate is NOT isolated !!!"<<std::endl;
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

    // charge veto using CDC
    TVector3 Pos_CDH;
    confMan->GetGeomMapManager()->GetPos( CID_CDH, ncdhhit->seg(), Pos_CDH );
    if(Verbosity) {
      std::cerr<<"CDH candidate = "<<ncdhhit->seg()<<" -> "<<Pos_CDH.Phi()/TwoPi*360.<<" deg"<<std::endl;
    }
    const double PhiMin = -15.0/360.*TwoPi; // rad
    const double PhiMax =  15.0/360.*TwoPi; // rad
    if(Verbosity) {
      std::cerr<<"Min/Max = "<<PhiMin/TwoPi*360.<<"/"<<PhiMax/TwoPi*360.<<" deg"<<std::endl;
    }
    int nCDCforVeto = 0;
    for( int ilr=14; ilr<16; ilr++ ) { // charge veto using layer 15, 16
      for( int icdchit=0; icdchit<cdsMan->nCDC(ilr); icdchit++ ) {
        CDCHit *cdc=cdsMan->CDC(ilr,icdchit);
        TVector3 Pos_CDC = cdc->wpos();
        Pos_CDC.SetZ(0); // only xy pos is used
        double angle = Pos_CDC.Angle(Pos_CDH); // rad
        if(Verbosity) {
          std::cerr<<"CDC "<<ilr<<" "<<icdchit<<" "<<cdc->wire()<<" -> "<<Pos_CDC.Phi()/TwoPi*360.
                   <<" deg :: diff = "<<angle/TwoPi*360<<" deg"<<std::endl;
        }
        Tools::Fill1D( Form("diff_CDH_CDC"), angle/TwoPi*360 );
        if( PhiMin<angle && angle<PhiMax ) nCDCforVeto++;
      }//icdchit
    }//ilr
    if(Verbosity) {
      std::cerr<<"# of CDC hits for nCDH candidate = "<<nCDCforVeto<<std::endl;
    }
    Pos_CDH.SetZ(-1*ncdhhit->hitpos()); // (-1*) is correct in data analysis [20170926]


    //** neutral particle in CDH **//
    if( !nCDCforVeto ) {
      if(NeutralCDHseg.size()!=1) {
        std::cout << "L." << __LINE__ << " # of seg for neutral hits " << NeutralCDHseg.size() << std::endl;
      } else {
        if(AddQAplots) Tools::Fill1D(Form("CDHNeutralSeg"),NeutralCDHseg.at(0));
      }


      CDSTrack *track_pip = trackMan->Track( pip_ID.at(0) ); // only 1 track
      CDSTrack *track_pim = trackMan->Track( pim_ID.at(0) ); // only 1 track

      //deuteron target
      TVector3 vtx_react;
      TVector3 vtx_dis;
      TVector3 vtx_beam_wpip;//vertex(beam-pip) on beam
      TVector3 vtx_pip;//vertex(beam-pip) on beam
      TVector3 vtx_beam;
      track_pip->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wpip, vtx_pip );
      TVector3 vtx_beam_wpim;//vertex(beam-pim) on beam
      TVector3 vtx_pim;//vertex(beam-pim) on beam
      track_pim->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_beam_wpim, vtx_pim );

      double dcapipvtx =  (vtx_pip-vtx_beam_wpip).Mag();
      double dcapimvtx =  (vtx_pim-vtx_beam_wpim).Mag();


      //reaction vertex is determined from beam and nearest vtx
      if(dcapipvtx < dcapimvtx) {
        //follows sakuma/sada's way , avg. of scattered particle ana beam particle [20180829]
        vtx_react = 0.5*(vtx_pip+vtx_beam_wpip);
        if(cdscuts::useclosestpi) vtx_dis  = vtx_pip;
        else              vtx_dis  = vtx_pim;
        vtx_beam = vtx_beam_wpip;
      } else {
        vtx_react = 0.5*(vtx_pim+vtx_beam_wpim);
        if(cdscuts::useclosestpi) vtx_dis = vtx_pim;
        else             vtx_dis = vtx_pip;
        vtx_beam = vtx_beam_wpim;
      }
      //vertex position from pi+/pi-
      TVector3 vtx1,vtx2;
      bool vtx_flag=TrackTools::Calc2HelixVertex(track_pip, track_pim, vtx1, vtx2);
      double dcapippim=-9999.;
      if(vtx_flag) dcapippim = (vtx2-vtx1).Mag();



      //** beam kaon tof **//
      TVector3 Pos_T0;
      confMan->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
      double beamtof=0;
      double momout=0;
      double z_pos = Pos_T0.Z();;
      //dE correction of beam
      ELossTools::CalcElossBeamTGeo( bpctrack->GetPosatZ(z_pos), vtx_react,
                                     LVec_beambf.Vect().Mag(), kpMass, momout, beamtof );
      LVec_beam.SetVectM( momout*LVec_beambf.Vect().Unit(), kpMass );
      double ntof = ncdhhit->ctmean()-ctmT0-beamtof;
      double nlen = (Pos_CDH-vtx_react).Mag();
      if(Verbosity>10) std::cout << "L." << __LINE__ << " flight length " << nlen << std::endl;
      NeutralBetaCDH = nlen/ntof/(Const*100.);
      double tmp_mom = NeutralBetaCDH<1. ? nMass*NeutralBetaCDH/sqrt(1.-NeutralBetaCDH*NeutralBetaCDH) : 0;
      if(Verbosity) {
        std::cerr<<"L. " << __LINE__ ;
        std::cerr<<" NeutralBetaCDH = "<<NeutralBetaCDH<<" mom_n = "<<tmp_mom<<std::endl; //" "<<1/sqrt(1+nMass*nMass)<<std::endl;
      }
      //** reconstructoin of missing neutorn **//
      TVector3 P_pim; // Momentum(pi-)
      TVector3 P_pip; // Momentum(pi+)
      TVector3 P_n;   // Momentum(n) (CDS)

      TLorentzVector LVec_pim; // 4-Momentum(pi-)
      TLorentzVector LVec_pip; // 4-Momentum(pi+)
      TLorentzVector LVec_n;   // 4-Momentum(n) (CDS)
      TLorentzVector LVec_nmiss; // 4-Momentum(n_miss)

      if( !track_pip->GetMomentum( vtx_pip, P_pip, true, true ) ) {
        std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }

      if( !track_pim->GetMomentum( vtx_pim, P_pim, true, true ) ) {
        std::cerr<<"L." << __LINE__ << " !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }
      P_n = tmp_mom*((Pos_CDH-vtx_react).Unit());

      LVec_pim.SetVectM( P_pim, piMass );
      LVec_pip.SetVectM( P_pip, piMass );
      LVec_n.SetVectM(   P_n,   nMass );//CDS n

      double cdhphi = Pos_CDH.Phi();
      double pimphi = P_pim.Phi();
      double pipphi = P_pip.Phi();
      Tools::Fill1D(Form("npimangle"),pimphi-cdhphi);
      Tools::Fill1D(Form("npipangle"),pipphi-cdhphi);

      double mm_mass   = (LVec_target+LVec_beam-LVec_pim-LVec_pip-LVec_n).M();
      TVector3 P_missn = (LVec_target+LVec_beam-LVec_pim-LVec_pip-LVec_n).Vect();
      LVec_nmiss.SetVectM( P_missn, nMass );
      if(Verbosity>10)std::cerr<<"  missing mass = "<<mm_mass<<std::endl;

      TVector3 boost = (LVec_target+LVec_beam).BoostVector();
      TLorentzVector LVec_nmiss_CM = LVec_nmiss;
      TLorentzVector LVec_beam_CM = LVec_beam;
      LVec_nmiss_CM.Boost(-boost);
      LVec_beam_CM.Boost(-boost);
      //cos in CM frame
      double cos_n = LVec_nmiss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_nmiss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
      if(Verbosity>10)std::cerr<<"  missing mom | cos_CM = "<<cos_n<<std::endl;


      //** + + + + + + + + + + + + + **//
      //**  fill histograms & tree   **//
      //** + + + + + + + + + + + + + **//
      kf_flag = -1;
      bool K0rejectFlag=false;
      bool MissNFlag=false;
      bool SigmaPFlag=false;
      bool SigmaMFlag=false;
      bool NBetaOK=false;
      bool NdEOK=false;

      Tools::Fill2D( Form("dE_betainv"), 1./NeutralBetaCDH, ncdhhit->emean() );
      Tools::Fill2D( Form("MMom_MMass"), mm_mass, P_missn.Mag() );

      if(AddQAplots) {
        Tools::Fill2D(Form("Vtx_ZX_nofid"),vtx_beam.Z(),vtx_beam.X());
        Tools::Fill2D(Form("Vtx_ZY_nofid"),vtx_beam.Z(),vtx_beam.Y());
        Tools::Fill2D(Form("Vtx_XY_nofid"),vtx_beam.X(),vtx_beam.Y());
      }
      //Fiducial cuts OK
      if( GeomTools::GetID(vtx_beam)==CID_Fiducial ) {
        if(AddQAplots) {
          Tools::Fill2D(Form("Vtx_ZX_fid"),vtx_beam.Z(),vtx_beam.X());
          Tools::Fill2D(Form("Vtx_ZY_fid"),vtx_beam.Z(),vtx_beam.Y());
          Tools::Fill2D(Form("Vtx_XY_fid"),vtx_beam.X(),vtx_beam.Y());
          Tools::Fill2D(Form("NeutraltimeEnergy"),ncdhhit->ctmean()-ctmT0-beamtof,ncdhhit->emean());
        }
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

            Tools::Fill2D( Form("IMnpim_IMnpip_dE_woK0"), (LVec_n+LVec_pip).M(), (LVec_n+LVec_pim).M() );
            Tools::Fill2D( Form("dE_MMom_fid_beta_woK0"), P_missn.Mag(), ncdhhit->emean() );
            Tools::Fill2D( Form("dE_MMass_fid_beta_woK0"), mm_mass, ncdhhit->emean() );
          } else {
            // K0 selection
            Tools::Fill2D( Form("dE_betainv_fid_beta_dE_wK0"), 1./NeutralBetaCDH, ncdhhit->emean() );
            Tools::Fill2D( Form("MMom_MMass_fid_beta_dE_wK0"), mm_mass, P_missn.Mag() );
          }

          if(K0rejectFlag && MissNFlag) {
            Tools::Fill2D( Form("dE_betainv_fid_beta_dE_woK0_n"), 1./NeutralBetaCDH, ncdhhit->emean() );
            Tools::Fill2D( Form("MMom_MMass_fid_beta_dE_woK0_n"), mm_mass, P_missn.Mag() );
            Tools::Fill2D( Form("NMom_NMom_fid_beta_dE_woK0_n"), P_n.Mag(), P_missn.Mag() );
            Tools::Fill2D( Form("IMnpim_IMnpip_dE_woK0_n"), (LVec_n+LVec_pip).M(), (LVec_n+LVec_pim).M() );
          }

          if(K0rejectFlag && (SigmaPFlag || SigmaMFlag)) {
            Tools::Fill2D(Form("MMom_MMass_fid_beta_dE_woK0_wSid"),mm_mass, P_missn.Mag());
          }

          if( MissNFlag ) {
            Tools::Fill1D( Form("IMnpipi_n"), (LVec_n+LVec_pim+LVec_pip).M() );
            Tools::Fill2D( Form("MMnmiss_IMnpipi_n"),(LVec_n+LVec_pim+LVec_pip).M(), P_missn.Mag());
          }

          if( MissNFlag && (SigmaPFlag || SigmaMFlag)) {
            Tools::Fill1D( Form("IMnpipi_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M() );
            Tools::Fill2D( Form("MMnmiss_IMnpipi_wSid_n"),(LVec_n+LVec_pim+LVec_pip).M(), P_missn.Mag());
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
            Tools::Fill2D( Form("MMnmiss_IMnpipi_woK0_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M(), P_missn.Mag());
            Tools::Fill2D( Form("nmom_IMnpipi_woK0_wSid_n"), (LVec_n+LVec_pim+LVec_pip).M(), LVec_n.P());

            //momentum transfer
            Tools::Fill2D( Form("q_IMnpipi_woK0_wSid_n"),(LVec_n+LVec_pim+LVec_pip).M(), (LVec_beam.Vect()-LVec_nmiss.Vect()).Mag());

            Tools::Fill1D( Form("DCA_pip"), dcapipvtx );
            Tools::Fill1D( Form("DCA_pim"), dcapimvtx );
            Tools::Fill1D( Form("DCA_pippim"), dcapippim);
          }

          if(K0rejectFlag && MissNFlag && SigmaPFlag) {
            Tools::Fill1D( Form("DCA_pip_SigmaP"),dcapipvtx);
            Tools::Fill1D( Form("DCA_pim_SigmaP"),dcapimvtx);
          }

          if(K0rejectFlag && MissNFlag && SigmaMFlag) {
            Tools::Fill1D( Form("DCA_pip_SigmaM"),dcapipvtx);
            Tools::Fill1D( Form("DCA_pim_SigmaM"),dcapimvtx);
          }

          if(DoKinFit) {
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
            // %%% Kinematical Fit using KinFitter %%% //
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
            //--- set TLorentzVector ---//
            // beam_K(K+), pi-/+, Sigma+/-,  missn, n from S, pi+/- from S
            //  = 1) TLorentzVector LVec_beam, LVec_pim, (LVec_n+LVec_pip), LVec_nmiss, LVec_n, LVec_pip = for pi- Sigma+
            TLorentzVector TL_meas_Spmode[kin::npart]; // measured
            TL_meas_Spmode[kin::kmbeam] = LVec_beam;
            TL_meas_Spmode[kin::pim_g1] = LVec_pim;
            TL_meas_Spmode[kin::Sp] = (LVec_n+LVec_pip);
            TL_meas_Spmode[kin::nmiss] = LVec_nmiss;
            TL_meas_Spmode[kin::ncds] = LVec_n;
            TL_meas_Spmode[kin::pip_g2] = LVec_pip;
            //  = 2) TLorentzVector LVec_beam, LVec_pip, (LVec_n+LVec_pim), LVec_nmiss, LVec_n, LVec_pim = for pi+ Sigma-
            TLorentzVector TL_meas_Smmode[kin::npart]; // measured
            TL_meas_Smmode[kin::kmbeam] = LVec_beam;
            TL_meas_Smmode[kin::pip_g1] = LVec_pip;
            TL_meas_Smmode[kin::Sm] = (LVec_n+LVec_pim);
            TL_meas_Smmode[kin::nmiss] = LVec_nmiss;
            TL_meas_Smmode[kin::ncds] = LVec_n;
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
            int PDG_Spmode[kin::npart] = {-321, -211, 3222, 2112, 2112,  211}; // pi-Sigma+
            //                    K-    pi+  S-     n     n     pi-
            int PDG_Smmode[kin::npart] = {-321,  211, 3112, 2112, 2112, -211}; // pi+Sigma-


            //--- KinFitter :: initialization ---//
            //  = 1) TLorentzVector LVec_beam, LVec_pim, (LVec_n+LVec_pip), LVec_nmiss, LVec_n, LVec_pip = for pi- Sigma+
            //  = 2) TLorentzVector LVec_beam, LVec_pip, (LVec_n+LVec_pim), LVec_nmiss, LVec_n, LVec_pim = for pi+ Sigma-
            //*** definition of fit particles in cartesian coordinates ***//
            TString str_particle_Spmode[kin::npart] = {"LVec_beam", "LVec_pim", "LVec_Sp", "LVec_mn", "LVec_n", "LVec_pip"};
            TString str_particle_Smmode[kin::npart] = {"LVec_beam", "LVec_pip", "LVec_Sm", "LVec_mn", "LVec_n", "LVec_pim"};
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
            TString str_constEp_Spmode[4]  = {"Px", "Py", "Pz", "E"};
            TString str_constEp_Smmode[4]  = {"Px", "Py", "Pz", "E"};
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
        mom_target = LVec_target; // 4-momentum(target)
        mom_pip = LVec_pip;        // 4-momentum(pi+)
        mom_pim = LVec_pim;        // 4-momentum(pi-)
        mom_n = LVec_n;            // 4-momentum(neutron)
        dE = ncdhhit->emean();
        // beta is already filled
        vtx_reaction = vtx_react; // vertex(reaction)
        run_num   = confMan->GetRunNumber(); // run number
        event_num = Event_Number;            // event number
        block_num = Block_Event_Number;      // block number

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
  // geneneral informantion **//
  Tools::newTH1F( Form("Time"), 3000, -0.5, 2999.5 );
  Tools::newTH1F( Form("nGoodTrack"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("EventCheck"), 20, 0, 20 );
  // Event Reduction Check memo
  // 1.  # of all events w/o cosmic trigger
  // 2.  after CDH cuts
  // 3.  after the # of CDS good tracks cut
  // 4.  after beam line analysis
  // 5.  after beam PID (kaon)
  // 6.  after BLC1-D5-BLC2 analysis
  // 7.  filled if CDS chi2 is bad
  // 8.  filled if a CDC track have no CDH hits
  // 9.  filled if CDC tracks share a CDH segments
  // 10. filled if FindMass2 is failed
  // 11. filled if FindMass2 is failed after Retiming
  // 12. filled if Energy loss calculation failed
  // 13. filled if CDSChargedAna is OK.
  // 14. filled if CDH isolation is OK.


  Tools::newTH1F( Form("Scaler"), 41, -0.5, 40.5 );

  // CDC and CDH information from CDC-trackig file **//
  Tools::newTH1F( Form("mul_CDH"),Form("CDH multiplicity"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_CDH_assoc"), 11, -0.5, 10.5 );
  if(AddQAplots) {
    Tools::newTH2F( Form("CDHtime"),36,0.5,36.5,4000,0,200);
    Tools::newTH2F( Form("NeutraltimeEnergy"),500,0,100,200,0,50);
    //Tools::newTH3F( Form("CDHtimeEnergy"),36,0.5,36.5,4000,0,200,100,0,100);
    Tools::newTH1F( Form("CDHNeutralSeg"),36, 0.5, 36.5);
  }
  Tools::newTH1F( Form("npimangle"),628, 0, 2*3.14);
  Tools::newTH1F( Form("npipangle"),628, 0, 2*3.14);

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
  Tools::newTH2F( Form("PID_CDS_beta"), 2000, 0, 10., 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS"), 1000, -0.6, 5, 1000, -1.2, 1.2 );
  Tools::newTH1F( Form("ntrack_CDS"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_pi_plus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_proton"), 6, -0.5, 5.5 );
  //Tools::newTH1F( Form("ntrack_deuteron"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_pi_minus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_K_minus"), 6, -0.5, 5.5 );

  if(AddQAplots) {
    Tools::newTH2F( Form("Vtx_ZX"),1000,-25,25,500,-12.5,12.5);
    Tools::newTH2F( Form("Vtx_ZY"),1000,-25,25,500,-12.5,12.5);
    Tools::newTH2F( Form("Vtx_XY"),500,-12.5,12.5,500,-12.5,12.5);

    Tools::newTH2F( Form("Vtx_ZX_nofid"),1000,-25,25,500,-12.5,12.5);
    Tools::newTH2F( Form("Vtx_ZY_nofid"),1000,-25,25,500,-12.5,12.5);
    Tools::newTH2F( Form("Vtx_XY_nofid"),500,-12.5,12.5,500,-12.5,12.5);
    Tools::newTH2F( Form("Vtx_ZX_fid"),1000,-25,25,500,-12.5,12.5);
    Tools::newTH2F( Form("Vtx_ZY_fid"),1000,-25,25,500,-12.5,12.5);
    Tools::newTH2F( Form("Vtx_XY_fid"),500,-12.5,12.5,500,-12.5,12.5);
  }
  //** forward counters **//
  Tools::newTH1F( Form("mul_BVC"), 9, -0.5, 8.5 );
  Tools::newTH1F( Form("mul_CVC"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_PC"), 11, -0.5, 10.5 );

  //** pi+ pi- X event **//
  Tools::newTH1F( Form("diff_CDH"), 73, -36.5, 36.5 );
  Tools::newTH1F( Form("diff_CDH_CDC"), 181, 0, 181 );
  Tools::newTH2F( Form("dE_betainv"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta_dE"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta_dE_woK0"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta_dE_wK0"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta_dE_woK0_n"), 200, 0, 10, 200, 0, 50);
  Tools::newTH2F( Form("dE_MMom_fid_beta_woK0"), 100, 0, 1.5, 200, 0, 50);
  Tools::newTH2F( Form("dE_MMass_fid_beta_woK0"), 140, 0.4, 1.8, 200, 0, 50);

  Tools::newTH2F( Form("MMom_MMass"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE_woK0"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE_woK0_wSid"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE_wK0"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE_woK0_n"), 140, 0.4, 1.8, 100, 0, 1.5 );



  Tools::newTH1F( Form("IMpipi_dE"), 200, 0.4, 0.6 );
  Tools::newTH2F( Form("IMpipi_NMom_dE"),100,0, 1.5, 200, 0.4,0.6);

  Tools::newTH2F( Form("NMom_NMom_fid_beta_dE_woK0_n"), 100, 0, 1.5, 100, 0, 1.5 );
  Tools::newTH2F( Form("IMnpim_IMnpip_dE_woK0"), 140, 1, 1.7, 140, 1, 1.7 );
  Tools::newTH2F( Form("IMnpim_IMnpip_dE_woK0_n"), 140, 1, 1.7, 140, 1, 1.7 );
  Tools::newTH2F( Form("IMmnpim_IMmnpip_woK0_wSid_n"), 70, 1, 1.7, 70, 1, 1.7 );
  Tools::newTH2F( Form("MMnpip_MMnpim_woK0_wSid_n"), 70, 1, 1.7, 70, 1, 1.7 );
  Tools::newTH2F( Form("Cosn_IMnpipi_woK0_wSid_n"), 100, 1, 2, 50, -1, 1 );
  Tools::newTH1F( Form("IMnpipi_n"), 100, 1, 2 );
  Tools::newTH1F( Form("IMnpipi_wSid_n"), 100, 1, 2 );
  Tools::newTH1F( Form("IMnpipi_woK0_wSid_n"), 100, 1, 2 );
  Tools::newTH2F( Form("dE_IMnpipi_woK0_wSid_n"), 100, 1, 2, 200, 0, 50);
  Tools::newTH2F( Form("MMnmiss_IMnpipi_n"),100,1,2,100,0,1.5);
  Tools::newTH2F( Form("MMnmiss_IMnpipi_wSid_n"),100,1,2,100,0,1.5);
  Tools::newTH2F( Form("MMnmiss_IMnpipi_woK0_wSid_n"),100,1,2,100,0,1.5);
  Tools::newTH2F( Form("nmom_IMnpipi_woK0_wSid_n"),100,1,2,100,0,1.0);
  Tools::newTH2F( Form("q_IMnpipi_woK0_wSid_n"),100,1,2,300,0,1.5);
  Tools::newTH1F( Form("DCA_pip"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pim"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pip_SigmaP"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pim_SigmaP"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pip_SigmaM"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pim_SigmaM"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pippim"), 500, 0, 5);

  Tools::newTH2F( Form("KFchi2_vs"),100,0,100,100,0,100);
  Tools::newTH1F( Form("KF_decision"), 2, -0.5, 1.5 );//TODO implement
  return;
}



bool EventAnalysis::EveSelectCDHMul()
{
  //** # of CDH-hits cut **//
  int nCDH = 0;
  for( int i=0; i<cdsMan->nCDH(); i++ ) {
    if(AddQAplots)Tools::Fill2D(Form("CDHtime"),cdsMan->CDH(i)->seg(),cdsMan->CDH(i)->ctmean());
    //if( cdsMan->CDH(i)->CheckRange() ) nCDH++; //** only requirement of TDC **//
    if( cdsMan->CDH(i)->CheckRange() && cdsMan->CDH(i)->ctmean()<cdscuts::tdc_cdh_max ) {
      nCDH++;
    }
  }
  Tools::Fill1D( Form("mul_CDH"), nCDH );
  if(Verbosity>100) {
    std::cout << __FILE__ <<  " L." << __LINE__ ;
    std::cout << " nCDH:" << nCDH << std::endl;
  }
  if( nCDH != cdscuts::cdhmulti  ) {
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
  for( int i=0; i<blMan->nBHD(); i++ ) {
    if( blMan->BHD(i)->CheckRange() ) nBHD++;
  }
  int nT0 = 0;
  for( int i=0; i<blMan->nT0(); i++ ) {
    if( blMan->T0(i)->CheckRange() ) nT0++;
  }
  Tools::Fill1D( Form("mul_BHD"), nBHD );
  Tools::Fill1D( Form("mul_T0"),  nT0 );

  //** T0 = 1hit selection **//
  if( nT0!=1 ) { //!! sada-D p.72 !!//
    Clear( nAbort_nT0 );
    return false;
  }

  //** Beam PID by T0-BHD TOF **//
  TVector3 vtxT0;
  for( int i=0; i<blMan->nT0(); i++ ) {
    if( blMan->T0(i)->CheckRange() ) {
      ctmT0 = blMan->T0(i)->ctmean();
      confMan->GetGeomMapManager()->GetGPos(CID_T0, blMan->T0(i)->seg(), vtxT0);
    }
  }
  double ctmBHD;
  PIDBeam = -1; // 0:pi 1:K 3:else
  for( int i=0; i<blMan->nBHD(); i++ ) {
    if( blMan->BHD(i)->CheckRange() ) {
      ctmBHD = blMan->BHD(i)->ctmean();
      double tofBHDT0 = ctmT0-ctmBHD;
      Tools::Fill1D( Form("tof_T0BHD"), tofBHDT0 );
      if( header->kaon() && blcuts::beam_tof_k_min <tofBHDT0 && tofBHDT0<blcuts::beam_tof_k_max  )
        PIDBeam = Beam_Kaon;
      else if( header->pion() && blcuts::beam_tof_pi_min<tofBHDT0 && tofBHDT0<blcuts::beam_tof_pi_max )
        PIDBeam = Beam_Pion;
    }
  }
  Tools::Fill1D(Form("PID_beam"), PIDBeam);
  if( PIDBeam== -1  ) { //** unidentified particle is discarded (other than pi/K) **//
    Clear( nAbort_pid_beam );
    return false;
  }

  unsigned int nblc1 = 0;
  unsigned int nblc2 = 0;
  blc1GoodTrackID = -1;
  blc2GoodTrackID = -1;
  //** timing selection of BLC1/BLC2 **//
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ) {
    LocalTrack *blc1 = bltrackMan->trackBLC1(i);
    Tools::Fill1D( Form("tracktime_BLC1"), blc1->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BLC1"), blc1->chi2all() );
    if( blc1->CheckRange(blcuts::blc1_time_window_min, blcuts::blc1_time_window_max)) {
      nblc1++;
      if( blc1->CheckRange(blcuts::blc1_time_min, blcuts::blc1_time_max) &&
          bltrackMan->trackBLC1(i)->chi2all()<blcuts::blc1_chi2_max ) blc1GoodTrackID = i;
    }
  }
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ) {
    LocalTrack *blc2 = bltrackMan->trackBLC2(i);
    Tools::Fill1D( Form("tracktime_BLC2"), blc2->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BLC2"), blc2->chi2all() );
    if( blc2->CheckRange(blcuts::blc2_time_window_min, blcuts::blc2_time_window_max) ) {
      nblc2++;
      if( blc2->CheckRange(blcuts::blc2_time_min, blcuts::blc2_time_max) &&
          bltrackMan->trackBLC2(i)->chi2all()<blcuts::blc2_chi2_max ) blc2GoodTrackID = i;
    }
  }
  Tools::Fill1D( Form("ntrack_BLC1"), nblc1 );
  Tools::Fill1D( Form("ntrack_BLC2"), nblc2 );

  //** single track selection in each BLC **//
  if( !(nblc1==1 && blc1GoodTrackID!=-1 && nblc2==1 && blc2GoodTrackID!=-1) ) { //** multi-good-beams event is discarded  **//
    Clear( nAbort_singleBLtrack );
    return false;
  }

  //** BPC track selection **//
  int nbpc = 0;
  bpcGoodTrackID = -1;
  double chibpc = 999;
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ) {
    Tools::Fill1D( Form("tracktime_BPC"), bltrackMan->trackBPC(i)->GetTrackTime() );
    if( bltrackMan->trackBPC(i)->CheckRange(blcuts::bpc_time_window_min,blcuts::bpc_time_window_max) ) {
      nbpc++;
      bpcGoodTrackID = i;
      chibpc = bltrackMan->trackBPC(i)->chi2all();
      Tools::Fill1D( Form("trackchi2_BPC"),chibpc);
    }
  }

  Tools::Fill1D( Form("ntrack_BPC"), nbpc );

  if( nbpc!=1 ) {
    Clear( nAbort_nbpc );
    return false;
  }

  LocalTrack *bpctrack = bltrackMan->trackBPC(bpcGoodTrackID);
  if( !(bpctrack->CheckRange(blcuts::bpc_time_min,blcuts::bpc_time_max))
      || bpctrack->chi2all()>blcuts::bpc_chi2_max) {
    Clear( nAbort_bpctrack );
    return false;
  }

  //** vertex calculation by CDS goodtrack and BPC tracks**/
  for( int it1=0; it1<trackMan->nGoodTrack(); it1++ ) {
    trackMan->CalcVertex_beam( trackMan->GoodTrackID(it1), bltrackMan, confMan );
  }

  //** BLC2-BPC track matching **//
  bool fblc2bpc = false;
  for( int ii=0; ii<bltrackMan->ntrackBLC2(); ii++ ) {
    if( ii!=blc2GoodTrackID ) continue;
    LocalTrack *blc2 = bltrackMan->trackBLC2(ii);
    double xblc2bpc[2]= {0,0};
    double yblc2bpc[2]= {0,0};
    double xpos[2]= {0,0};
    double ypos[2]= {0,0};

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
  if( !fblc2bpc ) {
    Clear( nAbort_fblc2bpc );
    return false;
  }

  return true;
}

//BLC1-D5-BLC2 momentum analysis
double EventAnalysis::AnaBeamSpec(ConfMan *confMan) {

  BeamSpectrometer *beamsp = new BeamSpectrometer( confMan );
  LocalTrack *blc1 = bltrackMan->trackBLC1(blc1GoodTrackID);
  LocalTrack *blc2 = bltrackMan->trackBLC2(blc2GoodTrackID);
  beamsp->TMinuitFit( blc1, blc2, confMan );
  double beammom = beamsp->mom();
  double bchi = beamsp->chisquare();
  delete beamsp;

  Tools::Fill1D( Form("trackchi2_beam"), bchi );
  if( bchi>blcuts::d5_chi2_max ) {
    Clear( nAbort_flagbmom );
    return -9999.;
  }
  Tools::Fill1D( Form("momentum_beam"), beammom );

  return beammom;
}

//

int EventAnalysis::CDSChargedAna(LocalTrack* bpctrack,
                                 const TLorentzVector LVec_beam,
                                 std::vector <int> &cdhseg,
                                 std::vector <int> &pimid,
                                 std::vector <int> &pipid,
                                 std::vector <int> &kmid,
                                 std::vector <int> &protonid)
{
  int CDHseg=-1;
  bool chi2OK = true;
  bool CDHseg1hitOK = true;
  bool CDHsharecheckOK = true;
  bool FindMass2OK1 = true;
  bool FindMass2OK2 = true;
  bool EnergyLossOK = true;

  for( int it=0; it<trackMan->nGoodTrack(); it++ ) {
    CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );

    // chi2 cut can be applied in CDSTrackingMan with MaxChi in CDSFittingParam_posi.param
    Tools::Fill1D( Form("trackchi2_CDC"), track->Chi() );
    if( track->Chi()>cdscuts::cds_chi2_max ) {
      chi2OK = false;
    }
    //asano memo
    //checking if there is CDH hits at the projected position of the track
    if( !track->CDHFlag() ) {
      CDHseg1hitOK = false;
    }
    double mom = track->Momentum();
    TVector3 vtxbline, vtxbhelix; //,vtxb;
    track->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtxbline, vtxbhelix );
    track->SetPID(-1);
    Tools::Fill2D(Form("Vtx_ZX"),vtxbline.Z(),vtxbline.X());
    Tools::Fill2D(Form("Vtx_ZY"),vtxbline.Z(),vtxbline.Y());
    Tools::Fill2D(Form("Vtx_XY"),vtxbline.X(),vtxbline.Y());
    //vtxb = (vtxbline+vtxbhelix)*0.5;//not used, so far

    double tof = 999.;
    double mass2 = -999.;
    int nCDHass = track->nCDHHit();
    Tools::Fill1D(Form("mul_CDH_assoc"),nCDHass);
    if(nCDHass>1) {
      CDHseg1hitOK = false;
    }
    for( int icdh=0; icdh<track->nCDHHit(); icdh++ ) {
      HodoscopeLikeHit *cdhhit = track->CDHHit( cdsMan, icdh );
      double tmptof = cdhhit->ctmean()-ctmT0;
      if( tmptof<tof || tof==999. ) {
        tof = tmptof;
        CDHseg = cdhhit->seg();
      }
    }

    //asano memo
    //check if the segment has been already used in the analysis.
    bool CDHflag = true;
    for( int icdhseg=0; icdhseg<(int)cdhseg.size(); icdhseg++) {
      if( CDHseg==cdhseg[icdhseg] ) CDHflag = false;
    }
    if( !CDHflag ) {
      Tools::Fill1D( Form("EventCheck"), 9 );//put here to get event based info.
      std::cerr << " The CDH segments is used by an another track !!! " << std::endl;
      CDHsharecheckOK = false;
      continue;//go to next CDStrack
    }
    cdhseg.push_back( CDHseg );

    // calculation of beta and squared-mass //
    double tmptof=0;
    double beta_calc=0;
    if( !TrackTools::FindMass2( track, bpctrack, tof, LVec_beam.Vect().Mag(),
                                PIDBeam, beta_calc, mass2, tmptof ) ) {
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      FindMass2OK1 = false;
      continue;
    }

    // Retiming of CDC track by CDH info. //
    if(DoCDCRetiming) {
      track->Retiming( cdsMan, confMan, beta_calc, true );
      //asano memo
      //why need this ?
      for( int m=0; m<5; m++ ) {
        track->HelixFitting( cdsMan ); // 5 times iteration //
      }
      track->Calc( confMan );
    }


    // finalize PID //
    if( !TrackTools::FindMass2( track, bpctrack, tof, LVec_beam.Vect().Mag(),
                                PIDBeam, beta_calc, mass2, tmptof ) ) { // not FindMass2C() [20170622] //
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      FindMass2OK2 = false;
      continue;
    }

    //RUN68 inoue's param.
    int pid = -1;
    pid = TrackTools::PIDcorr_wide(mom,mass2);

    track->SetPID( pid );
    Tools::Fill2D( "PID_CDS_beta", 1/beta_calc, mom );
    Tools::Fill2D( "PID_CDS", mass2, mom );


    // Energy loss calculation //
    double tmpl=0;
    TVector3 vtx_beam, vtx_cds;
    if( !track->CalcVertexTimeLength( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), track->Mass(),
                                      vtx_beam, vtx_cds, tmptof, tmpl, true ) ) {
      std::cerr<<" !!! failure in energy loss calculation [CalcVertexTimeLength()] !!! "<<std::endl;
      EnergyLossOK = false;
      continue;
    }

    if( pid==CDS_PiMinus ) {
      pimid.push_back( trackMan->GoodTrackID(it) );
    } else if( pid==CDS_PiPlus ) {
      pipid.push_back( trackMan->GoodTrackID(it) );
    } else if( pid==CDS_Proton ) {
      protonid.push_back( trackMan->GoodTrackID(it) );
    } else if( pid==CDS_Kaon ) {
      kmid.push_back( trackMan->GoodTrackID(it) );
    }
  } // for( int it=0; it<trackMan->nGoodTrack(); it++ )
  // end of PID (except for neutron) //

  if(!chi2OK) {
    Tools::Fill1D( Form("EventCheck"), 7 );
    return -1;
  }

  if(!CDHseg1hitOK) {
    Tools::Fill1D( Form("EventCheck"), 8 );
    return -1;
  }

  if(!CDHsharecheckOK) {
    Tools::Fill1D( Form("EventCheck"), 9 );
    return -1;
  }

  if(!FindMass2OK1) {
    Tools::Fill1D( Form("EventCheck"), 10 );
    return -1;
  }

  if(!FindMass2OK2) {
    Tools::Fill1D( Form("EventCheck"), 11 );
    return -1;
  }

  if(!EnergyLossOK) {
    Tools::Fill1D( Form("EventCheck"), 12 );//put here to get event based info.
    return -1;
  }

  return pimid.size()+pipid.size()+protonid.size()+kmid.size();
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



bool EventAnalysis::IsForwardCharge() {
  //** charge veto with BVC, CVC (TOF=CVC), & PC **//
  int nBVC = 0;
  int nCVC = 0;
  int nPC  = 0;
  for( int i=0; i<blMan->nBVC(); i++ ) {
    if( blMan->BVC(i)->CheckRange() ) nBVC++;
  }
  for( int i=0; i<blMan->nTOF(); i++ ) {
    if( blMan->TOF(i)->CheckRange() ) nCVC++;
  }
  for( int i=0; i<blMan->nPC(); i++ ) {
    if( blMan->PC(i)->CheckRange() ) nPC++;
  }
  Tools::Fill1D( Form("mul_BVC"), nBVC );
  Tools::Fill1D( Form("mul_CVC"), nCVC );
  Tools::Fill1D( Form("mul_PC"),  nPC );
  if( nBVC || nCVC || nPC ) return true;
  else return false;
}

//TODO implement ?
//purpose : check if there is a neutron at NC
bool EventAnalysis::IsForwardNeutron() {

  if(IsForwardCharge()) return false;
  /*
  int nNeutron=0;
  for(int i=0; i<blMan->nNC();i++){
    HodoscopeLikeHit* hit = blMan->NC(i);
    if(hit->CheckRange()){
      //pNC is defined in Particle.h
      int hitseg = hit->seg();
      TVector3 hitpos;
      hitpos.SetXYZ(hit->pos().X(),hit->hitpos(),hit->pos().Z());
      double meantime = hit->ctmean();
      double meanenergy = hit->cmean();


    }
  }

  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particle->nNC(); it++){
    pNC* nc = particle->nc(it);
    nc->CalcMom(beam,vertex);
    double seg = (nc->seg()-1)%16+1;
    double lay = (nc->seg()-1)/16+1;
    FillHist("FWDN_Overbeta",1.0/nc->beta());
    FillHist("FWDN_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("NCHitPosition_XY",nc->hitpos().X(),nc->hitpos().Y());
    FillHist("NCHitSegment",seg,lay);
    if(nc->energy()>8.0 && (time>9998||time>nc->time())){
      time = nc->time();
      fnc = it;
      FillHist("FWDN_Overbeta_withdECut",1.0/nc->beta());
      FillHist("FWDN_OverbetavsEnergy_withdECut",1.0/nc->beta(),nc->energy());
    }
  }
  if(fnc==-1) return false;
  pNC*  nc =0;
  if(fnc!=-1) nc = particle->nc(fnc);
  FillHist("FWDN_Overbeta_Selected",1.0/nc->beta());
  if(nc->pid()!=F_Neutron){ return false; }

  */
  return true;
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
  PIDBeam=-1;
  blc1GoodTrackID=-1;
  blc2GoodTrackID=-1;
  bpcGoodTrackID=-1;
  ctmT0 = 0;
  //CDS
  NeutralBetaCDH=9999.;
  dE=-9999.;

  return;
}


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
EventTemp *EventAlloc::EventAllocator()
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
