//----------------------------------------------------------------//
// ===== EventAnalysis_pipipnn_sakuma.cpp =====
// *** originated from EventAnalysis_Lpn_sakuma.cpp ***
// pi+pi-pnn reconstruction program using CDC-tracking-fil
// generated from EventAnalysis_Tracking.cpp (evtracking).
//----------------------------------------------------------------//
//  exe-file: evpipipnn_sakuma
//      [./evlpn_sakuma $(ConfFile) $(runnumber) $(OutFile) $(InFile) $(CDC-tracking-file)]
//  input : raw data, conf-file, & CDC-tracking-file
//  output: when $(OutFile) is "tmp.root", the following 3 files are generated.
//     "tmp.root":      histogram file
//     "tmp_CDC.root": condensed event file from the CDC-tracking-file with p/p/pi selection
//           <- including class EventHeader and CDSTrackingMan
//     "tmp_pipipnn.root": basic information of pi+pi-pnn event is listed up in TTree
//----------------------------------------------------------------//
//  updated by F.S, 2017 8/1
//----------------------------------------------------------------//

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
#include <KinFitter/TKinFitter.h>
#include <KinFitter/TFitParticlePxPyPz.h>
#include <KinFitter/TFitConstraintM.h>
#include <KinFitter/TFitConstraintEp.h>
#include <Math/ProbFuncMathCore.h>

#define VTX_SIGMA 1 // 1: Sigma reconstruction,  vertex = K- & p
                    // 0: Lambda reconstruction, vertex = K- & pi+

#define DEBUG 0

#define KFDEBUG 0 // verbose level of the KinFitter
// 0: quiet, 1: print result, 2: print iterations, 3: print also matrices

//-- set run# --//
//const std::string GRUN = "49c";
const std::string GRUN = "65";
//-------------//
//****** Sada analysis *******//
// data: /w/e15/data/Run49c/
// conf: conf/Run49c/analyzer_sada.conf
// CDC:  /w/e15/sada/root/Run49c/evtracking_fin2_*.root
//----------------------------//
//----------------------------//
//***** Yamaga analysis *****//
// data: /w/e15/data/Run65/
// conf: conf/Run65/analyzer_run65.conf
// CDC:  /w/e15/yamaga/cdctracking/Run65/run65_*_v07.root
//---------------------------//

//-----------------------------------------//
//--- covariance matrices for KinFitter ---//
//-----------------------------------------//
// ### obtained from (p_meas[j]-p_gene[j])*(p_meas[k]-p_gene[k])
// ###  using G4-data with TH1F(Form("cov_%d_%d_%d", i, j, k), 100, -cov_MAX, cov_MAX);
//   evaluated using "Air" Dora MC
// 1) TLorentzVector L3_beam, L_pim, (L_n+L_pip), L_p, L_nmiss, L_n, L_pip = for pi- Sigma+
const double covVal1[7][16] = {
    { 2.03675e-05, 0, 0, 0,
      0, 1.72317e-05, 0, 0,
      0, 0, 3.63879e-06, 0,
      0, 0, 0, 3.13463e-06 },
    { 9.61619e-06, 0, 0, 0,
      0, 9.71896e-06, 0, 0,
      0, 0, 1.95299e-05, 0,
      0, 0, 0, 1.23025e-05 },
    { 0.000344907, 0, 0, 0,
      0, 0.000344784, 0, 0,
      0, 0, 6.33015e-05, 0,
      0, 0, 0, 4.95287e-05 },
    { 9.63118e-05, 0, 0, 0,
      0, 9.49883e-05, 0, 0,
      0, 0, 0.000100833, 0,
      0, 0, 0, 6.21131e-05 },
    { 0.000624282, 0, 0, 0,
      0, 0.000664143, 0, 0,
      0, 0, 0.00028286, 0,
      0, 0, 0, 0.000121639 },
    { 0.000253317, 0, 0, 0,
      0, 0.000234635, 0, 0,
      0, 0, 1.96443e-05, 0,
      0, 0, 0, 2.26387e-05 },
    { 9.82579e-06, 0, 0, 0,
      0, 1.00648e-05, 0, 0,
      0, 0, 1.50993e-05, 0,
      0, 0, 0, 7.68311e-06 }
};
// 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
double covVal2[7][16];
// 3) TLorentzVector L3_beam, L_pip, (L_p+L_pim), L_n, L_nmiss, L_p, L_pim = for pi+ Lambda
const double covVal3[7][16] = { // temporaly same as piSigma - have to be evaluated
    { 1.92133e-05, 0, 0, 0,
      0, 1.84496e-05, 0, 0,
      0, 0, 3.86757e-06, 0,
      0, 0, 0, 3.2775e-06 },
    { 1.23926e-05, 0, 0, 0,
      0, 1.19661e-05, 0, 0,
      0, 0, 2.22957e-05, 0,
      0, 0, 0, 2.20691e-05 },
    { 0.000121839, 0, 0, 0,
      0, 0.000130653, 0, 0,
      0, 0, 9.66432e-05, 0,
      0, 0, 0, 4.40921e-05 },
    { 0.000345995, 0, 0, 0,
      0, 0.00038387, 0, 0,
      0, 0, 2.89356e-05, 0,
      0, 0, 0, 3.34355e-05 },
    { 0.000640011, 0, 0, 0,
      0, 0.000674321, 0, 0,
      0, 0, 0.000210562, 0,
      0, 0, 0, 0.000122066 },
    { 0.00014628, 0, 0, 0,
      0, 0.000132961, 0, 0,
      0, 0, 7.55614e-05, 0,
      0, 0, 0, 3.42671e-05 },
    { 2.10616e-05, 0, 0, 0,
      0, 1.72002e-05, 0, 0,
      0, 0, 1.12709e-05, 0,
      0, 0, 0, 2.56098e-06 }
};
TMatrixD *covZero;
TMatrixD *covParticle1[7];
TMatrixD *covParticle2[7];
TMatrixD *covParticle3[7];

const int MaxTreeSize = 1000000000;

const double TDC_CDH_MAX = 25; // ns

const int PDG1[7] = {321, -211, 3222, 2212, 2112, 2112,  211}; // pi-Sigma+
const int PDG2[7] = {321,  211, 3112, 2212, 2112, 2112, -211}; // pi+Sigma-
const int PDG3[7] = {321,  211, 3122, 2112, 2112, 2212, -211}; // pi+Lambda

// ================================================== //
//                                                    //
// class EventAnalysis delivered from class EventTemp //
//                                                    //
// ================================================== //
class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TFile *rtFile2;
  TFile *rtFile3;
  TFile *cdcFile;
  TTree *cdcTree;
  TTree *evTree;
  TTree *pipipnnTree;

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

  //** counters for filing **//
  int nFill_pppi;
  int nFill_pipipnn;
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

  //** cut parameters **//
  //@@ only different cut parameters btw Sada & Yamaga are listed @@//
  double PARA_tof_K_MIN;
  double PARA_tof_K_MAX;
  double PARA_tof_pi_MIN;
  double PARA_tof_pi_MAX;
  double PARA_time_range_BLC;
  double PARA_blc2bpc_dx_MIN;
  double PARA_blc2bpc_dx_MAX;
  double PARA_blc2bpc_dy_MIN;
  double PARA_blc2bpc_dy_MAX;
  double PARA_blc2bpc_dxdz_MIN;
  double PARA_blc2bpc_dxdz_MAX;
  double PARA_blc2bpc_dydz_MIN;
  double PARA_blc2bpc_dydz_MAX;
  double PARA_lnL_MAX;

  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector mom_beam;   // 4-momentum(beam)
  TLorentzVector mom_target; // 4-momentum(target)
  TLorentzVector mom_pip;    // 4-momentum(pi+)
  TLorentzVector mom_pim;    // 4-momentum(pi-)
  TLorentzVector mom_p;      // 4-momentum(proton)
  TLorentzVector mom_n;      // 4-momentum(neutron)
  double beta; // veracity of neutral particle on CDH
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
  TLorentzVector kf3mom_beam;   // 4-momentum(beam) after kinematical refit for pi+ Lambda
  TLorentzVector kf3mom_pip;    // 4-momentum(pi+) after kinematical refit for pi+ Lambda
  TLorentzVector kf3mom_pim;    // 4-momentum(pi-) after kinematical refit for pi+ Lambda
  TLorentzVector kf3mom_p;      // 4-momentum(proton) after kinematical refit for pi+ Lambda
  TLorentzVector kf3mom_n;      // 4-momentum(neutron) after kinematical refit for pi+ Lambda
  double kf3_chi2;   // chi2 of kinematical refit
  double kf3_NDF;    // NDF of kinematical refit
  double kf3_status; // status of kinematical refit -> details can be found in this code
  double kf3_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc
  //= = = = pipipnn final-sample tree = = = =//

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


//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
void EventAnalysis::Initialize( ConfMan *conf )
//**--**--**--**--**--**--**--**--**--**--**--**--**--**//
{
  std::cout << " *** Enter EventAnalysis::Initialize " << std::endl;

  INIT = true;
  //  spillinit = spillfini = -1;
  for( int i=0; i<40; i++ ){
    SCAOVERFLOW[i] = false;
    scaend[i] = 0;
  }

  //** Conf file open **//
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

  //** output file 2 : condensed event with p/p/pi selection **//
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

  //** output file 3 : pipipnn final-sample tree **//
  std::string outfile3 = confMan->GetOutFileName();
  outfile3.insert( outfile3.size()-5, "_pipipnn" );
  std::cout<<"pipipnn file "<<outfile3<<std::endl;
  rtFile3 = new TFile( outfile3.c_str(), "recreate" );
  rtFile3->cd();
  pipipnnTree = new TTree( "EventTree", "EventTree" );
  pipipnnTree->Branch( "mom_beam",   &mom_beam );
  pipipnnTree->Branch( "mom_target", &mom_target );
  pipipnnTree->Branch( "mom_pip", &mom_pip );
  pipipnnTree->Branch( "mom_pim", &mom_pim );
  pipipnnTree->Branch( "mom_p", &mom_p );
  pipipnnTree->Branch( "mom_n", &mom_n );
  pipipnnTree->Branch( "beta", &beta );
  pipipnnTree->Branch( "dE", &dE );
  pipipnnTree->Branch( "vtx_reaction", &vtx_reaction );
  pipipnnTree->Branch( "run_num", &run_num );
  pipipnnTree->Branch( "event_num", &event_num );
  pipipnnTree->Branch( "block_num", &block_num );
  pipipnnTree->Branch( "kf1mom_beam",   &kf1mom_beam );
  pipipnnTree->Branch( "kf1mom_pip", &kf1mom_pip );
  pipipnnTree->Branch( "kf1mom_pim", &kf1mom_pim );
  pipipnnTree->Branch( "kf1mom_p", &kf1mom_p );
  pipipnnTree->Branch( "kf1mom_n", &kf1mom_n );
  pipipnnTree->Branch( "kf1_chi2", &kf1_chi2 );
  pipipnnTree->Branch( "kf1_NDF", &kf1_NDF );
  pipipnnTree->Branch( "kf1_status", &kf1_status );
  pipipnnTree->Branch( "kf1_pvalue", &kf1_pvalue );
  pipipnnTree->Branch( "kf2mom_beam",   &kf2mom_beam );
  pipipnnTree->Branch( "kf2mom_pip", &kf2mom_pip );
  pipipnnTree->Branch( "kf2mom_pim", &kf2mom_pim );
  pipipnnTree->Branch( "kf2mom_p", &kf2mom_p );
  pipipnnTree->Branch( "kf2mom_n", &kf2mom_n );
  pipipnnTree->Branch( "kf2_chi2", &kf2_chi2 );
  pipipnnTree->Branch( "kf2_NDF", &kf2_NDF );
  pipipnnTree->Branch( "kf2_status", &kf2_status );
  pipipnnTree->Branch( "kf2_pvalue", &kf2_pvalue );
  pipipnnTree->Branch( "kf3mom_beam",   &kf3mom_beam );
  pipipnnTree->Branch( "kf3mom_pip", &kf3mom_pip );
  pipipnnTree->Branch( "kf3mom_pim", &kf3mom_pim );
  pipipnnTree->Branch( "kf3mom_p", &kf3mom_p );
  pipipnnTree->Branch( "kf3mom_n", &kf3mom_n );
  pipipnnTree->Branch( "kf3_chi2", &kf3_chi2 );
  pipipnnTree->Branch( "kf3_NDF", &kf3_NDF );
  pipipnnTree->Branch( "kf3_status", &kf3_status );
  pipipnnTree->Branch( "kf3_pvalue", &kf3_pvalue );
  pipipnnTree->Branch( "kf_flag", &kf_flag );

  //** making classes **//
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
  AllGoodTrack = 0;
  nTrack = 0;
  CDC_Event_Number = 0;

  nFill_pppi = 0;
  nFill_pipipnn  = 0;

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

  //confMan->GetSlewingMapManager()->PrintMapCDS();

  //** set cut parameters **//
  if( GRUN=="49c" ){
    //** from EventAnalysis_Lpn_sada.cpp **// //!! sada-D p.73 !!//
    PARA_tof_K_MIN = 28.176;
    PARA_tof_K_MAX = 29.400;
    PARA_tof_pi_MIN = 24.0;
    PARA_tof_pi_MAX = 28.176;
    //** from EventAnalysis_Lpn_sada.cpp **// //!! sada-D p.73 !!//
    PARA_time_range_BLC = 5;
    //** from EventAnalysis_Lpn_sada.cpp **// //!! sada-D p.83 !!//
    PARA_blc2bpc_dx_MIN = -0.795;
    PARA_blc2bpc_dx_MAX = 0.822;
    PARA_blc2bpc_dy_MIN = -0.865;
    PARA_blc2bpc_dy_MAX = 0.871;
    PARA_blc2bpc_dxdz_MIN = -0.0240;
    PARA_blc2bpc_dxdz_MAX = 0.0250;
    PARA_blc2bpc_dydz_MIN = -0.02481;
    PARA_blc2bpc_dydz_MAX = 0.02489;
    //** Sada PTEP value **//
    PARA_lnL_MAX = 6.0;
  }else if ( GRUN=="65" ){
    //** from Yamaga:MyAnalysisBL.cpp **//
    PARA_tof_K_MIN = 27.9637;
    PARA_tof_K_MAX = 29.4648;
    PARA_tof_pi_MIN = 25.0;
    PARA_tof_pi_MAX = 27.0;
    //!! Yamaga-D v03 !!//
    PARA_time_range_BLC = 10;
    //** from Yamaga:MyAnalysisBL.cpp **// //!! Yamaga-D v03 !!//
    PARA_blc2bpc_dx_MIN = -0.75;
    PARA_blc2bpc_dx_MAX = 0.75;
    PARA_blc2bpc_dy_MIN = -0.75;
    PARA_blc2bpc_dy_MAX = 0.75;
    PARA_blc2bpc_dxdz_MIN = -0.02;
    PARA_blc2bpc_dxdz_MAX = 0.02;
    PARA_blc2bpc_dydz_MIN = -0.02;
    PARA_blc2bpc_dydz_MAX = 0.02;
    //!! Yamaga-D v03 !!//
    PARA_lnL_MAX = 9.95;
  }else{
    std::cerr<<" !!! pleas set paramter::GRUN in EventAnalysis_Lpn_sakuma.cpp => 49c or 65 !!!"<<std::endl;
    exit( false );
  }

  cerr<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
  if( VTX_SIGMA==1 )
    cerr<<"  Sigma reconstruction mode,  vertex = K- & p "  <<endl;
  else if( VTX_SIGMA==0 )
    cerr<<"  Lambda reconstruction mode, vertex = K- & pi+ "<<endl;
  cerr<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;


  //-----------------------------------------//
  //--- covariance matrices for KinFitter ---//
  //-----------------------------------------//
  for(int i=0; i<7; i++){
    for(int j=0; j<16; j++){
      covVal2[i][j] = covVal1[i][j];
    }
  }
  covZero = new TMatrixD(4, 4);
  covZero->Zero();
  covZero->ResizeTo(3, 3); // resize from 4x4 to 3x3
  for( int i=0; i<7; i++ ){
    covParticle1[i] = new TMatrixD(4, 4);
    covParticle2[i] = new TMatrixD(4, 4);
    covParticle3[i] = new TMatrixD(4, 4);
    int n = 0;
    for( int j=0; j<4; j++ ){
      for( int k=0; k<4; k++ ){
	if( j==k ){
	  (*covParticle1[i])[j][k] = covVal1[i][n]; // only diagonal elements
	  (*covParticle2[i])[j][k] = covVal2[i][n]; // only diagonal elements
	  (*covParticle3[i])[j][k] = covVal3[i][n]; // only diagonal elements
	} else{
	  (*covParticle1[i])[j][k] = 0;
	  (*covParticle2[i])[j][k] = 0;
	  (*covParticle3[i])[j][k] = 0;
	}
	n++;
      }
    }
    covParticle1[i]->ResizeTo(3, 3); // resize from 4x4 to 3x3
    covParticle2[i]->ResizeTo(3, 3); // resize from 4x4 to 3x3
    covParticle3[i]->ResizeTo(3, 3); // resize from 4x4 to 3x3
    covParticle1[i]->Print(); // Print all
    covParticle2[i]->Print(); // Print all
    covParticle3[i]->Print(); // Print all
  }
  //-----------------------------------------//
  //--- covariance matrices for KinFitter ---//
  //-----------------------------------------//

  pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");

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

  TH1F *h1;
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
  h1 = (TH1F*)gFile->Get( "Time" );
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

  //** # of tracks cut **//
  if( nGoodTrack<3 ){
    Clear( nAbort_nGoodTrack );
    return true;
  }

  //** # of CDH-hits cut **//
  int nCDH = 0;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    //if( cdsMan->CDH(i)->CheckRange() ) nCDH++; //** only requirement of TDC **//
    if( cdsMan->CDH(i)->CheckRange() && cdsMan->CDH(i)->ctmean()<TDC_CDH_MAX )
      nCDH++;
  }
  Tools::Fill1D( Form("mul_CDH"), nCDH );
  if( nCDH!=4 ){ //** only 4 hits events **//
    Clear( nAbort_nCDH );
    return true;
  }

  //** + + + + + + + + + + + + **//
  //**  beamline analysis      **//
  //** + + + + + + + + + + + + **//

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
    return true;      
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
  int pid_beam = 3; // 0:pi 1:K 3:else
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      ctmBHD = blMan->BHD(i)->ctmean();
      double tofBHDT0 = ctmT0-ctmBHD;
      Tools::Fill1D( Form("tof_T0BHD"), tofBHDT0 );      
      if( header->kaon() && PARA_tof_K_MIN<tofBHDT0 && tofBHDT0<PARA_tof_K_MAX )
	pid_beam = Beam_Kaon; 
      else if( header->pion() && PARA_tof_pi_MIN<tofBHDT0 && tofBHDT0<PARA_tof_pi_MAX )
	pid_beam = Beam_Pion;
    }
  }
  if( pid_beam==3 ){ //** unidentified particle is discarded (other than pi/K) **//
    Clear( nAbort_pid_beam );
    return true;
  }

  //** BPC track selection **//
  int nbpc = 0;
  int bpcid = -1;
  double chibpc = 999;
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    Tools::Fill1D( Form("tracktime_BPC"), bltrackMan->trackBPC(i)->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BPC"), bltrackMan->trackBPC(i)->chi2all() );
    if( bltrackMan->trackBPC(i)->CheckRange(-30,100) ){
      nbpc++;
      bpcid = i;
      chibpc = bltrackMan->trackBPC(i)->chi2all();
    }
  }
  Tools::Fill1D( Form("ntrack_BPC"), nbpc );
  if( nbpc!=1 ){
    Clear( nAbort_nbpc );
    return true;
  }
  LocalTrack *bpctrack = bltrackMan->trackBPC(bpcid);
  if( !(bpctrack->CheckRange(-10,10)) || bpctrack->chi2all()>10 ){ //!! sada-D p.74 !!// //!! Yamaga-D v03 !!//
    Clear( nAbort_bpctrack );
    return true;
  }

  //** vertex calculation **//
  for( int it1=0; it1<trackMan->nGoodTrack(); it1++ ){
    trackMan->CalcVertex_beam( trackMan->GoodTrackID(it1), bltrackMan, confMan );
  }

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

  bool flagbmom = false;
  TLorentzVector Lpipdef, Lpimdef;
  TVector3 vtx_react;

  int nblc1 = 0;
  int nblc2 = 0;
  int blc1id = -1;
  int blc2id = -1;

  //** timing selection of BLC1/BLC2 **//
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    LocalTrack *blc1 = bltrackMan->trackBLC1(i);
    Tools::Fill1D( Form("tracktime_BLC1"), blc1->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BLC1"), blc1->chi2all() );
    if( blc1->CheckRange(-30,100) ){
      nblc1++;
      if( blc1->CheckRange(-PARA_time_range_BLC, PARA_time_range_BLC) &&
	  bltrackMan->trackBLC1(i)->chi2all()<10 ) blc1id = i;
    }	
  }
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    LocalTrack *blc2 = bltrackMan->trackBLC2(i);
    Tools::Fill1D( Form("tracktime_BLC2"), blc2->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BLC2"), blc2->chi2all() );
    if( blc2->CheckRange(-30,100) ){
      nblc2++;
      if( blc2->CheckRange(-PARA_time_range_BLC, PARA_time_range_BLC) &&
	  bltrackMan->trackBLC2(i)->chi2all()<10 ) blc2id = i;
    }	
  }
  Tools::Fill1D( Form("ntrack_BLC1"), nblc1 );
  Tools::Fill1D( Form("ntrack_BLC2"), nblc2 );

  //** sngle track selection in each BLC **//
  if( !(nblc1==1 && blc1id!=-1 && nblc2==1 && blc2id!=-1) ){ //** multi-good-beams event is ignored (20170614) **//
    Clear( nAbort_singleBLtrack );
    return true;
  }

  //** BLC2-BPC track matching **//
  bool fblc2bpc = false;
  for( int ii=0; ii<bltrackMan->ntrackBLC2(); ii++ ){
    if( ii!=blc2id ) continue;
    LocalTrack *blc2 = bltrackMan->trackBLC2(ii);
    double xblc2bpc[2], yblc2bpc[2];	 
    double xmom[2], ymom[2];

    TVector3 Pos_BPC, Pos_BLC2, tmp;
    confMan->GetBLDCWireMapManager()->GetGParam( CID_BPC, Pos_BPC, tmp );
    confMan->GetBLDCWireMapManager()->GetGParam( CID_BLC2a, Pos_BLC2, tmp );
    double zPos_BPC = Pos_BPC.Z();
    double zPos_BLC2 = Pos_BLC2.Z();
    double zPos_BPC_BLC2 = (Pos_BPC.Z()+Pos_BLC2.Z())/2;

    bpctrack->XYPosatZ( zPos_BPC_BLC2, xblc2bpc[0], yblc2bpc[0] );
    bpctrack->XYPosatZ( zPos_BPC, xmom[0], ymom[0] );
    blc2->XYPosatZ( zPos_BPC_BLC2, xblc2bpc[1], yblc2bpc[1]);
    blc2->XYPosatZ( zPos_BLC2, xmom[1], ymom[1]);
    double dxdz[2], dydz[2];
    dxdz[0] = (xmom[0]-xblc2bpc[0]) / (zPos_BPC-zPos_BPC_BLC2);
    dxdz[1] = (xmom[1]-xblc2bpc[1]) / (zPos_BLC2-zPos_BPC_BLC2);
    dydz[0] = (ymom[0]-yblc2bpc[0]) / (zPos_BPC-zPos_BPC_BLC2);
    dydz[1] = (ymom[1]-yblc2bpc[1]) / (zPos_BLC2-zPos_BPC_BLC2);

    if(      (xblc2bpc[1]-xblc2bpc[0])<PARA_blc2bpc_dx_MIN ||
	     (xblc2bpc[1]-xblc2bpc[0])>PARA_blc2bpc_dx_MAX ) fblc2bpc = false;
    else if( (yblc2bpc[1]-yblc2bpc[0])<PARA_blc2bpc_dy_MIN ||
	     (yblc2bpc[1]-yblc2bpc[0])>PARA_blc2bpc_dy_MAX ) fblc2bpc = false;
    else if( (dxdz[1]-dxdz[0])<PARA_blc2bpc_dxdz_MIN ||
	     (dxdz[1]-dxdz[0])>PARA_blc2bpc_dxdz_MAX ) fblc2bpc = false;
    else if( (dydz[1]-dydz[0])<PARA_blc2bpc_dydz_MIN ||
	     (dydz[1]-dydz[0])>PARA_blc2bpc_dydz_MAX ) fblc2bpc = false;
    else fblc2bpc = true;
    
    Tools::Fill2D( Form("dydx_BLC2BPC"), xblc2bpc[1]-xblc2bpc[0], yblc2bpc[1]-yblc2bpc[0] );
    Tools::Fill2D( Form("dydzdxdz_BLC2BPC"), dxdz[1]-dxdz[0], dydz[1]-dydz[0] );
  }
  if( !fblc2bpc ){
    Clear( nAbort_fblc2bpc );
    return true;
  }

  //** beam momentum calculation **//
  TLorentzVector L3_beambf;  // 4-Momentum(beam) in LAB
  TLorentzVector L3_beam;    // 4-Momentum(beam) in LAB with dE correcion
  TLorentzVector L3_target;  // 4-Momentum(He3-target) in LAB
  TLorentzVector L3_targetP; // 4-Momentum(p-target) in LAB
  TLorentzVector L3_beambfCM;  // 4-Momentum(beam) in CM
  TLorentzVector L3_beamCM;    // 4-Momentum(beam) in CM with dE correcion
  TLorentzVector L3_targetCM;  // 4-Momentum(He3-target) in CM
  TLorentzVector L3_targetPCM; // 4-Momentum(p-target) in CM

  double bchi = 999;
  for( int ii=0; ii<bltrackMan->ntrackBLC2(); ii++ ){
    for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
      if( nblc2!=1 || ii!=blc2id ) continue;
      if( nblc1!=1 || i!=blc1id ) continue;
      LocalTrack *blc1 = bltrackMan->trackBLC1(i);
      LocalTrack *blc2 = bltrackMan->trackBLC2(ii);
      TVector3 Pos_T0;
      confMan->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
      double zPos_T0 = Pos_T0.Z();
      
      TVector3 blc2t0 = blc2->GetPosatZ( zPos_T0 );
      TVector3 bpct0 = bpctrack->GetPosatZ( zPos_T0 );

      BeamSpectrometer *beamsp = new BeamSpectrometer( confMan );
      beamsp->TMinuitFit( blc1, blc2, confMan );
      double beammom = beamsp->mom();
      double bchitmp = beamsp->chisquare(); 

      if( pid_beam==Beam_Kaon && bchitmp<bchi ){
	bchi = bchitmp;
	double x1, y1, x2, y2;
	double z1 = 0, z2 = 20;
	bpctrack->XYPosatZ( z1, x1, y1 );		  
	bpctrack->XYPosatZ( z2, x2, y2 );
	TVector3 lp; lp.SetXYZ( x1, y1, z1 );
	TVector3 ls; ls.SetXYZ( x2-x1, y2-y1, z2-z1); ls = ls.Unit();	       
	TVector3 Pp_beam = beammom*ls; 
	TVector3 Pp_target; Pp_target.SetXYZ( 0, 0, 0 ); 

	L3_beambf.SetVectM( Pp_beam , kpMass );
	L3_target.SetVectM( Pp_target, ThreeHeMass );
	L3_targetP.SetVectM( Pp_target, pMass );
	L3_beam = L3_beambf;
	TVector3 boost = (L3_target+L3_beam).BoostVector();
	L3_beambfCM = L3_beam;
	L3_targetCM = L3_target;
	L3_targetPCM = L3_targetP;
	L3_beambfCM.Boost( -1*boost );
	L3_targetCM.Boost( -1*boost );
	L3_targetPCM.Boost( -1*boost );
	
	if( bchi<20 ) flagbmom = true; //!! sada-D p.80 !!// //!! Yamaga-D v03 !!//
      }
      delete beamsp;
      
    } // for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
  } // for( int ii=0; ii<bltrackMan->ntrackBLC2(); ii++ ){

  Tools::Fill1D( Form("trackchi2_beam"), bchi );
  if( !flagbmom ){
    Clear( nAbort_flagbmom );
    return true;
  }
  Tools::Fill1D( Form("momentum_beam"), L3_beambf.P() );


#if 0 //** not need to cut [20170626] **//
  //** Fiducial volume cut for the beam @ z=0 **//
  bool ftarget = false;
  double z = 0;
  double x, y;
  bpctrack->XYPosatZ( z, x, y );
  TVector3 vtx_tar;
  vtx_tar.SetXYZ( x, y, z );
  if( GeomTools::GetID(vtx_tar)==CID_Fiducial ){
    ftarget = true;
  } 
  if( !ftarget ){
    Clear( nAbort_ftarget );
    return true;
  }
#endif



  //** + + + + + + + + + + + + **//
  //**  PID in CDS             **//
  //** + + + + + + + + + + + + **//
  
  int CDHseg;

  //** PID of CDS tracks **//
  for( int it=0; it<trackMan->nGoodTrack(); it++ ){
    CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );

    //** chi2 cut can be applied in CDSTrackingMan with MaxChi in CDSFittingParam_posi.param **//
    Tools::Fill1D( Form("trackchi2_CDC"), track->Chi() );
    if( track->Chi()>30 ) continue; //!! sada-D p.90 !!//

    if( !track->CDHFlag() ) continue;

    double mom = track->Momentum();
    TVector3 vtxb1, vtxb2, vtxb;
    track->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtxb1, vtxb2 );
    track->SetPID(-1);
    vtxb = (vtxb1+vtxb2)*0.5;

    double tof = 999.;
    double mass2 = -999.;
    for( int icdh=0; icdh<track->nCDHHit(); icdh++ ){
      HodoscopeLikeHit *cdhhit = track->CDHHit( cdsMan, icdh );
      double tmptof = cdhhit->ctmean()-ctmT0;
      if( tmptof<tof || tof==999. ){
	tof = tmptof;
	CDHseg = cdhhit->seg();
      }
    }

    bool CDHflag = true;
    for( int m=0; m<(int)vCDHseg.size(); m++ ){
      if( CDHseg==vCDHseg[m] ) CDHflag = false;
    }
    if( !CDHflag ) continue;
    vCDHseg.push_back( CDHseg );

    //** calculation of beta and squared-mass **//
    double tmptof, beta_calc;
    if( !TrackTools::FindMass2( track, bpctrack, tof, L3_beam.Vect().Mag(),
				pid_beam, beta_calc, mass2, tmptof ) ){
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      continue;
    }

    //** Retiming of CDC track by CDH info. **//
    track->Retiming( cdsMan, confMan, beta_calc, true );
    for( int m=0; m<5; m++ ){
      track->HelixFitting( cdsMan ); //** 5 times iteration **//
    }
    track->Calc( confMan );

    //** finalize PID **//
    if( !TrackTools::FindMass2( track, bpctrack, tof, L3_beam.Vect().Mag(),
				pid_beam, beta_calc, mass2, tmptof ) ){ //** not FindMass2C() [20170622] **//
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      continue;
    }
    int pid = -1;
    if( GRUN=="49c" ){
      pid = TrackTools::PIDcorr( mom, mass2 ); //!! sada-D p.92 !!//
    }else if ( GRUN=="65" ){
      pid = TrackTools::PIDcorr3( mom, mass2 ); //** from Yamaga library **// //!! Yamaga-D v03 !!//
    }
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

    if( pid==CDS_PiMinus )
      pim_ID.push_back( trackMan->GoodTrackID(it) );
    else if( pid==CDS_PiPlus )
      pip_ID.push_back( trackMan->GoodTrackID(it) );
    else if( pid==CDS_Proton )
      p_ID.push_back( trackMan->GoodTrackID(it) );
    else if( pid==CDS_Deuteron )
      d_ID.push_back( trackMan->GoodTrackID(it) );
    else if( pid==CDS_Kaon )
      km_ID.push_back( trackMan->GoodTrackID(it) );

  } // for( int it=0; it<trackMan->nGoodTrack(); it++ ){
  //** end of PID **//

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
  if( GRUN=="49c" ){
    if( nBVC ) chargedhit = true;
  }else if( GRUN=="65" ){
    if( nBVC || nCVC || nPC ) chargedhit = true;
  }

#if 0
  std::cout<<Event_Number<<":"<<Block_Event_Number<<":"<<CDC_Event_Number<<" | "
	   <<p_ID.size()<<","<<pim_ID.size()<<","<<trackMan->nGoodTrack()<<" |"
	   <<flagbmom<<","<<chargedhit<<std::endl;
#endif


  //** + + + + + + + + + + + **//
  //**  p pi+ pi- X event  **//
  //** + + + + + + + + + + + **//

  if( flagbmom && p_ID.size()==1 && pim_ID.size()==1 &&pip_ID.size()==1 &&
      trackMan->nGoodTrack()==3 && !chargedhit ){

    //=== condense p pi+ pi- X candidates ===//
    rtFile2->cd();
    evTree->Fill();
    rtFile->cd();
    std::cout<<"### filled: Event_Number, Block_Event_Number, CDC_Event_Number = "
	     <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
    nFill_pppi++;
    //=== condense p pi+ pi- X event ===//
    
    //** find CDH hit from neutral particles **//
    std::vector <int> nCDHseg;
    std::vector <int> CDHhit_list;
    for( int n=0; n<cdsMan->nCDH(); n++ ){
      //if( cdsMan->CDH(n)->CheckRange() )
      if( cdsMan->CDH(n)->CheckRange() && cdsMan->CDH(n)->ctmean()<TDC_CDH_MAX )
	CDHhit_list.push_back( cdsMan->CDH(n)->seg() );
    }
    std::sort(vCDHseg.begin(), vCDHseg.end());
    std::sort(CDHhit_list.begin(), CDHhit_list.end());
    std::set_difference( CDHhit_list.begin(), CDHhit_list.end(),
			 vCDHseg.begin(), vCDHseg.end(),
			 std::back_inserter(nCDHseg) );

    if( nCDHseg.size()!=1 ){
      std::cerr<<" CDH neutral hit is not 1 :: "<<nCDHseg.size()<<std::endl;
    }

    std::cerr<<"# of diff = "<<nCDHseg.size()<<std::endl;
    std::cerr<<"CDH hits =   ";
    for( int n=0; n<(int)CDHhit_list.size(); n++ ){
      std::cerr<<CDHhit_list[n]<<" ";
    } std::cerr<<std::endl;
    std::cerr<<"track hits = ";
    for( int n=0; n<(int)vCDHseg.size(); n++ ){
      std::cerr<<vCDHseg[n]<<" ";
    } std::cerr<<std::endl;
    std::cerr<<"diff hits =  ";
    for( int n=0; n<(int)nCDHseg.size(); n++ ){
      std::cerr<<nCDHseg[n]<<" ";
    } std::cerr<<std::endl;
    
    //** isolation cut **//
    int flag_isolation = 0;
    for( int l=0; l<(int)nCDHseg.size(); l++ ){
      for( int m=0; m<(int)CDHhit_list.size(); m++ ){
	if( nCDHseg[l]-CDHhit_list[m] ) Tools::Fill1D( Form("diff_CDH"), nCDHseg[l]-CDHhit_list[m] );
	if( abs(nCDHseg[l]-CDHhit_list[m])==1 || abs(nCDHseg[l]-CDHhit_list[m])==35 )
	  flag_isolation++;
      }
    }
    if( flag_isolation ){
      std::cerr<<"CDH hit candidate is NOT isolated !!!"<<std::endl;
      Clear( nAbort_CDHiso );
      return true;
    }

    //** copy neutral CDH hit candidate **//
    int icdh = -1;
    for( int n=0; n<cdsMan->nCDH(); n++ ){
      if( cdsMan->CDH(n)->seg()==nCDHseg[0] ) icdh = n;
    }
    HodoscopeLikeHit *ncdhhit = cdsMan->CDH(icdh);

    //** charge veto using CDC **//
    TVector3 Pos_CDH;
    confMan->GetGeomMapManager()->GetPos( CID_CDH, ncdhhit->seg(), Pos_CDH );
    std::cerr<<"CDH candidate = "<<ncdhhit->seg()<<" -> "<<Pos_CDH.Phi()/TwoPi*360<<" deg"<<std::endl;

    const double PhiMin = -15.0/360*TwoPi; // rad
    const double PhiMax =  15.0/360*TwoPi; // rad
    std::cerr<<"Min/Max = "<<PhiMin/TwoPi*360<<"/"<<PhiMax/TwoPi*360<<" deg"<<std::endl;

    int nCDC = 0;
    for( int l=14; l<16; l++ ){ // charge veto using layer 14, 15
      for( int m=0; m<cdsMan->nCDC(l); m++ ){
	CDCHit *cdc=cdsMan->CDC(l,m);
	TVector3 Pos_CDC = cdc->wpos();
	Pos_CDC.SetZ(0); // only xy pos is used
	double angle = Pos_CDC.Angle(Pos_CDH); // rad
	std::cerr<<"CDC "<<l<<" "<<m<<" "<<cdc->wire()<<" -> "<<Pos_CDC.Phi()/TwoPi*360
		 <<" deg :: diff = "<<angle/TwoPi*360<<" deg"<<std::endl;
	Tools::Fill1D( Form("diff_CDH_CDC"), angle/TwoPi*360 );
	if( PhiMin<angle && angle<PhiMax ) nCDC++;
      }
    }
    std::cerr<<"# of CDC hits for nCDH candidate = "<<nCDC<<std::endl;

    Pos_CDH.SetZ(-1*ncdhhit->hitpos()); // (-1*) is correct in data analysis [20170926]
    //Pos_CDH.SetZ(ncdhhit->hitpos());


    //** neutral particle in CDH **//
    if( !nCDC ){
      CDSTrack *track_p   = trackMan->Track( p_ID[0] );   // only 1 track
      CDSTrack *track_pip = trackMan->Track( pip_ID[0] ); // only 1 track
      CDSTrack *track_pim = trackMan->Track( pim_ID[0] ); // only 1 track

      TVector3 vtx_b; // Vertex(baem-particle)_on_beam
      TVector3 vtx_p; // Vertex(baem-particle)_on_particle
      if( VTX_SIGMA==1 ){ // 1: Sigma reconstruction,  vertex = K- & p
	track_p->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      }
      else if( VTX_SIGMA==0 ){ // 0: Lambda reconstruction, vertex = K- & pi+
	track_pip->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      }
      vtx_react = 0.5*(vtx_b+vtx_p); // reaction vertex

      // from here
      double tof = 999.;
      for( int icdh=0; icdh<track_p->nCDHHit(); icdh++ ){
	HodoscopeLikeHit *cdhhit = track_p->CDHHit( cdsMan, icdh );
	double tmptof = cdhhit->ctmean()-ctmT0;
	if( tmptof<tof || tof==999. ){
	  tof = tmptof;
	  CDHseg = cdhhit->seg();
	}
      }
      // to here, meaningless??? [20180524]

      //** beam kaon tof **//
      TVector3 Pos_T0;
      confMan->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
      double beamtof, momout;
      double z_pos = Pos_T0.Z();;
      ELossTools::CalcElossBeamTGeo( bpctrack->GetPosatZ(z_pos), vtx_react,
				     L3_beambf.Vect().Mag(), kpMass, momout, beamtof );
      L3_beam.SetVectM( momout*L3_beambf.Vect().Unit(), kpMass );
      double ntof = ncdhhit->ctmean()-ctmT0-beamtof;
      double nlen = (Pos_CDH-vtx_react).Mag();
      beta = nlen/ntof/(Const*100);
      double tmp_mom = beta<1 ? nMass*beta/sqrt(1-beta*beta) : 0;
      std::cerr<<"$$$ beta = "<<beta<<" mom_n = "<<tmp_mom<<std::endl; //" "<<1/sqrt(1+nMass*nMass)<<std::endl;

      //** reconstructoin of missing neutorn **//
      TVector3 P_p;   // Momentum(p)
      TVector3 P_pim; // Momentum(pi-)
      TVector3 P_pip; // Momentum(pi+)
      TVector3 P_n;   // Momentum(n)
      
      TLorentzVector L_p;   // 4-Momentum(p)
      TLorentzVector L_pim; // 4-Momentum(pi-)
      TLorentzVector L_pip; // 4-Momentum(pi+)
      TLorentzVector L_n;   // 4-Momentum(n)
      TLorentzVector L_nmiss; // 4-Momentum(n_miss)

      track_p->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      double dca_p  = (vtx_b-vtx_p).Mag(); // DCA(beam-p)
      if( !track_p->GetMomentum( vtx_p, P_p, true, true ) ){
	std::cerr<<" !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }
      track_pip->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      double dca_pip  = (vtx_b-vtx_p).Mag(); // DCA(beam-pip)
      if( !track_pip->GetMomentum( vtx_p, P_pip, true, true ) ){
	std::cerr<<" !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }
      track_pim->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtx_b, vtx_p );
      double dca_pim  = (vtx_b-vtx_p).Mag(); // DCA(beam-pim)
      if( !track_pim->GetMomentum( vtx_p, P_pim, true, true ) ){
	std::cerr<<" !!! failure in momentum calculation [GetMomentum()] !!! "<<std::endl;
      }
      P_n = tmp_mom*(Pos_CDH-vtx_react).Unit();
      
      L_p.SetVectM(   P_p,   pMass );
      L_pim.SetVectM( P_pim, piMass );
      L_pip.SetVectM( P_pip, piMass );
      L_n.SetVectM(   P_n,   nMass );
      
      double mm_mass   = (L3_target+L3_beam-L_p-L_pim-L_pip-L_n).M();
      TVector3 P_missn = (L3_target+L3_beam-L_p-L_pim-L_pip-L_n).Vect();
      L_nmiss.SetVectM( P_missn, nMass );
      std::cerr<<"  missing mass = "<<mm_mass<<std::endl;

      TVector3 boost = (L3_target+L3_beam).BoostVector();
      TLorentzVector L_nmiss_CM = L_nmiss;
      TLorentzVector L3_beam_CM = L3_beam;
      L_nmiss_CM.Boost(-boost);
      L3_beam_CM.Boost(-boost);
      double cos_n = L_nmiss_CM.Vect().Dot(L3_beam_CM.Vect())/(L_nmiss_CM.Vect().Mag()*L3_beam_CM.Vect().Mag());
      std::cerr<<"  missing mom | cos_CM = "<<cos_n<<std::endl;


      //** + + + + + + + + + + + + + **//
      //**  fill histograms & tree   **//
      //** + + + + + + + + + + + + + **//

      const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
      //const double dE_MIN = 5.0; // 8.0MeVee * 3cm / 5cm;
      const double dE_MIN = 0.0;

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

      Tools::Fill2D( Form("dE_betainv"), 1/beta, ncdhhit->emean() );
      Tools::Fill2D( Form("MMom_MMass"), mm_mass, P_missn.Mag() );
      
      if( GeomTools::GetID(vtx_react)==CID_Fiducial ){
	Tools::Fill2D( Form("dE_betainv_fiducial"), 1/beta, ncdhhit->emean() );
	Tools::Fill2D( Form("MMom_MMass_fiducial"), mm_mass, P_missn.Mag() );

	if(  beta<beta_MAX ){
	  Tools::Fill2D( Form("dE_betainv_fiducial_beta"), 1/beta, ncdhhit->emean() );
	  Tools::Fill2D( Form("MMom_MMass_fiducial_beta"), mm_mass, P_missn.Mag() );

	  if( dE_MIN<ncdhhit->emean() ){
	    Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE"), 1/beta, ncdhhit->emean() );
	    Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE"), mm_mass, P_missn.Mag() );

	    Tools::Fill1D( Form("IMpipi"), (L_pim+L_pip).M() );
	    Tools::Fill1D( Form("IMppi"), (L_p+L_pim).M() );

	    // ********************** //	      
	    // *** pi Sigma mode *** //
	    // ********************** //
#if 0
	    if( VTX_SIGMA==1 &&
		((L_pim+L_pip).M()<pipi_MIN || pipi_MAX<(L_pim+L_pip).M()) &&
		((L_p+L_pim).M()<ppi_MIN || ppi_MAX<(L_p+L_pim).M()) ){ // K0 & Lambda subtraction
#else
	    if( VTX_SIGMA==1 ){
#endif
	      Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE_res"), 1/beta, ncdhhit->emean() );
	      Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE_res"), mm_mass, P_missn.Mag() );

	      if( neutron_MIN<mm_mass && mm_mass<neutron_MAX ){ // missing n selection
		Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE_res_n"), 1/beta, ncdhhit->emean() );
		Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE_res_n"), mm_mass, P_missn.Mag() );

		Tools::Fill2D( Form("MMom_NMom"), P_n.Mag(), P_missn.Mag() );
		Tools::Fill2D( Form("IMnpim_IMnpip"), (L_n+L_pip).M(), (L_n+L_pim).M() );
		
		if( (Sigmap_MIN<(L_n+L_pip).M() && (L_n+L_pip).M()<Sigmap_MAX) ||
		    (Sigmam_MIN<(L_n+L_pim).M() && (L_n+L_pim).M()<Sigmam_MAX) ){ // Sigma selection
		  Tools::Fill2D( Form("IMmnpim_IMmnpip"), (L_nmiss+L_pip).M(), (L_nmiss+L_pim).M() );
		  Tools::Fill2D( Form("MMnppip_MMnppim"), (L3_target+L3_beam-L_p-L_pim-L_n).M(),
				 (L3_target+L3_beam-L_p-L_pip-L_n).M() );
		  
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
	      //  = 1) TLorentzVector L3_beam, L_pim, (L_n+L_pip), L_p, L_nmiss, L_n, L_pip = for pi- Sigma+
	      TLorentzVector TL_meas1[7]; // measured
	      TL_meas1[0] = L3_beam;
	      TL_meas1[1] = L_pim;
	      TL_meas1[2] = (L_n+L_pip);
	      TL_meas1[3] = L_p;
	      TL_meas1[4] = L_nmiss;
	      TL_meas1[5] = L_n;
	      TL_meas1[6] = L_pip;
	      //  = 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
	      TLorentzVector TL_meas2[7]; // measured
	      TL_meas2[0] = L3_beam;
	      TL_meas2[1] = L_pip;
	      TL_meas2[2] = (L_n+L_pim);
	      TL_meas2[3] = L_p;
	      TL_meas2[4] = L_nmiss;
	      TL_meas2[5] = L_n;
	      TL_meas2[6] = L_pim;
	      TLorentzVector TL_kfit1[7]; // kinematical fitted
	      TLorentzVector TL_kfit2[7]; // kinematical fitted
	      // L3_target is defined as (0, 0, 0, M_3He)
	      TVector3 TV_target = L3_target.Vect();
	      TVector3 TV_meas1[7];
	      TVector3 TV_meas2[7];
	      for( int i=0; i<7; i++ ){
		TV_meas1[i] = TL_meas1[i].Vect();
		TV_meas2[i] = TL_meas2[i].Vect();
	      }
	      
	      //--- KinFitter :: initialization ---//
	      //  = 1) TLorentzVector L3_beam, L_pim, (L_n+L_pip), L_p, L_nmiss, L_n, L_pip = for pi- Sigma+
	      //  = 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
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
	      //*** definition of the fitter ***//
	      TKinFitter kinfitter1;
	      TKinFitter kinfitter2;
	      // add measured particles
	      kinfitter1.addMeasParticles(&Particle1[0], &Particle1[1], &Particle1[3], &Particle1[5], &Particle1[6]); // K, pi-, p, n, pi+
	      kinfitter2.addMeasParticles(&Particle2[0], &Particle2[1], &Particle2[3], &Particle2[5], &Particle2[6]); // K, pi+, p, n, pi-
	      kinfitter1.addUnmeasParticles(&Particle1[4]); // missing-n
	      kinfitter2.addUnmeasParticles(&Particle2[4]); // missing-n
	      // add constraints
	      kinfitter1.addConstraint(&ConstMS1); // mass of Sigma+
	      kinfitter2.addConstraint(&ConstMS2); // mass of Sigma-
	      for( int i=0; i<4; i++ ){
		kinfitter1.addConstraint(&ConstEp1[i]); // 4-momentum conservation
		kinfitter2.addConstraint(&ConstEp2[i]); // 4-momentum conservation
	      }
	      //*** perform the fit ***//
	      kinfitter1.setMaxNbIter(50);       // max number of iterations
	      kinfitter2.setMaxNbIter(50);       // max number of iterations
	      kinfitter1.setMaxDeltaS(5e-5);     // max delta chi2
	      kinfitter2.setMaxDeltaS(5e-5);     // max delta chi2
	      kinfitter1.setMaxF(1e-4);          // max sum of constraints
	      kinfitter2.setMaxF(1e-4);          // max sum of constraints
	      kinfitter1.setVerbosity(KFDEBUG);  // verbosity level
	      kinfitter2.setVerbosity(KFDEBUG);  // verbosity level
	      kinfitter1.fit();
	      kinfitter2.fit();
	      //*** copy fit results ***//
	      for( int i=0; i<7; i++ ){
		TL_kfit1[i] = (*Particle1[i].getCurr4Vec());
		TL_kfit2[i] = (*Particle2[i].getCurr4Vec());
	      }
	      TL_kfit1[2] = TL_kfit1[5]+TL_kfit1[6];
	      TL_kfit2[2] = TL_kfit2[5]+TL_kfit2[6];
	      
	      
	      Tools::Fill2D( Form("KFchi2_vs"), kinfitter1.getS()/kinfitter1.getNDF(),
			     kinfitter2.getS()/kinfitter2.getNDF() );
	      
	      std::cerr<<"pi- S+ : status = "<<kinfitter1.getStatus()<<", chi2/NDF = "<<kinfitter1.getS()<<"/"<<kinfitter1.getNDF()<<std::endl;
	      std::cerr<<"pi+ S- : status = "<<kinfitter2.getStatus()<<", chi2/NDF = "<<kinfitter2.getS()<<"/"<<kinfitter2.getNDF()<<std::endl;
	      
	      //** fill tree **//
	      kf1mom_beam   = TL_kfit1[0];
	      kf1mom_pip    = TL_kfit1[6];
	      kf1mom_pim    = TL_kfit1[1];
	      kf1mom_p      = TL_kfit1[3];
	      kf1mom_n      = TL_kfit1[5];
	      kf1_chi2      = kinfitter1.getS();
	      kf1_NDF       = kinfitter1.getNDF();
	      kf1_status    = kinfitter1.getStatus();
	      kf1_pvalue    = ROOT::Math::chisquared_cdf_c(kinfitter1.getS(), kinfitter1.getNDF());
	      kf2mom_beam   = TL_kfit2[0];
	      kf2mom_pip    = TL_kfit2[1];
	      kf2mom_pim    = TL_kfit2[6];
	      kf2mom_p      = TL_kfit2[3];
	      kf2mom_n      = TL_kfit2[5];
	      kf2_chi2      = kinfitter2.getS();
	      kf2_NDF       = kinfitter2.getNDF();
	      kf2_status    = kinfitter2.getStatus();
	      kf2_pvalue    = ROOT::Math::chisquared_cdf_c(kinfitter2.getS(), kinfitter2.getNDF());
	      kf_flag       = 1;

	      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	      // %%% Kinematical Fit using KinFitter %%% //
	      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //

	    } // if( ((L_pim+L_pip).M()<pipi_MIN || pipi_MAX<(L_pim+L_pip).M()) &&


	    // ********************** //
	    // *** pi Lambda mode *** //
	    // ********************** //
#if 0
	    else if( VTX_SIGMA==0 &&
		     ((L_pim+L_pip).M()<pipi_MIN || pipi_MAX<(L_pim+L_pip).M()) ){ // K0 subtraction
#else
	    else if( VTX_SIGMA==0 ){
#endif
	      Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE_Lambda"), 1/beta, ncdhhit->emean() );
	      Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE_Lambda"), mm_mass, P_missn.Mag() );
	      
	      if ( ppi_MIN<(L_p+L_pim).M() && (L_p+L_pim).M()<ppi_MAX ){ // Lambda selection
		if( neutron_MIN<mm_mass && mm_mass<neutron_MAX ){ // missing n selection
		  Tools::Fill2D( Form("dE_betainv_fiducial_beta_dE_Lambda_n"), 1/beta, ncdhhit->emean() );
		  Tools::Fill2D( Form("MMom_MMass_fiducial_beta_dE_Lambda_n"), mm_mass, P_missn.Mag() );
		  
		  Tools::Fill2D( Form("MMom_NMom_Lambda"), P_n.Mag(), P_missn.Mag() );
		  Tools::Fill2D( Form("Cosn_IMnppipi_Lambda"), (L_n+L_p+L_pim+L_pip).M(), cos_n );
		  Tools::Fill2D( Form("Cosn_IMppipi_Lambda"), (L_p+L_pim+L_pip).M(), cos_n );
		  Tools::Fill2D( Form("IMppipi_IMnppipi_Lambda"), (L_n+L_p+L_pim+L_pip).M(), (L_p+L_pim+L_pip).M() );
		  
		  Tools::Fill1D( Form("DCA_p_Lambda"), dca_p );
		  Tools::Fill1D( Form("DCA_pip_Lambda"), dca_pip );
		  Tools::Fill1D( Form("DCA_pim_Lambda"), dca_pim );
		}
	      } // if ( ppi_MIN<(L_p+L_pim).M() && (L_p+L_pim).M()<ppi_MAX ){

	      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	      // %%% Kinematical Fit using KinFitter %%% //
	      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	      //--- set TLorentzVector ---//
	      // beam_K(K+), pi+, Lambda, n, n, p from L, pi- from L 
	      //  = 3) TLorentzVector L3_beam, L_pip, (L_p+L_pim), L_n, L_nmiss, L_p, L_pim = for pi+ Lambda
	      TLorentzVector TL_meas3[7]; // measured
	      TL_meas3[0] = L3_beam;
	      TL_meas3[1] = L_pip;
	      TL_meas3[2] = (L_p+L_pim);
	      TL_meas3[3] = L_n;
	      TL_meas3[4] = L_nmiss;
	      TL_meas3[5] = L_p;
	      TL_meas3[6] = L_pim;
	      TLorentzVector TL_kfit3[7]; // kinematical fitted
	      // L3_target is defined as (0, 0, 0, M_3He)
	      TVector3 TV_target = L3_target.Vect();
	      TVector3 TV_meas3[7];
	      for( int i=0; i<7; i++ ){
		TV_meas3[i] = TL_meas3[i].Vect();
	      }
	      
	      //--- KinFitter :: initialization ---//
	      //  = 3) TLorentzVector L3_beam, L_pip, (L_p+L_pim), L_n, L_nmiss, L_p, L_pim = for pi+ Lambda
	      //*** definition of fit particles in cartesian coordinates ***//
	      TString str_particle3[7] = {"L_beam", "L_pip", "L_Lam", "L_n", "L_mn", "L_p", "L_pim"};
	      TFitParticlePxPyPz ParticleTgt = TFitParticlePxPyPz("target", "target", &TV_target,
								  pdg->GetParticle("He3")->Mass(), covZero);
	      TFitParticlePxPyPz Particle3[7];
	      for( int i=0; i<7; i++ ){
		Particle3[i] = TFitParticlePxPyPz(str_particle3[i], str_particle3[i], &TV_meas3[i],
						  pdg->GetParticle(PDG3[i])->Mass(), covParticle3[i]);
	      }
	      //*** definition of constraints ***//
	      // constraint :: mass of Sigma
	      TFitConstraintM ConstMS3 = TFitConstraintM("M_Lam", "M_Lam", 0, 0, pdg->GetParticle(PDG3[2])->Mass());
	      ConstMS3.addParticles1(&Particle3[5], &Particle3[6]);
	      // constraint :: 4-momentum conservation
	      TFitConstraintEp ConstEp3[4];
	      TString str_constEp3[4]  = {"Px", "Py", "Pz", "E"};
	      for( int i=0; i<4; i++ ){
		ConstEp3[i] = TFitConstraintEp(str_constEp3[i], str_constEp3[i], 0, TFitConstraintEp::component(i), 0);
		ConstEp3[i].addParticles1(&ParticleTgt, &Particle3[0]);
		ConstEp3[i].addParticles2(&Particle3[1], &Particle3[3], &Particle3[4], &Particle3[5], &Particle3[6]);
	      }
	      
	      //--- KinFitter :: execution ---//
	      //*** definition of the fitter ***//
	      TKinFitter kinfitter3;
	      // add measured particles
	      kinfitter3.addMeasParticles(&Particle3[0], &Particle3[1], &Particle3[3], &Particle3[5], &Particle3[6]); // K, pi+, n, p, pi-
	      kinfitter3.addUnmeasParticles(&Particle3[4]); // missing-n
	      // add constraints
	      kinfitter3.addConstraint(&ConstMS3); // mass of Lambda
	      for( int i=0; i<4; i++ ){
		kinfitter3.addConstraint(&ConstEp3[i]); // 4-momentum conservation
	      }
	      //*** perform the fit ***//
	      kinfitter3.setMaxNbIter(50);       // max number of iterations
	      kinfitter3.setMaxDeltaS(5e-5);     // max delta chi2
	      kinfitter3.setMaxF(1e-4);          // max sum of constraints
	      kinfitter3.setVerbosity(KFDEBUG);  // verbosity level
	      kinfitter3.fit();
	      //*** copy fit results ***//
	      for( int i=0; i<7; i++ ){
		TL_kfit3[i] = (*Particle3[i].getCurr4Vec());
	      }
	      TL_kfit3[2] = TL_kfit3[5]+TL_kfit3[6];
	      
	      
	      Tools::Fill2D( Form("KFchi2_vs"), kinfitter3.getS()/kinfitter3.getNDF(), 1 );
	      
	      std::cerr<<"pi+ L : status = "<<kinfitter3.getStatus()<<", chi2/NDF = "<<kinfitter3.getS()<<"/"<<kinfitter3.getNDF()<<std::endl;
	      
	      //** fill tree **//
	      kf3mom_beam   = TL_kfit3[0];
	      kf3mom_pip    = TL_kfit3[1];
	      kf3mom_pim    = TL_kfit3[6];
	      kf3mom_p      = TL_kfit3[5];
	      kf3mom_n      = TL_kfit3[3];
	      kf3_chi2      = kinfitter3.getS();
	      kf3_NDF       = kinfitter3.getNDF();
	      kf3_status    = kinfitter3.getStatus();
	      kf3_pvalue    = ROOT::Math::chisquared_cdf_c(kinfitter3.getS(), kinfitter3.getNDF());
	      kf_flag       = 1;

	      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	      // %%% Kinematical Fit using KinFitter %%% //
	      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //

	    } // else if( (L_pim+L_pip).M()<pipi_MIN || pipi_MAX<(L_pim+L_pip).M() ){ 

	  } // if( dE_MIN<ncdhhit->emean() ){
	} // if(  beta<beta_MAX ){

	//** fill tree **//
	mom_beam   = L3_beam;   // 4-momentum(beam)
	mom_target = L3_target; // 4-momentum(target)
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

	std::cout<<"%%% pipipnn event: Event_Number, Block_Event_Number, CDC_Event_Number = "
		 <<Event_Number<<" , "<<Block_Event_Number<<" , "<<CDC_Event_Number<<std::endl;
	rtFile3->cd();
	pipipnnTree->Fill();
	rtFile->cd();
	nFill_pipipnn++;
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
  std::cout<<"*** # of pi+ pi- p n n events = "<<nFill_pipipnn<<" ***"<<std::endl;

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
  //** gneneral informantion **//
  Tools::newTH1F( Form("Time"), 3000, -0.5, 2999.5 );
  Tools::newTH1F( Form("EventCheck"), 20, 0, 20 );
  Tools::newTH1F( Form("Scaler"), 41, -0.5, 40.5 );

  //** CDC and CDH information from CDC-trackig file **//
  Tools::newTH1F( Form("nGoodTrack"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_CDH"), 11, -0.5, 10.5 );

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

  Tools::newTH2F( Form("KFchi2_vs"), 100, 0, 100, 100, 0, 100 );

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
