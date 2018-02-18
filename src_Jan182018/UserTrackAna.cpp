

#include <iostream>
#include <string>
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
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
#include <TLorentzVector.h>

#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
#include "Display3D.h"

  //#######Need knucl lib file!! Set your dir!! ###
#include "/home/ysada/knucl3/include/KnuclRootData.h"
  //#######################################


/* ####2011/07/21########
const double MassThPPi = 0.8; // GeV/c2
const double TOFOffset = 3; // nsec
*/
/* #####2011/07/25 slew2######
const double MassThPPi = 0.7; // GeV/c2
const double TOFOffset = -2.5; // nsec
*/

const double MassThPPi = 0.7; // GeV/c2
const double TOFOffset = -3.5; // nsec


const double piMass = 0.13957;
const double pMass = 0.938272;
const double nMass = 0.939565;
const double lMass = 1.115683;
const double kpMass = 0.4936;
const double Const=29.97;// [cm/ns]
#define SIM 0

int PID( double tof, double mom ) // 0:unknown, 1:pi+, 2;pi-, 3:proton
{
  //TF1 *fun1 = new TF1("fun1",Form("%lf/sqrt( pow((x+%lf),2)-1)",MassThPPi,TOFOffset),-2,14);
  int val=0;
  if( tof<TOFOffset ){
  }
  if( 0<mom ){
    if( mom < MassThPPi/sqrt(pow(tof+TOFOffset,2)-1) ){
      val = 1;
    }
    else{
      val = 3;
    }
  }
  else{
    if( fabs(mom) < MassThPPi/sqrt(pow(tof+TOFOffset,2)-1) ){
      val = 2;
    }
    else{
      val = 0;
    }
  }

  //std::cout << " tof:" << tof << " mom:" << mom << " type:" << val << std::endl;
  return val;
}

int main( int argc, char **argv )
{
 
  std::cout << " argc:" << argc << std::endl;
  if(argc != 4 )
    std::cout << "Plese set Conffile Outputfile Inputfile  "<< std::endl;

  for(int i=0; i<argc; i++ ){
    std::cout << "  " << argv[i] << std::endl;
  }

  std::string confFile,outFile,inFile;
  if( argv[1][0] != '!' ) confFile = argv[1];
  if( argv[2][0] != '!' ) outFile = argv[2];
  if( argv[3][0] != '!' ) inFile = argv[3];


 
  TRint *theApp = new TRint( "theApp", &argc, argv );
  //TROOT root( "GUI", "GUI" );
  //TApplication theApp( "App", &argc, argv );
  //gROOT->SetStyle( "Plain" );
  gROOT->cd();

  int runnum=99999;
  if(argc>1) runnum=atoi(argv[1]);

  
  gSystem->Load("libPhysics.so");
  gSystem->Load( "./lib/libAll.so" );
  //#######Need knucl lib file!! Set your dir!! ###
  gSystem->Load("~/work/geant/knucl3/libKnuclRootData.so");
  //#######################################


  
  //  std::string conffile;
  //  conffile="conf/Oct2010/analyzerBLC.conf";
//   if(3057<=runnum && runnum<=3060)
//     conffile="conf/Oct2010/analyzer3057-3060.conf";
//   else  if(3061<=runnum && runnum<=3070)
//     conffile="conf/Oct2010/analyzer3061-3070.conf";
//   else  if(3071<=runnum && runnum<=3080)
//     conffile="conf/Oct2010/analyzer3071-3080.conf";
//   else  if(3081<=runnum && runnum<=3090)
//     conffile="conf/Oct2010/analyzer3081-3090.conf";
//   else  if(3091<=runnum && runnum<=3108)
//     conffile="conf/Oct2010/analyzer3091-3100.conf";
//   else 
//     conffile="conf/Oct2010/analyzer3057-3060.conf";

  ConfMan *conf = new ConfMan(confFile);


  conf->Initialize();


  TFile *f;
  TTree *evtree;

  //f[i] = new TFile( Form( "./root/track_%d_%d.root", runnum,run[i] ) );
  //       f[i] = new TFile( Form( "./root/track_re_%d.root", runnum ) );
// #if SIM
//   f = new TFile( Form( "./simdata/sim2_pi018.root" ) );
//   //  f = new TFile( Form( "./simout.root" ) );
// #else
//   if(argc==1) f = new TFile( Form( "./root/track_3057_3107.root", runnum ) );
//   // if(argc==1) f = new TFile( Form( "./root/track_piall.root", runnum ) );
// //   else f = new TFile( Form( "./root/track_%s_newBLDC.root", argv[1] ) );
// //  if(argc==1) f = new TFile( Form( "./root/track_3057_3107.root", runnum ) );
//   else f = new TFile( Form( "./root/track_%s.root", argv[1] ) );
// #endif
  f = new TFile( inFile.c_str());
  evtree = (TTree*)f->Get( "EventTree" );

  

  
  TFile *fout;


  //#if SIM
//   //  fout = new TFile( Form("./root/trout_sim_pi2.root",runnum), "recreate" );
//   fout = new TFile( Form("./trout_sim.root",runnum), "recreate" );
// #else 
//   //  if(argc==1) fout = new TFile( Form("./root/trout_piall.root",runnum), "recreate" );
//   if(argc==1) fout = new TFile( Form("./root/trout_3057_3107.root",runnum), "recreate" );
//   else fout = new TFile( Form("./root/trout_tmp_%s.root",argv[1]), "recreate" );
// #endif
  fout = new TFile( outFile.c_str(), "recreate" );
  
  new TH1F( "nTrack", "nTrack", 200, 0, 100 );
  new TH1F( "chi", "chi", 300, 0, 30 );
  for(int seg=0;seg<5;seg++)  new TH1F( Form("timeT0_%d",seg), Form("time of T0 %d",seg), 100, -5, 5 );
  new TH1F( "mom", "mom", 500, -1, 1 );
  new TH2F( "pid", "pid", 200, -5, 15, 300, -1.5, 1.5 );
  new TH1F( "phi", "phi of track", 200, 0, 360 );
  new TH1F( "momp", "mom proton", 500, -1, 1 );
  new TH1F( "mompip", "mom piplus", 500, -1, 1 );
  new TH1F( "mompim", "mom piminus", 500, -1, 1 );
  new TH1F( "vx", "vx", 500, -25, 25 );
  new TH1F( "vy", "vy", 500, -25, 25 );
  new TH1F( "vz", "vz", 1000, -50, 50 );
  new TH2F( "vxy", "vxy", 500, -25, 25 , 500, -25, 25 );
  new TH1F( "vx1", "vx1", 500, -25, 25 );
  new TH1F( "vy1", "vy1", 500, -25, 25 );
  new TH2F( "vxy1", "vxy1", 500, -25, 25 , 500, -25, 25 );
  new TH1F( "vx2", "vx2", 500, -25, 25 );
  new TH1F( "vy2", "vy2", 500, -25, 25 );
  new TH2F( "vxy2", "vxy2", 500, -25, 25 , 500, -25, 25 );
  new TH1F( "vx3", "vx3", 500, -25, 25 );
  new TH1F( "vy3", "vy3", 500, -25, 25 );
  new TH2F( "vxy3", "vxy3", 500, -25, 25 , 500, -25, 25 );
  new TH1F( "vrx_dis", "vrx_dis", 20000, 0, 4 );


  new TH1F( "chi_BLC", "chi_BLC", 3000, 0, 30 );
  new TH1F( "chixz_BLC", "chixz_BLC", 300, 0, 30 );
  new TH1F( "chiyz_BLC", "chiyz_BLC", 300, 0, 30 );
  new TH1F( "Beamx", "Beamx", 6000, -30, 30 );
  new TH1F( "Beamy", "Beamy", 6000, -30, 30 );
  new TH1F( "Beamdx", "Beamdx", 10000, -0.05, 0.05 );
  new TH1F( "Beamdy", "Beamdy", 10000, -0.05, 0.05 );
  new TH1F( "chiyz_BLC", "chiyz_BLC", 300, 0, 30 );
  new TH1F( "chiyz_BLC", "chiyz_BLC", 300, 0, 30 );
  new TH1F( "tofBHDT0", "tof BHD and T0", 200, -5, 35 );



  new TH1F( "vvrx_dis", "vvrx_dis", 2000, 0, 2 );
  new TH1F( "vvrx_cdc", "vvrx_cdc", 2000, 0, 2 );
  new TH1F( "vvrx_cdcx", "vvrx_cdcx", 2000, -1, 1 );
  new TH1F( "vvrx_cdcy", "vvrx_cdcy", 2000, -1, 1 );
  new TH1F( "vvrx_cdcz", "vvrx_cdcz", 2000, -1, 1 );
  new TH1F( "vvrx_cdc", "vvrx_cdc", 2000, 0, 2 );
  new TH1F( "vvrx_blc", "vvrx_blc", 2000, 0, 2 );
  new TH1F( "vvrx_blc2", "vvrx_blc2", 2000, 0, 2 );


  new TH1F( "vbx2", "vbx2", 500, -25, 25 );
  new TH1F( "vby2", "vby2", 500, -25, 25 );
  new TH1F( "vbz2", "vbz2", 1000, -50, 50 );
  new TH2F( "vbxy2", "vbxy2", 500, -25, 25 , 500, -25, 25 );
  new TH1F( "vrxb_dis2", "vrxb_dis2", 2000, 0, 10 );
  new TH2F( "bxy2", "bxy2", 1000, -10, 10 , 1000, -10, 10 );
  new TH1F( "v_vb_dis", "v_vb_dis", 2000, 0, 10 );


  new TH2F( "vx_vbx1", "vx_vbx1", 300, -15, 15, 300, -15, 15 );
  new TH2F( "vy_vby1", "vy_vby1", 300, -15, 15, 300, -15, 15 );
  new TH2F( "vx_vbx2", "vx_vbx2", 300, -15, 15, 300, -15, 15 );
  new TH2F( "vy_vby2", "vy_vby2", 300, -15, 15, 300, -15, 15 );
  new TH2F( "vx_vbx3", "vx_vbx3", 300, -15, 15, 300, -15, 15 );
  new TH2F( "vy_vby3", "vy_vby3", 300, -15, 15, 300, -15, 15 );


  new TH1F( "vbx", "vbx", 500, -25, 25 );
  new TH1F( "vby", "vby", 500, -25, 25 );
  new TH1F( "vbz", "vbz", 1000, -50, 50 );
  new TH2F( "vbxy", "vbxy", 500, -25, 25 , 500, -25, 25 );
  new TH1F( "vrxb_dis", "vrxb_dis", 5000, 0, 10 );
  new TH2F( "vrxb_dis_chi", "vrxb_dis_chi", 2000, 0, 10,1000,0,10 );
  for(int n=0;n<3;n++)
    {
      new TH1F(Form( "vrxb_disx%d",n), Form( "vrxb_disx%d",n), 1200, -6, 6 );
      new TH1F(Form( "vrxb_disy%d",n),Form( "vrxb_disy%d",n), 1200, -6, 6 );
      new TH1F(Form( "vrxb_disz%d",n),Form( "vrxb_disz%d",n), 1200, -6, 6 );
    }

  new TH2F( "vbx_dis", "vbx_dis", 100, -10.0, 10.0 , 100, -10, 10 );
  new TH2F( "vby_dis", "vby_dis", 100, -10.0, 10.0 , 100, -10, 10 );

  new TH2F( "BLCvtx_z_x", "BLCvtx_z_x", 400, -50.0, 50.0 , 80, -10, 10 );
  new TH2F( "BLCvtx_z_y", "BLCvtx_z_y", 400, -50.0, 50.0 , 80, -10, 10 );
  new TH2F( "BLCvtx_r_x", "BLCvtx_r_x", 90, 0.0, 30.0 , 80, -10, 10 );
  new TH2F( "BLCvtx_r_y", "BLCvtx_r_y", 90, 0.0, 30.0 , 80, -10, 10 );
  new TH2F( "BLCvtx_r_theta", "BLCvtx_r_theta", 90, 0.0, 30.0 , 80, -1.0, 1.0 );
  for(int n=0;n<3;n++)
    {
      new TH1F( Form("BLCvtx_x%d",n+1), Form("BLCvtx_x%d",n+1), 80, -10, 10 );
      new TH1F( Form("BLCvtx_y%d",n+1), Form("BLCvtx_y%d",n+1), 80, -10, 10 );
    }

  new TH2F( "vb_tanx", "vb_tanx", 100, -0.04, 0.04 , 100, -10, 10 );
  new TH2F( "vb_tany", "vb_tany", 100, -0.04, 0.04 , 100, -10, 10 );
  new TH2F( "bxy", "bxy", 1000, -10, 10 , 1000, -10, 10 );
  new TH2F( "bdxy", "bdxy", 1000, -0.1, 0.1 , 1000, -0.1, 0.1 );

  for(int n=0;n<5;n++)
  new TH1F( Form("BLCpos_T0_%d",n), Form("BLCpos_T0_%d",n) ,200,-10,10);

  new TH2F( "mass_mom", "mass_mom", 360, 0.0 , 1.2, 240,-1.2,1.2);
  new TH2F( "mass2_mom", "mass2_mom", 630, -0.1 , 2.0, 240,-1.2,1.2);


  new TH2F( "CDH_calc_diff_proton", "CDH_calc_diff_proton", 400, 0.0 , 20, 120,-3,3);
  new TH2F( "CDH_E_diff_proton", "CDH_E_diff_proton", 200, 0.0 , 100, 120,-3,3);
  
  new TH2F( Form("CDH_sqrtE_diff"), Form("CDH_sqrtE_diff"), 200, 0.0 , 1, 400,-10,10);
  new TH2F( Form("CDH_sqrtE_diff2"), Form("CDH_sqrtE_diff2"), 200, 0.0 , 1, 400,-10,10);

  //#####Momentum####
  for(int p=1;p<=4;p++)
    {
      new TH2F( Form("CDH_sqrtE_diff_p%d",p), Form("CDH_sqrtE_diff_p%d",p), 200, 0.0 , 1, 400,-10,10);
      new TH2F( Form("CDH_sqrtE_diff2_p%d",p), Form("CDH_sqrtE_diff2_p%d",p), 200, 0.0 , 1, 400,-10,10);
    }
  

  for(int seg=1;seg<=36;seg++)
    { 
  new TH2F( Form("CDHu%d_ADC_mom",seg), Form("CDHu%d_ADC_mom",seg), 4000,0,4000 , 200, 0.0 , 1);
  new TH2F( Form("CDHd%d_ADC_mom",seg), Form("CDHd%d_ADC_mom",seg), 4000,0,4000 , 200, 0.0 , 1);
  new TH2F( Form("CDH%d_ADCmean_mom_pi",seg), Form("CDH%d_ADCmean_mom_pi",seg), 1000,-5,45 , 200, 0.0 , 1);
  new TH2F( Form("CDH%d_ADCmean_mom_p",seg), Form("CDH%d_ADCmean_mom_p",seg), 1000,-5,45, 200, 0.0, 1.0);
  new TH2F( Form("CDHu%d_sqrtE_diff",seg), Form("CDHu%d_sqrtE_diff",seg), 200, 0.0 , 1, 400,-10,10);
  new TH2F( Form("CDHd%d_sqrtE_diff",seg), Form("CDHd%d_sqrtE_diff",seg), 200, 0.0 , 1, 400,-10,10);

  new TH2F( Form("CDHu%d_sqrtE_diff2",seg), Form("CDHu%d_sqrtE_diff2",seg), 200, 0.0 , 1, 400,-10,10);
  new TH2F( Form("CDHd%d_sqrtE_diff2",seg), Form("CDHd%d_sqrtE_diff2",seg), 200, 0.0 , 1, 400,-10,10);

  //#####Momentum####
  for(int p=1;p<=4;p++)
    {
      new TH2F( Form("CDHu%d_sqrtE_diff_p%d",seg,p), Form("CDHu%d_sqrtE_diff_p%d",seg,p), 200, 0.0 , 1, 400,-10,10);
      new TH2F( Form("CDHd%d_sqrtE_diff_p%d",seg,p), Form("CDHd%d_sqrtE_diff_p%d",seg,p), 200, 0.0 , 1, 400,-10,10);      
      new TH2F( Form("CDHu%d_sqrtE_diff2_p%d",seg,p), Form("CDHu%d_sqrtE_diff2_p%d",seg,p), 200, 0.0 , 1, 400,-10,10);
      new TH2F( Form("CDHd%d_sqrtE_diff2_p%d",seg,p), Form("CDHd%d_sqrtE_diff2_p%d",seg,p), 200, 0.0 , 1, 400,-10,10);
    }

  //###############

      new TH2F(Form("CDH%d_calc_diff_proton",seg),Form("CDH%d_calc_diff_proton",seg),400,0,20,400,-10,10);
  new TH2F( Form("CDH%d_E_diff_proton",seg), Form("CDH%d_E_diff_proton",seg), 200, 0.0 , 100, 400,-10,10);
  new TH2F( Form("CDHu%d_E_diff_proton",seg), Form("CDHu%d_E_diff_proton",seg), 200, 0.0 , 100, 400,-10,10);
  new TH2F( Form("CDHd%d_E_diff_proton",seg), Form("CDHd%d_E_diff_proton",seg), 200, 0.0 , 100, 400,-10,10);

    }

  //#########t0 #####
  for(int seg=0;seg<5;seg++)
    {
      new TH2F( Form("T0%d_sqrtE_diff",seg), Form("T0%d_E_diff",seg), 200, 0.0 , 1.0, 400,-10,10);
      new TH2F( Form("T0u%d_sqrtE_diff",seg), Form("T0u%d_E_diff",seg), 200, 0.0 , 1.0, 400,-10,10);
      new TH2F( Form("T0d%d_sqrtE_diff",seg), Form("T0d%d_E_diff",seg), 200, 0.0 , 1.0, 400,-10,10);
    }

  new TH1F( "vxppi", "vxppi", 500, -25, 25 );
  new TH1F( "vyppi", "vyppi", 500, -25, 25 );
  new TH1F( "vzppi", "vzppi", 1000, -50, 50 );
  new TH3F( "vppi",  "vppi", 200, -25, 25, 200, -25, 25, 200, -25, 25 );
  new TH1F( "IMppi", "IM proton piminus", 500, 1.0, 2 );
  new TH1F( "IMppi_m", "IM proton piminus p*0.99", 500, 1.0, 2 );
  new TH1F( "IMppi_p", "IM proton piminus p*1.01", 500, 1.0, 2 );
  new TH1F( "IMppiKselected", "IM proton piminus", 500, 1.0, 2 );
  new TH1F( "cosOAppi", "cosOA between proton and #pi^-", 200, -1, 1 );
  new TH2F( "IMcosOAppi", "IM and cosOA for ppi", 200, 1, 2, 200, -1, 1);
  new TH1F( "costppi", "cost of ppi system", 200, -1, 1 );
  new TH1F( "phippi", "phi of ppi system", 200, -180, 180 );
  new TH2F( "IMcostppi", "IM and cost of ppi system", 500, 1, 2, 200, -1, 1 );


  //   new TH3F( "vLambda", "vLambda", 1000, -25, 25,1000, -25, 25,1000, -25, 25 );
   new TH1F( "v_disLambda", "v_disLambda", 500, 0, 5 );
   new TH1F( "pLambda", "pLambda", 1000, 0, 1.0 );
   new TH1F( "disLambda", "disLambda", 2000, 0, 10 );
   new TH1F( "timeLambda", "timeLambda", 1000, 0, 1.0 );
   new TH1F( "missingMassLambda", "missingMassLambda", 1500, -0.5, 1.0 );

  new TH1F( "vxpipi", "vxpipi", 500, -25, 25 );
  new TH1F( "vypipi", "vypipi", 500, -25, 25 );
  new TH1F( "vzpipi", "vzpipi", 1000, -50, 50 );



  new TH3F( "vpipi",  "vpipi", 200, -25, 25, 200, -25, 25, 200, -25, 25 );
  new TH1F( "IMpipi", "IM piplus piminus", 500, 0, 1.0 );
  new TH1F( "IMpipi_m", "IM piplus piminus p*0.99", 500, 0, 1.0 );
  new TH1F( "IMpipi_p", "IM piplus piminus p*1.01", 500, 0, 1.0 );
  new TH1F( "IMpipi_cut", "IM piplus piminus cut OA", 500, 0, 1.0 );
  new TH1F( "IMpipi_nopid", "IM piplus piminus nopid", 500, 0, 1.0 );
  new TH1F( "cosOApipi", "cosOA between #pi^+ and #pi^-", 200, -1, 1 );
  new TH2F( "IMcosOApipi", "IM and cosOA for pipi", 500, 0, 1, 200, -1, 1);
  new TH1F( "costpipi", "cost of pipi system", 200, -1, 1);
  new TH1F( "phipipi", "phi of pipi system", 200, -180, 180 );
  new TH2F( "IMcostpipi", "IM and cost of pipi system", 500, 0, 1, 200, -1, 1 );
  

  TH1F *h1;
  TH2F *h2;
  TH3F *h3;


  CDSHitMan *cdsMan = 0;
  BeamLineHitMan *blMan = 0;
  BeamLineTrackMan *blTrackMan = 0;
  EventHeader *head = 0;
  CDSTrackingMan *trackMan = 0;
  MCData *mcdata=0;

  bool AllData=false;
  bool mcflag=false;

  TObjArray *objarr=evtree->GetListOfBranches();
  if(objarr->FindObject("CDSHitMan")!=0 && objarr->FindObject("CDSTrackingMan")!=0)
    {
      evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
      evtree->SetBranchAddress( "CDSTrackingMan", &trackMan );
      if(objarr->FindObject("BeamLineHitMan")!=0 && 
	 objarr->FindObject("BeamLineTrackMan")!=0 )
      {
	evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
	evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
	evtree->SetBranchAddress( "BeamLineTrackMan", &blTrackMan );
	evtree->SetBranchAddress( "CDSTrackingMan", &trackMan );

      }
      if(objarr->FindObject("EventHeader")!=0)
	{
	  evtree->SetBranchAddress( "EventHeader", &head );
	  AllData=true;
	}
      else if(objarr->FindObject("MCData")!=0) 
	{
	  evtree->SetBranchAddress( "MCData", &mcdata );
	  mcflag=true;
	}
    }
  else 
    {
      std::cout<<"error of EventTree branch"<<std::endl;
      return 0;
    }


  int nev = evtree->GetEntries();
  std::cout <<" entry:" << nev << std::endl;

  //###test for trackout###
  ofstream of;
  //  std::string OutFileName=Form("trackdata.txt");
  of.open("trackdata.txt");
  //  if( of.open("trackdata.txt")==0 )
  //{
      //  std::cerr << " File open fail. [" << OutFileName << "]" << std::endl;
  //  return 0;
  // }


  for( int iev=0; iev<nev; iev++ ){
    //    if(iev>10000) continue;

    if( iev%2000 == 0 )
      std::cout << " Event : " << iev << std::endl;
    evtree->GetEvent(iev);


    //    blMan->CheckContainerSize();
    //  blMan->Calc(conf);
    //blTrackMan->
    //blTrackMan->DoTracking(blMan);
    //    cdsMan->Calc(conf);
    //cdsMan->CheckContainerSize();
    //  
    /*
    std::cout<<"nGoodTrack "<<trackMan->nGoodTrack()<<std::endl;
   for(int n=0;n<trackMan->nGoodTrack();n++)
      {
	trackMan->CalcVertex_beam(trackMan->GoodTrackID(n),blTrackMan,conf);
      }
    */


    double ctmT0=0,euT0,edT0,emT0;

    int pid_beam;//0:pi 1:K 3:else
    int segT0;
    if(AllData)
	{
	  int nT0=0;
	  for( int i=0; i<blMan->nT0(); i++ ){
	    if( blMan->T0(i)->CheckRange() ) nT0++;
	  }
	  if( nT0!=1 ) ctmT0=0;//continue;


	  for( int i=0; i<blMan->nT0(); i++ ){
	    if( blMan->T0(i)->CheckRange() ){
	      //  blMan->T0(i)->Calc(conf);
#if SIM
	      ctmT0 = blMan->T0(i)->tmean();
#else
	      ctmT0 = blMan->T0(i)->ctmean();
#endif
	
	      euT0 = blMan->T0(i)->eu();
	      edT0 = blMan->T0(i)->ed();
	      emT0 = blMan->T0(i)->emean();
	      segT0=i;
	    }
	  }

	  int GoodTrack=0;
	  for( int it=0; it<trackMan->nTrack(); it++ ){
	    CDSTrack *track = trackMan->Track(it);
	    double chi = track->Chi();
	    if(chi<20) GoodTrack++;
	  }
	  
	  int nBHD=0;
	  for( int i=0; i<blMan->nBHD(); i++ ){
	    if( blMan->BHD(i)->CheckRange() ) nBHD++;
	  }
	  //       if( nBHD!=1 ) continue;
	  
	  double ctmBHD=0;
	  for( int i=0; i<blMan->nBHD(); i++ ){
	    if( blMan->BHD(i)->CheckRange() ){
#if SIM
	      ctmBHD = blMan->BHD(i)->tmean();
#else
	      ctmBHD = blMan->BHD(i)->ctmean();
#endif

	    }
	  }

// 	  double tofBHDT0 = ctmT0-ctmBHD;
	  
// 	  h1 = (TH1F*)gFile->Get("tofBHDT0");      
// 	  h1->Fill( tofBHDT0 );
// 	  h1 = (TH1F*)gFile->Get(Form("timeT0_%d",segT0));h1->Fill(ctmT0);
	  
       
// 	  if(24<tofBHDT0 && tofBHDT0<28) pid_beam=0;
// 	  else if(28<tofBHDT0 && tofBHDT0<32) pid_beam=1;
// 	  else pid_beam=3;
// 	  //if(pid_beam!=0) continue;



	}
    else
      {
       pid_beam=0;//0:pi 1:K 3:else
       segT0=1; 
       for( int i=0; i<blMan->nT0(); i++ ){
	 if( blMan->T0(i)->CheckRange() ){
	   //  blMan->T0(i)->Calc(conf);
#if SIM
	      ctmT0 = blMan->T0(i)->tmean();
#else
	      ctmT0 = blMan->T0(i)->ctmean();
#endif
	      
	      euT0 = blMan->T0(i)->eu();
	      edT0 = blMan->T0(i)->ed();
	      emT0 = blMan->T0(i)->emean();
	      segT0=i;
	      h1 = (TH1F*)gFile->Get(Form("timeT0_%d",segT0));h1->Fill(ctmT0);
	 }
       }

      }


    //###BLC Chi and pos T0#############

    //    if(blTrackMan->ntrackBLC2()==1)
      {
	for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
	  {
	    LocalTrack *blc=blTrackMan->trackBLC2(itr);
	    double x1,y1,x2,y2;
	    double gx,gy,gz,dx,dy,dz;
	    conf->GetGeomMapManager()->GetGParam(CID_T0,gx,gy,gz,dx,dy,dz);
	    blc->XYPosatZ(gz,x1,y1);		  
	    h1 = (TH1F*)gFile->Get(Form("BLCpos_T0_%d",segT0) ); h1->Fill(x1);
	    h1 = (TH1F*)gFile->Get(Form("chi_BLC") ); h1->Fill(blc->chi2all());	
	    h1 = (TH1F*)gFile->Get(Form("chixz_BLC") ); h1->Fill(blc->chi2xz());	
	    h1 = (TH1F*)gFile->Get(Form("chiyz_BLC") ); h1->Fill(blc->chi2yz());	

	    blc->XYPosatZ(0,x1,y1);
	    blc->XYPosatZ(10,x2,y2);
	    TVector3 dir;
	    dir.SetXYZ(x2-x1,y2-y1,10);
	    dir=dir.Unit();

	    //	    if(pid_beam==1)//0:pi 1:K 3:else
	      {
		h1 = (TH1F*)gFile->Get(Form("Beamx") ); h1->Fill(x1);		  
		h1 = (TH1F*)gFile->Get(Form("Beamy") ); h1->Fill(y1);
		h1 = (TH1F*)gFile->Get(Form("Beamdx") ); h1->Fill(asin(dir.x()) );		  
		h1 = (TH1F*)gFile->Get(Form("Beamdy") ); h1->Fill(asin(dir.y()) );		  		  
	      }
	  }
      }
      h1 = (TH1F*)gFile->Get("nTrack"); h1->Fill( trackMan->nGoodTrack() );

	//#####BLC#########//

	for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
	  {
	    LocalTrack *blc2=blTrackMan->trackBLC2(itr);
	    double x1,y1,x2,y2;
	    double z=0,z2=20;
	    blc2->XYPosatZ(z,x1,y1);		  
	    blc2->XYPosatZ(z2,x2,y2);
	    double len=sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z)*(z2-z) );
	    double dx=(x2-x1)/len;   double dy=(y2-y1)/len;		 
	    h2 = (TH2F*)gFile->Get("bxy"); h2->Fill( x1,y1 );
	    h2 = (TH2F*)gFile->Get("bdxy"); h2->Fill( dx,dy );
	  }




      //##########################
      //####CDSTrack Single#######
      //#########################
      int trackdatanum=0;
      for( int it=0; it<trackMan->nGoodTrack(); it++ ){
	CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
	double chi = track->Chi();
	h1 = (TH1F*)gFile->Get("chi"); h1->Fill( chi );
	HodoscopeLikeHit *cdh = track->CDHHit();
	if(!track->CDHFlag() ) continue;
#if SIM
	double tof = cdh->tmean()-ctmT0;
#else
	double tof = cdh->ctmean()-ctmT0;
#endif

	double param[5];
	track->GetParameters(param);
	//	if(param[2]==0 ||param[4]==0   ) continue;
	double drho=param[0], phi0=param[1], rho=track->Rho(), dz=param[3], tlam=param[4];
	double mom =track->Momentum();
	double x0,y0,z0;
	z0=dz;
	track->XYatZ(x0,y0,z0);

	h1 = (TH1F*)gFile->Get("mom"); h1->Fill( mom );
	h2 = (TH2F*)gFile->Get("pid"); h2->Fill( tof, mom );
	int ptype =-1;
	if(mom!=0) ptype = PID( tof, mom );
	if( ptype==1 ){ h1 = (TH1F*)gFile->Get("mompip"); h1->Fill( mom ); }
	if( ptype==2 ){ h1 = (TH1F*)gFile->Get("mompim"); h1->Fill( mom ); }
	if( ptype==3 ){ h1 = (TH1F*)gFile->Get("momp"); h1->Fill( mom ); }

	double phi=track->Theta();
	h1 = (TH1F*)gFile->Get("phi"); h1->Fill( phi);
	int pdgID;
	if(ptype==1) pdgID=211;
	else if(ptype==2) pdgID=-211;
	else if(ptype==3) pdgID=2212;
	else continue;
	TVector3 pos0,mom0;
	pos0.SetXYZ(x0,y0,z0);
	if(!track->GetMomentum(pos0,mom0)) continue;
	//	of<<pdgID<<" "<<0<<" "<<pos0.x()<<" "<<pos0.y()<<" "<<pos0.z()<<" "<<mom0.x()<<" "<<mom0.y()<<" "<<mom0.z()<<std::endl;
	trackdatanum++;
      }

      //####test for trackdata###
      of<<iev<<" "<<trackdatanum<<std::endl;
      for( int it=0; it<trackMan->nGoodTrack(); it++ ){
	CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
	double chi = track->Chi();
	h1 = (TH1F*)gFile->Get("chi"); h1->Fill( chi );
	HodoscopeLikeHit *cdh = track->CDHHit();
	if(!track->CDHFlag() ) continue;
#if SIM
	double tof = cdh->tmean()-ctmT0;
#else
	double tof = cdh->ctmean()-ctmT0;
#endif

	double param[5];
	track->GetParameters(param);
	//	if(param[2]==0 ||param[4]==0   ) continue;
	double drho=param[0], phi0=param[1], rho=track->Rho(), dz=param[3], tlam=param[4];
	double mom =track->Momentum();
	double x0,y0,z0;
	z0=dz;
	track->XYatZ(x0,y0,z0);

	h1 = (TH1F*)gFile->Get("mom"); h1->Fill( mom );
	h2 = (TH2F*)gFile->Get("pid"); h2->Fill( tof, mom );
	int ptype = -1;
	if(mom!=0) ptype = PID( tof, mom );
	double phi=track->Theta();
	int pdgID;
	if(ptype==1) pdgID=211;
	else if(ptype==2) pdgID=-211;
	else if(ptype==3) pdgID=2212;
	else continue;
	TVector3 pos0,mom0;
	pos0.SetXYZ(x0,y0,z0);
	if(!track->GetMomentum(pos0,mom0)) continue;
	of<<pdgID<<" "<<0<<" "<<pos0.x()<<" "<<pos0.y()<<" "<<pos0.z()<<" "<<mom0.x()<<" "<<mom0.y()<<" "<<mom0.z()<<std::endl;
      }



      //###############
      //#CDS Track 2###
      //###############

      for( int it1=0; it1<trackMan->nGoodTrack(); it1++ ){

	CDSTrack *track1 = trackMan->Track(trackMan->GoodTrackID(it1));
	bool cdhflag=true;;
	if(!track1->CDHFlag()) cdhflag=false;
	HodoscopeLikeHit *cdh1;
	cdh1=track1->CDHHit();
	//	if(cdhflag) cdh1->Calc(conf);
	//	if(cdh1->seg()==-1) cdhflag=false;
#if SIM
	double tof1 = cdh1->tmean()-ctmT0;
	double tof1_u = cdh1->tu()-ctmT0;
	double tof1_d = cdh1->td()-ctmT0;
#else
	double tof1 = cdh1->ctmean()-ctmT0;
	double tof1_u = cdh1->ctu()-ctmT0;
	double tof1_d = cdh1->ctd()-ctmT0;
#endif
	double param1[5];

	track1->GetParameters(param1);
	//	if(param1[2]==0  ) continue;
	double drho1=param1[0], phi01=param1[1], rho1=track1->Rho(), dz1=param1[3], tlam1=param1[4];
	double mom1 = track1->Momentum();
	int ptype1 = -1;
	if(mom1!=0 )ptype1 = PID( tof1, mom1 );
	if(cdhflag)
	  {
	h2 = (TH2F*)gFile->Get(Form("CDHu%d_ADC_mom",cdh1->seg() )); h2->Fill(cdh1->adcu(),fabs(mom1) );
	h2 = (TH2F*)gFile->Get(Form("CDHd%d_ADC_mom",cdh1->seg() )); h2->Fill(cdh1->adcd(),fabs(mom1) );
	  }
	if(ptype1==0 ) cdhflag=false;
	if(cdhflag &&( ptype1==1 || ptype1==2))
	  {
	h2 = (TH2F*)gFile->Get(Form("CDH%d_ADCmean_mom_pi",cdh1->seg() )); h2->Fill(cdh1->emean(),fabs(mom1) );
	  }
	if(cdhflag && ptype1==3 )
	  {
	h2 = (TH2F*)gFile->Get(Form("CDH%d_ADCmean_mom_p",cdh1->seg() )); h2->Fill(cdh1->emean(),mom1 );
	  }
	//#######vertex beam#####	  
	TrackVertex vertexb;
	TVector3 vtxb1, vtxb2, vtxb;

	//	if(blTrackMan->ntrackBLC()!=1 ) continue;
	if(trackMan->GetVertex_beam(trackMan->GoodTrackID(it1),0,vertexb) )
	  {	 
	    vtxb=vertexb.GetVertexPos_mean();
	    vtxb1=vertexb.GetVertexPos1();
	    vtxb2=vertexb.GetVertexPos2();
	    h1 = (TH1F*)gFile->Get("vbx"); h1->Fill( vtxb.X() );
	    h1 = (TH1F*)gFile->Get("vby"); h1->Fill( vtxb.Y() );
	    h1 = (TH1F*)gFile->Get("vbxy"); h1->Fill( vtxb.X(),vtxb.Y() );
	    h1 = (TH1F*)gFile->Get("vbz"); h1->Fill( vtxb.Z() );
	    TVector3 tmpbv=vtxb2-vtxb1;
	    double dis=tmpbv.Mag();
	    //	    h1 = (TH1F*)gFile->Get("vrxb_dis"); h1->Fill( dis );
	    h2 = (TH2F*)gFile->Get("vrxb_dis_chi"); h2->Fill( dis,track1->Chi() );
	    if(track1->Chi()<0.5)
	      h1 = (TH1F*)gFile->Get("vrxb_dis"); h1->Fill( dis );

	    //############MC ##############//
	    if(mcflag)
	      {
		std::cout<<"mctrack :"<<mcdata->trackSize()<<std::endl;
		for(int vt=0;vt<mcdata->trackSize();vt++)
		  {
		    Track *vtrack=mcdata->track(vt);
		    //		if(vtrack->momentum()==900) continue;

		    TVector3 vvertex =vtrack->vertex();
		    vvertex=0.1*vvertex;
		    //   vvertex.SetXYZ(0,0,83.88/2.0);
		    //		    if(vt==0) vvertex.SetXYZ(0,0,0);
		    //		std::cout<<"Mom :vertex : "<<std::endl;
		    std::cout<<"mom1 :Mom :vertex : "<<mom1<<" "<<vtrack->momentum().Mag()<<" "<<vvertex.x()<<" "<< vvertex.y()<<" "<< vvertex.z()<<" "<<std::endl; 
		    TVector3 tmp=vtxb-vvertex;
		    h1 = (TH1F*)gFile->Get("vvrx_dis"); h1->Fill( tmp.Mag() );
		    if(vt==0)
		      {
 			TVector3 fitpos;
 			double cdcdis;
// 			track1->PointToHelix(vvertex,param1,fitpos,cdcdis);
// 			double paramtmp[5],xest,yest,distmp;
// 			track1->GetTmpParameters(2,paramtmp);

// 			track1->PointToCircle(vvertex.x(),vvertex.y(),
// 					      fabs(1./paramtmp[2]),(paramtmp[0]+1./paramtmp[2])*cos(paramtmp[1]),(paramtmp[0]+1./paramtmp[2])*sin(paramtmp[1]),
// 												  distmp,xest,yest);
			if(param1[2]==0  )
			  {
			    TVector3 dline,oline;
			    oline=track1->LineOrigin();    dline=track1->LineDirection();
			    track1->PointToLine(vvertex,oline,dline,cdcdis,fitpos);
			  }
			else 
			  {
			    track1->PointToHelix(vvertex,param1,fitpos,cdcdis);
			  }
			h1 = (TH1F*)gFile->Get("vvrx_cdc"); h1->Fill( cdcdis );
			h1 = (TH1F*)gFile->Get("vvrx_cdcx"); h1->Fill( fitpos.x()-vvertex.x() );
			h1 = (TH1F*)gFile->Get("vvrx_cdcy"); h1->Fill( fitpos.y()-vvertex.y() );
			h1 = (TH1F*)gFile->Get("vvrx_cdcz"); h1->Fill( fitpos.z()-vvertex.z() );
		      }
		    else if(vt==1)
		      {
			LocalTrack *blctrack=blTrackMan->trackBLC(0);
			double bx[2],by[2];
			blctrack->XYPosatZ(vvertex.z() ,bx[1],by[1]);
			blctrack->XYPosatZ(vvertex.z()-10 ,bx[0],by[0]);
			TVector3 blcline,blcp;
			blcp.SetXYZ(bx[1],by[1],vvertex.z());
			blcline.SetXYZ(bx[1]-bx[0],by[1]-by[0],10);
			blcline=blcline.Unit();
			TVector3 nest;
			double blcdis;
			track1->PointToLine(vvertex,blcp,blcline,blcdis,nest);
			h1 = (TH1F*)gFile->Get("vvrx_blc"); h1->Fill( blcdis );
			double blcdis2=sqrt( (bx[1]-vvertex.x() )*(bx[1]-vvertex.x() )+(by[1]-vvertex.y() )*(by[1]-vvertex.y() ) );
			h1 = (TH1F*)gFile->Get("vvrx_blc2"); h1->Fill( blcdis2 );
		      }
		  }
	      }
	    //##############################//
	    double dx=(vtxb2.x()-vtxb1.x());
	    double dy=(vtxb2.y()-vtxb1.y());
	    h2 = (TH2F*)gFile->Get("vbx_dis"); h2->Fill( vtxb2.x(),tmpbv.x() );	   
	    h2 = (TH2F*)gFile->Get("vby_dis"); h2->Fill( vtxb2.y(),tmpbv.y() );
	    double tanx,tany;
	    if(dis>2.0) cdhflag=false;
 	    for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
	      {
		LocalTrack *blc2=blTrackMan->trackBLC(itr);
		double x1,x2,y1,y2;

		blc2->XYPosatZ(10,x2,y2);
		blc2->XYPosatZ(0,x1,y1);		  
		tanx=(x2-x1)/10.; tany=(y2-y1)/10.;
	      }

	    h2 = (TH2F*)gFile->Get("vb_tanx"); h2->Fill( tanx,dx );
	    h2 = (TH2F*)gFile->Get("vb_tany"); h2->Fill( tany,dy );

	    //	    h1 = (TH1F*)gFile->Get("v_vb_dis"); h1->Fill( tmpbv.Mag() );	    
	    if(vtxb.Z()<-32 && vtxb.Z()>-47 ) 
	      {
		h1 = (TH1F*)gFile->Get("vrxb_disx0"); h1->Fill(tmpbv.x() );
		h1 = (TH1F*)gFile->Get("vrxb_disy0"); h1->Fill(tmpbv.y() );
		h1 = (TH1F*)gFile->Get("vrxb_disz0"); h1->Fill( tmpbv.z() );
		for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
		  {
		    LocalTrack *blc2=blTrackMan->trackBLC(itr);
		    double x1,y1;
		    double z=-42;
		    blc2->XYPosatZ(z,x1,y1);		  
		    double cdc_x,cdc_y;
		    track1->XYatZ(cdc_x,cdc_y,z);

		    //h2 = (TH2F*)gFile->Get("vy_vby1"); h2->Fill( cdc_y,y1 );
		  }		
	      }

	    else if(vtxb.Z()<-4 && vtxb.Z()>-6) 
	      {
		h1 = (TH1F*)gFile->Get("vrxb_disx1"); h1->Fill(tmpbv.x() );
		h1 = (TH1F*)gFile->Get("vrxb_disy1"); h1->Fill(tmpbv.y() );
		h1 = (TH1F*)gFile->Get("vrxb_disz1"); h1->Fill( tmpbv.z() );
		for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
		  {
		    LocalTrack *blc2=blTrackMan->trackBLC(itr);
		    double x1,y1;
		    double z=-5;
		    blc2->XYPosatZ(z,x1,y1);		  
		    double cdc_x,cdc_y;
		    track1->XYatZ(cdc_x,cdc_y,z);
		    //		    h2 = (TH2F*)gFile->Get("vx_vbx2"); h2->Fill( cdc_x,x1 );
		    //h2 = (TH2F*)gFile->Get("vy_vby2"); h2->Fill( cdc_y,y1 );
		  }		
	      }
	    else if(vtxb.Z()>4 && vtxb.Z()<6 ) 
	      {
		h1 = (TH1F*)gFile->Get("vrxb_disx2"); h1->Fill(tmpbv.x() );
		h1 = (TH1F*)gFile->Get("vrxb_disy2"); h1->Fill(tmpbv.y() );
		h1 = (TH1F*)gFile->Get("vrxb_disz2"); h1->Fill( tmpbv.z() );
		for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
		  {
		    LocalTrack *blc2=blTrackMan->trackBLC(itr);
		    double x1,y1;
		    double z=5;
		    blc2->XYPosatZ(z,x1,y1);		  
		    double cdc_x,cdc_y;
		    track1->XYatZ(cdc_x,cdc_y,z);
		    //h2 = (TH2F*)gFile->Get("vx_vbx3"); h2->Fill( cdc_x,x1 );
		    // h2 = (TH2F*)gFile->Get("vy_vby3"); h2->Fill( cdc_y,y1 );
		  }		
	      }
	    

	  }

	//############## CDH #####################//
	
	if(cdhflag && blTrackMan->ntrackBLC()==1 /*&& pid_beam==1*/) //CDH hit && vrx_b &&dis<2.0 && pid OK
	  {

	    //#####beam dis #######
	    // if( blTrackMan->ntrackBLC2()!=1) std::cout<<"BLC error!"<<std::endl; 
	    LocalTrack *blc2=blTrackMan->trackBLC(0);
	    double x1,y1,x2,y2;
	    double gx,gy,gz,dx,dy,dz;
	    //   if(blc2->chi2all()>1) continue;
	    conf->GetGeomMapManager()->GetGParam(CID_T0,gx,gy,gz,dx,dy,dz);
	    //	    std::cout<<"gz "<<gz<<std::endl; 
	    double z_t0=gz,z_vtxb=vtxb.Z();
	    blc2->XYPosatZ(z_t0,x1,y1);blc2->XYPosatZ(z_vtxb,x2,y2);
	    double beam_dis=sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z_t0-z_vtxb)*(z_t0-z_vtxb));

	    double beta_b;
	    if(AllData)
	      {
		if(pid_beam==0)  beta_b=0.9/sqrt(0.9*0.9+piMass*piMass);
		else if(pid_beam==1) beta_b=0.9/sqrt(0.9*0.9+kpMass*kpMass);
	      }
	    else    beta_b=0.9/sqrt(0.9*0.9+piMass*piMass);

	    double time_beam=beam_dis/beta_b/Const;
	    //double time_beam=0;//beam_dis/beta_b/Const;
	    //#####CDC dis#######		  

	    double mass1=0;	    
	    TVector3 cdhvtx=track1->CDHVertex();
	    double phi1,phi2;
	    phi1=track1->CalcHelixPhi(vtxb.x(),vtxb.y());
	    phi2=track1->CalcHelixPhi(cdhvtx.x(),cdhvtx.y());
	    TVector3 tmp=cdhvtx-vtxb;
	    //	    double cdc_dis=tmp.Mag();
	    double cdc_dis=rho1*fabs(phi2-phi1)*sqrt(1+tlam1*tlam1);
	    double beta_calc=cdc_dis/(tof1-time_beam)/Const;
	    double calc_mass1=fabs(mom1)*sqrt(1/(beta_calc*beta_calc)-1);
	    double calc_mass2=mom1*mom1*(1/(beta_calc*beta_calc)-1);
	    //  if(ptype1==1 || ptype1==2) mass1= piMass; // 1,2:pi ,3: proton
	    // else if(ptype1==3 ) mass1= pMass; // 1,2:pi ,3: proton
	    if(calc_mass1>0.07 && calc_mass1<0.22 ) mass1=piMass;
	    else  if(calc_mass1>0.6 && calc_mass1<1.2 ) mass1=pMass;

	    double beta=fabs(mom1)/sqrt(mass1*mass1+mom1*mom1);
	    double time_cdc=cdc_dis/beta/Const;


	    //####Calc time#####
	    double calc_time=time_cdc+time_beam;
	    double calc_time2=cdc_dis/Const+time_beam;

	    double ctsub=cdh1->ctsub();
	    double calc_time_u=time_cdc+time_beam+ctsub;
	    double calc_time_d=time_cdc+time_beam-ctsub;
	    //	    std::cout<<"CDH seg= "<<cdh1->seg()<<std::endl;	    
			      

	    h2 = (TH2F*)gFile->Get(Form("mass2_mom" )); h2->Fill(calc_mass2,mom1 );
	    h2 = (TH2F*)gFile->Get(Form("mass_mom" )); h2->Fill(calc_mass1,mom1 );

	    //   if( (calc_time-tof1)>2 && cdh1->emean()<10  ) continue;
// 	    std::cout<<"TOF : clac_time : calc_mas1 : "<<tof1<<" "<<calc_time<<" "<<calc_mass1<<std::endl;
// 	    std::cout<<"beam_time : beta_cdc : "<<time_beam<<" "<<beta_calc<<std::endl;
// 	    std::cout<<"Vertex  : "<<vtxb.x()<<" "<<vtxb.y()<<" "<<vtxb.z()<<" "<<std::endl;
	    if(mass1==piMass )
	      {

		if(fabs(mom1)>0.1)
		  { 				
		h2 = (TH2F*)gFile->Get(Form("CDH_sqrtE_diff" )); h2->Fill(1./sqrt(cdh1->emean()),calc_time-tof1 );
		h2 = (TH2F*)gFile->Get(Form("CDH_sqrtE_diff2" )); h2->Fill(1./sqrt(cdh1->emean()),calc_time2-tof1 );

		//up
		h2 = (TH2F*)gFile->Get(Form("CDHu%d_sqrtE_diff",cdh1->seg() )); h2->Fill(1./sqrt(cdh1->eu()),calc_time-tof1 );
		h2 = (TH2F*)gFile->Get(Form("CDHu%d_sqrtE_diff2",cdh1->seg() )); h2->Fill(1./sqrt(cdh1->eu()),calc_time2-tof1 );
		
		//down		
		h2 = (TH2F*)gFile->Get(Form("CDHd%d_sqrtE_diff",cdh1->seg() )); h2->Fill(1./sqrt(cdh1->ed()),calc_time-tof1 );						
		h2 = (TH2F*)gFile->Get(Form("CDHd%d_sqrtE_diff2",cdh1->seg() )); h2->Fill(1./sqrt(cdh1->ed()),calc_time2-tof1 ); 
		  }		
		//##Momentum#####
		int pnum;
		if(fabs(mom1)<0.1) pnum=1; 
		else if(fabs(mom1)<0.2) pnum=2;
		else if(fabs(mom1)<0.3) pnum=3;
		else                    pnum=4;

		h2 = (TH2F*)gFile->Get(Form("CDH_sqrtE_diff_p%d",pnum )); h2->Fill(1./sqrt(cdh1->emean()),calc_time-tof1 );
		h2 = (TH2F*)gFile->Get(Form("CDH_sqrtE_diff2_p%d",pnum )); h2->Fill(1./sqrt(cdh1->emean()),calc_time2-tof1 );
 
		h2 = (TH2F*)gFile->Get(Form("CDHu%d_sqrtE_diff_p%d",cdh1->seg(),pnum )); h2->Fill(1./sqrt(cdh1->eu()),calc_time-tof1 );
		h2 = (TH2F*)gFile->Get(Form("CDHd%d_sqrtE_diff_p%d",cdh1->seg(),pnum )); h2->Fill(1./sqrt(cdh1->ed()),calc_time-tof1 );		

		h2 = (TH2F*)gFile->Get(Form("CDHu%d_sqrtE_diff2_p%d",cdh1->seg(),pnum )); h2->Fill(1./sqrt(cdh1->eu()),calc_time2-tof1 );
		h2 = (TH2F*)gFile->Get(Form("CDHd%d_sqrtE_diff2_p%d",cdh1->seg(),pnum )); h2->Fill(1./sqrt(cdh1->ed()),calc_time2-tof1 );		
		
		//###T0####
		h2 = (TH2F*)gFile->Get(Form("T0u%d_sqrtE_diff",segT0 )); h2->Fill(1./sqrt(euT0),calc_time-tof1 );
		h2 = (TH2F*)gFile->Get(Form("T0d%d_sqrtE_diff",segT0 )); h2->Fill(1./sqrt(edT0),calc_time-tof1 );
		h2 = (TH2F*)gFile->Get(Form("T0%d_sqrtE_diff",segT0 )); h2->Fill(1./sqrt(emT0),calc_time-tof1 );
		
	      }
	    if(mass1=pMass	)//proton
	      {

		h2 = (TH2F*)gFile->Get(Form("CDH_calc_diff_proton" )); h2->Fill(calc_time,calc_time-tof1 );
		h2 = (TH2F*)gFile->Get(Form("CDH_E_diff_proton" )); h2->Fill(cdh1->emean(),calc_time-tof1 );
		
		h2 = (TH2F*)gFile->Get(Form("CDH%d_calc_diff_proton",cdh1->seg() )); h2->Fill(calc_time,calc_time-tof1 );
		h2 = (TH2F*)gFile->Get(Form("CDH%d_E_diff_proton",cdh1->seg() )); h2->Fill(cdh1->emean(),calc_time-tof1 );		
		h2 = (TH2F*)gFile->Get(Form("CDHu%d_E_diff_proton",cdh1->seg() )); h2->Fill(cdh1->eu(),calc_time_u-tof1_u );	
		h2 = (TH2F*)gFile->Get(Form("CDHd%d_E_diff_proton",cdh1->seg() )); h2->Fill(cdh1->ed(),calc_time_d-tof1_d );
	      }
	      
	  }
	
	
	//######2nd CDSTrack##########
	for( int it2=it1+1; it2<trackMan->nGoodTrack(); it2++ ){
	  CDSTrack *track2 = trackMan->Track(trackMan->GoodTrackID(it2));
	  //	if(!track2->CDHFlag()) continue;
	  HodoscopeLikeHit *cdh2 = track2->CDHHit();
#if SIM
	  double tof2 = cdh2->tmean()-ctmT0;
#else
	  double tof2 = cdh2->ctmean()-ctmT0;
#endif

	  double param2[5];
	  track2->GetParameters(param2);

	  //	  if(param2[2]==0 ||param2[4]==0   ) continue;

	  double drho2=param2[0], phi02=param2[1], rho2=track2->Rho(), dz2=param2[3], tlam2=param2[4];
	  double mom2 = track2->Momentum();
	  int ptype2 =-1;
	  if(mom2!=0) ptype2 = PID( tof2, mom2 );
	  //#######vertex#####	  
	  TrackVertex vertex;
	  TVector3 vtx1, vtx2, vtx;
	  if(trackMan->GetVertex(trackMan->GoodTrackID(it1),trackMan->GoodTrackID(it2),vertex) )
	    {	 
	      vtx=vertex.GetVertexPos_mean();
	      vtx1=vertex.GetVertexPos1();
	      vtx2=vertex.GetVertexPos2();

	      h1 = (TH1F*)gFile->Get("vx"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vy"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vxy"); h1->Fill( vtx.X(),vtx.Y() );
	       
	      if(vtx.Z()<-10)
		{
	      h1 = (TH1F*)gFile->Get("vx1"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vy1"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vxy1"); h1->Fill( vtx.X(),vtx.Y() );
		}

	      if(vtx.Z()>-10 && vtx.Z()<0)
		{
	      h1 = (TH1F*)gFile->Get("vx2"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vy2"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vxy2"); h1->Fill( vtx.X(),vtx.Y() );
		}
	      if(vtx.Z()>0 && vtx.Z()<10)
		{
	      h1 = (TH1F*)gFile->Get("vx3"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vy3"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vxy3"); h1->Fill( vtx.X(),vtx.Y() );
		}
	      h1 = (TH1F*)gFile->Get("vz"); h1->Fill( vtx.Z() );
	      TVector3 tmpv=vtx1-vtx2;
	      double dis=tmpv.Mag();
	      h1 = (TH1F*)gFile->Get("vrx_dis"); h1->Fill( dis );

	      ////#######vertex beam#####	  
	      //if( sqrt(vtx.X()*vtx.X()+vtx.Y()*vtx.Y() )>5.0 || fabs(vtx.z()) <3 || fabs(vtx.z()) >7) continue;
 	    for(int itr=0;itr<blTrackMan->ntrackBLC();itr++)
	      {
		LocalTrack *blc2=blTrackMan->trackBLC(itr);
		double x1,y1;
		blc2->XYPosatZ(vtx.z(),x1,y1);		  
		double chi_blc=blc2->chi2all();
		//		if(chi_blc>10) continue;

		h2 = (TH2F*)gFile->Get("BLCvtx_z_x"); h2->Fill( vtx.Z(),x1-vtx.X() );
		h2 = (TH2F*)gFile->Get("BLCvtx_z_y"); h2->Fill( vtx.Z(),y1-vtx.Y() );
		if(-44<vtx.z() && vtx.z()<-40)
		  {
		    h1 = (TH1F*)gFile->Get("BLCvtx_x1"); h1->Fill( x1-vtx.X() );
		    h1 = (TH1F*)gFile->Get("BLCvtx_y1"); h1->Fill( y1-vtx.Y() );
		    h2 = (TH2F*)gFile->Get("vx_vbx1"); h2->Fill( vtx.x(),x1 );
		    h2 = (TH2F*)gFile->Get("vy_vby1"); h2->Fill( vtx.y(),y1 );
		  }
		else if(-7<vtx.z() && vtx.z()<-3)
		  {
		    h1 = (TH1F*)gFile->Get("BLCvtx_x2"); h1->Fill( x1-vtx.X() );
		    h1 = (TH1F*)gFile->Get("BLCvtx_y2"); h1->Fill( y1-vtx.Y() );
		    h2 = (TH2F*)gFile->Get("vx_vbx2"); h2->Fill( vtx.x(),x1 );
		    h2 = (TH2F*)gFile->Get("vy_vby2"); h2->Fill( vtx.y(),y1 );
		    h2 = (TH2F*)gFile->Get("BLCvtx_r_x"); h2->Fill( sqrt(vtx.X()*vtx.X()+vtx.Y()*vtx.Y()),x1-vtx.X() );
		    h2 = (TH2F*)gFile->Get("BLCvtx_r_y"); h2->Fill( sqrt(vtx.X()*vtx.X()+vtx.Y()*vtx.Y()),y1-vtx.Y() );
		    MathTools *tool=new MathTools();
		    double theta1=tool->CalcDeg(vtx.x(),vtx.y() )*TMath::DegToRad();
		    double theta2=tool->CalcDeg(x1,y1 )*TMath::DegToRad();

		    h2 = (TH2F*)gFile->Get("BLCvtx_r_theta"); h2->Fill( sqrt(vtx.X()*vtx.X()+vtx.Y()*vtx.Y()),theta2-theta1 );
		    delete tool;
		  }
		else if(3<vtx.z() && vtx.z()<7)
		  {
		    h1 = (TH1F*)gFile->Get("BLCvtx_x3"); h1->Fill( x1-vtx.X() );
		    h1 = (TH1F*)gFile->Get("BLCvtx_y3"); h1->Fill( y1-vtx.Y() );
		    h2 = (TH2F*)gFile->Get("vx_vbx3"); h2->Fill( vtx.x(),x1 );
		    h2 = (TH2F*)gFile->Get("vy_vby3"); h2->Fill( vtx.y(),y1 );
		  }

	      }



	      TrackVertex vertexb;
	      TVector3 vtxb1, vtxb2, vtxb;
	      //	      if(blTrackMan->ntrackBLC()!=1) continue;
	      if(trackMan->GetVertex_beam(trackMan->GoodTrackID(it1),0,vertexb) )
		{	 
		  vtxb=vertexb.GetVertexPos_mean();
		  vtxb1=vertexb.GetVertexPos1();
		  vtxb2=vertexb.GetVertexPos2();
		  h1 = (TH1F*)gFile->Get("vbx2"); h1->Fill( vtxb.X() );
		  h1 = (TH1F*)gFile->Get("vby2"); h1->Fill( vtxb.Y() );
		  h1 = (TH1F*)gFile->Get("vbxy2"); h1->Fill( vtxb.X(),vtxb.Y() );
		  h1 = (TH1F*)gFile->Get("vbz2"); h1->Fill( vtxb.Z() );
		  TVector3 tmpbv=vtxb1-vtxb2;
		  double dis=tmpbv.Mag();
		  TVector3 tmpv_vb=vtxb-vtx;
		  h1 = (TH1F*)gFile->Get("vrxb_dis2"); h1->Fill( dis );
		  h1 = (TH1F*)gFile->Get("v_vb_dis"); h1->Fill( tmpv_vb.Mag() );
		}
	      //	      if(blTrackMan->ntrackBLC()!=1) continue;
	      if(trackMan->GetVertex_beam(trackMan->GoodTrackID(it2),0,vertexb) )
		{	 
		  vtxb=vertexb.GetVertexPos_mean();
		  vtxb1=vertexb.GetVertexPos1();
		  vtxb2=vertexb.GetVertexPos2();
		  h1 = (TH1F*)gFile->Get("vbx2"); h1->Fill( vtxb.X() );
		  h1 = (TH1F*)gFile->Get("vby2"); h1->Fill( vtxb.Y() );
		  h1 = (TH1F*)gFile->Get("vbxy2"); h1->Fill( vtxb.X(),vtxb.Y() );
		  h1 = (TH1F*)gFile->Get("vbz2"); h1->Fill( vtxb.Z() );
		  TVector3 tmpbv=vtxb1-vtxb2;
		  double dis=tmpbv.Mag();
		  TVector3 tmpv_vb=vtxb-vtx;
		  h1 = (TH1F*)gFile->Get("vrxb_dis2"); h1->Fill( dis );
		  h1 = (TH1F*)gFile->Get("v_vb_dis"); h1->Fill( tmpv_vb.Mag() );
		}



	    }
	  



	  if( (ptype1!=3 || ptype1!=0)  && (ptype2!=0 ||ptype2!=3) ){// no use PID
	    CDSTrack *tracktmp; 
	    HodoscopeLikeHit *cdhtmp;
	    double paramtmp[5];
	    
	    double mass1, mass2;
	    mass1 = piMass; mass2 = piMass; 

	    MathTools *mtool = new MathTools();
	    TVector3 Pp1 = mtool->CalcHelixMom(  param1, vtx.Z() );
	    TVector3 Pp2 = mtool->CalcHelixMom(  param2, vtx.Z() );
	    delete mtool;
	    TVector3 Pp = Pp1 + Pp2;
	    TLorentzVector L_p1; L_p1.SetVectM( Pp1, mass1 );
	    TLorentzVector L_p2; L_p2.SetVectM( Pp2, mass2);

	    double cosOA = (Pp1.Dot(Pp2))/Pp1.Mag()/Pp2.Mag();
	    double im = (L_p1+L_p2).M();
	      h1 = (TH1F*)gFile->Get("IMpipi_nopid"); h1->Fill( im );

	  }
	  

	  if( ((ptype1==3&&ptype2==2) || (ptype1==2&&ptype2==3)) ||// select p pi-
		((ptype1==1&&ptype2==2) || (ptype1==2&&ptype2==1)) ){// select pi+ pi-
	    CDSTrack *tracktmp; 
	    HodoscopeLikeHit *cdhtmp;
	    double paramtmp[5];
	    
	    if( (ptype1==3&&ptype2==2) || (ptype1==1&&ptype2==2) ){}
	    else{
	      tracktmp = track1; track1 = track2; track2 = tracktmp;
	      cdhtmp = cdh1; cdh1 = cdh2; cdh2 = cdhtmp;
	      for( int i=0; i<5; i++ ){ paramtmp[i] = param1[i]; param1[i] = param2[i]; param2[i] = paramtmp[i]; }
	    }
	    
	    double mass1, mass2;
	    if( ptype1==3 ){ mass1 = pMass;  mass2 = piMass; }
	    else           { mass1 = piMass; mass2 = piMass; }
	    MathTools *mtool = new MathTools();
	    TVector3 Pp1 = mtool->CalcHelixMom(  param1, vtx1.Z() );
	    TVector3 Pp2 = mtool->CalcHelixMom(  param2, vtx2.Z() );
	    delete mtool;
	    TVector3 Pp = Pp1 + Pp2;
	    TLorentzVector L_p1; L_p1.SetVectM( Pp1, mass1 );
	    TLorentzVector L_p2; L_p2.SetVectM( Pp2, mass2);

	    TLorentzVector L_p1_p; L_p1_p.SetVectM( Pp1*1.01, mass1 );
	    TLorentzVector L_p2_p; L_p2_p.SetVectM( Pp2*1.01, mass2);

	    TLorentzVector L_p1_m; L_p1_m.SetVectM( Pp1*0.99, mass1 );
	    TLorentzVector L_p2_m; L_p2_m.SetVectM( Pp2*0.99, mass2);

	    double cosOA = (Pp1.Dot(Pp2))/Pp1.Mag()/Pp2.Mag();
	    double im = (L_p1+L_p2).M();
	    double im_p = (L_p1_p+L_p2_p).M();
	    double im_m = (L_p1_m+L_p2_m).M();
	    if( ptype1==3 ){ // ppi-
	      h1 = (TH1F*)gFile->Get("vxppi"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vyppi"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vzppi"); h1->Fill( vtx.Z() );
	      h3 = (TH3F*)gFile->Get("vppi");  h3->Fill( vtx.X(), vtx.Y(), vtx.Z() );
	      h1 = (TH1F*)gFile->Get("cosOAppi"); h1->Fill( cosOA );
	      h1 = (TH1F*)gFile->Get("IMppi"); h1->Fill( im );
	      h1 = (TH1F*)gFile->Get("IMppi_m"); h1->Fill( im_m );
	      h1 = (TH1F*)gFile->Get("IMppi_p"); h1->Fill( im_p );
	      h2 = (TH2F*)gFile->Get("IMcosOAppi"); h2->Fill( im, cosOA );
 	      h1 = (TH1F*)gFile->Get("costppi"); h1->Fill( Pp.CosTheta() );
 	      h1 = (TH1F*)gFile->Get("phippi"); h1->Fill( Pp.Phi()*TMath::RadToDeg() );
 	      h2 = (TH2F*)gFile->Get("IMcostppi"); h2->Fill( im, Pp.CosTheta() );
	      //   std::cout<<"Phi "<<Pp.Phi()*TMath::RadToDeg()<<std::endl;
// 	      if( 2<tofBHDT0 ){
// 		h1 = (TH1F*)gFile->Get("IMppiKselected"); h1->Fill( im );
// 	      }

//###########Study Lambda############
	      if( (1.14-0.12)<im && im<(1.14+0.12 ))
		{
		  
		  //		  h3 = (TH3F*)gFile->Get("vLambda"); h3->Fill( vtx.X(),vtx.Y(),vtx.Z() );
		  TVector3 tmpv=vtx1-vtx2;
		  double dis=tmpv.Mag();
		  h1 = (TH1F*)gFile->Get("v_disLambda"); h1->Fill( dis );
		  
		  h1 = (TH1F*)gFile->Get("pLambda"); h1->Fill( Pp.Mag() );
		  double target_pos;
		  if( (Pp.z()>0 &&vtx.z()>5.0) || (Pp.z()<0 && vtx.z()<5.0 && vtx.z()>-5.0) )  target_pos=5.0;
		  else 	if( (Pp.z()<0 &&vtx.z()<-5.0) || (Pp.z()>0 && vtx.z()<5.0 && vtx.z()>-5.0) )  target_pos=-5.0;
		  TVector3 tmp=(1/fabs(Pp.z()) )*Pp;
		  tmp=(vtx.z()-target_pos)*tmp;
		  double tof_dis=tmp.Mag();
		  double beta=Pp.Mag()/sqrt(Pp.Mag()*Pp.Mag()+lMass*lMass);
		  double tof=tof_dis/beta/Const;
		  h1 = (TH1F*)gFile->Get("timeLambda"); h1->Fill( tof );
		  h1 = (TH1F*)gFile->Get("disLambda"); h1->Fill( tof_dis );
		  TVector3 Pk;Pk.SetXYZ(0,0,0.9);//GeV/c
		  TVector3 Ppr;Ppr.SetXYZ(0,0,0);//GeV/c
		  TLorentzVector L_K; L_K.SetVectM( Pk, kpMass );
		  TLorentzVector L_proton; L_proton.SetVectM( Ppr, pMass);
		  TLorentzVector L_L; L_L.SetVectM( Pp, lMass);
		  TLorentzVector L_x; L_x.SetVect( Pk-Pp );L_x.SetE( L_K.E()+L_proton.E()-L_L.E() );
		  h1 = (TH1F*)gFile->Get("missingMassLambda"); h1->Fill( L_x.M() );
		}

	    }
	    else if( ptype1==1){ // pi+pi-
	      
	      h1 = (TH1F*)gFile->Get("vxpipi"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vypipi"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vzpipi"); h1->Fill( vtx.Z() );
	      h3 = (TH3F*)gFile->Get("vpipi");  h3->Fill( vtx.X(), vtx.Y(), vtx.Z() );
	      h1 = (TH1F*)gFile->Get("cosOApipi"); h1->Fill( cosOA );
	      h1 = (TH1F*)gFile->Get("IMpipi"); h1->Fill( im );
	      h1 = (TH1F*)gFile->Get("IMpipi_m"); h1->Fill( im_m );
	      h1 = (TH1F*)gFile->Get("IMpipi_p"); h1->Fill( im_p );

	      h2 = (TH2F*)gFile->Get("IMcosOApipi"); h2->Fill( im, cosOA );
	      h1 = (TH1F*)gFile->Get("costpipi"); h1->Fill( Pp.CosTheta() );
	      h1 = (TH1F*)gFile->Get("phipipi"); h1->Fill( Pp.Phi()*TMath::RadToDeg() );
 	      h2 = (TH2F*)gFile->Get("IMcostpipi"); h2->Fill( im, Pp.CosTheta() );
	      if(  (cosOA < 0.6 && cosOA > -0.95)) h1 = (TH1F*)gFile->Get("IMpipi_cut"); h1->Fill( im );
	    }
	    
	  }
	}
      }
      
  }
  
  
  gFile->Write();
  gFile->Close();
  
  gSystem->Exit(1);
  gROOT->GetApplication()->Run();
  //  theApp.Run();

  return 0;
}

