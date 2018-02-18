

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

//const double MassThPPi = 0.7; // GeV/c2
//const double TOFOffset = -3.5; // nsec Oct2010
//const double MassThPPi = 0.7; // GeV/c2
//const double TOFOffset = 7.0; // nsec Oct2010
//const double MassThPPi = 1.1; // GeV/c2
//const double TOFOffset = -2.5; // nsec Oct2010 //2012/3/22
//const double MassThPPi = 0.7; // GeV/c2
//const double TOFOffset = 6.2; // nsec 2012/03/07
//const double TOFOffset = 5.5; // nsec 2012/03/09

const double MassThPPi = 0.83; // GeV/c2 /2012/6/22
const double TOFOffset = -4.2; // nsec 2012/06/22

const double piMass = 0.13957;
const double pMass = 0.938272;
const double nMass = 0.939565;
const double lMass = 1.115683;
const double dMass = 1.87561;
const double kpMass = 0.4936;
const double k0sMass = 0.497614;
const double ThreeHeMass = 2.80839;


const double ELV=14.0;// [cm/ns]

const double Const=29.97;// [cm/ns]
#define SIM 0
#define CDHHIST 1


int PID( double tof, double mom ) // 0:unknown, 1:pi+, 2;pi-, 3:proton
{


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

  return val;
}




void InitializeHist()
{


  new TH1F( "tofBHDT0", "tof BHD and T0", 200, 0, 40 );
  for(int seg=0;seg<5;seg++)  new TH1F( Form("timeT0_%d",seg), Form("time of T0 %d",seg), 100, -5, 5 );
  new TTree("Tvtx","vertex tree");

  
  new TH1F( "nTrack", "nTrack", 200, 0, 100 );
  new TH1F( "chi", "chi", 300, 0, 30 );

  new TH1F( "mom", "mom", 500, -1, 1 );
  new TH2F( "pid", "pid", 300, -15, 15, 300, -1.5, 1.5 );
  new TH2F( "pid_cut", "pid_cut", 300, -10, 20, 300, -1.5, 1.5 );
  new TH1F( "phi", "phi of track", 200, 0, 360 );
  new TH1F( "momp", "mom proton", 500, -1, 1 );
  new TH1F( "mompip", "mom piplus", 500, -1, 1 );
  new TH1F( "mompim", "mom piminus", 500, -1, 1 );

  for(int layer=1;layer<=15;layer++)
    {
      new TH1F(Form("HitPatCDC%d",layer),Form("HitPatCDC%d",layer),150,0,150);
    }


  new TH1F( "chi_BLC", "chi_BLC", 3000, 0, 30 );
  new TH1F( "chixz_BLC", "chixz_BLC", 300, 0, 30 );
  new TH1F( "chiyz_BLC", "chiyz_BLC", 300, 0, 30 );
  new TH1F( "Beamx", "Beamx", 6000, -30, 30 );
  new TH1F( "Beamy", "Beamy", 6000, -30, 30 );
  new TH1F( "Beamdx", "Beamdx", 10000, -0.05, 0.05 );
  new TH1F( "Beamdy", "Beamdy", 10000, -0.05, 0.05 );




  new TH1F( "vvrx_dis", "vvrx_dis", 2000, 0, 2 );
  new TH1F( "vvrx_cdc", "vvrx_cdc", 2000, 0, 2 );
  new TH1F( "vvrx_cdcx", "vvrx_cdcx", 2000, -1, 1 );
  new TH1F( "vvrx_cdcy", "vvrx_cdcy", 2000, -1, 1 );
  new TH1F( "vvrx_cdcz", "vvrx_cdcz", 2000, -1, 1 );
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
  new TH2F( "vbzy", "vbzy", 600, -30, 30 , 500, -25, 25 );
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

  new TH2F( Form("CDHu_ADC_mom"), Form("CDHu_ADC_mom"), 4000,0,4000 , 200, 0.0 , 1);


//   for(int seg=1;seg<=36;seg++)
//     { 
//       new TH2F( Form("CDHu%d_ADC_mom",seg), Form("CDHu%d_ADC_mom",seg), 4000,0,4000 , 200, 0.0 , 1);
//       new TH2F( Form("CDHd%d_ADC_mom",seg), Form("CDHd%d_ADC_mom",seg), 4000,0,4000 , 200, 0.0 , 1);
//       new TH2F( Form("CDH%d_ADCmean_mom_pi",seg), Form("CDH%d_ADCmean_mom_pi",seg), 1000,-5,45 , 200, 0.0 , 1);
//       new TH2F( Form("CDH%d_ADCmean_mom_p",seg), Form("CDH%d_ADCmean_mom_p",seg), 1000,-5,45, 200, 0.0, 1.0);
//     }

  //CDH Slewing Hist

  new TH2F( "mass_mom", "mass_mom", 750, 0.0 , 2.5, 240,-1.2,1.2);
  new TH2F( "mass2_mom", "mass2_mom", 1020, -0.1 , 5.0, 240,-1.2,1.2);
  new TH2F( "mom_overbeta", "mom_overbeta", 1000, 0 , 10.0, 240,-1.2,1.2);

  new TH2F( "mass_mom2", "mass_mom2", 750, 0.0 , 2.5, 240,-1.2,1.2);
  new TH2F( "mass2_mom2", "mass2_mom2", 1020, -0.1 , 5.0, 240,-1.2,1.2);
  new TH2F( "mom_overbeta2", "mom_overbeta2", 1000, 0 , 10.0, 240,-1.2,1.2);

  new TH2F( "mass_mom3", "mass_mom3", 750, 0.0 , 2.5, 240,-1.2,1.2);
  new TH2F( "mass2_mom3", "mass2_mom3", 1020, -0.1 , 5.0, 240,-1.2,1.2);
  new TH2F( "mom_overbeta3", "mom_overbeta3", 1000, 0 , 10.0, 240,-1.2,1.2);

  new TH2F( "massp_mom3", "massp_mom3", 750, 0.0 , 2.5, 240,-1.2,1.2);
  new TH2F( "mass2p_mom3", "mass2p_mom3", 1020, -0.1 , 5.0, 240,-1.2,1.2);

  new TH2F( "massn_mom3", "massn_mom3", 750, 0.0 , 2.5, 240,-1.2,1.2);
  new TH2F( "mass2n_mom3", "mass2n_mom3", 1020, -0.1 , 5.0, 240,-1.2,1.2);


  new TH2F( "mass_mom_tmp", "mass_mom_tmp", 750, 0.0 , 2.5, 240,-1.2,1.2);
  new TH2F( "mass2_mom_tmp", "mass2_mom_tmp", 1020, -0.1 , 5.0, 240,-1.2,1.2);
  new TH2F( "mom_overbeta_tmp", "mom_overbeta_tmp", 1000, 0 , 10.0, 240,-1.2,1.2);

  //deuteron
  new TH1F( "ntrack_d", "ntrack_d", 10, 0, 10 );
  new TH1F( "mom_d", "mom_d", 2000, 0, 2 );
  new TH2F( "pid_pat_d", "pid_pat_d", 10, 0, 10,10,0,10 );
  new TH1F( "cos_d", "cos_d", 200, -1, 1 );
  new TH2F( "mom_cos_d", "cos_d", 2000, 0, 2,200, -1, 1 );
  new TH1F( "Mmd", "Missing mass 3He(K,d)X", 1500, 0.0, 3.0 );

  //kaon
  new TH1F( "ntrack_k", "ntrack_k", 10, 0, 10 );
  new TH1F( "mom_k", "mom_k", 2000, 0, 2 );
  new TH2F( "pid_pat_k", "pid_pat_k", 10, 0, 10,10,0,10 );
  new TH1F( "cos_k", "cos_k", 200, -1, 1 );
  new TH2F( "mom_cos_k", "cos_k", 2000, 0, 2,200, -1, 1 );

  new TH1F( "MmK", "Missing mass 3He(K,K)X", 1500, 0.0, 3.0 );
  new TH1F( "MmK_1", "Missing mass 3He(K,K)X track1", 1500, 0.0, 3.0 );

  new TH1F( "nK_cos", "cos of n @Missing mass N(K,K)N", 200, -1,1 );
  new TH2F( "nK_dxy", "nK_dxy", 600, -0.3, 0.3,600, -0.3, 0.3 );
  new TH1F( "nK_mom", "mom of n @Missing mass N(K,K)N", 2500, -1,1.5 );

  //Lambda
  new TH1F( "ntrack_L", "ntrack_L", 10, 0, 10 );
  //  new TH1F( "mom_L", "mom_L", 2000, 0, 2 );
  new TH2F( "pid_pat_L", "pid_pat_L", 10, 0, 10,10,0,10 );
  //  new TH1F( "cos_L", "cos_L", 200, -1, 1 );
  //  new TH2F( "mom_cos_L", "mom_cos_L", 2000, 0, 2,200, -1, 1 );
  new TH1F( "pLambda", "pLambda", 1000, 0, 1.0 );
  new TH1F( "cosLambda", "cosLambda", 200, -1, 1 );
  new TH2F( "mom_cos_L", "mom_cos_L", 2000, 0, 2,200, -1, 1 );
  new TH1F( "nK0s_cos", "cos of n @Missing mass 3He(K,K0s)n", 200, -1,1 );
  new TH2F( "nK0s_dxy", "nK0s_dxy", 600, -0.3, 0.3,600, -0.3, 0.3 );
  new TH1F( "nK0s_mom", "mom of n @Missing mass (K,K0s)X", 2500, -1,1.5 );
  //K0s
  new TH1F( "ntrack_k0s", "ntrack_k0s", 10, 0, 10 );
  //  new TH1F( "mom_L", "mom_L", 2000, 0, 2 );
  new TH2F( "pid_pat_k0s", "pid_pat_k0s", 10, 0, 10,10,0,10 );
  //  new TH1F( "cos_L", "cos_L", 200, -1, 1 );
  //  new TH2F( "mom_cos_L", "mom_cos_L", 2000, 0, 2,200, -1, 1 );
  new TH1F( "pK0s", "pK0s", 1000, 0, 1.0 );
  new TH1F( "cosK0s", "cosK0s", 200, -1, 1 );
  new TH2F( "mom_cos_k0s", "mom_cos_k0s", 2000, 0, 2,200, -1, 1 );


  // IM ppi
  new TH1F( "vxppi", "vxppi", 500, -25, 25 );
  new TH1F( "vyppi", "vyppi", 500, -25, 25 );
  new TH1F( "vzppi", "vzppi", 1000, -50, 50 );
  new TH3F( "vppi",  "vppi", 200, -25, 25, 200, -25, 25, 200, -25, 25 );
  new TH1F( "IMppi", "IM proton piminus", 2000, 1.0, 2 );
  new TH1F( "IMppi_bg", "IM proton piminus background", 2000, 1.0, 2 );
  new TH1F( "IMppi_bg2", "IM proton piminus background2", 2000, 1.0, 2 );
  new TH1F( "IMppi_cut", "IM proton piminus with cut", 2000, 1.0, 2 );
  new TH1F( "IMppi_target", "IM proton piminus from target", 2000, 1.0, 2 );
  new TH2F( "IMppi_dis", "IM proton piminus vs DC dis", 500, 1.0, 2,150,0,15 );
  new TH1F( "IMppiKselected", "IM proton piminus", 500, 1.0, 2 );
  new TH1F( "cosOAppi", "cosOA between proton and #pi^-", 200, -1, 1 );
  new TH2F( "IMcosOAppi", "IM and cosOA for ppi", 200, 1, 2, 200, -1, 1);
  new TH1F( "costppi", "cost of ppi system", 200, -1, 1 );
  new TH1F( "phippi", "phi of ppi system", 200, -180, 180 );
  new TH2F( "IMcostppi", "IM and cost of ppi system", 500, 1, 2, 200, -1, 1 );
  new TH1F( "IMppi_vdis", "IMppi_vdis", 500, 1, 2 );

//   new TH1F( "IMppi_test", "IMppi_test", 500, 1.0, 2 );
//   new TH1F( "IMppi_test_xy", "IMppi_test_xy", 500, 1.0, 2 );
//   new TH1F( "IMppi_test_z", "IMppi_test_z", 500, 1.0, 2 ); 
//   new TH1F( "IMppi_test_beam", "IMppi_test_beam", 500, 1.0, 2 );

  new TH2F( "IMppi_xest_yest", "IMppi_xest_yest", 500, -25, 25,500, -25, 25  );
  new TH2F( "IMppi_yest_zest", "IMppi_yest_zest", 500, -25, 25,500, -25, 25  );

  //   new TH3F( "vLambda", "vLambda", 1000, -25, 25,1000, -25, 25,1000, -25, 25 );
   new TH1F( "v_disLambda", "v_disLambda", 500, 0, 5 );

   // 3 part 
   new TH1F( "pL_proton_m", "pL_proton_m", 1000, 0, 1.0 );
   new TH1F( "pL_pi_m", "pL_pi_m", 1000, 0, 1.0 );
   new TH1F( "pLambda_m", "pLambda_m", 1000, 0, 1.0 );
   new TH1F( "LifetimeLambda_m", "LifetimeLambda_m", 1000, 0, 1.0 );
   new TH1F( "MmL", "Missing mass 3He(K,L)X", 1000, 0.0, 2.0 );
   new TH1F( "MmL_cut", "Missing mass 3He(K,L)X", 1000, 0.0, 2.0 );

   new TH1F( "pL_proton_p", "pL_proton_p", 1000, 0, 1.0 );
   new TH1F( "pL_pi_p", "pL_pi_p", 1000, 0, 1.0 );
   new TH1F( "pLambda_p", "pLambda_p", 1000, 0, 1.0 );
   new TH1F( "LifetimeLambda_p", "LifetimeLambda_p", 1000, 0, 1.0 );

   new TH1F( "pL_proton", "pL_proton", 1000, 0, 1.0 );
   new TH1F( "pL_pi", "pL_pi", 1000, 0, 1.0 );
  
   new TH1F( "LifetimeLambda", "LifetimeLambda", 1000, 0, 1.0 );


   new TH1F( "disLambda", "disLambda", 2000, 0, 10 );
   new TH1F( "missingMassLambda", "missingMassLambda", 1500, -0.5, 1.0 );

  new TH1F( "vxpipi", "vxpipi", 500, -25, 25 );
  new TH1F( "vypipi", "vypipi", 500, -25, 25 );
  new TH1F( "vzpipi", "vzpipi", 1000, -50, 50 );


  //IM pi pi 
  new TH3F( "vpipi",  "vpipi", 200, -25, 25, 200, -25, 25, 200, -25, 25 );
  new TH1F( "IMpipi", "IM piplus piminus", 500, 0, 1.0 );
  new TH1F( "IMpipi_cut", "IM piplus piminus wih cut", 500, 0, 1.0 );
  new TH1F( "IMpipi_target", "IM piplus piminus from traget", 500, 0, 1.0 );
  new TH2F( "IMpipi_dis", "IM piplus piminus vs DC dis", 500, 0, 1.0,150,0,15 );
  new TH1F( "IMpipi_m", "IM piplus piminus p*0.99", 500, 0, 1.0 );
  new TH1F( "IMpipi_p", "IM piplus piminus p*1.01", 500, 0, 1.0 );
  //  new TH1F( "IMpipi_cut", "IM piplus piminus cut OA", 500, 0, 1.0 );
  new TH1F( "IMpipi_nopid", "IM piplus piminus nopid", 500, 0, 1.0 );
  new TH1F( "cosOApipi", "cosOA between #pi^+ and #pi^-", 200, -1, 1 );
  new TH2F( "IMcosOApipi", "IM and cosOA for pipi", 500, 0, 1, 200, -1, 1);
  new TH1F( "costpipi", "cost of pipi system", 200, -1, 1);
  new TH1F( "phipipi", "phi of pipi system", 200, -180, 180 );
  new TH2F( "IMcostpipi", "IM and cost of pipi system", 500, 0, 1, 200, -1, 1 );
   new TH1F( "IMpipi_vdis", "IMpipi_vdis", 500, 0, 5 );
   new TH1F( "v_disK0s", "v_disK0s", 500, 0, 5 );

   new TH1F( "pK0s_pip_m", "pK0s_pip_m", 1000, 0, 1.0 );
   new TH1F( "pK0s_pim_m", "pK0s_pim_m", 1000, 0, 1.0 );
   new TH1F( "pK0s_m", "pK0s_m", 1000, 0, 1.0 );
   new TH1F( "LifetimeK0s_m", "LifetimeK0s_m", 1000, 0, 1.0 );
   new TH1F( "MmK0s_m", "Missing mass 3He(K,K0s)X", 1000, 0.0, 2.0 );

   new TH1F( "pK0s_pip_p", "pK0s_pip_p", 1000, 0, 1.0 );
   new TH1F( "pK0s_pim_p", "pK0s_pim_p", 1000, 0, 1.0 );
   new TH1F( "pK0s_p", "pK0s_p", 1000, 0, 1.0 );
   new TH1F( "LifetimeK0s_p", "LifetimeK0s_p", 1000, 0, 1.0 );
   new TH1F( "MmK0s_p", "Missing mass 3He(K,K0s)X", 1000, 0.0, 2.0 );

   new TH1F( "pK0s_pip", "pK0s_pip", 1000, 0, 1.0 );
   new TH1F( "pK0s_pim", "pK0s_pim", 1000, 0, 1.0 );
   new TH1F( "disK0s", "disK0s", 2000, 0, 10 );
   new TH1F( "LifetimeK0s", "LifetimeK0s", 1000, 0, 1.0 );
   new TH1F( "MmK0s", "Missing mass 3He(K,K0s)X", 1000, 0.0, 2.0 );
   new TH1F( "MmK0s_cut", "Missing mass 3He(K,K0s)X", 1000, 0.0, 2.0 );
   //   new TH1F( "missingMassLambda", "missingMassLambda", 1500, -0.5, 1.0 );

   //IM Lambda p
  new TH1F( "IMkpp", "IM Lambda proton", 500, 2.0, 3 );
  new TH1F( "MmLp", "Missing mass Lambda proton", 1000, 0.0, 2.0 );
  new TH1F( "vxkpp", "vxkpp", 500, -25, 25 );
  new TH1F( "vykpp", "vykpp", 500, -25, 25 );
  new TH1F( "vzkpp", "vzkpp", 1000, -50, 50 );
  new TH1F( "pMmLp", "Missing mass Lambda proton", 1000, 0.0, 2.0 );

  new TH1F( "IMkpp_m", "IM Lambda proton m", 500, 2.0, 3 );
  new TH1F( "MmLp_m", "Missing mass Lambda proton", 1000, 0.0, 2.0 );
  new TH1F( "pMmLp_m", "Missing mass Lambda proton", 1000, 0.0, 2.0 );
  new TH1F( "IMkpp_btb_m", "IM Lambda proton cosOA <-0.8", 500, 2.0, 3 );
  new TH2F( "IMcosOAkpp_m", "IM and cosOA for kpp", 500, 2.0, 3, 200, -1, 1);

  new TH1F( "IMkpp_p", "IM Lambda proton p", 500, 2.0, 3 );
  new TH1F( "MmLp_p", "Missing mass Lambda proton", 1000, 0.0, 2.0 );
  new TH1F( "pMmLp_p", "Missing mass Lambda proton", 1000, 0.0, 2.0 );
  new TH1F( "IMkpp_btb_p", "IM Lambda proton cosOA <-0.8", 500, 2.0, 3 );
  new TH2F( "IMcosOAkpp_p", "IM and cosOA for kpp", 500, 2.0, 3, 200, -1, 1);

  new TH3F( "vkpp",  "vkpp", 200, -25, 25, 200, -25, 25, 200, -25, 25 );
  new TH1F( "cosOAkpp", "cosOA between #Lambda+ and p", 200, -1, 1 );
  new TH1F( "CMcosOAkpp", "CMcosOA between #Lambda+ and p", 200, -1, 1 );
  new TH1F( "IMkpp_mn", "IM Lambda proton mm neutron cut", 500, 2.0, 3 );
  new TH1F( "IMkpp_btb", "IM Lambda proton CMcosOA <-0.9", 500, 2.0, 3 );
  new TH2F( "IMcosOAkpp", "IM and cosOA for kpp", 500, 2.0, 3, 200, -1, 1);
  new TH1F( "costkpp", "cost of kpp system", 200, -1, 1);
  new TH1F( "phikpp", "phi of kpp system", 200, -180, 180 );
  new TH2F( "IMcostkpp", "IM and cost of kpp system", 500, 2.0, 3, 200, -1, 1 );
  //NC
  new TH1F( "tofNC", "tofNC", 1020, -10, 500 );
  new TH1F( "tofNC_vtx", "tofNC_vtx", 1020, -10, 500 );
  new TH1F( "tofNC_L", "tofNC_L", 1020, -10, 500 );
  new TH1F( "tofNC_Lp", "tofNC_Lp", 1020, -10, 500 );
  new TH1F( "obetaNC_vtx", "obetaNC_vtx", 1100, -1, 10 );
  new TH1F( "obetaNC_L", "obetaNC_L", 1100, -1, 10 );
  new TH1F( "obetaNC_Lp", "obetaNC_Lp", 1100, -1, 10 );

  new TH1F( "NChit_nK", "NChit_nK", 2000, -100, 100 );
  new TH1F( "NChit_nK0s", "NChit_nK0s", 2000, -100, 100 );
  new TH2F( "NChitX_nK0s", "NChitX_nK0s", 2000, -200, 200,2000, -200, 200 );
  new TH2F( "NChitY_nK0s", "NChitY_nK0s", 2000, -200, 200,2000, -200, 200 );
  new TH1F( "NCeff_nK0s", "NCeff_nK0s", 16, -160, 160 );
  new TH1F( "NCeff_nK0s2", "NCeff_nK0s2", 16, -160, 160 );
  new TH2F( "NCxy_nK", "NCxy_nK", 3000, -300, 300,3000,-300,300 );
  new TH2F( "NCxy_nK0s", "NCxy_nK0s", 3000, -300, 300,3000,-300,300 );

  new TH1F( "momNC_vtx", "momNC_vtx", 3100, 0.1, 3 );
  new TH1F( "momNC_L", "momNC_L", 3100, 0.1, 3 );
  new TH1F( "momNC_Lp", "momNC_Lp", 3100, 0.1, 3 );
  new TH1F( "Q_Lpn", "Q_Lpn", 300, -3, 3 );
  new TH1F( "Mm_Lpn", "Mm_Lpn", 4000, -1, 3 );
  new TH2F( "dalitz_Lpn", "dalitz_Lpn", 1200, -0.6,0.6, 1000,0,1 );
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



  ConfMan *conf = new ConfMan(confFile);
  conf->Initialize();


  TFile *f;
  TTree *evtree;
  f = new TFile( inFile.c_str());
  evtree = (TTree*)f->Get( "EventTree" );
  
  TFile *fout;
  fout = new TFile( outFile.c_str(), "recreate" );

  InitializeHist();


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


  Float_t beamvx;
  Float_t beamvy;
  Float_t beamvz;
  Float_t beamvdis;

  Float_t cdcvx;
  Float_t cdcvy;
  Float_t cdcvz;
  Float_t cdcvdis;
  Float_t tpid_beam;
  TTree *tmptree;
  tmptree = (TTree*)gFile->Get("Tvtx");      
  tmptree->Branch("beamvx",&beamvx);
  tmptree->Branch("beamvy",&beamvy);
  tmptree->Branch("beamvz",&beamvz);
  tmptree->Branch("beamvdis",&beamvdis);

  tmptree->Branch("cdcvx",&cdcvx);
  tmptree->Branch("cdcvy",&cdcvy);
  tmptree->Branch("cdcvz",&cdcvz);
  tmptree->Branch("cdcvdis",&cdcvdis);
  tmptree->Branch("pid_beam",&tpid_beam);

  
  std::vector <TLorentzVector> vLv_p;
  std::vector <TLorentzVector> vLv_pim;


  int nev = evtree->GetEntries();
  std::cout <<" entry:" << nev << std::endl;

  //###test for trackout###
  //  ofstream of;
  //  std::string OutFileName=Form("trackdata.txt");
  //  of.open("trackdata.txt");

  std::cout << " AllEvent : " << nev << std::endl;
  std::cout << "|0%                  |50%                |100% "<<inFile     
	    <<   std::endl;
  int moniter=0;  
  for( int iev=0; iev<nev; iev++ ){

    //    if(iev>100000) continue;

    if( iev==0 )
      std::cout << "|";
    if( iev%(nev/40)==0 )
      std::cout << "*"<<std::flush;
    if( iev%(100)==0 )
      {
	if(moniter==0) std::cout << "\\\b"<<std::flush;
	else if(moniter==1) std::cout << "-\b"<<std::flush;
	else if(moniter==2) std::cout << "/\b"<<std::flush;
	else if(moniter==3) std::cout << "|\b"<<std::flush;
	moniter++;
	if(moniter==4) moniter=0;
      }
    if( iev==nev-1 )
	  std::cout << "| fin!" << std::endl;
    //    if(iev<15000) continue;
    //    if(iev>30000 && iev<90000) continue;
    //    if(iev>10000) continue;
    evtree->GetEvent(iev);



    beamvx=-999;beamvy=-999;beamvz=-999;beamvdis=-999;
    cdcvx=-999;cdcvy=-999;cdcvz=-999;cdcvdis=-999;
    tpid_beam=-999;


    bool secflag=false;

    if( head->pattern(12) >0 || head->pattern(13) >0 ||
	head->pattern(14) >0 || head->pattern(16) >0) secflag=true;
    if(!secflag) continue;

    //    blMan->CheckContainerSize();
    //    blMan->Calc(conf);
    //blTrackMan->
    //blTrackMan->DoTracking(blMan);
    //    cdsMan->Calc(conf);
    //    trackMan->Calc(conf);
    //cdsMan->CheckContainerSize();
    //  

    double ctmT0=0,euT0,edT0,emT0;

    int pid_beam;//0:pi 1:K 3:else
    int segT0;
    if(AllData)
	{
	  int nT0=0;
	  for( int i=0; i<blMan->nT0(); i++ ){
	    if( blMan->T0(i)->CheckRange() ) nT0++;
	  }
	  if( nT0<1 ) ctmT0=0;//continue;


	  for( int i=0; i<blMan->nT0(); i++ ){
	    if( blMan->T0(i)->CheckRange() ){
	      //	      blMan->T0(i)->Calc(conf);

	      ctmT0 = blMan->T0(i)->ctmean();	
	      euT0 = blMan->T0(i)->eu();
	      edT0 = blMan->T0(i)->ed();
	      emT0 = blMan->T0(i)->emean();
	      segT0=i;
	    }
	  }

	  int nBHD=0;
	  for( int i=0; i<blMan->nBHD(); i++ ){
	    if( blMan->BHD(i)->CheckRange() ) nBHD++;
	  }
	  //	  if( nBHD!=1 ) continue;
	  
	  double ctmBHD=0;
	  for( int i=0; i<blMan->nBHD(); i++ ){
	    if( blMan->BHD(i)->CheckRange() ){
	      //	      blMan->BHD(i)->Calc(conf);
	      ctmBHD = blMan->BHD(i)->ctmean();
	    }
	  }

	  double tofBHDT0 = ctmT0-ctmBHD;
	  
	  h1 = (TH1F*)gFile->Get("tofBHDT0");      
	  h1->Fill( tofBHDT0 );
	  h1 = (TH1F*)gFile->Get(Form("timeT0_%d",segT0));h1->Fill(ctmT0);
	  
       
	  if(0<tofBHDT0 && tofBHDT0<27.5) pid_beam=0;
	  else if(27.5<tofBHDT0 && tofBHDT0<32) pid_beam=1;
	  else pid_beam=3;
	  //if(pid_beam!=0) continue;

	  tpid_beam=pid_beam;

	}
    else
      {
       pid_beam=0;//0:pi 1:K 3:else
       segT0=1; 
       for( int i=0; i<blMan->nT0(); i++ ){
	 if( blMan->T0(i)->CheckRange() ){
	   //	   blMan->T0(i)->Calc(conf);

	      ctmT0 = blMan->T0(i)->tmean();	      
	      euT0 = blMan->T0(i)->eu();
	      edT0 = blMan->T0(i)->ed();
	      emT0 = blMan->T0(i)->emean();
	      segT0=i;
	      h1 = (TH1F*)gFile->Get(Form("timeT0_%d",segT0));h1->Fill(ctmT0);
	 }
       }

      }



      h1 = (TH1F*)gFile->Get("nTrack"); h1->Fill( trackMan->nGoodTrack() );


// 	for(int itr=0;itr<blTrackMan->ntrackBPC();itr++)
// 	  {
// 	    LocalTrack *bpc=blTrackMan->trackBPC(itr);

// 	    h1 = (TH1F*)gFile->Get(Form("BPCpos_T0_%d",segT0) ); h1->Fill(x1);
// 	    h1 = (TH1F*)gFile->Get(Form("chi_BPC") ); h1->Fill(bpc->chi2all());	
// 	    h1 = (TH1F*)gFile->Get(Form("chixz_BPC") ); h1->Fill(bpc->chi2xz());	
// 	    h1 = (TH1F*)gFile->Get(Form("chiyz_BPC") ); h1->Fill(bpc->chi2yz());	

// 	    double x1,y1,x2,y2;
// 	    double z=0,z2=20;
// 	    blc2->XYPosatZ(z,x1,y1);		  
// 	    blc2->XYPosatZ(z2,x2,y2);
// 	    double len=sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z)*(z2-z) );
// 	    double dx=(x2-x1)/len;   double dy=(y2-y1)/len;		 
// 	    h2 = (TH2F*)gFile->Get("bxy"); h2->Fill( x1,y1 );
// 	    h2 = (TH2F*)gFile->Get("bdxy"); h2->Fill( dx,dy );
// 	  }


      

      std::vector <int> cdsID;
      std::vector <double> mass;
      std::vector <int> pip_ID;
      std::vector <int> pim_ID;
      std::vector <int> km_ID;
      std::vector <int> p_ID;
      std::vector <int> d_ID;
      bool flaglambda=false;
      bool flagLp=false;
      bool flagCDC=false;
      bool flagnK=false;
      bool flagnK0s=false;
      double beta_b=0;
      TVector3 vtxNC;
      TVector3 vtxT0;
      TVector3 boost;
      TVector3 Pnk;
      TLorentzVector L3_p;
      TLorentzVector L3_L;
      TLorentzVector L3_target;
      TLorentzVector L3_beam;
      double TL;
      double Tp;
      double tQ;
      if(blTrackMan->ntrackBPC()!=1 ) continue;
      for( int it1=0; it1<trackMan->nGoodTrack(); it1++ )
	{
	  CDSTrack *track=trackMan->Track(trackMan->GoodTrackID(it1));	
	  bool cdhflag=true;
	  if(!track->CDHFlag()) {cdhflag=false;continue;}	  
	  HodoscopeLikeHit *cdh1;
	  cdh1=track->CDHHit();
	  //	  if(cdhflag) cdh1->Calc(conf);
	  //	  track->Calc(conf);
	  //	if(cdh1->seg()==-1) cdhflag=false;
	  
	  double tof1 = cdh1->ctmean()-ctmT0;
	
	  double param1[5];
	  
	  track->GetParameters(param1);
	//	if(param1[2]==0  ) continue;
	  double drho1=param1[0], phi01=param1[1], rho1=track->Rho(), dz1=param1[3], tlam1=param1[4];
	  double mom1 = track->Momentum();
	  int ptype1 = -1;
	  if(mom1!=0 )ptype1 = PID( tof1, mom1 );
	  
	  double mass1=0;	
	  if(ptype1==1 || ptype1==2) mass1=piMass;
	  else  if(ptype1==3) mass1=pMass;
	  //#######vertex beam#####	  
	  TrackVertex vertexb;
	  TVector3 vtxb1, vtxb2, vtxb;
	
	  if(trackMan->GetVertex_beam(trackMan->GoodTrackID(it1),0,vertexb))
	    {	 	      
	      vtxb=vertexb.GetVertexPos_mean();
	      vtxb1=vertexb.GetVertexPos1();
	      vtxb2=vertexb.GetVertexPos2();
	      beamvx=vtxb.x();
	      beamvy=vtxb.y();
	      beamvz=vtxb.z();
	      TVector3 tmpbv=vtxb2-vtxb1;
	      double dis=tmpbv.Mag();
	      beamvdis=dis;
	      h1 = (TH1F*)gFile->Get("vrxb_dis"); h1->Fill( dis );
	      if(dis<2&& sqrt(vtxb.x()*vtxb.x()+vtxb.y()*vtxb.y())<3 && fabs(vtxb.z())<6)
		{
		  vtxNC=vtxb;
		  flagCDC=true;
		}      
	    }
	  
	//############## CDH #####################//
	  if(cdhflag && blTrackMan->ntrackBPC()==1 ) 
	    {
	      
	      LocalTrack *bpc=blTrackMan->trackBPC(0);
	      double x1,y1,x2,y2;
	      double gx,gy,gz,dx,dy,dz;
	      //   if(blc2->chi2all()>1) continue;
	      conf->GetGeomMapManager()->GetGParam(CID_T0,gx,gy,gz,dx,dy,dz);
	      //	    std::cout<<"gz "<<gz<<std::endl; 
	      double z_t0=gz,z_vtxb=vtxb.Z();
	      bpc->XYPosatZ(z_t0,x1,y1);bpc->XYPosatZ(z_vtxb,x2,y2);
	      double beam_dis=sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z_t0-z_vtxb)*(z_t0-z_vtxb));
	      vtxT0.SetXYZ(x1,y1,z_t0);
	      if(AllData)
		{
		  if(pid_beam==0)  beta_b=1.0/sqrt(1.0*1.0+piMass*piMass);
		  else if(pid_beam==1) beta_b=1.0/sqrt(1.0*1.0+kpMass*kpMass);
		}
	      else    beta_b=1.0/sqrt(1.0*1.0+piMass*piMass);
	      double time_beam=beam_dis/beta_b/Const;
	      
	      //#####CDC dis#######		  

	      TVector3 cdhvtx=track->CDHVertex();
	      double phi1,phi2;
	      phi1=track->CalcHelixPhi(vtxb.x(),vtxb.y());
	      phi2=track->CalcHelixPhi(cdhvtx.x(),cdhvtx.y());
	      TVector3 tmp=cdhvtx-vtxb;
	      //	    double cdc_dis=tmp.Mag();
	      double cdc_dis=rho1*fabs(phi2-phi1)*sqrt(1+tlam1*tlam1);
	      double beta_calc=cdc_dis/(tof1-time_beam)/Const;
	      double calc_mass1=fabs(mom1)*sqrt(1/(beta_calc*beta_calc)-1);
	      double calc_mass2=mom1*mom1*(1/(beta_calc*beta_calc)-1);
	      
	      //	    if(ptype1==1 || ptype1==2) mass1= piMass; // 1,2:pi ,3: proton
	      //	    else if(ptype1==3 ) mass1= pMass; // 1,2:pi ,3: proton
	      if(calc_mass1>0.04 && calc_mass1<0.36 ) 
		{
		  mass1=piMass;
		  if(mom1<0)
		    {pim_ID.push_back(trackMan->GoodTrackID(it1));
		    track->SetPID(1);}
		  else 
		  {  pip_ID.push_back(trackMan->GoodTrackID(it1));
		  track->SetPID(2);}
		}
	      else  if(calc_mass1>0.36 && calc_mass1<0.61 ) 
		{
		  mass1=kpMass;
		  if(mom1<0)
		   { km_ID.push_back(trackMan->GoodTrackID(it1));
		   track->SetPID(3);
		   h1 = (TH1F*)gFile->Get(Form("mom_k" )); h1->Fill(mom1);
		   }
		  else
		    {
		      //     pip_ID.push_back(trackMan->GoodTrackID(it1));
		      track->SetPID(4);}
		}
	      else  if(calc_mass1>0.61 && calc_mass1<1.2 )
		{
		  mass1=pMass;
		  if(mom1<0)
		    {		     
		      //     pip_ID.push_back(trackMan->GoodTrackID(it1));
		      track->SetPID(6);}
		  else
		  {  p_ID.push_back(trackMan->GoodTrackID(it1));
		  track->SetPID(5);}
		}
	      else  if(calc_mass1>1.2 && calc_mass1<4.0 ) 
		{
		  mass1=dMass;
		  if(mom1<0)
		    {
		      track->SetPID(7);
		    }
		  else
		    {  d_ID.push_back(trackMan->GoodTrackID(it1));
		  track->SetPID(8);}
		  h1 = (TH1F*)gFile->Get(Form("mom_d" )); h1->Fill(mom1);
		}
	      else if(calc_mass1<0.04)
		{
		  track->SetPID(0);
		}
	      else
		{
		  track->SetPID(9);
		}
	      double beta=fabs(mom1)/sqrt(mass1*mass1+mom1*mom1);
	      double time_cdc=cdc_dis/beta/Const;
	      	      
	      //####Calc time#####
	      double calc_time=time_cdc+time_beam;
	      double calc_time2=cdc_dis/Const+time_beam;

	      double ctsub=cdh1->ctsub();
	      double calc_time_u=time_cdc+time_beam+ctsub;
	      double calc_time_d=time_cdc+time_beam-ctsub;
	      //	    std::cout<<"CDH seg= "<<cdh1->seg()<<std::endl;	    	      
	      if( pid_beam==0)
		{
		  h2 = (TH2F*)gFile->Get(Form("mom_overbeta" )); h2->Fill(1./beta_calc,mom1 );
		  h2 = (TH2F*)gFile->Get(Form("mass2_mom" )); h2->Fill(calc_mass2,mom1 );
		  h2 = (TH2F*)gFile->Get(Form("mass_mom" )); h2->Fill(calc_mass1,mom1 );
		}
	      if( pid_beam==1)
		{
		  h2 = (TH2F*)gFile->Get(Form("mom_overbeta2" )); h2->Fill(1./beta_calc,mom1 );
		  h2 = (TH2F*)gFile->Get(Form("mass2_mom2" )); h2->Fill(calc_mass2,mom1 );
		  h2 = (TH2F*)gFile->Get(Form("mass_mom2" )); h2->Fill(calc_mass1,mom1 );
		  if(cdh1->seg()==3)
		    {
		      h2 = (TH2F*)gFile->Get(Form("mom_overbeta_tmp" )); h2->Fill(1./beta_calc,mom1 );
		      h2 = (TH2F*)gFile->Get(Form("mass2_mom_tmp" )); h2->Fill(calc_mass2,mom1 );
		      h2 = (TH2F*)gFile->Get(Form("mass_mom_tmp" )); h2->Fill(calc_mass1,mom1 );
		    }

		}

	      if(sqrt(vtxb.x()*vtxb.x()+vtxb.y()*vtxb.y() )<3 && fabs(vtxb.z())<6)
		{

		  if(pid_beam==1)
		    {
		      h2 = (TH2F*)gFile->Get(Form("mom_overbeta3" )); h2->Fill(1./beta_calc,mom1 );
		      h2 = (TH2F*)gFile->Get(Form("mass2_mom3" )); h2->Fill(calc_mass2,mom1 );
		      h2 = (TH2F*)gFile->Get(Form("mass_mom3" )); h2->Fill(calc_mass1,mom1 );
		      if(mom1>0)
			{
		      h2 = (TH2F*)gFile->Get(Form("mass2p_mom3" )); h2->Fill(calc_mass2,mom1 );
		      h2 = (TH2F*)gFile->Get(Form("massp_mom3" )); h2->Fill(calc_mass1,mom1 );
			}
		      if(mom1<0)
			{
		      h2 = (TH2F*)gFile->Get(Form("mass2n_mom3" )); h2->Fill(calc_mass2,mom1 );
		      h2 = (TH2F*)gFile->Get(Form("massn_mom3" )); h2->Fill(calc_mass1,mom1 );
			}
		    }
		}
	      
	      
	    }//CDH
	  mass.push_back(mass1);
	  cdsID.push_back(trackMan->GoodTrackID(it1));
	}//trackMan


      //####deuteron
      if((int)d_ID.size()>0)
	{
	  LocalTrack *bpc=blTrackMan->trackBPC(0);
	  h1 = (TH1F*)gFile->Get(Form("ntrack_d" )); h1->Fill(trackMan->nGoodTrack());
	  for( int it=0; it<trackMan->nGoodTrack(); it++ )
	    {
	      CDSTrack *track=trackMan->Track(trackMan->GoodTrackID(it));
	      double mom=track->Momentum();
	      h2 = (TH2F*)gFile->Get(Form("pid_pat_d" )); h2->Fill(trackMan->nGoodTrack(),track->PID());
	      if(track->PID()==8)
		{

		  //#######vertex beam#####	  
		  TrackVertex vertexb;
		  TVector3 vtxb1, vtxb2, vtxb;
		  
		  if(trackMan->GetVertex_beam(trackMan->GoodTrackID(it),0,vertexb))
		    {	 	      
		      vtxb=vertexb.GetVertexPos_mean();
		      vtxb1=vertexb.GetVertexPos1();
		      vtxb2=vertexb.GetVertexPos2();
		      beamvx=vtxb.x();
		      beamvy=vtxb.y();
		      beamvz=vtxb.z();
		      TVector3 tmpbv=vtxb2-vtxb1;
		      double dis=tmpbv.Mag();
		      beamvdis=dis;
		    }

		  TVector3 Pd;
		  if( !track->GetMomentum(vtxb1, Pd) )
		    std::cout<<"error of GetMomentum() track"<<std::endl;;
		  
		  double x1,y1,x2,y2;
		  double z1=0,z2=20;
		  bpc->XYPosatZ(z1,x1,y1);		  
		  bpc->XYPosatZ(z2,x2,y2);
		  TVector3 lp;lp.SetXYZ(x1,y1,z1);
		  TVector3 ls;ls.SetXYZ(x2-x1,y2-y1,z2-z1);ls=ls.Unit();

		  TVector3 Pp_beam=1.0*ls; 
		  TVector3 Pp_target;Pp_target.SetXYZ(0,0,0); 
		  TLorentzVector L_beam; L_beam.SetVectM(Pp_beam , kpMass );
		  TLorentzVector L_target; L_target.SetVectM( Pp_target, ThreeHeMass);
		  boost=(L_beam+L_target).BoostVector();

		  double cos_d=Pd.Dot(Pp_beam)/Pd.Mag()/Pp_beam.Mag();
		  h1 = (TH1F*)gFile->Get(Form("cos_d" )); h1->Fill(cos_d);
		  h2 = (TH2F*)gFile->Get(Form("mom_cos_d" )); h2->Fill(mom,cos_d);
		  TLorentzVector L_d; L_d.SetVectM(Pd,dMass);
		  TLorentzVector L_x; L_x=L_beam+L_target-L_d;
		  h1 = (TH1F*)gFile->Get("Mmd"); h1->Fill( L_x.M() );

		}
	    }

	}



      //####kaon
      if((int)km_ID.size()>0)
	{
	  LocalTrack *bpc=blTrackMan->trackBPC(0);
	  h1 = (TH1F*)gFile->Get(Form("ntrack_k" )); h1->Fill(trackMan->nGoodTrack());
	  for( int it=0; it<trackMan->nGoodTrack(); it++ )
	    {
	      CDSTrack *track=trackMan->Track(trackMan->GoodTrackID(it));
	      double mom=fabs(track->Momentum());
	      h2 = (TH2F*)gFile->Get(Form("pid_pat_k" )); h2->Fill(trackMan->nGoodTrack(),track->PID());
	      if(track->PID()==3)
		{

		  //#######vertex beam#####	  
		  TrackVertex vertexb;
		  TVector3 vtxb1, vtxb2, vtxb;
		  
		  if(trackMan->GetVertex_beam(trackMan->GoodTrackID(it),0,vertexb))
		    {	 	      
		      vtxb=vertexb.GetVertexPos_mean();
		      vtxb1=vertexb.GetVertexPos1();
		      vtxb2=vertexb.GetVertexPos2();
		      beamvx=vtxb.x();
		      beamvy=vtxb.y();
		      beamvz=vtxb.z();
		      TVector3 tmpbv=vtxb2-vtxb1;
		      double dis=tmpbv.Mag();
		      beamvdis=dis;
		    }

		  TVector3 Pk;
		  if( !track->GetMomentum(vtxb1, Pk) )
		    std::cout<<"error of GetMomentum() track"<<std::endl;;
		  
		  double x1,y1,x2,y2;
		  double z1=0,z2=20;
		  bpc->XYPosatZ(z1,x1,y1);		  
		  bpc->XYPosatZ(z2,x2,y2);
		  TVector3 lp;lp.SetXYZ(x1,y1,z1);
		  TVector3 ls;ls.SetXYZ(x2-x1,y2-y1,z2-z1);ls=ls.Unit();

		  TVector3 Pp_beam=1.0*ls; 
		  TVector3 Pp_target;Pp_target.SetXYZ(0,0,0); 
		  TLorentzVector L_beam; L_beam.SetVectM(Pp_beam , kpMass );
		  TLorentzVector L_target; L_target.SetVectM( Pp_target, ThreeHeMass);
		   boost=(L_beam+L_target).BoostVector();

		  double cos_k=Pk.Dot(Pp_beam)/Pk.Mag()/Pp_beam.Mag();
		  h1 = (TH1F*)gFile->Get(Form("cos_k" )); h1->Fill(cos_k);
		  h2 = (TH2F*)gFile->Get(Form("mom_cos_k" )); h2->Fill(mom,cos_k);

		  TLorentzVector L_k; L_k.SetVectM(Pk , kpMass );
		  TVector3 Ppr;Ppr.SetXYZ(0,0,0);//GeV/c
		  TLorentzVector L_proton; L_proton.SetVectM( Ppr, pMass);
		  TLorentzVector L_x; L_x=L_beam+L_proton-L_k;
		  h1 = (TH1F*)gFile->Get("MmK"); h1->Fill( L_x.M() );
		  if(trackMan->nGoodTrack()==1) {h1 = (TH1F*)gFile->Get("MmK_1"); h1->Fill( L_x.M() );}
		  if(L_x.M()>0.85 && L_x.M()<1.1 && trackMan->nGoodTrack()==1)
		    {
		      h1 = (TH1F*)gFile->Get("nK_cos"); h1->Fill( L_x.CosTheta() );
		      h1 = (TH1F*)gFile->Get("nK_dxy"); h1->Fill( L_x.Px()/L_x.Pz(),L_x.Py()/L_x.Pz() );
		      h1 = (TH1F*)gFile->Get("nK_mom"); h1->Fill( L_x.P() );
		      Pnk=L_x.Vect();
		      double ztmp=1483-vtxNC.Z();
		      double xtmp=ztmp*Pnk.X()/Pnk.Z()+vtxNC.X();
		      double ytmp=ztmp*Pnk.Y()/Pnk.Z()+vtxNC.Y();
		      //		      std::cout<<"K x "<<xtmp<<std::endl;
		      h2 = (TH2F*)gFile->Get("NCxy_nK"); h2->Fill( xtmp,ytmp );
		      if(fabs(xtmp)<160 && fabs(ytmp)<75 && flagCDC)
			{
			  flagnK=true;
			}
		    }

		}
	    }

	}



      //#########Invariant mass############
      //###########ppi#######
      bool lambdaflag=false;
      for(int np=0;np<(int)p_ID.size();np++)
	{
	  for(int npi=0;npi<(int)pim_ID.size();npi++)
	    {
	      CDSTrack *track_p=trackMan->Track(p_ID[np] );	
	      CDSTrack *track_pim=trackMan->Track(pim_ID[npi]);	
	      //#######vertex#####	  
	      TrackVertex vertex;
	      TVector3 vtx_p, vtx_pim, vtx;
	      if(trackMan->GetVertex(pim_ID[npi],p_ID[np],vertex)  )
		{	 
		  vtx=vertex.GetVertexPos_mean();
		  if( p_ID[np]<pim_ID[npi] )
		    {		 
		      vtx_p=vertex.GetVertexPos1();
		      vtx_pim=vertex.GetVertexPos2();
		    }
		  else
		    {		 
		      vtx_p=vertex.GetVertexPos2();
		      vtx_pim=vertex.GetVertexPos1();
		    }

		}
	      //###################

	      TVector3 Pp_p;
	      if( !track_p->GetMomentum(vtx_p, Pp_p) )
		std::cout<<"error of GetMomentum() track proton"<<std::endl;;
	      TVector3 Pp_pim;
	      if( !track_pim->GetMomentum(vtx_pim, Pp_pim) )
		std::cout<<"error of GetMomentum() track pim "<<std::endl;;
	      
	      TVector3 Pp = Pp_p + Pp_pim;
	      TLorentzVector L_p; L_p.SetVectM( Pp_p, pMass );
	      TLorentzVector L_pim; L_pim.SetVectM( Pp_pim, piMass);
	      
	      double cosOA = (Pp_p.Dot(Pp_pim))/Pp_p.Mag()/Pp_pim.Mag();
	      double im = (L_p+L_pim).M();

	      h1 = (TH1F*)gFile->Get("cosOAppi"); h1->Fill( cosOA );
	      h1 = (TH1F*)gFile->Get("IMppi"); h1->Fill( im );
	      h2 = (TH2F*)gFile->Get("IMcosOAppi"); h2->Fill( im, cosOA );
 	      h1 = (TH1F*)gFile->Get("costppi"); h1->Fill( Pp.CosTheta() );
 	      h1 = (TH1F*)gFile->Get("phippi"); h1->Fill( Pp.Phi()*TMath::RadToDeg() );
 	      h2 = (TH2F*)gFile->Get("IMcostppi"); h2->Fill( im, Pp.CosTheta() );
	      //#####lambda background Cut######
		
	      LocalTrack *bpc=blTrackMan->trackBPC(0);
	      double x1,y1,x2,y2;
	      double z1=0,z2=20;
	      bpc->XYPosatZ(z1,x1,y1);		  
	      bpc->XYPosatZ(z2,x2,y2);
	      TVector3 lp;lp.SetXYZ(x1,y1,z1);
	      TVector3 ls;ls.SetXYZ(x2-x1,y2-y1,z2-z1);ls=ls.Unit();

	      TVector3 Pp_beam=1.0*ls; 
	      //		  TVector3 Pp_beam;Pp_beam.SetXYZ(0,0,1.0); 
	      TVector3 Pp_target;Pp_target.SetXYZ(0,0,0); 
	      TLorentzVector L_beam; L_beam.SetVectM(Pp_beam , kpMass );
	      TLorentzVector L_target; L_target.SetVectM( Pp_target, ThreeHeMass);
	      boost=(L_beam+L_target).BoostVector();
	      L3_target=L_target;
	      L3_beam=L_beam;
	      double dis,dist,dltmp=0;TVector3 xest,nest;
	      //		  track1->PointToLine( vtx, lp, ls,dis,xest );
	      track_p->LineToLine( vtx,Pp.Unit(),lp, ls,dltmp,dist,xest,nest );
	      vtxNC=xest;
	      dis=(xest-vtx).Mag();//dis/(sin(Pp.Theta() ));
	      h1 = (TH1F*)gFile->Get("IMppi_vdis"); h1->Fill( dist );
	      h2 = (TH2F*)gFile->Get("IMppi_dis"); h2->Fill( im,dis );
	      h2 = (TH2F*)gFile->Get("IMppi_xest_yest"); h2->Fill( xest.x(),xest.y() );
	      h2 = (TH2F*)gFile->Get("IMppi_yest_zest"); h2->Fill( xest.z(),xest.y() );
		  

	      if(sqrt( pow( (xest.x()),2)+pow( (xest.y()),2) )<3.5 && fabs(xest.z())<6  && pid_beam==1  )
		{

		  //		      if(dis>2 /*ct=7.89cm */  )
		  if(dis>2.5){
		    h1 = (TH1F*)gFile->Get("IMppi_cut"); h1->Fill( im );
		  }

		  h1 = (TH1F*)gFile->Get("IMppi_target"); h1->Fill( im );
		  if(im>1.12 || im<1.108)
		    {
		      if((int)vLv_p.size()>20 && dis>2)
			{
			  std::vector<TLorentzVector>::iterator itp;
			  std::vector<TLorentzVector>::iterator itpim;
			  itp=vLv_p.begin();
			  itpim=vLv_pim.begin();
			  vLv_p.erase(itp);
			  vLv_pim.erase(itpim);
			}
		      vLv_p.push_back(L_p);
		      vLv_pim.push_back(L_pim);
		    }

		  else //lambdapeak
		    {
		      flaglambda=true;
		      TLorentzVector L_L; L_L.SetVectM( Pp, lMass);  
		      double tof_dis=dis;//tmp.Mag();
		      double beta=Pp.Mag()/sqrt(Pp.Mag()*Pp.Mag()+lMass*lMass);
		      
		      double gamma=L_L.Gamma();
		      double tof=tof_dis/beta/Const/gamma;
		      
		      h1 = (TH1F*)gFile->Get(Form("ntrack_L" )); h1->Fill(trackMan->nGoodTrack());
		      for( int it=0; it<trackMan->nGoodTrack(); it++ )
			{
			  if(trackMan->GoodTrackID(it)==p_ID[np] || trackMan->GoodTrackID(it)==pim_ID[npi] ) continue;
			  CDSTrack *track=trackMan->Track(trackMan->GoodTrackID(it));
			  h2 = (TH2F*)gFile->Get(Form("pid_pat_L" )); h2->Fill(trackMan->nGoodTrack()-2,track->PID());			  
			}

		      double cos_L=Pp_beam.Dot(Pp)/Pp_beam.Mag()/Pp.Mag();
		      h1 = (TH1F*)gFile->Get("cosLambda"); h1->Fill( cos_L );
		      h1 = (TH1F*)gFile->Get("pLambda"); h1->Fill( Pp.Mag() );
		      h2 = (TH2F*)gFile->Get("mom_cos_L"); h2->Fill(Pp.Mag(), cos_L  );

		      h1 = (TH1F*)gFile->Get("LifetimeLambda"); h1->Fill( tof );


		      h1 = (TH1F*)gFile->Get("pL_proton"); h1->Fill( Pp_p.Mag() );
		      h1 = (TH1F*)gFile->Get("pL_pi"); h1->Fill( Pp_pim.Mag() ); 


		      TVector3 Ppr;Ppr.SetXYZ(0,0,0);//GeV/c
		      TLorentzVector L_proton; L_proton.SetVectM( Ppr, pMass);
		      TLorentzVector L_x; L_x=L_beam+L_proton-L_L;
		      h1 = (TH1F*)gFile->Get("MmL"); h1->Fill( L_x.M() );      
		      if(dis>2){		      h1 = (TH1F*)gFile->Get("MmL_cut"); h1->Fill( L_x.M() );      }
		      //		      TLorentzVector L_x; L_x.SetVect( Pp_beam-Pp );L_x.SetE( L_beam.E()+L_proton.E()-L_L.E() );

		      if((int)p_ID.size()==2 && trackMan->nGoodTrack()==3)
			{
			  int np2=0;
			  if(np==0) np2=1;
			  else np2=0;
			  CDSTrack *track3 = trackMan->Track(p_ID[np2]);
			  //			  track3->Calc(conf);
			  if(!track3->CDHFlag()) continue;
			  HodoscopeLikeHit *cdh3 = track3->CDHHit();
			  //			  cdh3->Calc(conf);
#if SIM
			  double tof3 = cdh3->tmean()-ctmT0;
#else
			  double tof3 = cdh3->ctmean()-ctmT0;
#endif		      
			  double param3[5];
			  track3->GetParameters(param3);		      
			  double drho3=param3[0], phi03=param3[1], rho3=track3->Rho(), dz3=param3[3], tlam3=param3[4];
			  double mom3 = track3->Momentum();
			  TVector3 lnest,hnest;
			  double dcdis;
			  track3->LineToHelix(vtx,Pp,param3,lnest,hnest,dcdis);
			  TVector3 Pp3;
			  TVector3 tmpline=vtx-hnest;
			  if(dcdis<1.0)
			    { 				      
			      if( !track3->GetMomentum(hnest, Pp3) )
				std::cout<<"error of GetMomentum() track3 "<<iev<<std::endl;
			      TLorentzVector L_p3; L_p3.SetVectM( Pp3, pMass );
			      double im_kpp = (L_L+L_p3).M();	
			      double cosOAkpp = (Pp.Dot(Pp3))/Pp.Mag()/Pp3.Mag();
			      //   TVector3 boost2=-1*(L_L+L_p3).BoostVector();
			      //     L_L.Boost(boost2); L_p3.Boost(boost2);


			      double CMcosOAkpp = ( L_L.Vect().Dot( L_p3.Vect() ) )/L_L.Vect().Mag()/L_p3.Vect().Mag();
			      TLorentzVector L_x; L_x=L_beam+L_target-L_L-L_p3;
			      
			      h1 = (TH1F*)gFile->Get("IMkpp"); h1->Fill( im_kpp );
			      h1 = (TH1F*)gFile->Get("MmLp"); h1->Fill( L_x.M() );
			      h2 = (TH2F*)gFile->Get("IMcosOAkpp"); h2->Fill( im_kpp, cosOAkpp );

			      if(cosOAkpp<-0.8)
				{h1 = (TH1F*)gFile->Get("IMkpp_btb"); h1->Fill( im_kpp );}
			      if(L_x.M()<1.1)
				{h1 = (TH1F*)gFile->Get("IMkpp_mn"); h1->Fill( im_kpp );}
			      double tmppz=L_p3.Pz();
			      //			      std::cout<<"pz pzcm "<<tmppz<<" "<<L_p3.Pz()<<std::endl;
			      L3_p=L_p3;
			      L3_L=L_L;
			      flagLp=true;

			    }//vertex dis <1
			}//2nd proton
		    }//lambdapeak
		  

		}
	    }
	}
      //background study 
      //      std::cout<<"bg "<<vLv_p.size()<<" "<<vLv_pim.size()<<std::endl;
      if((int)vLv_p.size()>15)
	{
	  for(int n=0;n<2;n++)
	    {
	      int itlp=0;
	      itlp=rand()/(double)RAND_MAX*( (int)vLv_p.size()-1);
	      int itlpim=0;
	      //	      do
		{
		  itlpim=rand()/(double)RAND_MAX*( (int)vLv_pim.size()-1);
		}
		//	      while(itlpim==itlp);
		//      std::cout<<"bg test1 "<<itlp<<" "<<itlpim<<std::endl;
	      //      double cosOAbc = (Pp_p.Dot(Pp_pim))/Pp_p.Mag()/Pp_pim.Mag();
	      double imbg = (vLv_p[itlp]+vLv_pim[itlpim]).M();
	      h1 = (TH1F*)gFile->Get("IMppi_bg2"); h1->Fill( imbg );
	    }
	}


      //###############pipi
      for(int np=0;np<(int)pip_ID.size();np++)
	{
	  for(int nm=0;nm<(int)pim_ID.size();nm++)
	    {
	      CDSTrack *track_p=trackMan->Track(pip_ID[np] );	
	      CDSTrack *track_m=trackMan->Track(pim_ID[nm]);	
	      //#######vertex#####	  
	      TrackVertex vertex;
	      TVector3 vtx_p, vtx_m, vtx;
	      if(trackMan->GetVertex(pim_ID[nm],pip_ID[np],vertex)  )
		{	 
		  vtx=vertex.GetVertexPos_mean();
		  if( pip_ID[np]<pim_ID[nm] )
		    {		 
		      vtx_p=vertex.GetVertexPos1();
		      vtx_m=vertex.GetVertexPos2();
		    }
		  else
		    {		 
		      vtx_p=vertex.GetVertexPos2();
		      vtx_m=vertex.GetVertexPos1();
		    }

		}
	      //###################

	      TVector3 Pp_p;
	      if( !track_p->GetMomentum(vtx_p, Pp_p) )
		std::cout<<"error of GetMomentum() track pip"<<std::endl;;
	      TVector3 Pp_m;
	      if( !track_m->GetMomentum(vtx_m, Pp_m) )
		std::cout<<"error of GetMomentum() track pim "<<std::endl;;
	      
	      TVector3 Pp = Pp_p + Pp_m;
	      TLorentzVector L_p; L_p.SetVectM( Pp_p, piMass );
	      TLorentzVector L_m; L_m.SetVectM( Pp_m, piMass);
	      
	      double cosOA = (Pp_p.Dot(Pp_m))/Pp_p.Mag()/Pp_m.Mag();
	      double im = (L_p+L_m).M();

	      h1 = (TH1F*)gFile->Get("cosOApipi"); h1->Fill( cosOA );
	      h1 = (TH1F*)gFile->Get("IMpipi"); h1->Fill( im );
	      h2 = (TH2F*)gFile->Get("IMcosOApipi"); h2->Fill( im, cosOA );
 	      h1 = (TH1F*)gFile->Get("costpipi"); h1->Fill( Pp.CosTheta() );
 	      h1 = (TH1F*)gFile->Get("phipipi"); h1->Fill( Pp.Phi()*TMath::RadToDeg() );
 	      h2 = (TH2F*)gFile->Get("IMcostpipi"); h2->Fill( im, Pp.CosTheta() );
	      //#####K0s background Cut######
		
	      LocalTrack *bpc=blTrackMan->trackBPC(0);
	      double x1,y1,x2,y2;
	      double z1=0,z2=20;
	      bpc->XYPosatZ(z1,x1,y1);		  
	      bpc->XYPosatZ(z2,x2,y2);
	      TVector3 lp;lp.SetXYZ(x1,y1,z1);
	      TVector3 ls;ls.SetXYZ(x2-x1,y2-y1,z2-z1);ls=ls.Unit();

	      TVector3 Pp_beam=1.0*ls; 
	      //		  TVector3 Pp_beam;Pp_beam.SetXYZ(0,0,1.0); 
	      TVector3 Pp_target;Pp_target.SetXYZ(0,0,0); 
	      TLorentzVector L_beam; L_beam.SetVectM(Pp_beam , kpMass );
	      TLorentzVector L_target; L_target.SetVectM( Pp_target, ThreeHeMass);
	      boost=(L_beam+L_target).BoostVector();

	      double dis,dist,dltmp=0;TVector3 xest,nest;
	      //		  track1->PointToLine( vtx, lp, ls,dis,xest );
	      track_p->LineToLine( vtx,Pp.Unit(),lp, ls,dltmp,dist,xest,nest );
	      vtxNC=xest;
	      dis=(xest-vtx).Mag();//dis/(sin(Pp.Theta() ));
	      h1 = (TH1F*)gFile->Get("IMpipi_vdis"); h1->Fill( dist );
	      h2 = (TH2F*)gFile->Get("IMpipi_dis"); h2->Fill( im,dis );

	      if(sqrt( pow( (xest.x()),2)+pow( (xest.y()),2) )<3.5 && fabs(xest.z())<6  && pid_beam==1  )
		{
		  h1 = (TH1F*)gFile->Get("IMpipi_target"); h1->Fill( im );

		  if(dis>1){
		    h1 = (TH1F*)gFile->Get("IMpipi_cut"); h1->Fill( im );
		  }// 		  if(im>1.12 || im<1.108)
// 		    {
// 		      if((int)vLv_p.size()>20)
// 			{
// 			  std::vector<TLorentzVector>::iterator itp;
// 			  std::vector<TLorentzVector>::iterator itpim;
// 			  itp=vLv_p.begin();
// 			  itpim=vLv_pim.begin();
// 			  vLv_p.erase(itp);
// 			  vLv_pim.erase(itpim);
// 			}
// 		      vLv_p.push_back(L_p);
// 		      vLv_pim.push_back(L_pim);
//		    }
	
		  //K0s peak
	      if( (0.49-0.02)<im && im<(0.49+0.02 ))
		{

		  TLorentzVector L_K0s; L_K0s.SetVectM( Pp, k0sMass);  
		  double tof_dis=dis;//tmp.Mag();
		  double beta=Pp.Mag()/sqrt(Pp.Mag()*Pp.Mag()+k0sMass*k0sMass);
		  double gamma=L_K0s.Gamma();
		  double tof=tof_dis/beta/Const/gamma;		      
		  h1 = (TH1F*)gFile->Get("ntrack_k0s"); h1->Fill(trackMan->nGoodTrack()  );
		  for( int it=0; it<trackMan->nGoodTrack(); it++ )
		    {
		      if(trackMan->GoodTrackID(it)==pip_ID[np] || trackMan->GoodTrackID(it)==pim_ID[nm] ) continue;
		      CDSTrack *track=trackMan->Track(trackMan->GoodTrackID(it));
		      h2 = (TH2F*)gFile->Get(Form("pid_pat_k0s" )); h2->Fill(trackMan->nGoodTrack()-2,track->PID());			  
		    }

		  double cos_k0s=Pp_beam.Dot(Pp)/Pp_beam.Mag()/Pp.Mag();
		  h1 = (TH1F*)gFile->Get("cosK0s"); h1->Fill( cos_k0s );
		  h2 = (TH2F*)gFile->Get("mom_cos_k0s"); h2->Fill(Pp.Mag(), cos_k0s  );
		  
		  h1 = (TH1F*)gFile->Get("pK0s_pip"); h1->Fill( Pp_p.Mag() );
		  h1 = (TH1F*)gFile->Get("pK0s_pim"); h1->Fill( Pp_m.Mag() );
		  
		  h1 = (TH1F*)gFile->Get("pK0s"); h1->Fill( Pp.Mag() );
		  h1 = (TH1F*)gFile->Get("LifetimeK0s"); h1->Fill( tof );
		  h1 = (TH1F*)gFile->Get("disK0s"); h1->Fill( tof_dis );
		  
		  TVector3 Ppr;Ppr.SetXYZ(0,0,0);//GeV/c
		  TLorentzVector L_proton; L_proton.SetVectM( Ppr, pMass);
		  TLorentzVector L_x; L_x=L_beam+L_proton-L_K0s;
		  h1 = (TH1F*)gFile->Get("MmK0s"); h1->Fill( L_x.M() );
		  if(dis>1 && trackMan->nGoodTrack()){	
		    h1 = (TH1F*)gFile->Get("MmK0s_cut"); h1->Fill( L_x.M() );
		    
		    if(L_x.M()>0.85 && L_x.M()<1.1 && trackMan->nGoodTrack()==2)
		      {

			h1 = (TH1F*)gFile->Get("nK0s_cos"); h1->Fill( L_x.CosTheta() );
			h1 = (TH1F*)gFile->Get("nK0s_dxy"); h1->Fill( L_x.Px()/L_x.Pz(),L_x.Py()/L_x.Pz() );
			h1 = (TH1F*)gFile->Get("nK0s_mom"); h1->Fill( L_x.P() );

			Pnk=L_x.Vect();
			double ztmp=1483-vtxNC.Z();
			double xtmp=ztmp*Pnk.X()/Pnk.Z()+vtxNC.X();
			double ytmp=ztmp*Pnk.Y()/Pnk.Z()+vtxNC.Y();
			//			std::cout<<"K0s x y "<<xtmp<<" "<<ytmp<<std::endl;
			h2 = (TH2F*)gFile->Get("NCxy_nK0s"); h2->Fill( xtmp,ytmp );

			if(fabs(xtmp)<160 && fabs(ytmp)<75 )
			{
			  flagnK0s=true;
			  h1 = (TH1F*)gFile->Get("NCeff_nK0s"); h1->Fill( xtmp );
			}
		      }
		  }
		  
		}
	      
		}
	    }
	}
      // #######NC spectrrum##############


      bool chargedhit=false;
      for(int i=0;i<blMan->nTOF();i++)
	{
	  HodoscopeLikeHit *hit=blMan->TOF(i);
	  if(hit->CheckRange() )
	    {
	      chargedhit=true;
	    }
	}
      for(int i=0;i<blMan->nCV();i++)
	{
	  HodoscopeLikeHit *hit=blMan->CV(i);
	  if(hit->CheckRange() )
	    {
	      chargedhit=true;
	    }
	}


      for( int i=0; i<blMan->nT0(); i++ ){
	if( blMan->T0(i)->CheckRange() ){
	  blMan->T0(i)->Calc(conf);	  
	  ctmT0 = blMan->T0(i)->ctmean();	
	  euT0 = blMan->T0(i)->eu();
	  edT0 = blMan->T0(i)->ed();
	  emT0 = blMan->T0(i)->emean();
	  segT0=i;
	}
      }

      int segNC=-1;
      int layerNC;
      double ctmNC=999;
      double ctsubNC=999;
      double euNC=999;
      double edNC=999;
      double emNC=999;
      double gx,gz;
      double gy=0;
      TVector3 gNC;
      if(!chargedhit)
	{
	  
	  for(int inc=0;inc<blMan->nNC();inc++)
	    {
	      HodoscopeLikeHit *hit=blMan->NC(inc);
	      
	      if(hit->CheckRange() )
		{		  
		  
		  hit->Calc(conf);
		  if( hit->ctmean() <ctmNC)
		    {
		      segNC=hit->seg();
		      if(segNC%16==0) layerNC=segNC/16;
		      else layerNC=(segNC+16-segNC%16)/16;

		      ctmNC = hit->ctmean();			 
		      ctsubNC = hit->ctsub();			 
		      euNC = hit->eu();
		      edNC = hit->ed();
		      emNC = hit->emean();
		      gx=150-((segNC-(layerNC-1)*16)-1)*20;
		      gz=1504.2+(layerNC-4)*7;
		      gy=ELV*ctsubNC;
		      gNC.SetXYZ(gx,gy,gz);
		      //		      std::cout<<"NC xyz "<<gx<<" "<<gy<<" "<<gz<<std::endl;
		    }
		}
	    }
	  if(segNC!=-1)
	    { 
	      
	      h1 = (TH1F*)gFile->Get("tofNC"); h1->Fill( ctmNC-ctmT0 );
	      //	  if(pid_beam) {h1 = (TH1F*)gFile->Get("tofNC"); h1->Fill( ctmNC-ctmT0 );}
	      //	  if(pid_beam) {h1 = (TH1F*)gFile->Get("tofNC"); h1->Fill( ctmNC-ctmT0 );}
	      
	      if(flagCDC&&pid_beam==1)
		{
		  double fdis=(vtxNC-gNC).Mag();
		  double beam_time= (vtxT0-vtxNC).Mag()/beta_b/Const;
		  double ntime=ctmNC-ctmT0- beam_time;
		  double beta_nc=fdis/ntime/Const;
		  //	      std::cout<<"fdis beamt ntime bn "<<fdis<<" "<<beam_time<<" "<<ntime<<" "<<beta_nc<<std::endl;
		  h1 = (TH1F*)gFile->Get("tofNC_vtx"); h1->Fill( ctmNC-ctmT0 );
		  h1 = (TH1F*)gFile->Get("obetaNC_vtx"); h1->Fill(1./beta_nc );
		  if(1./beta_nc>1.1)
		    {	      
		      double  mom_n=beta_nc*nMass/sqrt(1-beta_nc*beta_nc);
		      h1 = (TH1F*)gFile->Get("momNC_vtx"); h1->Fill( mom_n );
		    }
		}
	      
	      if(flagnK)
		{
		  double fdis=(vtxNC-gNC).Mag();
		  double beam_time= (vtxT0-vtxNC).Mag()/beta_b/Const;
		  double ntime=ctmNC-ctmT0- beam_time;
		  double beta_nc=fdis/ntime/Const;
		  //	      std::cout<<"fdis beamt ntime bn "<<fdis<<" "<<beam_time<<" "<<ntime<<" "<<beta_nc<<std::endl;
		  //		  h1 = (TH1F*)gFile->Get("tofNC_vtx"); h1->Fill( ctmNC-ctmT0 );
		  //		  h1 = (TH1F*)gFile->Get("obetaNC_vtx"); h1->Fill(1./beta_nc );
		  double ztmp2=(gNC-vtxNC).Z();
		  double xtmp2=ztmp2*Pnk.X()/Pnk.Z()+vtxNC.X();
		  std::cout<<"K x2 "<<xtmp2<<std::endl;
		  if( 1./beta_nc>1.1)
		    {
		      h1 = (TH1F*)gFile->Get("NChit_nK"); h1->Fill( xtmp2-gNC.X() );
		    }
		  
		  double  mom_n=beta_nc*nMass/sqrt(1-beta_nc*beta_nc);

		      //		      h1 = (TH1F*)gFile->Get("momNC_vtx"); h1->Fill( mom_n );
		}
	    
	      
	      if(flagnK0s){
		double fdis=(vtxNC-gNC).Mag();
		double beam_time= (vtxT0-vtxNC).Mag()/beta_b/Const;
		double ntime=ctmNC-ctmT0-beam_time;
		double beta_nc=fdis/ntime/Const;
		double ztmp2=(gNC-vtxNC).Z();
		double xtmp2=ztmp2*Pnk.X()/Pnk.Z()+vtxNC.X();		
		double ytmp2=ztmp2*Pnk.Y()/Pnk.Z()+vtxNC.Y();		
		std::cout<<"K0s x2 "<<xtmp2<<" "<<gNC.X()<<std::endl;
		std::cout<<"K0s y2 "<<ytmp2<<" "<<gNC.Y()<<std::endl;
		if( 1./beta_nc>1.1)
		  {
		    h1 = (TH1F*)gFile->Get("NChit_nK0s"); h1->Fill( xtmp2-gNC.X() );
		    h2 = (TH2F*)gFile->Get("NChitX_nK0s"); h2->Fill( gNC.X(),xtmp2 );
		    h2 = (TH2F*)gFile->Get("NChitY_nK0s"); h2->Fill( gNC.Y(),ytmp2 );
		    if(fabs(xtmp2-gNC.X())<10 )h1 = (TH1F*)gFile->Get("NCeff_nK0s2"); h1->Fill( gNC.X() );
		  }
		  
		double mom_n=beta_nc*nMass/sqrt(1-beta_nc*beta_nc);
		
	      }

	      if(flaglambda){
		double fdis=(vtxNC-gNC).Mag();
		double beam_time= (vtxT0-vtxNC).Mag()/beta_b/Const;
		double ntime=ctmNC-ctmT0-beam_time;
		double beta_nc=fdis/ntime/Const;
		h1 = (TH1F*)gFile->Get("tofNC_L"); h1->Fill( ctmNC-ctmT0 );
		h1 = (TH1F*)gFile->Get("obetaNC_L"); h1->Fill( 1./beta_nc );
		if(1./beta_nc>1.1)
		  {	      
		    double mom_n=beta_nc*nMass/sqrt(1-beta_nc*beta_nc);
		    h1 = (TH1F*)gFile->Get("momNC_L"); h1->Fill( mom_n );
	      }
		
	      }
	      
	      if(flagLp){
		double fdis=(vtxNC-gNC).Mag();
		double beam_time= (vtxT0-vtxNC).Mag()/beta_b/Const;
		double ntime=ctmNC-ctmT0-beam_time;
		double beta_nc=fdis/ntime/Const;
		h1 = (TH1F*)gFile->Get("tofNC_Lp"); h1->Fill( ctmNC-ctmT0 );
		h1 = (TH1F*)gFile->Get("obetaNC_Lp"); h1->Fill( 1./beta_nc );
		if(1./beta_nc>1.1)
		  {	      
		    double mom_n=beta_nc*nMass/sqrt(1-beta_nc*beta_nc);
		    h1 = (TH1F*)gFile->Get("momNC_Lp"); h1->Fill( mom_n );
		    TLorentzVector L_n;
		    TVector3 Pn;Pn=(gNC-vtxNC).Unit()*mom_n;
		    L_n.SetVectM(Pn,nMass);
		    TVector3 tP=(L3_L+L3_p-L3_target-L3_beam).P();
		    double tE=(L3_L+L3_p-L3_target-L3_beam).E();
		    //		h1 = (TH1F*)gFile->Get("tE_Lpn"); h1->Fill( Tp+Tn+TL-tQ );
		    //		h1 = (TH1F*)gFile->Get("tPz_tP_Lpn"); h1->Fill( Tp+Tn+TL-tQ );
		    //		h1 = (TH1F*)gFile->Get("tPz_Lpn"); h1->Fill( Tp+Tn+TL-tQ );
		    TVector3 boost1=(L3_target+L3_beam).BoostVector();
		    
		    double pzn=L_n.Pz();
		    L3_target.Boost(-1*boost1);
		    L3_beam.Boost(-1*boost1);
		    L3_p.Boost(-1*boost1);
		    L3_L.Boost(-1*boost1);
		    L_n.Boost(-1*boost1);
		    std::cout<<"missMass : "<<(L3_target+L3_beam-L3_p-L3_L-L_n).M()
			     <<" CMiniP "<<(L3_target+L3_beam).P()
			     <<" Pz Pzcmn " <<pzn<<" "<<L_n.Pz()<<std::endl;
		    std::cout<<" P " <<(L3_p+L3_L+L_n).P()<<std::endl;
		    double Tp=L3_p.E()-pMass;
		    double TL=L3_L.E()-lMass;
		    double Tn=L_n.E()-nMass;
		    double tQ=(L3_beam+L3_target).M()-lMass-pMass-nMass;
		    std::cout<<"S-Sf : "<<(L3_target+L3_beam).M()-(L3_p+L3_L+L_n).M()<<std::endl;
		    std::cout<<" tQ E " <<tQ<<" "<<(L3_beam+L3_target).E()<<std::endl;
		    std::cout<<" TL/q : Tp-Tn/s3Q " <<TL<<" "<<Tp<<" "<<Tn<<std::endl;
		    std::cout<<" iev" <<iev<<std::endl;
		    //		h1 = (TH1F*)gFile->Get("Q_Lpn"); h1->Fill( Tp+Tn+TL-tQ );
		    h1 = (TH1F*)gFile->Get("Mm_Lpn"); h1->Fill( (L3_beam+L3_target-L3_p-L3_L-L_n).M() );
		    h2 = (TH2F*)gFile->Get("dalitz_Lpn"); h2->Fill( (Tp-Tn)/sqrt(3)/tQ,TL/tQ );
		  }
	    
		
	      }
	      
	    }
	}
      



      tmptree->Fill();
      
      //#####################//
 

  }//iev
  
  
  gFile->Write();
  gFile->Close();
  
  gSystem->Exit(1);
  gROOT->GetApplication()->Run();
  //  theApp.Run();

  return 0;
}

