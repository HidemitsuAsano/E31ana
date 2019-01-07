//IMAnaParPiSigma
//H. Asano
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
  const int cds_ngoodtrack = 2;
  const int cdhmulti = 3;
  const double tdc_cdh_max = 40; // ns 
  const double cds_chi2_max = 30;
  const bool useclosestpi = false;
}

namespace anacuts{
  const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
  const double dE_MIN = 1.8; 

  const double pipi_MIN = 0.485;
  const double pipi_MAX = 0.510;
  const double ppi_MIN = 1.1075;
  const double ppi_MAX = 1.1225;

  const double neutron_MIN = 0.85;
  const double neutron_MAX = 1.03;

  const double Sigmap_MIN = 1.17;
  const double Sigmap_MAX = 1.21;
  const double Sigmam_MIN = 1.18;
  const double Sigmam_MAX = 1.22;
}



namespace kin{
  const int npart=6;
  const int kmbeam=0;
  const int pim_g1=1;//pi- 1st generation
  const int Sp=2;
  const int nmiss=3;
  const int ncds=4;
  const int pip_g2=5;//pi+ 2nd generation
  const int pip_g1=1;//pi- 1st generation
  const int Sm=2;
  const int pim_g2=5;//pi- 2nd generation
  
  const int maxitr=50;
  const double maxdchi2=5e-5;
  const double maxsumconst=1e-4;

  const double chi2cut=6;

  const double covValSpmode[kin::npart][16] = {
    { 1.81201e-05, 0, 0, 0,
      0, 1.83048e-05, 0, 0,
      0, 0, 4.01283e-06, 0,
      0, 0, 0, 3.21263e-06 },
    { 4.6e-05, 0, 0, 0,
      0, 4.41782e-05, 0, 0,
      0, 0, 0.000106602, 0,
      0, 0, 0, 0.00011834 },
    { 0.000448925, 0, 0, 0,
      0, 0.000471427, 0, 0,
      0, 0, 0.000188364, 0,
      0, 0, 0, 9.50713e-05 },
    { 0.000829595, 0, 0, 0,
      0, 0.000672813, 0, 0,
      0, 0, 0.000509478, 0,
      0, 0, 0, 0.000246371 },
    { 0.000231783, 0, 0, 0,
      0, 0.000234654, 0, 0,
      0, 0, 4.53459e-05, 0,
      0, 0, 0, 4.81177e-06 },
    { 1.66653e-05, 0, 0, 0,
      0, 1.68917e-05, 0, 0,
      0, 0, 3.96755e-05, 0,
      0, 0, 0, 2.14282e-05 },
  };
// 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
  const double covValSmmode[kin::npart][16] = {
    { 1.83012e-05, 0, 0, 0,
      0, 1.84406e-05, 0, 0,
      0, 0, 4.02792e-06, 0,
      0, 0, 0, 3.23206e-06 },
    { 4.15989e-05, 0, 0, 0,
      0, 4.00911e-05, 0, 0,
      0, 0, 9.65422e-05, 0,
      0, 0, 0, 0.000103395 },
    { 0.000505624, 0, 0, 0,
      0, 0.000572862, 0, 0,
      0, 0, 0.000219836, 0,
      0, 0, 0, 0.00012269 },
    { 0.000789841, 0, 0, 0,
      0, 0.000853735, 0, 0,
      0, 0, 0.000553476, 0,
      0, 0, 0, 0.000240124 },
    { 0.00025318, 0, 0, 0,
      0, 0.000252536, 0, 0,
      0, 0, 5.01769e-05, 0,
      0, 0, 0, 4.67395e-06 },
    { 2.1251e-05, 0, 0, 0,
      0, 2.14052e-05, 0, 0,
      0, 0, 4.74657e-05, 0,
      0, 0, 0, 2.9605e-05 },
    };
}

namespace gen{
  const int reactionID_Spmode = 1725; 
  const int reactionID_Smmode = 1525; 
}
