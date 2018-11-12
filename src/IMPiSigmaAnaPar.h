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
  const double tdc_cdh_max = 9999; // ns 
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

  const double covVal1[7][16] = {
    { 1.90578e-05, 0, 0, 0,
      0, 1.73984e-05, 0, 0,
      0, 0, 4.67911e-06, 0,
      0, 0, 0, 3.88385e-06 },
    { 1.23591e-05, 0, 0, 0,
      0, 1.08675e-05, 0, 0,
      0, 0, 2.21353e-05, 0,
      0, 0, 0, 7.2569e-06 },
    { 0.000304799, 0, 0, 0,
      0, 0.000331023, 0, 0,
      0, 0, 8.24242e-05, 0,
      0, 0, 0, 4.10355e-05 },
    { 0.000105953, 0, 0, 0,
      0, 7.71763e-05, 0, 0,
      0, 0, 0.000120953, 0,
      0, 0, 0, 3.90285e-05 },
    { 0.000468473, 0, 0, 0,
      0, 0.000831477, 0, 0,
      0, 0, 0.000354281, 0,
      0, 0, 0, 0.000117394 },
    { 0.000139189, 0, 0, 0,
      0, 0.000220018, 0, 0,
      0, 0, 2.41655e-05, 0,
      0, 0, 0, 5.73053e-06 },
    { 1.11052e-05, 0, 0, 0,
      0, 1.15129e-05, 0, 0,
      0, 0, 1.7214e-05, 0,
      0, 0, 0, 8.1515e-06 }
  };
// 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
  const double covVal2[7][16] = {
    { 1.90578e-05, 0, 0, 0,
      0, 1.73984e-05, 0, 0,
      0, 0, 4.67911e-06, 0,
      0, 0, 0, 3.88385e-06 },
    { 1.23591e-05, 0, 0, 0,
      0, 1.08675e-05, 0, 0,
      0, 0, 2.21353e-05, 0,
      0, 0, 0, 7.2569e-06 },
    { 0.000304799, 0, 0, 0,
      0, 0.000331023, 0, 0,
      0, 0, 8.24242e-05, 0,
      0, 0, 0, 4.10355e-05 },
    { 0.000105953, 0, 0, 0,
      0, 7.71763e-05, 0, 0,
      0, 0, 0.000120953, 0,
      0, 0, 0, 3.90285e-05 },
    { 0.000468473, 0, 0, 0,
      0, 0.000831477, 0, 0,
      0, 0, 0.000354281, 0,
      0, 0, 0, 0.000117394 },
    { 0.000139189, 0, 0, 0,
      0, 0.000220018, 0, 0,
      0, 0, 2.41655e-05, 0,
      0, 0, 0, 5.73053e-06 },
    { 1.11052e-05, 0, 0, 0,
      0, 1.15129e-05, 0, 0,
      0, 0, 1.7214e-05, 0,
    0, 0, 0, 8.1515e-06 }
  };
}



namespace gen{
  const int reactionID_Spmode = 1725; 
  const int reactionID_Smmode = 1525; 
}
