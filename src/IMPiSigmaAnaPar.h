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
   { 1.77541e-05, 0, 0, 0,
     0, 1.83282e-05, 0, 0,
     0, 0, 4.04348e-06, 0,
     0, 0, 0, 3.24158e-06 },
   { 4.17938e-05, 0, 0, 0,
     0, 4.08315e-05, 0, 0,
     0, 0, 9.80125e-05, 0,
     0, 0, 0, 0.000103327 },
   { 0.000332859, 0, 0, 0,
     0, 0.000326395, 0, 0,
     0, 0, 0.000122881, 0,
     0, 0, 0, 5.94971e-05 },
   { 0.000473442, 0, 0, 0,
     0, 0.000478175, 0, 0,
     0, 0, 0.000321085, 0,
     0, 0, 0, 0.000150869 },
   { 0.000187591, 0, 0, 0,
     0, 0.000181608, 0, 0,
     0, 0, 3.57942e-05, 0,
     0, 0, 0, 3.36171e-06 },
   { 1.31063e-05, 0, 0, 0,
     0, 1.26118e-05, 0, 0,
     0, 0, 3.50466e-05, 0,
     0, 0, 0, 1.51075e-05 },
   };
// 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
  const double covValSmmode[kin::npart][16] = {
    { 1.81631e-05, 0, 0, 0,
      0, 1.80253e-05, 0, 0,
      0, 0, 4.06498e-06, 0,
      0, 0, 0, 3.26103e-06 },
    { 3.81424e-05, 0, 0, 0,
      0, 3.66517e-05, 0, 0,
      0, 0, 8.88447e-05, 0,
      0, 0, 0, 8.88273e-05 },
    { 0.000356573, 0, 0, 0,
      0, 0.000378745, 0, 0,
      0, 0, 0.00014202, 0,
      0, 0, 0, 7.45554e-05 },
    { 0.000517849, 0, 0, 0,
      0, 0.000520009, 0, 0,
      0, 0, 0.000332048, 0,
      0, 0, 0, 0.000147404 },
    { 0.000192652, 0, 0, 0,
      0, 0.000184853, 0, 0,
      0, 0, 4.00797e-05, 0,
      0, 0, 0, 3.52531e-06 },
    { 1.60146e-05, 0, 0, 0,
      0, 1.61009e-05, 0, 0,
      0, 0, 4.13413e-05, 0,
      0, 0, 0, 2.07916e-05 },
    };
}

namespace gen{
  const int reactionID_Spmode = 1725; 
  const int reactionID_Smmode = 1525; 
}
