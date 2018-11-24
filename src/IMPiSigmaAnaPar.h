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

  const double chi2cut=6;

  const double covValSpmode[kin::npart][16] = {
    { 1.89064e-05, 0, 0, 0,
      0, 1.89622e-05, 0, 0,
      0, 0, 4.0112e-06, 0,
      0, 0, 0, 3.22341e-06 },
    { 0.000147362, 0, 0, 0,
      0, 0.000151568, 0, 0,
      0, 0, 0.000187105, 0,
      0, 0, 0, 0.000302456 },
    { 0.0226731, 0, 0, 0,
      0, 0.0227568, 0, 0,
      0, 0, 0.0152559, 0,
      0, 0, 0, 0.0109434 },
    { 0.0210347, 0, 0, 0,
      0, 0.0213206, 0, 0,
      0, 0, 0.0147438, 0,
      0, 0, 0, 0.00531022 },
    { 0.0206777, 0, 0, 0,
      0, 0.0209008, 0, 0,
      0, 0, 0.0146246, 0,
      0, 0, 0, 0.00461358 },
    { 2.96805e-05, 0, 0, 0,
      0, 2.972e-05, 0, 0,
      0, 0, 5.40635e-05, 0,
      0, 0, 0, 3.58992e-05 },
  };
// 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
  const double covValSmmode[kin::npart][16] = {
    { 1.81671e-05, 0, 0, 0,
      0, 1.83872e-05, 0, 0,
      0, 0, 4.03316e-06, 0,
      0, 0, 0, 3.24089e-06 },
    { 4.25399e-05, 0, 0, 0,
      0, 4.18085e-05, 0, 0,
      0, 0, 9.90933e-05, 0,
      0, 0, 0, 0.000107376 },
    { 0.000695792, 0, 0, 0,
      0, 0.000800459, 0, 0,
      0, 0, 0.000225677, 0,
      0, 0, 0, 0.000110685 },
    { 0.00110024, 0, 0, 0,
      0, 0.00105452, 0, 0,
      0, 0, 0.000646375, 0,
      0, 0, 0, 0.000296747 },
    { 0.000273254, 0, 0, 0,
      0, 0.000278707, 0, 0,
      0, 0, 4.94206e-05, 0,
      0, 0, 0, 2.54488e-05 },
    { 2.25895e-05, 0, 0, 0,
      0, 2.26291e-05, 0, 0,
      0, 0, 4.86585e-05, 0,
      0, 0, 0, 3.0968e-05 },
  };
}



namespace gen{
  const int reactionID_Spmode = 1725; 
  const int reactionID_Smmode = 1525; 
}
