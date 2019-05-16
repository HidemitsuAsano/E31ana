//IMAnaParPiSigma
//H. Asano

#ifndef IMPISIGMAANAPAR_HH
#define IMPISIGMAANAPAR_HH 1

namespace blcuts {

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

namespace cdscuts {
  const int cds_ngoodtrack = 2;
  const int cdhmulti = 3;
  //const double tdc_cdh_max = 40; // ns v25
  const double tdc_cdh_max = 50; // ns v32-v40
  //const double tdc_cdh_max = 9999; // ns v32-v40
  const double tdc_simoffset = -18;
  const double cds_chi2_max = 30;
  const bool useclosestpi = false;
  const double chargevetoangle = 15.0;
}

namespace cdscuts_lpim{
  const int cds_ngoodtrack = 3;
  const int cdhmulti = 3;
  //const double tdc_cdh_max = 40; // ns v25
  const double tdc_cdh_max = 50; // ns v32-v40
  //const double tdc_cdh_max = 9999; // ns v32-v40
  const double tdc_simoffset = -18;
  const double cds_chi2_max = 30;
  const bool useclosestpi = false;
  const double chargevetoangle = 15.0;
}

namespace anacuts {
  const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
  const double dE_MIN = 2.0;

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

namespace anacuts_lpim {
  const double ppi_MIN = 1.1075;
  const double ppi_MAX = 1.1225;

  const double proton_MIN = 0.83;
  const double proton_MAX = 1.04;
}


namespace kin {
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
    { 1.79228e-05, 0, 0, 0,
      0, 1.82956e-05, 0, 0,
      0, 0, 4.02874e-06, 0,
      0, 0, 0, 3.23335e-06 },
    { 2.05071e-05, 0, 0, 0,
      0, 1.9884e-05, 0, 0,
      0, 0, 7.15267e-05, 0,
      0, 0, 0, 4.39928e-05 },
    { 0.000165845, 0, 0, 0,
      0, 0.000166991, 0, 0,
      0, 0, 0.000115236, 0,
      0, 0, 0, 4.90095e-05 },
    { 0.000370377, 0, 0, 0,
      0, 0.000368532, 0, 0,
      0, 0, 0.000336299, 0,
      0, 0, 0, 0.000174533 },
    { 0.000141756, 0, 0, 0,
      0, 0.000139519, 0, 0,
      0, 0, 6.62114e-05, 0,
      0, 0, 0, 2.7891e-06 },
    { 1.23728e-05, 0, 0, 0,
      0, 1.27203e-05, 0, 0,
      0, 0, 3.29329e-05, 0,
      0, 0, 0, 1.50133e-05 },
  };
  // 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
  const double covValSmmode[kin::npart][16] = {
    { 1.81413e-05, 0, 0, 0,
      0, 1.81383e-05, 0, 0,
      0, 0, 4.00711e-06, 0,
      0, 0, 0, 3.21603e-06 },
    { 1.8442e-05, 0, 0, 0,
      0, 1.84139e-05, 0, 0,
      0, 0, 6.35269e-05, 0,
      0, 0, 0, 4.40656e-05 },
    { 0.0001729, 0, 0, 0,
      0, 0.000171629, 0, 0,
      0, 0, 0.000120585, 0,
      0, 0, 0, 5.8488e-05 },
    { 0.000397142, 0, 0, 0,
      0, 0.000367854, 0, 0,
      0, 0, 0.000327771, 0,
      0, 0, 0, 0.000168868 },
    { 0.000143633, 0, 0, 0,
      0, 0.000140488, 0, 0,
      0, 0, 6.66998e-05, 0,
      0, 0, 0, 2.84004e-06 },
    { 1.54122e-05, 0, 0, 0,
      0, 1.5112e-05, 0, 0,
      0, 0, 3.76225e-05, 0,
      0, 0, 0, 1.92646e-05 },
  };
}

namespace kinpLpim {
  const int npart=6;
  const int kmbeam=0;
  const int pim_g1=1;//pi- 1st generation
  const int L=2;
  const int pmiss=3;
  const int pcds=4;
  const int pim_g2=5;//pi- 2nd generation

  const int maxitr=50;
  const double maxdchi2=5e-5;
  const double maxsumconst=1e-4;

  const double chi2cut=6;

  const double covValpLpim[kin::npart][16] = {
    { 1.79228e-05, 0, 0, 0,
      0, 1.82956e-05, 0, 0,
      0, 0, 4.02874e-06, 0,
      0, 0, 0, 3.23335e-06 },
    { 2.05071e-05, 0, 0, 0,
      0, 1.9884e-05, 0, 0,
      0, 0, 7.15267e-05, 0,
      0, 0, 0, 4.39928e-05 },
    { 0.000165845, 0, 0, 0,
      0, 0.000166991, 0, 0,
      0, 0, 0.000115236, 0,
      0, 0, 0, 4.90095e-05 },
    { 0.000370377, 0, 0, 0,
      0, 0.000368532, 0, 0,
      0, 0, 0.000336299, 0,
      0, 0, 0, 0.000174533 },
    { 0.000141756, 0, 0, 0,
      0, 0.000139519, 0, 0,
      0, 0, 6.62114e-05, 0,
      0, 0, 0, 2.7891e-06 },
    { 1.23728e-05, 0, 0, 0,
      0, 1.27203e-05, 0, 0,
      0, 0, 3.29329e-05, 0,
      0, 0, 0, 1.50133e-05 },
  };
}

namespace gen {
  const int reactionID_Spmode = 1725;
  const int reactionID_Smmode = 1525;
  const int reactionID_pLpim  = 1600;
}

#endif 
