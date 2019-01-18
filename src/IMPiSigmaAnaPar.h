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
  const double tdc_cdh_max = 9999; // ns v26
  const double cds_chi2_max = 30;
  const bool useclosestpi = false;
}

namespace anacuts {
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
    { 1.78541e-05, 6.94256e-310, 0, 0,
      1.03615e-316, 1.76636e-05, 1.64242e-311, 0,
      0, 0, 4.01777e-06, 0,
      4.94066e-324, 4.94066e-324, 1.03615e-316, 3.21566e-06 },
    { 2.19618e-05, 7.8726e-312, 0, 0,
      0, 2.20191e-05, 0, 4.94066e-324,
      4.94066e-324, 1.03615e-316, 6.78153e-05, 0,
      7.8726e-312, 0, 0, 5.21473e-05 },
    { 0.000254424, 0, 4.94066e-324, 4.94066e-324,
      1.03615e-316, 0.000215581, 0, 7.8726e-312,
      0, 0, 0.000112185, 0,
      0, 4.94066e-324, 4.94066e-324, 5.51874e-05 },
    { 0.000394157, 0, 7.8726e-312, 0,
      0, 0.000348927, 0, 0,
      2.122e-314, 0, 0.00028936, 1.03615e-316,
      1.03615e-316, 1.64242e-311, 0, 0.000143852 },
    { 0.000138964, 0, 0, 4.94066e-324,
      4.94066e-324, 0.000137699, 1.03615e-316, 0,
      1.64242e-311, 0, 3.50564e-05, 0,
      0, 0, 4.94066e-324, 3.17828e-06 },
    { 1.27274e-05, 1.03615e-316, 0, 7.8726e-312,
      0, 1.28568e-05, 0, 0,
      0, 4.94066e-324, 3.53656e-05, 1.03615e-316,
      1.03615e-316, 0, 7.8726e-312, 1.52959e-05 },
  };
  // 2) TLorentzVector L3_beam, L_pip, (L_n+L_pim), L_p, L_nmiss, L_n, L_pim = for pi+ Sigma-
  const double covValSmmode[kin::npart][16] = {
    { 1.79699e-05, 6.93331e-310, 0, 0,
      1.03615e-316, 1.81322e-05, 1.64242e-311, 0,
      0, 0, 3.95469e-06, 0,
      4.94066e-324, 4.94066e-324, 1.03615e-316, 3.17793e-06
    },
    { 1.91565e-05, 7.8726e-312, 0, 0,
      0, 1.8974e-05, 0, 4.94066e-324,
      4.94066e-324, 1.03615e-316, 6.37481e-05, 0,
      7.8726e-312, 0, 0, 3.88994e-05
    },
    { 0.000267432, 0, 4.94066e-324, 4.94066e-324,
      1.03615e-316, 0.000268238, 0, 7.8726e-312,
      0, 0, 0.000130654, 0,
      0, 4.94066e-324, 4.94066e-324, 7.2603e-05
    },
    { 0.000427818, 0, 7.8726e-312, 0,
      0, 0.000412773, 0, 0,
      2.122e-314, 0, 0.000281576, 1.03615e-316,
      1.03615e-316, 1.64242e-311, 0, 0.000136006
    },
    { 0.000143547, 0, 0, 4.94066e-324,
      4.94066e-324, 0.000144877, 1.03615e-316, 0,
      1.64242e-311, 0, 3.86582e-05, 0,
      0, 0, 4.94066e-324, 3.39109e-06
    },
    { 1.59753e-05, 1.03615e-316, 0, 7.8726e-312,
      0, 1.61775e-05, 0, 0,
      0, 4.94066e-324, 4.22882e-05, 1.03615e-316,
      1.03615e-316, 0, 7.8726e-312, 2.06409e-05
    },
  };
}

namespace gen {
  const int reactionID_Spmode = 1725;
  const int reactionID_Smmode = 1525;
}

#endif 
