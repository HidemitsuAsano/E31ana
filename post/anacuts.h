namespace anacuts {
  const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
  //const double beta_MAX = 0.628786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
  const double beta_MIN = 0.000000001;
  //const double dE_MIN = 2.0/0.88;
  const double dE_MIN = 2.0;

  const double K0_NSigmacut=3.0;
  const double K0_NSigmacut_narrow=2.0;
  const double K0_center = 0.497855;
  const double K0_sigma = 0.0074845;
  const double pipi_MIN = K0_center - K0_NSigmacut*K0_sigma;
  const double pipi_MAX = K0_center + K0_NSigmacut*K0_sigma;
  const double pipi_MIN_narrow = K0_center - K0_NSigmacut_narrow*K0_sigma;
  const double pipi_MAX_narrow = K0_center + K0_NSigmacut_narrow*K0_sigma;
  
  //const double pipi_MIN = 0.475;
  //const double pipi_MAX = 0.522;
  //const double pipi_MIN = 0.487;// (almost 1.5 sigma)
  //const double pipi_MAX = 0.510;// (almost 1.5 sigma)
  
  const double neutron_NSigmacut = 2.0;//+/- N sigma
  const double neutron_NSigmacut_wide = 3.0;//+/- N sigma
  const double neutron_center = 0.939161;
  const double neutron_sigma = 0.0308562;
  const double neutron_MIN = neutron_center-neutron_NSigmacut*neutron_sigma;
  const double neutron_MAX = neutron_center+neutron_NSigmacut*neutron_sigma;
  const double neutron_MIN_wide = neutron_center-neutron_NSigmacut_wide*neutron_sigma; 
  const double neutron_MAX_wide = neutron_center+neutron_NSigmacut_wide*neutron_sigma; 

  //neutron cut of CDS, this is necesarry due to the range of CDH3 trigger
  const double nmomcut = 0.14;
  //const double nmomcut = 0.20;

  const double Sigma_NSigmacut = 2.0;
  const double Sigmap_center = 1.18911;
  const double Sigmap_sigma = 0.00540844;
  const double Sigmap_MIN = Sigmap_center-Sigma_NSigmacut*Sigmap_sigma;   
  const double Sigmap_MAX = Sigmap_center+Sigma_NSigmacut*Sigmap_sigma;
  const double Sigmam_center = 1.19723;
  const double Sigmam_sigma = 0.00601265;
  const double Sigmam_MIN = Sigmam_center-Sigma_NSigmacut*Sigmam_sigma;
  const double Sigmam_MAX = Sigmam_center+Sigma_NSigmacut*Sigmam_sigma;  
  
  const double Sigma_wideNSigmacut = 3.0;
  const double Sigmap_MIN_wide = Sigmap_center-Sigma_wideNSigmacut*Sigmap_sigma;   
  const double Sigmap_MAX_wide = Sigmap_center+Sigma_wideNSigmacut*Sigmap_sigma;   
  const double Sigmam_MIN_wide = Sigmam_center-Sigma_wideNSigmacut*Sigmam_sigma;
  const double Sigmam_MAX_wide = Sigmam_center+Sigma_wideNSigmacut*Sigmam_sigma;  
  
  const double SigmaPMomCut = 0.04;
  const double SigmaMMomCut = 0.04;

  //side band selection
  //const double Sigmap_sidelow_MIN = Sigmap_center-6.0*Sigmap_sigma;   
  //const double Sigmap_sidelow_MAX = Sigmap_center-3.0*Sigmap_sigma;   
  const double Sigmap_sidelow_MIN = Sigmap_center-4.0*Sigmap_sigma;   
  const double Sigmap_sidelow_MAX = Sigmap_center-2.0*Sigmap_sigma;   
  const double Sigmap_sidelow_center = (Sigmap_sidelow_MAX + Sigmap_sidelow_MIN)/2.0;
  const double Sigmap_sidehigh_MIN = Sigmap_center+2.0*Sigmap_sigma;   
  const double Sigmap_sidehigh_MAX = Sigmap_center+4.0*Sigmap_sigma;   
  const double Sigmap_sidehigh_center = (Sigmap_sidehigh_MAX + Sigmap_sidehigh_MIN)/2.0;
  
  const double Sigmam_sidelow_MIN = Sigmam_center-4.0*Sigmam_sigma;
  const double Sigmam_sidelow_MAX = Sigmam_center-2.0*Sigmam_sigma;  
  const double Sigmam_sidelow_center = (Sigmam_sidelow_MAX + Sigmam_sidelow_MIN)/2.0;
  const double Sigmam_sidehigh_MIN = Sigmam_center+2.0*Sigmam_sigma;
  const double Sigmam_sidehigh_MAX = Sigmam_center+4.0*Sigmam_sigma;  
  const double Sigmam_sidehigh_center = (Sigmam_sidehigh_MAX + Sigmam_sidehigh_MIN)/2.0;


  const double qvalcut = 0.35; 

  //isolation cut parameters
  const double Isonpip_shift = -0.05;//cm
  const double Isonpip_phicut = 0.40;//radian
  const double Isonpip_zcut = 25.0;//cm
  const double Isonpip_phicutwide = 0.45;//radian
  const double Isonpip_zcutwide = 28.0;//cm

  const double Isonpim_shift = 0.05;//cm
  const double Isonpim_phicut = 0.60;//radian
  const double Isonpim_zcut = 25.0;//cm
  const double Isonpim_phicutwide = 0.65;//radian
  const double Isonpim_zcutwide = 28.0;//cm
  
  const double CDHwidthphi = 0.12;//radian

  
//  const double Lambda_MIN = 1.1075;//not tuned
//  const double Lambda_MAX = 1.1225;//not tuned
  const double Lambda_MIN = 1.1067;//3 sigma
  const double Lambda_MAX = 1.1250;//3 sigma
  const double Lambda_MIN_narrow = 1.1113;//1.5 sigma
  const double Lambda_MAX_narrow = 1.1204;//1.5 sigma
  const double Proton_MIN = 0.8303;//wide 
  const double Proton_MAX = 1.0279;//wide
  const double Proton_MIN_narrow = 0.8825;
  const double Proton_MAX_narrow = 0.9794;
  //const double Proton_MIN = 0.8825;
  //const double Proton_MAX = 0.9794;

}

