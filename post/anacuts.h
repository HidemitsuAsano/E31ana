namespace anacuts {
  const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
  const double beta_MIN = 0.000000001;
  const double dE_MIN = 2.0;

  //const double pipi_MIN = 0.480;
  //const double pipi_MAX = 0.518;
  const double pipi_MIN = 0.475;
  const double pipi_MAX = 0.522;
  //const double neutron_MIN = 0.85;
  //const double neutron_MAX = 1.03;
  const double neutron_MIN = 0.939161-2.5*0.0308562;
  const double neutron_MAX = 0.939161+2.5*0.0308562;
  //const double neutron_MIN = 0.939161-1.5*0.0308562;
  //const double neutron_MAX = 0.939161+1.5*0.0308562;
  
  //const double Sigmap_MIN = 1.17;
  //const double Sigmap_MAX = 1.21;
  //const double Sigmam_MIN = 1.18;
  //const double Sigmam_MAX = 1.22;
  const double Sigmap_center = 1.18911;
  const double Sigmap_MIN = 1.18911-3.0*0.00540844;   
  const double Sigmap_MAX = 1.18911+3.0*0.00540844;
  const double Sigmam_center = 1.19723;
  const double Sigmam_MIN = 1.19723-3.0*0.00601265;
  const double Sigmam_MAX = 1.19723+3.0*0.00601265;  
  
  const double Sigmap_MIN_wide = 1.18911-5.0*0.00540844;   
  const double Sigmap_MAX_wide = 1.18911+5.0*0.00540844;   
  const double Sigmam_MIN_wide = 1.19723-5.0*0.00601265;
  const double Sigmam_MAX_wide = 1.19723+5.0*0.00601265;  
  
  //side band selection
  const double Sigmap_sigma = 0.00540844;
  const double Sigmap_sidelow_MIN = 1.18911-6.0*0.00540844;   
  const double Sigmap_sidelow_MAX = 1.18911-3.0*0.00540844;   
  //const double Sigmap_sidelow_MIN = 1.18911-11.0*0.00540844;   
  //const double Sigmap_sidelow_MAX = 1.18911-8.0*0.00540844;   
  const double Sigmap_sidelow_center = (Sigmap_sidelow_MAX + Sigmap_sidelow_MIN)/2.0;
  const double Sigmap_sidehigh_MIN = 1.18911+3.0*0.00540844;   
  const double Sigmap_sidehigh_MAX = 1.18911+6.0*0.00540844;   
  //const double Sigmap_sidehigh_MIN = 1.18911+8.0*0.00540844;   
  //const double Sigmap_sidehigh_MAX = 1.18911+11.0*0.00540844;   
  const double Sigmap_sidehigh_center = (Sigmap_sidehigh_MAX + Sigmap_sidehigh_MIN)/2.0;
  
  const double Sigmam_sigma = 0.00601265;
  const double Sigmam_sidelow_MIN = 1.19723-6.0*0.00601265;
  const double Sigmam_sidelow_MAX = 1.19723-3.0*0.00601265;  
  //const double Sigmam_sidelow_MIN = 1.19723-11.0*0.00601265;
  //const double Sigmam_sidelow_MAX = 1.19723-8.0*0.00601265;  
  const double Sigmam_sidelow_center = (Sigmam_sidelow_MAX + Sigmam_sidelow_MIN)/2.0;
  const double Sigmam_sidehigh_MIN = 1.19723+3.0*0.00601265;
  const double Sigmam_sidehigh_MAX = 1.19723+6.0*0.00601265;  
  //const double Sigmam_sidehigh_MIN = 1.19723+8.0*0.00601265;
  //const double Sigmam_sidehigh_MAX = 1.19723+11.0*0.00601265;  
  const double Sigmam_sidehigh_center = (Sigmam_sidehigh_MAX + Sigmam_sidehigh_MIN)/2.0;

  /*
  const double Sigmap_MIN_wide = 1.18911-4.0*0.00540844;   
  const double Sigmap_MAX_wide = 1.18911+4.0*0.00540844;   
  const double Sigmam_MIN_wide = 1.19723-4.0*0.00601265;
  const double Sigmam_MAX_wide = 1.19723+4.0*0.00601265;  
  
  //side band selection
  const double Sigmap_sidelow_MIN = 1.18911-6.0*0.00540844;   
  const double Sigmap_sidelow_MAX = 1.18911-4.0*0.00540844;   
  const double Sigmap_sidehigh_MIN = 1.18911+4.0*0.00540844;   
  const double Sigmap_sidehigh_MAX = 1.18911+6.0*0.00540844;   

  const double Sigmam_sidelow_MIN = 1.19723-6.0*0.00601265;
  const double Sigmam_sidelow_MAX = 1.19723-4.0*0.00601265;  
  const double Sigmam_sidehigh_MIN = 1.19723+4.0*0.00601265;
  const double Sigmam_sidehigh_MAX = 1.19723+6.0*0.00601265;  
  */

  const double qvalcut = 0.35; 
}

