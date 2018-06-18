#ifndef CS_PARAM
#define CS_PARAM 1

#include <iostream>

const double barn = 1.0e24;
const double micro = 1.0e6;

const double num_of_beam = 6.835e9;
const double err_of_beam = 8.5e7;

const double num_of_target=10*0.169*6.02e23/2.01;

const double eff_trigN=0.946*0.998;
const double err_trigN=sqrt(0.021*0.021+0.002*0.002);
const double eff_trigC=0.946*0.955;
const double err_trigC=sqrt(0.021*0.021+0.022*0.022);

const double eff_DAQ=0.818;
const double err_DAQ=0.01;

const double eff_CDC=0.977;
const double err_CDC=0.004;

const double omega_NC = 3.2*1.5/(14.9476*14.9476);

const double eff_NC=0.292;
const double err_NC=0.018;

const double eff_FDC1=0.966714;
const double err_FDC1=0.007;

/***************************/
/* old parameter 2017/2/22 */
/***************************/
/* const double barn = 1.0e24; */
/* const double micro = 1.0e6; */

/* const double num_of_beam = 6.95e9; */
/* const double num_of_target=10*0.169*6.02e23/2.01; */

/* const double eff_DAQ=0.818; */
/* const double eff_CDC=0.977; */

/* const double omega_NC = 3.2*1.5/(14.9476*14.9476); */
/* const double eff_NC=0.289; */

/* const double eff_FDC1=0.966714; */

#endif
