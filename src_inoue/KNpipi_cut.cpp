#include "KNpipi_cut.h"

const double K0_peak=0.497656;
const double K0_sigma=0.005958;
const double Sm_peak=1.197868;
const double Sm_sigma=0.006502;
const double Sp_peak=1.189276;
const double Sp_sigma=0.007325;

// const double K0_peak  = 0.497689;
// const double K0_sigma = 0.00536856;
// const double Sm_peak  = 1.19742;
// const double Sm_sigma = 0.0048656;
// const double Sp_peak  = 1.18901;
// const double Sp_sigma = 0.00380568;

const double range  = 3.0;
const double K0_min = K0_peak-range*K0_sigma;
const double K0_max = K0_peak+range*K0_sigma;
const double Sm_min = Sm_peak-range*Sm_sigma;
const double Sm_max = Sm_peak+range*Sm_sigma;
const double Sp_min = Sp_peak-range*Sp_sigma;
const double Sp_max = Sp_peak+range*Sp_sigma;

bool isK0(const double im){ return K0_min<im && im<K0_max; }
bool isIM_Sm(const double im){ return Sm_min<im && im<Sm_max; }
bool isIM_Sp(const double im){ return Sp_min<im && im<Sp_max; }
bool isSignal(const double pipi, const double npim, const double npip){ return !isK0(pipi) && !isIM_Sm(npim) && !isIM_Sp(npip); }

void printParam(TMacro *log){
  log-> AddLine(Form("K0_cut  %lf  %lf", K0_min, K0_max));
  log-> AddLine(Form("Sm_cut  %lf  %lf", Sm_min, Sm_max));
  log-> AddLine(Form("Sp_cut  %lf  %lf", Sp_min, Sp_max));
}
