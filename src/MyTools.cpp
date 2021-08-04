//v2//
#include "MyTools.h"

using namespace std;



int MyTools::PIDcorr_wide(double mom,double mass2)
{
  double pi_sigma = sqrt(4.*piMass*piMass*mom*mom*PID_Param[0]+4.*piMass*piMass*piMass*piMass*PID_Param[1]*(1+piMass*piMass/(mom*mom))+4.*mom*mom*(piMass*piMass+mom*mom)*PID_Param[2]);
  double k_sigma = sqrt(4.*kpMass*kpMass*mom*mom*PID_Param[0]+4.*kpMass*kpMass*kpMass*kpMass*PID_Param[1]*(1+kpMass*kpMass/(mom*mom))+4.*mom*mom*(kpMass*kpMass+mom*mom)*PID_Param[2]);
  double p_sigma = sqrt(4.*pMass*pMass*mom*mom*PID_Param[0]+4.*pMass*pMass*pMass*pMass*PID_Param[1]*(1+pMass*pMass/(mom*mom))+4.*mom*mom*(pMass*pMass+mom*mom)*PID_Param[2]);

  bool pim_flag=false, pip_flag=false, km_flag=false, p_flag=false;
  if( mom>0 ){
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pip_flag=true;
    if( pMass*pMass-2.5*p_sigma<mass2 && mass2<pMass*pMass+2.5*p_sigma && mom>0.1 ) p_flag=true;

    if( pip_flag && p_flag ){
      if( mass2<Ppi_mid_mass2 ) return CDS_PiPlus;
      else return CDS_Proton;
    }
    else if( pip_flag ) return CDS_PiPlus;
    else if( p_flag ) return CDS_Proton;
  }
  else{
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pim_flag=true;
    if( kpMass*kpMass-2.5*k_sigma<mass2 && mass2<kpMass*kpMass+2.5*k_sigma && mom<-0.03 ) km_flag=true;

    if( pim_flag ) return CDS_PiMinus;
    else if( km_flag ) return CDS_Kaon;
  }
  return CDS_Other;
}


