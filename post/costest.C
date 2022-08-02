#include "../src/GlobalVariables.h"

void costest()
{
  TLorentzVector mom_target;
  TVector3 P_target(0,0,0);
  mom_target.SetVectM(P_target,pMass);
  std::cout << mom_target.M() << std::endl;

  TVector3 P_beam(0,0,1.018);
  TLorentzVector mom_beam;
  mom_beam.SetVectM(P_beam,kpMass);

  TVector3 P_Sigma(0,0.8,0.0);
  TLorentzVector mom_S;
  mom_S.SetVectM(P_Sigma,spMass);

  TLorentzVector mom_misspi;
  mom_misspi = mom_target+mom_beam-mom_S;
  TVector3 boost = (mom_target+mom_beam).BoostVector();
  mom_misspi.Boost(-boost);
  std::cout << mom_misspi.CosTheta() << std::endl;


}
