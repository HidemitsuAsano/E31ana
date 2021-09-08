void rootSE31()
{

  const double kpMass = 0.4936;
  const double dMass = 1.87561;


  const TVector3 Pp_beam(0,0,1.0);
  TLorentzVector LVec_beam;
  LVec_beam.SetVectM(Pp_beam, kpMass);
  
  const TVector3 Pp_target(0,0,0.0);
  TLorentzVector LVec_target;
  LVec_target.SetVectM(Pp_target,dMass);

  TLorentzVector LVec_total = LVec_target + LVec_beam;
  std::cout << LVec_total.M() << std::endl;

};
