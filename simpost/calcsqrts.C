void calcsqrts()
{

  TLorentzVector km;
  //km.SetPz(1.0);
  //km.SetM(0.493);
  TVector3 beamom;
  beamom.SetXYZ(0,0,1.0);
  km.SetVectM(beamom,0.493);
  km.M();
  km.E();
  TLorentzVector target;
  TVector3 tmom;
  tmom.SetXYZ(0,0,0);
  target.SetVectM(tmom,1.87561);
  TLorentzVector total;
  total = target + km;
  std::cout << total.M() << std::endl;
  std::cout << total.E() << std::endl;







}
