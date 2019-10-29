void test_nmomshift()
{
  TLorentzVector tln ;
  tln.SetXYZM(100,0,0,940);

  TLorentzVector tlp ;
  tlp.SetXYZM(0,-300,0,140);

  TLorentzVector tls = tln+tlp;
  std::cout << tls.M() << std::endl;

}
