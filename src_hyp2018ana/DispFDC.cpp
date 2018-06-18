#include "DispFDC.h"

DispFDC::DispFDC(ConfMan *confMan) : fConfMan(confMan)
{
  fCanvas = new TCanvas("canvas", "canvas", 600, 600);
  fCanvas-> Divide(1, 3);

  fH2_U = new TH2F("U_Plane", "U Plane", 100, -20, 20, 100, -5, 5);
  fH2_U-> SetStats(0);
  fH2_U-> GetXaxis()-> CenterTitle();
  fH2_U-> GetYaxis()-> CenterTitle();
  fH2_U-> GetXaxis()-> SetTitleSize(0.05);
  fH2_U-> GetYaxis()-> SetTitleSize(0.05);
  fH2_U-> GetXaxis()-> SetTitle("Local Pos X[cm]");
  fH2_U-> GetXaxis()-> SetTitle("Local Pos Z[cm]");

  fH2_X = new TH2F("X_Plane", "X Plane", 100, -20, 20, 100, -5, 5);
  fH2_X-> SetStats(0);
  fH2_X-> GetXaxis()-> CenterTitle();
  fH2_X-> GetYaxis()-> CenterTitle();
  fH2_X-> GetXaxis()-> SetTitleSize(0.05);
  fH2_X-> GetYaxis()-> SetTitleSize(0.05);
  fH2_X-> GetXaxis()-> SetTitle("Local Pos X[cm]");
  fH2_X-> GetXaxis()-> SetTitle("Local Pos Z[cm]");

  fH2_V = new TH2F("V_Plane", "V Plane", 100, -20, 20, 100, -5, 5);
  fH2_V-> SetStats(0);
  fH2_V-> GetXaxis()-> CenterTitle();
  fH2_V-> GetYaxis()-> CenterTitle();
  fH2_V-> GetXaxis()-> SetTitleSize(0.05);
  fH2_V-> GetYaxis()-> SetTitleSize(0.05);
  fH2_V-> GetXaxis()-> SetTitle("Local Pos X[cm]");
  fH2_V-> GetXaxis()-> SetTitle("Local Pos Z[cm]");

}

bool DispFDC::draw(BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan)
{
  fCanvas-> Update();
  std::cout<<"Please input any word.  q: quit  s: save"<<std::endl;
  char in;
  std::cin>>in;
  if( in=='q' ) return false;
  if( in=='s' ){
    fCanvas-> Print("tmp.png");
  }
  else return true;;
}
