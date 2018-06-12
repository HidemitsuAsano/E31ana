#include "DisplayBLC.h"

DisplayBLC::DisplayBLC()
{
  fCanvas = new TCanvas("c1", "BLC2 Event Display", 900, 500);
  fCanvas-> Divide(2, 2);

  fBLC2a_xz = new TH2F("BLC2a_XZ", "BLC2 a XZ Plane", 100, -10, 10, 100, -17, -10);
  fBLC2a_xz-> SetStats(0);
  fBLC2a_xz-> GetXaxis()-> CenterTitle();
  fBLC2a_xz-> GetYaxis()-> CenterTitle();
  fBLC2a_xz-> GetXaxis()-> SetTitleSize(0.05);
  fBLC2a_xz-> GetYaxis()-> SetTitleSize(0.05);
  fBLC2a_xz-> GetXaxis()-> SetTitle("Pos X[cm]");
  fBLC2a_xz-> GetYaxis()-> SetTitle("Pos Z[cm]");

  fBLC2a_yz = new TH2F("BLC2a_YZ", "BLC2 a YZ Plane", 100, -10, 10, 100, -17, -10);
  fBLC2a_yz-> SetStats(0);
  fBLC2a_yz-> GetXaxis()-> CenterTitle();
  fBLC2a_yz-> GetYaxis()-> CenterTitle();
  fBLC2a_yz-> GetXaxis()-> SetTitleSize(0.05);
  fBLC2a_yz-> GetYaxis()-> SetTitleSize(0.05);
  fBLC2a_yz-> GetXaxis()-> SetTitle("Pos Y[cm]");
  fBLC2a_yz-> GetYaxis()-> SetTitle("Pos Z[cm]");

  fBLC2b_xz = new TH2F("BLC2b_XZ", "BLC2 b XZ Plane", 100, -10, 10, 100, 10, 17);
  fBLC2b_xz-> SetStats(0);
  fBLC2b_xz-> GetXaxis()-> CenterTitle();
  fBLC2b_xz-> GetYaxis()-> CenterTitle();
  fBLC2b_xz-> GetXaxis()-> SetTitleSize(0.05);
  fBLC2b_xz-> GetYaxis()-> SetTitleSize(0.05);
  fBLC2b_xz-> GetXaxis()-> SetTitle("Pos X[cm]");
  fBLC2b_xz-> GetYaxis()-> SetTitle("Pos Z[cm]");

  fBLC2b_yz = new TH2F("BLC2b_YZ", "BLC2 b YZ Plane", 100, -10, 10, 100, 10, 17);
  fBLC2b_yz-> SetStats(0);
  fBLC2b_yz-> GetXaxis()-> CenterTitle();
  fBLC2b_yz-> GetYaxis()-> CenterTitle();
  fBLC2b_yz-> GetXaxis()-> SetTitleSize(0.05);
  fBLC2b_yz-> GetYaxis()-> SetTitleSize(0.05);
  fBLC2b_yz-> GetXaxis()-> SetTitle("Pos Y[cm]");
  fBLC2b_yz-> GetYaxis()-> SetTitle("Pos Z[cm]");
}

void DisplayBLC::Draw(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan)
{
  fCanvas-> cd(3);
  fBLC2a_xz-> Draw();

  fCanvas-> cd(4);
  fBLC2a_yz-> Draw();

  fCanvas-> cd(1);
  fBLC2b_xz-> Draw();

  fCanvas-> cd(2);
  fBLC2b_yz-> Draw();

  TMarker mark;

  int nlayer=0;
  mark.SetMarkerColor(kBlack);
  mark.SetMarkerStyle(3);
  nlayer=DetectorList::GetInstance()->GetNlayers(CID_BLC2a);
  for( int lay=1; lay<=nlayer; lay++ ){
    int nwire, xy;
    double z, xy0, dxy, wl, tilt, ra;
    conf-> GetBLDCWireMapManager()-> GetParam(CID_BLC2a, lay, nwire, z, xy, xy0, dxy, wl, tilt, ra);
    if( xy==0 ) fCanvas-> cd(3);
    else fCanvas-> cd(4);
    for( int wire=0; wire<nwire; wire++ ) mark.DrawMarker(xy0+wire*dxy, z);
  }
  mark.SetMarkerStyle(8);
  mark.SetMarkerColor(kRed);
  for( int lay=1; lay<=nlayer; lay++ ){
    for( int i=0; i<blMan->nBLDC(CID_BLC2a, lay); i++ ){
      ChamberLikeHit *hit = blMan->BLDC(CID_BLC2a, lay, i);
      if( hit->xy()==0 ){
	fCanvas-> cd(3);
	mark.DrawMarker(hit->wx(), hit->wz());
      }
      else{
	fCanvas-> cd(4);
	mark.DrawMarker(hit->wy(), hit->wz());
      }
    }
  }

  mark.SetMarkerColor(kBlack);
  mark.SetMarkerStyle(3);
  nlayer=DetectorList::GetInstance()->GetNlayers(CID_BLC2b);
  for( int lay=1; lay<=nlayer; lay++ ){
    int nwire, xy;
    double z, xy0, dxy, wl, tilt, ra;
    conf-> GetBLDCWireMapManager()-> GetParam(CID_BLC2b, lay, nwire, z, xy, xy0, dxy, wl, tilt, ra);
    if( xy==0 ) fCanvas-> cd(1);
    else fCanvas-> cd(2);
    for( int wire=0; wire<nwire; wire++ ) mark.DrawMarker(xy0+wire*dxy, z);
  }

  mark.SetMarkerStyle(8);
  mark.SetMarkerColor(kBlue);
  for( int lay=1; lay<=nlayer; lay++ ){
    for( int i=0; i<blMan->nBLDC(CID_BLC2b, lay); i++ ){
      ChamberLikeHit *hit = blMan->BLDC(CID_BLC2b, lay, i);
      if( hit->xy()==0 ){
	fCanvas-> cd(1);
	mark.DrawMarker(hit->wx(), hit->wz());
      }
      else{
	fCanvas-> cd(2);
	mark.DrawMarker(hit->wy(), hit->wz());
      }
    }
  }

  TLine line;
  line.SetLineColor(kRed);
  for( int i=0; i<bltrackMan-> ntrackBLDC(CID_BLC1a); i++ ){
    LocalTrack *track = bltrackMan-> trackBLDC(CID_BLC1a, i);
    double z0=-17, z1=-10, z2=10, z3=17;
    double x0, y0, x1, y1, x2, y2, x3, y3;
    track-> XYLocalPosatZ(z0, x0, y0);
    track-> XYLocalPosatZ(z1, x1, y1);
    track-> XYLocalPosatZ(z2, x2, y2);
    track-> XYLocalPosatZ(z3, x3, y3);

    line.SetLineStyle(1);
    fCanvas-> cd(3);
    line.DrawLine(x0, z0, x1, z1);

    fCanvas-> cd(4);
    line.DrawLine(y0, z0, y1, z1);

    line.SetLineStyle(8);
    fCanvas-> cd(1);
    line.DrawLine(x2, z2, x3, z3);

    fCanvas-> cd(2);
    line.DrawLine(y2, z2, y3, z3);
  }

  for( int i=0; i<bltrackMan-> ntrackBLDC(CID_BLC1b); i++ ){
    LocalTrack *track = bltrackMan-> trackBLDC(CID_BLC1b, i);
    double z0=-17, z1=-10, z2=10, z3=17;
    double x0, y0, x1, y1, x2, y2, x3, y3;
    track-> XYLocalPosatZ(z0, x0, y0);
    track-> XYLocalPosatZ(z1, x1, y1);
    track-> XYLocalPosatZ(z2, x2, y2);
    track-> XYLocalPosatZ(z3, x3, y3);

    line.SetLineStyle(8);
    fCanvas-> cd(3);
    line.DrawLine(x0, z0, x1, z1);

    fCanvas-> cd(4);
    line.DrawLine(y0, z0, y1, z1);

    line.SetLineStyle(1);
    fCanvas-> cd(1);
    line.DrawLine(x2, z2, x3, z3);

    fCanvas-> cd(2);
    line.DrawLine(y2, z2, y3, z3);
  }
}

bool DisplayBLC::Wait()
{
  fCanvas-> Update();
  std::cout<<" Plase input any word. q:quit"<<std::endl;
  std::string input;
  std::cin>>input;
  if( input=="q" ) return false;
  return true;
}
