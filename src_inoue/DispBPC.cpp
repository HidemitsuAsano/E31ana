#include "DispBPC.h"

DispBPC::DispBPC(ConfMan *confMan) : fConfMan(confMan)
{
  fCanvas = new TCanvas("canvas", "canvas", 900, 450);
  fCanvas-> Divide(2, 1);

  fH2_xz = new TH2F("xz", "XZ Plane", 100, -7.5, 7.5, 100, -4.5, 4.5);
  fH2_xz-> SetStats(0);
  fH2_xz-> GetXaxis()-> CenterTitle();
  fH2_xz-> GetYaxis()-> CenterTitle();
  fH2_xz-> GetXaxis()-> SetTitleSize(0.05);
  fH2_xz-> GetYaxis()-> SetTitleSize(0.05);
  fH2_xz-> GetXaxis()-> SetTitle("Local Pos X[cm]");
  fH2_xz-> GetYaxis()-> SetTitle("Local Pos Y[cm]");

  fH2_yz = new TH2F("yz", "YZ Plane", 100, -7.5, 7.5, 100, -4.5, 4.5);
  fH2_yz-> SetStats(0);
  fH2_yz-> GetXaxis()->CenterTitle();
  fH2_yz-> GetYaxis()->CenterTitle();
  fH2_yz-> GetXaxis()-> SetTitleSize(0.05);
  fH2_yz-> GetYaxis()-> SetTitleSize(0.05);
  fH2_yz-> GetXaxis()-> SetTitle("Local Pos X[cm]");
  fH2_yz-> GetYaxis()-> SetTitle("Local Pos Y[cm]");
}

bool DispBPC::draw(BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BPCTrackMan &BPCTrack)
{
  fCanvas-> cd(1);
  fH2_xz-> Draw();

  fCanvas-> cd(2);
  fH2_yz-> Draw();

  TMarker mark;
  TArc arc;
  arc.SetFillColor(0);
  arc.SetFillStyle(0);
  arc.SetLineColor(kRed);

  BLDCWireMapMan *mapMan = fConfMan-> GetBLDCWireMapManager();
  for( int lay=1; lay<=8; lay++ ){
    BLDCWireMap *map = mapMan->GetWireMap(CID_BPC, lay);
    int xy = map-> GetXY();
    int nwire = map-> GetNWire();
    if( xy==0 ) fCanvas-> cd(1);
    else        fCanvas-> cd(2);
    double z = map-> GetZ();
    double xy0 = map-> GetXY0();
    double dxy = map-> GetdXY();

    mark.SetMarkerColor(kBlack);
    mark.SetMarkerStyle(3);
    for( int i=0; i<nwire; i++ ) mark.DrawMarker(xy0+i*dxy, z);

    mark.SetMarkerColor(kRed);
    mark.SetMarkerStyle(8);
    for( int i=0; i<blMan->nBPC(lay); i++ ){
      ChamberLikeHit *hit = blMan-> BPC(lay, i);
      double pos;
      if( xy==0 ) pos=hit-> wx();
      else pos=hit->wy();
      z = hit->wz();

      mark.DrawMarker(pos, z);
      arc.DrawArc(pos, z, hit->dl());
    }
  }
  TLine line;
  line.SetLineColor(kRed);
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    LocalTrack *track = bltrackMan-> trackBPC(i);
    double z0=-3, z1=3;
    double x0, x1, y0, y1;
    track-> XYLocalPosatZ(z0, x0, y0);
    track-> XYLocalPosatZ(z1, x1, y1);

    fCanvas-> cd(1);
    line.DrawLine(x0, z0, x1, z1);
    fCanvas-> cd(2);
    line.DrawLine(y0, z0, y1, z1);
  }
  fCanvas-> cd(1);
  for( int i=0; i<BPCTrack.nXTrack(); i++ ){
    XYTrack *track = BPCTrack.xTrack(i);
    double z0=-3, z1=3;
    double x0 = track->posatZ(z0);
    double x1 = track->posatZ(z1);
    line.SetLineColor(kBlue);
    line.DrawLine(x0, z0, x1, z1);
  }

  fCanvas-> cd(2);
  for( int i=0; i<BPCTrack.nYTrack(); i++ ){
    XYTrack *track = BPCTrack.yTrack(i);
    double z0=-3, z1=3;
    double y0 = track->posatZ(z0);
    double y1 = track->posatZ(z1);
    line.SetLineColor(kBlue);
    line.DrawLine(y0, z0, y1, z1);
  }


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
