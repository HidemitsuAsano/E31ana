#include "src/GlobalVariables.h"

static const double PI = 6*asin(0.5);

static const double CDS_z = 117;
static const double CDS_r =  59;

static const double Ushiwaka_posz = 250;
static const double Ushiwaka_z    = 140;
static const double Ushiwaka_x    = 240;

static const double FDC1_box_z = 14.6;
static const double FDC1_box_x = 57.8;

static const std::string particlename = "proton";

static const int NData=13;
static const double MOM_MIN=0.8;

void DrawTrack(double dmom=0.1)
{
  gSystem-> Load("src/lib/libAll.so");
  ConfMan *conf = new ConfMan("conf/Run47/analyzer.conf");
  conf-> Initialize();

  UshiwakaFieldMapMan *fieldMapMan = new UshiwakaFieldMapMan();
  fieldMapMan-> SetFileName("param/Run47/inoue/UshiwakaFieldMap.param");

  cout<<"=== USWK Simulation by Runge-Kutta Method ==="<<endl;
  fieldMapMan-> Initialize();
  fieldMapMan-> DumpStatus();
  cout<<"    Input position (x,y)=(0,0) fixed"<<endl;
  cout<<"    momentum direction is z fixed"<<endl;

  TCanvas *can = new TCanvas("can", "can", 600, 300);
  TH2F *h2_xz = new TH2F("h2_xz", "h2_xz", 100, 0, 1600, 100, -500, 200);
  h2_xz-> SetStats(false);
  h2_xz-> Draw();

  TBox box;
  box.SetFillStyle(0);
  double USWK_gpos[3];
  double USWK_range[3];

  //*** Draw CDS ***//
  double z1 = -CDS_z/2;
  double x1 = -CDS_r;
  double z2 = CDS_z/2;
  double x2 = CDS_r;
  box.SetLineColor(kBlue);
  box.DrawBox(z1,x1,z2,x2);

  //*** Draw Ushiwaka ***//
  fieldMapMan-> GetGPOS(USWK_gpos[0], USWK_gpos[1], USWK_gpos[2]);
  fieldMapMan-> GetRange(USWK_range[0], USWK_range[1], USWK_range[2]);

  box.SetLineColor(3);
  box.DrawBox(USWK_gpos[2]+USWK_range[2], USWK_gpos[0]+USWK_range[0],USWK_gpos[2]-USWK_range[2], USWK_gpos[0]-USWK_range[0]);

  z1 = Ushiwaka_posz - Ushiwaka_z/2;
  x1 = -Ushiwaka_x/2;
  z2 = Ushiwaka_posz + Ushiwaka_z/2;
  x2 = Ushiwaka_x/2;

  box.SetLineColor(kYellow);
  box.DrawBox(z1,x1,z2,x2);

  box.SetLineColor(kBlack);
  double gx, gy ,gz, dgx, dgy, dgz;
  double x, y, z, dx, dy, dz, len, wid, th;

  //*** Draw CVC***//
  conf-> GetGeomMapManager()-> GetGParam( CID_CVC, gx, gy, gz, dgx, dgy, dgz);
  double CVC_gposx = gx;
  double CVC_gposy = gy;
  double CVC_gposz = gz;
  for( int seg=1; seg<=NumOfCVCSegments; seg++ ){
    conf-> GetGeomMapManager()-> GetParam( CID_CVC, seg, x, y, z, dx, dy, dz, len, wid, th);
    z1 = CVC_gposz + z - th/2;
    x1 = CVC_gposx + x - wid/2;
    z2 = CVC_gposz + z + th/2;
    x2 = CVC_gposx + x + wid/2;
    box.DrawBox(z1,x1,z2,x2);
  }
  //*** Draw NC ***//
  conf-> GetGeomMapManager()-> GetGParam( CID_NC, gx, gy, gz, dgx, dgy, dgz);
  double NC_gposx = gx;
  double NC_gposy = gy;
  double NC_gposz = gz;
  for( int seg=1; seg<=NumOfNCLayer*NumOfNCSegmentsInLayer; seg++ ){
    conf-> GetGeomMapManager()-> GetParam( CID_NC, seg, x, y, z, dx, dy, dz, len, wid, th);
    z1 = NC_gposz + z - th/2;
    x1 = NC_gposx + x - wid/2;
    z2 = NC_gposz + z + th/2;
    x2 = NC_gposx + x + wid/2;
    box.DrawBox(z1,x1,z2,x2);
  }
  //*** Draw PC ***//
  conf-> GetGeomMapManager()-> GetGParam( CID_PC, gx, gy, gz, dgx, dgy, dgz);
  double PC_gposx = gx;
  double PC_gposy = gy;
  double PC_gposz = gz;

  double PC_z1;
  double PC_x1;
  double PC_z2;
  double PC_x2;

  for( int seg=1; seg<=NumOfPCSegments; seg++ ){
    conf-> GetGeomMapManager()-> GetParam( CID_PC, seg, x, y, z, dx, dy, dz, len, wid, th);
    x1 = PC_gposx + cos(PI*dgy/180.)*(x-wid/2) + sin(PI*dgy/180.)*(z-th/2);
    x2 = PC_gposx + cos(PI*dgy/180.)*(x+wid/2) + sin(PI*dgy/180.)*(z+th/2);
    z1 = PC_gposz - sin(PI*dgy/180.)*(x-wid/2) + cos(PI*dgy/180.)*(z-th/2);
    z2 = PC_gposz - sin(PI*dgy/180.)*(x+wid/2) + cos(PI*dgy/180.)*(z+th/2);
    box.DrawBox(z1,x1,z2,x2);

    if( seg==1 ){
      PC_z1 = z1;
      PC_x1 = x1;
    }
    if( seg==NumOfPCSegments ){
      PC_z2 = z1;
      PC_x2 = x1;
    }
  }
  double PC_court_A = (PC_x2-PC_x1)/(PC_z2-PC_z1);
  double PC_court_B = PC_x1 - PC_court_A*PC_z1;

  TLine line;
  line.SetLineColor(4);
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(2);

  x1 = -500;
  x2 = 0;
  z1 = (x1-PC_court_B)/PC_court_A;
  z2 = (x2-PC_court_B)/PC_court_A;
  line.DrawLine(z1,x1,z2,x2);

  int color_id =1;
  for( int i=0; i<NData; i++ ){
    color_id++;
    if( color_id==10) color_id=2;
    mark.SetMarkerColor(color_id);
    line.SetLineColor(color_id);
    double XA = 0, XB=0, YA=0, YB=0;
    double InPosZ = USWK_gpos[2]-USWK_range[2];
    double InPosX = XA*InPosZ+XB;
    double InPosY = YA*InPosZ+YB;

    TVector3 in_pos(InPosX, InPosY, InPosZ);
    double AngX = -atan(YA);
    double AngY = atan(XA);

    double mom = MOM_MIN + i*dmom;
    TVector3 in_mom(0.0, 0.0, mom);
    in_mom.RotateX(AngX);
    in_mom.RotateY(AngY);

    fieldMapMan-> RungeKutta(particlename, in_pos, in_mom);

    USWKTrack *uswk_track = fieldMapMan-> GetUSWKTrack(0);
    double mass = uswk_track-> mass();
    double out_time = uswk_track-> outtime();
    TVector3 out_pos = uswk_track-> outposition();
    TVector3 out_mom = uswk_track-> outmomentum();

    double in_dxdz = in_mom.X()/in_mom.Z();
    double in_dydz = in_mom.Y()/in_mom.Z();
    double in_x0 = in_pos.X()-in_dxdz*in_pos.Z();
    double in_y0 = in_pos.Y()-in_dydz*in_pos.Z();
    double in_abs_mom = in_mom.Mag();
    double in_v = 100.*Const*in_abs_mom/sqrt(mass*mass+in_abs_mom*in_abs_mom);
    double in_ang = atan2(in_mom.X(),in_mom.Z());
    double zz0 = 0.;
    double xx0 = in_pos.X()-in_dxdz*in_pos.Z();
    double yy0 = in_pos.Y()-in_dydz*in_pos.Z();

    double out_dxdz = out_mom.X()/out_mom.Z();
    double out_dydz = out_mom.Y()/out_mom.Z();
    double out_x0 = out_pos.X()-out_dxdz*out_pos.Z();
    double out_y0 = out_pos.Y()-out_dydz*out_pos.Z();
    double out_abs_mom = out_mom.Mag();
    double out_v = 100.*Const*out_abs_mom/sqrt(mass*mass+out_abs_mom*out_abs_mom);
    double out_ang = atan2(out_mom.X(),out_mom.Z());

    double bending_ang = 180.*(in_ang-out_ang)/PI;
    double PC_hitz = (PC_court_B-out_x0)/(out_dxdz-PC_court_A);
    double PC_hity = out_dydz*PC_hitz+out_y0; 
    double PC_hitx = out_dxdz*PC_hitz+out_x0; 

    TVector3 vtx_pos(xx0, yy0, zz0);
    double in_fl = (in_pos-vtx_pos).Mag();
    double in_tof = in_fl/in_v;

    double uswk_fl = out_v*out_time;

    TVector3 PC_hitpos(PC_hitx, PC_hity, PC_hitz);
    double out_fl = (out_pos-PC_hitpos).Mag();
    double out_tof = out_fl/out_v;

    double fl = in_fl+uswk_fl+out_fl;
    double tof = in_tof+out_time+out_tof;

    //*** Check Method ***//
    if( fabs(in_abs_mom-out_abs_mom)>10e-6 ){
      cout<<"not match input momentum & output momentum !!!"<<endl;
    }

    mark.DrawMarker(zz0, xx0);
    mark.DrawMarker(in_pos.Z(), in_pos.X());
    mark.DrawMarker(out_pos.Z(), out_pos.X());
    mark.DrawMarker(PC_hitz, PC_hitx);

    cout<<"  > Input momentum : "<<mom<<" [GeV/c] ==="<<endl;
    cout<<"  > Input particle : "<<uswk_track-> particlename()<<" mass : "<<mass<<endl;
    cout<<"    Num Of Point : "<<uswk_track->NPoint()<<endl;
    cout<<"    angle:"<<bending_ang<<"[degree]"<<endl;
    cout<<"    fligth length:"<<fl<<"[cm]"<<endl;
    cout<<"    TOF:"<<tof<<"[ns]"<<endl;

    TVector3 pos1;
    line.DrawLine(zz0, xx0, in_pos.Z(), in_pos.X());
    TVector3 pos = uswk_track-> position(0);
    line.DrawLine(in_pos.Z(), in_pos.X(), pos.Z(), pos.X());
    for( int i_p=1; i_p<uswk_track->NPoint(); i_p++ ){
      pos = uswk_track-> position(i_p);
      pos1 = uswk_track-> position(i_p-1);
      line.DrawLine(pos1.Z(), pos1.X(), pos.Z(), pos.X());
    }
    double XA2 = out_mom.X()/out_mom.Z();
    double zz = 1500;
    double xx = out_pos.X()+XA2*(1500-out_pos.Z());
    line.DrawLine(out_pos.Z(), out_pos.X(), zz, xx);

    fieldMapMan-> Clear();
  }
}
