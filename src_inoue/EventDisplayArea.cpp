#include <TCanvas.h>
#include <TApplication.h>
#include <TRint.h>
#include <TBox.h>

#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#include "DataContainer.h"
#include "AnalysisMan.h"

#include "UshiwakaFieldMapMan.h"
#include "UshiwakaTableMan.h"

static const double PI = 6*asin(0.5);

//*** for BLC2 Parameter ***//
static const double BLC2_box_z = 4.00;
static const double BLC2_box_x = 16.5;

//*** for CDC Paramterer ***//                      
static const double CDS_z = 117;
static const double CDS_r =  59;

//**** for Ushiwaka Parameter ***//                 
static const double Ushiwaka_posz = 250;
static const double Ushiwaka_z    = 140;
static const double Ushiwaka_x    = 240;

static const double FieldMin = 150;
static const double FieldMax = 350;

//*** for FDC1 Parameter ***//                      
static const double FDC1_box_z = 14.6;
static const double FDC1_box_x = 57.8;

//*** for BHD-T0 TOF parameter ***//
static const double BHDT0_PI_MIN = 25.156;
static const double BHDT0_PI_MAX = 26.8522;

static const double BHDT0_K_MIN = 27.3686;
static const double BHDT0_K_MAX = 28.6376;

static const double BHDT0_P_MIN = 33.964;
static const double BHDT0_P_MAX = 36.572;

//*** for Forward Charged Particle Parameter ***//
static const double P_MASS_MIN = 0.70;
static const double P_MASS_MAX = 1.30;

static const double D_MASS_MIN = 1.30;
static const double D_MASS_MAX = 2.30;

class EventAnalysisCheck: public EventTemp
{
public:
  EventAnalysisCheck();
  ~EventAnalysisCheck();
private:
  CDSHitMan        *cdsMan;
  BeamLineHitMan   *blMan;
  EventHeader      *header;
  ScalerMan        *scaMan;

  DataContainer    *data;

  UshiwakaFieldMapMan *fieldMan;
  UshiwakaTableMan *tableMan;

  TCanvas *can_area_xz;
  TCanvas *can_area_yz;
  TH2F *h2_area_xz;
  TH2F *h2_area_yz;

  // PC xz plane define x = PC_A*z + PC_B
  double PC_A;
  double PC_B;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  // for Draw Method
  void DrawArea();
  bool DrawHit();
};

EventAnalysisCheck::EventAnalysisCheck()
  : EventTemp()
{
}

EventAnalysisCheck::~EventAnalysisCheck()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisCheck::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisCheck::Initialize " << std::endl;
#endif
  int argc2;
  char **argv2;
  TRint *theApp = new TRint( "theApp", &argc2, argv2 );
  confMan = conf;
  
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  can_area_xz = new TCanvas("can_area_xz", "can_area_xz", 850, 350);
  can_area_yz = new TCanvas("can_area_yz", "can_area_yz", 850, 200);
  h2_area_xz = new TH2F("h2_area_xz", "Area XZ Plane", 100, -150, 1550, 100, -500, 200);
  h2_area_xz-> SetStats(false);
  h2_area_yz = new TH2F("h2_area_yz", "Area YZ Plane", 100, -150, 1550, 100, -100, 100);
  h2_area_yz-> SetStats(false);

  fieldMan = new UshiwakaFieldMapMan();
  fieldMan-> Initialize();
  tableMan = new UshiwakaTableMan();
  tableMan-> Initialize();

  data = new DataContainer();
}

void EventAnalysisCheck::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisCheck::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  std::cout << nsca<< std::endl;
  for( int i=0; i<nsca; i++ ){
   std::cout << "  " << sca[i];
  }
  std::cout << std::endl;
#endif
  scaMan->Clear();
}

bool EventAnalysisCheck::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisCheck::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }
    
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

  data-> SetHodoscope( blMan, cdsMan );
  if( data-> nCVC()>0 || data-> nPC()>0 ){
    data-> DoTrackingBL( blMan, confMan );
    data-> ExecuteCDS( cdsMan, confMan );
    data-> CalcVertex_beam( confMan );

    std::cout << "=== Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " ===" << std::endl;
    DrawArea();
    if( DrawHit() ){
      can_area_xz-> Update();
      can_area_yz-> Update();
      char dispflag;
      char dispin[100]="";
      std::cout << " Input any word or return" << std::endl;
      std::cout << " (q:quit)" << std::endl;
      fgets(dispin,100,stdin);
      if(sscanf(dispin,"%c",&dispflag)==1){
	if( dispflag=='q' ){
	  return false;
	}
	if( dispflag=='s' ){

	}
      }
    }
  }

  data-> Clear();

  header-> Clear();
  blMan-> Clear();
  cdsMan-> Clear();

  return true;
}

bool EventAnalysisCheck::DrawHit()
{
  can_area_xz-> cd();

  BeamLineTrackMan *bltrackMan = data-> GetBLTrackMan();
  CDSTrackingMan *trackingMan = data-> GetCDSTrackingMan();

  TMarker mark;
  mark.SetMarkerStyle(7);
  mark.SetMarkerColor(2);
  mark.SetMarkerSize(10);

  TLine line;
  line.SetLineColor(2);

  std::cout<<" BPC track is "<<bltrackMan->ntrackBPC()<<std::endl;
  std::cout<<" CDS Good Track is "<<trackingMan-> nGoodTrack()<<std::endl;
  std::cout<<" FDC1 track is "<<bltrackMan->ntrackFDC1()<<std::endl;
  std::cout<<" T0  hit is "<<data->nT0()<<std::endl;
  std::cout<<" PC  hit is "<<data->nPC()<<std::endl;
  std::cout<<" CVC hit is "<<data->nCVC()<<std::endl;

  bool PC1hit = false;
  bool CVC1hit = false;
  if( data->nPC()==1 && data->nCVC()==0 ) PC1hit = true;
  if( data->nCVC()==1 && data->nPC()==0 ) CVC1hit = true;

  TVector3 pos, dir, gpos, gdir;
  TVector3 hit_pos;
  for( int i=0; i<data->nCVC(); i++ ){
    int seg = data-> CVC(i)-> seg();
    confMan-> GetGeomMapManager()-> GetGParam(CID_CVC, gpos, gdir);
    confMan-> GetGeomMapManager()-> GetParam(CID_CVC, seg, pos, dir);
    mark.DrawMarker( gpos.z()+pos.z(), gpos.x()+pos.x() );

    hit_pos.SetXYZ( gpos.x()+pos.x(), gpos.y()+pos.y(), gpos.z()+pos.z() );
  }

  for( int i=0; i<data->nPC(); i++ ){
    int seg = data-> PC(i)-> seg();
    confMan-> GetGeomMapManager()-> GetGParam(CID_PC, gpos, gdir);
    confMan-> GetGeomMapManager()-> GetParam(CID_PC, seg, pos, dir);

    double PCz = gpos.z() + cos(PI*gdir.y()/180.)*pos.z() - sin(PI*gdir.y()/180.)*pos.x();
    double PCx = gpos.x() + sin(PI*gdir.y()/180.)*pos.z() + cos(PI*gdir.y()/180.)*pos.x();
    mark.DrawMarker( PCz, PCx );

    hit_pos.SetXYZ( PCx, gpos.y()+pos.y(), PCz );
  }

  for( int i=0; i<data->nNC(); i++ ){
    int seg = data-> NC(i)-> seg();
    confMan-> GetGeomMapManager()-> GetGParam(CID_NC, gpos, gdir);
    confMan-> GetGeomMapManager()-> GetParam(CID_NC, seg, pos, dir);
    mark.DrawMarker( gpos.z()+pos.z(), gpos.x()+pos.x() );
  }

  for( int bpc_id=0; bpc_id<bltrackMan->ntrackBPC(); bpc_id++ ){
    bool vtx_flag = false;
    TrackVertex vtx_beam;
    double tmp_dis = 9999;
    for( int cdc_id=0; cdc_id<trackingMan->nGoodTrack(); cdc_id++ ){
      int good_id = trackingMan->GoodTrackID(cdc_id);
      TrackVertex tmp_vtx;
      if( trackingMan-> GetVertex_beam(bpc_id, good_id, tmp_vtx) ){
	TVector3 pos1 = tmp_vtx.GetVertexPos1();
	TVector3 pos2 = tmp_vtx.GetVertexPos2();
	double tmp = (pos1-pos2).Mag();
	if( tmp<tmp_dis ){
	  vtx_flag = true;
	  tmp_dis = tmp;
	  vtx_beam = tmp_vtx;
	}
      }
    }
    if( vtx_flag ){
      TVector3 vtx_pos = vtx_beam.GetVertexPos_mean();
      double T0z = -110;
      double T0x, T0y;
      LocalTrack *trackBPC = bltrackMan-> trackBPC(bpc_id);
      trackBPC-> XYPosatZ(T0z, T0x, T0y);
      line.DrawLine(T0z, T0x, vtx_pos.z(), vtx_pos.x());
      TVector3 T0_pos(T0x, T0y, T0z);
      double t0vtx_fl = (vtx_pos-T0_pos).Mag();

      if( bltrackMan-> ntrackFDC1()==1 ){
	LocalTrack *trackFDC1 = bltrackMan-> trackFDC1(0);
	double FDC1z = trackFDC1-> hit(0)->gz();
	double FDC1x, FDC1y;
	trackFDC1-> XYPosatZ(FDC1z, FDC1x, FDC1y);
	mark.DrawMarker(FDC1z, FDC1x);
	TVector3 FDC1_pos(FDC1x, FDC1y, FDC1z);

	double in_dxdz = (vtx_pos.x()-FDC1x)/(vtx_pos.z()-FDC1z);
	double in_dydz = (vtx_pos.y()-FDC1y)/(vtx_pos.z()-FDC1z);
	double in_x0 = vtx_pos.x()-in_dxdz*vtx_pos.z();
	double in_y0 = vtx_pos.y()-in_dydz*vtx_pos.z();

	double USWK_z = Ushiwaka_posz;
	double USWK_x = in_dxdz*Ushiwaka_posz + in_x0;
	double USWK_y = in_dydz*Ushiwaka_posz + in_y0;

	line.DrawLine(vtx_pos.z(), vtx_pos.x(), USWK_z, USWK_x);

	if( data->nT0()==1 ){
	  double bhdt0_tof = -999;
	  double beam_v = 0;
	  double t0vtx_tof = -999;
	  if( data->nBHD()==1 ) bhdt0_tof = data->T0(0)->ctmean()-data->BHD(0)->ctmean();
	  bool bhdt0_pi = false;
	  bool bhdt0_k = false;
	  if( BHDT0_PI_MIN<bhdt0_tof && bhdt0_tof<BHDT0_PI_MAX ){
	    bhdt0_pi=true;
	    beam_v = 100.*Const*1.0/sqrt(piMass*piMass+1.0);
	    t0vtx_tof = t0vtx_fl/beam_v;
	    std::cout<<" pi beam"<<std::endl;
	  }
	  else if( BHDT0_K_MIN<bhdt0_tof && bhdt0_tof<BHDT0_K_MAX ){
	    bhdt0_k=true;
	    beam_v = 100.*Const*1.0/sqrt(kpMass*kpMass+1.0);
	    t0vtx_tof = t0vtx_fl/beam_v;
	    std::cout<<" k beam"<<std::endl;
	  }
	  else{
	    return false;
	    std::cout<<" unknown beam"<<std::endl;
	  }

	  double tof = -999;
	  if( PC1hit ){
	    tof = data->PC(0)->ctmean()  - data->T0(0)->ctmean();
	    std::cout<<" T0-PC  : "<<tof<<"[ns]"<<std::endl;
	  }
	  if( CVC1hit ){
	    tof = data->CVC(0)->ctmean() - data->T0(0)->ctmean();
	    std::cout<<" T0-CVC : "<<tof<<"[ns]"<<std::endl;
	  }
	  if( PC1hit || CVC1hit ){
	    double in_ang  = atan( in_dxdz );
	    double out_ang = atan( (hit_pos.x()-USWK_x)/(hit_pos.z()-USWK_z) );
	    double bending_angle = in_ang - out_ang;
	    line.DrawLine( USWK_z, USWK_x, hit_pos.z(), hit_pos.x() );

	    double in_x = in_dxdz*FieldMin + in_x0;
	    double in_y = in_dxdz*FieldMin + in_y0;

	    double mom, fl;
	    if( tableMan-> GetParam( in_x, in_y, in_dxdz, in_dydz, bending_angle, mom, fl) ){ 
	      std::cout<<" bending angle : "<<180*bending_angle/PI<<"[degree]"<<std::endl;
	      std::cout<<" momentum by angle : "<<mom<<"[GeV/c]"<<std::endl;

	      TVector3 in_pos(in_x, in_y, FieldMin);
	      TVector3 in_mom(in_dxdz, in_dydz, 1.0);
	      in_mom.SetMag(mom);
	      double vtx_in_fl = (in_pos-vtx_pos).Mag();

	      std::string particlename = "proton";
	      fieldMan-> RungeKutta(particlename, in_pos, in_mom);
	      USWKTrack *uswk_track = fieldMan-> GetUSWKTrack(0);

	      line.SetLineColor(3);
	      for( int i=1; i<uswk_track->NPoint(); i++ ){
		TVector3 pos1 = uswk_track-> position(i);
		TVector3 pos2 = uswk_track-> position(i-1);
		line.DrawLine(pos1.z(), pos1.x(), pos2.z(), pos2.x());
	      }
	      TVector3 out_pos = uswk_track-> outposition();
	      TVector3 out_mom = uswk_track-> outmomentum();
	      double out_time = uswk_track-> outtime();
	      if( fabs(in_mom.Mag()-out_mom.Mag())>0.0001 ) std::cout<<" differ in momentum & out momentum !!!"<<std::endl;
	      double velocity = 100.*Const*mom/sqrt(pMass*pMass+mom*mom);
	      double uswk_fl = out_time*velocity;

	      double out_dxdz = out_mom.x()/out_mom.z();
	      double out_x0 = out_pos.x()-out_dxdz*out_pos.z();
	      double out_dydz = out_mom.y()/out_mom.z();
	      double out_y0 = out_pos.y()-out_dydz*out_pos.z();

	      double calc_hitpos_z = (out_x0-PC_B)/(PC_A-out_dxdz);
	      double calc_hitpos_x = out_dxdz*calc_hitpos_z+out_x0;
	      double calc_hitpos_y = out_dydz*calc_hitpos_z+out_y0;
	      TVector3 calc_hitpos(calc_hitpos_x, calc_hitpos_y, calc_hitpos_z);
	      double out_fl = (calc_hitpos-out_pos).Mag();

	      line.DrawLine(out_pos.z(), out_pos.x(), calc_hitpos.z(), calc_hitpos.x());

	      if( bhdt0_pi || bhdt0_k ){
		double fl = vtx_in_fl + uswk_fl + out_fl;
		tof -= t0vtx_tof;
		double calc_beta = fl/(tof*100.*Const);
		double calc_mass = mom*sqrt(1-calc_beta*calc_beta)/calc_beta;

		std::cout<<" calc mass : "<<calc_mass<<"[GeV]"<<std::endl;
		double calc_mom = -999;
		if( P_MASS_MIN<calc_mass && calc_mass<P_MASS_MAX ){
		  calc_mom = pMass*calc_beta/sqrt(1-calc_beta*calc_beta);
		  std::cout<<" forward charged particle 'proton' "<<std::endl;
		  std::cout<<" momentum by TOF : "<<calc_mom<<"[GeV/c]"<<std::endl;
		}
		if( D_MASS_MIN<calc_mass && calc_mass<D_MASS_MAX ){
		  calc_mom = dMass*calc_beta/sqrt(1-calc_beta*calc_beta);
		  std::cout<<" forward charged particle 'deuteron' "<<std::endl;
		  std::cout<<" momentum by TOF : "<<calc_mom<<"[GeV/c]"<<std::endl;
		}
	      }
	      fieldMan-> Clear();
	      return true;
	    }
	    else std::cout<<" fault UshiwakaTableMan::GetParam("<<in_x<<","<<in_y<<","<<in_dxdz<<","<<in_dydz<<","
			  <<bending_angle<<","<<mom<<","<<fl<<")"<<std::endl;
	  }
	}
      }
    }
    else std::cout<<" not dicide vertex"<<std::endl;
  }
  return false;
}

void EventAnalysisCheck::DrawArea()
{
  can_area_xz-> cd();
  h2_area_xz-> Draw();
  // for CDS
  TBox box;
  box.SetFillStyle(0);

  double z1 = -CDS_z/2;
  double x1 = -CDS_r;
  double z2 = CDS_z/2;
  double x2 = CDS_r;
  box.SetLineColor(kBlue);
  box.DrawBox(z1,x1,z2,x2);
  // for Ushiwaka
  z1 = Ushiwaka_posz - Ushiwaka_z/2;
  x1 = -Ushiwaka_x/2;
  z2 = Ushiwaka_posz + Ushiwaka_z/2;
  x2 = Ushiwaka_x/2;

  box.SetLineColor(kYellow);
  box.DrawBox(z1,x1,z2,x2);
  // for FieldMap
  double field_x, field_y, field_z;
  double range_x, range_y, range_z;
  fieldMan-> GetGPOS(field_x, field_y, field_z);
  fieldMan-> GetRange(range_x, range_y, range_z);

  box.SetLineColor(3);
  box.DrawBox( field_z-range_z, field_x-range_x, field_z+range_z, field_x+range_x );

  TVector3 pos, dir, gpos, gdir;
  double len, wid, th, lv;

  box.SetLineColor(kBlack);
  // for CVC
  confMan-> GetGeomMapManager()-> GetGParam( CID_CVC, gpos, gdir);
  for( int seg=1; seg<=NumOfCVCSegments; seg++ ){
    confMan-> GetGeomMapManager()-> GetParam( CID_CVC, seg, pos, dir, len, wid, th, lv );
    z1 = gpos.z() + pos.z() - th/2;
    x1 = gpos.x() + pos.x() - wid/2;
    z2 = gpos.z() + pos.z() + th/2;
    x2 = gpos.x() + pos.x() + wid/2;
    box.DrawBox(z1,x1,z2,x2);
  }
  // for NC
  confMan-> GetGeomMapManager()-> GetGParam( CID_NC, gpos, gdir);
  for( int seg=1; seg<=NumOfNCSegments; seg++ ){
    confMan-> GetGeomMapManager()-> GetParam( CID_NC, seg, pos, dir, len, wid, th, lv );
    z1 = gpos.z() + pos.z() - th/2;
    x1 = gpos.x() + pos.x() - wid/2;
    z2 = gpos.z() + pos.z() + th/2;
    x2 = gpos.x() + pos.x() + wid/2;
    box.DrawBox(z1,x1,z2,x2);
  }
  // for PC
  double PC_x1, PC_x2, PC_z1, PC_z2;
  confMan-> GetGeomMapManager()-> GetGParam( CID_PC, gpos, gdir);
  for( int seg=1; seg<=NumOfPCSegments; seg++ ){
    confMan-> GetGeomMapManager()-> GetParam( CID_PC, seg, pos, dir, len, wid, th, lv );
    x1 = gpos.x() + cos(PI*gdir.y()/180.)*(pos.x()-wid/2) + sin(PI*gdir.y()/180.)*(pos.z()-th/2);
    x2 = gpos.x() + cos(PI*gdir.y()/180.)*(pos.x()+wid/2) + sin(PI*gdir.y()/180.)*(pos.z()+th/2);
    z1 = gpos.z() - sin(PI*gdir.y()/180.)*(pos.x()-wid/2) + cos(PI*gdir.y()/180.)*(pos.z()-th/2);
    z2 = gpos.z() - sin(PI*gdir.y()/180.)*(pos.x()+wid/2) + cos(PI*gdir.y()/180.)*(pos.z()+th/2);
    box.DrawBox(z1,x1,z2,x2);
    if( seg==1 ){
      PC_x1 = x1;
      PC_z1 = z1;
    }
    if( seg==27 ){
      PC_x2 = x1;
      PC_z2 = z1;
    }
  }
  TLine line;
  line.SetLineColor(4);
  PC_A = (PC_x1-PC_x2)/(PC_z1-PC_z2);
  PC_B = PC_x1-PC_A*PC_z1;
  double PC_xx1= -100., PC_xx2 = -500;
  double PC_zz1 = (PC_xx1-PC_B)/PC_A;
  double PC_zz2 = (PC_xx2-PC_B)/PC_A;
  line.DrawLine(PC_zz1, PC_xx1, PC_zz2, PC_xx2);

  // for FDC1
  confMan-> GetBLDCWireMapManager()-> GetGParam( CID_FDC1, gpos, gdir);
  x1 = gpos.x() - FDC1_box_x/2;
  z1 = gpos.z() - FDC1_box_z/2;
  x2 = gpos.x() + FDC1_box_x/2;
  z2 = gpos.z() + FDC1_box_z/2;
  box.SetLineColor(14);
  box.DrawBox(z1,x1,z2,x2);
}

void EventAnalysisCheck::Finalize()
{
  std::cout << " Enter EventAnalysisCheck::Finalize " << std::endl;

  delete data;

  delete blMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisCheck *event = new EventAnalysisCheck();
  return (EventTemp*)event;
}
