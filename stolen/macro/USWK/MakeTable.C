#include "src/GlobalVariables.h"

static const double PI = 6.*asin(0.5);

static const int NumOfX = 3;
static const double MaxX = 20;

static const int NumOfY = 3;
static const double MaxY = 20;

static const int NumOfDX = 3;
static const double MaxDX = 0.01;

static const int NumOfDY = 3;
static const double MaxDY = 0.01;

static const int NumOfMom = 3;
static const double MinMom = 1.0;
static const double MaxMom = 2.0;

void MakeTable(char *table_name = "tmpUshiwakaTable.param")
{
  gSystem-> Load("src/lib/libAll.so");

  ConfMan *conf = new ConfMan("conf/Run47/analyzer.conf");
  conf-> Initialize();
  long n_calc = NumOfX*NumOfY*NumOfDX*NumOfDY*NumOfMom;
  long i_calc = 0;

  UshiwakaFieldMapMan *fieldMapMan = new UshiwakaFieldMapMan();
  fieldMapMan-> SetFileName("param/Run47/inoue/UshiwakaFieldMap.param");

  fieldMapMan-> Initialize();
  fieldMapMan-> DumpStatus();

  double USWK_gpos[3];
  double USWK_range[3];

  fieldMapMan-> GetGPOS(USWK_gpos[0], USWK_gpos[1], USWK_gpos[2]);
  fieldMapMan-> GetRange(USWK_range[0], USWK_range[1], USWK_range[2]);
  double dt = fieldMapMan-> dt();

  const double dX = 2.*MaxX/(NumOfX-1);
  const double dDX = 2.*MaxDX/(NumOfDX-1);

  const double dY = 2.*MaxX/(NumOfY-1);
  const double dDY = 2.*MaxDX/(NumOfDY-1);

  const double dMom = (MaxMom-MinMom)/(NumOfMom-1);

  cout<<"##### Create Uhiwaka Table name : "<<table_name<<" #####"<<endl;
  cout<<"  === Num Of Calclation : "<<n_calc<<endl;
  ofstream ofs(table_name);
  ofs<<"###"<<endl;
  ofs<<"#   Ushiwaka Table"<<endl;
  ofs<<"#"<<endl;
  ofs<<"# use fieldMap : "<<fieldMapMan-> GetFileName()<<endl;
  ofs<<"#    [ns]      [cm]   [cm]"<<endl;       
  ofs<<"dt: "<<dt<<" range:   "<<USWK_gpos[2]-USWK_range[2]<<"  "<<USWK_gpos[2]+USWK_range[2]<<endl;
  ofs<<"#"<<endl;
  ofs.setf(ios::fixed, ios::floatfield);
  ofs<<"#   nx    x_max      x_min      dx"<<endl;
  ofs<<"X:   "<<NumOfX<<"  "<<-MaxX<<"  "<<MaxX<<"  "<<dX<<endl;
  ofs<<"#   ny    y_max      y_min      dy"<<endl;
  ofs<<"Y:   "<<NumOfY<<"  "<<-MaxY<<"  "<<MaxY<<"  "<<dY<<endl;
  ofs<<"#  ndx   dx_max      dx_min     d(dx)"<<endl;
  ofs<<"dX:  "<<NumOfDX<<"  "<<-MaxDX<<"  "<<MaxDX<<"  "<<dDX<<endl;
  ofs<<"#  ndy   dy_may      dy_min     d(dy)"<<endl;
  ofs<<"dY:  "<<NumOfDY<<"  "<<-MaxDY<<"  "<<MaxDY<<"  "<<dDY<<endl;
  ofs<<"# n(mom) mom_max    mom_min     d(mom)"<<endl;
  ofs<<"Mom: "<<NumOfMom<<"   "<<MaxMom<<"  "<<MinMom<<"  "<<dMom<<endl; 
  ofs<<"#"<<endl;
  ofs<<"###"<<endl;
  ofs<<">  pos_x ]cm]   pos_y [cm]    dx/dz        dy/dz         angle [rad]  mom [GeV/c]  flight length [cm]"<<endl;

  for( int i_x=0; i_x<NumOfX; i_x++ ){
    for( int i_y=0; i_y<NumOfY; i_y++ ){
      for( int i_dx=0; i_dx<NumOfDX; i_dx++ ){
	for( int i_dy=0; i_dy<NumOfDY; i_dy++ ){
	  for( int i_mom=0; i_mom<NumOfMom; i_mom++ ){
	    i_calc++;
	    if( i_calc%10000==0 ){
	      cout<<"    > "<<i_calc<<" finish"<<endl;
	    }

	    double in_x = -MaxX+i_x*dX;
	    double in_y = -MaxY+i_y*dY;
	    double in_z = USWK_gpos[2]-USWK_range[2];
	    TVector3 in_pos(in_x, in_y, in_z);

	    double dx = -MaxDX+i_dx*dDX;
	    double dy = -MaxDY+i_dy*dDY;
	    double mom = MinMom+i_mom*dMom;

	    TVector3 in_mom(0.0, 0.0, mom);
	    double ang_x = -atan(dy);
	    double ang_y = atan(dx);
	    in_mom.RotateX(ang_x);
	    in_mom.RotateY(ang_y);
	    double in_dxdz = in_mom.X()/in_mom.Z();
	    double in_dydz = in_mom.Y()/in_mom.Z();

	    fieldMapMan-> RungeKutta("proton", in_pos, in_mom);

	    USWKTrack *uswk_track = fieldMapMan-> GetUSWKTrack(0);
	    double out_time = uswk_track-> outtime();
	    TVector3 out_pos = uswk_track-> outposition();
	    TVector3 out_mom = uswk_track-> outmomentum();

	    double in_ang = atan2(in_mom.X(),in_mom.Z());
	    double out_ang = atan2(out_mom.X(),out_mom.Z());
	    double bending_ang = in_ang-out_ang;

	    if( out_pos.Z()<USWK_gpos[2] ){
	      cout<<"Error out pos"<<endl;
	    }
	    if( fabs(in_mom.Mag() - out_mom.Mag()) > 0.0001 ){
	      cout<<"Error not match in momentum & out momentum"<<fabs(in_mom.Mag()-out_mom.Mag())<<endl;
	    }

	    double out_dxdz = out_mom.X()/out_mom.Z();
	    double out_dydz = out_mom.Y()/out_mom.Z();
	    double out_z = USWK_gpos[2]+USWK_range[2];
	    double out_x = out_pos.X()+out_dxdz*(out_z-out_pos.Z());
	    double out_y = out_pos.Y()+out_dydz*(out_z-out_pos.Z());
	    double dx = out_x-out_pos.X();
	    double dy = out_y-out_pos.Y();
	    double dz = out_z-out_pos.Z();
	    double dl = sqrt(dx*dx+dy*dy+dz*dz);
	    double out_abs_mom = out_mom.Mag();
	    double out_v = 100.*Const*out_abs_mom/sqrt(pMass*pMass+out_abs_mom*out_abs_mom);
	    out_time -=dl/out_v;
	    double flight_length = out_v*out_time;

	    ofs<<setw(13)<<in_pos.X()<<setw(13)<<in_pos.Y()<<setw(13)<<in_dxdz<<setw(13)<<in_dydz<<setw(13)
	       <<bending_ang<<setw(13)<<mom<<setw(13)<<flight_length<<endl;

#if 0
	    cout.setf(ios::fixed, ios::floatfield);
	    cout<<setw(11)<<in_pos.X()<<setw(11)<<in_pos.Y()<<setw(11)<<in_mom.X()/in_mom.Z()<<setw(11)<<in_mom.Y()/in_mom.Z()
		<<setw(11)<<bending_ang<<setw(11)<<out_mom.Mag()<<setw(14)<<flight_length<<endl;
#endif
	    fieldMapMan-> Clear();
	  }
	}
      }
    }
  }
}
