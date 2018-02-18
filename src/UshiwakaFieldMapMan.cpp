//UshiwakaFieldMapMan.cpp
#include "UshiwakaFieldMapMan.h"

ClassImp(UshiwakaFieldMapMan);

UshiwakaFieldMapMan::UshiwakaFieldMapMan() : FileName(DefaultFieldMapName), DT(-999), SR(0), 
					     GX(0), GY(0), GZ(0), X_MAX(-999), Y_MAX(-999), Z_MAX(-999),
					     DX(-999), DY(-999), DZ(-999), NX(-1), NY(-1), NZ(-1)
{
}

UshiwakaFieldMapMan::UshiwakaFieldMapMan(const std::string &filename) : FileName(filename), DT(-999), SR(0), 
									GX(0), GY(0), GZ(0), X_MAX(-999), Y_MAX(-999), Z_MAX(-999),
									DX(-999), DY(-999), DZ(-999), NX(-1), NY(-1), NZ(-1)
{
}

UshiwakaFieldMapMan::~UshiwakaFieldMapMan()
{
}

const int MAXCHAR = 144;

bool UshiwakaFieldMapMan::Initialize()
{
  static const std::string funcname = "UshiwakaFieldMapMan::Initialize";
  std::cout<<"["<<funcname<<"] Initialization start ..."<<std::endl;

  if( FileName != DefaultFileName ){
    //int nd;
    double dt, sr;
    double x, y, z;
    double Bx, By, Bz;
    FILE *fp;
    char str[MAXCHAR];

    if( (fp=fopen(FileName.c_str(),"r"))==0 ){
      std::cerr << " File open fail. [" << FileName << "]" << std::endl;
      exit(-1);
    }

    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      if( sscanf(str, "dt: %lf",&dt)==1 ){
	DT = dt;
      }
      if( sscanf(str, "sr: %lf",&sr)==1 ){
	SR = (int)sr;
      }
      if( sscanf(str, "GPOS: %lf %lf %lf",&x,&y,&z)==3 ){
	GX = x, GY=y, GZ = z;
      }
      if( sscanf(str, "MAX: %lf %lf %lf",&x,&y,&z)==3 ){
        X_MAX = x, Y_MAX = y, Z_MAX = z;
      }
      if( sscanf(str, "dl: %lf %lf %lf",&x,&y,&z)==3 ){
        DX = x, DY = y, DZ = z;
      }
      if( strncmp(str,"Magnetic",8)==0 ){
	if( X_MAX>0 ){
	  NX = (int)(2*X_MAX/DX);
	  NY = (int)(2*Y_MAX/DY);
	  NZ = (int)(2*Z_MAX/DZ);
	}
	else return false;
      }
      if( sscanf(str, " %lf %lf %lf %lf %lf %lf",&x, &y, &z, &Bx, &By, &Bz)==6 ){
	TVector3 field(Bx, By, Bz);
	int index_x = (int)((x+X_MAX)/DX);
	int index_y = (int)((y+Y_MAX)/DY);
	int index_z = (int)((z+Z_MAX)/DZ);
	unsigned long key = index_x + NX*index_y + NX*NY*index_z;
	//	std::cout<<"  key:"<<key<<"  pos("<<x<<","<<y<<","<<z<<") field("<<Bx<<","<<By<<","<<Bz<<")"<<std::endl;
	fieldMap[key] = field;
      }
    }
  }
  return true;
}

TVector3 UshiwakaFieldMapMan::gvalue(const double gx, const double gy, const double gz)
{
  double x = gx-GX;
  double y = gy-GY;
  double z = gz-GZ;
  return value(x,y,z);
}

TVector3 UshiwakaFieldMapMan::value(const double x, const double y, const double z)
{
  TVector3 field(0.0, 0.0, 0.0);

  double int_x, int_y, int_z;
  double rate_x = modf((x+X_MAX)/DX, &int_x);
  double rate_y = modf((y+Y_MAX)/DY, &int_y);
  double rate_z = modf((z+Z_MAX)/DZ, &int_z);
  int ix = (int)int_x;
  int iy = (int)int_y;
  int iz = (int)int_z;

  if( ix<0 || NX-1<=ix ) return field;
  if( iy<0 || NY-1<=iy ) return field;
  if( iz<0 || NZ-1<=iz ) return field;

  unsigned long key1 = ix   + NX*iy     + NX*NY*iz;
  unsigned long key2 = ix+1 + NX*iy     + NX*NY*iz;
  unsigned long key3 = ix+1 + NX*(iy+1) + NX*NY*iz;
  unsigned long key4 = ix   + NX*(iy+1) + NX*NY*iz;
  unsigned long key5 = ix   + NX*iy     + NX*NY*(iz+1);
  unsigned long key6 = ix+1 + NX*iy     + NX*NY*(iz+1);
  unsigned long key7 = ix+1 + NX*(iy+1) + NX*NY*(iz+1);
  unsigned long key8 = ix   + NX*(iy+1) + NX*NY*(iz+1);

  TVector3 field1 = fieldMap[key1];
  TVector3 field2 = fieldMap[key2];
  TVector3 field3 = fieldMap[key3];
  TVector3 field4 = fieldMap[key4];
  TVector3 field5 = fieldMap[key5];
  TVector3 field6 = fieldMap[key6];
  TVector3 field7 = fieldMap[key7];
  TVector3 field8 = fieldMap[key8];

  TVector3 ave_x1 = rate_x*field1+(1-rate_x)*field2;
  TVector3 ave_x2 = rate_x*field4+(1-rate_x)*field3;
  TVector3 ave_x3 = rate_x*field5+(1-rate_x)*field6;
  TVector3 ave_x4 = rate_x*field8+(1-rate_x)*field7;

  TVector3 ave_y1 = rate_y*ave_x1+(1-rate_y)*ave_x2;
  TVector3 ave_y2 = rate_y*ave_x3+(1-rate_y)*ave_x4;

  field = rate_z*ave_y1+(1-rate_z)*ave_y2;

  return 0.0001*field;
}

bool UshiwakaFieldMapMan::RungeKutta(const char* name, TVector3 &inpos, TVector3 &inmom)
{
  std::string particlename(name);
  return RungeKutta(particlename, inpos, inmom);
}

bool UshiwakaFieldMapMan::RungeKutta(std::string &name, TVector3 &inpos, TVector3 &inmom)
{
  USWKTrack track;
  double time = 0.0;

  track.SetParticleName(name);
  track.SetInTime(time);
  track.SetInPosition(inpos);
  track.SetInMomentum(inmom);

  double mass = track.mass();

  TVector3 pos = inpos;
  TVector3 mom= inmom;
  int counter=0;
  while(InField(pos)){
    TVector3 v0 = (100.*Const/sqrt(mom.Mag2()+mass*mass))*mom;
    TVector3 kmom1 = 0.01*Const*v0.Cross(gvalue(pos.X(),pos.Y(),pos.Z()));
    TVector3 kpos1 = v0;

    TVector3 mom1 = mom+0.5*DT*kmom1;
    TVector3 pos1 = pos+0.5*DT*kpos1;
    TVector3 v1 = (100.*Const/sqrt(mom1.Mag2()+mass*mass))*mom1;
    TVector3 kmom2 = 0.01*Const*v1.Cross(gvalue(pos1.X(),pos1.Y(),pos1.Z()));
    TVector3 kpos2 = v1;

    TVector3 mom2 = mom+0.5*DT*kmom2;
    TVector3 pos2 = pos+0.5*DT*kpos2;
    TVector3 v2 = (100.*Const/sqrt(mom2.Mag2()+mass*mass))*mom2;
    TVector3 kmom3 = 0.01*Const*v2.Cross(gvalue(pos2.X(),pos2.Y(),pos2.Z()));
    TVector3 kpos3 = v2;

    TVector3 mom3 = mom+DT*kmom3;
    TVector3 pos3 = pos+DT*kpos3;
    TVector3 v3 = (100.*Const/sqrt(mom3.Mag2()+mass*mass))*mom3;
    TVector3 kmom4 = 0.01*Const*v3.Cross(gvalue(pos3.X(),pos3.Y(),pos3.Z()));
    TVector3 kpos4 = v3;

    pos += (kpos1+2.*kpos2+2.*kpos3+kpos4)*(DT/6.);
    mom += (kmom1+2.*kmom2+2.*kmom3+kmom4)*(DT/6.);

    counter++;
    time += DT;
    if( SR>0 ){
      if( counter%SR==0 ){
	track.SetTime(time);
	track.SetPosition(pos);
	track.SetMomentum(mom);
      }
    }
  }
  track.SetOutTime(time);
  track.SetOutPosition(pos);
  track.SetOutMomentum(mom);

  Tracks.push_back(track);

  return true;
}

bool UshiwakaFieldMapMan::InField(const TVector3 &pos)
{
  if( fabs(pos.X()-GX)>X_MAX+0.001 ) return false;
  if( fabs(pos.Y()-GY)>Y_MAX+0.001 ) return false;
  if( fabs(pos.Z()-GZ)>Z_MAX+0.001 ) return false;

  return true;
}

void UshiwakaFieldMapMan::Clear()
{
  Tracks.clear();
}

void UshiwakaFieldMapMan::DumpStatus()
{
  std::cout<<"===== Ushiwaka FieldMap Status ====="<<std::endl;
  std::cout<<"  === dt : "<<DT<<"[ns] Sample rate : 1/"<<SR<<std::endl;
  std::cout<<"  === GPOS  : "<<GX<<", "<<GY<<", "<<GZ<<"[cm] "<<std::endl;
  std::cout<<"  === Range : "<<X_MAX<<", "<<Y_MAX<<", "<<Z_MAX<<"[cm] "<<std::endl;
  std::cout<<"  === DL    : "<<DX<<", "<<DY<<", "<<DZ<<"[cm] "<<std::endl;
  std::cout<<"  === N Data: "<<NX<<", "<<NY<<", "<<NZ<<std::endl;

  std::cout<<"All Data ="<< NX*NY*NZ<<"  field Container size = "<<fieldMap.size()<<std::endl;
}

void UshiwakaFieldMapMan::DrawXGraph(TVirtualPad *pad, const int direction, const int iy, const int iz)
{
  pad-> cd();
  double pos[NX];
  double data[NX];

  double y = -Y_MAX + DY*iy;
  double z = -Z_MAX + DZ*iz;
  std::cout<<"=== Draw Ushiwaka Magnetic Field X distribution (y,z) = ("<<y<<","<<z<<") ";
  if( direction==0 ) std::cout<<"B_x ==="<<std::endl;
  else if( direction==1 ) std::cout<<"B_y ==="<<std::endl;
  else if( direction==2 ) std::cout<<"B_z ==="<<std::endl;

  for( int ix=0; ix<NX; ix++ ){
    pos[ix] = -X_MAX+ix*DX;
    unsigned long key = ix + NX*iy + NX*NY*iz;
    TVector3 field = fieldMap[key];
    if( direction==0 ) data[ix] = field.X();
    else if( direction==1 ) data[ix] = field.Y();
    else if( direction==2 ) data[ix] = field.Z();
  }

  TGraph *gra = new TGraph(NX, pos, data);
  gra-> Draw("APL");
  delete gra;
}

void UshiwakaFieldMapMan::DrawYGraph(TVirtualPad *pad, const int direction, const int ix, const int iz)
{
  pad-> cd();
  double pos[NY];
  double data[NY];

  double x = -X_MAX + DX*ix;
  double z = -Z_MAX + DZ*iz;
  std::cout<<"=== Draw Ushiwaka Magnetic Field X distribution (x,z) = ("<<x<<","<<z<<") ";
  if( direction==0 ) std::cout<<"B_x ==="<<std::endl;
  else if( direction==1 ) std::cout<<"B_y ==="<<std::endl;
  else if( direction==2 ) std::cout<<"B_z ==="<<std::endl;

  for( int iy=0; iy<NY; iy++ ){
    pos[iy] = -Y_MAX+iy*DY;
    unsigned long key = ix + NX*iy + NX*NY*iz;
    TVector3 field = fieldMap[key];
    if( direction==0 ) data[iy] = field.X();
    else if( direction==1 ) data[iy] = field.Y();
    else if( direction==2 ) data[iy] = field.Z();
  }

  TGraph *gra = new TGraph(NY, pos, data);
  gra-> Draw("APL");
  delete gra;
}

void UshiwakaFieldMapMan::DrawZGraph(TVirtualPad *pad, const int direction, const int ix, const int iy)
{
  pad-> cd();
  double pos[NZ];
  double data[NZ];

  double x = -X_MAX + DX*ix;
  double y = -Y_MAX + DY*iy;
  std::cout<<"=== Draw Ushiwaka Magnetic Field Z distribution (x,y) = ("<<x<<","<<y<<") ";
  if( direction==0 ) std::cout<<"B_x ==="<<std::endl;
  else if( direction==1 ) std::cout<<"B_y ==="<<std::endl;
  else if( direction==2 ) std::cout<<"B_z ==="<<std::endl;

  for( int iz=0; iz<NZ; iz++ ){
    pos[iz] = -Z_MAX+iz*DZ;
    unsigned long key = ix + NX*iy + NX*NY*iz;
    TVector3 field = fieldMap[key];
    if( direction==0 ) data[iz] = field.X();
    else if( direction==1 ) data[iz] = field.Y();
    else if( direction==2 ) data[iz] = field.Z();
  }

  TGraph *gra = new TGraph(NZ, pos, data);
  gra-> Draw("APL");
  delete gra;
}

bool UshiwakaFieldMapMan::DrawZGraph(TVirtualPad *pad, const int direction, const double x, const double y, const int n)
{
  pad-> cd();
  double pos[n];
  double data[n];

  if( fabs(x)>X_MAX ) return false;
  if( fabs(y)>Y_MAX ) return false;

  std::cout<<"=== Draw Ushiwaka Magnetic Field Z distribution (x,y) = ("<<x<<","<<y<<") ";
  if( direction==0 ) std::cout<<"B_x ==="<<std::endl;
  else if( direction==1 ) std::cout<<"B_y ==="<<std::endl;
  else if( direction==2 ) std::cout<<"B_z ==="<<std::endl;

  double dz = 2*Z_MAX/n;

  for( int iz=0; iz<n; iz++ ){
    pos[iz] = -Z_MAX+iz*dz;

    TVector3 field = value(x,y,pos[iz]);
    if( direction==0 ) data[iz] = field.X();
    else if( direction==1 ) data[iz] = field.Y();
    else if( direction==2 ) data[iz] = field.Z();
  }

  TGraph *gra = new TGraph(n, pos, data);
  gra-> Draw("APL");
  delete gra;
  return 0;
}
