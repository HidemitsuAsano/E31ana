#include "UshiwakaTableMan.h"

UshiwakaTable::UshiwakaTable() : X(-999), Y(-999), dX(-999), dY(-999), Mom(0), Ang(-999), FL(-999)
{
}

UshiwakaTable::UshiwakaTable(const double &x, const double &y, const double &dx, const double &dy, const double &mom, const double &ang, const double &fl)
  : X(x), Y(y), dX(dx), dY(dy), Mom(mom), Ang(ang), FL(fl)
{
}

UshiwakaTableMan::UshiwakaTableMan()
  : FileName(DefaultTableName), DT(-999), Zmax(-999), Zmin(-999),
    NX(-1), Xmax(-999), Xmin(-999), dX(-999), NDX(-1), DXmax(-999), DXmin(-999), dDX(-999),
    NY(-1), Ymax(-999), Ymin(-999), dY(-999), NDY(-1), DYmax(-999), DYmin(-999), dDY(-999),
    NMom(-1), MomMin(-999), MomMax(-999), dMom(-999)
{
}

UshiwakaTableMan::UshiwakaTableMan(const std::string &filename)
  : FileName(filename), DT(-999), Zmax(-999), Zmin(-999),
    NX(-1), Xmax(-999), Xmin(-999), dX(-999), NDX(-1), DXmax(-999), DXmin(-999), dDX(-999),
    NY(-1), Ymax(-999), Ymin(-999), dY(-999), NDY(-1), DYmax(-999), DYmin(-999), dDY(-999),
    NMom(-1), MomMin(-999), MomMax(-999), dMom(-999)
{
}

bool UshiwakaTableMan::Initialize()
{
  static const std::string funcname = "UshiwakaTableMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ..."<< std::endl;

  std::ifstream ifs(FileName.c_str());

  std::string str;
  while( !ifs.eof() ){
    getline(ifs, str);
    if( str.find("#")!=std::string::npos ) continue;
    int n;
    double max, min, dpara;
    double dt, zmax, zmin;

    if( sscanf(str.c_str(), "dt: %lf range: %lf %lf", &dt, &zmin, &zmax)==3 ){
      DT=dt, Zmax=zmax, Zmin=zmin;
      continue;
    }

    if( sscanf(str.c_str(), "X: %d %lf  %lf %lf", &n, &min, &max, &dpara)==4 ){
      NX=n, Xmax=max, Xmin=min, dX=dpara;
      continue;
    }
    if( sscanf(str.c_str(), "dX: %d %lf  %lf %lf", &n, &min, &max, &dpara)==4 ){
      NDX=n, DXmax=max, DXmin=min, dDX=dpara;
      continue;
    }
    if( sscanf(str.c_str(), "Y: %d %lf  %lf %lf", &n, &min, &max, &dpara)==4 ){
      NY=n, Ymax=max, Ymin=min, dY=dpara;
      continue;
    }
    if( sscanf(str.c_str(), "dY: %d %lf  %lf %lf", &n, &min, &max, &dpara)==4 ){
      NDY=n, DYmax=max, DYmin=min, dDY=dpara;
      continue;
    }
    if( sscanf(str.c_str(), "Mom: %d %lf  %lf %lf", &n, &min, &max, &dpara)==4 ){
      NMom=n, MomMax=max, MomMin=min, dMom=dpara;
      continue;
    }
    if( str.find(">")==0 ) break;

    std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
    std::cerr << std::string(str) << std::endl;
  }

  if( NX<0 ){
    std::cerr<<" Error Num of X parameter is "<<NX<<std::endl;
    return false;
  }
  if( NY<0 ){
    std::cerr<<" Error Num of Y parameter is "<<NY<<std::endl;
    return false;
  }
  if( NDX<0 ){
    std::cerr<<" Error Num of dX parameter is "<<NDX<<std::endl;
    return false;
  }
  if( NDY<0 ){
    std::cerr<<" Error Num of dY parameter is "<<NDY<<std::endl;
    return false;
  }
  if( NMom<0 ){
    std::cerr<<" Error Num of momentum parameter is "<<NMom<<std::endl;
    return false;
  }

  unsigned long key = 0;
  while(!ifs.eof()){
    getline(ifs, str);
    if( str.find("#")!=std::string::npos ) continue;
    double x, y, dx, dy, mom, ang, fl;
    if( sscanf(str.c_str(), " %lf %lf %lf %lf %lf %lf %lf", &x, &y, &dx, &dy, &ang, &mom, &fl)==7 ){
      UshiwakaTable table(x, y, dx, dy, mom, ang, fl);
      tableMap[key] = table;
      key++;
      continue;
    }
  }
  return true;
}

UshiwakaTable* UshiwakaTableMan::table(int ix, int iy, int idx, int idy, int iang)
{
  if( ix>=NX ){
    std::cerr<<" Error ix > Num Of X parameter"<<std::endl;
    return 0;
  }
  if( iy>=NY ){
    std::cerr<<" Error iy > Num Of Y parameter"<<std::endl;
    return 0;
  }
  if( idx>=NDX ){
    std::cerr<<" Error idx > Num Of dX parameter"<<std::endl;
    return 0;
  }
  if( idy>=NDY ){
    std::cerr<<" Error idy > Num Of dY parameter"<<std::endl;
    return 0;
  }
  if( iang>=NMom ){
    std::cerr<<" Error iang > Num Of angle parameter"<<std::endl;
    return 0;
  }
  if( ix<0 ) ix = (NX-1)/2;
  if( iy<0 ) iy = (NY-1)/2;
  if( idx<0 ) idx = (NDX-1)/2;
  if( idy<0 ) idy = (NDY-1)/2;
  if( iang<0 ) iang = (NMom-1)/2;

  unsigned long key = iang + NMom*idy + NMom*NDY*idx + NMom*NDY*NDX*iy + NMom*NDY*NDX*NY*ix;
  return &tableMap[key];
}

bool UshiwakaTableMan::GetParam(const double &x, const double &y, const double &dx, const double &dy, const double &ang, double &mom, double &fl)
{
  if( x<Xmin || Xmax<x ){
    //    std::cerr<<" UshiwakaTableMan::GetParam Error x is out of range ! "<<x<<std::endl;
    return false;
  }
  if( y<Ymin || Ymax<y ){
    //    std::cerr<<" UshiwakaTableMan::GetParam Error y is out of range ! "<<y<<std::endl;
    return false;
  }
  if( dx<DXmin || DXmax<dx ){
    //    std::cerr<<" UshiwakaTableMan::GetParam Error dx is out of range ! "<<dx<<std::endl;
    return false;
  }
  if( dy<DYmin || DYmax<dy ){
    //    std::cerr<<" UshiwakaTableMan::GetParam Error dy is out of range ! "<<dy<<std::endl;
    return false;
  }
  double int_x, int_y, int_dx, int_dy;
  double ratio_x = modf((x-Xmin)/dX, &int_x);
  double ratio_y = modf((y-Ymin)/dY, &int_y);
  double ratio_dx = modf((dx-DXmin)/dDX, &int_dx);
  double ratio_dy = modf((dy-DYmin)/dDY, &int_dy);

  int ix = (int)int_x;
  int iy = (int)int_y;
  int idx = (int)int_dx;
  int idy = (int)int_dy;

  int iang=-1;
  for( int i=0; i<NMom; i++ ){
    double get_ang = table(ix,iy,idx,idy,i)-> ang();
    if( get_ang-ang<0 ){
      iang = i-1;
      break;
    }
  }
  if( iang==-1 ){
    //    std::cerr<<" UshiwakaTableMan::GetParam Error angle is out of range ! "<<ang<<std::endl;
    return false;
  }

  double get_ang1 = table(ix,iy,idx,idy,iang)-> ang();
  double get_ang2 = table(ix,iy,idx,idy,iang+1)-> ang();
  double ratio_ang = (ang-get_ang1)/(get_ang2-get_ang1);

  //get_xxx[ix][iy][idx][idy][iang]
  double get_mom[2][2][2][2][2];
  double get_fl[2][2][2][2][2];

  for( int cou_ix=0; cou_ix<2; cou_ix++ ){
    for( int cou_iy=0; cou_iy<2; cou_iy++ ){
      for( int cou_idx=0; cou_idx<2; cou_idx++ ){
	for( int cou_idy=0; cou_idy<2; cou_idy++ ){
	  for( int cou_iang=0; cou_iang<2; cou_iang++ ){
	    get_mom[cou_ix][cou_iy][cou_idx][cou_idy][cou_iang] = table(ix+cou_ix,iy+cou_iy,idx+cou_idx,idy+cou_idy,iang+cou_iang)-> mom();
	    get_fl[cou_ix][cou_iy][cou_idx][cou_idy][cou_iang] = table(ix+cou_ix,iy+cou_iy,idx+cou_idx,idy+cou_idy,iang+cou_iang)-> fl();
	  }
	}
      }
    }
  }

  double mom_dy[2][2][2][2];
  double fl_dy[2][2][2][2];
  for( int cou_x=0; cou_x<2; cou_x++ ){
    for( int cou_y=0; cou_y<2; cou_y++ ){
      for( int cou_dx=0; cou_dx<2; cou_dx++ ){
	for( int cou_ang=0; cou_ang<2; cou_ang++ ){
	  mom_dy[cou_x][cou_y][cou_dx][cou_ang] = ratio_dy*get_mom[cou_x][cou_y][cou_dx][1][cou_ang] + (1.0-ratio_dy)*get_mom[cou_x][cou_y][cou_dx][0][cou_ang];
	  fl_dy[cou_x][cou_y][cou_dx][cou_ang] = ratio_dy*get_fl[cou_x][cou_y][cou_dx][1][cou_ang] + (1.0-ratio_dy)*get_fl[cou_x][cou_y][cou_dx][0][cou_ang];
	}
      }
    }
  }

  double mom_y[2][2][2];
  double fl_y[2][2][2];
  for( int cou_x=0; cou_x<2; cou_x++ ){
    for( int cou_dx=0; cou_dx<2; cou_dx++ ){
      for( int cou_ang=0; cou_ang<2; cou_ang++ ){
	mom_y[cou_x][cou_dx][cou_ang] = ratio_y*mom_dy[cou_x][1][cou_dx][cou_ang] + (1.0-ratio_y)*mom_dy[cou_x][0][cou_dx][cou_ang];
	fl_y[cou_x][cou_dx][cou_ang] = ratio_y*fl_dy[cou_x][1][cou_dx][cou_ang] + (1.0-ratio_y)*fl_dy[cou_x][0][cou_dx][cou_ang];
      }
    }
  }

  double mom_dx[2][2];
  double fl_dx[2][2];
  for( int cou_x=0; cou_x<2; cou_x++ ){
    for( int cou_ang=0; cou_ang<2; cou_ang++ ){
      mom_dx[cou_x][cou_ang] = ratio_dx*mom_y[cou_x][1][cou_ang] + (1.0-ratio_dx)*mom_y[cou_x][0][cou_ang];
      fl_dx[cou_x][cou_ang] = ratio_dx*fl_y[cou_x][1][cou_ang] + (1.0-ratio_dx)*fl_y[cou_x][0][cou_ang];
    }
  }

  double mom_x[2];
  double fl_x[2];
  for( int cou_ang=0; cou_ang<2; cou_ang++ ){
    mom_x[cou_ang] = ratio_x*mom_dx[1][cou_ang] + (1.0-ratio_x)*mom_dx[0][cou_ang];
    fl_x[cou_ang] = ratio_x*fl_dx[1][cou_ang] + (1.0-ratio_x)*fl_dx[0][cou_ang];
  }

  mom = ratio_ang*mom_x[1] + (1.0-ratio_ang)*mom_x[0];
  fl = ratio_ang*fl_x[1] + (1.0-ratio_ang)*fl_x[0];
  //  std::cout<<"UshiwakaTableMan::GetParam("<<x<<","<<y<<","<<dx<<","<<dy<<","<<ang<<","<<mom<<","<<fl<<")"<<std::endl;
  return true;
}

bool UshiwakaTableMan::DumpStatus()
{
  std::cout<<"###"<<std::endl;
  std::cout<<"#   Ushiwaka Table"<<std::endl;
  std::cout<<"#"<<std::endl;
  std::cout<<"#    [ns]      [cm]   [cm]"<<std::endl;
  std::cout<<"dt: "<<DT<<" range:   "<<Zmin<<"  "<<Zmax<<std::endl;
  std::cout<<"#"<<std::endl;
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout<<"#   nx    x_min      x_max      dx"<<std::endl;
  std::cout<<"X:   "<<NX<<"  "<<Xmin<<"  "<<Xmax<<"  "<<dX<<std::endl;
  std::cout<<"#   ny    y_min      y_max      dy"<<std::endl;
  std::cout<<"Y:   "<<NY<<"  "<<Ymin<<"  "<<Ymax<<"  "<<dY<<std::endl;
  std::cout<<"#  ndx   dx_min      dx_max     d(dx)"<<std::endl;
  std::cout<<"dX:  "<<NDX<<"  "<<DXmin<<"  "<<DXmax<<"  "<<dDX<<std::endl;
  std::cout<<"#  ndy   dy_min      dy_max     d(dy)"<<std::endl;
  std::cout<<"dY:  "<<NDY<<"  "<<DYmin<<"  "<<DYmax<<"  "<<dDY<<std::endl;
  std::cout<<"# n(mom) mom_min    mom_max     d(mom)"<<std::endl;
  std::cout<<"Mom: "<<NMom<<"   "<<MomMin<<"  "<<MomMax<<"  "<<dMom<<std::endl;
  std::cout<<"#"<<std::endl;
  std::cout<<"###"<<std::endl;

  unsigned long n_para = NX*NY*NDX*NDY*NMom;
  if( n_para!=tableMap.size() ){
    std::cerr<<" not macth Data Size !!!"<<std::endl;
    return false;
  }
  else{
    std::cout<<n_para<<" parameters"<<std::endl;
  }
  return true;
}
