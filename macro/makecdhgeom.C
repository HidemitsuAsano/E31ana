#include <iomanip>
void makecdhgeom(){
  int cid = 1;
  int nseg = 36;
  double rad = 55.;
  double phi0 = 5.;
  double dphi = 10.;
  double len = 70.;
  double wid = 9.5;
  double thick=3.0;
  
//   std::cout << "##" << std::endl  << "#" << std::endl<< "#" << std::endl
// 	    << "GPOS: " << cid << " 0 0 0 0 0 0" << std::endl
// 	    << "#"
// 	    << std::setw(11) << "CID"
// 	    << std::setw(12) << "Seg"
// 	    << std::setw(12) << "X"
// 	    << std::setw(12) << "Y"
// 	    << std::setw(12) << "Z"
// 	    << std::setw(12) << "dX"
// 	    << std::setw(12) << "dY"
// 	    << std::setw(12) << "dZ"
// 	    << std::setw(12) << "Length"
// 	    << std::setw(12) << "Width"
// 	    << std::setw(12) << "Thick"
// 	    << std::endl
// 	    << "#"
// 	    << std::setw(11) << " "
// 	    << std::setw(12) << " "
// 	    << std::setw(12) << "[cm]"
// 	    << std::setw(12) << "[cm]"
// 	    << std::setw(12) << "[cm]"
// 	    << std::setw(12) << " "
// 	    << std::setw(12) << " "
// 	    << std::setw(12) << " "
// 	    << std::setw(12) << "[cm]"
// 	    << std::setw(12) << "[cm]"
// 	    << std::setw(12) << "[cm]"
// 	    << std::endl;
  outheader(cid);
  
  for( int seg=1; seg<=nseg; seg++ ){
    double x = rad*cos( (phi0+dphi*(seg-1))*TMath::DegToRad() );
    double y = rad*sin( (phi0+dphi*(seg-1))*TMath::DegToRad() );
    double z = 0.;
    double dx = cos( (phi0+dphi*(seg-1))*TMath::DegToRad() );
    double dy = sin( (phi0+dphi*(seg-1))*TMath::DegToRad() );
    double dz =0;
    std::cout << std::setprecision(6) << setiosflags(ios::showpoint) 
	      << std::setw(12) << cid
	      << std::setw(12) << seg 
	      << std::setw(12) << x 
	      << std::setw(12) << y
	      << std::setw(12) << z
	      << std::setw(12) << dx
	      << std::setw(12) << dy
	      << std::setw(12) << dz
	      << std::setw(12) << len
	      << std::setw(12) << wid
	      << std::setw(12) << thick
	      << std::endl;
  }
  
}

void outheader(int cid){
  std::cout << "##" << std::endl  << "#" << std::endl<< "#" << std::endl
	    << "GPOS: " << cid << " 0.0 0.0 0.0 0.0 0.0 0.0" << std::endl
	    << "#"
	    << std::setw(11) << "CID"
	    << std::setw(12) << "Seg"
	    << std::setw(12) << "X"
	    << std::setw(12) << "Y"
	    << std::setw(12) << "Z"
	    << std::setw(12) << "dX"
	    << std::setw(12) << "dY"
	    << std::setw(12) << "dZ"
	    << std::setw(12) << "Length"
	    << std::setw(12) << "Width"
	    << std::setw(12) << "Thick"
	    << std::endl
	    << "#"
	    << std::setw(11) << " "
	    << std::setw(12) << " "
	    << std::setw(12) << "[cm]"
	    << std::setw(12) << "[cm]"
	    << std::setw(12) << "[cm]"
	    << std::setw(12) << " "
	    << std::setw(12) << " "
	    << std::setw(12) << " "
	    << std::setw(12) << "[cm]"
	    << std::setw(12) << "[cm]"
	    << std::setw(12) << "[cm]"
	    << std::endl;
}
