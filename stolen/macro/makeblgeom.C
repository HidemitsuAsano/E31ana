#include <iomanip>
void makeblgeom(){
  // BHD
  int cid = 2;
  int nseg = 20;
  double len = 10.;
  double wid = 2.;
  double thi = 1.;
  out(cid,nseg,len,wid,thi,"BHD");

  // PA
  cid = 3;
  nseg = 8;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"PA");

  // T0
  cid = 4;
  nseg = 5;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"T0");

  // E0
  cid = 5;
  nseg = 3;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"E0");

  // B1
  cid = 6;
  nseg = 1;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"B1");

  // B2
  cid = 13;
  nseg = 1;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"B2");

  // LC1
  cid = 7;
  nseg = 1;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"LC1");

  // LC2
  cid = 8;
  nseg = 2;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"LC2");

  // AC
  cid = 9;
  nseg = 1;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"AC");

  // WC
  cid = 10;
  nseg = 1;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"WC");

  // GC
  cid = 11;
  nseg = 1;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"GC");

  // Range
  cid = 12;
  nseg = 10;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"Range");

  // TOF
  cid = 14;
  nseg = 32;
  len = 10.;
  wid = 2.;
  thi = 1.;
  out(cid,nseg,len,wid,thi,"TOFstop");



}

void out( int cid, int nseg, double len, double wid, double thi, TString name="AAA" ){
  outheader(cid,name);

  for( int seg=1; seg<=nseg; seg++ ){
    double x = 0.;
    double y = 0.;
    double z = 0.;
    double dx = 0.;
    double dy = 0.;
    double dz = 0;
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
	      << std::setw(12) << thi
	      << std::endl;
  }

}

void outheader(int cid, TString name){
  std::cout << "##" << std::endl  << "#" << std::endl<< "# " << name << std::endl
	    << "GPOS:" << std::setw(4) << cid
	    << setiosflags(ios::showpoint) << "  " << 0.0 << "  " << 0.0 << "  " << 0.0 << "  " << 0.0 << "  " << 0.0 << "  " << 0.0 << std::endl
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
