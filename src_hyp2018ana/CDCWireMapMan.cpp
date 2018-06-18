// CDCWireMapMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>
#include <cctype>
#include <map>
#include <utility>

#include "MathTools.h"
#include "CDCWireMapMan.h"
#define OFFSET 1

ClassImp(CDCWireMap);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMap::CDCWireMap()
{
  Layer = Wire = -1;
  SuperLayer = ASDNum = TransType = ASDch = -1;
  Radius = Phi = Tilt = 0.0;
  for(int i=0;i<3;i++){
    Pos[i]= -999.;
    Posp[i]= -999.;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMap::CDCWireMap( int layer, int wire )
{
  Layer = layer; Wire = wire;
  SuperLayer = ASDNum = TransType = ASDch = -1;
  Radius = Phi = Tilt = 0.0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMap::CDCWireMap( int layer, int wire,
			int slayer, int asdnum, int ttype, int asdch,
			double rad, double phi, double tilt, double zlen )

{
  Layer = layer; Wire = wire;
  SuperLayer = slayer; ASDNum = asdnum; TransType = ttype; ASDch = asdch;
  Radius = rad; Phi = phi; Tilt = tilt;

  Pos[0]=Radius*cos(Phi*TMath::DegToRad());
  Pos[1]=Radius*sin(Phi*TMath::DegToRad());
  Pos[2]=zlen/2.;

#if OFFSET
  Posp[0]=Radius*cos((Phi+Tilt)*TMath::DegToRad());
  Posp[1]=Radius*sin((Phi+Tilt)*TMath::DegToRad());
#else
  double S = zlen*tan(Tilt*TMath::DegToRad());
  double rad2 = Radius*Radius;
  double cost = 1.-S*S/(2.*rad2);
  double sint = S/(2.*rad2)*sqrt(4.*rad2-S*S);
  if( 0<TMath::Abs(tilt) ){
    Posp[0]=Pos[0]*cost-Pos[1]*sint;
    Posp[1]=Pos[0]*sint+Pos[1]*cost;  
  }else{
    Posp[0]=Pos[0];
    Posp[1]=Pos[1];
  }
#endif
  Posp[2]=-zlen/2.;
#if 0
  if(wire==1){
    std::cout<<"layer "<<layer<<"    wire "<<wire;
    std::cout<<"    , wpos   "<< Pos[0]<<"  "<< Pos[1]<<"  "<< Pos[2]<<"  "<< pos().Phi()*TMath::RadToDeg();
    std::cout<<"    ,wposp  "<<Posp[0]<<"  "<<Posp[1]<<"  "<<Posp[2]<<"  "<< posp().Phi()*TMath::RadToDeg()<<std::endl;
  }
#endif
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //

ClassImp(RCDCWireMap);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
RCDCWireMap::RCDCWireMap()
{
  Crate = Slot = Channel = -1;
  SuperLayer = ASDNum = TransType = ASDch = -1;
  Radius = Phi = Tilt = 0.0;
  for(int i=0;i<3;i++){
    Pos[i]= -999.;
    Posp[i]= -999.;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
RCDCWireMap::RCDCWireMap( int cr, int sl, int ch )
{
  Crate = cr; Slot = sl; Channel = ch;
  SuperLayer = ASDNum = TransType = ASDch = -1;
  Radius = Phi = Tilt = 0.0;
  for(int i=0;i<3;i++){
    Pos[i]= -999.;
    Posp[i]= -999.;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
RCDCWireMap::RCDCWireMap( int cr, int sl, int ch,
			  int slayer, int asdnum, int ttype, int asdch,
			  double rad, double phi, double tilt,double zlen )
{
  Crate = cr; Slot = sl; Channel = ch;
  SuperLayer = slayer; ASDNum = asdnum; TransType = ttype; ASDch = asdch;
  Radius = rad; Phi = phi; Tilt = tilt;

  Pos[0]=Radius*cos(Phi*TMath::DegToRad());
  Pos[1]=Radius*sin(Phi*TMath::DegToRad());
  Pos[2]=zlen/2.;
#if OFFSET
  Posp[0]=Radius*cos((Phi+Tilt)*TMath::DegToRad());
  Posp[1]=Radius*sin((Phi+Tilt)*TMath::DegToRad());
#else
  double S = zlen*tan(Tilt*TMath::DegToRad());
  double rad2 = Radius*Radius;
  double cost = 1.-S*S/(2.*rad2);
  double sint = S/(2.*rad2)*sqrt(4.*rad2-S*S);
  if( 0<TMath::Abs(Tilt) ){
    Posp[0]=Pos[0]*cost-Pos[1]*sint;
    Posp[1]=Pos[0]*sint+Pos[1]*cost;  
  }else{
    Posp[0]=Pos[0];
    Posp[1]=Pos[1];
  }
#endif
  Posp[2]=-zlen/2.;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //

ClassImp( CDCWireMapMan )
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMapMan::CDCWireMapMan()
  : WireMapFileName()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMapMan::CDCWireMapMan( const std::string &filename1, const std::string &filename2, const std::string &filename3 )
  : ASDMapFileName(filename1), ChannelMapFileName(filename2), GeometryFileName(filename3)
{
  GX = GY = GZ = dGX = dGY = dGZ = -999.;
  for( int i=0; i<=NumOfCDCLayers; i++ ){
    Radius[i] = Phi0[i] = dPhi[i] = Tilt[i] = 0.;
    NWires[i] = 0;
  }
  ZLengthOfWire = InnerRadius = OuterRadius = 0.;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMapMan::~CDCWireMapMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CDCWireMapMan::Clear()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CDCWireMapMan::SetFileName( const std::string &filename1,
				 const std::string &filename2,
				 const std::string &filename3)
{
  ASDMapFileName = filename1;
  ChannelMapFileName = filename2;
  GeometryFileName = filename3;
}

const int CHMASK=0x003F;  /* 7 (0-63) */
const int SLMASK=0x001F;  /* 5 (0-31) */
const int CRMASK=0x000F;  /* 4 (0-15) */
const int CHSHIFT= 0;
const int SLSHIFT=10;
const int CRSHIFT=20;

#define KEY(cr,sl,ch) \
((((cr)&CRMASK)<<CRSHIFT) | (((sl)&SLMASK)<<SLSHIFT) | \
 (((ch)&CHMASK)<<CHSHIFT))

const int WIREMASK =0x00FF;  /* 8 (0-255) */
const int LAYERMASK=0x001F;  /* 5 (0- 31) */
const int LRMASK   =0x0007;  /* 3 (0-  7) */
const int WIRESHIFT = 0;
const int LAYERSHIFT=10;
const int LRSHIFT   =20;

#define RKEY(layer,wire) \
((((layer)&LAYERMASK)<<LAYERSHIFT) | (((wire)&WIREMASK)<<WIRESHIFT))

const int MAXCHAR = 144;

// used at once
const int MASK  = 0x00FF;
const int SHIFT = 16;
#define ASDKEY(x,y) \
((((x)&MASK)<<SHIFT) | ((y)&MASK))

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::Initialize()
{
  static const std::string funcname = "CDCWireMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";

  int cid;
  double x,y,z;
  double dx,dy,dz;
  int n,cr,sl,ch,ud,layer,wire,layer0,wire0,nwires;
  int slayer, asdnum, ttype, asdch;
  double radius,phi,dphi,tilt;
  double zlen, r1, r2;
  unsigned int key, rkey;
  char str[MAXCHAR];

  FILE *fp0, *fp1, *fp2;
  if( (fp0=fopen(GeometryFileName.c_str(), "r"))==0 ){
    std::cerr << funcname << " File open fail : " << GeometryFileName << std::endl;
    exit(-1);
  }
  if( (fp1=fopen(ChannelMapFileName.c_str(), "r"))==0 ){
    std::cerr << funcname << " File open fail : " << ChannelMapFileName << std::endl;
    exit(-1);
  }
  if( (fp2=fopen(ASDMapFileName.c_str(), "r"))==0 ){
    std::cerr << funcname << " File open fail : " << ASDMapFileName << std::endl;
    exit(-1);
  }

  while( fgets(str, MAXCHAR, fp0)!=0 ){
    if( str[0]=='#' ) continue;
    if( (n=sscanf(str, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &cid, &layer, &x, &y, &z, &dx, &dy, &dz, &r1, &r2, &phi, &zlen, &dphi ))==13 ){
      if(cid!=CID_CDC||layer!=0) continue;
      GX=x; GY=y; GZ=z; dGX=dx; dGY=dy; dGZ=dz; InnerRadius=r1; OuterRadius=r2; ZLengthOfWire=zlen; MotherVolume=(int)dphi;
      RotationAngle = dz; 
    }
    else if( (n=sscanf(str, "%d %d %d %lf %lf %lf %lf", &cid, &layer, &nwires, &radius, &phi, &dphi, &tilt ))==7 ){
      if( layer<=NumOfCDCLayers ){
	NWires[layer] = nwires;
	Radius[layer] = radius;
	Phi0[layer] = phi;
	dPhi[layer] = dphi;
	Tilt[layer] = tilt;
      }
    }
  }
  fclose(fp0);

#if 0
  for( int i=1; i<=NumOfCDCLayers; i++ ){
    std::cout << i << "  " << NWires[i] << "  " Radius[i] << "  " << Phi0[i] << "  " << dPhi[i] << "  " << Tilt[i] << std::endl;
  }
#endif  

  CDCWireMap  *fMap = 0;
  RCDCWireMap *bMap = 0;
  
  wContainer.clear();
  rwContainer.clear();

  std::map <int,int> chmap;
  typedef std::map <int,int>::const_iterator chmap_iterator;

  while( fgets(str, MAXCHAR, fp1)!=0 ){
    if( str[0]=='#' ) continue;
    if( (n=sscanf(str, "%d %d %d %d", &ttype, &asdch, &layer, &wire ))==4 ){
      chmap[ASDKEY(ttype,asdch)] = ASDKEY(layer,wire);
    }
  }
  fclose(fp1);

#if 0
  for( int t=1; t<=15; t++ ){
    for( int c=0; c<16; c++ ){
      std::cout << t << "  " << c << "  "
		<< chmap[ASDKEY(t,c)] << "  "
		<< (((chmap[ASDKEY(t,c)])>>SHIFT)&MASK) << "  "
		<< (((chmap[ASDKEY(t,c)])       )&MASK) << std::endl;
    }
  }	
#endif

  int slayer_old = 1;
  int layer_old=0;
  int wire_old[3];
  int maxlayer=0;
  int maxwire[3]; 
  for(int i=0;i<3;i++){
    wire_old[i]=0;
    maxwire[i]=0;
  }

  while( fgets(str, MAXCHAR, fp2)!=0 ){
    if( str[0]=='#' ) continue;
    if( (n=sscanf(str, "%d %d %d %d %d %d", &slayer, &asdnum, &ttype, &cr, &sl, &ud))==6 ){
      if( slayer!=slayer_old ){
	layer_old = maxlayer;
	for(int i=0;i<3;i++){
	  wire_old[i] = 0;
	  maxwire[i] = 0;
	}
      }
#if 0
      std::cout << "slayer:" << slayer << " asdnum:" << asdnum << " ttype:" << ttype
		<< " cr:" << cr << " sl:" << sl << " ud:" << ud 
		<< " layer_old:" << layer_old
		<< " wire_old:" << wire_old[0] << " " << wire_old[1] << "  " << wire_old[2]
		<< std::endl;
#endif
      for( asdch=0; asdch<16; asdch++ ){
	int lw;
	lw = chmap[ASDKEY(ttype,asdch)];
	layer0 = (((lw)>>SHIFT)&MASK);
	wire0  = ( (lw)        &MASK);
	if( layer0==0 || wire0==0 ) continue;
	
	wire = wire_old[layer0-1]+wire0;
	layer = layer_old+layer0;
	radius = Radius[layer];
	phi = Phi0[layer]+dPhi[layer]*(wire-1);
	tilt = Tilt[layer];
	ch = asdch+ud*16;
	key = KEY(cr,sl,ch);
	if( (fMap = new CDCWireMap(layer,wire,slayer,asdnum,ttype,asdch,radius,phi,tilt,ZLengthOfWire)) ){
	  wContainer[key] = *fMap;
	  delete fMap;
	}
	else{
	  std::cerr << "[" << funcname << "] : new fail " << std::endl;
	  exit(-1);
	}
	rkey = RKEY(layer,wire);
	if( (bMap = new RCDCWireMap(cr,sl,ch,slayer,asdnum,ttype,asdch,radius,phi,tilt,ZLengthOfWire)) ){
	  rwContainer[rkey] = *bMap;
	  delete bMap;
	}
	else{
	  std::cerr << "[" << funcname << "] : new fail " << std::endl;
	  exit(-1);
	}
	
	if(maxlayer<layer) maxlayer=layer;
	if(maxwire[layer0-1]<wire) maxwire[layer0-1]=wire;
#if 0	
	std::cout << "layer:" << layer << " wire:" << wire
		  << " slayer:" << slayer << " asdnum:" << asdnum
		  << " ttype:" << ttype << " asdch:" << asdch
		  << " rad:" << radius << " phi:" << phi << " tilt:" << tilt
		  << " cr:" << cr << " sl:" << sl << " ch:" << ch
		  << std::endl;
	  
#endif
      }
      for(int i=0;i<3;i++) wire_old[i] = maxwire[i];
      slayer_old = slayer;
    
    }
  }
  fclose(fp2);
  //  std::cout << "[" << funcname << "] Initialization finish." << std::endl;
  std::cout << " finish." << std::endl;

  //PrintWireMap();
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetWire( const int &cr, const int &sl, const int &ch,
			  int &slayer, int &asdnum, int &ttype, int &asdch,
			  double &rad, double &phi, double &tilt,
			  int &layer, int &wire )
{
  unsigned int key = KEY(cr,sl,ch);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }

  if(map){
    slayer = map->slayer();
    asdnum = map->asdnum();
    ttype  = map->transtype();
    asdch  = map->asdch();
    rad    = map->radius();
    phi    = map->phi();
    tilt   = map->tilt();
    layer  = map->layer();
    wire   = map->wire();
    return true;
  }
  else{
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TVector3 CDCWireMapMan::GetWirePos( const int &cr, const int &sl, const int &ch)
{
  unsigned int key = KEY(cr,sl,ch);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }
  if(map){
    return map->pos();
  }
  else{
    return TVector3(DEFVECT);
  }
}

TVector3 CDCWireMapMan::GetWirePosp( const int &cr, const int &sl, const int &ch)
{
  unsigned int key = KEY(cr,sl,ch);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }
  if(map){
    return map->posp();
  }
  else{
    return TVector3(DEFVECT);
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetWire( const int &cr, const int &sl, const int &ch,
			  int &layer, int &wire )
{
  unsigned int key = KEY(cr,sl,ch);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }

  if(map){
    layer  = map->layer();
    wire   = map->wire();
    return true;
  }
  else{
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetWire( const int &cr, const int &sl, const int &ch,
			    double &rad, double &phi, double &tilt )
{
  unsigned int key = KEY(cr,sl,ch);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }

  if(map){
    rad    = map->radius();
    phi    = map->phi();
    tilt   = map->tilt();
    return true;
  }
  else{
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetWire( const int &layer, const int &wire,
			    int &slayer, int &asdnum, int &ttype, int &asdch,
			    double &rad, double &phi, double &tilt,
			    int &cr, int &sl, int &ch )
{
  unsigned int rkey = RKEY(layer,wire);
  std::map <unsigned int, RCDCWireMap>::const_iterator ri = rwContainer.find(rkey);
  const RCDCWireMap *map=0;
  if( ri != rwContainer.end() ){
    map = &(ri->second);
  }

  if(map){
    slayer = map->slayer();
    asdnum = map->asdnum();
    ttype  = map->transtype();
    asdch  = map->asdch();
    rad    = map->radius();
    phi    = map->phi();
    tilt   = map->tilt();
    cr     = map->cr();
    sl     = map->sl();
    ch     = map->ch();
    return true;
  }
  else{
    return false;
  }
}    

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetWire( const int &layer, const int &wire,
			    double &rad, double &phi, double &tilt )
{
  unsigned int rkey = RKEY(layer,wire);
  std::map <unsigned int, RCDCWireMap>::const_iterator ri = rwContainer.find(rkey);
  const RCDCWireMap *map=0;
  if( ri != rwContainer.end() ){
    map = &(ri->second);
  }

  if(map){
    rad    = map->radius();
    phi    = map->phi();
    tilt   = map->tilt();
    return true;
  }
  else{
    return false;
  }
}    

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TVector3 CDCWireMapMan::GetWirePos(  int layer,  int wire )
{
  unsigned int rkey = RKEY(layer,wire);
  std::map <unsigned int, RCDCWireMap>::const_iterator ri = rwContainer.find(rkey);
  const RCDCWireMap *map=0;
  if( ri != rwContainer.end() ){
    map = &(ri->second);
  }

  if(map){
    return map->pos();
  }
  else{
    return TVector3(DEFVECT);
  }
}    
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetWirePosDir(  int layer,  int wire , TVector3 &pos, TVector3 &dir)
{
  unsigned int rkey = RKEY(layer,wire);
  std::map <unsigned int, RCDCWireMap>::const_iterator ri = rwContainer.find(rkey);
  const RCDCWireMap *map=0;
  if( ri != rwContainer.end() ){
    map = &(ri->second);
  }

  if(map){
    pos=map->posp();
    dir=map->pos()-map->posp();
    return true;
  }
  else{
    return false;
  }
}    
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetGWirePosDir(  int layer,  int wire , TVector3 &pos, TVector3 &dir)
{
  unsigned int rkey = RKEY(layer,wire);
  std::map <unsigned int, RCDCWireMap>::const_iterator ri = rwContainer.find(rkey);
  const RCDCWireMap *map=0;
  if( ri != rwContainer.end() ){
    map = &(ri->second);
  }

  if(map){
    TVector3 tmppos=map->posp();
    TVector3 tmpposp=map->pos();
    TVector3 posp;
    LocalToGlobal(tmpposp,pos);
    LocalToGlobal(tmppos ,posp);
    dir=posp-pos;
    return true;
  }
  else{
    return false;
  }
}    


// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TVector3 CDCWireMapMan::GetWirePosp(  int layer,  int wire )
{
  unsigned int rkey = RKEY(layer,wire);
  std::map <unsigned int, RCDCWireMap>::const_iterator ri = rwContainer.find(rkey);
  const RCDCWireMap *map=0;
  if( ri != rwContainer.end() ){
    map = &(ri->second);
  }

  if(map){
    return map->posp();
  }
  else{
    return TVector3(DEFVECT);
  }
}    

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::PrintSimpleWireMap( std::ostream &p_out )
{
  int slayer, asdnum, ttype, asdch, cr, sl, ch;
  double rad, phi, tilt;
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    for( int wire=1; wire<=200; wire++ ){
      if( GetWire( layer, wire, slayer, asdnum, ttype, asdch, rad, phi, tilt, cr, sl, ch ) ){
	p_out.setf(std::ios::showpoint);
	p_out << std::setw(5) << cr
	      << std::setw(5) << sl
	      << std::setw(5) << ch
	      << std::setw(5) << CID_CDC
	      << std::setw(5) << layer
	      << std::setw(5) << wire
	      << std::setw(5) << asdnum
	      << std::setw(5) << ttype
	      << std::setw(5) << asdch
	      << std::endl;
      }
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::PrintWireMap()
{
  int slayer, asdnum, ttype, asdch, cr, sl, ch;
  double rad, phi, tilt;
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    for( int wire=1; wire<=200; wire++ ){
      if( GetWire( layer, wire, slayer, asdnum, ttype, asdch, rad, phi, tilt, cr, sl, ch ) ){
	std::cout << "layer:" << std::setw(2) << layer << " wire:" << std::setw(3) << wire
		  << " slayer:" << std::setw(1) << slayer << " asdnum:" << std::setw(2) << asdnum
		  << " ttype:" << std::setw(2) << ttype << " asdch:" << std::setw(2) << asdch
		  << " rad:" << std::setw(6) << rad << " phi:" << std::setw(7) << phi 
		  << " tilt:" << std::setw(7) << tilt
		  << " cr:" << std::setw(2) << cr << " sl:" << std::setw(2) << sl 
		  << " ch:" << std::setw(2) << ch
		  << std::endl;
      }
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetGeom( const int &layer, double &radius, double &phi0, double &dphi, double &tilt )
{
  if( 0<layer && layer<=NumOfCDCLayers ){
    radius = Radius[layer]; phi0 = Phi0[layer]; dphi = dPhi[layer]; tilt = Tilt[layer];
    return true;
  }
  else{
    radius = -999.; phi0 = -999.; dphi = -999.; tilt = -999.;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetGPOS( double &gx, double &gy, double &gz,
			     double &dgx, double &dgy, double &dgz )
{
  gx = GX; gy = GY; gz = GZ; dgx = dGX; dgy = dGY; dgz = dGZ;
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetGParam( double *par )
{
  par[0] = GX; par[1] = GY; par[2] = GZ; par[3] = dGX; par[4] = dGY; par[5] = dGZ;
  par[6] = InnerRadius; par[7]=OuterRadius; par[8]=360; par[9]=ZLengthOfWire; par[10]=MotherVolume; 
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetFrame( double &zlen, double &rin, double &rout )
{
  zlen = ZLengthOfWire; rin = InnerRadius; rout = OuterRadius;
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::HelixLocalToGlobal( double *local , double *global, int charge)
{
  double tmp[5];
  for(int i=0;i<5;i++) tmp[i]=local[i];
  if(tmp[2]>0)     tmp[1]+=RotationAngle*Deg2Rad;
  else if(tmp[2]<0)     tmp[1]-=RotationAngle*Deg2Rad;
  for(int i=0;i<5;i++) global[i]=tmp[i];
  //  std::cout<<"GX,GY,GZ=  "<<GX<<"  "<<GY<<"  "<<GZ<<std::endl;
  MathTools::ChangePivot(TVector3(0,0,0),TVector3(-GX,-GY,-GZ),tmp,global,charge);
  return true;
}  
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::LocalToGlobal( const TVector3 &local , TVector3 &global)
{
  global=local;
  global.RotateZ(RotationAngle*Deg2Rad);
  global+=TVector3(GX,GY,GZ);
  return true;
}  
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GlobalToLocal( const TVector3 &global , TVector3 &local)
{
  local=global-TVector3(GX,GY,GZ);;
  local.RotateZ(-RotationAngle*Deg2Rad);
  return true;
}  
