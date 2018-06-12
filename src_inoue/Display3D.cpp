// Display3D.cpp

#include <cmath>

#include "Display3D.h"

ClassImp(Display3D);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Display3D::Display3D() : TObject()
{
  gSystem->Load("libPhysics.so");
  pad = 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Display3D::Display3D( TVirtualPad *p ) : TObject()
{
  gSystem->Load("libPhysics.so");
  pad = p;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Display3D::~Display3D()
{
  ShapesCon.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDCFrame( ConfMan *conf, Int_t transparency )
{
  int color = kGray;
  SetTube( "__CDCEndUp",   color, 0, 0, -40,  2./2., 15., 53., 0, 360 );
  SetTube( "__CDCEndDown", color, 0, 0,  40,  2./2., 15., 53., 0, 360 );
  SetTube( "__CDCCFRP",    color, 0, 0,   0, 80./2., 15.-0.01, 15.+0.01, 0, 360 );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDH( ConfMan *conf, Int_t transparency )
{
  int color = kGreen;
  double x,y,z,dx,dy,dz;
  double len, wid, th;
  double rad;
  for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
    if( conf->GetGeomMapManager() ){
      conf->GetGeomMapManager()->GetParam( CID_CDH, seg, x, y, z, dx, dy, dz, len, wid, th );
      rad = sqrt(x*x+y*y);
#if 0
      std::cout << " seg:" << seg << " x:" << x << " y:" << y << " z:" << z << " rad:" << rad
		<< " len:" << len << " wid:" << wid << " th:" << th << std::endl;
#endif
      SetTube( Form("__CDHBody%d",seg), color, 0, 0, 0, len/2., rad-th/2., rad+th/2., (seg-1)*10+0.5, seg*10-0.5, transparency );
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::ReSetCDH( ConfMan *conf, Int_t transparency )
{
  int color = kGreen;
  std::map <std::string, Shape*>::iterator is;
  for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
    is = ShapesCon.find( Form("__CDHBody%d",seg) );
    if( is != ShapesCon.end() ){
      (is->second)->SetColor(color);
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDHColor( ConfMan *conf, Int_t seg, Int_t color, Int_t transparency )
{
  std::map <std::string, Shape*>::iterator is;
  is = ShapesCon.find( Form("__CDHBody%d",seg) );
  if( is != ShapesCon.end() ){
    (is->second)->SetColor(color);
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDSSolenoidBody( ConfMan *conf, Double_t phimin, Double_t phimax, Int_t transparency )
{
  int color = kCyan;
  SetTube( "__CDSSolenoidBody", color, 0, 0, 0, 150./2., 58., 80., phimin, phimax, transparency );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDSSolenoidCap( ConfMan *conf, Int_t pattern,  Int_t transparency )
{ // pattern
  // 0001:UpTop, 0010:UpBottom, 0100:DownTop, 1000:DownBottom
  // default : pattern=0b1111
  std::cout << " Pattern:" << pattern << std::endl;
  int color = kCyan;
  if( (pattern&1)==1 )
    SetTube( "__CDSSolenoidCapUpstreamTop",      color, 0, 0, -80, 10./2., 15., 80.,   0., 180., transparency );
  if( (pattern&2)==2 )
    SetTube( "__CDSSolenoidCapUpstreamBottom",   color, 0, 0, -80, 10./2., 15., 80., 180., 360., transparency );
  if( (pattern&4)==4 )
    SetTube( "__CDSSolenoidCapDownstreamTop",    color, 0, 0,  80, 10./2., 15., 80.,   0., 180., transparency );
  if( (pattern&8)==8 )
    SetTube( "__CDSSolenoidCapDownstreamBottom", color, 0, 0,  80, 10./2., 15., 80., 180., 360., transparency );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetTargetOct2010( ConfMan *conf, Int_t transparency )
{
  int color;
  color = kRed;
  SetBox( "__CarbonTarget", color, 0., 0., -5., 5, 5., 0.25, transparency );
  SetBox( "__CopperTarget", color, 0., 0.,  5., 5, 5., 0.05, transparency );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetTarget( ConfMan *conf, Int_t transparency )
{
  int color = kRed;
  SetTube( "__HeliumTarget", color, 0., 0., 0., 2.5, 0., 2.5, 0., 360., transparency );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDCHit( ConfMan *conf, CDSHitMan *cds )
{
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
#if 0
    std::cout << " Layer:" << layer << " nCDC:" << cds->nCDC(layer) << std::endl;
#endif
    for( int i=0; i<cds->nCDC(layer); i++ ){
      int wire = cds->CDC(layer,i)->wire();
      double x,y,z;
      x = cds->CDC(layer,i)->x();
      y = cds->CDC(layer,i)->y();
      z = cds->CDC(layer,i)->z();
      SetPMark( Form("__CDCHit_L%dW%d",layer,wire), 1, &x, &y, &z, kBlue, 2 );
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDCHitWire( ConfMan *conf, CDSHitMan *cds )
{
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
#if 0
    std::cout << " Layer:" << layer << " nCDC:" << cds->nCDC(layer) << std::endl;
#endif
    for( int i=0; i<cds->nCDC(layer); i++ ){
      int wire = cds->CDC(layer,i)->wire();
      double wx[2],wy[2],wz[2];
      wx[0] = cds->CDC(layer,i)->wx();
      wy[0] = cds->CDC(layer,i)->wy();
      wz[0] = cds->CDC(layer,i)->wz();
      wx[1] = cds->CDC(layer,i)->wxp();
      wy[1] = cds->CDC(layer,i)->wyp();
      wz[1] = cds->CDC(layer,i)->wzp();
      SetPLine( Form("__CDCLINE_L%dW%d",layer,wire), 2, wx, wy, wz, kBlack, 2 );
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDHHit( ConfMan *conf, CDSHitMan *cds )
{
  for( int i=0; i<cds->nCDH(); i++ ){
    int seg = cds->CDH(i)->seg();
    if( cds->CDH(i)->CheckRange() ){
      SetCDHColor(conf,seg, kRed, 20 );
      double x,y,z;
      cds->CDH(i)->pos(x,y,z);
      double hitpos = cds->CDH(i)->hitpos();
      z += hitpos;
      SetPMark( Form("__CDHHit%d",seg), 1, &x, &y, &z, kBlue, 5 );
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDCVertexPoint( ConfMan *conf, CDSTrackingMan *trackMan )
{
  /*
  TVector3 vert = trackMan->GetVertex();
  if( 500<vert.Mag() ) return;
  double x = vert.X(), y = vert.Y(), z = vert.Z();
  SetPMark( "__CDCVertex", 1, &x, &y, &z, kBlue, 2 );
  */
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDCTrackHit( ConfMan *conf, CDSTrackingMan *trackMan )
{
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      for( int ih=0; ih<trackMan->Track(itr)->nTrackHit(layer); ih++ ){
	int wire = trackMan->Track(itr)->TrackHit(layer,ih)->wire();//cds->CDC(layer,i)->wire();
	double x,y,z;
	x = trackMan->Track(itr)->TrackHit(layer,ih)->x();
	y = trackMan->Track(itr)->TrackHit(layer,ih)->y();
	z = trackMan->Track(itr)->TrackHit(layer,ih)->z();
	if( z<-500 ) z=0;
	SetPMark( Form("__CDCTrackHit_L%dW%d",layer,wire), 1, &x, &y, &z, kRed, 2 );
      }
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCDCTrack( ConfMan *conf, CDSTrackingMan *trackMan )
{
  double param[5];
  double zlen = conf->GetCDCWireMapManager()->zlen();
  double rin = conf->GetCDCWireMapManager()->rin();
  double rout = conf->GetCDCWireMapManager()->rout();
  double dz = zlen/500.;
  double xx[1000], yy[1000], zz[1000];
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    trackMan->Track(itr)->GetParameters(param);
    int count=0;
    for( int i=0; i<=500; i++ ){
      double z = -zlen/2. + dz*i;
      double x,y;
      trackMan->Track(itr)->XYatZ(x,y,z);
      double r = sqrt(x*x+y*y);
      if( r<rout ){
	xx[count]=x;
	yy[count]=y;
	zz[count]=z; 
	//SetPMark( Form("__CDCTrackHitInTrack_%d_%d",itr,count), 1, &x, &y, &z, kGreen, 1 );
	count++;
      }
    }
    SetPLine( Form("__CDCTrackLine_%d",itr), count, xx, yy, zz, kRed, 2 );
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::DeleteCDSSolenoidBody()
{
  Delete( "__CDSSolenoidBody" );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::DeleteCDSSolenoidCap()
{
  Delete( "__CDSSolenoidCapUpstreamTop" );
  Delete( "__CDSSolenoidCapUpstreamBottom" );
  Delete( "__CDSSolenoidCapDownstreamTop" );
  Delete( "__CDSSolenoidCapDownstreamBottom" );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::DeleteCDH( Int_t seg )
{
  Delete( Form("__CDHBody%d",seg) );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::DeleteCDH()
{
  for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
    Delete( Form("__CDHBody%d",seg) );
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCounterFrame( ConfMan *conf, Int_t cid, Bool_t fLocal, Int_t transparency )
{
  if( cid==CID_TOFstop || cid==CID_T0 || cid==CID_B1 || cid==CID_B2 ||
      cid==CID_PA || cid==CID_BHD ){
    int numseg;
    int hitdirection; // x:0, y:1, z:2
    if( cid==CID_TOFstop  ){ numseg = NumOfTOFstopSegments; hitdirection=1; }
    else if( cid==CID_T0  ){ numseg = NumOfT0Segments;      hitdirection=1; }
    else if( cid==CID_B1  ){ numseg = NumOfB1Segments;      hitdirection=0; }
    else if( cid==CID_B2  ){ numseg = NumOfB2Segments;      hitdirection=0; }
    else if( cid==CID_PA  ){ numseg = NumOfPASegments;      hitdirection=1; }
    else if( cid==CID_BHD ){ numseg = NumOfBHDSegments;     hitdirection=1; }
    
    //pad->cd();
    //TPolyLine pl;
    for( int seg=1; seg<=numseg; seg++ ){
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      double gxc,gyc,gzc,dgxc,dgyc,dgzc;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
      conf->GetGeomMapManager()->GetGParam(cid,gxc,gyc,gzc,dgxc,dgyc,dgzc);
      double x,y,z;
      if( fLocal==true ){ x = xc; y = yc; z =zc; }
      else{               x = xc + gxc; y = yc + gyc; z = zc + gzc; }
#if 0
      std::cout << " seg:" << seg << " x:" << x << " y:" << y << " z:" << z
		<< " dx:" << dx << " dy:" << dy << " dz:" << dz
		<< " len:" << len << " wid:" << wid << " thi:" << th
		<< " hitdirection:" << hitdirection
		<< std::endl;
#endif
      double wx,wy,wz;
      if( hitdirection==0 ){
	wx = len/2.; wy = wid/2.; wz = th/2.;
      }
      else if( hitdirection==1 ){
	wx = wid/2.; wy = len/2.; wz = th/2.;
      }
      else{ // unnatural ... 
	wx = wid/2.; wy = th/2.; wz = len/2.;
      }
      SetBox( Form( "__Counter_ID%dSeg%d",cid,seg), kGreen, x, y, z,
	      wx, wy,wz, transparency );
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetCounterHit( ConfMan *conf, BeamLineHitMan *bl, Int_t cid, Bool_t fLocal, Int_t transparency )
{
  if( cid==CID_TOFstop || cid==CID_T0 || cid==CID_B1 || cid==CID_B2 ||
      cid==CID_PA || cid==CID_BHD ){
    int numseg;
    int hitdirection; // x:0, y:1, z:2
    if( cid==CID_TOFstop  ){ numseg = NumOfTOFstopSegments; hitdirection=1; }
    else if( cid==CID_T0  ){ numseg = NumOfT0Segments;      hitdirection=1; }
    else if( cid==CID_B1  ){ numseg = NumOfB1Segments;      hitdirection=0; }
    else if( cid==CID_B2  ){ numseg = NumOfB2Segments;      hitdirection=0; }
    else if( cid==CID_PA  ){ numseg = NumOfPASegments;      hitdirection=1; }
    else if( cid==CID_BHD ){ numseg = NumOfBHDSegments;     hitdirection=1; }

    std::map < std::string, Shape* >::iterator is;
    for( int seg=1; seg<=numseg; seg++ ){
      HodoscopeLikeHit *hod = bl->Hodo( cid, seg );
      if( hod==0 ) continue;
      if( !hod->CheckRange() ) continue;

      // change a color of segment
      is = ShapesCon.find( Form( "__Counter_ID%dSeg%d",cid,seg) );
      if( is != ShapesCon.end() ){
	(is->second)->SetColor( kRed );
      }

      // mark a hit point
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      double gxc,gyc,gzc,dgxc,dgyc,dgzc;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
      conf->GetGeomMapManager()->GetGParam(cid,gxc,gyc,gzc,dgxc,dgyc,dgzc);
      double hitpos = hod->hitpos();
      double x,y,z;
      if( fLocal==true ){ x = xc;       y = yc;       z = zc; }
      else{               x = xc + gxc; y = yc + gyc; z = zc + gzc; }

      if( hitdirection==0 ){
	x += hitpos;
      }
      else if( hitdirection==1 ){
	y += hitpos;
      }
      else{ // unnatural ... 
	z += hitpos;
      }

      SetPMark( Form("__CounterHit_ID%dSeg%d",cid,seg), 1, &x, &y, &z, kBlue, 2 );
#if 1
      std::cout << " cid:" << cid << " x:" << x << " y:" << y << " z:" << z
		<< " ctmean:" << hod->ctmean() << " emean:" << hod->emean()
		<< std::endl;
#endif
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetBox( const std::string &name, const Int_t &color,
			const Double_t  &x, const Double_t  &y, const Double_t  &z,
			const Double_t &dx, const Double_t &dy, const Double_t &dz,
			const Int_t &transparency )
{
  Shape *aShape;
  aShape = new Box( color, x, y, z, dx, dy, dz );
  aShape->SetTransparency( transparency );
  aShape->SetName(name);
  ShapesCon[name] = aShape;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetSphere( const std::string &name, const Int_t &color,
			   const Double_t  &x, const Double_t  &y, const Double_t  &z,
			   const Double_t &r,  const Int_t &transparency )
{
  Shape *aShape;
  aShape = new Sphere( color, x, y, z, r );
  aShape->SetTransparency( transparency );
  aShape->SetName(name);
  ShapesCon[name] = aShape;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetTube( const std::string &name, const Int_t &color,
			 const Double_t &x, const Double_t &y, const Double_t &z,
			 const Double_t &len, const Double_t &irad, const Double_t &orad,
			 const Double_t &phi_start, const Double_t &phi_end, const Int_t &transparency )
{
  Shape *aShape;
  aShape = new Tube( color, x, y, z, len, irad, orad, phi_start, phi_end );
  aShape->SetTransparency( transparency );
  aShape->SetName(name);
  ShapesCon[name] = aShape;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetPLine( const std::string &name, const Int_t &npoints, Double_t *x, Double_t *y, Double_t *z,
			  const Int_t &color, const Int_t &width )
{
  Line *line = new Line();
  line->SetPoints( name, npoints, x, y, z, color, width );
  line->SetName(name);
  PLineCon[name] = line;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::SetPMark( const std::string &name, const Int_t &npoints, Double_t *x, Double_t *y, Double_t *z,
			  const Int_t &color, const Int_t &size )
{
  Mark *mark = new Mark();
  mark->SetPoints( name, npoints, x, y, z, color, size );
  mark->SetName(name);
  PMarkCon[name] = mark;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::Draw( Option_t *option )
{
  TGLViewer *viewer;
  std::cout << " pad in Draw:" << pad << " name:" << pad->GetName() << std::endl;
  if( pad==0 )
    viewer = (TGLViewer*)gPad->GetViewer3D(option);
  else
    viewer = (TGLViewer*)pad->GetViewer3D(option);

  //viewer->SetCurrentCamera( TGLViewer::kCameraPerspXOZ );
  //viewer->SetCurrentCamera( TGLViewer::kCameraOrthoXOY );
  double ref[3] = {0,0,0};
  viewer->SetGuideState(2,true,false,ref);
  
  TObject::Draw(option);

  std::map < std::string, Line * >::const_iterator ip = PLineCon.begin();
  Line *line ;
  while( ip != PLineCon.end() ){
    line = ip->second;
    TPolyLine3D *pl = new TPolyLine3D( line->np() );
    for( int j=0; j<line->np(); j++ ){
      TVector3 vec = line->vec(j);
      pl->SetPoint( j, vec.X(), vec.Y(), vec.Z() );
    }
    pl->SetLineColor( line->color() );
    pl->SetLineWidth( line->width() );
    pl->Draw(option);
    ip++;
  }

  std::map < std::string, Mark * >::const_iterator im = PMarkCon.begin();
  Mark *mark;
  while( im != PMarkCon.end() ){
    mark = im->second;
    TPolyMarker3D *pl = new TPolyMarker3D( mark->np() );
    for( int j=0; j<mark->np(); j++ ){
      TVector3 vec = mark->vec(j);
      pl->SetPoint( j, vec.X(), vec.Y(), vec.Z() );
    }
    pl->SetMarkerStyle(20);
    pl->SetMarkerColor( mark->color() );
    pl->SetMarkerSize( mark->size() );
    pl->Draw(option);
    im++;
  }

}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::Paint( Option_t * )
{
  //TVirtualViewer3D *viewer = gPad->GetViewer3D();
  TGLViewer *viewer; 
  if( pad == 0 )
    viewer = (TGLViewer*)gPad->GetViewer3D();
  else
    viewer = (TGLViewer*)pad->GetViewer3D();
  //viewer->Clear();

  Shape *shape;

  std::map < std::string, Shape * >::const_iterator is = ShapesCon.begin();
  while( is != ShapesCon.end() ){
    shape = is->second;
    UInt_t reqSections = TBuffer3D::kCore | TBuffer3D::kBoundingBox | TBuffer3D::kShapeSpecific;
    TBuffer3D & buffer = shape->GetBuffer3D(reqSections);
    reqSections = viewer->AddObject(buffer);
    
    if( reqSections != TBuffer3D::kNone ){
      shape->GetBuffer3D(reqSections);
      viewer->AddObject(buffer);
    }
    is++;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::RotationX( const std::string &name,  const double &theta )
{
  std::map <std::string,Shape*>::const_iterator is;
  is = ShapesCon.find(name);
  if( is== ShapesCon.end() ){ return; }
  (is->second)->RotationX( theta );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::RotationY( const std::string &name,  const double &theta )
{
  std::map <std::string,Shape*>::const_iterator is;
  is = ShapesCon.find(name);
  if( is== ShapesCon.end() ){ return; }
  (is->second)->RotationY( theta );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::RotationZ( const std::string &name,  const double &theta )
{
  std::map <std::string,Shape*>::const_iterator is;
  is = ShapesCon.find(name);
  if( is == ShapesCon.end() ){ return; }
  (is->second)->RotationZ( theta );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::Scale( const std::string &name,  const double &sx, const double &sy, const double &sz )
{
  std::map <std::string,Shape*>::const_iterator is;
  is = ShapesCon.find(name);
  if( is== ShapesCon.end() ){ return; }
  (is->second)->Scale( sx, sy, sz );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::Scale( const std::string &name,  const double &sca )
{
  std::map <std::string,Shape*>::const_iterator is;
  is = ShapesCon.find(name);
  if( is== ShapesCon.end() ){ return; }
  (is->second)->Scale( sca );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Display3D::PrintObj()
{
  std::map <std::string, Shape*>::iterator is = ShapesCon.begin();
  while( is!=ShapesCon.end() ){
    std::cout << "Name:" << is->first << std::endl;
    is++;
  }
  std::map <std::string, Line*>::iterator ip = PLineCon.begin();
  while( ip!=PLineCon.end() ){
    std::cout << "Name:" << ip->first << std::endl;
    ip++;
  }
  std::map <std::string, Mark*>::iterator im = PMarkCon.begin();
  while( im!=PMarkCon.end() ){
    std::cout << "Name:" << im->first << std::endl;
    im++;
  }

}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool Display3D::DeleteAll()
{
  ShapesCon.clear();
  PLineCon.clear();
  PMarkCon.clear();
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool Display3D::Delete( const std::string &name )
{
  std::map <std::string, Shape*>::iterator is;
  is = ShapesCon.find(name);
  if( is != ShapesCon.end() ){
    ShapesCon.erase(is);
    return true;
  }

  std::map <std::string, Line*>::iterator ip;
  ip = PLineCon.find(name);
  if( ip != PLineCon.end() ){
    PLineCon.erase(ip);
    return true;
  }

  std::map <std::string, Mark*>::iterator im;
  im = PMarkCon.find(name);
  if( im != PMarkCon.end() ){
    PMarkCon.erase(im);
    return true;
  }

  return false;
}

ClassImp(Line);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Line::Line()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Line::SetPoints( std::string name, Int_t npoints, Double_t *x, Double_t *y, Double_t *z,
		      Int_t color, Int_t width )
{
  Name = name;
  Color = color;
  Width = width;
  for( int i=0; i<npoints; i++ ){
    pCon.push_back( TVector3(x[i],y[i],z[i]) );
  }
}

ClassImp(Mark);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Mark::Mark()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Mark::SetPoints( std::string name, Int_t npoints, Double_t *x, Double_t *y, Double_t *z,
		      Int_t color, Int_t size )
{
  Name = name;
  Color = color;
  Size = size;
  for( int i=0; i<npoints; i++ ){
    mCon.push_back( TVector3(x[i],y[i],z[i]) );
  }
}

ClassImp(Shape);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Shape::Shape( Int_t color, Double_t x, Double_t y, Double_t z ) :
  fX(x), fY(y), fZ(z), fColor(color)
{
  for( int i=0; i<16; i++ ) TransMat[i] = 0;
  TransMat[0] = TransMat[5] = TransMat[10] = TransMat[15] = 1.;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Shape::SetTranslation( Double_t *mat )
{ // mat = TransMat x mat
#if 0
  std::cout << "before " << std::endl;
  for( int i=0; i<16; i++ ){
    std::cout << " mat[" << i << "]=" << mat[i] << std::endl;
  }
  for( int i=0; i<16; i++ ){
    std::cout << " TransMat[" << i << "]=" << TransMat[i] << std::endl;
  }
#endif  
  Double_t tmp[16];
  for( int i=0; i<16; i++ ) tmp[i] = mat[i];
  GLTranslationCalc( TransMat, tmp, mat );
#if 0
  std::cout << "after " << std::endl;
  for( int i=0; i<16; i++ ){
    std::cout << " mat[" << i << "]=" << mat[i] << std::endl;
  }
#endif  
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Shape::GLTranslationCalc( Double_t *mat )
{ // TransMat = mat x TransMat
#if 0
  std::cout << "before [in Calc]" << std::endl;
  for( int i=0; i<16; i++ ){
    std::cout << " TransMat[" << i << "]=" << TransMat[i] << std::endl;
  }
  for( int i=0; i<16; i++ ){
    std::cout << " mat[" << i << "]=" << mat[i] << std::endl;
  }
#endif  
  Double_t tmp[16];
  for( int i=0; i<16; i++ ) tmp[i] = TransMat[i];
  GLTranslationCalc( mat, tmp, TransMat );
#if 0
  std::cout << "after [in Calc]" << std::endl;
  for( int i=0; i<16; i++ ){
    std::cout << " TransMat[" << i << "]=" << TransMat[i] << std::endl;
  }
#endif  
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Shape::GLTranslationCalc( Double_t *mat1, Double_t *mat2, Double_t *mat3 )
{
  // mat3 = mat1 x mat2
  for( int i=0; i<4; i++ ){
    mat3[i   ] = mat1[+i]*mat2[0];
    mat3[i+ 4] = mat1[+i]*mat2[4];
    mat3[i+ 8] = mat1[+i]*mat2[8];
    mat3[i+12] = mat1[+i]*mat2[12];
    for( int j=1; j<4; j++ ){
      mat3[i   ] += mat1[j*4+i]*mat2[j];
      mat3[i+ 4] += mat1[j*4+i]*mat2[j+4];
      mat3[i+ 8] += mat1[j*4+i]*mat2[j+8];
      mat3[i+12] += mat1[j*4+i]*mat2[j+12];
    }
  }
      
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Shape::RotationX( const Double_t &ang )
{
  Double_t ang2 = ang*TMath::DegToRad();
  Double_t ct = cos(ang2), st = sin(ang2);
  Double_t mat[16];
  for( Int_t i=0; i<16; i++ ) mat[i] = 0;
  mat[5] =  ct; mat[6] = st;
  mat[9] = -st; mat[10] = ct;
  mat[0] = mat[15] = 1;
  GLTranslationCalc(mat);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Shape::RotationY( const Double_t &ang )
{
  Double_t ang2 = ang*TMath::DegToRad();
  Double_t ct = cos(ang2), st = sin(ang2);
  Double_t mat[16];
  for( Int_t i=0; i<16; i++ ) mat[i] = 0;
  mat[0] = ct; mat[2] = -st;
  mat[8] = st; mat[10] = ct;
  mat[5] = mat[15] = 1;
  GLTranslationCalc(mat);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Shape::RotationZ( const Double_t &ang )
{
  Double_t ang2 = ang*TMath::DegToRad();
  Double_t ct = cos(ang2), st = sin(ang2);
  Double_t mat[16];
  for( Int_t i=0; i<16; i++ ) mat[i] = 0;
  mat[0] =  ct; mat[1] = st;
  mat[4] = -st; mat[5] = ct;
  mat[10] = mat[15] = 1;
  GLTranslationCalc(mat);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Shape::Scale( const Double_t &sx, const Double_t &sy, const Double_t &sz )
{
  Double_t mat[16];
  for( Int_t i=0; i<16; i++ ) mat[i] = 0;
  mat[0] =  sx; mat[5] = sy; mat[10] = sz; mat[15] = 1;
  GLTranslationCalc(mat);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void Shape::Scale( const Double_t &sca )
{
  Double_t mat[16];
  for( Int_t i=0; i<16; i++ ) mat[i] = 0;
  mat[0] =  sca; mat[5] = sca; mat[10] = sca; mat[15] = 1;
  GLTranslationCalc(mat);
}

ClassImp(Sphere);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Sphere::Sphere( Int_t color,  Double_t x, Double_t y, Double_t z, Double_t r ) :
  Shape(color,x,y,z), fRadius(r)
{}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TBuffer3D & Sphere::GetBuffer3D(UInt_t reqSections)
{
   static TBuffer3DSphere buffer;

   // Complete kCore section - this could be moved to Shape base class
   if (reqSections & TBuffer3D::kCore) { 
      buffer.ClearSectionsValid();
      buffer.fID = this; 
      buffer.fColor = fColor;       // Color index - see gROOT->GetColor()
      buffer.fTransparency = Transparency;//0;     // Transparency 0 (opaque) - 100 (fully transparent)

      // Complete local/master transformation matrix - simple x/y/z 
      // translation. Easiest way to set identity then override the 
      // translation components
      buffer.SetLocalMasterIdentity();
      buffer.fLocalMaster[12] = fX;
      buffer.fLocalMaster[13] = fY;
      buffer.fLocalMaster[14] = fZ;
      SetTranslation( buffer.fLocalMaster );

      buffer.fLocalFrame = kTRUE;  // Local frame

      buffer.fReflection = kFALSE;
      buffer.SetSectionsValid(TBuffer3D::kCore);
   }
   // Complete kBoundingBox section
   if (reqSections & TBuffer3D::kBoundingBox) {
     //Double_t origin[3] = { 0.0, 0.0, 0.0 };
     Double_t origin[3] = { fX, fY, fZ };
      Double_t halfLength[3] = { fRadius, fRadius, fRadius };
      buffer.SetAABoundingBox(origin, halfLength);
      buffer.SetSectionsValid(TBuffer3D::kBoundingBox);
   }
   // Complete kShapeSpecific section
   if (reqSections & TBuffer3D::kShapeSpecific) {
      buffer.fRadiusOuter = fRadius;
      buffer.fRadiusInner = 0.0;
      buffer.fThetaMin    = 0.0;
      buffer.fThetaMax    = 180.0;
      buffer.fPhiMin    = 0.0;
      buffer.fPhiMax    = 360.0;
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   // We don't implement raw tesselation of sphere - hence this will 
   // not appear in viewers which don't support directly (non-OpenGL)
   // Complete kRawSizes section
   if (reqSections & TBuffer3D::kRawSizes) {
      //buffer.SetSectionsValid(TBuffer3D::kRawSizes);
   }
   // Complete kRaw section
   if (reqSections & TBuffer3D::kRaw) {
      //buffer.SetSectionsValid(TBuffer3D::kRaw);
   }

   //printf( "np:%d\n",buffer.NbPnts() ); 

   return buffer;
}

ClassImp(Box);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Box::Box( Int_t color,  Double_t x, Double_t y, Double_t z,
	  Double_t dx, Double_t dy, Double_t dz ) :
  Shape(color,x,y,z), fDX(dx), fDY(dy), fDZ(dz)
{}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TBuffer3D & Box::GetBuffer3D( UInt_t reqSections )
{
  static TBuffer3D buffer( TBuffer3DTypes::kGeneric );

  if( reqSections & TBuffer3D::kCore ){
    buffer.ClearSectionsValid();
    buffer.fID = this;
    buffer.fColor = fColor;
    buffer.fTransparency = Transparency;//0;

    buffer.SetLocalMasterIdentity();
    buffer.fLocalMaster[12] = fX;
    buffer.fLocalMaster[13] = fY;
    buffer.fLocalMaster[14] = fZ;

    SetTranslation( buffer.fLocalMaster );

    buffer.fLocalFrame = kTRUE;

    buffer.fReflection = kFALSE;
    buffer.SetSectionsValid(TBuffer3D::kCore);
  }
  // Complete kBoundingBox section
  if( reqSections & TBuffer3D::kBoundingBox ){
    Double_t origin[3] = {fX,fY,fZ};
    //Double_t origin[3] = {0,0,0};
    Double_t halfLength[3] = { fDX, fDY, fDZ };
    buffer.SetAABoundingBox(origin,halfLength);
    buffer.SetSectionsValid(TBuffer3D::kBoundingBox);
  }
#if 0
  printf( "X,Y,Z=(%lf,%lf,%lf)\n", fX, fY, fZ );
  printf( "dX,dY,dZ=(%lf,%lf,%lf)\n", fDX, fDY, fDZ );
#endif  
  // Complete kRawSizes section
  if( reqSections & TBuffer3D::kRawSizes ){
    buffer.SetRawSizes( 8, 3*8, 12, 3*12, 6, 6*6 );
    buffer.SetSectionsValid( TBuffer3D::kRawSizes );
  }

  // Complete kRaw section
  if( reqSections & TBuffer3D::kRaw ){
      // Points (8)
      // 3 components: x,y,z
//       buffer.fPnts[ 0] = fX - fDX; buffer.fPnts[ 1] = fY - fDY; buffer.fPnts[ 2] = fZ - fDZ; // 0
//       buffer.fPnts[ 3] = fX + fDX; buffer.fPnts[ 4] = fY - fDY; buffer.fPnts[ 5] = fZ - fDZ; // 1
//       buffer.fPnts[ 6] = fX + fDX; buffer.fPnts[ 7] = fY + fDY; buffer.fPnts[ 8] = fZ - fDZ; // 2
//       buffer.fPnts[ 9] = fX - fDX; buffer.fPnts[10] = fY + fDY; buffer.fPnts[11] = fZ - fDZ; // 3
//       buffer.fPnts[12] = fX - fDX; buffer.fPnts[13] = fY - fDY; buffer.fPnts[14] = fZ + fDZ; // 4
//       buffer.fPnts[15] = fX + fDX; buffer.fPnts[16] = fY - fDY; buffer.fPnts[17] = fZ + fDZ; // 5
//       buffer.fPnts[18] = fX + fDX; buffer.fPnts[19] = fY + fDY; buffer.fPnts[20] = fZ + fDZ; // 6
//       buffer.fPnts[21] = fX - fDX; buffer.fPnts[22] = fY + fDY; buffer.fPnts[23] = fZ + fDZ; // 7
      buffer.fPnts[ 0] = - fDX; buffer.fPnts[ 1] = - fDY; buffer.fPnts[ 2] = - fDZ; // 0
      buffer.fPnts[ 3] = + fDX; buffer.fPnts[ 4] = - fDY; buffer.fPnts[ 5] = - fDZ; // 1
      buffer.fPnts[ 6] = + fDX; buffer.fPnts[ 7] = + fDY; buffer.fPnts[ 8] = - fDZ; // 2
      buffer.fPnts[ 9] = - fDX; buffer.fPnts[10] = + fDY; buffer.fPnts[11] = - fDZ; // 3
      buffer.fPnts[12] = - fDX; buffer.fPnts[13] = - fDY; buffer.fPnts[14] = + fDZ; // 4
      buffer.fPnts[15] = + fDX; buffer.fPnts[16] = - fDY; buffer.fPnts[17] = + fDZ; // 5
      buffer.fPnts[18] = + fDX; buffer.fPnts[19] = + fDY; buffer.fPnts[20] = + fDZ; // 6
      buffer.fPnts[21] = - fDX; buffer.fPnts[22] = + fDY; buffer.fPnts[23] = + fDZ; // 7

      // Segments (12)
      // 3 components: segment color(ignored), start point index, end point index
      // Indexes reference the above points
      buffer.fSegs[ 0] = fColor   ; buffer.fSegs[ 1] = 0   ; buffer.fSegs[ 2] = 1   ; // 0
      buffer.fSegs[ 3] = fColor   ; buffer.fSegs[ 4] = 1   ; buffer.fSegs[ 5] = 2   ; // 1
      buffer.fSegs[ 6] = fColor   ; buffer.fSegs[ 7] = 2   ; buffer.fSegs[ 8] = 3   ; // 2
      buffer.fSegs[ 9] = fColor   ; buffer.fSegs[10] = 3   ; buffer.fSegs[11] = 0   ; // 3
      buffer.fSegs[12] = fColor   ; buffer.fSegs[13] = 4   ; buffer.fSegs[14] = 5   ; // 4
      buffer.fSegs[15] = fColor   ; buffer.fSegs[16] = 5   ; buffer.fSegs[17] = 6   ; // 5
      buffer.fSegs[18] = fColor   ; buffer.fSegs[19] = 6   ; buffer.fSegs[20] = 7   ; // 6
      buffer.fSegs[21] = fColor   ; buffer.fSegs[22] = 7   ; buffer.fSegs[23] = 4   ; // 7
      buffer.fSegs[24] = fColor   ; buffer.fSegs[25] = 0   ; buffer.fSegs[26] = 4   ; // 8
      buffer.fSegs[27] = fColor   ; buffer.fSegs[28] = 1   ; buffer.fSegs[29] = 5   ; // 9
      buffer.fSegs[30] = fColor   ; buffer.fSegs[31] = 2   ; buffer.fSegs[32] = 6   ; // 10
      buffer.fSegs[33] = fColor   ; buffer.fSegs[34] = 3   ; buffer.fSegs[35] = 7   ; // 11
      
      // Polygons (6)
      // 5+ (2+n) components: polygon color (ignored), segment count(n=3+),
      // seg1, seg2 .... segn index
      // Segments indexes refer to the above 12 segments
      // Here n=4 - each polygon defines a rectangle - 4 sides.
      buffer.fPols[ 0] = fColor   ; buffer.fPols[ 1] = 4   ;  buffer.fPols[ 2] = 8  ; // 0
      buffer.fPols[ 3] = 4        ; buffer.fPols[ 4] = 9   ;  buffer.fPols[ 5] = 0  ;
      buffer.fPols[ 6] = fColor   ; buffer.fPols[ 7] = 4   ;  buffer.fPols[ 8] = 9  ; // 1
      buffer.fPols[ 9] = 5        ; buffer.fPols[10] = 10  ;  buffer.fPols[11] = 1  ;
      buffer.fPols[12] = fColor   ; buffer.fPols[13] = 4   ;  buffer.fPols[14] = 10  ; // 2
      buffer.fPols[15] = 6        ; buffer.fPols[16] = 11  ;  buffer.fPols[17] = 2  ;
      buffer.fPols[18] = fColor   ; buffer.fPols[19] = 4   ;  buffer.fPols[20] = 11 ; // 3
      buffer.fPols[21] = 7        ; buffer.fPols[22] = 8   ;  buffer.fPols[23] = 3 ;
      buffer.fPols[24] = fColor   ; buffer.fPols[25] = 4   ;  buffer.fPols[26] = 1  ; // 4
      buffer.fPols[27] = 2        ; buffer.fPols[28] = 3   ;  buffer.fPols[29] = 0  ;
      buffer.fPols[30] = fColor   ; buffer.fPols[31] = 4   ;  buffer.fPols[32] = 7  ; // 5
      buffer.fPols[33] = 6        ; buffer.fPols[34] = 5   ;  buffer.fPols[35] = 4  ;
      
      buffer.SetSectionsValid(TBuffer3D::kRaw);

  }

  return buffer;
}

ClassImp(Tube);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Tube::Tube( Int_t color,  Double_t x, Double_t y, Double_t z,
	    Double_t len, Double_t irad, Double_t orad,
	    Double_t phi_start, Double_t phi_end ) :
  Shape(color,x,y,z), fLength(len), fInnerRadius(irad), fOuterRadius(orad), fPhiStart(phi_start), fPhiEnd(phi_end)
{}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TBuffer3D & Tube::GetBuffer3D(UInt_t reqSections)
{
   static TBuffer3DTubeSeg buffer;

   // Complete kCore section - this could be moved to Shape base class
   if (reqSections & TBuffer3D::kCore) { 
      buffer.ClearSectionsValid();
      buffer.fID = this; 
      buffer.fColor = fColor;       // Color index - see gROOT->GetColor()
      buffer.fTransparency = Transparency;//0;     // Transparency 0 (opaque) - 100 (fully transparent)

      // Complete local/master transformation matrix - simple x/y/z 
      // translation. Easiest way to set identity then override the 
      // translation components
      buffer.SetLocalMasterIdentity();
      buffer.fLocalMaster[12] = fX;
      buffer.fLocalMaster[13] = fY;
      buffer.fLocalMaster[14] = fZ;

      SetTranslation( buffer.fLocalMaster );

      buffer.fLocalFrame = kTRUE;  // Local frame

      buffer.fReflection = kFALSE;
      buffer.SetSectionsValid(TBuffer3D::kCore);
   }
   // Complete kBoundingBox section
   if (reqSections & TBuffer3D::kBoundingBox) {
     //Double_t origin[3] = { 0.0, 0.0, 0.0 };
     Double_t origin[3] = { fX, fY, fZ };
      Double_t halfLength[3] = { fOuterRadius, fOuterRadius, fLength };
      buffer.SetAABoundingBox(origin, halfLength);
      buffer.SetSectionsValid(TBuffer3D::kBoundingBox);
   }
   // Complete kShapeSpecific section
   if (reqSections & TBuffer3D::kShapeSpecific) {
      buffer.fRadiusOuter = fOuterRadius;
      buffer.fRadiusInner = fInnerRadius;
      //buffer.fThetaMin    = 0.0;
      //buffer.fThetaMax    = 180.0;
      buffer.fHalfLength = fLength;
      buffer.fPhiMin    = fPhiStart;
      buffer.fPhiMax    = fPhiEnd;
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   // We don't implement raw tesselation of sphere - hence this will 
   // not appear in viewers which don't support directly (non-OpenGL)
   // Complete kRawSizes section
   if (reqSections & TBuffer3D::kRawSizes) {
     //std::cout << " kRawSize ... in Tube" << std::endl;
      //buffer.SetSectionsValid(TBuffer3D::kRawSizes);
   }
   // Complete kRaw section
   if (reqSections & TBuffer3D::kRaw) {
     //std::cout << " kRaw ... in Tube" << std::endl;
      //buffer.SetSectionsValid(TBuffer3D::kRaw);
   }

   //printf( "np:%d\n",buffer.NbPnts() ); 

   return buffer;
}

