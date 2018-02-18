#include <cmath>

#include "Display.h"

ClassImp(Display);

Display::Display() : TObject()
{
  frameCDSXY = frameCDSYZ = frameCDSXZ = 0;
  frameBLXZ =frameBLDCXZ =frameBLDCYZ = 0;
  //pl = 0;
}
bool Display::Wait()
{
  char dispflag;
  char dispin[100]="";

  std::cout << " Input any word or return" << std::endl;
  std::cout << " (q:quit)" << std::endl;
  fgets(dispin,100,stdin);
  if(sscanf(dispin,"%c",&dispflag)==1){
    if( dispflag=='q' ){
      return false;
      //gSystem->Exit(1);
    }
  }
  return true;
}

void Display::SetCDSFrameXY( double xmin, double xmax, double ymin, double ymax )
{
  frameCDSXY = new TH2F( "frameCDSXY", "", 100, xmin, xmax, 100, ymin, ymax );
  frameCDSXY->SetStats(0);
}

void Display::SetCDSFrameYZ( double ymin, double ymax, double zmin, double zmax )
{
  frameCDSYZ = new TH2F( "frameCDSYZ", "", 100, ymin, ymax, 100, zmin, zmax );
  frameCDSYZ->SetStats(0);
}

void Display::SetCDSFrameXZ( double xmin, double xmax, double zmin, double zmax )
{
  frameCDSXZ = new TH2F( "frameCDSXZ", "", 100, xmin, xmax, 100, zmin, zmax );
  frameCDSXZ->SetStats(0);
}

bool Display::DrawCDSFrameXY( TVirtualPad *pad ){
  if( frameCDSXY!=0 ){
    pad->cd();
    frameCDSXY->Draw();
    return true;
  }
  return false;
}

bool Display::DrawCDSFrameYZ( TVirtualPad *pad ){
  if( frameCDSYZ!=0 ){
    pad->cd();
    frameCDSYZ->Draw();
    return true;
  }
  return false;
}

bool Display::DrawCDSFrameXZ( TVirtualPad *pad ){
  if( frameCDSXZ!=0 ){
    pad->cd();
    frameCDSXZ->Draw();
    return true;
  }
  return false;
}

bool Display::DrawCDCLayersXY( TVirtualPad *pad, ConfMan *conf )
{
  pad->cd();
  TArc arc;
  Double_t rad,phi0,dphi,tilt;
  for( Int_t lay=1; lay<=15; lay++ ){
    conf->GetCDCWireMapManager()->GetGeom( lay, rad, phi0, dphi, tilt );
    // std::cout<<"conf "<<lay<<rad<<tilt<<std::endl;
    if( 0<tilt ) arc.SetLineColor(4);
    else if( tilt<0 ) arc.SetLineColor(6);
    else arc.SetLineColor(3);
    arc.SetFillStyle(0);
    arc.DrawArc( 0., 0., rad, 0., 360. );
#if 0
    Double_t cell_size = rad*dphi*3.141592/180.;
    Double_t phi;
    for( Int_t wire=1; wire<=200; wire++ ){
      if( !conf->GetCDCWireMapManager()->GetWire( lay, wire, rad, phi, tilt ) ) continue;
      arc.SetLineColor(4);
      arc.DrawArc( rad*cos(phi*3.141592/180.), rad*sin(phi*3.141592/180.),
                   cell_size/2., 0, 360 );
    }
#endif
  }
  arc.SetFillStyle(0);
  arc.SetLineWidth(2);
  arc.SetLineColor(1);
  arc.DrawArc( 0., 0., 15., 0., 360. );
  arc.DrawArc( 0., 0., 53., 0., 360. );
  arc.SetLineWidth(1);
  return true;
}

bool Display::DrawCDSHitXY( TVirtualPad *pad, ConfMan *conf, CDSHitMan *cds, int cid )
{
  
  if( cid==CID_CDH ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    for( int i=0; i<cds->nCDH(); i++ ){
      int seg = cds->CDH(i)->seg();
      if( !cds->CDH(i)->CheckRange() ) continue;
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      std::cout<<"seg, ctmean, adcu, adcd"<<seg<<"\t"<<cds->CDH(i)->ctmean()
	       <<"\t"<<cds->CDH(i)->adcu()<<"\t"<<cds->CDH(i)->adcd()<<std::endl;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
      mark.DrawMarker(xc,yc);
    }
    return true;
  }
  else if( cid==CID_CDC ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    mark.SetMarkerSize(0.5);
    TArc arc;
    arc.SetFillStyle(0);
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      for( int i=0; i<cds->nCDC(layer); i++ ){
	double wx = cds->CDC(layer,i)->wx();
	double wy = cds->CDC(layer,i)->wy();
	//	double rad = cds->CDC(layer,i)->dl();
	mark.DrawMarker(wx,wy);
	//arc.DrawArc( wx, wy, rad, 0., 360. );
      }
    }
  }

  return false;
}

bool Display::DrawCDSHitYZ( TVirtualPad *pad, ConfMan *conf, CDSHitMan *cds, int cid )
{
  if( cid==CID_CDH ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    for( int i=0; i<cds->nCDH(); i++ ){
      int seg = cds->CDH(i)->seg();
      if( !cds->CDH(i)->CheckRange() ) continue;
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
      double hitl=cds->CDH(i)->hitpos();
      // std::cout << " seg:" << seg << " y:" << yc << " z:" << zc << std::endl;
      mark.DrawMarker(zc+hitl,yc);
    }
    return true;
  }
  else if( cid==CID_CDC ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    mark.SetMarkerSize(0.5);
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      for( int i=0; i<cds->nCDC(layer); i++ ){
	//	double x = cds->CDC(layer,i)->x();
	double y = cds->CDC(layer,i)->y();
	double z = cds->CDC(layer,i)->z();
	mark.DrawMarker(z,y);
	//	std::cout << " layer:" << layer << " y:" << y << " z:" << z << std::endl;
      }
    }
  }

  return false;
}

bool Display::DrawCDSHitXZ( TVirtualPad *pad, ConfMan *conf, CDSHitMan *cds, int cid )
{
  if( cid==CID_CDH ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    for( int i=0; i<cds->nCDH(); i++ ){
      int seg = cds->CDH(i)->seg();
      if( !cds->CDH(i)->CheckRange() ) continue;
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
      // std::cout << " seg:" << seg << " y:" << yc << " z:" << zc << std::endl;
      mark.DrawMarker(zc,xc);
    }
    return true;
  }
  else if( cid==CID_CDC ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    mark.SetMarkerSize(0.5);
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      for( int i=0; i<cds->nCDC(layer); i++ ){
	double x = cds->CDC(layer,i)->x();
	//	double y = cds->CDC(layer,i)->y();
	double z = cds->CDC(layer,i)->z();
	mark.DrawMarker(z,x);
	//	std::cout << " layer:" << layer << " y:" << y << " z:" << z << std::endl;
      }
    }
  }

  return false;
}

bool Display::DrawCDCTrackHitXY( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(2);
  mark.SetMarkerSize(0.5);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    for( int layer=1; layer<NumOfCDCLayers; layer++ ){
      for( int ih=0; ih<trackMan->Track(itr)->nTrackHit(layer); ih++ ){
	double x,y,z;
	x = trackMan->Track(itr)->TrackHit(layer,ih)->x();
	y = trackMan->Track(itr)->TrackHit(layer,ih)->y();
	z = trackMan->Track(itr)->TrackHit(layer,ih)->z();
	if( z<-500 ) z=0;
	mark.DrawMarker(x,y);
      }
    }
  }
  return true;
}

bool Display::DrawCDCTrackHitYZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(2);
  mark.SetMarkerSize(0.5);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    for( int layer=1; layer<NumOfCDCLayers; layer++ ){
      for( int ih=0; ih<trackMan->Track(itr)->nTrackHit(layer); ih++ ){
	double x,y,z;
	x = trackMan->Track(itr)->TrackHit(layer,ih)->x();
	y = trackMan->Track(itr)->TrackHit(layer,ih)->y();
	z = trackMan->Track(itr)->TrackHit(layer,ih)->z();
	if( z<-500 ) z=0;
	mark.DrawMarker(z,y);
      }
    }
  }
  return true;
}

bool Display::DrawCDCTrackHitXZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(2);
  mark.SetMarkerSize(0.5);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    for( int layer=1; layer<NumOfCDCLayers; layer++ ){
      for( int ih=0; ih<trackMan->Track(itr)->nTrackHit(layer); ih++ ){
	double x,y,z;
	x = trackMan->Track(itr)->TrackHit(layer,ih)->x();
	y = trackMan->Track(itr)->TrackHit(layer,ih)->y();
	z = trackMan->Track(itr)->TrackHit(layer,ih)->z();
	if( z<-500 ) z=0;
	mark.DrawMarker(z,x);
      }
    }
  }
  return true;
}

bool Display::DrawCDCTrackXY( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  double zlen = conf->GetCDCWireMapManager()->zlen();
  //  double rin = conf->GetCDCWireMapManager()->rin();
  double rout = conf->GetCDCWireMapManager()->rout();
  double dz = zlen/500.;
  double xx[1000], yy[1000], zz[1000];
  TPolyLine pline;
  pline.SetLineColor(2);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
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
	count++;
      }
    }
    pline.DrawPolyLine(count,xx,yy);
  }
  return true;
}

bool Display::DrawCDCTrackYZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  double zlen = conf->GetCDCWireMapManager()->zlen();
  // double rin = conf->GetCDCWireMapManager()->rin();
  double rout = conf->GetCDCWireMapManager()->rout();
  double dz = zlen/500.;
  double xx[1000], yy[1000], zz[1000];
  TPolyLine pline;
  pline.SetLineColor(2);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
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
	count++;
      }
    }
    pline.DrawPolyLine(count,zz,yy);
  }
  return true;
}

bool Display::DrawCDCTrackXZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  double zlen = conf->GetCDCWireMapManager()->zlen();
  //  double rin = conf->GetCDCWireMapManager()->rin();
  double rout = conf->GetCDCWireMapManager()->rout();
  double dz = zlen/500.;
  double xx[1000], yy[1000], zz[1000];
  TPolyLine pline;
  pline.SetLineColor(2);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
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
	count++;
      }
    }
    pline.DrawPolyLine(count,zz,xx);
  }
  return true;
}

bool Display::DrawSegmentsXY( TVirtualPad *pad, ConfMan *conf, int cid )
{
  if( cid==CID_CDH ){
    pad->cd();
    TPolyLine pl;
    for( int seg=1; seg<=36; seg++ ){
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
#if 0
      std::cout << " seg:" << seg << " x:" << xc << " y:" << yc << " z:" << zc
		<< " dx:" << dx << " dy:" << dy << " dz:" << dz
		<< " len:" << len << " wid:" << wid << " thi:" << th
		<< std::endl;
#endif
      double x[5],y[5];
      double nor = atan2(dy,dx);
      x[0] = xc - th/2.*cos(nor) - wid/2.*sin(nor);
      x[1] = xc + th/2.*cos(nor) - wid/2.*sin(nor);
      x[2] = xc + th/2.*cos(nor) + wid/2.*sin(nor);
      x[3] = xc - th/2.*cos(nor) + wid/2.*sin(nor);
      x[4] = xc - th/2.*cos(nor) - wid/2.*sin(nor);
      y[0] = yc - th/2.*sin(nor) + wid/2.*cos(nor);
      y[1] = yc + th/2.*sin(nor) + wid/2.*cos(nor);
      y[2] = yc + th/2.*sin(nor) - wid/2.*cos(nor);
      y[3] = yc - th/2.*sin(nor) - wid/2.*cos(nor);
      y[4] = yc - th/2.*sin(nor) + wid/2.*cos(nor);
      pl.DrawPolyLine(5,x,y);
    }
  }
  return true;
}

bool Display::DrawSegmentsXZ( TVirtualPad *pad, ConfMan *conf, int cid )
{
  if( cid==CID_TOFstop || cid==CID_T0 || cid==CID_B1 || cid==CID_B2 ||
      cid==CID_PA || cid==CID_BHD ){
    int numseg=0;
    if( cid==CID_TOFstop ) numseg = NumOfTOFstopSegments;
    else if( cid==CID_T0 ) numseg = NumOfT0Segments;
    else if( cid==CID_B1 ) numseg = NumOfB1Segments;
    else if( cid==CID_B2 ) numseg = NumOfB2Segments;
    else if( cid==CID_PA ) numseg = NumOfPASegments;
    else if( cid==CID_BHD ) numseg = NumOfBHDSegments;
    
    pad->cd();
    TPolyLine pl;
    for( int seg=1; seg<=numseg; seg++ ){
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      double gxc,gyc,gzc,dgxc,dgyc,dgzc;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
      conf->GetGeomMapManager()->GetGParam(cid,gxc,gyc,gzc,dgxc,dgyc,dgzc);
#if 0
      std::cout << " seg:" << seg << " x:" << xc << " y:" << yc << " z:" << zc
		<< " dx:" << dx << " dy:" << dy << " dz:" << dz
		<< " len:" << len << " wid:" << wid << " thi:" << th
		<< std::endl;
#endif
      double x[5],z[5];
      x[0] = xc - wid/2. + gxc;
      x[1] = xc - wid/2. + gxc;
      x[2] = xc + wid/2. + gxc;
      x[3] = xc + wid/2. + gxc;
      x[4] = xc - wid/2. + gxc;
      z[0] = zc - th/2. + gzc;
      z[1] = zc + th/2. + gzc;
      z[2] = zc + th/2. + gzc;
      z[3] = zc - th/2. + gzc;
      z[4] = zc - th/2. + gzc;
      pl.DrawPolyLine(5,z,x);
    }
  }
  return true;
}

bool Display::DrawSegmentsYZ( TVirtualPad *pad, ConfMan *conf, int cid )
{
  return true;
}

bool Display::DrawCDCLayersYZ( TVirtualPad *pad, ConfMan *conf )
{
  pad->cd();
  TPolyLine pline;
  
  pline.SetLineWidth(2);
  Double_t xcdc[5], ycdc[5];
  xcdc[0] = -85./2.   ; ycdc[0] =  106./2.;
  xcdc[1] = -85./2.   ; ycdc[1] = -106./2.;
  xcdc[2] = -85./2.-5.; ycdc[2] = -106./2.;
  xcdc[3] = -85./2.-5.; ycdc[3] =  106./2.;
  xcdc[4] = -85./2.   ; ycdc[4] =  106./2.;
  pline.DrawPolyLine(5,xcdc,ycdc);
  xcdc[0] = 85./2.   ; ycdc[0] =  106./2.;
  xcdc[1] = 85./2.   ; ycdc[1] = -106./2.;
  xcdc[2] = 85./2.+5.; ycdc[2] = -106./2.;
  xcdc[3] = 85./2.+5.; ycdc[3] =  106./2.;
  xcdc[4] = 85./2.   ; ycdc[4] =  106./2.;
  pline.DrawPolyLine(5,xcdc,ycdc);
  return true;
}

bool Display::DrawCDCLayersXZ( TVirtualPad *pad, ConfMan *conf )
{
  pad->cd();
  TPolyLine pline;
  
  pline.SetLineWidth(2);
  Double_t xcdc[5], ycdc[5];
  xcdc[0] = -85./2.   ; ycdc[0] =  106./2.;
  xcdc[1] = -85./2.   ; ycdc[1] = -106./2.;
  xcdc[2] = -85./2.-5.; ycdc[2] = -106./2.;
  xcdc[3] = -85./2.-5.; ycdc[3] =  106./2.;
  xcdc[4] = -85./2.   ; ycdc[4] =  106./2.;
  pline.DrawPolyLine(5,xcdc,ycdc);
  xcdc[0] = 85./2.   ; ycdc[0] =  106./2.;
  xcdc[1] = 85./2.   ; ycdc[1] = -106./2.;
  xcdc[2] = 85./2.+5.; ycdc[2] = -106./2.;
  xcdc[3] = 85./2.+5.; ycdc[3] =  106./2.;
  xcdc[4] = 85./2.   ; ycdc[4] =  106./2.;
  pline.DrawPolyLine(5,xcdc,ycdc);
  return true;
}


// ====================================================
void Display::SetBLFrameXZ( double xmin, double xmax, double zmin, double zmax )
{
  frameBLXZ = new TH2F( "frameBLXZ", "", 100, xmin, xmax, 100, zmin, zmax );
  frameBLXZ->SetStats(0);
}

bool Display::DrawBLFrameXZ( TVirtualPad *pad ){
  if( frameBLXZ!=0 ){
    pad->cd();
    frameBLXZ->Draw();
    return true;
  }
  return false;
}

bool Display::DrawBLHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid )
{
  
  if( cid==CID_TOFstop || cid==CID_T0 || cid==CID_B1 || cid==CID_B2 ||
      cid==CID_PA || cid==CID_BHD ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    int numseg=0;
    if( cid==CID_TOFstop ) numseg = NumOfTOFstopSegments;
    else if( cid==CID_T0 ) numseg = NumOfT0Segments;
    else if( cid==CID_B1 ) numseg = NumOfB1Segments;
    else if( cid==CID_B2 ) numseg = NumOfB2Segments;
    else if( cid==CID_PA ) numseg = NumOfPASegments;
    else if( cid==CID_BHD ) numseg = NumOfBHDSegments;

    for( int seg=1; seg<=numseg; seg++ ){
      HodoscopeLikeHit *hod = bl->Hodo( cid, seg );
      if( hod==0 ) continue;
      if( !hod->CheckRange() ) continue;
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      double gxc,gyc,gzc,dgxc,dgyc,dgzc;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
      conf->GetGeomMapManager()->GetGParam(cid,gxc,gyc,gzc,dgxc,dgyc,dgzc);
      std::cout << " x:" << xc+gxc << " z:" << zc+gzc << std::endl;
      mark.DrawMarker(zc+gzc,xc+gxc);
    }
    return true;
  }

  return false;
}

bool Display::DrawBLHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid )
{
  
  if( cid==CID_TOFstop || cid==CID_T0 || cid==CID_B1 || cid==CID_B2 ||
      cid==CID_PA || cid==CID_BHD ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    int numseg=0;
    if( cid==CID_TOFstop ) numseg = NumOfTOFstopSegments;
    else if( cid==CID_T0 ) numseg = NumOfT0Segments;
    else if( cid==CID_B1 ) numseg = NumOfB1Segments;
    else if( cid==CID_B2 ) numseg = NumOfB2Segments;
    else if( cid==CID_PA ) numseg = NumOfPASegments;
    else if( cid==CID_BHD ) numseg = NumOfBHDSegments;

    for( int seg=1; seg<=numseg; seg++ ){
      HodoscopeLikeHit *hod = bl->Hodo( cid, seg );
      if( hod==0 ) continue;
      if( !hod->CheckRange() ) continue;
      double xc,yc,zc,dx,dy,dz,len,wid,th;
      double gxc,gyc,gzc,dgxc,dgyc,dgzc;
      conf->GetGeomMapManager()->GetParam(cid,seg,xc,yc,zc,dx,dy,dz,len,wid,th);
      conf->GetGeomMapManager()->GetGParam(cid,gxc,gyc,gzc,dgxc,dgyc,dgzc);
      std::cout << " y:" << hod->hitpos()+gyc << " z:" << zc+gzc << std::endl;
      mark.DrawMarker(zc+gzc,hod->hitpos()+gyc);
    }
    return true;
  }

  return false;
}



//#######for BLDC  ########

void Display::SetBLDCFrameXZ( double xmin, double xmax, double ymin, double ymax, char* title )
{
  if(frameBLDCXZ) delete frameBLDCXZ;
  frameBLDCXZ = new TH2F( "frameBLDCXZ", title, 100, xmin, xmax, 100, ymin, ymax );
  frameBLDCXZ->SetStats(0);
}


void Display::SetBLDCFrameYZ( double xmin, double xmax, double ymin, double ymax,char* title )
{
  if(frameBLDCYZ) delete frameBLDCYZ;
  frameBLDCYZ = new TH2F( "frameBLDCYZ", title, 100, xmin, xmax, 100, ymin, ymax );
  frameBLDCYZ->SetStats(0);
}


bool Display::DrawBLDCFrameXZ( TVirtualPad *pad ){
  if( frameBLDCXZ!=0 ){
    pad->cd();
    frameBLDCXZ->DrawCopy();
    return true;
  }
  return false;
}


bool Display::DrawBLDCFrameYZ( TVirtualPad *pad ){
  if( frameBLDCYZ!=0 ){
    pad->cd();
    frameBLDCYZ->DrawCopy();
    return true;
  }
  return false;
}


bool Display::DrawBLDCLayersXZ( TVirtualPad *pad, ConfMan *conf,int cid )
{
  pad->cd();
  TMarker mark;
  int nw,xy;
  double z,xy0,dxy,wl,tilt,ra;
  for( int layer=1; layer<=8; layer++ ){
    conf->GetBLDCWireMapManager()->GetParam( cid, layer,
		 nw, z, xy, xy0, dxy, wl, tilt, ra );
    // std::cout<<"conf "<<lay<<rad<<tilt<<std::endl;
    if( xy!=0 ) continue;
    mark.SetMarkerStyle(7);

    for(int wire=0;wire<nw;wire++)
      {
	mark.DrawMarker( z,xy0+wire*dxy );
      }
  }
  return true;
}

bool Display::DrawBLDCLayersYZ( TVirtualPad *pad, ConfMan *conf,int cid )
{
  pad->cd();
  TMarker mark;
  int nw,xy;
  double z,xy0,dxy,wl,tilt,ra;
  for( int layer=1; layer<=8; layer++ ){
    conf->GetBLDCWireMapManager()->GetParam( cid, layer,
		 nw, z, xy, xy0, dxy, wl, tilt, ra );
    // std::cout<<"conf "<<lay<<rad<<tilt<<std::endl;
    if( xy!=1 ) continue;
    mark.SetMarkerStyle(7);
    for(int wire=0;wire<nw;wire++)
      {
	mark.DrawMarker( z,xy0+wire*dxy );
      }
  }
  return true;
}

bool Display::DrawBLDCHit( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *blMan, int cid,int xy )
{
  //xz plane :xy=0 
  //yz plane :xy=1
  if( !(xy==0 || xy==1) ) return false;
  if( !( cid==CID_BLC1a || cid==CID_BLC1b || cid==CID_BLC1 || cid==CID_BLC2 ||cid==CID_BLC2a || cid==CID_BLC2b || cid==CID_BPC ) ) return false; 

  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerSize(1.);
  for( int layer=1; layer<=8; layer++ )
    {
      int nw,xytmp;
      double z,xy0,dxy,wl,tilt,ra;
      conf->GetBLDCWireMapManager()->GetParam( cid, layer,
					       nw, z, xytmp, xy0, 
					       dxy, wl, tilt, ra );
      
      if(xytmp!=xy)	  continue;
      for( int i=0; i<blMan->nBLDC(cid,layer); i++ )
	{
	  double x,z;
	  ChamberLikeHit *hit=blMan->BLDC(cid,layer,i);
	  std::cout<<"lay,wire,dl= "<<layer<<"\t"<<hit->wire()<<"\t"<<hit->dl()<<"\t"<<hit->CheckRange()<<std::endl;
	  if(hit->CheckRange())   mark.SetMarkerColor(2);
	  else   mark.SetMarkerColor(3);
	  if(xy==0)   x=hit->wx();
	  else    x=hit->wy();
	  z=hit->wz();
	  mark.DrawMarker(z,x);
	}
    }
  
  return true;
}

bool Display::DrawBLDCTrackHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id)
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(4);
  mark.SetMarkerSize(0.5);
  
  int ntr=trackMan->ntrackBLDC(id);
  std::cout<<id<<"\t"<<ntr<<" track in X plane"<<std::endl;
  for( int itr=0; itr<ntr; itr++ ){
    std::cout<<"track "<<itr<<"\t chi2xz= "<<trackMan->trackBLDC(id,itr)->chi2xz()<<std::endl;
    for( int ih=0; ih<trackMan->trackBLDC(id,itr)->nhitxz(); ih++ ){
      double x,y,z;
      ChamberLikeHit *hit=trackMan->trackBLDC(id,itr)->hitxz(ih);
      if(hit->leftright())
	x = hit->wx()-hit->dl();
      else
	x = hit->wx()+hit->dl();
      y = hit->y();
      z = hit->z();
      mark.DrawMarker(z,x);
    }
  }  
  return true;
}

bool Display::DrawBLDCTrackXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id, const int &col)
{
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  
  TPolyLine pline;
  pline.SetLineColor(col);
  for( int itr=0; itr<trackMan->ntrackBLDC(id); itr++ ){
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[0],x[0],y[0]);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[1],x[1],y[1]);
    pline.DrawPolyLine(2,z,x);
  }
  return true;
}

bool Display::DrawBLDCLTrackXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id, const int &col)
{
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  
  TPolyLine pline;
  pline.SetLineColor(col);
  for( int itr=0; itr<trackMan->nltrackBLDC(id); itr++ ){
    trackMan->ltrackBLDC(id,itr)->XYLocalPosatZ(z[0],x[0],y[0]);
    trackMan->ltrackBLDC(id,itr)->XYLocalPosatZ(z[1],x[1],y[1]);
    pline.DrawPolyLine(2,z,x);
  }
  return true;
}

bool Display::DrawBLDCTrackHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id)
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(4);
  mark.SetMarkerSize(0.5);

  int ntr=trackMan->ntrackBLDC(id);
  std::cout<<id<<"\t"<<ntr<<" track in Y plane"<<std::endl;
  
  for( int itr=0; itr<trackMan->ntrackBLDC(id); itr++ ){
    std::cout<<"track "<<itr<<"\t chi2yz= "<<trackMan->trackBLDC(id,itr)->chi2yz()<<std::endl;
    for( int ih=0; ih<trackMan->trackBLDC(id,itr)->nhityz(); ih++ ){
      double x,y,z;
      ChamberLikeHit *hit=trackMan->trackBLDC(id,itr)->hityz(ih);
      x = hit->x();
      if(hit->leftright())
	y = hit->wy()-hit->dl();
      else
	y = hit->wy()+hit->dl();
      z = hit->z();
      mark.DrawMarker(z,y);
    }
  }  
  return true;
}

bool Display::DrawBLDCTrackYZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id, const int &col)
{
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  
  TPolyLine pline;
  pline.SetLineColor(col);
  for( int itr=0; itr<trackMan->ntrackBLDC(id); itr++ ){
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[0],x[0],y[0]);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[1],x[1],y[1]);
    pline.DrawPolyLine(2,z,y);
  }
  return true;
}

bool Display::DrawBLDCLTrackYZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id, const int &col)
{
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  
  TPolyLine pline;
  pline.SetLineColor(col);
  for( int itr=0; itr<trackMan->nltrackBLDC(id); itr++ ){
    trackMan->ltrackBLDC(id,itr)->XYLocalPosatZ(z[0],x[0],y[0]);
    trackMan->ltrackBLDC(id,itr)->XYLocalPosatZ(z[1],x[1],y[1]);
    pline.DrawPolyLine(2,z,y);
  }
  return true;
}

bool Display::DrawBLC2TrackfromBLC1YZ( TVirtualPad *pad, ConfMan *conf, LinearTrack *track, const double &mom,const int &col)
{
  /*
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  double a,b,c,d,e,f;
  track->gabc(a,b,c);
  track->gdef(d,e,f);

  double parblc1[6];
  parblc1[0]=c; //cm
  parblc1[1]=TMath::ATan(b)*1000;//mrad
  parblc1[2]=f; //cm
  parblc1[3]=TMath::ATan(e)*1000;//mrad
  parblc1[4]=0.;
  parblc1[5]=(mom-1000)/1000.*100.;
  double mat[36];

  conf->GetTransferMatrixManager()->GetD5Matrix(mat);
  double parblc[6];
  for(int i=0;i<6;i++){
    parblc[i]=0;
    for(int j=0;j<6;j++){
      parblc[i]+=parblc1[j]*mat[6*i+j];
    }
  }
  TVector2 pos(parblc[0],parblc[2]);
  TVector2 dir(parblc[1],parblc[3]);

  TVector2 pos2=pos.Rotate(-135./180.*TMath::Pi());
  TVector2 dir2=dir.Rotate(-135./180.*TMath::Pi());

  y[0]=pos2.Y()+z[0]*TMath::Tan(dir2.Y()/1000.);  
  y[1]=pos2.Y()+z[1]*TMath::Tan(dir2.Y()/1000.);  
  TPolyLine pline;
  pline.SetLineColor(col);
  pline.DrawPolyLine(2,z,y);
  */
  return true;
}

bool Display::DrawBLC2TrackfromBLC1XZ( TVirtualPad *pad, ConfMan *conf, LinearTrack *track, const double &mom,const int &col)
{
  /*
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  double a,b,c,d,e,f;
  track->gabc(a,b,c);
  track->gdef(d,e,f);

  double parblc1[6];
  parblc1[0]=c; //cm
  parblc1[1]=TMath::ATan(b)*1000;//mrad
  parblc1[2]=f; //cm
  parblc1[3]=TMath::ATan(e)*1000;//mrad
  parblc1[4]=0.;
  parblc1[5]=(mom-1000)/1000.*100.;
  double mat[36];

  conf->GetTransferMatrixManager()->GetD5Matrix(mat);
  double parblc[6];
  for(int i=0;i<6;i++){
    parblc[i]=0;
    for(int j=0;j<6;j++){
      parblc[i]+=parblc1[j]*mat[6*i+j];
    }
  }
  TVector2 pos(parblc[0],parblc[2]);
  TVector2 dir(parblc[1],parblc[3]);

  TVector2 pos2=pos.Rotate(-135./180.*TMath::Pi());
  TVector2 dir2=dir.Rotate(-135./180.*TMath::Pi());

  x[0]=pos2.X()+z[0]*TMath::Tan(dir2.X()/1000.);  
  x[1]=pos2.X()+z[1]*TMath::Tan(dir2.X()/1000.);  
  TPolyLine pline;
  pline.SetLineColor(col);
  pline.DrawPolyLine(2,z,x);
  */
  return true;
}

bool Display::DrawBLDCXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
			  BeamLineTrackMan *track,int cid)
{
  DrawBLDCLayersXZ(pad,conf,cid);
  DrawBLDCHit(pad,conf,bl,cid,0);
  DrawBLDCTrackHitXZ(pad,conf,track,cid);
  DrawBLDCTrackXZ(pad,conf,track,cid);
  return true;
}

bool Display::DrawTrackBLDCXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
			  BeamLineTrackMan *track,int cid)
{
  DrawBLDCTrackHitXZ(pad,conf,track,cid);
  DrawBLDCTrackXZ(pad,conf,track,cid,3);
  return true;
}

bool Display::DrawBLDCYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
			  BeamLineTrackMan *track,int cid)
{
  DrawBLDCLayersYZ(pad,conf,cid);
  DrawBLDCHit(pad,conf,bl,cid,1);
  DrawBLDCTrackHitYZ(pad,conf,track,cid);
  DrawBLDCTrackYZ(pad,conf,track,cid);
  return true;
}

bool Display::DrawTrackBLDCYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
			  BeamLineTrackMan *track,int cid)
{
  DrawBLDCTrackHitYZ(pad,conf,track,cid);
  DrawBLDCTrackYZ(pad,conf,track,cid,3);
  return true;
}
