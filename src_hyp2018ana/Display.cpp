#include <cmath>

#include "Display.h"
#include "DetectorList.h"

ClassImp(Display);

Display::Display() : TObject()
{
  frameCDSXY = frameCDSYZ = frameCDSXZ = 0;
  frameBLXZ =frameBLDCXZ =frameBLDCYZ = 0;
  //pl = 0;
}

Display::~Display()
{
  if(frameCDSXY) delete frameCDSXY;
  if(frameCDSYZ) delete frameCDSYZ;
  if(frameCDSXZ) delete frameCDSXZ;
  if(frameBLXZ) delete frameBLXZ;

  if(frameBLDCXZ) delete frameBLDCXZ;
  if(frameBLDCYZ) delete frameBLDCYZ;
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
  if(frameCDSXY) delete frameCDSXY;
  frameCDSXY = new TH2F( "frameCDSXY", "", 100, xmin, xmax, 100, ymin, ymax );
  frameCDSXY->SetStats(0);
}

void Display::SetCDSFrameYZ( double ymin, double ymax, double zmin, double zmax )
{
  if(frameCDSYZ) delete frameCDSYZ;
  frameCDSYZ = new TH2F( "frameCDSYZ", "", 100, ymin, ymax, 100, zmin, zmax );
  frameCDSYZ->SetStats(0);
}

void Display::SetCDSFrameXZ( double xmin, double xmax, double zmin, double zmax )
{
  if(frameCDSXZ) delete frameCDSXZ;
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
  
  if( cid==CID_CDH || cid==CID_IH){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    for( int i=0; i<cds->nHodo(cid); i++ ){
      HodoscopeLikeHit *hit=cds->Hodoi(cid,i);
      int seg = hit->seg();
      TVector3 pos;
      conf->GetGeomMapManager()->GetGPos(cid,seg,pos);
#if 0
      std::cout<<"seg, ctmean, adcu, adcd, x, y  "<<seg<<"\t"<<hit->ctmean()
	       <<"\t"<<hit->adc(0)<<"\t"<<hit->adc(1)
	       <<"\t"<<pos.X()<<"\t"<<pos.Y()<<"  "<<hit->nsensor()<<hit->CheckRange()<<std::endl;
#endif 
      if( !hit->CheckRange() ) continue;
      mark.DrawMarker(pos.X(),pos.Y());
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
      TVector3 pos;
      conf->GetGeomMapManager()->GetPos(cid,seg,pos);
      double hitl=cds->CDH(i)->hitpos();
      // std::cout << " seg:" << seg << " y:" << yc << " z:" << zc << std::endl;
      mark.DrawMarker(pos.Z()+hitl,pos.Y());
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
      TVector3 pos;
      conf->GetGeomMapManager()->GetPos(cid,seg,pos);
      // std::cout << " seg:" << seg << " y:" << yc << " z:" << zc << std::endl;
      mark.DrawMarker(pos.X(),pos.Y());
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


bool Display::DrawCDCTrackHitXY( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan,CDSHitMan *cdsMan )
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerSize(0.5);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    mark.SetMarkerColor(2+itr);
    for( int layer=1; layer<NumOfCDCLayers; layer++ ){
      for( int ih=0; ih<trackMan->Track(itr)->nTrackHit(layer); ih++ ){
	double x,y,z;
	x = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->x();
	y = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->y();
	z = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->z();
	if( z<-500 ) z=0;
	mark.DrawMarker(x,y);
      }
    }
  }
  return true;
}


bool Display::DrawCDCTrackHitYZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan,CDSHitMan *cdsMan )
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerSize(0.5);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    mark.SetMarkerColor(2+itr);
    for( int layer=1; layer<NumOfCDCLayers; layer++ ){
      for( int ih=0; ih<trackMan->Track(itr)->nTrackHit(layer); ih++ ){
	double x,y,z;
	x = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->x();
	y = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->y();
	z = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->z();
	if( z<-500 ) z=0;
	mark.DrawMarker(z,y);
      }
    }
  }
  return true;
}


bool Display::DrawCDCTrackHitXZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan,CDSHitMan *cdsMan )
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerSize(0.5);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    mark.SetMarkerColor(2+itr);
    for( int layer=1; layer<NumOfCDCLayers; layer++ ){
      for( int ih=0; ih<trackMan->Track(itr)->nTrackHit(layer); ih++ ){
	double x,y,z;

	x = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->x();
	y = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->y();
	z = trackMan->Track(itr)->TrackHit(cdsMan,layer,ih)->z();
	if( z<-500 ) z=0;
	mark.DrawMarker(z,x);
      }
    }
  }
  return true;
}

bool Display::DrawCDCTrackXY( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan,CDSHitMan *cdsMan )
{
  pad->cd(); 
  for(int i=0;i<trackMan->nTrack();i++)
    {
      CDSTrack *track=trackMan->Track(i);
      double aparam[5];
      track->GetParameters(aparam);      
      TVector3 ihpos;
      for(int n=0;n<track->nIHHit();n++)
	{

	  ihpos.SetXYZ(track->IHHit(cdsMan,n)->x(),
		       track->IHHit(cdsMan,n)->y(),
		       track->IHHit(cdsMan,n)->z() );
	  //	  std::cout<<"cdh?? "<<cdhpos.x()<<" "<<cdhpos.y()<<" "<<track->CDHHit(cdsMan)->seg()<<" "<<track->CDHHit(cdsMan)->ctmean()<<std::endl;
	  TMarker mihpos;
	  mihpos.SetMarkerStyle(20);
	  //	  mcdhpos.SetMarkerStyle(5);
	  mihpos.SetMarkerColor(2+i);
	  mihpos.SetMarkerSize(0.7);
	  mihpos.DrawMarker(ihpos.x(),ihpos.y());
	}
      TVector3 cdhpos;
      for(int n=0;n<track->nCDHHit();n++)
	{

	  cdhpos.SetXYZ(track->CDHHit(cdsMan,n)->x(),
			track->CDHHit(cdsMan,n)->y(),
			track->CDHHit(cdsMan,n)->z() );
	  std::cout<<"cdh?? "<<cdhpos.x()<<" "<<cdhpos.y()<<" "<<track->CDHHit(cdsMan,n)->seg()<<" "<<track->CDHHit(cdsMan,n)->ctmean()<<std::endl;
	  TMarker mcdhpos;
	  mcdhpos.SetMarkerStyle(20);
	  //	  mcdhpos.SetMarkerStyle(5);
	  mcdhpos.SetMarkerColor(2+i);
	  mcdhpos.SetMarkerSize(0.7);
	  mcdhpos.DrawMarker(cdhpos.x(),cdhpos.y());
	}
	if(aparam[2]==0)
	  {
	    TLine lxy,lzy;
	    lxy.SetLineColor(2+i);
	    double a,b,c;
	    a=track->A(); b=track->B(); c=track->C();

	    double x1,y1,z1,x2,y2,z2;
	    
	    x1=aparam[0]*cos(aparam[1])-(0)*sin(aparam[1]);
	    y1=aparam[0]*sin(aparam[1])+(0)*cos(aparam[1]);
	    z1=aparam[3]-(0)*aparam[4];
	    x2=aparam[0]*cos(aparam[1])-60*sin(aparam[1]);
	    y2=aparam[0]*sin(aparam[1])+60*cos(aparam[1]);
	    z2=aparam[3]-60*aparam[4];

	    lxy.DrawLine(x1,y1,x2,y2);
	  }
	else
	  {    
	    double arho=fabs(1./aparam[2]);
	    double ax_c=(aparam[0]+1./aparam[2])*cos(aparam[1]);
	    double ay_c=(aparam[0]+1./aparam[2])*sin(aparam[1]);
	    
	    TArc arc;
	    arc.SetFillStyle(0);
	    arc.SetLineColor(2+i);
	    arc.DrawArc(ax_c,ay_c,arho,0,360);
	    
	  }
    }
  return true;
}

bool Display::DrawCDCTrackYZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  for(int i=0;i<trackMan->nTrack();i++)
    {
      CDSTrack *track=trackMan->Track(i);
      double aparam[5];
      track->GetParameters(aparam);
      
      if(aparam[2]==0)
	{
	  TLine lzy;
	  double a,b,c;
	  a=track->A(); b=track->B(); c=track->C();	  
	  double x1,y1,z1,x2,y2,z2;
	  
	  x1=aparam[0]*cos(aparam[1])-(0)*sin(aparam[1]);
	  y1=aparam[0]*sin(aparam[1])+(0)*cos(aparam[1]);
	  z1=aparam[3]-(0)*aparam[4];
	  x2=aparam[0]*cos(aparam[1])-60*sin(aparam[1]);
	  y2=aparam[0]*sin(aparam[1])+60*cos(aparam[1]);
	  z2=aparam[3]-60*aparam[4];
	  
	  lzy.SetLineColor(2+i);	    
	  lzy.DrawLine(z1,y1,z2,y2);	        
	}
      else
	{    
	  
	  TF1 *func_y=new TF1("func_y","([0])*sin([1])+1./[2]*(sin([1])-sin( [1]+([3]-x)/([4]*1./[2]) ) )",aparam[3]-1./aparam[2]*aparam[4]*(-TMath::Pi()/2. ),aparam[3]-1./aparam[2]*aparam[4]*(TMath::Pi()/2. ) );
	  func_y->SetParameters(aparam[0],aparam[1],aparam[2],aparam[3],aparam[4]);
	  func_y->SetLineColor(2+i);
	  func_y->SetLineWidth(1);
	  
	  func_y->Draw("LPSAME");
	}
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
  if( cid==CID_CDH || cid == CID_IH ){
    pad->cd();
    TPolyLine pline;
    int nseg=36;
    if(cid==CID_IH) nseg=24;
    for( int seg=1; seg<=nseg; seg++ ){
      conf->GetGeomMapManager()->GetXYCDS(cid,seg,pline);
    }
  }
  return true;
}

bool Display::DrawSegmentsXZ( TVirtualPad *pad, ConfMan *conf, int cid, bool GLOBAL )
{
  TPolyLine pl;
  pad->cd();
  if( cid==CID_CVC || cid==CID_T0 || cid==CID_BHD ||
      cid==CID_PC || cid==CID_NC
      ){
    int numseg=DetectorList::GetInstance()->GetNsegs(cid);
    for( int seg=1; seg<=numseg; seg++ ){
      conf->GetGeomMapManager()->GetZX(cid,seg,pl,GLOBAL);
    }
  }else if(cid==CID_CDH){
    double xc=55.9,wid=3.,gxc=0.;
    double zc=0.,th=79.0,gzc=0.;
    double rot=0;
    MakePLine(zc+gzc,xc+gxc,th,wid,rot,pl);
    pl.Draw();
    MakePLine(zc+gzc,-xc+gxc,th,wid,rot,pl);
    pl.Draw();
  }else if(cid==CID_IH){
    double xc=13.9,wid=0.3,gxc=0.;
    double zc=0.,th=60.0,gzc=0.;
    double rot=0;
    MakePLine(zc+gzc,xc+gxc,th,wid,rot,pl);
    pl.Draw();
    MakePLine(zc+gzc,-xc+gxc,th,wid,rot,pl);
    pl.Draw();
  }
  return true;
}
// ====================================================
bool Display::DrawSegmentsYZ( TVirtualPad *pad, ConfMan *conf, int cid )
{
  TPolyLine pl;
  pad->cd();
  if(cid==CID_CDH){
    double xc=55.9,wid=3.,gxc=0.;
    double zc=0.,th=79.0,gzc=0.;
    double rot=0.;
    MakePLine(zc+gzc,xc+gxc,th,wid,rot,pl);
    MakePLine(zc+gzc,-xc+gxc,th,wid,rot,pl);
  }else if(cid==CID_IH){
    double xc=13.9,wid=0.3,gxc=0.;
    double zc=0.,th=60.0,gzc=0.;
    double rot=0.;
    MakePLine(zc+gzc,xc+gxc,th,wid,rot,pl);
    MakePLine(zc+gzc,-xc+gxc,th,wid,rot,pl);
  }
  return true;
}
// ====================================================
bool Display::MakePLine(const double &xc,const double &yc, const double &dx, const double &dy, const double &rot, TPolyLine &pline)
{
  double x[5],y[5];
  TVector3 tmppos[5];
  tmppos[0]=TVector3( dx/2, dy/2,0);
  tmppos[1]=TVector3(-dx/2, dy/2,0);
  tmppos[2]=TVector3(-dx/2,-dy/2,0);
  tmppos[3]=TVector3( dx/2,-dy/2,0);
  tmppos[4]=TVector3( dx/2, dy/2,0);
  for(int i=0;i<5;i++){
    tmppos[i].RotateZ(rot);
    x[i]=tmppos[i].X()+xc;
    y[i]=tmppos[i].Y()+yc;
  }
  pline.DrawPolyLine(5,x,y);
  return true;
}
// ====================================================
bool Display::DrawCDCLayersYZ( TVirtualPad *pad, ConfMan *conf )
{
  pad->cd();
  TPolyLine pline;  
  pline.SetLineWidth(2);
  double x=90./2.,dx=5;
  double y=0.,dy=106;
  double rot=0;
  MakePLine(x,y,dx,dy,rot,pline);
  MakePLine(-x,y,dx,dy,rot,pline);
  return true;
}
// ====================================================
bool Display::DrawCDCLayersXZ( TVirtualPad *pad, ConfMan *conf )
{
  return DrawCDCLayersYZ(pad,conf);
}
// ====================================================
void Display::SetBLFrameXZ( double xmin, double xmax, double zmin, double zmax )
{
  if(frameBLXZ) delete frameBLXZ;
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
// ====================================================
bool Display::DrawBLHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid )
{
  
  if( cid==CID_CVC || cid==CID_T0 || cid==CID_BHD  ||
      cid==CID_PC || cid==CID_NC
      ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(2);
    int numseg=DetectorList::GetInstance()->GetNsegs(cid);

    for( int seg=1; seg<=numseg; seg++ ){
      HodoscopeLikeHit *hod = bl->Hodo( cid, seg );
      if( hod==0 ) continue;
      if( !hod->CheckRange() ) continue;
      TVector3 pos=hod->pos();
      mark.DrawMarker(pos.Z(),pos.X());
    }
    return true;
  }

  return false;
}

bool Display::DrawBLHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid )
{
  
  if( cid==CID_CVC || cid==CID_T0 || cid==CID_BHD
      ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    int numseg=DetectorList::GetInstance()->GetNsegs(cid);
    for( int seg=1; seg<=numseg; seg++ ){
      HodoscopeLikeHit *hod = bl->Hodo( cid, seg );
      if( hod==0 ) continue;
      if( !hod->CheckRange() ) continue;
      TVector3 pos=hod->pos();
      mark.DrawMarker(pos.Z(),pos.Y());
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
    if(!conf->GetBLDCWireMapManager()->GetParam( cid, layer,
						 nw, z, xy, xy0, dxy, wl, tilt, ra )) continue;
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
  mark.SetMarkerSize(0.5);
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
	  //	  std::cout<<"lay,wire,dl= "<<layer<<"\t"<<hit->wire()<<"\t"<<hit->dl()<<"\t"<<hit->CheckRange()<<std::endl;
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


bool Display::DrawBLDCTrack( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *blTrackMan, int cid,int xy )
{
  //xz plane :xy=0 
  //yz plane :xy=1
  if( !(xy==0 || xy==1) ) return false;
  if( !( cid==CID_BLC1a || cid==CID_BLC1b || cid==CID_BLC1 || cid==CID_BLC2 ||cid==CID_BLC2a || cid==CID_BLC2b || cid==CID_BPC ) ) return false; 

  pad->cd();
  for(int itr=0;itr<blTrackMan->ntrackBLDC(cid);itr++)
    {
      LocalTrack *bldc=blTrackMan->trackBLDC(cid,itr);
      // BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
      // TVector3 gpos,gdir;
      // BLwireman->GetGParam(CID_BLC2a,gpos,gdir);
      double z1,z2,x1,x2,y1,y2;
      z1=-25;z2=25;
      //      bldc->XYLocalPosatZ(z1+gpos.z(),x1,y1);
      //      bldc->XYLocalPosatZ(z2+gpos.z(),x2,y2);
      bldc->XYLocalPosatZ(z1,x1,y1);
      bldc->XYLocalPosatZ(z2,x2,y2);

      TLine line;
      line.SetLineStyle(2);
      line.SetLineColor(itr+2);
      if(xy==0)       line.DrawLine(z1,x1,z2,x2);
      else if(xy==1)  line.DrawLine(z1,y1,z2,y2);
      double a,b,c;
      bldc->abc(a,b,c);

      for(int ih=0;ih<bldc->nhit();ih++)
   	{
   	  ChamberLikeHit *hit=bldc->hit(ih);
   	  TMarker blt_m;	
   	  double hpos,wpos,dltrack;
   	  if(hit->xy()==0) 
   	    {
   	      hpos=hit->x(); wpos=hit->wx();dltrack=fabs(hpos-wpos);
	      if(xy==0)
		{
		  blt_m.SetMarkerStyle(20);
		  blt_m.SetMarkerSize(0.2);
		  blt_m.SetMarkerColor(3+itr );
		  blt_m.DrawMarker(hit->wz(),hit->wx() );
		}
   	    }  
   	  else if(hit->xy()==1) 
   	    {
   	      hpos=hit->y(); wpos=hit->wy();dltrack=fabs(hpos-wpos);
	      if(xy==1)
		{
		  blt_m.SetMarkerStyle(20);
		  blt_m.SetMarkerSize(0.2);
		  blt_m.SetMarkerColor(3+itr );
		  blt_m.DrawMarker(hit->wz(),hit->wy() );
		}
   	    }  
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
  std::cout<<id<<"\ttracking status: "<<trackMan->status(id)<<"\t"<<ntr<<" track in X plane"<<std::endl;
  for( int itr=0; itr<ntr; itr++ ){
    LocalTrack *track=trackMan->trackBLDC(id,itr);
    std::cout<<"track "<<itr<<"\t chi2xz= "<<track->chi2xz()
	     <<"\t chi2= "<<track->chi2all()
	     <<"\t tracktimeX= "<<track->GetTrackTime(0)
	     <<"\t tracktime= "<<track->GetTrackTime()
	     <<"\t RMS= "<<track->GetTrackTimeRMS()<<std::endl;
    for( int ih=0; ih<track->nhit(0); ih++ ){
      double x,y,z;
      ChamberLikeHit *hit=track->hit(0,ih);
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
  int ntr=trackMan->ntrackBLDC(id);
  std::cout<<id<<"\ttracking status: "<<trackMan->status(id)<<"\t"<<ntr<<" track in X plane"<<std::endl;
  
  TPolyLine pline;
  pline.SetLineColor(col);
  for( int itr=0; itr<ntr; itr++ ){
    if(col==-1)
      pline.SetLineColor(itr+2);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[0],x[0],y[0]);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[1],x[1],y[1]);
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
  std::cout<<id<<"\ttracking status: "<<trackMan->status(id)<<"\t"<<ntr<<" track in Y plane"<<std::endl;
  
  for( int itr=0; itr<ntr; itr++ ){
    LocalTrack *track=trackMan->trackBLDC(id,itr);
    std::cout<<"track "<<itr<<"\t chi2yz= "<<track->chi2yz()
	     <<"\t chi2= "<<track->chi2all()
	     <<"\t tracktimeY= "<<track->GetTrackTime(1)
	     <<"\t tracktime= "<<track->GetTrackTime()<<std::endl;
    for( int ih=0; ih<track->nhit(1); ih++ ){
      double x,y,z;
      ChamberLikeHit *hit=track->hit(1,ih);
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
  int ntr=trackMan->ntrackBLDC(id);
  std::cout<<id<<"\ttracking status: "<<trackMan->status(id)<<"\t"<<ntr<<" track in Y plane"<<std::endl;
  TPolyLine pline;
  pline.SetLineColor(col);
  for( int itr=0; itr<ntr; itr++ ){
    if(col==-1)
      pline.SetLineColor(itr+2);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[0],x[0],y[0]);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[1],x[1],y[1]);
    pline.DrawPolyLine(2,z,y);
  }
  return true;
}

bool Display::DrawClusterTime( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id, const int &col)
{
  pad->cd();
  double x,y;
  int ntr=trackMan->ntrackBLDC(id);
  
  TMarker mark;
  mark.SetMarkerSize(0.5);
  mark.SetMarkerStyle(20);

  for( int itr=0; itr<ntr; itr++ ){
    LocalTrack *track=trackMan->trackBLDC(id,itr);
    mark.SetMarkerColor(itr+2);
    for(int xy=0;xy<2;xy++)
      for(int i=0;i<track->ncluster(xy);i++){
	x=track->cluster(xy,i)->hit(0)->z();
	y=track->cluster(xy,i)->GetCTime();
	mark.DrawMarker(x,y);
      }
  }
  return true;
}

bool Display::DrawBLC2TrackfromBLC1YZ( TVirtualPad *pad, ConfMan *conf, LocalTrack *track, const double &mom,const int &col)
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

bool Display::DrawBLC2TrackfromBLC1XZ( TVirtualPad *pad, ConfMan *conf, LocalTrack *track, const double &mom,const int &col)
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

bool Display::DrawFDC(TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl )
{
  DrawFDCWire( pad, conf );
  DrawFDCHitWire( pad, conf ,bl);
  return true;
}

bool Display::DrawFDCWire( TVirtualPad *pad, ConfMan *conf )
{
  pad->cd();
  BLDCWireMapMan *wireman=conf->GetBLDCWireMapManager();

  TLine line;
  line.SetLineStyle(2);
  int nw,xy;
  double z,xy0,dxy,wl,tilt,ra;

  for(int layer=0;layer<3;layer++){
    wireman->GetParam(CID_FDC1,layer*2+1,
		      nw,z,xy,xy0,dxy,wl,tilt,ra);
    TVector2 tmp(0,wl/2);
    TVector2 rot=tmp.Rotate(TMath::Pi()*ra/180.);
    for(int wire=1;wire<=nw;wire++){
      TVector2 zero(xy0+dxy*(wire-1),0);
      TVector2 rot2=zero.Rotate(TMath::Pi()*ra/180.);

#if 0
      line.DrawLine( (zero+rot).X(), (zero+rot).Y(),
		     (zero-rot).X(), (zero-rot).Y() );
#else
      line.DrawLine( (rot2+rot).X(), (rot2+rot).Y(),
		     (rot2-rot).X(), (rot2-rot).Y() );
#endif
    }
  }
  return true;
}

bool Display::DrawFDCHitWire( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl )
{
  pad->cd();
  BLDCWireMapMan *wireman=conf->GetBLDCWireMapManager();

  TLine line;
  line.SetLineWidth(2);
  int nw,xy;
  double z,xy0,dxy,wl,tilt,ra;

  for(int layer=1;layer<=6;layer++){
    line.SetLineColor(layer+1);
    wireman->GetParam(CID_FDC1,layer,
		      nw,z,xy,xy0,dxy,wl,tilt,ra);
    TVector2 tmp(0,wl/2);
    TVector2 rot=tmp.Rotate(TMath::Pi()*ra/180.);    
    for(int i=0;i<bl->nBLDC(CID_FDC1,layer);i++){
      int wire=bl->BLDC(CID_FDC1,layer,i)->wire();
      TVector2 zero(xy0+dxy*(wire-1),0);
      line.DrawLine( (zero+rot).X(), (zero+rot).Y(),
		     (zero-rot).X(), (zero-rot).Y() );
    }
  }
  return true;
}
