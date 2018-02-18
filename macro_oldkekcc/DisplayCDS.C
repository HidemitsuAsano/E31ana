#include <fstream.h>
#include <iostream.h>
#include <vector.h>
#include <math.h>
#include "TMinuit.h"
#include "TMath.h"

#define DISPLAY 1

#define SIM 1

class ConfMan;
class Display;
class CDSHitMan;
class CDCHit;
class CDSTrackingMan;
class CDSTrack;
class CircleFit;
class HelixFit;

const double piMass     = 0.13957018; //GeV
const double pMass      = 0.93827203; //GeV
const double lambdaMass = 1.115683; //GeV

void DisplayCDS()
{

  /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
#if SIM
  //TFile *f = new TFile("simout.root");
  //TFile *f = new TFile("K0production_simout.root");
  TFile *f = new TFile("Lambdaproduction_simout.root");
  TTree *evtree = (TTree*)f->Get("CDSSimTree");
#else
  //TFile *f = new TFile("3030.root");
  TFile *f = new TFile("3045.root");
  TTree *evtree = (TTree*)f->Get("EventTree");
#endif

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Oct2010/analyzer.conf");
  conf->Initialize();

  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;
  CDSHitMan *cdsMan = 0;
  EventHeader *head = 0;
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
#if !SIM
  //evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "EventHeader", &head );
#endif  

  CDSTrackingMan *trackMan=new CDSTrackingMan();

  /*** output file/tree ***/
  TFile *of=new TFile("wtrack.root","recreate");
  TTree *otree = new TTree( "CDSTrackTree", "CDSTrackTree" );
  otree->Branch( "CDSHitMan",  &cdsMan );
  otree->Branch( "CDSTrack", &trackMan );
#if !SIM
  //otree->Branch( "TKOHitCol", &tko );
  otree->Branch( "EventHeader", &head );
#endif

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

#if DISPLAY
  Display *disp = new Display();
  TCanvas *c_disp = new TCanvas( "c_disp", "c_disp", 0, 0, 1200, 600 );
  c_disp->Divide(2,1);
  disp->SetCDSFrameXY(-70,70,-70,70);
  disp->SetCDSFrameYZ(-70,70,-70,70);
#endif

  /*                */
  /* event analysis */
  /*                */
  int nev = evtree->GetEntries();
  std::cout << " AllEvent : " << nev << std::endl;
  for( int iev=0; iev<nev; iev++ ){
    cerr<<"$$$ event # = "<<iev<<endl;
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    //###### event cut######
    int numCDH=0;
    for( int i=0; i<cdsMan->nCDH(); i++ )
      {
        if( 0<cdsMan->CDH(i)->tdcu() && cdsMan->CDH(i)->tdcu()<4000
            && 0<cdsMan->CDH(i)->tdcd() && cdsMan->CDH(i)->tdcd()<4000)
          numCDH++;
      }
    int numCDC=0;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      for( int i=0; i<cdsMan->nCDC(layer); i++ ){
	int tdc = cdsMan->CDC(layer,i)->tdc();
	if( 400<tdc && tdc<900) numCDC++;
      }      
    }
    //if( !(10< numCDC &&numCDC<40 && numCDH>=2 ) ) {continue;}
    if( !( numCDH>=2 ) ) {continue;}
#if !SIM
    int cdh_flag = false;
    for( int i=0; i<cdsMan->nCDH(); i++ ){
      for( int j=i+1; j<cdsMan->nCDH(); j++ ){
	if( 0<cdsMan->CDH(i)->tdcu() && cdsMan->CDH(i)->tdcu()<4000
            && 0<cdsMan->CDH(i)->tdcd() && cdsMan->CDH(i)->tdcd()<4000
	    && 0<cdsMan->CDH(j)->tdcu() && cdsMan->CDH(j)->tdcu()<4000
            && 0<cdsMan->CDH(j)->tdcd() && cdsMan->CDH(j)->tdcd()<4000 ){
	  TVector3 cdh1 = TVector3(cdsMan->CDH(i)->x(), cdsMan->CDH(i)->y(), 0);
	  TVector3 cdh2 = TVector3(cdsMan->CDH(j)->x(), cdsMan->CDH(j)->y(), 0);
	  double ang = cdh1.Angle(cdh2);
	  if( ang>TMath::Pi()*3/4 ) cdh_flag = true;
	}
      }
    }
    if( !cdh_flag ) continue;
#endif
    //###### event cut######

#if 0
    numCDC=0;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      for( int i=0; i<cdsMan->nCDC(layer); i++ ){
	double x = cdsMan->CDC(layer,i)->x();
	double y = cdsMan->CDC(layer,i)->y();
	double z = cdsMan->CDC(layer,i)->z();
	double wx = cdsMan->CDC(layer,i)->wx();
	double wy = cdsMan->CDC(layer,i)->wy();
	double wz = cdsMan->CDC(layer,i)->wz();
	int tdc = cdsMan->CDC(layer,i)->tdc();
	int wire = cdsMan->CDC(layer,i)->wire();
	double dl = cdsMan->CDC(layer,i)->dl();
	double t = cdsMan->CDC(layer,i)->dt();
	double phi = cdsMan->CDC(layer,i)->phi();
	cerr<<"layer, wid, dx, tdc, t : "<<layer<<", "<<wire<<", "<<dl<<", "<<tdc<<", "<<t<<endl;
      }      
    }

    for( int i=0; i<cdsMan->nCDH(); i++ ){
      int seg = cdsMan->CDH(i)->seg();
      int tdcu = cdsMan->CDH(i)->tdcu();
      int tdcd = cdsMan->CDH(i)->tdcd();
      cerr<<"CDH->seg, tdcu, tdcd : "<<seg<<" "<<tdcu<<" "<<tdcd<<endl;
    }
#endif

#if DISPLAY
    TVirtualPad *pad;
    pad = c_disp->GetPad(1);
    disp->DrawCDSFrameXY( pad );
    disp->DrawSegmentsXY( pad, conf, CID_CDH );
    disp->DrawCDCLayersXY( pad, conf );
    disp->DrawCDSHitXY( pad, conf, cdsMan, CID_CDH );
    disp->DrawCDSHitXY( pad, conf, cdsMan, CID_CDC );
    
    pad = c_disp->GetPad(2);
    disp->DrawCDSFrameYZ( pad );
    disp->DrawCDCLayersYZ( pad, conf );
    disp->DrawCDSHitYZ( pad, conf, cdsMan, CID_CDH );
    disp->DrawCDSHitYZ( pad, conf, cdsMan, CID_CDC );

    c_disp->Update();

#endif    

    //-----------------------//
    // tracking & fitting
    //-----------------------//
    int fit_flag[10];
    int nfit = 0;
    for(i=0; i<10; i++) fit_flag[i] = false;

    double arho[10],ax_c[10],ay_c[10],aPt[10];
    double aparam[10][5];

    trackMan->Clear();
    trackMan->SearchCluster(cdsMan);
    trackMan->FindingTrack(cdsMan);

    cerr<<"nTrack = "<<trackMan->nTrack()<<endl;

    if(trackMan->nTrack()>10) {continue;}

    double param[5];
    for(int i=0; i<trackMan->nTrack(); i++){
      cerr<<"=== track # "<<i<<" ==="<<endl;

      CDSTrack *track=trackMan->Track(i);
      
      double a=track->A();
      double b=track->B();
      double c=track->C();
      
      double x_mid=track->MidX();
      double y_mid=track->MidY();
      
      double ap=b/a;
      double bp=-1;
      double cp=y_mid-ap*x_mid;
      
      double jitter=track->Jitter();
      
      if(fabs(jitter)<1e-100 ) continue;
      double rho=fabs( (30*30+jitter*jitter)/(2*jitter) );
      
      double xp=rho/sqrt(1+b*b);
      
      double x_c;
      if(jitter>0 && y_mid>0 )       x_c=x_mid-xp;
      else if(jitter<0 && y_mid>0 )  x_c=x_mid+xp;
      else if(jitter>0 && y_mid<0 )  x_c=x_mid-xp;
      else if(jitter<0 && y_mid<0 )  x_c=x_mid+xp;
      double y_c=b*x_c+cp;
      
      double d_c=fabs( sqrt(x_c*x_c+y_c*y_c) );
      double cos_c=-x_c/d_c;
      double sin_c=-y_c/d_c;

      double x_o=x_c+rho*cos_c;
      double y_o=y_c+rho*sin_c;
      param[0]=sqrt(x_o*x_o+y_o*y_o );

      rho=fabs(rho);
      if(rho<d_c) param[0]=param[0];
      else if(rho>d_c) param[0]=-param[0];

      if(cos_c>0 && sin_c>0)   param[1]=atan(sin_c/cos_c);
      else if(cos_c<0 && sin_c>0)   param[1]=TMath::Pi()+atan(sin_c/cos_c);
      else if(cos_c<0 && sin_c<0)   param[1]=TMath::Pi()+atan(sin_c/cos_c);
      else if(cos_c>0 && sin_c<0)   param[1]=2*TMath::Pi()+atan(sin_c/cos_c);

      if(jitter>0 && y_mid>0 )
	{param[2]=1./rho;param[0]=param[0];param[1]-=TMath::Pi();}
      else if(jitter<0 && y_mid>0 )
	{param[2]=-1./rho;param[0]=-1*param[0];param[1]=param[1];}
      else if(jitter>0 && y_mid<0 )
	{param[2]=-1./rho;param[0]=-1*param[0];param[1]=param[1];}
      else if(jitter<0 && y_mid<0 )
	{param[2]=1./rho;param[0]=param[0];param[1]-=TMath::Pi();}

      param[3]=0;
      param[4]=0;

      //--- circle fit ---//
      track->SetParameters(param);
      if(!track->FirstTrackFitting() ) continue;
      cerr<<"FirstTrackFitting() chi2 =  "<<track->Chi()<<endl;

      bool flag=true;
      //for(int trial=0;trial<10;trial++) { if(!track->TrackFitting() ) flag=false;}
      //if(!flag) continue;
      if(! track->TrackFitting() ) continue;
      cerr<<"TrackFitting() chi2 =  "<<track->Chi()<<endl;
      if( track->Chi()>100 ) continue;

      //--- helix fit ---//      
      if(! track->FirstHelixFitting() ) continue;
      cerr<<"FirstHelixFitting() chi2 =  "<<track->Chi()<<endl;
      if( track->Chi()>100 ) continue;

      //cerr<<track->SecondHelixFitting()<<endl;
      if(! track->SecondHelixFitting() ) continue;
      cerr<<"SecondHelixFitting() chi2 =  "<<track->Chi()<<endl;;
      track->SetHitPos();

      //for(int trial=0;trial<5;trial++) { if(!track->HelixFitting() ) {flag=false;break;}}
      //if(!flag) continue;
      if(! track->HelixFitting() ) continue;
      cerr<<"first iteretaion HelixFitting() chi2 =  "<<track->Chi()<<endl;;

      track->DeleteExtraHit();

#if 0
      int numt=0;
      double Chi=999,Chitmp=999;
      bool hflag1=true;
      while(  numt<20 )	{
	Chi=Chitmp;
	for(int trial=0;trial<10;trial++)
	  { if(!track->HelixFitting() ){hflag1=false; break;} }
	Chitmp=track->Chi();
	if( (Chi-Chitmp)<0.1 ) break;
	numt++;
      }
      if(!hflag1) {std::cout<<"skip event"<<std::endl;continue;}
#endif

      std::cout<<" $$$$$ helix chi2 = "<<track->Chi()<<std::endl;
      fit_flag[i] = true;
      nfit++;

#if 0
      for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
	for(int j=0; j<track->nTrackHit(layer); j++){
	  cerr<<"layer="<<layer<<" wid="<<track->TrackHit(layer,j)->wire()
	      <<" dx="<<track->TrackHit(layer,j)->dl()<<endl;
	}
      }
#endif


      track->GetParameters(aparam[i]);
      std::cout<<"param= "
	       <<aparam[i][0]<<" "
	       <<aparam[i][1]<<" "
	       <<aparam[i][2]<<" "
	       <<aparam[i][3]<<" "
	       <<aparam[i][4]<<" "
	       <<std::endl;
      arho[i]=fabs(1./aparam[i][2]);
      ax_c[i]=(aparam[i][0]+1./aparam[i][2])*cos(aparam[i][1]);
      ay_c[i]=(aparam[i][0]+1./aparam[i][2])*sin(aparam[i][1]);
      aPt[i]=0.3*0.5*arho[i]/100.;
      std::cout<<"Pt"<<i<<"= "<<aPt[i]<<"GeV/c"<<std::endl;

    } // for(int i=0; i<trackMan->nTrack(); i++){      

#if DISPLAY 
    pad = c_disp->GetPad(1);
    pad->cd();
    TArc arc[10];
    
    for(int i=0; i<trackMan->nTrack(); i++){
      if( !fit_flag[i] ) continue;
      arc[i].SetFillStyle(0);
      arc[i].SetLineColor(1+i);
      arc[i].DrawArc(ax_c[i],ay_c[i],arho[i],0,360);
    }

    c_disp->cd(2);     
    TF1 *func_y[10];
    for(int i=0; i<trackMan->nTrack(); i++){
      if( !fit_flag[i] ) continue;
      func_y[i]=new TF1("func_y","([0])*sin([1])+1./[2]*(sin([1])-sin( [1]+([3]-x)/([4]*1./[2]) ) )",aparam[i][3]-1./aparam[i][2]*aparam[i][4]*(-TMath::Pi()/2. ),aparam[i][3]-1./aparam[i][2]*aparam[i][4]*(TMath::Pi()/2. ) );
      func_y[i]->SetParameters(aparam[i][0],aparam[i][1],aparam[i][2],aparam[i][3],aparam[i][4]);
      func_y[i]->SetLineColor(1+i);
      func_y[i]->SetLineWidth(1);
      func_y[i]->Draw("LPSAME");
    }
 
    c_disp->Update();
#endif

    //-----------------------//
    // vertexing
    //-----------------------//
    //if( nfit>9 || nfit<2 ) continue;
    double tpar[10][5];
    if( nfit>2 ){
      for(int i=0; i<trackMan->nTrack(); i++){
	if( !fit_flag[i] ) continue;
	trackMan->Track(i)->GetParameters(tpar[i]);
      }
      TVector3 vtx_tmp[100];
      int nn = 0;
      for(int ii=0; ii<trackMan->nTrack(); ii++){
	for(int jj=ii+1; jj<trackMan->nTrack(); jj++){
	  if( !(fit_flag[ii] && fit_flag[jj]) ) continue;
	  TVector3 vtx1;
	  TVector3 vtx2;
	  TVector3 vtx;
	  CalcHelixDCA(tpar[ii], tpar[jj], vtx1, vtx2, vtx);
	  vtx_tmp[nn] = vtx;
	  nn++;
	  //cerr<<"vtx1 = "<<vtx1.x()<<" "<<vtx1.y()<<" "<<vtx1.z()<<endl;
	  //cerr<<"vtx2 = "<<vtx2.x()<<" "<<vtx2.y()<<" "<<vtx2.z()<<endl;
	  //cerr<<"vtx = "<<vtx.x()<<" "<<vtx.y()<<" "<<vtx.z()<<endl;
	}
      }
      TVector3 vertex;
      for(int ii=0; ii<nn; ii++){
	vertex += vtx_tmp[ii];
      }
      vertex *= 1/double(nn);
      //cerr<<"vertex = "<<vertex.x()<<" "<<vertex.y()<<" "<<vertex.z()<<endl;

      nn = 0;
      TVector3 cvtx[10];
      double dist[10];
      for(int ii=0; ii<trackMan->nTrack(); ii++){
	if( !fit_flag[ii] ) continue;
	double param[5];
	trackMan->Track(ii)->GetParameters(param);
        trackMan->Track(ii)->PointToHelix(vertex, param, cvtx[nn], dist[nn]);
	//cerr<<"cvtx"<<nn<<" = "<<cvtx[nn].x()<<" "<<cvtx[nn].y()<<" "<<cvtx[nn].z()<<endl;
	nn++;
      }
    }// if
    else if( nfit==2 ){    //if( nfit!=2 ) continue;
      int j = 0;
      for(int i=0; i<trackMan->nTrack(); i++){
	if( !fit_flag[i] ) continue;
	trackMan->Track(i)->GetParameters(tpar[j]);
	j++;
      }
      TVector3 vtx1;
      TVector3 vtx2;
      TVector3 vtx;
      CalcHelixDCA(tpar[0], tpar[1], vtx1, vtx2, vtx);
      //cerr<<"vtx1 = "<<vtx1.x()<<" "<<vtx1.y()<<" "<<vtx1.z()<<endl;
      //cerr<<"vtx2 = "<<vtx2.x()<<" "<<vtx2.y()<<" "<<vtx2.z()<<endl;
      //cerr<<"vtx = "<<vtx.x()<<" "<<vtx.y()<<" "<<vtx.z()<<endl;
      
      TVector3 mom1 = CalcHelixMom(tpar[0], vtx1.z());
      TVector3 mom2 = CalcHelixMom(tpar[1], vtx2.z());
      TVector3 mom  = mom1+mom2;
      //cerr<<"mom1 = "<<mom1.x()<<" "<<mom1.y()<<" "<<mom1.z()<<" : Pt ="<<mom1.Pt()<<" GeV/c"<<endl;
      //cerr<<"mom2 = "<<mom2.x()<<" "<<mom2.y()<<" "<<mom2.z()<<" : Pt ="<<mom2.Pt()<<" GeV/c"<<endl;
      //cerr<<"mom = "<<mom.x()<<" "<<mom.y()<<" "<<mom.z()<<" : P ="<<mom.Mag()<<" GeV/c"<<endl;
    }      
//    //-----------------------//
//    // vertexing
//    //-----------------------//
//    if( nfit!=2 ) continue;
//    double tpar[2][5];
//    double tarho[2];
//    double tPt[2];
//    int j = 0;
//    for(int i=0; i<trackMan->nTrack(); i++){
//      if( !fit_flag[i] ) continue;
//      trackMan->Track(i)->GetParameters(tpar[j]);
//      j++;
//    }
//    TVector3 vtx1;
//    TVector3 vtx2;
//    TVector3 vtx;
//    CalcHelixDCA(tpar[0], tpar[1], vtx1, vtx2, vtx);
//    cerr<<"vtx1 = "<<vtx1.x()<<" "<<vtx1.y()<<" "<<vtx1.z()<<endl;
//    cerr<<"vtx2 = "<<vtx2.x()<<" "<<vtx2.y()<<" "<<vtx2.z()<<endl;
//    cerr<<"vtx = "<<vtx.x()<<" "<<vtx.y()<<" "<<vtx.z()<<endl;
//
//    TVector3 mom1 = CalcHelixMom(tpar[0], vtx1.z());
//    TVector3 mom2 =CalcHelixMom(tpar[1], vtx2.z());
//    cerr<<"mom1 = "<<mom1.x()<<" "<<mom1.y()<<" "<<mom1.z()<<" : Pt ="<<mom1.Pt()<<" GeV/c"<<endl;
//    cerr<<"mom2 = "<<mom2.x()<<" "<<mom2.y()<<" "<<mom2.z()<<" : Pt ="<<mom2.Pt()<<" GeV/c"<<endl;
//
//    //-----------------------//
//    // calc. invariant mass
//    //-----------------------//
//    double mass1 = CalcInvMass(mom1, mom2);
//    double mass2 = CalcInvMass(mom2, mom1);
//    double mass = fabs(lambdaMass-mass1) < fabs(lambdaMass-mass2) ? mass1 : mass2;
//    cerr<<"inv.mass = "<<mass1<<" or "<<mass2<<" ---> "<<mass<<endl;

#if DISPLAY
    bool status = disp->Wait();
    if( !status ) return;
#endif    

    otree->Fill();
    trackMan->Clear();

  } //for( int iev=0; iev<nev; iev++ ){

  //-----------------------//
  // file close
  //-----------------------//
  of->Write();
  //of->Print();
  of->Close();

}

TVector3 CalcHelixPos(const double par[5], double z)
{
  double phi2 = par[2]/par[4]*(par[3]-z);
  double x    = (par[0]+1./par[2] )*cos(par[1])-1./par[2]*cos(par[1]+phi2);
  double y    = (par[0]+1./par[2] )*sin(par[1])-1./par[2]*sin(par[1]+phi2);
  return TVector3(x,y,z);
}

double CalcHelixDCA(const double par1[5], const double par2[5],
		    TVector3& vtx1, TVector3& vtx2, TVector3& vtx)
{
  const int NUM = 10;
  const int npoint = 10;
  const double initz = -20; //cm
  double region = 20.0; //+/-cm
  TVector3 pos1[npoint+1];
  TVector3 pos2[npoint+1];
  int num = 1;
  TVector3 now_vtx1;
  TVector3 now_vtx2;
  TVector3 now_vtx;
  double nowz = initz;
  double minl = 999.0;
  while(num<=NUM){
    for(int i=0; i<npoint+1; i++){
      double z=nowz+(2*region/npoint)*(-npoint/2+i);
      pos1[i] = CalcHelixPos(par1,z);
      pos2[i] = CalcHelixPos(par2,z);
    }
    for(int i=0; i<npoint+1; i++){
      for(int j=0; j<npoint+1; j++){
	TVector3 diff = pos1[i]-pos2[j];
        double l = diff.Mag();
        if(l<minl){
          minl = l;
          now_vtx1 = pos1[i];
          now_vtx2 = pos2[j];
          now_vtx = now_vtx1+now_vtx2;
	  now_vtx *=0.5;
          nowz = now_vtx.z();
        }
      }
    }
    num++;
    region = 2*region/10;
  }
  vtx1 = now_vtx1;
  vtx2 = now_vtx2;
  vtx = now_vtx;
  return minl;
}

TVector3 CalcHelixMom(const double par[5], double z)
{
  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = -1*0.5; //T, "-1" is needed.
  double phi2 = (par[3]-z)*par[2]/par[4];
  double pt = fabs(1/par[2])*(Const*dMagneticField)/100.;
  double px = pt*(-1*sin(par[1]+phi2));
  double py = pt*(cos(par[1]+phi2));
  double pz = pt*(par[4]);
  return TVector3(px,py,pz);
}

double CalcInvMass(const TVector3 & piMom, const TVector3& pMom)
{
  double piE = sqrt(piMass*piMass+piMom.Mag2());
  double pE  = sqrt(pMass*pMass+pMom.Mag2());
  double E   = piE+pE;
  TVector3 P = piMom+pMom;
  double M   = E*E-P.Mag2();
  if(M>0) return sqrt(M);
  else    return 0;
}
