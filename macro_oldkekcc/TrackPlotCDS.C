#include <fstream.h>
#include <iostream.h>
#include <vector.h>
#include <math.h>
#include "TMinuit.h"
#include "TMath.h"

#define DISPLAY 0

#define SIM 1
#define LAMBDA 0

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

void TrackPlotCDS()
{

  /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  //TFile *f = new TFile("wtrack.root");
  //TFile *f = new TFile("wtrack_test5.root");
  //TFile *f = new TFile("wtrack_3045.root");
  //TFile *f = new TFile("tmp.root");
  //TFile *f = new TFile("wtrack_3092.root");
#if LAMBDA
  TFile *f = new TFile("wtrack_Lambdaproduction_simout.root");
#else
  TFile *f = new TFile("wtrack_K0production_simout.root");
#endif
#if SIM
  TTree *evtree = (TTree*)f->Get("CDSTrackTree"); // from macro
#else
  TTree *evtree = (TTree*)f->Get("EventTree"); // from cc
#endif

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Oct2010/analyzer.conf");
  conf->Initialize();

  /*** declaration of classes ***/
  CDSHitMan *cdsMan = 0;
  CDSTrackingMan *trackMan = 0;
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
#if SIM
  evtree->SetBranchAddress( "CDSTrack", &trackMan ); // from macro
#else
  evtree->SetBranchAddress( "CDSTrackingMan", &trackMan ); // from cc
#endif

  //*** hist definition ***//
  TH1F *h_vtxx    = new TH1F("h_vtxx", "vertex x", 100, -20, 20);
  TH1F *h_vtxy    = new TH1F("h_vtxy", "vertex y", 100, -20, 20);
  TH1F *h_vtxz    = new TH1F("h_vtxz", "vertex z", 100, -60, 60);
  TH1F *h_vtxr    = new TH1F("h_vtxr", "vertex r", 100, 0, 20);
  TH1F *h_chi2    = new TH1F("h_chi2", "chi2",     100, 0, 10);
  TH1F *h_mompix  = new TH1F("h_mompix", "mom #pi x", 100, -1, 1);
  TH1F *h_mompiy  = new TH1F("h_mompiy", "mom #pi y", 100, -1, 1);
  TH1F *h_mompiz  = new TH1F("h_mompiz", "mom #pi z", 100, -1, 1);
  TH1F *h_mompi   = new TH1F("h_mompi", "mom #pi", 100, 0, 1);  
  TH1F *h_mompx   = new TH1F("h_mompx", "mom p x", 100, -1, 1);
  TH1F *h_mompy   = new TH1F("h_mompy", "mom p y", 100, -1, 1);
  TH1F *h_mompz   = new TH1F("h_mompz", "mom p z", 100, -1, 1);
  TH1F *h_momp    = new TH1F("h_momp", "mom p", 100, 0, 1);  
  TH1F *h_momlamx = new TH1F("h_momlamx", "mom #Lambda x", 100, -1, 1);
  TH1F *h_momlamy = new TH1F("h_momlamy", "mom #Lambda y", 100, -1, 1);
  TH1F *h_momlamz = new TH1F("h_momlamz", "mom #Lambda z", 100, -1, 1);
  TH1F *h_momlam  = new TH1F("h_momlam", "mom #Lambda", 100, 0, 1);  
#if LAMBDA
  TH1F *h_mass    = new TH1F("h_mass", "invariant mass", 50, 1.05, 1.25);
  TH1F *h_mass1    = new TH1F("h_mass1", "invariant mass1", 50, 1.05, 1.25);
  TH1F *h_mass2    = new TH1F("h_mass2", "invariant mass2", 50, 1.05, 1.25);
#else
  TH1F *h_mass    = new TH1F("h_mass", "invariant mass", 50, 0.35 ,0.65);
  TH1F *h_mass1    = new TH1F("h_mass1", "invariant mass1", 50, 0.35, 0.65);
  TH1F *h_mass2    = new TH1F("h_mass2", "invariant mass2", 50, 0.35, 0.65);
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
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    //if(trackMan->nTrack()) cerr<<"nTrack = "<<trackMan->nTrack()<<endl;
#if !DISPLAY
    if( !trackMan->nTrack() ) continue;
#endif
    //-----------------------//
    // check tracks
    //-----------------------//
    int fit_flag[10];
    int nfit = 0;
    for( int i=0; i<10; i++ ) fit_flag[i] = false;

    for(int i=0; i<trackMan->nTrack(); i++){
      CDSTrack *track=trackMan->Track(i);
      if( track->Chi()<100 ){
	fit_flag[i] = true;
	nfit++;
      }
      else continue;
      h_chi2->Fill(track->Chi());
    }

    //if( nfit ) cerr<<"nfit = "<<nfit<<endl;
    //if( nfit<3 ) continue;

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

    pad = c_disp->GetPad(1);
    pad->cd();
   
    double arho[10],ax_c[10],ay_c[10],aPt[10];
    double aparam[10][5];
    TArc arc[10];

    for(int i=0; i<trackMan->nTrack(); i++){
      if( !fit_flag[i] ) continue;
      CDSTrack *track=trackMan->Track(i);
      //cerr<<"=== track # "<<i<<" ==="<<endl;
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

      arc[i].SetFillStyle(0);
      arc[i].SetLineColor(1+i);
      arc[i].DrawArc(ax_c[i],ay_c[i],arho[i],0,360);
    }
    c_disp->Update();

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
    if( nfit<2 || nfit>9 ) continue;
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
      vertex = vertex/double(nn);
      //cerr<<"vertex = "<<vertex.x()<<" "<<vertex.y()<<" "<<vertex.z()<<endl;
      h_vtxx->Fill(vertex.x());
      h_vtxy->Fill(vertex.y());
      h_vtxz->Fill(vertex.z());
      h_vtxr->Fill(vertex.Mag());


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
    else{    //if( nfit!=2 ) continue;
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
      h_vtxx->Fill(vtx.x());
      h_vtxy->Fill(vtx.y());
      h_vtxz->Fill(vtx.z());
      h_vtxr->Fill(vtx.Mag());
      
      TVector3 mom1 = CalcHelixMom(tpar[0], vtx1.z());
      TVector3 mom2 = CalcHelixMom(tpar[1], vtx2.z());
      TVector3 mom  = mom1+mom2;
      //cerr<<"mom1 = "<<mom1.x()<<" "<<mom1.y()<<" "<<mom1.z()<<" : Pt ="<<mom1.Pt()<<" GeV/c"<<endl;
      //cerr<<"mom2 = "<<mom2.x()<<" "<<mom2.y()<<" "<<mom2.z()<<" : Pt ="<<mom2.Pt()<<" GeV/c"<<endl;
      //cerr<<"mom = "<<mom.x()<<" "<<mom.y()<<" "<<mom.z()<<" : P ="<<mom.Mag()<<" GeV/c"<<endl;
      h_mompix->Fill(mom1.x());
      h_mompiy->Fill(mom1.y());
      h_mompiz->Fill(mom1.z());
      h_mompi->Fill(mom1.Mag());
      h_mompx->Fill(mom2.x());
      h_mompy->Fill(mom2.y());
      h_mompz->Fill(mom2.z());
      h_momp->Fill(mom2.Mag());
      h_momlamx->Fill(mom.x());
      h_momlamy->Fill(mom.y());
      h_momlamz->Fill(mom.z());
      h_momlam->Fill(mom.Mag());
      
      //-----------------------//
      // calc. invariant mass
      //-----------------------//
      //--- charge sign cut ---//
      if( tpar[0][2]*tpar[1][2]>0 ) continue;
      //--- charge sign cut ---//
      double mass1 = CalcInvMass(mom1, mom2);
      double mass2 = CalcInvMass(mom2, mom1);
      double mass = fabs(lambdaMass-mass1) < fabs(lambdaMass-mass2) ? mass1 : mass2;
      //cerr<<"inv.mass = "<<mass1<<" or "<<mass2<<" ---> "<<mass<<endl;
      //h_mass->Fill(mass);
      h_mass->Fill(mass1);
      h_mass->Fill(mass2);
      h_mass1->Fill(mass1);
      h_mass2->Fill(mass2);
    } // else

#if DISPLAY
    bool status = disp->Wait();
    if( !status ) return;
#endif    

  } //for( int iev=0; iev<nev; iev++ ){

#if 1
  //-----------------------//
  // plot hist
  //-----------------------//
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  gROOT->cd();

  TCanvas *c1 = new TCanvas("c1", "");
  h_mass->Draw();
#if LAMBDA
  h_mass->Fit("gaus","ev","",1.1,1.14);
#else
  h_mass->Fit("gaus","ev","",0.45,0.55);
#endif
  h_mass->GetFunction("gaus")->SetLineColor(4);
  h_mass->SetXTitle("p#pi^{-} invariant mass [GeV/c^{2}]");
  c1->Print("tmp1.pdf");

  TCanvas *c2 = new TCanvas("c2", "");
  c2->Divide(2,2);
  c2->cd(1); h_vtxx->Draw();
  c2->cd(2); h_vtxy->Draw();
  c2->cd(3); h_vtxz->Draw();
  c2->cd(4); h_vtxr->Draw();
  c2->Print("tmp2.pdf");

  TCanvas *c3 = new TCanvas("c3", "");
  c3->cd(1); h_chi2->Draw();
  c3->Print("tmp3.pdf");

  TCanvas *c4 = new TCanvas("c4", "");
  c4->Divide(2,2);
  c4->cd(1); h_mompix->Draw();
  c4->cd(2); h_mompiy->Draw();
  c4->cd(3); h_mompiz->Draw();
  c4->cd(4); h_mompi->Draw();
  c4->Print("tmp4.pdf");

  TCanvas *c5 = new TCanvas("c5", "");
  c5->Divide(2,2);
  c5->cd(1); h_mompx->Draw();
  c5->cd(2); h_mompy->Draw();
  c5->cd(3); h_mompz->Draw();
  c5->cd(4); h_momp->Draw();
  c5->Print("tmp5.pdf");

  TCanvas *c6 = new TCanvas("c6", "");
  c6->Divide(2,2);
  c6->cd(1); h_momlamx->Draw();
  c6->cd(2); h_momlamy->Draw();
  c6->cd(3); h_momlamz->Draw();
  c6->cd(4); h_momlam->Draw();
  c6->Print("tmp6.pdf");

  TCanvas *c7 = new TCanvas("c7", "");
  h_vtxz->Draw();
  c7->Print("tmp7.pdf");

  TCanvas *c8 = new TCanvas("c8", "");
  c8->Divide(1,2);
  c8->cd(1); h_mass1->Draw();
  c8->cd(2); h_mass2->Draw();
  c8->Print("tmp8.pdf");
#endif

}

TVector3 CalcHelixPos(const double par[5], double z)
{
  if( par[2] == 0 || par[4] == 0 ) return TVector3(-999,-999,-999);
  double phi2 = par[2]/par[4]*(par[3]-z);
  double x    = (par[0]+1./par[2] )*cos(par[1])-1./par[2]*cos(par[1]+phi2);
  double y    = (par[0]+1./par[2] )*sin(par[1])-1./par[2]*sin(par[1]+phi2);
  return TVector3(x,y,z);
}

double CalcHelixDCA(const double par1[5], const double par2[5],
		    TVector3& vtx1, TVector3& vtx2, TVector3& vtx)
{
  const int NUM = 2;
  const int npoint = 100;
  const double initz = 0; //cm
  double region = 60.0; //+/-cm
  TVector3 pos1[npoint+1];
  TVector3 pos2[npoint+1];
  int num = 1;
  TVector3 now_vtx1;
  TVector3 now_vtx2;
  TVector3 now_vtx;
  double nowz = initz;
  double minl = 999.0;
  while(num<=NUM){
    //cerr<<"---"<<num<<endl;
    for(int i=0; i<npoint+1; i++){
      double z=nowz+(2*region/npoint)*(-npoint/2+i);
      //cerr<<z<<endl;
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
  if( par[2] == 0 || par[4] == 0 ) return TVector3(-999,-999,-999);
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
#if LAMBDA
  double pE  = sqrt(pMass*pMass+pMom.Mag2());
#else
  double pE  = sqrt(piMass*piMass+pMom.Mag2());
#endif
  double E   = piE+pE;
  TVector3 P = piMom+pMom;
  double M   = E*E-P.Mag2();
  if(M>0) return sqrt(M);
  else    return 0;
}
