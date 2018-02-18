#include <fstream.h>
#include <iostream.h>
#include <vector.h>
#include <math.h>
#include "TMinuit.h"
#include "TMath.h"

#define DISPLAY 1
#define DISPLAY2 1

class ConfMan;
class Display;
class CDSHitMan;
class CDCHit;
class CDSTrackingMan;
class CDSTrack;
class CircleFit;
class HelixFit;

void DispTrack()
{

  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  gSystem->Load("libPhysics.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("tmp.root");
  TTree *evtree = (TTree*)f->Get("EventTree");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/run-1980/analyzer.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */




#if DISPLAY
  Display *disp = new Display();
  TCanvas *c_disp = new TCanvas( "c_disp", "c_disp", 0, 0, 800, 400 );
  c_disp->Divide(2,1);

    disp->SetCDSFrameXY(-70,70,-70,70);
    disp->SetCDSFrameYZ(-70,70,-70,70);
#endif

  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;
  CDSHitMan *cdsMan = 0;
  EventHeader *head = 0;
  CDSTrackingMan *trackMan=0;
  //  evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "EventHeader", &head );
  evtree->SetBranchAddress( "CDSTrackingMan", &trackMan );
  
  TH1F *h1;
  char dispflag;
  char dispin[100]="";
  TFile *of=new TFile("disout_tmp.root","recreate");
  TH1F *h_jitter= new TH1F("h_jitter","jitter",1000,-5,5);
  TH1F *h_dltheta= new TH1F("h_dltheta","dltheta",1000,-50,50);
  TH1F *h_pt= new TH1F("h_pt","Pt",100,0.,10);
  TH1F *h_chi= new TH1F("h_chi","chi",300,0.,30);
  TH2F *h_dPt= new TH2F("dPt","deltaPt",20,0,2,240,0,1.2);

  /*                */
  /* event analysis */
  /*                */
  int nev = evtree->GetEntries();
  int ntrack=0;
  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);

    if( iev%100==0 )
      std::cout << " Event : " << iev << std::endl;
    //    std::cout << " AllEvent : " << nev << std::endl;
    if(iev%100==0 ) c_disp->Update();    
    //    if(iev>10000 ) break;    
    if(trackMan->nTrack()!=2) continue;
   
    /*** if some parameters are newed, you can re-calc all quantities as below, ***/
    cdsMan->Calc(conf);


#if DISPLAY2
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
    // disp->DrawCDSHitYZ( pad, conf, cdsMan, CID_CDC );

    c_disp->Update();

#endif    


    //####### each fit###############
 
    double arho[2],ax_c[2],ay_c[2],aPt[2];
    double aparam[2][5];

    for(int i=0;i<trackMan->nTrack();i++)
      {
	CDSTrack *track=trackMan->Track(i);
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
	//    std::cout<<"rho= "<<rho<<" x_c= "<<x_c<<" y_c= "<<y_c<<std::endl;
	aPt[i]=0.3*0.5*arho[i]/100.;
	std::cout<<"Pt"<<i<<"= "<<aPt[i]<<"GeV/c"<<std::endl;
	if(aPt[i]>20){flag=false;continue;}
	double Mass =0.105658;
	double beta=sqrt(aPt[i]*aPt[i]/(aPt[i]*aPt[i]+Mass*Mass));	  
	double Chi=track->Chi();
	if(Chi>9999) continue;
	h_chi->Fill(Chi);

	int nhit[15]={0};
	double jitter,theta,y_mid,chi2;
	for(int numlayer=1;numlayer<=15;numlayer++)
	  {
	    for(int n=0;n<track->nTrackHit(numlayer);n++)
	      {
		TVector3 hpos,wpos,wposp;
		double dl=track->TrackHit(numlayer,n)->dl();		
		jitter=track->Jitter();		
		y_mid=track->MidY();		
		theta=track->Theta();		
		chi2=track->Chi();		
		
		hpos.SetXYZ(track->TrackHit(numlayer,n)->x(),
			    track->TrackHit(numlayer,n)->y(),
			    track->TrackHit(numlayer,n)->z() );
		TMarker mhpos;
		mhpos.SetMarkerStyle(i+20);
		mhpos.SetMarkerColor(3+i);
		mhpos.SetMarkerSize(0.5);
		c_disp->cd(2);
		mhpos.DrawMarker(hpos.z(),hpos.y());
		nhit[numlayer-1]++;
	      }

	  }
	/*
	cout<<"TrackHit"<<i<<"  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15"<<endl;
	cout<<"         ";
	for(int n=0;n<15;n++)
	  {
	    cout<<"  "<<nhit[n];
	  }
	cout<<endl;
	*/
	ntrack++;
	cout<<"Track Num "<<ntrack<<" jitter= "<<jitter<<" theta= "<<theta<<" y_mid= "<<y_mid<<" rchi2 "<<chi2<<endl;
      }
 
    double Pt=(aPt[0]+aPt[1])/2.;
    double dPt=fabs(aPt[0]-aPt[1]);
    h_dPt->Fill(Pt,sqrt(2)*dPt/Pt);


    //    if( (aPt[0]>1.0||aPt[1]>1.0) || dPt>0.02) continue;
    //    c_disp->cd(3);
    // h_dPt->Draw("colz");


    /*
    h_dPt->Fill( (aPt[0]+aPt[1])/2. , sqrt(2)*fabs(aPt[0]-aPt[1])/(aPt[0]+aPt[1]));
    h_dPt->Draw("colz");
    */
  
  
#if DISPLAY2 
    pad = c_disp->GetPad(1);
    pad->cd();
    TArc arc1,arc2,arc3;


    arc1.SetFillStyle(0);
    arc1.SetLineColor(2);

    arc2.SetFillStyle(0);
    arc2.SetLineColor(3);

    arc3.SetFillStyle(0);
    arc3.SetLineColor(4);

    //    arc1.DrawArc(x_c,y_c,rho,0,360);
    arc2.DrawArc(ax_c[0],ay_c[0],arho[0],0,360);
    arc3.DrawArc(ax_c[1],ay_c[1],arho[1],0,360);


	TF1 *func_y[3];
	//func_y[0]=new TF1("func_y0","([0])*sin([1])+[2]*(sin([1])-sin( [1]+([3]-x)/([4]*[2]) ) )",param[3]-param[2]*param[4]*(-TMath::Pi()/2. ),param[3]-param[2]*param[4]*(TMath::Pi()/2. ) );
	func_y[1]=new TF1("func_y1","([0])*sin([1])+1./[2]*(sin([1])-sin( [1]+([3]-x)/([4]*1./[2]) ) )",aparam[0][3]-1./aparam[0][2]*aparam[0][4]*(-TMath::Pi()/2. ),aparam[0][3]-1./aparam[0][2]*aparam[0][4]*(TMath::Pi()/2. ) );
	func_y[2]=new TF1("func_y2","([0])*sin([1])+1./[2]*(sin([1])-sin( [1]+([3]-x)/([4]*1./[2]) ) )",aparam[1][3]-1./aparam[1][2]*aparam[1][4]*(-TMath::Pi()/2. ),aparam[1][3]-1./aparam[1][2]*aparam[1][4]*(TMath::Pi()/2. ) );

	//func_y[0]->SetParameters(param[0],param[1],param[2],param[3],param[4]);
	func_y[1]->SetParameters(aparam[0][0],aparam[0][1],aparam[0][2],aparam[0][3],aparam[0][4]);
	func_y[2]->SetParameters(aparam[1][0],aparam[1][1],aparam[1][2],aparam[1][3],aparam[1][4]);
	//func_y[0]->SetLineColor(2);
	func_y[1]->SetLineColor(3);
	func_y[2]->SetLineColor(4);
	func_y[1]->SetLineWidth(1);
	func_y[2]->SetLineWidth(1);
	c_disp->cd(2);
	//func_y[0]->Draw("LPSAME");
	func_y[1]->Draw("LPSAME");
	func_y[2]->Draw("LPSAME");


    c_disp->Update();    


    bool status = disp->Wait();
    if( !status ) return;
#endif



  }
    c_disp->Update();    
    of->Write();
    of->Close();
}

