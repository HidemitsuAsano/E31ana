#include <iostream>
#include <fstream>

#include "TMinuit.h"

#define DISPLAY 0
#define DISPLAY2 0


class ConfMan;
class Display;
class CDSHitMan;
class CDCHit;
class CDSTrackingMan;
class CDSTrack;
class CircleFit;
class HelixFit;

//######## Calc XT #### ///

void CalcXT( )
{

  
  int runnum = 3060;
  int run[] = {1,2,3,4};
  //int run[] = {1,2,3,4,5,6};
  int nrun = sizeof(run)/sizeof(int);



  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");


  /*** conf file for new parameters ***/
  //  ConfMan *conf = new ConfMan("conf/Oct2010/analyzer.conf");
  ConfMan *conf = new ConfMan("conf/Oct2010/analyzer.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */


  TFile *f[10];
  TTree *evtree[10];
  for( int i=0; i<nrun; i++ ){
    f[i] = new TFile( Form( "./root/track_%d_%d.root", runnum, run[i] ) );
    evtree[i] = (TTree*)f[i]->Get( "EventTree" );
  }


  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;
  CDSHitMan *cdsMan = 0;
  EventHeader *head = 0;
  CDSTrackingMan *trackMan=0;

  TFile *rootof=new TFile("reslout3060.root","recreate");  
  TH2F *h_dt_dl[15];
  for(int i=0;i<15;i++)
    h_dt_dl[i]=new TH2F(Form("h_dt_dl%d",i+1),Form("dt vs dl in layer%d",i+1) ,700,-50.,300.0,120,-0.2,1.0);
  TH2F *h_dt_resi[15];
  for(int i=0;i<15;i++)
    h_dt_resi[i]=new TH2F(Form("h_dt_resi%d",i+1),Form("dt vs resi in layer%d",i+1) ,700,-50.,300.0,440,-0.22,0.22);
  TH1F *h_resi[15];
  for(int i=0;i<15;i++)
    h_resi[i]=new TH1F(Form("h_resi%d",i+1),Form("resi in layer%d",i+1) ,3000,
		       -0.15,0.15);


    double eff[15]={0};
    int ntrack=0;


  for( int irun=0; irun<nrun; irun++ ){
    evtree[irun]->SetBranchAddress( "CDSHitMan", &cdsMan );
    //    evtree[irun]->SetBranchAddress( "BeamLineHitMan", &blMan );
    evtree[irun]->SetBranchAddress( "EventHeader", &head );
    evtree[irun]->SetBranchAddress( "CDSTrackingMan", &trackMan );
    int nev = evtree[irun]->GetEntries();
    std::cout << " irun:" << irun << " entry:" << nev << std::endl;
    for( int iev=0; iev<nev; iev++ ){
      if( iev%100 == 0 )
	std::cout << " Event[" << irun << "] : " << iev << std::endl;
      evtree[irun]->GetEvent(iev);


    //    if(iev>10000 ) break;    
    /*** if some parameters are newed, you can re-calc all quantities as below, ***/
    //cdsMan->Calc(conf);

    double param[5];
    for(int i=0;i<trackMan->nTrack();i++)
      {
	CDSTrack *track=trackMan->Track(i);
	track->GetParameters(param);
	if(fabs(param[2])<1e-100 ) continue;	
	if(track->Chi()>100 ) continue;	
	double rho=fabs(1./param[2]);
	double x_c=(param[0]+1./param[2])*cos(param[1]);
	double y_c=(param[0]+1./param[2])*sin(param[1]);

	for(int numlayer=1;numlayer<=15;numlayer++)
	  {
	    for(int n=0;n<track->nTrackHit(numlayer);n++)
	      {
		double dt=track->TrackHit(numlayer,n)->dt();
		TVector3 hpos,wpos,wposp;
		double dl_c=track->TrackHit(numlayer,n)->dl();		
		
		hpos.SetXYZ(track->TrackHit(numlayer,n)->x(),
			track->TrackHit(numlayer,n)->y(),
			    track->TrackHit(numlayer,n)->z() );
		wpos.SetXYZ(track->TrackHit(numlayer,n)->wx(),
			    track->TrackHit(numlayer,n)->wy(),
			    track->TrackHit(numlayer,n)->wz() );
		wposp.SetXYZ(track->TrackHit(numlayer,n)->wxp(),
			     track->TrackHit(numlayer,n)->wyp(),
			     track->TrackHit(numlayer,n)->wzp() );
		
		double dl_r,resi,xest,yest,distmp;
		if( (1<=numlayer && numlayer<=3) || (8<=numlayer && numlayer<=9) || (14<=numlayer && numlayer<=15) ) 
		  {
		    track->PointToCircle(hpos.x(),hpos.y(),rho,x_c,y_c,resi,xest,yest);
		    track->PointToCircle(wpos.x(),wpos.y(),rho,x_c,y_c,dl_r,xest,yest);
		    
		    int sign2=0;
		    if(dl_c<dl_r) sign2=1;
		    else sign2=-1;
		    resi=sign2*resi;
		    h_resi[numlayer-1]->Fill(resi);
		    if(dl_r < 1. && n==0) eff[numlayer-1]++;
		  }	
		else if( (4<=numlayer && numlayer<=7) || (10<=numlayer && numlayer<=13) ) 
		  {
		    
		    TVector3 hnest;
		    track->PointToHelix(hpos,param,hnest,resi);	    
		    TVector3 whpos,dline;
		    dline=wposp-wpos;
		    dline=dline.Unit();
		    double k=hpos*dline-wpos*dline;
		    whpos=wpos+k*dline;
		    TVector3 lnest;
		    if(!track->LineToHelix(whpos,dline,param,lnest,hnest,dl_r) )
		      {std::cout<<"Miss LineToHelix!!"<<std::endl;continue;}
		    resi=dl_r-dl_c;
		    h_resi[numlayer-1]->Fill(resi);
		    if(dl_r < 1. && n==0) eff[numlayer-1]++;
		  }	
		h_dt_dl[numlayer-1]->Fill(dt,dl_r);
		h_dt_resi[numlayer-1]->Fill(dt,resi);
	      }

	  }
	ntrack++;
      }
  }
  }
  
  cout<<"###########eff#########"<<endl;
  for(int layer=1;layer<=15;layer++) 
    cout<<"layer "<<layer<<eff[layer-1]/ntrack<<" "<<eff[layer-1]<<" "<<ntrack<<endl;

  //########  Calc XT ###############//
  
  int itanum=0;
  
  TF1 *resl_ita=new TF1("resl_ita","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x**4+[5]*x**5");
  //   TF1 *resl_ita=new TF1("resl_ita","pol4",-100,800.);
  resl_ita->SetLineColor(2);
    
  //  resl_ita->FixParameter(4,0);
  //resl_ita->FixParameter(4,0);
  
  TCanvas *c1=new TCanvas();
  c1->Divide(5,3);
  ofstream of;
  of.open(Form("XTCurveCDSita%d.param",itanum));
  
  of<<"#Resl itaration resl=p0+p1*dt+p2*dt^2+p3*dt^3+p4*dt^4+p5*dt^5"
    <<std::endl;
  of<<"#cid layer wire np  p0   p1   p2   p3   p4  p5"<<std::endl;
  std::cout<<"#Resl itaration resl=p0+p1*dt+p2*dt^2+p3*dt^3+p4*dt^4+p5*dt^5"
	   <<std::endl;
  std::cout<<"#cid layer wire np  p0   p1   p2   p3   p4 p5"<<std::endl;
    
  for(int i=0;i<15;i++)
    {
      resl_ita->SetParameters(0,0,0,0,0);
      //      resl_ita->FixParameter(3,0);
      //  resl_ita->SetParLimits(0,-0.05,0.05);
      c1->cd(i+1);
      
      //	  TH2F *h2d=(TH2F*)f->Get(Form("h_resl_l%dw%d",i+1,wire+1) );
      //std::cout<<"ok"<<std::endl;  
      
      h_dt_dl[i]->Draw("colz");
      // h2d->RebinX(3);
      TProfile *prof=h_dt_dl[i]->ProfileX();
      
      //prof->GetXaxis()->SetRangeUser(-50,300);
      //prof->GetYaxis()->SetRangeUser(-0.2,0.2);
      prof->Fit("resl_ita","W","",0,240);
      //      if(h_dt_dl[i]->GetEntries()>200 )prof->Fit("resl_ita","I","",0,240);
      double par[6];
      resl_ita->GetParameters(par);
      //pol5->GetParameters(par);
      double chi2,dof;
      dof=resl_ita->GetNDF();
      chi2=resl_ita->GetChisquare();
      
      of<<"0 "<<i+1<<" 0"<<" "<<"6 "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<par[4]<<" "<<par[5]<<std::endl;
      std::cout<<"0 "<<i+1<<" 0" <<" "<<"6 "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<par[4]<<" "<<par[5]<<std::endl;
      
      //      delete h2d;
      //      delete prof;
    }
  
  c1->SaveAs(Form("Reslita%d.eps",itanum));
  
  rootof->Write();
  rootof->Close();  
  
}

