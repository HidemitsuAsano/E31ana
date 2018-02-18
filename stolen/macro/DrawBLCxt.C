#include <iostream>
#include <fstream>

#include "TMinuit.h"

#define DISPLAY 0
#define DISPLAY2 0

class GlobalVariables;
class ConfMan;
class Display;
class CDSHitMan;
class CDCHit;
class CDSTrackingMan;
class CDSTrack;
class CircleFit;
class HelixFit;

//######## Calc XT #### ///

void DrawBLCxt()
{

  //#####################

  int itanum=3;  
  int runnum = 3036;

  //#####################


  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");


  /*** conf file for new parameters ***/
  //  ConfMan *conf = new ConfMan("conf/Oct2010/analyzer.conf");
  ConfMan *conf = new ConfMan("conf/Oct2010/analyzerBLC.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */



  TFile *f;

  
  f = new TFile( Form( "./root/xtout_blc%d_%d.root",itanum, runnum ) );

  //########  Calc XT ###############//
  
  int layer=5;
  int ud=1;

  TCanvas *c1=new TCanvas();  
  c1->Divide(4,4);
  for(int blcnum=1;blcnum<=2;blcnum++)
    {
      for(int layer=1;layer<=8;layer++)
	{
	  if(blcnum==1 &&layer>5) continue;
	  for(int ud=0;ud<2;ud++)
	    {	      
	      for(int n=0;n<16;n++)
		{
		  c1->cd(n+1);
		  int wire=1+n+ud*16;
		  TH1F *h1=(TH1F*)f->Get(Form("BLC%ddt%d_%d",blcnum,layer,wire) );
		  h1->Draw();
		}
	      if(blcnum==1 && layer==1 &&ud==0)     c1->SaveAs("BLCdt.ps(");
	      else if(blcnum==2 && layer==8 &&ud==1)     c1->SaveAs("BLCdt.ps)");
	      else    c1->SaveAs("BLCdt.ps");
	    }
	}
    }
  //######dt_dl2#####

  for(int blc=1;blc<=2;blc++)
    {
      for(int layer=1;layer<=8;layer++)
	{
	  if(blc==1 &&layer>5 ) continue;
	  for(int ud=0;ud<2;ud++)
	    {
	      
	      for(int n=0;n<16;n++)
		{
		  c1->cd(n+1);
		  int wire=1+n+ud*16;
		  TH2F *h2=(TH2F*)f->Get(Form("BLC%ddt_dl2_%d_%d",blc,layer,wire) );
		  h2->Draw("colz");
		}
	      if(blc==1 && layer==1 &&ud==0)     c1->SaveAs("BLCdtdl.ps(");
	      else if(blc==2&&layer==8 &&ud==1)     c1->SaveAs("BLCdtdl.ps)");
	      else    c1->SaveAs("BLCdtdl.ps");
	    }
	}
    }

  //######dt_resid#####

  for(int blc=1;blc<=2;blc++)
    {
      for(int layer=1;layer<=8;layer++)
	{
	  if(blc==1 &&layer>5 ) continue;
	  for(int ud=0;ud<2;ud++)
	    {
	      
	      for(int n=0;n<16;n++)
		{
		  c1->cd(n+1);
		  int wire=1+n+ud*16;
		  TH2F *h2=(TH2F*)f->Get(Form("BLC%ddt_resid%d_%d",blc,layer,wire) );
		  h2->Draw("colz");
		}
	      if(blc==1 && layer==1 &&ud==0)     c1->SaveAs("BLCdtresid.ps(");
	      else if(blc==2&&layer==8 &&ud==1)     c1->SaveAs("BLCdtresid.ps)");
	      else    c1->SaveAs("BLCdtresid.ps");
	    }
	}
    }


  //######dt_resid#####

  TF1 *dgaus=new TF1("dgaus","gaus(0)+gaus(3)");
  TF1 *gaus1=new TF1("gaus1","gaus(0)");
  gaus1->SetParameters(100,0,0.01);
  dgaus->SetLineColor(2);
  double nlayer[16]={0};
  double resid[16]={0};
  double ex[16]={0};
  double ey[16]={0}; 

  for(int blc=1;blc<=2;blc++)
    {
      for(int layer=1;layer<=8;layer++)
	{
	  int n=(blc-1)*8+layer;
	  if(blc==1 &&layer>5 )
	    {
	      nlayer[n-1]=n;
	      resid[n-1]=0.01;
	      ey[n-1]=0.0;
	      continue;
	    }
	  c1->cd(n);
	  TH1F *h1=(TH1F*)f->Get(Form("BLC%dresid%d",blc,layer) );
	  h1->Fit("gaus1");
	  double par[3],apar[6]; 
	  gaus1->GetParameters(par);
	  dgaus->SetParameters(par[0],par[1],par[2]/2.0,par[0]/20.0,par[1],par[2]*2);
 	  h1->Fit("dgaus");
	  dgaus->GetParameters(apar);
	  nlayer[n-1]=n;
	  resid[n-1]=apar[2];
	  ey[n-1]=dgaus->GetParError(2);

	}

    }
  c1->SaveAs("BLCresid.ps(");
  TGraphErrors *g1=new TGraphErrors(16,nlayer,resid,ex,ey);
  g1->SetLineColor(2);

  TCanvas *c2=new TCanvas();
  g1->Draw("ALP");
  c2->SaveAs("BLCresid.ps)");
  
}
