#include <iostream>
#include <fstream>
#include <iomanip.h>
#include "TMinuit.h"

#define DISPLAY 0
#define DISPLAY2 0

class GlobalVariables;
class ConfMan;

#define CDH 0
#define T0 0
#define NC 1

#define SlewType 1

void PlotSlewNC3(std::string InFileName)
{

  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");
  //########parameter file
  gStyle->SetOptFit();
#define MAXCHAR 256   
  ifstream fp;
  //  string InFileName="SlewNC.param";
  char str[MAXCHAR],str1[MAXCHAR];

  int runnum = 83;
  int Layer[2]={6,7};
  int sseg[2][4]={9,10,11,12,9,10,11,12};
  int ita_ini=15;
  int ita_fin=16;
  string conffile;
  string slewfile;
  int writeparam=0;
  std::cout << InFileName <<endl;
  if( fp.open(InFileName.c_str() )==0 )
    {
      std::cerr << " File open fail. [" << InFileName << "]" << std::endl;
      return;
    }

  while( fp.getline(str,MAXCHAR) )
    {

      if( str[0]=='#' )
	{
	  //	  ofp<<str<<std::endl;
	} 

      if( sscanf(str,"Run: %d", &runnum)==1 ) 
	std::cout<<"#Run= "<<runnum<<std::endl;
      else if( sscanf(str,"Layer: %d %d", &Layer[0],&Layer[1])==2 ) 
      	std::cout<<"Layer[0]= "<<Layer[0]<<" Layer[1]= "<<Layer[1]<<std::endl;
      else if( sscanf(str,"segment[0]: %d %d %d %d", &sseg[0][0], &sseg[0][1], &sseg[0][2], &sseg[0][3])==4 )
      	std::cout<<"sseg[0]= "<<sseg[0][0]<<" "<<sseg[0][1]<<" "<<sseg[0][2]<<" "<<sseg[0][3]<<std::endl;
      else if( sscanf(str,"segment[1]: %d %d %d %d", &sseg[1][0], &sseg[1][1], &sseg[1][2], &sseg[1][3])==4 )
      	std::cout<<"sseg[1]= "<<sseg[1][0]<<" "<<sseg[1][1]<<" "<<sseg[1][2]<<" "<<sseg[1][3]<<std::endl;
      else if( sscanf(str,"itanum: %d %d", &ita_ini, &ita_fin)==2 )
	std::cout<<"itanum= "<<ita_ini<<" "<<ita_fin<<std::endl;
      else if( sscanf(str,"Conffile: %s", str1)==1 )
	{
	  conffile=str1;	
	  std::cout<<"conffile= "<<conffile<<std::endl;
	}
      else if( sscanf(str,"Slewfile: %s", str1)==1 )
	{
	  slewfile=str1;	
	  std::cout<<"slewfile= "<<slewfile<<std::endl;
	}
      else if( sscanf(str,"WriteParam: %d", &writeparam)==1 )
	std::cout<<"writeparam= "<<writeparam<<std::endl;
    }


  int rseg[2][4];
  for(int iset=0;iset<2;iset++){
    for(int iseg=0;iseg<4;iseg++){
      rseg[iset][iseg]=(Layer[iset]-1)*16+sseg[iset][iseg];
      cout<<"rseg "<<iset <<" "<<iseg<<" "<<rseg[iset][iseg]<<endl;;
      }
    }

  int itanum=ita_fin;
  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  TFile *f;
  f = new TFile( Form( "./root/NCcalib/NCr000%d_ita%d.root",runnum,itanum) );


  TH1F *h1;
  TH2F *h2;
  TProfile *prof;
  TCanvas *c1;
  for(int iset=0;iset<2;iset++){
    int iseg12=0;
    for(int iseg1=0;iseg1<4;iseg1++){
      for(int iseg2=iseg1+1;iseg2<4;iseg2++){      
	c1=new TCanvas();
	c1->Divide(3,2);
	iseg12++;
	int seg1=sseg[iset][iseg1]+(Layer[iset]-1)*16;
	int seg2=sseg[iset][iseg2]+(Layer[iset]-1)*16;		
	for(int ud=0;ud<2;ud++){      
	  c1->cd(1+ud);
	  if(ud==0)
	    h2 = (TH2F*)f->Get( Form("NCU%d_E_TOF%d_%d",seg1,seg1,seg2) );
	  else 
	    h2 = (TH2F*)f->Get( Form("NCD%d_E_TOF%d_%d",seg1,seg1,seg2) );
	  h2->GetXaxis()->SetRangeUser(-50,300);
	  h2->GetYaxis()->SetRangeUser(-2,2);
	  h2->Draw("colz");
	  if(ud==0)
	    prof = (TProfile*)f->Get( Form("NCU%d_E_TOF%d_%d_pfx",seg1,seg1,seg2) );
	  else 
	    prof = (TProfile*)f->Get( Form("NCD%d_E_TOF%d_%d_pfx",seg1,seg1,seg2) );
	  prof->GetFunction("slewing")->Draw("same");

	  
	  if(ud==0)
	    h2 = (TH2F*)f->Get( Form("NCU%d_E_TOF%d_%d",seg2,seg1,seg2) );
	  else 
	    h2 = (TH2F*)f->Get( Form("NCD%d_E_TOF%d_%d",seg2,seg1,seg2) );
	  c1->cd(4+ud);
	  h2->GetXaxis()->SetRangeUser(-50,300);
	  h2->GetYaxis()->SetRangeUser(-2,2);
	  h2->Draw("colz");
	  if(ud==0)
	    prof = (TProfile*)f->Get( Form("NCU%d_E_TOF%d_%d_pfx",seg2,seg1,seg2) );
	  else 
	    prof = (TProfile*)f->Get( Form("NCD%d_E_TOF%d_%d_pfx",seg2,seg1,seg2) );
	  prof->GetFunction("slewing")->Draw("same");

	  if(ud==0) continue;
	  c1->cd(3);
	  h1 = (TH1F*)f->Get( Form("NCD%d_E_TOF%d_%d_py",seg2,seg1,seg2) );
	  h1->GetXaxis()->SetRangeUser(-2,2);
	  h1->Draw();
	  h1->GetFunction("gaus")->Draw("same");
	}
	if(iseg12==1)c1->Print(Form("NC_slew_run%d_ita%d.ps(",runnum,itanum));
	else c1->Print(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
      }
    }
  }
  c1=new TCanvas();
  c1->Divide(4,2);
  
  int isetseg=0;
  TGraphErrors *g1;
  TF1* f1;
  for(int iset=0;iset<2;iset++)
    for(int iseg=0;iseg<4;iseg++)
      {
	isetseg++;

	g1=(TGraphErrors*)f->Get(Form("layer%dSeg%d",Layer[iset],sseg[iset][iseg]));
	c1->cd(isetseg);
	g1->GetYaxis()->SetRangeUser(50,150);
	g1->Draw("AP");
	f1=g1->GetFunction("pol0");
	f1->SetLineColor(2);
	f1->Draw("same");
	double fitsigma=f1->GetParameter(0);
	double err_sigma=f1->GetParError(0);
	TLatex latex;
	latex.SetTextColor(2);
	latex.SetTextSize(0.08);
	latex.DrawLatex( 1, 130, Form("Layer%d Seg%d ", Layer[iset],sseg[iset][iseg]) );
	latex.DrawLatex( 1, 120, Form("#sigma = %2.2lf +- %1.2lf", fitsigma,err_sigma) ); 
      }

  c1->Print(Form("NC_slew_run%d_ita%d.ps)",runnum,itanum));
}
