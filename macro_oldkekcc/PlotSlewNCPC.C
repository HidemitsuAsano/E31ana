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
ofstream ofs;
#define SlewType 1
void PlotSlewNCPC(std::string InFileName);
void PlotSlewNCPC(){
  ofs.open("resol.txt");
  /*
  PlotSlewNCPC("SlewNC34.param");
  PlotSlewNCPC("SlewNC61.param");
  PlotSlewNCPC("SlewNC67.param");
  PlotSlewNCPC("SlewNC73.param");
  PlotSlewNCPC("SlewNC82.param");
  PlotSlewNCPC("SlewNC90.param");
  PlotSlewNCPC("SlewNC103.param");
  PlotSlewNCPC("SlewNC104.param");
  PlotSlewNCPC("SlewNC105.param");
  PlotSlewNCPC("SlewNC118.param");
  PlotSlewNCPC("SlewNC124.param");
  PlotSlewNCPC("SlewNC136.param");
  PlotSlewNCPC("SlewNC142.param");
  PlotSlewNCPC("SlewNC148.param");
  PlotSlewNCPC("SlewNC165.param");
  PlotSlewNCPC("SlewNC176.param");
  */
  PlotSlewNCPC("SlewPC170.param");
  PlotSlewNCPC("SlewPC176.param");
  PlotSlewNCPC("SlewPC182.param");
  PlotSlewNCPC("SlewPC188.param");
  PlotSlewNCPC("SlewPC197.param");
  ofs.close();
}
void PlotSlewNCPC(std::string InFileName)
{  
  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");
  gStyle->SetOptFit();
#define MAXCHAR 256   
  ifstream fp;
  char str[MAXCHAR],str1[MAXCHAR];

  int runnum = 83;
  int Layer[2]={6,7};
  int sseg[2][4]={9,10,11,12,9,10,11,12};
  char* CName[2]={"NC","NC"};
  int CID[2]={CID_NC,CID_NC};
  int ita_ini=15;
  int ita_fin=16;
  int nset=2;
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
      if( str[0]=='#' )continue;
      if( sscanf(str,"Run: %d", &runnum)==1 ) 
	std::cout<<"#Run= "<<runnum<<std::endl;
      else if( sscanf(str,"NumberOfSet: %d", &nset)==1 ) 
      	std::cout<<"nset= "<<nset<<std::endl;
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

  for(int i=0;i<2;i++)
    if(Layer[i]<0){
      CID[i]=CID_PC; 
      CName[i]="PC";
    }
  
  int rseg[2][4];
  for(int iset=0;iset<nset;iset++){
    if(CID[iset]==CID_NC)
      for(int iseg=0;iseg<4;iseg++){
	rseg[iset][iseg]=(Layer[iset]-1)*16+sseg[iset][iseg];
	std::cout<<"rseg "<<iset <<" "<<iseg<<" "<<CName[iset]<<" "<<rseg[iset][iseg]<<std::endl;;
      }
    else
      for(int iseg=0;iseg<4;iseg++){
	rseg[iset][iseg]=sseg[iset][iseg];
	std::cout<<"rseg "<<iset <<" "<<iseg<<" "<<CName[iset]<<" "<<rseg[iset][iseg]<<std::endl;;
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
  for(int iset=0;iset<nset;iset++){
    c1=new TCanvas();
    c1->Divide(4,4);
    for(int iseg1=0;iseg1<4;iseg1++){
      int seg1=rseg[iset][iseg1];
      c1->cd(iseg1+1);
      //      gPad->SetLogz();
      char* opt="col";
      std::cout<<iset<<" "<<iseg1<<" "<<Form("%sU%dADCEvNum",CName[iset],seg1)<<std::endl;
      h2=(TH2F*)f->Get( Form("%sU%dADCEvNum",CName[iset],seg1) );
      h2->GetYaxis()->SetRangeUser(0,1000);
      h2->Draw(opt);
      c1->cd(iseg1+5);
      h2=(TH2F*)f->Get( Form("%sD%dADCEvNum",CName[iset],seg1) );
      h2->GetYaxis()->SetRangeUser(0,1000);
      h2->Draw(opt);
      c1->cd(iseg1+9);
      h2=(TH2F*)f->Get( Form("%sU%dTDCEvNum",CName[iset],seg1) );
      double mean=h2->GetMean(2);
      h2->GetYaxis()->SetRangeUser(mean-200,mean+200);
      h2->Draw(opt);
      c1->cd(iseg1+13);
      h2=(TH2F*)f->Get( Form("%sD%dTDCEvNum",CName[iset],seg1) );
      mean=h2->GetMean(2);
      h2->GetYaxis()->SetRangeUser(mean-200,mean+200);
      h2->Draw(opt);
    }
    if(iset==0)c1->Print(Form("NC_slew_run%d_ita%d.ps(",runnum,itanum));
    else c1->Print(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  }

  for(int iset=0;iset<nset;iset++){
    int iseg12=0;
    for(int iseg1=0;iseg1<4;iseg1++){
      for(int iseg2=iseg1+1;iseg2<4;iseg2++){      
	c1=new TCanvas();
	c1->Divide(3,2);
	iseg12++;
	int seg1=rseg[iset][iseg1];
	int seg2=rseg[iset][iseg2];
	for(int ud=0;ud<2;ud++){      
	  c1->cd(1+ud);
	  if(ud==0)
	    h2 = (TH2F*)f->Get( Form("%sU%d_E_TOF%d_%d",CName[iset],seg1,seg1,seg2) );
	  else 
	    h2 = (TH2F*)f->Get( Form("%sD%d_E_TOF%d_%d",CName[iset],seg1,seg1,seg2) );
	  h2->GetXaxis()->SetRangeUser(0,250);
	  h2->GetYaxis()->SetRangeUser(-1,1);
	  h2->Draw("colz");
	  if(ud==0)
	    prof = (TProfile*)f->Get( Form("%sU%d_E_TOF%d_%d_pfx",CName[iset],seg1,seg1,seg2) );
	  else 
	    prof = (TProfile*)f->Get( Form("%sD%d_E_TOF%d_%d_pfx",CName[iset],seg1,seg1,seg2) );
	  prof->GetFunction("slewing")->Draw("same");

	  
	  if(ud==0)
	    h2 = (TH2F*)f->Get( Form("%sU%d_E_TOF%d_%d",CName[iset],seg2,seg1,seg2) );
	  else 
	    h2 = (TH2F*)f->Get( Form("%sD%d_E_TOF%d_%d",CName[iset],seg2,seg1,seg2) );
	  c1->cd(4+ud);
	  h2->GetXaxis()->SetRangeUser(0,250);
	  h2->GetYaxis()->SetRangeUser(-1,1);
	  h2->Draw("colz");
	  if(ud==0)
	    prof = (TProfile*)f->Get( Form("%sU%d_E_TOF%d_%d_pfx",CName[iset],seg2,seg1,seg2) );
	  else 
	    prof = (TProfile*)f->Get( Form("%sD%d_E_TOF%d_%d_pfx",CName[iset],seg2,seg1,seg2) );
	  prof->GetFunction("slewing")->Draw("same");

	  if(ud==0) continue;
	  c1->cd(3);
	  h1 = (TH1F*)f->Get( Form("%s_TOF%d_%d",CName[iset],seg1,seg2) );
	  h1->GetXaxis()->SetRangeUser(-2,2);
	  h1->Draw();
	  h1->GetFunction("gaus")->Draw("same");
	}
	c1->Print(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
      }
    }
  }
  c1=new TCanvas();
  c1->Divide(4,2);
  
  int isetseg=0;
  TGraphErrors *g1;
  TF1* f1;
  for(int iset=0;iset<nset;iset++)
    for(int iseg=0;iseg<4;iseg++)
      {
	isetseg++;

	g1=(TGraphErrors*)f->Get(Form("%slayer%dSeg%d",CName[iset],Layer[iset],sseg[iset][iseg]));
	c1->cd(isetseg);
	g1->GetYaxis()->SetRangeUser(50,150);
	g1->Draw("AP");
	f1=g1->GetFunction("pol0");
	f1->SetLineColor(2);
	f1->Draw("same");
	double fitsigma=f1->GetParameter(0);
	double err_sigma=f1->GetParError(0);
	if(ofs)	ofs<<"Run"<<runnum<<" "<<CName[iset]<<" Layer"<<Layer[iset]<<" Seg"<<sseg[iset][iseg]<<"  : "<<fitsigma<<" +- "<<err_sigma<<std::endl;
	std::cout<<"Run"<<runnum<<" "<<CName[iset]<<" Layer"<<Layer[iset]<<" Seg"<<sseg[iset][iseg]<<"  : "<<fitsigma<<" +- "<<err_sigma<<std::endl;
	TLatex latex;
	latex.SetTextColor(2);
	latex.SetTextSize(0.08);
	latex.DrawLatex( 1, 130, Form("%s Layer%d Seg%d ",CName[iset], Layer[iset],sseg[iset][iseg]) );
	latex.DrawLatex( 1, 120, Form("#sigma = %2.2lf +- %1.2lf", fitsigma,err_sigma) ); 
      }

  c1->Print(Form("NC_slew_run%d_ita%d.ps)",runnum,itanum));
}
