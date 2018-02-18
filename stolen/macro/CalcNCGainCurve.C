#include <iostream>
#include <fstream>
#include <iomanip.h>
#include "TMinuit.h"

#define DISPLAY 0
#define DISPLAY2 0

class GlobalVariables;
class ConfMan;



void CalcNCGainCurve()
{

  //#####################

  int Layer[2]={7,7};
  int sseg[2][4]={1,2,3,4,5,6,7,8};
  double hv[2][4][2]={1510,1410,1440,1400,1480,1450,1420,1480,
		      1520,1560,1510,1440,1480,1540,1550,1460};

  //#####################

  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");


  /*** conf file for new parameters ***/
  //  ConfMan *conf = new ConfMan("conf/Oct2010/analyzer.conf");
  ConfMan *conf = new ConfMan("conf/Run40/analyzerNCcalib2.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */
  double  gain[2][4][2][5];//layer,subseg,ud,hvnum
  double  gain_err[2][4][2][5];//layer,subseg,ud,hvnum
  double ihv[2][4][2][5];

  for(int ilay=0;ilay<2;ilay++){
    for(int iseg=0;iseg<4;iseg++){
      for(int ud=0;ud<2;ud++){
	ihv[ilay][iseg][ud][0]=hv[ilay][iseg][ud]-80;
	ihv[ilay][iseg][ud][1]=hv[ilay][iseg][ud]-20;
	ihv[ilay][iseg][ud][2]=hv[ilay][iseg][ud];
	ihv[ilay][iseg][ud][3]=hv[ilay][iseg][ud]+20;
	ihv[ilay][iseg][ud][4]=hv[ilay][iseg][ud]+80;
      }
    }
  }
  ifstream fp[5];
  ofstream ofp;
  string InFileName[5];
#define MAXCHAR 256   
  char str[MAXCHAR];

  for(int i=0;i<5;i++)
    {
    InFileName[i]=Form("./NCgain_lay%d_fseg%d_lay%d_fseg%d_hv%d.dat",Layer[0],sseg[0][0],Layer[1],sseg[1][0],i);

    std::cout << InFileName[i] <<endl;
    //  cout<<setprecision(6);
    if( fp[i].open(InFileName[i].c_str() )==0 )
      {
	std::cerr << " File open fail. [" << InFileName << "]" << std::endl;
	return;
      }

    while( fp[i].getline(str,MAXCHAR) )
      {
	int ilay=-1;
	int iseg=-1;
	int lay,seg,ud;
	double gaintmp,gainerrtmp, miptmp,miperrtmp, pedtmp,pederrtmp;
	if( str[0]=='#' )
	  {
	    ofp<<str<<std::endl;
	  } 
	else if( (nd=sscanf(str,"%d %d %d %lf %lf %lf %lf %lf %lf", &lay, &seg, &ud, &gaintmp, &gainerrtmp, &miptmp, &miperrtmp, &pedtmp, &pederrtmp) ) ==9 ) 
	{
	  if(lay==Layer[0] ) ilay=0;
	  else  if(lay==Layer[1] ) ilay=1;
	  if(ilay==-1) 
	    {
	      std::cout<<"error! ilay";
		return;
	    }
	  if(seg==sseg[0][0]) iseg=0;
	  else if(seg==sseg[0][1]) iseg=1;
	  else if(seg==sseg[0][2]) iseg=2;
	  else if(seg==sseg[0][3]) iseg=3;
	  if(seg==sseg[1][0]) {ilay=1; iseg=0;} 
	  else if(seg==sseg[1][1]) {ilay=1; iseg=1;} 
	  else if(seg==sseg[1][2]) {ilay=1; iseg=2;} 
	  else if(seg==sseg[1][3]) {ilay=1; iseg=3;} 
	  if(iseg==-1) 
	    {
	      std::cout<<"error! iseg";
		return;
	    }
	  gain[ilay][iseg][ud][i]=gaintmp;
	  gain_err[ilay][iseg][ud][i]=gainerrtmp;
	}
      }
    }

    std::cout <<"test1"<<endl;
  TFile *f;
  f = new TFile( Form("./root/NCcalib/NCgain_lay%d_fseg%d_lay%d_fseg%d.root",Layer[0],sseg[0][0],Layer[1],sseg[1][0]),"recreate" );
    std::cout <<"test2"<<endl;
    TCanvas *c1=new TCanvas();
    c1->Divide(4,4);
  double ex[5]={0,0,0,0,0};
  TGraphErrors *ge[2][4][2];
  double rhv[2][4][2];
  int tmp=0;
  for(int ilay=0;ilay<2;ilay++){
    for(int iseg=0;iseg<4;iseg++){
      for(int ud=0;ud<2;ud++){
	tmp++;
	ge[ilay][iseg][ud]=new TGraphErrors(5,ihv[ilay][iseg][ud],gain[ilay][iseg][ud],ex,gain_err[ilay][iseg][ud]);
	ge[ilay][iseg][ud]->SetLineColor(2);
	//  g1->Draw("ALP");
	ge[ilay][iseg][ud]->SetName(Form("NC_lay%d_seg%d_%d_Gain",Layer[ilay],sseg[ilay][iseg],ud));
	TF1 *f1=new TF1("f1","pol0(0)+expo(1)");
	f1->SetLineWidth(0.2);
	ge[ilay][iseg][ud]->Fit("f1","0q","",1000,2000);
	rhv[ilay][iseg][ud]=ge[ilay][iseg][ud]->GetFunction("f1")->GetX(65);
	double chi=ge[ilay][iseg][ud]->GetFunction("f1")->GetChisquare()/ge[ilay][iseg][ud]->GetFunction("f1")->GetNDF();
	c1->cd(tmp);
	ge[ilay][iseg][ud]->Draw("AP");
	ge[ilay][iseg][ud]->GetFunction("f1")->Draw("SAME");

	//	ge[ilay][iseg][ud]->Fit("f2","q","",ihv[ilay][iseg][ud][1]-10,ihv[ilay][iseg][ud][3]+10);
	ge[ilay][iseg][ud]->Write();

	//	double rhvtmp=ge[ilay][iseg][ud]->GetFunction("f2")->GetX(65);
	std::cout<<"Lay : seg : ud  "<<Layer[ilay]<<" "<<sseg[ilay][iseg]<<" "<<ud<<" HV@65[pC] = "<<	rhv[ilay][iseg][ud]<<" chi "<<chi<<std::endl;

	TLatex latex;
	latex.SetTextColor(2);
	latex.SetTextSize(0.1)
;	latex.DrawLatex( ihv[ilay][iseg][ud][2]-60, 80, Form("HV@65[pC] =%4.2lf ", rhv[ilay][iseg][ud]) );
	latex.DrawLatex( ihv[ilay][iseg][ud][2]-70, 90, Form("Layer%d Seg%d ud=%d ", Layer[ilay],sseg[ilay][iseg],ud) );


	//	std::cout<<"Lay : seg : ud  "<<Layer[ilay]<<" "<<sseg[ilay][iseg]<<" "<<ud<<" HV@65[pC] = "<<	rhvtmp<<" chi "<<chi<<std::endl;
      }
    }
  }
  
  c1->SaveAs( Form("./NCgain_lay%d_fseg%d_lay%d_fseg%d.ps",Layer[0],sseg[0][0],Layer[1],sseg[1][0]) );
  for(int i=0;i<5;i++) fp[i].close();  
  cout<<"test fin"<<endl;
  return;
}
