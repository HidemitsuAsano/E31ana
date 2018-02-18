#include <fstream.h>
#include <iostream.h>
#include <stdio.h>
#include <vector.h>
#include <math.h>
#include "TMinuit.h"
#include "TMath.h"

#define Crate_PDC 0
#define Crate_BLC 1
#define Crate_BLHodo 2
#define Crate_CDC1 3
#define Crate_CDC2 4
#define Crate_CDC3 5
#define Crate_CDSHodo 6

const int Layer[2]={6,6};
//int sseg[2][4]={9,10,11,12,13,14,15,16};
//const int sseg[2][4]={1,2,3,4,5,6,7,8};
const int sseg[2][4]={1,2,3,4,9,10,11,12};
//layer:7,7,seg:1-8
//const double hv[2][4][2]={1510,1410,1440,1400,1480,1450,1420,1480,
//			  1520,1560,1510,1440,1480,1540,1550,1460};
//layer:7,7,seg:9-16
//const double hv[2][4][2]={1470,1750,1580,1510,1330,1460,1450,1450,
//			  1450,1530,1620,1410,1350,1470,1710,1550};
//layer:6,6,seg:1-4,9-13
const double hv[2][4][2]={1630,1360,1530,1520,1400,1600,1430,1590,
			  1660,1590,1420,1480,1380,1510,1320,1440};

const int hv2[5]={0,20,-20,80,-80};
const int startrun=47;

void NCGainCorr(int runnum,  int hvnum);// -80, -20, 0 ,20, 80 )
void CalcNCGainCurve();
void NCHVcalib()
{
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");
  for(int i=0;i<5;i++)
    NCGainCorr(startrun+i,hv2[i]);
  CalcNCGainCurve();
}

void NCGainCorr(int runnum,  int hvnum)// -80, -20, 0 ,20, 80 )
{
  /*** assign input file & call tree ***/
  TFile *f = new TFile(Form("root/NCcalib/NCr000%d.root",runnum));
  int cr = 7; // crate number
  int slall[] = {23};
  int ud=1;//ud=0 ch0~15 ud=1 ch16~17

  int nsl = sizeof(slall)/sizeof(int);

  int seg[2][4];
  for(int i=0;i<4;i++)
    {
      seg[0][i]=(Layer[0]-1)*16+sseg[0][i];
      seg[1][i]=(Layer[1]-1)*16+sseg[1][i];
    }
  // TFile *f = new TFile( filename );

  TH1F *h1;
  TH2F *h2;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  TCanvas *c2 = new TCanvas( "c2", "c2", 800, 600 );
  c1->Divide(4,2);
  c2->Divide(4,2);
  double ped[2][4][2];
  double ped_err[2][4][2];
  double mip[2][4][2];
  double mip_err[2][4][2];

  for( int iset=1;iset<=2; iset++ ){
  for( int iseg=0; iseg<4; iseg++ ){
   if(iset==1) c1->cd(iseg+1);
   else if(iset==2) c2->cd(iseg+1);   
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ENCU%d",seg[iset-1][iseg]) );
    //   h1->SetAxisRange(-100,600,"X");
    h1->Rebin(4);
    h1->Draw();
    h1->GetXaxis()->SetRangeUser(-100,600);
    //    if(iset==1 && iseg==0)    h1->Fit("gaus","q","SAME",-25,20);

    h1->Fit("gaus","q","SAME",-25,20);
    ped[iset-1][iseg][0] = h1->GetFunction( "gaus" )->GetParameter(1);
    ped_err[iset-1][iseg][0] = h1->GetFunction( "gaus" )->GetParError(1);
    double rchi = h1->GetFunction( "gaus" )->GetChisquare()/h1->GetFunction( "gaus" )->GetNDF();;
    if(rchi >300 || ped_err[iset-1][iseg][0] >5)
      {
	std::cout<<"Error!! ped set seg ud "<<iset<<" "<<iseg+1<<" "<<0<<" rchi "<<rchi<<" err "<<ped_err[iset-1][iseg][0]<<std::endl;
      }

    if(iset==1 && iseg==0)      h1->Fit("landau","q","SAME",20,180);
    landau->SetLineColor(2);
    h1->Fit("landau","q","SAME",20,180);
    mip[iset-1][iseg][0] = h1->GetFunction( "landau" )->GetParameter(1);
    mip_err[iset-1][iseg][0] = h1->GetFunction( "landau" )->GetParError(1);
    rchi = h1->GetFunction( "landau" )->GetChisquare()/h1->GetFunction( "landau" )->GetNDF();;
    if(rchi >10 || mip_err[iset-1][iseg][0] >5)
      {
	std::cout<<"Error!! landau set seg ud "<<iset<<" "<<iseg+1<<" "<<0<<" rchi "<<rchi<<std::endl;
      }
    gaus->SetLineColor(4);
    gaus->SetRange(-25,25);
    gaus->Draw("SAME");
    
    if(iset==1) c1->cd(iseg+5);
    else if(iset==2) c2->cd(iseg+5);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ENCD%d",seg[iset-1][iseg]) );
    //    h1->SetAxisRange(0,600,"X");
    h1->Rebin(4);
    h1->Draw();
    h1->GetXaxis()->SetRangeUser(-100,600);
    h1->Fit("gaus","q","",-25,20);
    ped[iset-1][iseg][1] = h1->GetFunction( "gaus" )->GetParameter(1);
    ped_err[iset-1][iseg][1] = h1->GetFunction( "gaus" )->GetParError(1);
    rchi = h1->GetFunction( "gaus" )->GetChisquare()/h1->GetFunction( "gaus" )->GetNDF();;
    if(rchi >300 || ped_err[iset-1][iseg][1] >5)
      {
	std::cout<<"Error!! ped set seg ud "<<iset<<" "<<iseg+1<<" "<<1<<" rchi "<<rchi<<" err "<<ped_err[iset-1][iseg][1]<<std::endl;
      }

    h1->Fit("landau","q","SAME",20,180);
    mip[iset-1][iseg][1] = h1->GetFunction( "landau" )->GetParameter(1);
    mip_err[iset-1][iseg][1] = h1->GetFunction( "landau" )->GetParError(1);
    rchi = h1->GetFunction( "landau" )->GetChisquare()/h1->GetFunction( "landau" )->GetNDF();;
    if(rchi >10 || mip_err[iset-1][iseg][1] >5)
      {
	std::cout<<"Error!! landau set seg ud "<<iset<<" "<<iseg+1<<" "<<1<<" rchi "<<rchi<<std::endl;
      }
    gaus->SetLineColor(4);
    gaus->SetRange(-25,25);
    gaus->Draw("SAME");


  }
  }
  
  c1->SaveAs(Form("NCGainCorr_r000%d.ps(",runnum));
  c2->SaveAs(Form("NCGainCorr_r000%d.ps)",runnum));
  double gain[2][4][2];
  double gain_err[2][4][2];

  ofstream ofs(Form("NCgain_lay%d_fseg%d_lay%d_fseg%d_hv%d.dat",Layer[0],sseg[0][0],Layer[1],sseg[1][0],hvnum ));

  for( int iset=1;iset<=2; iset++ ){
    for( int iseg=0; iseg<4; iseg++ ){
      gain[iset-1][iseg][0]=mip[iset-1][iseg][0]-ped[iset-1][iseg][0];
      gain[iset-1][iseg][1]=mip[iset-1][iseg][1]-ped[iset-1][iseg][1];
      gain_err[iset-1][iseg][0]=sqrt( pow(mip_err[iset-1][iseg][0],2)+pow(ped_err[iset-1][iseg][0],2) );
      gain_err[iset-1][iseg][1]=sqrt( pow(mip_err[iset-1][iseg][1],2)+pow(ped_err[iset-1][iseg][1],2) );


      std::cout<<"Set seg "<<iset<<" "<<iseg<<" u gain "<<gain[iset-1][iseg][0]<<" +- "<<gain_err[iset-1][iseg][0]<<" ( "<<  mip[iset-1][iseg][0]<<" - " <<ped[iset-1][iseg][0]<<" )"<<std::endl;
      std::cout<<"Set seg "<<iset<<" "<<iseg<<" d gain "<<gain[iset-1][iseg][1]<<" +- "<<gain_err[iset-1][iseg][1]<<" ( "<<  mip[iset-1][iseg][1]<<" - " <<ped[iset-1][iseg][1]<<" )"<<std::endl;
      ofs<<Layer[iset-1]<<" "<<sseg[iset-1][iseg]<<" "<<0<<" "
	 <<gain[iset-1][iseg][0]<<" "<<gain_err[iset-1][iseg][0]<<" "
	 << mip[iset-1][iseg][0]<<" "<< mip_err[iset-1][iseg][0]<<" " 
	 <<ped[iset-1][iseg][0]<<"  "<<ped_err[iset-1][iseg][0]<<std::endl;
      ofs<<Layer[iset-1]<<" "<<sseg[iset-1][iseg]<<" "<<1<<" "
	 <<gain[iset-1][iseg][1]<<" "<<gain_err[iset-1][iseg][1]<<" "
	 <<  mip[iset-1][iseg][1]<<" "<<  mip_err[iset-1][iseg][1]<<" " 
	 <<ped[iset-1][iseg][1]	 <<" "<<ped_err[iset-1][iseg][1]<<std::endl;
    }
  }



  //  std::cout<<"gain u"<<gain[iset][iseg][0]<<" d "<<  gain[iset][iseg][01]<<std::endl;
  ofs.close();
  //gFile->Write();
  //gFile->Close();
}
void CalcNCGainCurve()
{
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
 
   /* --------------- */
  /* ---- start ---- */
  /* --------------- */
  double  gain[2][4][2][5];//layer,subseg,ud,hvnum
  double  gain_err[2][4][2][5];//layer,subseg,ud,hvnum
  double ihv[2][4][2][5];

  for(int ilay=0;ilay<2;ilay++){
    for(int iseg=0;iseg<4;iseg++){
      for(int ud=0;ud<2;ud++){
	for(int i=0;i<5;i++)
	  ihv[ilay][iseg][ud][i]=hv[ilay][iseg][ud]+hv2[i];
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
    InFileName[i]=Form("./NCgain_lay%d_fseg%d_lay%d_fseg%d_hv%d.dat",Layer[0],sseg[0][0],Layer[1],sseg[1][0],hv2[i]);

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
	  std::cout<<"ilay,iseg,ud,i,gain,= "<<ilay<<","<<iseg<<","<<ud<<","<<i<<","<<gain[ilay][iseg][ud][i]<<std::endl;
	}
      }
    }

    std::cout <<"test1"<<endl;
    TFile *f;
    //    f = new TFile( Form("./root/NCcalib/NCgain_lay%d_fseg%d_lay%d_fseg%d.root",Layer[0],sseg[0][0],Layer[1],sseg[1][0]),"recreate" );
    f = new TFile( Form("./NCgain_lay%d_fseg%d_lay%d_fseg%d.root",Layer[0],sseg[0][0],Layer[1],sseg[1][0]),"recreate" );
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
	  ge[ilay][iseg][ud]->SetTitle(Form("NC_lay%d_seg%d_%d_Gain",Layer[ilay],sseg[ilay][iseg],ud));
	  TF1 *f1=new TF1("f1","pol0(0)+expo(1)");
	  f1->SetLineWidth(0.2);
	  ge[ilay][iseg][ud]->Fit("f1","0q","",1000,2000);
	  rhv[ilay][iseg][ud]=ge[ilay][iseg][ud]->GetFunction("f1")->GetX(65);
	  //	double chi=ge[ilay][iseg][ud]->GetFunction("f1")->GetChisquare()/ge[ilay][iseg][ud]->GetFunction("f1")->GetNDF();
	  c1->cd(tmp);
	  ge[ilay][iseg][ud]->Draw("AP");
	  ge[ilay][iseg][ud]->GetFunction("f1")->Draw("SAME");
	  
	  //	ge[ilay][iseg][ud]->Fit("f2","q","",ihv[ilay][iseg][ud][1]-10,ihv[ilay][iseg][ud][3]+10);
	  ge[ilay][iseg][ud]->Write();
	  
	  //	double rhvtmp=ge[ilay][iseg][ud]->GetFunction("f2")->GetX(65);
	  std::cout<<"Lay : seg : ud  "<<Layer[ilay]<<" "<<sseg[ilay][iseg]<<" "<<ud<<" HV@65[pC] = "<<	rhv[ilay][iseg][ud]<<std::endl;//" chi "<<chi<<std::endl;
	  
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
