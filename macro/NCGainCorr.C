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

void NCGainCorr()
{


  //Layer 
  int Layer[2]={7,7};
  int sseg[2][4]={1,2,3,4,5,6,7,8};
  int hvnum=0;// -80, -20, 0 ,20, 80 
  /*** assign input file & call tree ***/
  int runnum=39;
  TFile *f = new TFile(Form("root/NCcalib/NCr000%d.root",runnum));
  int cr = 7; // crate number
  int slall[] = {23};
  int ud=1;//ud=0 ch0~15 ud=1 ch16~17

  int nsl = sizeof(slall)/sizeof(int);


  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Run40/analyzerNCcalib2.conf");
  conf->Initialize();


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
