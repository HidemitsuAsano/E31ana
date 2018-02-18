#include "TMinuit.h"

void drawCDC()
{
  /*** load library ***/
  //  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("tmphist.root");
  //  TTree *evtree = (TTree*)f->Get("EventTree");

  /*** conf file for new parameters ***/
  //  ConfMan *conf = new ConfMan("conf/Oct2010/analyzer.conf");
  // conf->Initialize();

  /*** assign output file ***/
  //TFile *of = new TFile("root/out.root","recreate");
  //  TFile *of = new TFile("root/out1367.root","recreate");

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  TH1F *h1;
  TCanvas *c1=new TCanvas();
  c1->Divide(5,3);
  for(int i=1;i<=15;i++)
    {
      c1->cd(i);
      h1 = (TH1F*)gFile->Get( Form("CDCHitPattern%d",i) );
      h1->Draw();
    }

  TCanvas *c2=new TCanvas();
  c2->Divide(5,3);
  for(int i=1;i<=15;i++)
    {
      c2->cd(i);
      h1 = (TH1F*)gFile->Get( Form("CDCTDC%d",i) );
      h1->Draw();
    }
  c1->SaveAs("cdc.ps(");
  c2->SaveAs("cdc.ps)");
  //  delete c1;


}

void drawCDH()
{

 TFile *f = new TFile("tmphist.root");
 
  TCanvas *c1=new TCanvas();
  c1->Divide(6,3);
  for(int i=1;i<=18;i++)
    {
      c1->cd(i);
      gPad->SetLogy();
      h1 = (TH1F*)gFile->Get( Form("CDHADCU%d",i) );
      //      h1->GetXaxis()->SetRangeUser();
      h1->Draw();
      h1 = (TH1F*)gFile->Get( Form("CDHADCUWT%d",i) );
      h1->SetLineColor(2);
      h1->Draw("same");
    }
  c1->SaveAs("cdh.ps(");
  delete c1;

  TCanvas *c2=new TCanvas();
  c2->Divide(6,3);
  for(int i=19;i<=36;i++)
    {
      c2->cd(i-18);
      gPad->SetLogy();
      h1 = (TH1F*)gFile->Get( Form("CDHADCU%d",i) );
      h1->Draw();
      h1 = (TH1F*)gFile->Get( Form("CDHADCUWT%d",i) );
      h1->SetLineColor(2);
      h1->Draw("same");
    }
  c2->SaveAs("cdh.ps");
  delete c2;


  TCanvas *c3=new TCanvas();
  c3->Divide(6,3);
  for(int i=1;i<=18;i++)
    {
      c3->cd(i);
      gPad->SetLogy();
      h1 = (TH1F*)gFile->Get( Form("CDHADCD%d",i) );
      h1->Draw();
      h1 = (TH1F*)gFile->Get( Form("CDHADCDWT%d",i) );
      h1->SetLineColor(2);
      h1->Draw("same");
    }
  c3->SaveAs("cdh.ps");
  delete c3;


  TCanvas *c4=new TCanvas();
  c4->Divide(6,3);
  for(int i=19;i<=36;i++)
    {
      c4->cd(i-18);
      gPad->SetLogy();
      h1 = (TH1F*)gFile->Get( Form("CDHADCD%d",i) );
      h1->Draw();
      h1 = (TH1F*)gFile->Get( Form("CDHADCDWT%d",i) );
      h1->SetLineColor(2);
      h1->Draw("same");
    }
  c4->SaveAs("cdh.ps");
  delete c4;



  TCanvas *c5=new TCanvas();
  c5->Divide(6,3);
  for(int i=1;i<=18;i++)
    {
      c5->cd(i);
      gPad->SetLogy();
      h1 = (TH1F*)gFile->Get( Form("CDHTDCU%d",i) );
      h1->Draw();
    }
  c5->SaveAs("cdh.ps");
  delete c5;


  TCanvas *c6=new TCanvas();
  c6->Divide(6,3);
  for(int i=19;i<=36;i++)
    {
      c6->cd(i-18);
      gPad->SetLogy();
      h1 = (TH1F*)gFile->Get( Form("CDHTDCU%d",i) );
      h1->Draw();
    }
  c6->SaveAs("cdh.ps");
  delete c6;


  TCanvas *c7=new TCanvas();
  c7->Divide(6,3);
  for(int i=1;i<=18;i++)
    {
      c7->cd(i);
      gPad->SetLogy();
      h1 = (TH1F*)gFile->Get( Form("CDHTDCD%d",i) );
      h1->Draw();
    }
  c7->SaveAs("cdh.ps");
  delete c7;


  TCanvas *c8=new TCanvas();
  c8->Divide(6,3);
  for(int i=19;i<=36;i++)
    {
      c8->cd(i-18);
      gPad->SetLogy();
      h1 = (TH1F*)gFile->Get( Form("CDHTDCD%d",i) );
      h1->Draw();
    }
  c8->SaveAs("cdh.ps)");
  delete c8;


}
