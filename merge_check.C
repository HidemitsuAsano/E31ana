#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TGClient.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TCanvas.h>

#define NV1290 2

void merge_check(TString tree_file, TString hist_file)
{
  gStyle->SetOptStat(1111111);

  //To summarize efficiency: 
  //1, merge eff, 
  //2, MTDC eff compared to TKO
  double mtdc_t0t_eff[5];
  double mtdc_t0b_eff[5];
  double mtdc_bhd_eff[20];
  double eff[100]={0.0};

  TFile *f_tree = new TFile(tree_file);
  TTree *tree_mtdc = (TTree*)f_tree->Get("tree_mtdc");
  double total_nevents = tree_mtdc->GetEntries();

  TFile *f_hist = new TFile(hist_file);
  TH2F *h_gpio_matching=(TH2F*)f_hist->Get("gpio_matching");
  double matched_nevents = h_gpio_matching->GetEntries();

  //  cout<<"total_nevents = "<<total_nevents<<" matched_nevents = "<<matched_nevents<<endl;
  double merge_eff = matched_nevents/total_nevents;
  //  cout<<"merge_eff = "<<merge_eff<<endl;
  eff[0]=merge_eff;

  int t0t_bin_l[5]={16300, 16350, 16240, 16200, 16350};
  int t0t_bin_h[5]={16500, 16550, 16440, 16400, 16550};
  int t0b_bin_l[5]={16280, 16350, 16250, 16220, 16340};
  int t0b_bin_h[5]={16480, 16550, 16450, 16420, 16540};
  int bhd_bin_l[20]={14000, 14000, 13950, 13900, 13900, 13800, 13850, 13900, 13800, 13850,
		     13800, 13800, 13850, 13900, 13850, 13800, 13850, 13800, 13900, 13900};
  int bhd_bin_h[20]={14400, 14400, 14350, 14300, 14300, 14200, 14250, 14300, 14200, 14250,
		     14200, 14200, 14250, 14300, 14250, 14200, 14250, 14200, 14300, 14300};

  TH1F *h_hodo_nevents=(TH1F*)f_hist->Get("h_hodo_nevents");

  TH1F *h_tmp;
  for(int i=1; i<6; i++){
    h_tmp=(TH1F*)f_hist->Get(Form("TKO_MTDC_T0t%d_eff",i));
    int bin_l = h_tmp->GetXaxis()->FindBin(t0t_bin_l[i-1]);
    int bin_h = h_tmp->GetXaxis()->FindBin(t0t_bin_h[i-1]);
    double counts = 0;
    for(int j=bin_l; j<bin_h+1; j++){
      counts = counts+h_tmp->GetBinContent(j);
      //      cout<<"counts = "<<counts<<endl;
    }
    int bin_hodo = h_hodo_nevents->GetXaxis()->FindBin(i);
    mtdc_t0t_eff[i-1]=counts/(h_hodo_nevents->GetBinContent(bin_hodo));
    eff[i+10]=counts/(h_hodo_nevents->GetBinContent(bin_hodo));
    // cout<<"counts2 = "<<counts<<endl;
    // cout<<"h_hodo_nevents = "<<h_hodo_nevents->GetBinContent(bin_hodo)<<endl;
    // cout<<"mtdc_t0t_eff = "<<mtdc_t0t_eff[i-1]<<endl;

  }

  for(int i=1; i<6; i++){
    h_tmp=(TH1F*)f_hist->Get(Form("TKO_MTDC_T0b%d_eff",i));
    int bin_l = h_tmp->GetXaxis()->FindBin(t0b_bin_l[i-1]);
    int bin_h = h_tmp->GetXaxis()->FindBin(t0b_bin_h[i-1]);
    double counts = 0;
    for(int j=bin_l; j<bin_h+1; j++){
      counts = counts+h_tmp->GetBinContent(j);
      //      cout<<"counts = "<<counts<<endl;
    }
    int bin_hodo = h_hodo_nevents->GetXaxis()->FindBin(i+8);
    mtdc_t0b_eff[i-1]=counts/(h_hodo_nevents->GetBinContent(bin_hodo));
    eff[i+20]=counts/(h_hodo_nevents->GetBinContent(bin_hodo));
    // cout<<"counts2 = "<<counts<<endl;
    // cout<<"h_hodo_nevents = "<<h_hodo_nevents->GetBinContent(bin_hodo)<<endl;
    // cout<<"mtdc_t0b_eff = "<<mtdc_t0b_eff[i-1]<<endl;

  }

  for(int i=1; i<21; i++){
    h_tmp=(TH1F*)f_hist->Get(Form("TKO_MTDC_BHD%d_eff",i));
    int bin_l = h_tmp->GetXaxis()->FindBin(bhd_bin_l[i-1]);
    int bin_h = h_tmp->GetXaxis()->FindBin(bhd_bin_h[i-1]);
    double counts = 0;
    for(int j=bin_l; j<bin_h+1; j++){
      counts = counts+h_tmp->GetBinContent(j);
      //      cout<<"counts = "<<counts<<endl;
    }
    int bin_hodo = h_hodo_nevents->GetXaxis()->FindBin(i+100);
    mtdc_bhd_eff[i-1]=counts/(h_hodo_nevents->GetBinContent(bin_hodo));
    eff[i+30]=counts/(h_hodo_nevents->GetBinContent(bin_hodo));
    // cout<<"counts2 = "<<counts<<endl;
    // cout<<"h_hodo_nevents = "<<h_hodo_nevents->GetBinContent(bin_hodo)<<endl;
    // cout<<"mtdc_t0t_eff = "<<mtdc_bhd_eff[i-1]<<endl;

  }

  TH1F *h_eff=new TH1F("h_eff", "h_eff", 50, -0.5, 49.5);
  int bin_l = h_eff->GetXaxis()->FindBin(0);
  int bin_h = h_eff->GetXaxis()->FindBin(149);
  int eff_bin=0;
  for(int i=bin_l; i<bin_h+1; i++){
    h_eff->SetBinContent(i, eff[eff_bin]);
    //    cout<<"eff="<<eff[eff_bin]<<endl;
    eff_bin++;
  }
  h_eff->SetMaximum(2.); 

  printf("create canvas \n");
  c0 = new TCanvas("c0", "c0", 800, 800);
  h_eff->Draw("TEXT");

  c2 = new TCanvas("c2", "c2", 800, 800);
  c2->Divide(3,4);
  c2->cd(1);
  TKO_MTDC_T0t1->Draw("colz");
  c2->cd(2);
  TKO_MTDC_T0t2->Draw("colz");
  c2->cd(3);
  TKO_MTDC_T0t3->Draw("colz");
  c2->cd(4);
  TKO_MTDC_T0t4->Draw("colz");
  c2->cd(5);
  TKO_MTDC_T0t5->Draw("colz");
  c2->cd(6);

  c2->cd(7);
  TKO_MTDC_T0b1->Draw("colz");
  c2->cd(8);
  TKO_MTDC_T0b2->Draw("colz");
  c2->cd(9);
  TKO_MTDC_T0b3->Draw("colz");
  c2->cd(10);
  TKO_MTDC_T0b4->Draw("colz");
  c2->cd(11);
  TKO_MTDC_T0b5->Draw("colz");
  c2->cd(12);


  c3 = new TCanvas("c3", "c3", 800, 800);
  c3->Divide(4,5);
  c3->cd(1);
  TKO_MTDC_BHD1->Draw("colz");
  c3->cd(2);
  TKO_MTDC_BHD2->Draw("colz");
  c3->cd(3);
  TKO_MTDC_BHD3->Draw("colz");
  c3->cd(4);
  TKO_MTDC_BHD4->Draw("colz");
  c3->cd(5);
  TKO_MTDC_BHD5->Draw("colz");
  c3->cd(6);
  TKO_MTDC_BHD6->Draw("colz");
  c3->cd(7);
  TKO_MTDC_BHD7->Draw("colz");
  c3->cd(8);
  TKO_MTDC_BHD8->Draw("colz");
  c3->cd(9);
  TKO_MTDC_BHD9->Draw("colz");
  c3->cd(10);
  TKO_MTDC_BHD10->Draw("colz");
  c3->cd(11);
  TKO_MTDC_BHD11->Draw("colz");
  c3->cd(12);
  TKO_MTDC_BHD12->Draw("colz");
  c3->cd(13);
  TKO_MTDC_BHD13->Draw("colz");
  c3->cd(14);
  TKO_MTDC_BHD14->Draw("colz");
  c3->cd(15);
  TKO_MTDC_BHD15->Draw("colz");
  c3->cd(16);
  TKO_MTDC_BHD16->Draw("colz");
  c3->cd(17);
  TKO_MTDC_BHD17->Draw("colz");
  c3->cd(18);
  TKO_MTDC_BHD18->Draw("colz");
  c3->cd(19);
  TKO_MTDC_BHD19->Draw("colz");
  c3->cd(20);
  TKO_MTDC_BHD20->Draw("colz");

  c6 = new TCanvas("c6", "c6", 800, 800);
  c6->Divide(4,8);
  c6->cd(1);
  gPad->SetLogy();
  nhitMTDC0ch0->Draw();
  c6->cd(2);
  gPad->SetLogy();
  nhitMTDC0ch1->Draw();
  c6->cd(3);
  gPad->SetLogy();
  nhitMTDC0ch2->Draw();
  c6->cd(4);
  gPad->SetLogy();
  nhitMTDC0ch3->Draw();
  c6->cd(5);
  gPad->SetLogy();
  nhitMTDC0ch4->Draw();
  c6->cd(6);
  gPad->SetLogy();
  nhitMTDC0ch5->Draw();
  c6->cd(7);
  gPad->SetLogy();
  nhitMTDC0ch6->Draw();
  c6->cd(8);
  gPad->SetLogy();
  nhitMTDC0ch7->Draw();
  c6->cd(9);
  gPad->SetLogy();
  nhitMTDC0ch8->Draw();
  c6->cd(10);
  gPad->SetLogy();
  nhitMTDC0ch9->Draw();
  c6->cd(11);
  gPad->SetLogy();
  nhitMTDC0ch10->Draw();
  c6->cd(12);
  gPad->SetLogy();
  nhitMTDC0ch11->Draw();
  c6->cd(13);
  gPad->SetLogy();
  nhitMTDC0ch12->Draw();
  c6->cd(14);
  gPad->SetLogy();
  nhitMTDC0ch13->Draw();
  c6->cd(15);
  gPad->SetLogy();
  nhitMTDC0ch14->Draw();
  c6->cd(16);
  gPad->SetLogy();
  nhitMTDC0ch15->Draw();
  c6->cd(17);
  gPad->SetLogy();
  nhitMTDC0ch16->Draw();
  c6->cd(18);
  gPad->SetLogy();
  nhitMTDC0ch17->Draw();
  c6->cd(19);
  gPad->SetLogy();
  nhitMTDC0ch18->Draw();
  c6->cd(20);
  gPad->SetLogy();
  nhitMTDC0ch19->Draw();
  c6->cd(21);
  gPad->SetLogy();
  nhitMTDC0ch20->Draw();
  c6->cd(22);
  gPad->SetLogy();
  nhitMTDC0ch21->Draw();
  c6->cd(23);
  gPad->SetLogy();
  nhitMTDC0ch22->Draw();
  c6->cd(24);
  gPad->SetLogy();
  nhitMTDC0ch23->Draw();
  c6->cd(25);
  gPad->SetLogy();
  nhitMTDC0ch24->Draw();
  c6->cd(26);
  gPad->SetLogy();
  nhitMTDC0ch25->Draw();
  c6->cd(27);
  gPad->SetLogy();
  nhitMTDC0ch26->Draw();
  c6->cd(28);
  gPad->SetLogy();
  nhitMTDC0ch27->Draw();
  c6->cd(29);
  gPad->SetLogy();
  nhitMTDC0ch28->Draw();
  c6->cd(30);
  gPad->SetLogy();
  nhitMTDC0ch29->Draw();
  c6->cd(31);
  gPad->SetLogy();
  nhitMTDC0ch30->Draw();
  c6->cd(32);
  gPad->SetLogy();
  nhitMTDC0ch31->Draw();

  c7 = new TCanvas("c7", "c7", 800, 800);
  c7->Divide(4,8);
  c7->cd(1);
  gPad->SetLogy();
  nhitMTDC1ch0->Draw();
  c7->cd(2);
  gPad->SetLogy();
  nhitMTDC1ch1->Draw();
  c7->cd(3);
  gPad->SetLogy();
  nhitMTDC1ch2->Draw();
  c7->cd(4);
  gPad->SetLogy();
  nhitMTDC1ch3->Draw();
  c7->cd(5);
  gPad->SetLogy();
  nhitMTDC1ch4->Draw();
  c7->cd(6);
  gPad->SetLogy();
  nhitMTDC1ch5->Draw();
  c7->cd(7);
  gPad->SetLogy();
  nhitMTDC1ch6->Draw();
  c7->cd(8);
  gPad->SetLogy();
  nhitMTDC1ch7->Draw();
  c7->cd(9);
  gPad->SetLogy();
  nhitMTDC1ch8->Draw();
  c7->cd(10);
  gPad->SetLogy();
  nhitMTDC1ch9->Draw();
  c7->cd(11);
  gPad->SetLogy();
  nhitMTDC1ch10->Draw();
  c7->cd(12);
  gPad->SetLogy();
  nhitMTDC1ch11->Draw();
  c7->cd(13);
  gPad->SetLogy();
  nhitMTDC1ch12->Draw();
  c7->cd(14);
  gPad->SetLogy();
  nhitMTDC1ch13->Draw();
  c7->cd(15);
  gPad->SetLogy();
  nhitMTDC1ch14->Draw();
  c7->cd(16);
  gPad->SetLogy();
  nhitMTDC1ch15->Draw();
  c7->cd(17);
  gPad->SetLogy();
  nhitMTDC1ch16->Draw();
  c7->cd(18);
  gPad->SetLogy();
  nhitMTDC1ch17->Draw();
  c7->cd(19);
  gPad->SetLogy();
  nhitMTDC1ch18->Draw();
  c7->cd(20);
  gPad->SetLogy();
  nhitMTDC1ch19->Draw();
  c7->cd(21);
  gPad->SetLogy();
  nhitMTDC1ch20->Draw();
  c7->cd(22);
  gPad->SetLogy();
  nhitMTDC1ch21->Draw();
  c7->cd(23);
  gPad->SetLogy();
  nhitMTDC1ch22->Draw();
  c7->cd(24);
  gPad->SetLogy();
  nhitMTDC1ch23->Draw();
  c7->cd(25);
  gPad->SetLogy();
  nhitMTDC1ch24->Draw();
  c7->cd(26);
  gPad->SetLogy();
  nhitMTDC1ch25->Draw();
  c7->cd(27);
  gPad->SetLogy();
  nhitMTDC1ch26->Draw();
  c7->cd(28);
  gPad->SetLogy();
  nhitMTDC1ch27->Draw();
  c7->cd(29);
  gPad->SetLogy();
  nhitMTDC1ch28->Draw();
  c7->cd(30);
  gPad->SetLogy();
  nhitMTDC1ch29->Draw();
  c7->cd(31);
  gPad->SetLogy();
  nhitMTDC1ch30->Draw();
  c7->cd(32);
  gPad->SetLogy();
  nhitMTDC1ch31->Draw();

  // c0->Print("plot.pdf(");
  // c2->Print("plot.pdf");
  // c3->Print("plot.pdf");
  // c6->Print("plot.pdf");
  // c7->Print("plot.pdf)");

  c0->Print(hist_file+".pdf(");
  c2->Print(hist_file+".pdf");
  c3->Print(hist_file+".pdf");
  c6->Print(hist_file+".pdf");
  c7->Print(hist_file+".pdf)");

}
