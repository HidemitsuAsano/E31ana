#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h>
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TColor.h>
#include <TProfile.h>
#include <TPaletteAxis.h>
#include <TComplex.h>

#define LOOP 1

int main()
{

  TComplex K_mass   = TComplex(493.677, 0);
  TComplex p_mass   = TComplex(938.272046, 0);
  TComplex pi_mass  = TComplex(139.57, 0);
  TComplex Sp_mass  = TComplex(1189.37, 0);
  TComplex Sm_mass  = TComplex(1197.44, 0);
  TComplex S_mass   = 0.5*(Sp_mass+Sm_mass);

  double Kp_mass     = K_mass.Re()+p_mass.Re();
 
  // definition :: grid of (Re(E), Im(E))
  const double Re_MIN = 1340;
  const double Re_MAX = 1500;
  const int    Re_BIN = Re_MAX-Re_MIN; // 1MeV/BIN
  const double Im_MIN = -50;
  const double Im_MAX = 50;
  const int    Im_BIN = 100; // 0.5MeV/BIN
  
  // PDF
  TH2F *PDF = new TH2F("PDF", "PDF", Re_BIN, Re_MIN, Re_MAX, Im_BIN, Im_MIN, Im_MAX);
  TH2F *K1 = new TH2F("K1", "K1", Re_BIN, Re_MIN, Re_MAX, Im_BIN, Im_MIN, Im_MAX);
  TH2F *K2 = new TH2F("K2", "K2", Re_BIN, Re_MIN, Re_MAX, Im_BIN, Im_MIN, Im_MAX);
  
  
  // Flatte's parameters obtained by the fitting
  TComplex mR = TComplex(1411.6, 0); // (mass_R)
  TComplex g1 = TComplex(0.29, 0); // (g_piS)
  TComplex g2 = TComplex(0.10, 0); // (g_KN)

  TComplex mR_e = TComplex(sqrt(12.4*12.4+3.9*3.9), 0); // (mass_R)
  TComplex g1_e = TComplex(sqrt(0.14*0.14+0.05*0.05), 0); // (g_piS)
  TComplex g2_e = TComplex(sqrt(1.28*1.28+0.10*0.10), 0); // (g_KN)

  // definition :: grid of (mR, g1, g2)
  const double mR_MIN = mR.Re()-mR_e.Re();
  const double mR_MAX = mR.Re()+mR_e.Re();
  const int    mR_BIN = 10;
  const double g1_MIN = g1.Re()-g1_e.Re();
  const double g1_MAX = g1.Re()+g1_e.Re();
  const int    g1_BIN = 10;
  const double g2_MIN = 0; // g2.Re()-g2_e.Re();
  const double g2_MAX = g2.Re()+g2_e.Re();
  const int    g2_BIN = 10;

  TComplex iunit = TComplex(0,1);
  TComplex E, A;
  TComplex k1, k2, G1, G2;

  std::vector<double> pole;
  std::vector<double> width;

  TH2F *POLE = new TH2F("POLE", "POLE", Re_BIN, Re_MIN, Re_MAX, Im_BIN, Im_MIN, Im_MAX);

#if LOOP
  for( int i=0; i<mR_BIN; i++ ){
    for( int j=0; j<g1_BIN; j++ ){
      for( int k=0; k<g2_BIN; k++ ){

	double ii = mR_MIN+(i+0.5)*(mR_MAX-mR_MIN)/mR_BIN;
	double jj = g1_MIN+(j+0.5)*(g1_MAX-g1_MIN)/g1_BIN;
	double kk = g2_MIN+(k+0.5)*(g2_MAX-g2_MIN)/g2_BIN;

	TComplex mR = TComplex(ii, 0); // (mass_R)
	TComplex g1 = TComplex(jj, 0); // (g_piS)
	TComplex g2 = TComplex(kk, 0); // (g_KN)
#endif

	
  for( int x=0; x<Re_BIN; x++ ){
    for( int y=0; y<Im_BIN; y++ ){
      double xx = Re_MIN+(x+0.5)*(Re_MAX-Re_MIN)/Re_BIN;
      double yy = Im_MIN+(y+0.5)*(Im_MAX-Im_MIN)/Im_BIN;

      E = TComplex(xx,yy);

      // decay momentum at the CM flame
      k1 = 0.5/E*TComplex::Sqrt((E*E-(pi_mass+S_mass)*(pi_mass+S_mass))*(E*E-(pi_mass-S_mass)*(pi_mass-S_mass)));
      G1 = g1*k1;
      if( Kp_mass<E.Re() ){ // above the KN threshold
	k2 = 0.5/E*TComplex::Sqrt((E*E-(K_mass+p_mass)*(K_mass+p_mass))*(E*E-(K_mass-p_mass)*(K_mass-p_mass)));
	G2 = g2*k2;
      }else{ // below the KN threshold 
	k2 = 0.5/E*TComplex::Sqrt(-(E*E-(K_mass+p_mass)*(K_mass+p_mass))*(E*E-(K_mass-p_mass)*(K_mass-p_mass)));
	G2 = iunit*g2*k2;
      }
      
      A =  mR*TComplex::Sqrt(G1) / ( mR*mR - E*E - iunit*mR*(G1+G2) );

      PDF->Fill(xx,yy,A.Rho2());
      K1->Fill(xx,yy,k1.Im());
      K2->Fill(xx,yy,k2.Im());
    }
  }
  
  pole.push_back( PDF->ProjectionX()->GetMean() );
  width.push_back( PDF->ProjectionY()->GetMean() );
  POLE->Fill(PDF->ProjectionX()->GetMean(), PDF->ProjectionY()->GetMean());
  //std::cerr<<PDF->ProjectionX()->GetMean()<<" "<<PDF->ProjectionY()->GetMean()<<std::endl;

#if LOOP
      }
    }
  }
#endif

  double center_p = 1415.52;
  double center_w = -18.7373;
  
  std::sort(pole.begin(), pole.end());
  std::vector<double>::iterator itr;
  itr = pole.begin(); // min
  std::cerr<<*itr<<" "<<*itr-center_p<<std::endl;
  itr = pole.end(); --itr; // max
  std::cerr<<*itr<<" "<<*itr-center_p<<std::endl;

  std::sort(width.begin(), width.end());
  itr = width.begin(); // min
  std::cerr<<*itr<<" "<<*itr-center_w<<std::endl;
  itr = width.end(); --itr; // max
  std::cerr<<*itr<<" "<<*itr-center_w<<std::endl;
  
  
  TCanvas *y0;
  y0 = new TCanvas("y0", "", 600, 600);
  y0->Divide(2,2);
  y0->cd(1);
  gPad->SetLogz();
  PDF->Draw("colz");
  y0->cd(2);
  PDF->ProjectionY()->Draw();
  y0->cd(3);
  PDF->ProjectionX()->Draw();
  y0->Print("tmp.pdf");

  TCanvas *y1;
  y1 = new TCanvas("y1", "", 1200, 600);
  y1->Divide(2,1);
  y1->cd(1);
  gPad->SetRightMargin(0.2);
  //gPad->SetLogz();
  K1->Draw("colz");
  y1->cd(2);
  //gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  K2->Draw("colz");
  y1->Print("tmp2.pdf");

  TCanvas *y2;
  y2 = new TCanvas("y2", "", 600, 600);
  y2->Divide(2,2);
  y2->cd(1);
  gPad->SetLogz();
  POLE->Draw("colz");
  y2->cd(2);
  POLE->ProjectionY()->Draw();
  y2->cd(3);
  POLE->ProjectionX()->Draw();
  y2->Print("tmp3.pdf");

  
  //*** for PRC ***//
  TLatex *tex = new TLatex();
  tex->SetTextSize( 0.04 );
  TCanvas *prc = new TCanvas("prc", "", 600, 600);
  gPad->SetLogz();
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gStyle->SetPalette(56);
  TH2F *his = (TH2F*) PDF;
  his->Draw("colz");
  his->GetYaxis()->SetTitleOffset(1.3);
  his->SetXTitle("Re(z) [MeV/c^{2}]");
  his->SetYTitle("Im(z) [MeV/c^{2}]");
  his->SetTitle(0);
  his->SetStats(0);
  his->SetMinimum(0.5);
  his->GetXaxis()->CenterTitle();
  his->GetYaxis()->CenterTitle();
  his->GetXaxis()->SetNdivisions(505);
  his->GetXaxis()->SetTitleSize(0.05);
  his->GetXaxis()->SetLabelSize(0.05);
  his->GetYaxis()->SetNdivisions(507);
  his->GetYaxis()->SetTitleSize(0.05);
  his->GetYaxis()->SetLabelSize(0.05);
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)his->GetListOfFunctions()->FindObject("palette");
  palette->SetY1NDC(0.15);
  palette->SetY2NDC(0.9);
  palette->SetX1NDC(0.9);
  palette->SetX2NDC(0.93);
  TLine *line = new TLine(Kp_mass, Im_MIN, Kp_mass, Im_MAX);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();

  tex->SetTextSize(0.03);
  tex->SetTextAngle(90);
  tex->DrawLatex( Kp_mass-1, Im_MIN+1, "M(K^{-}p)" );

  prc->Print("pole.pdf");
  prc->Print("pole.eps");

  //*** for PRC ***//

  return 0;
}
