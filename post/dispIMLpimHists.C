#include <anacuts.h>
#include <TGeoManager.h>
#include "../src/GlobalVariables.h"

void dispIMLpimHists()
{


  TFile *file = TFile::Open("evanaIMLambdaPim_ppimpim_v15_out.root","READ");
  //TFile *file = TFile::Open("../simpost/simIMLpim_ppimpim_v17_out.root","READ");
  bool SimMode = (std::string(file->GetName()).find("sim")!= std::string::npos);
  
  TCanvas *cBeamMom = new TCanvas("cBeamMom","cBeamMom",1000,800);
  TH1D* BeamMom = (TH1D*)file->Get("BeamMom");
  BeamMom->Draw("HE");
  
  TCanvas *cVtx_XY = new TCanvas("cVtx_XY","cVtx_XY",1000,800);
  TH2D* Vtx_XY_Lcan = (TH2D*)file->Get("Vtx_XY_Lcan");
  Vtx_XY_Lcan->Draw("colz");

  TH2D* Vtx_XY_Lcan_fid = (TH2D*)file->Get("Vtx_XY_Lcan_fid");
  Vtx_XY_Lcan_fid->RebinX(2);
  Vtx_XY_Lcan_fid->RebinY(2);
  Vtx_XY_Lcan_fid->Draw("boxsame");


  TCanvas *cVtx_ZX = new TCanvas("cVtx_ZX","cVtx_ZX",1000,800);
  TH2D* Vtx_ZX_Lcan = (TH2D*)file->Get("Vtx_ZX_Lcan");
  Vtx_ZX_Lcan->Draw("colz");

  TH2D* Vtx_ZX_Lcan_fid = (TH2D*)file->Get("Vtx_ZX_Lcan_fid");
  Vtx_ZX_Lcan_fid->RebinX(2);
  Vtx_ZX_Lcan_fid->RebinY(2);
  Vtx_ZX_Lcan_fid->Draw("boxsame");

  
  TCanvas *cVtx_ZY = new TCanvas("cVtx_ZY","cVtx_ZY",1000,800);
  TH2D* Vtx_ZY_Lcan = (TH2D*)file->Get("Vtx_ZY_Lcan");
  Vtx_ZY_Lcan->Draw("colz");
  TH2D* Vtx_ZY_Lcan_fid = (TH2D*)file->Get("Vtx_ZY_Lcan_fid");
  Vtx_ZY_Lcan_fid->RebinX(2);
  Vtx_ZY_Lcan_fid->RebinY(2);
  Vtx_ZY_Lcan_fid->Draw("boxsame");

  TCanvas *cpcos = new TCanvas("cpcos","cpcos",1000,800);
  TH1D* pcos =  (TH1D*)file->Get("pcos");
  TH1D* pmisscos = (TH1D*)file->Get("pmisscos");
  TH1D* p2cos = (TH1D*)file->Get("p2cos");
  pcos->Draw("HE");
  pmisscos->SetLineColor(2);
  pmisscos->Draw("HEsame");
  p2cos->SetLineColor(3);
  p2cos->Draw("HEsame");

  
  TCanvas *chp2flag = new TCanvas("chp2flag","chp2flag");
  TH1I *hp2flag = (TH1I*)file->Get("hp2flag");
  hp2flag->Draw("HE");

  
  TCanvas *cIMppim1_IMppim2 = new TCanvas("cIMppim1_IMppim2","cIMppim1_IMppim2",1000,800);
  TH2F* IMppim1_IMppim2 = (TH2F*)file->Get("IMppim1_IMppim2");
  IMppim1_IMppim2->Draw("colz");

  TCanvas *cIMppim1_IMp2pim1 = new TCanvas("cIMppim1_IMp2pim1","cIMppim1_IMp2pim1",1000,800);
  TH2F* IMppim1_IMp2pim1 = (TH2F*)file->Get("IMppim1_IMp2pim1");
  IMppim1_IMp2pim1->Draw("colz");

  TCanvas *cIMppim2_IMp2pim2 = new TCanvas("cIMppim2_IMp2pim2","cIMppim2_IMp2pim2",1000,800);
  TH2F* IMppim2_IMp2pim2 = (TH2F*)file->Get("IMppim2_IMp2pim2");
  IMppim2_IMp2pim2->Draw("colz");

  TCanvas *cIMp2pim1_IMp2pim2 = new TCanvas("cIMp2pim1_IMp2pim2","cIMp2pim1_IMp2pim2",1000,800);
  TH2F* IMp2pim1_IMp2pim2 = (TH2F*)file->Get("IMp2pim1_IMp2pim2");
  IMp2pim1_IMp2pim2->Draw("colz");

  TCanvas *cMMom_PMom = new TCanvas("cMMom_PMom","cMMom_PMom",1000,800);
  TH2F* MMom_PMom = (TH2F*)file->Get("MMom_PMom");
  MMom_PMom->Draw("colz");

  TCanvas *cMMom_PMom_2 = new TCanvas("cMMom_PMom_2","cMMom_PMom_2",1000,800);
  TH2F* MMom_PMom_2 = (TH2F*)file->Get("MMom_PMom_2");
  MMom_PMom_2->Draw("colz");

  TCanvas *cq_PMom = new TCanvas("cq_PMom","cq_PMom",1000,800);
  TH2F* q_PMom = (TH2F*)file->Get("q_PMom");
  q_PMom->Draw("colz");

  TCanvas *cPmom = new TCanvas("cPmom","cPmom",1000,800);
  TH1D* PMom = (TH1D*)q_PMom->ProjectionX("PMom");
  PMom->Draw("HE");

  TCanvas *cq_P2Mom = new TCanvas("cq_P2Mom","cq_P2Mom",1000,800);
  TH2F* q_P2Mom = (TH2F*)file->Get("q_P2Mom");
  q_P2Mom->Draw("colz");

  TCanvas *cP2mom = new TCanvas("cP2mom","cP2mom",1000,800);
  TH1D* P2Mom = (TH1D*)q_P2Mom->ProjectionX("P2Mom");
  P2Mom->Draw("HE");

  TCanvas *cMiss = new TCanvas("cMiss","cMiss",1000,800);
  TH2F* MMom_MMass = (TH2F*)file->Get("MMom_MMass");
  TH2F* MMom_MMass_2 = (TH2F*)file->Get("MMom_MMass_2");
  TH1D* MMass = (TH1D*)MMom_MMass->ProjectionX("MMass");
  TH1D* MMass_2 = (TH1D*)MMom_MMass_2->ProjectionX("MMass_2");
  MMass->Draw("HE");
  MMass_2->SetLineColor(2);
  MMass_2->Draw("HEsame");
  double mmax = MMass->GetMaximum();
  TLine *plow = new TLine(anacuts::Proton_MIN,0,anacuts::Proton_MIN,mmax);
  plow->SetLineColor(4);
  plow->SetLineWidth(2.0);
  plow->SetLineStyle(10);
  plow->Draw();
  TLine *phigh = new TLine(anacuts::Proton_MAX,0,anacuts::Proton_MAX,mmax);
  phigh->SetLineColor(4);
  phigh->SetLineWidth(2.0);
  phigh->SetLineStyle(10);
  phigh->Draw();

  TCanvas *cMiss_sum = new TCanvas("cMiss_sum","cMiss_sum",1000,800);
  TH1D* MMass_sum = (TH1D*)MMass->Clone("MMass_sum");
  MMass_sum->Add(MMass_2);
  MMass_sum->Draw("HE");

  TH2F* MMom_MMass_p = (TH2F*)file->Get("MMom_MMass_p");
  TH2F* MMom_MMass_p2 = (TH2F*)file->Get("MMom_MMass_p2");
 
  TH1D* MMass_p = (TH1D*)MMom_MMass_p->ProjectionX("MMass_p");
  TH1D* MMass_p2 = (TH1D*)MMom_MMass_p2->ProjectionX("MMass_p2");

  TH1D* MMass_p_sum = (TH1D*)MMass_p->Clone("MMass_p_sum");
  MMass_p_sum->Add(MMass_p2);
  MMass_p_sum->SetLineColor(2);
  MMass_p_sum->SetFillColor(46);
  MMass_p_sum->Draw("HEsame");
  

  TCanvas *cLambda  = new TCanvas("cLambda","cLambda",1000,800);
  TH1D* IMppim1 = (TH1D*)IMppim1_IMppim2->ProjectionY("IMppim1");
  TH1D* IMppim2 = (TH1D*)IMppim1_IMppim2->ProjectionX("IMppim2");
  TH1D* IMp2pim1 = (TH1D*)IMp2pim1_IMp2pim2->ProjectionY("IMp2pim1");
  TH1D* IMp2pim2 = (TH1D*)IMp2pim1_IMp2pim2->ProjectionX("IMp2pim2");
  IMppim1->GetXaxis()->SetRangeUser(1,1.3);
  IMppim1->Draw("HE");
  IMppim2->SetLineColor(2);
  IMppim2->Draw("HEsame");
  IMp2pim1->SetLineColor(3);
  IMp2pim1->Draw("HEsame");
  IMp2pim2->SetLineColor(4);
  IMp2pim2->Draw("HEsame");
  double imax = IMppim1->GetMaximum();
  TLine *llow = new TLine(anacuts::Lambda_MIN,0,anacuts::Lambda_MIN,imax);
  llow->SetLineColor(5);
  llow->SetLineWidth(2.0);
  llow->SetLineStyle(10);
  llow->Draw();
  TLine *lhigh = new TLine(anacuts::Lambda_MAX,0,anacuts::Lambda_MAX,imax);
  lhigh->SetLineColor(5);
  lhigh->SetLineWidth(2.0);
  lhigh->SetLineStyle(10);
  lhigh->Draw();


  TCanvas *cMMass_wL_or = new TCanvas("cMMass_wL_or","cMMass_wL_or",1000,800);
  TH1D* MMass_wL_or = (TH1D*)file->Get("MMass_wL_or");
  MMass_wL_or->Draw("HE");
  plow->Draw();
  phigh->Draw();
  TLine *phigh_narrow = new TLine(anacuts::Proton_MAX_narrow,0,anacuts::Proton_MAX_narrow,mmax);
  phigh_narrow->SetLineColor(5);
  phigh_narrow->SetLineWidth(2.0);
  phigh_narrow->SetLineStyle(10);
  phigh_narrow->Draw();
  TLine *plow_narrow = new TLine(anacuts::Proton_MIN_narrow,0,anacuts::Proton_MIN_narrow,mmax);
  plow_narrow->SetLineColor(5);
  plow_narrow->SetLineWidth(2.0);
  plow_narrow->SetLineStyle(10);
  plow_narrow->Draw();
  

  TCanvas *cMMass_IMppipi_wL_sum = new TCanvas("cMMass_IMppipi_wL_sum","cMMass_IMppipi_wL_sum",1000,800);
  TH2D* MMass_IMppipi_wL_sum = (TH2D*)file->Get("MMass_IMppipi_wL_sum");
  MMass_IMppipi_wL_sum->Draw("colz");


  TCanvas *cq_IMppipi_p_wL = new TCanvas("cq_IMppipi_p_wL","cq_IMppipi_p_wL",1000,800);
  TH2F* q_IMppipi_p_wL = (TH2F*)file->Get("q_IMppipi_p_wL");
  q_IMppipi_p_wL->Draw("colz");

  TCanvas *cq_IMp2pipi_p2_wL = new TCanvas("cq_IMp2pipi_p2_wL","cq_IMp2pipi_p2_wL",1000,800);
  TH2F* q_IMp2pipi_p2_wL = (TH2F*)file->Get("q_IMp2pipi_p2_wL");
  q_IMp2pipi_p2_wL->Draw("colz");


  TCanvas *cq_IMppipi_p_wL_sum  = new TCanvas("cq_IMppipi_p_wL_sum","cq_IMppipi_p_wL_sum",1000,800);
  TH2F* q_IMppipi_p_wL_sum = (TH2F*)file->Get("q_IMppipi_p_wL_sum");
  q_IMppipi_p_wL_sum->Draw("colz");
  
  TCanvas *cq_IMppipi_p_wL_sum_forward  = new TCanvas("cq_IMppipi_p_wL_sum_f","cq_IMppipi_p_wL_sum_f",1000,800);
  TH2F* q_IMppipi_p_wL_sum_forward = (TH2F*)file->Get("q_IMppipi_p_wL_sum_forward");
  q_IMppipi_p_wL_sum_forward->Draw("colz");
  
  TCanvas *cq_IMppipi_p_wL_sum_fp  = new TCanvas("cq_IMppipi_p_wL_sum_fp","cq_IMppipi_p_wL_sum_fp",1000,800);
  TH2F* q_IMppipi_p_wL_sum_fp = (TH2F*)file->Get("q_IMppipi_p_wL_sum_fp");
  q_IMppipi_p_wL_sum_fp->Draw("colz");

  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname;
  pdfname= "Lpim_data.pdf";
  if(SimMode) pdfname= "Lpim_sim.pdf";
   

  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    c->Modified();
    c->Update();
    std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }
  
  TFile *fbeam = new TFile("fbeam_lpim.root","RECREATE");
  fbeam->cd();
  BeamMom->Write();
  fbeam->Close();

}

