const bool gridon=true;
const bool staton=true;
const bool labelon=true;

#include <iostream>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLegend.h>


#include "../src/IMPiSigmaAnaPar.h"

void QAbeamline(TFile *f);
void QACDS(TFile *f);
void QAForward(TFile *f);
void PhysicsPlots(TFile *f);

void plothists(const char *filename="evanaIMpisigma_all_v23.root")
{
  std::cout << "filename " << filename << std::endl;
  gStyle->SetPadGridX(gridon);
  gStyle->SetPadGridY(gridon);
  gStyle->SetTitleYOffset(1.6);
  if(staton)gStyle->SetOptStat("emruo");
  else      gStyle->SetOptStat(0);
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);
  //gStyle->SetPalette(56);
  TFile *f = new TFile(filename,"READ"); 
  
  std::cout << "infile " << filename <<std::endl;
  TString pdfname = std::string(filename);
  pdfname.Replace(std::string(filename).size()-4,5,"pdf");
  std::cout << "pdfname: " << pdfname << std::endl;
  
  TCanvas *cEventCheck = new TCanvas("cEventCheck","EventCheck");
  cEventCheck->cd();
  TH1F *h1_EventCheck = (TH1F*)f->Get("EventCheck");
  h1_EventCheck->SetXTitle("Event tag"); 
  h1_EventCheck->Draw();

  QAbeamline(f);
   
  //CDS QA 
  QACDS(f);

  //Forward QA
  //QAForward(f);

  if(staton)gStyle->SetOptStat("emruo");
  else      gStyle->SetOptStat(0);
  
  //PhysicsPlots(f);
  
  //centering title of all histograms 
  TIter nexthist(gDirectory->GetList());
  TH1F *h1 = nullptr;
  TH1D *h1d = nullptr;
  TH2F *h2 = nullptr;
  TObject *obj = nullptr;
  while( ((obj = (TObject*)nexthist())!=nullptr) && labelon  ){
    if(obj->InheritsFrom("TH1F")){
      h1 = (TH1F*) obj;
      //h1->SetFillStyle(3004);
      h1->SetFillStyle(3004);
      h1->GetXaxis()->CenterTitle();
    }
    if(obj->InheritsFrom("TH1D")){
      h1d = (TH1D*) obj;
      h1d->SetFillStyle(3004);
      h1d->GetXaxis()->CenterTitle();
    }
    if(obj->InheritsFrom("TH2")){
      h2 = (TH2F*) obj;
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
    }
  }

  TCanvas *c = nullptr;
  //TPDF *pdf = new TPDF(pdfname);
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  // while((c= (TCanvas*)next())){
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    TPaveText *pt = new TPaveText(.82,0.90,0.98,0.99,"NDC");
    pt->AddText("Real Data");
    pt->SetFillColor(kCyan-9);
    pt->SetBorderSize(1);
    pt->Draw();
    c->Modified();
    c->Update();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("Title:%s",c->GetTitle())); 
  }
  //pdf->Close();
  //std::cout << "closing pdf " << std::endl;


  return;
}

void QAForward(TFile *f){
  
  gStyle->SetTitleYOffset(1.6);
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.03);
  
  TCanvas *cNC_overbeta_dE = new TCanvas("NC_overbeta_dE","NC_overbeta_dE");
  cNC_overbeta_dE->cd();
  cNC_overbeta_dE->SetLogz();
  TH2F* h2_NC_overbeta_dE = (TH2F*) f->Get("NC_overbeta_dE");
  h2_NC_overbeta_dE->SetXTitle("1/#beta");
  h2_NC_overbeta_dE->SetYTitle("dE [MeV]");
  h2_NC_overbeta_dE->Draw("colz");
  
  TCanvas *cNC_overbeta = new TCanvas("NC_overbeta","NC_overbeta");
  cNC_overbeta->cd();
  cNC_overbeta->SetLogz();
  TH1D* h1_NC_overbeta = h2_NC_overbeta_dE->ProjectionX("px",9,200);
  h1_NC_overbeta->SetMinimum(50);
  h1_NC_overbeta->GetXaxis()->SetRangeUser(0.5,4.5);
  h1_NC_overbeta->Draw();

  return;
}



void QACDS(TFile *f){

  gStyle->SetTitleYOffset(1.6);
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);
  
  //CDH time distribution
  TCanvas *cCDHtime = new TCanvas("cCDHtime","CDHtime");
  cCDHtime->cd();
  TH2F* h2_CDHtime = (TH2F*)f->Get("CDHtime");
  h2_CDHtime->SetXTitle("CDH seg#");
  h2_CDHtime->SetYTitle("time [nsec]");
  h2_CDHtime->Draw("colz");

  TCanvas *cCDHNeutralSeg = new TCanvas("cCDHNeutralSeg","CDHNeutralSeg");
  cCDHNeutralSeg->cd();
  TH1F* h1_CDHNeutralSeg = (TH1F*)f->Get("CDHNeutralSeg");
  h1_CDHNeutralSeg->SetXTitle("CDH seg#");
  h1_CDHNeutralSeg->Draw();

  TCanvas *cdE_CDHtime_2track = new TCanvas("cdE_CDHtime_2track","cdE_CDHtime_2track");
  cdE_CDHtime_2track->SetLogz();
  TH2F* dE_CDHtime_2track = (TH2F*)f->Get("dE_CDHtime_2track");
  dE_CDHtime_2track->SetXTitle("CDH time [nsec]");
  dE_CDHtime_2track->SetYTitle("dE [MeVee]");
  dE_CDHtime_2track->Draw("colz");
  
  
  TCanvas *cdE_CDHtime_pippimn = new TCanvas("cdE_CDHtime_pippimn","cdE_CDHtime_pippimn");
  cdE_CDHtime_pippimn->SetLogz();
  TH2F* dE_CDHtime_pippimn = (TH2F*)f->Get("dE_CDHtime_pippimn");
  dE_CDHtime_pippimn->SetXTitle("CDH time [nsec]");
  dE_CDHtime_pippimn->SetYTitle("dE [MeVee]");
  dE_CDHtime_pippimn->Draw("colz");
  
  //CDC track chi2 distribution
  TCanvas *ctrackchi2_CDC = new TCanvas("ctrackchi2_CDC","trackchi2_CDC");
  ctrackchi2_CDC->SetLogy();
  TH1F* h1_trackchi2_CDC = (TH1F*)f->Get("trackchi2_CDC");
  h1_trackchi2_CDC->SetXTitle("CDC track chi2/ndf");
  h1_trackchi2_CDC->Draw();
  TH1F* h1_trackchi2_CDC_clone = (TH1F*) h1_trackchi2_CDC->Clone();
  h1_trackchi2_CDC_clone->GetXaxis()->SetRangeUser(0,cdscuts::cds_chi2_max);
  h1_trackchi2_CDC_clone->SetFillColor(kGreen);
  h1_trackchi2_CDC_clone->Draw("same");

  TCanvas *cntrack = new TCanvas("cntrack_CDS","ntrack_CDS");
  TH1F* h1_ntrack_CDS = (TH1F*)f->Get("ntrack_CDS");
  //h1_ntrack_CDS->SetXTitle("# of tracks");
  //h1_ntrack_CDS->Draw();
  TH1F* h1_ntrack_pi_plus = (TH1F*) f->Get("ntrack_pi_plus");
  h1_ntrack_pi_plus->SetTitle("");
  h1_ntrack_pi_plus->SetLineColor(1);
  h1_ntrack_pi_plus->SetXTitle("# of tracks per event");
  h1_ntrack_pi_plus->Draw("");
  TH1F* h1_ntrack_pi_minus = (TH1F*)f->Get("ntrack_pi_minus");
  h1_ntrack_pi_minus->SetLineColor(2);
  h1_ntrack_pi_minus->Draw("same");
  TH1F* h1_ntrack_proton = (TH1F*)f->Get("ntrack_proton");
  h1_ntrack_proton->SetLineColor(3);
  h1_ntrack_proton->Draw("same");
  TH1F* h1_ntrack_K_minus = (TH1F*)f->Get("ntrack_K_minus");
  h1_ntrack_K_minus->SetLineColor(4);
  h1_ntrack_K_minus->Draw("same");
  
  TLegend *legend = new TLegend(0.6,0.7,0.9,0.9);
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  //legend->AddEntry(h1_ntrack_CDS,"# of good CDS track per event","L");
  legend->AddEntry(h1_ntrack_pi_plus,"# of #pi^{+}","L");
  legend->AddEntry(h1_ntrack_pi_minus,"# of #pi^{-}","L");
  legend->AddEntry(h1_ntrack_proton,"# of proton","L");
  legend->AddEntry(h1_ntrack_K_minus,"# of K^{-}","L");
  legend->Draw();


  TCanvas *cCDS_mom_betainv  = new TCanvas("CDS_mom_betainv","CDS_mom_betainv");
  cCDS_mom_betainv->cd();
  cCDS_mom_betainv->SetLogz();
  TH2F *h2_CDS_mom_betainv = (TH2F*)f->Get("PID_CDS_beta");
  //h2_CDS_mom_betainv->GetXaxis()->SetRangeUser(-0.3,1.5);
  h2_CDS_mom_betainv->GetXaxis()->SetTitle("mass^{2} [(GeV/c^{2})^{2}]");
  //h2_CDS_mom_betainv->GetYaxis()->SetRangeUser(-1.0,1.0);
  h2_CDS_mom_betainv->GetYaxis()->SetTitle("1/#beta");
  h2_CDS_mom_betainv->Draw("colz");
  gStyle->SetOptStat(0);
  cCDS_mom_betainv->Update();
  cCDS_mom_betainv->Modified();

   
  TCanvas *cCDS_mass2_mom  = new TCanvas("CDS_mass2_mom","CDS_mass2_mom");
  cCDS_mass2_mom->cd();
  cCDS_mass2_mom->SetLogz();
  TH2F *h2_CDS_mass2_mom = (TH2F*)f->Get("PID_CDS");
  h2_CDS_mass2_mom->GetXaxis()->SetRangeUser(-0.3,1.5);
  h2_CDS_mass2_mom->GetXaxis()->SetTitle("mass^{2} [(GeV/c^{2})^{2}]");
  h2_CDS_mass2_mom->GetYaxis()->SetRangeUser(-1.0,1.0);
  h2_CDS_mass2_mom->GetYaxis()->SetTitle("charge X momentum [GeV/c]");
  h2_CDS_mass2_mom->Draw("colz");
  gStyle->SetOptStat(0);
  cCDS_mass2_mom->Update();
  cCDS_mass2_mom->Modified();
  
  
  TCanvas *cVtx_ZX  = new TCanvas("Vtx_ZX","Vtx_ZX");
  TH2F *h2_Vtx_ZX = (TH2F*)f->Get("Vtx_ZX");
  h2_Vtx_ZX->SetTitle("CDS Vertex");
  h2_Vtx_ZX->GetXaxis()->SetTitle("Z Vertex [cm]");
  h2_Vtx_ZX->GetYaxis()->SetTitle("X Vertex [cm]");
  h2_Vtx_ZX->Draw("colz");
  double maxZX = h2_Vtx_ZX->GetMaximum();
  TH2F *h2_Vtx_ZX_fid = (TH2F*)f->Get("Vtx_ZX_fid");
  h2_Vtx_ZX_fid->RebinX(10);
  h2_Vtx_ZX_fid->RebinY(10);
  h2_Vtx_ZX_fid->SetMaximum(maxZX);
  h2_Vtx_ZX_fid->Draw("box same");

  TCanvas *cVtx_ZY  = new TCanvas("Vtx_ZY","Vtx_ZY");
  TH2F *h2_Vtx_ZY = (TH2F*)f->Get("Vtx_ZY");
  h2_Vtx_ZY->SetTitle("CDS Vertex");
  h2_Vtx_ZY->GetXaxis()->SetTitle("Z Vertex [cm]");
  h2_Vtx_ZY->GetYaxis()->SetTitle("Y Vertex [cm]");
  h2_Vtx_ZY->Draw("colz");
  double maxZY = h2_Vtx_ZY->GetMaximum();
  TH2F *h2_Vtx_ZY_fid = (TH2F*)f->Get("Vtx_ZY_fid");
  h2_Vtx_ZY_fid->SetMaximum(maxZY);
  h2_Vtx_ZY_fid->RebinX(10);
  h2_Vtx_ZY_fid->RebinY(10);
  h2_Vtx_ZY_fid->Draw("box same");
  
  TCanvas *cVtx_XY  = new TCanvas("Vtx_XY","Vtx_XY");
  TH2F *h2_Vtx_XY = (TH2F*)f->Get("Vtx_XY");
  h2_Vtx_XY->SetTitle("CDS Vertex");
  h2_Vtx_XY->GetXaxis()->SetTitle("X Vertex [cm]");
  h2_Vtx_XY->GetYaxis()->SetTitle("Y Vertex [cm]");
  h2_Vtx_XY->Draw("colz");
  double maxXY = h2_Vtx_XY->GetMaximum();
  TH2F *h2_Vtx_XY_fid = (TH2F*) f->Get("Vtx_XY_fid");
  h2_Vtx_XY_fid->SetMaximum(maxXY);
  h2_Vtx_XY_fid->RebinX(10);
  h2_Vtx_XY_fid->RebinY(10);
  h2_Vtx_XY_fid->Draw("box same");
    
  /*
  //fiducial cuts
  TCanvas *cVtx_ZX_fid  = new TCanvas("Vtx_ZX_fid","Vtx_ZX_fid");
  TH2F *h2_Vtx_ZX_fid = (TH2F*)f->Get("Vtx_ZX_fid");
  h2_Vtx_ZX_fid->SetTitle("CDS Vertex");
  h2_Vtx_ZX_fid->GetXaxis()->SetTitle("Z Vertex [cm]");
  h2_Vtx_ZX_fid->GetYaxis()->SetTitle("X Vertex [cm]");
  h2_Vtx_ZX_fid->SetMaximum(maxZX);
  h2_Vtx_ZX_fid->Draw("colz");
  
  TCanvas *cVtx_ZY_fid  = new TCanvas("Vtx_ZY_fid","Vtx_ZY_fid");
  TH2F *h2_Vtx_ZY_fid = (TH2F*)f->Get("Vtx_ZY_fid");
  h2_Vtx_ZY_fid->SetTitle("CDS Vertex");
  h2_Vtx_ZY_fid->GetXaxis()->SetTitle("Z Vertex [cm]");
  h2_Vtx_ZY_fid->GetYaxis()->SetTitle("Y Vertex [cm]");
  h2_Vtx_ZY_fid->SetMaximum(maxZY);
  h2_Vtx_ZY_fid->Draw("colz");
  
  TCanvas *cVtx_XY_fid  = new TCanvas("Vtx_XY_fid","Vtx_XY_fid");
  TH2F *h2_Vtx_XY_fid = (TH2F*) f->Get("Vtx_XY_fid");
  h2_Vtx_XY_fid->SetTitle("CDS Vertex");
  h2_Vtx_XY_fid->GetXaxis()->SetTitle("X Vertex [cm]");
  h2_Vtx_XY_fid->GetYaxis()->SetTitle("Y Vertex [cm]");
  h2_Vtx_XY_fid->SetMaximum(maxXY);
  h2_Vtx_XY_fid->Draw("colz");
  */

  /*
  TCanvas *cVtx_X_center  = new TCanvas("Vtx_X_center","Vtx_X_center");
  TH1F *h1_Vtx_X_center = f->Get("Vtx_X_center");
  gPad->SetLogz(0);
  h1_Vtx_X_center->SetTitle("CDS Vertex");
  h1_Vtx_X_center->GetXaxis()->SetTitle("X Vertex [cm]");
  h1_Vtx_X_center->GetYaxis()->SetTitle("Counts / 25.0 [#mum]");
  h1_Vtx_X_center->Draw("");
  
  TCanvas *cVtx_Z_center  = new TCanvas("Vtx_Z_center","Vtx_Z_center");
  TH1F *h1_Vtx_Z_center = f->Get("Vtx_Z_center");
  h1_Vtx_Z_center->SetTitle("CDS Vertex");
  h1_Vtx_Z_center->GetXaxis()->SetTitle("Z Vertex [cm]");
  h1_Vtx_Z_center->GetYaxis()->SetTitle("Counts / 25.0 [#mum]");
  h1_Vtx_Z_center->Draw("");
  */

  TCanvas *cmul_CDH = new TCanvas("cmul_CDH","mul_CDH");
  TH1F* h1_mul_CDH = (TH1F*)f->Get("mul_CDH");
  h1_mul_CDH->SetTitle("CDH multiplicity");
  h1_mul_CDH->SetXTitle("# of CDH segment");
  h1_mul_CDH->Draw();
  TH1F* h1_mul_CDH_clone = (TH1F*)h1_mul_CDH->Clone();
  h1_mul_CDH_clone->GetXaxis()->SetRangeUser(2.5,3.5);
  h1_mul_CDH_clone->SetFillColor(2);
  h1_mul_CDH_clone->Draw("same");

  TCanvas *cmul_CDH_assoc = new TCanvas("cmul_CDH_assoc","mul_CDH_assoc");
  TH1F* h1_mul_CDH_assoc = (TH1F*)f->Get("mul_CDH_assoc");
  h1_mul_CDH_assoc->SetTitle("# CDH seg. associated with CDC track");
  h1_mul_CDH_assoc->SetXTitle("# of CDH segment");
  h1_mul_CDH_assoc->Draw();

  TCanvas *cdiff_CDH = new TCanvas("cdiff_CDH","diff_CDH");
  TH1F* h1_diff_CDH = (TH1F*)f->Get("diff_CDH");
  h1_diff_CDH->SetTitle("Not used CDH segment # - CDH seg# used for charged");
  h1_diff_CDH->SetXTitle("diff of CDH seg#");
  h1_diff_CDH->Draw();
  
  TCanvas *cdiff_CDH_CDC = new TCanvas("cdiff_CDH_CDC","diff_CDH_CDC");
  TH1F* h1_cdiff_CDH_CDC = (TH1F*)f->Get("diff_CDH_CDC");
  h1_cdiff_CDH_CDC->SetXTitle("angle [deg]");
  h1_cdiff_CDH_CDC->Draw();
  
  
  TCanvas *cdE_betainv_fiducial = new TCanvas("cdE_betainv_fiducial","dE_betainv_fid");
  TH2F* h2_dE_betainv = (TH2F*)f->Get("dE_betainv_fid");
  h2_dE_betainv->SetXTitle("1/#beta");
  h2_dE_betainv->SetYTitle("dE [MeV]");
  h2_dE_betainv->Draw("colz");

  TCanvas *cdE_betainv_fiducial_px = new TCanvas("cdE_betainv_fiducial_px","dE_betainv_fiducial_px");
  int bin2mev = h2_dE_betainv->GetYaxis()->FindBin(2.0);
  int bin4mev = h2_dE_betainv->GetYaxis()->FindBin(4.0);
  int bin6mev = h2_dE_betainv->GetYaxis()->FindBin(6.0);
  
  TH1D* h1_nocut = h2_dE_betainv->ProjectionX("px");
  TH1D* h1_2mevcut = h2_dE_betainv->ProjectionX("px1",bin2mev,-1);
  TH1D* h1_4mevcut = h2_dE_betainv->ProjectionX("px2",bin4mev,-1);
  TH1D* h1_6mevcut = h2_dE_betainv->ProjectionX("px3",bin6mev,-1);
  
  h1_nocut->Draw();
  h1_2mevcut->SetLineColor(2);
  h1_2mevcut->Draw("same");
  h1_4mevcut->SetLineColor(3);
  h1_4mevcut->Draw("same");
  h1_6mevcut->SetLineColor(4);
  h1_6mevcut->Draw("same");
  
  TCanvas *c_dE_fiducial = new TCanvas("c_dE_fiducial","dE_betainv_fid_beta");
  TH2F* h2_dE_betainv_fid_beta = (TH2F*)f->Get("dE_betainv_fid_beta");
  TH1D* h1_dE_beta = h2_dE_betainv_fid_beta->ProjectionY("py");
  h1_dE_beta->Draw();

  TCanvas *cDCA_pippim = new TCanvas("cDCA_pippim","DCA_pippim");
  TH1F* h1_DCA_pippim = new (TH2F*)f->Get("DCA_pippim");
  h1_DCA_pipim->Draw();




  return;
}

void PhysicsPlots(TFile *f){

  //physics plots
  //IM pi+pi- when nCDH==3, fudicial OK,
  TCanvas *c_IMpipi = new TCanvas("c_IMpipi","c_IMpipi");
  TH1F* h1_IMpipi = (TH1F*)f->Get("IMpipi_dE");
  c_IMpipi->cd();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.03);
  h1_IMpipi->SetTitle("Invariant Mass #pi^{+}#pi^{-}");
  h1_IMpipi->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  h1_IMpipi->SetYTitle("Counts/ 1 MeV");
  h1_IMpipi->GetYaxis()->SetTitleOffset(1.4);
  h1_IMpipi->Draw();
  TH1F* h1_IMpipi_clone = (TH1F*) h1_IMpipi->Clone();
  h1_IMpipi_clone->SetFillColor(2);
  h1_IMpipi_clone->GetXaxis()->SetRangeUser(anacuts::pipi_MIN,anacuts::pipi_MAX);
  h1_IMpipi_clone->Draw("same");
  

  TCanvas *cMMom_MMass_fid_beta_dE = new TCanvas("cMMom_MMass_fid_beta_dE","cMMom_MMass_fid_beta_dE");
  TH2F* h2_MMom_MMass_fid_beta_dE = (TH2F*) f->Get("MMom_MMass_fid_beta_dE");
  h2_MMom_MMass_fid_beta_dE->SetXTitle("missing mass [GeV/c^{2}]");
  h2_MMom_MMass_fid_beta_dE->SetYTitle("missing momentum [GeV/c]");
  h2_MMom_MMass_fid_beta_dE->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_dE_px = new TCanvas("cMMom_MMass_fid_beta_dE_px","cMMom_MMass_fid_beta_dE_px");
  TH1D* proMMass = h2_MMom_MMass_fid_beta_dE->ProjectionX();
  proMMass->Draw();
  
  TCanvas *cMMom_MMass_fid_beta = new TCanvas("cMMom_MMass_fid_beta","cMMom_MMass_fid_beta");
  TH2F* h2_MMom_MMass_fid_beta = (TH2F*)f->Get("MMom_MMass_fid_beta_dE");
  h2_MMom_MMass_fid_beta->SetXTitle("missing mass [GeV/c^{2}]");
  h2_MMom_MMass_fid_beta->SetYTitle("missing momentum [GeV/c]");
  h2_MMom_MMass_fid_beta->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_px = new TCanvas("cMMom_MMass_fid_beta_px","cMMom_MMass_fid_beta_px");
  TH1D* proMMass2 = h2_MMom_MMass_fid_beta->ProjectionX();
  proMMass2->Draw();

 
  TCanvas *cMMom_MMass_fid_beta_wSid = new TCanvas("cMMom_MMass_fid_beta_wSid","cMMom_MMass_fid_beta_wSid");
  TH2F* h2_MMom_MMass_fid_beta_wSid = (TH2F*) f->Get("MMom_MMass_fid_beta_dE_woK0_wSid");
  h2_MMom_MMass_fid_beta_wSid->SetXTitle("missing mass [GeV/c^{2}]");
  h2_MMom_MMass_fid_beta_wSid->SetYTitle("missing momentum [GeV/c]");
  h2_MMom_MMass_fid_beta_wSid->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_wSid_px = new TCanvas("cMMom_MMass_fid_beta_wSid_px","cMMom_MMass_fid_beta_wSid_px");
  TH1D* proMMass2_wSid = h2_MMom_MMass_fid_beta_wSid->ProjectionX();
  proMMass2->Draw();
  proMMass2_wSid->SetLineColor(2);
  proMMass2_wSid->Draw("same");

  
  // Sigma+/- selection w/o missing N ID after removing K0
  TCanvas *cIMnpim_IMnpip = new TCanvas("cIMnpim_IMnpip_dE_woK0","cIMnpim_IMnpip_dE_woK0");
  TH2F* h2_IMnpim_IMnpip_dE_woK0 = (TH2F*)f->Get("IMnpim_IMnpip_dE_woK0");
  h2_IMnpim_IMnpip_dE_woK0->SetXTitle("IM(#pi^{+}n) [GeV/c^{2}]");
  h2_IMnpim_IMnpip_dE_woK0->SetYTitle("IM(#pi^{-}n) [GeV/c^{2}]");
  h2_IMnpim_IMnpip_dE_woK0->Draw("colz");
  

  TCanvas *cIMnpim_IMnpip_px = new TCanvas("cIMnpim_IMnpip_dE_woK0_px","cIMnpim_IMnpip_dE_woK0_px");
  TH1D* h2_IMnpim_IMnpip_dE_woK0_px = h2_IMnpim_IMnpip_dE_woK0->ProjectionX();
  h2_IMnpim_IMnpip_dE_woK0_px->Draw();
  TH1D* h2_IMnpim_IMnpip_dE_woK0_px_clone = (TH1D*)h2_IMnpim_IMnpip_dE_woK0_px->Clone();
  h2_IMnpim_IMnpip_dE_woK0_px_clone->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  h2_IMnpim_IMnpip_dE_woK0_px_clone->SetFillColor(3);
  h2_IMnpim_IMnpip_dE_woK0_px_clone->Draw("same");

  TCanvas *cIMnpim_IMnpip_py = new TCanvas("cIMnpim_IMnpip_dE_woK0_py","cIMnpim_IMnpip_dE_woK0_py");
  TH1D* h2_IMnpim_IMnpip_dE_woK0_py = h2_IMnpim_IMnpip_dE_woK0->ProjectionY();
  h2_IMnpim_IMnpip_dE_woK0_py->Draw();
  TH1D* h2_IMnpim_IMnpip_dE_woK0_py_clone = (TH1D*)h2_IMnpim_IMnpip_dE_woK0_py->Clone();
  h2_IMnpim_IMnpip_dE_woK0_py_clone->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  h2_IMnpim_IMnpip_dE_woK0_py_clone->SetFillColor(4);
  h2_IMnpim_IMnpip_dE_woK0_py_clone->Draw("same");
  
  /*
  // Sigma+/- selection w/ missing N ID after removing K0
  TCanvas *cIMnpim_IMnpip_wmn = new TCanvas("cIMnpim_IMnpip_wmn","cIMnpim_IMnpip_wmn");
  TH2F* h2_IMnpim_IMnpip_wmn = f->Get("IMnpim_IMnpip_wmn_dE");
  h2_IMnpim_IMnpip_wmn->SetXTitle("IM(#pi^{+}n) [GeV/c^{2}]");
  h2_IMnpim_IMnpip_wmn->SetYTitle("IM(#pi^{-}n) [GeV/c^{2}]");
  h2_IMnpim_IMnpip_wmn->Draw("colz");
  

  TCanvas *cIMnpim_IMnpip_wmn_px = new TCanvas("cIMnpim_IMnpip_wmn_px","cIMnpim_IMnpip_wmn_px");
  TH1D* h2_IMnpim_IMnpip_wmn_px = h2_IMnpim_IMnpip_wmn->ProjectionX();
  h2_IMnpim_IMnpip_wmn_px->Draw();
  TH1D* h2_IMnpim_IMnpip_wmn_px_clone = h2_IMnpim_IMnpip_wmn_px->Clone();
  h2_IMnpim_IMnpip_wmn_px_clone->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  h2_IMnpim_IMnpip_wmn_px_clone->SetFillColor(3);
  h2_IMnpim_IMnpip_wmn_px_clone->Draw("same");

  TCanvas *cIMnpim_IMnpip_wmn_py = new TCanvas("cIMnpim_IMnpip_wmn_py","cIMnpim_IMnpip_wmn_py");
  TH1D* h2_IMnpim_IMnpip_wmn_py = h2_IMnpim_IMnpip_wmn->ProjectionY();
  h2_IMnpim_IMnpip_wmn_py->Draw();
  TH1D* h2_IMnpim_IMnpip_wmn_py_clone = h2_IMnpim_IMnpip_wmn_py->Clone();
  h2_IMnpim_IMnpip_wmn_py_clone->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  h2_IMnpim_IMnpip_wmn_py_clone->SetFillColor(4);
  h2_IMnpim_IMnpip_wmn_py_clone->Draw("same");
  

  
  TCanvas *cIMmnpim_IMmnpip = new TCanvas("cIMmnpim_IMmnpip","cIMmnpim_IMmnpip");
  TH2F* h2_IMmnpim_IMmnpip = f->Get("IMmnpim_IMmnpip");
  h2_IMmnpim_IMmnpip->SetXTitle("IM(#pi^{+}n) [GeV/c^{2}]");
  h2_IMmnpim_IMmnpip->SetYTitle("IM(#pi^{-}n) [GeV/c^{2}]");
  h2_IMmnpim_IMmnpip->Draw("colz");
  
  TCanvas *cIMmnpim_IMmnpip_px = new TCanvas("cIMmnpim_IMmnpip_px","cIMmnpim_IMmnpip_px");
  TH1D* h2_IMmnpim_IMmnpip_px = h2_IMmnpim_IMmnpip->ProjectionX();
  h2_IMmnpim_IMmnpip_px->Draw();
  //TH1D* h2_IMmnpim_IMmnpip_px_clone = h2_IMmnpim_IMmnpip_px->Clone();
  //h2_IMmnpim_IMmnpip_px_clone->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  //h2_IMmnpim_IMmnpip_px_clone->SetFillColor(3);
  //h2_IMmnpim_IMmnpip_px_clone->Draw("same");

  TCanvas *cIMmnpim_IMmnpip_py = new TCanvas("cIMmnpim_IMmnpip_py","cIMmnpim_IMmnpip_py");
  TH1D* h2_IMmnpim_IMmnpip_py = h2_IMmnpim_IMmnpip->ProjectionY();
  h2_IMmnpim_IMmnpip_py->Draw();
  //TH1D* h2_IMmnpim_IMmnpip_py_clone = h2_IMmnpim_IMmnpip_py->Clone();
  //h2_IMmnpim_IMmnpip_py_clone->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  //h2_IMmnpim_IMmnpip_py_clone->SetFillColor(4);
  //h2_IMmnpim_IMmnpip_py_clone->Draw("same");

  TCanvas *cIMnpipi_wSid = new TCanvas("cIMnpipi_wSid","cIMnpipi_wSid");
  //TH1F* h1_IMnpipi_wSid = f->Get("IMnpipi_wSid");
  //h1_IMnpipi_wSid->SetXTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  //h1_IMnpipi_wSid->SetYTitle("Counts per 10 MeV/c^{2}");
  //h1_IMnpipi_wSid->Draw("HE");
  //TH1F* h1_IMnpipi_woK0_wSid = f->Get("IMnpipi_woK0_wSid");
  //h1_IMnpipi_woK0_wSid->SetLineColor(2);
  //h1_IMnpipi_woK0_wSid->Draw("HEsame");
  
  TH1F* h1_IMnpipi_woK0_wSid = f->Get("IMnpipi_woK0_wSid");
  h1_IMnpipi_woK0_wSid->SetXTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  h1_IMnpipi_woK0_wSid->SetYTitle("Counts per 10 MeV/c^{2}");
  h1_IMnpipi_woK0_wSid->Draw("HE");
  

  //dE dep. check
  TCanvas *cdE_MMom_fid_beta_woK0 = new TCanvas("cdE_MMom_fid_beta_woK0","cdE_MMom_fid_beta_woK0");
  TH2F* h2_dE_MMom_fid_beta_woK0 = f->Get("dE_MMom_fid_beta_woK0");
  h2_dE_MMom_fid_beta_woK0->SetXTitle("Missing Momentum [GeV/c]");
  h2_dE_MMom_fid_beta_woK0->SetYTitle("Energy Dep. [MeVee]");
  h2_dE_MMom_fid_beta_woK0->Draw("colz");

  TCanvas *cdE_MMass_fid_beta_woK0 = new TCanvas("cdE_MMass_fid_beta_woK0","cdE_MMass_fid_beta_woK0");
  TH2F* h2_dE_MMass_fid_beta_woK0 = f->Get("dE_MMass_fid_beta_woK0");
  h2_dE_MMass_fid_beta_woK0->SetXTitle("Missing Mass [GeV/c^{2}]");
  h2_dE_MMass_fid_beta_woK0->SetYTitle("Energy Dep. [MeVee]");
  h2_dE_MMass_fid_beta_woK0->Draw("colz");

  //
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid","cMMnmiss_IMnpipi_woK0_wSid");
  TH2F* h2_MMnmiss_IMnpipi_woK0_wSid = f->Get("MMnmiss_IMnpipi_woK0_wSid");
  h2_MMnmiss_IMnpipi_woK0_wSid->SetXTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  h2_MMnmiss_IMnpipi_woK0_wSid->SetYTitle("Missing Momentum [GeV/c]");
  h2_MMnmiss_IMnpipi_woK0_wSid->Draw("colz");
  */ 
  TCanvas *cCosn_IMnpipi = new TCanvas("cCosn_IMnpipi","cCosn_IMnpipi");
  TH2F* h2_cCosn_IMnpipi = (TH2F*) f->Get("Cosn_IMnpipi_woK0_wSid_n");
  h2_cCosn_IMnpipi->SetXTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  h2_cCosn_IMnpipi->SetYTitle("cos#theta_{CM}");
  h2_cCosn_IMnpipi->Draw("colz");
  
  TCanvas *cq_IMnpipi_woK0_wSid = new TCanvas("cq_IMnpipi_woK0_wSid","cq_IMnpipi_woK0_wSid");
  TH2F* h2_q_IMnpipi_woK0_wSid_n = (TH2F*)f->Get("q_IMnpipi_woK0_wSid_n");
  h2_q_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  h2_q_IMnpipi_woK0_wSid_n->SetYTitle("mom. transfer [GeV/c]");
  h2_q_IMnpipi_woK0_wSid_n->RebinX(2);
  h2_q_IMnpipi_woK0_wSid_n->RebinY(6);
  h2_q_IMnpipi_woK0_wSid_n->GetXaxis()->SetRangeUser(1.2,2.0);
  //h2_q_IMnpipi_woK0_wSid_n->SetTitle("");
  h2_q_IMnpipi_woK0_wSid_n->Draw("colz");
  
  //find Kp threshold
  int kpbin = h2_q_IMnpipi_woK0_wSid_n->GetXaxis()->FindBin(1.432);
  //find 350 MeV/c bin
  int qcut350 = h2_q_IMnpipi_woK0_wSid_n->GetYaxis()->FindBin(0.350);  

  TCanvas *cq_IMnpipi_woK0_wSid_px = new TCanvas("cq_IMnpipi_woK0_wSid_px","");   
  TH1D* h2_q_IMnpipi_woK0_wSid_n_px = h2_q_IMnpipi_woK0_wSid_n->ProjectionX("px");
  h2_q_IMnpipi_woK0_wSid_n_px->GetXaxis()->SetRangeUser(1.2,2.0);
  //h2_q_IMnpipi_woK0_wSid_n_px->SetTitle("");
  //h2_q_IMnpipi_woK0_wSid_n_px->GetXaxis()->SetRange(0,kpbin);
  //h2_q_IMnpipi_woK0_wSid_n_px->SetMaximum(1700);
  h2_q_IMnpipi_woK0_wSid_n_px->Draw("E");

  TH1D* h2_q_IMnpipi_woK0_wSid_n_px_qcut350 = h2_q_IMnpipi_woK0_wSid_n->ProjectionX("pxcut",qcut350,300);
  h2_q_IMnpipi_woK0_wSid_n_px_qcut350->SetLineColor(2);
  //h2_q_IMnpipi_woK0_wSid_n_px_qcut350->Draw("Esame");
  //h2_q_IMnpipi_woK0_wSid_n_px->SetTitle("");



  std::cout << "kpbin:"<< kpbin << std::endl;
  std::cout << "qcut350:" << qcut350 << std::endl;

  TCanvas *cq_IMnpipi_woK0_wSid_py_nocut = new TCanvas("cq_IMnpipi_woK0_wSid_py_nocut","cq_IMnpipi_woK0_wSid_py_nocut");   
  TH1D* h2_q_IMnpipi_woK0_wSid_n_py_nocut = h2_q_IMnpipi_woK0_wSid_n->ProjectionY("pycut");
  //h2_q_IMnpipi_woK0_wSid_n_py_nocut->SetTitle("");
  //h2_q_IMnpipi_woK0_wSid_n_py_nocut->Rebin(4);
  h2_q_IMnpipi_woK0_wSid_n_py_nocut->Draw("E");
  
  
  //TCanvas *cq_IMnpipi_woK0_wSid_py = new TCanvas("cq_IMnpipi_woK0_wSid_py","cq_IMnpipi_woK0_wSid_py");   
  TH1D* h2_q_IMnpipi_woK0_wSid_n_py = h2_q_IMnpipi_woK0_wSid_n->ProjectionY("py",0,kpbin-1);
  h2_q_IMnpipi_woK0_wSid_n_py->SetTitle("");
  //h2_q_IMnpipi_woK0_wSid_n_py->Rebin(4);
  h2_q_IMnpipi_woK0_wSid_n_py->SetLineColor(2);
  h2_q_IMnpipi_woK0_wSid_n_py->Draw("Esame");

  //TCanvas *cq_IMnpipi_woK0_wSid_py2 = new TCanvas("cq_IMnpipi_woK0_wSid_py2","cq_IMnpipi_woK0_wSid_py2");   
  TH1D* h2_q_IMnpipi_woK0_wSid_n_py2 = h2_q_IMnpipi_woK0_wSid_n->ProjectionY("py2",kpbin,200);
  h2_q_IMnpipi_woK0_wSid_n_py2->SetTitle("");
  //h2_q_IMnpipi_woK0_wSid_n_py2->Rebin(4);
  h2_q_IMnpipi_woK0_wSid_n_py2->SetLineColor(3);
  h2_q_IMnpipi_woK0_wSid_n_py2->Draw("Esame");

}


void QAbeamline(TFile *f){
  
  gStyle->SetTitleYOffset(1.6);
  gStyle->SetOptStat("emruo");
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);

  //QA plots
  TH1F *h1_nBHD = (TH1F*) f->Get("mul_BHD");
  TCanvas *cnBHD =  new TCanvas("nBHD","mul_BHD");
  cnBHD->cd();
  h1_nBHD->GetXaxis()->SetTitle("BHD multiplicity");
  h1_nBHD->Draw();
  TH1F *h1_nT0 = (TH1F*) f->Get("mul_T0");
  TCanvas *cnT0 =  new TCanvas("nT0","mul_T0");
  cnT0->cd();
  h1_nT0->GetXaxis()->SetTitle("T0 multiplicity");
  h1_nT0->Draw();

  TH1F *h1_T0BHDtof = (TH1F*)f->Get("tof_T0BHD");
  TCanvas *ctof = new TCanvas("tof","tof_T0BHD");
  ctof->cd();
  h1_T0BHDtof->GetXaxis()->SetRangeUser(24,38);
  h1_T0BHDtof->GetXaxis()->SetTitle("T0-BHD tof [nsec.]");
  h1_T0BHDtof->Draw();

  TH1F *h1_T0BHDtof_cut = (TH1F*) h1_T0BHDtof->Clone();
  h1_T0BHDtof_cut->GetXaxis()->SetRangeUser(blcuts::beam_tof_k_min,blcuts::beam_tof_k_max);
  h1_T0BHDtof_cut->SetFillColor(kGreen);
  h1_T0BHDtof_cut->Draw("same");

  //BLC1
  TH1F *h1_nBLC1 = (TH1F*) f->Get("ntrack_BLC1");
  TCanvas *cnBLC1 = new TCanvas("nBLC1","ntrack_BLC1");
  cnBLC1->cd();
  h1_nBLC1->GetXaxis()->SetTitle("BLC1 # of tracks");
  h1_nBLC1->Draw();
  
  TH1F *h1_BLC1time = (TH1F*) f->Get("tracktime_BLC1");
  TCanvas *cBLC1time = new TCanvas("BLC1time","tracktime_BLC1");
  cBLC1time->cd();
  cBLC1time->SetLogy();
  h1_BLC1time->GetXaxis()->SetRangeUser(-50,50);
  h1_BLC1time->GetXaxis()->SetTitle("BLC1 track time [nsec]");
  h1_BLC1time->Draw();
  TH1F *h1_BLC1time_cut = (TH1F*)h1_BLC1time->Clone();
  h1_BLC1time_cut->GetXaxis()->SetRangeUser(blcuts::blc1_time_min,blcuts::blc1_time_max);
  h1_BLC1time_cut->SetFillColor(kGreen);
  h1_BLC1time_cut->Draw("same");

  TH1F *h1_BLC1chi2 = (TH1F*)f->Get("trackchi2_BLC1");
  TCanvas *cBLC1chi2 = new TCanvas("BLC1chi2","trackchi2_BLC1");
  cBLC1chi2->cd();
  cBLC1chi2->SetLogy();
  h1_BLC1chi2->GetXaxis()->SetTitle("BLC1 chi2/ndf");
  h1_BLC1chi2->Draw();
  TH1F *h1_BLC1chi2_cut = (TH1F*)h1_BLC1chi2->Clone();
  h1_BLC1chi2_cut->GetXaxis()->SetRangeUser(0,blcuts::blc1_chi2_max);
  h1_BLC1chi2_cut->SetFillColor(kGreen);
  h1_BLC1chi2_cut->Draw("same");


  //BLC2
  TH1F *h1_nBLC2 = (TH1F*)f->Get("ntrack_BLC2");
  TCanvas *cnBLC2 = new TCanvas("nBLC2","ntrack_BLC2");
  cnBLC2->cd();
  h1_nBLC2->GetXaxis()->SetTitle("BLC2 # of tracks");
  h1_nBLC2->Draw();

  TH1F *h1_BLC2time = (TH1F*)f->Get("tracktime_BLC2");
  TCanvas *cBLC2time = new TCanvas("BLC2time","tracktime_BLC2");
  cBLC2time->cd();
  cBLC2time->SetLogy();
  h1_BLC2time->GetXaxis()->SetTitle("BLC2 track time [nsec.]");
  h1_BLC2time->GetXaxis()->SetRangeUser(-30,30);
  h1_BLC2time->Draw();
  TH1F *h1_BLC2time_cut = (TH1F*)h1_BLC2time->Clone();
  h1_BLC2time_cut->GetXaxis()->SetRangeUser(blcuts::blc2_time_min,blcuts::blc2_time_max);
  h1_BLC2time_cut->SetFillColor(kGreen);
  h1_BLC2time_cut->Draw("same");

  TH1F *h1_BLC2chi2 = (TH1F*)f->Get("trackchi2_BLC2");
  TCanvas *cBLC2chi2 = new TCanvas("BLC2chi2","trackchi2_BLC2");
  cBLC2chi2->cd();
  cBLC2chi2->SetLogy();
  h1_BLC2chi2->GetXaxis()->SetTitle("BLC2 chi2/ndf");
  h1_BLC2chi2->Draw();
  TH1F* h1_BLC2chi2_cut = (TH1F*)h1_BLC2chi2->Clone();
  h1_BLC2chi2_cut->GetXaxis()->SetRangeUser(0,blcuts::blc2_chi2_max);
  h1_BLC2chi2_cut->SetFillColor(kGreen);
  h1_BLC2chi2_cut->Draw("same");
  


  //BPC
  TH1F *h1_nBPC = (TH1F*)f->Get("ntrack_BPC");
  TCanvas *cnBPC = new TCanvas("nBPC","ntrack_BPC");
  cnBPC->cd();
  h1_nBPC->GetXaxis()->SetTitle("BPC # of tracks");
  h1_nBPC->Draw();

  TH1F *h1_BPCtime = (TH1F*)f->Get("tracktime_BPC");
  TCanvas *cBPCtime = new TCanvas("BPCtime","tracktime_BPC");
  cBPCtime->cd();
  cBPCtime->SetLogy();
  h1_BPCtime->GetXaxis()->SetTitle("BPC track time [nsec.]");
  h1_BPCtime->GetXaxis()->SetRangeUser(-30,30);
  h1_BPCtime->Draw();
  TH1F *h1_BPCtime_cut = (TH1F*)h1_BPCtime->Clone();
  h1_BPCtime_cut->GetXaxis()->SetRangeUser(blcuts::bpc_time_min,blcuts::bpc_time_max);
  h1_BPCtime_cut->SetFillColor(kGreen);
  h1_BPCtime_cut->Draw("same");

  TH1F *h1_BPCchi2 = (TH1F*)f->Get("trackchi2_BPC");
  TCanvas *cBPCchi2 = new TCanvas("BPCchi2","trackchi2_BPC");
  cBPCchi2->cd();
  cBPCchi2->SetLogy();
  h1_BPCchi2->GetXaxis()->SetTitle("BPC chi2/ndf");
  h1_BPCchi2->Draw();
  TH1F *h1_BPCchi2_cut = (TH1F*)h1_BPCchi2->Clone();
  h1_BPCchi2_cut->GetXaxis()->SetRangeUser(0,blcuts::bpc_chi2_max);
  h1_BPCchi2_cut->SetFillColor(kGreen);
  h1_BPCchi2_cut->Draw("same");


  //D5
  TH1F *h1_D5chi2 = (TH1F*)f->Get("trackchi2_beam");
  TCanvas *cD5chi2 = new TCanvas("D5chi2","trackchi2_beam");
  cD5chi2->cd();
  cD5chi2->SetLogy();
  h1_D5chi2->GetXaxis()->SetTitle("D5 chi2/ndf");
  h1_D5chi2->Draw();
  TH1F *h1_D5chi2_cut = (TH1F*)h1_D5chi2->Clone();
  h1_D5chi2_cut->GetXaxis()->SetRangeUser(0,blcuts::d5_chi2_max);
  h1_D5chi2_cut->SetFillColor(kGreen);
  h1_D5chi2_cut->Draw("same");

  TH1F *h1_D5mom = (TH1F*)f->Get("momentum_beam");
  TCanvas *cDmom = new TCanvas("D5mom","momentum_beam");
  cDmom->cd();
  h1_D5mom->GetXaxis()->SetTitle("beam mom. [GeV/c]");
  h1_D5mom->Draw();
  
  //BLC2-BPC matching
  TH2F* h2_BLC2_BPC_posdiff = (TH2F*)f->Get("dydx_BLC2BPC");
  TCanvas *cBLC2_BPC_posdiff = new TCanvas("cdydx_BLC2BPC","dydx_BLC2BPC");
  h2_BLC2_BPC_posdiff->SetXTitle("BLC2BPC track dx [cm]");
  h2_BLC2_BPC_posdiff->SetYTitle("BLC2BPC track dy [cm]");
  //gPad->SetLeftMargin(0.13);
  //gPad->SetRightMargin(0.03);
  cBLC2_BPC_posdiff->SetLeftMargin(0.13);   
  cBLC2_BPC_posdiff->SetRightMargin(0.03);  
  cBLC2_BPC_posdiff->cd();
  h2_BLC2_BPC_posdiff->Draw("colz");

  TCanvas *cBLC2_BPC_posdiff_px = new TCanvas("cBLC2_BPC_posdiff_px","BLC2_BPC_posdiff_px");
  TH1D* h2_BLC2_BPC_posdiff_px =  h2_BLC2_BPC_posdiff->ProjectionX();
  h2_BLC2_BPC_posdiff_px->Draw();
  TH1D* h2_BLC2_BPC_posdiff_px_cut = (TH1D*) h2_BLC2_BPC_posdiff_px->Clone();
  h2_BLC2_BPC_posdiff_px_cut->GetXaxis()->SetRangeUser(blcuts::blc2bpc_x_min,blcuts::blc2bpc_x_max);
  h2_BLC2_BPC_posdiff_px_cut->SetFillColor(kGreen);
  h2_BLC2_BPC_posdiff_px_cut->Draw("same");

  TCanvas *cBLC2_BPC_posdiff_py = new TCanvas("cBLC2_BPC_posdiff_py","BLC2_BPC_posdiff_py");
  TH1D* h2_BLC2_BPC_posdiff_py =  h2_BLC2_BPC_posdiff->ProjectionY();
  h2_BLC2_BPC_posdiff_py->Draw();
  TH1D* h2_BLC2_BPC_posdiff_py_cut = (TH1D*) h2_BLC2_BPC_posdiff_py->Clone();
  h2_BLC2_BPC_posdiff_py_cut->GetXaxis()->SetRangeUser(blcuts::blc2bpc_y_min,blcuts::blc2bpc_y_max);
  h2_BLC2_BPC_posdiff_py_cut->SetFillColor(kGreen);
  h2_BLC2_BPC_posdiff_py_cut->Draw("same");



  TH2F* h2_BLC2_BPC_dirdiff = (TH2F*)f->Get("dydzdxdz_BLC2BPC");
  TCanvas *cBLC2_BPC_dirdiff = new TCanvas("cBLC2_BPC_dirdiff","BLC2_BPC_dirdiff");
  cBLC2_BPC_dirdiff->cd();
  h2_BLC2_BPC_dirdiff->SetXTitle("BCL2BPC track dx/dz");
  h2_BLC2_BPC_dirdiff->SetYTitle("BCL2BPC track dy/dz");
  h2_BLC2_BPC_dirdiff->Draw("colz");
   
  TCanvas *cBLC2_BPC_dirdiff_px = new TCanvas("cBLC2_BPC_dirdiff_px","BLC2_BPC_dirdiff_px");
  cBLC2_BPC_dirdiff_px->cd();
  TH1D* h2_BLC2_BPC_dirdiff_px = h2_BLC2_BPC_dirdiff->ProjectionX();
  h2_BLC2_BPC_dirdiff_px->Draw();
  TH1D* h2_BLC2_BPC_dirdiff_px_cut = (TH1D*) h2_BLC2_BPC_dirdiff_px->Clone();
  h2_BLC2_BPC_dirdiff_px_cut->GetXaxis()->SetRangeUser(blcuts::blc2bpc_dx_min,blcuts::blc2bpc_dx_max);
  h2_BLC2_BPC_dirdiff_px_cut->SetFillColor(kGreen);
  h2_BLC2_BPC_dirdiff_px_cut->Draw("same");
  
  TCanvas *cBLC2_BPC_dirdiff_py = new TCanvas("cBLC2_BPC_dirdiff_py","BLC2_BPC_dirdiff_py");
  cBLC2_BPC_dirdiff_py->cd();
  TH1D* h2_BLC2_BPC_dirdiff_py = h2_BLC2_BPC_dirdiff->ProjectionY();
  h2_BLC2_BPC_dirdiff_py->Draw();
  TH1D* h2_BLC2_BPC_dirdiff_py_cut = (TH1D*) h2_BLC2_BPC_dirdiff_py->Clone();
  h2_BLC2_BPC_dirdiff_py_cut->GetXaxis()->SetRangeUser(blcuts::blc2bpc_dy_min,blcuts::blc2bpc_dy_max);
  h2_BLC2_BPC_dirdiff_py_cut->SetFillColor(kGreen);
  h2_BLC2_BPC_dirdiff_py_cut->Draw("same");

  return;
}
