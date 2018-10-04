const bool gridon=true;
const bool staton=true;

namespace anacuts{
  const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
  const double dE_MIN = 2.0; // 8.0MeVee * 3cm / 5cm;

  const double pipi_MIN = 0.485;
  const double pipi_MAX = 0.510;
  const double ppi_MIN = 1.1075;
  const double ppi_MAX = 1.1225;

  const double neutron_MIN = 0.85;
  const double neutron_MAX = 1.03;

  const double Sigmap_MIN = 1.17;
  const double Sigmap_MAX = 1.21;
  const double Sigmam_MIN = 1.18;
  const double Sigmam_MAX = 1.22;
}

void plothists_v15(const char *filename="evanaIMpisigma_all_v15.root")
{
  cout << "filename " << filename << endl;
  gStyle->SetPadGridX(gridon);
  gStyle->SetPadGridY(gridon);
  gStyle->SetTitleYOffset(1.6);
  if(staton)gStyle->SetOptStat("emruo");
  else      gStyle->SetOptStat(0);
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);
  
  //QAbeamline(filename);

  //CDS QA 
  QACDS(filename);

  //Forward QA
  //QAForward(filename);

  if(staton)gStyle->SetOptStat("emruo");
  else      gStyle->SetOptStat(0);
  
  //PhysicsPlots(filename);

  return;
}

void QAForward(const char* filename){
  
  TFile *f = new TFile(filename,"READ"); 
  gStyle->SetTitleYOffset(1.6);
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.03);
  
  TCanvas *cNC_overbeta_dE = new TCanvas("NC_overbeta_dE","NC_overbeta_dE");
  cNC_overbeta_dE->cd();
  TH2F* h2_NC_overbeta_dE = f->Get("NC_overbeta_dE");
  gPad->SetLogz();
  h2_NC_overbeta_dE->SetXTitle("1/#beta");
  h2_NC_overbeta_dE->GetXaxis()->CenterTitle();
  h2_NC_overbeta_dE->SetYTitle("dE [MeV]");
  h2_NC_overbeta_dE->GetYaxis()->CenterTitle();
  h2_NC_overbeta_dE->Draw("colz");
  
  TCanvas *cNC_overbeta = new TCanvas("NC_overbeta","NC_overbeta");
  cNC_overbeta->cd();
  gPad->SetLogz(0);
  gPad->SetLogy();
  TH1D* h1_NC_overbeta = h2_NC_overbeta_dE->ProjectionX("px",9,200);
  h1_NC_overbeta->SetMinimum(50);
  h1_NC_overbeta->GetXaxis()->SetRangeUser(0.5,4.5);
  h1_NC_overbeta->Draw();

  return;
}



void QACDS(const char *filename){

  TFile *f = new TFile(filename,"READ"); 
  gStyle->SetTitleYOffset(1.6);
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);
  
  //CDH time distribution
  TCanvas *cCDHtime = new TCanvas("cCDHtime","CDHtime");
  cCDHtime->cd();
  TH2F* h2_CDHtime = f->Get("CDHtime");
  h2_CDHtime->SetXTitle("CDH seg#");
  h2_CDHtime->GetXaxis()->CenterTitle();
  h2_CDHtime->SetYTitle("time [nsec]");
  h2_CDHtime->GetYaxis()->CenterTitle();
  h2_CDHtime->Draw("colz");

  TCanvas *cCDHNeutralSeg = new TCanvas("cCDHNeutralSeg","cCDHNeutralSeg");
  cCDHNeutralSeg->cd();
  TH1F* h1_CDHNeutralSeg = f->Get("CDHNeutralSeg");
  h1_CDHNeutralSeg->SetXTitle("CDH seg#");
  h1_CDHNeutralSeg->GetXaxis()->CenterTitle();
  h1_CDHNeutralSeg->Draw();

  
  //CDC track chi2 distribution
  TCanvas *ctrackchi2_CDC = new TCanvas("ctrackchi2_CDC","ctrackchi2_CDC");
  TH1F* h1_trackchi2_CDC = f->Get("trackchi2_CDC");
  h1_trackchi2_CDC->SetXTitle("CDC track chi2/ndf");
  h1_trackchi2_CDC->GetXaxis()->CenterTitle();
  h1_trackchi2_CDC->Draw();


  TCanvas *cntrack = new TCanvas("cntrac","cntrack");
  TH1F* h1_ntrack_CDS = f->Get("ntrack_CDS");
  //h1_ntrack_CDS->SetXTitle("# of tracks");
  //h1_ntrack_CDS->GetXaxis()->CenterTitle();
  //h1_ntrack_CDS->Draw();
  TH1F* h1_ntrack_pi_plus = f->Get("ntrack_pi_plus");
  h1_ntrack_pi_plus->SetTitle("");
  h1_ntrack_pi_plus->SetLineColor(1);
  h1_ntrack_pi_plus->SetXTitle("# of tracks per event");
  h1_ntrack_pi_plus->GetXaxis()->CenterTitle();
  h1_ntrack_pi_plus->Draw("");
  TH1F* h1_ntrack_pi_minus = f->Get("ntrack_pi_minus");
  h1_ntrack_pi_minus->SetLineColor(2);
  h1_ntrack_pi_minus->Draw("same");
  TH1F* h1_ntrack_proton = f->Get("ntrack_proton");
  h1_ntrack_proton->SetLineColor(3);
  h1_ntrack_proton->Draw("same");
  TH1F* h1_ntrack_K_minus = f->Get("ntrack_K_minus");
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
  TH2F *h2_CDS_mom_betainv = f->Get("PID_CDS_beta");
  //h2_CDS_mom_betainv->GetXaxis()->SetRangeUser(-0.3,1.5);
  h2_CDS_mom_betainv->GetXaxis()->SetTitle("mass^{2} [(GeV/c^{2})^{2}]");
  h2_CDS_mom_betainv->GetXaxis()->CenterTitle();
  //h2_CDS_mom_betainv->GetYaxis()->SetRangeUser(-1.0,1.0);
  h2_CDS_mom_betainv->GetYaxis()->SetTitle("1/#beta");
  h2_CDS_mom_betainv->GetYaxis()->CenterTitle();
  gPad->SetLogz();
  h2_CDS_mom_betainv->Draw("colz");
  gStyle->SetOptStat(0);
  cCDS_mom_betainv->Update();
  cCDS_mom_betainv->Modified();

   
  TCanvas *cCDS_mass2_mom  = new TCanvas("CDS_mass2_mom","CDS_mass2_mom");
  TH2F *h2_CDS_mass2_mom = f->Get("PID_CDS");
  h2_CDS_mass2_mom->GetXaxis()->SetRangeUser(-0.3,1.5);
  h2_CDS_mass2_mom->GetXaxis()->SetTitle("mass^{2} [(GeV/c^{2})^{2}]");
  h2_CDS_mass2_mom->GetXaxis()->CenterTitle();
  h2_CDS_mass2_mom->GetYaxis()->SetRangeUser(-1.0,1.0);
  h2_CDS_mass2_mom->GetYaxis()->SetTitle("charge X momentum [GeV/c]");
  h2_CDS_mass2_mom->GetYaxis()->CenterTitle();
  gPad->SetLogz();
  h2_CDS_mass2_mom->Draw("colz");
  gStyle->SetOptStat(0);
  cCDS_mass2_mom->Update();
  cCDS_mass2_mom->Modified();
  
  
  TCanvas *cVtx_ZX  = new TCanvas("Vtx_ZX","Vtx_ZX");
  TH2F *h2_Vtx_ZX = f->Get("Vtx_ZX");
  h2_Vtx_ZX->SetTitle("CDS Vertex");
  h2_Vtx_ZX->GetXaxis()->SetTitle("Z Vertex [cm]");
  h2_Vtx_ZX->GetXaxis()->CenterTitle();
  h2_Vtx_ZX->GetYaxis()->SetTitle("X Vertex [cm]");
  h2_Vtx_ZX->GetYaxis()->CenterTitle();
  h2_Vtx_ZX->Draw("colz");
  double maxZX = h2_Vtx_ZX->GetMaximum();

  TCanvas *cVtx_ZY  = new TCanvas("Vtx_ZY","Vtx_ZY");
  TH2F *h2_Vtx_ZY = f->Get("Vtx_ZY");
  h2_Vtx_ZY->SetTitle("CDS Vertex");
  h2_Vtx_ZY->GetXaxis()->SetTitle("Z Vertex [cm]");
  h2_Vtx_ZY->GetXaxis()->CenterTitle();
  h2_Vtx_ZY->GetYaxis()->SetTitle("Y Vertex [cm]");
  h2_Vtx_ZY->GetYaxis()->CenterTitle();
  h2_Vtx_ZY->Draw("colz");
  double maxZY = h2_Vtx_ZY->GetMaximum();
  
  TCanvas *cVtx_XY  = new TCanvas("Vtx_XY","Vtx_XY");
  TH2F *h2_Vtx_XY = f->Get("Vtx_XY");
  h2_Vtx_XY->SetTitle("CDS Vertex");
  h2_Vtx_XY->GetXaxis()->SetTitle("X Vertex [cm]");
  h2_Vtx_XY->GetXaxis()->CenterTitle();
  h2_Vtx_XY->GetYaxis()->SetTitle("Y Vertex [cm]");
  h2_Vtx_XY->GetYaxis()->CenterTitle();
  h2_Vtx_XY->Draw("colz");
  double maxXY = h2_Vtx_XY->GetMaximum();
    
  //fiducial cuts
  TCanvas *cVtx_ZX_fid  = new TCanvas("Vtx_ZX_fid","Vtx_ZX_fid");
  TH2F *h2_Vtx_ZX_fid = f->Get("Vtx_ZX_fid");
  h2_Vtx_ZX_fid->SetTitle("CDS Vertex");
  h2_Vtx_ZX_fid->GetXaxis()->SetTitle("Z Vertex [cm]");
  h2_Vtx_ZX_fid->GetXaxis()->CenterTitle();
  h2_Vtx_ZX_fid->GetYaxis()->SetTitle("X Vertex [cm]");
  h2_Vtx_ZX_fid->GetYaxis()->CenterTitle();
  h2_Vtx_ZX_fid->SetMaximum(maxZX);
  h2_Vtx_ZX_fid->Draw("colz");
  
  TCanvas *cVtx_ZY_fid  = new TCanvas("Vtx_ZY_fid","Vtx_ZY_fid");
  TH2F *h2_Vtx_ZY_fid = f->Get("Vtx_ZY_fid");
  h2_Vtx_ZY_fid->SetTitle("CDS Vertex");
  h2_Vtx_ZY_fid->GetXaxis()->SetTitle("Z Vertex [cm]");
  h2_Vtx_ZY_fid->GetXaxis()->CenterTitle();
  h2_Vtx_ZY_fid->GetYaxis()->SetTitle("Y Vertex [cm]");
  h2_Vtx_ZY_fid->GetYaxis()->CenterTitle();
  h2_Vtx_ZY_fid->SetMaximum(maxZY);
  h2_Vtx_ZY_fid->Draw("colz");
  
  TCanvas *cVtx_XY_fid  = new TCanvas("Vtx_XY_fid","Vtx_XY_fid");
  TH2F *h2_Vtx_XY_fid = f->Get("Vtx_XY_fid");
  h2_Vtx_XY_fid->SetTitle("CDS Vertex");
  h2_Vtx_XY_fid->GetXaxis()->SetTitle("X Vertex [cm]");
  h2_Vtx_XY_fid->GetXaxis()->CenterTitle();
  h2_Vtx_XY_fid->GetYaxis()->SetTitle("Y Vertex [cm]");
  h2_Vtx_XY_fid->GetYaxis()->CenterTitle();
  h2_Vtx_XY_fid->SetMaximum(maxXY);
  h2_Vtx_XY_fid->Draw("colz");
  
  /*
  TCanvas *cVtx_X_center  = new TCanvas("Vtx_X_center","Vtx_X_center");
  TH1F *h1_Vtx_X_center = f->Get("Vtx_X_center");
  gPad->SetLogz(0);
  h1_Vtx_X_center->SetTitle("CDS Vertex");
  h1_Vtx_X_center->GetXaxis()->SetTitle("X Vertex [cm]");
  h1_Vtx_X_center->GetXaxis()->CenterTitle();
  h1_Vtx_X_center->GetYaxis()->SetTitle("Counts / 25.0 [#mum]");
  h1_Vtx_X_center->GetYaxis()->CenterTitle();
  h1_Vtx_X_center->Draw("");
  
  TCanvas *cVtx_Z_center  = new TCanvas("Vtx_Z_center","Vtx_Z_center");
  TH1F *h1_Vtx_Z_center = f->Get("Vtx_Z_center");
  h1_Vtx_Z_center->SetTitle("CDS Vertex");
  h1_Vtx_Z_center->GetXaxis()->SetTitle("Z Vertex [cm]");
  h1_Vtx_Z_center->GetXaxis()->CenterTitle();
  h1_Vtx_Z_center->GetYaxis()->SetTitle("Counts / 25.0 [#mum]");
  h1_Vtx_Z_center->GetYaxis()->CenterTitle();
  h1_Vtx_Z_center->Draw("");
  */

  TCanvas *cmul_CDH = new TCanvas("cmul_CDH","cmul_CDH");
  TH1F* h1_mul_CDH = f->Get("mul_CDH");
  h1_mul_CDH->SetTitle("CDH multiplicity");
  h1_mul_CDH->SetXTitle("# of CDH segment");
  h1_mul_CDH->GetXaxis()->CenterTitle();
  h1_mul_CDH->Draw();
  TH1F* h1_mul_CDH_clone = (TH1F*)h1_mul_CDH->Clone();
  h1_mul_CDH_clone->GetXaxis()->SetRangeUser(2.5,3.5);
  h1_mul_CDH_clone->SetFillColor(2);
  h1_mul_CDH_clone->Draw("same");

  TCanvas *cmul_CDH_assoc = new TCanvas("cmul_CDH_assoc","cmul_CDH_assoc");
  TH1F* h1_mul_CDH_assoc = f->Get("mul_CDH_assoc");
  h1_mul_CDH_assoc->SetTitle("# CDH seg. associated with CDC track");
  h1_mul_CDH_assoc->SetXTitle("# of CDH segment");
  h1_mul_CDH_assoc->GetXaxis()->CenterTitle();
  h1_mul_CDH_assoc->Draw();

  TCanvas *cdiff_CDH = new TCanvas("cdiff_CDH","cdiff_CDH");
  TH1F* h1_diff_CDH = f->Get("diff_CDH");
  h1_diff_CDH->SetTitle("Not used CDH segment # - CDH seg# used for charged");
  h1_diff_CDH->SetXTitle("diff of CDH seg#");
  h1_diff_CDH->Draw();
  
  TCanvas *cdiff_CDH_CDC = new TCanvas("cdiff_CDH_CDC","cdiff_CDH_CDC");
  TH1F* h1_cdiff_CDH_CDC = f->Get("diff_CDH_CDC");
  h1_cdiff_CDH_CDC->SetXTitle("angle [deg]");
  h1_cdiff_CDH_CDC->GetXaxis()->CenterTitle();
  h1_cdiff_CDH_CDC->Draw();

  TCanvas *cdE_betainv_fiducial = new TCanvas("cdE_betainv_fiducial","cdE_betainv_fiducial");
  TH2F* h2_dE_betainv = f->Get("dE_betainv_fid");
  h2_dE_betainv->SetXTitle("1/#beta");
  h2_dE_betainv->GetXaxis()->CenterTitle();
  h2_dE_betainv->SetYTitle("dE [MeV]");
  h2_dE_betainv->GetYaxis()->CenterTitle();
  h2_dE_betainv->Draw("colz");

  TCanvas *c_betainv_fiducial = new TCanvas("c_betainv_fiducial","c_betainv_fiducial");
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
  
  
  //TCanvas *c


  return;
}

void PhysicsPlots(const char *filename){

  TFile *f = new TFile(filename,"READ"); 
  //physics plots
  //IM pi+pi- when nCDH==3, fudicial OK,
  TCanvas *c_IMpipi = new TCanvas("c_IMpipi","c_IMpipi");
  TH1F* h1_IMpipi = f->Get("IMpipi_dE");
  c_IMpipi->cd();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.03);
  h1_IMpipi->SetTitle("Invariant Mass #pi^{+}#pi^{-}");
  h1_IMpipi->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  h1_IMpipi->GetXaxis()->CenterTitle() ;
  h1_IMpipi->SetYTitle("Counts/ 1 MeV");
  h1_IMpipi->GetYaxis()->CenterTitle() ;
  h1_IMpipi->GetYaxis()->SetTitleOffset(1.4);
  h1_IMpipi->Draw();
  TH1F* h1_IMpipi_clone = h1_IMpipi->Clone();
  h1_IMpipi_clone->SetFillColor(2);
  h1_IMpipi_clone->GetXaxis()->SetRangeUser(anacuts::pipi_MIN,anacuts::pipi_MAX);
  h1_IMpipi_clone->Draw("same");
  

  TCanvas *cMMom_MMass_fid_beta_dE = new TCanvas("cMMom_MMass_fid_beta_dE","");
  TH2F* h2_MMom_MMass_fid_beta_dE = f->Get("MMom_MMass_fid_beta_dE");
  h2_MMom_MMass_fid_beta_dE->SetXTitle("missing mass [GeV/c^{2}]");
  h2_MMom_MMass_fid_beta_dE->GetXaxis()->CenterTitle();
  h2_MMom_MMass_fid_beta_dE->SetYTitle("missing momentum [GeV/c]");
  h2_MMom_MMass_fid_beta_dE->GetYaxis()->CenterTitle();
  h2_MMom_MMass_fid_beta_dE->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_dE_px = new TCanvas("cMMom_MMass_fid_beta_dE_px","");
  TH1D* proMMass = h2_MMom_MMass_fid_beta_dE->ProjectionX();
  proMMass->Draw();
  
  TCanvas *cMMom_MMass_fid_beta = new TCanvas("cMMom_MMass_fid_beta","");
  TH2F* h2_MMom_MMass_fid_beta = f->Get("MMom_MMass_fid_beta");
  h2_MMom_MMass_fid_beta->SetXTitle("missing mass [GeV/c^{2}]");
  h2_MMom_MMass_fid_beta->GetXaxis()->CenterTitle();
  h2_MMom_MMass_fid_beta->SetYTitle("missing momentum [GeV/c]");
  h2_MMom_MMass_fid_beta->GetYaxis()->CenterTitle();
  h2_MMom_MMass_fid_beta->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_px = new TCanvas("cMMom_MMass_fid_beta_px","");
  TH1D* proMMass2 = h2_MMom_MMass_fid_beta->ProjectionX();
  proMMass2->Draw();

 
  TCanvas *cMMom_MMass_fid_beta_wSid = new TCanvas("cMMom_MMass_fid_beta_wSid","");
  TH2F* h2_MMom_MMass_fid_beta_wSid = f->Get("MMom_MMass_fid_beta_dE_woK0_wSid");
  h2_MMom_MMass_fid_beta_wSid->SetXTitle("missing mass [GeV/c^{2}]");
  h2_MMom_MMass_fid_beta_wSid->GetXaxis()->CenterTitle();
  h2_MMom_MMass_fid_beta_wSid->SetYTitle("missing momentum [GeV/c]");
  h2_MMom_MMass_fid_beta_wSid->GetYaxis()->CenterTitle();
  h2_MMom_MMass_fid_beta_wSid->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_wSid_px = new TCanvas("cMMom_MMass_fid_beta__wSid_px","");
  TH1D* proMMass2_wSid = h2_MMom_MMass_fid_beta_wSid->ProjectionX();
  proMMass2_wSid->Draw();


  
  TCanvas *cIMnpim_IMnpip = new TCanvas("cIMnpim_IMnpip","");
  TH2F* h2_IMnpim_IMnpip = f->Get("IMnpim_IMnpip_dE");
  h2_IMnpim_IMnpip->SetXTitle("IM(#pi^{+}n) [GeV/c^{2}]");
  h2_IMnpim_IMnpip->GetXaxis()->CenterTitle();
  h2_IMnpim_IMnpip->SetYTitle("IM(#pi^{-}n) [GeV/c^{2}]");
  h2_IMnpim_IMnpip->GetYaxis()->CenterTitle();
  h2_IMnpim_IMnpip->Draw("colz");
  

  TCanvas *cIMnpim_IMnpip_px = new TCanvas("cIMnpim_IMnpip_px","");
  TH1D* h2_IMnpim_IMnpip_px = h2_IMnpim_IMnpip->ProjectionX();
  h2_IMnpim_IMnpip_px->Draw();
  TH1D* h2_IMnpim_IMnpip_px_clone = h2_IMnpim_IMnpip_px->Clone();
  h2_IMnpim_IMnpip_px_clone->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  h2_IMnpim_IMnpip_px_clone->SetFillColor(3);
  h2_IMnpim_IMnpip_px_clone->Draw("same");

  TCanvas *cIMnpim_IMnpip_py = new TCanvas("cIMnpim_IMnpip_py","");
  TH1D* h2_IMnpim_IMnpip_py = h2_IMnpim_IMnpip->ProjectionY();
  h2_IMnpim_IMnpip_py->Draw();
  TH1D* h2_IMnpim_IMnpip_py_clone = h2_IMnpim_IMnpip_py->Clone();
  h2_IMnpim_IMnpip_py_clone->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  h2_IMnpim_IMnpip_py_clone->SetFillColor(4);
  h2_IMnpim_IMnpip_py_clone->Draw("same");


  //
  TCanvas *cIMnpim_IMnpip_wmn = new TCanvas("cIMnpim_IMnpip_wmn","");
  TH2F* h2_IMnpim_IMnpip_wmn = f->Get("IMnpim_IMnpip_wmn_dE");
  h2_IMnpim_IMnpip_wmn->SetXTitle("IM(#pi^{+}n) [GeV/c^{2}]");
  h2_IMnpim_IMnpip_wmn->GetXaxis()->CenterTitle();
  h2_IMnpim_IMnpip_wmn->SetYTitle("IM(#pi^{-}n) [GeV/c^{2}]");
  h2_IMnpim_IMnpip_wmn->GetYaxis()->CenterTitle();
  h2_IMnpim_IMnpip_wmn->Draw("colz");
  

  TCanvas *cIMnpim_IMnpip_wmn_px = new TCanvas("cIMnpim_IMnpip_wmn_px","");
  TH1D* h2_IMnpim_IMnpip_wmn_px = h2_IMnpim_IMnpip_wmn->ProjectionX();
  h2_IMnpim_IMnpip_wmn_px->Draw();
  TH1D* h2_IMnpim_IMnpip_wmn_px_clone = h2_IMnpim_IMnpip_wmn_px->Clone();
  h2_IMnpim_IMnpip_wmn_px_clone->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  h2_IMnpim_IMnpip_wmn_px_clone->SetFillColor(3);
  h2_IMnpim_IMnpip_wmn_px_clone->Draw("same");

  TCanvas *cIMnpim_IMnpip_wmn_py = new TCanvas("cIMnpim_IMnpip_wmn_py","");
  TH1D* h2_IMnpim_IMnpip_wmn_py = h2_IMnpim_IMnpip_wmn->ProjectionY();
  h2_IMnpim_IMnpip_wmn_py->Draw();
  TH1D* h2_IMnpim_IMnpip_wmn_py_clone = h2_IMnpim_IMnpip_wmn_py->Clone();
  h2_IMnpim_IMnpip_wmn_py_clone->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  h2_IMnpim_IMnpip_wmn_py_clone->SetFillColor(4);
  h2_IMnpim_IMnpip_wmn_py_clone->Draw("same");


  
  TCanvas *cIMmnpim_IMmnpip = new TCanvas("cIMmnpim_IMmnpip","");
  TH2F* h2_IMmnpim_IMmnpip = f->Get("IMmnpim_IMmnpip");
  h2_IMmnpim_IMmnpip->SetXTitle("IM(#pi^{+}n) [GeV/c^{2}]");
  h2_IMmnpim_IMmnpip->GetXaxis()->CenterTitle();
  h2_IMmnpim_IMmnpip->SetYTitle("IM(#pi^{-}n) [GeV/c^{2}]");
  h2_IMmnpim_IMmnpip->GetYaxis()->CenterTitle();
  h2_IMmnpim_IMmnpip->Draw("colz");
  
  TCanvas *cIMmnpim_IMmnpip_px = new TCanvas("cIMmnpim_IMmnpip_px","");
  TH1D* h2_IMmnpim_IMmnpip_px = h2_IMmnpim_IMmnpip->ProjectionX();
  h2_IMmnpim_IMmnpip_px->Draw();
  //TH1D* h2_IMmnpim_IMmnpip_px_clone = h2_IMmnpim_IMmnpip_px->Clone();
  //h2_IMmnpim_IMmnpip_px_clone->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  //h2_IMmnpim_IMmnpip_px_clone->SetFillColor(3);
  //h2_IMmnpim_IMmnpip_px_clone->Draw("same");

  TCanvas *cIMmnpim_IMmnpip_py = new TCanvas("cIMmnpim_IMmnpip_py","");
  TH1D* h2_IMmnpim_IMmnpip_py = h2_IMmnpim_IMmnpip->ProjectionY();
  h2_IMmnpim_IMmnpip_py->Draw();
  //TH1D* h2_IMmnpim_IMmnpip_py_clone = h2_IMmnpim_IMmnpip_py->Clone();
  //h2_IMmnpim_IMmnpip_py_clone->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  //h2_IMmnpim_IMmnpip_py_clone->SetFillColor(4);
  //h2_IMmnpim_IMmnpip_py_clone->Draw("same");

  TCanvas *cIMnpipi_wSid = new TCanvas("cIMnpipi_wSid","");
  //TH1F* h1_IMnpipi_wSid = f->Get("IMnpipi_wSid");
  //h1_IMnpipi_wSid->SetXTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  //h1_IMnpipi_wSid->GetXaxis()->CenterTitle();
  //h1_IMnpipi_wSid->SetYTitle("Counts per 10 MeV/c^{2}");
  //h1_IMnpipi_wSid->GetYaxis()->CenterTitle();
  //h1_IMnpipi_wSid->Draw("HE");
  //TH1F* h1_IMnpipi_woK0_wSid = f->Get("IMnpipi_woK0_wSid");
  //h1_IMnpipi_woK0_wSid->SetLineColor(2);
  //h1_IMnpipi_woK0_wSid->Draw("HEsame");
  
  TH1F* h1_IMnpipi_woK0_wSid = f->Get("IMnpipi_woK0_wSid");
  h1_IMnpipi_woK0_wSid->SetXTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  h1_IMnpipi_woK0_wSid->GetXaxis()->CenterTitle();
  h1_IMnpipi_woK0_wSid->SetYTitle("Counts per 10 MeV/c^{2}");
  h1_IMnpipi_woK0_wSid->GetYaxis()->CenterTitle();
  h1_IMnpipi_woK0_wSid->Draw("HE");




  TCanvas *cIMnpipi_MMnmiss_woK0_wSid = new TCanvas("cIMnpipi_MMnmiss_woK0_wSid","");
  TH2F* h2_IMnpipi_MMnmiss_woK0_wSid = f->Get("IMnpipi_MMnmiss_woK0_wSid");
  h2_IMnpipi_MMnmiss_woK0_wSid->SetXTitle("Missing Momentum [GeV/c]");
  h2_IMnpipi_MMnmiss_woK0_wSid->GetXaxis()->CenterTitle();
  h2_IMnpipi_MMnmiss_woK0_wSid->SetYTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  h2_IMnpipi_MMnmiss_woK0_wSid->GetYaxis()->CenterTitle();
  h2_IMnpipi_MMnmiss_woK0_wSid->Draw("colz");
  
  TCanvas *cCosn_IMnpipi = new TCanvas("cCosn_IMnpipi","");
  TH2F* h2_cCosn_IMnpipi = f->Get("Cosn_IMnpipi");
  h2_cCosn_IMnpipi->SetXTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  h2_cCosn_IMnpipi->GetXaxis()->CenterTitle();
  h2_cCosn_IMnpipi->SetYTitle("cos#theta_{CM}");
  h2_cCosn_IMnpipi->GetYaxis()->CenterTitle();
  h2_cCosn_IMnpipi->Draw("colz");
  
  
  TCanvas *cIMnpipi_q_woK0_wSid = new TCanvas("cIMnpipi_q_woK0_wSid","");
  TH2F* h2_IMnpipi_q_woK0_wSid = f->Get("IMnpipi_q_woK0_wSid");
  h2_IMnpipi_q_woK0_wSid->SetXTitle("mom transfer [GeV/c]");
  h2_IMnpipi_q_woK0_wSid->SetYTitle("IM(n#pi^{-}#pi^{+}) [GeV/c^{2}]");
  h2_IMnpipi_q_woK0_wSid->GetXaxis()->CenterTitle();
  h2_IMnpipi_q_woK0_wSid->GetYaxis()->CenterTitle();
  h2_IMnpipi_q_woK0_wSid->Draw("colz");








}


void QAbeamline(const char *filename){
  
  TFile *f = new TFile(filename,"READ"); 
  gStyle->SetTitleYOffset(1.6);
  gStyle->SetOptStat("emruo");
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);

  //QA plots
  TH1F *h1_nBHD = f->Get("mul_BHD");
  TCanvas *cnBHD =  new TCanvas("nBHD","nBHD");
  cnBHD->cd();
  h1_nBHD->Draw();
  
  TH1F *h1_nT0 = f->Get("mul_T0");
  TCanvas *cnT0 =  new TCanvas("nT0","nT0");
  cnT0->cd();
  h1_nT0->Draw();

  TH1F *h1_T0BHDtof = f->Get("tof_T0BHD");
  TCanvas *ctof = new TCanvas("tof","tof");
  ctof->cd();
  h1_T0BHDtof->GetXaxis()->SetRangeUser(24,38);
  h1_T0BHDtof->Draw();
  TH1F *h1_T0BHDtof_cut = (TH1F*) h1_T0BHDtof->Clone();
  h1_T0BHDtof_cut->GetXaxis()->SetRangeUser(27.8588,29.5663);
  h1_T0BHDtof_cut->SetFillColor(30);
  h1_T0BHDtof_cut->Draw("same");

  //BLC1
  TH1F *h1_nBLC1 = f->Get("ntrack_BLC1");
  TCanvas *cnBLC1 = new TCanvas("nBLC1","nBLC1");
  cnBLC1->cd();
  h1_nBLC1->Draw();
  
  TH1F *h1_BLC1time = f->Get("tracktime_BLC1");
  TCanvas *cBLC1time = new TCanvas("BLC1time","BLC1time");
  cBLC1time->cd();
  h1_BLC1time->Draw();

  TH1F *h1_BLC1chi2 = f->Get("trackchi2_BLC1");
  TCanvas *cBLC1chi2 = new TCanvas("BLC1chi2","BLC1chi2");
  cBLC1chi2->cd();
  h1_BLC1chi2->Draw();


  //BLC2
  TH1F *h1_nBLC2 = f->Get("ntrack_BLC2");
  TCanvas *cnBLC2 = new TCanvas("nBLC2","nBLC2");
  cnBLC2->cd();
  h1_nBLC2->Draw();

  TH1F *h1_BLC2time = f->Get("tracktime_BLC2");
  TCanvas *cBLC2time = new TCanvas("BLC2time","BLC2time");
  cBLC2time->cd();
  h1_BLC2time->Draw();

  TH1F *h1_BLC2chi2 = f->Get("trackchi2_BLC2");
  TCanvas *cBLC2chi2 = new TCanvas("BLC2chi2","BLC2chi2");
  cBLC2chi2->cd();
  h1_BLC2chi2->Draw();

  //BPC
  TH1F *h1_nBPC = f->Get("ntrack_BPC");
  TCanvas *cnBPC = new TCanvas("nBPC","nBPC");
  cnBPC->cd();
  h1_nBPC->Draw();

  TH1F *h1_BPCtime = f->Get("tracktime_BPC");
  TCanvas *cBPCtime = new TCanvas("BPCtime","BPCtime");
  cBPCtime->cd();
  h1_BPCtime->Draw();

  TH1F *h1_BPCchi2 = f->Get("trackchi2_BPC");
  TCanvas *cBPCchi2 = new TCanvas("BPCchi2","BPCchi2");
  cBPCchi2->cd();
  h1_BPCchi2->Draw();
  
  //D5
  TH1F *h1_D5chi2 = f->Get("trackchi2_beam");
  TCanvas *cD5chi2 = new TCanvas("D5chi2","D5chi2");
  cD5chi2->cd();
  h1_D5chi2->Draw();

  TH1F *h1_D5mom = f->Get("momentum_beam");
  TCanvas *cDmom = new TCanvas("D5mom","D5mom");
  cDmom->cd();
  h1_D5mom->Draw();
  
  //BLC2-BPD matching
  TH2F* h2_BLC2_BPC_posdiff = f->Get("dydx_BLC2BPC");
  TCanvas *cBLC2_BPC_posdiff = new TCanvas("cBLC2_BPC_posdiff","cBLC2_BPC_posdiff");
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.03);
  cBLC2_BPC_posdiff->cd();
  h2_BLC2_BPC_posdiff->Draw("colz");

  TH2F* h2_BLC2_BPC_dirdiff = f->Get("dydzdxdz_BLC2BPC");
  TCanvas *cBLC2_BPC_dirdiff = new TCanvas("cBLC2_BPC_dirdiff","cBLC2_BPC_dirdiff");
  cBLC2_BPC_dirdiff->cd();
  h2_BLC2_BPC_dirdiff->Draw("colz");

  return;
}
