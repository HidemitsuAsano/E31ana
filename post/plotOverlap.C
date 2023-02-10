//2D fit for K0nn modeling 
//x[0]: IM(npi+)
//x[1]: IM(npi-)
Double_t K0fit2dNoconvert(Double_t *x,Double_t *par)
{
  double xx = 1.0/sqrt(2.0)*(x[0]-x[1]);
  double yy = 1.0/sqrt(2.0)*(x[0]+x[1]);
  //double yy2 = yy-(sqrt(6.76*xx*xx+2.725)-1.0);
  double yy2 = yy-(sqrt(6.76*xx*xx+2.765)-1.0);//slightly tuned from original value by looking final fitting result

  Double_t r1 = (xx-par[1])/par[2];
  Double_t r2 = yy2*par[3];
  const Double_t th = 1.00;
  const Double_t a  = 0.002;
  const Double_t thend = 1.95;
  const Double_t aend = 0.05;
  double ret = par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(-1.0*r2*r2)/(1.0+TMath::Exp((-yy2+th)/a))/(1.0+TMath::Exp((yy-thend)/aend));    
  //double ret = par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(-1.0*r2*r2)/(1.0+TMath::Exp((-yy2+th)/a))/(1.0+TMath::Exp((yy-thend)/aend));    
  return ret;
  /*
  //if(yy2>1.0){
  if(yy2>1.0){
    return  ret;
   }else{
    return 0;
  }*/
}
const int Version = 245;

void plotOverlap()
{
  TFile *fr = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE2_iso_nostop_sys0_sub.root",Version),"READ");
  TFile *flo = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE2_iso_qlo_nostop_sys0_sub.root",Version),"READ");
  TFile *fhi = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE2_iso_qhi_nostop_sys0_sub.root",Version),"READ");
  TFile *ffitlo = TFile::Open(Form("fout_qlo_v%d_dE2_sys0.root",Version),"READ");
  TFile *ffithi = TFile::Open(Form("fout_qhi_v%d_dE2_sys0.root",Version),"READ");

  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  gStyle->SetOptStat(0);   
  IMnpim_IMnpip_dE_wK0orwSid_n->Rebin2D(4,4);
  TF2* f3widelo = (TF2*)ffitlo->Get("f3wide"); 
  TF2* f3widehi = (TF2*)ffithi->Get("f3wide"); 
  TH2F* hf3widelo = (TH2F*)f3widelo->GetHistogram();
  TH2F* hf3widehi = (TH2F*)f3widehi->GetHistogram();

  TH2F* hf3add = (TH2F*)hf3widelo->Clone("hf3add");
  hf3add->Add(hf3widehi);
  hf3add->Print("base");  
  IMnpim_IMnpip_dE_wK0orwSid_n->Print("base");
  const int nx = IMnpim_IMnpip_dE_wK0orwSid_n->GetNbinsX();
  const int ny = IMnpim_IMnpip_dE_wK0orwSid_n->GetNbinsY();
  const double xmin = IMnpim_IMnpip_dE_wK0orwSid_n->GetXaxis()->GetXmin();
  const double xmax = IMnpim_IMnpip_dE_wK0orwSid_n->GetXaxis()->GetXmax();
  const double ymin = IMnpim_IMnpip_dE_wK0orwSid_n->GetYaxis()->GetXmin();
  const double ymax = IMnpim_IMnpip_dE_wK0orwSid_n->GetYaxis()->GetXmax();
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_mev =  new TH2F("IMnpim_IMnpip_dE_wK0orwSid_n_mev","",nx,xmin*1000.,xmax*1000.,ny,ymin*1000.,ymax*1000.);
  
  for(int ix=1;ix<nx;ix++){
    for(int iy=1;iy<ny;iy++){
      double cont = IMnpim_IMnpip_dE_wK0orwSid_n->GetBinContent(ix,iy);
      IMnpim_IMnpip_dE_wK0orwSid_n_mev->SetBinContent(ix,iy,cont);
    }
  }

  const int fnx = hf3add->GetNbinsX();
  const int fny = hf3add->GetNbinsY();
  const double fxmin = hf3add->GetXaxis()->GetXmin();
  const double fxmax = hf3add->GetXaxis()->GetXmax();
  const double fymin = hf3add->GetYaxis()->GetXmin();
  const double fymax = hf3add->GetYaxis()->GetXmax();
  TH2F* hf3add_mev =  new TH2F("hf3add_mev","",fnx,fxmin*1000.,fxmax*1000.,fny,fymin*1000.,fymax*1000.);
  for(int ix=1;ix<fnx;ix++){
    for(int iy=1;iy<fny;iy++){
      double cont = hf3add->GetBinContent(ix,iy);
      hf3add_mev->SetBinContent(ix,iy,cont);
    }
  }

  const bool gridon=false;
  gStyle->SetPadGridX(gridon);
  gStyle->SetPadGridY(gridon);
  gStyle->SetTitleYOffset(1.6);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetNdivisions(505,"x");
  gStyle->SetNdivisions(505,"y");
  gStyle->SetNdivisions(505,"z");
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->SetTitle("");
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetXaxis()->SetTitle("IM(n#pi^{+}) [MeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetXaxis()->CenterTitle();
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetYaxis()->SetTitle("IM(n#pi^{-}) [MeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetYaxis()->CenterTitle();
  //IMnpim_IMnpip_dE_wK0orwSid_n->GetXaxis()->SetRangeUser(0.9,1.7);
  //IMnpim_IMnpip_dE_wK0orwSid_n->GetYaxis()->SetRangeUser(0.9,1.7);
  //IMnpim_IMnpip_dE_wK0orwSid_n->GetYaxis()->SetRangeUser(0.9,1.7);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetXaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetYaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetXaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetYaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetXaxis()->SetTitleOffset(1.31);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->GetYaxis()->SetTitleOffset(1.42);
  IMnpim_IMnpip_dE_wK0orwSid_n_mev->Draw("colz");
  hf3add_mev->SetLineColor(2);
  //hf3add_mev->SetLineWidth(1);
  hf3add_mev->SetMinimum(-1.);
  hf3add_mev->Draw("cont2same");
  c1->SaveAs("K0orwSid.pdf");

  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_qlo = (TH2F*)flo->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->Rebin2D(2,2);
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->SetTitle("");
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetXaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetYaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetXaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetYaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetXaxis()->SetTitleOffset(1.31);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetYaxis()->SetTitleOffset(1.42);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->Draw("colz");
  hf3widelo->Draw("cont2same");
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_qhi = (TH2F*)fhi->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->Rebin2D(2,2);
  TCanvas *c3 = new TCanvas("c3","",1000,800);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->SetTitle("");
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetXaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetYaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetXaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetYaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetXaxis()->SetTitleOffset(1.31);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetYaxis()->SetTitleOffset(1.42);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->Draw("colz");
  hf3widehi->Draw("cont2same");
}
