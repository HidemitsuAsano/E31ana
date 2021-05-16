#include "anacuts.h"

void SpSmDecoError(const int qcut=2)
{
  TFile *f = NULL;
  if(qcut==1){
    f = TFile::Open("fout_qlo.root","READ");
  }else if(qcut==2){
    f = TFile::Open("fout_qhi.root","READ");
  }else{
    std::cout << "no file" << std::endl;
    return;
  }

  TH1D* IMnpip_K0sub_woSp = (TH1D*)f->Get("IMnpip_K0sub_woSp"); 
  TH1D* IMnpim_K0sub_woSm = (TH1D*)f->Get("IMnpim_K0sub_woSm");
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  c1->Divide(2,1);
  c1->cd(1);
  IMnpip_K0sub_woSp->Draw("HE");
  const int Ntry = 10000;
  TH1D* IMnpip_K0sub_woSp_est = (TH1D*)IMnpip_K0sub_woSp->Clone("IMnpip_K0sub_woSp_est");
  TH1D* IMnpim_K0sub_woSm_est = (TH1D*)IMnpim_K0sub_woSm->Clone("IMnpim_K0sub_woSm_est");
  TGraph *gIMnpip_all = new TGraph();
  TGraph *gIMnpim_all = new TGraph();
  const int nbinSp = IMnpip_K0sub_woSp_est->GetNbinsX();
  const double xminSp = IMnpip_K0sub_woSp_est->GetXaxis()->GetXmin();
  const double xmaxSp = IMnpip_K0sub_woSp_est->GetXaxis()->GetXmax();
  TH2D* Est_IMnpip_woSp_pol1 = new TH2D("Est_IMnpip_woSp_pol1","Est_IMnpip_woSp_pol1",nbinSp,xminSp,xmaxSp,220,0,220);
  TH2D* Est_IMnpip_woSp_3rd = new TH2D("Est_IMnpip_woSp_3rd","Est_IMnpip_woSp_3rd",nbinSp,xminSp,xmaxSp,220,0,220);
  const int nbinSm = IMnpim_K0sub_woSm_est->GetNbinsX();
  const double xminSm = IMnpim_K0sub_woSm_est->GetXaxis()->GetXmin();
  const double xmaxSm = IMnpim_K0sub_woSm_est->GetXaxis()->GetXmax();
  TH2D* Est_IMnpim_woSm_pol1 = new TH2D("Est_IMnpim_woSm_pol1","Est_IMnpim_woSm_pol1",nbinSm,xminSm,xmaxSm,220,0,220);
  TH2D* Est_IMnpim_woSm_3rd = new TH2D("Est_IMnpim_woSm_3rd","Est_IMnpim_woSm_3rd",nbinSm,xminSm,xmaxSm,220,0,220);
  

  const int Spbin = IMnpip_K0sub_woSp->GetXaxis()->FindBin(anacuts::Sigmap_center);
  const int Smbin = IMnpim_K0sub_woSm->GetXaxis()->FindBin(anacuts::Sigmam_center);
  for(int it=0;it<Ntry;it++){
    if(it%1000==0) std::cout << "itry " << it << std::endl;
    IMnpip_K0sub_woSp_est->Reset();
    IMnpim_K0sub_woSm_est->Reset();
    TGraph *gIMnpip = new TGraph();
    TGraph *gIMnpim = new TGraph();
    for(int ibin=0;ibin<IMnpip_K0sub_woSp->GetNbinsX();ibin++){
      double cont = IMnpip_K0sub_woSp->GetBinContent(ibin);
      double bincent = IMnpip_K0sub_woSp->GetBinCenter(ibin);
      double err = IMnpip_K0sub_woSp->GetBinError(ibin);
      double gen = gRandom->Gaus(cont,err);
      if(cont<0.0001) continue;
      IMnpip_K0sub_woSp_est->SetBinContent(ibin,gen);
      IMnpip_K0sub_woSp_est->SetBinError(ibin,err);
      gIMnpip_all->AddPoint(bincent,gen);
      gIMnpip->AddPoint(bincent,gen);
    }
    TSpline3 *sIMnpip = new TSpline3("snpip",gIMnpip);
    const double Splowbincen = IMnpip_K0sub_woSp->GetBinCenter(Spbin-1);
    const double Sphighbincen = IMnpip_K0sub_woSp->GetBinCenter(Spbin+1);
    TF1 *fSp = new TF1("fSp","pol1",Splowbincen,Sphighbincen);
    IMnpip_K0sub_woSp_est->Fit("fSp","q","",Splowbincen,Sphighbincen);
    double est_Sp_pol1 = fSp->Eval(anacuts::Sigmap_center);
    Est_IMnpip_woSp_pol1->Fill(anacuts::Sigmap_center,est_Sp_pol1);
    double est_Sp_3rd  = sIMnpip->Eval(anacuts::Sigmap_center); 
    Est_IMnpip_woSp_3rd->Fill(anacuts::Sigmap_center,est_Sp_3rd);

    //IMnpip_K0sub_woSp_est->Draw("H");
    //sIMnpip->Draw("same");
    for(int ibin=0;ibin<IMnpim_K0sub_woSm->GetNbinsX();ibin++){
      double cont = IMnpim_K0sub_woSm->GetBinContent(ibin);
      double bincent = IMnpim_K0sub_woSm->GetBinCenter(ibin);
      double err = IMnpim_K0sub_woSm->GetBinError(ibin);
      double gen = gRandom->Gaus(cont,err);
      //std::cout << bincent << " " << gen << " " << cont <<  std::endl;
      if(cont<0.0001) continue;
      IMnpim_K0sub_woSm_est->SetBinContent(ibin,gen);
      IMnpim_K0sub_woSm_est->SetBinError(ibin,err);
      gIMnpim_all->AddPoint(bincent,gen);
      gIMnpim->AddPoint(bincent,gen);
    }
    TSpline3 *sIMnpim = new TSpline3("snpim",gIMnpim);
    const double Smlowbincen = IMnpim_K0sub_woSm->GetBinCenter(Smbin-1);
    const double Smhighbincen = IMnpim_K0sub_woSm->GetBinCenter(Smbin+1);
    TF1 *fSm = new TF1("fSm","pol1",Smlowbincen,Smhighbincen);
    IMnpim_K0sub_woSm_est->Fit("fSm","q","",Smlowbincen,Smhighbincen);
    double est_Sm_pol1 = fSm->Eval(anacuts::Sigmam_center);
    Est_IMnpim_woSm_pol1->Fill(anacuts::Sigmam_center,est_Sm_pol1);
    double est_Sm_3rd  = sIMnpim->Eval(anacuts::Sigmam_center); 
    Est_IMnpim_woSm_3rd->Fill(anacuts::Sigmam_center,est_Sm_3rd);
    //IMnpim_K0sub_woSm_est->Draw("H");
    //sIMnpim->Draw("same");
    //gIMnpim->Print();
    //sIMnpim->Print();
    //break;
  }
  
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(2,1);
  c2->cd(1);
  gIMnpip_all->SetMarkerStyle(20);
  gIMnpip_all->Draw("ap");
  c2->cd(2);
  gIMnpim_all->SetMarkerStyle(20);
  gIMnpim_all->Draw("ap");
  
  TCanvas *c3 = new TCanvas("c3","c3");
  Est_IMnpip_woSp_pol1->Draw("colz");
  double stddev_Sp_pol1 = Est_IMnpip_woSp_pol1->GetStdDev(2);
  std::cout  << Est_IMnpip_woSp_pol1->GetStdDev(2) << std::endl;

  TCanvas *c4 = new TCanvas("c4","c4");
  Est_IMnpip_woSp_3rd->Draw("colz");
  double stddev_Sp_3rd = Est_IMnpip_woSp_3rd->GetStdDev(2);
  std::cout << Est_IMnpip_woSp_3rd->GetStdDev(2) << std::endl;
  
  TCanvas *c5 = new TCanvas("c5","c5");
  Est_IMnpim_woSm_pol1->Draw("colz");
  double stddev_Sm_pol1 = Est_IMnpim_woSm_pol1->GetStdDev(2);
  std::cout << Est_IMnpim_woSm_pol1->GetStdDev(2) << std::endl;

  TCanvas *c6 = new TCanvas("c6","c6");
  Est_IMnpim_woSm_3rd->Draw("colz");
  double stddev_Sm_3rd = Est_IMnpim_woSm_3rd->GetStdDev(2);
  std::cout << Est_IMnpim_woSm_3rd->GetStdDev(2) << std::endl;

  TCanvas *c7 = new TCanvas("c7","c7");
  IMnpip_K0sub_woSp->Draw("H");

  TCanvas *c8 = new TCanvas("c8","c8");
  IMnpim_K0sub_woSm->Draw("H");
}
