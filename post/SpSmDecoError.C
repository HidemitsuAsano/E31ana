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
  TH1D* IMnpip_K0sub_woSp_est = IMnpip_K0sub_woSp->Clone("IMnpip_K0sub_woSp_est");
  TH1D* IMnpim_K0sub_woSm_est = IMnpim_K0sub_woSm->Clone("IMnpim_K0sub_woSm_est");

  for(int it=0;it<Ntry;it++){
    IMnpip_K0sub_woSp_est->Reset();
    IMnpim_K0sub_woSm_est->Reset();
    for(int ibin=0;ibin<IMnpip_K0sub_woSp->GetBbinsX();ibin++){
      double cont = IMnpip_K0sub_woSp->GetBinContent(ibin);
      double err = IMnpip_K0sub_woSp->GetBinError(ibin);
      double gen = gRandom->Gaus(cont,err);
    }

    c1->cd(2);
    IMnpim_K0sub_woSm->Draw("HE");
    for(int ibin=0;ibin<IMnpim_K0sub_woSm->GetBbinsX();ibin++){
      double cont = IMnpim_K0sub_woSm->GetBinContent(ibin);
      double err = IMnpim_K0sub_woSm->GetBinError(ibin);
      double gen = gRandom->Gaus(cont,err);
    }
  }

}
