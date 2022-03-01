void calc_purity_H2()
{
  const int version = 9;
  const int sysud=0;
  int dEcut[3]={2,4,6};
  TFile *f[3];
  TFile *fr[3];
  gStyle->SetOptStat(0);

  for(int iEcut=0;iEcut<3;iEcut++){
    f[iEcut] = TFile::Open(Form("evanaIMsigma_npi_h2_v%d_out_dE%d_iso_nostop_sub_sys%d.root",version,dEcut[iEcut],sysud));
    fr[iEcut] = TFile::Open(Form("evanaIMsigma_npi_h2_v%d_out_dE%d_iso_nostop.root",version,dEcut[iEcut]));
  }

  TCanvas *ccos_purity_Sp = new TCanvas("ccos_purity_Sp","ccos_purity_Sp");

  TH1F* Cospicm_pi_Sp_S[3];
  TH1F* Cospicm_pi_Sp_SN[3];
  TH1F* Cospicm_purity_Sp[3];
  for(int ie=0;ie<3;ie++){
    Cospicm_pi_Sp_S[ie]  = (TH1F*)f[ie]->Get("Cospicm_pi_Sp");
    Cospicm_pi_Sp_S[ie]->RebinX(5);
    Cospicm_pi_Sp_SN[ie] = (TH1F*)fr[ie]->Get("Cospicm_pi_Sp");
    Cospicm_pi_Sp_SN[ie]->RebinX(5);
    Cospicm_purity_Sp[ie] = (TH1F*)Cospicm_pi_Sp_S[ie]->Clone(Form("Cospicm_putiry%d",ie));
    Cospicm_purity_Sp[ie]->SetTitle("purity S/(S+N) #Sigma^{+}#pi^{-}");
    Cospicm_purity_Sp[ie]->Divide(Cospicm_pi_Sp_S[ie],Cospicm_pi_Sp_SN[ie],1.,1.,"b");
  }
  Cospicm_purity_Sp[0]->SetMarkerStyle(20);
  Cospicm_purity_Sp[0]->SetMaximum(1);
  Cospicm_purity_Sp[0]->SetMinimum(0);
  Cospicm_purity_Sp[0]->Draw();
  Cospicm_purity_Sp[1]->SetLineColor(2);
  Cospicm_purity_Sp[1]->SetMarkerStyle(21);
  Cospicm_purity_Sp[1]->SetMarkerColor(2);
  Cospicm_purity_Sp[1]->Draw("same");
  Cospicm_purity_Sp[2]->SetMarkerStyle(22);
  Cospicm_purity_Sp[2]->SetMarkerColor(3);
  Cospicm_purity_Sp[2]->SetLineColor(3);
  Cospicm_purity_Sp[2]->Draw("same");


  TCanvas *ccos_purity_Sm = new TCanvas("ccos_purity_Sm","ccos_purity_Sm");

  TH1F* Cospicm_pi_Sm_S[3];
  TH1F* Cospicm_pi_Sm_SN[3];
  TH1F* Cospicm_purity_Sm[3];
  for(int ie=0;ie<3;ie++){
    Cospicm_pi_Sm_S[ie]  = (TH1F*)f[ie]->Get("Cospicm_pi_Sm");
    Cospicm_pi_Sm_S[ie]->RebinX(5);
    Cospicm_pi_Sm_SN[ie] = (TH1F*)fr[ie]->Get("Cospicm_pi_Sm");
    Cospicm_pi_Sm_SN[ie]->RebinX(5);
    Cospicm_purity_Sm[ie] = (TH1F*)Cospicm_pi_Sm_S[ie]->Clone(Form("Cospicm_putiry%d",ie));
    Cospicm_purity_Sm[ie]->SetTitle("purity S/(S+N) #Sigma^{-}#pi^{+}");
    Cospicm_purity_Sm[ie]->Divide(Cospicm_pi_Sm_S[ie],Cospicm_pi_Sm_SN[ie],1.,1.,"b");
  }
  Cospicm_purity_Sm[0]->SetMarkerStyle(20);
  Cospicm_purity_Sm[0]->SetMaximum(1);
  Cospicm_purity_Sm[0]->SetMinimum(0);
  Cospicm_purity_Sm[0]->Draw();
  Cospicm_purity_Sm[1]->SetLineColor(2);
  Cospicm_purity_Sm[1]->SetMarkerStyle(21);
  Cospicm_purity_Sm[1]->SetMarkerColor(2);
  Cospicm_purity_Sm[1]->Draw("same");
  Cospicm_purity_Sm[2]->SetMarkerStyle(22);
  Cospicm_purity_Sm[2]->SetMarkerColor(3);
  Cospicm_purity_Sm[2]->SetLineColor(3);
  Cospicm_purity_Sm[2]->Draw("same");







}
