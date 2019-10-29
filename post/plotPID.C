void plotPID(const char* filename=""){
  gStyle->SetPalette(56);
  TFile *infile = new TFile(filename,"READ");
  TCanvas *cPID_CDS = new TCanvas("cPID_CDS","PID_CDS");
  //TH2F* PID_CDS = (TH2F*)infile->Get("PID_CDS");
  //PID_CDS->RebinX(20);
  //PID_CDS->RebinY(20);
  //PID_CDS->Draw("colz");
  
  TH2F* PID_CDS_PIM = (TH2F*)infile->Get("PID_CDS_PIM");
  PID_CDS_PIM->RebinX(20);
  PID_CDS_PIM->RebinY(20);
  PID_CDS_PIM->SetFillColor(2);
  PID_CDS_PIM->SetMarkerColor(2);
  PID_CDS_PIM->SetMarkerSize(1.4);
  PID_CDS_PIM->Draw("scat same");

  TH2F* PID_CDS_PIP = (TH2F*)infile->Get("PID_CDS_PIP");
  PID_CDS_PIP->RebinX(20);
  PID_CDS_PIP->RebinY(20);
  PID_CDS_PIP->SetFillColor(3);
  PID_CDS_PIP->SetMarkerColor(3);
  PID_CDS_PIP->SetMarkerSize(1.4);
  PID_CDS_PIP->Draw("scat same");
  
  TH2F* PID_CDS_Proton = (TH2F*)infile->Get("PID_CDS_Proton");
  PID_CDS_Proton->RebinX(20);
  PID_CDS_Proton->RebinY(20);
  PID_CDS_Proton->SetFillColor(4);
  PID_CDS_Proton->SetMarkerColor(4);
  PID_CDS_Proton->SetMarkerSize(1.4);
  PID_CDS_Proton->Draw("scat same");

  TH2F* PID_CDS_Kaon = (TH2F*)infile->Get("PID_CDS_Kaon");
  PID_CDS_Kaon->RebinX(20);
  PID_CDS_Kaon->RebinY(20);
  PID_CDS_Kaon->SetFillColor(5);
  PID_CDS_Kaon->SetMarkerColor(5);
  PID_CDS_Kaon->SetMarkerSize(1.4);
  PID_CDS_Kaon->Draw("scat same");
}
