void Plot_Kine800and1000(){
  

  TFile *file800 = TFile::Open("NumericalRootFinder_fine20_Sp800MeV.root");
  TMultiGraph *mg800 = (TMultiGraph*)file800->Get("mg");  
  mg800->GetYaxis()->SetRangeUser(0,0.65);
  mg800->GetXaxis()->SetTitle("IM(#pi#Sigma) [GeV/c^2]");
  mg800->GetXaxis()->CenterTitle();
  mg800->GetYaxis()->SetTitle("mom. transfer [GeV/c]");
  mg800->GetYaxis()->CenterTitle();
  mg800->Draw("ac");
  TFile *file1 = TFile::Open("NumericalRootFinder_fine20_Sp1GeV.root");
  TMultiGraph *mg1 = (TMultiGraph*)file1->Get("mg");  
  mg1->Draw("c");








}
