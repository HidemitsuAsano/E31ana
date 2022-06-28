void CountKaonv2()
{
  //v225, unbiased kaon selected file 
  TFile *_file0 = TFile::Open("evanaIMpisigma_v244.root");
  TH1D* Scaler = (TH1D*)_file0->Get("Scaler");
  TCanvas *c1 = new TCanvas("c1","c1");
  Scaler->Draw("HIST");
  const int ScCH_K = Scaler->FindBin(10);
  double nK = Scaler->GetBinContent(ScCH_K);
  std::cout << "kaon total " << (float) nK << std::endl;


  ifstream ifs("../goodrunlist/goodrun_list");
  std::string str;
  
  int runnum[1024];
 
  int iline=0;
  while(getline(ifs,str)){
    runnum[iline] = atoi(str.c_str());
   // std::cout << runnum[iline] << std::endl;
    iline++;
  }

  TH1D* hnK = new TH1D("hnK","Scaler Kaon #",1000,0,1000);
  for(int irun=0;irun<iline;irun++){
    //TFile *_file = TFile::Open(Form("/gpfs/group/had/knucl/e15/asano/Run78/IMpisigmav225/evanaIMpisigma_0%03d.root",runnum[irun]),"READ");
    TFile *_file = TFile::Open(Form("/gpfs/group/had/knucl/e15/asano/Run78/IMpisigmav244/evanaIMpisigma_0%03d.root",runnum[irun]),"READ");
    ///TH1D* hscaler = (TH1D*)_file->Get("Scaler");
    TH1D* SCA11 = (TH1D*)_file->Get("SCA11");
    //double nK = hscaler->GetBinContent(ScCH_K);
    int nbinsX = SCA11 ->GetNbinsX();
    nK=0;
    for(int iev=0;iev<nbinsX-1;iev++){
      int val = SCA11->GetBinContent(iev);
      if(0< val && val<1e6) nK += val; 
    }
    
    if(runnum[irun]>=233) nK*=10;
    //if(runnum[irun]>=217) nK*=10;
    hnK->SetBinContent(runnum[irun],nK);
    hnK->SetBinError(runnum[irun],sqrt(nK));

    _file->Close();
  }

  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  //gdaqeff->SetMarkerStyle(20);
  //gdaqeff->SetMarkerSize(2);
  gStyle->SetOptStat("i");
  hnK->SetMarkerStyle(20);
  hnK->SetXTitle("RUN#");
  hnK->GetXaxis()->CenterTitle();
  hnK->Draw("E");
  std::cout << hnK->Integral() << std::endl;

  TFile *file = new TFile("lumi.root","RECREATE");
  hnK->Write();
  file->Close();

}
