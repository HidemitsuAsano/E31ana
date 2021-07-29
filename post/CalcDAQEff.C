void CalcDAQEff()
{
  TFile *_file0 = TFile::Open("evanaIMpisigma_v227.root");
  TH1D* Scaler = (TH1D*)_file0->Get("Scaler");
  TCanvas *c1 = new TCanvas("c1","c1");
  Scaler->Draw("HIST");
  const int ScCH_firstTrig = 4;
  const int ScCH_firstCompVeto = 5;
  const int ScCH_Cosmic = 17;
  double firstTrig = Scaler->GetBinContent(ScCH_firstTrig);
  double firstCompVeto = Scaler->GetBinContent(ScCH_firstCompVeto);
  std::cout << "DAQ acc. " << firstTrig<< std::endl;
  std::cout << "DAQ req. " << firstCompVeto << std::endl;
  std::cout << "DAQ eff. " << (float) (firstCompVeto)/firstTrig << std::endl;
  std::cout << "err    . " << sqrt(firstCompVeto*(firstTrig-firstCompVeto)/firstTrig)/firstTrig << std::endl;

  
  ifstream ifs("../goodrunlist/goodrun_list");
  std::string str;
  
  int runnum[1024];
 
  int iline=0;
  while(getline(ifs,str)){
    runnum[iline] = atoi(str.c_str());
    //std::cout << runnum[iline] << std::endl;
    iline++;
  }
  
  TGraphErrors *gdaqeff = new TGraphErrors();
  TH1D* hacc = new TH1D("hacc","hacc",1000,0,1000);
  hacc->Sumw2();
  TH1D* hreq = new TH1D("hreq","hreq",1000,0,1000);
  hreq->Sumw2();
  for(int irun=0;irun<iline;irun++){
    TFile *_file = TFile::Open(Form("/gpfs/group/had/knucl/e15/asano/Run78/IMpisigmav227/evanaIMpisigma_0%03d.root",runnum[irun]),"READ");
    TH1D* hscaler = (TH1D*)_file->Get("Scaler");
    double firstTrig = hscaler->GetBinContent(ScCH_firstTrig);
    double firstCompVeto = hscaler->GetBinContent(ScCH_firstCompVeto);
    //hreq->SetBinContent(runnum[irun],firstTrig);
    //hreq->SetBinError(runnum[irun],sqrt(firstTrig));
    //hacc->SetBinContent(runnum[irun],firstCompVeto);
    //hacc->SetBinError(runnum[irun],sqrt(firstCompVeto));
    //gdaqeff->SetPoint(gdaqeff->GetN(),runnum[irun],firstTrig/(firstCompVeto+firstTrig)) ;
    TH1D* SCA3 = (TH1D*)_file->Get("SCA3");
    TH1D* SCA4 = (TH1D*)_file->Get("SCA4");
    TH1D* SCA16 = (TH1D*)_file->Get("SCA16");
    double req =0.0;
    double acc =0.0;
    for(int ie=0;ie<SCA3->GetNbinsX();ie++){
      double cosmic = SCA16->GetBinContent(16);
      if(cosmic<1){
        req +=SCA3->GetBinContent(ie);
        acc +=SCA4->GetBinContent(ie);
      }
    }
    hreq->SetBinContent(runnum[irun],req);
    hreq->SetBinError(runnum[irun],sqrt(req));
    hacc->SetBinContent(runnum[irun],acc);
    hacc->SetBinError(runnum[irun],sqrt(acc));
   
    _file->Close();
  }
  TH1D* heff = (TH1D*)hacc->Clone("heff");
  heff->Sumw2();
  heff->Divide(hacc,hreq,1.,1.,"B");

  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  //gdaqeff->SetMarkerStyle(20);
  //gdaqeff->SetMarkerSize(2);
  heff->SetMarkerStyle(20);
  heff->Draw("E");
  //gdaqeff->Draw("AP");
  
  TFile *file = new TFile("daqeff.root","RECREATE");
  heff->Write();
  file->Close();

}
