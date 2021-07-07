void CalcDAQEff()
{


  TFile *_file0 = TFile::Open("evanaIMpisigma_v202.root");
  TH1D* Scaler = (TH1D*)_file0->Get("Scaler");
  TCanvas *c1 = new TCanvas("c1","c1");
  Scaler->Draw("HIST");
  const int ScCH_daqacc = 5;
  const int ScCH_daqrej = 6;
  double daqacc = Scaler->GetBinContent(ScCH_daqacc);
  double daqrej = Scaler->GetBinContent(ScCH_daqrej);
  std::cout << "DAQ acc. " << daqacc<< std::endl;
  std::cout << "DAQ rej. " << daqrej << std::endl;
  std::cout << "DAQ eff. " << daqacc/(daqrej+daqacc) << std::endl;


  
  ifstream ifs("../goodrunlist/goodrun_list");
  std::string str;
  
  int runnum[1024];
 
  int iline=0;
  while(getline(ifs,str)){
    runnum[iline] = atoi(str.c_str());
    std::cout << runnum[iline] << std::endl;


    iline++;
  }
  
  TGraphErrors *gdaqeff = new TGraphErrors();
  TH1D* hacc = new TH1D("hacc","hacc",1000,0,1000);
  hacc->Sumw2();
  TH1D* hreq = new TH1D("hreq","hreq",1000,0,1000);
  hreq->Sumw2();
  for(int irun=0;irun<iline;irun++){
    TFile *_file = TFile::Open(Form("/gpfs/group/had/knucl/e15/asano/Run78/IMpisigmav207/evanaIMpisigma_0%03d.root",runnum[irun]),"READ");
    TH1D* hscaler = (TH1D*)_file->Get("Scaler");
    double daqacc = hscaler->GetBinContent(ScCH_daqacc);
    double daqrej = hscaler->GetBinContent(ScCH_daqrej);
    hacc->SetBinContent(runnum[irun],daqacc);
    hacc->SetBinError(runnum[irun],sqrt(daqacc));
    hreq->SetBinContent(runnum[irun],daqacc+daqrej);
    hreq->SetBinError(runnum[irun],sqrt(daqacc+daqrej));
    //gdaqeff->SetPoint(gdaqeff->GetN(),runnum[irun],daqacc/(daqrej+daqacc)) ;

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

}
