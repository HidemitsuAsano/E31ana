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
  double reqall = 0.0;
  double accall = 0.0;
  for(int irun=0;irun<iline;irun++){
    std::cout << irun << std::endl;
    TFile *_file = TFile::Open(Form("/gpfs/group/had/knucl/e15/asano/Run78/IMpisigmav227/evanaIMpisigma_0%03d.root",runnum[irun]),"READ");
    //TFile *_file = TFile::Open(Form("~/v207tmp/evanaIMpisigma_0%03d.root",runnum[irun]),"READ");
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
    for(int ie=1;ie<SCA3->GetNbinsX();ie++){
      double cosmic = SCA16->GetBinContent(ie);
      double s3c = SCA3->GetBinContent(ie);
      double s4c = SCA4->GetBinContent(ie);
      if( (cosmic<1) && (s3c>0) && (s4c>0) ){
        req += s3c;
        acc += s4c;
        reqall += s3c;
        accall += s4c;
      }

      //test
      //if(irun==0 && (ie==(SCA3->GetNbinsX()-1)) ){
      //  std::cout << "ie " << ie << std::endl;
      //  std::cout << "reqall " << reqall << std::endl;
      //  std::cout << "accall " << accall << std::endl;
      //  std::cout << "inte   " << SCA3->Integral() << std::endl;
      //}
    }
    if(req>acc){
      hreq->SetBinContent(runnum[irun],req);
      hreq->SetBinError(runnum[irun],sqrt(req));
      hacc->SetBinContent(runnum[irun],acc);
      hacc->SetBinError(runnum[irun],sqrt(acc));
    }
    _file->Close();
  }
  gStyle->SetOptStat(0);
  gPad->SetTopMargin(0.20);
  gStyle->SetTitleY(0.96);
  TH1D* heff = (TH1D*)hacc->Clone("heff");
  heff->Sumw2();
  heff->Divide(hacc,hreq,1.,1.,"B");
  heff->SetTitle("#splitline{DAQ efficiency}{without Cosmic trigger}");
  heff->SetXTitle("run#");
  heff->GetXaxis()->CenterTitle();
  
  //correction
  for(int ix=0;ix<heff->GetNbinsX();ix++){
    double val = heff->GetBinContent(ix);
    if(0<val && val <0.70){
      heff->SetBinContent(ix,0.772601);
    }
  }


  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  //gdaqeff->SetMarkerStyle(20);
  //gdaqeff->SetMarkerSize(2);
  heff->SetMarkerStyle(20);
  heff->Draw("E");
  //gdaqeff->Draw("AP");
  std::cout << "DAQ acc. " << reqall << std::endl;
  std::cout << "DAQ req. " << accall << std::endl;
  std::cout << "DAQ eff. " << (float) accall/reqall << std::endl;
  std::cout << "err    . " << sqrt(accall*(reqall-accall)/reqall)/reqall << std::endl;
  
  TFile *file = new TFile("daqeff.root","RECREATE");
  heff->Write();
  file->Close();

}
