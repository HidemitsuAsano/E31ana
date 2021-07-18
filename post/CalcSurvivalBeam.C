void CalcSurvivalBeam()
{

  TFile *_file0 = TFile::Open("evanaIMpisigma_v222.root");
   
  TH1F *EventCheck = (TH1F*)_file0->Get("EventCheck");
  TCanvas *c1 = new TCanvas("c1","c1");
  EventCheck->Draw("HIST");
  //total detected kaon
  double kaonEvent = EventCheck->GetBinContent(2);
  std::cout << "kaon evt. " << kaonEvent << std::endl;
  //T0 multiplicity selection 
  double T0selection = EventCheck->GetBinContent(16);
  std::cout << "T0 OK evt.  " << T0selection << std::endl;
  std::cout << "T0 survival rate  " << T0selection/kaonEvent << std::endl;
  
  //Beam PID
  double beamPID = EventCheck->GetBinContent(5);
  std::cout << "beam PID OK evt " << beamPID << std::endl;
  std::cout << "PID survival rate " << beamPID/T0selection << std::endl;
  
  //BLC single track
  double BLCsingleTrack = EventCheck->GetBinContent(17);
  std::cout << "BLC1,2 killed evt " << BLCsingleTrack << std::endl;
  std::cout << "BLC1,2 survival rate" << (BLCsingleTrack)/beamPID << std::endl;
  
  //
  double nBPCtrack = EventCheck->GetBinContent(18);
  std::cout << "nBPC 1 track evt " << nBPCtrack << std::endl;
  std::cout << "nBPC 1 track survival rate " << nBPCtrack/BLCsingleTrack << std::endl;

  double BPCBLCmatch = EventCheck->GetBinContent(20);
  std::cout << "BLC-BPC matching " << BPCBLCmatch << std::endl;
  std::cout << "BLC-BPC matching survival rate " << BPCBLCmatch/(nBPCtrack) << std::endl;

  double D5AnaOK = EventCheck->GetBinContent(7);
  std::cout << "D5AnaOK " << D5AnaOK << std::endl;
  std::cout << "D5Ana survival rate" << D5AnaOK/(BPCBLCmatch) << std::endl;
  
  std::cout << "Survival before fiducial cut" << D5AnaOK/kaonEvent << std::endl;

  TCanvas *c2 = new TCanvas("c2","c2");
  TH2F* bpcVtx_nofid = (TH2F*)_file0->Get("bpcVtx_nofid");
  TH2F* bpcVtx_fid = (TH2F*)_file0->Get("bpcVtx_fid");
  bpcVtx_nofid->Draw("colz");
  bpcVtx_fid->Draw("boxsame");

  double Nnofid = bpcVtx_nofid->Integral();
  double Nfid = bpcVtx_fid->Integral();
  
  std::cout << "bpc Vtx cut" << Nfid/Nnofid << std::endl;

  std::cout << "Total survival  "<< D5AnaOK*Nfid/Nnofid/kaonEvent << std::endl;

  ifstream ifs("../goodrunlist/goodrun_list");
  std::string str;
  
  int runnum[1024];
 
  int nline=0;
  while(getline(ifs,str)){
    runnum[nline] = atoi(str.c_str());
    std::cout << runnum[nline] << std::endl;
    nline++;
  }

  TH1D* hD5AnaOK = new TH1D("hD5AnaOK","hD5AnaOK",1000,0,1000);
  hD5AnaOK->Sumw2();
  TH1D* hbeam = new TH1D("hbeam","hbeam",1000,0,1000);
  hbeam->Sumw2();
  for(int irun=0;irun<nline;irun++){
    TFile *_file = TFile::Open(Form("/gpfs/group/had/knucl/e15/asano/Run78/IMpisigmav222/evanaIMpisigma_0%03d.root",runnum[irun]),"READ");
    TH1D* hEventCheck = (TH1D*)_file->Get("EventCheck");
    double nKaonEvt = hEventCheck->GetBinContent(2);
    double nD5AnaOK = hEventCheck->GetBinContent(7);
    hbeam->SetBinContent(runnum[irun],nKaonEvt);
    hbeam->SetBinError(runnum[irun],sqrt(nKaonEvt));
    hD5AnaOK->SetBinContent(runnum[irun],nD5AnaOK);
    hD5AnaOK->SetBinError(runnum[irun],sqrt(nD5AnaOK));
  }
  TH1D* hBeamSurvivalR = (TH1D*)hbeam->Clone("hBeamSurvivalR");

  hBeamSurvivalR->Divide(hD5AnaOK,hbeam,1.,1.,"B");
  TCanvas *c3 = new TCanvas("c3","c3");
  hBeamSurvivalR->Scale(Nfid/Nnofid);
  hBeamSurvivalR->SetMarkerStyle(20);
  hBeamSurvivalR->SetXTitle("RUN #");
  hBeamSurvivalR->SetYTitle("Survival Rate");
  hBeamSurvivalR->GetXaxis()->CenterTitle();
  hBeamSurvivalR->GetYaxis()->CenterTitle();
  hBeamSurvivalR->GetYaxis()->SetRangeUser(0,1);
  hBeamSurvivalR->Draw("E");

  TFile *file = new TFile("beamSurvival.root","RECREATE");
  hBeamSurvivalR->Write();
  file->Close();

}
