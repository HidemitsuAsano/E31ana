void CalcSurvivalBeam()
{

  //TFile *_file0 = TFile::Open("evanaIMpisigma_v209.root");
  TFile *_file0 = TFile::Open("evanaIMpisigma_v211.root");
   
  TH1F *EventCheck = (TH1F*)_file0->Get("EventCheck");
  TCanvas *c1 = new TCanvas("c1","c1");
  EventCheck->Draw("HIST");
  //total detected kaon
  double kaonEvent = EventCheck->GetBinContent(2);
  std::cout << "kaon evt. " << kaonEvent << std::endl;
  //T0 multiplicity selection 
  double T0selection = EventCheck->GetBinContent(16);
  std::cout << "T0 killed evt.  " << T0selection << std::endl;
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
  
  TCanvas *c2 = new TCanvas("c2","c2");
  TH2F* bpcVtx_nofid = (TH2F*)_file->Get("bpcVtx_nofid");
  TH2F* bpcVtx_fid = (TH2F*)_file->Get("bpcVtx_fid");
  bpcVtx_nofid->Draw("colz");
  bpcVtx_fid->Draw("boxsame");

  double Nnofid = bpcVtx_nofid->Integral();
  double Nfid = bpcVtx_fid->Integral();
  
  std::cout << "bpc Vtx cut" << Nnofid/Nfid << std::endl;

  std::cout << "Total survival  "<< D5AnaOK*Nnofid/Nfid/kaonEvent << std::endl;
}
