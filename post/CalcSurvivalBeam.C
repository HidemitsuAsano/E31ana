void CalcSurvivalBeam()
{

  TFile *_file0 = TFile::Open("evanaIMpisigma_v209.root");
   
  TH1F *EventCheck = (TH1F*)_file0->Get("EventCheck");
  TCanvas *c1 = new TCanvas("c1","c1");
  EventCheck->Draw("HIST");
  //total detected kaon
  double kaonEvent = EventCheck->GetBinContent(2);
  std::cout << "kaon evt. " << kaonEvent << std::endl;
  //T0 multiplicity selection 
  double T0selection = EventCheck->GetBinContent(16);
  std::cout << "T0 kill evt.  " << T0selection << std::endl;
  std::cout << "T0 survival rate  " << 1.0 - T0selection/kaonEvent << std::endl;
  
  //Beam PID
  double beamPID = EventCheck->GetBinContent
  
}
