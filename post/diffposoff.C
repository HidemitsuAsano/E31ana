void diffposoff(){

  TFile *_file0 = TFile::Open("evanaIMpisigma_v75.root");
  //TFile *_file0 = TFile::Open("evanaIMpisigma_v68.root");
  _file0->cd();
  char hname[256];
  TH1F *h[37];
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(6,6);
  _file0->Print();
  for(int iseg=1;iseg<=36;iseg++){
    sprintf(hname,"CDH_diffpos_pi_z_seg%d",iseg);
    c1->cd(iseg);
    h[iseg] = (TH1F*)_file0->Get(hname);
    h[iseg]->Draw();
    //std::cout << iseg << "  " << h[iseg]->GetMean() << " " << h[iseg]->GetRMS() << std::endl;
    char word[256];
    sprintf(word,"if(Seg==%d) HitPosition -= %f", iseg, h[iseg]->GetMean());
    cout << word << endl;
  }







}
