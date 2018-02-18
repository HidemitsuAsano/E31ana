void DispTKO(const char* filename){
  TFile *f = new TFile(filename,"READ");
  
  const int nTKOcrate=9;
  const int nTKOslot;
  TCanvas *c[nTKOcrate];
  for(int icr=0;icr < nTKOcrate; icr++){
    c[icr] = new TCanvas(Form("c%d",icr), "TKO entry", 600, 600 );
    c[icr]->Divide(4,8);
    for(int ich=0;ich<=31; ich++){
      c[isl]->cd(ich+1);
      TH1F* htko = (TH1F*)f->Get( Form( "TKOc%dn%da%d", cr,isl,ich ) ); 
      if(!htko){
        std::cout << " !!! " << cr <<" "<< sl <<" "<< ich <<std::endl;
        continue;
      }
      htko->Draw();
    }
    char fname[256];
    sprintf(fname,"cr%d.png",isl);
    c[isl]->SaveAs(fname,"PNG");
  }
}
