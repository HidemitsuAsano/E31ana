void DispTKO(const char* filename, int cr=1){
  TFile *f = new TFile(filename,"READ");
  
  TCanvas *c[24] ;
  for(int isl=1;isl<22;isl++){
    c[isl] = new TCanvas(Form("c%d",isl), "TKO entry", 600, 600 );
    c[isl]->Divide(4,8);
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
