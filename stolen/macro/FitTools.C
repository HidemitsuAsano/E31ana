test(int aaa){
  std::cout<<"test input: "<<aaa<<std::endl;
}
bool SearchPedPeak(TH1F* h1,double par[4],double width=5){
  //  return false;
  if(!h1){ std::cout<<h1->GetName()<<" !!!"<<std::endl; return false; } // just check
  if(h1->GetEntries()<100){
    std::cout<<h1->GetName()<<" !!!"<<std::endl; 
    return false;
  }
  int max=h1->GetBinCenter(h1->GetMaximumBin());
  int fitStatus=h1->Fit("gaus","0q","L",max-width,max+width);
  //  h1->GetFunction("gaus")->SetLineStyle(2);
  //  h1->GetFunction("gaus")->Draw("same");
  if(fitStatus){
    std::cout<<h1->GetName()<<" !!!"<<std::endl; 
    return false;
  }
  if(0){
    double offs=h1->GetFunction("gaus")->GetParameter(1);
    double sigma=h1->GetFunction("gaus")->GetParameter(2);
    int fitStatus=h1->Fit("gaus","0q","",offs-5*fabs(sigma),offs+5*fabs(sigma));
    std::cout<<max<<"  "<<offs<<"  "<<sigma<<std::endl;
    if(fitStatus) return false;
    h1->GetFunction("gaus")->Draw("same");
  }
  h1->GetFunction("gaus")->GetParameters(par);
  par[3]=h1->GetFunction("gaus")->GetChisquare();
  return true;
}

bool SearchMIPPeak(TH1F* h1,double par[4],double width){
  //  return false;
  if(h1==0){ std::cerr<<h1->GetName()<<std::endl; return false; } // just check
  if(h1->GetEntries()<100) return false;
  int max=h1->GetBinCenter(h1->GetMaximumBin());
  //  int max=h1->GetMean();
  int fitStatus=h1->Fit("gaus","0q","L",max-width,max+width);
  if(fitStatus){
    std::cout<<h1->GetName()<<"  gaus1 failed"<<std::endl;
    return false;
  }
  double offs=h1->GetFunction("gaus")->GetParameter(1);
  double sigma=h1->GetFunction("gaus")->GetParameter(2);
  //  std::cout<<max<<"  "<<offs<<"  "<<sigma<<std::endl;
  int fitStatus=h1->Fit("gaus","0q","L",offs-3*sigma,offs+3*sigma);
  if(fitStatus){
    std::cout<<h1->GetName()<<"  gaus2 failed"<<std::endl;
    return false;
  }
  offs=h1->GetFunction("gaus")->GetParameter(1);
  sigma=h1->GetFunction("gaus")->GetParameter(2);
  //  std::cout<<max<<"  "<<offs<<"  "<<sigma<<std::endl;
  int fitStatus=h1->Fit("landau","0q","L",offs-1*sigma,offs+2*sigma);
  if(fitStatus){
    std::cout<<h1->GetName()<<"  landau failed"<<std::endl;
    return false;
  }
  TF1 *f1=h1->GetFunction("landau");
  f1->SetLineColor(2);
  f1->Draw("same");
  f1->GetParameters(par);
  par[3]=f1->GetChisquare();
  return true;
}
