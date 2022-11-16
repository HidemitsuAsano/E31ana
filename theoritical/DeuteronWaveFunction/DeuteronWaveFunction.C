//dueteron wave function analytics formula taken from
//PHYSICAL REVIEW C, VOLUME 63, 024001
//(C21) (only s-wave)

//const double hbarc=0.1973269804; // GeV*fm
const double hbarc=197.3269804; // MeV*fm
const double lv = 2.99792458e27; //light velocity
//C(5)
const double gammaMass=0.2315380; // fm^(-1)
const double m0=0.9; // fm^(-1)

//TABLE XX. Coefficients for the parametrized deuteron wave
//functions (n=11).
const double C_swave[]={ 0.88472985e0, -0.26408759e0, -0.44114404e-1, -0.14397512e2, 0.85591256e2, -0.31876761e3, 0.70336701e3, -0.90049586e3, 0.66145441e3, -0.25958894e3 };
const double C_dwave[]={ 0.22623762e-1,-0.50471056e0, 0.56278897e0,-0.16079764e2,0.11126803e3,-0.44667490e3,0.10985907e4,-0.16114995e4};

//S+D wave 
//const double normF = 0.00506526;
const double normF = 0.0048196;
//const double normF = 0.00405278;
//const double normF = 1.0;

double massMn2(int j);
void DeuteronWaveFunction(){
  double sum_C_swave=0.;
  std::vector<double> vec_coeS;
  for( int i=0; i<sizeof(C_swave)/sizeof(C_swave[0]); i++ ){
    std::cout<<"C_S"<<i<<":"<<C_swave[i]<<std::endl;
    vec_coeS.push_back(C_swave[i]);
    sum_C_swave+=C_swave[i];
  }
  std::cout<<-sum_C_swave<<std::endl;
  vec_coeS.push_back(-sum_C_swave);

  sum_C_swave=0;
  for( int i=0; i<vec_coeS.size(); i++ ){
    sum_C_swave+=vec_coeS[i];
  }
  //check sum should be 0
  std::cout<<sum_C_swave<<std::endl;
   
  //calc C_dwave, j=9,10,11
  std::vector<double> vec_coeD;
  
  double sum_C_dwave_uptoj8divmj2=0.0;
  double sum_C_dwave_uptoj8=0.0;
  double sum_C_dwave_uptoj8mulmj2=0.0;
  for(int i=0;i<8;i++){
    std::cout <<"C_D" << i << ":" << C_dwave[i] << std::endl;
    vec_coeD.push_back(C_dwave[i]);
    sum_C_dwave_uptoj8divmj2 += C_dwave[i]/massMn2(i+1);
    sum_C_dwave_uptoj8 += C_dwave[i];
    sum_C_dwave_uptoj8mulmj2 += C_dwave[i]*massMn2(i+1);
    std::cout << i << " "  << massMn2(i+1) << std::endl;
    std::cout << i << " "  << sum_C_dwave_uptoj8divmj2     << std::endl;
  }
  std::cout << "sum_C_dwave_uptoj8divmj2:" << sum_C_dwave_uptoj8divmj2 << std::endl;
  std::cout << "sum_C_dwave_uptoj8:" << sum_C_dwave_uptoj8 << std::endl;
  std::cout << "sum_C_dwave_uptoj8mulmj2:" << sum_C_dwave_uptoj8mulmj2 << std::endl;
   
  //circular permutation of n-2,n-1,n (n=11); 
  double cdj9 = 0.0;
  cdj9 = massMn2(9)/((massMn2(11)-massMn2(9))*(massMn2(10)-massMn2(9)));
  std::cout << "C_D 9" <<  ":" << cdj9 << std::endl; 
  cdj9 = cdj9*(-1.0*massMn2(10)*massMn2(11)*sum_C_dwave_uptoj8divmj2+ (massMn2(10)+massMn2(11))*sum_C_dwave_uptoj8- sum_C_dwave_uptoj8mulmj2   )  ;
  std::cout << "C_D 9" <<  ":" << cdj9 << std::endl; 
  vec_coeD.push_back(cdj9);
  
  double cdj10 = 0.0;
  cdj10 = massMn2(10)/((massMn2(9)-massMn2(10))*(massMn2(11)-massMn2(10)));
  std::cout << "C_D 10" <<  ":" << cdj10 << std::endl; 
  cdj10 = cdj10*(-1.0*massMn2(11)*massMn2(9)*sum_C_dwave_uptoj8divmj2+ (massMn2(11)+massMn2(9))*sum_C_dwave_uptoj8- sum_C_dwave_uptoj8mulmj2   )  ;
  vec_coeD.push_back(cdj10);
  std::cout << "C_D 10" <<  ":" << cdj10 << std::endl; 

  double cdj11 = 0.0;
  cdj11 = massMn2(11)/((massMn2(10)-massMn2(11))*(massMn2(9)-massMn2(11)));
  std::cout << "C_D 11" <<  ":" << cdj11 << std::endl; 
  cdj11 = cdj11*(-1.0*massMn2(9)*massMn2(10)*sum_C_dwave_uptoj8divmj2+ (massMn2(9)+massMn2(10))*sum_C_dwave_uptoj8- sum_C_dwave_uptoj8mulmj2   )  ;
  vec_coeD.push_back(cdj11);
  std::cout << "C_D 11" <<  ":" << cdj11 << std::endl; 

  std::vector<double> x, y,y2p2;
  //  for( double p=0.1; p<1.0; p+=0.1 ){
  //for( double p=0; p<1.0; p+=0.001 ){
  for( double p=0; p<1000; p+=1 ){
    double val=0;
    for( int i=0; i<vec_coeS.size(); i++ ){
      double m2=massMn2(i+1)*hbarc*hbarc;
      val+=vec_coeS[i]/(p*p+m2);
      //      std::cout<<i<<"  "<<m<<std::endl;
    }
    x.push_back(p);
    val = val*sqrt(2./TMath::Pi()/normF);
    //    std::cout<<p<<"   "<<val<<std::endl;
    y.push_back(val);
    y2p2.push_back(p*p*val*val/4./TMath::Pi() );
  }

  std::vector<double> xd,yd,y2p2d;
  for( double p=0; p<1000; p+=1 ){
    double val=0;
    for( int i=0; i<vec_coeD.size(); i++ ){
      double m2=massMn2(i+1)*hbarc*hbarc;
      val+=vec_coeD[i]/(p*p+m2);
      //      std::cout<<i<<"  "<<m<<std::endl;
    }
    xd.push_back(p);
    val = val*sqrt(2./TMath::Pi()/normF);
    //    std::cout<<p<<"   "<<val<<std::endl;
    yd.push_back(val);
    y2p2d.push_back(p*p*val*val/4./TMath::Pi() );
  }

  std::vector<double> xsum,ysum,y2p2sum;
  for( double p=0; p<1000; p+=1 ){
    double valS=0;
    double valD=0;
    for( int i=0; i<vec_coeS.size(); i++ ){
      double m2=massMn2(i+1)*hbarc*hbarc;
      valS+=vec_coeS[i]/(p*p+m2);
      //      std::cout<<i<<"  "<<m<<std::endl;
    }
    for( int i=0; i<vec_coeD.size(); i++ ){
      double m2=massMn2(i+1)*hbarc*hbarc;
      valD+=vec_coeD[i]/(p*p+m2);
      //      std::cout<<i<<"  "<<m<<std::endl;
    }
    xsum.push_back(p);
    valS = valS*sqrt(2./TMath::Pi()/normF);
    valD = valD*sqrt(2./TMath::Pi()/normF);
    //    std::cout<<p<<"   "<<val<<std::endl;
    ysum.push_back(valS+valD);
    y2p2sum.push_back(p*p*(valS*valS+valD*valD)/4./TMath::Pi() );
  }
  

  
  TCanvas *cwave = new TCanvas("cwave","cwave");
  TGraph *grw=new TGraph(x.size(), &x[0], &y[0]);
  grw->SetLineColor(2);
  grw->SetMarkerColor(2);
  grw->GetXaxis()->SetRangeUser(0,700);
  grw->GetXaxis()->SetTitle("Nucleon mom. (CM) [MeV/c]");
  grw->GetXaxis()->CenterTitle();
  grw->GetYaxis()->SetTitle("#Phi(p)");
  grw->GetYaxis()->CenterTitle();
  grw->Draw("APL");

  TCanvas *cwaveq2 = new TCanvas("cwaveq2","cwaveq2");
  TGraph *grw2 = new TGraph(x.size(), &x[0],&y2p2[0]);
  double intw2 = grw2->Integral(0,1000);
  grw2->GetXaxis()->SetRangeUser(0,700);
  grw2->GetXaxis()->SetTitle("Nucleon mom. (CM) [MeV/c]");
  grw2->GetXaxis()->CenterTitle();
  grw2->GetYaxis()->SetTitle("|#Phi(p)|^{2}");
  grw2->GetYaxis()->CenterTitle();
  grw2->SetLineColor(2);
  grw2->SetMarkerColor(2);
  std::cout <<"S-wave norm " <<   intw2*4.0*TMath::Pi() << std::endl;
  grw2->Draw("apl");

  TCanvas *cdwave = new TCanvas("cdwave","cdwave");
  TGraph *grwd=new TGraph(xd.size(), &xd[0], &yd[0]);
  grwd->SetLineColor(2);
  grwd->SetMarkerColor(2);
  grwd->GetXaxis()->SetRangeUser(0,700);
  grwd->GetXaxis()->SetTitle("Nucleon mom. (CM) [MeV/c]");
  grwd->GetXaxis()->CenterTitle();
  grwd->GetYaxis()->SetTitle("#Phi(p)");
  grwd->GetYaxis()->CenterTitle();
  grwd->Draw("APL");
  
  TCanvas *cdwaveq2 = new TCanvas("cdwaveq2","cdwaveq2");
  TGraph *grwd2 = new TGraph(xd.size(), &xd[0],&y2p2d[0]);
  grwd2->GetXaxis()->SetRangeUser(0,700);
  grwd2->GetXaxis()->SetTitle("Nucleon mom. (CM) [MeV/c]");
  grwd2->GetXaxis()->CenterTitle();
  grwd2->GetYaxis()->SetTitle("|#Phi(p)|^{2}");
  grwd2->GetYaxis()->CenterTitle();
  grwd2->SetLineColor(2);
  grwd2->SetMarkerColor(2);
  grwd2->Draw("apl");

  TCanvas *cwavesum = new TCanvas("cwavesum","cdwavesum");
  TGraph *grwsum = new TGraph(xsum.size(), &xsum[0],&ysum[0]);
  grwsum->GetXaxis()->SetRangeUser(0,700);
  grwsum->GetXaxis()->SetTitle("Nucleon mom. (CM) [MeV/c]");
  grwsum->GetXaxis()->CenterTitle();
  grwsum->GetYaxis()->SetTitle("#Phi(p)");
  grwsum->GetYaxis()->CenterTitle();
  grwsum->SetLineColor(2);
  grwsum->SetMarkerColor(2);
  grwsum->Draw("apl");
  
  TCanvas *cwavesumq2 = new TCanvas("cwavesumq2","cdwavesumq2");
  TGraph *grwsum2 = new TGraph(xsum.size(), &xsum[0],&y2p2sum[0]);
  grwsum2->GetXaxis()->SetRangeUser(0,700);
  grwsum2->GetXaxis()->SetTitle("Nucleon mom. (CM) [MeV/c]");
  grwsum2->GetXaxis()->CenterTitle();
  grwsum2->GetYaxis()->SetTitle("|#Phi(p)|^{2}");
  grwsum2->GetYaxis()->CenterTitle();
  grwsum2->SetLineColor(2);
  grwsum2->SetMarkerColor(2);
  grwsum2->Draw("apl");

  double intwsum2 = grwsum2->Integral(0,1000);
  std::cout << "S+D wave norm" <<  intwsum2*4.0*TMath::Pi() << std::endl;

}

double massMn2(int j)
{
  double m=(gammaMass+m0*(j-1));
  return m*m;
}
