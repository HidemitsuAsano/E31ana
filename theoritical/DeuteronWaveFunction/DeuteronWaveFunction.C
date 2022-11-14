//dueteron wave function analytics formula taken from
//PHYSICAL REVIEW C, VOLUME 63, 024001
//(C21) (only s-wave)

const double hbarc=0.1973269804; // GeV*fm

//C(5)
const double gammaMass=0.2315380; // fm^(-1)
const double m0=0.9; // fm^(-1)
const double C_swave[]={ 0.88472985, -0.26408759, -0.44114404e-1, -0.14397512e2, 0.85591256e2, -0.31876761e3, 0.70336701e3, -0.90049586e3, 0.66145441e3, -0.25958894e3 };
const double C_dwave[]={                                                };


void DeuteronWaveFunction(){
  double sum_C_swave=0;
  std::vector<double> conv_C;
  for( int i=0; i<sizeof(C_swave)/sizeof(C_swave[0]); i++ ){
    std::cout<<"C_"<<i<<"="<<C_swave[i]<<std::endl;
    conv_C.push_back(C_swave[i]);
    sum_C_swave+=C_swave[i];
  }
  std::cout<<-sum_C_swave<<std::endl;
  conv_C.push_back(-sum_C_swave);

  sum_C_swave=0;
  for( int i=0; i<conv_C.size(); i++ ){
    sum_C_swave+=conv_C[i];
  }
  //check sum should be 0
  std::cout<<sum_C_swave<<std::endl;

  std::vector<double> x, y;
  //  for( double p=0.1; p<1.0; p+=0.1 ){
  for( double p=0; p<1.0; p+=0.001 ){
    double val=0;
    for( int i=0; i<conv_C.size(); i++ ){
      double m=(gammaMass+m0*i)*hbarc;
      val+=conv_C[i]/(p*p+m*m);
      //      std::cout<<i<<"  "<<m<<std::endl;
    }
    x.push_back(p);
    //    std::cout<<p<<"   "<<val<<std::endl;
    y.push_back(val);
  }

  TGraph *gra=new TGraph(x.size(), &x[0], &y[0]);
  gra-> Draw("APL");
}
