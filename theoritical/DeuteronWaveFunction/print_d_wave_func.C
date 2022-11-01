const double Const_GeVFM=0.1973269602; // GeV*fm

const double gammaMass=0.2315380; // fm^(-1)
const double m0=0.9; // fm^(-1)
const double C[]={ 0.88472985, -0.26408759, -0.44114404e-1, -0.14397512e2, 0.85591256e2, -0.31876761e3, 0.70336701e3, -0.90049586e3, 0.66145441e3, -0.25958894e3 };

void print_d_wave_func(){
  double sum=0;
  std::vector<double> conv_C;
  for( int i=0; i<sizeof(C)/sizeof(C[0]); i++ ){
    std::cout<<"C_"<<i<<"="<<C[i]<<std::endl;
    conv_C.push_back(C[i]);
    sum+=C[i];
  }
  std::cout<<-sum<<std::endl;
  conv_C.push_back(-sum);

  sum=0;
  for( int i=0; i<conv_C.size(); i++ ){
    sum+=conv_C[i];
  }
  std::cout<<sum<<std::endl;

  std::vector<double> x, y;
  //  for( double p=0.1; p<1.0; p+=0.1 ){
  for( double p=0; p<1.0; p+=0.001 ){
    double val=0;
    for( int i=0; i<conv_C.size(); i++ ){
      double m=(gammaMass+m0*i)*Const_GeVFM;
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
