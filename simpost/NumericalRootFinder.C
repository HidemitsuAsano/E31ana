//----------------------------------------//
// from Iwasaki-san slide on 2019, Jan.30 //
//----------------------------------------//
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>

#include <TDatabasePDG.h>

#if 1
const double d_mass  = 1.87561;
const double K_mass    = 0.493677;
const double K0_mass   = 0.497611;
const double n_mass    = 0.939565;
const double K0n_mass  = K0_mass + n_mass;
const double piSp_mass = 1.18937+0.13957;
const double piSm_mass = 1.197449+0.13957;
#endif

//const double pK = 1.05; //GeV/c = maximum ~ 1.018*1.025
const double pK = 1.0; //default value of simulation
const double EK = sqrt(K_mass*K_mass+pK*pK);

double func_EM(double *x, double *par);
double func_cq(double *x, double *par);
double func_thetalab(double q, double M);
double func_C(double *x, double *par);
double func_S2(double *x, double *par);
double func(double *x, double *par);
double funclab(double *x, double *par);


int NumericalRootFinder()
{
  double target_mass = d_mass;
  TVector3 target_mom(0.0, 0.0, 0.0);
  TLorentzVector target(target_mom, sqrt(target_mom.Mag2()+target_mass*target_mass));
  
  double beam_mass = K_mass;
  TVector3 beam_mom(0.0, 0.0, pK);
  TLorentzVector beam(beam_mom, EK);
  
  TLorentzVector sys = beam+target;

  const double COS_MIN = -1;
  const double COS_MAX = 1;
  const int    COS_BIN = 200;

  //const double M_MIN = piSm_mass;
  const double M_MIN = piSp_mass;
  //const double M_MIN = K0n_mass;
  const double M_MAX = sys.M()-n_mass+0.00001;
  const int    M_BIN = 100;
  cerr<<M_MIN<<" "<<M_MAX<<endl;
  
  double value[2][COS_BIN+1][M_BIN]; //[m,q][][]
  double value_th[2][2]; //[m,q][]
  
  // parameter = m, cos
  // variable  = q_Kn = x
  
  for( int i=0; i<COS_BIN+1; i++ ){
    double cos = COS_MIN+i*(COS_MAX-COS_MIN)/COS_BIN; // cos(theta_n^CM)

    for( int j=0; j<M_BIN; j++ ){
      double m = M_MIN+j*(M_MAX-M_MIN)/M_BIN; // IM(piSp)

      // Create the function and wrap it
      TF1 f( "func", func, 0, 2, 2 );
      //TF1 f( "funclab", funclab, 0, 2, 2 );
      f.SetParameter( 0, m ); // m
      f.SetParameter( 1, cos ); // cos
      ROOT::Math::WrappedTF1 wf1( f );
      
      // Create the Integrator
      //ROOT::Math::BrentRootFinder brf; // does not work?
      ROOT::Math::Roots::Brent brf;
      
      // Set parameters of the method
      brf.SetFunction( wf1, 0.0000001, 2 );
      brf.Solve();

      value[0][i][j] = m;
      value[1][i][j] = brf.Root(); // q
      //cout << i << " " << j << " " << m << " " << cos << " " << brf.Root() << endl;
      
    } // for( int j=0; j<M_BIN; j++ ){
  } // for( int i=0; i<COS_BIN+1; i++ ){

  value_th[0][0] = M_MIN;
  value_th[0][1] = M_MIN;
  value_th[1][0] = value[1][0][0];
  value_th[1][1] = value[1][COS_BIN][0];

  char com[64];
  TGraph *gr[COS_BIN+1];
  TMultiGraph *mg = new TMultiGraph("mg","kinematic limits");
  for( int i=0; i<COS_BIN+1; i++ ){
    gr[i] = new TGraph(M_BIN, value[0][i], value[1][i]);
    sprintf(com, "gr_%d", i);
    gr[i]->SetName(com);
    gr[i]->SetLineColor(1);
    gr[i]->SetLineWidth(3);
    gr[i]->SetLineStyle(2);
    mg->Add(gr[i]);
  }

  TGraph *gr_th;
  gr_th = new TGraph(2, value_th[0], value_th[1]);
  gr_th->SetName("th");
  gr_th->SetLineColor(1);
  gr_th->SetLineWidth(4);
  gr_th->SetLineStyle(10);
  mg->Add(gr_th);
  TCanvas *c1;
  c1 = new TCanvas("c1", "", 600, 600);
  TH2F *his = new TH2F("his","his",100,1,2,120,-0.2,2);
  his->Draw();
  for( int i=0; i<COS_BIN+1; i++ ){
    gr[i]->Draw("same");
  }
  gr_th->Draw("same");

  TFile *out = new TFile("NumericalRootFinder_fine200_Sm1GeV.root", "recreate");
  for( int i=0; i<COS_BIN+1; i++ ){
    gr[i]->Write();
  }
  gr_th->Write();
  mg->Write();
  out->Close();

  return 0;
}

double func_EM(double *x, double *par)
// x      = q
// par[0] = m
{
  double q = x[0];
  double f = TMath::Sqrt(par[0]*par[0]+q*q);
  return f;
}

double func_cq(double *x, double *par)
// x      = q
// par[0] = m
{
  double q = x[0];
  double f = (n_mass*n_mass+pK*pK+q*q-TMath::Power(EK+d_mass-func_EM(x,par),2))/(2*q*pK);
  return f;
}

double func_thetalab(double q, double M)
// x      = q
// par[0] = m
{ 
  double f = (n_mass*n_mass+pK*pK+q*q-TMath::Power(EK+d_mass-TMath::Sqrt(M*M+q*q),2.))/(2.*q*pK);
  double tannlab = -q*(sqrt(1.0-f*f))/(pK-q*f);
  return -1.0*atan(tannlab)*180./3.1415926535;
}

double func_C(double *x, double *par)
// x      = q
// par[0] = m
{
  double q = x[0];
  double f = (EK+d_mass)*(pK-q*func_cq(x,par))-pK*(EK+d_mass-func_EM(x,par));
  return f;
}

double func_S2(double *x, double *par)
// x      = q
// par[0] = m
{
  double q = x[0];
  double f = q*q*(1-func_cq(x,par)*func_cq(x,par))*(d_mass*d_mass+2*d_mass*EK+K_mass*K_mass);
  return f;
}

double func(double *x, double *par)
// x      = q
// par[0] = m
// par[1] = cos
{
  double q = x[0];
  double f = par[1]-func_C(x,par)/TMath::Sqrt(func_C(x,par)*func_C(x,par)+func_S2(x,par));
  //cerr<<q<<" "<<f<<endl;
  return f;
}

double functheta(double q, double M)
{
  double x[1];
  x[0]=q;
  double par[1];
  par[0]=M;
  double costhetaCM = func_C(x,par)/TMath::Sqrt(func_C(x,par)*func_C(x,par)+func_S2(x,par));
  return acos(costhetaCM)*180./3.1415926535;
}

double funcos(double q, double M)
{
  double x[1];
  x[0]=q;
  double par[1];
  par[0]=M;
  double costhetaCM = func_C(x,par)/TMath::Sqrt(func_C(x,par)*func_C(x,par)+func_S2(x,par));
  return costhetaCM;

}

