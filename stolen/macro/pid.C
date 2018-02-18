#include "Unit.C"

const double mmu = 106.5*MeV;
//const double mmu = 940*MeV;
const double toflength = 1.*m;
const double tofresl = 0.1*nsec;

pid()
{

  double P[10000],Beta[10000],TOFwR[10000], BetawR[10000], BetaInv[10000], PwR2[10000];
  TH1F *h_pwr2 = new TH1F("h_pwr2","",500,0,1);
  for(int i=0; i<10000; i++ ){
    double pmu = gRandom->Uniform()*GeV;//i*MeV;
    double beta = pmu/sqrt(pmu*pmu+mmu*mmu);
    double tof = toflength/(beta*LightVelocity);
    double tof_wresl = tof + tofresl*gRandom->Gaus();
    double beta_wresl = (toflength/tof_wresl)/LightVelocity;
    double pmu_wresl2 = mmu*mmu*beta_wresl*beta_wresl/(1-beta_wresl*beta_wresl);
    std::cout << " pmu:" << pmu
	      << " beta:" << beta
	      << " tof:" << tof
	      << " tof_wresl:" << tof_wresl
	      << " beta_wresl:" << beta_wresl
	      << " pmu_wresl2:" << pmu_wresl2
	      << std::endl;

    P[i] = pmu/GeV;
    Beta[i] = beta;
    TOFwR[i] = tof_wresl/nsec;
    BetawR[i] = beta_wresl;
    BetaInv[i] = 1./beta;
    h_pwr2->Fill(pmu_wresl2/GeV/GeV);
  }

  TGraph *gra0 = new TGraph( 10000, BetaInv, P );
  TGraph *gra1 = new TGraph( 10000, Beta, P );
  TGraph *gra2 = new TGraph( 10000, BetawR, P );
  
  TCanvas *c0 = new TCanvas( "c0", "c0", 600, 600 );
  gra0->Draw("apw");
  TCanvas *c1 = new TCanvas( "c1", "c1", 600, 600 );
  gra1->Draw("apw");
  gra2->SetMarkerColor(2);
  gra2->Draw("p");

  TCanvas *c2 = new TCanvas( "c2", "c2", 600, 600 );
  h_pwr2->Draw();



}
