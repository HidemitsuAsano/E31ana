#include "MyHistT0.h"

void initHistT0()
{
  const int nseg=5;
  new TH1F("hitpatT0", "hitpatT0", nseg, 0.5, nseg+0.5);
  new TH1F("mulT0", "mulT0", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("AT0_u%d", seg), Form("T0 seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("AT0_d%d", seg), Form("T0 seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("AT0_u%d_wt", seg), Form("T0 seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("AT0_d%d_wt", seg), Form("T0 seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eT0_u%d", seg), Form("T0 seg%d up dE", seg), 50, -0.5, 4.5);
    new TH1F(Form("eT0_d%d", seg), Form("T0 seg%d down dE", seg), 50, -0.5, 4.5);
    new TH1F(Form("eT0_%d",  seg), Form("T0 seg%d dE", seg), 50, -0.5, 4.5);

    new TH1F(Form("eT0_u%d_wt", seg), Form("T0 seg%d up dE", seg), 50, -0.5, 4.5);
    new TH1F(Form("eT0_d%d_wt", seg), Form("T0 seg%d down dE", seg), 50, -0.5, 4.5);
    new TH1F(Form("eT0_%d_wt",  seg), Form("T0 seg%d dE", seg), 50, -0.5, 4.5);
  }
}

void fillT0_ADC(BeamLineHitMan *blMan)
{
  int nT0=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    HodoscopeLikeHit *hit= blMan->T0(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("AT0_u%d", seg), au);
    MyHistTools::fillTH(Form("AT0_d%d", seg), ad);

    MyHistTools::fillTH(Form("eT0_u%d", seg), eu);
    MyHistTools::fillTH(Form("eT0_d%d", seg), ed);
    MyHistTools::fillTH(Form("eT0_%d", seg), emean);

    if( hit->CheckRange() ){
      nT0++;
      MyHistTools::fillTH("hitpatT0", seg);
      MyHistTools::fillTH(Form("AT0_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("AT0_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eT0_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eT0_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eT0_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulT0", nT0);
}
