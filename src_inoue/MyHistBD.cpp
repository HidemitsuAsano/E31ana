#include "MyHistBD.h"

void initHistBD()
{
  const int nseg=27;
  new TH1F("hitpatBD", "hitpatBD", nseg, 0.5, nseg+0.5);
  new TH1F("mulBD", "mulBD", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("TBD_u%d", seg), Form("BD seg%d up TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("TBD_d%d", seg), Form("BD seg%d down TDC",seg), 4000, 0, 4000);  

    new TH1F(Form("ABD_u%d", seg), Form("BD seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("ABD_d%d", seg), Form("BD seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("ABD_u%d_wt", seg), Form("BD seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ABD_d%d_wt", seg), Form("BD seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eBD_u%d", seg), Form("BD seg%d up dE", seg),      550, -10, 100);
    new TH1F(Form("eBD_d%d", seg), Form("BD seg%d down dE", seg),    550, -10, 100);
    new TH1F(Form("eBD_%d",  seg), Form("BD seg%d dE", seg),         550, -10, 100);

    new TH1F(Form("eBD_u%d_wt", seg), Form("BD seg%d up dE", seg),   550, -10, 100);
    new TH1F(Form("eBD_d%d_wt", seg), Form("BD seg%d down dE", seg), 550, -10, 100);
    new TH1F(Form("eBD_%d_wt",  seg), Form("BD seg%d dE", seg),      550, -10, 100);
  }
}

void fillBD_ADC(BeamLineHitMan *blMan)
{
  int nBD=0;
  for( int i=0; i<blMan->nBD(); i++ ){
    HodoscopeLikeHit *hit= blMan->BD(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("TBD_u%d", seg), tu);
    MyHistTools::fillTH(Form("TBD_d%d", seg), td);
    MyHistTools::fillTH(Form("ABD_u%d", seg), au);
    MyHistTools::fillTH(Form("ABD_d%d", seg), ad);

    MyHistTools::fillTH(Form("eBD_u%d", seg), eu);
    MyHistTools::fillTH(Form("eBD_d%d", seg), ed);
    MyHistTools::fillTH(Form("eBD_%d", seg), emean);

    if( hit->CheckRange() ){
      nBD++;
      MyHistTools::fillTH("hitpatBD", seg);
      MyHistTools::fillTH(Form("ABD_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("ABD_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eBD_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eBD_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eBD_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulBD", nBD);
}
