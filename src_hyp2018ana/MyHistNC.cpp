#include "MyHistNC.h"

void initHistNC()
{
  const int nseg=112;
  new TH1F("hitpatNC", "hitpatNC", nseg, 0.5, nseg+0.5);
  new TH1F("mulNC", "mulNC", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("ANC_u%d", seg), Form("NC seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("ANC_d%d", seg), Form("NC seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("ANC_u%d_wt", seg), Form("NC seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ANC_d%d_wt", seg), Form("NC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eNC_u%d", seg), Form("NC seg%d up dE", seg),      550, -10, 100);
    new TH1F(Form("eNC_d%d", seg), Form("NC seg%d down dE", seg),    550, -10, 100);
    new TH1F(Form("eNC_%d",  seg), Form("NC seg%d dE", seg),         550, -10, 100);

    new TH1F(Form("eNC_u%d_wt", seg), Form("NC seg%d up dE", seg),   550, -10, 100);
    new TH1F(Form("eNC_d%d_wt", seg), Form("NC seg%d down dE", seg), 550, -10, 100);
    new TH1F(Form("eNC_%d_wt",  seg), Form("NC seg%d dE", seg),      550, -10, 100);
  }
}

void fillNC_ADC(BeamLineHitMan *blMan)
{
  int nNC=0;
  for( int i=0; i<blMan->nNC(); i++ ){
    HodoscopeLikeHit *hit= blMan->NC(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("ANC_u%d", seg), au);
    MyHistTools::fillTH(Form("ANC_d%d", seg), ad);

    MyHistTools::fillTH(Form("eNC_u%d", seg), eu);
    MyHistTools::fillTH(Form("eNC_d%d", seg), ed);
    MyHistTools::fillTH(Form("eNC_%d", seg), emean);

    if( hit->CheckRange() ){
      nNC++;
      MyHistTools::fillTH("hitpatNC", seg);
      MyHistTools::fillTH(Form("ANC_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("ANC_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eNC_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eNC_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eNC_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulNC", nNC);
}
