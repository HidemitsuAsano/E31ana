#include "MyHistDEF.h"

void initHistDEF()
{
  const int nseg=8;
  new TH1F("hitpatDEF", "hitpatDEF", nseg, 0.5, nseg+0.5);
  new TH1F("mulDEF", "mulDEF", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("ADEF_u%d", seg), Form("DEF seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("ADEF_d%d", seg), Form("DEF seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("ADEF_u%d_wt", seg), Form("DEF seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ADEF_d%d_wt", seg), Form("DEF seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eDEF_u%d", seg), Form("DEF seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eDEF_d%d", seg), Form("DEF seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eDEF_%d",  seg), Form("DEF seg%d dE", seg), 500, -0.5, 4.5);

    new TH1F(Form("eDEF_u%d_wt", seg), Form("DEF seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eDEF_d%d_wt", seg), Form("DEF seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eDEF_%d_wt",  seg), Form("DEF seg%d dE", seg), 500, -0.5, 4.5);
  }
}

void fillDEF_ADC(BeamLineHitMan *blMan)
{
  int nDEF=0;
  for( int i=0; i<blMan->nDEF(); i++ ){
    HodoscopeLikeHit *hit= blMan->DEF(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("ADEF_u%d", seg), au);
    MyHistTools::fillTH(Form("ADEF_d%d", seg), ad);

    MyHistTools::fillTH(Form("eDEF_u%d", seg), eu);
    MyHistTools::fillTH(Form("eDEF_d%d", seg), ed);
    MyHistTools::fillTH(Form("eDEF_%d", seg), emean);

    if( hit->CheckRange() ){
      nDEF++;
      MyHistTools::fillTH("hitpatDEF", seg);
      MyHistTools::fillTH(Form("ADEF_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("ADEF_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eDEF_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eDEF_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eDEF_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulDEF", nDEF);
}
