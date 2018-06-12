#include "MyHistPC.h"

void initHistPC()
{
  const int nseg=27;
  new TH1F("hitpatPC", "hitpatPC", nseg, 0.5, nseg+0.5);
  new TH1F("mulPC", "mulPC", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("APC_u%d", seg), Form("PC seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("APC_d%d", seg), Form("PC seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("APC_u%d_wt", seg), Form("PC seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("APC_d%d_wt", seg), Form("PC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("ePC_u%d", seg), Form("PC seg%d up dE", seg),      550, -10, 100);
    new TH1F(Form("ePC_d%d", seg), Form("PC seg%d down dE", seg),    550, -10, 100);
    new TH1F(Form("ePC_%d",  seg), Form("PC seg%d dE", seg),         550, -10, 100);

    new TH1F(Form("ePC_u%d_wt", seg), Form("PC seg%d up dE", seg),   550, -10, 100);
    new TH1F(Form("ePC_d%d_wt", seg), Form("PC seg%d down dE", seg), 550, -10, 100);
    new TH1F(Form("ePC_%d_wt",  seg), Form("PC seg%d dE", seg),      550, -10, 100);
  }
}

void fillPC_ADC(BeamLineHitMan *blMan)
{
  int nPC=0;
  for( int i=0; i<blMan->nPC(); i++ ){
    HodoscopeLikeHit *hit= blMan->PC(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("APC_u%d", seg), au);
    MyHistTools::fillTH(Form("APC_d%d", seg), ad);

    MyHistTools::fillTH(Form("ePC_u%d", seg), eu);
    MyHistTools::fillTH(Form("ePC_d%d", seg), ed);
    MyHistTools::fillTH(Form("ePC_%d", seg), emean);

    if( hit->CheckRange() ){
      nPC++;
      MyHistTools::fillTH("hitpatPC", seg);
      MyHistTools::fillTH(Form("APC_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("APC_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("ePC_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("ePC_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("ePC_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulPC", nPC);
}
