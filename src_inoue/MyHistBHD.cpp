#include "MyHistBHD.h"

void initHistBHD()
{
  const int nseg=20;
  new TH1F("hitpatBHD", "hitpatBHD", nseg, 0.5, nseg+0.5);
  new TH1F("mulBHD", "mulBHD", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("ABHD_u%d", seg), Form("BHD seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("ABHD_d%d", seg), Form("BHD seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("ABHD_u%d_wt", seg), Form("BHD seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ABHD_d%d_wt", seg), Form("BHD seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eBHD_u%d", seg), Form("BHD seg%d up dE", seg), 50, -0.5, 4.5);
    new TH1F(Form("eBHD_d%d", seg), Form("BHD seg%d down dE", seg), 50, -0.5, 4.5);
    new TH1F(Form("eBHD_%d",  seg), Form("BHD seg%d dE", seg), 50, -0.5, 4.5);

    new TH1F(Form("eBHD_u%d_wt", seg), Form("BHD seg%d up dE", seg), 50, -0.5, 4.5);
    new TH1F(Form("eBHD_d%d_wt", seg), Form("BHD seg%d down dE", seg), 50, -0.5, 4.5);
    new TH1F(Form("eBHD_%d_wt",  seg), Form("BHD seg%d dE", seg), 50, -0.5, 4.5);
  }
}

void fillBHD_ADC(BeamLineHitMan *blMan)
{
  int nBHD=0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    HodoscopeLikeHit *hit= blMan->BHD(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("ABHD_u%d", seg), au);
    MyHistTools::fillTH(Form("ABHD_d%d", seg), ad);

    MyHistTools::fillTH(Form("eBHD_u%d", seg), eu);
    MyHistTools::fillTH(Form("eBHD_d%d", seg), ed);
    MyHistTools::fillTH(Form("eBHD_%d", seg), emean);

    if( hit->CheckRange() ){
      nBHD++;
      MyHistTools::fillTH("hitpatBHD", seg);
      MyHistTools::fillTH(Form("ABHD_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("ABHD_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eBHD_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eBHD_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eBHD_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulBHD", nBHD);
}
