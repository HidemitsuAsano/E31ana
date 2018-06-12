#include "MyHistBVC.h"

void initHistBVC()
{
  const int nseg=8;
  new TH1F("hitpatBVC", "hitpatBVC", nseg, 0.5, nseg+0.5);
  new TH1F("mulBVC", "mulBVC", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("ABVC_u%d", seg), Form("BVC seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("ABVC_d%d", seg), Form("BVC seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("ABVC_u%d_wt", seg), Form("BVC seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ABVC_d%d_wt", seg), Form("BVC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eBVC_u%d", seg), Form("BVC seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eBVC_d%d", seg), Form("BVC seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eBVC_%d",  seg), Form("BVC seg%d dE", seg), 500, -0.5, 4.5);

    new TH1F(Form("eBVC_u%d_wt", seg), Form("BVC seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eBVC_d%d_wt", seg), Form("BVC seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eBVC_%d_wt",  seg), Form("BVC seg%d dE", seg), 500, -0.5, 4.5);
  }
}

void fillBVC_ADC(BeamLineHitMan *blMan)
{
  int nBVC=0;
  for( int i=0; i<blMan->nBVC(); i++ ){
    HodoscopeLikeHit *hit= blMan->BVC(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("ABVC_u%d", seg), au);
    MyHistTools::fillTH(Form("ABVC_d%d", seg), ad);

    MyHistTools::fillTH(Form("eBVC_u%d", seg), eu);
    MyHistTools::fillTH(Form("eBVC_d%d", seg), ed);
    MyHistTools::fillTH(Form("eBVC_%d", seg), emean);

    if( hit->CheckRange() ){
      nBVC++;
      MyHistTools::fillTH("hitpatBVC", seg);
      MyHistTools::fillTH(Form("ABVC_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("ABVC_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eBVC_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eBVC_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eBVC_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulBVC", nBVC);
}
