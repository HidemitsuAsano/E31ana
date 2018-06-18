#include "MyHistCVC.h"

void initHistCVC()
{
  const int nseg=34;
  new TH1F("hitpatCVC", "hitpatCVC", nseg, 0.5, nseg+0.5);
  new TH1F("mulCVC", "mulCVC", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("ACVC_u%d", seg), Form("CVC seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("ACVC_d%d", seg), Form("CVC seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("ACVC_u%d_wt", seg), Form("CVC seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ACVC_d%d_wt", seg), Form("CVC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eCVC_u%d", seg), Form("CVC seg%d up dE", seg),      500, -10, 100);
    new TH1F(Form("eCVC_d%d", seg), Form("CVC seg%d down dE", seg),    500, -10, 100);
    new TH1F(Form("eCVC_%d",  seg), Form("CVC seg%d dE", seg),         500, -10, 100);

    new TH1F(Form("eCVC_u%d_wt", seg), Form("CVC seg%d up dE", seg),   550, -10, 100);
    new TH1F(Form("eCVC_d%d_wt", seg), Form("CVC seg%d down dE", seg), 550, -10, 100);
    new TH1F(Form("eCVC_%d_wt",  seg), Form("CVC seg%d dE", seg),      550, -10, 100);
  }
}

void fillCVC_ADC(BeamLineHitMan *blMan)
{
  int nCVC=0;
  for( int i=0; i<blMan->nCVC(); i++ ){
    HodoscopeLikeHit *hit= blMan->CVC(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("ACVC_u%d", seg), au);
    MyHistTools::fillTH(Form("ACVC_d%d", seg), ad);

    MyHistTools::fillTH(Form("eCVC_u%d", seg), eu);
    MyHistTools::fillTH(Form("eCVC_d%d", seg), ed);
    MyHistTools::fillTH(Form("eCVC_%d", seg), emean);

    if( hit->CheckRange() ){
      nCVC++;
      MyHistTools::fillTH("hitpatCVC", seg);
      MyHistTools::fillTH(Form("ACVC_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("ACVC_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eCVC_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eCVC_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eCVC_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulCVC", nCVC);
}
