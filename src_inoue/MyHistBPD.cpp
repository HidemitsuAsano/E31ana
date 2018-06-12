#include "MyHistBPD.h"

void initHistBPD()
{
  const int nseg=70;
  new TH1F("hitpatBPD", "hitpatBPD", nseg, 0.5, nseg+0.5);
  new TH1F("mulBPD", "mulBPD", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("ABPD_u%d", seg), Form("BPD seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("ABPD_d%d", seg), Form("BPD seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("ABPD_u%d_wt", seg), Form("BPD seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ABPD_d%d_wt", seg), Form("BPD seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eBPD_u%d", seg), Form("BPD seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eBPD_d%d", seg), Form("BPD seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eBPD_%d",  seg), Form("BPD seg%d dE", seg), 500, -0.5, 4.5);

    new TH1F(Form("eBPD_u%d_wt", seg), Form("BPD seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eBPD_d%d_wt", seg), Form("BPD seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("eBPD_%d_wt",  seg), Form("BPD seg%d dE", seg), 500, -0.5, 4.5);
  }
}

void fillBPD_ADC(BeamLineHitMan *blMan)
{
  int nBPD=0;
  for( int i=0; i<blMan->nBPD(); i++ ){
    HodoscopeLikeHit *hit= blMan->BPD(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("ABPD_u%d", seg), au);
    MyHistTools::fillTH(Form("ABPD_d%d", seg), ad);

    MyHistTools::fillTH(Form("eBPD_u%d", seg), eu);
    MyHistTools::fillTH(Form("eBPD_d%d", seg), ed);
    MyHistTools::fillTH(Form("eBPD_%d", seg), emean);

    if( hit->CheckRange() ){
      nBPD++;
      MyHistTools::fillTH("hitpatBPD", seg);
      MyHistTools::fillTH(Form("ABPD_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("ABPD_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eBPD_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eBPD_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eBPD_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulBPD", nBPD);
}
