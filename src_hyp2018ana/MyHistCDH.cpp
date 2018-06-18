#include "MyHistCDH.h"

using namespace std;

void initHistCDH()
{
  const int nseg=36;
  new TH1F("hitpatCDH", "hitpatCDH", nseg, 0.5, nseg+0.5);
  new TH1F("mulCDH", "mulCDH", nseg+1, -0.5, nseg+0.5);

  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("ACDH_u%d", seg), Form("CDH seg%d up ADC",seg),   4000, 0, 4000);  
    new TH1F(Form("ACDH_d%d", seg), Form("CDH seg%d down ADC",seg), 4000, 0, 4000);  
    new TH1F(Form("ACDH_u%d_wt", seg), Form("CDH seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ACDH_d%d_wt", seg), Form("CDH seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eCDH_u%d", seg), Form("CDH seg%d up dE", seg), 500, -0.5, 95.5);
    new TH1F(Form("eCDH_d%d", seg), Form("CDH seg%d down dE", seg), 500, -0.5, 95.5);
    new TH1F(Form("eCDH_%d",  seg), Form("CDH seg%d dE", seg), 500, -0.5, 95.5);

    new TH1F(Form("eCDH_u%d_wt", seg), Form("CDH seg%d up dE", seg), 500, -0.5, 95.5);
    new TH1F(Form("eCDH_d%d_wt", seg), Form("CDH seg%d down dE", seg), 500, -0.5, 95.5);
    new TH1F(Form("eCDH_%d_wt",  seg), Form("CDH seg%d dE", seg), 500, -0.5, 95.5);

    new TH1F(Form("ACDH_u%d_wt_pi", seg), Form("CDH seg%d up ADC w/ TDC",seg),   4000, 0, 4000);  
    new TH1F(Form("ACDH_d%d_wt_pi", seg), Form("CDH seg%d down ADC w/ TDC", seg), 4000, 0, 4000);  

    new TH1F(Form("eCDH_u%d_wt_pi", seg), Form("CDH seg%d up dE", seg), 500, -0.5, 95.5);
    new TH1F(Form("eCDH_d%d_wt_pi", seg), Form("CDH seg%d down dE", seg), 500, -0.5, 95.5);
    new TH1F(Form("eCDH_%d_wt_pi",  seg), Form("CDH seg%d dE", seg), 500, -0.5, 95.5);

    new TH1F(Form("CDH%d_offset",     seg), Form("CDH seg%d offset", seg), 1000, -5, 5);
    new TH1F(Form("CDH%d_offset_pim", seg), Form("CDH seg%d offset", seg), 1000, -5, 5);
    new TH1F(Form("CDH%d_offset_pip", seg), Form("CDH seg%d offset", seg), 1000, -5, 5);
    new TH1F(Form("CDH%d_offset_km",  seg), Form("CDH seg%d offset", seg), 1000, -5, 5);
    new TH1F(Form("CDH%d_offset_p",   seg), Form("CDH seg%d offset", seg), 1000, -5, 5);
    for( int seg2=1; seg2<=5; seg2++ ){
      new TH1F(Form("CDH%d_offset_T0%d",     seg, seg2), Form("CDH seg%d offset", seg), 1000, -5, 5);
      new TH1F(Form("CDH%d_offset_T0%d_pim", seg, seg2), Form("CDH seg%d offset", seg), 1000, -5, 5);
      new TH1F(Form("CDH%d_offset_T0%d_pip", seg, seg2), Form("CDH seg%d offset", seg), 1000, -5, 5);
      new TH1F(Form("CDH%d_offset_T0%d_km",  seg, seg2), Form("CDH seg%d offset", seg), 1000, -5, 5);
      new TH1F(Form("CDH%d_offset_T0%d_p",   seg, seg2), Form("CDH seg%d offset", seg), 1000, -5, 5);
    }
    new TNtuple(Form("CDH%d_slewing_info", seg), Form("CDH seg%d slewing info", seg), "eu:ed:tu:td:ctu:ctd:ctm");
  }

  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("T0%d_offset_CDH", seg), Form("T0 seg%d offset by CDH", seg),1000, -5, 5);
    new TH1F(Form("T0%d_offset_CDH_pim", seg), Form("T0 seg%d offset by CDH", seg), 1000, -5, 5);
    new TH1F(Form("T0%d_offset_CDH_pip", seg), Form("T0 seg%d offset by CDH", seg), 1000, -5, 5);
    new TH1F(Form("T0%d_offset_CDH_km", seg), Form("T0 seg%d offset by CDH", seg), 1000, -5, 5);
    new TH1F(Form("T0%d_offset_CDH_p", seg), Form("T0 seg%d offset by CDH", seg), 1000, -5, 5);
  }
}

void fillCDH_ADC(CDSHitMan *cdsMan)
{
  int nCDH=0;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    HodoscopeLikeHit *hit= cdsMan->CDH(i);
    int seg=hit-> seg();
    int tu=hit->tdcu(), td=hit->tdcd(), au=hit->adcu(), ad=hit->adcd();
    double eu=hit->eu(), ed=hit->ed(), emean=hit->emean();

    MyHistTools::fillTH(Form("ACDH_u%d", seg), au);
    MyHistTools::fillTH(Form("ACDH_d%d", seg), ad);

    MyHistTools::fillTH(Form("eCDH_u%d", seg), eu);
    MyHistTools::fillTH(Form("eCDH_d%d", seg), ed);
    MyHistTools::fillTH(Form("eCDH_%d", seg), emean);

    if( hit->CheckRange() ){
      nCDH++;
      MyHistTools::fillTH("hitpatCDH", seg);
      MyHistTools::fillTH(Form("ACDH_u%d_wt", seg), au);
      MyHistTools::fillTH(Form("ACDH_d%d_wt", seg), ad);

      MyHistTools::fillTH(Form("eCDH_u%d_wt", seg), eu);
      MyHistTools::fillTH(Form("eCDH_d%d_wt", seg), ed);
      MyHistTools::fillTH(Form("eCDH_%d_wt", seg), emean);
    }
  }
  MyHistTools::fillTH("mulCDH", nCDH);
}

void fillCDH(BeamInfo *beam, CDSInfo *info, CDSHitMan *cdsMan)
{

  if( info->flag() && (info->pid()==CDS_PiMinus || info->pid()==CDS_PiPlus) && fabs(info->mom())>0.3 ){
    HodoscopeLikeHit *hit = info->CDH(cdsMan);
    MyHistTools::fillTH(Form("ACDH_u%d_wt_pi", hit->seg()), hit->adcu());
    MyHistTools::fillTH(Form("ACDH_d%d_wt_pi", hit->seg()), hit->adcd());

    MyHistTools::fillTH(Form("eCDH_u%d_wt_pi", hit->seg()), hit->eu());
    MyHistTools::fillTH(Form("eCDH_d%d_wt_pi", hit->seg()), hit->ed());
    MyHistTools::fillTH(Form("eCDH_%d_wt_pi", hit->seg()), hit->emean());
  }

  if( info->flag() && info-> dca()<1.0 ){
    int T0seg=beam->T0seg();
    HodoscopeLikeHit *hit=info->CDH(cdsMan);
    double eu=hit->eu(), ed=hit->ed();
    double tu=hit->tu()-beam->T0time()-info->offset();
    double td=hit->td()-beam->T0time()-info->offset();
    double ctu=hit->tu()-beam->T0time()-info->offset();
    double ctd=hit->td()-beam->T0time()-info->offset();
    double offset=hit->ctmean()-beam->T0time()-info->offset();
    int seg=hit->seg();

    TNtuple *tup=(TNtuple*)gFile->Get(Form("CDH%d_slewing_info", seg));
    tup->Fill(eu, ed, tu, td, ctu, ctd, offset);
    //    cout<<" CDH seg"<<seg<<"  offset : "<<offset<<endl;
    MyHistTools::fillTH(Form("CDH%d_offset", seg), offset);
    MyHistTools::fillTH(Form("CDH%d_offset_T0%d", seg, T0seg), offset);
    MyHistTools::fillTH(Form("T0%d_offset_CDH", T0seg), -offset);
    if( info->pid()==CDS_PiMinus ){
      MyHistTools::fillTH(Form("CDH%d_offset_pim", seg), offset);
      MyHistTools::fillTH(Form("CDH%d_offset_T0%d_pim", seg, T0seg), offset);
      MyHistTools::fillTH(Form("T0%d_offset_CDH_pim", T0seg), -offset);
    }
    if( info->pid()==CDS_PiPlus  ){
      MyHistTools::fillTH(Form("CDH%d_offset_pip", seg), offset);
      MyHistTools::fillTH(Form("CDH%d_offset_T0%d_pip", seg, T0seg), offset);
      MyHistTools::fillTH(Form("T0%d_offset_CDH_pip", T0seg), -offset);
    }
    if( info->pid()==CDS_Kaon    ){
      MyHistTools::fillTH(Form("CDH%d_offset_km",  seg), offset);
      MyHistTools::fillTH(Form("CDH%d_offset_T0%d_km", seg, T0seg), offset);
      MyHistTools::fillTH(Form("T0%d_offset_CDH_km", T0seg), -offset);
    }
    if( info->pid()==CDS_Proton  ){
      MyHistTools::fillTH(Form("CDH%d_offset_p",   seg), offset);
      MyHistTools::fillTH(Form("CDH%d_offset_T0%d_p", seg, T0seg), offset);
      MyHistTools::fillTH(Form("T0%d_offset_CDH_p", T0seg), -offset);
    }
  }
}
