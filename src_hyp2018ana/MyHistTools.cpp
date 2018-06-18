#include "MyHistTools.h"
#include "MyParam.h"

using namespace std;

void MyHistTools::fillTH(TString name, double val)
{
  TH1F *h1 = dynamic_cast<TH1F*>(gFile->Get(name.Data()));
  if( !h1 ){
    cout<<"  !!! TH1F "<<name<<" not found Break !!!"<<endl;
    exit(1);
  }

  h1-> Fill(val);
}


void MyHistTools::fillTH(TString name, double val1, double val2)
{
  TH2F *h2 = dynamic_cast<TH2F*>(gFile->Get(name.Data()));
  if( !h2 ){
    cout<<"  !!! TH2F "<<name<<" not found Break !!!"<<endl;
    exit(1);
  }

  h2-> Fill(val1, val2);
}

void MyHistTools::initFDC1()
{
  for( int lay=1; lay<=6; lay++ ){
    new TH1F(Form("hitpatFDC1_%d",lay), Form("FDC1 lay%d hit pattern", lay), 64, 0.5, 64.5);
    new TH1F(Form("nFDC1_%d",lay), Form("FDC1 lay%d multiplicity", lay), 64, 0.5, 64.5);
    new TH1F(Form("FDC1_%d_res", lay),  Form("FDC1 lay%d residual", lay), 500, -1.0, 1.0);
    for( int wire=1; wire<=64; wire++ ){
      new TH1F(Form("FDC1_%d_%d_tdc",lay, wire), Form("FDC1 TDC lay%d wire%d", lay, wire), 4000, 0, 4000);
      new TH1F(Form("FDC1_%d_%d_dt",lay, wire),  Form("FDC1 dt lay%d wire%d", lay, wire), 500, -100, 400);
      new TH1F(Form("FDC1_%d_%d_tdc_track",lay, wire), Form("FDC1 TDC lay%d wire%d", lay, wire), 4000, 0, 4000);
      new TH1F(Form("FDC1_%d_%d_dt_track",lay, wire),  Form("FDC1 dt lay%d wire%d", lay, wire), 500, -100, 400);
    }

    if( lay%2==0 ){
      new TH2F(Form("FDC1_%d_tsub_tmean", lay), Form("FDC1 lay%d tsub vs tmean", lay),       400, -200, 200, 400, -100, 300);
      new TH2F(Form("FDC1_%d_tsub_ctm",   lay), Form("FDC1 lay%d tsub vs corr. tmean", lay), 400, -200, 200, 400, -100, 300);
      for( int wire=1; wire<=64; wire++ ){
	new TH2F(Form("FDC1_%d_%d_tsub_tmean", lay, wire), Form("FDC1 l%d w%d tsub vs tmean", lay, wire),       400, -200, 200, 400, -100, 300);
	new TH2F(Form("FDC1_%d_%d_tsub_ctm",   lay, wire), Form("FDC1 l%d w%d tsub vs corr. tmean", lay, wire), 400, -200, 200, 400, -100, 300);
      }
    }
  }
}

void MyHistTools::fillFDC1(BeamLineHitMan *blMan)
{
  TH1F *h1;
  for( int lay=1; lay<=6; lay++ ){
    h1=(TH1F*)gFile->Get(Form("nFDC1_%d", lay)), h1->Fill(blMan->nFDC1(lay));
    for( int i=0; i<blMan->nFDC1(lay); i++ ){
      ChamberLikeHit *hit=blMan-> FDC1(lay, i);
      int lay=hit->layer(), wire=hit->wire(), tdc=hit->tdc();
      double dt=hit->dt();
      h1=(TH1F*)gFile->Get(Form("hitpatFDC1_%d", lay)), h1->Fill(wire);
      h1=(TH1F*)gFile->Get(Form("FDC1_%d_%d_tdc", lay, wire)), h1-> Fill(tdc);
      h1=(TH1F*)gFile->Get(Form("FDC1_%d_%d_dt", lay, wire)),  h1-> Fill(dt);
    }
  }
}

void MyHistTools::fillFDC1(LocalTrack *track)
{
  TH1F *h1;
  for( int i=0; i<track->nhit(); i++ ){
    ChamberLikeHit *hit = track->hit(i);
    int lay=hit->layer(), wire=hit->wire(), tdc=hit->tdc();
    double res=hit->resl(), dt=hit->dt();
    h1=(TH1F*)gFile->Get(Form("FDC1_%d_%d_tdc_track", lay, wire)), h1-> Fill(tdc);
    h1=(TH1F*)gFile->Get(Form("FDC1_%d_%d_dt_track", lay, wire)),  h1-> Fill(dt);
    h1=(TH1F*)gFile->Get(Form("FDC1_%d_res", lay)),  h1-> Fill(res);
  }

  for( int i=0; i<track->ncluster(0); i++ ){
    BLDCCluster *cluster = track->cluster(0, i);
    MyHistTools::fillTH(Form("FDC1_%d_tsub_tmean", cluster->hit(0)->layer()), cluster->GetTimeSub(), cluster->GetTimeMean());
    MyHistTools::fillTH(Form("FDC1_%d_tsub_ctm", cluster->hit(0)->layer()), cluster->GetTimeSub(), cluster->GetCTime());
  }
}

void MyHistTools::initCVC()
{
  new TH1F("nCVC", "CVC multiplicity", 6, -0.5, 5.5);
  new TH1F("hitpatCVC",   "CVC hit pattern",    35, 0.5, 34.5);
  new TH1F("hitpatCVCPC", "CVC PC hit pattern", 61, 0.5, 61.5);
  for( int seg=1; seg<=34; seg++ ){
    new TNtuple(Form("CVC_%d_slewing_info", seg), Form("CVC seg%d slewing info", seg), "eu:tu:ed:td"); 

    new TH1F(Form("CVC_%d_ua", seg), Form("CVC seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("CVC_%d_da", seg), Form("CVC seg%d down ADC", seg), 4000, 0, 4000);

    new TH1F(Form("CVC_%d_ua_wt", seg), Form("CVC seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("CVC_%d_da_wt", seg), Form("CVC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH2F(Form("CVC_%d_eu_ut", seg), Form("CVC seg%d dE_{up} vd tu", seg),         1000, 0, 100, 1000, -10, 10);
    new TH2F(Form("CVC_%d_ed_dt", seg), Form("CVC seg%d dE_{down} vs dt", seg),       1000, 0, 100, 1000, -10, 10); 
    new TH2F(Form("CVC_%d_eu_offset", seg), Form("CVC seg%d dE_{up} vs time", seg),   1000, 0, 100, 1000, -10, 10);
    new TH2F(Form("CVC_%d_ed_offset", seg), Form("CVC seg%d dE_{down} vs time", seg), 1000, 0, 100, 1000, -10, 10); 
    new TH2F(Form("CVC_%d_dE_offset", seg), Form("CVC seg%d dE vs time", seg),        1000, 0, 100, 1000, -10, 10);
  }
}

void MyHistTools::fillCVC(BeamLineHitMan *blMan)
{
  TH1F *h1;
  int nCVC=0;
  for( int i=0; i<blMan->nCVC(); i++ ){
    HodoscopeLikeHit *hit = blMan-> CVC(i);
    int seg = hit->seg();
    int ua = hit->adcu(), ud=hit->adcd();
    h1 = (TH1F*)gFile-> Get(Form("CVC_%d_ua", seg)), h1-> Fill(ua);
    h1 = (TH1F*)gFile-> Get(Form("CVC_%d_da", seg)), h1-> Fill(ud);

    if( hit-> CheckRange() ){
      nCVC++;
      h1 = (TH1F*)gFile-> Get("hitpatCVC"), h1-> Fill(seg);
      h1 = (TH1F*)gFile-> Get("hitpatCVCPC"), h1-> Fill(seg);
      h1 = (TH1F*)gFile-> Get(Form("CVC_%d_ua_wt", seg)), h1-> Fill(ua);
      h1 = (TH1F*)gFile-> Get(Form("CVC_%d_da_wt", seg)), h1-> Fill(ud);
    }
  }
  h1 = (TH1F*)gFile-> Get("nCVC"), h1-> Fill(nCVC);
}

void MyHistTools::initPC()
{
  new TH1F("nPC",      "PC multiplicity", 28, -0.5, 27.5);
  new TH1F("hitpatPC", "PC hit pattern",  61,  0.5, 61.5);
  for( int seg=1; seg<=27; seg++ ){
    new TNtuple(Form("PC_%d_slewing_info", seg), Form("PC seg%d slewing info", seg), "eu:tu:ed:td"); 

    new TH1F(Form("PC_%d_ua", seg), Form("PC seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("PC_%d_da", seg), Form("PC seg%d down ADC", seg), 4000, 0, 4000);

    new TH1F(Form("PC_%d_ua_wt", seg), Form("PC seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("PC_%d_da_wt", seg), Form("PC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH2F(Form("PC_%d_eu_ut", seg), Form("PC seg%d dE_{up} vs ut", seg),         1000, 0, 100, 1000, -10, 10);
    new TH2F(Form("PC_%d_ed_dt", seg), Form("PC seg%d dE_{down} vs dt", seg),       1000, 0, 100, 1000, -10, 10); 
    new TH2F(Form("PC_%d_eu_offset", seg), Form("PC seg%d dE_{up} vs time", seg),   1000, 0, 100, 1000, -10, 10);
    new TH2F(Form("PC_%d_ed_offset", seg), Form("PC seg%d dE_{down} vs time", seg), 1000, 0, 100, 1000, -10, 10); 
    new TH2F(Form("PC_%d_dE_offset", seg), Form("PC seg%d dE vs time", seg),        1000, 0, 100, 1000, -10, 10);
  }
}

void MyHistTools::fillPC(BeamLineHitMan *blMan)
{
  TH1F *h1;
  int nPC=0;
  for( int i=0; i<blMan->nPC(); i++ ){
    HodoscopeLikeHit *hit = blMan-> PC(i);
    int seg = hit->seg();
    int ua = hit->adcu(), ud=hit->adcd();
    h1 = (TH1F*)gFile-> Get(Form("PC_%d_ua", seg)), h1-> Fill(ua);
    h1 = (TH1F*)gFile-> Get(Form("PC_%d_da", seg)), h1-> Fill(ud);
    if( hit-> CheckRange() ){
      nPC++;
      h1 = (TH1F*)gFile-> Get("hitpatPC"), h1-> Fill(seg);
      h1 = (TH1F*)gFile-> Get("hitpatCVCPC"), h1-> Fill(34+seg);
      h1 = (TH1F*)gFile-> Get(Form("PC_%d_ua_wt", seg)), h1-> Fill(ua);
      h1 = (TH1F*)gFile-> Get(Form("PC_%d_da_wt", seg)), h1-> Fill(ud);
    }
  }
  h1 = (TH1F*)gFile-> Get("nPC"), h1-> Fill(nPC);
}

void MyHistTools::fillT0PC(const BeamInfo &beam, const ForwardChargeInfo &forward, BeamLineHitMan *blMan)
{
  TH2F *h2;
  TNtuple *tup;
  HodoscopeLikeHit *T0hit=beam.T0(blMan);
  HodoscopeLikeHit *PChit=forward.hodo(blMan);
  if( PChit->cid()!=CID_PC ){
    cout<<"  !!! PC not found !!!"<<endl;
    return;
  }
  int PCseg=PChit->seg();
  double tof=PChit->ctmean()-T0hit->ctmean();
  double ctu=PChit->ctu()-T0hit->ctmean();
  double ctd=PChit->ctd()-T0hit->ctmean();
  double tu=PChit->tu()-T0hit->ctmean();
  double td=PChit->td()-T0hit->ctmean();
  double calc_tof=forward.calc_tof();
  double eu=PChit->eu(), ed=PChit->ed(), dE=PChit->emean();

  tup = (TNtuple*)gFile->Get(Form("PC_%d_slewing_info", PCseg)), tup-> Fill(eu, tu-calc_tof, ed, td-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("PC_%d_eu_ut", PCseg)),     h2-> Fill(eu, ctu-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("PC_%d_ed_dt", PCseg)),     h2-> Fill(ed, ctd-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("PC_%d_eu_offset", PCseg)), h2-> Fill(eu, tof-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("PC_%d_ed_offset", PCseg)), h2-> Fill(ed, tof-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("PC_%d_dE_offset", PCseg)), h2-> Fill(dE, tof-calc_tof);
}

void MyHistTools::fillT0CVC(const BeamInfo &beam, const ForwardChargeInfo &forward, BeamLineHitMan *blMan)
{
  TH2F *h2;
  TNtuple *tup;
  HodoscopeLikeHit *T0hit=beam.T0(blMan);
  HodoscopeLikeHit *CVChit=forward.hodo(blMan);
  if( CVChit->cid()!=CID_CVC ){
    cout<<"  !!! CVC not found !!!"<<endl;
    return;
  }
  int CVCseg=CVChit->seg();
  double tof=CVChit->ctmean()-T0hit->ctmean();
  double ctu=CVChit->ctu()-T0hit->ctmean();
  double ctd=CVChit->ctd()-T0hit->ctmean();
  double tu=CVChit->tu()-T0hit->ctmean();
  double td=CVChit->td()-T0hit->ctmean();
  double calc_tof=forward.calc_tof();
  double eu=CVChit->eu(), ed=CVChit->ed(), dE=CVChit->emean();

  tup = (TNtuple*)gFile->Get(Form("CVC_%d_slewing_info", CVCseg)), tup-> Fill(eu, tu-calc_tof, ed, td-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("CVC_%d_eu_ut", CVCseg)),     h2-> Fill(eu, ctu-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("CVC_%d_ed_dt", CVCseg)),     h2-> Fill(ed, ctd-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("CVC_%d_eu_offset", CVCseg)), h2-> Fill(eu, tof-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("CVC_%d_ed_offset", CVCseg)), h2-> Fill(ed, tof-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("CVC_%d_dE_offset", CVCseg)), h2-> Fill(dE, tof-calc_tof);
}

void MyHistTools::initT0()
{
  new TH1F("nT0",      "T0 multiplicity", 6, -0.5, 5.5);
  new TH1F("hitpatT0", "T0 hit pattern",  5,  0.5, 5.5);
  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("T0_%d_ua", seg), Form("T0 seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("T0_%d_da", seg), Form("T0 seg%d down ADC", seg), 4000, 0, 4000);

    new TH1F(Form("T0_%d_ua_wt", seg), Form("T0 seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("T0_%d_da_wt", seg), Form("T0 seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH2F(Form("T0_%d_eu_tu", seg), Form("T0 seg%d dE_{up} vs ut", seg),         100, 0, 10, 1000, -10, 10);
    new TH2F(Form("T0_%d_ed_td", seg), Form("T0 seg%d dE_{down} vs dt", seg),       100, 0, 10, 1000, -10, 10); 
    new TH2F(Form("T0_%d_eu_offset", seg), Form("T0 seg%d dE_{up} vs time", seg),   100, 0, 10, 1000, -10, 10);
    new TH2F(Form("T0_%d_ed_offset", seg), Form("T0 seg%d dE_{down} vs time", seg), 100, 0, 10, 1000, -10, 10); 
    new TH2F(Form("T0_%d_dE_offset", seg), Form("T0 seg%d dE vs time", seg),        100, 0, 10, 1000, -10, 10);
  }
}

void MyHistTools::fillT0(BeamLineHitMan *blMan)
{
  TH1F *h1;
  int nT0=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    HodoscopeLikeHit *hit = blMan-> T0(i);
    int seg = hit->seg();
    int ua = hit->adcu(), ud=hit->adcd();
    h1 = (TH1F*)gFile-> Get(Form("T0_%d_ua", seg)), h1-> Fill(ua);
    h1 = (TH1F*)gFile-> Get(Form("T0_%d_da", seg)), h1-> Fill(ud);
    if( hit-> CheckRange() ){
      nT0++;
      h1 = (TH1F*)gFile-> Get("hitpatT0"), h1-> Fill(seg);
      h1 = (TH1F*)gFile-> Get(Form("T0_%d_ua_wt", seg)), h1-> Fill(ua);
      h1 = (TH1F*)gFile-> Get(Form("T0_%d_da_wt", seg)), h1-> Fill(ud);
    }
  }
  h1 = (TH1F*)gFile-> Get("nT0"), h1-> Fill(nT0);
}

void MyHistTools::initBHD()
{
  new TH1F("nBHD",      "BHD multiplicity", 21, -0.5, 20.5);
  new TH1F("hitpatBHD", "BHD hit pattern",  20,  0.5, 20.5);
  for( int seg=1; seg<=20; seg++ ){
    new TH1F(Form("BHD_%d_ua", seg), Form("BHD seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("BHD_%d_da", seg), Form("BHD seg%d down ADC", seg), 4000, 0, 4000);

    new TH1F(Form("BHD_%d_ua_wt", seg), Form("BHD seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("BHD_%d_da_wt", seg), Form("BHD seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH2F(Form("BHD_%d_eu_tu", seg), Form("BHD seg%d dE_{up} vs ut", seg),         200, 0, 10, 1000, -10, 10);
    new TH2F(Form("BHD_%d_ed_td", seg), Form("BHD seg%d dE_{down} vs dt", seg),       200, 0, 10, 1000, -10, 10); 
    new TH2F(Form("BHD_%d_eu_offset", seg), Form("BHD seg%d dE_{up} vs time", seg),   200, 0, 10, 1000, -10, 10);
    new TH2F(Form("BHD_%d_ed_offset", seg), Form("BHD seg%d dE_{down} vs time", seg), 200, 0, 10, 1000, -10, 10); 
    new TH2F(Form("BHD_%d_dE_offset", seg), Form("BHD seg%d dE vs time", seg),        200, 0, 10, 1000, -10, 10);
  }
}

void MyHistTools::fillBHDT0(double tof, EventHeader *header)
{
  TH1F *h1;
  h1 = (TH1F*)gFile-> Get("BHDT0"), h1-> Fill(tof);
  if( header->kaon() )  h1 = (TH1F*)gFile-> Get("BHDT0_K"), h1-> Fill(tof);
  if( header->pion() )  h1 = (TH1F*)gFile-> Get("BHDT0_pi"), h1-> Fill(tof);
  if( header->proton() )  h1 = (TH1F*)gFile-> Get("BHDT0_p"), h1-> Fill(tof);
  
  if( header->IsTrig(Trig_Kf) )  h1 = (TH1F*)gFile-> Get("BHDT0_Kf"), h1-> Fill(tof);
}

void MyHistTools::initBHDT0()
{
  new TH1F("BHDT0",    "BHD-T0 TOF", 10000, 0, 50);
  new TH1F("BHDT0_K",  "BHD-T0 TOF trig K",  10000, 0, 50);
  new TH1F("BHDT0_pi", "BHD-T0 TOF trig pi", 10000, 0, 50);
  new TH1F("BHDT0_p",  "BHD-T0 TOF trig p",  10000, 0, 50);
  
  new TH1F("BHDT0_Kf", "BHD-T0 TOF trig K/f",  10000, 0, 50);
}

void MyHistTools::fillBHD(BeamLineHitMan *blMan)
{
  TH1F *h1;
  int nBHD=0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    HodoscopeLikeHit *hit = blMan-> BHD(i);
    int seg = hit->seg();
    int ua = hit->adcu(), ud=hit->adcd();
    h1 = (TH1F*)gFile-> Get(Form("BHD_%d_ua", seg)), h1-> Fill(ua);
    h1 = (TH1F*)gFile-> Get(Form("BHD_%d_da", seg)), h1-> Fill(ud);
    if( hit-> CheckRange() ){
      nBHD++;
      h1 = (TH1F*)gFile-> Get("hitpatBHD"), h1-> Fill(seg);
      h1 = (TH1F*)gFile-> Get(Form("BHD_%d_ua_wt", seg)), h1-> Fill(ua);
      h1 = (TH1F*)gFile-> Get(Form("BHD_%d_da_wt", seg)), h1-> Fill(ud);
    }
  }
  h1 = (TH1F*)gFile-> Get("nBHD"), h1-> Fill(nBHD);
}

void MyHistTools::fillBHDT0(const BeamInfo &beam, BeamLineHitMan *blMan)
{
  TH2F *h2;
  HodoscopeLikeHit *T0hit=beam.T0(blMan);
  HodoscopeLikeHit *BHDhit=beam.BHD(blMan);
  double calc_tof=beam.calc_tof();
  int T0seg=T0hit->seg(), BHDseg=BHDhit->seg();
  double T0ctu=T0hit->ctu(), T0ctd=T0hit->ctd(), T0time=T0hit->ctmean();
  double T0eu=T0hit->eu(),  T0ed=T0hit->ed(),  T0dE=T0hit->emean();
  double BHDctu=BHDhit->ctu(), BHDctd=BHDhit->ctd(), BHDtime=BHDhit->ctmean();
  double BHDeu=BHDhit->eu(),  BHDed=BHDhit->ed(),  BHDdE=BHDhit->emean();

  h2 = (TH2F*)gFile->Get(Form("BHD_%d_eu_tu", BHDseg)),     h2-> Fill(BHDeu, BHDctu-T0time+calc_tof);
  h2 = (TH2F*)gFile->Get(Form("BHD_%d_ed_td", BHDseg)),     h2-> Fill(BHDed, BHDctd-T0time+calc_tof);
  h2 = (TH2F*)gFile->Get(Form("BHD_%d_eu_offset", BHDseg)), h2-> Fill(BHDeu, BHDtime-T0time+calc_tof);
  h2 = (TH2F*)gFile->Get(Form("BHD_%d_ed_offset", BHDseg)), h2-> Fill(BHDed, BHDtime-T0time+calc_tof);
  h2 = (TH2F*)gFile->Get(Form("BHD_%d_dE_offset", BHDseg)), h2-> Fill(BHDdE, BHDtime-T0time+calc_tof);

  h2 = (TH2F*)gFile->Get(Form("T0_%d_eu_tu", T0seg)),     h2-> Fill(T0eu, T0ctu-BHDtime-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("T0_%d_ed_td", T0seg)),     h2-> Fill(T0ed, T0ctd-BHDtime-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("T0_%d_eu_offset", T0seg)), h2-> Fill(T0eu, T0time-BHDtime-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("T0_%d_ed_offset", T0seg)), h2-> Fill(T0ed, T0time-BHDtime-calc_tof);
  h2 = (TH2F*)gFile->Get(Form("T0_%d_dE_offset", T0seg)), h2-> Fill(T0dE, T0time-BHDtime-calc_tof);
}

