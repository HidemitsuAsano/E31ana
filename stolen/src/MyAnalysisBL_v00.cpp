// MyAnalysisBL.cpp

#include "MyAnalysisBL.h"

MyAnalysisBL::MyAnalysisBL(TFile* rtFile, ConfMan* conf)
{
  D5 = new BeamSpectrometer(conf);
  Initialize(rtFile, conf);
  Clear();
}

MyAnalysisBL::~MyAnalysisBL()
{
  if(D5) delete D5;
}

bool MyAnalysisBL::Clear()
{
  BeamPID = Beam_Other;
  BeamMom = 0.0; 
  D5 -> Clear();
}

bool MyAnalysisBL::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{
  DetectorList *dlist=DetectorList::GetInstance();

  // T0 timing
  double T0Timing[5]; for(int i=0;i<5;i++){T0Timing[i]=-999.;}
  for(int i=0; i<blMan->nT0(); i++){
    HodoscopeLikeHit* hit = blMan->T0(i);
    int seg = hit->seg();
    if(hit->CheckRange()) T0Timing[seg-1] = hit->ctmean();
  }

  // Hodoscope
  const int nhodo=4;
  int hodoid[nhodo]={CID_BHD,CID_T0,CID_BPD,CID_DEF};
  int hodonHit[nhodo]; for(int i=0; i<nhodo; i++) hodonHit[i]=0;
  int hodoseg[nhodo][100];
  double hodoemean[nhodo][100];
  double hodoctmean[nhodo][100];
  double hodohitpos[nhodo][100];
  TVector3 hodopos[nhodo][100];
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist->GetName(cid).data();
    const int nsegs = dlist->GetNsegs(cid);
    int ihit = 0;
    for(int i=0; i<blMan->nHodo(cid); i++){
      HodoscopeLikeHit* hit = blMan->Hodoi(cid,i);
      if(hit->CheckRange()){
        hodoseg[ihodo][ihit] = hit->seg();
        hodoemean[ihodo][ihit] = hit->emean();
        hodoctmean[ihodo][ihit] = hit->ctmean();
        hodohitpos[ihodo][ihit] = hit->hitpos();
        FillHist(Form("%s_HitPat",name),hodoseg[ihodo][ihit]);
        FillHist(Form("%s%d_dEMean",name,hodoseg[ihodo][ihit]),hodoemean[ihodo][ihit]);
        FillHist(Form("%s%d_CTimeMean",name,hodoseg[ihodo][ihit]),hodoctmean[ihodo][ihit]);
        FillHist(Form("%s%d_Position",name,hodoseg[ihodo][ihit]),hodohitpos[ihodo][ihit]);
        FillHist(Form("%s%d_dEMeanvCTimeMean",name,hodoseg[ihodo][ihit]),hodoemean[ihodo][ihit],hodoctmean[ihodo][ihit]);
        FillHist(Form("%s_SegmentvPosition",name),hodoseg[ihodo][ihit],hodohitpos[ihodo][ihit]);
        for(int t0seg=0; t0seg<5; t0seg++){
          if(T0Timing[t0seg]>-900){
            double tof = hodoctmean[ihodo][ihit]-T0Timing[t0seg]; if(ihodo==0){ tof=tof*-1.; }
            FillHist(Form("T0%d%s_TOF",t0seg+1,name),hodoseg[ihodo][ihit],tof);
            FillHist(Form("T0%s_TOF",name),hodoseg[ihodo][ihit],tof);
            char* ud[2] = {"u","d"};
            for(int iud=0; iud<2; iud++){
              FillHist(Form("%s%s%d_dEvTOF",name,ud[iud],hodoseg[ihodo][ihit]),hit->ene(iud),tof);
            }
            FillHist(Form("%s%d_dEMeanvTOF",name,hodoseg[ihodo][i]),hodoemean[ihodo][i],tof);
            FillHist(Form("T0%d%s_SegmentvPosition",t0seg+1,name),hodoseg[ihodo][i],hodohitpos[ihodo][i]);
          } 
        }
        hodonHit[ihodo]++;
        ihit++;
      }
    }
    FillHist(Form("%s_Mult",name),hodonHit[ihodo]);
  }

  // Chamber
  const int ndc = 5;
  int dccid[ndc] = {CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_BPC};
  for(int idc=0;idc<ndc;idc++){
    const int cid = dccid[idc];
    const char* name = dlist->GetName(cid).data();
    const int nlays = dlist->GetNlayers(cid);
    for( int layer=1; layer<=nlays; layer++ ){
      FillHist(Form("%s_Layer%d_Mult",name,layer),blMan->nBLDC(cid,layer)); 
      for( int i=0; i<blMan->nBLDC(cid,layer); i++ ){
        ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
        double wire = hit->wire();
        double dt = hit->dt();
        double dl = hit->dl();
        FillHist(Form("%s_Layer%d_Time",name,layer),dt);
        FillHist(Form("%s_Layer%d_Length",name,layer),dl);
        FillHist(Form("%s_Layer%d_HitPat",name,layer),wire);
      }
    }
  }

  // Cherenkov
  const int nchere=1;
  int chereid[nchere]={CID_AC};
  char* cherename[nchere]={"AC"};
  double chereadc[nchere]; 
  for(int ichere=0; ichere<nchere; ichere++){
    double adcsum = 0;
    double tdc = 0;
    for(int i=0; i<blMan->nChere(chereid[ichere]); i++){
      CherenkovLikeHit* hit = blMan->Cherei(chereid[ichere],i);
      tdc = hit->tdc(1);
      FillHist(Form("%s_TDC",cherename[ichere]),tdc);
      for(int ich=1; ich<=4; ich++){
        adcsum += hit->adc(ich);
      }
    }
    FillHist(Form("%s_ADC",cherename[ichere]),adcsum);
  }

  // Beam Momentum
  if(bltrackMan->ntrackBLC1()==1 && bltrackMan->ntrackBLC2()==1){
    D5 -> TMinuitFit(bltrackMan->trackBLC1(0), bltrackMan->trackBLC2(0),conf);
    BeamMom = D5 -> mom();
    FillHist(Form("BeamMom"), BeamMom);
    FillHist(Form("BeamMom_Chi2"), D5->chisquare());
    for(int ihodo=0; ihodo<nhodo; ihodo++){
      for(int i=0; i<hodonHit[ihodo]; i++){
        const char* name = dlist->GetName(hodoid[ihodo]).data();
        FillHist(Form("%sSegvBeamMom",name),hodoseg[ihodo][i],BeamMom); 
      }
    }
  }

  BeamPID = Beam_Kaon;

  return true;

}

bool MyAnalysisBL::FillHist(TString name, double val1)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1){
    h1 -> Fill(val1);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBL::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){
    h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBL::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  // Hodoscope
  const int nhodo = 4;
  int NumOfSegments[nhodo] = {20,5,70,8};
  TString CounterName[nhodo] = {"BHD","T0","BPD","DEF"};
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    std::cout << "Define Histograms for " << CounterName[ihodo].Data() << std::endl;
    new TH1F( Form("%s_HitPat",CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
    new TH2F( Form("%s_SegmentvPosition",CounterName[ihodo].Data()), Form("Hit Position %s;%s segment;Position (cm)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1, 201, -50.25, 50.25 );
    for(int seg=1; seg<=NumOfSegments[1]; seg++){
      new TH1F( Form("T0%d%s_HitPat",seg,CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
      new TH2F( Form("T0%d%s_SegmentvPosition",seg,CounterName[ihodo].Data()), Form("Hit Position %s;%s segment;Position (cm)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1, 201, -50.25, 50.25 );
    }
    for( int seg=1; seg<=NumOfSegments[ihodo]; seg++ ){
      new TH1F( Form("%s%d_dEMean",CounterName[ihodo].Data(),seg),   Form("Mean dE %s%d;dE (MeV);Counts",CounterName[ihodo].Data(),seg),    400,    0, 40 );
      new TH1F( Form("%s%d_CTimeMean",CounterName[ihodo].Data(),seg),   Form("Mean CTime %s%d;Time (ns);Counts",CounterName[ihodo].Data(),seg),    4000,    -100, 100 );
      new TH1F( Form("%s%d_Position",CounterName[ihodo].Data(),seg),   Form("Hit position %s%d;Position (cm);Counts",CounterName[ihodo].Data(),seg),    500,    -50, 50 );
      new TH2F( Form("%s%d_dEMeanvCTimeMean",CounterName[ihodo].Data(),seg),   Form("dEMean CTimeMean corr. %s%d;dE (MeV);Time (ns)",CounterName[ihodo].Data(),seg),     400,    0, 40,  4000,    -100, 100 );
      new TH2F( Form("%su%d_dEvTOF",CounterName[ihodo].Data(),seg),   Form("dE TOF corr. %su%d;dE (MeV);TOF (ns)",CounterName[ihodo].Data(),seg),     400,    0, 40,  4000,    -100, 100 );
      new TH2F( Form("%sd%d_dEvTOF",CounterName[ihodo].Data(),seg),   Form("dE TOF corr. %sd%d;dE (MeV);TOF (ns)",CounterName[ihodo].Data(),seg),     400,    0, 40,  4000,    -100, 100 );
      new TH2F( Form("%s%d_dEMeanvTOF",CounterName[ihodo].Data(),seg),   Form("dEMean TOF corr. %s%d;dE (MeV);TOF (ns)",CounterName[ihodo].Data(),seg),     400,    0, 40,  4000,    -100, 100 );
    }
    new TH1F( Form("%s_Mult",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
  }
  // Chamber
  const int ndc=5;
  int dccid[ndc]={CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_BPC};
  TString dcname[ndc]={"BLC1a","BLC1b","BLC2a","BLC2b","BPC"};
  int NumOfLayers[ndc]={NumOfBLC1Layers,NumOfBLC1Layers,NumOfBLC2Layers,NumOfBLC2Layers,NumOfBPCLayers};

  for(int idc=0;idc<ndc;idc++){
    std::cout << "Define Histgram for " << dcname[idc] << std::endl;
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire( dccid[idc], layer );
      new TH1F( Form("%s_Layer%d_Mult",dcname[idc].Data(),layer), Form("Multiplicity %s Layer%d;Multiplicity;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_HitPat",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d;Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_TDC",dcname[idc].Data(),layer), Form("TDC %s Layer%d;TDC ch.;Counts",dcname[idc].Data(),layer), 4000, 0, 4000 );
      new TH1F( Form("%s_Layer%d_Time",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift time (ns);Counts",dcname[idc].Data(),layer), 3000, -500, 1000 );
      new TH1F( Form("%s_Layer%d_Length",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift length (mm);Counts",dcname[idc].Data(),layer), 1000, -1.0, 9.0 );
      new TH1F( Form("%s_Layer%d_Residual",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift length (mm);Counts",dcname[idc].Data(),layer), 1000, -1.0, 9.0 );
    }
  }
  // Cherenkov
  const int nchere = 1;
  TString ChereName[nchere] = {"AC"};
  for(int ichere=0; ichere<nchere; ichere++){
    std::cout << "Define Histograms for " << ChereName[ichere].Data() << std::endl;
    new TH1F( Form("%s_ADC",ChereName[ichere].Data()), Form("ADC Sum %s;ADC ch.;Counts",ChereName[ichere].Data()),4000, 0, 16000 );
    new TH1F( Form("%s_TDC",ChereName[ichere].Data()), Form("TDC %s;TDC ch,;Counts",ChereName[ichere].Data()),4000, 0, 4000 );
  }
  // TOF
  std::cout << "Define Histograms for TOF" << std::endl;
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    new TH2F( Form("T0%s_TOF",CounterName[ihodo].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo], 1, NumOfSegments[ihodo]+1, 4000, -50, 50 );
    for(int seg=1;seg<=5;seg++){
      new TH2F( Form("T0%d%s_TOF",seg,CounterName[ihodo].Data()), Form("TOF T0%d-%s;%s segment;TOF (ns)",seg,CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo], 1, NumOfSegments[ihodo]+1, 4000, -50, 50 );
    }
  }
  // Beam Momentum
  std::cout << "Define Histograms for Beam momentum" << std::endl;
  new TH1F( Form("BeamMom"), Form("Beam momentum;Beam momentum (GeV/c);Counts"),300, 0.0, 1.5 );
  new TH1F( Form("BeamMom_Chi2"), Form("Beam mom #chi^{2};#chi^{2};Counts"),101, 0, 101 );
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    new TH2F( Form("%sSegvBeamMom",CounterName[ihodo].Data()), Form("%s sement vs. Beam momentum;%s segment;Beam momentum (GeV/c)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo], 1, NumOfSegments[ihodo]+1, 300, 0.0, 1.5 );
  }
}


