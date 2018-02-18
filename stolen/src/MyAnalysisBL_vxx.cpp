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
  BeamMom = -9999.0; 
  BeamTOF = -9999.0;
  trackblc1 = 0;
  trackblc2 = 0;
  trackbpc  = 0;
  D5 -> Clear();
}

bool MyAnalysisBL::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{
  DetectorList *dlist=DetectorList::GetInstance();

  // ################################# //
  // T0 timing set //
  // ################################# //
  double T0Timing = -9999.; 
  int MulT0 = 0;
  for(int i=0; i<blMan->nT0(); i++){
    HodoscopeLikeHit* hit = blMan->T0(i);
    int seg = hit->seg();
    if(hit->CheckRange()){
      MulT0++;
      T0Timing = hit->ctmean();
    }
  }
  // ################ //
  // Chamber          //
  // ################ //
  const int ndc = 3;
  int dccid[ndc] = {CID_BLC1,CID_BLC2,CID_BPC};
  double ll[ndc] = {-20.0,-20.0,-20.0};
  double ul[ndc] = { 20.0, 20.0, 20.0};
  double chiul[ndc] = {30.0, 30.0, 30.0};
  int nhitdc[ndc] = { 0, 0, 0};
  int trackid[ndc] = { -1, -1, -1};
  for(int idc=0;idc<ndc;idc++){
    const int cid = dccid[idc];
    const char* name = dlist->GetName(cid).data();
    for(int i=0; i<bltrackMan->ntrackBLDC(cid); i++){
      LocalTrack* track = bltrackMan->trackBLDC(cid,i);
      double tracktime = track -> GetTrackTime();
      double chi2 = track -> chi2all();
        FillHist(Form("%s_TrackTime",name),tracktime);
        FillHist(Form("%s_Chi2",name),chi2);
      if(ll[idc]<=tracktime&&tracktime<=ul[idc] && chi2<=chiul[idc]){
        nhitdc[idc]++;        /* Num. of hit in cut region */
        trackid[idc] = i;     /* Selected track ID */
        FillHist(Form("%s_TrackTimewCut",name),tracktime);
        FillHist(Form("%s_Chi2wCut",name),chi2);
      }
    }
    FillHist(Form("%s_nTrack",name),nhitdc[idc]);
  }
  // ################ //
  // Hodoscope        //
  // ################ //
  const int nhodo=4;
  int hodoid[nhodo]={CID_BHD,CID_T0,CID_BPD,CID_DEF};
  int hodonHit[nhodo]; for(int i=0; i<nhodo; i++) hodonHit[i]=0;
  int hodoseg[nhodo][100];
  double hodoemean[nhodo][100];
  double hodoctmean[nhodo][100];
  double hodohitpos[nhodo][100];
  double hodotof[nhodo][100];
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
        FillHist(Form("%s_SegmentvsPosition",name),hodoseg[ihodo][ihit],hodohitpos[ihodo][ihit]);

        // === TOF calculation === //
        if(T0Timing>-900){
          hodotof[ihodo][i] = hodoctmean[ihodo][ihit]-T0Timing;
          if(ihodo==0){ hodotof[ihodo][i]=hodotof[ihodo][i]*-1.; }
          FillHist(Form("T0%s_TOF",name),hodoseg[ihodo][ihit],hodotof[ihodo][i]);
        }
        // === TOF calculation === //

        hodonHit[ihodo]++;
        ihit++;
      }
    }
    FillHist(Form("%s_Mult",name),hodonHit[ihodo]);
  }
  // ################ //
  // Beam analysis    //
  // ################ //
  double tofll = 24.0;
  double toful = 40.0;
  if(nhitdc[0]!=1||nhitdc[1]!=1||nhitdc[2]!=1) return false;
  trackblc1 = bltrackMan->trackBLC1(trackid[0]);
  trackblc2 = bltrackMan->trackBLC2(trackid[1]);
  trackbpc  = bltrackMan->trackBPC(trackid[2]);
  D5 -> TMinuitFit(trackblc1, trackblc2,conf);
  BeamMom = D5 -> mom();
  for(int i=0; i<blMan->nBHD(); i++){
    BeamTOF = hodotof[0][i];
    if(tofll<=BeamTOF&&BeamTOF<=toful) break;
  }
  FillHist(Form("BeamMomentum"), BeamMom);
  FillHist(Form("BeamTOF"), BeamTOF);
  FillHist(Form("Beam_Chi2"), D5->chisquare());
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    for(int i=0; i<hodonHit[ihodo]; i++){
      const char* name = dlist->GetName(hodoid[ihodo]).data();
      FillHist(Form("%sSegvsBeamMomentum",name),hodoseg[ihodo][i],BeamMom); 
    }
  }
  BeamP = trackbpc->GetMomDir();
  BeamP.SetMag(BeamMom);
  // === Beam PID === //
  double tofpi[2] = {24.0, 27.0};
  double tofk[2] = {27.0, 31.0};
  if(tofpi[0]<=BeamTOF&&BeamTOF<=tofpi[1]){
    BeamPID = Beam_Pion;
    BeamMass = piMass;
  }
  else if(tofk[0]<=BeamTOF&&BeamTOF<=tofk[1]){
    BeamPID = Beam_Kaon;
    BeamMass = kpMass;
  }
  else {
    BeamPID = Beam_Other;
    BeamMass = -1;
  }
  // === Beam PID === //
  BeamL.SetVectM( BeamP, BeamMass );

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
  std::cout << "Define Histograms for Hodoscope" << std::endl;
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    new TH1F( Form("%s_HitPat",CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
    new TH2F( Form("%s_SegmentvsPosition",CounterName[ihodo].Data()), Form("Hit Position %s;%s segment;Position (cm)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1, 201, -50.25, 50.25 );
    new TH1F( Form("%s_Mult",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
  }
  // Chamber
  const int ndc=3;
  int dccid[ndc]={CID_BLC1,CID_BLC2,CID_BPC};
  TString dcname[ndc]={"BLC1","BLC2","BPC"};
  std::cout << "Define Histgram for DC" << std::endl;
  for(int idc=0;idc<ndc;idc++){
    new TH1F( Form("%s_nTrack",dcname[idc].Data()), Form("Num. of tracks %s;Time (ns);Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
    new TH1F( Form("%s_TrackTime",dcname[idc].Data()), Form("Track time %s;Time (ns);Counts",dcname[idc].Data()), 600, -200.0, 400.0 );
    new TH1F( Form("%s_Chi2",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2};Counts",dcname[idc].Data()), 101, 0.0, 101.0 );
    new TH1F( Form("%s_TrackTimewCut",dcname[idc].Data()), Form("Track time %s;Time (ns);Counts",dcname[idc].Data()), 600, -200.0, 400.0 );
    new TH1F( Form("%s_Chi2wCut",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2};Counts",dcname[idc].Data()), 101, 0.0, 101.0 );
  }
  // TOF
  std::cout << "Define Histograms for TOF" << std::endl;
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    new TH2F( Form("T0%s_TOF",CounterName[ihodo].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo], 1, NumOfSegments[ihodo]+1, 4000, -50, 50 );
  }
  // Beam Momentum
  std::cout << "Define Histograms for Beam momentum" << std::endl;
  new TH1F( Form("BeamMomentum"), Form("Beam momentum;Beam momentum (GeV/c);Counts"),300, 0.0, 1.5 );
  new TH1F( Form("BeamTOF"), Form("Beam TOF;TOF between BHD and T0 (ns);Counts"),400, 20.0, 40.0 );
  new TH1F( Form("Beam_Chi2"), Form("Beam mom #chi^{2};#chi^{2};Counts"),101, 0, 101 );
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    new TH2F( Form("%sSegvsBeamMomentum",CounterName[ihodo].Data()), Form("%s sement vs. Beam momentum;%s segment;Beam momentum (GeV/c)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo], 1, NumOfSegments[ihodo]+1, 300, 0.0, 1.5 );
  }
}

