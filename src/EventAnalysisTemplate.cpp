#include "EventAnalysisTemplate.h"

using namespace std;

// Please return true  if UAna() return false, program stop.
bool EventAnalysis::UAna()
{
  // Please write user program w/o CDSTrackingMan
  int nT0=0;
  int segT0=-1;
  double ctmT0=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
      ctmT0 = blMan->T0(i)->ctmean();
      segT0=blMan->T0(i)->seg();
    }
  }
  if( nT0!=1 ) return true;

  int nBHD=0;
  double ctmBHD=0;
  int pid_beam=2;//0:pi 1:K 2:else
  double tofBHDT0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      nBHD++;
      ctmBHD = blMan->BHD(i)->ctmean();
      tofBHDT0 = ctmT0-ctmBHD;
      if(header->kaon()&&tofBHDT0>27.5&&tofBHDT0<31.) pid_beam=Beam_Kaon;
      else if(header->pion()&&tofBHDT0<27.5&&tofBHDT0>24.5) pid_beam=Beam_Pion;
    }
  }
  if(pid_beam!=Beam_Kaon)  return true;

  int nblc1=0;
  int blc1id=-1;
  for(int i=0;i<bltrackMan->ntrackBLC1();i++)
    if(bltrackMan->trackBLC1(i)->CheckRange(-30,100)){
      nblc1++;
      blc1id=i;
    }
  if(nblc1!=1) return true;
  LocalTrack *blc1track=bltrackMan->trackBLC1(blc1id);

  int nblc2=0;
  int blc2id=-1;
  for(int i=0;i<bltrackMan->ntrackBLC2();i++)
    if(bltrackMan->trackBLC2(i)->CheckRange(-30,100)){
      nblc2++;
      blc2id=i;
    }
  if(nblc2!=1) return true;
  LocalTrack *blc2track=bltrackMan->trackBLC2(blc2id);

  beamSpec->TMinuitFit(blc1track, blc2track, confMan);
  double D5mom=beamSpec->mom();

  int nbpc=0;
  int bpcid=-1;
  for(int i=0;i<bltrackMan->ntrackBPC();i++)
    if(bltrackMan->trackBPC(i)->CheckRange(-30,100)){
      nbpc++;
      bpcid=i;
    }
  if(nbpc!=1) return true;
  LocalTrack *bpctrack=bltrackMan->trackBPC(bpcid);

  if( cdsFileMatching() ){
    // Please wite user program using CDSTrackingMan
    int n_pim = 0;
    int n_pip = 0;
    CDSTrack *pimTrack = 0;
    CDSTrack *pipTrack = 0;

    for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
      CDSTrack *track = cdstrackMan->GoodTrack(i);
      if( track->Chi()>30 ) continue;
      if( !track->CDHFlag() ) continue;
      int CDHseg;
      double CDHtime;
      track->GetCDHHit(cdsMan, CDHseg, CDHtime);
    
      TVector3 vtxBeam, vtxCDS;
      if( !track->GetVertex(bpctrack->GetPosatZ(110-0.5), bpctrack->GetMomDir(),
                            vtxBeam, vtxCDS) ) continue;
      double beam_out, beam_tof;

      ELossTools::CalcElossBeamTGeo(bpctrack->GetPosatZ(110-0.5), vtxBeam, D5mom,
                                    kpMass, beam_out, beam_tof);

      double mom = track->Momentum();
      TVector3 vtxCDH = track->CDHVertex();
      double par[5];
      track->GetParameters(par);
      double tof = CDHtime-ctmT0;
      double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, vtxCDS);
      double beta = cdc_dis/(tof-beam_tof)/(100.*Const);
      double mass2 = mom*mom*(1./(beta*beta)-1);
      double tofvtxcdc;
    
      if( !TrackTools::FindMass2(track, bpctrack, tof, D5mom, Beam_Kaon, beta,
                                 mass2, tofvtxcdc) ) continue;
    
      track->Retiming(cdsMan, confMan, beta, true); 
      track->Calc(confMan); 
      if( !TrackTools::FindMass2(track, bpctrack, tof, D5mom, Beam_Kaon, beta, mass2, tofvtxcdc) ) continue;

      int pid = confMan->GetCDSFittingParamManager()->PID(mom, mass2);

      TH2F *h21 = (TH2F*)gFile->Get("CDS_PID");
      TH2F *h22 = (TH2F*)gFile->Get("CDS_PID_pi");
      h21->Fill(mass2, mom);
      if( pid==CDS_PiMinus || pid==CDS_PiPlus ) h22->Fill(mass2, mom);

      double tmpfl;
      track-> CalcVertexTimeLength(bpctrack->GetPosatZ(110-0.5), bpctrack->GetMomDir(),
                                   cdsMass[pid], vtxBeam, vtxCDS, tof, tmpfl, true);

      if( pid==CDS_PiMinus ){
        n_pim++;
        pimTrack = track;
      }
      if( pid==CDS_PiPlus ){
        n_pip++;
        pipTrack = track;
      }
    }

    if( n_pim==1 && n_pip==1 ){
      TVector3 vtx_pim, vtx_pip;
      if( !TrackTools::Calc2HelixVertex(pimTrack, pipTrack, vtx_pim, vtx_pip) ) return true;
      TVector3 p_pim, p_pip;
      if( !pimTrack->GetMomentum(vtx_pim, p_pim, true, true) ) return true;
      if( !pipTrack->GetMomentum(vtx_pip, p_pip, true, true) ) return true;

      TLorentzVector pim_lmom, pip_lmom;
      pim_lmom.SetVectM(p_pim, piMass);
      pip_lmom.SetVectM(p_pip, piMass);

      TH1F *h1 = (TH1F*)gFile->Get("CDS_IM_pipi");
      h1-> Fill((pim_lmom+pip_lmom).M());
    }
  }
  return true;
}

void EventAnalysis::InitializeHistogram()
{
  new TH1F("CDS_IM_pipi", "CDS IM #pi^{+} #pi^{-}", 200, 0.4, 0.6);
  new TH2F("CDS_PID",     "CDS PID",     100, -0.5, 2.0, 100, -1.5, 1.5);
  new TH2F("CDS_PID_pi",  "CDS PID #pi", 100, -0.5, 2.0, 100, -1.5, 1.5);
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
