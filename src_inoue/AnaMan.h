#ifndef AnaMan_h
#define AaaMan_h 1

#include <iostream>
#include <vector>

#include "ConfMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "BeamSpectrometer.h"
#include "HodoscopeLikeHit.h"
#include "LocalTrack.h"

#include "CDS1Data.h"
#include "CDS2Data.h"

class AnaMan
{
 public:
  AnaMan(ConfMan *conf);
  ~AnaMan(){};

 private:
  BeamSpectrometer *fBeamSpec;

  int fStatus; // This is event reduction status
               // 1 : T0 1hit
               // 2 : TOF Kaon
               // 3 : BLC1 1Track w/ Time Window & chi2 cut
               // 4 : BLC2 1Track w/ Time Window & chi2 cut
               // 5 : Beam Momentum reconstruction (BLC1 1track & BLC2 1track required)
               // 6 : BPC 1track
               // 7 : CDS 1track

  int fBeamPID;
  double fT0time;
  TVector3 fT0pos;
  double fD5mom;
  std::vector<HodoscopeLikeHit*> fBHD_hit;
  std::vector<HodoscopeLikeHit*> fT0_hit;
  std::vector<HodoscopeLikeHit*> fBPD_hit;
  std::vector<HodoscopeLikeHit*> fDEF_hit;
  std::vector<HodoscopeLikeHit*> fCDH_hit;
  std::vector<HodoscopeLikeHit*> fBVC_hit;
  std::vector<HodoscopeLikeHit*> fCVC_hit;
  std::vector<HodoscopeLikeHit*> fPC_hit;
  std::vector<HodoscopeLikeHit*> fNC_hit[8];
  std::vector<HodoscopeLikeHit*> fBD_hit;
  std::vector<HodoscopeLikeHit*> fLB_hit;
  std::vector<HodoscopeLikeHit*> fWVC_hit;

  std::vector<LocalTrack*> fTrackBLC1;
  std::vector<LocalTrack*> fTrackBLC2;
  std::vector<LocalTrack*> fTrackBPC;

  //*** -50, 200 ***//
  std::vector<LocalTrack*> fTrackBLC1_1;
  std::vector<LocalTrack*> fTrackBLC2_1;
  std::vector<LocalTrack*> fTrackBPC_1;

  //*** -30, 100 ***//
  std::vector<LocalTrack*> fTrackBLC1_2;
  std::vector<LocalTrack*> fTrackBLC2_2;
  std::vector<LocalTrack*> fTrackBPC_2;

  int fVtxID;
  std::vector<CDS1Data> fCDS1;
  std::vector<CDS2Data> fCDS2;

  int fFPID;
  TVector3 fFMom;
  TVector3 fFHitPos;
  double fFdE;
  double fFBeta;
  double fFtime;

 public:
  int status() const { return fStatus; };
  int beam_pid() const { return fBeamPID; };
  double D5mom() const { return fD5mom; };

  int nBHD() const { return fBHD_hit.size(); };
  int nT0()  const { return fT0_hit.size(); };
  int nBPD() const { return fBPD_hit.size(); };
  int nDEF() const { return fDEF_hit.size(); };
  int nCDH() const { return fCDH_hit.size(); };
  int nBVC() const { return fBVC_hit.size(); };
  int nCVC() const { return fCVC_hit.size(); };
  int nPC()  const { return fPC_hit.size(); };
  int nNC()  const { 
    int n=0;
    for( int i=0; i<8; i++ ){
      n += fNC_hit[i].size();
    }
    return n;
  }
  int nNC(const int &lay) const { return fNC_hit[lay-1].size(); }; // lay is 1 initial;
  int nBD() const { return fBD_hit.size(); };
  int nLB() const { return fLB_hit.size(); };
  int nWVC() const { return fWVC_hit.size(); };

  HodoscopeLikeHit* BHD(const int &i){ return fBHD_hit[i]; };
  HodoscopeLikeHit* T0(const int &i) { return fT0_hit[i]; };
  HodoscopeLikeHit* BPD(const int &i){ return fBPD_hit[i]; };
  HodoscopeLikeHit* DEF(const int &i){ return fDEF_hit[i]; };
  HodoscopeLikeHit* CDH(const int &i){ return fCDH_hit[i]; };
  HodoscopeLikeHit* BVC(const int &i){ return fBVC_hit[i]; };
  HodoscopeLikeHit* CVC(const int &i){ return fCVC_hit[i]; };
  HodoscopeLikeHit* PC(const int &i) { return fPC_hit[i]; };
  HodoscopeLikeHit* NC(const int &lay, const int &i){ return fNC_hit[lay-1][i]; };
  HodoscopeLikeHit* NC(const int &i){
    int n=0;
    for( int lay=0; lay<8; lay++ ){
      if( i<n+fNC_hit[lay].size() ) return fNC_hit[lay][i-n];
      n += fNC_hit[i].size();
    }
    return 0;
  }
  HodoscopeLikeHit* BD(const int &i) { return fBD_hit[i]; };
  HodoscopeLikeHit* LB(const int &i) { return fLB_hit[i]; };
  HodoscopeLikeHit* WVC(const int &i) { return fWVC_hit[i]; };

  int ntrackBLC1(){ return fTrackBLC1.size(); };
  int ntrackBLC2(){ return fTrackBLC2.size(); };
  int ntrackBPC(){ return fTrackBPC.size(); };
  LocalTrack* trackBLC1(const int &i){ return fTrackBLC1[i]; };
  LocalTrack* trackBLC2(const int &i){ return fTrackBLC2[i]; };
  LocalTrack* trackBPC(const int &i){ return fTrackBPC[i]; };

  int ntrackBLC1_1(){ return fTrackBLC1_1.size(); };
  int ntrackBLC2_1(){ return fTrackBLC2_1.size(); };
  int ntrackBPC_1(){ return fTrackBPC_1.size(); };
  LocalTrack* trackBLC1_1(const int &i){ return fTrackBLC1_1[i]; };
  LocalTrack* trackBLC2_1(const int &i){ return fTrackBLC2_1[i]; };
  LocalTrack* trackBPC_1(const int &i){ return fTrackBPC_1[i]; };

  int ntrackBLC1_2(){ return fTrackBLC1_2.size(); };
  int ntrackBLC2_2(){ return fTrackBLC2_2.size(); };
  int ntrackBPC_2(){ return fTrackBPC_2.size(); };
  LocalTrack* trackBLC1_2(const int &i){ return fTrackBLC1_2[i]; };
  LocalTrack* trackBLC2_2(const int &i){ return fTrackBLC2_2[i]; };
  LocalTrack* trackBPC_2(const int &i){ return fTrackBPC_2[i]; };

  void set(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, CDSHitMan *cdsMan);
  void set(ConfMan *conf, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackingMan);

  void clear();
};
#endif
