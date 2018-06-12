#ifndef HistManwMC_h
#define HistManwMC_h 1

#include <iostream>
#include <string>

#include "TTree.h"
#include "TNtuple.h"

#include "KnuclRootData.h"
#include "BeamLineHitMan.h"
#include "CDSHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"
#include "TrackTools.h"
#include "ELossTools.h"
#include "MyTools.h"
#include "EventHeader.h"
#include "SimDataReader.h"
#include "NpipiData.h"
#include "MyParam.h"

#include "KinFitMan.h"
#include "HistTools.h"
#include "HistTools.h"

//#define CALIB 1

class HistManwMC
{
 public:
  HistManwMC(TFile *f, ConfMan *conf);
  virtual ~HistManwMC(){};

 private:
  // Output TFile pointer;
  TFile *rtFile;
  ConfMan *confMan;
  // Data Manager for knucl
  SimDataReader *simReader;
  // Data Manager for k18ana
  BeamLineHitMan *blMan;
  CDSHitMan *cdsMan;
  BeamLineTrackMan *bltrackMan;
  CDSTrackingMan *cdstrackMan;
  // for Analysis Data
  int fBeamPID;
  bool fBeamPion;
  bool fBeamKaon;
  double fT0time;
  double fD5mom;
  LocalTrack *fTrackBPC;
  TLorentzVector fBeamLmom;
  TVector3 fT0pos;

  std::vector<HodoscopeLikeHit*> fT0_hit;
  std::vector<HodoscopeLikeHit*> fBHD_hit;
  std::vector<HodoscopeLikeHit*> fBPD_hit;
  std::vector<HodoscopeLikeHit*> fCDH_hit;
  std::vector<HodoscopeLikeHit*> fBVC_hit;
  std::vector<HodoscopeLikeHit*> fCVC_hit;
  std::vector<HodoscopeLikeHit*> fPC_hit;
  std::vector<HodoscopeLikeHit*> fNC_hit[8];
  std::vector<HodoscopeLikeHit*> fBD_hit;

  //********************//
  //*** for CDS Data ***//
  //********************//
  std::vector<int> fCDSPID;
  std::vector<double> fCDSbeta;
  std::vector<double> fCDSmass2;
  std::vector<double> fCDSmom;
  std::vector<CDSTrack*> fTrackPim;
  std::vector<CDSTrack*> fTrackKm;
  std::vector<CDSTrack*> fTrackPip;
  std::vector<CDSTrack*> fTrackP;
  std::vector<CDSTrack*> fTrackD;
  TVector3 fVtxCDS; // Best Vertex Helix
  TVector3 fVtxBeam; // Best Vertex Beam

  //*******************//
  //*** for NC Data ***//
  //*******************//
  int fFPID;
  HodoscopeLikeHit *fNC_eff_hit;
  double fNCbeta;
  TLorentzVector fFLmom;
  int fNCseg;
  double fNCtime;
  double fNCdE;
  TVector3 fNCpos;

  // for Forward Charge
  bool fFCflag;
  HodoscopeLikeHit *fFChit;
  LocalTrack *fFDC1track;
  TVector3 fFC_start;
  TVector3 fFC_FDC1pos;
  TVector3 fFC_hitpos;
  double fFC_Angle;
  double fFC_mass2;
  double fFC_beta;
  double fFC_Mom_USWK;
  double fFC_Mom_TOF;

  //***************************//
  //*** for n pi+ pi- event ***//
  //***************************//
  //  NpipiData *fNpipiData;
  // for IM
  bool fKNpipi_K0_flag;
  bool fKNpipi_K0_SB0_flag;
  bool fKNpipi_K0_SB1_flag;
  bool fKNpipi_K0_SB2_flag;
  bool fKNpipi_K0_SB3_flag;
  bool fKNpipi_K0_SB4_flag;
  bool fKNpipi_K0_SB5_flag;

  bool fKNpipi_Sm_SB0_flag;
  bool fKNpipi_Sm_SB1_flag;
  bool fKNpipi_Sm_SB2_flag;
  bool fKNpipi_Sm_SB3_flag;
  bool fKNpipi_Sm_SB4_flag;
  bool fKNpipi_Sm_SB5_flag;

  bool fKNpipi_Sp_SB0_flag;
  bool fKNpipi_Sp_SB1_flag;
  bool fKNpipi_Sp_SB2_flag;
  bool fKNpipi_Sp_SB3_flag;
  bool fKNpipi_Sp_SB4_flag;
  bool fKNpipi_Sp_SB5_flag;

  bool fKNpipi_Sp_flag;
  bool fKNpipi_Sm_flag;
  bool fKNpipi_K0_flag25;
  bool fKNpipi_Sp_flag25;
  bool fKNpipi_Sm_flag25;
  bool fKNpipi_K0_flag3;
  bool fKNpipi_Sp_flag3;
  bool fKNpipi_Sm_flag3;
  bool fKNpipi_K0_flag4;
  bool fKNpipi_Sp_flag4;
  bool fKNpipi_Sm_flag4;
  // for MM
  bool fKNpipi_N_flag;
  bool fKNpipi_N0_flag;
  bool fKNpipi_N1_flag;
  bool fKNpipi_N2_flag;
  bool fKNpim_Sp_flag;
  bool fKNpip_Sm_flag;
  // "n" pi-> S
  bool fSm_mass_rc_flag;
  bool fSp_mass_rc_flag;
  bool fSm_mass_rc_flag2;
  bool fSp_mass_rc_flag2;
  bool fSm_mass_rc_flag3;
  bool fSp_mass_rc_flag3;
  bool fSm_mass_rc_flag4;
  bool fSp_mass_rc_flag4;
  bool fSm_mass_rc_flag5;
  bool fSp_mass_rc_flag5;

 public:
  void setSimReader(SimDataReader *reader){ simReader=reader; };
  void setK18ana(BeamLineHitMan *blman, CDSHitMan *cdsman, BeamLineTrackMan *bltrackman, CDSTrackingMan *cdstrackman)
  {
    blMan=blman, cdsMan=cdsman, bltrackMan=bltrackman, cdstrackMan=cdstrackman;
  }

  void setBeamPID(const int &pid){ fBeamPID=pid; };
  void setBeamKaon(const bool &flag=true){ fBeamKaon=flag; };
  void setBeamPion(const bool &flag=true){ fBeamKaon=flag; };
  void setT0time(const double &time){ fT0time=time; };
  void setD5mom(const double &mom){ fD5mom=mom; };
  void setTrackBPC(LocalTrack *track){ fTrackBPC=track; };

  int beamPID() const { return fBeamPID; };		 
  double T0time() const { return fT0time; };
  double D5mom() const { return fD5mom; };
  LocalTrack* trakcBPC() const { return fTrackBPC; };

  int nT0() const { return fT0_hit.size(); };
  int nBPD() const { return fT0_hit.size(); };
  int nCDH() const { return fT0_hit.size(); };
  int nBVC() const { return fT0_hit.size(); };
  int nCVC() const { return fT0_hit.size(); };
  int nPC() const { return fT0_hit.size(); };
  int nNC() const {
    int nNC = 0;
    for( int i=0; i<8; i++ ) nNC += fNC_hit[i].size();
    return nNC;
  }
  int nNC(const int &lay){ return fNC_hit[lay-1].size(); };
  int nBD() const { return fBD_hit.size(); };

  HodoscopeLikeHit *T0(const int &i){ return fT0_hit[i]; };
  HodoscopeLikeHit *BPD(const int &i){ return fBPD_hit[i]; };
  HodoscopeLikeHit *CDH(const int &i){ return fCDH_hit[i]; };
  HodoscopeLikeHit *BVC(const int &i){ return fBVC_hit[i]; };
  HodoscopeLikeHit *CVC(const int &i){ return fCVC_hit[i]; };
  HodoscopeLikeHit *PC(const int &i){ return fPC_hit[i]; };
  HodoscopeLikeHit *NC(const int &lay, const int &i){ return fNC_hit[lay-1][i]; };
  HodoscopeLikeHit *BD(const int &i){ return fBD_hit[i]; };

  void initHistCom();
  void initHist();
  void initCalib();
  void fillCalib(EventHeader* header);
  //  void anaMC_NC();
  bool ana(EventHeader *header=0); // return false if not find CDS GoodTrack w/ CDH hit
  void fill(EventHeader *header=0);
  void fill_pppim(EventHeader *header=0);
  void finit();

  bool Clear(const bool flag=true);

 public:
  void anaNC(EventHeader *header);
  void anaFC(EventHeader *header);
  void anaNC_MC();

  void printHit();
};
#endif
