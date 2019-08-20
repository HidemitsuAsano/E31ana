// TrackTools.h
#ifndef TrackTools_h
#define TrackTools_h 1

#include <string>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TMath.h"
#include "GlobalVariables.h"
#include "GeomTools.h"
#include "MathTools.h"
#include "CDSTrack.h"
#include "CDSTrackingMan.h"
#include "Particle.h"
#include "LocalTrack.h"

namespace TrackTools
{
  bool CalcBetaMass2TOF(CDSTrack *track, LocalTrack* bpc,double tof,double mass,double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof);
  bool FindMass2(CDSTrack *track, LocalTrack* bpc,double tof, double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof);
  bool FindMass2C(CDSTrack *track, LocalTrack* bpc,double tof, double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof);
  bool Calc2HelixVertex(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2);

  bool Calc2HelixVertexWresl(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2);


  bool Calc2HelixVertex2(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2);
  bool CalcLineHelixVertex(LocalTrack *track,CDSTrack *cds, TVector3 &vtx1, TVector3 &vtx2, double &dis);
  bool CalcLineHelixVertex(pBeam *beam,CDSTrack *cds,TVector3 &vtx1,TVector3 &vtx2,double &dis, bool ELOSS=false);
  void CalcBetaMass(TVector3 vertex,LocalTrack *beam,CDSTrack* cdc, 
		     ConfMan *conf,int pid_beam,double tof, double &beta,double &mass2);
  void CalcBetaMass2(TVector3 vertex,LocalTrack *beam,CDSTrack* cdc, 
		     ConfMan *conf,int pid_beam,double beammom,double tof, double &beta,double &mass2);

  int PID(double mom,double mass2);
  int PIDcorr(double mom,double mass2);
  int PIDcorr3(double mom,double mass2);
  int PIDcorr_wide(double mom, double mass2);

  pCDS *CalcSingleAll(pBeam *beam,CDSTrack* cdc,CDSHitMan *cdsMan, bool ELOSS=false);
  pCDS *Calc2HelixAll(pBeam *beam, pCDS* cds1, pCDS* cds2,CDSTrackingMan *trackMan, bool ELOSS=false);
};

#endif
