#ifndef MyTools_h
#define MyTools_h 1

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
#include "TrackTools.h"

namespace MyTools
{
  int FindMass2(CDSTrack *track, LocalTrack *bpc, const double tof, const double beammom, const int pid_beam, double &beta_calc, double &mass2, double &tmptof,
		TVector3 &vtx_beam, TVector3 &vtx_cds, TVector3 &p, double &out_beam, double &beam_tof);

  bool IsTarget(const TVector3 &pos);

};
#endif
