//v2//
#ifndef MYTOOLS_H
#define MYTOOLS_H 1

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
//#include "KnuclRootData.h"
#include "MyParam.h"

//#include "AnaInfo.h"
#include "EventHeader.h"

namespace MyTools
{
  int PIDcorr_wide(double mom, double mass2);

};
#endif

