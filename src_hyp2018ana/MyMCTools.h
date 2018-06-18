#ifndef MYHISTTOOLS_HH
#define MYHISTTOOLS_HH 1

#include "GlobalVariables.h"
#include "AnaInfo.h"
#include "KnuclRootData.h"

#include "MyParam.h"

namespace MyMCTools{
  bool effNC(DetectorData *detData);

  std::vector<DetectorHit*> NChits(DetectorData *detData);
  bool goodBeam(MCData *mcData, DetectorData *detData);
  bool goodPi_beam(MCData *mcData, DetectorData *detData);
  bool goodP_beam(MCData *mcData, DetectorData *detData);
  Track *initTrack(MCData *mcData, int pdg);
  Track *track(MCData *mcData, int id);
  Track *trackByPDG(MCData *mcData, int pdg);
  double lmom(int pdg, ReactionData *reacData);
  std::vector<Track*> getInitTracks(MCData *mcData);
  std::vector<Track*> getDaughter(Track *parent, MCData *mcData);
  int nHit(int cid, DetectorData *detData);
  bool isHit(int cid, Track *track, DetectorData *detData);

  bool effFC(DetectorData *detData, MCData *mcData, BeamLineHitMan *blMan);

  DetectorHit *getHit(DetectorData *detData, HodoscopeLikeHit *hit);

  double Ystar_mass(ReactionData *reacData);
  double S1385mass(ReactionData *reacData);
  bool pimS0all_hit(MCData *mcData, DetectorData *detData);
  CDSInfo *getCDSInfo(DetectorHit *hit, AnaInfo *anaInfo, CDSHitMan *cdsMan);

};
#endif
