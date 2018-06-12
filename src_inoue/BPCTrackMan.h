#ifndef BPCTrackMan_h
#define BPCTrackMan_h 1

#include "ConfMan.h"
#include "BeamLineTrackMan.h"
#include "BLDCTrack.h"
#include "XYTrack.h"

class BPCTrackMan
{
 public:
  BPCTrackMan(ConfMan *conf, BeamLineHitMan *blMan);
  virtual ~BPCTrackMan(){};

 private:
  ConfMan *fConfMan;
  BeamLineHitMan *fBLMan;

  std::vector<ChamberLikeHit*> fXZHits[4];
  std::vector<ChamberLikeHit*> fYZHits[4];

  std::vector<XYTrack> fXTracks;
  std::vector<XYTrack> fYTracks;

  void set();
  int nXCluster(const int &i);
  int nYCluster(const int &i);
  bool setXCluster(const int &i1, const int &i2, XYTrack &track);
  bool setYCluster(const int &i1, const int &i2, XYTrack &track);
  void deleteHit(XYTrack &track);

 public:
  void DoTracking();
  int nXTrack() const { return fXTracks.size(); };
  int nYTrack() const { return fYTracks.size(); };
  XYTrack* xTrack(const int &i){ return &fXTracks[i]; };
  XYTrack* yTrack(const int &i){ return &fYTracks[i]; };

  void dumpHits();
  void dumpTrack();
};
#endif
