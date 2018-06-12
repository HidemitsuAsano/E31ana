#ifndef DCEffMan_h
#define DCEffMan_h 1

#include <vector>
#include <iostream>

#include "BeamLineTrackMan.h"

class DCEffMan
{
 public:
  DCEffMan(const int &cid, ConfMan *confMan, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan);
  virtual ~DCEffMan(){};

 private:
  ConfMan *fConfMan;
  BeamLineHitMan *fBLMan;
  BeamLineTrackMan *fBLTrackMan;
  int fCID;
  std::vector<std::vector<ChamberLikeHit*> > fHits;

 public:
  void get();
  void dump();
  bool trig18();
  bool trig18_2();
  bool trig(const int &lay, const int &max=99)const;
  bool eff(const int &lay, const int &max=3)const;
  bool all1hit();
  LocalTrack* getTrack();
  bool trigTrack(const int &lay);
  bool effTrack(const int &lay);

 private:
  int maxWire() const;
  void clear();
};
#endif
