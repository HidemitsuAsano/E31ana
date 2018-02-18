/* BeamLineChamber.h */
#ifndef BeamLineChamber_h
#define BeamLineChamber_h 1
#define BeamLineChamberCheck 5

#include "Common.h"

class BeamLineChamber : public TObject
{
 private:
  BeamLineChamber();

 public:
  BeamLineChamber(const int &cid, ConfMan *conf);

 private:
  int CID;

 public:
  std::string GetName();

  ClassDef( BeamLineChamber, 1 );
};

#endif
