#ifndef SimDataMan_h
#define SimDataMan_h 1

#include <vector>
#include <string>

#include "TObject.h"
#include "TTree.h"

#include "ConfMan.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "HodoscopeLikeHit.h"
#include "ChamberLikeHit.h"
#include "CDCHit.h"

#include "ComCrossSectionTable.h"
#include "KnuclRootData.h"

class SimDataMan : public TObject
{
 public:
  SimDataMan();
  virtual ~SimDataMan(){};

  void Convert(DetectorData*, ConfMan*, BeamLineHitMan*, CDSHitMan*);

  ClassDef(SimDataMan, 1);
};

#endif

