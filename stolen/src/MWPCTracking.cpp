// MWPCTracking.cpp

#include "MWPCTracking.h"

MWPCTracking::MWPCTracking(ConfMan* c,BeamLineHitMan* bl)
{
	Initialize(c,bl);
	Clear();
}

MWPCTracking::~MWPCTracking()
{
	Clear();
}

bool MWPCTracking::Initialize(ConfMan* c, BeamLineHitMan* bl)
{
  conf = c;
  blMan = bl;
	dlist=DetectorList::GetInstance();

  return true;
}

void MWPCTracking::Clear()
{
	return;
}

void MWPCTracking::DoTracking()
{
}

void MWPCTracking::DoTracking()
{
	DetectorList *dlist=DetectorList::GetInstance();
}

