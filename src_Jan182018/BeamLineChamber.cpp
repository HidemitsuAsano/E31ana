/* BeamLineChamber.cpp */
#include "BeamLineChamber.h"

ClassImp(BeamLineChamber);

BeamLineChamber::BeamLineChamber()
{
}

BeamLineChamber::BeamLineChamber(const int &cid, ConfMan *conf) : CID(cid)
{


}

std::string BeamLineChamber::GetName()
{
  if( CID == CID_PDC1 ) return "PDC1";
  if( CID == CID_PDC2 ) return "PDC2";
  if( CID == CID_BLC1 ) return "BLC1";
  if( CID == CID_BLC2 ) return "BLC2";
  if( CID == CID_BPC  ) return "BPC ";
  return "unkown Beam Line Chamber !!!";
}
