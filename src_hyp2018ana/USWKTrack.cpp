#include "USWKTrack.h"

USWKTrack::USWKTrack()
  : ParticleName("unknown"), Mass(-999), Charged(-999),
    InTime(-999), InPosition(TVector3(-999,-999,-999)), InMomentum(TVector3(-999,-999,-999)),
    OutTime(-999), OutPosition(TVector3(-999,-999,-999)), OutMomentum(TVector3(-999,-999,-999))
{
}

USWKTrack::USWKTrack(const USWKTrack &right)
  : ParticleName(right.ParticleName), Mass(right.Mass), Charged(right.Charged),
    InTime(right.InTime), InPosition(right.InPosition), InMomentum(right.InMomentum),
    Time(right.Time), Position(right.Position), Momentum(right.Momentum),
    OutTime(right.OutTime), OutPosition(right.OutPosition), OutMomentum(right.OutMomentum)
{
}

void USWKTrack::SetParticleName(std::string &name)
{
  ParticleName = name;
  if( name=="proton" ){
    Mass = pMass;
    Charged = 1.0;
  }
  else if( name=="deuteron" ){
    Mass = dMass;
    Charged = 1.0;
  }
  else if( name=="pi+" ){
    Mass = piMass;
    Charged = 1.0;
  }
  else{
    std::cout<<" USWKTrack set unknown particle ! "<<name<<std::endl;
  }
}

void USWKTrack::Clear()
{
  InTime = -999;
  InPosition.SetXYZ(-999, -999, -999);
  InMomentum.SetXYZ(-999, -999, -999);

  Time.clear();
  Position.clear();
  Momentum.clear();

  OutTime = -999;
  OutPosition.SetXYZ(-999, -999, -999);
  OutMomentum.SetXYZ(-999, -999, -999);
}
