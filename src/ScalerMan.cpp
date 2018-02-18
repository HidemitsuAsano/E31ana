#include "ScalerMan.h"

ClassImp(Scaler);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
Scaler::Scaler() : TObject()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void Scaler::Clear()
{
  Name="";
  Value=-1;
}

ClassImp(ScalerMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
ScalerMan::ScalerMan() : TObject()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void ScalerMan::AddHit( const int &val, const std::string &name )
{
  Scaler *sca = new Scaler();
  sca->SetName( name );
  sca->SetValue(val);
  ScaContainer.push_back(*sca);
  delete sca;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void ScalerMan::Clear()
{
  ScaContainer.clear();
}
