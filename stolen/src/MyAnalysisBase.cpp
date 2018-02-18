// MyAnalysisBase.cpp
// 2014.12.22
// T. Yamaga

#include "MyAnalysisBase.h"

MyAnalysisBase::MyAnalysisBase()
{
	confMan = 0;
  cdsTrackMan = 0;
  blTrackMan = 0;
  rtFile = 0;
}

void MyAnalysisBase::SetTrackMan(CDSTrackingMan* cds, BeamLineTrackMan* bl)
{
  cdsTrackMan = cds;
  blTrackMan = bl;
}

void MyAnalysisBase::SetTrackMan(CDSTrackingMan* cds)
{
  cdsTrackMan = cds;
}

void MyAnalysisBase::SetTrackMan(BeamLineTrackMan* bl)
{
  blTrackMan = bl;
}

void MyAnalysisBase::Initialize(ConfMan* conf)
{
  confMan = conf;
  rtFile = new TFile(confMan->GetOutFileName().c_str(),"RECREATE"); 
}

void MyAnalysisBase::Finalize()
{
  if(!rtFile) return;

  rtFile -> Write();
  rtFile -> Close();
  delete rtFile;
}

