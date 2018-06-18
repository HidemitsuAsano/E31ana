#include "CDS2Data.h"

CDS2Data::CDS2Data() : fID1(-1), fID2(-1), fPID1(CDS_Other), fPID2(CDS_Other), fDCA(DBL_MIN)
{
}

CDS2Data::CDS2Data(const int &id1, const int &id2) : fID1(id1), fID2(id2), fPID1(CDS_Other), fPID2(CDS_Other), fDCA(DBL_MIN)
{
}


int CDS2Data::set(const int &beam_pid, const double &beam_mom, LocalTrack *beam, CDSTrack *cds1, CDSTrack *cds2)
{
  fBeamPID = beam_pid;
  TVector3 T0pos = beam->GetPosatZ(-110.5);
  fPID1 = cds1-> PID();
  fPID2 = cds2-> PID();

  if( !TrackTools::Calc2HelixVertex(cds1, cds2, fVtxCDS1, fVtxCDS2) ) return 1;

  if( !cds1->GetMomentum(fVtxCDS1, fCDSMom1, true, true) ) return 2;
  if( !cds2->GetMomentum(fVtxCDS2, fCDSMom2, true, true) ) return 3;

  TVector3 vtx_mean = 0.5*(fVtxCDS1+ fVtxCDS2);
  TVector3 mom_sum = fCDSMom1+fCDSMom2;

  double dltmp, dist;
  TVector3 vtx;
  MathTools::LineToLine(vtx_mean, mom_sum.Unit(), beam->GetPosatZ(-20), beam->GetMomDir(), dltmp, fDCA, fVtxCDS, fVtxBeam);

  double beam_out, beam_tof;
  ELossTools::CalcElossBeamTGeo(T0pos, fVtxBeam, beam_mom, parMass[fBeamPID], beam_out, beam_tof);
  fBeamMom = beam->GetMomDir();
  fBeamMom.SetMag(beam_out);

  return 0;
}
