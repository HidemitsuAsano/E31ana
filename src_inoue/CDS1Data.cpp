#include "CDS1Data.h"

CDS1Data::CDS1Data() : fID(-1), fPID(CDS_Other), fBeamPID(Beam_Other), fCDHseg(-1), fTOF(DBL_MIN), fMass2(DBL_MIN), fMom(0), fBeta(DBL_MIN), fdE(0.0)
{
}

CDS1Data::CDS1Data(const int &id) : fID(id), fPID(CDS_Other), fBeamPID(Beam_Other), fCDHseg(-1), fTOF(DBL_MIN), fMass2(DBL_MIN), fMom(0), fBeta(DBL_MIN), fdE(0.0)
{
}

int CDS1Data::set(const int &beam_pid, const double &t0time, const double &beam_mom, LocalTrack *beam, CDSTrack *cdc, CDSHitMan *cdsMan){
  fBeamPID = beam_pid;
  TVector3 T0pos = beam-> GetPosatZ(-110.5);
  double CDHtime, dis;
  if( !cdc-> GetCDHHit(cdsMan, fCDHseg, CDHtime) ) return 1;
  for( int i=0; i<cdc->nCDHHit(); i++ ){
    fdE += cdc-> CDHHit(cdsMan, i)->emean();
  }
  double par[5];
  cdc-> GetParameters(CID_CDC, par, fVtxCDS);
  fMom = cdc-> Momentum();
  TrackTools::CalcLineHelixVertex(beam, cdc, fVtxBeam, fVtxCDS, dis);
  double beam_out, beam_tof;
  ELossTools::CalcElossBeamTGeo(T0pos, fVtxBeam, beam_mom, parMass[fBeamPID], beam_out, beam_tof);

  TVector3 vtxCDH = cdc-> CDHVertex();
  double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, fVtxCDS);
  fBeta = cdc_dis/(CDHtime-t0time)/(100.*Const);
  fMass2 = fMom*fMom*(1./(fBeta*fBeta)-1);

  double beta_calc, tofvtxcdc;
  if( !TrackTools::FindMass2C(cdc, beam, beam_tof, beam_mom, parMass[fBeamPID], fBeta, fMass2, tofvtxcdc) ) return 2;
  fPID = TrackTools::PID(fMom, fMass2);
  cdc-> SetPID(fPID);

  double tmpl;
  TVector3 vtxb, vtxcds;
  if( !cdc-> CalcVertexTimeLength(T0pos, beam->GetMomDir(), cdsMass[fPID], vtxb, vtxcds, fTOF, tmpl, true) ) return 3;
  //  if( !cdc-> CalcVertexTimeLength(T0pos, beam->GetMomDir(), cdsMass[fPID], fVtxBeam, fVtxCDS, fTOF, tmpl) ) return 3;

  if( !cdc-> GetMomentum(fVtxCDS, fCDSMom, true, true) ) return 4;

  fBeamMom = beam->GetMomDir();
  fBeamMom.SetMag(beam_out);

  return 0;
}

