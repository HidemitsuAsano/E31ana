//AnalysisMan2.h
#ifndef AnalysisMan2_h
#define AnalysisMan2_h 1

#include "GlobalVariables.h"
#include "ConfMan.h"
#include "ScalerMan.h"
#include "EventHeader.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "CDSTrackingMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "UshiwakaFieldMapMan.h"
#include "UshiwakaTableMan.h"

#include "InvariantMass.h"
#include "MissingMass.h"

class AnalysisMan
{
 public:
  AnalysisMan();
  ~AnalysisMan();

//********** AnalysisFlag **********//
 private:
  std::string ParamFileName;
  int DumpLevel;
  int CDSAnaMode;
  bool D5AnaFlag;
  bool FCAnaFlag;

 public:
  void SetParamFileName(const std::string &name){ ParamFileName = name; };
  void SetDumpLevel(const int &i){ DumpLevel=i; };
  void SetCDSAnaMode(const int &i){ CDSAnaMode=i; };
  void OffD5Ana(){ D5AnaFlag=false; };
  void OffFCAna(){ FCAnaFlag=false; };

//********** Classes for Analysis **********//
 private:
  //*** These classes manage out of AnalysisMan ***//
  ConfMan *confMan;
  ScalerMan *scaMan;
  EventHeader *header;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  CDSTrackingMan *cdstrackMan;
  BeamLineTrackMan *bltrackMan;
  //*** These classes manage in AnalysisMan ***//
  BeamSpectrometer *beamSpec;
  UshiwakaFieldMapMan *fieldMan;
  UshiwakaTableMan *tableMan;

 public:
  void SetConfMan(ConfMan *conf){ confMan = conf; };
  void SetScaMan(ScalerMan *sca){ scaMan = sca; };
  void SetHeader(EventHeader *head){ header = head; };
  void SetCDSHitMan(CDSHitMan *cds){ cdsMan = cds; };
  void SetBeamLineHitMan(BeamLineHitMan *blman){ blMan = blman; };
  void SetCDSTrackingMan(CDSTrackingMan *cdstrackman){ cdstrackMan = cdstrackman; };
  void SetBeamLineTrackMan(BeamLineTrackMan *bltrackman){ bltrackMan = bltrackman; };

  ConfMan *GetConfMan(){ return confMan; };
  ScalerMan *GetScaMan(){ return scaMan; };
  EventHeader *GetHeader(){ return header; };
  CDSHitMan *GetCDSHitMan(){ return cdsMan; };
  BeamLineHitMan *GetBeamLineHitMan(){ return blMan; };
  CDSTrackingMan *GetCDSTrackingMan(){ return cdstrackMan; };
  BeamLineTrackMan *GetBeamLineTrackingMan(){ return bltrackMan; };

  BeamSpectrometer *GetBeamSpectrometer(){ return beamSpec; };
  UshiwakaFieldMapMan *GetFieldMapMan(){ return fieldMan; };
  UshiwakaTableMan *GetTableMan(){ return tableMan; };

//********** Analysis Parameter **********//
 private:
  double BHDT0_Param[3][2];
  double Vtx_ZX_Param[2];
  double Vtx_ZY_Param[2];
  double Vtx_R2_Param;
  double Vtx_Z_Param[2];
  double CDS_PID_Param[4][2];
  double FC_PID_Param[3][2];
  double FC_Cut_Param[2];
  double FC_ADC_Param[4];
  double CDS_Lambda_Param[2];
  double CDS_K0_Param[2];
  double FC_Lambda_Param[2];
  double FC_Lambda1520_Param[2];
  double PC_ZX_Param[2];
  double USWK_Z;
  double FieldInZ;

  //***** for Event Flag *****//
 private:
  bool goodNCflag;
  bool Neutronflag;
  bool Gammaflag;
  bool Chargedflag;
  bool goodPCflag;
  bool goodCVCflag;

  bool VtxCutflag;
  bool FCADCCutflag;

  bool Npipflag;
  bool Npimflag;
  bool Nkflag;

  bool Ppimflag;
  bool Ppipflag;
  bool Pkflag;

  bool Dkflag;
  bool Dpimflag;
  bool Dpipflag;

 public:
  bool GoodNC() const { return goodNCflag; };
  bool Neutron() const { return Neutronflag; };
  bool Gamma() const { return Gammaflag; };
  bool ForwardCharged() const { return Chargedflag; };
  bool GoodPC() const { return goodPCflag; };
  bool GoodCVC() const { return goodCVCflag; };

  bool VtxCut() const { return VtxCutflag; };
  bool FCADCCut() const { return FCADCCutflag; };

  bool Npip() const { return Npipflag; };
  bool Npim() const { return Npimflag; };
  bool Nk() const { return Nkflag; };

  bool Ppip() const { return Ppipflag; };
  bool Ppim() const { return Ppimflag; };
  bool Pk() const { return Pkflag; };

  bool Dpip() const { return Dpipflag; };
  bool Dpim() const { return Dpimflag; };
  bool Dk() const { return Dkflag; };

  //***** Analysis Data *****//
 private:
  int beamPID;
  TLorentzVector beam_lmom;

  int blc1ID;
  int blc2ID;
  int bpcID;

  HodoscopeLikeHit *BHDHit;
  HodoscopeLikeHit *T0Hit;
  HodoscopeLikeHit *NC1st;
  HodoscopeLikeHit *ChargedHit;

  double D5Mom;
  double T0Vtx_tof;

 public:
  bool D5() const { return (D5Mom>0) ? true : false; };
  bool BHDT0() const { return (BHDHit!=0 && T0Hit!=0) ? true : false; };
  bool T0Vtx() const { return (T0Vtx_tof>-100) ? true : false; };

  int BeamPID() const { return beamPID; };
  TLorentzVector BeamLmom() const { return beam_lmom; };
  double D5mom() const { return D5Mom; };

  int BLC1_ID() const { return blc1ID; };
  int BLC2_ID() const { return blc2ID; };
  int BPC_ID() const { return bpcID; };

  double BHDT0_TOF() const { return (BHDHit!=0 && T0Hit!=0) ? T0Hit->ctmean()-BHDHit->ctmean() : -999.; }; 
  double T0Vtx_TOF() const { return T0Vtx_tof; };

  HodoscopeLikeHit* GetBHD(){ return BHDHit; };
  HodoscopeLikeHit* GetT0(){ return T0Hit; };
  HodoscopeLikeHit* GetNC1st(){ return NC1st; };
  HodoscopeLikeHit* GetFC(){ return ChargedHit; };

  //***** for Hit Position *****//
 private:
  TVector3 T0HitPosBLC;
  TVector3 T0HitPosBPC;

  TVector3 BPDHitPosBLC;
  TVector3 BPDHitPosBPC;

  TVector3 BPCHitPos;
  TVector3 VertexPos;
  TVector3 FDC1HitPos;
  TVector3 USWKHitPos;

  TVector3 NC1stHitPos;
  std::vector<TVector3> NCHitPos;
  TVector3 FCHitPos;
  TVector3 CalcFCHitPos;

 public:
  bool IsT0HitPosBLC() const { return (T0HitPosBLC.Mag()<1000) ? true : false; };
  bool IsT0HitPosBPC() const { return (T0HitPosBPC.Mag()<1000) ? true : false; };
  bool IsBPDHitPosBLC() const { return (T0HitPosBLC.Mag()<1000) ? true : false; };
  bool IsBPDHitPosBPC() const { return (T0HitPosBPC.Mag()<1000) ? true : false; };
  bool IsVertexPos() const { return (VertexPos.Mag()<1000) ? true : false; };
  bool IsFDC1HitPos() const { return (FDC1HitPos.Mag()<1000) ? true : false; };
  bool IsUSWKHitPos() const { return (USWKHitPos.Mag()<1000) ? true : false; };

  TVector3 GetT0HitPosByBLC() const { return T0HitPosBLC; };
  TVector3 GetT0HitPosByBPC() const { return T0HitPosBPC; };

  TVector3 GetBPDHitPosByBLC() const { return BPDHitPosBLC; };
  TVector3 GetBPDHitPosByBPC() const { return BPDHitPosBPC; };

  TVector3 GetBPCHitPos() const { return BPCHitPos; };
  TVector3 GetVertexPos() const { return VertexPos; };
  TVector3 GetFDC1HitPos() const { return FDC1HitPos; };
  TVector3 GetUSWKHitPos() const { return USWKHitPos; };

  TVector3 GetNC1stHitPos() const { return NC1stHitPos; };
  TVector3 GetNCHitPos(const int &i){ return NCHitPos[i]; };
  TVector3 GetFCHitPos() const { return FCHitPos; };

  //***** for Forward Counter ******//
 private:
  double T0PC_tof;
  double T0CVC_tof;
  std::vector<double> T0NC_tof;

  double VtxNC1st_tof;
  double NC1st_beta;
  double Neutron_mom;
  std::vector<double> VtxNC_tof;
  std::vector<double> NC_beta;
  std::vector<double> NC_offset;

  int FC_pid;
  double FCInParam[5];
  double FC_fl;
  double FC_mom;
  double VtxFC_tof;
  double FC_beta;
  double FC_mass2;
  double FC_mom2; // momentum by TOF
  double FC_offset;

  MissingMass BeamN_mm;
  MissingMass BeamP_mm;
  MissingMass BeamD_mm;

  InvariantMass Npip_im;
  InvariantMass Npim_im;
  InvariantMass Nk_im;

  InvariantMass Ppip_im;
  InvariantMass Ppim_im;
  InvariantMass Pk_im;

  InvariantMass Dpip_im;
  InvariantMass Dpim_im;
  InvariantMass Dk_im;

 public:
  bool T0PC() const { return (T0PC_tof>-100 ) ? true : false; };
  bool T0CVC() const { return (T0CVC_tof>-100 ) ? true : false; };
  double T0PC_TOF() const { return T0PC_tof; };
  double T0CVC_TOF() const { return T0CVC_tof; };

  double T0NC_TOF(const int i){ return T0NC_tof[i]; };
  double VtxNC1st_TOF() const { return VtxNC1st_tof; };
  double NC1st_Beta() const { return NC1st_beta; };
  double NeutronMom() const { return Neutron_mom; };
  double VtxNC_TOF(const int &i) const { return VtxNC_tof[i]; };
  double NC_Beta(const int &i) const { return NC_beta[i]; };
  double NC_Offset(const int &i) const{ return NC_offset[i]; };

  int FC_PID() const { return FC_pid; };
  bool GetFCInParam(double param[5]){ param[0] = FCInParam[0], param[1] = FCInParam[1], param[2] = FCInParam[2], param[3] = FCInParam[3], param[4] = FCInParam[4]; };
  double FCFL() const { return FC_fl; };
  double FCMom() const { return FC_mom; };
  double VtxFC_TOF() const { return VtxFC_tof; };
  double FC_Beta() const { return FC_beta; };
  double FC_Mass2() const { return FC_mass2; };
  double FCMom2() const { return FC_mom2; };
  double FC_Offset() const { return FC_offset; };

  MissingMass* BeamN_MM(){ return &BeamN_mm; };
  MissingMass* BeamP_MM(){ return &BeamP_mm; };
  MissingMass* BeamD_MM(){ return &BeamD_mm; };

  InvariantMass* NPip_IM(){ return &Npip_im; };
  InvariantMass* NPim_IM(){ return &Npim_im; };
  InvariantMass* NK_IM(){ return &Nk_im; };

  InvariantMass* PPip_IM(){ return &Ppip_im; };
  InvariantMass* PPim_IM(){ return &Ppim_im; };
  InvariantMass* PK_IM(){ return &Pk_im; };

  InvariantMass* DPip_IM(){ return &Dpip_im; };
  InvariantMass* DPim_IM(){ return &Dpim_im; };
  InvariantMass* DK_IM(){ return &Dk_im; };

  //***** for CDS Data *****//
 private:
  TrackVertex vertexBeam;

  std::vector<int> GoodTrack_ID;
  std::vector<double> CDS_beta;
  std::vector<double> CDS_mass2;
  std::vector<double> CDH_offset;

  std::vector<int> cdsPim_ID;
  std::vector<int> cdsPip_ID;
  std::vector<int> cdsKp_ID;
  std::vector<int> cdsKm_ID;
  std::vector<int> cdsP_ID;
  std::vector<int> cdsD_ID;
  std::vector<int> cdsLow_ID;
  std::vector<int> cdsHigh_ID;

  std::vector<InvariantMass> cdsPiPiIM;
  std::vector<InvariantMass> cdsPiPIM;
  std::vector<InvariantMass> cdsKPIM;

 public:
  int nGoodTrack()  const { return GoodTrack_ID.size(); };
  int nCDSpim()     const { return cdsPim_ID.size(); };
  int nCDSpip()     const { return cdsPip_ID.size(); };
  int nCDSkm()      const { return cdsKm_ID.size(); };
  int nCDSkp()      const { return cdsKp_ID.size(); };
  int nCDSp()       const { return cdsP_ID.size(); };
  int nCDSd()       const { return cdsD_ID.size(); };
  int nCDSlow()     const { return cdsLow_ID.size(); };
  int nCDShigh()    const { return cdsHigh_ID.size(); };

  TrackVertex *VertexBeam(){ return &vertexBeam; };

  CDSTrack* CDSGoodTrack(const int &i){ return cdstrackMan-> Track(GoodTrack_ID[i]); };
  HodoscopeLikeHit* CDHHit(const int &i){ return cdstrackMan-> Track(GoodTrack_ID[i])-> CDHHit(cdsMan); };
  double CDS_Beta(const int &i) const { return CDS_beta[i]; };
  double CDS_Mass2(const int &i) const { return CDS_mass2[i]; };
  double CDS_Mass(const int &i) const { return sqrt(CDS_mass2[i]); };
  double CDH_Offset(const int &i) const{ return CDH_offset[i]; };

  CDSTrack* CDSpim(const int &i){ return cdstrackMan-> Track(cdsPim_ID[i]); };
  CDSTrack* CDSpip(const int &i){ return cdstrackMan-> Track(cdsPip_ID[i]); };
  CDSTrack* CDSkm(const int &i){ return cdstrackMan-> Track(cdsKm_ID[i]); };
  CDSTrack* CDSkp(const int &i){ return cdstrackMan-> Track(cdsKp_ID[i]); };
  CDSTrack* CDSp(const int &i){ return cdstrackMan-> Track(cdsP_ID[i]); };
  CDSTrack* CDSd(const int &i){ return cdstrackMan-> Track(cdsD_ID[i]); };
  CDSTrack* CDSlow(const int &i){ return cdstrackMan-> Track(cdsLow_ID[i]); };
  CDSTrack* CDShigh(const int &i){ return cdstrackMan-> Track(cdsHigh_ID[i]); };

  int nCDSPiPi() const { return cdsPiPiIM.size(); };
  int nCDSPiP() const { return cdsPiPIM.size(); };
  int nCDSKP() const { return cdsKPIM.size(); };

  InvariantMass* CDSPiPiIM(const int &i){ return &cdsPiPiIM[i]; };
  InvariantMass* CDSPiPIM(const int &i){ return &cdsPiPIM[i]; };
  InvariantMass* CDSKPIM(const int &i){ return &cdsKPIM[i]; };

//********** Main Routine **********//
 public:
  bool Initialize();
  void Execute();
  void Clear();
  bool Finalize();

  void SetHodo();
  void AnaBeam();
  void AnaCDS();
  void CalcCDSIM();
  void AnaNC();
  void AnaFC();
  void CalcMM();
  void CalcForwardIM();
  bool CalcCDSIM(const int &trackID1, const int &trackID2, InvariantMass &im);

//********** for Hodoscopes **********//
 private:
  std::vector<int> BHD_ID;
  std::vector<int> BHDpost_ID;
  std::vector<int> T0_ID;
  std::vector<int> T0pre_ID;
  std::vector<int> T0post_ID;
  std::vector<int> E0_ID;
  std::vector<int> CVC_ID;
  std::vector<int> BPD_ID;
  std::vector<int> BVC_ID;
  std::vector<int> PC_ID;
  std::vector<int> BD_ID;
  std::vector<int> LB_ID;
  std::vector<int> WVC_ID;
  std::vector<int> HVC1_ID;
  std::vector<int> HVC2_ID;
  std::vector<int> Temp1_ID;
  std::vector<int> Temp2_ID;

  std::vector<int> NC_ID[NumOfNCLayer];

  std::vector<int> CDH_ID;
  std::vector<int> IH_ID;

 public:
  int nBHD()     const { return BHD_ID.size(); };
  int nBHDpost() const { return BHDpost_ID.size(); };
  int nT0()      const { return T0_ID.size(); };
  int nT0pre()   const { return T0pre_ID.size(); };
  int nT0post()  const { return T0post_ID.size(); };
  int nE0()      const { return E0_ID.size(); };
  int nDEF()     const { return E0_ID.size(); };
  int nCVC()     const { return CVC_ID.size(); };
  int nBPD()     const { return BPD_ID.size(); };
  int nBVC()     const { return BVC_ID.size(); };
  int nPC()      const { return PC_ID.size(); };
  int nBD()      const { return BD_ID.size(); };
  int nLB()      const { return LB_ID.size(); };
  int nWVC()     const { return WVC_ID.size(); };
  int nHVC1()    const { return HVC1_ID.size(); };
  int nHVC2()    const { return HVC2_ID.size(); };
  int nTemp1()   const { return Temp1_ID.size(); };
  int nTemp2()   const { return Temp2_ID.size(); };

  int nNC() const;
  int nNC(const int &layer) const { return NC_ID[layer-1].size(); };

  int nCDH() const { return CDH_ID.size(); };
  int nIH()  const { return IH_ID.size(); };

  HodoscopeLikeHit* BHD(const int &i)     { return blMan->BHD(BHD_ID[i]); };
  HodoscopeLikeHit* BHDpost(const int &i) { return blMan->BHDpost(BHD_ID[i]); };
  HodoscopeLikeHit* T0(const int &i)      { return blMan->T0(T0_ID[i]); };
  HodoscopeLikeHit* T0pre(const int &i)   { return blMan->T0pre(T0pre_ID[i]); };
  HodoscopeLikeHit* T0post(const int &i)  { return blMan->T0post(T0post_ID[i]); };
  HodoscopeLikeHit* E0(const int &i)      { return blMan->E0(E0_ID[i]); };
  HodoscopeLikeHit* DEF(const int &i)     { return blMan->E0(E0_ID[i]); };
  HodoscopeLikeHit* CVC(const int &i)     { return blMan->CVC(CVC_ID[i]); };
  HodoscopeLikeHit* BPD(const int &i)     { return blMan->BPD(BPD_ID[i]); };
  HodoscopeLikeHit* BVC(const int &i)     { return blMan->BVC(BVC_ID[i]); };
  HodoscopeLikeHit* PC(const int &i)      { return blMan->PC(PC_ID[i]); };
  HodoscopeLikeHit* BD(const int &i)      { return blMan->BD(BD_ID[i]); };
  HodoscopeLikeHit* LB(const int &i)      { return blMan->LB(LB_ID[i]); };
  HodoscopeLikeHit* WVC(const int &i)     { return blMan->WVC(WVC_ID[i]); };
  HodoscopeLikeHit* HVC1(const int &i)    { return blMan->HVC1(HVC1_ID[i]); };
  HodoscopeLikeHit* HVC2(const int &i)    { return blMan->HVC2(HVC2_ID[i]); };
  HodoscopeLikeHit* Temp1(const int &i)   { return blMan->Temp1(Temp1_ID[i]); };
  HodoscopeLikeHit* Temp2(const int &i)   { return blMan->Temp2(Temp2_ID[i]); };

  HodoscopeLikeHit* NC(const int &layer, const int &i) { return blMan->NC(NC_ID[layer-1][i]); };
  HodoscopeLikeHit* NC(const int &i);

  HodoscopeLikeHit *CDH(const int &i) { return cdsMan->CDH(CDH_ID[i]); };
  HodoscopeLikeHit *IH(const int &i)  { return cdsMan->IH(IH_ID[i]); };

//********** for Dumpping **********//
  void DumpParam();
  void DumpStatus();
  void DumpHit();
  void DumpBLTrack();
  void DumpCDSTrack();
  void DumpHitPos();
  void DumpForward();

  void Wait();
};
#endif
