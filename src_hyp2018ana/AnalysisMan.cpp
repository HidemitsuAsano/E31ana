#include "AnalysisMan.h"

static const double PI = 6*asin(0.5);

static const std::string DefaultParamFileName = "param/Analysis.param";
static const int MAXCHAR = 144;

AnalysisMan::AnalysisMan()
  : DumpLevel(0), CDSAnaMode(1), D5AnaFlag(true), FCAnaFlag(true), ParamFileName(DefaultParamFileName),
    confMan(0), scaMan(0), header(0), cdsMan(0), blMan(0), cdstrackMan(0), bltrackMan(0), beamSpec(0), fieldMan(0), tableMan(0)
{
  Clear();
}

AnalysisMan::~AnalysisMan()
{
  if( beamSpec!=0 ) delete beamSpec;
  if( fieldMan!=0 ) delete fieldMan;
  if( tableMan!=0 ) delete tableMan;
}

void AnalysisMan::Execute()
{
  SetHodo();
  bltrackMan-> DoTracking(blMan, confMan, true, true);
  AnaBeam();
  if( CDSAnaMode==1 ) AnaCDS();
  else if( CDSAnaMode==2 ){
    if( nNC()>0 || nCVC()>0 || nPC()>0 ) AnaCDS();
  }
  else if( CDSAnaMode==3 ){
    if( nNC()>0 ) AnaCDS();
  }
  else if( CDSAnaMode==4 ){
    if( nCVC()>0 || nPC()>0 ) AnaCDS();
  }
  else if( CDSAnaMode>4 ){
    if( nT0()==1 && bpcID>=0 ){
      if( CDSAnaMode==5 ) AnaCDS();
      else if( CDSAnaMode==6 ){
	if( nNC()>0 || nCVC()>0 || nPC()>0 ) AnaCDS();
      }
      else if( CDSAnaMode==7 ){
	if( nNC()>0 ) AnaCDS();
      }
      else if( CDSAnaMode==8 ){
	if( nCVC()>0 || nPC()>0 ) AnaCDS();
      }
      else if( CDSAnaMode==9 ){
	if( nBVC()==0 && nCVC()==0 && nPC()>0 ) AnaCDS();
      }
      else if( CDSAnaMode==10 ){
	if( nBVC()==1 && bltrackMan->ntrackFDC1()==1 ){
	  if( nCVC()>0 || nPC()>0 ) AnaCDS();
	}
      }
    }
  }

  AnaNC();
  if( FCAnaFlag ){
    AnaFC();
  }
  CalcMM();
  CalcForwardIM();
}

void AnalysisMan::AnaBeam()
{
  if( nBHD()==1 ) BHDHit = BHD(0);
  if( nT0()==1  ) T0Hit  = T0(0);
  if( BHDT0() ){
    double tof = BHDT0_TOF();
    if( BHDT0_Param[0][0]<tof && tof<BHDT0_Param[0][1] ) beamPID = Beam_Pion;
    else if( BHDT0_Param[1][0]<tof && tof<BHDT0_Param[1][1] ) beamPID = Beam_Kaon;
    else if( BHDT0_Param[2][0]<tof && tof<BHDT0_Param[2][1] ) beamPID = Beam_Proton;
    else beamPID = Beam_Other;
  }

  if( bltrackMan->ntrackBLC1()==1 ) blc1ID=0;
  if( bltrackMan->ntrackBLC2()==1 ) blc2ID=0;

  if( D5AnaFlag ){
    if( blc1ID>=0 && blc2ID>=0 ){
      beamSpec-> TMinuitFit( bltrackMan->trackBLC1(blc1ID), bltrackMan->trackBLC2(blc2ID), confMan );
      D5Mom = beamSpec-> mom();
    }
  }

  if( bltrackMan->ntrackBLC2()==1 ){
    TVector3 pos, dir, gpos, gdir;
    double wid, len, th, lv;

    LocalTrack *trackBLC2 = bltrackMan->trackBLC2(0);

    confMan-> GetGeomMapManager()-> GetGParam(CID_T0, gpos, gdir);
    confMan-> GetGeomMapManager()-> GetParam(CID_T0, 1, pos, dir, wid, len, th, lv);
    double T0z = gpos.z()+pos.z()-0.5*th;
    double T0x, T0y;
    trackBLC2-> XYPosatZ(T0z, T0x, T0y);
    T0HitPosBLC.SetXYZ(T0x, T0y, T0z);

    confMan-> GetGeomMapManager()-> GetGParam(CID_BPD, gpos, gdir);
    confMan-> GetGeomMapManager()-> GetParam(CID_BPD, 1, pos, dir, wid, len, th, lv);
    double BPDz = gpos.z()-0.5*th;
    double BPDx, BPDy;
    trackBLC2-> XYPosatZ(BPDz, BPDx, BPDy);
    BPDHitPosBLC.SetXYZ(BPDx, BPDy, BPDz);

    if( bltrackMan->ntrackBPC()>0 ){
      double BLC2_dx = trackBLC2-> gdx();
      double BLC2_dy = trackBLC2-> gdy();

      double tmp_dxdz_diff=10e10;
      for( int bpc_id=0; bpc_id<bltrackMan->ntrackBPC(); bpc_id++ ){
	LocalTrack *trackBPC = bltrackMan->trackBPC(bpc_id);

	double BPC_dx = trackBPC-> gdx();
	double BPC_dy = trackBPC-> gdy();

	double dxdz_diff = (BLC2_dx-BPC_dx)*(BLC2_dx-BPC_dx)+(BLC2_dy-BPC_dy)*(BLC2_dy-BPC_dy);
	if( dxdz_diff<tmp_dxdz_diff ){
	  tmp_dxdz_diff = dxdz_diff;
	  bpcID = bpc_id;
	}
      }

      LocalTrack *trackBPC = bltrackMan->trackBPC(bpcID);
      trackBPC-> XYPosatZ(T0z, T0x, T0y);
      T0HitPosBPC.SetXYZ(T0x, T0y, T0z);

      trackBPC-> XYPosatZ(BPDz, BPDx, BPDy);
      BPDHitPosBPC.SetXYZ(BPDx, BPDy, BPDz);
      double BPCz = trackBPC-> hit(0,0)-> gz();
      double BPCx, BPCy;
      trackBPC-> XYPosatZ(BPCz, BPCx, BPCy);
      BPCHitPos.SetXYZ(BPCx, BPCy, BPCz);
    }
  }

  if( bpcID>=0 ){
    LocalTrack *trackBPC = bltrackMan->trackBPC(bpcID);
    double dx = trackBPC->gdx();
    double dy = trackBPC->gdy();
    TVector3 mom(dx, dy, 1.0);
    if( D5() )  mom.SetMag(D5mom());
    else        mom.SetMag(1.0);

    if( beamPID==Beam_Pion )        beam_lmom.SetVectM(mom, piMass);
    else if( beamPID==Beam_Proton ) beam_lmom.SetVectM(mom, pMass);
    else if( beamPID==Beam_Kaon   ) beam_lmom.SetVectM(mom, kpMass);
    else                            beam_lmom.SetVectM(mom, kpMass);
  }

  if( bltrackMan-> ntrackFDC1()==1 ){
    LocalTrack *trackFDC1 = bltrackMan-> trackFDC1(0);
    double FDC1z = trackFDC1-> hit(0,0)-> gz();
    double FDC1x, FDC1y;
    trackFDC1-> XYPosatZ(FDC1z, FDC1x, FDC1y);
    FDC1HitPos.SetXYZ(FDC1x, FDC1y, FDC1z);
  }

  if( nT0()==1 && nPC()==1 )   T0PC_tof = PC(0)->ctmean() - T0(0)->ctmean();
  if( nT0()==1 && nCVC()==1 )  T0CVC_tof = CVC(0)->ctmean() - T0(0)->ctmean();
}

void AnalysisMan::AnaCDS()
{
  cdstrackMan-> Execute(cdsMan, confMan);
  cdstrackMan-> Calc(cdsMan, confMan, true);

  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    int goodtrack_id = cdstrackMan->GoodTrackID(i);
    cdstrackMan-> CalcVertex_beam(goodtrack_id, bltrackMan, confMan);
  }

  // Calc T0-Vertex TOF etc.
  if( bpcID>=0 && nT0()==1 ){
    double tmp_dis = 10e10;
    for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
      int good_track_id = cdstrackMan-> GoodTrackID(id);
      CDSTrack *cdsTrack = cdstrackMan-> Track(good_track_id);
      double param[5];
      cdsTrack-> GetParameters(param);

      if( cdsTrack-> CDHFlag() ){
	TrackVertex vertex_beam;
	if( cdstrackMan-> GetVertex_beam(good_track_id, bpcID, vertex_beam) ){
	  TVector3 vtxPos1 = vertex_beam.GetVertexPos1();
	  TVector3 vtxPos2 = vertex_beam.GetVertexPos2();
	  TVector3 vtxPosMean = vertex_beam.GetVertexPos_mean();

	  LocalTrack *trackBPC = bltrackMan->trackBPC(bpcID);
	  double vtx_x, vtx_y;
	  trackBPC->XYPosatZ(vtxPosMean.z(), vtx_x, vtx_y);
	  TVector3 vtxPosBeam(vtx_x, vtx_y, vtxPosMean.z());

	  double T0CDH_tof = cdsTrack->CDHHit(cdsMan)->ctmean() - T0(0)->ctmean();
	  double T0Vtx_dis = (T0HitPosBPC-vtxPosBeam).Mag();

	  double beam_v, beam_mom;
	  double beam_dis = (VertexPos-T0HitPosBPC).Mag();

	  if( D5() ) beam_mom=D5Mom;
	  else beam_mom=1.0;

	  if( beamPID==Beam_Proton )    beam_v = 100.*Const*beam_mom/sqrt(beam_mom*beam_mom+pMass*pMass);
	  else if( beamPID==Beam_Pion ) beam_v = 100.*Const*beam_mom/sqrt(beam_mom*beam_mom+piMass*piMass);
	  else if( beamPID==Beam_Kaon ) beam_v = 100.*Const*beam_mom/sqrt(beam_mom*beam_mom+kpMass*kpMass);
	  else                          beam_v = 100.*Const*beam_mom/sqrt(beam_mom*beam_mom+kpMass*kpMass);
	  double vtxCDH_tof = T0CDH_tof - T0Vtx_dis/beam_v;

	  TVector3 CDHvertex = cdsTrack-> CDHVertex();
	  double phi1 = MathTools::CalcHelixPhi(vtx_x, vtx_y, param);
	  double phi2 = MathTools::CalcHelixPhi(CDHvertex.x(), CDHvertex.y(), param);

	  double rho  = cdsTrack->Rho();
	  double tlam = param[4];

	  double cds_dis   = rho*fabs(phi1-phi2)*sqrt(1.0+tlam*tlam);
	  double cds_beta  = cds_dis/(100.*Const*vtxCDH_tof);
	  double cds_mom   = cdsTrack-> Momentum();
	  double cds_mass2 = cds_mom*cds_mom*(1./(cds_beta*cds_beta)-1.);
	  double cds_mass  = sqrt(cds_mass2);

          double calc_cds_pi_beta = cds_mom/sqrt(cds_mom*cds_mom+piMass*piMass);
          double calc_cds_tof = cds_dis/(100.*Const*calc_cds_pi_beta);
          double offset = calc_cds_tof-vtxCDH_tof;
          CDH_offset.push_back(offset);

	  GoodTrack_ID.push_back(good_track_id);
	  CDS_beta.push_back(cds_beta);
	  CDS_mass2.push_back(cds_mass2);

	  // CDS PID
          if( cds_mass<CDS_PID_Param[0][0] ){
            cdsTrack-> SetPID(CDS_Other);
            cdsLow_ID.push_back(good_track_id);
          }
          else if( CDS_PID_Param[0][0]<cds_mass && cds_mass<CDS_PID_Param[0][1] ){
            if( cds_mom<0 ){
              cdsTrack-> SetPID(CDS_PiMinus);
              cdsPim_ID.push_back(good_track_id);
            }
	    else{
              cdsTrack-> SetPID(CDS_PiPlus);
              cdsPip_ID.push_back(good_track_id);
            }
          }
          else if( CDS_PID_Param[1][0]<cds_mass && cds_mass<CDS_PID_Param[1][1] ){
            if( cds_mom<0 ){
              cdsTrack-> SetPID(CDS_Kaon);
              cdsKm_ID.push_back(good_track_id);
            }
	    else{
              cdsTrack-> SetPID(CDS_Kaon);
              cdsKp_ID.push_back(good_track_id);
            }
          }
          else if( CDS_PID_Param[2][0]<cds_mass && cds_mass<CDS_PID_Param[2][1] ){
            if( cds_mom<0 ){
              cdsTrack-> SetPID(CDS_Other);
              cdsHigh_ID.push_back(good_track_id);
            }
	    else{
              cdsTrack-> SetPID(CDS_Proton);
              cdsP_ID.push_back(good_track_id);
            }
          }
          else if( CDS_PID_Param[3][0]<cds_mass && cds_mass<CDS_PID_Param[3][1] ){
            if( cds_mom<0 ){
              cdsTrack-> SetPID(CDS_Other);
              cdsHigh_ID.push_back(good_track_id);
            }
	    else{
              cdsTrack-> SetPID(CDS_Deuteron);
              cdsD_ID.push_back(good_track_id);
            }
          }
          else{
            cdsTrack-> SetPID(CDS_Other);
            cdsHigh_ID.push_back(good_track_id);
          }

	  double dis = (vtxPos1-vtxPos2).Mag();
	  if( dis<tmp_dis ){
	    tmp_dis = dis;
	    VertexPos = vtxPos2;
	    vertexBeam = vertex_beam;
	  }
	}
	else std::cout<<" !!! fault CDSTrackMan-> GetVertex_beam("<<good_track_id<<","<<bpcID<<")"<<std::endl;
      }
    }
    if( IsVertexPos() && IsT0HitPosBPC() ){
      double beam_v, beam_mom;
      double beam_dis = (VertexPos-T0HitPosBPC).Mag();

      if( D5() ) beam_mom=D5Mom;
      else beam_mom=1.0;

      if( beamPID==Beam_Proton )    beam_v = 100.*Const*beam_mom/sqrt(beam_mom*beam_mom+pMass*pMass);
      else if( beamPID==Beam_Pion ) beam_v = 100.*Const*beam_mom/sqrt(beam_mom*beam_mom+piMass*piMass);
      else if( beamPID==Beam_Kaon ) beam_v = 100.*Const*beam_mom/sqrt(beam_mom*beam_mom+kpMass*kpMass);
      else                          beam_v = 100.*Const*beam_mom/sqrt(beam_mom*beam_mom+kpMass*kpMass);

      T0Vtx_tof = beam_dis/beam_v;

      if( VertexPos.z()<Vtx_Z_Param[0] || Vtx_Z_Param[1]<VertexPos.z() ) VtxCutflag = true;
      else{
	double center_x = Vtx_ZX_Param[1]*VertexPos.z()+Vtx_ZX_Param[0];
	double center_y = Vtx_ZY_Param[1]*VertexPos.z()+Vtx_ZY_Param[0];
	double vtx_r2 = (VertexPos.x()-center_x)*(VertexPos.x()-center_x) + (VertexPos.y()-center_y)*(VertexPos.y()-center_y);
	if( vtx_r2>Vtx_R2_Param ) VtxCutflag = true;
      }
    }
  }

  CalcCDSIM();
}

void AnalysisMan::CalcCDSIM()
{
  for( int i=0; i<nCDSpim(); i++ ){
    for( int j=0; j<nCDSpip(); j++ ){
      InvariantMass im;
      if( CalcCDSIM( cdsPim_ID[i], cdsPip_ID[j], im) ){
	cdsPiPiIM.push_back(im);
      }
    }
  }

  for( int i=0; i<nCDSpim(); i++ ){
    for( int j=0; j<nCDSp(); j++ ){
      InvariantMass im;
      if( CalcCDSIM( cdsPim_ID[i], cdsP_ID[j], im) ){
	cdsPiPIM.push_back(im);
      }
    }
  }

  for( int i=0; i<nCDSkm(); i++ ){
    for( int j=0; j<nCDSp(); j++ ){
      InvariantMass im;
      if( CalcCDSIM( cdsKm_ID[i], cdsP_ID[j], im) ){
	cdsKPIM.push_back(im);
      }
    }
  }
}

bool AnalysisMan::CalcCDSIM(const int &trackID1, const int &trackID2, InvariantMass &im)
{
  CDSTrack *track1 = cdstrackMan-> Track(trackID1);
  CDSTrack *track2 = cdstrackMan-> Track(trackID2);

  int pid1 = track1-> PID();
  int pid2 = track2-> PID();

  double mass1, mass2;
  if( pid1==CDS_PiPlus || pid1==CDS_PiMinus ) mass1=piMass;
  else if( pid1==CDS_Kaon )     mass1=kpMass;
  else if( pid1==CDS_Proton )   mass1=pMass;
  else if( pid1==CDS_Deuteron ) mass1=dMass;
  else{
    std::cout<<" !!! AnalysisMan::CalcCDSIM() track1 is unknown track return false"<<std::endl;
    return false;
  }
  if( pid2==CDS_PiPlus || pid2==CDS_PiMinus ) mass2=piMass;
  else if( pid2==CDS_Kaon )     mass2=kpMass;
  else if( pid2==CDS_Proton )   mass2=pMass;
  else if( pid2==CDS_Deuteron ) mass2=dMass;
  else{
    std::cout<<" !!! AnalysisMan::CalcCDSIM() track2 is unknown track return false"<<std::endl;
    return false;
  }

  TrackVertex vertex;
  TVector3 vtx_pos0, vtx_pos1, vtx_pos2;
  if( cdstrackMan->GetVertex(trackID1, trackID2, vertex) ){
    vtx_pos0 = vertex.GetVertexPos_mean();
    if( trackID1<trackID2 ){
      vtx_pos1 = vertex.GetVertexPos1();
      vtx_pos2 = vertex.GetVertexPos2();
    }
    else{
      vtx_pos1 = vertex.GetVertexPos2();
      vtx_pos2 = vertex.GetVertexPos1();
    }
  }
  else{
    std::cout<<" !!! AnalysisMan::CalcCDSIM() GetVertex fault TrackID1:"<<trackID1<<" TrackID2:"<<trackID2<<" return false"<<std::endl;
    return false;
  }

  TVector3 mom1, mom2;
  if( !track1-> GetMomentum(vtx_pos1, mom1) ){
    std::cout<<" !!! AnalysisMan::CalcCDSIM() GetMomentum fault return false"<<std::endl;
    return false;
  }
  if( !track2-> GetMomentum(vtx_pos2, mom2) ){
    std::cout<<" !!! AnalysisMan::CalcCDSIM() GetMomentum fault return false"<<std::endl;
    return false;
  }

  TLorentzVector lmom1;
  lmom1.SetVectM(mom1, mass1);
  TLorentzVector lmom2;
  lmom2.SetVectM(mom2, mass2);

  im.SetParentLmom(lmom1+lmom2);
  im.SetSonPID(pid1);
  im.SetSonLmom(lmom1);
  im.SetSonPID(pid2);
  im.SetSonLmom(lmom2);

  return true;
}

void AnalysisMan::AnaNC()
{
  if( nT0()==1 &&  nCVC()==0 && nNC()>0 ){
    for( int layer=1; layer<=NumOfNCLayer; layer++ ){
      if( nNC(layer)>0 ){
	double tmp_dE=-999;
	for( int i=0; i<nNC(layer); i++ ){
	  if( NC(layer,i)->emean()>tmp_dE ){
	    tmp_dE = NC(layer,i)->emean();
	    NC1st = NC(layer,i);
	  }
	}
	if( tmp_dE>0 ) break;
      }
    }

    for( int i=0; i<nNC(); i++ ){
      double tof = NC(i)->ctmean() - T0(0)->ctmean();
      T0NC_tof.push_back(tof);
    }

    if( IsVertexPos() && NC1st!=0 ){
      goodNCflag = true;

      TVector3 pos, dir, gpos, gdir;
      double wid, len, th, lv;
      int seg = NC1st-> seg();
      confMan-> GetGeomMapManager()-> GetGParam(CID_NC, gpos, gdir);
      confMan-> GetGeomMapManager()-> GetParam(CID_NC, seg, pos, dir, wid, len, th, lv);
      NC1stHitPos.SetXYZ(gpos.x()+pos.x(), gpos.y()+pos.y(), gpos.z()+pos.z()-0.5*th);
      double vtxNC1st_fl = (NC1stHitPos-VertexPos).Mag();
      VtxNC1st_tof = NC1st->ctmean() - T0(0)-> ctmean() - T0Vtx_tof;
      NC1st_beta = vtxNC1st_fl/(100.*Const*VtxNC1st_tof);

      if( NC1st_beta>0.95 ) Gammaflag = true;
      else{
	Neutronflag = true;
	Neutron_mom = nMass*NC1st_beta/sqrt(1.0-NC1st_beta*NC1st_beta);
      }
      for( int i=0; i<nNC(); i++ ){
	int NCseg = NC(i)-> seg();
	confMan-> GetGeomMapManager()-> GetGParam(CID_NC, gpos, gdir);
	confMan-> GetGeomMapManager()-> GetParam(CID_NC, NCseg, pos, dir, wid, len, th, lv);
	TVector3 NChitpos(gpos.x()+pos.x(), gpos.y()+pos.y(), gpos.z()+gpos.z());

	double vtxNC_fl = (NChitpos-VertexPos).Mag();
	double vtxNC_tof = NC(i)->ctmean() - T0(0)-> ctmean() - T0Vtx_tof;
	double beta = vtxNC_fl/(100.*Const*vtxNC_tof);

	NCHitPos.push_back(NChitpos);
	VtxNC_tof.push_back(vtxNC_tof);
	NC_beta.push_back(beta);

        double Calc_VtxNC_tof = vtxNC_fl/(100.*Const);
        double offset = Calc_VtxNC_tof-vtxNC_tof;
	NC_offset.push_back(offset);
      }
    }
  }
}

void AnalysisMan::AnaFC()
{
  if( nT0()==1 && IsVertexPos() && IsFDC1HitPos() && nBVC()>0 ){
    TVector3 pos, dir, gpos, gdir;
    double len, wid, th, lv;

    double in_dx = (FDC1HitPos.x()-VertexPos.x())/(FDC1HitPos.z()-VertexPos.z());
    double in_dy = (FDC1HitPos.y()-VertexPos.y())/(FDC1HitPos.z()-VertexPos.z());
    double in_x0 = VertexPos.x() - in_dx*VertexPos.z();
    double in_y0 = VertexPos.y() - in_dy*VertexPos.z();

    double FieldInX = in_dx*FieldInZ + in_x0;
    double FieldInY = in_dy*FieldInZ + in_y0;

    double USWK_X = in_dx*USWK_Z + in_x0;
    double USWK_Y = in_dy*USWK_Z + in_y0;
    USWKHitPos.SetXYZ(USWK_X, USWK_Y, USWK_Z);

    if( nPC()==1 && nCVC()==0 ){
      Chargedflag = true, goodPCflag = true;
      ChargedHit = PC(0);

      int seg = PC(0)-> seg();
      confMan-> GetGeomMapManager()-> GetGParam( CID_PC, gpos, gdir );
      confMan-> GetGeomMapManager()-> GetParam( CID_PC, seg, pos, dir, len, wid, th, lv );
      double xx = gpos.x() + cos(PI*gdir.y()/180.)*pos.x() + sin(PI*gdir.y()/180.)*(pos.z()-0.5*th);
      double zz = gpos.z() - sin(PI*gdir.y()/180.)*pos.x() + cos(PI*gdir.y()/180.)*(pos.z()-0.5*th);
      double yy = gpos.y() + pos.y();
      FCHitPos.SetXYZ(xx,yy,zz);
    }

    if( nPC()==0 && nCVC()==1 ){
      Chargedflag = true, goodCVCflag = true;
      ChargedHit = CVC(0);

      int seg = CVC(0)-> seg();
      confMan-> GetGeomMapManager()-> GetGParam( CID_CVC, gpos, gdir );
      confMan-> GetGeomMapManager()-> GetParam( CID_CVC, seg, pos, dir, len, wid, th, lv );
      double xx = gpos.x() + pos.x();
      double yy = gpos.y() + pos.y();
      double zz = gpos.z() + pos.z()-0.5*th;
      FCHitPos.SetXYZ(xx,yy,zz);
    }

    if( Chargedflag ){
      double out_dx = (FCHitPos.x()-USWKHitPos.x())/(FCHitPos.z()-USWKHitPos.z());
      double out_dy = (FCHitPos.y()-USWKHitPos.y())/(FCHitPos.z()-USWKHitPos.z());
      double in_ang = atan(in_dx);
      double out_ang = atan(out_dx);
      FCInParam[0] = FieldInX, FCInParam[1] = FieldInY, FCInParam[2] = in_dx, FCInParam[3] = in_dy, FCInParam[4] = in_ang - out_ang;

      double uswk_fl;
      if( tableMan-> GetParam(FCInParam[0], FCInParam[1], FCInParam[2], FCInParam[3], FCInParam[4], FC_mom, uswk_fl) ){
	TVector3 fieldInPos(FieldInX, FieldInY, FieldInZ);
	TVector3 inMom(in_dx, in_dy, 1.0);
	inMom.SetMag(FC_mom);

	fieldMan-> RungeKutta( "proton", fieldInPos, inMom );
	USWKTrack *uswk_track = fieldMan-> GetUSWKTrack(0);
	double out_time = uswk_track-> outtime();
	TVector3 out_pos = uswk_track-> outposition();
	TVector3 out_mom = uswk_track-> outmomentum();

	double calc_out_dx = out_mom.x()/out_mom.z();
	double calc_out_dy = out_mom.y()/out_mom.z();
	double calc_out_x0 = out_pos.x() - calc_out_dx*out_pos.z();
	double calc_out_y0 = out_pos.y() - calc_out_dy*out_pos.z();

	if( goodPCflag ){
	  double calc_hitpos_z = ( calc_out_x0-PC_ZX_Param[0] )/( PC_ZX_Param[1]-calc_out_dx );
	  double calc_hitpos_x = calc_out_dx*calc_hitpos_z + calc_out_x0;
	  double calc_hitpos_y = calc_out_dy*calc_hitpos_z + calc_out_y0;
	  CalcFCHitPos.SetXYZ(calc_hitpos_x, calc_hitpos_y, calc_hitpos_z);
	}
	if( goodCVCflag ){
	  int seg = ChargedHit-> seg();
	  confMan-> GetGeomMapManager()-> GetGParam( CID_CVC, gpos, gdir );
	  confMan-> GetGeomMapManager()-> GetParam( CID_CVC, seg, pos, dir, len, wid, th, lv );
	  double calc_hitpos_z = gpos.z() + pos.z() - 0.5*th;
	  double calc_hitpos_x = calc_out_dx*calc_hitpos_z + calc_out_x0;
	  double calc_hitpos_y = calc_out_dy*calc_hitpos_z + calc_out_y0;
	  CalcFCHitPos.SetXYZ(calc_hitpos_x, calc_hitpos_y, calc_hitpos_z);
	}

	double in_fl = (fieldInPos-VertexPos).Mag();
	double uswk_v = 100.*Const*FC_mom/sqrt(pMass*pMass+FC_mom*FC_mom);
	double uswk_fl = out_time*uswk_v;
	double out_fl = (CalcFCHitPos-out_pos).Mag();
	FC_fl = in_fl + uswk_fl + out_fl;
	if( DumpLevel>20 ){
	  std::cout<<" FC fl : "<<in_fl<<"[cm]  "<<uswk_fl<<"[cm]  "<<out_fl<<"[cm]"<<std::endl;
	}

	double FCdE = ChargedHit-> emean();
	double T0FC_tof = ChargedHit->ctmean() - T0(0)->ctmean();
	VtxFC_tof = T0FC_tof - T0Vtx_tof;
	FC_beta = FC_fl/(100.*Const*VtxFC_tof);
	FC_mass2 = FC_mom*FC_mom*(1./(FC_beta*FC_beta)-1.);

        double FC_calc_p_beta = FC_mom/sqrt(FC_mom*FC_mom+pMass*pMass);
        double Calc_VtxFC_tof = FC_fl/(100.*Const*FC_calc_p_beta);
        FC_offset = Calc_VtxFC_tof - VtxFC_tof;

	double FC_mass = sqrt(FC_mass2);
	if( VtxFC_tof<FC_ADC_Param[2] ){
	  if( FCdE < FC_Cut_Param[1]*VtxFC_tof + FC_Cut_Param[0] ) FCADCCutflag =true;
	}
	else{
	  if( FCdE < FC_ADC_Param[3] ) FCADCCutflag =true;
	}

	bool FC_pid_flag = true;
	if( FC_PID_Param[0][0]<FC_mass && FC_mass<FC_PID_Param[0][1] ){
	  FC_pid = CDS_PiPlus;
	  FC_mom2 = piMass*FC_beta/sqrt(1.0-FC_beta*FC_beta);
	}
	else if( FC_PID_Param[1][0]<FC_mass && FC_mass<FC_PID_Param[1][1] ){
	  FC_pid = CDS_Proton;
	  FC_mom2 = pMass*FC_beta/sqrt(1.0-FC_beta*FC_beta);
	}
	else if( FC_PID_Param[2][0]<FC_mass && FC_mass<FC_PID_Param[2][1] ){
	  FC_pid = CDS_Deuteron;
	  FC_mom2 = dMass*FC_beta/sqrt(1.0-FC_beta*FC_beta);
	}
	else{
	  FC_pid = CDS_Other;
	  FC_pid_flag = false;
	}

      }
      else{
	Chargedflag = false;
      }
    }
  }
}

void AnalysisMan::CalcMM()
{
  if( beamPID==Beam_Pion || beamPID==Beam_Kaon ){
    int beam_pid;
    if( beamPID==Beam_Pion ) beam_pid = CDS_PiMinus;
    else if( beamPID==Beam_Kaon ) beam_pid = CDS_Kaon;

    TLorentzVector target_lmom(0., 0., 0., ThreeHeMass);

    if( Neutronflag ){
      double dx = (NC1stHitPos-VertexPos).x();
      double dy = (NC1stHitPos-VertexPos).y();
      double dz = (NC1stHitPos-VertexPos).z();
      TVector3 n_mom(dx/dz, dy/dz, 1.0);
      n_mom.SetMag(Neutron_mom);
      TLorentzVector N_lmom;
      N_lmom.SetVectM(n_mom, nMass);

      TLorentzVector missing_lmom = target_lmom+beam_lmom-N_lmom;
      BeamN_mm.SetMissingLmom(missing_lmom);
      BeamN_mm.SetParentLmom(beam_lmom);
      BeamN_mm.SetParentPID(beam_pid);
      BeamN_mm.SetParentLmom(target_lmom);
      BeamN_mm.SetParentPID(CDS_Helium3);
      BeamN_mm.SetSonLmom(N_lmom);
      BeamN_mm.SetSonPID(CDS_Neutron);
    }

    if( Chargedflag ){
      TVector3 fc_mom(FCInParam[2], FCInParam[3], 1.0);
      fc_mom.SetMag(FC_mom2);

      if( FC_pid==CDS_Proton ){
	TLorentzVector FP_lmom;
	FP_lmom.SetVectM(fc_mom, pMass);
	TLorentzVector missing_lmom = target_lmom+beam_lmom-FP_lmom;

	BeamP_mm.SetMissingLmom(missing_lmom);
	BeamP_mm.SetParentLmom(beam_lmom);
	BeamP_mm.SetParentPID(beam_pid);
	BeamP_mm.SetParentLmom(target_lmom);
	BeamP_mm.SetParentPID(CDS_Helium3);
	BeamP_mm.SetSonLmom(FP_lmom);
	BeamP_mm.SetSonPID(CDS_Proton);
      }
      if( FC_pid==CDS_Deuteron ){
	TLorentzVector FD_lmom;
	FD_lmom.SetVectM(fc_mom, dMass);
	TLorentzVector missing_lmom = target_lmom+beam_lmom-FD_lmom;

	BeamD_mm.SetMissingLmom(missing_lmom);
	BeamD_mm.SetParentLmom(beam_lmom);
	BeamD_mm.SetParentPID(beam_pid);
	BeamD_mm.SetParentLmom(target_lmom);
	BeamD_mm.SetParentPID(CDS_Helium3);
	BeamD_mm.SetSonLmom(FD_lmom);
	BeamD_mm.SetSonPID(CDS_Deuteron);
      }
    }
  }
}

void AnalysisMan::CalcForwardIM()
{
  if( Neutronflag ){
    double dx = (NC1stHitPos-VertexPos).x();
    double dy = (NC1stHitPos-VertexPos).y();
    double dz = (NC1stHitPos-VertexPos).z();
    TVector3 n_mom(dx/dz, dy/dz, 1.0);
    n_mom.SetMag(Neutron_mom);
    TLorentzVector N_lmom;
    N_lmom.SetVectM(n_mom, nMass);

    if( nCDSpip()==1 ){
      TrackVertex vertex_beam;
      if( cdstrackMan-> GetVertex_beam(cdsPip_ID[0], bpcID, vertex_beam) ){
	TVector3 vtx1 = vertex_beam.GetVertexPos1();
	CDSTrack *CDStrack = cdstrackMan-> Track(cdsPip_ID[0]);
	TVector3 cds_mom;
	if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	  Npipflag = true;
	  TLorentzVector cds_lmom;
	  cds_lmom.SetVectM(cds_mom, piMass);
	  TLorentzVector im_lmom = N_lmom + cds_lmom;
	  Npip_im.SetParentLmom(im_lmom);
	  Npip_im.SetSonPID(CDS_Neutron);
	  Npip_im.SetSonLmom(N_lmom);
	  Npip_im.SetSonPID(CDS_PiPlus);
	  Npip_im.SetSonLmom(cds_lmom);
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
      }
      else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
    }

    if( nCDSpim()==1 ){
      TrackVertex vertex_beam;
      if( cdstrackMan-> GetVertex_beam(cdsPim_ID[0], bpcID, vertex_beam) ){
	TVector3 vtx1 = vertex_beam.GetVertexPos1();
	CDSTrack *CDStrack = cdstrackMan-> Track(cdsPim_ID[0]);
	TVector3 cds_mom;
	if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	  Npimflag = true;
	  TLorentzVector cds_lmom;
	  cds_lmom.SetVectM(cds_mom, piMass);
	  TLorentzVector im_lmom = N_lmom + cds_lmom;
	  Npim_im.SetParentLmom(im_lmom);
	  Npim_im.SetSonPID(CDS_Neutron);
	  Npim_im.SetSonLmom(N_lmom);
	  Npim_im.SetSonPID(CDS_PiMinus);
	  Npim_im.SetSonLmom(cds_lmom);
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
      }
      else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
    }

    if( nCDSkm()==1 ){
      TrackVertex vertex_beam;
      if( cdstrackMan-> GetVertex_beam(cdsKm_ID[0], bpcID, vertex_beam) ){
	TVector3 vtx1 = vertex_beam.GetVertexPos1();
	CDSTrack *CDStrack = cdstrackMan-> Track(cdsKm_ID[0]);
	TVector3 cds_mom;
	if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	  Nkflag = true;
	  TLorentzVector cds_lmom;
	  cds_lmom.SetVectM(cds_mom, kpMass);
	  TLorentzVector im_lmom = N_lmom + cds_lmom;
	  Nk_im.SetParentLmom(im_lmom);
	  Nk_im.SetSonPID(CDS_Neutron);
	  Nk_im.SetSonLmom(N_lmom);
	  Nk_im.SetSonPID(CDS_Kaon);
	  Nk_im.SetSonLmom(cds_lmom);
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
      }
      else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
    }
  }

  if( Chargedflag ){
    TVector3 fc_mom(FCInParam[2], FCInParam[3], 1.0);
    fc_mom.SetMag(FC_mom2);
    if( FC_pid==CDS_Proton ){
      TLorentzVector FC_lmom;
      FC_lmom.SetVectM(fc_mom, pMass);

      if( nCDSpip()==1 ){
	TrackVertex vertex_beam;
	if( cdstrackMan-> GetVertex_beam(cdsPip_ID[0], bpcID, vertex_beam) ){
	  TVector3 vtx1 = vertex_beam.GetVertexPos1();
	  CDSTrack *CDStrack = cdstrackMan-> Track(cdsPip_ID[0]);
	  TVector3 cds_mom;
	  if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	    Ppipflag = true;
	    TLorentzVector cds_lmom;
	    cds_lmom.SetVectM(cds_mom, piMass);
	    TLorentzVector im_lmom = FC_lmom + cds_lmom;
	    Ppip_im.SetParentLmom(im_lmom);
	    Ppip_im.SetSonPID(CDS_Proton);
	    Ppip_im.SetSonLmom(FC_lmom);
	    Ppip_im.SetSonPID(CDS_PiPlus);
	    Ppip_im.SetSonLmom(cds_lmom);
	  }
	  else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
      }

      if( nCDSpim()==1 ){
	TrackVertex vertex_beam;
	if( cdstrackMan-> GetVertex_beam(cdsPim_ID[0], bpcID, vertex_beam) ){
	  TVector3 vtx1 = vertex_beam.GetVertexPos1();
	  CDSTrack *CDStrack = cdstrackMan-> Track(cdsPim_ID[0]);
	  TVector3 cds_mom;
	  if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	    Ppimflag = true;
	    TLorentzVector cds_lmom;
	    cds_lmom.SetVectM(cds_mom, piMass);
	    TLorentzVector im_lmom = FC_lmom + cds_lmom;
	    Ppim_im.SetParentLmom(im_lmom);
	    Ppim_im.SetSonPID(CDS_Proton);
	    Ppim_im.SetSonLmom(FC_lmom);
	    Ppim_im.SetSonPID(CDS_PiMinus);
	    Ppim_im.SetSonLmom(cds_lmom);
	  }
	  else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
      }

      if( nCDSkm()==1 ){
	TrackVertex vertex_beam;
	if( cdstrackMan-> GetVertex_beam(cdsKm_ID[0], bpcID, vertex_beam) ){
	  TVector3 vtx1 = vertex_beam.GetVertexPos1();
	  CDSTrack *CDStrack = cdstrackMan-> Track(cdsKm_ID[0]);
	  TVector3 cds_mom;
	  if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	    Pkflag = true;
	    TLorentzVector cds_lmom;
	    cds_lmom.SetVectM(cds_mom, kpMass);
	    TLorentzVector im_lmom = FC_lmom + cds_lmom;
	    Pk_im.SetParentLmom(im_lmom);
	    Pk_im.SetSonPID(CDS_Proton);
	    Pk_im.SetSonLmom(FC_lmom);
	    Pk_im.SetSonPID(CDS_Kaon);
	    Pk_im.SetSonLmom(cds_lmom);
	  }
	  else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
      }
    }

    if( FC_pid==CDS_Deuteron ){
      TLorentzVector FC_lmom;
      FC_lmom.SetVectM(fc_mom, dMass);

      if( nCDSpip()==1 ){
	TrackVertex vertex_beam;
	if( cdstrackMan-> GetVertex_beam(cdsPip_ID[0], bpcID, vertex_beam) ){
	  TVector3 vtx1 = vertex_beam.GetVertexPos1();
	  CDSTrack *CDStrack = cdstrackMan-> Track(cdsPip_ID[0]);
	  TVector3 cds_mom;
	  if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	    Dpipflag = true;
	    TLorentzVector cds_lmom;
	    cds_lmom.SetVectM(cds_mom, piMass);
	    TLorentzVector im_lmom = FC_lmom + cds_lmom;
	    Dpip_im.SetParentLmom(im_lmom);
	    Dpip_im.SetSonPID(CDS_Deuteron);
	    Dpip_im.SetSonLmom(FC_lmom);
	    Dpip_im.SetSonPID(CDS_PiPlus);
	    Dpip_im.SetSonLmom(cds_lmom);
	  }
	  else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
      }

      if( nCDSpim()==1 ){
	TrackVertex vertex_beam;
	if( cdstrackMan-> GetVertex_beam(cdsPim_ID[0], bpcID, vertex_beam) ){
	  TVector3 vtx1 = vertex_beam.GetVertexPos1();
	  CDSTrack *CDStrack = cdstrackMan-> Track(cdsPim_ID[0]);
	  TVector3 cds_mom;
	  if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	    Dpimflag = true;
	    TLorentzVector cds_lmom;
	    cds_lmom.SetVectM(cds_mom, piMass);
	    TLorentzVector im_lmom = FC_lmom + cds_lmom;
	    Dpim_im.SetParentLmom(im_lmom);
	    Dpim_im.SetSonPID(CDS_Deuteron);
	    Dpim_im.SetSonLmom(FC_lmom);
	    Dpim_im.SetSonPID(CDS_PiMinus);
	    Dpim_im.SetSonLmom(cds_lmom);
	  }
	  else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
      }

      if( nCDSkm()==1 ){
	TrackVertex vertex_beam;
	if( cdstrackMan-> GetVertex_beam(cdsKm_ID[0], bpcID, vertex_beam) ){
	  TVector3 vtx1 = vertex_beam.GetVertexPos1();
	  CDSTrack *CDStrack = cdstrackMan-> Track(cdsKm_ID[0]);
	  TVector3 cds_mom;
	  if( CDStrack->GetMomentum(vtx1, cds_mom) ){
	    Dkflag = true;
	    TLorentzVector cds_lmom;
	    cds_lmom.SetVectM(cds_mom, kpMass);
	    TLorentzVector im_lmom = FC_lmom + cds_lmom;
	    Dk_im.SetParentLmom(im_lmom);
	    Dk_im.SetSonPID(CDS_Deuteron);
	    Dk_im.SetSonLmom(FC_lmom);
	    Dk_im.SetSonPID(CDS_Kaon);
	    Dk_im.SetSonLmom(cds_lmom);
	  }
	  else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetMomentum fault n & pi+"<<std::endl;
	}
	else std::cout<<" !!! AnalysisMan::CalcForwardIM() GetVertex_beam fault n &pi+"<<std::endl;
      }
    }
  }
}    
  

void AnalysisMan::Clear()
{
//***** for Clear Classes ******//
  if( beamSpec!=0 )  beamSpec-> Clear();
  if( fieldMan!=0 )  fieldMan-> Clear();

//***** for Event Flag ******//
  goodNCflag = false;
  Neutronflag = false;
  Gammaflag = false;
  Chargedflag = false;
  goodPCflag = false;
  goodCVCflag = false;

  VtxCutflag = false;
  FCADCCutflag = false;

  Npipflag = false;
  Npimflag = false;
  Nkflag = false;

  Ppimflag = false;
  Ppipflag = false;
  Pkflag = false;

  Dpimflag = false;
  Dpipflag = false;
  Dkflag = false;

//***** for Analysis Data *****//
  beamPID = -999;
  beam_lmom.SetXYZM(0.,0.,0.,0.);

  blc1ID = -999;
  blc2ID = -999;
  bpcID  = -999;

  BHDHit = 0;
  T0Hit = 0;
  NC1st = 0;
  ChargedHit = 0;

  D5Mom = -999.;
  T0Vtx_tof = -999.;

//***** for Hit Position *****//
  T0HitPosBLC.SetXYZ(-9999., -9999., -9999.);
  T0HitPosBPC.SetXYZ(-9999., -9999., -9999.);

  BPDHitPosBLC.SetXYZ(-9999., -9999., -9999.);
  BPDHitPosBPC.SetXYZ(-9999., -9999., -9999.);

  BPCHitPos.SetXYZ(-9999., -9999., -9999);
  VertexPos.SetXYZ(-9999., -9999., -9999.);
  FDC1HitPos.SetXYZ(-9999., -9999., -9999);
  USWKHitPos.SetXYZ(-9999., -9999., -9999);

  NC1stHitPos.SetXYZ(-9999., -9999., -9999.);
  NCHitPos.clear();
  FCHitPos.SetXYZ(-9999., -9999., -9999.);
  CalcFCHitPos.SetXYZ(-9999., -9999., -9999.);

//***** for Forward Counter *****//
  FC_pid = -999;
  T0NC_tof.clear();

  VtxNC1st_tof = -999.;
  NC1st_beta = -999.;
  Neutron_mom = -999;
  VtxNC_tof.clear();
  NC_beta.clear();
  NC_offset.clear();

  for( int i=0; i<5; i++ ){
    FCInParam[i] = -999.;
  }
  FC_fl = -999;
  FC_mom = -999.;
  FC_mom2 = -999.;
  FC_offset = -9999.;

  T0PC_tof = -999.;
  T0CVC_tof = -999.;
  VtxFC_tof = -999.;

  BeamN_mm.Clear();
  BeamP_mm.Clear();
  BeamD_mm.Clear();

  Npip_im.Clear();
  Npim_im.Clear();
  Nk_im.Clear();

  Ppip_im.Clear();
  Ppim_im.Clear();
  Pk_im.Clear();

  Dpip_im.Clear();
  Dpim_im.Clear();
  Dk_im.Clear();

//***** for CDS data *****//
  vertexBeam.Clear();
  GoodTrack_ID.clear();
  CDS_beta.clear();
  CDS_mass2.clear();
  CDH_offset.clear();

  cdsPim_ID.clear();
  cdsPip_ID.clear();
  cdsKm_ID.clear();
  cdsKp_ID.clear();
  cdsP_ID.clear();
  cdsD_ID.clear();
  cdsLow_ID.clear();
  cdsHigh_ID.clear();

  cdsPiPiIM.clear();
  cdsPiPIM.clear();
  cdsKPIM.clear();

//***** for Hodoscope *****//
  BHD_ID.clear();
  BHDpost_ID.clear();
  T0_ID.clear();
  T0pre_ID.clear();
  T0post_ID.clear();
  E0_ID.clear();
  CVC_ID.clear();
  BPD_ID.clear();
  BVC_ID.clear();
  PC_ID.clear();
  BD_ID.clear();
  LB_ID.clear();
  WVC_ID.clear();
  HVC1_ID.clear();
  HVC2_ID.clear();
  Temp1_ID.clear();
  Temp2_ID.clear();
  for( int i=0; i<NumOfNCLayer; i++ ){
    NC_ID[i].clear();
  }

  CDH_ID.clear();
  IH_ID.clear();
}


bool AnalysisMan::Initialize()
{
  std::cout<<"====================================="<<std::endl;
  std::cout<<"=== AnalysisMan::Initialize start ==="<<std::endl;
  std::cout<<"====================================="<<std::endl;

  if( confMan==0 ) exit(-1);
  if( scaMan==0 ) exit(-1);
  if( header==0 ) exit(-1);
  if( cdsMan==0 ) exit(-1);
  if( blMan==0 ) exit(-1);
  if( bltrackMan==0 ) exit(-1);

  if( D5AnaFlag ){
    std::cout<<" > D5 momentum Analysis is 'On'"<<std::endl;
    beamSpec = new BeamSpectrometer(confMan);
  }
  else std::cout<<" > D5 momentum Analysis is 'Off'"<<std::endl;

  if( CDSAnaMode>0 ){
    if( cdstrackMan==0 ) exit(-1);
    if( CDSAnaMode==1 )   std::cout<<" > CDS Tracking is 'All Event'"<<std::endl; 
    if( CDSAnaMode==2 )   std::cout<<" > CDS Tracking is 'Forward Hit Event'"<<std::endl; 
    if( CDSAnaMode==3 )   std::cout<<" > CDS Tracking is 'NC hit Event'"<<std::endl; 
    if( CDSAnaMode==4 )   std::cout<<" > CDS Tracking is 'PC or CVC Event'"<<std::endl; 
    if( CDSAnaMode==5 )   std::cout<<" > CDS Tracking is 'T0 1hit & BPC track Event'"<<std::endl; 
    if( CDSAnaMode==6 )   std::cout<<" > CDS Tracking is 'T0 1hit & BPC track & Forward hit Event'"<<std::endl; 
    if( CDSAnaMode==7 )   std::cout<<" > CDS Tracking is 'T0 1hit & BPC track & NC hit Event'"<<std::endl; 
    if( CDSAnaMode==8 )   std::cout<<" > CDS Tracking is 'T0 1hit & BPC track & CVC or PC hit Event'"<<std::endl; 
    if( CDSAnaMode==9 )   std::cout<<" > CDS Tracking is 'T0 1hit & BPC track & Neutral Event'"<<std::endl; 
    if( CDSAnaMode==10 )  std::cout<<" > CDS Tracking is 'T0 1hit & BPC track & Charged Event'"<<std::endl; 
    if( CDSAnaMode>10  ){
      std::cout<<" > CDS Tracking is 'unknown mode'"<<std::endl;;
      exit(-1);
    }
  }
  else{
    std::cout<<" > CDS Tracking is 'Off'"<<std::endl;
  }

  if( FCAnaFlag ){
    std::cout<<" > Forward Charged Analysis is 'On'"<<std::endl;

    fieldMan = new UshiwakaFieldMapMan();
    fieldMan-> Initialize();
    tableMan = new UshiwakaTableMan();
    tableMan-> Initialize();

    double PCz1, PCx1, PCz2, PCx2;
    TVector3 pos, dir, gpos, gdir;
    double len, wid, th, lv;
    confMan-> GetGeomMapManager()-> GetGParam( CID_PC, gpos, gdir );
    confMan-> GetGeomMapManager()-> GetParam(CID_PC, 1, pos, dir, len, wid ,th, lv);
    PCx1 = gpos.x() + cos(PI*gdir.y()/180.)*pos.x() + sin(PI*gdir.y()/180.)*(pos.z()-0.5*th);
    PCz1 = gpos.z() - sin(PI*gdir.y()/180.)*pos.x() + cos(PI*gdir.y()/180.)*(pos.z()-0.5*th);

    confMan-> GetGeomMapManager()-> GetParam(CID_PC, NumOfPCSegments, pos, dir, len, wid ,th, lv);
    PCx2 = gpos.x() + cos(PI*gdir.y()/180.)*pos.x() + sin(PI*gdir.y()/180.)*(pos.z()-0.5*th);
    PCz2 = gpos.z() - sin(PI*gdir.y()/180.)*pos.x() + cos(PI*gdir.y()/180.)*(pos.z()-0.5*th);

    PC_ZX_Param[1] = (PCx1-PCx2)/(PCz1-PCz2);
    PC_ZX_Param[0] = PCx1 - PC_ZX_Param[1]*PCz1;

    double center[3];
    double range[3];
    fieldMan-> GetGPOS(center[0], center[1], center[2]);
    fieldMan-> GetRange(range[0], range[1], range[2]);
    USWK_Z = center[2];
    FieldInZ = center[2]-range[2];
  }
  else std::cout<<" > Forward Charged Analysis is 'Off'"<<std::endl;

  FILE *fp;
  char str[MAXCHAR];
  double param1, param2, param3, param4;
  if( (fp=fopen(ParamFileName.c_str(),"r"))==0 ){
    std::cerr << " File open fail. [" << ParamFileName << "]" << std::endl;
    exit(-1);
  }

  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;

    if( sscanf(str, "BHDT0 pi: %lf %lf",&param1, &param2)==2 ){
      BHDT0_Param[0][0] = param1;
      BHDT0_Param[0][1] = param2;
    }
    if( sscanf(str, "BHDT0 k : %lf %lf",&param1, &param2)==2 ){
      BHDT0_Param[1][0] = param1;
      BHDT0_Param[1][1] = param2;
    }
    if( sscanf(str, "BHDT0 p : %lf %lf",&param1, &param2)==2 ){
      BHDT0_Param[2][0] = param1;
      BHDT0_Param[2][1] = param2;
    }

    if( sscanf(str, "Vtx ZX: %lf %lf",&param1, &param2)==2 ){
      Vtx_ZX_Param[0] = param1;
      Vtx_ZX_Param[1] = param2;
    }
    if( sscanf(str, "Vtx ZY: %lf %lf",&param1, &param2)==2 ){
      Vtx_ZY_Param[0] = param1;
      Vtx_ZY_Param[1] = param2;
    }
    if( sscanf(str, "Vtx R2: %lf",&param1)==1 ){
      Vtx_R2_Param = param1;
    }
    if( sscanf(str, "Vtx Z : %lf %lf",&param1, &param2)==2 ){
      Vtx_Z_Param[0] = param1;
      Vtx_Z_Param[1] = param2;
    }

    if( sscanf(str, "CDS pi: %lf %lf",&param1, &param2)==2 ){
      CDS_PID_Param[0][0] = param1;
      CDS_PID_Param[0][1] = param2;
    }
    if( sscanf(str, "CDS k : %lf %lf",&param1, &param2)==2 ){
      CDS_PID_Param[1][0] = param1;
      CDS_PID_Param[1][1] = param2;
    }
    if( sscanf(str, "CDS p : %lf %lf",&param1, &param2)==2 ){
      CDS_PID_Param[2][0] = param1;
      CDS_PID_Param[2][1] = param2;
    }
    if( sscanf(str, "CDS d : %lf %lf",&param1, &param2)==2 ){
      CDS_PID_Param[3][0] = param1;
      CDS_PID_Param[3][1] = param2;
    }

    if( sscanf(str, "FC pi: %lf %lf",&param1, &param2)==2 ){
      FC_PID_Param[0][0] = param1;
      FC_PID_Param[0][1] = param2;
    }
    if( sscanf(str, "FC p : %lf %lf",&param1, &param2)==2 ){
      FC_PID_Param[1][0] = param1;
      FC_PID_Param[1][1] = param2;
    }
    if( sscanf(str, "FC d : %lf %lf",&param1, &param2)==2 ){
      FC_PID_Param[2][0] = param1;
      FC_PID_Param[2][1] = param2;
    }

    if( sscanf(str, "FC cut: %lf %lf %lf %lf",&param1, &param2, &param3, &param4)==4 ){
      FC_ADC_Param[0] = param1, FC_ADC_Param[1] = param2, FC_ADC_Param[2] = param3, FC_ADC_Param[3] = param4;
      FC_Cut_Param[1] = (param2-param4)/(param1-param3);
      FC_Cut_Param[0] = param2 - FC_Cut_Param[1]*param1;
    }

    if( sscanf(str, "CDS lambda: %lf %lf",&param1, &param2)==2 ){
      CDS_Lambda_Param[0] = param1;
      CDS_Lambda_Param[1] = param2;
    }
    if( sscanf(str, "CDS k0: %lf %lf",&param1, &param2)==2 ){
      CDS_K0_Param[0] = param1;
      CDS_K0_Param[1] = param2;
    }
    if( sscanf(str, "FC lambda: %lf %lf",&param1, &param2)==2 ){
      FC_Lambda_Param[0] = param1;
      FC_Lambda_Param[1] = param2;
    }
    if( sscanf(str, "FC lambda(1520): %lf %lf",&param1, &param2)==2 ){
      FC_Lambda1520_Param[0] = param1;
      FC_Lambda1520_Param[1] = param2;
    }
  }

  DumpParam();
  std::cout<<"==================================="<<std::endl;
  std::cout<<"=== AnalysisMan::Initialize end ==="<<std::endl;
  std::cout<<"==================================="<<std::endl;
}

bool AnalysisMan::Finalize()
{
}

int AnalysisMan::nNC() const
{
  int nhit = 0;
  for( int layer=1; layer<=NumOfNCLayer; layer++ ){
    nhit += nNC(layer);
  }
  return nhit;
}

HodoscopeLikeHit* AnalysisMan::NC(const int &i)
{
  int tmp = i;
  for( int layer=1; layer<=NumOfNCLayer; layer++ ){
    if( tmp<nNC(layer) ){
      return blMan->NC(NC_ID[layer-1][tmp]);
    }
    else{
      tmp -= nNC(layer);
    }
  }
  std::cout<<" AnalysisMan::NC("<<i<<") out of size  return 0 !!!"<<std::endl;
  return 0;
}

void AnalysisMan::SetHodo()
{
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      BHD_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nBHDpost(); i++ ){
    if( blMan->BHDpost(i)->CheckRange() ){
      BHDpost_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      T0_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nT0pre(); i++ ){
    if( blMan->T0pre(i)->CheckRange() ){
      T0pre_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nT0post(); i++ ){
    if( blMan->T0post(i)->CheckRange() ){
      T0post_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nE0(); i++ ){
    if( blMan->E0(i)->CheckRange() ){
      E0_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nCVC(); i++ ){
    if( blMan->CVC(i)->CheckRange() ){
      CVC_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nBPD(); i++ ){
    if( blMan->BPD(i)->CheckRange() ){
      BPD_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nBVC(); i++ ){
    if( blMan->BVC(i)->CheckRange() ){
      BVC_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nPC(); i++ ){
    if( blMan->PC(i)->CheckRange() ){
      PC_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nBD(); i++ ){
    if( blMan->BD(i)->CheckRange() ){
      BD_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nLB(); i++ ){
    if( blMan->LB(i)->CheckRange() ){
      LB_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nHVC1(); i++ ){
    if( blMan->HVC1(i)->CheckRange() ){
      HVC1_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nHVC2(); i++ ){
    if( blMan->HVC2(i)->CheckRange() ){
      HVC2_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nTemp1(); i++ ){
    if( blMan->Temp1(i)->CheckRange() ){
      Temp1_ID.push_back(i);
    }
  }
  for( int i=0; i<blMan->nTemp2(); i++ ){
    if( blMan->Temp2(i)->CheckRange() ){
      Temp2_ID.push_back(i);
    }
  }

  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->CheckRange() ){
      int seg = blMan->NC(i)->seg();
      int id = (seg-1)/NumOfNCSegmentsInLayer;
      NC_ID[id].push_back(i);
    }
  }
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    if( cdsMan->CDH(i)-> CheckRange() ){
      CDH_ID.push_back(i);
    }
  }

  // for Single readout
  for( int i=0; i<blMan->nWVC(); i++ ){
    int tdc = blMan->WVC(i)-> tdcu();
    if( 0<tdc && tdc<4000 ){
      WVC_ID.push_back(i);
    }
  }

  for( int i=0; i<cdsMan->nIH(); i++ ){
    int tdc = cdsMan->IH(i)-> tdcu();
    if( 0<tdc && tdc<4000 ){
      IH_ID.push_back(i);
    }
  }
}

void AnalysisMan::DumpStatus()
{
  int event_number = header->ev();
  int block_event_number = header->blev();
  std::cout<<"=== Event Number : "<<event_number<<" Dumping Dump Level : "<<DumpLevel<<" ==="<<std::endl;

  if( DumpLevel>4 ){
    DumpHit();
    DumpBLTrack();
    DumpCDSTrack();
    DumpForward();
  }

  if( BHDT0() ){
    std::cout<<" > BHD-T0 TOF : "<<BHDT0_TOF()<<" [ns]";
    if( beamPID==Beam_Pion )   std::cout<<" Pion Beam "<<std::endl;
    if( beamPID==Beam_Kaon )   std::cout<<" Kaon Beam "<<std::endl;
    if( beamPID==Beam_Proton ) std::cout<<" Proton Beam "<<std::endl;
    if( beamPID==Beam_Other )  std::cout<<" nuknown Beam "<<std::endl;
  }
  else{
    std::cout<<" > BHD : "<<nBHD()<<"  &  T0 : "<<nT0()<<std::endl;
  }

  if( D5() ){
    std::cout<<" > D5 Momentum : "<<D5Mom<<"[GeV/c]  Chi-Square : "<<beamSpec-> chisquare()<<std::endl;
  }
  else{
    std::cout<<" > BLC1 : "<<bltrackMan->ntrackBLC1()<<"  &  BLC2 : "<<bltrackMan->ntrackBLC2()<<std::endl;
  }
}

void AnalysisMan::DumpHit()
{
  std::cout<<" === Beam Line Counter Hit "<<std::endl;
  std::cout<<"   > BHD : "<<nBHD()<<std::endl;
  std::cout<<"   > T0  : "<<nT0()<<std::endl;
  std::cout<<"   > DEF : "<<nDEF()<<std::endl;
  std::cout<<"   > BPD : "<<nBPD()<<std::endl;
  std::cout<<std::endl;
  std::cout<<" === CDS Counter Hit "<<std::endl;
  std::cout<<"   > IH  ; "<<nIH()<<std::endl;
  std::cout<<"   > CDH ; "<<nCDH()<<std::endl;
  std::cout<<std::endl;
  std::cout<<" === Forward Counter Hit "<<std::endl;
  std::cout<<"   > BVC : "<<nBVC()<<std::endl;
  std::cout<<"   > CVC : "<<nCVC()<<std::endl;
  std::cout<<"   > PC  : "<<nPC()<<std::endl;
  std::cout<<"   > NC  : "<<nNC()<<std::endl;
  for( int layer=1; layer<=NumOfNCLayer; layer++ )
    std::cout<<"     > NC layer"<<layer<<" : "<<nNC(layer)<<std::endl;
  std::cout<<"   > LB  : "<<nLB()<<std::endl;
  std::cout<<"   > BD  : "<<nBD()<<std::endl;
  std::cout<<std::endl;
  std::cout<<" === Veto Counter Hit "<<std::endl;
  std::cout<<"   > HVC1: "<<nHVC1()<<std::endl;
  std::cout<<"   > HVC2: "<<nHVC2()<<std::endl;
  std::cout<<"   > WVC : "<<nWVC()<<std::endl;
  std::cout<<std::endl;
  std::cout<<" === pile up Hit "<<std::endl;
  std::cout<<"   > BHD post : "<<nBHDpost()<<std::endl;
  std::cout<<"   > T0 pre   : "<<nT0pre()<<std::endl;
  std::cout<<"   > T0 post  : "<<nT0post()<<std::endl;
}

void AnalysisMan::DumpBLTrack()
{
  std::cout<<" === BLC1 ==="<<std::endl;
  std::cout<<"   > BLC1  : "<<bltrackMan->ntrackBLC1()<<std::endl;
  std::cout<<"    > BLC1a : "<<bltrackMan->ntrackBLC1a()<<std::endl;
  std::cout<<"    > BLC1b : "<<bltrackMan->ntrackBLC1b()<<std::endl;
  std::cout<<" === BLC2 ==="<<std::endl;
  std::cout<<"   > BLC2  : "<<bltrackMan->ntrackBLC2()<<std::endl;
  std::cout<<"    > BLC2a : "<<bltrackMan->ntrackBLC2a()<<std::endl;
  std::cout<<"    > BLC2b : "<<bltrackMan->ntrackBLC2b()<<std::endl;
  std::cout<<" === BPC & FDC1 === "<<std::endl;
  std::cout<<"   > BPC   : "<<bltrackMan->ntrackBPC()<<std::endl;
  std::cout<<"   > FDC1  : "<<bltrackMan->ntrackFDC1()<<std::endl;
}

void AnalysisMan::DumpHitPos()
{
  if( IsT0HitPosBLC() ) std::cout<<"   > T0 Hit Position by BLC2 : ("<<T0HitPosBLC.x()<<" ,"<<T0HitPosBLC.y()<<" ,"<<T0HitPosBLC.z()<<") [cm]"<<std::endl;
  else std::cout<<"   > T0 Hit Position by BLC2 not dicided"<<std::endl;
  if( IsT0HitPosBLC() ) std::cout<<"   > T0 Hit Position by BPC  : ("<<T0HitPosBPC.x()<<" ,"<<T0HitPosBPC.y()<<" ,"<<T0HitPosBPC.z()<<") [cm]"<<std::endl;
  else std::cout<<"   > T0 Hit Position by BPC  not dicided"<<std::endl;

  if( IsBPDHitPosBLC() ) std::cout<<"   > BPD Hit Position by BLC2 : ("<<BPDHitPosBLC.x()<<" ,"<<BPDHitPosBLC.y()<<" ,"<<BPDHitPosBLC.z()<<") [cm]"<<std::endl;
  else std::cout<<"   > BPD Hit Position by BLC2 not dicided"<<std::endl;
  if( IsBPDHitPosBLC() ) std::cout<<"   > BPD Hit Position by BPC  : ("<<BPDHitPosBPC.x()<<" ,"<<BPDHitPosBPC.y()<<" ,"<<BPDHitPosBPC.z()<<") [cm]"<<std::endl;
  else std::cout<<"   > BPD Hit Position by BPC  not dicided"<<std::endl;

  if( bltrackMan->ntrackBPC()>0 ) std::cout<<"   > BPC Hit Position : ("<<BPCHitPos.x()<<" ,"<<BPCHitPos.y()<<" ,"<<BPCHitPos.z()<<") [cm]"<<std::endl;
  else std::cout<<"   > BPC Hit Position not dicided"<<std::endl;

  if( IsVertexPos() ) std::cout<<"   > Vertex Position : ("<<VertexPos.x()<<" ,"<<VertexPos.y()<<" ,"<<VertexPos.z()<<") [cm]"<<std::endl;
  else std::cout<<"   > Vertex Position not decided BPC : "<<bltrackMan->ntrackBPC()<<" track & CDS "<<nGoodTrack()<<" track"<<std::endl;

  if( IsFDC1HitPos() ) std::cout<<"   > FDC1 Hit Position : ("<<FDC1HitPos.x()<<" ,"<<FDC1HitPos.y()<<" ,"<<FDC1HitPos.z()<<") [cm]"<<std::endl;
  else std::cout<<"   > FDC1 Hit Position not decided"<<std::endl;

  if( goodNCflag ) std::cout<<"   > NC1st Hit Position : ("<<NC1stHitPos.x()<<" ,"<<NC1stHitPos.y()<<" ,"<<NC1stHitPos.z()<<") [cm]"<<std::endl;
  if( Chargedflag ) std::cout<<"   > Charged Hit Position : ("<<FCHitPos.x()<<" ,"<<FCHitPos.y()<<" ,"<<FCHitPos.z()<<") [cm]"<<std::endl;
}

void AnalysisMan::DumpCDSTrack()
{
  std::cout<<" === CDS Track Dumping ==="<<std::endl;
  std::cout<<"   > CDS Good Track : "<<cdstrackMan->nGoodTrack()<<std::endl;
  if( nT0()==1 && bpcID>=0 ){
    std::cout<<"   > CDS Good Track w CDH : "<<nGoodTrack()<<std::endl;
    for( int i=0; i<nGoodTrack(); i++ ){
      CDSTrack *cdsTrack = CDSGoodTrack(i);
      std::cout<<"    > CDH seg : "<<cdsTrack-> CDHHit(cdsMan)-> seg()<<"  mom : "<<cdsTrack->Momentum()<<"[GeV/c]";
      if( cdsTrack-> PID()==CDS_PiPlus )   std::cout<<"  pi+"<<std::endl;
      if( cdsTrack-> PID()==CDS_Proton )   std::cout<<"  proton"<<std::endl;
      if( cdsTrack-> PID()==CDS_Deuteron ) std::cout<<"  deuteron"<<std::endl;
      if( cdsTrack-> PID()==CDS_Triton )   std::cout<<"  triton"<<std::endl;
      if( cdsTrack-> PID()==CDS_PiMinus )  std::cout<<"  pi-"<<std::endl;
      if( cdsTrack-> PID()==CDS_Kaon )     std::cout<<"  K-"<<std::endl;
      if( cdsTrack-> PID()==CDS_Other )    std::cout<<"  nuknown"<<std::endl;
    }
  }
  else{
    std::cout<<"   > nT0 : "<<nT0()<<"  BPC ID : "<<bpcID<<std::endl;
  }
}

void AnalysisMan::DumpForward()
{
  if( goodNCflag ){
    std::cout<<" === Good Nenutron Counter Hit Event  ==="<<std::endl;
    std::cout<<"   > Vtx-NC 1st hit  TOF : "<<VtxNC1st_tof<<"[ns]  beta : "<<NC1st_beta;
    if( Neutronflag ) std::cout<<"  Neutron event"<<std::endl;
    if( Gammaflag ) std::cout<<"  Gamma event"<<std::endl;
    for( int i=0; i<VtxNC_tof.size(); i++ ){
      std::cout<<"     > Vtx-NC TOF : "<<VtxNC_tof[i]<<"[ns]  beta : "<<NC_beta[i]<<std::endl;
    }
  }

  if( Chargedflag ){
    std::cout<<" === Forward Charged Event ==="<<std::endl;
    std::cout<<"   > Field In PAram (  InPos x  ,  InPos y  ,  In dx/dz  ,  In dy/dz  ,  Bending Angle  )"<<std::endl;
    std::cout<<"                    ("<<FCInParam[0]<<" ,"<<FCInParam[1]<<" ,"<<FCInParam[2]<<" ,"<<FCInParam[3]<<" ,"<<360.*FCInParam[4]/PI<<")"<<std::endl;
    std::cout<<"                    momentum : "<<FC_mom<<"[GeV/c]  flight length : "<<FC_fl<<"[cm]"<<std::endl;
    std::cout<<"   > Vtx-Charged TOF : "<<VtxFC_tof<<"[ns]  Forward Charged beta : "<<FC_beta<<std::endl;
    std::cout<<"     > mass : "<<sqrt(FC_mass2);
    if( FC_pid==CDS_PiPlus ) std::cout<<"  pi+"<<std::endl;
    if( FC_pid==CDS_Proton ) std::cout<<"  proton  mom : "<<FC_mom2<<"[GeV/c]"<<std::endl;
    if( FC_pid==CDS_Deuteron ) std::cout<<"  deuteron  mom : "<<FC_mom2<<"[GeV/c]"<<std::endl;
    if( FC_pid==CDS_Other ) std::cout<<"  unknown"<<std::endl;
    if( DumpLevel>20 ){
      double x_diff = FCHitPos.x()-CalcFCHitPos.x();
      double z_diff = FCHitPos.z()-CalcFCHitPos.z();
      double pos_diff = sqrt( x_diff*x_diff + z_diff*z_diff );
      if( FCHitPos.x()>CalcFCHitPos.x() )   pos_diff *= -1.0;

      std::cout<<"   > FC   Hit Pos ("<<FCHitPos.x()<<", "<<FCHitPos.y()<<", "<<FCHitPos.z()<<")"<<std::endl;
      std::cout<<"   > Calc Hit Pos ("<<CalcFCHitPos.x()<<", "<<CalcFCHitPos.y()<<", "<<CalcFCHitPos.z()<<")"<<std::endl;
      std::cout<<"     > XZ plane diff : "<<pos_diff<<"[cm]"<<std::endl;
    }
  }
}

void AnalysisMan::DumpParam()
{
  std::ios::fmtflags flagsSaved = std::cout.flags();
  std::cout<<"##### AnalysisMan::DumpParameter #####"<<std::endl;
  std::cout<<"# BHDT0-TOF Param"<<std::endl;
  std::cout<<"> pi ("<<BHDT0_Param[0][0]<<" : "<<BHDT0_Param[0][1]<<") [ns]"<<std::endl;
  std::cout<<"> k  ("<<BHDT0_Param[1][0]<<" : "<<BHDT0_Param[1][1]<<") [ns]"<<std::endl;
  std::cout<<"> p  ("<<BHDT0_Param[2][0]<<" : "<<BHDT0_Param[2][1]<<") [ns]"<<std::endl;
  std::cout<<"# Vtx Cut Param"<<std::endl;
  std::cout<<"> X center : "<<Vtx_ZX_Param[1]<<"*z"<<std::showpos<<Vtx_ZX_Param[0]<<" [cm]"<<std::endl;
  std::cout.flags(flagsSaved);
  std::cout<<"> Y center : "<<Vtx_ZY_Param[1]<<"*z"<<std::showpos<<Vtx_ZY_Param[0]<<" [cm]"<<std::endl;
  std::cout.flags(flagsSaved);
  std::cout<<"> Target Z ("<<Vtx_Z_Param[0]<<" : "<<Vtx_Z_Param[1]<<") [cm]"<<std::endl;
  std::cout<<"> radius : "<<sqrt(Vtx_R2_Param)<<" [cm]"<<std::endl;
  std::cout<<"# CDS pid parameter"<<std::endl;
  std::cout<<"> pi ("<<CDS_PID_Param[0][0]<<" : "<<CDS_PID_Param[0][1]<<") [GeV]"<<std::endl;
  std::cout<<"> k  ("<<CDS_PID_Param[1][0]<<" : "<<CDS_PID_Param[1][1]<<") [GeV]"<<std::endl;
  std::cout<<"> p  ("<<CDS_PID_Param[2][0]<<" : "<<CDS_PID_Param[2][1]<<") [GeV]"<<std::endl;
  std::cout<<"> d  ("<<CDS_PID_Param[3][0]<<" : "<<CDS_PID_Param[3][1]<<") [GeV]"<<std::endl;
  std::cout<<"# FC pid parameter"<<std::endl;
  std::cout<<"> pi ("<<FC_PID_Param[0][0]<<" : "<<FC_PID_Param[0][1]<<") [GeV]"<<std::endl;
  std::cout<<"> p  ("<<FC_PID_Param[1][0]<<" : "<<FC_PID_Param[1][1]<<") [GeV]"<<std::endl;
  std::cout<<"> d  ("<<FC_PID_Param[2][0]<<" : "<<FC_PID_Param[2][1]<<") [GeV]"<<std::endl;
  std::cout<<"# FC Cut parameter"<<std::endl;
  std::cout<<"> dE = "<<FC_Cut_Param[1]<<"*tof"<<std::showpos<<FC_Cut_Param[0]<<" [MeV]"<<std::endl;
  std::cout<<"# FC XZ Parameter"<<std::endl;
  std::cout<<"> X = "<<PC_ZX_Param[1]<<"*z"<<std::showpos<<PC_ZX_Param[0]<<std::endl;
  std::cout.flags(flagsSaved);
  std::cout<<"> (tof[ns], dE[MeV]) : ("<<FC_ADC_Param[0]<<" ,"<<FC_ADC_Param[1]<<")  ("<<FC_ADC_Param[2]<<" ,"<<FC_ADC_Param[3]<<")"<<std::endl;
  std::cout<<"# Invariant Mass Cut"<<std::endl;
  std::cout<<"> CDS Lambda ("<<CDS_Lambda_Param[0]<<" : "<<CDS_Lambda_Param[1]<<") [GeV]"<<std::endl;
  std::cout<<"> CDS K0     ("<<CDS_K0_Param[0]<<" : "<<CDS_K0_Param[1]<<") [GeV]"<<std::endl;
  std::cout<<"> FC Lambda ("<<FC_Lambda_Param[0]<<" : "<<FC_Lambda_Param[1]<<") [GeV]"<<std::endl;
  std::cout<<"> FC Lambda(1520) ("<<FC_Lambda1520_Param[0]<<" : "<<FC_Lambda1520_Param[1]<<") [GeV]"<<std::endl;
}

void AnalysisMan::Wait()
{
  char dispflag;
  char dispin[100]="";
  std::cout<<" Input any word or return"<<std::endl;
  std::cout<<" (q:quite)"<<std::endl;
  fgets(dispin,100,stdin);
  if(sscanf(dispin,"%c",&dispflag)==1){
    if( dispflag=='q' ){
      exit(-1);
    }
  }
}
