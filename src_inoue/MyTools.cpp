#include "MyTools.h"

using namespace std;

double MyTools::cosOP(CDS2Info *info){
  TVector3 p1=info->momentum1(), p2=info->momentum2();
  return p1.Dot(p2)/(p1.Mag()*p2.Mag());
}

double MyTools::calc_offset(const double &mm, const TLorentzVector tot_lmom, const TLorentzVector lmom, AnaInfo *anaInfo)
{
  TLorentzVector fit_lmom=lmom;
  double fit_mm=(tot_lmom-lmom).M();
  double delta=fit_mm-mm;

  int count=0;
  while( true ){
    TVector3 mom=fit_lmom.Vect();
    mom.SetMag(mom.Mag()+delta);
    fit_lmom.SetVectM(mom, lmom.M());
    fit_mm=(tot_lmom-fit_lmom).M();
    delta=fit_mm-mm;

    count++;
    if( count>100 ) return -99999;
    if( fabs(delta)<0.0001 ) break;
  }

  TVector3 vertex=anaInfo->minDCA()->vertexBeam();
  TVector3 FDC1pos=anaInfo->forwardCharge(0)->FDC1().pos();
  double fl=anaInfo->forwardCharge(0)-> fl();
  double tmp_mom, tof1, tof2;
  ELossTools::CalcElossForwardTGeo(vertex, FDC1pos, fl, lmom.Vect().Mag(), lmom.M(), tmp_mom, tof1);
  ELossTools::CalcElossForwardTGeo(vertex, FDC1pos, fl, fit_lmom.Vect().Mag(), fit_lmom.M(), tmp_mom, tof2);

  double mom_out, tmp_tof;

  return tof1-tof2;
}

bool MyTools::searchCDHHit(CDSTrack *track, CDSHitMan *cdsMan, bool z_match)
{
  double param[5];
  track-> ClearCDH();
  track->GetParameters(param);
  TVector3 vtxCDH=MathTools::CalcHelixPosatR(param, 54.4);
  track->SetCDHVertex(vtxCDH); 

  vector<HodoscopeLikeHit*> CDHhits=MyTools::getCDH(cdsMan);
  double min_dphi=DBL_MAX;
  HodoscopeLikeHit *cand=0;
  for( int i=0; i<(int)CDHhits.size(); i++ ){
    double dphi=fabs(dphiCDH(track, CDHhits[i]));
    if( dphi<min_dphi ){
      min_dphi=dphi;
      cand=CDHhits[i];
    }
  }
  if( min_dphi>7.0 ) return false;
  if( !z_match ){
    track-> SetHodoHit(*cand, vtxCDH);
    return true;
  }
  
  double min_z_diff=fabs(vtxCDH.Z()-cand->hitpos());
  int seg=cand->seg();
  for( int i=0; i<(int)CDHhits.size(); i++ ){
    if( fabs(seg-CDHhits[i]->seg() )==1 ){
      double z_diff=fabs(vtxCDH.Z()-CDHhits[i]->hitpos());
      if( z_diff<min_z_diff ){
	min_z_diff=z_diff;
	cand=CDHhits[i];
      }
    }
  }

  track-> SetHodoHit(*cand, vtxCDH);
  return true;
}

bool MyTools::searchIHHit(CDSTrack *track, CDSHitMan *cdsMan)
{
  double param[5];
  track-> ClearIH();
  track->GetParameters(param);
  TVector3 vtxIH=MathTools::CalcHelixPosatR(param, 13.7);
  track->SetIHVertex(vtxIH); 

  vector<HodoscopeLikeHit*> IHhits=MyTools::getIH(cdsMan);
  double min_dphi=DBL_MAX;
  HodoscopeLikeHit *cand=0;
  for( int i=0; i<(int)IHhits.size(); i++ ){
    double dphi=fabs(dphiIH(track, IHhits[i]));
    if( dphi<min_dphi ){
      min_dphi=dphi;
      cand=IHhits[i];
    }
  }
  if( 10.0<min_dphi ) return false;
  track-> SetHodoHit(*cand, vtxIH);
  return true;
}

double MyTools::dphiCDH(CDSTrack *track, HodoscopeLikeHit *hit)
{
  TVector3 vtxCDH=track->CDHVertex();
  if( vtxCDH.X()==0 && vtxCDH.Y()==0 && vtxCDH.Z()==0 ){
    double param[5];
    track->GetParameters(param);
    vtxCDH=MathTools::CalcHelixPosatR(param, 54.4);
    track->SetCDHVertex(vtxCDH);
  }
  if( !hit || hit->cid()!=CID_CDH ) return DBL_MAX;

  double dphi = vtxCDH.Phi()-hit->pos().Phi();
  if( dphi>TMath::Pi() )       dphi -= 2*TMath::Pi();
  else if( dphi<-TMath::Pi() ) dphi += 2*TMath::Pi();

  return dphi*TMath::RadToDeg();
}

double MyTools::dphiIH(CDSTrack *track, HodoscopeLikeHit *hit)
{
  TVector3 vtxIH=track->IHVertex();
  if( vtxIH.X()==0 && vtxIH.Y()==0 && vtxIH.Z()==0 ){
    double param[5];
    track->GetParameters(param);
    vtxIH=MathTools::CalcHelixPosatR(param, 13.7);
    track->SetCDHVertex(vtxIH);
  }
  if( !hit || hit->cid()!=CID_IH ) return DBL_MAX;

  double dphi = vtxIH.Phi()-hit->pos().Phi();
  if( dphi>TMath::Pi() )       dphi -= 2*TMath::Pi();
  else if( dphi<-TMath::Pi() ) dphi += 2*TMath::Pi();

  return dphi*TMath::RadToDeg();
}

bool MyTools::sameCDHCluster(HodoscopeLikeHit *hit1, HodoscopeLikeHit *hit2)
{
  if( abs(hit1->seg()-hit2->seg())>1 ) return false; 
  if( fabs(hit1->ctmean()-hit2->ctmean())>0.2 ) return false;
  if( fabs(hit1->hitpos()-hit2->hitpos())>6.0 ) return false;
  return true;
}

bool MyTools::reanaFC(ConfMan *conf, BeamLineHitMan *blMan, AnaInfo *anaInfo, bool refit)
{
  //  cout<<"===== MyTools::MakeFC START ====="<<endl;
  ForwardChargeInfo fcInfo;
  vector<HodoscopeLikeHit*> PC_hits=getHodo(blMan, CID_PC);
  vector<HodoscopeLikeHit*> CVC_hits=getHodo(blMan, CID_CVC);
  vector<HodoscopeLikeHit*> BVC_hits=getHodo(blMan, CID_BVC);

  // cout<<"nCVC : "<<CVC_hits.size()<<endl;
  // cout<<"nPC  : "<<PC_hits.size()<<endl;
  // cout<<"nBVC : "<<BVC_hits.size()<<endl;
  // cout<<"nFDC1 : "<<anaInfo->beam(0)->nFDC1()<<endl;

  HodoscopeLikeHit *fc_hit=0;
  if( CVC_hits.size()==0 && PC_hits.size()==1 ) fc_hit=PC_hits[0];
  if( PC_hits.size()==0 && CVC_hits.size()==1 ) fc_hit=CVC_hits[0];

  if( !fc_hit ) return false;
  if( anaInfo->beam(0)->nFDC1()!=1 ) return false;
  if( !anaInfo->minDCA() ) return false;

  if( anaInfo->nFCharge()!=1 ){
    cout<<"!!!!! nFCharge="<<anaInfo->nFCharge()<<" !!!!!"<<endl;
    return false;
  }
  if( !anaInfo->minDCA() ) return false;

  anaInfo->forwardCharge(0)->SetHodo(fc_hit);
  anaInfo->forwardCharge(0)->fit_forward(blMan, anaInfo->beam(0), anaInfo->minDCA(), conf, refit);

  return true;
}

ForwardChargeInfo MyTools::makeFC(ConfMan *conf, BeamLineHitMan *blMan, AnaInfo *anaInfo)
{
  //  cout<<"===== MyTools::MakeFC START ====="<<endl;
  ForwardChargeInfo fcInfo;
  vector<HodoscopeLikeHit*> PC_hits=getHodo(blMan, CID_PC);
  vector<HodoscopeLikeHit*> CVC_hits=getHodo(blMan, CID_CVC);
  vector<HodoscopeLikeHit*> BVC_hits=getHodo(blMan, CID_BVC);

  // cout<<"nCVC : "<<CVC_hits.size()<<endl;
  // cout<<"nPC  : "<<PC_hits.size()<<endl;
  // cout<<"nBVC : "<<BVC_hits.size()<<endl;
  // cout<<"nFDC1 : "<<anaInfo->beam(0)->nFDC1()<<endl;
  // cout<<"minDCA : "<<anaInfo->minDCA()<<endl;

  HodoscopeLikeHit *fc_hit=0;
  if( CVC_hits.size()==0 && PC_hits.size()==1 ) fc_hit=PC_hits[0];
  if( PC_hits.size()==0 && CVC_hits.size()==1 ) fc_hit=CVC_hits[0];

  if( !fc_hit ) return fcInfo;
  if( anaInfo->beam(0)->nFDC1()!=1 ) return fcInfo;
  if( !anaInfo->minDCA() ) return fcInfo;

  fcInfo.SetVertex(anaInfo->minDCA()->vertexBeam());
  fcInfo.SetHodo(fc_hit);
  fcInfo.SetFDC1(anaInfo->beam(0)->FDC1(0));

  fcInfo.fit_forward(blMan, anaInfo->beam(0), anaInfo->minDCA(), conf);
  //  fcInfo.dump();
  //  cout<<"===== MyTools::MakeFC FINISH ====="<<endl;
  return fcInfo;
}

CDS2Info MyTools::makeCDS2Info(CDSTrackingMan *cdstrackMan, int id1, int id2, AnaInfo *anaInfo)
{
  CDS2Info info;
  if( anaInfo->nBeam()!=1 && !anaInfo->beam(0)->flag() ) return info;
  if( !anaInfo->CDSbyID(id1) && !anaInfo->CDSbyID(id2) ) return info;
  BeamInfo *beam=anaInfo->beam(0);

  CDSInfo *cds1 = anaInfo->CDSbyID(id1);
  CDSInfo *cds2 = anaInfo->CDSbyID(id2);
  CDSTrack *track1=cdstrackMan->Track(id1);
  CDSTrack *track2=cdstrackMan->Track(id2);
  TVector3 vtx1=DEFVECT, vtx2=DEFVECT;
  bool vtx_flag=TrackTools::Calc2HelixVertex(track1, track2, vtx1, vtx2);
  TVector3 p1=DEFVECT, p2=DEFVECT;
  bool mom_flag1=false, mom_flag2=false;
  if( vtx_flag ){
    mom_flag1=track1-> GetMomentum(vtx1, p1, true, true);
    mom_flag2=track2-> GetMomentum(vtx2, p2, true, true);
  }

  TVector3 vtx=0.5*(vtx1+vtx2);
  TVector3 sum_p=p1+p2;
  double dist, dltmp=0;
  TVector3 vtxCDS, vtxBeam;
  MathTools::LineToLine(vtx, sum_p, beam->BPCpos(), beam->BPCdir(), dltmp, dist, vtxCDS, vtxBeam);
  info.SetTrackID1(id1);
  info.SetTrackID2(id2);
  info.SetPID1(cds1->pid());
  info.SetPID2(cds2->pid());
  info.SetVertex1(vtx1);
  info.SetVertex2(vtx2);
  info.SetVertexBeam(vtxBeam);
  info.SetMomentum1(p1);
  info.SetMomentum2(p2);
  if( vtx_flag && mom_flag1 && mom_flag2 && cds1->flag() && cds2->flag() ) info.SetFlag(true);

  return info;
}

ForwardNeutralInfo MyTools::makeFN(BeamLineHitMan *blMan, AnaInfo *anaInfo, const double &thre, bool sim)
{
  ForwardNeutralInfo info;
  if( anaInfo->nBeam()!=1 ){ return info; }
  if( !anaInfo->beam(0)->flag() ){ return info; }

  CDSInfo *cdsInfo =  anaInfo->minDCA();
  if( !cdsInfo || !cdsInfo->flag() ){ return info; }

  vector<HodoscopeLikeHit*> CVC_hits=getHodo(blMan, CID_CVC);
  vector<HodoscopeLikeHit*> PC_hits=getHodo(blMan, CID_PC);
  vector<vector<HodoscopeLikeHit*> > NC_hits=getNChits(blMan);
  vector<HodoscopeLikeHit*> NC_hits2=getHodo(blMan, CID_NC);

  HodoscopeLikeHit *nc_hit=0;
  for( int lay=0; lay<(int)NC_hits.size(); lay++ ){
    //    cout<<"NC layer"<<lay<<" search "<<endl;
    double time=DBL_MAX;
    for( int i=0; i<(int)NC_hits[lay].size(); i++ ){
      if( NC_hits[lay][i]->emean()>thre && NC_hits[lay][i]->ctmean()<time ){
        time=NC_hits[lay][i]->ctmean();
        nc_hit=NC_hits[lay][i];
      }
    }
    if( nc_hit ) break;
  }

  if( !nc_hit ) return info;
  //    cout<<"find NC hit"<<endl;
  vector<int> clus_seg;
  bool isAdd=true;
  clus_seg.push_back(nc_hit->seg());
  while(isAdd){
    isAdd=false;
    for( int i=0; i<(int)NC_hits2.size(); i++ ){
      int seg=NC_hits2[i]->seg();
      int lay=(seg-1)/16;
      int seg2=seg-16*lay;
      for( int j=0; j<(int)clus_seg.size(); j++ ){
	int cl_lay=(clus_seg[j]-1)/16;
	int cl_seg2=clus_seg[j]-16*cl_lay;
	if( seg==clus_seg[j] ){
	  isAdd=false;
	  break;
	}
	if( abs(lay-cl_lay)==1 && cl_seg2==seg2 ) isAdd=true;
	if( lay==cl_lay && abs(cl_seg2-seg2)==1 ) isAdd=true;
      }
      if( isAdd ){
	clus_seg.push_back(seg);
	break;
      }
    }
  }
  //    cout<<"make fnInfo"<<endl;
  info.SetVertex(cdsInfo->vertexBeam());
  info.SetHodo(nc_hit);
  for( int i=0; i<(int)clus_seg.size(); i++ ) info.SetCluster(clus_seg[i]);
  info.calc(anaInfo->beam(0), sim);
  //    cout<<"return true"<<endl;
  return info;
}

BeamInfo MyTools::makeBeamInfo(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamSpectrometer *spec)
{
  BeamInfo beam;
  vector<HodoscopeLikeHit*> T0hits=getHodo(blMan, CID_T0);
  if( T0hits.size()!=1 ) return beam;
  HodoscopeLikeHit *T0hit=T0hits[0];

  int trig_par = header->trigparticle();
  double expected_tof;
  if( trig_par==Beam_Kaon || trig_par==Beam_Other ) expected_tof=28.6429;
  if( trig_par==Beam_Pion   ) expected_tof=25.9334;
  if( trig_par==Beam_Proton ) expected_tof=35.22;

  double tof_diff=DBL_MAX;
  HodoscopeLikeHit *BHDhit=0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    HodoscopeLikeHit *hit=blMan->BHD(i);
    if( hit-> CheckRange() ){
      double tof=T0hit->ctmean()-hit->ctmean();
      if( fabs(tof-expected_tof)<tof_diff ){
        BHDhit=hit;
        tof_diff=fabs(tof-expected_tof);
      }
    }
  }

  int BLC1id=-1, BLC2id=-1, BPCid=-1;
  int ntrackBLC1=0, ntrackBLC2=0, ntrackBPC=0;  
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    LocalTrack *track = bltrackMan->trackBLC1(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      BLC1id=i;
      ntrackBLC1++;
    }
  }
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    LocalTrack *track = bltrackMan->trackBLC2(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      BLC2id=i;
      ntrackBLC2++;
    }
  }
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    LocalTrack *track = bltrackMan->trackBPC(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      BPCid=i;
      ntrackBPC++;
    }
  }

  beam.SetPID(trig_par);
  beam.SetBLMan(blMan);
  beam.SetBLDC(bltrackMan);
  if( T0hit ) beam.SetT0(T0hit);
  if( BHDhit ) beam.SetBHD(BHDhit);
  if( ntrackBLC1==1 ) beam.SetBLC1id(BLC1id);
  if( ntrackBLC2==1 ) beam.SetBLC2id(BLC2id);
  if( ntrackBPC==1  ) beam.SetBPCid(BPCid, conf);
  if( spec->stat()>0 ) beam.SetD5(spec);
  if( T0hit && BHDhit && ntrackBLC1==1 && ntrackBLC2==1 && ntrackBPC==1 && spec->stat()>0 ){
    beam.SetFlag(true);
  }
  return beam;
}

int MyTools::CDSanaLevel(CDSTrack *track, CDSHitMan *cdsMan, BeamInfo *beam, BeamLineTrackMan *bltrackMan, ConfMan *conf, bool sim)
{
  double mom=track->Momentum();
  TVector3 vtxCDH=track->CDHVertex();
  if( vtxCDH.X()==0 && vtxCDH.Y()==0 && vtxCDH.Z()==0 ){
    double param[5];
    track->GetParameters(param);
    vtxCDH=MathTools::CalcHelixPosatR(param, 54.4);
    track->SetCDHVertex(vtxCDH);
  }
  if( !track->CDHFlag() ) searchCDHHit(track, cdsMan);
  if( !track->IHFlag() ) searchIHHit(track, cdsMan);
  HodoscopeLikeHit *CDH_hit=0;
  if( track->CDHFlag() ) CDH_hit=track->CDHHit(cdsMan, 0);
  HodoscopeLikeHit *IH_hit=0;
  if( track->IHFlag() ) IH_hit=track->IHHit(cdsMan, 0);

  TVector3 vtxBeam, vtxCDS, p=DEFVECT;
  bool vtx_flag=track->GetVertex(beam->T0pos(), beam->BPCdir(), vtxBeam, vtxCDS);

  bool find_mass2_flag=false, calc_tv_flag=false, mom_flag=false;
  double mass2=DEFAULTD, calc_beta=DEFAULTD, tofvtxcdc=DEFAULTD, offset=DEFAULTD;
  int pid=CDS_Other;

  if( !CDH_hit || !vtx_flag ) return 0;
  if( CDH_hit && vtx_flag ){
    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(beam->T0pos(), vtxBeam, beam->D5mom(), beam->mass(), beam_out, beam_tof);

    double par[5];
    track-> GetParameters(par);
    double tof=CDH_hit->ctmean()-beam->T0time();
    double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, vtxCDS);
    double beta = cdc_dis/(tof-beam_tof)/(100.*Const);
    mass2 = mom*mom*(1./(beta*beta)-1);
    cout<<beam->BPC(bltrackMan)<<"  "<<track<<endl;
    if( !beam->BPC(bltrackMan) ) return 1;
    if( !track ) return 1;
    if( !TrackTools::FindMass2(track, beam->BPC(bltrackMan), tof, beam->D5mom(), beam->pid(), calc_beta, mass2, tofvtxcdc) ) return 1;

    if( find_mass2_flag ){
      if( conf && !sim ){
        track->Retiming(cdsMan, conf, calc_beta, true);
        track->Calc(conf);
        if( !TrackTools::FindMass2(track, beam->BPC(bltrackMan), tof, beam->D5mom(), beam->pid(), calc_beta, mass2, tofvtxcdc) ) return 1;
      }
      pid=PIDcorr_wide(mom, mass2);
      // if( conf->GetCDSFittingParamManager()-> PID(CDS_PiMinus, mom, mass2) ) pid=CDS_PiMinus;
      // else if( conf->GetCDSFittingParamManager()-> PID(CDS_PiPlus, mom, mass2) ) pid=CDS_PiPlus;
      // else if( conf->GetCDSFittingParamManager()-> PID(CDS_Proton, mom, mass2) ) pid=CDS_Proton;
      // else if( conf->GetCDSFittingParamManager()-> PID(CDS_Kaon, mom, mass2)   ) pid=CDS_Kaon;
      // else if( conf->GetCDSFittingParamManager()-> PID(CDS_Deuteron, mom, mass2) ) pid=CDS_Deuteron;

      double tmpfl;
      if( !track-> CalcVertexTimeLength(beam->T0pos(), beam->BPCdir(), cdsMass[pid], vtxBeam, vtxCDS, tof, tmpfl, false) ) return 2;
      if( !track->GetMomentum(vtxCDS, p, true, true) ) return 3;
      if( pid==CDS_PiMinus || pid==CDS_PiPlus || pid==CDS_Kaon || pid==CDS_Proton || pid==CDS_Deuteron ){
	ELossTools::CalcElossBeamTGeo(beam->T0pos(), vtxBeam, beam->D5mom(), beam->mass(), beam_out, beam_tof);
	TVector3 vtxCFRP;
	double param[5];
	track-> GetParameters(CID_CDCCFRP, param, vtxCFRP);
	double fl = MathTools::CalcHelixArc(param, vtxCFRP, vtxCDH);
	double mom2 = fabs(track-> Momentum(CID_CDCCFRP));
	double mass = cdsMass[pid];
	double momout;
	double tmptof;
	ELossTools::CalcdE(mom2, mass, fl, "CDCGas", momout, -1, tmptof);
	offset = beam_tof+tof+tmptof;
      }
    }
  }
  return 10;
}


CDSInfo MyTools::makeCDSInfo(int id, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, BeamInfo *beam, BeamLineTrackMan *bltrackMan, ConfMan *conf, bool sim)
{
  CDSTrack *track = cdstrackMan->Track(id);
  double mom=track->Momentum();
  TVector3 vtxCDH=track->CDHVertex();
  if( vtxCDH.X()==0 && vtxCDH.Y()==0 && vtxCDH.Z()==0 ){
    double param[5];
    track->GetParameters(param);
    vtxCDH=MathTools::CalcHelixPosatR(param, 54.4);
    track->SetCDHVertex(vtxCDH);
  }
  if( !track->CDHFlag() ) searchCDHHit(track, cdsMan);
  if( !track->IHFlag() ) searchIHHit(track, cdsMan);
  HodoscopeLikeHit *CDH_hit=0;
  if( track->CDHFlag() ) CDH_hit=track->CDHHit(cdsMan, 0);
  HodoscopeLikeHit *IH_hit=0;
  if( track->IHFlag() ) IH_hit=track->IHHit(cdsMan, 0);

  TVector3 vtxBeam, vtxCDS, p=DEFVECT;
  bool vtx_flag=track->GetVertex(beam->T0pos(), beam->BPCdir(), vtxBeam, vtxCDS);

  bool find_mass2_flag=false, calc_tv_flag=false, mom_flag=false;
  double mass2=DEFAULTD, calc_beta=DEFAULTD, tofvtxcdc=DEFAULTD, offset=DEFAULTD;
  int pid=CDS_Other;

  if( CDH_hit && vtx_flag ){
    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(beam->T0pos(), vtxBeam, beam->D5mom(), beam->mass(), beam_out, beam_tof);

    double par[5];
    track-> GetParameters(par);
    double tof=CDH_hit->ctmean()-beam->T0time();
    double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, vtxCDS);
    double beta = cdc_dis/(tof-beam_tof)/(100.*Const);
    mass2 = mom*mom*(1./(beta*beta)-1);
    find_mass2_flag=TrackTools::FindMass2(track, beam->BPC(bltrackMan), tof, beam->D5mom(), beam->pid(), calc_beta, mass2, tofvtxcdc);

    if( find_mass2_flag ){
      if( conf && !sim ){
        track->Retiming(cdsMan, conf, calc_beta, true);
        track->Calc(conf);
        find_mass2_flag=TrackTools::FindMass2(track, beam->BPC(bltrackMan), tof, beam->D5mom(), beam->pid(), calc_beta, mass2, tofvtxcdc);
      }
      pid=PIDcorr_wide(mom, mass2);
      // if( conf->GetCDSFittingParamManager()-> PID(CDS_PiMinus, mom, mass2) ) pid=CDS_PiMinus;
      // else if( conf->GetCDSFittingParamManager()-> PID(CDS_PiPlus, mom, mass2) ) pid=CDS_PiPlus;
      // else if( conf->GetCDSFittingParamManager()-> PID(CDS_Proton, mom, mass2) ) pid=CDS_Proton;
      // else if( conf->GetCDSFittingParamManager()-> PID(CDS_Kaon, mom, mass2)   ) pid=CDS_Kaon;
      // else if( conf->GetCDSFittingParamManager()-> PID(CDS_Deuteron, mom, mass2) ) pid=CDS_Deuteron;

      double tmpfl;
      calc_tv_flag=track-> CalcVertexTimeLength(beam->T0pos(), beam->BPCdir(), cdsMass[pid], vtxBeam, vtxCDS, tof, tmpfl, true);
      if( calc_tv_flag ){
	mom_flag=track->GetMomentum(vtxCDS, p, true, true);
	if( pid==CDS_PiMinus || pid==CDS_PiPlus || pid==CDS_Kaon || pid==CDS_Proton || pid==CDS_Deuteron ){
	  ELossTools::CalcElossBeamTGeo(beam->T0pos(), vtxBeam, beam->D5mom(), beam->mass(), beam_out, beam_tof);
	  TVector3 vtxCFRP;
	  double param[5];
	  track-> GetParameters(CID_CDCCFRP, param, vtxCFRP);
	  double fl = MathTools::CalcHelixArc(param, vtxCFRP, vtxCDH);
	  double mom2 = fabs(track-> Momentum(CID_CDCCFRP));
	  double mass = cdsMass[pid];
	  double momout;
	  double tmptof;
	  ELossTools::CalcdE(mom2, mass, fl, "CDCGas", momout, -1, tmptof);
	  offset = beam_tof+tof+tmptof;
	}
      }
    }
  }

  CDSInfo cds;
  cds.SetTrackID(id);
  cds.SetPID(pid);
  if( CDH_hit ) cds.SetCDHseg(CDH_hit->seg());
  if( IH_hit  ) cds.SetIHseg(IH_hit->seg());
  cds.SetMom(mom);
  cds.SetMass2(mass2);
  cds.SetBeta(calc_beta);
  cds.SetVertexCDS(vtxCDS);
  cds.SetVertexBeam(vtxBeam);
  cds.SetMomentum(p);
  cds.SetOffset(offset);
  if( find_mass2_flag && calc_tv_flag && mom_flag && track->Chi()<30 ){
    cds.SetFlag(true);
  }

  return cds;
}

CDSInfo MyTools::makeCDSInfo_z40
(int id, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, BeamInfo *beam, BeamLineTrackMan *bltrackMan, ConfMan *conf, bool sim)
{
  CDSTrack *track = cdstrackMan->Track(id);
  double mom=track->Momentum();
  TVector3 vtxCDH=track->CDHVertex();
  if( vtxCDH.X()==0 && vtxCDH.Y()==0 && vtxCDH.Z()==0 ){
    double param[5];
    track->GetParameters(param);
    vtxCDH=MathTools::CalcHelixPosatR(param, 54.4);
    track->SetCDHVertex(vtxCDH);
  }
  if( !track->CDHFlag() ) searchCDHHit(track, cdsMan);
  if( !track->IHFlag() ) searchIHHit(track, cdsMan);
  HodoscopeLikeHit *CDH_hit=0;
  if( track->CDHFlag() ) CDH_hit=track->CDHHit(cdsMan, 0);
  HodoscopeLikeHit *IH_hit=0;
  if( track->IHFlag() ) IH_hit=track->IHHit(cdsMan, 0);

  TVector3 vtxBeam, vtxCDS, p=DEFVECT;
  bool vtx_flag=track->GetVertex(beam->T0pos(), beam->BPCdir(), vtxBeam, vtxCDS);

  bool find_mass2_flag=false, calc_tv_flag=false, mom_flag=false;
  double mass2=DEFAULTD, calc_beta=DEFAULTD, tofvtxcdc=DEFAULTD, offset=DEFAULTD;
  int pid=CDS_Other;

  if( CDH_hit && vtx_flag ){
    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(beam->T0pos(), vtxBeam, beam->D5mom(), beam->mass(), beam_out, beam_tof);

    double par[5];
    track-> GetParameters(par);
    double tof=CDH_hit->ctmean()-beam->T0time();
    double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, vtxCDS);
    double beta = cdc_dis/(tof-beam_tof)/(100.*Const);
    mass2 = mom*mom*(1./(beta*beta)-1);
    find_mass2_flag=TrackTools::FindMass2(track, beam->BPC(bltrackMan), tof, beam->D5mom(), beam->pid(), calc_beta, mass2, tofvtxcdc);

    if( find_mass2_flag ){
      if( conf && !sim ){
        track->Retiming(cdsMan, conf, calc_beta, true);
        track->Calc(conf);
        find_mass2_flag=TrackTools::FindMass2(track, beam->BPC(bltrackMan), tof, beam->D5mom(), beam->pid(), calc_beta, mass2, tofvtxcdc);
      }
      pid=PIDcorr_wide(mom, mass2);

      double tmpfl;
      calc_tv_flag=track-> CalcVertexTimeLength(beam->T0pos(), beam->BPCdir(), cdsMass[pid], vtxBeam, vtxCDS, tof, tmpfl, true);
      if( calc_tv_flag ){
	mom_flag=track->GetMomentum(vtxCDS, p, true, true);
	if( pid==CDS_PiMinus || pid==CDS_PiPlus || pid==CDS_Kaon || pid==CDS_Proton || pid==CDS_Deuteron ){
	  ELossTools::CalcElossBeamTGeo(beam->T0pos(), vtxBeam, beam->D5mom(), beam->mass(), beam_out, beam_tof);
	  TVector3 vtxCFRP;
	  double param[5];
	  track-> GetParameters(CID_CDCCFRP, param, vtxCFRP);
	  double fl = MathTools::CalcHelixArc(param, vtxCFRP, vtxCDH);
	  double mom2 = fabs(track-> Momentum(CID_CDCCFRP));
	  double mass = cdsMass[pid];
	  double momout;
	  double tmptof;
	  ELossTools::CalcdE(mom2, mass, fl, "CDCGas", momout, -1, tmptof);
	  offset = beam_tof+tof+tmptof;
	}
      }
    }
  }

  CDSInfo cds;
  cds.SetTrackID(id);
  cds.SetPID(pid);
  if( CDH_hit ) cds.SetCDHseg(CDH_hit->seg());
  if( IH_hit  ) cds.SetIHseg(IH_hit->seg());
  cds.SetMom(mom);
  cds.SetMass2(mass2);
  cds.SetBeta(calc_beta);
  cds.SetVertexCDS(vtxCDS);
  cds.SetVertexBeam(vtxBeam);
  cds.SetMomentum(p);
  cds.SetOffset(offset);
  if( find_mass2_flag && calc_tv_flag && mom_flag && track->Chi()<30 && fabs(vtxCDH.Z())<40 ){
    cds.SetFlag(true);
  }

  return cds;
}

HodoscopeLikeHit* MyTools::getCDH(CDSTrack *track, CDSHitMan *cdsMan)
{
  vector<HodoscopeLikeHit*> CDH_hit=getCDH(cdsMan);

  double min_dphiCDH=DBL_MAX;
  HodoscopeLikeHit *nearest_CDH=0;
  for( int j=0; j<(int)CDH_hit.size(); j++ ){
    double dphi = dphiCDH(track, CDH_hit[j]);
    if( fabs(dphi)<min_dphiCDH ){
      min_dphiCDH= fabs(dphi);
      nearest_CDH = CDH_hit[j];
    }
  }

  return nearest_CDH;
}

HodoscopeLikeHit* MyTools::getCDH2(CDSTrack *track, CDSHitMan *cdsMan)
{
  vector<HodoscopeLikeHit*> CDH_hit=getCDH(cdsMan);

  double min_dphiCDH=DBL_MAX;
  HodoscopeLikeHit *nearest_CDH=0;
  for( int j=0; j<(int)CDH_hit.size(); j++ ){
    double dphi = dphiCDH(track, CDH_hit[j]);

    if( fabs(dphi)<min_dphiCDH ){
      min_dphiCDH= fabs(dphi);
      nearest_CDH = CDH_hit[j];
    }
  }
  if( !nearest_CDH ) return 0;
  if( min_dphiCDH>7.0 ) return 0;

  if( fabs(nearest_CDH->hitpos()-track->CDHVertex().Z())<6.0 ){
    //    cout<<" return nearest hit  dphi="<<min_dphiCDH<<endl;
    return nearest_CDH;
  }

  HodoscopeLikeHit *cand_CDH1=0;
  HodoscopeLikeHit *cand_CDH2=0;
  for( int j=0; j<(int)CDH_hit.size(); j++ ){
    if( nearest_CDH->seg()==CDH_hit[j]->seg()-1 ) cand_CDH1=CDH_hit[j];
    if( nearest_CDH->seg()==CDH_hit[j]->seg()+1 ) cand_CDH2=CDH_hit[j];
  }

  double z_diff0=fabs(track->CDHVertex().Z()-nearest_CDH->hitpos());
  double z_diff1=DBL_MAX; double z_diff2=DBL_MAX;
  if( cand_CDH1 ) z_diff1=fabs(track->CDHVertex().Z()-cand_CDH1->hitpos());
  if( cand_CDH2 ) z_diff2=fabs(track->CDHVertex().Z()-cand_CDH2->hitpos());
  if( z_diff1<z_diff0 && z_diff1<z_diff2 ){
    cout<<" return side hit z_dfii="<<z_diff1<<"  z_diff0="<<z_diff0<<endl;
    return cand_CDH1;
  }
  if( z_diff2<z_diff0 && z_diff2<z_diff1 ){
    cout<<" return side hit z_dfii="<<z_diff2<<"  z_diff0="<<z_diff0<<endl;
    return cand_CDH2;
  }
  //  cout<<" return nearest hit    dphi="<<min_dphiCDH<<endl;
  return nearest_CDH;
}

HodoscopeLikeHit* MyTools::getIH(CDSTrack *track, CDSHitMan *cdsMan)
{
  vector<HodoscopeLikeHit*> IH_hit=getIH(cdsMan);

  TVector3 vtxIH=track->IHVertex();
  double min_dphiIH=DBL_MAX;
  HodoscopeLikeHit *nearest_IH=0;
  for( int j=0; j<(int)IH_hit.size(); j++ ){
    double dphi = dphiIH(track, IH_hit[j]);

    if( fabs(dphi)<min_dphiIH ){
      min_dphiIH= fabs(dphi);
      nearest_IH = IH_hit[j];
    }
  }

  return nearest_IH;
}

bool MyTools::isBeamFiducial(BeamLineTrackMan *bltrackMan)
{
  LocalTrack *BPC=trackBPC(bltrackMan);
  if( !BPC ) return false;

  TVector3 beamFF=BPC->GetPosatZ(0);
  if( GeomTools::GetID(beamFF)==CID_Fiducial ) return true;
  return false;
}

bool MyTools::isBeamMom(BeamSpectrometer *beam)
{
  if( beam->stat()>0 ){
    if( beam->chisquare()<30 ) return true;
  }

  return false;
}

LocalTrack *MyTools::trackBLC1(BeamLineTrackMan *bltrackMan)
{
  int ntrack=0;
  LocalTrack *track=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    double time=bltrackMan->trackBLC1(i)->GetTrackTime();
    if( -30<time && time<100 ){
      ntrack++;
      track=bltrackMan->trackBLC1(i);
    }
  }
  if( ntrack==1 ) return track;

  return 0;
}

LocalTrack *MyTools::trackBLC2(BeamLineTrackMan *bltrackMan)
{
  int ntrack=0;
  LocalTrack *track=0;
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    double time=bltrackMan->trackBLC2(i)->GetTrackTime();
    if( -30<time && time<100 ){
      ntrack++;
      track=bltrackMan->trackBLC2(i);
    }
  }
  if( ntrack==1 ) return track;

  return 0;
}

LocalTrack *MyTools::trackBPC(BeamLineTrackMan *bltrackMan)
{
  int ntrack=0;
  LocalTrack *track=0;
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    double time=bltrackMan->trackBPC(i)->GetTrackTime();
    if( -30<time && time<100 ){
      ntrack++;
      track=bltrackMan->trackBPC(i);
    }
  }
  if( ntrack==1 ) return track;

  return 0;
}

bool MyTools::CDCall1hit(CDSHitMan* cdsMan)
{
  bool flag=true;
  for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
    if( cdsMan->nCDC(lay)!=1 ) flag=false;
  }
  return flag;
}

vector<vector<HodoscopeLikeHit*> > MyTools::getNChits(BeamLineHitMan *blMan)
{
  vector<vector<HodoscopeLikeHit*> > cluster(7);
  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->CheckRange() ){
      int seg=blMan->NC(i)->seg();
      int lay=(seg-1)/16;
      if( 7<lay ) cout<<"  !!! NC layer>7    "<<lay<<" !!!"<<endl;
      cluster[lay].push_back(blMan->NC(i));
    }
  }
  return cluster;
}

vector<HodoscopeLikeHit*> MyTools::getCDHCluster(HodoscopeLikeHit *seed, CDSHitMan *cdsMan)
{
  vector<HodoscopeLikeHit*> cluster;
  cluster.push_back(seed);
  while( true ){
    HodoscopeLikeHit *next=0;
    for( int i=0; i<cdsMan->nCDH(); i++ ){
      if( cdsMan->CDH(i)->CheckRange() ){
	bool add=false; bool same=false;
	for( int j=0; j<(int)cluster.size(); j++ ){
	  if( cdsMan->CDH(i)->seg()==cluster[j]->seg() ) same=true;
	  if( abs(cdsMan->CDH(i)->seg()-cluster[j]->seg())==1 ) add=true;
	}
	if( same ) continue;
	if( add ){
	  next=cdsMan->CDH(i);
	  break;
	}
      }
    }

    if( next ) cluster.push_back(next);
    if( !next ) break;
  }
  return cluster;
}

vector<HodoscopeLikeHit*> MyTools::getCDH(CDSHitMan *cdsMan)
{
  vector<HodoscopeLikeHit*> hits;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    if( cdsMan->CDH(i)->CheckRange() ){
      hits.push_back(cdsMan->CDH(i));
    }
  }
  return hits;
}

vector<HodoscopeLikeHit*> MyTools::getIH(CDSHitMan *cdsMan)
{
  vector<HodoscopeLikeHit*> hits;
  for( int i=0; i<cdsMan->nIH(); i++ ){
    if( 0<cdsMan->IH(i)->tdcu() && cdsMan->IH(i)->tdcu()<4000 ){
      hits.push_back(cdsMan->IH(i));
    }
  }
  return hits;
}

vector<HodoscopeLikeHit*> MyTools::getHodo(BeamLineHitMan *blMan, const int &cid)
{
  vector<HodoscopeLikeHit*> hits;
  for( int i=0; i<blMan->nHodo(cid); i++ ){
    HodoscopeLikeHit *hit=blMan->Hodoi(cid, i);
    //    cout<<"hit="<<hit<<endl;
    if( hit-> CheckRange() ){
      hits.push_back(hit);
    }
  }
  return hits;
}

template <class T>
T* get(const TString &name, TFile *f)
{
  T* ptr = dynamic_cast<T*>(f->Get(name));
  if( !ptr ){
    cout<<"  !!! "<<name<<" not found !!!"<<endl;
  }
  return ptr;
}

bool MyTools::isFiducial(AnaInfo *info)
{
  CDSInfo *min_cds=info->minDCA();
  if( !min_cds ) return false;

  return isFiducial(min_cds);
}

bool MyTools::isFiducial(CDSInfo *info)
{
  if( GeomTools::GetID(info->vertexBeam())==CID_Fiducial ) return true;
  return false;
}

int MyTools::checkArea(double mom, double mass2)
{
  double pi_sigma = sqrt(4.*piMass*piMass*mom*mom*PID_Param[0]+4.*piMass*piMass*piMass*piMass*PID_Param[1]*(1+piMass*piMass/(mom*mom))+4.*mom*mom*(piMass*piMass+mom*mom)*PID_Param[2]);
  double k_sigma = sqrt(4.*kpMass*kpMass*mom*mom*PID_Param[0]+4.*kpMass*kpMass*kpMass*kpMass*PID_Param[1]*(1+kpMass*kpMass/(mom*mom))+4.*mom*mom*(kpMass*kpMass+mom*mom)*PID_Param[2]);
  double p_sigma = sqrt(4.*pMass*pMass*mom*mom*PID_Param[0]+4.*pMass*pMass*pMass*pMass*PID_Param[1]*(1+pMass*pMass/(mom*mom))+4.*mom*mom*(pMass*pMass+mom*mom)*PID_Param[2]);

  if( mom<0 ){
    if( piMass*piMass+2.5*pi_sigma<mass2 && mass2<kpMass*kpMass-2.5*k_sigma ) return -1;
  }
  else{
    if( piMass*piMass+2.5*pi_sigma<mass2 && mass2<pMass*pMass-2.5*p_sigma ) return 1;
  }
  return 0;
}

bool MyTools::isShareCDH(CDSHitMan *cdsMan, CDSTrackingMan *trackMan, const int &id)
{
  CDSTrack *track = trackMan->GoodTrack(id);
  std::vector<int> segments;
  for( int i=0; i<track->nCDHHit(); i++ ){ segments.push_back(track->CDHHit(cdsMan, i)->seg()); }

  TString tstr=" seg : ";
  for( int i=0; i<(int)segments.size(); i++ ){ tstr += Form("%d", segments[i]); }

  bool result = false;
  for( int i=0; i<trackMan->nGoodTrack(); i++ ){
    if( id==i ) continue;
    CDSTrack *track2 = trackMan->GoodTrack(i);
    for( int j=0; j<track2->nCDHHit(); j++ ){
      for( int k=0; k<(int)segments.size(); k++ ){
	if( track2->CDHHit(cdsMan, j)->seg()==segments[k] ) result=true;
      }
    }
  }

  // if( result ){
  //   cout<<"MyTool::isShareCDH id="<<id<<endl;
  //   cout<<tstr<<endl;
  //   cout<<" nGoodTrack : "<<trackMan->nGoodTrack()<<endl;
  // }
  //  if( result ) cout<<" This track share CDH Hit"<<endl;

  return result;
}

int MyTools::PID(double mom,double mass2)
{
  int ptype=-1;
  if(mom>0) {
    if( mass2<-0.02 ) ptype=CDS_Other;
    else if( mass2<0.06) ptype=CDS_PiPlus;
    else if( mass2<0.3) ptype=CDS_Other;
    else if(mass2<2.0) ptype=CDS_Proton;
    else if(mass2<5.) ptype=CDS_Deuteron;
    else if(mass2<9.) ptype=CDS_Triton;
    else ptype=CDS_Other;
  }
  else{
    if( mass2<-0.02 ) ptype=CDS_Other;
    else if( mass2<0.06) ptype=CDS_PiMinus;
    else if( mass2<0.12) ptype=CDS_Other;
    else if(mass2<0.5) ptype=CDS_Kaon;
    else ptype=CDS_Other;
  }
  return ptype;
}

int MyTools::PIDcorr(double mom,double mass2)
{
  double pi_sigma = sqrt(4.*piMass*piMass*mom*mom*PID_Param[0]+4.*piMass*piMass*piMass*piMass*PID_Param[1]*(1+piMass*piMass/(mom*mom))+4.*mom*mom*(piMass*piMass+mom*mom)*PID_Param[2]);
  double k_sigma = sqrt(4.*kpMass*kpMass*mom*mom*PID_Param[0]+4.*kpMass*kpMass*kpMass*kpMass*PID_Param[1]*(1+kpMass*kpMass/(mom*mom))+4.*mom*mom*(kpMass*kpMass+mom*mom)*PID_Param[2]);
  double p_sigma = sqrt(4.*pMass*pMass*mom*mom*PID_Param[0]+4.*pMass*pMass*pMass*pMass*PID_Param[1]*(1+pMass*pMass/(mom*mom))+4.*mom*mom*(pMass*pMass+mom*mom)*PID_Param[2]);

  bool pim_flag=false, pip_flag=false, km_flag=false, p_flag=false;
  if( mom>0 ){
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pip_flag=true;
    if( pMass*pMass-2.5*p_sigma<mass2 && mass2<pMass*pMass+2.5*p_sigma ) p_flag=true;

    if( pip_flag ) return CDS_PiPlus;
    else if( p_flag ) return CDS_Proton;
  }
  else{
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pim_flag=true;
    if( kpMass*kpMass-2.5*k_sigma<mass2 && mass2<kpMass*kpMass+2.5*k_sigma ) km_flag=true;

    if( pim_flag ) return CDS_PiMinus;
    else if( km_flag ) return CDS_Kaon;
  }
  return CDS_Other;
}

int MyTools::PIDcorr2(double mom,double mass2)
{
  double pi_sigma = sqrt(4.*piMass*piMass*mom*mom*PID_Param[0]+4.*piMass*piMass*piMass*piMass*PID_Param[1]*(1+piMass*piMass/(mom*mom))+4.*mom*mom*(piMass*piMass+mom*mom)*PID_Param[2]);
  double k_sigma = sqrt(4.*kpMass*kpMass*mom*mom*PID_Param[0]+4.*kpMass*kpMass*kpMass*kpMass*PID_Param[1]*(1+kpMass*kpMass/(mom*mom))+4.*mom*mom*(kpMass*kpMass+mom*mom)*PID_Param[2]);
  double p_sigma = sqrt(4.*pMass*pMass*mom*mom*PID_Param[0]+4.*pMass*pMass*pMass*pMass*PID_Param[1]*(1+pMass*pMass/(mom*mom))+4.*mom*mom*(pMass*pMass+mom*mom)*PID_Param[2]);

  bool pim_flag=false, pip_flag=false, km_flag=false, p_flag=false;
  if( mom>0 ){
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pip_flag=true;
    if( pMass*pMass-2.5*p_sigma<mass2 && mass2<pMass*pMass+2.5*p_sigma && mom>0.1 ) p_flag=true;

    if( pip_flag && p_flag ){
      if( mass2<Ppi_mid_mass2 ) return CDS_PiPlus;
      else return CDS_Proton;
    }
    else if( pip_flag ) return CDS_PiPlus;
    else if( p_flag ) return CDS_Proton;
  }
  else{
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pim_flag=true;
    if( kpMass*kpMass-2.5*k_sigma<mass2 && mass2<kpMass*kpMass+2.5*k_sigma && mom<-0.03 ) km_flag=true;

    if( pim_flag && km_flag ){
      if( mass2<Kpi_mid_mass2 ) return CDS_PiMinus;
      else return CDS_Kaon;
    }

    if( pim_flag ) return CDS_PiMinus;
    else if( km_flag ) return CDS_Kaon;
  }
  return CDS_Other;
}

int MyTools::PIDcorr_wide(double mom,double mass2)
{
  double pi_sigma = sqrt(4.*piMass*piMass*mom*mom*PID_Param[0]+4.*piMass*piMass*piMass*piMass*PID_Param[1]*(1+piMass*piMass/(mom*mom))+4.*mom*mom*(piMass*piMass+mom*mom)*PID_Param[2]);
  double k_sigma = sqrt(4.*kpMass*kpMass*mom*mom*PID_Param[0]+4.*kpMass*kpMass*kpMass*kpMass*PID_Param[1]*(1+kpMass*kpMass/(mom*mom))+4.*mom*mom*(kpMass*kpMass+mom*mom)*PID_Param[2]);
  double p_sigma = sqrt(4.*pMass*pMass*mom*mom*PID_Param[0]+4.*pMass*pMass*pMass*pMass*PID_Param[1]*(1+pMass*pMass/(mom*mom))+4.*mom*mom*(pMass*pMass+mom*mom)*PID_Param[2]);

  bool pim_flag=false, pip_flag=false, km_flag=false, p_flag=false;
  if( mom>0 ){
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pip_flag=true;
    if( pMass*pMass-2.5*p_sigma<mass2 && mass2<pMass*pMass+2.5*p_sigma && mom>0.1 ) p_flag=true;

    if( pip_flag && p_flag ){
      if( mass2<Ppi_mid_mass2 ) return CDS_PiPlus;
      else return CDS_Proton;
    }
    else if( pip_flag ) return CDS_PiPlus;
    else if( p_flag ) return CDS_Proton;
  }
  else{
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pim_flag=true;
    if( kpMass*kpMass-2.5*k_sigma<mass2 && mass2<kpMass*kpMass+2.5*k_sigma && mom<-0.03 ) km_flag=true;

    if( pim_flag ) return CDS_PiMinus;
    else if( km_flag ) return CDS_Kaon;
  }
  return CDS_Other;
}

int MyTools::PIDcorr_tigt(double mom,double mass2)
{
  double pi_sigma = sqrt(4.*piMass*piMass*mom*mom*PID_Param[0]+4.*piMass*piMass*piMass*piMass*PID_Param[1]*(1+piMass*piMass/(mom*mom))+4.*mom*mom*(piMass*piMass+mom*mom)*PID_Param[2]);
  double k_sigma = sqrt(4.*kpMass*kpMass*mom*mom*PID_Param[0]+4.*kpMass*kpMass*kpMass*kpMass*PID_Param[1]*(1+kpMass*kpMass/(mom*mom))+4.*mom*mom*(kpMass*kpMass+mom*mom)*PID_Param[2]);
  double p_sigma = sqrt(4.*pMass*pMass*mom*mom*PID_Param[0]+4.*pMass*pMass*pMass*pMass*PID_Param[1]*(1+pMass*pMass/(mom*mom))+4.*mom*mom*(pMass*pMass+mom*mom)*PID_Param[2]);

  bool pim_flag=false, pip_flag=false, km_flag=false, p_flag=false;
  if( mom>0 ){
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pip_flag=true;
    if( pMass*pMass-2.5*p_sigma<mass2 && mass2<pMass*pMass+2.5*p_sigma && mom>0.1 ) p_flag=true;

    if( pip_flag && p_flag ) return CDS_Other;
    else if( pip_flag ) return CDS_PiPlus;
    else if( p_flag ) return CDS_Proton;
  }
  else{
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pim_flag=true;
    if( kpMass*kpMass-2.5*k_sigma<mass2 && mass2<kpMass*kpMass+2.5*k_sigma && mom<-0.03 ) km_flag=true;

    if( pim_flag && km_flag ) return CDS_Other;
    else if( pim_flag ) return CDS_PiMinus;
    else if( km_flag ) return CDS_Kaon;
  }
  return CDS_Other;
}

bool MyTools::IsElectron(const double &beta, const double &mom)
{
  double x=1./beta;
  if( fabs(mom)>0.15/(x*x) ) return false;
  return true;
}

bool MyTools::IsTarget0(const TVector3 &pos, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
  double z = pos.Z()-fiducial_pos.Z();

  if( -6.0<z && z<6.0 && r<3.5 ) return true;

  return false;
}

bool MyTools::IsTarget1(const TVector3 &pos, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
  double z = pos.Z()-fiducial_pos.Z();

  if( -4.5<z && z<4.5 && r<2.5 ) return true;

  return false;
}

bool MyTools::IsTarget2(const TVector3 &pos, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
  double z = pos.Z()-fiducial_pos.Z();

  if( -4<z && z<4 && r<2. ) return true;

  return false;
}

bool MyTools::IsTarget3(const TVector3 &pos, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
  double z = pos.Z()-fiducial_pos.Z();

  if( -3.<z && z<3. && r<2. ) return true;

  return false;
}

bool MyTools::IsTarget4(const TVector3 &pos, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
  double z = pos.Z()-fiducial_pos.Z();

  if( -2.5<z && z<2.5 && r<1.5 ) return true;

  return false;
}

bool MyTools::IsTarget5(const TVector3 &pos, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
  double z = pos.Z()-fiducial_pos.Z();

  if( -2.0<z && z<2.0 && r<1.0 ) return true;

  return false;
}

bool MyTools::IsTarget0(const std::vector<TVector3> &pos_vec, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  for( int i=0; i<(int)pos_vec.size(); i++ ){
    TVector3 pos = pos_vec[i];
    double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
    double z = pos.Z()-fiducial_pos.Z();

    if( z<-6.0 || 6.0<z || 3.5<r ) return false;
  }
  return true;
}


bool MyTools::IsTarget1(const std::vector<TVector3> &pos_vec, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  for( int i=0; i<(int)pos_vec.size(); i++ ){
    TVector3 pos = pos_vec[i];
    double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
    double z = pos.Z()-fiducial_pos.Z();

    if( z<-4.5 || 4.5<z || 2.5<r ) return false;
  }

  return true;
}

bool MyTools::IsTarget2(const std::vector<TVector3> &pos_vec, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  for( int i=0; i<(int)pos_vec.size(); i++ ){
    TVector3 pos = pos_vec[i];
    double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
    double z = pos.Z()-fiducial_pos.Z();

    if( z<-4.0 || 4.0<z || 2.0<r ) return false;
  }

  return true;
}

bool MyTools::IsTarget3(const std::vector<TVector3> &pos_vec, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  for( int i=0; i<(int)pos_vec.size(); i++ ){
    TVector3 pos = pos_vec[i];
    double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
    double z = pos.Z()-fiducial_pos.Z();

    if( z<-3.0 || 3.0<z || 2.0<r ) return false;
  }

  return true;
}

bool MyTools::IsTarget4(const std::vector<TVector3> &pos_vec, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  for( int i=0; i<(int)pos_vec.size(); i++ ){
    TVector3 pos = pos_vec[i];
    double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
    double z = pos.Z()-fiducial_pos.Z();

    if( z<-2.5 ||2.5<z || 1.5<r ) return false;
  }

  return true;
}


bool MyTools::IsTarget5(const std::vector<TVector3> &pos_vec, ConfMan *conf)
{
  TVector3 fiducial_pos;
  conf->GetGeomMapManager()->GetGPos(CID_Fiducial, 0, fiducial_pos);
  for( int i=0; i<(int)pos_vec.size(); i++ ){
    TVector3 pos = pos_vec[i];
    double r = sqrt((fiducial_pos.X()-pos.X())*(fiducial_pos.X()-pos.X())+(fiducial_pos.Y()-pos.Y())*(fiducial_pos.Y()-pos.Y()));
    double z = pos.Z()-fiducial_pos.Z();

    if( z<-2.0 || 2.0<z || 1.0<r ) return false;
  }

  return true;
}

int MyTools::generation(const int &parentID, MCData *mcData)
{
  if( parentID==0 ){
    return 1;
  }
  int parent=parentID;
  int gen=1;
  while(true){
    bool next=false;
    for(int i=0; i<mcData->trackSize(); i++ ){
      if( mcData->track(i)->trackID()==parent ){
	gen++;
	parent = mcData->track(i)->parentTrackID();
	if( mcData->track(i)->parentTrackID()==0 ){
	  return gen;
	}
	else{
	  next = true;
	  break;
	}
      }
    }
    if( !next ) break;
  }
  return 0;
}

int MyTools::trackStatus(const int &trackID, MCData *mcData)
{
  int track_id=trackID;
  while(true){
    bool next=false;
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData-> track(i);
      if( track->trackID()==track_id ){
	next = true;
	track_id = track->parentTrackID();
	std::cout<<"  Track : pdgID="<<track->pdgID()<<" parentID : "<<track->parentTrackID()<<std::endl;
	if( track->parentTrackID()==0 ) return 0;
	if( track->pdgID()==3222 ) return 3222;
	if( track->pdgID()==3112 ) return 3112;
	if( track->pdgID()==3212 ) return 3212;
	if( track->pdgID()==3122 ) return 3122;
	if( track->pdgID()==3124 ) return 3124;
      }
    }
    if( !next ){
      std::cout<<"  !!! MyTools::isInit not find track id="<<track_id<<" !!!"<<std::endl;
      return -1;
    }
  }
}

int MyTools::nStatus(const int &trackID, MCData *mcData)
{
  int track_id=trackID;
  while(true){
    bool next=false;
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData-> track(i);
      if( track->trackID()==track_id ){
	if( track->parentTrackID()==0 ){
	  //	  std::cout<<" from Initial"<<std::endl;
	  return 0;
	}
	//	std::cout<<"  Track : pdgID="<<track->pdgID()<<" parentID : "<<track->parentTrackID()<<std::endl;

	if( track->pdgID()==2112 ){
	  next = true;
	  track_id = track->parentTrackID();
	}
	else{
	  if( track->pdgID()==3222 || track->pdgID()==3112 || track->pdgID()==3212 ){
	    //	    std::cout<<" decay from Sigma"<<std::endl;
	    return 1;
	  }
	  if( track->pdgID()==3122 ){
	    //	    std::cout<<" decay from Lambda"<<std::endl;
	    return 2;
	  }
	  if( track->pdgID()==3124 ){
	    //	    std::cout<<" decay from Lambda(1520)"<<std::endl;
	    return 3;
	  }
	  //	  std::cout<<" accidental"<<std::endl;
	  return 999; //return error:
	}
      }
    }
    if( !next ){
      std::cout<<"  !!! MyTools::isInit not find track id="<<track_id<<" !!!"<<std::endl;
      return -1;
    }
  }
}

int MyTools::nStatus(const int &seg, DetectorData *detData, MCData *mcData)
{
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData-> detectorHit(i);
    if( hit->detectorID()==CID_NC && hit->channelID()+1==seg ){
      return MyTools::nStatus(hit->trackID(), mcData);
    }
  }
  std::cout<<"  !!! MyTools::nStatus not find DetectorHit !!!"<<std::endl;
  return -2;
}

int MyTools::gammaStatus(const int &trackID, MCData *mcData)
{
  int track_id=trackID;
  while(true){
    bool next=false;
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData-> track(i);
      if( track->trackID()==track_id ){
	if( track->parentTrackID()==0 ){
	  //	  std::cout<<" from Initial"<<std::endl;
	  return 0;
	}
	//	std::cout<<"  Track : pdgID="<<track->pdgID()<<" parentID : "<<track->parentTrackID()<<std::endl;

	if( track->pdgID()==22 ){
	  next = true;
	  track_id = track->parentTrackID();
	}
	else{
	  if( track->pdgID()==11 ){
	    //	    std::cout<<" from e-"<<std::endl;
	    return 1;
	  }
	  if( track->pdgID()==-11 ){
	    //	    std::cout<<" from e+"<<std::endl;
	    return 2;
	  }
	  if( track->pdgID()==3212 ){
	    //	    std::cout<<" decay from Sigma0"<<std::endl;
	    return 3;
	  }
	  if( track->pdgID()==111 ){
	    //	    std::cout<<" decay from pi0"<<std::endl;
	    return 4;
	  }

	  //	  std::cout<<" accidental"<<std::endl;
	  return 999; //return error:
	}
      }
    }
    if( !next ){
      std::cout<<"  !!! MyTools::isInit not find track id="<<track_id<<" !!!"<<std::endl;
      return -1;
    }
  }
}


std::string MyTools::getDetName(const TVector3 &pos)
{
  int id = GeomTools::GetID(pos);
  if( id==160 ) return "Fiducial";
  else if( id==150 ) return "Target";
  else return Form("%d", id);
}

void MyTools::printCDSParams(CDSTrack *track)
{
  int pid=track->PID();
  std::cout<<"===== Print CDS Parameter ";
  if( pid==CDS_PiMinus ) std::cout<<" pi- ";
  else if( pid==CDS_PiPlus   ) std::cout<<" pi+ ";
  else if( pid==CDS_Proton   ) std::cout<<" p   ";
  else if( pid==CDS_Deuteron ) std::cout<<" d   ";
  else if( pid==CDS_Triton   ) std::cout<<" t   ";
  else if( pid==CDS_Helium3  ) std::cout<<" He3 ";
  else if( pid==CDS_Kaon     ) std::cout<<" K-  ";
  else std::cout<<" unknown ";
  std::cout<<"====="<<std::endl;
  int n_param=track->nParamSets();
  for( int i=0; i<n_param; i++ ){
    TVector3 vtx;
    int id;
    double param[5];
    track-> GetNthParameters(i, id, param, vtx);
    std::cout<<"> "<<i<<"  id : "<<id<<"   vtx r="<<vtx.Perp()<<" phi="<<vtx.Phi()<<" z="<<vtx.Z()<<std::endl;
  }
}

