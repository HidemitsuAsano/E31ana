#include "ForwardChargeInfo.h"

#define RUNGE_KUTTA 1

using namespace std;
static const double FC_P_MIN=0.5;
static const double FC_P_MAX=1.5;

ForwardChargeInfo::ForwardChargeInfo() 
  : fVertex(DEFVECT), fFitVtx(DEFVECT), fFitFDC1(DEFVECT), fDiffCounter(DEFAULTD),
    fSeg(-1), fTime(DEFAULTD), fStep(-1), fPCflag(false), fPID(F_Other), 
    fMass2ByAng(DEFAULTD), fMass2ByRK(DEFAULTD), fFlightLengthByArc(DEFAULTD), fFlightLengthByRK(DEFAULTD), fBeta(DEFAULTD), 
    fMomByTOF(DEFAULTD), fMomByAng(DEFAULTD), fMomByRK(DEFAULTD), fCalcTOF(DEFAULTD), fOffset(DEFAULTD), fHitPos(DEFVECT), fMomentum(DEFVECT)
{
}

bool ForwardChargeInfo::calc_forward(BeamLineHitMan *blMan, BeamInfo *beam, CDSInfo *cds, ConfMan *conf)
{
  //  cout<<"==== calc Forward Charge particle by Arc ====="<<endl;

  double beam_tof, beam_out;
  ELossTools::CalcElossBeamTGeo(beam->T0pos(), cds->vertexBeam(), beam->D5mom(), beam->mass(), beam_out, beam_tof);
  TVector3 FDC1pos=fFDC1info.pos();
  double tof=fTime-beam->T0time()-beam_tof;

  fPID=F_Other;
  TVector3 pdir=(fFDC1info.pos()-fVertex).Unit();
  TVector3 uswk_pos=fFDC1info.pos()+((250.-fFDC1info.pos().Z())/pdir.Z())*pdir;
  
  TVector3 counter_pos;
  if( fPCflag ) conf->GetGeomMapManager()->GetGPos(CID_PC, fSeg, counter_pos);
  else          conf->GetGeomMapManager()->GetGPos(CID_CVC, fSeg, counter_pos);
  fHitPos=counter_pos;  

  TVector3 pdir2=counter_pos-uswk_pos;
  double angle=pdir.Angle(pdir2);
  
  double leff=100.;
  double r=leff/sqrt(2*(1-cos(angle)));
  double larc=r*angle;
  double line=2*leff/sqrt(2*(1+cos(angle)));
  double diff=line-larc;

  fFlightLengthByArc=(uswk_pos-fVertex).Mag()+(counter_pos-uswk_pos).Mag()-diff-1.5;
  fMomByAng=r*0.9782*Const/100.+0.023;
  fBeta =fFlightLengthByArc/(tof*Const*100.);
  fMomentum=FDC1pos-fVertex;
  fMomentum.SetMag(fMomByAng);

  fMass2ByAng=findMass2(fFlightLengthByArc, fMomByAng, tof);

  return true;
}

double ForwardChargeInfo::findMass2(const double &fl, const double &mom, const double &tof)
{
  double beta=fl/(tof*Const*100.);
  double mass2=mom*mom*(1./(beta*beta)-1);
  //  cout<<beta<<"  "<<mass2<<endl;
  if( mass2<0 ) return mass2;

  double mom_out, tmp_tof, tmp_tof2;
  //  cout<<beta<<"  "<<sqrt(mass2)<<endl;
  double mass=sqrt(mass2);
  if( !ELossTools::CalcElossForwardTGeo(fVertex, fFDC1info.pos(), fl, mom, mass, mom_out, tmp_tof) ){
    // cout<<"!!!!! ELossTools::CalcElossForwardTGeo return false    Break !!!!!"<<endl;
    // exit(0);
    return mass2;
  }

  double tmp_mass2=mass2;
  while( true ){
    tmp_mass2*=0.9;
    mass=sqrt(tmp_mass2);
    if( !ELossTools::CalcElossForwardTGeo(fVertex, fFDC1info.pos(), fl, mom, mass, mom_out, tmp_tof2) ){
      cout<<"!!!!! ELossTools::CalcElossForwardTGeo return false !!!!!"<<endl;
      //      cout<<"!!!!! ELossTools::CalcElossForwardTGeo return false   Break !!!!!"<<endl;
      cout<<"mass2 : "<<mass2<<"  "<<tmp_mass2<<endl;
      cout<<"Step2"<<endl;
      return 0.5*(tmp_mass2+mass2);
      //      exit(0);
    }
    if( (tmp_tof-tof)*(tmp_tof2-tof)<0 ) break;
  }
  // cout<<"2nd step   TOF : "<<tof<<" "<<tmp_tof2<<endl;
  // cout<<"   mass2 : "<<mass2<<"  "<<tmp_mass2<<endl;

  while( true ){
    double mass2_mid=0.5*(mass2+tmp_mass2);
    mass=sqrt(mass2_mid);    
    double tmp_tof_mid;
    if( !ELossTools::CalcElossForwardTGeo(fVertex, fFDC1info.pos(), fl, mom, mass, mom_out, tmp_tof_mid) ){
      cout<<"!!!!! ELossTools::CalcElossForwardTGeo return false !!!!!"<<endl;
      cout<<"mass2 : "<<mass2<<"  "<<tmp_mass2<<endl;
      cout<<"Step2"<<endl;
      return 0.5*(tmp_mass2+mass2);
      // cout<<"!!!!! ELossTools::CalcElossForwardTGeo return false    Break !!!!!"<<endl;
      // cout<<"mass2 : "<<mass2<<"  "<<tmp_mass2<<endl;
      // cout<<"Find Bynari "<<endl;
      //      exit(0);
    }

    if( (tmp_tof_mid-tof)*(tmp_tof-tof)>0 ) mass2=mass2_mid;
    else tmp_mass2=mass2_mid;
 
    if( fabs(tmp_mass2-mass2)<1.0e-8 ){
      // cout<<"binary search finish"<<endl;
      // cout<<"   TOF   : "<<tof<<" "<<tmp_tof2<<endl;
      // cout<<"   mass2 : "<<mass2<<" "<<tmp_mass2<<endl;
      mass2=0.5*(tmp_mass2+mass2);
      break;
    }
  }

  return mass2;
}

bool ForwardChargeInfo::fit_forward(BeamLineHitMan *blMan, BeamInfo *beam, CDSInfo *cds, ConfMan *conf, bool refit)
{
  // cout<<"===== ForwardChargeInfo::fit_forward ====="<<endl;
  // cout<<"      > step : "<<fStep<<endl;
  double beam_tof, beam_out;
  ELossTools::CalcElossBeamTGeo(beam->T0pos(), beam->vertex(), beam->D5mom(), beam->mass(), beam_out, beam_tof);
  TVector3 FDC1pos=fFDC1info.pos();
  HodoscopeLikeHit *hit=hodo(blMan);
  //  if( hit ) fTime=hit->ctmean();
  if( hit ) fTime=hit->ctmean();
  double tof=fTime-beam->T0time()-beam_tof;
  if( fStep<0 ){
    if( !calc_forward(blMan, beam, cds, conf) ){ return false; }
    fStep=0;
  }
#if RUNGE_KUTTA
  if( fStep<1 ){
    // cout<<"===== ForwardChargeInfo::fit_forward ====="<<endl;
    // cout<<"      > Fit w/o PID (Eloss)"<<endl;

    TVector3 vtx=fVertex;
    TVector3 mom=FDC1pos-fVertex;
    mom.SetMag(fMomByAng);
    if( !ProtonArm::fit(0.0, vtx, mom, FDC1pos, hodo(blMan)) ){
      cout<<"  !!!  ProtonArm::fit fault !!!"<<endl;
    }
    else{
      ChargeParticle particle=ProtonArm::shootChargeParticle(particleMass[fPID], vtx, mom, conf);
      if( !particle.flag() ){
	cout<<"!!!!! particle.flag error !!!!!"<<endl;
	exit(0);
      }
      fFitFDC1=particle.FDC1pos();
      fDiffCounter=particle.diffCounter();

      //      particle.dump();
      fFlightLengthByRK=particle.fl();
      //      cout<<particleMass[fPID]<<"  "<<fFlightLengthByArc<<"  "<<fFlightLengthByRK<<endl;
      //      fBeta =fFlightLengthByRK/(tof*Const*100.);
      fFitVtx=vtx;
      fMomByRK=mom.Mag();
      //      fMass2ByRK=fMomByRK*fMomByRK*(1./(fBeta*fBeta)-1);
      fMass2ByRK=findMass2(fFlightLengthByRK, fMomByRK, tof);
      fMomentum=mom;
      fHitPos=particle.pos();

      if( -0.5<fMass2ByRK && fMass2ByRK<0.1 ) fPID=F_Pion;
      if( 0.1<fMass2ByRK && fMass2ByRK<2.5 ) fPID=F_Proton;
      if( 2.5<fMass2ByRK && fMass2ByRK<5.0 ) fPID=F_Deuteron;

      // cout<<"> Beta : "<<fBeta<<endl;      
      // cout<<Form("vtx (%lf, %lf, %lf)", vtx.x(), vtx.y(), vtx.z())<<endl;
      // cout<<Form("mom (%lf, %lf, %lf)", mom.x(), mom.y(), mom.z())<<endl;
      // cout<<"> mom by RK      : "<<fMomByRK<<" [GeV/c]"<<endl;
      // cout<<"> mass2 by RK    : "<<fMass2ByRK<<" [(GeV/c2)2]"<<endl;
      // cout<<"> fl : "<<fFlightLengthByRK<<" [cm]"<<endl;
      
      fStep=1;
    }
  }
#else
  if( -0.5<fMass2ByAng && fMass2ByAng<0.1 ) fPID=F_Pion;
  if( 0.1<fMass2ByAng && fMass2ByAng<2.5 ) fPID=F_Proton;
  if( 2.5<fMass2ByAng && fMass2ByAng<5.0 ) fPID=F_Deuteron;
#endif
  TString parName;
  if( fPID==F_Pion ) parName="pion";
  else if( fPID==F_Proton ) parName="proton";
  else if( fPID==F_Deuteron ) parName="proton";
  if( parName.Length()==0 ) return false;
#if 1
  if( fStep<2 ){
    // cout<<"===== ForwardChargeInfo::fit_forward ====="<<endl;
    // cout<<"      > Fit particle : "<<parName<<endl;

    TVector3 vtx=fVertex;
    TVector3 mom=FDC1pos-fVertex;
    mom.SetMag(fMomByRK);
    if( !ProtonArm::fit(particleMass[fPID], vtx, mom, FDC1pos, hodo(blMan)) ){
      //      cout<<"  !!!  ProtonArm::fit fault "<<parName<<" !!!"<<endl;
      return false;
    }
    else{
      ChargeParticle particle=ProtonArm::shootChargeParticle(particleMass[fPID], vtx, mom, conf);
      if( !particle.flag() ){
	cout<<"!!!!! particle.flag error !!!!!"<<endl;
	exit(0);
      }
      fFitFDC1=particle.FDC1pos();
      fDiffCounter=particle.diffCounter();

      fFlightLengthByRK=particle.fl();
      fCalcTOF=particle.time();
      //      cout<<particleMass[fPID]<<"  "<<fFlightLengthByArc<<"  "<<fFlightLengthByRK<<endl;
      //      fBeta =fFlightLengthByRK/(tof*Const*100.);
      // cout<<"Fit mom : "<<mom.X()<<", "<<mom.Y()<<", "<<mom.Z()<<endl;
      // cout<<"Particle mom : "<<particle.mom().X()<<", "<<particle.mom().Y()<<", "<<particle.mom().Z()<<endl;
      fFitVtx=vtx;
      fMomByRK=mom.Mag();
      //      fMass2ByRK=findMass2(fFlightLengthByRK, fMomByRK, tof);
      fMomentum=mom;
	
      // cout<<Form("vtx (%lf, %lf, %lf)", vtx.x(), vtx.y(), vtx.z())<<endl;
      // cout<<Form("mom (%lf, %lf, %lf)", mom.x(), mom.y(), mom.z())<<endl;
      // cout<<"> Particle : "<<parName<<endl;
      // cout<<"> mom by RK      : "<<fMomByRK<<" [GeV/c]"<<endl;
      // cout<<"> mass2 by RK    : "<<fMass2ByRK<<" [(GeV/c2)2]"<<endl;
      // cout<<"> fl : "<<fFlightLengthByRK<<" [cm]"<<endl;

      fStep=2;
    }
  }
#endif

  if( fPID==F_Proton || fPID==F_Deuteron ){
    if( refit ){
      // cout<<"===== ForwardChargeInfo::fit_forward ====="<<endl;
      // cout<<"      > reFit particle : "<<parName<<endl;

      TVector3 vtx=fVertex;
      TVector3 mom=FDC1pos-fVertex;

      //      if( !ProtonArm::fine_fit(particleMass[fPID], vtx, mom, FDC1pos, hodo(blMan)) ){
	//	cout<<"  !!!  ProtonArm::fine_fit fault "<<parName<<" !!!"<<endl;
      // 	return false;
      // }
      // std::cout<<"FL by Arc ; "<<fFlightLengthByArc<<std::endl;
      // std::cout<<"FL by RK  ; "<<fFlightLengthByRK<<std::endl;
      // std::cout<<"befor calc : "<<fMomByTOF<<std::endl;
    }

    calcMomByTOF(beam, cds);
    //    std::cout<<"after calc : "<<fMomByTOF<<std::endl;
  }
  // cout<<"===== fit Forward succeesfully finish ====="<<endl;
  // dump();

  return true;
}

bool ForwardChargeInfo::calcMomByTOF(BeamInfo *beam, CDSInfo *cds)
{
  double beam_tof, beam_out;
  TVector3 FDC1pos=fFDC1info.pos();
  ELossTools::CalcElossBeamTGeo(beam->T0pos(), beam->vertex(), beam->D5mom(), beam->mass(), beam_out, beam_tof);
  double tof=fTime-beam->T0time()-beam_tof;

  // std::cout<<"   T0 time : "<<beam->T0time()<<std::endl;
  // std::cout<<"   start time : "<<-beam->T0time()-beam_tof<<std::endl;
  // std::cout<<"   stop  time : "<<fTime<<std::endl;

  double fl=-999;
  if( fStep==0 ){
    fBeta =fFlightLengthByArc/(tof*Const*100.);
    fl=fFlightLengthByArc;
  }
  else if( fStep>0 ){
    fBeta =fFlightLengthByRK/(tof*Const*100.);
    fl=fFlightLengthByRK;
  }
  else return false;

  //  if( fFlightLengthByArc>0 ) fl=fFlightLengthByArc;

  double mom1=particleMass[fPID]/sqrt(1./(fBeta*fBeta)-1.);
  double mom2=mom1;
  double tmpmom, tmptof;
  //      ELossTools::CalcElossForwardTGeo(fFitVtx, FDC1pos, fFlightLength, mom1, particleMass[fPID], tmpmom,tmptof);
  ELossTools::CalcElossForwardTGeo(fVertex, FDC1pos, fl, mom1, particleMass[fPID], tmpmom, tmptof);
  double tmp_diff=tmptof-tof;
  double diff;

  for( int i=0; i<100; i++ ){
    mom2+=0.001;
    //	ELossTools::CalcElossForwardTGeo(fFitVtx, FDC1pos, fFlightLength, mom2, particleMass[fPID], tmpmom,tmptof);
    ELossTools::CalcElossForwardTGeo(fVertex, FDC1pos, fl, mom2, particleMass[fPID], tmpmom,tmptof);
    diff=tmptof-tof;
    
    //	cout<<i<<"  "<<tmptof<<"  "<<tof<<"  "<<tmptof*tof<<"  "<<mom1<<"  "<<mom2<<endl;
    if( i==99 ){
      cout<<"!!!!! calcMomByTOF iteration over !!!!!"<<endl;
    }
    if( diff*tmp_diff<0 ) break;
  }
  //  cout<<"calc_tof0 "<<tmptof<<"  "<<tof<<"  "<<tmptof*tof<<"  "<<mom1<<"  "<<mom2<<endl;
  
  for( int i=0; i<10000; i++ ){
    double mommid=0.5*(mom1+mom2);
    //	ELossTools::CalcElossForwardTGeo(fFitVtx, FDC1pos, fFlightLength, mommid, particleMass[fPID], tmpmom,tmptof);
    if( mom2<mom1 ){
      //      std::cout<<"     mom 1 : "<<mom1<<"  mom2 : "<<mom2<<std::endl;
      std::swap(mom1, mom2);
    }

    ELossTools::CalcElossForwardTGeo(fVertex, FDC1pos, fl, mommid, particleMass[fPID], tmpmom, tmptof);
    diff=tmptof-tof;
    
    if( diff>0 ) mom1=mommid;
    else mom2=mommid;

    //	cout<<i<<"  "<<tmptof<<"  "<<tof<<"  "<<fabs(mom1-mom2)<<"  "<<mom1<<"  "<<mom2<<endl;
    if( i==9999 ){
      cout<<"!!!!! calcMomByTOF iteration over !!!!!"<<endl;
    }
    if( fabs(mom1-mom2)<0.00001 && fabs(tof-tmptof)<0.01 ){
      //	  cout<<"calc_tof1 "<<i<<"  "<<tmptof<<"  "<<tof<<"  "<<tmptof*tof<<"  "<<mom1<<"  "<<mom2<<endl;
      break;
    }
  }
  // std::cout<<"   tof : "<<tof<<"  tmptof : "<<tmptof<<std::endl;
  // std::cout<<"   mom1 : "<<mom1<<"  mom2 : "<<mom2<<std::endl;

  //      fMomentum=FDC1pos-fFitVtx;
  ELossTools::CalcElossForwardTGeo(fVertex, FDC1pos, fl, fMomByRK, particleMass[fPID], tmpmom, fCalcTOF);
  fMomentum=FDC1pos-fVertex;
  fMomentum.SetMag(0.5*(mom1+mom2));
  fMomByTOF=0.5*(mom1+mom2);
  
  // cout<<"> mom by angle   : "<<fMomByAng<<" [GeV/c]"<<endl;
  // cout<<"> mom by RK      : "<<fMomByRK<<" [GeV/c]"<<endl;
  // cout<<"> mom by TOF     : "<<fMomByTOF<<" [GeV/c]"<<endl;
  // string input;
  // cin>>input;
  // if( input=="q" ) exit(-1);
  fStep=3;

  return true;
}

double ForwardChargeInfo::angle(ConfMan *conf) const
{
  TVector3 hitpos=fHitPos;
  if( conf ){
    if( fPCflag ) conf->GetGeomMapManager()->GetGPos(CID_PC, fSeg, hitpos);
    else          conf->GetGeomMapManager()->GetGPos(CID_CVC, fSeg, hitpos);
  }

  TVector3 pdir=(fFDC1info.pos()-fVertex).Unit();
  TVector3 uswk_pos=fFDC1info.pos()+((250.-fFDC1info.pos().Z())/pdir.Z())*pdir;
  TVector3 pdir2=hitpos-uswk_pos;
  double angle=pdir.Angle(pdir2);
  return angle;
}

void ForwardChargeInfo::SetHodo(HodoscopeLikeHit *hit)
{
  if( hit->cid()==CID_PC ) fPCflag=true;
  else fPCflag=false;
  fSeg=hit->seg();
  fTime=hit->ctmean();
}

HodoscopeLikeHit* ForwardChargeInfo::hodo(BeamLineHitMan *blMan) const
{
  if( fPCflag ){
    for( int i=0; i<blMan->nPC(); i++ ){
      if( fSeg==blMan->PC(i)->seg() ) return blMan->PC(i);
    }
  }
  else{
    for( int i=0; i<blMan->nCVC(); i++ ){
      if( fSeg==blMan->CVC(i)->seg() ) return blMan->CVC(i);
    }
  }
  return 0;
}

bool ForwardChargeInfo::calc_simple_beam_through(const BeamInfo &beam, ConfMan *conf, BeamLineHitMan *blMan)
{
  double mass=beam.mass();
  TVector3 T0pos=beam.T0pos();

  TVector3 pdir=(fFDC1info.pos()-T0pos).Unit();
  TVector3 uswk_pos=fFDC1info.pos()+((250.-fFDC1info.pos().Z())/pdir.Z())*pdir;
  TVector3 counter_pos;

  if( fPCflag ) conf->GetGeomMapManager()->GetGPos(CID_PC, fSeg, counter_pos);
  else          conf->GetGeomMapManager()->GetGPos(CID_CVC, fSeg, counter_pos);

  TVector3 pdir2=counter_pos-uswk_pos;
  double angle=pdir.Angle(pdir2);

  double leff=100.;
  double r=leff/sqrt(2*(1-cos(angle)));
  double larc=r*angle;
  double line=2*leff/sqrt(2*(1+cos(angle)));
  double diff=line-larc;
  fFlightLengthByArc=(uswk_pos-T0pos).Mag()+(counter_pos-uswk_pos).Mag()-diff;
  fFlightLengthByArc-=1.5;
  fHitPos.SetXYZ(counter_pos.X(), counter_pos.Y(), fFlightLengthByRK*pdir.Y()/pdir.Z());

  bool flag=ELossTools::CalcElossForwardTGeo(T0pos, fFDC1info.pos(), fFlightLengthByArc, beam.D5mom(), mass, fMomByAng, fCalcTOF);
  if( flag ) fStep=0;
  return flag;
}

bool ForwardChargeInfo::calc_simple_p_through(const BeamInfo &beam, ConfMan *conf, BeamLineHitMan *blMan)
{
  HodoscopeLikeHit *FChit=hodo(blMan);
  double mass=beam.mass();
  TVector3 T0pos=beam.T0pos();

  TVector3 pdir=(fFDC1info.pos()-T0pos).Unit();
  TVector3 uswk_pos=fFDC1info.pos()+((250.-fFDC1info.pos().Z())/pdir.Z())*pdir;
  TVector3 counter_pos;

  if( fPCflag ) conf->GetGeomMapManager()->GetGPos(CID_PC, fSeg, counter_pos);
  else          conf->GetGeomMapManager()->GetGPos(CID_CVC, fSeg, counter_pos);

  TVector3 pdir2=counter_pos-uswk_pos;
  double angle=pdir.Angle(pdir2);

  double leff=100.;
  double r=leff/sqrt(2*(1-cos(angle)));
  double larc=r*angle;
  double line=2*leff/sqrt(2*(1+cos(angle)));
  double diff=line-larc;
  fFlightLengthByArc=(uswk_pos-T0pos).Mag()+(counter_pos-uswk_pos).Mag()-diff;
  fFlightLengthByArc-=1.5;
  fHitPos.SetXYZ(counter_pos.X(), counter_pos.Y(), fFlightLengthByRK*pdir.Y()/pdir.Z());

  bool flag=ELossTools::CalcElossForwardTGeo(T0pos, fFDC1info.pos(), fFlightLengthByArc, beam.D5mom(), mass, fMomByAng, fCalcTOF);
  if( flag ) fStep=0;

  double tof=FChit->ctmean()-beam.T0time();

  fBeta=fFlightLengthByArc/(tof*100.*Const);

  double mom1=particleMass[fPID]/sqrt(1./(fBeta*fBeta)-1.);
  double mom2=mom1;
  double tmpmom, tmptof;

  TVector3 FDC1pos=fFDC1info.pos();
  ELossTools::CalcElossForwardTGeo(beam.T0pos(), FDC1pos, fFlightLengthByArc, mom1, particleMass[fPID], tmpmom, tmptof);
  double tmp_diff=tmptof-tof;

  // std::cout<<"                       tof : "<<tof<<std::endl;
  // std::cout<<"step0  mom : "<<mom1<<"  tof : "<<tmptof<<std::endl;
  for( int i=0; i<100; i++ ){
    mom2+=0.001;

    ELossTools::CalcElossForwardTGeo(beam.T0pos(), FDC1pos, fFlightLengthByArc, mom2, particleMass[fPID], tmpmom,tmptof);
    diff=tmptof-tof;

    if( i==99 ){
      cout<<"!!!!! calcMomByTOF iteration over !!!!!"<<endl;
    }
    if( diff*tmp_diff<0 ) break;
  }
  //  std::cout<<"step1  mom : "<<0.5*(mom1+mom2)<<"  tof : "<<tmptof<<std::endl;
  //  std::cout<<"mom 1 : "<<mom1<<"  mom2 : "<<mom2<<std::endl;

  for( int i=0; i<10000; i++ ){
    if( mom2<mom1 ){
      //      std::cout<<"mom 1 : "<<mom1<<"  mom2 : "<<mom2<<std::endl;
      std::swap(mom1, mom2);
    }
    double mommid=0.5*(mom1+mom2);
    ELossTools::CalcElossForwardTGeo(beam.T0pos(), FDC1pos, fFlightLengthByArc, mommid, particleMass[fPID], tmpmom,tmptof);
    diff=tmptof-tof;

    if( diff>0 ) mom1=mommid;
    else mom2=mommid;

    if( i==9999 ){
      cout<<"!!!!! calcMomByTOF iteration over !!!!!"<<endl;
    }
    //    if( fabs(mom1-mom2)<0.00001 && fabs(tof-tmptof)<0.01 ){
    if( fabs(mom1-mom2)<0.000001 ){
      break;
    }
  }
  //  std::cout<<"step2  mom : "<<0.5*(mom1+mom2)<<"  tof : "<<tmptof<<std::endl;

  // std::cout<<tof<<std::endl;
  // std::cout<<tmptof<<std::endl;

  fCalcTOF=tmptof;
  fMomentum=FDC1pos-beam.T0pos();
  fMomentum.SetMag(0.5*(mom1+mom2));
  fMomByTOF=0.5*(mom1+mom2);

  return flag;
}

vector<vector<double> > ForwardChargeInfo::calc_beam_through_wUSWK(const BeamInfo &beam, ConfMan *conf, BeamLineHitMan *blMan)
{
  HodoscopeLikeHit *FChit=hodo(blMan);
  std::cout<<Form("FDC1pos : (%lf, %lf, %lf)",  fFDC1info.pos().X(), fFDC1info.pos().Y(), fFDC1info.pos().Z())<<std::endl;
  TVector3 mom=fFDC1info.pos()-beam.BPCpos();
  mom.SetMag(beam.D5mom());
  TVector3 T0pos=beam.BPCpos()+((beam.T0pos().Z()-beam.BPCpos().Z())/mom.Z())*mom;
  vector<vector<double> > trajectory;
  ChargeParticle particle=ProtonArm::shootChargeParticle(particleMass[fPID], T0pos, mom, conf, hodo(blMan), trajectory);
  // if( !particle.flag() ){
  //   cout<<"  ForwardChargeInfo::calc_beam_throught_wUSWK"<<endl;
  //   cout<<"  !!! ProtonArm::shootChargeParticle  false !!!"<<endl;
  // }
  fFlightLengthByRK=particle.fl();
  fCalcTOF=particle.time();
  fHitPos=particle.pos();
  fMomByAng=beam.D5mom();

  double tof=FChit->ctmean()-beam.T0time();

  fBeta=fFlightLengthByRK/(tof*100.*Const);

  double mom1=particleMass[fPID]/sqrt(1./(fBeta*fBeta)-1.);
  double mom2=mom1;
  double tmpmom, tmptof;

  TVector3 FDC1pos=fFDC1info.pos();
  ELossTools::CalcElossForwardTGeo(beam.T0pos(), FDC1pos, fFlightLengthByRK, mom1, particleMass[fPID], tmpmom, tmptof);
  double tmp_diff=tmptof-tof;
  double diff;

  for( int i=0; i<100; i++ ){
    mom2+=0.001;

    ELossTools::CalcElossForwardTGeo(beam.T0pos(), FDC1pos, fFlightLengthByRK, mom2, particleMass[fPID], tmpmom,tmptof);
    diff=tmptof-tof;

    if( i==99 ){
      cout<<"!!!!! calcMomByTOF iteration over !!!!!"<<endl;
    }
    if( diff*tmp_diff<0 ) break;
  }

  for( int i=0; i<100; i++ ){
    double mommid=0.5*(mom1+mom2);
    ELossTools::CalcElossForwardTGeo(beam.T0pos(), FDC1pos, fFlightLengthByRK, mom2, particleMass[fPID], tmpmom,tmptof);
    diff=tmptof-tof;

    if( diff>0 ) mom1=mommid;
    else mom2=mommid;

    if( i==99 ){
      cout<<"!!!!! calcMomByTOF iteration over !!!!!"<<endl;
    }
    if( fabs(mom1-mom2)<0.00001 ){
      break;
    }
  }
  std::cout<<tof<<std::endl;
  std::cout<<tmptof<<std::endl;

  fMomentum=FDC1pos-beam.T0pos();
  fMomentum.SetMag(0.5*(mom1+mom2));
  fMomByTOF=0.5*(mom1+mom2);

  return trajectory;
}

void ForwardChargeInfo::dump() const
{
  cout<<"===== Forward Charge Dump ====="<<endl; 
  if( fPID==F_Kaon ) cout<<"> Kaon"<<endl;
  else if( fPID==F_Pion     ) cout<<"> Pion"<<endl;
  else if( fPID==F_Proton   ) cout<<"> Proton"<<endl;
  else if( fPID==F_Deuteron ) cout<<"> Deuteron"<<endl;
  else if( fPID==F_Neutron  ) cout<<"> Neutron"<<endl;
  else if( fPID==F_Gamma    ) cout<<"> Gamma"<<endl;
  else if( fPID==F_Other    ) cout<<"> Other"<<endl;
  else cout<<"> unknown"<<endl;

  if( fPCflag ) cout<<"PC seg"<<fSeg;
  else          cout<<"CVC seg"<<fSeg;
  cout<<" "<<fTime<<"[ns]"<<endl;
  cout<<"Vertex   : "<<Form("(%6.4lf, %6.4lf, %6.4lf)", fVertex.X(), fVertex.Y(), fVertex.Z())<<endl;
  cout<<"FDC1 pos : "<<Form("(%6.4lf, %6.4lf, %6.4lf)", fFDC1info.pos().X(), fFDC1info.pos().Y(), fFDC1info.pos().Z())<<endl;
  cout<<"FDC1 dir : "<<Form("(%6.4lf, %6.4lf, %6.4lf)", fFDC1info.dir().X(), fFDC1info.dir().Y(), fFDC1info.dir().Z())<<endl;
  cout<<"Flight Length  Arc : "<<fFlightLengthByArc<<" [cm]   RK : "<<fFlightLengthByRK<<endl;
  cout<<"beta : "<<fBeta<<" [cm/ns]"<<endl;
  cout<<"mass2 Ang : "<<fMass2ByAng<<" [(GeV/c2)2]"<<"  RK : "<<fMass2ByRK<<" [(GeV/c2)2]"<<endl;
  cout<<"mom Ang : "<<fMomByAng<<" [GeV/c]"<<"  RK : "<<fMomByRK<<" [GeV/c]  TOF : "<<fMomByTOF<<" [GeV/c]"<<endl;

  cout<<">>> Fit Info"<<endl;
  cout<<"> Vtx : "<<(fFitVtx-fVertex).Mag()<<endl;
  cout<<"> FDC1 : "<<(fFDC1info.pos()-fFitFDC1).Mag()<<endl;
  cout<<"> Counter : "<<fDiffCounter<<endl;
}

//******************************************//
//*                                        *//      
//*   for Check Position by Graphic        *//
//*                                        *//
//******************************************//
vector<vector<double> > ForwardChargeInfo::fit_beam_through_wPoints(const BeamInfo &beam, ConfMan *conf)
{
  vector<vector<double> > trajectory;
  TVector3 BPCpos=beam.BPCpos();
  TVector3 momentum=fFDC1info.pos()-BPCpos;
  TVector3 T0gpos;
  conf-> GetGeomMapManager()-> GetGPos(CID_T0, 1, T0gpos);
  momentum.SetMag(beam.D5mom());
  TVector3 T0pos=BPCpos+((T0gpos.Z()-0.5-BPCpos.Z())/momentum.Z())*momentum;

  if( !ProtonArm::fitBT(particleMass[fPID], beam.D5mom(), beam.BPCpos(), fFDC1info.pos(), conf) ) return trajectory;
  double par[5];
  double err[5];
  for( int i=0; i<gMinuit->GetNumPars(); i++ ){
    TString name; double bnd1, bnd2; int flag;
    gMinuit-> mnpout(i, name, par[i], err[i], bnd1, bnd2, flag);
  }

  // cout<<"> Initial parameter : "<<Form("%lf  %lf  %lf  %lf  %lf", T0pos.X(), T0pos.Y(), momentum.X(), momentum.Y(), momentum.Z())<<endl;
  // cout<<"> Fitted parameter  : "<<Form("%lf  %lf  %lf  %lf  %lf", par[0], par[1], par[2], par[3], par[4])<<endl;

  TVector3 fit_pos(par[0], par[1], T0pos.Z());
  TVector3 fit_mom(par[2], par[3], par[4]);
  ChargeParticle particle=ProtonArm::shootChargeParticle(particleMass[fPID], fit_pos, fit_mom, conf, trajectory);
  if( !particle.flag() ){
    cout<<"  !!! fit_beam_through !!!"<<endl;
    cout<<"  !!! ProtonArm::shootChargeParticle  false !!!"<<endl;
  }
  fFlightLengthByRK=particle.fl();
  fCalcTOF=particle.time();
  fHitPos=particle.pos();
  fMomByAng=fit_mom.Mag();

  return trajectory;
}
