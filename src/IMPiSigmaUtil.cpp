//author: Hidemitsu Asano
//email : hidemitsu.asano@riken.jp
//
//These are helper functions for EventAnalysisIMPiSigma.cpp and UserSimIMPiSigma.cpp to avoid many duplications of codes.
//All functions are defined in the namespace::Util


#include "IMPiSigmaUtil.h"
#include "Tools.h"
#include "TrackTools.h"



//analyze # of CDH hits within the time range and judge if it is cdscuts::cdhmulti defined in IMPiSigmaAnaPar.h
int Util::GetCDHMul(CDSHitMan *cdsman, const int ntrack, const bool MCFlag)
{
  //** # of CDH-hits cut **//
  int nCDH = 0;
  for( int i=0; i<cdsman->nCDH(); i++ ) {
    Tools::Fill2D(Form("CDHtime"),cdsman->CDH(i)->seg(),cdsman->CDH(i)->ctmean());
    Tools::Fill2D(Form("dE_CDHtime"), cdsman->CDH(i)->ctmean(), cdsman->CDH(i)->emean());
    if(ntrack == cdscuts::cds_ngoodtrack){
      Tools::Fill2D(Form("dE_CDHtime_2track"), cdsman->CDH(i)->ctmean(), cdsman->CDH(i)->emean());
    }
    //if( cdsman->CDH(i)->CheckRange() ) nCDH++; //** only requirement of TDC **//
    if(MCFlag){
      if( cdsman->CDH(i)->CheckRange() && cdsman->CDH(i)->ctmean()<(cdscuts::tdc_cdh_max+cdscuts::tdc_simoffset) ) nCDH++;
    }else{
      if( cdsman->CDH(i)->CheckRange() && cdsman->CDH(i)->ctmean()<cdscuts::tdc_cdh_max ) nCDH++;
    }
  }
  Tools::Fill1D( Form("mul_CDH"), nCDH );

  return nCDH;
}

bool Util::IsForwardCharge(BeamLineHitMan *blman)
{
  int nBVC = 0;
  int nCVC = 0;
  int nPC  = 0;
  for( int i=0; i<blman->nBVC(); i++ ) {
    if( blman->BVC(i)->CheckRange() ) nBVC++;
  }
  for( int i=0; i<blman->nTOF(); i++ ) {
    if( blman->TOF(i)->CheckRange() ) nCVC++;
  }
  for( int i=0; i<blman->nPC(); i++ ) {
    if( blman->PC(i)->CheckRange() ) nPC++;
  }
  Tools::Fill1D( Form("mul_BVC"), nBVC );
  Tools::Fill1D( Form("mul_CVC"), nCVC );
  Tools::Fill1D( Form("mul_PC"),  nPC );
  if( nBVC || nCVC || nPC ) return true;
  else return false;
}



//returns # of neighboring hits of CDH segments listed in vector cdhseg
//used for Isolation cuts of neutral particle candidates
int Util::GetCDHNeighboringNHits(const std::vector <int> &seg, const std::vector <int> &allhit )
{
  int NNeighboringHits=0;
  for( int ineuseg=0; ineuseg<(int)seg.size(); ineuseg++ ) {
    for( int ihit=0; ihit<(int)allhit.size(); ihit++ ) {
      if( seg[ineuseg]-allhit[ihit] ) {
        Tools::Fill1D( Form("diff_CDH"), seg[ineuseg]-allhit[ihit] );
      }
      //CDH has 36 segments. Requiring there is no hits on neighboring segments.
      if( abs(seg[ineuseg]-allhit[ihit])==1 || abs(seg[ineuseg]-allhit[ihit])==35 )
        NNeighboringHits++;
    }
  }

  return NNeighboringHits;
}


//returns # of cdc hits of layer 15 and 16 within -+ 15 degree of the CDH hit
int Util::GetNHitsCDCOuter(const TVector3 PosCDH, CDSHitMan *cdsman)
{
  int nCDChit = 0;
  const double PhiMin = -15.0/360.*TwoPi; // rad
  const double PhiMax =  15.0/360.*TwoPi; // rad
  for( int ilr=14; ilr<16; ilr++ ) { // charge veto using layer 15, 16
    for( int icdchit=0; icdchit<cdsman->nCDC(ilr); icdchit++ ) {
      CDCHit *cdc=cdsman->CDC(ilr,icdchit);
      TVector3 Pos_CDC = cdc->wpos();
      Pos_CDC.SetZ(0); // only xy pos is used
      const double angle = Pos_CDC.Angle(PosCDH); // rad
      //  std::cerr<<"CDC "<<ilr<<" "<<icdchit<<" "<<cdc->wire()<<" -> "<<Pos_CDC.Phi()/TwoPi*360.
      //    <<" deg :: diff = "<<angle/TwoPi*360<<" deg"<<std::endl;
      //
      Tools::Fill1D( Form("diff_CDH_CDC"), angle/TwoPi*360 );
      if( PhiMin<angle && angle<PhiMax ) nCDChit++;
    }//icdchit
  }//ilr
  
  return nCDChit;
}


//BLC1-D5-BLC2 momentum analysis
double Util::AnaBeamSpec(ConfMan *confman, BeamLineTrackMan *bltrackman,const int blc1id, const int blc2id)
{
  BeamSpectrometer *beamsp = new BeamSpectrometer( confman );
  LocalTrack *blc1 = bltrackman->trackBLC1(blc1id);
  LocalTrack *blc2 = bltrackman->trackBLC2(blc2id);
  beamsp->TMinuitFit( blc1, blc2, confman );
  const double beammom = beamsp->mom();
  const double bchi = beamsp->chisquare();
  delete beamsp;

  Tools::Fill1D( Form("trackchi2_beam"), bchi );
  if( bchi>blcuts::d5_chi2_max ) {
    return -9999.;
  }
  Tools::Fill1D( Form("momentum_beam"), beammom );

  return beammom;
}


//CDSChargedAna 
//a function for general analysis of charged particle in CDS. It is not allowded that two good CDC tracks share a CDH segment.
//beta and squared-mass is calculated and PID is done by TrackTools::PIDcorr_wide().
//
//input 
//bpctrack, CDSinfo(CDSHitMan & CDSTrackingMan), ConfMan, 
//TLorenzVector of beam (const), T0 hit time (const)
//
//output
//1. cdhseg: vector of CDH segment # which is used for charged particle analysys
//2. pimid,pipid, kmid, protonid : vector of track ID which is IDed to pi-,pii+,K-,proton respectively 
//

int Util::CDSChargedAna(const bool docdcretiming, 
                        LocalTrack* bpctrack,
                        CDSHitMan *cdsman,
                        CDSTrackingMan *trackman,
                        ConfMan *confman,
                        const TLorentzVector LVecbeam,
                        const double ctmt0, 
                        std::vector <int> &cdhseg,
                        std::vector <int> &pimid,
                        std::vector <int> &pipid,
                        std::vector <int> &kmid,
                        std::vector <int> &protonid,
                        const bool MCFlag)
{
  int CDHseg=-1;
  bool chi2OK = true;
  bool CDHseg1hitOK = true;
  bool CDHsharecheckOK = true;
  bool FindMass2OK1 = true;//before CDC retiming 
  bool FindMass2OK2 = true;//after CDC retiming
  bool EnergyLossOK = true;
  
  if(cdhseg.size()!=0 ){
    std::cout << __FILE__ << "input cdhseg size " << cdhseg.size() << std::endl;
  }
  if(pimid.size()!=0 ){
    std::cout << __FILE__ << "input pimid size " << pimid.size() << std::endl;
  }
  if(pipid.size()!=0 ){
    std::cout << __FILE__ << "input pipid size " << pipid.size() << std::endl;
  }
  if(kmid.size()!=0 ){
    std::cout << __FILE__ << "input kmid size " << kmid.size() << std::endl;
  }
  if(protonid.size()!=0 ){
    std::cout << __FILE__ << "input proton size " << protonid.size() << std::endl;
  }

  for( int it=0; it<trackman->nGoodTrack(); it++ ) {
    CDSTrack *track = trackman->Track( trackman->GoodTrackID(it) );

    // chi2 cut can be applied in CDSTrackingMan with MaxChi in CDSFittingParam_posi.param
    Tools::Fill1D( Form("trackchi2_CDC"), track->Chi() );
    if( track->Chi()>cdscuts::cds_chi2_max ) {
      chi2OK = false;
    }
    //asano memo
    //checking if there is CDH hits at the projected position of the track
    if( !track->CDHFlag() ) {
      CDHseg1hitOK = false;
    }
    double mom = track->Momentum();
    TVector3 vtxbline, vtxbhelix ,vtxb; //,vtxb;
    

    //TODO: why using different bpc track origin btw mc and real data ?
    if(MCFlag){
      TVector3 Pos_T0;
      confman->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
      const double zPos_T0 = Pos_T0.Z();
      track->GetVertex( bpctrack->GetPosatZ(zPos_T0), bpctrack->GetMomDir(), vtxbline, vtxbhelix );
    }else{
      track->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtxbline, vtxbhelix );
    }
    track->SetPID(-1);
    vtxb = 0.5*(vtxbline+vtxbhelix);
    Tools::Fill2D(Form("Vtx_ZX"),vtxb.Z(),vtxb.X());
    Tools::Fill2D(Form("Vtx_ZY"),vtxb.Z(),vtxb.Y());
    Tools::Fill2D(Form("Vtx_XY"),vtxb.X(),vtxb.Y());

    double tof = 999.;//cdh-T0
    double mass2 = -999.;
    int nCDHass = track->nCDHHit();
    Tools::Fill1D(Form("mul_CDH_assoc"),nCDHass);
    if(nCDHass>1) {
      CDHseg1hitOK = false;
    }

    for( int icdh=0; icdh<track->nCDHHit(); icdh++ ) {
      HodoscopeLikeHit *cdhhit = track->CDHHit( cdsman, icdh );
      double tmptof = cdhhit->ctmean()-ctmt0;
      if( tmptof<tof || tof==999. ) {
        tof = tmptof;
        CDHseg = cdhhit->seg();
      }
    }

    //asano memo
    //check if the segment has been already used in the analysis.
    bool CDHflag = true;
    for( int icdhseg=0; icdhseg<(int)cdhseg.size(); icdhseg++) {
      if( CDHseg==cdhseg[icdhseg] ) CDHflag = false;
    }
    if( !CDHflag ) {
      //if(Verbosity)std::cerr << " The CDH segments is used by an another track !!! " << std::endl;
      CDHsharecheckOK = false;
      continue;//go to next CDStrack
    }
    cdhseg.push_back( CDHseg );

    // calculation of beta and squared-mass //
    double tmptof=0;
    double beta_calc=0;
    if( !TrackTools::FindMass2( track, bpctrack, tof, LVecbeam.Vect().Mag(),
                                Beam_Kaon, beta_calc, mass2, tmptof ) ) {
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      FindMass2OK1 = false;
      continue;
    }

    // Retiming of CDC track by CDH info. //
    if(docdcretiming) {
      track->Retiming( cdsman, confman, beta_calc, true );
      //asano memo
      //why need this ?
      for( int m=0; m<5; m++ ) {
        track->HelixFitting( cdsman ); // 5 times iteration //
      }
      track->Calc( confman );
    }


    // finalize PID //
    if( !TrackTools::FindMass2( track, bpctrack, tof, LVecbeam.Vect().Mag(),
                                Beam_Kaon, beta_calc, mass2, tmptof ) ) { // not FindMass2C() [20170622] //
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      FindMass2OK2 = false;
      continue;
    }

    const int pid = TrackTools::PIDcorr_wide(mom,mass2);

    track->SetPID( pid );
    Tools::Fill2D( "PID_CDS_beta", 1/beta_calc, mom );
    Tools::Fill2D( "PID_CDS", mass2, mom );


    // Energy loss calculation //
    double tmpl=0;
    TVector3 vtx_beam, vtx_cds;
    if( !track->CalcVertexTimeLength( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), track->Mass(),
                                      vtx_beam, vtx_cds, tmptof, tmpl, true ) ) {
      std::cerr<< __FILE__ << " L. "<<  __LINE__ << 
      " !!! failure in energy loss calculation [CalcVertexTimeLength()] !!! "<<std::endl;
      EnergyLossOK = false;
      continue;
    }

    if(MCFlag &&  
       chi2OK && 
       CDHseg1hitOK && 
       CDHsharecheckOK && 
       FindMass2OK1 && 
       FindMass2OK2 && 
       EnergyLossOK && 
      ((pid==CDS_PiMinus) || (pid==CDS_PiPlus))){
      Util::AnaCDHHitPos(tof,beta_calc,bpctrack,LVecbeam,track,cdsman,confman);
    }

    if( pid==CDS_PiMinus ) {
      pimid.push_back( trackman->GoodTrackID(it) );
    } else if( pid==CDS_PiPlus ) {
      pipid.push_back( trackman->GoodTrackID(it) );
    } else if( pid==CDS_Proton ) {
      protonid.push_back( trackman->GoodTrackID(it) );
    } else if( pid==CDS_Kaon ) {
      kmid.push_back( trackman->GoodTrackID(it) );
    }
  } // for( int it=0; it<trackman->nGoodTrack(); it++ )
  // end of PID (except for neutron) //

  if(!chi2OK) {
    return -7;
  }

  if(!CDHseg1hitOK) {
    return -8;
  }

  if(!CDHsharecheckOK) {
    return -9;
  }

  if(!FindMass2OK1) {
    return -10;
  }

  if(!FindMass2OK2) {
    return -11;
  }

  if(!EnergyLossOK) {
    return -12;
  }

  return pimid.size()+pipid.size()+protonid.size()+kmid.size();
}

double Util::AnalyzeT0(BeamLineHitMan *blman,ConfMan *confman)
{
  //** BHD & T0 **//
  int nBHD = 0;
  for( int i=0; i<blman->nBHD(); i++ ) {
    if( blman->BHD(i)->CheckRange() ) nBHD++;
  }
  int nT0 = 0;
  for( int i=0; i<blman->nT0(); i++ ) {
    if( blman->T0(i)->CheckRange() ) nT0++;
  }
  Tools::Fill1D( Form("mul_BHD"), nBHD );
  Tools::Fill1D( Form("mul_T0"),  nT0 );

  //** T0 = 1hit selection **//
  if( nT0!=1 ) { //!! sada-D p.72 !!//
    return -9999;
  }

  //** Beam PID by T0-BHD TOF **//
  TVector3 vtxT0;
  double ctmt0=0;
  for( int i=0; i<blman->nT0(); i++ ) {
    if( blman->T0(i)->CheckRange() ) {
      ctmt0 = blman->T0(i)->ctmean();
      confman->GetGeomMapManager()->GetGPos(CID_T0, blman->T0(i)->seg(), vtxT0);
    }
  }

  return ctmt0;
}

int Util::BeamPID(EventHeader *header, const double ctmt0,BeamLineHitMan *blman)
{
  //  beamline analysis & event selection
  double ctmBHD=0;
  int PIDBeam = -1; // 0:pi 1:K 3:else
  for( int i=0; i<blman->nBHD(); i++ ) {
    if( blman->BHD(i)->CheckRange() ) {
      ctmBHD = blman->BHD(i)->ctmean();
      double tofBHDT0 = ctmt0-ctmBHD;
      Tools::Fill1D( Form("tof_T0BHD"), tofBHDT0 );
      if( header->kaon() && blcuts::beam_tof_k_min <tofBHDT0 && tofBHDT0<blcuts::beam_tof_k_max  )
        PIDBeam = Beam_Kaon;
      else if( header->pion() && blcuts::beam_tof_pi_min<tofBHDT0 && tofBHDT0<blcuts::beam_tof_pi_max )
        PIDBeam = Beam_Pion;
    }
  }
  Tools::Fill1D(Form("PID_beam"), PIDBeam);
  if( PIDBeam== -1  ) { //** unidentified particle is discarded (other than pi/K) **//
    return -1;
  }

  return PIDBeam;
}


//Beam PID + event selection of BLC1,2, BPC
//TODO: add comment
int Util::EveSelectBeamline(BeamLineTrackMan *bltrackman,
                            CDSTrackingMan *trackman,
                            ConfMan *confman, 
                            int &blc1id, 
                            int &blc2id, 
                            int &bpcid)
{

  unsigned int nblc1 = 0;
  unsigned int nblc2 = 0;
  blc1id = -1;
  blc2id = -1;
  //** timing selection of BLC1/BLC2 **//
  for( int i=0; i<bltrackman->ntrackBLC1(); i++ ) {
    LocalTrack *blc1 = bltrackman->trackBLC1(i);
    Tools::Fill1D( Form("tracktime_BLC1"), blc1->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BLC1"), blc1->chi2all() );
    if( blc1->CheckRange(blcuts::blc1_time_window_min, blcuts::blc1_time_window_max)) {
      nblc1++;
      if( blc1->CheckRange(blcuts::blc1_time_min, blcuts::blc1_time_max) &&
          bltrackman->trackBLC1(i)->chi2all()<blcuts::blc1_chi2_max ) blc1id = i;
    }
  }
  for( int i=0; i<bltrackman->ntrackBLC2(); i++ ) {
    LocalTrack *blc2 = bltrackman->trackBLC2(i);
    Tools::Fill1D( Form("tracktime_BLC2"), blc2->GetTrackTime() );
    Tools::Fill1D( Form("trackchi2_BLC2"), blc2->chi2all() );
    if( blc2->CheckRange(blcuts::blc2_time_window_min, blcuts::blc2_time_window_max) ) {
      nblc2++;
      if( blc2->CheckRange(blcuts::blc2_time_min, blcuts::blc2_time_max) &&
          bltrackman->trackBLC2(i)->chi2all()<blcuts::blc2_chi2_max ) blc2id = i;
    }
  }
  Tools::Fill1D( Form("ntrack_BLC1"), nblc1 );
  Tools::Fill1D( Form("ntrack_BLC2"), nblc2 );

  //** single track selection in each BLC **//
  if( !(nblc1==1 && blc1id!=-1 && nblc2==1 && blc2id!=-1) ) { //** multi-good-beams event is discarded  **//
    return -16;
  }

  //** BPC track selection **//
  int nbpc = 0;
  bpcid = -1;
  double chibpc = 999;
  for( int i=0; i<bltrackman->ntrackBPC(); i++ ) {
    Tools::Fill1D( Form("tracktime_BPC"), bltrackman->trackBPC(i)->GetTrackTime() );
    if( bltrackman->trackBPC(i)->CheckRange(blcuts::bpc_time_window_min,blcuts::bpc_time_window_max) ) {
      nbpc++;
      bpcid = i;
      chibpc = bltrackman->trackBPC(i)->chi2all();
      Tools::Fill1D( Form("trackchi2_BPC"),chibpc);
    }
  }

  Tools::Fill1D( Form("ntrack_BPC"), nbpc );

  if( nbpc!=1 ) {
    return -17;
  }

  LocalTrack *bpctrack = bltrackman->trackBPC(bpcid);
  if( !(bpctrack->CheckRange(blcuts::bpc_time_min,blcuts::bpc_time_max))
      || bpctrack->chi2all()>blcuts::bpc_chi2_max) {
    return -18;
  }

  //** vertex calculation by CDS goodtrack and BPC tracks**/
  for( int it1=0; it1<trackman->nGoodTrack(); it1++ ) {
    trackman->CalcVertex_beam( trackman->GoodTrackID(it1), bltrackman, confman );
  }

  //** BLC2-BPC track matching **//
  bool fblc2bpc = false;
  for( int ii=0; ii<bltrackman->ntrackBLC2(); ii++ ) {
    if( ii!=blc2id ) continue;
    LocalTrack *blc2 = bltrackman->trackBLC2(ii);
    double xblc2bpc[2]= {0,0};
    double yblc2bpc[2]= {0,0};
    double xpos[2]= {0,0};
    double ypos[2]= {0,0};

    TVector3 Pos_BPC, Pos_BLC2, rot;
    confman->GetBLDCWireMapManager()->GetGParam( CID_BPC, Pos_BPC, rot );
    confman->GetBLDCWireMapManager()->GetGParam( CID_BLC2a, Pos_BLC2, rot );
    const double zPos_BPC = Pos_BPC.Z();
    const double zPos_BLC2 = Pos_BLC2.Z();
    const double zPos_BPC_BLC2 = (Pos_BPC.Z()+Pos_BLC2.Z())/2;

    bpctrack->XYPosatZ( zPos_BPC_BLC2, xblc2bpc[0], yblc2bpc[0] );
    bpctrack->XYPosatZ( zPos_BPC, xpos[0], ypos[0] );
    blc2->XYPosatZ( zPos_BPC_BLC2, xblc2bpc[1], yblc2bpc[1]);
    blc2->XYPosatZ( zPos_BLC2, xpos[1], ypos[1]);
    double dxdz[2], dydz[2];
    dxdz[0] = (xpos[0]-xblc2bpc[0]) / (zPos_BPC-zPos_BPC_BLC2);
    dxdz[1] = (xpos[1]-xblc2bpc[1]) / (zPos_BLC2-zPos_BPC_BLC2);
    dydz[0] = (ypos[0]-yblc2bpc[0]) / (zPos_BPC-zPos_BPC_BLC2);
    dydz[1] = (ypos[1]-yblc2bpc[1]) / (zPos_BLC2-zPos_BPC_BLC2);

    if( (xblc2bpc[1]-xblc2bpc[0])< blcuts::blc2bpc_x_min ||
        (xblc2bpc[1]-xblc2bpc[0])> blcuts::blc2bpc_x_max ) fblc2bpc = false;
    else if( (yblc2bpc[1]-yblc2bpc[0])<blcuts::blc2bpc_y_min ||
             (yblc2bpc[1]-yblc2bpc[0])>blcuts::blc2bpc_y_max ) fblc2bpc = false;
    else if( (dxdz[1]-dxdz[0])<blcuts::blc2bpc_dx_min ||
             (dxdz[1]-dxdz[0])>blcuts::blc2bpc_dx_max ) fblc2bpc = false;
    else if( (dydz[1]-dydz[0])<blcuts::blc2bpc_dy_min ||
             (dydz[1]-dydz[0])>blcuts::blc2bpc_dy_max ) fblc2bpc = false;
    else fblc2bpc = true;

    Tools::Fill2D( Form("dydx_BLC2BPC"), xblc2bpc[1]-xblc2bpc[0], yblc2bpc[1]-yblc2bpc[0] );
    Tools::Fill2D( Form("dydzdxdz_BLC2BPC"), dxdz[1]-dxdz[0], dydz[1]-dydz[0] );
  }
  if( !fblc2bpc ) {
    return -19;
  }

  return 0;
}

//usage: use in the loop of CDSTrack
//
//input
//1.meas_tof (T0toCDH)
//2.beta_cal (beta of CDS track) 
//3.bpc track
//4.TLorenzVector of beam
//5.CDS track
//6.ConfMan *confman, 
//
//
//output:
//histograms 
void Util::AnaCDHHitPos(const double meas_tof, const double beta_calc, 
                 LocalTrack *bpc,
                 const TLorentzVector LVecbeam,
                 CDSTrack *track, 
                 CDSHitMan *cdsman,
                 ConfMan *confman)
{
  //projected position of the track at the CDH segment at the surface (r=54.4 cm)
  TVector3 track_pos = track->CDHVertex();
  TVector3 hit_pos;
  if( track->nCDHHit()!=1 ){
    std::cerr<<" #of CDH hit / track is not 1 !!! continue;"<<std::endl;
    return;
  }
  for( int icdh=0; icdh<track->nCDHHit(); icdh++ ){
    HodoscopeLikeHit *cdhhit=track->CDHHit(cdsman,icdh);
    confman->GetGeomMapManager()->GetPos( CID_CDH, cdhhit->seg(), hit_pos );
    hit_pos.SetZ(cdhhit->hitpos());//works in MC. 
  }
  TVector3 diff = track_pos-hit_pos;
  Tools::Fill2D( Form("CDH_mom_diffpos_pi_phi"), (track_pos.Phi()-hit_pos.Phi())/TwoPi*360., track->Momentum() );
  Tools::Fill2D( Form("CDH_mom_diffpos_pi_z"), diff.Z(), track->Momentum() );

  TVector3 Pos_T0;
  confman->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
  double beamtof=0;// T0-VTX
  double momout=0;
  double z_pos = Pos_T0.Z();
  TVector3 vtxb1,vtxb2,vtxb;
  track->GetVertex( bpc->GetPosatZ(z_pos), bpc->GetMomDir(), vtxb1, vtxb2 );
  vtxb = (vtxb1+vtxb2)*0.5;
  ELossTools::CalcElossBeamTGeo( bpc->GetPosatZ(z_pos), vtxb,
				     LVecbeam.Vect().Mag(), kpMass, momout, beamtof );
  double beam_len = (bpc->GetPosatZ(z_pos)-vtxb).Mag();//T0toVTX
  double beam_mom = momout;
  double part_tof = meas_tof-beamtof; // VTX-CDH, scattered particle tof
  double part_len = part_tof*(Const*100.)*beta_calc; //
  double part_mom = track->Momentum(); // from CDC
  double part_beta = sqrt(part_mom*part_mom/(part_mom*part_mom+piMass*piMass));
  double beam_beta = sqrt(beam_mom*beam_mom/(beam_mom*beam_mom+kpMass*kpMass));
  double calc_tof = beam_len/(Const*100.)/beam_beta+part_len/(Const*100.)/part_beta; // T0-VTX part is measured value
	Tools::Fill2D( Form("CDH_mom_TOF_pi"), meas_tof-calc_tof, track->Momentum() );
}
