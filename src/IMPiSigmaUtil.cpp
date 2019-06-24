//author: Hidemitsu Asano
//email : hidemitsu.asano@riken.jp
//
//These are helper functions for EventAnalysisIMPiSigma.cpp and UserSimIMPiSigma.cpp to avoid many duplications of codes.
//All functions are defined in the namespace::Util


#include "IMPiSigmaUtil.h"
#include "Tools.h"
#include "TrackTools.h"


//analyzes # of CDH hits within the TDC range and returns the CDH muliplicity
int Util::GetCDHMul(CDSHitMan *cdsman, const int ntrack, const bool MCFlag)
{
  //** # of CDH-hits cut **//
  int nCDH = 0;
  for( int i=0; i<cdsman->nCDH(); i++ ){
    Tools::Fill2D(Form("CDHtime"),cdsman->CDH(i)->seg(),cdsman->CDH(i)->ctmean());
    Tools::Fill2D(Form("dE_CDHtime"), cdsman->CDH(i)->ctmean(), cdsman->CDH(i)->emean());
    if(ntrack == cdscuts::cds_ngoodtrack){
      Tools::Fill2D(Form("dE_CDHtime_2track"), cdsman->CDH(i)->ctmean(), cdsman->CDH(i)->emean());
    }
    //if( cdsman->CDH(i)->CheckRange() ) nCDH++; //** only requirement of TDC **//
    if(MCFlag){
      if((cdsman->CDH(i)->CheckRange()) && (cdsman->CDH(i)->ctmean()<(cdscuts::tdc_cdh_max+cdscuts::tdc_simoffset))) nCDH++;
    }else{
      if((cdsman->CDH(i)->CheckRange()) && (cdsman->CDH(i)->ctmean()<cdscuts::tdc_cdh_max)){
        //std::cout << cdsman->CDH(i)->seg() << " " <<  cdsman->CDH(i)->ctmean() << std::endl;
        nCDH++;
      }
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
  int nNC  = 0;
  for( int i=0; i<blman->nBVC(); i++ ) {
    if( blman->BVC(i)->CheckRange() ) nBVC++;
  }
  for( int i=0; i<blman->nTOF(); i++ ) {
    if( blman->TOF(i)->CheckRange() ) nCVC++;
  }
  for( int i=0; i<blman->nPC(); i++ ) {
    if( blman->PC(i)->CheckRange() ) nPC++;
  }
  for( int i=0; i<blman->nNC(); i++ ) {
    if( blman->NC(i)->CheckRange() ) nNC++;
  }
  Tools::Fill1D( Form("mul_BVC"), nBVC );
  Tools::Fill1D( Form("mul_CVC"), nCVC );
  Tools::Fill1D( Form("mul_PC"),  nPC );
  Tools::Fill1D( Form("mul_NC"),  nNC );

  if( nBVC || nCVC || nPC ) return true;
  else return false;
}


//inputs
//1. seg, vector of CDH segment IDs to check # of neigboring hits
//2. segall, vector of CDH segment IDs of all CDH hits
//3. pippimhit, vector of CDH segment IDs (size == 2) of pi+ and pi- tracks
//4, cdsman
//returns # of neighboring hits of CDH segments listed in vector cdhseg
//used for Isolation cuts of neutral particle candidates
int Util::GetCDHNeighboringNHits(const std::vector <int> &seg, const std::vector <int> &segall, const std::vector <int> &pippimhit,CDSHitMan *cdsman )
{
  int NNeighboringHits=0;
  for( int ineuseg=0; ineuseg<(int)seg.size(); ineuseg++ ) {
    for( int iseg=0; iseg<(int)segall.size(); iseg++ ) {
      if( seg[ineuseg]-segall[iseg] ) {
        Tools::Fill1D( Form("diff_CDH"), seg[ineuseg]-segall[iseg] );
        int cdhnuhit=-1;
        for(int ihit=0; ihit<cdsman->nCDH(); ihit++) {
          if( cdsman->CDH(ihit)->seg()==seg[ineuseg] ) cdhnuhit = ihit;
        }
        int cdhpippimhit = -1;
        for(int ihit=0; ihit<cdsman->nCDH(); ihit++){
          if( cdsman->CDH(ihit)->seg()==segall[iseg] ) cdhpippimhit = ihit;
        }
        HodoscopeLikeHit *ncdhhit = cdsman->CDH(seg[ineuseg]);
        HodoscopeLikeHit *pippimcdhhit = cdsman->CDH(segall[iseg]);
        double ncdhz = -1.0*ncdhhit->hitpos();
        double pippimcdhz = -1.0*pippimcdhhit->hitpos();
        Tools::Fill2D( Form("diff2D_CDH"), seg[ineuseg]-pippimhit[iseg],ncdhz-pippimcdhz );
      }
      //CDH has 36 segments. # of hits on neghiboring cdh segments
      if( (abs(seg[ineuseg]-segall[iseg])==1) || (abs(seg[ineuseg]-segall[iseg])==35) )
        NNeighboringHits++;
    }
  }
  
  if(pippimhit.size()!=2){
    std::cout << __FILE__ << " L." << __LINE__ << "# of pion tracks !=2 " << std::endl; 
  }
  for( int ineuseg=0; ineuseg<(int)seg.size(); ineuseg++ ) {
    for( int iseg=0; iseg<(int)pippimhit.size(); iseg++ ) {
      if( seg[ineuseg]-pippimhit[iseg] ) {
        Tools::Fill1D( Form("diff_CDH_pippim"), seg[ineuseg]-pippimhit[iseg]);
        int cdhnuhit=-1;
        for(int ihit=0; ihit<cdsman->nCDH(); ihit++) {
          if( cdsman->CDH(ihit)->seg()==seg[ineuseg] ) cdhnuhit = ihit;
        }
        int cdhpippimhit = -1;
        for(int ihit=0; ihit<cdsman->nCDH(); ihit++){
          if( cdsman->CDH(ihit)->seg()==pippimhit[iseg] ) cdhpippimhit = ihit;
        }

        HodoscopeLikeHit *ncdhhit = cdsman->CDH(cdhnuhit);
        HodoscopeLikeHit *pippimcdhhit = cdsman->CDH(cdhpippimhit);
        double ncdhz = -1.0*ncdhhit->hitpos();
        double pippimcdhz = -1.0*pippimcdhhit->hitpos();
        Tools::Fill2D( Form("diff2D_CDH_pippim"), seg[ineuseg]-pippimhit[iseg],ncdhz-pippimcdhz );
      }
    }
  }

  return NNeighboringHits;
}

//returns # of CDH hits which is 2 segment away from 
//CDH segments listed in vector cdhseg
//used for Isolation cuts of neutral particle candidates
int Util::GetCDHTwoSegAwayNHits(const std::vector <int> &seg, const std::vector <int> &segall )
{
  int NTwoSegAwayHits=0;
  for( int ineuseg=0; ineuseg<(int)seg.size(); ineuseg++ ) {
    for( int iseg=0; iseg<(int)segall.size(); iseg++ ) {
      //CDH has 36 segments. 
      if( (abs(seg[ineuseg]-segall[iseg])==2) || (abs(seg[ineuseg]-segall[iseg])==34) )
        NTwoSegAwayHits++;
    }
  }
  //std::cout << "two seg" << NTwoSegAwayHits << std::endl;
  return NTwoSegAwayHits;
}


//returns # of cdc hits of layer 15 and 16 within -+ 15 degree of the CDH hit
int Util::GetNHitsCDCOuter(const TVector3 PosCDH, CDSHitMan *cdsman, const double rangedeg)
{
  int nCDChit = 0;
  static bool state=false;
  if(!state){
    std::cout << __FILE__ << " L." << __LINE__ << " CDC charge veto angle +/- " << rangedeg << std::endl; 
    state = true;
  }
  const double PhiMin = -rangedeg/360.*TwoPi; // rad
  const double PhiMax =  rangedeg/360.*TwoPi; // rad
  for( int ilr=14; ilr<16; ilr++ ) { // charge veto using layer 14, 15
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

void Util::AnaPipPimCDCCDH(const TVector3 PosCDH,const std::vector <int> &seg, const int pip_ID, const int pim_ID, CDSHitMan *cdsman,CDSTrackingMan *trackman)
{
  CDSTrack *track_pip = trackman->Track( pip_ID ); // only 1 track
  CDSTrack *track_pim = trackman->Track( pim_ID ); // only 1 track
  HodoscopeLikeHit *cdhhit_pip = track_pip->CDHHit( cdsman, 0 );
  HodoscopeLikeHit *cdhhit_pim = track_pim->CDHHit( cdsman, 0 );
  int CDHseg_pip = cdhhit_pip->seg();
  int CDHseg_pim = cdhhit_pim->seg();
  Tools::Fill1D( Form("diff_CDH_pip"), seg.at(0)-CDHseg_pip);
  Tools::Fill1D( Form("diff_CDH_pim"), seg.at(0)-CDHseg_pim);

  const int cdchitpipL14 = track_pip->nTrackHit(14);
  const int cdchitpipL15 = track_pip->nTrackHit(15);
  for(int ipip14=0;ipip14<cdchitpipL14;ipip14++){
    CDCHit *cdc = track_pip->hit(cdsman,14,ipip14);
    TVector3 Pos_CDC = cdc->wpos();
    Pos_CDC.SetZ(0);
    const double angle = Pos_CDC.Angle(PosCDH);
    Tools::Fill1D( Form("diff_CDH_CDC_pip"), angle/TwoPi*360. );
  }
  for(int ipip15=0;ipip15<cdchitpipL15;ipip15++){
    CDCHit *cdc = track_pip->hit(cdsman,15,ipip15);
    TVector3 Pos_CDC = cdc->wpos();
    Pos_CDC.SetZ(0);
    const double angle = Pos_CDC.Angle(PosCDH);
    Tools::Fill1D( Form("diff_CDH_CDC_pip"), angle/TwoPi*360. );
  }


  const int cdchitpimL14 = track_pim->nTrackHit(14);
  const int cdchitpimL15 = track_pim->nTrackHit(15);
  for(int ipim14=0;ipim14<cdchitpimL14;ipim14++){
    CDCHit *cdc = track_pim->hit(cdsman,14,ipim14);
    TVector3 Pos_CDC = cdc->wpos();
    Pos_CDC.SetZ(0);
    const double angle = Pos_CDC.Angle(PosCDH);
    Tools::Fill1D( Form("diff_CDH_CDC_pim"), angle/TwoPi*360. );
  }
  for(int ipim15=0;ipim15<cdchitpimL15;ipim15++){
    CDCHit *cdc = track_pim->hit(cdsman,15,ipim15);
    TVector3 Pos_CDC = cdc->wpos();
    Pos_CDC.SetZ(0);
    const double angle = Pos_CDC.Angle(PosCDH);
    Tools::Fill1D( Form("diff_CDH_CDC_pim"), angle/TwoPi*360. );
  }

  return;
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
//a function for general analysis of charged particle in CDS. It is not allowded that two good CDC tracks share one CDH segment.
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
                        BeamLineHitMan *blman,
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
    double mom = track->Momentum();//charge X momentum
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

    double tof = 999.;//TOF of CDH-T0 (slewing corrected)
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
    double correctedtof=0;//CDH-T0 (corrected by energy loss)
    double beta_calc=0;
    if( !TrackTools::FindMass2( track, bpctrack, tof, LVecbeam.Vect().Mag(),
                                Beam_Kaon, beta_calc, mass2, correctedtof ) ) {
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
                                Beam_Kaon, beta_calc, mass2, correctedtof ) ) { // not FindMass2C() [20170622] //
      std::cerr<<" !!! failure in PID_CDS [FindMass2()] !!! "<<std::endl;
      FindMass2OK2 = false;
      continue;
    }

    const int pid = TrackTools::PIDcorr_wide(mom,mass2);

    track->SetPID( pid );
    Tools::Fill2D( "PID_CDS_beta", 1/beta_calc, mom );
    Tools::Fill2D( "PID_CDS", mass2, mom );
    if(pid == CDS_PiMinus) Tools::Fill2D("PID_CDS_PIM",mass2,mom);
    else if(pid == CDS_PiPlus) Tools::Fill2D("PID_CDS_PIP",mass2,mom);
    else if(pid == CDS_Proton) Tools::Fill2D("PID_CDS_Proton",mass2,mom);
    else if(pid == CDS_Kaon) Tools::Fill2D("PID_CDS_Kaon",mass2,mom);

    // Energy loss calculation //
    double tmpl=0;
    TVector3 vtx_beam, vtx_cds;
    if( !track->CalcVertexTimeLength( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), track->Mass(),
                                      vtx_beam, vtx_cds, correctedtof, tmpl, true ) ) {
      std::cerr<< __FILE__ << " L. "<<  __LINE__ << 
      " !!! failure in energy loss calculation [CalcVertexTimeLength()] !!! "<<std::endl;
      EnergyLossOK = false;
      continue;
    }

    if(
       chi2OK && 
       CDHseg1hitOK && 
       CDHsharecheckOK && 
       FindMass2OK1 && 
       FindMass2OK2 && 
       EnergyLossOK && 
      ((pid==CDS_PiMinus) || (pid==CDS_PiPlus))){
        Util::AnaCDHHitPos(tof,beta_calc,bpctrack,LVecbeam,ctmt0,track,cdsman,confman,blman,correctedtof,MCFlag);
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
      if( header->kaon() && (blcuts::beam_tof_k_min <tofBHDT0) && (tofBHDT0<blcuts::beam_tof_k_max)  )
        PIDBeam = Beam_Kaon;
      else if( header->pion() && (blcuts::beam_tof_pi_min<tofBHDT0) && (tofBHDT0<blcuts::beam_tof_pi_max) )
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
                 const double ctmt0,
                 CDSTrack *track, 
                 CDSHitMan *cdsman,
                 ConfMan *confman,
                 BeamLineHitMan *blman,
                 const double correctedtof,
                 const bool MCFlag
                 )
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
    if(MCFlag)hit_pos.SetZ(cdhhit->hitpos());//works in MC. 
    else   hit_pos.SetZ(-1*cdhhit->hitpos());
  }
  TVector3 diff = track_pos-hit_pos;
  TVector3 Pos_T0;
  confman->GetGeomMapManager()->GetPos( CID_T0, 0, Pos_T0 );
  double beamtof=0;// T0-VTX
  double momout=0;
  double z_pos = Pos_T0.Z();
  TVector3 vtxb1,vtxb2,vtxb;
  track->GetVertex( bpc->GetPosatZ(z_pos), bpc->GetMomDir(), vtxb1, vtxb2 );
  vtxb = (vtxb1+vtxb2)*0.5;
  if(!(GeomTools::GetID(vtxb)==CID_Fiducial)) return;
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
  HodoscopeLikeHit *cdh=track->CDHHit(cdsman,0);
	Tools::Fill2D( Form("CDH_mom_TOF_pi") , meas_tof-calc_tof, track->Momentum() );
	Tools::Fill2D( Form("CDH%d_mom_TOF_pi",cdh->seg()) , meas_tof-calc_tof, track->Momentum() );
  
  Tools::Fill2D( Form("CDH_mom_diffpos_pi_phi"), (track_pos.Phi()-hit_pos.Phi())/TwoPi*360., track->Momentum() );
  Tools::Fill2D( Form("CDH_mom_diffpos_pi_z"), diff.Z(), track->Momentum() );
  Tools::H2( Form("CDH_diffpos_z_pi_z"), track_pos.Z(),diff.Z(),1000,-50,50,1000,-50,50);
  
  //return here as checking histogram is not necessary
  return;
  for( int icdh=0; icdh<track->nCDHHit(); icdh++ ){
    HodoscopeLikeHit *cdhhit=track->CDHHit(cdsman,icdh);
    HodoscopeLikeHit *t0hit = blman->T0(0);
    int t0seg = t0hit->seg();
    //Tools::H2( Form("CDH_diffpos_z_pi_z_seg%d",cdhhit->seg()), track_pos.Z(),diff.Z(),1000,-50,50,1000,-50,50);
    //Tools::H2( Form("CDH_diffpos_ctsub_pi_z_seg%d",cdhhit->seg()),cdhhit->ctsub(),diff.Z(),1000,-4,4,1000,-50,50);
    if(-5 < track_pos.Z() && track_pos.Z() <5 ){
      //Tools::H1(Form("CDH_diffpos_pi_z_seg%d",cdhhit->seg()),diff.Z(),1000,-50,50);
      //Tools::H1( Form("CDH_diffpos_pi_z_seg%d",cdhhit->seg()),diff.Z(),400,-20,20);
      //Tools::H1( Form("CTMean%d",cdhhit->seg()), cdhhit->ctmean(), 2000,-50,150 );
      //Tools::H1( Form("CTSub%d",cdhhit->seg()), cdhhit->ctsub(), 1000,-50,50 );
    }
    if(!MCFlag && (fabs(track->Momentum())>0.1) && t0seg==3// &&
    //if(!MCFlag  && t0seg==3 &&
      // && (-5 < track_pos.Z()) && (track_pos.Z() <5)
      ){
      double tof = (cdhhit->tmean()) - ctmt0 - correctedtof;
      double ctof = (cdhhit->ctmean()) - ctmt0 - correctedtof;
      if(track->charge()>0){
        Tools::H2( Form("dECDH_T0%d_CDHU%d_PIP",t0seg,cdhhit->seg()), cdhhit->eu(),tof,500, -0.5,49.5,300,0,15.);
        Tools::H2( Form("dECDH_T0%d_CDHD%d_PIP",t0seg,cdhhit->seg()), cdhhit->ed(),tof,500, -0.5,49.5,300,0,15.);
        Tools::H2( Form("dECDH_T0%d_CDHU%d_PIPc",t0seg,cdhhit->seg()), cdhhit->eu(),ctof,500, -0.5,49.5,300,0,15.);
        Tools::H2( Form("dECDH_T0%d_CDHD%d_PIPc",t0seg,cdhhit->seg()), cdhhit->ed(),ctof,500, -0.5,49.5,300,0,15.);
      }else{
        Tools::H2( Form("dECDH_T0%d_CDHU%d_PIM",t0seg,cdhhit->seg()), cdhhit->eu(),tof,500, -0.5,49.5,300,0,15.);
        Tools::H2( Form("dECDH_T0%d_CDHD%d_PIM",t0seg,cdhhit->seg()), cdhhit->ed(),tof,500, -0.5,49.5,300,0,15.);
        Tools::H2( Form("dECDH_T0%d_CDHU%d_PIMc",t0seg,cdhhit->seg()), cdhhit->eu(),ctof,500, -0.5,49.5,300,0,15.);
        Tools::H2( Form("dECDH_T0%d_CDHD%d_PIMc",t0seg,cdhhit->seg()), cdhhit->ed(),ctof,500, -0.5,49.5,300,0,15.);
      }
      //std::cout << "tmean " << cdhhit->tmean() << " ctmean " << cdhhit->ctmean() << std::endl;
      Tools::H2( Form("dECDH_T0%d_CDH%d_PI",t0seg,cdhhit->seg()), cdhhit->emean(),tof,500, -0.5,49.5,600,0,30.);
      Tools::H2( Form("dECDH_T0%d_CDH%d_PIc",t0seg,cdhhit->seg()), cdhhit->emean(),ctof,500, -0.5,49.5,600,0,30.);
      //std::cout << "t0seg " << t0hit->ctu() << "  " << t0hit->ctd() << std::endl;
      //Tools::H2(Form("ectT0U%d",t0seg),t0hit->eu(),t0hit->ctu(),200,-0.5,4.5,300,5,20);
      //Tools::H2(Form("ectT0D%d",t0seg),t0hit->ed(),t0hit->ctd(),200,-0.5,4.5,300,5,20);
    }//!MCFlag
  }

}

void Util::CorrectCDHz(CDSHitMan *cdsman){
   for(int i=0;i<cdsman->nCDH();i++){
     double HitPosition = cdsman->CDH(i)->hitpos();
     int CDHSeg = cdsman->CDH(i)->seg();
    {
      //CDCz <10
      if(CDHSeg==1) HitPosition -= 16.639714;
      else if(CDHSeg==2) HitPosition -= 14.798649;
      else if(CDHSeg==3) HitPosition -= 5.953724;
      else if(CDHSeg==4) HitPosition -= 12.213804;
      else if(CDHSeg==5) HitPosition -= 8.428154;
      else if(CDHSeg==6) HitPosition -= 15.883103;
      else if(CDHSeg==7) HitPosition -= 12.115939;
      else if(CDHSeg==8) HitPosition -= 16.009665;
      else if(CDHSeg==9) HitPosition -= 14.164350;
      else if(CDHSeg==10) HitPosition -= 13.751159;
      else if(CDHSeg==11) HitPosition -= 16.320796;
      else if(CDHSeg==12) HitPosition -= -17.983675;
      else if(CDHSeg==13) HitPosition -= 0.183959;
      else if(CDHSeg==14) HitPosition -= 6.421490;
      else if(CDHSeg==15) HitPosition -= 10.357269;
      else if(CDHSeg==16) HitPosition -= 4.211167;
      else if(CDHSeg==17) HitPosition -= 12.573929;
      else if(CDHSeg==18) HitPosition -= 13.799147;
      else if(CDHSeg==19) HitPosition -= 8.857060;
      else if(CDHSeg==20) HitPosition -= -11.439437;
      else if(CDHSeg==21) HitPosition -= 14.573345;
      else if(CDHSeg==22) HitPosition -= 15.119745;
      else if(CDHSeg==23) HitPosition -= 9.142359;
      else if(CDHSeg==24) HitPosition -= 12.666169;
      else if(CDHSeg==25) HitPosition -= 13.793150;
      else if(CDHSeg==26) HitPosition -= 0.132931;
      else if(CDHSeg==27) HitPosition -= 11.586085;
      else if(CDHSeg==28) HitPosition -= 25.751971;
      else if(CDHSeg==29) HitPosition -= 8.215086;
      else if(CDHSeg==30) HitPosition -= 7.648705;
      else if(CDHSeg==31) HitPosition -= 0.085094;
      else if(CDHSeg==32) HitPosition -= -0.417300;
      else if(CDHSeg==33) HitPosition -= -5.000405;
      else if(CDHSeg==34) HitPosition -= -11.931577;
      else if(CDHSeg==35) HitPosition -= -7.746246;
      else if(CDHSeg==36) HitPosition -= -9.454847;
      //CDCz <1 cm correction
      if(CDHSeg==1) HitPosition -= 0.044796;
      else if(CDHSeg==2) HitPosition -= 0.050329;
      else if(CDHSeg==3) HitPosition -= 0.047539;
      else if(CDHSeg==4) HitPosition -= 0.011800;
      else if(CDHSeg==5) HitPosition -= 0.025058;
      else if(CDHSeg==6) HitPosition -= 0.066763;
      else if(CDHSeg==7) HitPosition -= 0.048025;
      else if(CDHSeg==8) HitPosition -= 0.053559;
      else if(CDHSeg==9) HitPosition -= 0.041002;
      else if(CDHSeg==10) HitPosition -= 0.043818;
      else if(CDHSeg==11) HitPosition -= 0.016847;
      else if(CDHSeg==12) HitPosition -= 0.060654;
      else if(CDHSeg==13) HitPosition -= 0.021380;
      else if(CDHSeg==14) HitPosition -= 0.054790;
      else if(CDHSeg==15) HitPosition -= 0.027565;
      else if(CDHSeg==16) HitPosition -= 0.012008;
      else if(CDHSeg==17) HitPosition -= 0.047283;
      else if(CDHSeg==18) HitPosition -= 0.038955;
      else if(CDHSeg==19) HitPosition -= 0.028200;
      else if(CDHSeg==20) HitPosition -= 0.033359;
      else if(CDHSeg==21) HitPosition -= 0.055034;
      else if(CDHSeg==22) HitPosition -= -0.006987;
      else if(CDHSeg==23) HitPosition -= 0.024364;
      else if(CDHSeg==24) HitPosition -= 0.013297;
      else if(CDHSeg==25) HitPosition -= 0.014870;
      else if(CDHSeg==26) HitPosition -= 0.050980;
      else if(CDHSeg==27) HitPosition -= 0.020656;
      else if(CDHSeg==28) HitPosition -= 0.047146;
      else if(CDHSeg==29) HitPosition -= 0.064573;
      else if(CDHSeg==30) HitPosition -= 0.020890;
      else if(CDHSeg==31) HitPosition -= 0.026834;
      else if(CDHSeg==32) HitPosition -= 0.022567;
      else if(CDHSeg==33) HitPosition -= 0.064572;
      else if(CDHSeg==34) HitPosition -= -0.070272;
      else if(CDHSeg==35) HitPosition -= 0.015730;
      else if(CDHSeg==36) HitPosition -= 0.056956;
    }

     cdsman->CDH(i)->SetHitPosition(HitPosition);

   }
}



void Util::AnaReactionData( ReactionData *reactionData){
  int ndecay    = reactionData->NParticle(0);
  int nspec     = reactionData->NParticle(1);
  int nparticle = ndecay+nspec;
  if(nparticle!=3) std::cout << "nparticle is 3 ! "  << nparticle << std::endl;
  for(int ipart=0;ipart<nparticle;ipart++){
    //0: neutron, or proton
    //1: Sigma,  or Lambda
    //2: pi
    double coscm = reactionData->GetCMParticle(ipart).CosTheta();
    Tools::H1(Form("ReactCosCM_%d",ipart),coscm, 100,-1.0,1.0);
    double cosl = reactionData->GetParticle(ipart).CosTheta();
    Tools::H1(Form("ReactCos_%d",ipart),cosl, 100,-1.0,1.0);
    double momcm = (reactionData->GetCMParticle(ipart).P())/1000.;
    Tools::H1(Form("Reactmomcm_%d",ipart),momcm, 150,0,1.5);
    double moml = (reactionData->GetParticle(ipart).P())/1000.;
    Tools::H1(Form("Reactmom_%d",ipart),moml, 150,0,1.5);
  }
  TVector3 beammom(0,0,1000);//1 MeV/c 
  TLorentzVector TL_beam;
  TL_beam.SetVectM(beammom, 493.677);
  TLorentzVector TL_nmiss = reactionData->GetParticle(0);
  double q = (TL_beam.Vect()-TL_nmiss.Vect()).Mag()/1000.;
  TLorentzVector TL_Sigma = reactionData->GetParticle(1);
  TLorentzVector TL_piSigma = TL_Sigma + reactionData->GetParticle(2);
  double mass = TL_piSigma.M()/1000.;
  Tools::H2(Form("React_q_IMPiSigma"),mass,q,500,1,2,300,0,1.5);
}

