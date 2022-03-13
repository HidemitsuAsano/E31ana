//author: Hidemitsu Asano
//email : hidemitsu.asano@riken.jp
//
//These are helper functions for EventAnalysisIMPiSigma.cpp and UserSimIMPiSigma.cpp to avoid many duplications of codes.
//All functions are defined in the namespace::Util


#include "IMPiSigmaUtil.h"
#include "Tools.h"
#include "TrackTools.h"


//analyzes # of CDH hits within the TDC range and returns the CDH muliplicity
int Util::GetCDHMul(CDSHitMan *cdsman, const int ntrack,const bool CDH3trgfired,  const bool MCFlag)
{
  //** # of CDH-hits cut **//
  int nCDH = 0;
  double firsthittime=0;
  double lasthittime=0;
  std::vector <double> vtime;
  for( int i=0; i<cdsman->nCDH(); i++ ){
    Tools::Fill2D(Form("CDHtime"),cdsman->CDH(i)->seg(),cdsman->CDH(i)->ctmean());
    Tools::Fill2D(Form("CDHdE"),cdsman->CDH(i)->seg(),cdsman->CDH(i)->emean());
    Tools::Fill2D(Form("dE_CDHtime"), cdsman->CDH(i)->ctmean(), cdsman->CDH(i)->emean());
    if(ntrack == cdscuts::cds_ngoodtrack){
      Tools::Fill2D(Form("dE_CDHtime_2track"), cdsman->CDH(i)->ctmean(), cdsman->CDH(i)->emean());
    }
    if(MCFlag){
      if((cdsman->CDH(i)->CheckRange()) && (cdsman->CDH(i)->ctmean()<(cdscuts::tdc_cdh_max+cdscuts::tdc_simoffset))){
        int seg = cdsman->CDH(i)->seg();
        double emean = cdsman->CDH(i)->emean();
        Tools::Fill2D(Form("CDHdE_wt"),seg,emean);
        vtime.push_back(cdsman->CDH(i)->ctmean());
        nCDH++;
      }
    }else{
      if((cdsman->CDH(i)->CheckRange()) && (cdsman->CDH(i)->ctmean()<cdscuts::tdc_cdh_max)){
        Tools::Fill2D(Form("CDHdE_wt"),cdsman->CDH(i)->seg(),cdsman->CDH(i)->emean());
        //std::cout << cdsman->CDH(i)->seg() << " " <<  cdsman->CDH(i)->ctmean() << std::endl;
        vtime.push_back(cdsman->CDH(i)->ctmean());
        nCDH++;
      }
    }
  }
  Tools::Fill1D( Form("mul_CDH"), nCDH );
  
  if(vtime.size()>=2){
    std::sort(vtime.begin(),vtime.end());
    Tools::H2( Form("lasttime_mul_CDH2"), nCDH,vtime.at(1),11, -0.5, 10.5,1000,0,100);
  }
    
  if(vtime.size()>=3){
    std::sort(vtime.begin(),vtime.end());
    Tools::H2( Form("lasttime_mul_CDH"), nCDH,vtime.at(2),11, -0.5, 10.5,1000,0,100);
    if(CDH3trgfired){
      Tools::H2( Form("lasttime_mul_CDH_CDH3trg"), nCDH,vtime.at(2),11, -0.5, 10.5,1000,0,100);
    }
    if(nCDH==3){
      Tools::H1( Form("firstlasttimediff"), vtime.back()-vtime.front(),1000,0,100);
      if(CDH3trgfired){
        Tools::H1( Form("firstlasttimediff_CDH3trg"), vtime.back()-vtime.front(),1000,0,100);
      }
    }
  }
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
int Util::GetCDHNeighboringNHits(const std::vector <int> &seg, const std::vector <int> &segall, 
                                 const std::vector <int> &pippimhit,CDSHitMan *cdsman, bool MCFlag )
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
        HodoscopeLikeHit *ncdhhit = cdsman->CDH(cdhnuhit);
        HodoscopeLikeHit *pippimcdhhit = cdsman->CDH(cdhpippimhit);
        double ncdhz = -1.0*ncdhhit->hitpos();
        double pippimcdhz = -1.0*pippimcdhhit->hitpos();
        if(MCFlag){
          ncdhz *= -1.0;
          pippimcdhz *= -1.0;
        }
        Tools::Fill2D( Form("diff2D_CDH"), seg[ineuseg]-pippimhit[iseg],ncdhz-pippimcdhz );
      }
      //CDH has 36 segments. # of hits on neghiboring cdh segments
      if( (abs(seg[ineuseg]-segall[iseg])==1) || (abs(seg[ineuseg]-segall[iseg])==35) )
        NNeighboringHits++;
    }
  }
  
  if(pippimhit.size()!=2){
    //std::cout << __FILE__ << " L." << __LINE__ << "# of pion tracks !=2 " << std::endl; 
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
        if(MCFlag){
          ncdhz *= -1.0;
          pippimcdhz *= -1.0;
        }
        Tools::Fill2D( Form("diff2D_CDH_pippim"), seg[ineuseg]-pippimhit[iseg],ncdhz-pippimcdhz );
      }
    }
  }

  return NNeighboringHits;
}

//returns # of CDH hits which is 2 segment away from 
//CDH segments listed in vector cdhseg
//used for Isolation cuts of neutral particle candidates
int Util::GetCDHTwoSegAwayNHits(const std::vector <int> &seg, const std::vector <int> &segall, bool MCFlag )
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

int Util::GetNHitsCDCInner3Lay(CDSHitMan *cdsman)
{
  int nCDChit = 0;

  for( int ilr=0; ilr<3; ilr++ ) { // 
    for( int icdchit=0; icdchit<cdsman->nCDC(ilr); icdchit++ ) {
      CDCHit *cdc=cdsman->CDC(ilr,icdchit);
      nCDChit++;
    }//icdchit
  }//ilr
  
  return nCDChit;
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

//returns # of cdc hits of layer 15 and 16 within -+ 15 degree of the CDH hit,
//which are not associated with CDS tracks
int Util::GetNHitsCDCOuterNoAss(const TVector3 PosCDH, CDSHitMan *cdsman,CDSTrackingMan *trackman, const double rangedeg)
{
  int nCDChit = 0;
  static bool state=false;
  if(!state){
    std::cout << __FILE__ << " L." << __LINE__ << " CDC charge veto (No. Association) angle +/- " << rangedeg << std::endl; 
    state = true;
  }
  const double PhiMin = -rangedeg/360.*TwoPi; // rad
  const double PhiMax =  rangedeg/360.*TwoPi; // rad
  for( int ilr=14; ilr<16; ilr++ ) { // charge veto using layer 14, 15
    for( int icdchit=0; icdchit<cdsman->nCDC(ilr); icdchit++ ) {
      CDCHit *cdc=cdsman->CDC(ilr,icdchit);
      
      bool AssFlag = false;
      for( int it=0; it<trackman->nGoodTrack(); it++ ) {
        CDSTrack *track = trackman->Track( trackman->GoodTrackID(it) );
        int ntrackhit = track->nTrackHit(ilr);
        for(int ihit=0;ihit<ntrackhit;ihit++){
          CDCHit *cdcass =  track->hit(cdsman,ilr,ihit);
          if(cdcass == cdc) AssFlag = true;
          if(AssFlag){
            //std::cout << "cdc       angle " << cdc->wpos().Phi()  << std::endl;
            //std::cout << "cdchitass angle " << cdcass->wpos().Phi()  << std::endl;
          }
        }
      }

      /*
      int s7nCluster = trackman->nCluster(7);
      std::cout << "s7nCluster " << s7nCluster << std::endl;
      for(int iclss7=0;iclss7<s7nCluster;iclss7++){
        HitCluster *s7Cluster = trackman->Cluster(7,iclss7);
        int nhit = s7Cluster->nHit();
        std::cout << "nhit " << nhit << std::endl;
        for(int ihit=0;ihit<nhit;ihit++){
          CDCHit *cdchitass = s7Cluster->CDC(cdsman,ihit);
          if(cdc == cdchitass) AssFlag = true;
          std::cout << "cdc       angle " << cdc->wpos().Phi()  << std::endl;
          std::cout << "cdchitass angle " << cdchitass->wpos().Phi()  << std::endl;
        }
      }*/
     
      TVector3 Pos_CDC = cdc->wpos();
      Pos_CDC.SetZ(0); // only xy pos is used
      const double angle = Pos_CDC.Angle(PosCDH); // rad
      Tools::Fill1D( Form("diff_CDH_CDC"), angle/TwoPi*360 );
      if( !AssFlag && (PhiMin<angle && angle<PhiMax) ){ 
        nCDChit++;
        //std::cerr<<"CDC "<<ilr<<" "<<icdchit<<" "<<cdc->wire()<<" -> "<<Pos_CDC.Phi()/TwoPi*360.
        //  <<" deg :: diff = "<<angle/TwoPi*360<<" deg"<<std::endl;
        //std::cerr << "AssFlag " << AssFlag << std::endl;
      }
    }//icdchit
  }//ilr
  
  return nCDChit;
}


void Util::AnaPipPimCDCCDH(const TVector3 PosCDH,const HodoscopeLikeHit *nhit, const int pip_ID, const int pim_ID, CDSHitMan *cdsman,CDSTrackingMan *trackman)
{
  CDSTrack *track_pip = trackman->Track( pip_ID ); // only 1 track
  CDSTrack *track_pim = trackman->Track( pim_ID ); // only 1 track
  HodoscopeLikeHit *cdhhit_pip = track_pip->CDHHit( cdsman, 0 );
  HodoscopeLikeHit *cdhhit_pim = track_pim->CDHHit( cdsman, 0 );
  int CDHseg_pip = cdhhit_pip->seg();
  int CDHseg_pim = cdhhit_pim->seg();
  Tools::Fill1D( Form("diff_CDH_pip"), nhit->seg()-CDHseg_pip);
  Tools::Fill1D( Form("diff_CDH_pim"), nhit->seg()-CDHseg_pim);

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
  Tools::Fill2D( Form("D5chi2_mom"),beammom,bchi);
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
//3. pim_projected, pip_projected: CDC track projected position on the CDH segment

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
                        TVector3 &pim_projected,
                        TVector3 &pip_projected,
                        const bool MCFlag,
                        const int H2data)
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
    
    double dE = -999.;
    for( int icdh=0; icdh<track->nCDHHit(); icdh++ ) {
      HodoscopeLikeHit *cdhhit = track->CDHHit( cdsman, icdh );
      double tmptof = cdhhit->ctmean()-ctmt0;
      if( tmptof<tof || tof==999. ) {
        tof = tmptof;
        CDHseg = cdhhit->seg();
        dE = cdhhit->emean();
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

    int pid = -1;
    if(H2data==0){
      pid = TrackTools::PIDcorr_wide(mom,mass2);
    }else if(H2data==1){
      //yamaga param ?
      pid = TrackTools::PIDcorr3(mom,mass2);
    }

    track->SetPID( pid );
    Tools::Fill2D( "PID_CDS_beta", 1./beta_calc, mom );
    Tools::Fill2D( "PID_CDS", mass2, mom );
    if(pid == CDS_PiMinus){
      Tools::Fill2D("PID_CDS_PIM_beta",1./beta_calc,mom);
      Tools::Fill2D("PID_CDS_PIM",mass2,mom);
      Tools::H2("PIM_TOF_MOM",fabs(mom),correctedtof,100,0,2,10000,0,100);
      Tools::H2("dE_mom_pi",mom,dE,100,0,2,500,0,50);
    }else if(pid == CDS_PiPlus){
      Tools::Fill2D("PID_CDS_PIP_beta",1./beta_calc,mom);
      Tools::Fill2D("PID_CDS_PIP",mass2,mom);
      Tools::H2("PIP_TOF_MOM",mom,correctedtof,100,0,2,10000,0,100);
      Tools::H2("dE_mom_pi",mom,dE,100,0,2,500,0,50);
    }
    else if(pid == CDS_Proton){
      Tools::Fill2D("PID_CDS_Proton",mass2,mom);
      Tools::H2("dE_mom_proton",mom,dE,100,0,2,500,0,50);
    }else if(pid == CDS_Kaon) Tools::Fill2D("PID_CDS_Kaon",mass2,mom);

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
      ( (pid==CDS_PiMinus) || (pid==CDS_PiPlus) || (pid==CDS_Proton)) 
      ){
        //Util::AnaCDHHitPos(tof,beta_calc,bpctrack,LVecbeam,ctmt0,track,cdsman,confman,blman,correctedtof,pid,MCFlag);
    }

    TVector3 cdh_projected = track->CDHVertex();
    if(pid==CDS_PiMinus) pim_projected = cdh_projected;
    else if(pid==CDS_PiPlus) pip_projected = cdh_projected;
    
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

double Util::AnalyzeT0(BeamLineHitMan *blman,ConfMan *confman,int &t0seg)
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
      t0seg = blman->T0(i)->seg();
    }
  }
  Tools::Fill1D("T0time",ctmt0);

  return ctmt0;
}


int Util::BeamPID(EventHeader *header, const double ctmt0,BeamLineHitMan *blman,const int H2data)
{
  //  beamline analysis & event selection
  double ctmBHD=0;
  int PIDBeam = -1; // 0:pi 1:K 3:else
  static int state = 0;
  if(H2data==0){
    if(!state){
      std::cout << __FILE__ << " L." <<  __LINE__ << " d2 data " << std::endl;
      std::cout << "kaon beam tof range " << blcuts::beam_tof_k_min << " - " << blcuts::beam_tof_k_max << std::endl;
      state = 1;
    }
  }else if(H2data==1){
    if(!state){
      std::cout << __FILE__ << " L." <<  __LINE__ << " H2 data " << std::endl;
      std::cout << "kaon beam tof range " << blcuts::beam_tof_k_min_h2 << " - " << blcuts::beam_tof_k_max_h2 << std::endl;
      state = 1;
    }
  }
  for( int i=0; i<blman->nBHD(); i++ ) {
    if( blman->BHD(i)->CheckRange() ) {
      ctmBHD = blman->BHD(i)->ctmean();
      double tofBHDT0 = ctmt0-ctmBHD;
      Tools::Fill1D( Form("tof_T0BHD"), tofBHDT0 );
      if(!H2data){
        if( header->kaon() && (blcuts::beam_tof_k_min <tofBHDT0) && (tofBHDT0<blcuts::beam_tof_k_max)  )
          PIDBeam = Beam_Kaon;
        else if( header->pion() && (blcuts::beam_tof_pi_min<tofBHDT0) && (tofBHDT0<blcuts::beam_tof_pi_max) )
          PIDBeam = Beam_Pion;
      }else if(H2data){
        if( header->kaon() && (blcuts::beam_tof_k_min_h2 <tofBHDT0) && (tofBHDT0<blcuts::beam_tof_k_max_h2)  )
          PIDBeam = Beam_Kaon;
        else if( header->pion() && (blcuts::beam_tof_pi_min_h2<tofBHDT0) && (tofBHDT0<blcuts::beam_tof_pi_max_h2) )
          PIDBeam = Beam_Pion;
      }
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
    const double zPos_BPC_BLC2 = -63.15; //(Pos_BPC.Z()+Pos_BLC2.Z())/2;
    
    //std::cout << "zPos_BPC " << zPos_BPC << std::endl;
    //std::cout << "zPos_BLC2 " << zPos_BLC2 << std::endl;
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
                 const int pid,
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
  double Mass =0;
  if(pid==CDS_PiMinus || pid==CDS_PiPlus) Mass = piMass;
  if(pid==CDS_Proton) Mass = pMass;
  double part_beta = sqrt(part_mom*part_mom/(part_mom*part_mom+Mass*Mass));
  double beam_beta = sqrt(beam_mom*beam_mom/(beam_mom*beam_mom+kpMass*kpMass));
  double calc_tof = beam_len/(Const*100.)/beam_beta+part_len/(Const*100.)/part_beta; // T0-VTX part is measured value
  HodoscopeLikeHit *cdh=track->CDHHit(cdsman,0);
	if( (pid==CDS_PiMinus)|| (pid==CDS_PiPlus)){
    Tools::Fill2D( Form("CDH_mom_TOF_pi") , meas_tof-calc_tof, track->Momentum() );
    Tools::Fill2D( Form("CDH%d_mom_TOF_pi",cdh->seg()) , meas_tof-calc_tof, track->Momentum() );

    Tools::Fill2D( Form("CDH_mom_diffpos_pi_phi"), (track_pos.Phi()-hit_pos.Phi())/TwoPi*360., track->Momentum() );
    Tools::Fill2D( Form("CDH_mom_diffpos_pi_z"), diff.Z(), track->Momentum() );
    Tools::H2( Form("CDH_diffpos_z_pi_z"), track_pos.Z(),diff.Z(),1000,-50,50,1000,-50,50);
  }else if(pid==CDS_Proton){
    Tools::Fill2D( Form("CDH_mom_TOF_p") , meas_tof-calc_tof, track->Momentum() );
    Tools::Fill2D( Form("CDH%d_mom_TOF_p",cdh->seg()) , meas_tof-calc_tof, track->Momentum() );
  
    Tools::Fill2D( Form("CDH_mom_diffpos_p_phi"), (track_pos.Phi()-hit_pos.Phi())/TwoPi*360., track->Momentum() );
    Tools::Fill2D( Form("CDH_mom_diffpos_p_z"), diff.Z(), track->Momentum() );
    Tools::Fill2D( Form("CDH_diffpos_z_p_z"), track_pos.Z(),diff.Z());
  }
  
  for( int icdh=0; icdh<track->nCDHHit(); icdh++ ){
    HodoscopeLikeHit *cdhhit=track->CDHHit(cdsman,icdh);
    HodoscopeLikeHit *t0hit = blman->T0(0);
    int t0seg = t0hit->seg();
    if( (pid==CDS_PiMinus)|| (pid==CDS_PiPlus)){
      Tools::H2( Form("CDH_diffpos_z_pi_z_seg%d",cdhhit->seg()), track_pos.Z(),diff.Z(),1000,-50,50,1000,-50,50);
      Tools::H2( Form("CDH_diffpos_ctsub_pi_z_seg%d",cdhhit->seg()),cdhhit->ctsub(),diff.Z(),1000,-4,4,1000,-50,50);
    }else if(pid==CDS_Proton){
      Tools::H2( Form("CDH_diffpos_z_p_z_seg%d",cdhhit->seg()), track_pos.Z(),diff.Z(),1000,-50,50,1000,-50,50);
      Tools::H2( Form("CDH_diffpos_ctsub_p_z_seg%d",cdhhit->seg()),cdhhit->ctsub(),diff.Z(),1000,-4,4,1000,-50,50);
    }
    if(-5 < track_pos.Z() && track_pos.Z() <5 ){
      if( (pid==CDS_PiMinus)|| (pid==CDS_PiPlus)){
        Tools::H1(Form("CDH_diffpos_pi_z_seg%d",cdhhit->seg()),diff.Z(),1000,-50,50);
        Tools::H1( Form("CDH_diffpos_pi_z_seg%d",cdhhit->seg()),diff.Z(),400,-20,20);
        Tools::H1( Form("CTMean%d",cdhhit->seg()), cdhhit->ctmean(), 2000,-50,150 );
        Tools::H1( Form("CTSub%d",cdhhit->seg()), cdhhit->ctsub(), 1000,-50,50 );
      }else if(pid==CDS_Proton){
        Tools::H1(Form("CDH_diffpos_p_z_seg%d",cdhhit->seg()),diff.Z(),1000,-50,50);
        Tools::H1( Form("CDH_diffpos_p_z_seg%d",cdhhit->seg()),diff.Z(),400,-20,20);
        Tools::H1( Form("CTMean%d_p",cdhhit->seg()), cdhhit->ctmean(), 2000,-50,150 );
        Tools::H1( Form("CTSub%d_p",cdhhit->seg()), cdhhit->ctsub(), 1000,-50,50 );
      }
    }
    /*
    if( (fabs(track->Momentum())>0.1) && t0seg==3// &&
    if(!MCFlag  && t0seg==3 &&
       && (-5 < track_pos.Z()) && (track_pos.Z() <5)
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
    */ 
  }
}

HodoscopeLikeHit* Util::CDHsegToCDHhitIndex(const int cdhsegid, CDSHitMan *cdsman){
  int icdh = -1;
  for( int icdhhit=0; icdhhit<cdsman->nCDH(); icdhhit++ ){
    if( cdsman->CDH(icdhhit)->seg()==cdhsegid) icdh = icdhhit;
  }
  HodoscopeLikeHit *ncdhhitcan = cdsman->CDH(icdh);
  
  return ncdhhitcan;
}

void Util::AnaReactionData( ReactionData *reactionData){
  int ndecay    = reactionData->NParticle(0);
  int nspec     = reactionData->NParticle(1);
  int nparticle = ndecay+nspec;
  static bool isState=true;
  const int reactionID = reactionData->ReactionID();
  if(isState){
    std::cout << __func__ << " reactionID  "  << reactionID << std::endl;
    std::cout << __func__ << " nparticle is  "  << nparticle << std::endl;
    isState = false;
  }
  for(int ipart=0;ipart<nparticle;ipart++){
    //PiSigma (nmiss) | Lpim (pmiss) | npipi (Lmiss)
    //0: neutron      | proton       | neutron
    //1: Sigma        | Lambda       | pi+
    //2: pi+ or pi-   | pi-          | pi-
    //3:              |              | Lambda
    //
    double coscm = reactionData->GetCMParticle(ipart).CosTheta();
    Tools::H1(Form("ReactCosCM_%d",ipart),coscm, 100,-1.0,1.0);
    double cosl = reactionData->GetParticle(ipart).CosTheta();
    Tools::H1(Form("ReactCos_%d",ipart),cosl, 100,-1.0,1.0);
    double momcm = (reactionData->GetCMParticle(ipart).P())/1000.;
    Tools::H1(Form("Reactmomcm_%d",ipart),momcm, 150,0,1.5);
    double moml = (reactionData->GetParticle(ipart).P())/1000.;
    Tools::H1(Form("Reactmom_%d",ipart),moml, 150,0,1.5);
  }

  //TVector3 beammom(0,0,1000);//1 MeV/c 
  TLorentzVector TL_beam = reactionData->GetInitParticle(0);
  //TL_beam.SetVectM(beammom, 493.677);
  if(reactionID == gen::reactionID_Spmode || gen::reactionID_Smmode){
    TLorentzVector TL_nmiss = reactionData->GetParticle(0);
    double q = (TL_beam.Vect()-TL_nmiss.Vect()).Mag()/1000.;
    TLorentzVector TL_Sigma = reactionData->GetParticle(1);
    TLorentzVector TL_piSigma = TL_Sigma + reactionData->GetParticle(2);
    double mass = TL_piSigma.M()/1000.;
    Tools::H2(Form("React_q_IMPiSigma"),mass,q,900,1.2,2.1,300,0,1.5);
  }
  if(reactionID == gen::reactionID_pLpim || reactionID == gen::reactionID_pLpimpi0){
    TLorentzVector TL_pmiss = reactionData->GetParticle(2);
    double q = (TL_beam.Vect()-TL_pmiss.Vect()).Mag()/1000.;
    double costhetap = TL_pmiss.CosTheta();
    TLorentzVector TL_Lambda = reactionData->GetParticle(0);
    TLorentzVector TL_LambdaPim = TL_Lambda + reactionData->GetParticle(1);
    double mass = TL_LambdaPim.M()/1000.;
    Tools::H2(Form("React_q_IMLPim"),mass,q,900,1.2,2.1,300,0,1.5);
    Tools::H2(Form("React_costhetap_IMLPim"),mass,costhetap,900,1.2,2.1,2000,-1,1);
  }
  if(reactionID == gen::reactionID_pS0pim || reactionID == gen::reactionID_pS0pim_ps){
    TLorentzVector TL_pmiss = reactionData->GetParticle(2);
    double q = (TL_beam.Vect()-TL_pmiss.Vect()).Mag()/1000.;
    TLorentzVector TL_S0 = reactionData->GetParticle(0);
    TLorentzVector TL_S0Pim = TL_S0 + reactionData->GetParticle(1);
    double mass = TL_S0Pim.M()/1000.;
    Tools::H2(Form("React_q_IMS0Pim"),mass,q,900,1.2,2.1,300,0,1.5);
  }
  if(reactionID == gen::reactionID_npipiL){
    TLorentzVector TL_n = reactionData->GetParticle(0);
    TLorentzVector TL_pip = reactionData->GetParticle(1);
    TLorentzVector TL_pim = reactionData->GetParticle(2);
    TLorentzVector TL_Lmiss = reactionData->GetParticle(3);
    double q = (TL_beam.Vect()-TL_Lmiss.Vect()).Mag()/1000.;
    TLorentzVector TL_npipi = TL_n+TL_pip+TL_pim;
    double mass = TL_npipi.M()/1000.;
    Tools::H2(Form("React_q_IMnpipi"),mass,q,900,1.2,2.1,300,0,1.5);
  }
  if(reactionID == gen::reactionID_nK0n){
    TLorentzVector TL_K0 = reactionData->GetParticle(0);
    TLorentzVector TL_n = reactionData->GetParticle(1);
    TLorentzVector TL_nmiss = reactionData->GetParticle(2);
    double q = (TL_beam.Vect()-TL_nmiss.Vect()).Mag()/1000.;
    TLorentzVector TL_K0n = TL_K0+TL_n;
    double mass = TL_K0n.M()/1000.;
    Tools::H2(Form("React_q_IMnpipi"),mass,q,900,1.2,2.1,300,0,1.5);
  }
}

void Util::AnaMcData(MCData *mcdata, 
                     DetectorData  *detdata,
                     CDSHitMan *cdsman,
                     ReactionData *reactionData,
                     double &ncanvtxr,
                     double &ncanvtxz,
                     int &ncangeneration
                     )
{
  ncanvtxr=999.0;
  ncanvtxz=999.0;
  ncangeneration=999;

  for(int itr=0;itr<mcdata->trackSize();itr++){
    TVector3 vtx = mcdata->track(itr)->vertex();
    Tools::H1(Form("track_vtxr_cdh"),vtx.Perp()/10.,1200,0,120);
  }
  
  
  int EventType=0;
  struct pimInfo{
    int gen;
    int processID;
    int parentProcessID;
    double dE;
    double mom;
    TLorentzVector LVec;
    double time;
    pimInfo(){
      gen=-1;
      processID=-1;
      parentProcessID=-1;
      LVec.SetXYZM(999,999,999,99999);
      dE=0.0;
      mom=0.0;
      time=0.0;
    }
  };
  
  struct pipInfo{
    int gen;
    int processID;
    int parentProcessID;
    double dE;
    double mom;
    TLorentzVector LVec;
    double time;
    pipInfo(){
      gen=-1;
      processID=-1;
      parentProcessID=-1;
      LVec.SetXYZM(999,999,999,99999);
      dE=0.0;
      mom=0.0;
      time=0.0;
    }
  };
  
  //Neutron candidate info.
  //the actual particle which fires CDH is proton or electron
  //the parent should be the neutron less it is not fake
  struct NcanInfo{
    int pdg;
    int parentpdg;
    int gen;
    int processID;
    int parentProcessID;
    double dE;
    double mom;
    TLorentzVector LVec;//stores parent LVec.
    TVector3 calc_mom;
    double vtx_r;
    double vtx_r_g1parent;
    double vtx_r_g2parent;
    DetectorHit *dhitncan;
    double time;
    bool IsFromDecay;
    NcanInfo(){
      pdg=-9999;
      parentpdg=-9999;
      gen=-1;
      processID=-1;
      parentProcessID=-1;
      dE=0.0;
      mom=0.0;
      LVec.SetXYZM(999,999,999,99999);
      calc_mom.SetXYZ(999.,999.,999.);
      vtx_r=0.0;
      vtx_r_g1parent=0.0;
      vtx_r_g2parent=0.0;
      dhitncan=NULL;
      time=0.0;
      IsFromDecay=false;
    }
  };
  
  pimInfo piminfo;
  pipInfo pipinfo;
  NcanInfo ncaninfo;

  int npim=0;
  int npip=0;
  int trackIDcont[3]={-1,-1,-1};
  int iokhit=0;
  //check detectorhit at first
  //store pi+/-, neutron candidate info.
  unsigned int ncdhhit=0;
  for(int idethit=0;idethit<detdata->detectorHitSize(); idethit++){  
    DetectorHit *dhit = detdata->detectorHit(idethit);
    int cid    = dhit->detectorID();
    if(cid!=CID_CDH){
      continue;
    }
    ncdhhit++;
    int seg  = dhit->channelID();//0 origin
    double dEreco = -100.;
    int segReco = -1;
    bool isHitMatch = false;
    for(int icdhhit=0;icdhhit<cdsman->nCDH();icdhhit++){
      if(cdsman->CDH(icdhhit)->CheckRange() && cdsman->CDH(icdhhit)->ctmean()<(cdscuts::tdc_cdh_max+cdscuts::tdc_simoffset)){
        segReco = (cdsman->CDH(icdhhit)->seg())-1;//1 origin   
        if(seg==segReco){
          isHitMatch=true;
          dEreco = cdsman->CDH(icdhhit)->emean();
        }
      }//if
    }//for
    if(!isHitMatch) continue;
    //}else{
    int pdg    = dhit->pdg();
    int trackID = dhit->trackID();
    double time = dhit->time();
    double dEtrue = dhit->de();
    trackIDcont[iokhit]=trackID;
    iokhit++;
    Track *track_p  = Util::FindTrackFromMcIndex(mcdata,trackID);
    int parentID = dhit->parentID();
    Track *parenttrack_p = Util::FindTrackFromMcIndex(mcdata,parentID);
    int parentpdg = 0;
    int parentprocessID=-1;
    if(parenttrack_p != 0){
      parentpdg= parenttrack_p->pdgID();
      std::string parentProcessname = std::string(parenttrack_p->creatorProcess()); 
      parentprocessID = ProcessNameToProcessID(parentProcessname);
    }
    double truemom = (track_p->momentum()).Mag()/1000.;
    double parenttruemom = -999.0;
    TVector3 parenttruemomVec;
    if(parenttrack_p != 0){
      parenttruemom = (parenttrack_p->momentum()).Mag()/1000.;
      parenttruemomVec = parenttrack_p->momentum();//MeV/c
    }
    TVector3 truemomVec = track_p->momentum();//MeV/c
    std::string processname = std::string(track_p->creatorProcess());
    int processID =ProcessNameToProcessID(processname);
    //std::cout << "pdg dhit" << pdg << std::endl;//->OK
    //std::cout << "pdg track" << track_p->pdgID() << std::endl;//->OK
    //std::cout << "processname " << processname.c_str() << std::endl;
    //std::cout << "processID " << processID << std::endl;
    
    
   
    int generation = Util::CalcGeneration(mcdata,parentID);
    //std::cout << "gen " << generation << std::endl;

    Tools::H2(Form("CDHdE_processID"),processID, dEtrue,222,-0.5,221.5, 100,0,10);
    Tools::H2(Form("CDHdE_generation"),generation, dEtrue,10,0,10, 100,0,10);
    if(pdg>=8000){// std::cout << "large pdgID " << pdg << std::endl; 
      if(pdg==1000010020) pdg=4000;//deuteron
      else if(pdg==1000010030) pdg=4001;//triton
      else if(pdg==1000020030) pdg=4002;//3He
      else if(pdg==1000020040) pdg=4003;//alpha
      else if(pdg==1000030050) pdg=4004;//Li5
      else if(pdg==1000030060) pdg=4005;//Li6
      else if(pdg==1000030070) pdg=4006;//Li7
      else if(pdg==1000040070) pdg=4007;//Be7
      else if(pdg==1000040080) pdg=4008;//Be8
      else if(pdg==1000040090) pdg=4009;//Be9
      else if(pdg==1000040100) pdg=4010;//Be10
      else if(pdg==1000050090) pdg=4011;//B9
      else if(pdg==1000050100) pdg=4012;//B10
      else if(pdg==1000050110) pdg=4013;//B11
      else if(pdg==1000060110) pdg=4014;//C12
      else if(pdg==1000060120) pdg=4015;//C12
      else if(pdg==1000060130) pdg=4016;//C13
      else{
        std::cout << "large pdgID " << pdg << std::endl;  
      }
    } 

    Tools::H1(Form("pdgID"),pdg,16000,-8000,8000);
    //pi- track
    if(pdg== -211){
      npim++;
      piminfo.dE = dEtrue;
      piminfo.mom = truemom;
      piminfo.LVec.SetVectM(truemomVec,piMass*1000.0);
      piminfo.processID = processID;
      piminfo.parentProcessID = parentprocessID;
      piminfo.gen = generation;
      piminfo.time = time;
      Tools::H2(Form("mom_processID_pim"),processID, truemom ,222,-0.5,221.5, 200,0,1);
      Tools::H2(Form("mom_parentprocessID_pim"),parentprocessID, truemom ,222,-0.5,221.5, 200,0,1);
      if(generation>2){
        Tools::H2(Form("mom_processID_pim_over2g"),processID, truemom ,222,-0.5,221.5, 200,0,1);
        Tools::H2(Form("mom_parentprocessID_pim_over2g"),parentprocessID, truemom ,222,-0.5,221.5, 200,0,1);
      }
      Tools::H2(Form("mom_generation_pim"),generation, truemom,10,0,10, 200,0,1);
      Tools::H2(Form("CDHdE_processID_pim"),processID, dEtrue,222,-0.5,221.5, 100,0,10);
      Tools::H2(Form("CDHdE_generation_pim"),generation, dEtrue,10,0,10, 100,0,10);
    }
    //pi+ track
    else if(pdg== 211){
      npip++;
      pipinfo.dE = dEtrue;
      pipinfo.mom = truemom;
      pipinfo.LVec.SetVectM(truemomVec,piMass*1000.0);
      pipinfo.processID = processID;
      pipinfo.parentProcessID = parentprocessID;
      pipinfo.gen = generation;
      pipinfo.time = time;
      Tools::H2(Form("mom_processID_pip"),processID, truemom ,222,-0.5,221.5, 200,0,1);
      Tools::H2(Form("mom_parentprocessID_pip"),parentprocessID, truemom ,222,-0.5,221.5, 200,0,1);
      if(generation>2){
        Tools::H2(Form("mom_processID_pip_over2g"),processID, truemom ,222,-0.5,221.5, 200,0,1);
        Tools::H2(Form("mom_parentprocessID_pip_over2g"),parentprocessID, truemom ,222,-0.5,221.5, 200,0,1);
      }
      Tools::H2(Form("mom_generation_pip"),generation, truemom,10,0,10, 200,0,1);
      Tools::H2(Form("CDHdE_processID_pip"),processID, dEtrue,222,-0.5,221.5, 100,0,10);
      Tools::H2(Form("CDHdE_generation_pip"),generation, dEtrue,10,0,10, 100,0,10);
    }
    //neutron candidate 
    else{
      ncaninfo.pdg = pdg;
      ncaninfo.parentpdg = parentpdg;
      ncaninfo.dE = dEtrue;
      ncaninfo.mom = parenttruemom;
      ncaninfo.LVec.SetVectM(parenttruemomVec,nMass*1000.0);
      ncaninfo.processID = processID;
      ncaninfo.gen = generation;
      ncaninfo.dhitncan = dhit;
      ncaninfo.time = time;
      Tools::H1(Form("ncan_pdg"),pdg,16000,-8000,8000);
      Tools::H1(Form("ncan_parentpdg"),parentpdg,16000,-8000,8000);
      Tools::H2(Form("mom_processID_ncan"),processID, truemom ,222,-0.5,221.5, 200,0,2);
      Tools::H2(Form("mom_generation_ncan"),generation, truemom,10,0,10, 200,0,2);
      Tools::H2(Form("CDHdE_processID_ncan"),processID, dEtrue,222,-0.5,221.5, 100,0,10);
      Tools::H2(Form("CDHdE_generation_ncan"),generation, dEtrue,10,0,10, 100,0,10);
      Tools::H1("ncan_time",ncaninfo.time,1000,0,100);
    }
  }//for idethit
  //std::cout << __LINE__ << "ncdh hit " << ncdhhit << std::endl;

  if(npip==1 && npim==1) EventType = 1;
  //check track sharing 
  for(int i=0;i<3;i++){
    int tidi = trackIDcont[i];
    for(int j=0;j<3;j++){
      int tidj = trackIDcont[j];
      if(i!=j && tidi==tidj) EventType=2;
    }
  }
  //EventType
  //0:CDH 3 hits,pi+ or pi- does not hit CDH
  //1:CDH 3 hits,pi+,pi- hits, no CDH sharing by pi+ or pi-
  //2:CDH 3 hits,pi+,pi- hits, at least one sharing by pi+/- tracks
  Tools::H1(Form("EventType"),EventType,3,-0.5,2.5);
  if(EventType==1){
    Tools::H2(Form("mom_processID_pim_select"),piminfo.processID, piminfo.mom ,222,-0.5,221.5, 200,0,1);
    Tools::H1("time_pim_select",piminfo.time,1000,0,100);
    Tools::H2(Form("CDHdE_mom_pim_select"),piminfo.mom,piminfo.dE,100,0,2, 100,0,10);
    if(piminfo.gen>2){ 
      Tools::H2(Form("mom_processID_pim_select_over2g"),piminfo.processID, piminfo.mom ,222,-0.5,221.5, 200,0,1);
      Tools::H2(Form("mom_parentprocessID_pim_select_over2g"),piminfo.parentProcessID, piminfo.mom ,222,-0.5,221.5, 200,0,1);
      Tools::H1("time_pim_select_over2g",piminfo.time,1000,0,100);
    }
    Tools::H2(Form("mom_generation_pim_select"),piminfo.gen, piminfo.mom,10,0,10, 200,0,1);
    Tools::H2(Form("CDHdE_processID_pim_select"),piminfo.processID, piminfo.dE,222,-0.5,221.5, 100,0,10);
    Tools::H2(Form("CDHdE_generation_pim_select"),piminfo.gen, piminfo.dE,10,0,10, 100,0,10);
    
    Tools::H2(Form("mom_processID_pip_select"),pipinfo.processID, pipinfo.mom ,222,-0.5,221.5, 200,0,1);
    Tools::H1("time_pip_select",pipinfo.time,1000,0,100);
    Tools::H2(Form("CDHdE_mom_pip_select"),pipinfo.mom,pipinfo.dE,100,0,2, 100,0,10);
    if(pipinfo.gen>2){ 
      Tools::H2(Form("mom_processID_pip_select_over2g"),pipinfo.processID, pipinfo.mom ,222,-0.5,221.5, 200,0,1);
      Tools::H2(Form("mom_parentprocessID_pip_select_over2g"),pipinfo.parentProcessID, pipinfo.mom ,222,-0.5,221.5, 200,0,1);
      Tools::H1("time_pip_select_over2g",pipinfo.time,1000,0,100);
    }
    Tools::H2(Form("mom_generation_pip_select"),pipinfo.gen, pipinfo.mom,10,0,10, 200,0,1);
    Tools::H2(Form("CDHdE_processID_pip_select"),pipinfo.processID, pipinfo.dE,222,-0.5,221.5, 100,0,10);
    Tools::H2(Form("CDHdE_generation_pip_select"),pipinfo.gen, pipinfo.dE,10,0,10, 100,0,10);
    
    Tools::H1(Form("ncan_pdg_select"),ncaninfo.pdg,8000,-4000,4000);
    Tools::H1(Form("ncan_parentpdg_select"),ncaninfo.parentpdg,8000,-4000,4000);
    Tools::H2(Form("mom_processID_ncan_select"),ncaninfo.processID, ncaninfo.mom ,222,-0.5,221.5, 100,0,2);
    Tools::H2(Form("mom_generation_ncan_select"),ncaninfo.gen, ncaninfo.mom,10,0,10, 100,0,2);
    Tools::H2(Form("CDHdE_processID_ncan_select"),ncaninfo.processID, ncaninfo.dE,222,-0.5,221.5, 100,0,10);
    Tools::H2(Form("CDHdE_generation_ncan_select"),ncaninfo.gen, ncaninfo.dE,10,0,10, 100,0,10);
    Tools::H2(Form("CDHdE_mom_ncan_select"),ncaninfo.mom,ncaninfo.dE,200,0,2, 100,0,10);
    Tools::H1("ncan_time_select",ncaninfo.time,1000,0,100);
    Tools::H1("ncan_pim_time_select",ncaninfo.time-piminfo.time,1000,-100,100);
    Tools::H1("ncan_pip_time_select",ncaninfo.time-pipinfo.time,1000,-100,100);
    double parentvtxr=Util::FillAncestryVertexR(mcdata,ncaninfo.dhitncan,ncaninfo.dE);
    double parentvtxz=Util::FillAncestryVertexZ(mcdata,ncaninfo.dhitncan,ncaninfo.dE);
    if(ncaninfo.dE>2.0){
      //Tools::H2("ncan_mom_parentvtxr_select",ncaninfo.mom,parentvtxr,200,0,2,480,0,120.);
      //Tools::H2("ncan_mom_parentvtxz_select",ncaninfo.mom,parentvtxz,200,0,2,750,-750,750);
      if(ncaninfo.parentpdg==2112){
        //Tools::H2("ncan_mom_parentvtxr_select_nparent",ncaninfo.mom,parentvtxr,200,0,2,480,0,120.);
      }
    }
    if(Util::IsFromSigma(mcdata,ncaninfo.dhitncan)){
      ncanvtxr = parentvtxr;
      ncanvtxz = parentvtxz;
      ncangeneration = ncaninfo.gen;
      //Tools::H2(Form("CDHdE_generation_ncan_select_sigma"),ncaninfo.gen, ncaninfo.dE,10,0,10, 100,0,10);
      //Tools::H2(Form("vtxr_generation_ncan_select_sigma"),ncaninfo.gen, parentvtxr,10,0,10, 480,0,120.0);
      TLorentzVector LVec_n_pim_mc = ncaninfo.LVec + piminfo.LVec;
      TLorentzVector LVec_n_pip_mc = ncaninfo.LVec + pipinfo.LVec;
      TLorentzVector TL_Sigma = reactionData->GetParticle(1);
      TLorentzVector TL_diff_n_pim = LVec_n_pim_mc - TL_Sigma;
      TLorentzVector TL_diff_n_pip = LVec_n_pip_mc - TL_Sigma;
      //std::cout << "MC " << LVec_n_pim_mc.P() << std::endl;
      //std::cout << "react " << TL_Sigma.P() << std::endl;
      //std::cout << "MC_M " << LVec_n_pim_mc.M() << std::endl;
      //std::cout << "react_M " << TL_Sigma.M() << std::endl;
      double diffmom_npim = (LVec_n_pim_mc.P()-TL_Sigma.P())*0.001;
      double diffmom_npip = (LVec_n_pip_mc.P()-TL_Sigma.P())*0.001;
      //Tools::H2(Form("vtxr_diffmom_npim_ncan_select_sigma"),diffmom_npim, parentvtxr,200,-0.1,0.1, 480,0,120.0);
      //Tools::H2(Form("vtxr_diffmom_npip_ncan_select_sigma"),diffmom_npip, parentvtxr,200,-0.1,0.1, 480,0,120.0);
      if(ncaninfo.gen==3){
        //Tools::H2(Form("vtxr_diffmom_npim_ncan_select_sigma_3rdgen"),diffmom_npim, parentvtxr,200,-0.1,0.1, 480,0,120.0);
        //Tools::H2(Form("vtxr_diffmom_npip_ncan_select_sigma_3rdgen"),diffmom_npip, parentvtxr,200,-0.1,0.1, 480,0,120.0);
      }
      //Tools::H2(Form("generation_diffmom_npim_ncan_select_sigma"),diffmom_npim,ncaninfo.gen,200,-0.1,0.1, 10,0,10);
      //Tools::H2(Form("generation_diffmom_npip_ncan_select_sigma"),diffmom_npip,ncaninfo.gen,200,-0.1,0.1, 10,0,10);
      if(ncaninfo.dE>2.0){
        //Tools::H2(Form("vtxr_generation_ncan_select_sigma_dE"),ncaninfo.gen, parentvtxr,10,0,10, 480,0,120.0);
      }
    }
    if(parentvtxr<58.0){
      //Tools::H2("CDHdE_mom_ncan_select_vtxr",ncaninfo.mom,ncaninfo.dE,200,0,2, 100,0,10);
      //if(ncaninfo.dE>0.6 && Util::IsFromSigma(mcdata,ncaninfo.dhitncan) && ncaninfo.gen==3)Tools::H2("CDHdE_mom_ncan_select_vtxr_sigma",ncaninfo.mom,ncaninfo.dE,100,0,2, 100,0,10);
    }
  }//event Type 1


  return;
}


void Util::AnaMcData2(MCData *mcdata, 
                     DetectorData  *detdata,
                     const int CDHseg,
                     double &ncanvtxr,
                     double &ncanvtxz,
                     double &firstvtxr,
                     double &firstvtxz,
                     int &ncangeneration,
                     int &mcpattern,
                     int &nanc
                     )
{
  ncanvtxr=999.0;
  ncanvtxz=999.0;
  firstvtxr=999.0;
  firstvtxz=999.0;
  ncangeneration=999;
  mcpattern=999;
  nanc=999;

  
  int trackID = -1;
  int pdg = 0;
  int ncandidate = 0;
  int parentID = 0;
  int generation = 0;
  int generation_init = 0;
  Track *track_p=0;
  //bool isPip=false;
  //bool isPim=false;
  //TVector3 mom_pip;
  //TVector3 mom_pim;
  //TLorentzVector LVec_pip;
  //TLorentzVector LVec_pim;
  for(int idethit=0;idethit<detdata->detectorHitSize(); idethit++){  
    DetectorHit *dhit = detdata->detectorHit(idethit);
    int cid    = dhit->detectorID();
    if(cid!=CID_CDH){
      continue;
    }
    //require CDH hit seg matching
    int dhitseg  = dhit->channelID()+1;//0 origin
    if(dhitseg != CDHseg) continue;
    
    trackID = dhit->trackID();
    pdg    = dhit->pdg();
    track_p = Util::FindTrackFromMcIndex(mcdata, trackID);
    parentID = track_p->parentTrackID();
    /*
    if(pdg==-211){
      isPim=true;
      mom_pim = dhit->momentum();
      LVec_pim.SetVectM(mom_pim*0.001,piMass);
    }
    if(pdg== 211){
      isPip=true;
      mom_pip = dhit->momentum();
      LVec_pip.SetVectM(mom_pip*0.001,piMass);
    }*/
    //parentID of detector 
    if(dhit->parentID() != parentID) std::cout << __LINE__ << " somothing is wrong !!!! " << std::endl;
    generation = Util::CalcGeneration(mcdata,parentID);
    generation_init = generation;
    ncandidate++;
  }
  
  //this is 
  if(parentID==0){
    std::cout << __FILE__ << " L." << __LINE__ << " primary particle is detected: " << pdg << std::endl;
  }
  if(ncandidate!=1){
    std::cout << __FILE__ << " L." << __LINE__ << " !!! ncandidate =  " << ncandidate << std::endl;
  }
  //std::cout << __FILE__ << " L." << __LINE__ << std::endl;
  //std::cout << "parentID " << parentID << std::endl;
  //std::cout << "CDH 1st pdg " << pdg << std::endl;
  //std::cout << "generation " << generation << std::endl;

  if( (pdg== -211) || (pdg==211) ){
    //std::cout << __FILE__<<  " L." << __LINE__ << " pion hit return" << std::endl;
    return; 
  }
 

  //std::cout << __LINE__ << "  " << pdg << std::endl;
  Track *AncestorTr[32]={0};//array index <= generation
  int AncestorTrackID[32]={0};
  AncestorTrackID[0]=parentID;
  
  TVector3 AncestorVTX[32];
  int AncestorPDG[32]={0};
  unsigned int anc=0;
  bool isFromNeutron = false;
  bool isFromSigma = false;
  bool isSigmaNeutronChain = false;
  bool isFromPion = false; 
  double vtxRNeutron = 999.0;//cm
  double vtxZNeutron = 999.0;//cm
  unsigned int nNeutrons = 0;
  int genNeutron = -1;
  bool isWentCDHOutSide = false;
  bool isinFiducialORinCDH = true;

  while(true){
    //std::cout << __LINE__ << std::endl;
    if(generation == 1) break;
    else generation--;
    //In the first loop, 
    //get the Track of the parent particle which fired CDH
    AncestorTr[anc]=Util::FindTrackFromMcIndex(mcdata,AncestorTrackID[anc]);
    AncestorVTX[anc] = AncestorTr[anc]->vertex();
    AncestorVTX[anc] *=0.1;
    AncestorPDG[anc] = AncestorTr[anc]->pdgID();
    Track *parenttrack_p = Util::FindTrackFromMcIndex(mcdata,AncestorTrackID[anc]);
    
    //if(GeomTools::GetID(AncestorVTX[anc])==CID_CDH){
    //std::cout << __LINE__ << std::endl;
    //std::cout << "GeomID " << GeomTools::GetID(AncestorVTX[anc]) << std::endl;
    //std::cout << "vtxr " << AncestorVTX[anc].Perp() << std::endl;
    //std::cout << "vtxz " << AncestorVTX[anc].Z() << std::endl;
    //std::cout << std::endl;
    //}

    //check vertex pos.
    //Since Sigma has some flight length, it can go out side of target volume
    if(GeomTools::GetID(AncestorVTX[anc])==CID_TarChm || 
       //GeomTools::GetID(AncestorVTX[anc])==CID_TarSys ||
       GeomTools::GetID(AncestorVTX[anc])==CID_RadS ||
       GeomTools::GetID(AncestorVTX[anc])==CID_TarCFRP ||
       GeomTools::GetID(AncestorVTX[anc])==CID_TarCap// ||
       //GeomTools::GetID(AncestorVTX[anc])==CID_TarCell 
    ){
      isinFiducialORinCDH = false;
    }else if( (AncestorVTX[anc].Perp() > 10 ) &&
             // ((AncestorVTX[anc].Z() < -15.0) || (5 < AncestorVTX[anc].Z())) && 
              (GeomTools::GetID(AncestorVTX[anc])!=CID_CDH)
              ){
      isinFiducialORinCDH = false;
    }

    if(anc==0){
      Tools::H2(Form("vtxrz_cdhhitparenet"),AncestorVTX[anc].Z(),
          AncestorVTX[anc].Perp(),
          750,-75.,75.,480,0.,120.);
      Tools::H1(Form("gen_cdhhitparent"),generation,10,0,10);
      firstvtxr = AncestorVTX[anc].Perp();
      firstvtxz = AncestorVTX[anc].Z();
    }
    if( (AncestorVTX[anc].Perp()) > 58.0) isWentCDHOutSide = true;
    if( fabs(AncestorVTX[anc].Z()) > 40.0) isWentCDHOutSide = true; 
    
    //2112 :Neutron 
    if( AncestorPDG[anc]==2112){
      nNeutrons++;
      isFromNeutron = true;
      genNeutron = generation;
      if(nNeutrons==1){
        vtxRNeutron = AncestorVTX[anc].Perp();
        vtxZNeutron = AncestorVTX[anc].Z();
        Tools::H2(Form("vtxrz_n"),AncestorVTX[anc].Z(),
            AncestorVTX[anc].Perp(),
            750,-75.,75.,480,0.,120.);
      }
    }
    /*
    std::cout << __LINE__ << " gen " << generation << std::endl;
    std::cout << __LINE__ << " pdg " << AncestorPDG[anc] << std::endl;
    std::cout << __LINE__ << " parentTrackID " << AncestorTrackID[anc] << std::endl;
    std::cout << __LINE__ << " parentID " << AncestorTrackID[anc] << std::endl;
    std::cout << __LINE__ << " parentID_2 " << parenttrack_p->parentTrackID() << std::endl;
    std::cout << __LINE__ << " anc "  << anc << std::endl;
    std::cout << __LINE__ << " R "  << AncestorVTX[anc].Perp() << std::endl;
    std::cout << __LINE__ << " Z "  << AncestorVTX[anc].Z() << std::endl;
    std::cout <<           "Pim " << isPim << std::endl;
    std::cout <<           "Pip " << isPip << std::endl;
    std::cout <<         "CDHout " << isWentCDHOutSide << std::endl;
    */
    //check originated from sigma+/-
    if( (AncestorPDG[anc]==3112) || (AncestorPDG[anc]==3222) ){
      isFromSigma = true;
      if(generation==genNeutron-1) isSigmaNeutronChain = true;
    }
    //check also sigma -> pion -> .... -> neutron Path
    if( (AncestorPDG[anc]==211) || (AncestorPDG[anc]==-211)) isFromPion = true;
    if( (generation == 1)  && !isWentCDHOutSide){
      //std::cout << "gen " << generation << std::endl;
      //std::cout << "anc " << anc << std::endl;
      //std::cout << "pdg " << AncestorPDG[anc] << std::endl;
      //std::cout << "generation_init " << generation_init << std::endl;
      //std::cout << "isFromSigma " << isFromSigma << std::endl;
      //std::cout << "isFromPion " << isFromPion << std::endl;
    }
    AncestorTrackID[anc+1]=AncestorTr[anc]->parentTrackID();
    anc++;
  }

  //if(nNeutrons!=1) std::cout << __LINE__ << " # of neutrons " <<  nNeutrons << std::endl; 

  int pattern=0;
  if( isFromNeutron && !isFromSigma  && !isFromPion && isWentCDHOutSide ) pattern=1;//BG neutron
  else if( isFromNeutron && isFromSigma && !isFromPion && isSigmaNeutronChain && !isWentCDHOutSide && isinFiducialORinCDH ) pattern=2;//Signal 
  else if( isFromNeutron && isFromSigma && !isFromPion && isSigmaNeutronChain && isWentCDHOutSide ) pattern=3;//BG 
  else if( isFromNeutron && isFromSigma && !isFromPion && !isSigmaNeutronChain ) pattern=4;//BG
  else if( isFromNeutron && !isFromSigma && isFromPion ) pattern=5;//BG
  else if( isFromNeutron && isFromSigma  && isFromPion ) pattern=6;//BG
  else if( isFromNeutron && !isFromSigma && !isFromPion && !isWentCDHOutSide && isinFiducialORinCDH) pattern=7;//initial neutron
  //std::cout << "pattern = " << pattern << std::endl;
  Tools::H1(Form("NfakePattern"),pattern,10,-0.5,9.5);
  
  //std::cout << "  " << std::endl;
  /*
  if(pattern==7){
    std::cout << "gen " << genNeutron << std::endl;
    std::cout << "anc " << anc << std::endl;
    std::cout << "generation_init " << generation_init << std::endl;
    std::cout << "isFromSigma " << isFromSigma << std::endl;
    std::cout << "isFromPion " << isFromPion << std::endl;
    std::cout << "vtx R " << vtxRNeutron << std::endl;
    std::cout << "vtx Z " << vtxZNeutron << std::endl;
  }*/

  ncanvtxr = vtxRNeutron;
  ncanvtxz = vtxZNeutron;
  ncangeneration = genNeutron;
  mcpattern = pattern;
  nanc = anc;

}


int Util::ProcessNameToProcessID(const std::string &name)
{
  if( name=="" ) return 0; //Initial
  else if( name=="Decay"                 ) return 1;
  else if( name=="conv"                  ) return 2;
  else if( name=="Transportation"        ) return 3;
  else if( name=="phot"                  ) return 4;
  else if( name=="annihil"               ) return 5;
  else if( name=="compt"                 ) return 6;
  else if( name=="eBrem"                 ) return 7;
  else if( name=="hadElastic"            ) return 8;
  else if( name=="CoulombScat"           ) return 9;
  else if( name=="nKiller"               ) return 10;
  else if( name=="photonNuclear"         ) return 11;
  else if( name=="msc"                   ) return 12;
  else if( name=="pi-Inelastic"          ) return 100;
  else if( name=="pi+Inelastic"          ) return 101;
  else if( name=="kaon-Inelastic"        ) return 102;
  else if( name=="kaon+Inelastic"        ) return 103;
  else if( name=="kaon0LInelastic"       ) return 104;
  else if( name=="kaon0SInelastic"       ) return 105;
  else if( name=="lambdaInelastic"       ) return 106;
  else if( name=="sigma+Inelastic"       ) return 107;
  else if( name=="sigma-Inelastic"       ) return 108;
  else if( name=="sigma0Inelastic"       ) return 109;
  else if( name=="protonInelastic"       ) return 110;
  else if( name=="neutronInelastic"      ) return 111;
  else if( name=="dInelastic"            ) return 112;
  else if( name=="tInelastic"            ) return 113;
  else if( name=="He3Inelastic"          ) return 114;
  else if( name=="alphaInelastic"        ) return 115;
  else if( name.find("Inelastic")!=std::string::npos ) return 199;
  else if( name=="eIoni"                 ) return 200;
  else if( name=="hIoni"                 ) return 201;
  else if( name=="ionIoni"               ) return 202;
  else if( name=="muIoni"                ) return 203;
  else if( name=="hBertiniCaptureAtRest" ) return 204;
  else if( name=="nCapture"              ) return 205;
  else if( name=="muMinusCaptureAtRest"  ) return 206;
  else if( name=="unknown"               ) return -999;
  else return -1;
}

std::string Util::ProcessIDToProcessName(const int &id)
{
  if( id==0                 ) return "initial";
  else if( id==1            ) return "Decay";
  else if( id==2            ) return "conv";
  else if( id==3            ) return "Transportation";
  else if( id==4            ) return "phot";
  else if( id==5            ) return "annihil";
  else if( id==6            ) return "compt";
  else if( id==7            ) return "eBrem";
  else if( id==8            ) return "hadElastic";
  else if( id==9            ) return "CoulombScat";
  else if( id==10           ) return "nKiller";
  else if( id==11           ) return "photoNuclear";
  else if( id==12           ) return "msc";
  else if( id==100          ) return "pi-Inelastic";
  else if( id==101          ) return "pi+Inelastic";
  else if( id==102          ) return "kaon-Inelastic";
  else if( id==103          ) return "kaon+Inelastic";
  else if( id==104          ) return "kaon0Inelastic";
  else if( id==105          ) return "kaon0SInelastic";
  else if( id==106          ) return "lambdaInelastic";
  else if( id==107          ) return "sigma-Inelastic";
  else if( id==108          ) return "sigma+Inelastic";
  else if( id==109          ) return "sigma0Inelastic";
  else if( id==110          ) return "protonInelastic";
  else if( id==111          ) return "neutronInelastic";
  else if( id==112          ) return "dInelastic";
  else if( id==113          ) return "tInelastic";
  else if( id==114          ) return "He3Inelastic";
  else if( id==115          ) return "alphaInelastic";
  else if( 100<id && id<200 ) return "Inelastic";
  else if( id==200          ) return "eIoni";
  else if( id==201          ) return "hIoni";
  else if( id==202          ) return "ionIoni";
  else if( id==203          ) return "muIoni";
  else if( id==204          ) return "hBertiniCaptureAtRest";
  else if( id==205          ) return "nCapture";
  else if( id==206          ) return "muMinusCaptureAtRest";
  else if( id==-999         ) return "unknown";
  else return "error";
}

//parentID from detectorhit node = parenttrackID of Track node
int Util::CalcGeneration(MCData *mcdata,int parentID)
{
  int gen=1;
  Track *parentTr=Util::FindTrackFromMcIndex(mcdata, parentID);
  while( parentTr!=NULL){
    //std::cout << "parent ID " << parentTr->parentTrackID() << std::endl;
    parentTr=Util::FindTrackFromMcIndex(mcdata,parentTr->parentTrackID());
    gen++;
  }
  return gen;
}

//input 
//MCData node
//trackid 
//
//output
//pointer of the Track node
Track* Util::FindTrackFromMcIndex(MCData *mcdata, int trackid)
{
  for( int itr=0;itr<mcdata->trackSize();itr++){
    Track *track=mcdata->track(itr);
    if( track->trackID()==trackid) return track;
  }
  //std::cout << __FILE__ << "  L." <<  __LINE__ << " track NOT found !!! " << std::endl; 
  return NULL;
}


double Util::FillAncestryVertexR(MCData *mcdata,DetectorHit *dhit, double dE)
{
  //if(dE<0.06) return 999;
  Track *parentTr=Util::FindTrackFromMcIndex(mcdata, dhit->parentID());
  
  TVector3 vtxp = parentTr->vertex();
  Tools::H2(Form("dE_track_vtxr_ncan_1stg"),vtxp.Perp()/10.,dE,1200,0,120,100,0,10);
  int gen=0;
  while( parentTr!=0){
    TVector3 vtx = parentTr->vertex();
    Tools::H2(Form("dE_track_vtxr_ncan"),vtx.Perp()/10.,dE,1200,0,120,100,0,10);
    Tools::H2(Form("time_track_vtxr_ncan"),vtx.Perp()/10.,dhit->time(),1200,0,120,1000,0,100);
    Tools::H2(Form("gen_track_vtxr_ncan"),vtx.Perp()/10.,gen+1,1200,0,120,10,0,10);
    parentTr=Util::FindTrackFromMcIndex(mcdata,parentTr->parentTrackID());
    gen++;
  }
  return vtxp.Perp()/10.;
}

double Util::FillAncestryVertexZ(MCData *mcdata,DetectorHit *dhit, double dE)
{
  Track *parentTr=Util::FindTrackFromMcIndex(mcdata, dhit->parentID());
  TVector3 vtxp = parentTr->vertex();
  Tools::H2(Form("dE_track_vtxz_ncan_1stg"),vtxp.Z()/10.,dE,1200,-200,400,100,0,10);
  int gen=0;
  while( parentTr!=0){
    TVector3 vtx = parentTr->vertex();
    Tools::H2(Form("dE_track_vtxz_ncan"),vtx.Z()/10.,dE,1200,-200,400,100,0,10);
    Tools::H2(Form("time_track_vtxz_ncan"),vtx.Z()/10.,dhit->time(),1200,-200,400,1000,0,100);
    Tools::H2(Form("gen_track_vtxz_ncan"),vtx.Z()/10.,gen+1,1200,-200,400,10,0,10);
    parentTr=Util::FindTrackFromMcIndex(mcdata,parentTr->parentTrackID());
    gen++;
  }
  return vtxp.Z()/10.;
}


bool Util::IsFromSigma(MCData *mcdata,DetectorHit *dhit)
{
  bool isFromSigma=false;
  Track *parentTr=Util::FindTrackFromMcIndex(mcdata, dhit->parentID());
  if( (parentTr->pdgID()==3112) || (parentTr->pdgID()==3222)) isFromSigma = true;

  while(parentTr!=0){
    if( (parentTr->pdgID()==3112) || (parentTr->pdgID()==3222)) isFromSigma = true;
    parentTr=Util::FindTrackFromMcIndex(mcdata,parentTr->parentTrackID());
  }

  return isFromSigma;
}


/*
TLorentzVector* Util::GetForwardNeutralLVec(BeamLineHitMan *blman,const TVector3 vtxpos, const double t0time,const double beamtof,const double thre)
{
  TLorentzVector *fNLVec = NULL;
  
  double fTime=0.0;
  std::vector<std::vector<HodoscopeLikeHit*> > NChits=Util::getNChits(blman);
  std::vector<HodoscopeLikeHit*>  NChits2=Util::getHodo(blman);
  
  HodoscopeLikeHit *nc_hit=0;
  for( int lay=0; lay<(int)NChits.size(); lay++ ){
    //    cout<<"NC layer"<<lay<<" search "<<endl;
    double time=DBL_MAX;
    for( int i=0; i<(int)NChits[lay].size(); i++ ){
      if( NChits[lay][i]->emean()>thre && NChits[lay][i]->ctmean()<time ){
        time=NChits[lay][i]->ctmean();
        nc_hit=NChits[lay][i];
      }
    }
    if( nc_hit ) break;
  }

  if( !nc_hit ) return fNLVec;

  std::vector<int> nc_seg;
  bool isAdd=true;
  nc_seg.push_back(nc_hit->seg());
  while(isAdd){
    isAdd=false;
    for(int ihit=0;ihit<(int)NChits2.size();ihit++){
      int seg=NChits2[ihit]->seg();
      int lay=(seg-1)/16;
      int seg2=seg-16*lay;
      for(int jhit=0;jhit<(int)nc_seg.size();jhit++){
        int cl_lay=(nc_seg[jhit]-1)/16;
        int cl_seg2=nc_seg[jhit]-16*cl_lay;
        if(seg==nc_seg[jhit]){
          isAdd=false;
          break;
        }
        if( abs(lay-cl_lay)==1 && cl_seg2==seg2) isAdd=true;
        if( lay==cl_lay && abs(cl_seg2-seg2)==1) isAdd=true;
      }
      if(isAdd){
        nc_seg.push_back(seg);
        break;
      }
    }
  }

  TVector3 NChitPos = nc_hit->pos();
  double ForwardNTof = fTime-t0time-beamtof;
  double Forwardl = (NChitPos-vtxpos).Mag();
  if(MCFlag)Forwardl-=2.5;
  double ForwardBeta=Forwardl/(ForwardNTof*100.*Const);
  //Neutron
  Tools::H1("ForwardBetaInv",1/ForwardBeata);
  TVector3 momvec;
  double mom=nMass*ForwardBeta/sqrt(1.0-ForwardBeta*ForwardBeata);
  momvec = NhitPos-vtxpos;
  momvec.SetMag(mom);

  

  return fNLVec;
}


  
std::vector<std::vector<HodoscopeLikeHit*> > Util::getNChits(BeamLineHitMan *blman)
{
  std::vector<std::vector<HodoscopeLikeHit*> > nchit(7);
  for( int i=0; i<blman->nNC(); i++ ){
    if( blman->NC(i)->CheckRange() ){
      int seg=blman->NC(i)->seg();
      int lay=(seg-1)/16;
      if( 7<lay ) std::cout<<"  !!! NC layer>7    "<<lay<<" !!!"<<std::endl;
      nchit[lay].push_back(blman->NC(i));
    }
  }
  return nchit;
}

std::vector<HodoscopeLikeHit*> Util::getHodo(BeamLineHitMan *blman)
{
  std::vector<HodoscopeLikeHit*> hits;
  for( int i=0; i<blman->nHodo(CID_NC); i++ ){
    HodoscopeLikeHit *hit=blman->Hodoi(CID_NC, i);
    //    cout<<"hit="<<hit<<endl;
    if( hit-> CheckRange() ){
      hits.push_back(hit);
    }
  }
  return hits;
}*/


