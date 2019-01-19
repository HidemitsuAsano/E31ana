#include "IMPiSigmaUtil.h"
#include "Tools.h"
#include "TrackTools.h"

//analyze # of CDH hits within the time range and judge if it is cdscuts::cdhmulti defined in IMPiSigmaAnaPar.h
bool Util::EveSelectCDHMul(CDSHitMan *cdsman)
{
  //** # of CDH-hits cut **//
  int nCDH = 0;
  for( int i=0; i<cdsman->nCDH(); i++ ) {
    Tools::Fill2D(Form("CDHtime"),cdsman->CDH(i)->seg(),cdsman->CDH(i)->ctmean());
    //if( cdsman->CDH(i)->CheckRange() ) nCDH++; //** only requirement of TDC **//
    if( cdsman->CDH(i)->CheckRange() && cdsman->CDH(i)->ctmean()<cdscuts::tdc_cdh_max ) {
      nCDH++;
    }
  }
  Tools::Fill1D( Form("mul_CDH"), nCDH );
  if( nCDH != cdscuts::cdhmulti  ) {
    return false;
  }

  return true;
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
int Util::GetCDHNeighboringNHits(const std::vector <int> seg, const std::vector <int> allhit  )
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
      double angle = Pos_CDC.Angle(PosCDH); // rad
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
  double beammom = beamsp->mom();
  double bchi = beamsp->chisquare();
  delete beamsp;

  Tools::Fill1D( Form("trackchi2_beam"), bchi );
  if( bchi>blcuts::d5_chi2_max ) {
    return -9999.;
  }
  Tools::Fill1D( Form("momentum_beam"), beammom );

  return beammom;
}



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
                        std::vector <int> &protonid)
{
  int CDHseg=-1;
  bool chi2OK = true;
  bool CDHseg1hitOK = true;
  bool CDHsharecheckOK = true;
  bool FindMass2OK1 = true;
  bool FindMass2OK2 = true;
  bool EnergyLossOK = true;

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
    TVector3 vtxbline, vtxbhelix; //,vtxb;
    track->GetVertex( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), vtxbline, vtxbhelix );
    track->SetPID(-1);
    Tools::Fill2D(Form("Vtx_ZX"),vtxbline.Z(),vtxbline.X());
    Tools::Fill2D(Form("Vtx_ZY"),vtxbline.Z(),vtxbline.Y());
    Tools::Fill2D(Form("Vtx_XY"),vtxbline.X(),vtxbline.Y());
    //vtxb = (vtxbline+vtxbhelix)*0.5;//not used, so far

    double tof = 999.;
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

    int pid = -1;
    pid = TrackTools::PIDcorr_wide(mom,mass2);

    track->SetPID( pid );
    Tools::Fill2D( "PID_CDS_beta", 1/beta_calc, mom );
    Tools::Fill2D( "PID_CDS", mass2, mom );


    // Energy loss calculation //
    double tmpl=0;
    TVector3 vtx_beam, vtx_cds;
    if( !track->CalcVertexTimeLength( bpctrack->GetPosatZ(0), bpctrack->GetMomDir(), track->Mass(),
                                      vtx_beam, vtx_cds, tmptof, tmpl, true ) ) {
      std::cerr<<" !!! failure in energy loss calculation [CalcVertexTimeLength()] !!! "<<std::endl;
      EnergyLossOK = false;
      continue;
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

