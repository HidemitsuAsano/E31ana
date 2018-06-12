#include "HistManwMC.h"

static const TLorentzVector TGT_LMOM = TLorentzVector(0., 0., 0., dMass);

HistManwMC::HistManwMC(TFile *f, ConfMan *conf)
  : rtFile(f), confMan(conf), runHeaderMC(0), evHeaderMC(0), mcData(0), detData(0), reacData(0), blMan(0), cdsMan(0), bltrackMan(0), cdstrackMan(0)
{
  Clear();
}

void HistManwMC::initHist()
{
  rtFile-> cd();
  std::cout<<"===== HistManwMC::initHist START  ====="<<std::endl;
  new TH2F("CDS_mass2_mom", "CDS mass2 vs mom", 2000, -0.5, 9.5, 2000, -1.0, 1.0);

  new TH1F("NC_overbeta_8MeVee", "NC 1/#beta dE>8MeVee", 1000, 0.0, 10.0);
  new TH1F("NC_overbeta_8MeVee_wtar", "NC 1/#beta dE>8MeVee", 1000, 0.0, 10.0);

  new TH1F("KN_MM", "d(K^{-}, n)", 2000, 0.0, 2.0);
  std::cout<<"===== HistManwMC::initHist FINISH ====="<<std::endl;
}

bool HistManwMC::ana()
{
  TH2F *h2;
  TH1F *h1;
  if( fT0time==DBL_MIN ) return false;
  if( fD5mom==DBL_MIN ) return false;
  if( !fTrackBPC ) return false;

  for( int i=0; i<blMan->nT0(); i++   ) if( blMan->T0(i)->CheckRange()   ) fT0_hit.push_back(blMan->T0(i));
  for( int i=0; i<blMan->nBPD(); i++  ) if( blMan->BPD(i)->CheckRange()  ) fBPD_hit.push_back(blMan->BPD(i));
  for( int i=0; i<cdsMan->nCDH(); i++ ) if( cdsMan->CDH(i)->CheckRange() ) fCDH_hit.push_back(cdsMan->CDH(i));
  for( int i=0; i<blMan->nBVC(); i++  ) if( blMan->BVC(i)->CheckRange()  ) fBVC_hit.push_back(blMan->BVC(i));
  for( int i=0; i<blMan->nCVC(); i++  ) if( blMan->CVC(i)->CheckRange()  ) fCVC_hit.push_back(blMan->CVC(i));
  for( int i=0; i<blMan->nPC(); i++   ) if( blMan->PC(i)->CheckRange()   ) fPC_hit.push_back(blMan->PC(i));
  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->CheckRange() ){
      int seg = blMan->NC(i)->seg();
      int lay = (seg-1)/15;
      fNC_hit[lay].push_back(blMan->NC(i));
    }
  }

  fT0pos = fTrackBPC->GetPosatZ(-110.5);

  double good_dis=DBL_MAX;
  bool good_vtx =false;
  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *cdc = cdstrackMan->GoodTrack(i);
    double CDHtime, dis;
    double CDHdE=0;
    int CDHseg;
    if( !cdc-> GetCDHHit(cdsMan, CDHseg, CDHtime) ) continue;
    for( int i=0; i<cdc->nCDHHit(); i++ ) CDHdE += cdc->CDHHit(cdsMan, i)->emean();
    double par[5];
    TVector3 vtxCDS, vtxBeam;
    cdc-> GetParameters(CID_CDC, par, vtxCDS);
    double mom = cdc-> Momentum();
    TrackTools::CalcLineHelixVertex(fTrackBPC, cdc, vtxBeam, vtxCDS, dis);
    if( dis<good_dis ){
      good_vtx = true;
      good_dis=dis;
      fVtxBeam = vtxBeam;
      fVtxCDS = vtxCDS;
    }

    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, vtxBeam, fD5mom, kpMass, beam_out, beam_tof);

    TVector3 vtxCDH = cdc-> CDHVertex();
    double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, vtxCDS);
    double beta = cdc_dis/(CDHtime-fT0time)/(100.*Const);
    double mass2 = mom*mom*(1./(beta*beta)-1);

    double calc_beta, tofvtxcdc;
    if( !TrackTools::FindMass2C(cdc, fTrackBPC, CDHtime-fT0time, fD5mom, kpMass, calc_beta, mass2, tofvtxcdc) ) continue;
    int pid=TrackTools::PID(mom, mass2);
    cdc-> SetPID(pid);
    if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom"), h2-> Fill(mass2, mom);
    }

    TVector3 vtxCDS2, vtxBeam2;
    double tof, tmpl;
    if( !cdc->CalcVertexTimeLength(fT0pos, fTrackBPC->GetMomDir(), cdsMass[pid], vtxCDS2, vtxBeam2, tof, tmpl, true) ) continue;
    if( pid==CDS_PiMinus  ) fTrackPim.push_back(cdc);
    if( pid==CDS_Kaon     ) fTrackKm.push_back(cdc);
    if( pid==CDS_PiPlus   ) fTrackPip.push_back(cdc);
    if( pid==CDS_Proton   ) fTrackP.push_back(cdc);
    if( pid==CDS_Deuteron ) fTrackD.push_back(cdc);
  }
  if( good_vtx ){
    TVector3 beam_mom = fTrackBPC->GetMomDir();
    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, kpMass, beam_out, beam_tof);
    beam_mom.SetMag(beam_out);
    fBeamLmom.SetVectM(beam_mom, kpMass);
  }

  if( good_vtx && fCVC_hit.size()==0 && fBVC_hit.size()==0 && nNC()>0 ){
    bool NC_flag = false;
    for( int lay=0; lay<8; lay++ ){
      fNCtime = DBL_MAX;
      for( int i=0; i<fNC_hit[lay].size(); i++ ){
	if( fNC_hit[lay][i]->emean()>8.0 && fNC_hit[lay][i]->ctmean()<fNCtime){
	  NC_flag =true;
	  fNCtime =fNC_hit[lay][i]->ctmean();
	  fNCdE = fNC_hit[lay][i]->emean();
	  confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNC_hit[lay][i]->seg() , fNCpos);
	  //	  fNCpos.SetY(fNC_hit[lay][i]->hitpos());
	}
      }
      if( NC_flag ) break;
    }
    if( NC_flag ){
      double fl = (fNCpos-fVtxBeam).Mag();
      double beam_out, beam_tof;
      ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, kpMass, beam_out, beam_tof);
      fNCbeta = fl/((fNCtime-fT0time-beam_tof)*100.*Const);

      h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee"), h1-> Fill(1./fNCbeta);
      if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_wtar"), h1-> Fill(1./fNCbeta);

      if( fNCbeta>0.9 ) fFPID=F_Gamma;
      else{
	fFPID=F_Neutron;
	double NCmom = nMass*fNCbeta/sqrt(1-fNCbeta*fNCbeta);
	TVector3 n_mom = fNCpos-fVtxBeam;
	n_mom.SetMag(NCmom);
	fFLmom.SetVectM(n_mom, nMass);
	double mm = (fBeamLmom+TGT_LMOM-fFLmom).M();
	h1 = (TH1F*)rtFile-> Get("KN_MM"), h1-> Fill(mm);
      }
    }
  }
}

void HistManwMC::fill()
{
  rtFile-> cd();
  //  std::cout<<"===== HistManwMC::fill START  ====="<<std::endl;


  //  std::cout<<"===== HistManwMC::fill FINISH ====="<<std::endl;
}

//*** after this only test ***//
void HistManwMC::dump()
{
  std::cout<<" detData="<<detData<<std::endl;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==CID_T0 ){
      std::cout<<">>> Detector Hit  T0 seg"<<hit->channelID()+1<<std::endl;
      std::cout<<">>>    time : "<<hit->time()<<std::endl; 
      std::cout<<">>>    dE   : "<<hit->de()<<"[MeV/c]"<<std::endl; 
      std::cout<<">>>    hitZ : "<<hit->pos().Z()<<std::endl;
    }
  }

  for( int i=0; i<reacData->InitParticleSize(); i++ ){
    int pdgId = reacData-> InitPDG(i);
    TLorentzVector lmom = reacData-> GetInitParticle(i);
    if( pdgId==-321 ){
      std::cout<<">>> Init K- ====="<<std::endl;
      std::cout<<">>>   lmom ("<<lmom.X()<<", "<<lmom.Y()<<", "<<lmom.Z()<<", "<<lmom.T()<<")"<<std::endl;
    }
  }
}

void HistManwMC::finit()
{
  rtFile-> Write();
  confMan-> SaveParams();
  confMan-> SaveCode();
  rtFile-> Close();

}

bool HistManwMC::Clear(const bool flag)
{
  if( blMan ) blMan->Clear();
  if( cdsMan ) cdsMan->Clear();
  if( bltrackMan ) bltrackMan->Clear();
  if( cdstrackMan ) cdstrackMan->Clear();
  fT0time = DBL_MIN;
  fD5mom = DBL_MIN;
  fTrackBPC = 0;
  fNCtime = DBL_MAX;
  fNCdE = 0.0;
  fBeamLmom.SetXYZT(DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX);

  fT0pos.SetXYZ(DBL_MAX, DBL_MAX, DBL_MAX);
  fNCpos.SetXYZ(DBL_MAX, DBL_MAX, DBL_MAX);
  fVtxCDS.SetXYZ(DBL_MAX, DBL_MAX, DBL_MAX);
  fVtxBeam.SetXYZ(DBL_MAX, DBL_MAX, DBL_MAX);

  fT0_hit.clear();
  fBPD_hit.clear();
  fCDH_hit.clear();
  fBVC_hit.clear();
  fCVC_hit.clear();
  fPC_hit.clear();
  for( int i=0; i<8; i++ ) fNC_hit[i].clear();
  fBD_hit.clear();

  fTrackPim.clear();
  fTrackKm.clear();
  fTrackPip.clear();
  fTrackP.clear();
  fTrackD.clear();

  fFPID = F_Other;
  fNCbeta = DBL_MIN;
  fFLmom.SetXYZT(DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX);

  return flag;
}
