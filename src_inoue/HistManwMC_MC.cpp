#include "HistManwMC.h"
#include "MyParam.h"

void HistManwMC::anaMC_NC()
{
  TH1F *h1;

  DetectorHit *hit_1st=0;
  double first_time=DBL_MAX;
  int nNC_hit = 0;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit-> detectorID()==CID_NC ){
      if( hit-> pdg()==2112 || hit-> pdg()==22 ){
	if( hit-> tdc()<first_time ){
	  nNC_hit++;
	  hit_1st = hit;
	  first_time = hit-> tdc();
	}
      }
    }
  }

  if( !hit_1st ){
    fFPID=F_Other;
    return;
  }

  fNCseg = hit_1st-> channelID()+1;
  fNC_eff_hit = 0;
  HodoscopeLikeHit *hodo_hit = 0;
  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->seg()==fNCseg ) fNC_eff_hit=blMan->NC(i);
  }
  if( !fNC_eff_hit ){
    std::cout<<"  !!! HistManwwMC::anaMC_NC not find HodoscopeLikeHit !!!"<<std::endl;
    return;
  }

  fNCdE = fNC_eff_hit-> emean();
  fNCtime = fNC_eff_hit->ctmean();
  confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNCseg , fNCpos);

  double fl = (fNCpos-fVtxBeam).Mag()-2.5;
  double beam_out, beam_tof;
  ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
  fNCbeta = fl/((fNCtime-fT0time-beam_tof)*100.*Const);

  if( fBeamPID==Beam_Kaon ){
    h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee"), h1-> Fill(1./fNCbeta);
    if( hit_1st->pdg()==22 ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_gamma"), h1-> Fill(1./fNCbeta);
    h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_8MeVee_gamma_seg%d", fNCseg)), h1-> Fill(1./fNCbeta);
    if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
      h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_wtar"), h1-> Fill(1./fNCbeta);
      if( hit_1st->pdg()==22 ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_gamma_wtar"), h1-> Fill(1./fNCbeta);
      h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_8MeVee_gamma_wtar_seg%d", fNCseg)), h1-> Fill(1./fNCbeta);
      h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_gamma_wtar"), h1-> Fill(1./fNCbeta);
    }
  }
  //  std::cout<<" nNC hit : "<<nNC_hit<<"  beta="<<1./fNCbeta<<std::endl;

  if( hit_1st->pdg()==22 ){
    //    std::cout<<"===== NC gamma : "<<MyTools::generation(hit_1st->parentID(), mcData)<<std::endl;
    int status = MyTools::gammaStatus(hit_1st->trackID(), mcData);
    h1 = (TH1F*)rtFile-> Get("NC_overbeta_gamma_MC"), h1-> Fill(1./fNCbeta);

    fFPID=F_Gamma;
  }
  if( hit_1st->pdg()==2112 ){
    //    std::cout<<"===== NC Neutrion : "<<MyTools::generation(hit_1st->parentID(), mcData)<<std::endl;
    int status = MyTools::nStatus(hit_1st->trackID(), mcData);
    h1 = (TH1F*)rtFile-> Get("NC_overbeta_n_MC"), h1-> Fill(1./fNCbeta);

    fFPID=F_Neutron;
    double NCmom = nMass*fNCbeta/sqrt(1-fNCbeta*fNCbeta);
    TVector3 n_mom = fNCpos-fVtxBeam;
    n_mom.SetMag(NCmom);
    fFLmom.SetVectM(n_mom, nMass);
  }
  //  std::cout<<std::endl;
}

DetectorHit *HistManwMC::getDetectorHit(const int &cid, const int &seg)
{
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(i);
    if( hit->detectorID()==cid && hit->channelID()+1==seg ) return hit;
  }
  std::cout<<"  !!! HistManwMC::getDetectorHit not found !!!"<<std::endl;
  return 0;
}
