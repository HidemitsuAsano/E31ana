#include "KinFitMan.h"

KinFitMan::KinFitMan() : minuit(0), fStatus(-1)
{
}

double KinFitMan::pipiMMSA_n_spec(const TLorentzVector &beam_lmom, const TLorentzVector &tar_lmom, const TLorentzVector &pip_lmom, const TLorentzVector &pim_lmom, TLorentzVector &n_lmom, TFile *f, SimDataReader *simReader)
{
  //  std::cout<<"===== KinFitMan::pipiMMSA neutron spectator START  ====="<<std::endl;
  TH1F *h1;
  TH2F *h2;

  fBeamLmom = beam_lmom;
  fTarLmom = tar_lmom;
  fPID1 = CDS_PiPlus;
  fLmom1 = pip_lmom;
  fPID2 = CDS_PiMinus;
  fLmom2 = pim_lmom;

  TLorentzVector beam_lmom_MC;
  TLorentzVector pim_lmom_MC;
  TLorentzVector pip_lmom_MC;
  TLorentzVector spec_lmom_MC;
  TLorentzVector n_lmom_MC;
  if( simReader ){
    MCData *mcData = simReader->getMCData();
    ReactionData *reacData = simReader->getReactionData();
    if( reacData->ReactionID()==476 ){
      int sigma_ID=-999;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData->track(i);
	if( track-> parentTrackID()==0 ){
	  if( track->pdgID()==2112 ) spec_lmom_MC.SetVectM(0.001*track->momentum(), nMass);
	  if( track->pdgID()==-211 ) pim_lmom_MC.SetVectM(0.001*track->momentum(), piMass);
	  if( track->pdgID()== 321 ) beam_lmom_MC.SetVectM(-0.001*track->momentum(), kpMass);
	  if( track->pdgID()==3222 ) sigma_ID=track->trackID();
	}
      }
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData->track(i);
	if( track-> parentTrackID()==sigma_ID ){
	  if( track->pdgID()== 211 ) pip_lmom_MC.SetVectM(0.001*track->momentum(), piMass);
	  if( track->pdgID()==2112 ) n_lmom_MC.SetVectM(0.001*track->momentum(), nMass);
	}
      }
    }
    else if( reacData->ReactionID()==536 ){
      int sigma_ID=-999;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData->track(i);
	if( track-> parentTrackID()==0 ){
	  if( track->pdgID()==2112 ) spec_lmom_MC.SetVectM(0.001*track->momentum(), nMass);
	  if( track->pdgID()== 211 ) pip_lmom_MC.SetVectM(0.001*track->momentum(), piMass);
	  if( track->pdgID()== 321 ) beam_lmom_MC.SetVectM(-0.001*track->momentum(), kpMass);
	  if( track->pdgID()==3112 ) sigma_ID=track->trackID();
	}
      }
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData->track(i);
	if( track-> parentTrackID()==sigma_ID ){
	  if( track->pdgID()==-211 ) pim_lmom_MC.SetVectM(0.001*track->momentum(), piMass);
	  if( track->pdgID()==2112 ) n_lmom_MC.SetVectM(0.001*track->momentum(), nMass);
	}
      }
    }
    //    else {
      //      std::cout<<"!!! MCData not found Reaction : "<<reacData-> ReactionID()<<" !!!"<<std::endl;
    //    }
  }

  TLorentzVector missing_lmom = fBeamLmom+fTarLmom-pim_lmom-pip_lmom;
  TLorentzVector missing_lmom_MC = missing_lmom;
  TVector3 LabToCM = -(fBeamLmom+fTarLmom).BoostVector();
  missing_lmom_MC.Boost(LabToCM);
  double pn_mass = missing_lmom.M();
  if( pn_mass<pMass+nMass ) pn_mass=pMass+nMass;

  double mom_cm = sqrt((pn_mass+pMass+nMass)*(pn_mass-pMass+nMass)*(pn_mass+pMass-nMass)*(pn_mass-pMass-nMass))/(2*pn_mass);
  double beta = missing_lmom.Vect().Mag()/missing_lmom.E();
  double gamma = missing_lmom.E()/pn_mass;
  double Ecm = sqrt(nMass*nMass+mom_cm*mom_cm);
  fMinMom = gamma*(beta*Ecm-mom_cm);
  double reac_mom = gamma*(beta*Ecm+mom_cm);
  TVector3 mom = reac_mom*missing_lmom.Vect().Unit();
  fReacLmom.SetVectM(mom, nMass);
  h1 = (TH1F*)f-> Get("MinMom_MMSA_npipi"), h1-> Fill(fMinMom);

  if( simReader ){
    h1 = (TH1F*)f-> Get("n_mom_diff_MMSA_MC"), h1-> Fill(n_lmom.Vect().Mag()-fReacLmom.Vect().Mag());
    h2 = (TH2F*)f-> Get("n_mom_2D_MMSA_MC"), h2-> Fill(n_lmom.Vect().Mag(), fReacLmom.Vect().Mag());
  }
#if 0
  if( simReader ){
    std::cout<<"> fermi mon projection : "<<spec_lmom_MC.Vect().Dot(missing_lmom.Vect())<<std::endl;
    std::cout<<"> n_mom : "<<n_lmom_MC.Vect().Mag()<<std::endl;

    if( simReader ){
      ReactionData *reacData = simReader->getReactionData();
      if( reacData->ReactionID()==536 ){
	double im = (n_lmom+pim_lmom).M();
	TLorentzVector mm_n_lmom;
	mm_n_lmom.SetVectM((beam_lmom-pim_lmom-pip_lmom).M(), nMass);
	double im_mm = (mm_n_lmom+pim_lmom).M();
	double mmsa_im = (fReacLmom+pim_lmom).M();

	std::cout<<" im diff           : "<<im-smMass<<std::endl;
	std::cout<<" missing n im diff : "<<im_mm-smMass<<std::endl;
	std::cout<<" mmsa im diff      : "<<mmsa_im-smMass<<std::endl;
      }
    }
  }
  std::cout<<"> min mom : "<<fMinMom<<std::endl;
  std::cout<<"> reaction mom : "<<fReacLmom.Vect().Mag()<<std::endl;

  std::cout<<"Please input : q (quit)"<<std::endl;
  std::string input;
  std::cin>>input;
  if( input=="q" ) exit(-1);
#endif
  fStatus =0;
  //  std::cout<<"===== KinFitMan::pipiMMSA neutron spectator FINISH ====="<<std::endl;
  return fMinMom;
}
