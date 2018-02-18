// MyAnalysisBL.cpp

#include "MyAnalysisBL.h"
//#define DEBUG

MyAnalysisBL::MyAnalysisBL(TFile* rtFile, ConfMan* conf, bool simflag)
{
	D5 = new BeamSpectrometer(conf);
	Initialize(rtFile, conf);
	SIMULATION = simflag;
	Clear();
}

MyAnalysisBL::~MyAnalysisBL()
{
	if(D5) delete D5;
}

bool MyAnalysisBL::Clear()
{
	T0Timing = -9999.0;
	T0Segment = -1;
	BHDSegment = -1;
	BeamPID = Beam_Other;
	BeamMom = -9999.0; 
	BeamChi2 = -9999.0;
	BeamTOF = -9999.0;
	trackblc1 = 0;
	trackblc2 = 0;
	trackbpc  = 0;
	D5 -> Clear();
	FIDUCIAL = false;

	BLC1Chi2 = -1;
	BLC2Chi2 = -1;
	BPCChi2 = -1;

	BLC1Prob = -1;
	BLC2Prob = -1;
	BPCProb = -1;

	return true;
}

bool MyAnalysisBL::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, Particle* particle, double _mom)
{
	DetectorList *dlist=DetectorList::GetInstance();
	BLDCFittingParamMan* BLDCParam = conf->GetBLDCFittingParamManager();

	// ############# //
	// Trigger check //
	// ############# //
	/* ## Select K/f trigger ## */
	if(!header->IsTrig(Trig_Kaon)) return false;
	//if(!header->IsTrig(Trig_KCDH1)&&!header->IsTrig(Trig_PiCDH1)) return false;
	//if(!header->IsTrig(Trig_KCDH1)) return false;
	//if(!header->IsTrig(Trig_Kf)) return false;
	//if(!header->IsTrig(Trig_Kf)||!header->IsTrig(Trig_1stMix)) return false;
	//if(!header->IsTrig(Trig_KCDH1)||!header->IsTrig(Trig_1stMix)) return false;

	/* ## All events ## */
	FillHist("Beam_Count",0);

	// ############# //
	// T0 timing set //
	// ############# //
	{
		T0Timing = -9999.; 
		int MulT0 = 0;
		for(int i=0; i<blMan->nT0(); i++){
			if(blMan->T0(i)->CheckRange()){
				MulT0++;
				T0Timing = blMan->T0(i)->ctmean();
				T0Segment = blMan->T0(i)->seg();
			}
		}
		FillHist("T0_Mult",MulT0);
		/* ## T0 1hit request ## */
		if(MulT0!=1) return false;
		FillHist("Beam_Count",1);
		FillHist("T0_MultwCut",MulT0);
	}

	// ###################### //
	// TOF between BHD and T0 //
	// ###################### //
	{
		int MulBHD = 0;
		int seg[20]; for(int i=0; i<20; i++) {seg[i] = -1;}
		double tof[20]; for(int i=0; i<20; i++) {tof[i] = -9999.0;}
		for(int i=0; i<blMan->nBHD(); i++){
			if(blMan->BHD(i)->CheckRange()){
				if(blMan->BHD(i)->seg()==5){
					FillHist("T0BHD05_TOF",T0Timing-blMan->BHD(i)->ctmean());
				}
				if(blMan->BHD(i)->seg()==16){
					FillHist("T0BHD16_TOF",T0Timing-blMan->BHD(i)->ctmean());
				}
				if(blMan->BHD(i)->seg()<6||15<blMan->BHD(i)->seg()) continue;
				MulBHD++;
				seg[MulBHD-1] = blMan->BHD(i)->seg();
				tof[MulBHD-1] = T0Timing - blMan->BHD(i)->ctmean();;
				FillHist("T0BHD_TOF",seg[MulBHD-1],tof[MulBHD-1]);
				FillHist(Form("Beam_TOF"), tof[MulBHD-1]);
				FillHist(Form("Beam_TOF_Wide"), tof[MulBHD-1]);
			}
		}
		FillHist("BHD_Mult",MulBHD);
		if(MulBHD==1){
			FillHist(Form("Beam_TOF_Single"), tof[MulBHD-1]);
		}
		/* ## TOF kaon request ## */
		double tofk[2] = {27.9637,29.4648};
		//double tofk[2] = {27.9637,30.4648};
		//double tofpi[2] = {25.0, 27.0};
		bool TOFflag=false;
		for(int i=0; i<MulBHD; i++){
			if(tofk[0]<tof[i]&&tof[i]<tofk[1]) {
				TOFflag=true;
				BeamTOF = tof[i];
				BHDSegment = seg[i];
				TVector3 gpos;
				conf->GetGeomMapManager()->GetGPos(CID_BHD,BHDSegment,gpos);
				BHDX = gpos.X();
				break;
			}
		}
		if(!TOFflag) return false;
		FillHist("Beam_Count",2);
		FillHist("BHD_MultwCut",MulBHD);
		FillHist(Form("Beam_TOFwCut"), BeamTOF);
	}
	// #################### //
	// DC tracking          //
	// #################### //
	// ==== First Tracking of BLDC ==== //
	int ntrack_first[3] = {0,0,0};
#ifdef DEBUG
	std::cout << "First tracking" << std::endl;
#endif
	bltrackMan->DoTracking(blMan,conf,true,true);
	// === Check BLC1 === //
#ifdef DEBUG
	std::cout << "-----  BLC1  -----" << std::endl;
	std::cout << "ntrack = " << bltrackMan->ntrackBLC1() << std::endl;
#endif
	ntrack_first[0] = bltrackMan->ntrackBLC1();
	double maxresl = 0.000; int maxlayer = 0, maxcid = -1;
	for(int i=0; i<bltrackMan->ntrackBLC1(); i++){
		LocalTrack* track = bltrackMan->trackBLC1(i);
		for(int ihit=0; ihit<track->nhit(); ihit++){
			ChamberLikeHit* hit = track->hit(ihit);
			int layer = hit->layer();
			int cid = hit->cid();
			double resl = hit->resl();
			if(TMath::Abs(maxresl)<TMath::Abs(resl)){ maxresl=resl; maxlayer=layer; maxcid=cid; }
		}
		if(track->nhit(0)>4&&track->nhit(1)>4&&track->chi2all()<50){
			if(maxcid==CID_BLC1a){BLDCParam->SetLayerDeadBLC1a(maxlayer-1,true);}
			if(maxcid==CID_BLC1b){BLDCParam->SetLayerDeadBLC1b(maxlayer-1,true);}
		}
	}
#ifdef DEBUG
	std::cout << "------------------" << std::endl;
#endif
	// === Check BLC2 === //
#ifdef DEBUG
	std::cout << "-----  BLC2  -----" << std::endl;
	std::cout << "ntrack = " << bltrackMan->ntrackBLC2() << std::endl;
#endif
	ntrack_first[1] = bltrackMan->ntrackBLC2();
	maxresl = 0.000; maxlayer = 0; maxcid = -1;
	for(int i=0; i<bltrackMan->ntrackBLC2(); i++){
		LocalTrack* track = bltrackMan->trackBLC2(i);
		for(int ihit=0; ihit<track->nhit(); ihit++){
			ChamberLikeHit* hit = track->hit(ihit);
			int layer = hit->layer();
			int cid = hit->cid();
			double resl = hit->resl();
			if(TMath::Abs(maxresl)<TMath::Abs(resl)){ maxresl=resl; maxlayer=layer; maxcid=cid; }
		}
		if(track->nhit(0)>4&&track->nhit(1)>4&&track->chi2all()<50){
			if(maxcid==CID_BLC2a){BLDCParam->SetLayerDeadBLC2a(maxlayer-1,true);}
			if(maxcid==CID_BLC2b){BLDCParam->SetLayerDeadBLC2b(maxlayer-1,true);}
		}
	}
#ifdef DEBUG
	std::cout << "------------------" << std::endl;
#endif
	// === Check BPC === //
#ifdef DEBUG
	std::cout << "-----  BPC  -----" << std::endl;
	std::cout << "ntrack = " << bltrackMan->ntrackBPC() << std::endl;
#endif
	ntrack_first[2] = bltrackMan->ntrackBPC();
	maxresl = 0.000; maxlayer = 0; maxcid = -1;
	for(int i=0; i<bltrackMan->ntrackBPC(); i++){
		LocalTrack* track = bltrackMan->trackBPC(i);
		for(int ihit=0; ihit<track->nhit(); ihit++){
			ChamberLikeHit* hit = track->hit(ihit);
			int layer = hit->layer();
			int cid = hit->cid();
			double resl = hit->resl();
			if(TMath::Abs(maxresl)<TMath::Abs(resl)){ maxresl=resl; maxlayer=layer; maxcid=cid; }
		}
		if(track->nhit(0)>3&&track->nhit(1)>3&&track->chi2all()<50){BLDCParam->SetLayerDeadBPC(maxlayer-1,true);}
	}
#ifdef DEBUG
	std::cout << "------------------" << std::endl;
#endif
#ifdef DEBUG
	std::cout << "Second tracking" << std::endl;
#endif
	bltrackMan->DoTracking(blMan,conf,true,true);

	{
		const int ndc = 3;
		int dccid[ndc] = {CID_BLC1,CID_BLC2,CID_BPC};
		int nhitdc[ndc] = {0,0,0};
		int trackid[ndc][2] = {{-1,-1},{-1,-1},{-1,-1}};
		double ll[ndc] = {-30.0,-30.0,-30.0};
		double ul[ndc] = {100.0,100.0,100.0};
		double trigll[ndc] = {-10.0,-10.0,-10.0};
		double trigul[ndc] = { 10.0, 10.0, 10.0};
		double chiul[ndc] = {10.0, 10.0, 10.0};
		double probll[ndc] = {0.02, 0.02, 0.02};
		int idc;

		// ==== BLC1 and BLC2 analysis ==== // 
		for(idc=0;idc<2;idc++){
			const char* dcname = dlist->GetName(dccid[idc]).data();
#ifdef DEBUG
			std::cout << dcname << " ===== || # of track = " << bltrackMan->ntrackBLDC(dccid[idc]) << std::endl;
#endif
			FillHist(Form("%s_nTrack",dcname),bltrackMan->ntrackBLDC(dccid[idc]));
			for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
				LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
				for(int ihit=0; ihit<track->nhit(); ihit++){
					ChamberLikeHit* hit = track->hit(ihit);
					int layer = hit->layer();
					double resl = hit->resl();
					if(i<=ntrack_first[idc])FillHist(Form("%sl%d_Residual_inTrack",dcname,layer),resl,1);
				}

				double tracktime = track -> GetTrackTime();
				double chi2 = track -> chi2all();
				int dof = track->dofxz()+track->dofyz();
				double prob = TMath::Prob(chi2*(double)dof,dof);
#ifdef DEBUG
				std::cout << "track # = " << i <<  " --- " << chi2 << std::endl;
#endif
				if(i<=ntrack_first[idc])FillHist(Form("%s_TrackTime",dcname),tracktime);
				if(ll[idc]<tracktime&&tracktime<ul[idc]&&chi2<50){
					if(i<=ntrack_first[idc])FillHist(Form("%s_Chi2",dcname),chi2);
					if(i<=ntrack_first[idc])FillHist(Form("%s_Prob",dcname),prob);
					nhitdc[idc]++;
					if(nhitdc[idc]<3){
						trackid[idc][nhitdc[idc]-1] = i;
					}
				}
			}
			FillHist(Form("%s_nGoodTrack",dcname),nhitdc[idc]);
		}
		for(int i=0; i<8; i++){
			BLDCParam->SetLayerDeadBLC1a(i,false);
			BLDCParam->SetLayerDeadBLC1b(i,false);
			BLDCParam->SetLayerDeadBLC2a(i,false);
			BLDCParam->SetLayerDeadBLC2b(i,false);
		}

		bool BLCGoodflag = false;
		/* nBLC1=2 & nBLC2=2 */
		if(nhitdc[0]==2&&nhitdc[1]==2){
			LocalTrack* blc1_1 = bltrackMan->trackBLC1(trackid[0][0]); LocalTrack* blc1_2 = bltrackMan->trackBLC1(trackid[0][1]);
			LocalTrack* blc2_1 = bltrackMan->trackBLC2(trackid[1][0]); LocalTrack* blc2_2 = bltrackMan->trackBLC2(trackid[1][1]);
			if(blc1_1->chi2all()<blc1_2->chi2all()){ trackblc1 = blc1_1; } else{ trackblc1 = blc1_2; }
			if(blc2_1->chi2all()<blc2_2->chi2all()){ trackblc2 = blc2_1; } else{ trackblc2 = blc2_2; }
			double blc1time = trackblc1->GetTrackTime();
			double blc2time = trackblc2->GetTrackTime();
			BLC1Chi2 = trackblc1->chi2all();
			BLC1Prob = TMath::Prob(trackblc1->chi2all()*(double)(trackblc1->dofxz()+trackblc1->dofyz()),trackblc1->dofxz()+trackblc1->dofyz());
			BLC2Chi2 = trackblc2->chi2all();
			BLC2Prob = TMath::Prob(trackblc2->chi2all()*(double)(trackblc2->dofxz()+trackblc2->dofyz()),trackblc2->dofxz()+trackblc2->dofyz());
			FillHist("BLC1_Chi21_vs_Chi22",blc1_1->chi2all(),blc1_2->chi2all());
			FillHist("BLC2_Chi21_vs_Chi22",blc2_1->chi2all(),blc2_2->chi2all());
			if(trigll[0]<blc1time&&blc1time<trigul[0]
					&& trigll[1]<blc2time&&blc2time<trigul[1])
				BLCGoodflag = true;
		}
		/* nBLC1=1 & nBLC2=2 */
		if(nhitdc[0]==1&&trackid[0][0]<ntrack_first[0]&&nhitdc[1]==2){
			trackblc1 = bltrackMan->trackBLC1(trackid[0][0]);
			LocalTrack* blc2_1 = bltrackMan->trackBLC2(trackid[1][0]); LocalTrack* blc2_2 = bltrackMan->trackBLC2(trackid[1][1]);
			if(blc2_1->chi2all()<blc2_2->chi2all()){ trackblc2 = blc2_1; } else{ trackblc2 = blc2_2; }
			double blc1time = trackblc1->GetTrackTime();
			double blc2time = trackblc2->GetTrackTime();
			BLC1Chi2 = trackblc1->chi2all();
			BLC1Prob = TMath::Prob(trackblc1->chi2all()*(double)(trackblc1->dofxz()+trackblc1->dofyz()),trackblc1->dofxz()+trackblc1->dofyz());
			BLC2Chi2 = trackblc2->chi2all();
			BLC2Prob = TMath::Prob(trackblc2->chi2all()*(double)(trackblc2->dofxz()+trackblc2->dofyz()),trackblc2->dofxz()+trackblc2->dofyz());
			FillHist("BLC2_Chi21_vs_Chi22",blc2_1->chi2all(),blc2_2->chi2all());
			if(trigll[0]<blc1time&&blc1time<trigul[0]
					&& trigll[1]<blc2time&&blc2time<trigul[1])
				BLCGoodflag = true;
		}
		/* nBLC1=2 & nBLC2=1 */
		if(nhitdc[0]==2&&nhitdc[1]==1&&trackid[1][0]<ntrack_first[1]){
			LocalTrack* blc1_1 = bltrackMan->trackBLC1(trackid[0][0]); LocalTrack* blc1_2 = bltrackMan->trackBLC1(trackid[0][1]);
			trackblc2 = bltrackMan->trackBLC2(trackid[1][0]);
			if(blc1_1->chi2all()<blc1_2->chi2all()){ trackblc1 = blc1_1; } else{ trackblc1 = blc1_2; }
			double blc1time = trackblc1->GetTrackTime();
			double blc2time = trackblc2->GetTrackTime();
			BLC1Chi2 = trackblc1->chi2all();
			BLC1Prob = TMath::Prob(trackblc1->chi2all()*(double)(trackblc1->dofxz()+trackblc1->dofyz()),trackblc1->dofxz()+trackblc1->dofyz());
			BLC2Chi2 = trackblc2->chi2all();
			BLC2Prob = TMath::Prob(trackblc2->chi2all()*(double)(trackblc2->dofxz()+trackblc2->dofyz()),trackblc2->dofxz()+trackblc2->dofyz());
			FillHist("BLC1_Chi21_vs_Chi22",blc1_1->chi2all(),blc1_2->chi2all());
			if(trigll[0]<blc1time&&blc1time<trigul[0]
					&& trigll[1]<blc2time&&blc2time<trigul[1])
				BLCGoodflag = true;
		}
		/* nBLC1=1 & nBLC2=1 */
		if(nhitdc[0]==1&&trackid[0][0]<ntrack_first[0]&&nhitdc[1]==1&&trackid[1][0]<ntrack_first[1]){
			trackblc1 = bltrackMan->trackBLC1(trackid[0][0]);
			trackblc2 = bltrackMan->trackBLC2(trackid[1][0]);
			double blc1time = trackblc1->GetTrackTime();
			double blc2time = trackblc2->GetTrackTime();
			BLC1Chi2 = trackblc1->chi2all();
			BLC1Prob = TMath::Prob(trackblc1->chi2all()*(double)(trackblc1->dofxz()+trackblc1->dofyz()),trackblc1->dofxz()+trackblc1->dofyz());
			BLC2Chi2 = trackblc2->chi2all();
			BLC2Prob = TMath::Prob(trackblc2->chi2all()*(double)(trackblc2->dofxz()+trackblc2->dofyz()),trackblc2->dofxz()+trackblc2->dofyz());
			if(trigll[0]<blc1time&&blc1time<trigul[0]
					&& trigll[1]<blc2time&&blc2time<trigul[1])
				BLCGoodflag = true;
		}
		/* ## BLC1/2 1 good track request ## */
		if(!BLCGoodflag) return false;
		FillHist("Beam_Count",3);
		// ==== Beam momentum analysis ==== // 
		D5->TMinuitFit(trackblc1,trackblc2,conf);
		BeamChi2 = D5->chisquare();
		int beam_dof = D5->dof();
		double beam_prob = TMath::Prob(BeamChi2*(double)beam_dof,beam_dof);
		BeamMom = D5->mom();
		FillHist(Form("Beam_Momentum"), BeamMom);
		FillHist(Form("Beam_MomentumvsTOF"), BeamMom,BeamTOF);
		FillHist(Form("Beam_Chi2"), BeamChi2);
		FillHist(Form("Beam_Prob"), beam_prob);
		FillHist(Form("BHD_SegmentvsBeamMomentum"),BHDSegment,BeamMom); 
		/* ## Beam chisquare < 20 request ## */
		//if(!(BeamChi2<20)) return false;
		FillHist("Beam_Count",4);
		FillHist(Form("Beam_MomentumwCut"), BeamMom);
		FillHist(Form("Beam_MomentumvsTOFwCut"), BeamMom,BeamTOF);
		FillHist(Form("Beam_Chi2wCut"), BeamChi2);
		FillHist(Form("Beam_ProbwCut"), beam_prob);
		FillHist(Form("BHD_SegmentvsBeamMomentumwCut"),BHDSegment,BeamMom); 

		// ==== BPC analysis ==== // 
		idc=2;
		const char* dcname = dlist->GetName(dccid[idc]).data();
#ifdef DEBUG
		std::cout << dcname << " ===== || # of track = " << bltrackMan->ntrackBLDC(dccid[idc]) << std::endl;
#endif

		/* Second tracking */
		FillHist(Form("%s_nTrack",dcname),bltrackMan->ntrackBLDC(dccid[idc]));
		for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
			LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
			for(int ihit=0; ihit<track->nhit(); ihit++){
				ChamberLikeHit* hit = track->hit(ihit);
				int layer = hit->layer();
				double resl = hit->resl();
				if(i<=ntrack_first[idc])FillHist(Form("%sl%d_Residual_inTrack",dcname,layer),resl,1);
			}
			double tracktime = track -> GetTrackTime();
			double chi2 = track -> chi2all();
			int dof = track->dofxz()+track->dofyz();
			double prob = TMath::Prob(chi2*(double)dof,dof);
			if(i<=ntrack_first[idc])FillHist(Form("%s_TrackTime",dcname),tracktime);
			if(ll[idc]<tracktime&&tracktime<ul[idc]&&chi2<50){
				if(i<=ntrack_first[idc])FillHist(Form("%s_Chi2",dcname),chi2);
				if(i<=ntrack_first[idc])FillHist(Form("%s_Prob",dcname),prob);
				nhitdc[idc]++;
				if(nhitdc[idc]<3){
					trackid[idc][nhitdc[idc]-1] = i;
				}
			}
		}
		FillHist(Form("%s_nGoodTrack",dcname),nhitdc[idc]);
		for(int i=0; i<8; i++){
			BLDCParam->SetLayerDeadBPC(i,false);
		}

		bool BPCGoodflag = false;
		// BPC 1 track
		if(nhitdc[idc]==2){
			LocalTrack* bpc_1 = bltrackMan->trackBLDC(dccid[idc],trackid[idc][0]);
			LocalTrack* bpc_2 = bltrackMan->trackBLDC(dccid[idc],trackid[idc][1]);
			//if(bpc_1->chi2all()<bpc_2->chi2all()){ trackbpc = bpc_1; } else{ trackbpc = bpc_2; }
			trackbpc = bpc_1; /* Not use re-tracking */
			double tracktime = trackbpc -> GetTrackTime();
			BPCChi2 = trackbpc->chi2all();
			BPCProb = TMath::Prob(trackbpc->chi2all()*(double)(trackbpc->dofxz()+trackbpc->dofyz()),trackbpc->dofxz()+trackbpc->dofyz());
			FillHist("BPC_Chi21_vs_Chi22",bpc_1->chi2all(),bpc_2->chi2all());
			if(trigll[idc]<tracktime&&tracktime<trigul[idc])
				BPCGoodflag = true;
		}
		else if(nhitdc[idc]==1&&trackid[idc][0]<ntrack_first[idc]){
			trackbpc = bltrackMan->trackBLDC(dccid[idc],trackid[idc][0]);
			double tracktime = trackbpc -> GetTrackTime();
			BPCChi2 = trackbpc->chi2all();
			BPCProb = TMath::Prob(trackbpc->chi2all()*(double)(trackbpc->dofxz()+trackbpc->dofyz()),trackbpc->dofxz()+trackbpc->dofyz());
			if(trigll[idc]<tracktime&&tracktime<trigul[idc])
				BPCGoodflag = true;
		}
		/* ## BPC 1 good track request ## */
		if(!BPCGoodflag) return false;
		FillHist("Beam_Count",5);
		TVector3 gpos;
		conf->GetGeomMapManager()->GetGPos(CID_T0,T0Segment,gpos);
		T0Position = trackbpc->GetPosatZ(gpos.Z()-0.5);

		BeamP = trackbpc->GetMomDir();
		BeamP.SetMag(BeamMom);
		BeamPID = Beam_Kaon;
		BeamMass = kpMass;
		BeamL.SetVectM( BeamP, BeamMass );
	}
	// ==== BLC2 and BPC maching ==== // 
	{
		double xll  = -0.75; double xul  = 0.75;
		double yll  = -0.75; double yul  = 0.75;
		double dxll = -0.02; double dxul = 0.02;
		double dyll = -0.02; double dyul = 0.02;
		TVector3 posblc2 = trackblc2->GetPosatZ(-75.0);
		TVector3 posbpc  = trackbpc->GetPosatZ(-75.0);
		TVector3 dirblc2 = trackblc2->GetMomDir();
		TVector3 dirbpc  = trackbpc->GetMomDir();
		TVector3 dist = posblc2-posbpc;
		double dxdz = dirblc2.X()/dirblc2.Z() - dirbpc.X()/dirbpc.Z();
		double dydz = dirblc2.Y()/dirblc2.Z() - dirbpc.Y()/dirbpc.Z();
		FillHist("BLC2_75Position",posblc2.X(),posblc2.Y());
		FillHist("BPC_75Position",posbpc.X(),posbpc.Y());
		FillHist("BLC2_Direction",dirblc2.X()/dirblc2.Z(),dirblc2.Y()/dirblc2.Z());
		FillHist("BPC_Direction",dirbpc.X()/dirbpc.Z(),dirbpc.Y()/dirbpc.Z());
		FillHist("BLC2BPC_Distance",dist.X(),dist.Y());
		FillHist("BLC2BPC_Distance",dist.X(),dist.Y());
		FillHist("BLC2BPC_Angle",dxdz,dydz);
		TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
		TVector3 posffbpc  = trackbpc->GetPosatZ(0.0);
		FillHist("BLC2_FFPosition",posffblc2.X(),posffblc2.Y());
		FillHist("BPC_FFPosition",posffbpc.X(),posffbpc.Y());
		/* ## BLC2 and BPC track maching ## */
		{
			if( !((xll<dist.X()&&dist.X()<xul) && (yll<dist.Y()&&dist.Y()<yul)
						&& (dxll<dxdz&&dxdz<dxul) && (dyll<dydz&&dydz<dyul)) )
				return false;
			FillHist("Beam_Count",6);
			FillHist("BLC2BPC_DistancewCut",dist.X(),dist.Y());
			FillHist("BLC2BPC_AnglewCut",dxdz,dydz);
		}
	}
	// ==== Beam position at FF ==== // 
	{
		TVector3 posff  = trackbpc->GetPosatZ(0.0);
		TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
		double dxdz = BeamP.X()/BeamP.Z();
		double dydz = BeamP.Y()/BeamP.Z();
		FillHist("Beam_FFPosition",posff.X(),posff.Y());
		FillHist("Beam_Angle",dxdz,dydz);
		/* ## On target selection ## */
		if(GeomTools::GetID(posff)==CID_Fiducial){
			FIDUCIAL = true;
			FillHist("Beam_Count",7);
			FillHist("Beam_FFPositionwCut",posff.X(),posff.Y());
			FillHist("Beam_AnglewCut",dxdz,dydz);
		}
	}

	pBeam* beam = new pBeam();
	beam->SetPID(Beam_Kaon);
	beam->SetT0seg(T0Segment);
	beam->SetT0Time(T0Timing);
	beam->SetT0Pos(T0Position);
	beam->SetBHDseg(BHDSegment);
	beam->SetBHDT0TOF(BeamTOF);
	beam->SetBHDX(BHDX);
	beam->SetBPCPos(trackbpc->GetPosatZ(0.0));
	beam->SetBPCDir(trackbpc->GetMomDir());
	TVector3 dir=trackblc2->GetMomDir();
	beam->SetBLCPos(trackblc2->GetPosatZ(0.0));
	beam->SetBLCDir(dir);
	beam->SetMomentum(BeamMom);
	beam->SetBeamChi2(BeamChi2);
	beam->SetBLC1Chi2(BLC1Chi2);
	beam->SetBLC2Chi2(BLC2Chi2);
	beam->SetBPCChi2(BPCChi2);
	beam->SetBLC1Prob(BLC1Prob);
	beam->SetBLC2Prob(BLC2Prob);
	beam->SetBPCProb(BPCProb);

	particle->AddBeam(*beam);

	return true;
}

bool MyAnalysisBL::DoAnalysisSim(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, Particle* particle, double _mom)
{
	if(!SIMULATION) {return false;}

	DetectorList *dlist=DetectorList::GetInstance();

	/* ## All events ## */
	FillHist("Beam_Count",0);

	// ############# //
	// T0 timing set //
	// ############# //
	{
		T0Timing = -9999.; 
		int MulT0 = 0;
		for(int i=0; i<blMan->nT0(); i++){
			if(blMan->T0(i)->CheckRange()){
				MulT0++;
				T0Timing = blMan->T0(i)->ctmean();
				T0Segment = blMan->T0(i)->seg();
			}
		}
		FillHist("T0_Mult",MulT0);
		/* ## T0 1hit request ## */
		if(MulT0!=1) return false;
		FillHist("Beam_Count",1);
		FillHist("T0_MultwCut",MulT0);
	}

	// #################### //
	// DC tracking          //
	// #################### //
	{
		bltrackMan->DoTracking(blMan,conf,true,true);
		const int ndc = 3;
		int dccid[ndc] = {CID_BLC1,CID_BLC2,CID_BPC};
		int nhitdc[ndc] = {0,0,0};
		int trackid[ndc] = {-1,-1,-1};
		double ll[ndc] = {-30.0,-30.0,-30.0};
		double ul[ndc] = {100.0,100.0,100.0};
		double trigll[ndc] = {-10.0,-10.0,-10.0};
		double trigul[ndc] = { 10.0, 10.0, 10.0};
		double chiul[ndc] = {10.0, 10.0, 10.0};
		double probll[ndc] = {0.02, 0.02, 0.02};
		int idc;
		// ==== BLC1 and BLC2 analysis ==== // 
		for(idc=1;idc<2;idc++){
			const char* dcname = dlist->GetName(dccid[idc]).data();
			FillHist(Form("%s_nTrack",dcname),bltrackMan->ntrackBLDC(dccid[idc]));
			for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
				if(i>=1) break;
				LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
				for(int ihit=0; ihit<track->nhit(); ihit++){
					ChamberLikeHit* hit = track->hit(ihit);
					int layer = hit->layer();
					double resl = hit->resl();
					FillHist(Form("%sl%d_Residual_inTrack",dcname,layer),resl,1);
				}
				double chi2 = track -> chi2all();
				int dof = track->dofxz()+track->dofyz();
				double prob = TMath::Prob(chi2*(double)dof,dof);
				double tracktime = track -> GetTrackTime();
				if(idc==0){
					BLC1Chi2 = chi2;
					BLC1Prob = prob;
				}
				else if(idc==1){
					BLC2Chi2 = chi2;
					BLC2Prob = prob;
				}
				FillHist(Form("%s_TrackTime",dcname),tracktime);
				FillHist(Form("%s_Chi2",dcname),chi2);
				FillHist(Form("%s_Prob",dcname),prob);
				if(1){
					nhitdc[idc]++;
					trackid[idc] = i;
					FillHist(Form("%s_TrackTimewCut",dcname),tracktime);
					FillHist(Form("%s_Chi2wCut",dcname),chi2);
					FillHist(Form("%s_ProbwCut",dcname),prob);
				}
			}
			FillHist(Form("%s_nGoodTrack",dcname),nhitdc[idc]);
		}
		bool BLCGoodflag = false;
		if(nhitdc[1]==1){
			BLCGoodflag = true;
		}
		/* ## BLC1/2 1 good track request ## */
		if(!BLCGoodflag) return false;
		FillHist("Beam_Count",3);
		trackblc2 = bltrackMan->trackBLC2(trackid[1]);
		// ==== Beam momentum analysis ==== // 
		BeamChi2 = 0.0;
		BeamMom = _mom;
		FillHist(Form("Beam_Momentum"), BeamMom);
		FillHist(Form("Beam_MomentumvsTOF"), BeamMom,BeamTOF);
		FillHist(Form("Beam_Chi2"), BeamChi2);
		/* ## Beam chisquare < 20 request ## */
		//if(!(BeamChi2<20)) return false;
		FillHist("Beam_Count",4);
		FillHist(Form("Beam_MomentumwCut"), BeamMom);
		FillHist(Form("Beam_MomentumvsTOFwCut"), BeamMom,BeamTOF);
		FillHist(Form("Beam_Chi2wCut"), BeamChi2);

		// ==== BPC analysis ==== // 
		idc=2;
		const char* dcname = dlist->GetName(dccid[idc]).data();
		FillHist(Form("%s_nTrack",dcname),bltrackMan->ntrackBLDC(dccid[idc]));
		for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
			if(i>=1) break;
			LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
			for(int ihit=0; ihit<track->nhit(); ihit++){
				ChamberLikeHit* hit = track->hit(ihit);
				int layer = hit->layer();
				double resl = hit->resl();
				FillHist(Form("%sl%d_Residual_inTrack",dcname,layer),resl,1);
			}
			double chi2 = track -> chi2all();
			int dof = track->dofxz()+track->dofyz();
			double prob = TMath::Prob(chi2*(double)dof,dof);
			double tracktime = track -> GetTrackTime();
			if(idc==2){
				BPCChi2 = chi2;
				BPCProb = prob;
			}
			FillHist(Form("%s_TrackTime",dcname),tracktime);
			FillHist(Form("%s_Chi2",dcname),chi2);
			FillHist(Form("%s_Prob",dcname),prob);
			if(1){
				nhitdc[idc]++;
				trackid[idc] = i;
				FillHist(Form("%s_TrackTimewCut",dcname),tracktime);
				FillHist(Form("%s_Chi2wCut",dcname),chi2);
				FillHist(Form("%s_ProbwCut",dcname),prob);
			}
		}
		FillHist(Form("%s_nGoodTrack",dcname),nhitdc[idc]);
		bool BPCGoodflag = false;
		// BPC 1 track
		if(nhitdc[idc]==1){
			BPCGoodflag = true;
		}
		/* ## BPC 1 good track request ## */
		if(!BPCGoodflag) return false;
		FillHist("Beam_Count",5);
		trackbpc  = bltrackMan->trackBPC(trackid[2]);
		TVector3 gpos;
		conf->GetGeomMapManager()->GetGPos(CID_T0,T0Segment,gpos);
		T0Position = trackbpc->GetPosatZ(gpos.Z()-0.5);

		BeamP = trackbpc->GetMomDir();
		BeamP.SetMag(BeamMom);
		BeamPID = Beam_Kaon;
		BeamMass = kpMass;
		BeamL.SetVectM( BeamP, BeamMass );
	}
	// ==== BLC2 and BPC maching ==== // 
	{
		double xll  = -0.75; double xul  = 0.75;
		double yll  = -0.75; double yul  = 0.75;
		double dxll = -0.02; double dxul = 0.02;
		double dyll = -0.02; double dyul = 0.02;
		TVector3 posblc2 = trackblc2->GetPosatZ(-75.0);
		TVector3 posbpc  = trackbpc->GetPosatZ(-75.0);
		TVector3 dirblc2 = trackblc2->GetMomDir();
		TVector3 dirbpc  = trackbpc->GetMomDir();
		TVector3 dist = posblc2-posbpc;
		double dxdz = dirblc2.X()/dirblc2.Z() - dirbpc.X()/dirbpc.Z();
		double dydz = dirblc2.Y()/dirblc2.Z() - dirbpc.Y()/dirbpc.Z();
		FillHist("BLC2_75Position",posblc2.X(),posblc2.Y());
		FillHist("BPC_75Position",posbpc.X(),posbpc.Y());
		FillHist("BLC2_Direction",dirblc2.X()/dirblc2.Z(),dirblc2.Y()/dirblc2.Z());
		FillHist("BPC_Direction",dirbpc.X()/dirbpc.Z(),dirbpc.Y()/dirbpc.Z());
		FillHist("BLC2BPC_Distance",dist.X(),dist.Y());
		FillHist("BLC2BPC_Distance",dist.X(),dist.Y());
		FillHist("BLC2BPC_Angle",dxdz,dydz);
		TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
		TVector3 posffbpc  = trackbpc->GetPosatZ(0.0);
		FillHist("BLC2_FFPosition",posffblc2.X(),posffblc2.Y());
		FillHist("BPC_FFPosition",posffbpc.X(),posffbpc.Y());
		/* ## BLC2 and BPC track maching ## */
		{
			if( !((xll<dist.X()&&dist.X()<xul) && (yll<dist.Y()&&dist.Y()<yul)
						&& (dxll<dxdz&&dxdz<dxul) && (dyll<dydz&&dydz<dyul)) )
				return false;
			FillHist("Beam_Count",6);
			FillHist("BLC2BPC_DistancewCut",dist.X(),dist.Y());
			FillHist("BLC2BPC_AnglewCut",dxdz,dydz);
		}
	}
	// ==== Beam position at FF ==== // 
	{
		TVector3 posff  = trackbpc->GetPosatZ(0.0);
		TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
		double dxdz = BeamP.X()/BeamP.Z();
		double dydz = BeamP.Y()/BeamP.Z();
		FillHist("Beam_FFPosition",posff.X(),posff.Y());
		FillHist("Beam_Angle",dxdz,dydz);
		/* ## On target selection ## */
		if(GeomTools::GetID(posff)==CID_Fiducial){
			FIDUCIAL = true;
			FillHist("Beam_Count",7);
			FillHist("Beam_FFPositionwCut",posff.X(),posff.Y());
			FillHist("Beam_AnglewCut",dxdz,dydz);
		}
	}

	pBeam* beam = new pBeam();
	beam->SetPID(Beam_Kaon);
	beam->SetT0seg(T0Segment);
	beam->SetT0Time(T0Timing);
	beam->SetT0Pos(T0Position);
	beam->SetBPCPos(trackbpc->GetPosatZ(0.0));
	beam->SetBPCDir(trackbpc->GetMomDir());
	TVector3 dir=trackblc2->GetMomDir();
	beam->SetBLCPos(trackblc2->GetPosatZ(0.0));
	beam->SetBLCDir(dir);
	beam->SetMomentum(BeamMom);
	beam->SetBeamChi2(BeamChi2);
	beam->SetBLC1Chi2(BLC1Chi2);
	beam->SetBLC2Chi2(BLC2Chi2);
	beam->SetBPCChi2(BPCChi2);
	beam->SetBLC1Prob(BLC1Prob);
	beam->SetBLC2Prob(BLC2Prob);
	beam->SetBPCProb(BPCProb);

	particle->AddBeam(*beam);

	return true;
}

bool MyAnalysisBL::FillHist(TString name, double val1)
{
	TH1F* h1 = (TH1F*)gFile -> Get(name);
	if(h1){
		h1 -> Fill(val1);
		return true;
	}
	else {
		return false;
	}
}

bool MyAnalysisBL::FillHist(TString name, double val1, double val2)
{
	TH2F* h2 = (TH2F*)gFile -> Get(name);
	if(h2){
		h2 -> Fill(val1,val2);
		return true;
	}
	else {
		return false;
	}
}

bool MyAnalysisBL::Initialize(TFile* rtFile, ConfMan* confMan)
{
	rtFile -> cd();
	DetectorList *dlist=DetectorList::GetInstance();

	// Scaler
	new TH1F( "Scaler", "Scaler", 50, 0, 50 );

	// Hodoscope
	const int nhodo = 4;
	int NumOfSegments[nhodo] = {20,5,70,8};
	TString CounterName[nhodo] = {"BHD","T0","BPD","DEF"};
	std::cout << "Define Histograms for Hodoscope" << std::endl;
	for(int ihodo=0; ihodo<2; ihodo++){
		new TH1F( Form("%s_HitPat",CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
		new TH1F( Form("%s_Mult",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
		new TH1F( Form("%s_MultwCut",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
	}
	// Chamber
	const int ndc=3;
	int dccid[ndc]={CID_BLC1,CID_BLC2,CID_BPC};
	TString dcname[ndc]={"BLC1","BLC2","BPC"};
	std::cout << "Define Histgram for DC" << std::endl;
	for(int idc=0;idc<ndc;idc++){
		const int nlays= dlist->GetNlayers(dccid[idc]);
		for( int layer=1; layer<=nlays; layer++ ){
			new TH1F( Form("%sl%d_Residual_inTrack",dcname[idc].Data(),layer),Form("Residual %s-L%d;Residual (cm);Counts",dcname[idc].Data(),layer),500,-0.5,0.5);
		}
		new TH2F( Form("%s_Chi21_vs_Chi22",dcname[idc].Data()),Form("1st #chi^{2} vs. 2nd #chi^{2} %s;1st #chi^{2};2nd #chi^{2}",dcname[idc].Data()),500,0.0,50.0,500,0.0,50.0);
		new TH1F( Form("%s_nTrack",dcname[idc].Data()), Form("Num. of tracks %s;Number of Tracks;Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
		new TH1F( Form("%s_nGoodTrack",dcname[idc].Data()), Form("Num. of tracks %s;Number of Good Tracks;Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
		new TH1F( Form("%s_nBeamTrack",dcname[idc].Data()), Form("Num. of tracks %s;Number of Beam Tracks;Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
		new TH1F( Form("%s_TrackTime",dcname[idc].Data()), Form("Track time %s;Time (ns);Counts",dcname[idc].Data()), 1200, -200.0, 400.0 );
		new TH1F( Form("%s_TrackTimewCut",dcname[idc].Data()), Form("Track time %s;Time (ns);Counts",dcname[idc].Data()), 1200, -200.0, 400.0 );
		new TH1F( Form("%s_Chi2",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2}/ndf;Counts",dcname[idc].Data()), 1000, 0.0, 100.0 );
		new TH1F( Form("%s_Chi2wCut",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2}/ndf;Counts",dcname[idc].Data()), 1000, 0.0, 100.0 );
		new TH1F( Form("%s_Prob",dcname[idc].Data()), Form("probability %s;probability;Counts",dcname[idc].Data()), 1000, 0.0, 1.0 );
		new TH1F( Form("%s_ProbwCut",dcname[idc].Data()), Form("probabilith %s;probability;Counts",dcname[idc].Data()), 1000, 0.0, 1.0 );
	}
	// TOF
	std::cout << "Define Histograms for TOF" << std::endl;
	new TH2F( Form("T0%s_TOF",CounterName[0].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[0].Data(),CounterName[0].Data()), NumOfSegments[0], 1, NumOfSegments[0]+1, 4000, -50, 50 );
	new TH1F( "T0BHD05_TOF", "TOF T0-BHD05;TOF (ns);Counts",1000, 0.0, 100.0 );
	new TH1F( "T0BHD16_TOF", "TOF T0-BHD16;TOF (ns);Counts",1000, 0.0, 100.0 );
	// Beam analysis
	std::cout << "Define Histograms for Beam analysis" << std::endl;
	new TH1F( "Beam_Count", "Beam counts;Selection;Counts", 10, 0, 10 );
	new TH1F( "Beam_Momentum", "Beam momentum;Beam momentum (GeV/c);Counts",1000, 0.5, 1.5 );
	new TH1F( "Beam_MomentumwCut", "Beam momentum;Beam momentum (GeV/c);Counts",1000, 0.5, 1.5 );
	new TH1F( "Beam_TOF", "Beam TOF;TOF between BHD and T0 (ns);Counts",2000, 20.0, 40.0 );
	new TH1F( "Beam_TOF_Wide", "Beam TOF;TOF between BHD and T0 (ns);Counts",7000, -25.0, 45.0 );
	new TH1F( "Beam_TOF_Single", "Beam TOF;TOF between BHD and T0 (ns);Counts",7000, -25.0, 45.0 );
	new TH1F( "Beam_TOFwCut", "Beam TOF;TOF between BHD and T0 (ns);Counts",1000, 20.0, 40.0 );
	new TH2F( "Beam_MomentumvsTOF", "Beam Momentum vs. TOF;Beam Momentum (GeV/c);TOF between BHD and T0 (ns)",1000,0.5,1.5,1000, 20.0, 40.0 );
	new TH2F( "Beam_MomentumvsTOFwCut", "Beam Momentum vs. TOF;Beam Momentum (GeV/c);TOF between BHD and T0 (ns)",1000,0.5,1.5,1000, 20.0, 40.0 );
	new TH1F( "Beam_Chi2", "Beam mom #chi^{2};#chi^{2}/ndf;Counts",1000, 0.0, 100.0 );
	new TH1F( "Beam_Chi2wCut", "Beam mom #chi^{2};#chi^{2}/ndf;Counts",1000, 0.0, 100.0 );
	new TH1F( "Beam_Prob", "Beam mom probability;probability;Counts",1000, 0.0, 1.0 );
	new TH1F( "Beam_ProbwCut", "Beam mom probability;probability;Counts",1000, 0.0, 1.0 );
	new TH2F( "Beam_FFPosition", "Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
	new TH2F( "Beam_FFPositionwCut", "Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
	new TH2F( "Beam_CellRing1Position", "Position at CellRing1;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
	new TH2F( "Beam_CellRing2Position", "Position at CellRing2;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
	new TH2F( "Beam_Angle", "Beam Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
	new TH2F( "Beam_AnglewCut", "Beam Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
	new TH2F( "BHD_SegmentvsBeamMomentum", "BHD sement vs. Beam momentum;BHD segment;Beam momentum (GeV/c)", 20, 1, 21, 1000, 0.5, 1.5 );
	new TH2F( "BHD_SegmentvsBeamMomentumwCut", "BHD sement vs. Beam momentum;BHD segment;Beam momentum (GeV/c)", 20, 1, 21, 1000, 0.5, 1.5 );

	new TH2F( "BLC2_FFPosition", "BLC2 Track Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
	new TH2F( "BLC2_75Position", "BLC2 Track Position at z=-75.0 cm;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
	new TH2F( "BPC_FFPosition", "BPC Track Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
	new TH2F( "BPC_75Position", "BPC Track Position at z=-75.0 cm;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
	new TH2F( "BLC2_Direction", "BLC2 Track Direction;dX/dZ;dY/dZ",1000,-0.05,0.05,1000,-0.05,0.05);
	new TH2F( "BPC_Direction", "BPC Track Direction;dX/dZ;dY/dZ",1000,-0.05,0.05,1000,-0.05,0.05);
	new TH2F( "BLC2BPC_Distance", "Distance between BLC2 and BPC;X distance (cm);Y distance (cm)",1000,-5.0,5.0,1000,-5.0,5.0 );
	new TH2F( "BLC2BPC_DistancewCut", "Distance between BLC2 and BPC;X distance (cm);Y distance (cm)",1000,-5.0,5.0,1000,-5.0,5.0 );
	new TH2F( "BLC2BPC_Angle", "Angle between BLC2 and BPC;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
	new TH2F( "BLC2BPC_AnglewCut", "Angle between BLC2 and BPC;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);

	std::cout << "=== End of [MyAnalysisBL::Initialize] === " << std::endl;

	return true;
}

