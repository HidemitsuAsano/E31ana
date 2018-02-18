// MyAnalysisCheck.cpp

#include "MyAnalysisCheck.h"

MyAnalysisCheck::MyAnalysisCheck(TFile* rt, ConfMan* conf)
{
	Initialize(conf);
	CutCondition();
	Clear();
}

MyAnalysisCheck::~MyAnalysisCheck()
{
	Clear();
	rtFile->cd();
	rtFile->Write();
	rtFile->Close();
}

void MyAnalysisCheck::Clear()
{
}

bool MyAnalysisCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
	rtFile->cd();

	DetectorList *dlist=DetectorList::GetInstance();
	FillHist("CDC_nGoodTrack",cdstrackMan->nGoodTrack());
	if(particle->nBeam()!=1) return false;
	pBeam* beam = particle->beam(0);

	// ##################### //
	// CDH all hits analysis //
	// ##################### //
	int hit_seg[40]={0};
	double hit_de[40]={0};
	double hit_time[40]={0};
	double hit_pos[40]={0};
	int cdh_nhit = 0;
	for(int id=0; id<cdsMan->nCDH(); id++){
		HodoscopeLikeHit* cdh = cdsMan->CDH(id);
		if(!cdh->CheckRange()) continue;
		double time = cdh->ctmean();
		double de = cdh->emean();
		double hitpos = cdh->hitpos();
		double phi = cdh->pos().Phi();
		int seg = cdh->seg();
		hit_pos[seg-1] = hitpos;
		hit_seg[seg-1] = 1;
		hit_de[seg-1] = de;
		hit_time[seg-1] = time;
		cdh_nhit++;
	}
	if(cdh_nhit==0) return true;

#ifdef DEBUT
	std::cout << "####################################" << std::endl;
	std::cout << "------------------------------------" << std::endl;
	for(int i=0; i<36; i++){
		if(hit_seg[i]){ std::cout << "o"; }
		else{ std::cout << " "; }
	}
	std::cout << std::endl;
	std::cout << "------------------------------------" << std::endl;
#endif


	FillHist("CDH_nHit",cdh_nhit);
	int cdh_cluster[40][40] = {};
	int ncluster = 0;
	int end = 36;
	for(int i=0; i<end; i++){
		int nhit=0;
		if(!hit_seg[i]==1){continue;}
		nhit++;

		if(i==0){
			int start = 0;
			while(hit_seg[36-nhit]==1){
				start = 36-nhit;
				nhit++;
				end = start-1;
			}
			int nhit2=0;
			while(hit_seg[i+nhit2+1]==1){
				nhit2++;
			}
			nhit+=nhit2;
			cdh_cluster[ncluster][0]=nhit;
			for(int icluster=1; icluster<=nhit; icluster++){
				if(start+icluster<=36){
					cdh_cluster[ncluster][icluster]=start+icluster;
				}
				else{
					cdh_cluster[ncluster][icluster]=start+icluster-36;
				}
			}
			ncluster++;
			i+=nhit2;
		}	
		else{
			while(hit_seg[i+nhit]==1){
				nhit++;
			}
			cdh_cluster[ncluster][0]=nhit;
			for(int icluster=1; icluster<=nhit; icluster++){
				cdh_cluster[ncluster][icluster]=i+icluster;
			}
			ncluster++;
			i+=nhit;
		}
	}

#ifdef DEBUG
	std::cout << "ncluster = " << ncluster << std::endl;
# endif

	int hit_incdc_seg[40]={0};
	for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
		CDSTrack* cdc=cdstrackMan->GoodTrack(id);
		if(cdc->Chi()>50) continue;
		double cdhtime=0,dphi=0;
		int cdhseg=0;
		if(!cdc->GetCDHHit2(cdsMan,cdhseg,cdhtime,dphi)) continue;
		hit_incdc_seg[cdhseg-1]++;
		if(hit_incdc_seg[cdhseg-1]>1){return true;}
	}
	int nt[40]={0};
	for(int i=0; i<ncluster; i++){
		for(int j=1; j<=cdh_cluster[i][0]; j++){
			if(hit_incdc_seg[cdh_cluster[i][j]-1]){ nt[i]++; }
		}
	}

	double delta = 0.000000001;
	bool flag=false;

	FillHist("CDH_nCluster",ncluster);
	FillHist("CDS_nGoodTrack_vs_nHit",(double)cdstrackMan->nGoodTrack(),(double)cdh_nhit);
	FillHist("CDS_nGoodTrack_vs_nCluster",(double)cdstrackMan->nGoodTrack(),(double)ncluster);
	for(int i=0; i<ncluster; i++){
		FillHist("CDH_nHitinCluster",cdh_cluster[i][0]);
		double tmp_t[36] = {};
		double tmp_p[36] = {};
		double tmp_de[36] = {};
		for(int j=1; j<=cdh_cluster[i][0]; j++){
			tmp_t[j-1] = hit_time[cdh_cluster[i][j]-1];
			tmp_p[j-1] = hit_pos[cdh_cluster[i][j]-1];
			tmp_de[j-1] = hit_de[cdh_cluster[i][j]-1];
#ifdef DEBUG
			std::cout << tmp_t[j-1] << std::endl;
#endif
		}
		for(int j=1; j<=cdh_cluster[i][0]; j++){
			if(cdh_cluster[i][0]>j){
				if(j==1){
					if(fabs(tmp_t[j-1]-tmp_t[j])<1.0&&fabs(tmp_p[j-1]-tmp_p[j])>35.0) flag=true;
				}
				//FillHist(Form("CDHCluster_timediff_vs_posdiff_%d%d",j,j+1),fabs(tmp_t[j-1]-tmp_t[j]),fabs(tmp_p[j-1]-tmp_p[j]));
				FillHist(Form("CDHCluster_timediff_vs_posdiff_%d%d",cdh_cluster[i][j],cdh_cluster[i][j+1]),fabs(tmp_t[j-1]-tmp_t[j]),fabs(tmp_p[j-1]-tmp_p[j]));
			}
		}
	}
#ifdef DEBUG
	std::cout << "####################################" << std::endl;
#endif

	std::vector<CDHCluster> cdhcontainer;
	for(int i=0; i<ncluster; i++){
		CDHCluster* cdhcluster = new CDHCluster();
		double pre_time=0;
		double pre_pos=0;
		int nhitincluster=0;
		int ntrackincluster=0;
		for(int j=1; j<=cdh_cluster[i][0]; j++){
			if(hit_incdc_seg[cdh_cluster[i][j]-1]){ ntrackincluster++; }
			HodoscopeLikeHit* hit = 0;
			for(int id=0; id<cdsMan->nCDH(); id++){
				if((int)cdh_cluster[i][j]==(int)cdsMan->CDH(id)->seg()){
					hit=cdsMan->CDH(id);
				}
			}
			if(!hit){ continue; }
			if(nhitincluster==0){
				cdhcluster->AddHit(hit);
				nhitincluster++;
			}
			else{
				if((ntrackincluster<=1)&&(fabs(hit->ctmean()-pre_time)<2.0&&fabs(hit->hitpos()-pre_pos)<10)){
					cdhcluster->AddHit(hit);
					nhitincluster++;
				}
				else{
					cdhcluster->Calc();
					cdhcontainer.push_back(*cdhcluster);
					delete cdhcluster;
					cdhcluster = new CDHCluster();
					nhitincluster=0;
					if(ntrackincluster<=1){ ntrackincluster=0; }
					else{ ntrackincluster=1; }
					cdhcluster->AddHit(hit);
					nhitincluster++;
				}
			}
			pre_time = hit->ctmean();
			pre_pos = hit->hitpos();
		}
		if(nhitincluster){
			cdhcluster->Calc();
			cdhcontainer.push_back(*cdhcluster);
		}
	}

//#ifdef DEBUG
if(0){
	std::cout << "####################################" << std::endl;
	for(int i=0; i<36; i++){
		if(hit_seg[i]){ std::cout << "o"; }
		else{ std::cout << " "; }
	}
	std::cout << std::endl;
	for(int i=0; i<36; i++){
		if(hit_incdc_seg[i]){ std::cout << hit_incdc_seg[i]; }
		else{ std::cout << " "; }
	}
	std::cout << std::endl;
	std::cout << "-----------" << std::endl;
	for(int i=0; i<cdhcontainer.size(); i++){
		CDHCluster cdhcluster = cdhcontainer[i];
		for(int j=0; j<cdhcluster.nCDH(); j++){
			std::cout << cdhcluster.CDH(j)->seg() << ", ";
		}
		std::cout << std::endl;
	}
}
//#endif

	FillHist("CDS_nGoodTrack_vs_nCluster_Re",(double)cdstrackMan->nGoodTrack(),(double)cdhcontainer.size());
	FillHist("CDS_nCluster_vs_nCluster_Re",(double)ncluster,(double)cdhcontainer.size());
	for(int i=0; i<cdhcontainer.size(); i++){
		CDHCluster cdhcluster = cdhcontainer[i];
		FillHist("CDH_nHitinCluster_Re",cdhcluster.nCDH());
	}

	return true;





	// ###################### //
	// CDS all track analysis //
	// ###################### //
	for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
		int icut=0;
		FillHist("CDS_NumOfTrack",icut); icut++; /* Initial */
		CDSTrack* cdc=cdstrackMan->GoodTrack(id);
		CDHCluster* cdh_c = 0;
		double cdhtime=0,dis=0,dphi=0,de=0;
		int cdhseg=0;

		if(cdc->FittingLevel()!=5){
			if(!cdc->GetCDHHit2(cdsMan,cdhseg,cdhtime,dphi)) continue;
			FillHist("CDC_DPhi",dphi);

			for(int ic=0; ic<cdhcontainer.size(); ic++){
				CDHCluster cdhcluster =  cdhcontainer[ic];
				for(int ihit=0; ihit<cdhcluster.nCDH(); ihit++){
					if(cdhcluster.CDH(ihit)->seg()==cdhseg){ cdh_c = &cdhcontainer[ic]; }
				}
			}
			if(!cdh_c) continue;
			de=cdh_c->emean();

			FillHist("CDS_NumOfTrack",icut); icut++; /* Find CDH */

			TVector3 cdsvtx,beamvtx;
			double par[5];
			cdc->GetParameters(CID_CDC,par,cdsvtx);
			double mom=cdc->Momentum();
			if(!cdc->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
			dis = (beamvtx-cdsvtx).Mag();

			double beam_out=0, beam_tof=0;
			ELossTools::CalcElossBeamTGeo(beam->t0pos(),beamvtx,beam->mom(),kpMass,beam_out,beam_tof);
			TVector3 cdhvtx=cdc->CDHVertex();
			double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
			double beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
			double mass2=mom*mom*(1/(beta*beta)-1);

			double calc_beta,tofvtxcdc;
			if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
			beta = calc_beta;
#ifdef DEBUG
			std::cout << "First fitting -----" << std::endl;
			std::cout << "chi = " << cdc->Chi() << std::endl;
#endif
			int maxlayer=-1, maxwire=-1;
			double maxresl=0.0;
			for(int layer=1; layer<=15; layer++){
				for(int nhit=0; nhit<(int)cdc->nTrackHit(layer); nhit++){
					CDCHit* cdchit = cdc->hit(cdsMan,layer,nhit);
					if(TMath::Abs(maxresl)<TMath::Abs(cdchit->resl())){
						maxlayer=layer; maxwire=cdchit->wire();maxresl=cdchit->resl();
					}
				}
			}
#ifdef DEBUG
			std::cout << "max = l:" << maxlayer << ", w:" << maxwire << std::endl;
			std::cout << maxresl << std::endl;
#endif
			double tmpchi = cdc->Chi();
			cdc->Retiming(cdsMan,conf,beta,true,maxlayer,maxwire);
#ifdef DEBUG
			std::cout << "Second fitting -----" << std::endl;
			std::cout << "chi = " << cdc->Chi() << std::endl;
#endif
			double maxresl_first = maxresl;
			int maxlayer_first=maxlayer, maxwire_first=maxwire;
			maxlayer=-1; maxwire=-1;
			maxresl=0.0;
			for(int layer=1; layer<=15; layer++){
				for(int nhit=0; nhit<(int)cdc->nTrackHit(layer); nhit++){
					CDCHit* cdchit = cdc->hit(cdsMan,layer,nhit);
					if(layer==maxlayer_first&&cdchit->wire()==maxwire_first) continue;
					if(TMath::Abs(maxresl)<TMath::Abs(cdchit->resl())){
						maxlayer=layer; maxwire=cdchit->wire();maxresl=cdchit->resl();
					}
				}
			}
#ifdef DEBUG
			std::cout << "max = l:" << maxlayer << ", w:" << maxwire << std::endl;
			std::cout << maxresl << std::endl;
#endif
			FillHist("CDC_Chi21_vs_Chi22", tmpchi, cdc->Chi());
			if(tmpchi<cdc->Chi()){
				cdc->Retiming(cdsMan,conf,beta,true);
#ifdef DEBUG
				std::cout << "Third fitting -----" << std::endl;
				std::cout << "chi = " << cdc->Chi() << std::endl;
				std::cout << "max = l:" << maxlayer << ", w:" << maxwire << std::endl;
				std::cout << maxresl << std::endl;
#endif
			}
		}

		FillHist("CDS_NumOfTrack",icut); icut++; /* After Retiming */

		FillHist("CDC_Chi2",cdc->Chi());
		//if(cdc->Chi()>30) continue;

		TVector3 cdsvtx,beamvtx;
		double par[5];
		cdc->GetParameters(CID_CDC,par,cdsvtx);
		double mom=cdc->Momentum();
		if(!cdc->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
		dis = (beamvtx-cdsvtx).Mag();

		double beam_out=0, beam_tof=0;
		ELossTools::CalcElossBeamTGeo(beam->t0pos(),beamvtx,beam->mom(),kpMass,beam_out,beam_tof);
		TVector3 cdhvtx=cdc->CDHVertex();
		double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
		double beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
		double mass2=mom*mom*(1/(beta*beta)-1);
		FillHist("CDS_OverbetavsMomentum",1./beta,mom);
		FillHist("CDS_Mass2vsMomentum",mass2,mom);
		FillHist("CDS_OverbetavsdE",1./beta,de);
		FillHist("CDS_Mass2vsdE",mass2,de);

		double calc_beta,tofvtxcdc;
		if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
		beta = calc_beta;
		FillHist("CDS_OverbetavsMomentum_AfterFindMass",1./beta,mom);
		FillHist("CDS_Mass2vsMomentum_AfterFindMass",mass2,mom);
		FillHist("CDS_OverbetavsdE_AfterFindMass",1./beta,de);
		FillHist("CDS_Mass2vsdE_AfterFindMass",mass2,de);

		int pid=TrackTools::PID(mom,mass2);
		pid=TrackTools::PIDcorr3(mom,mass2);
		cdc->SetPID(pid);

		TVector3 cdsvtx2,beamvtx2;
		double tof=0,tmpl=0;
		if(!cdc->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],beamvtx2,cdsvtx2,tof,tmpl,true)) continue;
		TVector3 Pp;
		if(!cdc->GetMomentum(cdsvtx2,Pp,true,true)) continue;
		mom = mom/TMath::Abs(mom)*Pp.Mag();
		FillHist("CDS_OverbetavsMomentum_AfterELossCorr",1./beta,mom);
		FillHist("CDS_Mass2vsMomentum_AfterELossCorr",mass2,mom);
		FillHist("CDS_OverbetavsdE_AfterELossCorr",1./beta,de);
		FillHist("CDS_Mass2vsdE_AfterELossCorr",mass2,de);

		for(int layer=1;layer<=15;layer++){
			int mult = cdc->nTrackHit(layer);
			for(int i=0;i<mult;i++){
				CDCHit *hit = cdc->TrackHit(cdsMan,layer,i);
				double resl = hit->resl();
				FillHist(Form("CDCl%d_Residual",layer),resl,1);
				if(cdc->Chi()<30){
					FillHist(Form("CDCl%d_Residual_withChi2Cut",layer),resl,1);
				}
			}
		}

		FillHist("CDS_NumOfTrack",icut); icut++; /* After Analysis */

		/* Set pCDS */
		pCDS* cdstrack = new pCDS();
		cdstrack->SetTrackID(id);
		cdstrack->SetPID(pid);
		cdstrack->SetdE(de);
		cdstrack->SetMomentum(mom);
		cdstrack->SetMass(TMath::Sqrt(mass2));
		cdstrack->SetMass2(mass2);
		cdstrack->SetBeta(beta);
		cdstrack->SetTOF(cdhtime-beam->t0time());
		cdstrack->SetFL(tmpl);
		cdstrack->SetChi2(cdc->Chi());
		cdstrack->SetProb(TMath::Prob(cdc->Chi()*(double)cdc->Dof(),cdc->Dof()));
		cdstrack->SetVertexCDC(cdsvtx2);
		cdstrack->SetVertexBeam(beamvtx2);
		cdstrack->SetMomDir(Pp);
		cdstrack->SetVBDis((beamvtx2-cdsvtx2).Mag());
		for(int icdh=0; icdh<cdc->nCDHHit(); icdh++){
			cdstrack->SetCDHSeg(cdc->CDHHit(cdsMan,icdh)->seg());
		}
		for(int iih=0; iih<cdc->nIHHit(); iih++){
			cdstrack->SetIHSeg(cdc->IHHit(cdsMan,iih)->seg());
		}

		//particle->AddCDS(*cdstrack);
		delete cdstrack;
	}


	cdhcontainer.clear();
	return true;
}

bool MyAnalysisCheck::FillHist(TString name, double val1, int weight)
{
	TH1F* h1 = (TH1F*)gFile -> Get(name);
	if(h1&&weight>0){ 
		for(int i=0; i<weight; i++){
			h1 -> Fill(val1,1);
		}
		return true;
	}
	else {
		return false;
	}
}

bool MyAnalysisCheck::FillHist(TString name, TString val1, int weight)
{
	TH1F* h1 = (TH1F*)gFile -> Get(name);
	if(h1&&weight>0){
		for(int i=0; i<weight; i++){
			h1 -> Fill(val1,1);
		}
		return true;
	}
	else {
		return false;
	}
}

bool MyAnalysisCheck::FillHist(TString name, double val1, double val2, int weight)
{
	TH2F* h2 = (TH2F*)gFile -> Get(name);
	if(h2&&weight>0){
		for(int i=0; i<weight; i++){
			h2 -> Fill(val1,val2,1);
		}
		return true;
	}
	else {
		return false;
	}
}

bool MyAnalysisCheck::FillHist(TString name, TString val1, TString val2, int weight)
{
	TH2F* h2 = (TH2F*)gFile -> Get(name);
	if(h2&&weight>0){
		for(int i=0; i<weight; i++){
			h2 -> Fill(val1,val2,1);
		}
		return true;
	}
	else {
		return false;
	}
}

void MyAnalysisCheck::CutCondition()
{
}

bool MyAnalysisCheck::Initialize(ConfMan* confMan)
{
	std::cout << "### MyAnalysisCheck::Initialize ###" << std::endl;

	std::string ofname = confMan->GetOutFileName();
	ofname.insert(ofname.find(".root"),"_anaCheck");

	rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
	rtFile -> cd();

	std::cout << "== Finish Histogram Initialization ==" << std::endl;

	new TH1F( "CDH_nHit", "Number of Hit;Num. of hits;Counts", 10, 0, 10);
	new TH1F( "CDH_nCluster", "Number of Clusters;Num. of clusters;Counts", 10, 0, 10);
	new TH2F( "CDS_nGoodTrack_vs_nHit", "Number of Good track vs. Hit;Num. of track;Num. of hits", 10, 0, 10, 10, 0, 10);
	new TH2F( "CDS_nGoodTrack_vs_nCluster", "Number of Good track vs. Clusters;Num. of track;Num. of clusters", 10, 0, 10, 10, 0, 10);
	new TH2F( "CDS_nGoodTrack_vs_nCluster_Re", "Number of Good track vs. Clusters;Num. of track;Num. of clusters", 10, 0, 10, 10, 0, 10);
	new TH2F( "CDS_nCluster_vs_nCluster_Re", "Number of Cluster vs. Clusters;Num. of clusters;Num. of clusters", 10, 0, 10, 10, 0, 10);
	new TH1F( "CDH_nHitinCluster", "Number of Hit in a Cluster;Num. of hit;Counts", 10, 0, 10);
	new TH1F( "CDH_nHitinCluster_Re", "Number of Hit in a Cluster;Num. of hit;Counts", 10, 0, 10);
	for(int i=1; i<=35; i++){
		new TH2F( Form("CDHCluster_timediff_vs_posdiff_%d%d",i,i+1), "time vs. pos;time diff. (ns);pos. diff. (cm)", 200, 0, 20, 200, 0, 40 );
	}


	new TH1F( Form("CDS_NumOfTrack"), Form("Number of tracks in CDS;Cut Number;Counts"), 10, 0, 10 );

	// CDC
	const int ndc=1;
	TString dcname[ndc]={"CDC"};

	for(int idc=0;idc<ndc;idc++){
		std::cout << "Define Histgram for " << dcname[idc] << std::endl;
		new TH1F( Form("%s_NumOfCDH",dcname[idc].Data()), Form("N cdh hits in %s;Num. of cdh hits;Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
		new TH1F( Form("%s_nGoodTrack",dcname[idc].Data()), Form("N good track %s;# of good tracks;Counts",dcname[idc].Data()), 10, 0, 10 );
		new TH1F( Form("%s_Chi2",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2};Counts",dcname[idc].Data()), 2000, 0.0, 100.0 );
		new TH2F( Form("%s_Chi21_vs_Chi22",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2};Counts",dcname[idc].Data()), 1000, 0.0, 100.0, 1000, 0.0, 100.0 );
		new TH1F( Form("%s_DPhi",dcname[idc].Data()), Form("d#phi %s;#delta #phi;Counts",dcname[idc].Data()), 3000, -1.5, 1.5 );
	}
	for(int layer=1;layer<=15;layer++){
		int idc=0;
		new TH1F( Form("%sl%d_Residual",dcname[idc].Data(),layer), Form("Residual %s-Layer%d;Residual (cm);Counts",dcname[idc].Data(),layer), 4000, -0.2, 0.2 );
		new TH1F( Form("%sl%d_Residual_withChi2Cut",dcname[idc].Data(),layer), Form("Residual %s-Layer%d;Residual (cm);Counts",dcname[idc].Data(),layer), 4000, -0.2, 0.2 );
	}
	// CDS Particles
	std::cout << "Define Histograms for CDS particles" << std::endl;
	new TH2F( "CDS_Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
	new TH2F( "CDS_Mass2vsMomentum_AfterFindMass", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
	new TH2F( "CDS_Mass2vsMomentum_AfterELossCorr", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
	new TH2F( "CDS_Mass2vsMomentum_Fiducial", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
	new TH2F( "CDS_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
	new TH2F( "CDS_OverbetavsMomentum_AfterFindMass", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
	new TH2F( "CDS_OverbetavsMomentum_AfterELossCorr", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
	new TH2F( "CDS_OverbetavsMomentum_Fiducial", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
	new TH2F( "CDS_Mass2vsdE", "Mass^{2} vs. dE;Mass^{2} (GeV/c^{2})^{2};dE (MeVee)", 3000, -5, 10, 800, 0, 40 );
	new TH2F( "CDS_Mass2vsdE_AfterFindMass", "Mass^{2} vs. dE;Mass^{2} (GeV/c^{2})^{2};dE (MeVee)", 3000, -5, 10, 800, 0, 40 );
	new TH2F( "CDS_Mass2vsdE_AfterELossCorr", "Mass^{2} vs. dE;Mass^{2} (GeV/c^{2})^{2};dE (MeVee)", 3000, -5, 10, 800, 0, 40 );
	new TH2F( "CDS_Mass2vsdE_Fiducial", "Mass^{2} vs. dE;Mass^{2} (GeV/c^{2})^{2};dE (MeVee)", 3000, -5, 10, 800, 0, 40 );
	new TH2F( "CDS_OverbetavsdE", "1/#beta vs. dE;1/#beta;dE (MeVee)", 1000, 0, 10, 800, 0, 40 );
	new TH2F( "CDS_OverbetavsdE_AfterFindMass", "1/#beta vs. dE;1/#beta;dE (MeVee)", 1000, 0, 10, 800, 0, 40 );
	new TH2F( "CDS_OverbetavsdE_AfterELossCorr", "1/#beta vs. dE;1/#beta;dE (MeVee)", 1000, 0, 10, 800, 0, 40 );
	new TH2F( "CDS_OverbetavsdE_Fiducial", "1/#beta vs. dE;1/#beta;dE (MeVee)", 1000, 0, 10, 800, 0, 40 );
	TString pname[13]  = {"pip","p","d","t","he","pim","k","o","pipi","ppim","kpip","pk","proo"};
	TString pname2[13] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o","#pi^{+}#pi^{-}","p#pi^{-}","#pi^{+}K^{-}","pK^{-}","o"};
	for(int ip=0; ip<13; ip++){
		if(ip<7){
			new TH1F( Form("%s_VDIS",pname[ip].Data()), Form("VDIS of %s;VDIS (cm);Coutns",pname2[ip].Data()), 100, 0, 10);
			new TH2F( Form("%s_Vertex_XY",pname[ip].Data()), Form("Vertex XY plane (%s);x (cm);y (cm)",pname2[ip].Data()), 600, -15, 15, 600, -15, 15 );
			new TH2F( Form("%s_Vertex_ZX",pname[ip].Data()), Form("Vertex ZX plane (%s);z (cm);x (cm)",pname2[ip].Data()), 1200, -30, 30, 600, -15, 15 );
			new TH2F( Form("%s_Vertex_ZY",pname[ip].Data()), Form("Vertex ZY plane (%s);z (cm);y (cm)",pname2[ip].Data()), 1200, -30, 30, 600, -15, 15 );
			new TH1F( Form("%s_Vertex_Z",pname[ip].Data()), Form("Vertex Z axis (%s);z (cm);Counts",pname2[ip].Data()), 1200, -30, 30 );
		}
		new TH1F( Form("%s_Mass_CM",pname[ip].Data()), Form("Invariant Mass (%s);IM(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
		new TH1F( Form("%s_Momentum_CM",pname[ip].Data()), Form("Momentum;Momentum (%s) (GeV/c);Counts",pname2[ip].Data()), 5000, 0.0, 5.0);
		new TH1F( Form("%s_CosT_CM",pname[ip].Data()), Form("cos#theta (%s);cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5 );
		new TH1F( Form("%s_MMass_CM",pname[ip].Data()), Form("Missing mass (%s);Missing Mass(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
		new TH1F( Form("%s_MMomentum_CM",pname[ip].Data()), Form("Missing Momentum (%s);Missing Momentum(%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
		new TH1F( Form("%s_MCosT_CM",pname[ip].Data()), Form("Missing cos#theta (%s);Missing cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5);
	}
	new TH1F( Form("All_VDIS"), Form("VDIS of All;VDIS (cm);Coutns"), 100, 0, 10);
	new TH1F( "CDS_NumOfParticle", "Number of CDS tracks", 10, 0, 10);
	TH2F* h2;
	h2 = new TH2F("CDS_Particle","Detected particles by CDS;Particle;Number of all tracks",6,0,6,6,0,6);
	h2->GetYaxis()->SetBinLabel(1,"1 track");
	h2->GetYaxis()->SetBinLabel(2,"2 track");
	h2->GetYaxis()->SetBinLabel(3,"3 track");
	h2->GetYaxis()->SetBinLabel(4,"4 track");
	h2->GetYaxis()->SetBinLabel(5,"5 track");
	h2->GetYaxis()->SetBinLabel(6,"6 track");
	h2->GetXaxis()->SetBinLabel(1,"#pi^{+}");
	h2->GetXaxis()->SetBinLabel(2,"#pi^{-}");
	h2->GetXaxis()->SetBinLabel(3,"K^{-}");
	h2->GetXaxis()->SetBinLabel(4,"p");
	h2->GetXaxis()->SetBinLabel(5,"d");
	h2->GetXaxis()->SetBinLabel(6,"Other");

	// Vertex
	std::cout << "Define Histograms for Vertex" << std::endl;
	new TH2F( "Vertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH2F( "Vertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH2F( "Vertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH1F( "Vertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
	new TH2F( "Vertex_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH2F( "Vertex_ZX_YCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH2F( "Vertex_ZY_XCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH1F( "Vertex_Z_XYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
	new TH2F( "Vertex_XY_Fiducial", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH2F( "Vertex_ZX_Fiducial", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH2F( "Vertex_ZY_Fiducial", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH1F( "Vertex_Z_Fiducial", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
	new TH2F( "Vertex_XY_FiducialZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH2F( "Vertex_ZX_FiducialYCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH2F( "Vertex_ZY_FiducialXCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH1F( "Vertex_Z_FiducialXYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
	new TH2F( "Vertex_XY_TargetCell", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH2F( "Vertex_ZX_TargetCell", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH2F( "Vertex_ZY_TargetCell", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );

	new TH2F( "PVertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH2F( "PVertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH2F( "PVertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH1F( "PVertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
	new TH2F( "PVertex_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH1F( "PVertex_Z_XYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
	new TH2F( "PVertex_XY_Fiducial", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH2F( "PVertex_ZX_Fiducial", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH2F( "PVertex_ZY_Fiducial", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
	new TH1F( "PVertex_Z_Fiducial", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
	new TH2F( "PVertex_XY_FiducialZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
	new TH1F( "PVertex_Z_FiducialXYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );

	std::cout << "=== End of [MyAnalysisCDS::Initialize] === " << std::endl;
	return true;

}
