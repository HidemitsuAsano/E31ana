// MyAnalysisCDS.cpp

#include "MyAnalysisCDS.h"

MyAnalysisCDS::MyAnalysisCDS(TFile* rtFile, ConfMan* conf, bool simflag)
{
  Initialize(rtFile, conf);
  Clear();
  SIMULATION = simflag;
}

MyAnalysisCDS::~MyAnalysisCDS()
{
}

void MyAnalysisCDS::Clear()
{
}

bool MyAnalysisCDS::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  FillHist("CDC_nGoodTrack",cdstrackMan->nGoodTrack());
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  // ###################### //
  // CDS all track analysis //
  // ###################### //
  for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
    int icut = 0;
    FillHist("CDS_NumOfTrack",icut); icut++; /* All Track */

    CDSTrack* cdc=cdstrackMan->GoodTrack(id);
		FillHist("CDC_NumOfCDH",cdc->nCDHHit());
    if(cdc->nCDHHit()<1) continue;
    FillHist("CDS_NumOfTrack",icut); icut++; /* CDH 1 hit in Track */

		if(cdc->FittingLevel()!=5){
			double cdhtime=0,dis=0,dphi=0;
			int cdhseg=0;
			if(!cdc->GetCDHHit2(cdsMan,cdhseg,cdhtime,dphi)) continue;
			HodoscopeLikeHit* cdh=0;
			for(int i=0; i<cdsMan->nCDH(); i++){
				HodoscopeLikeHit* hit = cdsMan->CDH(i);
				if(hit->CheckRange()){
					if(hit->seg()==cdhseg) cdh=hit;
				}
			}
		FillHist("CDC_NumOfCDH",cdc->nCDHHit());
			double de = cdh->emean();

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
		double cdhtime=0,dis=0;
		int cdhseg=0;
		if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
		HodoscopeLikeHit* cdh=0;
		for(int i=0; i<cdsMan->nCDH(); i++){
			HodoscopeLikeHit* hit = cdsMan->CDH(i);
			if(hit->CheckRange()){
				if(hit->seg()==cdhseg) cdh=hit;
			}
		}
		double de = cdh->emean();

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
		FillHist("CDS_dEvsMomentum",de,mom);
		FillHist("CDS_OverbetavsdE",1./beta,de);

		double calc_beta,tofvtxcdc;
		if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
		beta = calc_beta;
		FillHist("CDS_OverbetavsMomentum_AfterFindMass",1./beta,mom);
		FillHist("CDS_Mass2vsMomentum_AfterFindMass",mass2,mom);
		FillHist("CDS_dEvsMomentum_AfterFindMass",de,mom);
		FillHist("CDS_OverbetavsdE_AfterFindMass",1./beta,de);

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
		FillHist("CDS_dEvsMomentum_AfterELossCorr",de,mom);
		FillHist("CDS_OverbetavsdE_AfterELossCorr",1./beta,de);

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

		particle->AddCDS(*cdstrack);
	}

	for(int it=0; it<particle->nCDS(); it++){
		pCDS* cdstrack = particle->cdsi(it);
		double mass2 = cdstrack->mass2();
		double beta = cdstrack->beta();
		double mom = cdstrack->mom();
		TVector3 vtxb = cdstrack->vbeam();
		TVector3 vtxh = cdstrack->vcdc();
		FillHist("Vertex_XY", vtxb.X(), vtxb.Y());
		FillHist("Vertex_ZX", vtxb.Z(), vtxb.X());
		FillHist("Vertex_ZY", vtxb.Z(), vtxb.Y());
		FillHist("Vertex_Z", vtxb.Z());
		if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
			FillHist("Vertex_XY_ZCut", vtxb.X(), vtxb.Y());
		}
		if(-0.5<vtxb.X()&&vtxb.X()<0.5){
			FillHist("Vertex_ZY_XCut", vtxb.Z(), vtxb.Y());
		}
		if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
			FillHist("Vertex_ZX_YCut", vtxb.Z(), vtxb.X());
		}
		if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
			FillHist("Vertex_Z_XYCut", vtxb.Z());
		}

		if(GeomTools::GetID(vtxb)==CID_Fiducial){
			FillHist("Vertex_XY_Fiducial", vtxb.X(), vtxb.Y());
			FillHist("Vertex_ZX_Fiducial", vtxb.Z(), vtxb.X());
			FillHist("Vertex_ZY_Fiducial", vtxb.Z(), vtxb.Y());
			FillHist("Vertex_Z_Fiducial", vtxb.Z());
			if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
				FillHist("Vertex_XY_FiducialZCut", vtxb.X(), vtxb.Y());
			}
			if(-0.5<vtxb.X()&&vtxb.X()<0.5){
				FillHist("Vertex_ZY_FiducialXCut", vtxb.Z(), vtxb.Y());
			}
			if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
				FillHist("Vertex_ZX_FiducialYCut", vtxb.Z(), vtxb.X());
			}
			if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
				FillHist("Vertex_Z_FiducialXYCut", vtxb.Z());
			}
			FillHist("CDS_OverbetavsMomentum_Fiducial",1./beta,mom);
			FillHist("CDS_Mass2vsMomentum_Fiducial",mass2,mom);
		}
		if(GeomTools::GetID(vtxb)==CID_TarCell){
			FillHist("Vertex_XY_TargetCell", vtxb.X(), vtxb.Y());
			FillHist("Vertex_ZX_TargetCell", vtxb.Z(), vtxb.X());
			FillHist("Vertex_ZY_TargetCell", vtxb.Z(), vtxb.Y());
		}
	}

	// ##################### //
	// CDS Particle analysis //
	// ##################### //
	int ncds[6] = {0,0,0,0,0,0};
	ncds[0] = particle->nPiplus();
	ncds[1] = particle->nPiminus();
	ncds[2] = particle->nKaon();
	ncds[3] = particle->nProton();
	ncds[4] = particle->nDeuteron();
	ncds[5] = particle->nTriton();
	ncds[5] += particle->nHelium3();
	ncds[5] += particle->nOther();

	FillHist("CDS_NumOfParticle",particle->nCDS());
	for(int it=0; it<particle->nCDS(); it++){
		pCDS* cds = particle->cdsi(it);
		if(GeomTools::GetID(cds->vbeam())==CID_Fiducial){
			double cds_vdis = cds->vdis();
			TVector3 vtxb = cds->vbeam();
			TVector3 ZeroV;
			TLorentzVector Ltgt;
			if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
				Ltgt.SetVectM(ZeroV,pMass);
			}
			if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
				Ltgt.SetVectM(ZeroV,dMass);
			}
			if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
				Ltgt.SetVectM(ZeroV,ThreeHeMass);
			}
			if(dlist->GetMaterial(CID_Target)=="Vacuum"){
				Ltgt.SetVectM(ZeroV,ThreeHeMass);
			}
			TLorentzVector Lbeam = beam->GetLorentzVector(vtxb);
			TVector3 boost = (Ltgt+Lbeam).BoostVector();
			TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
			TLorentzVector Lcds = cds->GetLorentzVector();
			TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
			TString frame[2] = {"Lab","CM"};
			TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
			double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
			double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
			double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
			double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
			double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
			double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};

			FillHist(Form("%s_Vertex_XY",pname[cds->pid()].Data()), vtxb.X(), vtxb.Y());
			FillHist(Form("%s_Vertex_ZX",pname[cds->pid()].Data()), vtxb.Z(), vtxb.X());
			FillHist(Form("%s_Vertex_ZY",pname[cds->pid()].Data()), vtxb.Z(), vtxb.Y());
			FillHist(Form("%s_Vertex_Z",pname[cds->pid()].Data()), vtxb.Z());
			FillHist(Form("%s_VDIS",pname[cds->pid()].Data()),cds_vdis);
			FillHist(Form("All_VDIS"),cds_vdis);
			for(int f=1; f<2; f++){
				FillHist(Form("%s_Mass_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mass[f]);
				FillHist(Form("%s_Momentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mom[f]);
				FillHist(Form("%s_CosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_cost[f]);
				FillHist(Form("%s_MMass_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f]);
				FillHist(Form("%s_MMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmom[f]);
				FillHist(Form("%s_MCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mcost[f]);
			}
		}
	}

	// ###################### //
	// CDS 2 track analysis //
	// ###################### //
	for( int id1=0; id1<particle->nCDS(); id1++ ){
		for( int id2=id1+1; id2<particle->nCDS(); id2++ ){
			pCDS* cds1 = particle->cdsi(id1);
			pCDS* cds2 = particle->cdsi(id2);
			CDSTrack* cdc1=cdstrackMan->GoodTrack(cds1->id());
			CDSTrack* cdc2=cdstrackMan->GoodTrack(cds2->id());
			TVector3 vtx=DEFVECT,vtx1=DEFVECT,vtx2=DEFVECT;
			if(!TrackTools::Calc2HelixVertex(cdc1,cdc2,vtx1,vtx2)) continue;
			vtx = (vtx1+vtx2)*0.5;

			TVector3 Pp1=DEFVECT,Pp2=DEFVECT;
			if(!cdc1->GetMomentum(vtx1,Pp1,true,true)) continue;
			if(!cdc2->GetMomentum(vtx2,Pp2,true,true)) continue;
			TVector3 Pp = Pp1+Pp2;
			TLorentzVector L1; L1.SetVectM(Pp1,cds1->pdgmass());
			TLorentzVector L2; L2.SetVectM(Pp2,cds2->pdgmass());
			double im = (L1+L2).M();
			double dist=0, dltmp=0;
			TVector3 xest=DEFVECT, nest=DEFVECT;
			MathTools::LineToLine(beam->bpcpos(),beam->bpcdir(),vtx,Pp.Unit(),dltmp,dist,xest,nest);

			pCDS* pro = new pCDS();
			pro->SetDaughterID1(cds1->id());
			pro->SetDaughterID2(cds2->id());
			pro->SetCombID((int)(pow(2,cds1->pid())+pow(2,cds2->pid())));
			pro->SetMomentum(Pp.Mag());
			pro->SetMass(im);
			pro->SetVertex(vtx);
			pro->SetVertexBeam(xest);
			pro->SetMomDir(Pp.Unit());
			pro->SetVDis((vtx1-vtx2).Mag());
			pro->SetVBDis((xest-vtx).Mag());
			pro->SetPBDCA(dist);
			pro->SetOA(Pp1.Angle(Pp2));
			particle->AddProduct(*pro);
		}
	}

	for(int it=0; it<particle->nProduct(); it++){
		pCDS* cds = particle->product(it);
		TVector3 vtxb = cds->vbeam();
		FillHist("PVertex_XY", vtxb.X(), vtxb.Y());
		FillHist("PVertex_ZX", vtxb.Z(), vtxb.X());
		FillHist("PVertex_ZY", vtxb.Z(), vtxb.Y());
		FillHist("PVertex_Z", vtxb.Z());
		if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
			FillHist("PVertex_XY_ZCut", vtxb.X(), vtxb.Y());
		}
		if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
			FillHist("PVertex_Z_XYCut", vtxb.Z());
		}
		if(GeomTools::GetID(cds->vbeam())==CID_Fiducial){
			TVector3 ZeroV;
			TLorentzVector Ltgt;
			if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
				Ltgt.SetVectM(ZeroV,pMass);
			}
			if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
				Ltgt.SetVectM(ZeroV,dMass);
			}
			if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
				Ltgt.SetVectM(ZeroV,ThreeHeMass);
			}
			if(dlist->GetMaterial(CID_Target)=="Vacuum"){
				Ltgt.SetVectM(ZeroV,ThreeHeMass);
			}
			TLorentzVector Lbeam = beam->GetLorentzVector(vtxb);
			TVector3 boost = (Ltgt+Lbeam).BoostVector();
			TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
			TLorentzVector Lcds = cds->GetLorentzVector();
			TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
			TString frame[2] = {"Lab","CM"};
			TString pname[5] = {"pipi","ppim","kpip","pk","pro_o"};
			double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
			double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
			double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
			double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
			double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
			double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};

			FillHist("PVertex_XY_Fiducial", vtxb.X(), vtxb.Y());
			FillHist("PVertex_ZX_Fiducial", vtxb.Z(), vtxb.X());
			FillHist("PVertex_ZY_Fiducial", vtxb.Z(), vtxb.Y());
			FillHist("PVertex_Z_Fiducial", vtxb.Z());
			if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
				FillHist("PVertex_XY_FiducialZCut", vtxb.X(), vtxb.Y());
			}
			if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
				FillHist("PVertex_Z_FiducialXYCut", vtxb.Z());
			}

			int combid = 4;
			int comb = cds->comb();
			if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
				combid = 0;
			}
			else if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
				combid = 1;
			}
			else if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_Kaon)){
				combid = 2;
			}
			else if(comb==pow(2,CDS_Proton)+pow(2,CDS_Kaon)){
				combid = 3;
			}

			for(int f=1; f<2; f++){
				FillHist(Form("%s_Mass_%s",pname[combid].Data(),frame[f].Data()),cds_mass[f]);
				FillHist(Form("%s_Momentum_%s",pname[combid].Data(),frame[f].Data()),cds_mom[f]);
				FillHist(Form("%s_CosT_%s",pname[combid].Data(),frame[f].Data()),cds_cost[f]);
				FillHist(Form("%s_MMass_%s",pname[combid].Data(),frame[f].Data()),cds_mmass[f]);
				FillHist(Form("%s_MMomentum_%s",pname[combid].Data(),frame[f].Data()),cds_mmom[f]);
				FillHist(Form("%s_MCosT_%s",pname[combid].Data(),frame[f].Data()),cds_mcost[f]);
			}
		}
	}

	return true;

}

bool MyAnalysisCDS::DoAnalysisSim(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
	if(!SIMULATION) return false;

	DetectorList *dlist=DetectorList::GetInstance();
	FillHist("CDC_nGoodTrack",cdstrackMan->nGoodTrack());
	if(particle->nBeam()!=1) return false;
	pBeam* beam = particle->beam(0);
	if(cdstrackMan->nGoodTrack()!=2&&cdstrackMan->nGoodTrack()!=3&&cdstrackMan->nGoodTrack()!=4) return false;

	// ###################### //
	// CDS all track analysis //
	// ###################### //
	for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
		int icut = 0;
		FillHist("CDS_NumOfTrack",icut); icut++; /* All Track */

		CDSTrack* cdc=cdstrackMan->GoodTrack(id);
		if(cdc->nCDHHit()!=1) continue;
		FillHist("CDS_NumOfTrack",icut); icut++; /* CDH 1 hit in Track */

		FillHist("CDC_Chi2",cdc->Chi());
		//if(cdc->Chi()>30) continue;
		double cdhtime=0,dis=0;
		int cdhseg=0;
		if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
		HodoscopeLikeHit* cdh=0;
		for(int i=0; i<cdsMan->nCDH(); i++){
			HodoscopeLikeHit* hit = cdsMan->CDH(i);
			if(hit->CheckRange()){
				if(hit->seg()==cdhseg) cdh=hit;
			}
		}
		double de = cdh->emean();

		TVector3 cdsvtx,beamvtx;
		double par[5];
		cdc->GetParameters(CID_CDC,par,cdsvtx);
		double mom=cdc->Momentum();
		if(!cdc->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
		dis = (beamvtx-cdsvtx).Mag();

		double beam_out=0, beam_tof=0;
		ELossTools::CalcElossBeamTGeo(beamvtx,beam->t0pos(),beam->mom(),kpMass,beam_out,beam_tof);
		TVector3 cdhvtx=cdc->CDHVertex();
		double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
		double beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
		double mass2=mom*mom*(1/(beta*beta)-1);
		FillHist("CDS_OverbetavsMomentum",1./beta,mom);
		FillHist("CDS_Mass2vsMomentum",mass2,mom);
		FillHist("CDS_dEvsMomentum",de,mom);

		double calc_beta,tofvtxcdc;
		if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
		beta = calc_beta;
		FillHist("CDS_OverbetavsMomentum_AfterFindMass",1./beta,mom);
		FillHist("CDS_Mass2vsMomentum_AfterFindMass",mass2,mom);
		FillHist("CDS_dEvsMomentum_AfterFindMass",de,mom);

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

		particle->AddCDS(*cdstrack);
	}

	for(int it=0; it<particle->nCDS(); it++){
		pCDS* cdstrack = particle->cdsi(it);
		double mass2 = cdstrack->mass2();
		double beta = cdstrack->beta();
		double mom = cdstrack->mom();
		TVector3 vtxb = cdstrack->vbeam();
		TVector3 vtxh = cdstrack->vcdc();
		FillHist("Vertex_XY", vtxb.X(), vtxb.Y());
		FillHist("Vertex_ZX", vtxb.Z(), vtxb.X());
		FillHist("Vertex_ZY", vtxb.Z(), vtxb.Y());
		FillHist("Vertex_Z", vtxb.Z());
		if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
			FillHist("Vertex_XY_ZCut", vtxb.X(), vtxb.Y());
		}
		if(-0.5<vtxb.X()&&vtxb.X()<0.5){
			FillHist("Vertex_ZY_XCut", vtxb.Z(), vtxb.Y());
		}
		if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
			FillHist("Vertex_ZX_YCut", vtxb.Z(), vtxb.X());
		}
		if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
			FillHist("Vertex_Z_XYCut", vtxb.Z());
		}

		if(GeomTools::GetID(vtxb)==CID_Fiducial){
			FillHist("Vertex_XY_Fiducial", vtxb.X(), vtxb.Y());
			FillHist("Vertex_ZX_Fiducial", vtxb.Z(), vtxb.X());
			FillHist("Vertex_ZY_Fiducial", vtxb.Z(), vtxb.Y());
			FillHist("Vertex_Z_Fiducial", vtxb.Z());
			if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
				FillHist("Vertex_XY_FiducialZCut", vtxb.X(), vtxb.Y());
			}
			if(-0.5<vtxb.X()&&vtxb.X()<0.5){
				FillHist("Vertex_ZY_FiducialXCut", vtxb.Z(), vtxb.Y());
			}
			if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
				FillHist("Vertex_ZX_FiducialYCut", vtxb.Z(), vtxb.X());
			}
			if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
				FillHist("Vertex_Z_FiducialXYCut", vtxb.Z());
			}
			FillHist("CDS_OverbetavsMomentum_Fiducial",1./beta,mom);
			FillHist("CDS_Mass2vsMomentum_Fiducial",mass2,mom);
		}
		if(GeomTools::GetID(vtxb)==CID_TarCell){
			FillHist("Vertex_XY_TargetCell", vtxb.X(), vtxb.Y());
			FillHist("Vertex_ZX_TargetCell", vtxb.Z(), vtxb.X());
			FillHist("Vertex_ZY_TargetCell", vtxb.Z(), vtxb.Y());
		}
	}

	// ##################### //
	// CDS Particle analysis //
	// ##################### //
	int ncds[6] = {0,0,0,0,0,0};
	ncds[0] = particle->nPiplus();
	ncds[1] = particle->nPiminus();
	ncds[2] = particle->nKaon();
	ncds[3] = particle->nProton();
	ncds[4] = particle->nDeuteron();
	ncds[5] = particle->nTriton();
	ncds[5] += particle->nHelium3();
	ncds[5] += particle->nOther();

	FillHist("CDS_NumOfParticle",particle->nCDS());
	for(int it=0; it<particle->nCDS(); it++){
		pCDS* cds = particle->cdsi(it);
		if(GeomTools::GetID(cds->vbeam())==CID_Fiducial){
			double cds_vdis = cds->vdis();
			TVector3 vtxb = cds->vbeam();
			TVector3 ZeroV;
			TLorentzVector Ltgt;
			if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
				Ltgt.SetVectM(ZeroV,pMass);
			}
			if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
				Ltgt.SetVectM(ZeroV,dMass);
			}
			if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
				Ltgt.SetVectM(ZeroV,ThreeHeMass);
			}
			if(dlist->GetMaterial(CID_Target)=="Vacuum"){
				Ltgt.SetVectM(ZeroV,ThreeHeMass);
			}
			TLorentzVector Lbeam = beam->GetLorentzVector(vtxb);
			TVector3 boost = (Ltgt+Lbeam).BoostVector();
			TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
			TLorentzVector Lcds = cds->GetLorentzVector();
			TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
			TString frame[2] = {"Lab","CM"};
			TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
			double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
			double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
			double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
			double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
			double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
			double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};

			FillHist(Form("%s_Vertex_XY",pname[cds->pid()].Data()), vtxb.X(), vtxb.Y());
			FillHist(Form("%s_Vertex_ZX",pname[cds->pid()].Data()), vtxb.Z(), vtxb.X());
			FillHist(Form("%s_Vertex_ZY",pname[cds->pid()].Data()), vtxb.Z(), vtxb.Y());
			FillHist(Form("%s_Vertex_Z",pname[cds->pid()].Data()), vtxb.Z());
			FillHist(Form("%s_VDIS",pname[cds->pid()].Data()),cds_vdis);
			FillHist(Form("All_VDIS"),cds_vdis);
			for(int f=1; f<2; f++){
				FillHist(Form("%s_Mass_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mass[f]);
				FillHist(Form("%s_Momentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mom[f]);
				FillHist(Form("%s_CosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_cost[f]);
				FillHist(Form("%s_MMass_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f]);
				FillHist(Form("%s_MMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmom[f]);
				FillHist(Form("%s_MCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mcost[f]);
			}
		}
	}

	// ###################### //
	// CDS 2 track analysis //
	// ###################### //
	for( int id1=0; id1<particle->nCDS(); id1++ ){
		for( int id2=id1+1; id2<particle->nCDS(); id2++ ){
			pCDS* cds1 = particle->cdsi(id1);
			pCDS* cds2 = particle->cdsi(id2);
			CDSTrack* cdc1=cdstrackMan->GoodTrack(cds1->id());
			CDSTrack* cdc2=cdstrackMan->GoodTrack(cds2->id());
			TVector3 vtx=DEFVECT,vtx1=DEFVECT,vtx2=DEFVECT;
			if(!TrackTools::Calc2HelixVertex(cdc1,cdc2,vtx1,vtx2)) continue;
			vtx = (vtx1+vtx2)*0.5;

			TVector3 Pp1=DEFVECT,Pp2=DEFVECT;
			if(!cdc1->GetMomentum(vtx1,Pp1,true,true)) continue;
			if(!cdc2->GetMomentum(vtx2,Pp2,true,true)) continue;
			TVector3 Pp = Pp1+Pp2;
			TLorentzVector L1; L1.SetVectM(Pp1,cds1->pdgmass());
			TLorentzVector L2; L2.SetVectM(Pp2,cds2->pdgmass());
			double im = (L1+L2).M();
			double dist=0, dltmp=0;
			TVector3 xest=DEFVECT, nest=DEFVECT;
			MathTools::LineToLine(beam->bpcpos(),beam->bpcdir(),vtx,Pp.Unit(),dltmp,dist,xest,nest);

			pCDS* pro = new pCDS();
			pro->SetDaughterID1(cds1->id());
			pro->SetDaughterID2(cds2->id());
			pro->SetCombID((int)(pow(2,cds1->pid())+pow(2,cds2->pid())));
			pro->SetMomentum(Pp.Mag());
			pro->SetMass(im);
			pro->SetVertex(vtx);
			pro->SetVertexBeam(xest);
			pro->SetMomDir(Pp.Unit());
			pro->SetVDis((vtx1-vtx2).Mag());
			pro->SetVBDis((xest-vtx).Mag());
			pro->SetPBDCA(dist);
			pro->SetOA(Pp1.Angle(Pp2));
			particle->AddProduct(*pro);
		}
	}

	for(int it=0; it<particle->nProduct(); it++){
		pCDS* cds = particle->product(it);
		TVector3 vtxb = cds->vbeam();
		FillHist("PVertex_XY", vtxb.X(), vtxb.Y());
		FillHist("PVertex_ZX", vtxb.Z(), vtxb.X());
		FillHist("PVertex_ZY", vtxb.Z(), vtxb.Y());
		FillHist("PVertex_Z", vtxb.Z());
		if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
			FillHist("PVertex_XY_ZCut", vtxb.X(), vtxb.Y());
		}
		if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
			FillHist("PVertex_Z_XYCut", vtxb.Z());
		}
		if(GeomTools::GetID(cds->vbeam())==CID_Fiducial){
			TVector3 ZeroV;
			TLorentzVector Ltgt;
			if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
				Ltgt.SetVectM(ZeroV,pMass);
			}
			if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
				Ltgt.SetVectM(ZeroV,dMass);
			}
			if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
				Ltgt.SetVectM(ZeroV,ThreeHeMass);
			}
			if(dlist->GetMaterial(CID_Target)=="Vacuum"){
				Ltgt.SetVectM(ZeroV,ThreeHeMass);
			}
			TLorentzVector Lbeam = beam->GetLorentzVector(vtxb);
			TVector3 boost = (Ltgt+Lbeam).BoostVector();
			TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
			TLorentzVector Lcds = cds->GetLorentzVector();
			TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
			TString frame[2] = {"Lab","CM"};
			TString pname[5] = {"pipi","ppim","kpip","pk","pro_o"};
			double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
			double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
			double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
			double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
			double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
			double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};

			FillHist("PVertex_XY_Fiducial", vtxb.X(), vtxb.Y());
			FillHist("PVertex_ZX_Fiducial", vtxb.Z(), vtxb.X());
			FillHist("PVertex_ZY_Fiducial", vtxb.Z(), vtxb.Y());
			FillHist("PVertex_Z_Fiducial", vtxb.Z());
			if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
				FillHist("PVertex_XY_FiducialZCut", vtxb.X(), vtxb.Y());
			}
			if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
				FillHist("PVertex_Z_FiducialXYCut", vtxb.Z());
			}

			int combid = 4;
			int comb = cds->comb();
			if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
				combid = 0;
			}
			else if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
				combid = 1;
			}
			else if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_Kaon)){
				combid = 2;
			}
			else if(comb==pow(2,CDS_Proton)+pow(2,CDS_Kaon)){
				combid = 3;
			}

			for(int f=1; f<2; f++){
				FillHist(Form("%s_Mass_%s",pname[combid].Data(),frame[f].Data()),cds_mass[f]);
				FillHist(Form("%s_Momentum_%s",pname[combid].Data(),frame[f].Data()),cds_mom[f]);
				FillHist(Form("%s_CosT_%s",pname[combid].Data(),frame[f].Data()),cds_cost[f]);
				FillHist(Form("%s_MMass_%s",pname[combid].Data(),frame[f].Data()),cds_mmass[f]);
				FillHist(Form("%s_MMomentum_%s",pname[combid].Data(),frame[f].Data()),cds_mmom[f]);
				FillHist(Form("%s_MCosT_%s",pname[combid].Data(),frame[f].Data()),cds_mcost[f]);
			}
		}
	}

	return true;

}

bool MyAnalysisCDS::FillHist(TString name, double val1)
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

bool MyAnalysisCDS::FillHist(TString name, double val1, double val2)
{
	TH2F* h2 = (TH2F*)gFile -> Get(name);
	if(h2){ h2 -> Fill(val1,val2);
		return true;
	}
	else {
		return false;
	}
}

bool MyAnalysisCDS::Initialize(TFile* rtFile, ConfMan* confMan)
{
	rtFile -> cd();

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
	//new TH2F( "CDS_OverbetavsMomentum_AfterFindMass", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
	//new TH2F( "CDS_OverbetavsMomentum_AfterELossCorr", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
	new TH2F( "CDS_OverbetavsMomentum_Fiducial", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
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
