// MyAnalysisBPD.cpp

#include "MyAnalysisBPD.h"

MyAnalysisBPD::MyAnalysisBPD(TFile* rt, ConfMan* conf)
{
	Initialize(conf);
	CutCondition();
	Clear();
}

MyAnalysisBPD::~MyAnalysisBPD()
{
	Clear();
	rtFile->cd();
	rtFile->Write();
	rtFile->Close();
}

void MyAnalysisBPD::Clear()
{
}

bool MyAnalysisBPD::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{

	rtFile->cd();

	DetectorList *dlist=DetectorList::GetInstance();
	GeomMapMan* gmapman=confFile->GetGeomMapManager();
	if(particle->nBeam()!=1) return false;
	pBeam* beam = particle->beam(0);

	/* Event selection */
	//if(!header->trigmode(Mode_KCDH2)) return true;
	//if(header->IsTrig(Trig_KCDH2)&&header->IsTrig(Trig_1stMix)){
	if(!header->IsTrig(Trig_Kaon)) return true;
	//}
	//else {
	//  return true;
	//}

	// Trigger Pattern
	FillHist("TriggerPattern",0);
	for( int i=0; i<20; i++ ){
		int val = header->pattern(i);
		if( 0<val ){
			FillHist("TriggerPattern",i);
		}
	}
	double bpcposX = -999, bpcposY = -999;
	double bpcdirX = -999, bpcdirY = -999;
	TVector3 bpcpos;
	gmapman->GetGPos(CID_BPD,0,bpcpos);
	beam -> bpcpos(bpcpos.Z(),bpcposX,bpcposY);
	beam -> bpcdir(bpcdirX,bpcdirY);
	bpcpos.SetX(bpcposX); bpcpos.SetY(bpcposY);
	FillHist("BPC_HitPosition",bpcpos.X(),bpcpos.Y());

	// Fill BPD
	int nBPD=0;
	for( int i=0; i<blMan->nBPD(); i++ ){
		HodoscopeLikeHit *hit = blMan->BPD(i);
		int seg = hit->seg();
		int au = hit->adcu(), ad = hit->adcd();
		int tu = hit->tdcu(), td = hit->tdcd();
		double hitpos = hit->hitpos();
		TVector3 pos = hit->pos();
		pos.SetY(hitpos);
		double timeu = hit->time(0), timed = hit->time(1);
		double ctimeu = hit->ctime(0), ctimed = hit->ctime(1);
		double eneu = hit->ene(0), ened = hit->ene(1);
		double emean = hit->emean();
		double tmean = hit->tmean();
		double ctmean = hit->ctmean();
		double tof = ctmean - beam->t0time();
		double tsub = hit->tsub();
		double ctsub = hit->ctsub();
		double timeoffs = CalcTimeOffs(hit, particle);
		double tsuboffs = CalcTimeSubOffs(hit, particle);

		FillHist(Form("BPDu%d_ADC",seg),au);
		FillHist(Form("BPDd%d_ADC",seg),ad);
		if(tu>0){
			FillHist(Form("BPDu%d_TDC",seg),tu);
			FillHist(Form("BPDu%d_Time",seg),timeu);
			FillHist(Form("BPDu%d_CTime",seg),ctimeu);
		}
		if(td>0){
			FillHist(Form("BPDd%d_TDC",seg),td);
			FillHist(Form("BPDd%d_Time",seg),timed);
			FillHist(Form("BPDd%d_CTime",seg),ctimed);
		}
		if( hit->CheckRange() ){
			FillHist(Form("BPDu%d_ADCwT",seg),au);
			FillHist(Form("BPDd%d_ADCwT",seg),ad);
			FillHist(Form("BPDu%d_dE",seg),eneu);
			FillHist(Form("BPDd%d_dE",seg),ened);
			FillHist(Form("BPD%d_dEMean",seg),emean);
			FillHist(Form("BPD%d_TimeMean",seg),tmean);
			FillHist(Form("BPD%d_TimeSub",seg),tsub);
			FillHist(Form("BPD%d_CTimeMean",seg),ctmean);
			FillHist(Form("BPD%d_CTimeSub",seg),ctsub);
			FillHist(Form("BPD%d_Position",seg),hitpos);
			FillHist(Form("BPDu%d_AvT",seg),au,tu);
			FillHist(Form("BPDd%d_AvT",seg),ad,td);
			FillHist(Form("BPDu%d_dEvTime",seg),eneu,timeu);
			FillHist(Form("BPDd%d_dEvTime",seg),ened,timed);
			FillHist(Form("BPDu%d_dEvCTime",seg),eneu,ctimeu);
			FillHist(Form("BPDd%d_dEvCTime",seg),ened,ctimed);
			FillHist(Form("BPD%d_TOFT0",seg),tof);
			FillHist(Form("BPD_HitPosition"),pos.X(),pos.Y());
			FillHist(Form("BPDvsBPC_HitPositionatBPDX"),pos.X(),bpcpos.X());
			FillHist(Form("BPDvsBPC_HitPositionatBPDY"),pos.Y(),bpcpos.Y());
			FillHist(Form("BPDvsBPC_HitPositionDifference"),pos.X()-bpcpos.X(),pos.Y()-bpcpos.Y());
			FillHist(Form("BPDvsBPC_HitPositionDifferenceXvsDirX"),pos.X()-bpcpos.X(),bpcdirX);
			FillHist(Form("BPDvsBPC_HitPositionDifferenceXvsDirY"),pos.X()-bpcpos.X(),bpcdirY);
			FillHist(Form("BPDvsBPC_HitPositionDifferenceYvsDirX"),pos.Y()-bpcpos.Y(),bpcdirX);
			FillHist(Form("BPDvsBPC_HitPositionDifferenceYvsDirY"),pos.Y()-bpcpos.Y(),bpcdirY);
			FillHist(Form("BPD%d_TimeOffset",seg),timeoffs);
			FillHist(Form("BPD%d_TimeSubOffset",seg),tsuboffs);
			FillHist("BPD_HitPat",seg);
			nBPD++;
		}
	}
	FillHist("BPD_Mult",nBPD);


}

bool MyAnalysisBPD::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisBPD::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisBPD::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisBPD::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisBPD::CutCondition()
{
}

double MyAnalysisBPD::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
	double tmpoffset = -9999.9;
	if(!hit->CheckRange()) return tmpoffset;
	pBeam* beam = particle->beam(0);
	double tof = hit->ctmean() - beam->t0time();

	TVector3 bpdpos; double bpdposX, bpdposY;
	confFile->GetGeomMapManager()->GetGPos(CID_BPD,0,bpdpos);
	beam -> bpcpos(bpdpos.Z(),bpdposX,bpdposY);
	bpdpos.SetX(bpdposX); bpdpos.SetY(bpdposY);
	TVector3 t0pos; double t0posX, t0posY;
	confFile->GetGeomMapManager()->GetGPos(CID_T0,0,t0pos);
	beam -> blcpos(t0pos.Z(),t0posX,t0posY);
	t0pos.SetX(t0posX); t0pos.SetY(t0posY);
	double calctof = (bpdpos-t0pos).Mag() / beam->beta() * Const;
	tmpoffset = tof - calctof;
	return tmpoffset;
}

double MyAnalysisBPD::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
	double tmpoffset = -9999.9;
	if(!hit->CheckRange()) return tmpoffset;
	pBeam* beam = particle->beam(0);
	double ctsub = hit->ctsub();

	TVector3 bpcpos; double bpcposX, bpcposY;
	confFile->GetGeomMapManager()->GetGPos(CID_BPD,0,bpcpos);
	beam -> bpcpos(bpcpos.Z(),bpcposX,bpcposY);
	bpcpos.SetX(bpcposX); bpcpos.SetY(bpcposY);
	double lv;
	if(confFile->GetGeomMapManager()->GetLightVelocity(CID_BPD,hit->seg(),lv)){
		double calctsub = -bpcpos.Y() / lv;
		tmpoffset = ctsub - calctsub;
	}
	return tmpoffset;
}

bool MyAnalysisBPD::Initialize(ConfMan* confMan)
{
	std::cout << "### MyAnalysisBPD::Initialize ###" << std::endl;

	confFile = confMan;

	std::string ofname = confMan->GetOutFileName();
	ofname.insert(ofname.find(".root"),"_anaBPD");

	rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
	rtFile -> cd();

	/* Trigger pattern */
	TH1F* h1 = new TH1F( "TriggerPattern", "Trigger Pattern", 18, -0.5, 17.5);
	h1->GetXaxis()->SetBinLabel(1,"All");
	h1->GetXaxis()->SetBinLabel(2,"Beam");
	h1->GetXaxis()->SetBinLabel(3,"Kaon");
	h1->GetXaxis()->SetBinLabel(4,"KCDH1f");
	h1->GetXaxis()->SetBinLabel(5,"Pion");
	h1->GetXaxis()->SetBinLabel(6,"Proton");
	h1->GetXaxis()->SetBinLabel(7,"KCDH1");
	h1->GetXaxis()->SetBinLabel(8,"KCDH2");
	h1->GetXaxis()->SetBinLabel(9,"PivBVC");
	h1->GetXaxis()->SetBinLabel(10,"PiCDH1");
	h1->GetXaxis()->SetBinLabel(11,"PiCDH2");
	h1->GetXaxis()->SetBinLabel(12,"Kf");
	h1->GetXaxis()->SetBinLabel(13,"1stMix");
	h1->GetXaxis()->SetBinLabel(14,"Charged");
	h1->GetXaxis()->SetBinLabel(15,"Neutral");
	h1->GetXaxis()->SetBinLabel(16,"Cosmic");
	h1->GetXaxis()->SetBinLabel(17,"Reject");
	h1->GetXaxis()->SetBinLabel(18,"SIM");

	/* BPD */ 
	std::cout << "Define Histograms for BPD" << std::endl;
	int NumOfBPDSegments = 70;
	for( int seg=1; seg<=NumOfBPDSegments; seg++ ){
		new TH1F( Form("BPDu%d_ADC",seg),   Form("ADC BPDU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("BPDd%d_ADC",seg),   Form("ADC BPDD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("BPDu%d_dE",seg),   Form("dE BPDU%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("BPDd%d_dE",seg),   Form("dE BPDD%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("BPD%d_dEMean",seg),   Form("Mean dE BPD%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("BPDu%d_TDC",seg),   Form("TDC BPDU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("BPDd%d_TDC",seg),   Form("TDC BPDD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("BPDu%d_Time",seg),   Form("Time BPDU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPDd%d_Time",seg),   Form("Time BPDD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPD%d_TimeMean",seg),   Form("Mean Time BPD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPD%d_TimeSub",seg),   Form("Subtract Time BPD%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPDu%d_CTime",seg),   Form("CTime BPDU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPDd%d_CTime",seg),   Form("CTime BPDD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPD%d_CTimeMean",seg),   Form("Mean CTime BPD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPD%d_CTimeSub",seg),   Form("Subtract CTime BPD%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPD%d_Position",seg),   Form("Hit position BPD%d;Position (cm);Counts",seg),    1000,    -50, 50 );
		new TH1F( Form("BPDu%d_ADCwT",seg), Form("ADC wTDC BPDU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
		new TH1F( Form("BPDd%d_ADCwT",seg), Form("ADC wTDC BPDD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
		new TH2F( Form("BPDu%d_AvT",seg),   Form("ADC TDC corr. BPDU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
		new TH2F( Form("BPDd%d_AvT",seg),   Form("ADC TDC corr. BPDD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
		new TH2F( Form("BPDu%d_dEvTime",seg),   Form("dE Time corr. BPDU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("BPDd%d_dEvTime",seg),   Form("dE Time corr. BPDD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("BPDu%d_dEvCTime",seg),   Form("dE CTime corr. BPDU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("BPDd%d_dEvCTime",seg),   Form("dE CTime corr. BPDD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH1F( Form("BPD%d_TOFT0",seg),   Form("TOF BPD%d - T0;TOF (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPD%d_TimeOffset",seg),   Form("Time offset BPD%d;Time offset (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("BPD%d_TimeSubOffset",seg),   Form("Subtract Time offset BPD%d;Subtract Time offset (ns);Counts",seg),    2000,    -100, 100 );
	}
	new TH1F( "BPD_HitPat", "Hit pattern BPD;BPD Segment;Counts", NumOfBPDSegments+1, 0, NumOfBPDSegments+1 );
	new TH1F( "BPD_Mult", "Multiplicity BPD;Multiplicity;Counts", NumOfBPDSegments+1, 0, NumOfBPDSegments+1 );
	new TH2F( Form("BPD_HitPosition"),   Form("Position BPD;X hit position (cm);Y hit position (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("BPC_HitPosition"),   Form("BPC track position at BPC;X track position at BPD (cm);Y track position at BPD (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("BPDvsBPC_HitPositionatBPDX"),   Form("Position at BPD;BPD X hit position (cm);BPC track X position at BPD (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("BPDvsBPC_HitPositionatBPDY"),   Form("Position at BPD;BPD Y hit position (cm);BPC track Y position at BPD (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("BPDvsBPC_HitPositionDifference"),   Form("Position difference between BPD and BPC;X position difference (cm);Y position difference (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("BPDvsBPC_HitPositionDifferenceXvsDirX"),   Form("Position difference(BPD-BPC) vs. X direction;X position difference (cm);X track direction"), 80, -20, 20, 100, -0.05, 0.05 );
	new TH2F( Form("BPDvsBPC_HitPositionDifferenceXvsDirY"),   Form("Position difference(BPD-BPC) vs. X direction;X position difference (cm);Y track direction"), 80, -20, 20, 100, -0.05, 0.05 );
	new TH2F( Form("BPDvsBPC_HitPositionDifferenceYvsDirX"),   Form("Position difference(BPD-BPC) vs. Y direction;Y position difference (cm);X track direction"), 80, -20, 20, 100, -0.05, 0.05 );
	new TH2F( Form("BPDvsBPC_HitPositionDifferenceYvsDirY"),   Form("Position difference(BPD-BPC) vs. Y direction;Y position difference (cm);Y track direction"), 80, -20, 20, 100, -0.05, 0.05 );

	return true;
}
