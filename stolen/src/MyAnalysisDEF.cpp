// MyAnalysisDEF.cpp

#include "MyAnalysisDEF.h"

MyAnalysisDEF::MyAnalysisDEF(TFile* rt, ConfMan* conf)
{
	Initialize(conf);
	CutCondition();
	Clear();
}

MyAnalysisDEF::~MyAnalysisDEF()
{
	Clear();
	rtFile->cd();
	rtFile->Write();
	rtFile->Close();
}

void MyAnalysisDEF::Clear()
{
}

bool MyAnalysisDEF::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
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
	gmapman->GetGPos(CID_DEF,0,bpcpos);
	beam -> bpcpos(bpcpos.Z(),bpcposX,bpcposY);
	beam -> bpcdir(bpcdirX,bpcdirY);
	bpcpos.SetX(bpcposX); bpcpos.SetY(bpcposY);
	FillHist("BPC_HitPosition",bpcpos.X(),bpcpos.Y());

	// Fill DEF
	int nDEF=0;
	for( int i=0; i<blMan->nDEF(); i++ ){
		HodoscopeLikeHit *hit = blMan->DEF(i);
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
		//double timeoffs = CalcTimeOffs(hit, particle);
		//double tsuboffs = CalcTimeSubOffs(hit, particle);

		FillHist(Form("DEFu%d_ADC",seg),au);
		FillHist(Form("DEFd%d_ADC",seg),ad);
		if(tu>0){
			FillHist(Form("DEFu%d_TDC",seg),tu);
			FillHist(Form("DEFu%d_Time",seg),timeu);
			FillHist(Form("DEFu%d_CTime",seg),ctimeu);
		}
		if(td>0){
			FillHist(Form("DEFd%d_TDC",seg),td);
			FillHist(Form("DEFd%d_Time",seg),timed);
			FillHist(Form("DEFd%d_CTime",seg),ctimed);
		}
		if( hit->CheckRange() ){
			FillHist(Form("DEFu%d_ADCwT",seg),au);
			FillHist(Form("DEFd%d_ADCwT",seg),ad);
			FillHist(Form("DEFu%d_dE",seg),eneu);
			FillHist(Form("DEFd%d_dE",seg),ened);
			FillHist(Form("DEF%d_dEMean",seg),emean);
			FillHist(Form("DEF%d_TimeMean",seg),tmean);
			FillHist(Form("DEF%d_TimeSub",seg),tsub);
			FillHist(Form("DEF%d_CTimeMean",seg),ctmean);
			FillHist(Form("DEF%d_CTimeSub",seg),ctsub);
			FillHist(Form("DEF%d_Position",seg),hitpos);
			FillHist(Form("DEFu%d_AvT",seg),au,tu);
			FillHist(Form("DEFd%d_AvT",seg),ad,td);
			FillHist(Form("DEFu%d_dEvTime",seg),eneu,timeu);
			FillHist(Form("DEFd%d_dEvTime",seg),ened,timed);
			FillHist(Form("DEFu%d_dEvCTime",seg),eneu,ctimeu);
			FillHist(Form("DEFd%d_dEvCTime",seg),ened,ctimed);
			FillHist(Form("DEF%d_TOFT0",seg),tof);
			FillHist(Form("DEF_HitPosition"),pos.X(),pos.Y());
			FillHist(Form("DEFvsBPC_HitPositionatDEFX"),pos.X(),bpcpos.X());
			FillHist(Form("DEFvsBPC_HitPositionatDEFY"),pos.Y(),bpcpos.Y());
			FillHist(Form("DEFvsBPC_HitPositionDifference"),pos.X()-bpcpos.X(),pos.Y()-bpcpos.Y());
			FillHist(Form("DEFvsBPC_HitPositionDifferenceXvsDirX"),pos.X()-bpcpos.X(),bpcdirX);
			FillHist(Form("DEFvsBPC_HitPositionDifferenceXvsDirY"),pos.X()-bpcpos.X(),bpcdirY);
			FillHist(Form("DEFvsBPC_HitPositionDifferenceYvsDirX"),pos.Y()-bpcpos.Y(),bpcdirX);
			FillHist(Form("DEFvsBPC_HitPositionDifferenceYvsDirY"),pos.Y()-bpcpos.Y(),bpcdirY);
			FillHist(Form("DEF%d_TimeOffset",seg),timeoffs);
			FillHist(Form("DEF%d_TimeSubOffset",seg),tsuboffs);
			FillHist("DEF_HitPat",seg);
			nDEF++;
		}
	}
	FillHist("DEF_Mult",nDEF);


}

bool MyAnalysisDEF::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisDEF::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisDEF::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisDEF::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisDEF::CutCondition()
{
}

double MyAnalysisDEF::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
	double tmpoffset = -9999.9;
	if(!hit->CheckRange()) return tmpoffset;
	pBeam* beam = particle->beam(0);
	double tof = hit->ctmean() - beam->t0time();

	TVector3 bpdpos; double bpdposX, bpdposY;
	confFile->GetGeomMapManager()->GetGPos(CID_DEF,0,bpdpos);
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

double MyAnalysisDEF::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
	double tmpoffset = -9999.9;
	if(!hit->CheckRange()) return tmpoffset;
	pBeam* beam = particle->beam(0);
	double ctsub = hit->ctsub();

	TVector3 bpcpos; double bpcposX, bpcposY;
	confFile->GetGeomMapManager()->GetGPos(CID_DEF,0,bpcpos);
	beam -> bpcpos(bpcpos.Z(),bpcposX,bpcposY);
	bpcpos.SetX(bpcposX); bpcpos.SetY(bpcposY);
	double lv;
	if(confFile->GetGeomMapManager()->GetLightVelocity(CID_DEF,hit->seg(),lv)){
		double calctsub = -bpcpos.Y() / lv;
		tmpoffset = ctsub - calctsub;
	}
	return tmpoffset;
}

bool MyAnalysisDEF::Initialize(ConfMan* confMan)
{
	std::cout << "### MyAnalysisDEF::Initialize ###" << std::endl;

	confFile = confMan;

	std::string ofname = confMan->GetOutFileName();
	ofname.insert(ofname.find(".root"),"_anaDEF");

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

	/* DEF */ 
	std::cout << "Define Histograms for DEF" << std::endl;
	int NumOfDEFSegments = 70;
	for( int seg=1; seg<=NumOfDEFSegments; seg++ ){
		new TH1F( Form("DEFu%d_ADC",seg),   Form("ADC DEFU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("DEFd%d_ADC",seg),   Form("ADC DEFD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("DEFu%d_dE",seg),   Form("dE DEFU%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("DEFd%d_dE",seg),   Form("dE DEFD%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("DEF%d_dEMean",seg),   Form("Mean dE DEF%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("DEFu%d_TDC",seg),   Form("TDC DEFU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("DEFd%d_TDC",seg),   Form("TDC DEFD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("DEFu%d_Time",seg),   Form("Time DEFU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEFd%d_Time",seg),   Form("Time DEFD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEF%d_TimeMean",seg),   Form("Mean Time DEF%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEF%d_TimeSub",seg),   Form("Subtract Time DEF%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEFu%d_CTime",seg),   Form("CTime DEFU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEFd%d_CTime",seg),   Form("CTime DEFD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEF%d_CTimeMean",seg),   Form("Mean CTime DEF%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEF%d_CTimeSub",seg),   Form("Subtract CTime DEF%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEF%d_Position",seg),   Form("Hit position DEF%d;Position (cm);Counts",seg),    1000,    -50, 50 );
		new TH1F( Form("DEFu%d_ADCwT",seg), Form("ADC wTDC DEFU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
		new TH1F( Form("DEFd%d_ADCwT",seg), Form("ADC wTDC DEFD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
		new TH2F( Form("DEFu%d_AvT",seg),   Form("ADC TDC corr. DEFU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
		new TH2F( Form("DEFd%d_AvT",seg),   Form("ADC TDC corr. DEFD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
		new TH2F( Form("DEFu%d_dEvTime",seg),   Form("dE Time corr. DEFU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("DEFd%d_dEvTime",seg),   Form("dE Time corr. DEFD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("DEFu%d_dEvCTime",seg),   Form("dE CTime corr. DEFU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("DEFd%d_dEvCTime",seg),   Form("dE CTime corr. DEFD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH1F( Form("DEF%d_TOFT0",seg),   Form("TOF DEF%d - T0;TOF (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEF%d_TimeOffset",seg),   Form("Time offset DEF%d;Time offset (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("DEF%d_TimeSubOffset",seg),   Form("Subtract Time offset DEF%d;Subtract Time offset (ns);Counts",seg),    2000,    -100, 100 );
	}
	new TH1F( "DEF_HitPat", "Hit pattern DEF;DEF Segment;Counts", NumOfDEFSegments+1, 0, NumOfDEFSegments+1 );
	new TH1F( "DEF_Mult", "Multiplicity DEF;Multiplicity;Counts", NumOfDEFSegments+1, 0, NumOfDEFSegments+1 );
	new TH2F( Form("DEF_HitPosition"),   Form("Position DEF;X hit position (cm);Y hit position (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("BPC_HitPosition"),   Form("BPC track position at BPC;X track position at DEF (cm);Y track position at DEF (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("DEFvsBPC_HitPositionatDEFX"),   Form("Position at DEF;DEF X hit position (cm);BPC track X position at DEF (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("DEFvsBPC_HitPositionatDEFY"),   Form("Position at DEF;DEF Y hit position (cm);BPC track Y position at DEF (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("DEFvsBPC_HitPositionDifference"),   Form("Position difference between DEF and BPC;X position difference (cm);Y position difference (cm)"), 80, -20, 20, 80, -20, 20 );
	new TH2F( Form("DEFvsBPC_HitPositionDifferenceXvsDirX"),   Form("Position difference(DEF-BPC) vs. X direction;X position difference (cm);X track direction"), 80, -20, 20, 100, -0.05, 0.05 );
	new TH2F( Form("DEFvsBPC_HitPositionDifferenceXvsDirY"),   Form("Position difference(DEF-BPC) vs. X direction;X position difference (cm);Y track direction"), 80, -20, 20, 100, -0.05, 0.05 );
	new TH2F( Form("DEFvsBPC_HitPositionDifferenceYvsDirX"),   Form("Position difference(DEF-BPC) vs. Y direction;Y position difference (cm);X track direction"), 80, -20, 20, 100, -0.05, 0.05 );
	new TH2F( Form("DEFvsBPC_HitPositionDifferenceYvsDirY"),   Form("Position difference(DEF-BPC) vs. Y direction;Y position difference (cm);Y track direction"), 80, -20, 20, 100, -0.05, 0.05 );

	return true;
}
