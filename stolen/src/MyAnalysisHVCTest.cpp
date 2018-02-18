// MyAnalysisHVCTest.cpp

#include "MyAnalysisHVCTest.h"

MyAnalysisHVCTest::MyAnalysisHVCTest(TFile* rt, ConfMan* conf)
{
	Initialize(conf);
	CutCondition();
	Clear();
}

MyAnalysisHVCTest::~MyAnalysisHVCTest()
{
	Clear();
	rtFile->cd();
	rtFile->Write();
	rtFile->Close();
}

void MyAnalysisHVCTest::Clear()
{
}

bool MyAnalysisHVCTest::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{

	rtFile->cd();

	/* Event selection */
  if(blMan->nHVC2()<1){
    return true;
  }

	// Trigger Pattern

	// Fill HVC
	int nHVC=0;
	for( int i=0; i<blMan->nHVC2(); i++ ){
		HodoscopeLikeHit *hit = blMan->HVC2(i);
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
		double tsub = hit->tsub();
		double ctsub = hit->ctsub();
		//double timeoffs = CalcTimeOffs(hit, particle);
		//double tsuboffs = CalcTimeSubOffs(hit, particle);

		FillHist(Form("HVCu%d_ADC",seg),au);
		FillHist(Form("HVCd%d_ADC",seg),ad);
		if(tu>0){
			FillHist(Form("HVCu%d_TDC",seg),tu);
			FillHist(Form("HVCu%d_Time",seg),timeu);
			FillHist(Form("HVCu%d_CTime",seg),ctimeu);
		}
		if(td>0){
			FillHist(Form("HVCd%d_TDC",seg),td);
			FillHist(Form("HVCd%d_Time",seg),timed);
			FillHist(Form("HVCd%d_CTime",seg),ctimed);
		}
		if( hit->CheckRange() ){
			FillHist(Form("HVCu%d_ADCwT",seg),au);
			FillHist(Form("HVCd%d_ADCwT",seg),ad);
			FillHist(Form("HVCu%d_dE",seg),eneu);
			FillHist(Form("HVCd%d_dE",seg),ened);
			FillHist(Form("HVC%d_dEMean",seg),emean);
			FillHist(Form("HVC%d_TimeMean",seg),tmean);
			FillHist(Form("HVC%d_TimeSub",seg),tsub);
			FillHist(Form("HVC%d_CTimeMean",seg),ctmean);
			FillHist(Form("HVC%d_CTimeSub",seg),ctsub);
			FillHist(Form("HVC%d_Position",seg),hitpos);
			FillHist(Form("HVCu%d_AvT",seg),au,tu);
			FillHist(Form("HVCd%d_AvT",seg),ad,td);
			FillHist(Form("HVCu%d_dEvTime",seg),eneu,timeu);
			FillHist(Form("HVCd%d_dEvTime",seg),ened,timed);
			FillHist(Form("HVCu%d_dEvCTime",seg),eneu,ctimeu);
			FillHist(Form("HVCd%d_dEvCTime",seg),ened,ctimed);
			FillHist(Form("HVC_HitPosition"),pos.X(),pos.Y());
			FillHist("HVC_HitPat",seg);
			nHVC++;
		}
	}
	FillHist("HVC_Mult",nHVC);


}

bool MyAnalysisHVCTest::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHVCTest::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHVCTest::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHVCTest::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHVCTest::CutCondition()
{
}

double MyAnalysisHVCTest::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
	double tmpoffset = -9999.9;
	if(!hit->CheckRange()) return tmpoffset;
	pBeam* beam = particle->beam(0);
	double tof = hit->ctmean() - beam->t0time();

	TVector3 bpdpos; double bpdposX, bpdposY;
	confFile->GetGeomMapManager()->GetGPos(CID_HVC2,0,bpdpos);
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

double MyAnalysisHVCTest::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
	double tmpoffset = -9999.9;
	if(!hit->CheckRange()) return tmpoffset;
	pBeam* beam = particle->beam(0);
	double ctsub = hit->ctsub();

	TVector3 bpcpos; double bpcposX, bpcposY;
	confFile->GetGeomMapManager()->GetGPos(CID_HVC2,0,bpcpos);
	beam -> bpcpos(bpcpos.Z(),bpcposX,bpcposY);
	bpcpos.SetX(bpcposX); bpcpos.SetY(bpcposY);
	double lv;
	if(confFile->GetGeomMapManager()->GetLightVelocity(CID_HVC2,hit->seg(),lv)){
		double calctsub = -bpcpos.Y() / lv;
		tmpoffset = ctsub - calctsub;
	}
	return tmpoffset;
}

bool MyAnalysisHVCTest::Initialize(ConfMan* confMan)
{
	std::cout << "### MyAnalysisHVCTest::Initialize ###" << std::endl;

	confFile = confMan;

	std::string ofname = confMan->GetOutFileName();
	ofname.insert(ofname.find(".root"),"_anaHVC");

	rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
	rtFile -> cd();

	/* HVC */ 
	std::cout << "Define Histograms for HVC" << std::endl;
	int NumOfHVCSegments = 8;
	for( int seg=1; seg<=NumOfHVCSegments; seg++ ){
		new TH1F( Form("HVCu%d_ADC",seg),   Form("ADC HVCU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("HVCd%d_ADC",seg),   Form("ADC HVCD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("HVCu%d_dE",seg),   Form("dE HVCU%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("HVCd%d_dE",seg),   Form("dE HVCD%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("HVC%d_dEMean",seg),   Form("Mean dE HVC%d;dE (MeV);Counts",seg),    200,    0, 20 );
		new TH1F( Form("HVCu%d_TDC",seg),   Form("TDC HVCU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("HVCd%d_TDC",seg),   Form("TDC HVCD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
		new TH1F( Form("HVCu%d_Time",seg),   Form("Time HVCU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("HVCd%d_Time",seg),   Form("Time HVCD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("HVC%d_TimeMean",seg),   Form("Mean Time HVC%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("HVC%d_TimeSub",seg),   Form("Subtract Time HVC%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("HVCu%d_CTime",seg),   Form("CTime HVCU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("HVCd%d_CTime",seg),   Form("CTime HVCD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("HVC%d_CTimeMean",seg),   Form("Mean CTime HVC%d;Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("HVC%d_CTimeSub",seg),   Form("Subtract CTime HVC%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
		new TH1F( Form("HVC%d_Position",seg),   Form("Hit position HVC%d;Position (cm);Counts",seg),    1000,    -50, 50 );
		new TH1F( Form("HVCu%d_ADCwT",seg), Form("ADC wTDC HVCU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
		new TH1F( Form("HVCd%d_ADCwT",seg), Form("ADC wTDC HVCD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
		new TH2F( Form("HVCu%d_AvT",seg),   Form("ADC TDC corr. HVCU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
		new TH2F( Form("HVCd%d_AvT",seg),   Form("ADC TDC corr. HVCD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
		new TH2F( Form("HVCu%d_dEvTime",seg),   Form("dE Time corr. HVCU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("HVCd%d_dEvTime",seg),   Form("dE Time corr. HVCD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("HVCu%d_dEvCTime",seg),   Form("dE CTime corr. HVCU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
		new TH2F( Form("HVCd%d_dEvCTime",seg),   Form("dE CTime corr. HVCD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
	}
	new TH1F( "HVC_HitPat", "Hit pattern HVC;HVC Segment;Counts", NumOfHVCSegments+1, 0, NumOfHVCSegments+1 );
	new TH1F( "HVC_Mult", "Multiplicity HVC;Multiplicity;Counts", NumOfHVCSegments+1, 0, NumOfHVCSegments+1 );
	new TH2F( Form("HVC_HitPosition"),   Form("Position HVC;X hit position (cm);Y hit position (cm)"), 80, -20, 20, 80, -20, 20 );

	return true;
}
