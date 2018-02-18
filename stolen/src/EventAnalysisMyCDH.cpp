#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#define Debug 0

class EventAnalysisMyCDH: public EventTemp
{
	public:
		EventAnalysisMyCDH();
		~EventAnalysisMyCDH();
	private:
		TFile *rtFile;
		TTree *evTree;
		TTree *scaTree;
		CDSHitMan *cdsMan;
		BeamLineHitMan *blMan;
		EventHeader *header;
		ScalerMan *scaMan;
	public:
		void Initialize( ConfMan *conf );
		void USca( int nsca, unsigned int *sca );
		bool UAna( TKOHitCollection *tko );
		void Finalize();

		void InitializeHistogram();
};

	EventAnalysisMyCDH::EventAnalysisMyCDH()
: EventTemp()
{
}

EventAnalysisMyCDH::~EventAnalysisMyCDH()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyCDH::Initialize( ConfMan *conf )
{
#if 0
	std::cout << " Enter EventAnalysisMyCDH::Initialize " << std::endl;
#endif
	confMan = conf;

	rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
	InitializeHistogram();

	header = new EventHeader();
	if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
	cdsMan = new CDSHitMan();
	if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
	blMan = new BeamLineHitMan();
	if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
	scaMan = new ScalerMan();
	if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
}

void EventAnalysisMyCDH::USca( int nsca, unsigned int *sca )
{
#if 0
	std::cout << " Enter EventAnalysisMyCDH::USca " << std::endl;
#endif
	Block_Event_Number++;
	header->SetBlockEventNumber( Block_Event_Number );
	scaMan->SetBlockEventNumber( Block_Event_Number );
	for( int i=0; i<nsca; i++ ){
		scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
	}
#if 0
	std::cout << nsca << std::endl;
	for( int i=0; i<nsca; i++ ){
		std::cout << "  " << sca[i];
	}
	std::cout << std::endl;
#endif

	TH1F *h1;
	for( int i=0; i<scaMan->nsca(); i++ ){
		int val = scaMan->sca(i)->val();
		TString name = scaMan->sca(i)->name();
		h1 = (TH1F*)gFile->Get("Scaler"); h1->Fill(i,val);
	}


	scaMan->Clear();
}

bool EventAnalysisMyCDH::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyCDH::UAna " << std::endl;
#endif

	Event_Number++;
	{ int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
		if( status==1 ) return true;
		if( status==2 ) return false; }

	if( Event_Number%5000==0 )
		std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;

	header->SetRunNumber(0);
	header->SetEventNumber(Event_Number);
	header->Convert( tko, confMan );
	blMan->Convert( tko, confMan );
	cdsMan->Convert( tko, confMan );

	rtFile->cd();
	TH1F *h1;
	TH2F *h2;

	// ======== //
	// Raw Data //
	// ======== //

	// Multihit of CDH and IH //
	int MulT0=0;	
	for(int i=0; i<blMan->nT0(); i++){
		HodoscopeLikeHit* hit = blMan->T0(i);
		if(hit->CheckRange()) MulT0++;
	}
	int MulCDH=0;	
	for(int i=0; i<cdsMan->nCDH(); i++){
		HodoscopeLikeHit* hit = cdsMan->CDH(i);
		if(hit->CheckRange()) MulCDH++;
	}
	int MulIH=0;	
	for(int i=0; i<cdsMan->nIH(); i++){
		HodoscopeLikeHit* hit = cdsMan->IH(i);
		if(hit->CheckRange()) MulIH++;
	}

	// Selection //
	//if(MulCDH==0){
	//	header->Clear();
	//	blMan->Clear();
	//	cdsMan->Clear();
	//	return true;
	//}

	// Trigger Pattern //
	for( int i=0; i<20; i++ ){
		int val = header->pattern(i);
		if( 0<val ){
			h1 = (TH1F*)gFile->Get("Pattern"); h1->Fill(i);
		}
	}
	///* ON beam */
	if(MulT0==1 && !header->IsTrig(Trig_Cosmic)){
		// CDH
		int nCDH=0;
		int HitSeg[2]={-1,-1};
		for( int i=0; i<cdsMan->nCDH(); i++ ){
			HodoscopeLikeHit *hit = cdsMan->CDH(i);
			int seg = hit->seg();
			int au = hit->adcu(), ad = hit->adcd();
			int tu = hit->tdcu(), td = hit->tdcd();
			double hitpos = hit->hitpos();
			double timeu = hit->time(0), timed = hit->time(1);
			double ctimeu = hit->ctime(0), ctimed = hit->ctime(1);
			double eneu = hit->ene(0), ened = hit->ene(1);
			double emean = hit->emean();
			double tmean = hit->tmean();
			double ctmean = hit->ctmean();

			h1 = (TH1F*)gFile->Get( Form("CDHu%d_ADC_OnBeam",seg) ); h1->Fill( au );
			h1 = (TH1F*)gFile->Get( Form("CDHd%d_ADC_OnBeam",seg) ); h1->Fill( ad );
			if(tu>0){
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_TDC_OnBeam",seg) ); h1->Fill( tu );
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_Time_OnBeam",seg) ); h1->Fill( timeu );
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_CTime_OnBeam",seg) ); h1->Fill( ctimeu );
			}
			if(td>0){
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_TDC_OnBeam",seg) ); h1->Fill( td );
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_Time_OnBeam",seg) ); h1->Fill( timed );
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_CTime_OnBeam",seg) ); h1->Fill( ctimed );
			}
			if( hit->CheckRange() ){
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_ADCwT_OnBeam",seg) ); h1->Fill( au );
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_ADCwT_OnBeam",seg) ); h1->Fill( ad );
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_dE_OnBeam",seg) ); h1->Fill( eneu );
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_dE_OnBeam",seg) ); h1->Fill( ened );
				h1 = (TH1F*)gFile->Get( Form("CDH%d_dEMean_OnBeam",seg) ); h1->Fill( emean );
				h1 = (TH1F*)gFile->Get( Form("CDH%d_TimeMean_OnBeam",seg) ); h1->Fill( tmean );
				h1 = (TH1F*)gFile->Get( Form("CDH%d_CTimeMean_OnBeam",seg) ); h1->Fill( ctmean );
				h1 = (TH1F*)gFile->Get( Form("CDH%d_Position_OnBeam",seg) ); h1->Fill( hitpos );
				//if((eneu<1.2||eneu>2.8)&&(ened<1.2||ened>2.8)){
				h2 = (TH2F*)gFile->Get( Form("CDHu%d_AvT_OnBeam",seg) );  h2->Fill( au, tu );
				h2 = (TH2F*)gFile->Get( Form("CDHd%d_AvT_OnBeam",seg) );  h2->Fill( ad, td );
				h2 = (TH2F*)gFile->Get( Form("CDHu%d_dEvTime_OnBeam",seg) );  h2->Fill( eneu, timeu );
				h2 = (TH2F*)gFile->Get( Form("CDHd%d_dEvTime_OnBeam",seg) );  h2->Fill( ened, timed );
				h2 = (TH2F*)gFile->Get( Form("CDHu%d_dEvCTime_OnBeam",seg) );  h2->Fill( eneu, ctimeu );
				h2 = (TH2F*)gFile->Get( Form("CDHd%d_dEvCTime_OnBeam",seg) );  h2->Fill( ened, ctimed );
				h2 = (TH2F*)gFile->Get( Form("CDHu%d_dEvCTimeMean_OnBeam",seg) );  h2->Fill( eneu, ctmean );
				h2 = (TH2F*)gFile->Get( Form("CDHd%d_dEvCTimeMean_OnBeam",seg) );  h2->Fill( ened, ctmean );
				//}
				h2 = (TH2F*)gFile->Get( Form("CDH%d_dEMeanvCTimeMean_OnBeam",seg) );  h2->Fill( emean, ctmean );
				h1 = (TH1F*)gFile->Get( "CDH_HitPat_OnBeam" ); h1->Fill( seg );
				nCDH++;
			}
		}
		h1 = (TH1F*)gFile->Get( "CDH_Mult_OnBeam" ); h1->Fill( nCDH );
	}
	/* On Beam */

	/* Off beam */
	if(header->IsTrig(Trig_Cosmic)){
		// CDH
		int nCDH=0;
		int HitSeg[2]={-1,-1};
		for( int i=0; i<cdsMan->nCDH(); i++ ){
			HodoscopeLikeHit *hit = cdsMan->CDH(i);
			int seg = hit->seg();
			int au = hit->adcu(), ad = hit->adcd();
			int tu = hit->tdcu(), td = hit->tdcd();
			double hitpos = hit->hitpos();
			double timeu = hit->time(0), timed = hit->time(1);
			double ctimeu = hit->ctime(0), ctimed = hit->ctime(1);
			double eneu = hit->ene(0), ened = hit->ene(1);
			double emean = hit->emean();
			double tmean = hit->tmean();
			double ctmean = hit->ctmean();

			h1 = (TH1F*)gFile->Get( Form("CDHu%d_ADC_OffBeam",seg) ); h1->Fill( au );
			h1 = (TH1F*)gFile->Get( Form("CDHd%d_ADC_OffBeam",seg) ); h1->Fill( ad );
			if(tu>0){
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_TDC_OffBeam",seg) ); h1->Fill( tu );
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_Time_OffBeam",seg) ); h1->Fill( timeu );
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_CTime_OffBeam",seg) ); h1->Fill( ctimeu );
			}
			if(td>0){
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_TDC_OffBeam",seg) ); h1->Fill( td );
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_Time_OffBeam",seg) ); h1->Fill( timed );
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_CTime_OffBeam",seg) ); h1->Fill( ctimed );
			}
			if( hit->CheckRange() ){
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_ADCwT_OffBeam",seg) ); h1->Fill( au );
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_ADCwT_OffBeam",seg) ); h1->Fill( ad );
				h1 = (TH1F*)gFile->Get( Form("CDHu%d_dE_OffBeam",seg) ); h1->Fill( eneu );
				h1 = (TH1F*)gFile->Get( Form("CDHd%d_dE_OffBeam",seg) ); h1->Fill( ened );
				h1 = (TH1F*)gFile->Get( Form("CDH%d_dEMean_OffBeam",seg) ); h1->Fill( emean );
				h1 = (TH1F*)gFile->Get( Form("CDH%d_TimeMean_OffBeam",seg) ); h1->Fill( tmean );
				h1 = (TH1F*)gFile->Get( Form("CDH%d_CTimeMean_OffBeam",seg) ); h1->Fill( ctmean );
				h1 = (TH1F*)gFile->Get( Form("CDH%d_Position_OffBeam",seg) ); h1->Fill( hitpos );
				if((eneu<1.2||eneu>2.8)&&(ened<1.2||ened>2.8)){
				h2 = (TH2F*)gFile->Get( Form("CDHu%d_AvT_OffBeam",seg) );  h2->Fill( au, tu );
				h2 = (TH2F*)gFile->Get( Form("CDHd%d_AvT_OffBeam",seg) );  h2->Fill( ad, td );
				h2 = (TH2F*)gFile->Get( Form("CDHu%d_dEvTime_OffBeam",seg) );  h2->Fill( eneu, timeu );
				h2 = (TH2F*)gFile->Get( Form("CDHd%d_dEvTime_OffBeam",seg) );  h2->Fill( ened, timed );
				h2 = (TH2F*)gFile->Get( Form("CDHu%d_dEvCTime_OffBeam",seg) );  h2->Fill( eneu, ctimeu );
				h2 = (TH2F*)gFile->Get( Form("CDHd%d_dEvCTime_OffBeam",seg) );  h2->Fill( ened, ctimed );
				h2 = (TH2F*)gFile->Get( Form("CDHu%d_dEvCTimeMean_OffBeam",seg) );  h2->Fill( eneu, ctmean );
				h2 = (TH2F*)gFile->Get( Form("CDHd%d_dEvCTimeMean_OffBeam",seg) );  h2->Fill( ened, ctmean );
				}
				h2 = (TH2F*)gFile->Get( Form("CDH%d_dEMeanvCTimeMean_OffBeam",seg) );  h2->Fill( emean, ctmean );
				h1 = (TH1F*)gFile->Get( "CDH_HitPat_OffBeam" ); h1->Fill( seg );
				nCDH++;
			}
		}
		h1 = (TH1F*)gFile->Get( "CDH_Mult_OffBeam" ); h1->Fill( nCDH );
	}
	/* Off Beam */

	header->Clear();
	blMan->Clear();
	cdsMan->Clear();
	return true;
}

void EventAnalysisMyCDH::Finalize()
{
	std::cout << " Enter EventAnalysisMyCDH::Finalize " << std::endl;

	rtFile->cd();
	gFile->Write();
	gFile->Close();

	delete blMan;
	delete cdsMan;
	delete header;
}

void EventAnalysisMyCDH::InitializeHistogram()
{
	Int_t NumOfCDHSegments = 36;

	rtFile->cd();

	// Scaler
	new TH1F( "Scaler", "Scaler", 50, 0, 50 );

	// Trigger Pattern
	new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

	// CDH
	std::cout << "Define Histograms for CDH" << std::endl;
	TString beam[2] = {"OnBeam","OffBeam"};
	for( int i=0; i<2; i++ ){
		for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
			new TH1F( Form("CDHu%d_ADC_%s",seg,beam[i].Data()),   Form("ADC CDHU%d (%s);ADC ch.;Counts",seg,beam[i].Data()),    4000,    0, 4000 );
			new TH1F( Form("CDHd%d_ADC_%s",seg,beam[i].Data()),   Form("ADC CDHD%d (%s);ADC ch.;Counts",seg,beam[i].Data()),    4000,    0, 4000 );
			new TH1F( Form("CDHu%d_dE_%s",seg,beam[i].Data()),   Form("dE CDHU%d (%s);dE (MeV);Counts",seg,beam[i].Data()),    400,    0, 40 );
			new TH1F( Form("CDHd%d_dE_%s",seg,beam[i].Data()),   Form("dE CDHD%d (%s);dE (MeV);Counts",seg,beam[i].Data()),    400,    0, 40 );
			new TH1F( Form("CDH%d_dEMean_%s",seg,beam[i].Data()),   Form("Mean dE CDH%d (%s);dE (MeV);Counts",seg,beam[i].Data()),    400,    0, 40 );
			new TH1F( Form("CDHu%d_TDC_%s",seg,beam[i].Data()),   Form("TDC CDHU%d (%s);TDC ch.;Counts",seg,beam[i].Data()),    4000,    0, 4000 );
			new TH1F( Form("CDHd%d_TDC_%s",seg,beam[i].Data()),   Form("TDC CDHD%d (%s);TDC ch.;Counts",seg,beam[i].Data()),    4000,    0, 4000 );
			new TH1F( Form("CDHu%d_Time_%s",seg,beam[i].Data()),   Form("Time CDHU%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("CDHd%d_Time_%s",seg,beam[i].Data()),   Form("Time CDHD%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("CDH%d_TimeMean_%s",seg,beam[i].Data()),   Form("Mean Time CDH%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("CDHu%d_CTime_%s",seg,beam[i].Data()),   Form("CTime CDHU%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("CDHd%d_CTime_%s",seg,beam[i].Data()),   Form("CTime CDHD%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("CDH%d_CTimeMean_%s",seg,beam[i].Data()),   Form("Mean CTime CDH%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("CDH%d_Position_%s",seg,beam[i].Data()),   Form("Hit position CDH%d (%s);Position (cm);Counts",seg,beam[i].Data()),    500,    -50, 50 );
			new TH1F( Form("CDHu%d_ADCwT_%s",seg,beam[i].Data()), Form("ADC wTDC CDHU%d (%s);ADC ch.;Counts",seg,beam[i].Data()),  4000,    0, 4000 );
			new TH1F( Form("CDHd%d_ADCwT_%s",seg,beam[i].Data()), Form("ADC wTDC CDHD%d (%s);ADC ch.;Counts",seg,beam[i].Data()),  4000,    0, 4000 );
			new TH2F( Form("CDHu%d_AvT_%s",seg,beam[i].Data()),   Form("ADC TDC corr. CDHU%d (%s);ADC ch.;TDC ch.",seg,beam[i].Data()),     200,    0, 4000,  200,    0, 4000 );
			new TH2F( Form("CDHd%d_AvT_%s",seg,beam[i].Data()),   Form("ADC TDC corr. CDHD%d (%s);ADC ch.;TDC ch.",seg,beam[i].Data()),     200,    0, 4000,  200,    0, 4000 );
			new TH2F( Form("CDHu%d_dEvTime_%s",seg,beam[i].Data()),   Form("dE Time corr. CDHU%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("CDHd%d_dEvTime_%s",seg,beam[i].Data()),   Form("dE Time corr. CDHD%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("CDHu%d_dEvCTime_%s",seg,beam[i].Data()),   Form("dE CTime corr. CDHU%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("CDHd%d_dEvCTime_%s",seg,beam[i].Data()),   Form("dE CTime corr. CDHD%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("CDHu%d_dEvCTimeMean_%s",seg,beam[i].Data()),   Form("dE CTimeMean corr. CDHU%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("CDHd%d_dEvCTimeMean_%s",seg,beam[i].Data()),   Form("dE CTimeMean corr. CDHD%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("CDH%d_dEMeanvCTimeMean_%s",seg,beam[i].Data()),   Form("dEMean CTimeMean corr. CDH%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
		}
		new TH1F( Form("CDH_Mult_%s",beam[i].Data()), Form("Multiplicity CDH (%s);Multiplicity;Counts",beam[i].Data()), NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
		new TH1F( Form("CDH_HitPat_%s",beam[i].Data()), Form("Hit Pattern CDH (%s);CDH hit segment;Counts",beam[i].Data()), NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
	}

}

EventTemp *EventAlloc::EventAllocator()
{
	EventAnalysisMyCDH *event = new EventAnalysisMyCDH();
	return (EventTemp*)event;
}
