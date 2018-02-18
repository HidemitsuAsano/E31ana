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

class EventAnalysisMyIH: public EventTemp
{
	public:
		EventAnalysisMyIH();
		~EventAnalysisMyIH();
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

	EventAnalysisMyIH::EventAnalysisMyIH()
: EventTemp()
{
}

EventAnalysisMyIH::~EventAnalysisMyIH()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyIH::Initialize( ConfMan *conf )
{
#if 0
	std::cout << " Enter EventAnalysisMyIH::Initialize " << std::endl;
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

void EventAnalysisMyIH::USca( int nsca, unsigned int *sca )
{
#if 0
	std::cout << " Enter EventAnalysisMyIH::USca " << std::endl;
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

bool EventAnalysisMyIH::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyIH::UAna " << std::endl;
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

	// Multihit of IH and IH //
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
	//if(MulIH==0){
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
		// IH
		int nIH=0;
		int HitSeg[2]={-1,-1};
		for( int i=0; i<cdsMan->nIH(); i++ ){
			HodoscopeLikeHit *hit = cdsMan->IH(i);
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

			h1 = (TH1F*)gFile->Get( Form("IHu%d_ADC_OnBeam",seg) ); h1->Fill( au );
			h1 = (TH1F*)gFile->Get( Form("IHd%d_ADC_OnBeam",seg) ); h1->Fill( ad );
			if(tu>0){
				h1 = (TH1F*)gFile->Get( Form("IHu%d_TDC_OnBeam",seg) ); h1->Fill( tu );
				h1 = (TH1F*)gFile->Get( Form("IHu%d_Time_OnBeam",seg) ); h1->Fill( timeu );
				h1 = (TH1F*)gFile->Get( Form("IHu%d_CTime_OnBeam",seg) ); h1->Fill( ctimeu );
			}
			if(td>0){
				h1 = (TH1F*)gFile->Get( Form("IHd%d_TDC_OnBeam",seg) ); h1->Fill( td );
				h1 = (TH1F*)gFile->Get( Form("IHd%d_Time_OnBeam",seg) ); h1->Fill( timed );
				h1 = (TH1F*)gFile->Get( Form("IHd%d_CTime_OnBeam",seg) ); h1->Fill( ctimed );
			}
			if( hit->CheckRange() ){
				h1 = (TH1F*)gFile->Get( Form("IHu%d_ADCwT_OnBeam",seg) ); h1->Fill( au );
				h1 = (TH1F*)gFile->Get( Form("IHd%d_ADCwT_OnBeam",seg) ); h1->Fill( ad );
				h1 = (TH1F*)gFile->Get( Form("IHu%d_dE_OnBeam",seg) ); h1->Fill( eneu );
				h1 = (TH1F*)gFile->Get( Form("IHd%d_dE_OnBeam",seg) ); h1->Fill( ened );
				h1 = (TH1F*)gFile->Get( Form("IH%d_dEMean_OnBeam",seg) ); h1->Fill( emean );
				h1 = (TH1F*)gFile->Get( Form("IH%d_TimeMean_OnBeam",seg) ); h1->Fill( tmean );
				h1 = (TH1F*)gFile->Get( Form("IH%d_CTimeMean_OnBeam",seg) ); h1->Fill( ctmean );
				h1 = (TH1F*)gFile->Get( Form("IH%d_Position_OnBeam",seg) ); h1->Fill( hitpos );
				//if((eneu<1.2||eneu>2.8)&&(ened<1.2||ened>2.8)){
				h2 = (TH2F*)gFile->Get( Form("IHu%d_AvT_OnBeam",seg) );  h2->Fill( au, tu );
				h2 = (TH2F*)gFile->Get( Form("IHd%d_AvT_OnBeam",seg) );  h2->Fill( ad, td );
				h2 = (TH2F*)gFile->Get( Form("IHu%d_dEvTime_OnBeam",seg) );  h2->Fill( eneu, timeu );
				h2 = (TH2F*)gFile->Get( Form("IHd%d_dEvTime_OnBeam",seg) );  h2->Fill( ened, timed );
				h2 = (TH2F*)gFile->Get( Form("IHu%d_dEvCTime_OnBeam",seg) );  h2->Fill( eneu, ctimeu );
				h2 = (TH2F*)gFile->Get( Form("IHd%d_dEvCTime_OnBeam",seg) );  h2->Fill( ened, ctimed );
				h2 = (TH2F*)gFile->Get( Form("IHu%d_dEvCTimeMean_OnBeam",seg) );  h2->Fill( eneu, ctmean );
				h2 = (TH2F*)gFile->Get( Form("IHd%d_dEvCTimeMean_OnBeam",seg) );  h2->Fill( ened, ctmean );
				//}
				h2 = (TH2F*)gFile->Get( Form("IH%d_dEMeanvCTimeMean_OnBeam",seg) );  h2->Fill( emean, ctmean );
				h1 = (TH1F*)gFile->Get( "IH_HitPat_OnBeam" ); h1->Fill( seg );
				nIH++;
			}
		}
		h1 = (TH1F*)gFile->Get( "IH_Mult_OnBeam" ); h1->Fill( nIH );
	}
	/* On Beam */

	/* Off beam */
	if(header->IsTrig(Trig_Cosmic)){
		// IH
		int nIH=0;
		int HitSeg[2]={-1,-1};
		for( int i=0; i<cdsMan->nIH(); i++ ){
			HodoscopeLikeHit *hit = cdsMan->IH(i);
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

			h1 = (TH1F*)gFile->Get( Form("IHu%d_ADC_OffBeam",seg) ); h1->Fill( au );
			h1 = (TH1F*)gFile->Get( Form("IHd%d_ADC_OffBeam",seg) ); h1->Fill( ad );
			if(tu>0){
				h1 = (TH1F*)gFile->Get( Form("IHu%d_TDC_OffBeam",seg) ); h1->Fill( tu );
				h1 = (TH1F*)gFile->Get( Form("IHu%d_Time_OffBeam",seg) ); h1->Fill( timeu );
				h1 = (TH1F*)gFile->Get( Form("IHu%d_CTime_OffBeam",seg) ); h1->Fill( ctimeu );
			}
			if(td>0){
				h1 = (TH1F*)gFile->Get( Form("IHd%d_TDC_OffBeam",seg) ); h1->Fill( td );
				h1 = (TH1F*)gFile->Get( Form("IHd%d_Time_OffBeam",seg) ); h1->Fill( timed );
				h1 = (TH1F*)gFile->Get( Form("IHd%d_CTime_OffBeam",seg) ); h1->Fill( ctimed );
			}
			if( hit->CheckRange() ){
				h1 = (TH1F*)gFile->Get( Form("IHu%d_ADCwT_OffBeam",seg) ); h1->Fill( au );
				h1 = (TH1F*)gFile->Get( Form("IHd%d_ADCwT_OffBeam",seg) ); h1->Fill( ad );
				h1 = (TH1F*)gFile->Get( Form("IHu%d_dE_OffBeam",seg) ); h1->Fill( eneu );
				h1 = (TH1F*)gFile->Get( Form("IHd%d_dE_OffBeam",seg) ); h1->Fill( ened );
				h1 = (TH1F*)gFile->Get( Form("IH%d_dEMean_OffBeam",seg) ); h1->Fill( emean );
				h1 = (TH1F*)gFile->Get( Form("IH%d_TimeMean_OffBeam",seg) ); h1->Fill( tmean );
				h1 = (TH1F*)gFile->Get( Form("IH%d_CTimeMean_OffBeam",seg) ); h1->Fill( ctmean );
				h1 = (TH1F*)gFile->Get( Form("IH%d_Position_OffBeam",seg) ); h1->Fill( hitpos );
				if((eneu<1.2||eneu>2.8)&&(ened<1.2||ened>2.8)){
				h2 = (TH2F*)gFile->Get( Form("IHu%d_AvT_OffBeam",seg) );  h2->Fill( au, tu );
				h2 = (TH2F*)gFile->Get( Form("IHd%d_AvT_OffBeam",seg) );  h2->Fill( ad, td );
				h2 = (TH2F*)gFile->Get( Form("IHu%d_dEvTime_OffBeam",seg) );  h2->Fill( eneu, timeu );
				h2 = (TH2F*)gFile->Get( Form("IHd%d_dEvTime_OffBeam",seg) );  h2->Fill( ened, timed );
				h2 = (TH2F*)gFile->Get( Form("IHu%d_dEvCTime_OffBeam",seg) );  h2->Fill( eneu, ctimeu );
				h2 = (TH2F*)gFile->Get( Form("IHd%d_dEvCTime_OffBeam",seg) );  h2->Fill( ened, ctimed );
				h2 = (TH2F*)gFile->Get( Form("IHu%d_dEvCTimeMean_OffBeam",seg) );  h2->Fill( eneu, ctmean );
				h2 = (TH2F*)gFile->Get( Form("IHd%d_dEvCTimeMean_OffBeam",seg) );  h2->Fill( ened, ctmean );
				}
				h2 = (TH2F*)gFile->Get( Form("IH%d_dEMeanvCTimeMean_OffBeam",seg) );  h2->Fill( emean, ctmean );
				h1 = (TH1F*)gFile->Get( "IH_HitPat_OffBeam" ); h1->Fill( seg );
				nIH++;
			}
		}
		h1 = (TH1F*)gFile->Get( "IH_Mult_OffBeam" ); h1->Fill( nIH );
	}
	/* Off Beam */

	header->Clear();
	blMan->Clear();
	cdsMan->Clear();
	return true;
}

void EventAnalysisMyIH::Finalize()
{
	std::cout << " Enter EventAnalysisMyIH::Finalize " << std::endl;

	rtFile->cd();
	gFile->Write();
	gFile->Close();

	delete blMan;
	delete cdsMan;
	delete header;
}

void EventAnalysisMyIH::InitializeHistogram()
{
	Int_t NumOfIHSegments = 36;

	rtFile->cd();

	// Scaler
	new TH1F( "Scaler", "Scaler", 50, 0, 50 );

	// Trigger Pattern
	new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

	// IH
	std::cout << "Define Histograms for IH" << std::endl;
	TString beam[2] = {"OnBeam","OffBeam"};
	for( int i=0; i<2; i++ ){
		for( int seg=1; seg<=NumOfIHSegments; seg++ ){
			new TH1F( Form("IHu%d_ADC_%s",seg,beam[i].Data()),   Form("ADC IHU%d (%s);ADC ch.;Counts",seg,beam[i].Data()),    4000,    0, 4000 );
			new TH1F( Form("IHd%d_ADC_%s",seg,beam[i].Data()),   Form("ADC IHD%d (%s);ADC ch.;Counts",seg,beam[i].Data()),    4000,    0, 4000 );
			new TH1F( Form("IHu%d_dE_%s",seg,beam[i].Data()),   Form("dE IHU%d (%s);dE (MeV);Counts",seg,beam[i].Data()),    400,    0, 40 );
			new TH1F( Form("IHd%d_dE_%s",seg,beam[i].Data()),   Form("dE IHD%d (%s);dE (MeV);Counts",seg,beam[i].Data()),    400,    0, 40 );
			new TH1F( Form("IH%d_dEMean_%s",seg,beam[i].Data()),   Form("Mean dE IH%d (%s);dE (MeV);Counts",seg,beam[i].Data()),    400,    0, 40 );
			new TH1F( Form("IHu%d_TDC_%s",seg,beam[i].Data()),   Form("TDC IHU%d (%s);TDC ch.;Counts",seg,beam[i].Data()),    4000,    0, 4000 );
			new TH1F( Form("IHd%d_TDC_%s",seg,beam[i].Data()),   Form("TDC IHD%d (%s);TDC ch.;Counts",seg,beam[i].Data()),    4000,    0, 4000 );
			new TH1F( Form("IHu%d_Time_%s",seg,beam[i].Data()),   Form("Time IHU%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("IHd%d_Time_%s",seg,beam[i].Data()),   Form("Time IHD%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("IH%d_TimeMean_%s",seg,beam[i].Data()),   Form("Mean Time IH%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("IHu%d_CTime_%s",seg,beam[i].Data()),   Form("CTime IHU%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("IHd%d_CTime_%s",seg,beam[i].Data()),   Form("CTime IHD%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("IH%d_CTimeMean_%s",seg,beam[i].Data()),   Form("Mean CTime IH%d (%s);Time (ns);Counts",seg,beam[i].Data()),    4000,    -100, 100 );
			new TH1F( Form("IH%d_Position_%s",seg,beam[i].Data()),   Form("Hit position IH%d (%s);Position (cm);Counts",seg,beam[i].Data()),    500,    -50, 50 );
			new TH1F( Form("IHu%d_ADCwT_%s",seg,beam[i].Data()), Form("ADC wTDC IHU%d (%s);ADC ch.;Counts",seg,beam[i].Data()),  4000,    0, 4000 );
			new TH1F( Form("IHd%d_ADCwT_%s",seg,beam[i].Data()), Form("ADC wTDC IHD%d (%s);ADC ch.;Counts",seg,beam[i].Data()),  4000,    0, 4000 );
			new TH2F( Form("IHu%d_AvT_%s",seg,beam[i].Data()),   Form("ADC TDC corr. IHU%d (%s);ADC ch.;TDC ch.",seg,beam[i].Data()),     200,    0, 4000,  200,    0, 4000 );
			new TH2F( Form("IHd%d_AvT_%s",seg,beam[i].Data()),   Form("ADC TDC corr. IHD%d (%s);ADC ch.;TDC ch.",seg,beam[i].Data()),     200,    0, 4000,  200,    0, 4000 );
			new TH2F( Form("IHu%d_dEvTime_%s",seg,beam[i].Data()),   Form("dE Time corr. IHU%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("IHd%d_dEvTime_%s",seg,beam[i].Data()),   Form("dE Time corr. IHD%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("IHu%d_dEvCTime_%s",seg,beam[i].Data()),   Form("dE CTime corr. IHU%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("IHd%d_dEvCTime_%s",seg,beam[i].Data()),   Form("dE CTime corr. IHD%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("IHu%d_dEvCTimeMean_%s",seg,beam[i].Data()),   Form("dE CTimeMean corr. IHU%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("IHd%d_dEvCTimeMean_%s",seg,beam[i].Data()),   Form("dE CTimeMean corr. IHD%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
			new TH2F( Form("IH%d_dEMeanvCTimeMean_%s",seg,beam[i].Data()),   Form("dEMean CTimeMean corr. IH%d (%s);dE (MeV);Time (ns)",seg,beam[i].Data()),     400,    0, 40,  4000,    -100, 100 );
		}
		new TH1F( Form("IH_Mult_%s",beam[i].Data()), Form("Multiplicity IH (%s);Multiplicity;Counts",beam[i].Data()), NumOfIHSegments+1, 0, NumOfIHSegments+1 );
		new TH1F( Form("IH_HitPat_%s",beam[i].Data()), Form("Hit Pattern IH (%s);IH hit segment;Counts",beam[i].Data()), NumOfIHSegments+1, 0, NumOfIHSegments+1 );
	}

}

EventTemp *EventAlloc::EventAllocator()
{
	EventAnalysisMyIH *event = new EventAnalysisMyIH();
	return (EventTemp*)event;
}
