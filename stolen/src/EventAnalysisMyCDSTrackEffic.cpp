#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "Particle.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "GlobalVariables.h"
#include "MyAnalysisBL.h"

#define RUN43 0
#define TKOHIS 0
class EventAnalysisMyCDC: public EventTemp
{
	public:
		EventAnalysisMyCDC();
		~EventAnalysisMyCDC();
	private:
		TFile *rtFile;
		TFile *cdcFile;
		TTree *cdcTree;
		TTree *evTree;
		TTree *scaTree;
		BeamLineHitMan *blMan;
		BeamLineTrackMan *bltrackMan;
		CDSHitMan *cdsMan;
		CDSTrackingMan *cdstrackMan;
		ScalerMan *scaMan;
		EventHeader *header;
		EventHeader *header2;
		MyAnalysisBL *blAna;

		int AllGoodTrack;
		int nTrack;
		int CDC_Event_Number;

	public:
		void Initialize( ConfMan *conf );
		void USca( int nsca, unsigned int *sca );
		bool UAna( TKOHitCollection *tko );
		void Finalize();
		void Clear();

		void InitializeHistogram();
};

	EventAnalysisMyCDC::EventAnalysisMyCDC()
: EventTemp()
{

}

EventAnalysisMyCDC::~EventAnalysisMyCDC()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysisMyCDC::Initialize( ConfMan *conf )
{
#if 0
	std::cout << " Enter EventAnalysisMyCDC::Initialize " << std::endl;
#endif
	confMan = conf;

	TString cdcfname=confMan->GetCDSTrackFileName();
	cdcFile = new TFile(cdcfname);
	if(!cdcFile->IsOpen()){
		std::cout<<" failed to open " <<cdcfname<< "  !!!"<<std::endl;
		exit(false);
	}
	std::cout<<" CDC Track File Successfully opend !!! " <<cdcfname<<std::endl;
	cdstrackMan   = new CDSTrackingMan();
	header2    = new EventHeader();
	cdcTree=(TTree*)cdcFile->Get("EventTree");
	cdcTree->SetBranchAddress( "CDSTrackingMan", &cdstrackMan );
	cdcTree->SetBranchAddress( "EventHeader" ,&header2);

	rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
	InitializeHistogram();

	header = new EventHeader();
	if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
	blMan = new BeamLineHitMan();
	if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
	bltrackMan = new BeamLineTrackMan();
	if( bltrackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
	cdsMan = new CDSHitMan();
	if( cdsMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
	cdsMan = new CDSHitMan();
	if( cdsMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }

	blAna = new MyAnalysisBL(rtFile, confMan);

	AllGoodTrack=0;
	nTrack=0;
	CDC_Event_Number=0;
}

void EventAnalysisMyCDC::USca( int nsca, unsigned int *sca )
{
#if 0
	std::cout << " Enter EventAnalysisMyCDC::USca " << std::endl;
#endif
	Block_Event_Number++;
	header->SetBlockEventNumber( Block_Event_Number );
}


bool EventAnalysisMyCDC::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyCDC::UAna " << std::endl;
#endif

	Event_Number++;
	{ int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
		if( status==1 ) return true;
		if( status==2 ) return false; }

	if( Event_Number%5000==0 )
		std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack  <<"  "<<CDC_Event_Number<<" / "<<cdcTree->GetEntries() << std::endl;

	header->SetRunNumber(0);
	header->SetEventNumber(Event_Number);

	TH1F *h1;
	TH2F *h2;

	header->Convert( tko, confMan );
	blMan->Convert( tko, confMan );
	cdsMan->Convert( tko, confMan );

	Particle* particle = new Particle();

	if(!blAna->DoAnalysis(confMan, header, blMan, bltrackMan, particle)){
		Clear();
		return true;
	}
	if(!blAna->Fiducial()){
		Clear();
		return true;
	}

	rtFile->cd();

	int MulBHD=0;
	for(int i=0; i<blMan->nBHD(); i++){
		HodoscopeLikeHit* hit = blMan->BHD(i);
		if(hit->CheckRange()) MulBHD++;
	}
	double tmptof = -999;
	int MulT0=0;
	for(int i=0; i<blMan->nT0(); i++){
		HodoscopeLikeHit* hit = blMan->T0(i);
		if(hit->CheckRange()){
			MulT0++;
			tmptof = -(hit->ctmean());
		}
	}
	int MulCDH=0;
	int CDHSeg=-1;
	bool CDHflag = false;
	double CDHazimang = -100.0;
	for(int i=0; i<cdsMan->nCDH(); i++){
		HodoscopeLikeHit* hit = cdsMan->CDH(i);
		if(hit->CheckRange()){
			MulCDH++;
			h1 = (TH1F*)gFile->Get( "CDH_dE" ); h1->Fill(hit->emean());
			if(5.5<hit->emean()&&hit->emean()<9.0){
				tmptof += hit->ctmean();
				CDHflag = true;
				CDHSeg = hit->seg();
				TVector3 pos = hit->pos();
				CDHazimang = atan(pos.Y()/pos.X());
				if(pos.Y()<0&&pos.X()<0){
					CDHazimang -= TMath::Pi();
				}
				if(pos.Y()<0&&pos.X()>0){
					CDHazimang += TMath::Pi();
				}
				h1 = (TH1F*)gFile->Get( "CDH_AzimAngle" ); h1->Fill(CDHazimang/TMath::Pi()*180.0 + 180.0);
			}
		}
	}
	int MulIH=0;
	bool IHflag = false;
	double IHazimang = -100.0;
	for(int i=0; i<cdsMan->nIH(); i++){
		HodoscopeLikeHit* hit = cdsMan->IH(i);
		if(hit->CheckRange()){
			MulIH++;
			h1 = (TH1F*)gFile->Get( "IH_dE" ); h1->Fill(hit->ene(0));
			if(0.3<hit->ene(0)&&hit->ene(0)<0.9){
				IHflag = true;
				TVector3 pos = hit->pos();
				IHazimang = atan(pos.Y()/pos.X());
				if(pos.Y()<0&&pos.X()<0){
					IHazimang -= TMath::Pi();
				}
				if(pos.Y()<0&&pos.X()>0){
					IHazimang += TMath::Pi();
				}
				h1 = (TH1F*)gFile->Get( "IH_AzimAngle" ); h1->Fill(IHazimang/TMath::Pi()*180.0 + 180.0);
			}
		}
	}
	// Selection //
	if(MulT0!=1 || MulCDH!=1 || MulIH!=1 || header->IsTrig(Trig_Cosmic)){
		Clear();
		return true;
	}
	if(!CDHflag || !IHflag){
		Clear();
		return true;
	}
	double CDHIH_angledif = TMath::Abs(CDHazimang-IHazimang);
	if(CDHIH_angledif>TMath::Pi()){
		CDHIH_angledif -= TMath::Pi();
	}
	h1 = (TH1F*)gFile->Get( "CDHIH_AzimAngleDifference" ); h1->Fill(CDHIH_angledif/TMath::Pi()*180.0);
	h2 = (TH2F*)gFile->Get( "CDHvsIH_AzimAngle" ); h2->Fill(CDHazimang/TMath::Pi()*180.0+180.0,IHazimang/TMath::Pi()*180.0+180.0);
	if(CDHIH_angledif/TMath::Pi()*180.0>40){
		Clear();
		return true;
	}
	h1 = (TH1F*)gFile->Get( "CDH_TOF" ); h1->Fill(tmptof);
	if(!(0.0<tmptof&&tmptof<10.0)){
		Clear();
		return true;
	}
	h1 = (TH1F*)gFile->Get( "CDC_TrackEfficiency" ); h1->Fill(0);


	/* CDC tracking event maching */
	if(CDC_Event_Number>=cdcTree->GetEntries()) return false;  
	cdcTree->GetEntry(CDC_Event_Number);
	if(header2->ev()!=Event_Number){
		int tmpnum=CDC_Event_Number;
		while(header2->ev()<Event_Number&&tmpnum<=cdcTree->GetEntries()){
			tmpnum++;    
			cdcTree->GetEntry(tmpnum);
		}
		if(header2->ev()!=Event_Number)  return true;
		else CDC_Event_Number=tmpnum;
	}
	CDC_Event_Number++;
	/* CDC tracking event maching */
	cdstrackMan->Calc( cdsMan, confMan );  
	int nGoodTrack=cdstrackMan->nGoodTrack();
	int nallTrack=  cdstrackMan->nTrack();
	AllGoodTrack+=nGoodTrack;
	nTrack+=nallTrack;
	h1 = (TH1F*)gFile->Get( "CDC_NTrack" ); h1->Fill( nallTrack );
	h1 = (TH1F*)gFile->Get( "CDC_NGTrack" ); h1->Fill( nGoodTrack );

	// Selection //
	if(nTrack<1){
		Clear();
		return true;
	}

	double T0_Timing = -999.0;
	double CDH_Timing[36];
	double CDH_Seg[36];
	for(int i=0; i<36; i++){
		CDH_Timing[i] = -9999.0;
		CDH_Seg[i] = -9999.0;
	}

	// T0
	int nT0=0;
	{
		int HitSeg[2]={-1,-1};
		for( int i=0; i<blMan->nT0(); i++ ){
			HodoscopeLikeHit *hit = blMan->T0(i);
			int seg = hit->seg();
			double emean = hit->emean();
			double ctmean = hit->ctmean();
#if Debug
			std::cout << "= trg =" << std::endl;
			std::cout << "seg = " << seg << std::endl;
#endif

			if( hit->CheckRange() ){
				h1 = (TH1F*)gFile->Get( Form("T0%d_dEMean",seg) ); h1->Fill( emean );
				h1 = (TH1F*)gFile->Get( Form("T0%d_CTimeMean",seg) ); h1->Fill( ctmean );
				h1 = (TH1F*)gFile->Get( "T0_HitPat" ); h1->Fill( seg );
				nT0++;
				HitSeg[i] = seg;
				T0_Timing = ctmean;
			}
		}
		h1 = (TH1F*)gFile->Get( "T0_Mult" ); h1->Fill( nT0 );
	}

	double momentum[2] = {0.0,0.0};
	for(int it=0; it<cdstrackMan->nGoodTrack(); it++){
		CDSTrack* track = cdstrackMan->Track(cdstrackMan->GoodTrackID(it));
		double mom = track->Momentum();
		double Chi =  track->Chi();
		h1 = (TH1F*)gFile->Get( "CDC_ChiSquare" ); h1->Fill(Chi); 
		h1 = (TH1F*)gFile->Get( Form("Momentum") ); h1->Fill(mom); 
		double dE_IH = 0.;
		double dE_CDH = 0.;
		momentum[it] = mom;

		// IH
		int nIH=0;
		h1 = (TH1F*)gFile->Get( "IH_Efficiency" ); h1->Fill( 0 );
		{
			int HitSeg[2]={-1,-1};
			for( int i=0; i<track->nIHHit(); i++ ){
				HodoscopeLikeHit *hit = track->IHHit(cdsMan,i);
				int seg = hit->seg();
				double emean = hit->emean();
				double ctmean = hit->ctmean();
				dE_IH = emean;
#if Debug
				std::cout << "= trg =" << std::endl;
				std::cout << "seg = " << seg << std::endl;
#endif

				if( hit->CheckRange() ){
					h1 = (TH1F*)gFile->Get( Form("IH%d_dEMean",seg) ); h1->Fill( emean );
					h1 = (TH1F*)gFile->Get( Form("IH%d_CTimeMean",seg) ); h1->Fill( ctmean );
					h1 = (TH1F*)gFile->Get( "IH_HitPat" ); h1->Fill( seg );
					nIH++;
					HitSeg[i] = seg;
				}
			}
			h1 = (TH1F*)gFile->Get( "IH_Mult" ); h1->Fill( nIH );
			if(track->SearchHodoHit(cdsMan,confMan,CID_IH,0)){ 
				h1 = (TH1F*)gFile->Get( "IH_Efficiency" ); h1->Fill( 1 );
			}
		}
		// CDH
		int nCDH=0;
		h1 = (TH1F*)gFile->Get( "CDH_Efficiency" ); h1->Fill( 0 );
		{
			int HitSeg[2]={-1,-1};
			for( int i=0; i<track->nCDHHit(); i++ ){
				HodoscopeLikeHit *hit = track->CDHHit(cdsMan,i);
				int seg = hit->seg();
				double emean = hit->emean();
				double ctmean = hit->ctmean();
				dE_CDH = emean;
#if Debug
				std::cout << "= trg =" << std::endl;
				std::cout << "seg = " << seg << std::endl;
#endif

				if( hit->CheckRange() ){
					h1 = (TH1F*)gFile->Get( Form("CDH%d_dEMean",seg) ); h1->Fill( emean );
					h1 = (TH1F*)gFile->Get( Form("CDH%d_CTimeMean",seg) ); h1->Fill( ctmean );
					h1 = (TH1F*)gFile->Get( "CDH_HitPat" ); h1->Fill( seg );
					nCDH++;
					CDH_Timing[i] = ctmean;
					HitSeg[i] = seg;
					CDH_Seg[i] = seg;
				}
			}
			h1 = (TH1F*)gFile->Get( "CDH_Mult" ); h1->Fill( nCDH );
			if(track->SearchHodoHit(cdsMan,confMan,CID_CDH,0)){ 
				h1 = (TH1F*)gFile->Get( "CDH_Efficiency" ); h1->Fill( 1 );
			}
		}
		// CDC
		if(nCDH==1 && CDH_Seg[0]==CDHSeg && Chi<30.0){
			h1 = (TH1F*)gFile->Get( "CDC_TrackEfficiency" ); h1->Fill(1);
			CDHSeg = -1;
		}
		{
			int nhit[NumOfCDCLayers];
			for(int i=0; i<NumOfCDCLayers; i++){
				nhit[i] = 0;
			}
			for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
				h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Mult",layer) ); 
				//h1->Fill( track->nTrackHit(layer) );
				h1->Fill( cdsMan->nCDC(layer) );
				nhit[layer-1]=track->nTrackHit(layer);
				for( int i=0; i<track->nTrackHit(layer); i++ ){
					CDCHit *hit = track->TrackHit(cdsMan,layer,i);
					int wire = hit->wire();
					double dt = hit->dt();
					double dl = hit->dl();
					double resi = hit->resl();
					h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Time",layer) ); h1->Fill( dt );
					h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Length",layer) ); h1->Fill( dl );
					h2 = (TH2F*)gFile->Get( Form("CDC_Layer%d_TimevLength",layer) ); h2->Fill( dt,dl );
					h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Residual",layer) ); h1->Fill( resi );
					h2 = (TH2F*)gFile->Get( Form("CDC_Layer%d_TimevResidual",layer) ); h2->Fill( dt,resi );
					h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_HitPat",layer) ); h1->Fill( wire );
				}
			}
		}
		for(int i=0; i<nCDH; i++){
			if(CDH_Timing[i]>-900&&T0_Timing>-900){
				double TOF = CDH_Timing[i] - T0_Timing;
				h2 = (TH2F*)gFile->Get("T0CDH_TOF"); h2->Fill(CDH_Seg[i],TOF);
				h2 = (TH2F*)gFile->Get("TOFvMomentum"); h2->Fill(TOF,mom);
			}
		}

		h2 = (TH2F*)gFile->Get( Form("dE_IHvMomentum") ); h2->Fill( dE_IH,mom );
		h2 = (TH2F*)gFile->Get( Form("dE_CDHvMomentum") ); h2->Fill( dE_CDH,mom );
	}

	Clear();
	return true;
}

void EventAnalysisMyCDC::Finalize()
{
	std::cout << " Enter EventAnalysisMyCDC::Finalize " << std::endl;

	rtFile->cd();
	gFile->Write();
	gFile->Close();

	delete cdsMan;
	delete blMan;
	delete header;
}

void EventAnalysisMyCDC::Clear()
{
	blAna->Clear();
	blMan->Clear();
	bltrackMan->Clear();
	cdsMan->Clear();
	cdstrackMan->Clear();
	header->Clear();
}


void EventAnalysisMyCDC::InitializeHistogram()
{
	Int_t NumOfT0Segments = 5;
	Int_t NumOfIHSegments = 24;
	Int_t NumOfCDHSegments = 36;
	const int ndc=1;
	int dccid[ndc]={CID_CDC};
	TString dcname[ndc]={"CDC"};
	int NumOfLayers[ndc]={NumOfCDCLayers};

	// TOF
	new TH2F("T0CDH_TOF","TOF",37,0,37,1000,0.,100.);

	new TH1F( "Mass","Mass (GeV/c2)", 200, 0.0, 1.0 );
	new TH1F( "Momentum","Momentum (GeV/c)", 300, -1.5, 1.5 );
	new TH2F("TOFvMomentum","TOF T0-CDH vs. Momentum;Time of flight (ns);Momentum (GeV/c)",750,-25,50,300,-1.5,1.5);
	new TH2F("dE_IHvMomentum","IH dE vs. Momentum;dE (MeV);Momentum (GeV/c)",200,0,20,300,-1.5,1.5);
	new TH2F("dE_CDHvMomentum","CDH dE vs. Momentum;dE (MeV);Momentum (GeV/c)",200,0,20,300,-1.5,1.5);

	// T0
	std::cout << "Define Histograms for T0" << std::endl;
	new TH1F( "T0_HitPat", "Hit Pattern T0", NumOfT0Segments+1, 0, NumOfT0Segments+1 );
	for( int seg=1; seg<=NumOfT0Segments; seg++ ){
		new TH1F( Form("T0%d_dEMean",seg),   Form("Mean dE T0%d;dE (MeV);Counts",seg),    400,    0, 40 );
		new TH1F( Form("T0%d_CTimeMean",seg),   Form("Mean CTime T0%d;Time (ns);Counts",seg),    4000,    -100, 100 );
	}
	new TH1F( "T0_Mult", "Multiplicity T0;Multiplicity;Counts", NumOfT0Segments+1, 0, NumOfT0Segments+1 );

	// IH
	std::cout << "Define Histograms for IH" << std::endl;
	new TH1F( "IH_HitPat", "Hit Pattern IH", NumOfIHSegments+1, 0, NumOfIHSegments+1 );
	for( int seg=1; seg<=NumOfIHSegments; seg++ ){
		new TH1F( Form("IH%d_dEMean",seg),   Form("Mean dE IH%d;dE (MeV);Counts",seg),    400,    0, 40 );
		new TH1F( Form("IH%d_CTimeMean",seg),   Form("Mean CTime IH%d;Time (ns);Counts",seg),    4000,    -100, 100 );
	}
	new TH1F( "IH_AzimAngle", "Azimuthal angle IH;Azimuthal angle;Counts", 15, 0, 360 );
	new TH1F( "IH_Mult", "Multiplicity IH;Multiplicity;Counts", NumOfIHSegments+1, 0, NumOfIHSegments+1 );
	new TH1F( "IH_Efficiency", "Efficiency IH;Efficiency;Counts", 2, 0, 2 );
	new TH1F( "IH_dE", "Energy deposit IH;dE (MeV);Counts", 400, 0.0, 40.0 );

	// CDH
	std::cout << "Define Histograms for CDH" << std::endl;
	new TH1F( "CDH_HitPat", "Hit Pattern CDH", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
	for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
		new TH1F( Form("CDH%d_dEMean",seg),   Form("Mean dE CDH%d;dE (MeV);Counts",seg),    400,    0, 40 );
		new TH1F( Form("CDH%d_CTimeMean",seg),   Form("Mean CTime CDH%d;Time (ns);Counts",seg),    4000,    -100, 100 );
	}
	new TH1F( "CDH_AzimAngle", "Azimuthal angle CDH;Azimuthal angle;Counts", 36, 0, 360 );
	new TH1F( "CDH_Mult", "Multiplicity CDH;Multiplicity;Counts", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
	new TH1F( "CDH_Efficiency", "Efficiency CDH;Efficiency;Counts", 2, 0, 2 );
	new TH1F( "CDH_TOF", "TOF CDH;TOF (ns);Counts", 750, -25.0, 50.0 );
	new TH1F( "CDH_dE", "Energy deposit CDH;dE (MeV);Counts", 400, 0.0, 40.0 );

	// CDC
	for(int idc=0;idc<ndc;idc++){
		std::cout << "Define Histgram for " << dcname[idc] << std::endl;
		new TH1F( Form("%s_ChiSquare",dcname[idc].Data()), Form("ChiSquare %s;ChiSquare;Counts",dcname[idc].Data()), 400, 0, 100 );
		new TH1F( Form("%s_NTrack",dcname[idc].Data()), Form("Number of tracks %s;# of Tracks;Counts",dcname[idc].Data()), 20, 0, 20 );
		new TH1F( Form("%s_NGTrack",dcname[idc].Data()), Form("Number of good tracks %s;# of Good tracks;Counts",dcname[idc].Data()), 20, 0, 20 );
		for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
			int nwire = NumOfCDCWiresInLayer[layer-1];
			new TH1F( Form("%s_Layer%d_Mult",dcname[idc].Data(),layer), Form("Multiplicity %s Layer%d;Multiplicity;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_HitPat",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d;Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_HitPat_mult1",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 1);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_HitPat_mult2",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 2);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_HitPat_mult3",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 3);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_Time",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift time (ns);Counts",dcname[idc].Data(),layer), 4000, -500, 1500 );
			new TH1F( Form("%s_Layer%d_Length",dcname[idc].Data(),layer), Form("Drift length %s Layer%d;Drift length (mm);Counts",dcname[idc].Data(),layer), 1100, -0.1, 1.0 );
			new TH1F( Form("%s_Layer%d_Residual",dcname[idc].Data(),layer), Form("Residual %s Layer%d;Residial (mm);Counts",dcname[idc].Data(),layer), 1000, -0.5, 0.5 );
			new TH2F( Form("%s_Layer%d_TimevLength",dcname[idc].Data(),layer), Form("Drift time vs. Drift length %s Layer%d;Drift time (ns);Drift length (mm)",dcname[idc].Data(),layer), 400, -500, 1500, 110, -0.1, 1.0 );
			new TH2F( Form("%s_Layer%d_TimevResidual",dcname[idc].Data(),layer), Form("Drift time vs. Drift length %s Layer%d;Drift time (ns);Residual (mm)",dcname[idc].Data(),layer), 800, -100, 300, 800, -0.4, 0.4 );
		}
	}
	new TH2F( "CDHvsIH_AzimAngle", "Azimuthal angle CDH vs. IH;Azimuthal angle CDH;Azimuthal angle IH", 36, 0, 360, 15, 0, 360 );
	new TH1F( "CDHIH_AzimAngleDifference", "Azimuthal angle difference between CDH and IH;Azimuthal angle difference;Counts", 18, 0, 180 );
	new TH1F( "CDC_Efficiency", "Efficiency CDC;CDC Layer;Counts", 16, 0, 16 );
	new TH1F( "CDC_TrackEfficiency", "Tracking efficiency CDC;CDC Layer;Counts", 16, 0, 16 );
}



EventTemp *EventAlloc::EventAllocator()
{
	EventAnalysisMyCDC *event = new EventAnalysisMyCDC();
	return (EventTemp*)event;
}
