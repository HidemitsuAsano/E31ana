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
#include "GlobalVariables.h"

#define RUN43 0
#define TKOHIS 0
class EventAnalysisMyBLC1Test: public EventTemp
{
  public:
    EventAnalysisMyBLC1Test();
    ~EventAnalysisMyBLC1Test();
  private:
    TFile *rtFile;
    TTree *evTree;
    TTree *scaTree;
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

  EventAnalysisMyBLC1Test::EventAnalysisMyBLC1Test()
: EventTemp()
{

}

EventAnalysisMyBLC1Test::~EventAnalysisMyBLC1Test()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysisMyBLC1Test::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisMyBLC1Test::Initialize " << std::endl;
#endif
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }

}

void EventAnalysisMyBLC1Test::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyBLC1Test::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}


bool EventAnalysisMyBLC1Test::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisMyBLC1Test::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

    if( Event_Number%5000==0 )
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
    header->SetRunNumber(0);
    header->SetEventNumber(Event_Number);

    TH1F *h1;
    TH2F *h2;

    header->Convert( tko, confMan );
    blMan->Convert( tko, confMan );
    rtFile->cd();

		int MulBLC1a[8]={0,0,0,0,0,0,0,0};
		for(int layer=1; layer<=8; layer++){
			for(int i=0; i<blMan->nBLC1a(layser); i++){
				HodoscopeLikeHit* hit = blMan->BLC1a(i);
				MulBLC1a[i]++;
			}
		}
		int MulBLC1b[8]={0,0,0,0,0,0,0,0};
		for(int layer=1; layer<=8; layer++){
			for(int i=0; i<blMan->nBLC1b(layser); i++){
				HodoscopeLikeHit* hit = blMan->BLC1b(i);
				MulBLC1b[i]++;
			}
		}

		// Selection //
		if(BLC1a[0]==0 || BLC1a[8]==0 || BLC1b[0]==0 || BLC1b[8]){
			

		// BLC1a
		{    
			int nhit[NumOfBLC2Layers];
			for(int i=0; i<NumOfBLC2Layers; i++){
				nhit[i] = 0;
			}
			for( int layer=1; layer<=NumOfBLC2Layers; layer++ ){
				h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_Mult",layer) ); 
				h1->Fill( blMan->nBLC2b(layer) );
				nhit[layer-1]=blMan->nBLC2b(layer);
				for( int i=0; i<blMan->nBLC2b(layer); i++ ){
					ChamberLikeHit *hit = blMan->BLC2b(layer,i);
					int wire = hit->wire();
					int tdc = hit->tdc();
					double dt = hit->dt(), cdt = hit->cdt();
					double dl = hit->dl(), cdl = hit->cdl();
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_TDC",layer) ); h1->Fill( tdc );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_Wire%d_TDC",layer,wire) ); h1->Fill( tdc );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_Time",layer) ); h1->Fill( dt );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_Wire%d_Time",layer,wire) ); h1->Fill( dt );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_CTime",layer) ); h1->Fill( cdt );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_Wire%d_CTime",layer,wire) ); h1->Fill( cdt );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_Length",layer) ); h1->Fill( dl );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_Wire%d_Length",layer,wire) ); h1->Fill( dl );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_CLength",layer) ); h1->Fill( cdl );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_Wire%d_CLength",layer,wire) ); h1->Fill( cdl );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_HitPat",layer) ); h1->Fill( wire );
					for(int mult=1; mult<=3; mult++){
						if(nhit[layer-1]==mult){
							h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%d_HitPat_mult%d",layer,mult) ); h1->Fill( wire );
						}
					}
				}
			}
			for(int mult=1; mult<=9; mult++){
				h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency_mult%d",mult) ); h1->Fill( 0 );
				h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency_multleq%d",mult) ); h1->Fill( 0 );
				if( nhit[0]==1 && nhit[NumOfBLC2Layers-1]==1){
					h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency_mult%d",mult) ); h1->Fill( 1 );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency_multleq%d",mult) ); h1->Fill( 1 );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency_mult%d",mult) ); h1->Fill( NumOfBLC2Layers );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency_multleq%d",mult) ); h1->Fill( NumOfBLC2Layers );
					for( int layer=2; layer<NumOfBLC2Layers; layer++ ){
						if(nhit[layer-1]==mult){
							h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency_mult%d",mult) ); h1->Fill( layer );
						}
						if(nhit[layer-1]<=mult && nhit[layer-1]!=0){
							h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency_multleq%d",mult) ); h1->Fill( layer );
						}
					}
				}
				bool hitflag1 = false;
				for(int layer=1; layer<=NumOfBLC2Layers; layer++){
					for(int i=1; i<=NumOfBLC2Layers; i++){
						if(i!=layer){
							if(nhit[i-1]==1){
								hitflag1 = true;
							}
							else{
								hitflag1 = false;
							}
							if(!hitflag1) break;
						}
					}
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 0 );
					h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 0 );
					if(hitflag1){
						h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 1 );
						if(nhit[layer-1]==mult){
							h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 2 );
						}
						h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 1 );
						if(nhit[layer-1]<=mult){
							h1 = (TH1F*)gFile->Get( Form("BLC2b_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 2 );
						}
					}
				}
			}
		}

		header->Clear();
		blMan->Clear();
		return true;
}

void EventAnalysisMyBLC1Test::Finalize()
{
	std::cout << " Enter EventAnalysisMyBLC1Test::Finalize " << std::endl;

	rtFile->cd();
	gFile->Write();
	gFile->Close();

	delete blMan;
	delete header;
}

void EventAnalysisMyBLC1Test::InitializeHistogram()
{
	const int ndc=4;
	int dccid[ndc]={CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b};
	TString dcname[ndc]={"BLC1a","BLC1b","BLC2a","BLC2b"};
	int NumOfLayers[ndc]={NumOfBLC1Layers,NumOfBLC1Layers,NumOfBLC2Layers,NumOfBLC2Layers};

	for(int idc=0;idc<ndc;idc++){
		std::cout << "Define Histgram for " << dcname[idc] << std::endl;
		for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
			int nwire = confMan->GetBLDCWireMapManager()->GetNWire( dccid[idc], layer );
			new TH1F( Form("%s_Layer%d_Mult",dcname[idc].Data(),layer), Form("Multiplicity %s Layer%d;Multiplicity;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_HitPat",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d;Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_HitPat_mult1",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 1);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_HitPat_mult2",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 2);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_HitPat_mult3",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 3);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
			new TH1F( Form("%s_Layer%d_TDC",dcname[idc].Data(),layer), Form("TDC %s Layer%d;TDC ch.;Counts",dcname[idc].Data(),layer), 4000, 0, 4000 );
			for( int wire=1; wire<=nwire; wire++ ){
				new TH1F( Form("%s_Layer%d_Wire%d_TDC",dcname[idc].Data(),layer,wire), Form("TDC %s Layer%d Wire%d;TDC ch.;Counts",dcname[idc].Data(),layer,wire), 4000, 0, 4000 );
			}
			new TH1F( Form("%s_Layer%d_Time",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift time (ns);Counts",dcname[idc].Data(),layer), 3000, -500, 1000 );
			new TH1F( Form("%s_Layer%d_CTime",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift time (ns);Counts",dcname[idc].Data(),layer), 3000, -500, 1000 );
			new TH1F( Form("%s_Layer%d_Length",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift length (mm);Counts",dcname[idc].Data(),layer), 1000, -1.0, 9.0 );
			new TH1F( Form("%s_Layer%d_CLength",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift length (mm);Counts",dcname[idc].Data(),layer), 1000, -1.0, 9.0 );
			for( int wire=1; wire<=nwire; wire++ ){
				new TH1F( Form("%s_Layer%d_Wire%d_Time",dcname[idc].Data(),layer,wire), Form("Drift time %s Layer%d Wire%d;Drift time (ns);Counts",dcname[idc].Data(),layer,wire), 3000, -500, 1000 );
				new TH1F( Form("%s_Layer%d_Wire%d_CTime",dcname[idc].Data(),layer,wire), Form("Drift time %s Layer%d Wire%d;Drift time (ns);Counts",dcname[idc].Data(),layer,wire), 3000, -500, 1000 );
				new TH1F( Form("%s_Layer%d_Wire%d_Length",dcname[idc].Data(),layer,wire), Form("Drift time %s Layer%d Wire%d;Drift Length (mm);Counts",dcname[idc].Data(),layer,wire), 1000, -1.0,9.0 );
				new TH1F( Form("%s_Layer%d_Wire%d_CLength",dcname[idc].Data(),layer,wire), Form("Drift time %s Layer%d Wire%d;Drift Length (mm);Counts",dcname[idc].Data(),layer,wire), 1000, -1.0,9.0 );
			}
		}
		for(int mult=1; mult<=9; mult++){
			new TH1F(Form("%s_LayerEfficiency_mult%d",dcname[idc].Data(),mult),Form("%s Layer Efficiency (mult=%d);Layer;Counts",dcname[idc].Data(),mult),NumOfLayers[idc]+1,0,NumOfLayers[idc]+1);
			new TH1F(Form("%s_LayerEfficiency_multleq%d",dcname[idc].Data(),mult),Form("%s Layer Efficiency (mult#leq%d);Layer;Counts",dcname[idc].Data(),mult),NumOfLayers[idc]+1,0,NumOfLayers[idc]+1);
			for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
				new TH1F(Form("%s_Layer%dEfficiency_mult%d",dcname[idc].Data(),layer,mult),Form("%s Layer %d Efficiency (mult=%d);;Counts",dcname[idc].Data(),layer,mult),3,0,3);
				new TH1F(Form("%s_Layer%dEfficiency_multleq%d",dcname[idc].Data(),layer,mult),Form("%s Layer %d Efficiency (mult#leq%d);;Counts",dcname[idc].Data(),layer,mult),3,0,3);
			}
		}
	}
}



EventTemp *EventAlloc::EventAllocator()
{
	EventAnalysisMyBLC1Test *event = new EventAnalysisMyBLC1Test();
	return (EventTemp*)event;
}
