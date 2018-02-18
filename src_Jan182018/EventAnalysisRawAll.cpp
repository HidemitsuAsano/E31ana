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
class EventAnalysisRawAll: public EventTemp
{
public:
  EventAnalysisRawAll();
  ~EventAnalysisRawAll();
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

EventAnalysisRawAll::EventAnalysisRawAll()
  : EventTemp()
{

}

EventAnalysisRawAll::~EventAnalysisRawAll()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysisRawAll::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::Initialize " << std::endl;
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

void EventAnalysisRawAll::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::USca " << std::endl;
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

#if 0
  for( int i=0; i<scaMan->nsca(); i++ ){
    int val = scaMan->sca(i)->val();
    TString name = scaMan->sca(i)->name();
    std::cout << std::setw(14) << name << std::setw(10) << val;
    if( (i+1)%8==0 ) std::cout << std::endl;
  }
  std::cout << std::endl;
#endif

  scaMan->Clear();
}


bool EventAnalysisRawAll::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::UAna " << std::endl;
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

  //TKO
#if TKOHIS  
  for(int ntko=0;ntko<tko->entries();ntko++)
    {
      TKOHit *tkohit=tko->hit(ntko);
      h1 = (TH1F*)gFile->Get(Form("TKOc%ds%d",tkohit->cr(),tkohit->sl() )); h1->Fill(tkohit->ch());
      h1 = (TH1F*)gFile->Get(Form("TKOc%ds%da%d",tkohit->cr(),tkohit->sl(),tkohit->ch() ) ); h1->Fill( tkohit->data() );
      
    }
#endif 
//#########

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  rtFile->cd();

  if(header->IsTrig(Trig_Cosmic)){
    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
  }
  // Trigger Pattern
  bool patternflag[4];
  patternflag[0]=false;  patternflag[1]=false;  patternflag[2]=false;
  patternflag[3]=false;
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    //std::cout << " i:" << i << " val:" << val << std::endl;
    if( 0<val ){
      h1 = (TH1F*)gFile->Get("Pattern"); h1->Fill(i);
      h1 = (TH1F*)gFile->Get(Form("TPattern%d",i)); h1->Fill(val);
    }
    if( (i==2 || i==4) && val > 0 ) // k /pi
      patternflag[0]=true;
    if( (6<=i && i<=11 ) && val > 0 ) 
      patternflag[1]=true;
    if( (i==12 || i==13 || i==14 || i==16 ) && val > 0 ) 
      patternflag[2]=true;
  }
  h1 = (TH1F*)gFile->Get("Pattern2"); 
  h1->Fill(0);
  if(patternflag[0] ) h1->Fill(1);
  if(patternflag[1] ) h1->Fill(2);
  if(patternflag[2] ) h1->Fill(3);
  if(patternflag[2] && patternflag[0] ) h1->Fill(4);
  // CDS

  // CDC
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    h1 = (TH1F*)gFile->Get( Form("MulCDC%d",layer) ); h1->Fill( cdsMan->nCDC(layer) );
    for( int i=0; i<cdsMan->nCDC(layer); i++ ){
      CDCHit *hit = cdsMan->CDC(layer,i);
      int wire = hit->wire();
      int tdc = hit->tdc();
      double dt = hit->dt();
      h1 = (TH1F*)gFile->Get( Form("TCDC%d",layer) ); h1->Fill( tdc );
      h1 = (TH1F*)gFile->Get( Form("dTCDC%d",layer) ); h1->Fill( dt );
      h1 = (TH1F*)gFile->Get( Form("TCDC%d_%d",layer,wire) ); h1->Fill( tdc );
      h1 = (TH1F*)gFile->Get( Form("HitPatCDC%d",layer) ); h1->Fill( wire );
    }
  }

  // BeamLine Chamber
  const int ndc=6;
  int dccid[ndc]={CID_BLC1a,CID_BLC1b,
		  CID_BLC2a,CID_BLC2b,
		  CID_BPC,
		  CID_FDC1
  };
  TString dcname[ndc]={"BLC1a","BLC1b",
		     "BLC2a","BLC2b",
		     "BPC",
		     "FDC1"
  };
  int NumOfLayers[ndc]={NumOfBLCLayers,NumOfBLCLayers,
			NumOfBLCLayers,NumOfBLCLayers,
			NumOfBPCLayers,
			NumOfFDC1Layers
  };

  for(int idc=0;idc<ndc;idc++){
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      h1 = (TH1F*)gFile->Get( Form("Mul%s_%d",dcname[idc].Data(),layer) ); 
      h1->Fill( blMan->nBLDC(dccid[idc],layer) );
      for( int i=0; i<blMan->nBLDC(dccid[idc],layer); i++ ){
	ChamberLikeHit *hit = blMan->BLDC(dccid[idc],layer,i);
	int wire = hit->wire();
	int tdc = hit->tdc();
	h1 = (TH1F*)gFile->Get( Form("T%s_%d",dcname[idc].Data(),layer) ); h1->Fill( tdc );
	h1 = (TH1F*)gFile->Get( Form("T%s_%d_%d",dcname[idc].Data(),layer,wire) ); h1->Fill( tdc );
	h1 = (TH1F*)gFile->Get( Form("HitPat%s_%d",dcname[idc].Data(),layer) ); h1->Fill( wire );
	if(layer%2==1)
	  {
	    for( int i2=0; i2<blMan->nBLDC(dccid[idc],layer+1); i2++ ){
	      ChamberLikeHit *hit2 = blMan->BLDC(dccid[idc],layer+1,i2);
	      int wire2 = hit2->wire();
	      //	      int tdc2 = hit2->tdc();
	      h2 = (TH2F*)gFile->Get( Form("WirePat%s_%d",dcname[idc].Data(),layer) ); h2->Fill( wire,wire2 );
	    }
	  }
      }
    }
  }

  // Hodoscope

  const int nhodo=13;
  int hodoid[nhodo]={CID_BHD,CID_T0,
		     CID_DEF,CID_BPD,
		     CID_BVC,CID_CVC,
		     CID_PC, CID_NC,
		     CID_BD, CID_LB,
		     CID_WVC,CID_HVC1,
		     CID_HVC2,
  };
  TString hodoname[nhodo]={"BHD","T0",
			   "DEF","BPD",
			   "BVC", "CVC",
			   "PC", "NC",
			   "Beamdump","Longbar",
			   "WVC","HVC1","HVC2"
  };

  bool CVC=false;
  for(int ihodo=0;ihodo<nhodo;ihodo++){
    int nHodo=0;
    for( int i=0; i<blMan->nHodo(hodoid[ihodo]); i++ ){
      HodoscopeLikeHit *hit = blMan->Hodoi(hodoid[ihodo],i);
      int seg = hit->seg();
      int au = hit->adcu(), ad = hit->adcd();
      int tu = hit->tdcu(), td = hit->tdcd();
      h1 = (TH1F*)gFile->Get( Form("A%sU%d",hodoname[ihodo].Data(),seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("A%sD%d",hodoname[ihodo].Data(),seg) ); h1->Fill( ad );
      h1 = (TH1F*)gFile->Get( Form("T%sU%d",hodoname[ihodo].Data(),seg) ); h1->Fill( tu );
      h1 = (TH1F*)gFile->Get( Form("T%sD%d",hodoname[ihodo].Data(),seg) ); h1->Fill( td );
      if(tu>0&&tu<4095){
	h1 = (TH1F*)gFile->Get( Form("AwT2%sU%d",hodoname[ihodo].Data(),seg) ); h1->Fill( au );
      }
      if(td>0&&td<4095){
	h1 = (TH1F*)gFile->Get( Form("AwT2%sD%d",hodoname[ihodo].Data(),seg) ); h1->Fill( ad );
      }
      if(
#if RUN43
	 hodoid[ihodo]==CID_BVC
#elseif RUN46
	 hodoid[ihodo]==CID_WVC||hodoid[ihodo]==CID_HVC1
#else
	 hodoid[ihodo]==CID_WVC
#endif
	 ){
	if( tu>0 ){
	  h1 = (TH1F*)gFile->Get( Form("AwT%sU%d",hodoname[ihodo].Data(),seg) ); h1->Fill( au );
	  h2 = (TH2F*)gFile->Get( Form("AT%sU%d",hodoname[ihodo].Data(),seg) );  h2->Fill( tu, au );
	  nHodo++;
	}else{
	  h1 = (TH1F*)gFile->Get( Form("AwoT%sU%d",hodoname[ihodo].Data(),seg) ); h1->Fill( au );
	}
      }else{      
	if( hit->CheckRange() ){
	  h1 = (TH1F*)gFile->Get( Form("AwT%sU%d",hodoname[ihodo].Data(),seg) ); h1->Fill( au );
	  h1 = (TH1F*)gFile->Get( Form("AwT%sD%d",hodoname[ihodo].Data(),seg) ); h1->Fill( ad );
	  h2 = (TH2F*)gFile->Get( Form("AT%sU%d",hodoname[ihodo].Data(),seg) );  h2->Fill( tu, au );
	  h2 = (TH2F*)gFile->Get( Form("AT%sD%d",hodoname[ihodo].Data(),seg) );  h2->Fill( td, ad );
	  h1 = (TH1F*)gFile->Get( Form("HitPat%s",hodoname[ihodo].Data()) ); h1->Fill( seg );
	  if(hodoid[ihodo] ==CID_BD){
	    h1 = (TH1F*)gFile->Get( "HitPatCharged" ); h1->Fill( -12+2*seg );
	  }
	  if(hodoid[ihodo] ==CID_CVC){
	    CVC=true;
	    h1 = (TH1F*)gFile->Get( "HitPatCharged" ); h1->Fill( seg );
	  }
	  if(hodoid[ihodo] ==CID_PC){
	    h1 = (TH1F*)gFile->Get( "HitPatCharged" ); h1->Fill( 34+seg );
	  }
	  nHodo++;
	  if(hodoid[ihodo]==CID_NC){
	    int layer=-1;
	    int sseg=-1;
	    layer=(seg-1)/16+1;
	    sseg=(seg-16*layer+16)%16+1;
	    if( CVC ) {h2 = (TH2F*)gFile->Get( "HitPatNC2DwC" ); h2->Fill(  sseg,layer );}
	    if(!CVC ) {h2 = (TH2F*)gFile->Get( "HitPatNC2DwoC" ); h2->Fill(  sseg,layer );}
	    h2 = (TH2F*)gFile->Get( "HitPatNC2D" ); h2->Fill(  sseg,layer);
	  }
	}else{
	  h1 = (TH1F*)gFile->Get( Form("AwoT%sU%d",hodoname[ihodo].Data(),seg) ); h1->Fill( au );
	  h1 = (TH1F*)gFile->Get( Form("AwoT%sD%d",hodoname[ihodo].Data(),seg) ); h1->Fill( ad );
	}
      }
    }
    h1 = (TH1F*)gFile->Get( Form("Mul%s",hodoname[ihodo].Data()) ); h1->Fill( nHodo );
  }
  
  for( int i=0; i<blMan->nHodo(CID_T0pre); i++ ){
    HodoscopeLikeHit *hit = blMan->Hodoi(CID_T0pre,i);
    int seg = hit->seg();
    int tu = hit->tdcu(), td = hit->tdcd();
    h1 = (TH1F*)gFile->Get( Form("TT0preU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("TT0preD%d",seg) ); h1->Fill( td );
  }
  for( int i=0; i<blMan->nHodo(CID_T0post); i++ ){
    HodoscopeLikeHit *hit = blMan->Hodoi(CID_T0post,i);
    int seg = hit->seg();
    int tu = hit->tdcu(), td = hit->tdcd();
    h1 = (TH1F*)gFile->Get( Form("TT0postU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("TT0postD%d",seg) ); h1->Fill( td );
    h2 = (TH2F*)gFile->Get( Form("TT0main_postU%d",seg) ); h2->Fill( blMan->Hodo(CID_T0,seg)->tdcu(),tu );
    h2 = (TH2F*)gFile->Get( Form("TT0main_postD%d",seg) ); h2->Fill( blMan->Hodo(CID_T0,seg)->tdcd(),td );
  }
  for( int i=0; i<blMan->nHodo(CID_BHDpost); i++ ){
    HodoscopeLikeHit *hit = blMan->Hodoi(CID_BHDpost,i);
    int seg = hit->seg();
    int tu = hit->tdcu();
    h1 = (TH1F*)gFile->Get( Form("TBHDpost%d",seg) ); h1->Fill( tu );
    h2 = (TH2F*)gFile->Get( Form("TBHDmain_post%d",seg) ); h2->Fill( blMan->Hodo(CID_BHD,seg)->tdcmean(),tu );
  }
    
  // CDS Hodoscope
  const int ncdshodo=2;
  int cdshodoid[ncdshodo]={CID_CDH, CID_IH };
  TString cdshodoname[ncdshodo]={"CDH","IH"};
  
  for(int ihodo=0;ihodo<ncdshodo;ihodo++){
    int nHodo=0;
    for( int i=0; i<cdsMan->nHodo(cdshodoid[ihodo]); i++ ){
      HodoscopeLikeHit *hit = cdsMan->Hodoi(cdshodoid[ihodo],i);
      int seg = hit->seg();
      int au = hit->adcu(), ad = hit->adcd();
      int tu = hit->tdcu(), td = hit->tdcd();
      h1 = (TH1F*)gFile->Get( Form("A%sU%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("A%sD%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( ad );
      h1 = (TH1F*)gFile->Get( Form("T%sU%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( tu );
      h1 = (TH1F*)gFile->Get( Form("T%sD%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( td );
      if( hit->CheckRange() ){
	h1 = (TH1F*)gFile->Get( Form("AwT%sU%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( au );
	h1 = (TH1F*)gFile->Get( Form("AwT%sD%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( ad );
	h2 = (TH2F*)gFile->Get( Form("AT%sU%d",cdshodoname[ihodo].Data(),seg) );  h2->Fill( tu, au );
	h2 = (TH2F*)gFile->Get( Form("AT%sD%d",cdshodoname[ihodo].Data(),seg) );  h2->Fill( td, ad );
	h1 = (TH1F*)gFile->Get( Form("HitPat%s",cdshodoname[ihodo].Data()) ); h1->Fill( seg );
	nHodo++;
      }else{
	h1 = (TH1F*)gFile->Get( Form("AwoT%sU%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( au );
	h1 = (TH1F*)gFile->Get( Form("AwoT%sD%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( ad );
      }
      if( cdshodoid[ihodo]==CID_IH && tu>0 && tu<4096 ){
	h1 = (TH1F*)gFile->Get( Form("AwT%sU%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( au );
	h2 = (TH2F*)gFile->Get( Form("AT%sU%d",cdshodoname[ihodo].Data(),seg) );  h2->Fill( tu, au );
	h1 = (TH1F*)gFile->Get( Form("HitPat%s",cdshodoname[ihodo].Data()) ); h1->Fill( seg );
	nHodo++;
      }else{
	h1 = (TH1F*)gFile->Get( Form("AwoT%sU%d",cdshodoname[ihodo].Data(),seg) ); h1->Fill( au );
      }
    }
    h1 = (TH1F*)gFile->Get( Form("Mul%s",cdshodoname[ihodo].Data()) ); h1->Fill( nHodo );
  }
  
  // Cherenkov
  const int nchere=1;
  int chereid[nchere]={CID_AC};
  TString cherename[nchere]={"AC"};
  int nch[nchere]={4};
  for(int ichere=0;ichere<nchere;ichere++){
    for( int i=0; i<blMan->nChere(chereid[ichere]); i++ ){
      CherenkovLikeHit *hit = blMan->Cherei(chereid[ichere],i);
      //      int seg = hit->seg();
      int adcsum=0;
      for(int ich=1;ich<=nch[ichere];ich++){
	int a = hit->adc(ich);
	int t = hit->tdc(ich);
	h1 = (TH1F*)gFile->Get( Form("A%s%d",cherename[ichere].Data(),ich) ); h1->Fill( a );
	h1 = (TH1F*)gFile->Get( Form("T%s%d",cherename[ichere].Data(),ich) ); h1->Fill( t );
	adcsum+=a;
	if( hit->CheckRange(1) ){
	  h1 = (TH1F*)gFile->Get( Form("AwT%s%d",cherename[ichere].Data(),ich) ); h1->Fill( a );
	  h2 = (TH2F*)gFile->Get( Form("AT%s%d",cherename[ichere].Data(),ich) );  h2->Fill( t, a );
	}
      }
      h1 = (TH1F*)gFile->Get( Form("A%sSUM",cherename[ichere].Data()) ); h1->Fill( adcsum );
      if(hit->CheckRange(1)) h1 = (TH1F*)gFile->Get( Form("AwT%sSUM",cherename[ichere].Data()) ); h1->Fill( adcsum );
    }
    //  h1 = (TH1F*)gFile->Get( Form("Mul%s",cherename[ichere]) ); h1->Fill( nChere[ichere] );
  }

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalysisRawAll::Finalize()
{
  std::cout << " Enter EventAnalysisRawAll::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisRawAll::InitializeHistogram()
{
  rtFile->cd();
  
  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );
  new TH1F( "Pattern2", "Trigger Pattern2", 5, 0, 5 );
  for(int i=0;i<20;i++)
    new TH1F( Form("TPattern%d",i), Form("Trigger Pattern %d TDC",i), 4000, 0, 4000 );
  
#if TKOHIS
  std::cout << "Define Histgram for TKO" << std::endl;
  for( int cr=0; cr<=6; cr++ ){
    for( int sl=1; sl<=22; sl++ ){
      new TH1F( Form("TKOc%ds%d",cr,sl), Form("TKO cr%d sl%d" ,cr,sl), 32, 0, 32+1 );
      for( int ch=0; ch<=32; ch++ )
	{  
	  new TH1F( Form("TKOc%ds%da%d",cr,sl,ch), Form("TKO cr%d sl%d ch%d" ,cr,sl,ch), 4000, 0, 4000 );
	}
    }
  }
#endif

  // CDS
  // CDC
  std::cout << "Define Histgram for CDC" <<  std::endl;
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    int nwire = confMan->GetCDCWireMapManager()->nw(layer);
    new TH1F( Form("MulCDC%d",layer), Form("Multiplicity CDC Layer%d",layer), nwire+1, -0.5, nwire+0.5 );
    new TH1F( Form("HitPatCDC%d",layer), Form("Hit Pattern CDC Layer%d",layer), nwire, 0.5, nwire+0.5 );
    new TH1F( Form("TCDC%d",layer), Form("TDC CDC Layer%d",layer), 2000, 0, 2000 );
    new TH1F( Form("dTCDC%d",layer), Form("dTDC CDC Layer%d",layer), 2200, -50, 500 );
    for( int wire=1; wire<=nwire; wire++ ){
      new TH1F( Form("TCDC%d_%d",layer,wire), Form("TDC CDC Layer%d Wire%d",layer,wire), 2000, 0, 2000 );
    }
  }
  // Hodoscope
  const int nhodo=15;
  TString hodoname[nhodo]={"BHD","T0",
			   "DEF","BPD",
			   "BVC", "CVC",
			   "PC", "NC",
			   "Beamdump","Longbar",
			   "CDH","IH",
			   "WVC","HVC1","HVC2"
  };
  int NumOfSegments[nhodo]={NumOfBHDSegments,
			    NumOfT0Segments,
			    NumOfDEFSegments,
			    NumOfBPDSegments,
			    NumOfBVCSegments,
			    NumOfCVCSegments,
			    NumOfPCSegments,
			    NumOfNCSegments,
			    NumOfBDSegments,
			    NumOfLBSegments,
			    NumOfCDHSegments,
			    NumOfIHSegments,
			    NumOfWVCSegments,
			    NumOfHVC1Segments,
			    NumOfHVC2Segments
  };
  for(int ihodo=0;ihodo<nhodo;ihodo++){
    std::cout << "Define Histgram for "<< hodoname[ihodo].Data() << std::endl;
    new TH1F( Form("Mul%s", hodoname[ihodo].Data()), Form("Multiplicity %s", hodoname[ihodo].Data()), NumOfSegments[ihodo]+1, -0.5, NumOfSegments[ihodo]+0.5 );
    new TH1F( Form("HitPat%s", hodoname[ihodo].Data()), Form("Hit Pattern %s", hodoname[ihodo].Data()), NumOfSegments[ihodo], 0.5, NumOfSegments[ihodo]+0.5 );
    for( int seg=1; seg<=NumOfSegments[ihodo]; seg++ ){
      new TH1F( Form("A%sU%d",hodoname[ihodo].Data(),seg),   Form("ADC %sU%d",hodoname[ihodo].Data(),seg),    4000,    0, 4000 );
      new TH1F( Form("A%sD%d",hodoname[ihodo].Data(),seg),   Form("ADC %sD%d",hodoname[ihodo].Data(),seg),    4000,    0, 4000 );
      new TH1F( Form("T%sU%d",hodoname[ihodo].Data(),seg),   Form("TDC %sU%d",hodoname[ihodo].Data(),seg),    4000,    0, 4000 );
      new TH1F( Form("T%sD%d",hodoname[ihodo].Data(),seg),   Form("TDC %sD%d",hodoname[ihodo].Data(),seg),    4000,    0, 4000 );
      new TH1F( Form("AwT%sU%d",hodoname[ihodo].Data(),seg), Form("ADC wTDC %sU%d",hodoname[ihodo].Data(),seg),  4000,    0, 4000 );
      new TH1F( Form("AwT%sD%d",hodoname[ihodo].Data(),seg), Form("ADC wTDC %sD%d",hodoname[ihodo].Data(),seg),  4000,    0, 4000 );
      new TH1F( Form("AwoT%sU%d",hodoname[ihodo].Data(),seg), Form("ADC wo TDC %sU%d",hodoname[ihodo].Data(),seg),  4000,    0, 4000 );
      new TH1F( Form("AwoT%sD%d",hodoname[ihodo].Data(),seg), Form("ADC wo TDC %sD%d",hodoname[ihodo].Data(),seg),  4000,    0, 4000 );
      new TH1F( Form("AwT2%sU%d",hodoname[ihodo].Data(),seg), Form("ADC wTDC %sU%d",hodoname[ihodo].Data(),seg),  4000,    0, 4000 );
      new TH1F( Form("AwT2%sD%d",hodoname[ihodo].Data(),seg), Form("ADC wTDC %sD%d",hodoname[ihodo].Data(),seg),  4000,    0, 4000 );
      new TH2F( Form("AT%sU%d",hodoname[ihodo].Data(),seg),   Form("ADC TDC corr. %sU%d",hodoname[ihodo].Data(),seg),     200,    0, 4000,  200,    0, 4000 );
      new TH2F( Form("AT%sD%d",hodoname[ihodo].Data(),seg),   Form("ADC TDC corr. %sD%d",hodoname[ihodo].Data(),seg),     200,    0, 4000,  200,    0, 4000 );
    }
  }

  for(int seg=1;seg<=5;seg++){
    new TH1F( Form("TT0preU%d",seg), Form("TT0preU%d",seg), 4000,0,4000);
    new TH1F( Form("TT0preD%d",seg), Form("TT0preD%d",seg), 4000,0,4000);
    new TH1F( Form("TT0postU%d",seg), Form("TT0postU%d",seg), 4000,0,4000);
    new TH1F( Form("TT0postD%d",seg), Form("TT0postD%d",seg), 4000,0,4000);
    new TH2F( Form("TT0main_postU%d",seg), Form("TT0main_postU%d",seg), 200,0,4000,200,0,4000);
    new TH2F( Form("TT0main_postD%d",seg), Form("TT0main_postD%d",seg), 200,0,4000,200,0,4000);
  }    
  for(int seg=1;seg<=20;seg++){
    new TH1F( Form("TBHDpost%d",seg), Form("TBHDpost%d",seg), 4000,0,4000);
    new TH2F( Form("TBHDmain_post%d",seg), Form("TBHDmain_post%d",seg), 200,0,4000,200,0,4000);
  }    
  // BeamLine
  // BLC
  const int ndc=6;
  int dccid[ndc]={CID_BLC1a,CID_BLC1b,
		  CID_BLC2a,CID_BLC2b,
		  CID_BPC,
		  CID_FDC1
  };
  TString dcname[ndc]={"BLC1a","BLC1b",
		     "BLC2a","BLC2b",
		     "BPC",
		     "FDC1"
  };
  int NumOfLayers[ndc]={NumOfBLCLayers,NumOfBLCLayers,
			NumOfBLCLayers,NumOfBLCLayers,
			NumOfBPCLayers,
			NumOfFDC1Layers
  };
  for(int idc=0;idc<ndc;idc++){
    std::cout << "Define Histgram for " << dcname[idc] << std::endl;
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire( dccid[idc], layer );
      new TH1F( Form("Mul%s_%d",dcname[idc].Data(),layer), Form("Multiplicity %s Layer%d",dcname[idc].Data(),layer), nwire+1, -0.5, nwire+0.5 );
      new TH1F( Form("HitPat%s_%d",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d",dcname[idc].Data(),layer), nwire, 0.5, nwire+0.5 );
      if(layer==1 ||layer==3 || layer==5 ||layer==7) 
	new TH2F( Form("WirePat%s_%d",dcname[idc].Data(),layer), Form("Wire Pattern %s Layer%d",dcname[idc].Data(),layer), nwire, 0.5, nwire+0.5, nwire, 0.5, nwire+0.5 );
      new TH1F( Form("T%s_%d",dcname[idc].Data(),layer), Form("TDC %s Layer%d",dcname[idc].Data(),layer), 2000, 0, 2000 );
      for( int wire=1; wire<=nwire; wire++ ){
	new TH1F( Form("T%s_%d_%d",dcname[idc].Data(),layer,wire), Form("TDC %s Layer%d Wire%d",dcname[idc].Data(),layer,wire), 2000, 0, 2000 );
      }
    }
  }
  // Cherenkov
  const int nchere=1;
  //  int chereid[nchere]={CID_AC};
  TString cherename[nchere]={"AC"};
  int nch[nchere]={4};
  for(int ichere=0;ichere<nchere;ichere++){
    std::cout << "Define Histgram for "<< cherename[ichere].Data() << std::endl;
    new TH1F( Form("Mul%s", cherename[ichere].Data()), Form("Multiplicity %s", cherename[ichere].Data()), NumOfSegments[ichere]+1, 0, NumOfSegments[ichere]+1 );
    new TH1F( Form("HitPat%s", cherename[ichere].Data()), Form("Hit Pattern %s", cherename[ichere].Data()), NumOfSegments[ichere]+1, 0, NumOfSegments[ichere]+1 );
    for( int seg=1; seg<=nch[ichere]; seg++ ){
      new TH1F( Form("A%s%d",cherename[ichere].Data(),seg),   Form("ADC %s%d",cherename[ichere].Data(),seg),    4000,    0, 4000 );
      new TH1F( Form("T%s%d",cherename[ichere].Data(),seg),   Form("TDC %s%d",cherename[ichere].Data(),seg),    4000,    0, 4000 );
      new TH1F( Form("AwT%s%d",cherename[ichere].Data(),seg), Form("ADC wTDC %s%d",cherename[ichere].Data(),seg),  4000,    0, 4000 );
      new TH1F( Form("AwoT%s%d",cherename[ichere].Data(),seg), Form("ADC wo TDC %s%d",cherename[ichere].Data(),seg),  4000,    0, 4000 );
      new TH2F( Form("AT%s%d",cherename[ichere].Data(),seg),   Form("ADC TDC corr. %s%d",cherename[ichere].Data(),seg),     200,    0, 4000,  200,    0, 4000 );
    }
    new TH1F( Form("A%sSUM",cherename[ichere].Data()),   Form("ADC %s SUM",cherename[ichere].Data()),    4000,    0, 4000 );
    new TH1F( Form("AwT%sSUM",cherename[ichere].Data()),   Form("ADC wTDC %s SUM",cherename[ichere].Data()),    4000,    0, 4000 );
  }
  // Foward ToF
  std::cout << "Define Histgram for Foward TOF" << std::endl;
  new TH1F( "HitPatCharged", "Hit Pattern Charged TOF", 74+1, -10, 64+1 );
  new TH2F( "HitPatNC2D", "Hit Pattern NC 2D", 16, 0.5, 16.5,7,0.5,7.5 );
  new TH2F( "HitPatNC2DwC", "Hit Pattern NC 2D with charged", 16, 0.5, 16.5,7,0.5,7.5 );
  new TH2F( "HitPatNC2DwoC", "Hit Pattern NC 2D without charged", 16, 0.5, 16.5,7,0.5,7.5 );
}



EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisRawAll *event = new EventAnalysisRawAll();
  return (EventTemp*)event;
}
