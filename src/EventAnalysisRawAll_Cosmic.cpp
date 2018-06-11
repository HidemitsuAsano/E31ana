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
#include "Tools.h"

#include <TParameter.h>

#define RUN43 0
#define TKOHIS 0
#define ENERGY 0

class EventAnalysisRawAll: public EventTemp
{
public:
  EventAnalysisRawAll();
  ~EventAnalysisRawAll();
private:
  TFile *rtFile;
  TTree *evTree;
  //TTree *scaTree;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;
  double scainit[40];
  double scaend[40];
  bool INIT;
  bool SCAOVERFLOW[40];
  Time_t Time;
  long t_init;
  time_t t0;
  int spillinit;
  int spillfini;
  TString tmpname;
  //flag if there is tdc calibration (i.e) GainMap.param
  bool is_TDCcalib;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  bool EndOfAnEvent(bool flag=true);
  void UTime( int time ){ Time = time; };
};

EventAnalysisRawAll::EventAnalysisRawAll()
  : EventTemp()
{
  is_TDCcalib = true;
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
  // if(gSystem->Exec(Form("test -d tmp")))
  //   gSystem->Exec(Form("mkdir tmp"));
  // int i=0;
  // tmpname=Form("tmp/tmpevanaraw%d_%d.root",i,gRandom->Integer(100));
  // while( !gSystem->Exec( Form("test -f %s",tmpname.Data()) ))
  //   tmpname=Form("tmp/tmpevanaraw%d_%d.root",++i,gRandom->Integer(100));
  // std::cout<<tmpname<<std::endl;
  // rtFile = new TFile( tmpname , "recreate" );
  rtFile = new TFile( conf->GetOutFileName().c_str() , "recreate" );
  INIT=true;
  spillinit=spillfini=-1;
  for(int i=0;i<40;i++){
    SCAOVERFLOW[i]=false;
    scaend[i]=0;
  }
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!! object of EventHeader is NULL " << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!! object of cdsMan is NULL" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!! object of BeamLineHitMan is NULL" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!! object of ScalerMan is NULL" << std::endl; return; }

}

void EventAnalysisRawAll::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::USca " << std::endl;
#endif
  rtFile->cd();
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  for( int i=0; i<8; i++ )
    std::cout<<sca[i]<<"\t";
  std::cout<<std::endl;
#endif
  //  std::cout<<"nsca: "<<nsca<<std::endl;
  if(INIT&&sca[0]<5){
    t_init=Time;    
    INIT=false;
    for( int i=0; i<scaMan->nsca(); i++ )   scainit[i] = scaMan->sca(i)->val();
  }
  if(INIT){
    scaMan->Clear();
    return;
  }
  for( int i=0; i<scaMan->nsca(); i++ ){
    //   int val = scaMan->sca(i)->val()-scaend[i];
    //   TString name = scaMan->sca(i)->name();
    //   Tools::Fill2D(name,Block_Event_Number,val);
    if(scaend[i]>9.9e+07&&scaend[i]>scaMan->sca(i)->val()) SCAOVERFLOW[i]=true;
    scaend[i] = scaMan->sca(i)->val();
    if(SCAOVERFLOW[i]) scaend[i]+=1.0e+08;
    Tools::H1(scaMan->sca(i)->name(),Block_Event_Number,3000,-0.5,2999.5,scaend[i]);
  }
  Tools::H1("Time",Block_Event_Number,3000,-0.5,2999.5,Time-t_init);
  scaMan->Clear();
}


bool EventAnalysisRawAll::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::UAna " << std::endl;
#endif
  Event_Number++;
  //???? what's going ? (asano memo)
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return EndOfAnEvent();
    if( status==2 ) return EndOfAnEvent(false);}
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number <<"  " <<Time<<"  "<<scaend[10]<<std::endl;

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  DetectorList *dlist=DetectorList::GetInstance();
  //TKO
#if TKOHIS  
  for(int ntko=0;ntko<tko->entries();ntko++)
    {
      TKOHit *tkohit=tko->hit(ntko);
      Tools::H1(Form("TKOc%ds%d",tkohit->cr(),tkohit->sl() ),tkohit->ch(), 32, -0.5,31.5 );
      Tools::SetXTitleH1(Form("TKOc%ds%d",tkohit->cr(),tkohit->sl() ),"TKO ch.");
      Tools::H1(Form("TKOc%ds%da%d",tkohit->cr(),tkohit->sl(),tkohit->ch() ) , tkohit->data(), 4196,-0.5,4195.5 );
      Tools::SetXTitleH1(Form("TKOc%ds%da%d",tkohit->cr(),tkohit->sl(),tkohit->ch() ) , "TKO ch.");
    }
#endif 
//#########

  header->Convert( tko, confMan );
  if(spillinit<0) spillinit=header->spillt(0);
  spillfini=header->spillt(0);
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  rtFile->cd();

  // Trigger Pattern
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      Tools::H1("Pattern",i,21,-0.5,20.5);
      Tools::H1(Form("TPattern%d",i),val,4196,-0.5,4195.5);
    }
  }

  //??
  for(int i=0;i<20;i++)
    if(header->trigmode(i))  Tools::H1("TrigMode",i, 21,-0.5,20.5); 

  int DAQFLAG=0;
  for(int i=0;i<20;i++)
    if(header->trigmode2(i)) //??
      {
	Tools::H1("DAQMode",i,21,-0.5,20.5); 
	DAQFLAG++;
      }
  if(DAQFLAG==0)  Tools::H1("DAQMode",18,21,-0.5,20.5); 
  if(DAQFLAG==1)  Tools::H1("DAQMode",19,21,-0.5,20.5); 
  if(DAQFLAG>1)   Tools::H1("DAQMode",20,21,-0.5,20.5); 
  
  if(!header->IsTrig(Trig_Cosmic)) return EndOfAnEvent();

  // CDC
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    int nwire = confMan->GetCDCWireMapManager()->nw(layer);
    Tools::H1( Form("MulCDC_%d",layer) , cdsMan->nCDC(layer), nwire+1, -0.5, nwire+0.5 );
    for( int i=0; i<cdsMan->nCDC(layer); i++ ){
      CDCHit *hit = cdsMan->CDC(layer,i);
      int wire = hit->wire();
      int tdc = hit->tdc();
      if(is_TDCcalib){
        double dt = hit->dt();
        Tools::H1( Form("dTCDC_%d",layer) , dt , 6000, -1500, 1500,"drift time [nsec.]" );
      }
      Tools::H1( Form("TCDC_%d",layer) , tdc, 2048,-0.6, 2047.5,"TDC Ch.");
      Tools::H1( Form("TCDC_%d_%d",layer,wire) , tdc , 2048,-0.6, 2047.5,"TDC Ch.");
      Tools::H1( Form("HitPatCDC_%d",layer) , wire , nwire, 0.5, nwire+0.5,"wire#" );
    }
  }
  
  // BeamLine Chamber
  const int ndc=6;
  int dccid[ndc]={CID_BLC1a,CID_BLC1b,
		  CID_BLC2a,CID_BLC2b,
		  CID_BPC,
		  CID_FDC1
  };

  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    for( int layer=1; layer<=nlays; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire( cid, layer );
      int tmpmul=blMan->nBLDC(cid,layer);
      Tools::H1( Form("Mul%s_%d",name,layer),tmpmul, nwire+1, -0.5, nwire+0.5 ); 
      Tools::SetXTitleH1( Form("Mul%s_%d",name,layer),"Multiplicity" ); 
      for( int i=0; i<tmpmul; i++ ){
	ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
	int wire = hit->wire();
	int tdc = hit->tdc();
	Tools::H1( Form("T%s_%d",name,layer) , tdc ,2048,-0.5,2047.5,"TDC Ch.");
	Tools::H1( Form("T%s_%d_%d",name,layer,wire) , tdc  ,2048,-0.5,2047.5,"TDC Ch.");
	Tools::H1( Form("HitPat%s_%d",name,layer) , wire , nwire, 0.5, nwire+0.5,"wire#");
  if(is_TDCcalib){
    double dt = hit->dt();
    Tools::H1( Form("d%s_%d_%d",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
  }
	if(layer%2==1)
	  {
	    for( int i2=0; i2< (blMan->nBLDC(cid,layer+1)); i2++ ){
	      ChamberLikeHit *hit2 = blMan->BLDC(cid,layer+1,i2);
	      int wire2 = hit2->wire();
	      //	      int tdc2 = hit2->tdc();
	      Tools::H2( Form("WirePat%s_%d",name,layer), wire,wire2, nwire, 0.5, nwire+0.5, nwire, 0.5, nwire+0.5 );
	      Tools::SetXTitleH2( Form("WirePat%s_%d",name,layer), Form("Layer%d wire#",layer));
	      Tools::SetYTitleH2( Form("WirePat%s_%d",name,layer), Form("Layer%d wire#",layer+1));
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
  bool CVC=false;
  int nbin=4096;
  int nbin2=256;
  double lbin=-0.5;
  double ubin=4095.5;
#if ENERGY
  int enbin=2000;
  double elbin=-5;
  double eubin=45;
#endif
  for(int ihodo=0;ihodo<nhodo;ihodo++){
    const int cid = hodoid[ihodo];
    const char* name= dlist->GetName(cid).data();
    const int nsegs= dlist->GetNsegs(cid);
    //    std::cout<<cid<<"  "<<name<<std::endl;
    int nHodo=0;
    for( int i=0; i<blMan->nHodo(cid); i++ ){
      HodoscopeLikeHit *hit = blMan->Hodoi(cid,i);
      int seg = hit->seg();
      int au = hit->adc(0), ad = hit->adc(1);
      int tu = hit->tdc(0), td = hit->tdc(1);

      Tools::H1( Form("A%sU%d",name,seg) , au , nbin, lbin, ubin, "ADC Ch.");
      Tools::H1( Form("A%sD%d",name,seg) , ad , nbin, lbin, ubin,"ADC Ch.");
      Tools::H1( Form("T%sU%d",name,seg) , tu , nbin, lbin, ubin, "TDC Ch.");
      Tools::H1( Form("T%sD%d",name,seg) , td , nbin, lbin, ubin, "TDC Ch.");
#if ENERGY
      double eu = hit->eu(), ed = hit->ed();
      Tools::H1( Form("E%sU%d",name,seg) , eu , enbin, elbin, eubin);
      Tools::H1( Form("E%sD%d",name,seg) , ed , enbin, elbin, eubin);
      if(tu>0&&tu<4095)
	Tools::H1( Form("EwT2%sU%d",name,seg) , eu , enbin, elbin, eubin, "Energy [MeV]");
      if(td>0&&td<4095)
	Tools::H1( Form("EwT2%sD%d",name,seg) , ed , enbin, elbin, eubin, "Energy [MeV]");
#endif
      if(tu>0&&tu<4095){
        Tools::H1( Form("AwT2%sU%d",name,seg) , au , nbin, lbin, ubin,"ADC Ch." );
      }
      if(td>0&&td<4095){
        Tools::H1( Form("AwT2%sD%d",name,seg) , ad , nbin, lbin, ubin, "ADC Ch.");
      }
      if(
#if RUN43
	 cid==CID_BVC
#elseif RUN46
	 cid==CID_WVC||cid==CID_HVC1
#else
	 cid==CID_WVC
#endif
	 ){
	if( tu>0 ){
	  Tools::H1( Form("AwT%sU%d",name,seg) , au , nbin, lbin, ubin,"ADC Ch.");
	  Tools::H2( Form("AT%sU%d",name,seg) , tu, au , nbin2, lbin, ubin, nbin2, lbin, ubin);
	  Tools::SetXTitleH2( Form("AT%sU%d",name,seg) , "ADC Ch.");
	  Tools::SetYTitleH2( Form("AT%sU%d",name,seg) , "TDC Ch.");
	  nHodo++;
	}else{
	  Tools::H1( Form("AwoT%sU%d",name,seg) , au , nbin, lbin, ubin,"ADC Ch.");
	}
      }else{      
	if( hit->CheckRange() ){
	  Tools::H1( Form("AwT%sU%d",name,seg) , au , nbin, lbin, ubin, "ADC Ch.");
	  Tools::H1( Form("AwT%sD%d",name,seg) , ad , nbin, lbin, ubin, "ADC Ch.");
	  Tools::H2( Form("AT%sU%d",name,seg) , tu, au , nbin2, lbin, ubin, nbin2, lbin, ubin);
	  Tools::SetXTitleH2( Form("AT%sU%d",name,seg) , "TDC Ch.");
	  Tools::SetYTitleH2( Form("AT%sU%d",name,seg) , "ADC Ch.");
	  Tools::H2( Form("AT%sD%d",name,seg) , td, ad , nbin2, lbin, ubin, nbin2, lbin, ubin);
	  Tools::SetXTitleH2( Form("AT%sD%d",name,seg) , "TDC Ch.");
	  Tools::SetYTitleH2( Form("AT%sD%d",name,seg) , "ADC Ch.");
#if ENERGY
	  Tools::H1( Form("EwT%sU%d",name,seg) , eu , enbin, elbin, eubin);
	  Tools::H1( Form("EwT%sD%d",name,seg) , ed , enbin, elbin, eubin);
#endif
	  Tools::H1( Form("HitPat%s",name) , seg , nsegs, 0.5, nsegs+0.5, "segment #"  );
	  int tmpseg=-999;
	  if(cid ==CID_BD)   tmpseg=-12+2*seg;
	  if(cid ==CID_PC)   tmpseg= 27+seg;
	  if(cid ==CID_CVC)  tmpseg=seg;
	  if(cid ==CID_CVC)  CVC=true;
	  if(tmpseg>-100){
	    Tools::H1( "HitPatCharged" , tmpseg , 74+1, -10, 64+1, "segment #"  );
	    Tools::H2( "Event_HitPatCharged" , Event_Number, tmpseg , 1000,0,1000000,74+1, -10, 64+1 );
	  }
	  nHodo++;
	  if(cid==CID_NC){
	    //int layer=-1;
	    //int sseg=-1;
	    int layer=(seg-1)/16+1;
	    int sseg=(seg-16*layer+16)%16+1;
	    if( CVC ){ 
         Tools::H2( "HitPatNC2DwC", sseg,layer, 16, 0.5, 16.5,7,0.5,7.5 );
         Tools::SetXTitleH2( "HitPatNC2DwC", "seg" );
         Tools::SetYTitleH2( "HitPatNC2DwC", "layer" );
	    }else{
        Tools::H2( "HitPatNC2DwoC", sseg,layer , 16, 0.5, 16.5,7,0.5,7.5 );
        Tools::SetXTitleH2( "HitPatNC2DwoC", "seg" );
        Tools::SetYTitleH2( "HitPatNC2DwoC", "layer" );
      }
	    Tools::H2( "HitPatNC2D" , sseg, layer, 16, 0.5, 16.5,7,0.5,7.5 );
	    Tools::SetXTitleH2( "HitPatNC2D" , "seg");
	    Tools::SetYTitleH2( "HitPatNC2D" , "layer");
	  }
	}else{
	  Tools::H1( Form("AwoT%sU%d",name,seg) , au , nbin, lbin, ubin, "ADC Ch.");
	  Tools::H1( Form("AwoT%sD%d",name,seg) , ad , nbin, lbin, ubin, "ADC Ch.");
#if ENERGY
	  Tools::H1( Form("EwoT%sU%d",name,seg) , eu , enbin, elbin, eubin);
	  Tools::H1( Form("EwoT%sD%d",name,seg) , ed , enbin, elbin, eubin);
#endif
	}
      }
    }
    Tools::H1( Form("Mul%s",name) , nHodo , nsegs+1, -0.5, nsegs+0.5, "segment #"  );
  }
  
  for( int i=0; i<blMan->nHodo(CID_T0pre); i++ ){
    HodoscopeLikeHit *hit = blMan->Hodoi(CID_T0pre,i);
    int seg = hit->seg();
    int tu = hit->tdcu(), td = hit->tdcd();
    Tools::H1( Form("TT0preU%d",seg) , tu , nbin, lbin, ubin);
    Tools::H1( Form("TT0preD%d",seg) , td , nbin, lbin, ubin);
  }
  for( int i=0; i<blMan->nHodo(CID_T0post); i++ ){
    HodoscopeLikeHit *hit = blMan->Hodoi(CID_T0post,i);
    int seg = hit->seg();
    int tu = hit->tdcu(), td = hit->tdcd();
    Tools::H1( Form("TT0postU%d",seg) , tu , nbin, lbin, ubin);
    Tools::H1( Form("TT0postD%d",seg) , td , nbin, lbin, ubin);
    HodoscopeLikeHit *hit2 = blMan->Hodo(CID_T0,seg);
    if(hit2){
      Tools::H2( Form("TT0main_postU%d",seg) , hit2->tdcu(),tu , nbin2, lbin, ubin, nbin2, lbin, ubin);
      Tools::H2( Form("TT0main_postD%d",seg) , hit2->tdcd(),td , nbin2, lbin, ubin, nbin2, lbin, ubin);
    }
  }
  for( int i=0; i<blMan->nHodo(CID_BHDpost); i++ ){
    HodoscopeLikeHit *hit = blMan->Hodoi(CID_BHDpost,i);
    int seg = hit->seg();
    int tu = hit->tdcu();
    Tools::H1( Form("TBHDpost%d",seg) , tu , nbin, lbin, ubin);
    Tools::H2( Form("TBHDmain_post%d",seg) , blMan->Hodo(CID_BHD,seg)->tdcmean(),tu , nbin2, lbin, ubin, nbin2, lbin, ubin);
  }
    
  // CDS Hodoscope
  const int ncdshodo=2;
  int cdshodoid[ncdshodo]={CID_CDH, CID_IH };
  
  for(int ihodo=0;ihodo<ncdshodo;ihodo++){
    const int cid = cdshodoid[ihodo];
    const char* name= dlist->GetName(cid).data();
    const int nsegs= dlist->GetNsegs(cid);
    int nHodo=0;
    for( int i=0; i<cdsMan->nHodo(cid); i++ ){
      HodoscopeLikeHit *hit = cdsMan->Hodoi(cid,i);
      int seg = hit->seg();
      int au = hit->adcu(), ad = hit->adcd();
      int tu = hit->tdcu(), td = hit->tdcd();
      Tools::H1( Form("A%sU%d",name,seg) , au , nbin, lbin, ubin);
      Tools::SetXTitleH1( Form("A%sU%d",name,seg) , "ADC Ch.");
      Tools::H1( Form("A%sD%d",name,seg) , ad , nbin, lbin, ubin);
      Tools::SetXTitleH1( Form("A%sD%d",name,seg) , "ADC Ch.");
      Tools::H1( Form("T%sU%d",name,seg) , tu , nbin, lbin, ubin);
      Tools::SetXTitleH1( Form("T%sU%d",name,seg) , "TDC Ch.");
      Tools::H1( Form("T%sD%d",name,seg) , td , nbin, lbin, ubin);
      Tools::SetXTitleH1( Form("T%sD%d",name,seg) , "TDC Ch.");
      if( cid!=CID_IH ){
	if( hit->CheckRange() ){
	  Tools::H1( Form("AwT%sU%d",name,seg) , au , nbin, lbin, ubin);
	  Tools::SetXTitleH1( Form("AwT%sU%d",name,seg) , "ADC Ch.");
	  Tools::H1( Form("AwT%sD%d",name,seg) , ad , nbin, lbin, ubin);
	  Tools::SetXTitleH1( Form("AwT%sD%d",name,seg) , "ADC Ch.");
	  Tools::H2( Form("AT%sU%d",name,seg) , tu, au  , nbin2, lbin, ubin, nbin2, lbin, ubin);
	  Tools::SetXTitleH2( Form("AT%sU%d",name,seg) , "TDC Ch." );
	  Tools::SetYTitleH2( Form("AT%sU%d",name,seg) , "ADC Ch." );
	  Tools::H2( Form("AT%sD%d",name,seg) , td, ad  , nbin2, lbin, ubin, nbin2, lbin, ubin);
	  Tools::SetXTitleH2( Form("AT%sD%d",name,seg) , "TDC Ch.") ;
	  Tools::SetYTitleH2( Form("AT%sD%d",name,seg) , "ADC Ch.") ;
	  Tools::H1( Form("HitPat%s",name) , seg , nsegs, 0.5, nsegs+0.5 );
	  Tools::SetXTitleH1( Form("HitPat%s",name) , "seg." );
	  nHodo++;
	}else{
	  Tools::H1( Form("AwoT%sU%d",name,seg) , au , nbin, lbin, ubin);
	  Tools::H1( Form("AwoT%sD%d",name,seg) , ad , nbin, lbin, ubin);
	}
      }else{ //IF IH, no longer used
	if( tu>0 && tu<4096 ){
	  Tools::H1( Form("AwT%sU%d",name,seg) , au , nbin, lbin, ubin);
	  Tools::H2( Form("AT%sU%d",name,seg) , tu, au , nbin2, lbin, ubin, nbin2, lbin, ubin);
	  Tools::H1( Form("HitPat%s",name) , seg , nsegs, 0.5, nsegs+0.5 );
	  nHodo++;
	}else{
	  Tools::H1( Form("AwoT%sU%d",name,seg) , au , nbin, lbin, ubin);
	}
      }
    }
    Tools::H1( Form("Mul%s",name) , nHodo , nsegs+1, -0.5, nsegs+0.5 );
  }

  // Cherenkov
  const int nchere=1;
  int chereid[nchere]={CID_AC};
  int nch[nchere]={4};
  for(int ichere=0;ichere<nchere;ichere++){
    int cid=chereid[ichere];
    const char* name=dlist->GetName(cid).data();
    for( int i=0; i<blMan->nChere(cid); i++ ){
      CherenkovLikeHit *hit = blMan->Cherei(cid,i);
      //      //      int seg = hit->seg();
      if(hit->seg()!=1) continue;
      int adcsum=0;
      for(int ich=1;ich<=nch[ichere];ich++){
	int a = hit->adc(ich);
	int t = hit->tdc(ich);
	Tools::H1( Form("A%s%d",name,ich) , a , nbin, lbin, ubin);
	Tools::SetXTitleH1( Form("A%s%d",name,ich) , "ADC Ch.");
	Tools::H1( Form("T%s%d",name,ich) , t , nbin, lbin, ubin);
	Tools::SetXTitleH1( Form("T%s%d",name,ich) , "TDC Ch.");
	adcsum+=a;
	if( hit->CheckRange(1) ){
	  Tools::H1( Form("AwT%s%d",name,ich) , a , nbin, lbin, ubin);
	  Tools::SetXTitleH1( Form("AwT%s%d",name,ich) , "ADC Ch.");
	  Tools::H2( Form("AT%s%d",name,ich) , t, a , nbin2, lbin, ubin, nbin2, lbin, ubin);
	  Tools::SetXTitleH2( Form("AT%s%d",name,ich) , "TDC Ch.");
	  Tools::SetYTitleH2( Form("AT%s%d",name,ich) , "ADC Ch.");
	}else{
	  Tools::H1( Form("AwoT%s%d",name,ich) , a , nbin, lbin, ubin);
	  Tools::SetXTitleH1( Form("AwoT%s%d",name,ich) ,"ADC Ch.");
  }
      }//for ich
      Tools::H1( Form("A%sSUM",name) , adcsum ,  nbin, lbin, ubin);
      Tools::SetXTitleH1( Form("A%sSUM",name) , "ADC Ch.");
      if(hit->CheckRange(1)){
        Tools::H1( Form("AwT%sSUM",name),  adcsum , nbin, lbin, ubin);   
        Tools::SetXTitleH1( Form("AwT%sSUM",name), "ADC Ch.");   
      }else{                   
        Tools::H1( Form("AwoT%sSUM",name), adcsum , nbin, lbin, ubin);   
        Tools::SetXTitleH1( Form("AwoT%sSUM",name), "ADC Ch.");   
      }
    }
  }

  for(int i=0;i<tko->entries();i++)
    {
      TKOHit *hit=tko->hit(i);
      int cr=hit->cr();
      int sl=hit->sl();
      int ch=hit->ch();
      if(cr==2&&sl==4&&ch==15)           Tools::H1("start1",hit->data(),nbin,lbin,ubin);
      else if(cr==2&&sl==7&&ch==11)      Tools::H1("start2",hit->data(),nbin,lbin,ubin);
      else if(cr==2&&sl==7&&ch==12)      Tools::H1("start3",hit->data(),nbin,lbin,ubin);
      else if(cr==2&&sl==7&&ch==13)      Tools::H1("start4",hit->data(),nbin,lbin,ubin);
      else if(cr==2&&sl==7&&ch==14)      Tools::H1("start5",hit->data(),nbin,lbin,ubin);
      else if(cr==2&&sl==7&&ch==15)      Tools::H1("start6",hit->data(),nbin,lbin,ubin);
      else if(cr==6&&sl==9&&ch==15)      Tools::H1("start7",hit->data(),nbin,lbin,ubin);
      else if(cr==8&&sl==1&&ch==15)      Tools::H1("start8",hit->data(),nbin,lbin,ubin);
      else if(cr==8&&sl==2&&ch==12)      Tools::H1("start9",hit->data(),nbin,lbin,ubin);
      else if(cr==8&&sl==2&&ch==15)      Tools::H1("start10",hit->data(),nbin,lbin,ubin);
      else if(cr==8&&sl==4&&ch==15)      Tools::H1("start11",hit->data(),nbin,lbin,ubin);
      else if(cr==8&&sl==5&&ch==12)      Tools::H1("start12",hit->data(),nbin,lbin,ubin);
      else if(cr==8&&sl==5&&ch==14)      Tools::H1("start13",hit->data(),nbin,lbin,ubin);
      else if(cr==8&&sl==5&&ch==15)      Tools::H1("start14",hit->data(),nbin,lbin,ubin);
      else if(cr==8&&sl==7&&ch==15)      Tools::H1("start15",hit->data(),nbin,lbin,ubin);
    }  
  return EndOfAnEvent();  
}

void EventAnalysisRawAll::Finalize()
{
  std::cout << " Enter EventAnalysisRawAll::Finalize " << std::endl;
  for(int i=0;i<40;i++){
    Tools::H1("Scaler",i,41,-0.5,40.5,scaend[i]-scainit[i]);
    std::cout<<i<<"  "<<scaend[i]<<std::endl;
  }
  TParameter<long> inittime("Time_init",t_init);
  inittime.Write();
  std::cout<<inittime.GetVal()<<std::endl;
  TParameter<long> finitime("Time_fini",Time);
  finitime.Write();
  std::cout<<finitime.GetVal()<<std::endl;
  TParameter<int> initspill("Spill_init",spillinit);
  initspill.Write();
  TParameter<int> finispill("Spill_fini",spillfini);
  finispill.Write();

  rtFile->cd();
  gFile->Write();
  gFile->Close();
  //  gSystem->Sleep(10000);
  //  gSystem->Exec(Form("mv %s %s",tmpname.Data(),confMan->GetOutFileName().c_str()));
  delete blMan;
  delete cdsMan;
  delete header;
}

bool EventAnalysisRawAll::EndOfAnEvent(bool flag){
  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return flag;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisRawAll *event = new EventAnalysisRawAll();
  return (EventTemp*)event;
}
