#include "mktree.h"

#include <string>
#include <iostream>
#include <new>
#include <signal.h>

#include "TFile.h"
#include "TTree.h"

#include "ConfMan.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "TKO.h"
#include "EventHeader.h"
#include "ScalerMan.h"

std::string conffilename = "conf/analyzer.conf";
std::string rootfilename = "tmp.root";
TFile *rtfile;
TTree *evtree;
TTree *scatree;
ConfMan *conf;
TKOHitCollection *tkoCol;
CDSHitMan *cdsMan;
BeamLineHitMan *blMan;
EventHeader *header;
ScalerMan *scaMan;

//const int MaxTreeSize = 1900000000000;
const Long64_t MaxTreeSize = 190000000000;

using namespace std;

extern "C" int cmain( int, char** );

int main( int argc, char **argv )
{
  if( 3<argc ){
    conffilename = argv[1];
    rootfilename = argv[2];
  }
  else{
    std::cout << " please type as follows," << std::endl
	      << " >a.out conffile out.root data.dat" << std::endl;
    return 0;
  }

  cmain(argc, argv );

  return 0;
}

int uopt( int i )
{
  return 0;
}

void uintr( int isig )
{
  cout << "signal handler" << endl;
  ufini();
  exit(0);
}

int uinit()
{
  cout << " Enter Initialize ... " << endl;

  (void)signal(SIGHUP,uintr);

  rtfile = new TFile( rootfilename.c_str(), "recreate" );
  //rtfile->SetBufferSize(1024*1000);
  evtree = new TTree( "EventTree", "EventTree" );
  scatree = new TTree( "ScalerTree", "ScalerTree" );

//   evtree->SetMaxTreeSize( MaxTreeSize );
//   scatree->SetMaxTreeSize( MaxTreeSize );

  conf = new ConfMan( conffilename.c_str() );
  conf->Initialize();

  tkoCol = new TKOHitCollection();
  if( tkoCol==NULL ){ cerr << "!!!!" << endl; return 0; }
  //evtree->Branch( "TKOHitCol", &tkoCol );

  header = new EventHeader();
  if( header==NULL ){ cerr << "!!!!" << endl; return 0; }
  evtree->Branch( "EventHeader", &header );

  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ) { cerr << "!!!!" << endl; return 0; }
  evtree->Branch( "CDSHitMan", &cdsMan );

  blMan = new BeamLineHitMan();
  if( blMan==NULL ){ cerr << "!!!!" << endl; return 0; }
  evtree->Branch( "BeamLineHitMan", &blMan );

  scaMan = new ScalerMan();
  if( scaMan==NULL ){ cerr << "!!!!" << endl; return 0; }
  scatree->Branch( "ScalerMan", &scaMan );

  return 0;
}

static int event_number = 0; // !! temporary !!
static int block_event_number = 0; // !! temporary !!

int usca( int nsca, unsigned int *sca )
{
  block_event_number++;
  header->SetBlockEventNumber( block_event_number );
  scaMan->SetBlockEventNumber( block_event_number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], conf->GetScalerMapManager()->GetName(i) );
  }
#if 0
  std::cout << nsca<< std::endl;
  for( int i=0; i<nsca; i++ ){
    std::cout << "  " << sca[i];
  }
  std::cout << std::endl;
#endif
  scatree->Fill();
  scaMan->Clear();
  return 0;
}

//#define CHKCHK 1

int uana( int num_word, s_raw *raw )
{
  event_number++;
  if( event_number%5000==0 )
    cout << " Event# : " << event_number << endl;

  //if( 350000<event_number ) return 0;

  header->SetRunNumber(0);
  header->SetEventNumber(event_number);

  for( int i=0; i<num_word; i++ ){
#ifdef CHKCHK
    int layer, seg; // wire or seg
    int cid, at, ud;
    conf->GetCounterMapManager()->GetInfo(raw[i].c,raw[i].m,raw[i].s,cid,layer,seg,at,ud);
    if( cid == CID_CDC )
      conf->GetCDCWireMapManager()->GetWire(raw[i].c,raw[i].m,raw[i].s,layer,seg);

    //if( raw[i].c == 2 && ( ( raw[i].m==7 && 11<=raw[i].s && raw[i].s <=15 ) || ( raw[i].m == 18 && 17<=raw[i].s && raw[i].s <=21 ) )  )
    //if( raw[i].c==2 && ( raw[i].m==7 || raw[i].m==18 ) )
    //if( raw[i].c>2 )
      std::cout << raw[i].c << "  " <<  raw[i].m << "  " << raw[i].s << "  " <<  raw[i].v << "  "
		<< layer << "  " << seg
		<< std::endl;
#else
    TKOHit *hit = new TKOHit( raw[i].c, raw[i].m, raw[i].s, raw[i].v );
    tkoCol->AddHit( *hit );
    delete hit;
#endif
  }
  
#if 0
  std::cout << " from now, converting ... " << std::endl;
#endif
#ifndef CHKCHK
  header->Convert( tkoCol, conf );
  blMan->Convert( tkoCol, conf );
  cdsMan->Convert( tkoCol, conf );
#endif
#if 0
  std::cout << " convert finish! " << std::endl;
#endif

  evtree->Fill();
  tkoCol->Clear();
  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return 0;
}

int ufini()
{
  cout << " Enter fInitialize ... " << endl;

  gFile->Write();
  gFile->Close();

  delete conf;  
  delete blMan;
  delete cdsMan;
  delete tkoCol;
  delete header;
  
  return 0;
}

