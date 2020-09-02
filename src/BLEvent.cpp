#include "BLEvent.h"

#include <string>
#include <iostream>
#include <new>
#include <signal.h>

#include "TFile.h"
#include "TTree.h"

ConfMan *conf;
TKOHitCollection *tkoCol;
EventAlloc *evAlloc;
EventTemp *evTemp;

using namespace std;

int BLEvent::BMain( int argc, char **argv, ConfMan *con )
{
  if( 3<argc ){
    conf = confMan = con;
  }
  else{
    std::cout << " please type as follows," << std::endl
	      << " >a.out conffile out.root data.dat (vme.root)" << std::endl;
    return 0;
  }
  int nsmp = con->GetCounterMapManager()->GetNumSMP();
  int nch_sca = con->GetCounterMapManager()->GetNumScaler();
  for(int i=0;i<nsmp;i++){
    for(int j=1;j<=23;j++){
      set_crate_type(i,j,con->GetCounterMapManager()->GetCrateType(i,j));  //sl 1 origin
    }
  }
  set_num_modules( nch_sca, nsmp );
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

  tkoCol = new TKOHitCollection();
  evAlloc = new EventAlloc();
  evTemp = evAlloc->EventAllocator();
  evTemp->Initialize(conf);

  return 0;
}

int usca( int nsca, unsigned int *sca )
{
  evTemp->USca( nsca, sca );

  return 0;
}

int utime( int time )
{
  evTemp->UTime( time );

  return 0;
}

//#define CHKCHK 1

int uana( int num_word, s_raw *raw )
{
  
  for( int i=0; i<num_word; i++ ){
#ifdef CHKCHK
    int layer, seg; // wire or seg
    int cid, at, ud;
    conf->GetCounterMapManager()->GetInfo(raw[i].c,raw[i].m,raw[i].s,cid,layer,seg,at,ud);
    if( cid == CID_CDC )
      conf->GetCDCWireMapManager()->GetWire(raw[i].c,raw[i].m,raw[i].s,layer,seg);

      std::cout << raw[i].c << "  " <<  raw[i].m << "  " << raw[i].s << "  " <<  raw[i].v << "  "
		<< layer << "  " << seg
		<< std::endl;
#else
    TKOHit *hit = new TKOHit( raw[i].c, raw[i].m, raw[i].s, raw[i].v );
    tkoCol->AddHit( *hit );
    delete hit;
#endif
  }

  int status;
  if( evTemp->UAna( tkoCol ) ) status = 1;
  else                         status = -1;

  tkoCol->Clear();
  return status;
}

int ufini()
{
  evTemp->Finalize();

  delete tkoCol;
  delete evTemp;
  delete evAlloc;
  return 0;
}

