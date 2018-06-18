// VEvent.cpp

#include <cstdio>
#include <string>
#include <new>
#include <cctype>
#include <cstddef>
#include <iostream>
#include <iomanip>

#include "VEvent.h"
// #include "DataForm.h"

bool VEvent::EventCheck()
{
  //bool retval = false;
  //unsigned int *ptr;
  
  Size = Buf[0];
  RecordType = Buf[1];
  RunNum = Buf[2];
  BlockEventNum = Buf[3];
  //EventNum = BlockEventNum;
  EventNum = EventNumLast;

  if( Size==0 ){
    std::cerr << "Size 0 !! data broken?" << std::endl;
    exit(-1);
  }

#if 0
  std::cout << std::endl << " Size:" << Size
	    << " RecordType:" << RecordType
	    << " RunNum:" << RunNum
	    << " BlockEventNum:" << BlockEventNum
	    << std::endl;
#endif
  return true;
}

void VEvent::SetBuffer( unsigned int *buf, int evnum )
{
  Buf = buf;
  EventNumLast = evnum;
}

bool VEvent::Processing3( ConfMan *confMan )
{
  bool retval = false;
  unsigned int *ptr;
  
  Size = Buf[0];
  RecordType = Buf[1];
  RunNum = Buf[2];
  BlockEventNum = Buf[3];
  //EventNum = BlockEventNum;
  EventNum = EventNumLast;

  if( Size==0 ){
    std::cerr << "Size 0 !! data broken?" << std::endl;
    exit(-1);
  }

#if 0
  std::cout << std::endl << " Size:" << Size
	    << " RecordType:" << RecordType
	    << " RunNum:" << RunNum
	    << " BlockEventNum:" << BlockEventNum
	    << std::endl;
#endif

  if( RecordType==UniRecTypeNormal ){
    ptr = &Buf[6];
    int length=6;
    int size_read_max=0;

    EventHeaderPointerContainer.clear();
    EventBufferPointerContainer.clear();

    if( (int)((*ptr)&Data_Header)!=Data_Header ){
      std::cerr << "Undefined Data Structure!!!"
		<< std::hex << *ptr << std::dec << std::endl;
      return 0;
    }
    int NumOfSMP = (*ptr)&Data_NumSMPMask; ptr++; length++;
    
    for( int i=0; i<NumOfSMP; i++ ){

      if( (int)((*ptr)&SMP_Header) != SMP_Header ){ // just check
	std::cerr << " something wrong ... " << std::endl;
      }

      EventHeaderPointerContainer.push_back(ptr);
      EventBufferPointerContainer.push_back(ptr+2);
      unsigned int *ptmp; ptmp = ptr;
      ptmp++;
      int nhit = (*ptmp)&SMP_NumDataMask;//SD_GetNumEv(*ptmp);
      size_read_max += nhit;
#if 0
      std::cout << " nhit : " << nhit
		<< std::hex << "  " << *ptr << "  " << *(ptr+1) << "  " << *(ptr+2) << std::dec
		<< std::endl;
      int nh = SD_GetNumEv( *(ptr+2) );
      std::cout << std::hex << *(ptr+2+nh) << "  " << *(ptr+2+nh+1) << std::dec << std::endl;
#endif
      ptr += nhit+3;
      if( (int)((*(ptr-1))&SMP_Footer) != SMP_Footer ){
	std::cout << " ################## " << std::endl;
	std::cout << " invalid : " << std::hex << *(ptr-1) << std::dec << std::endl;
	std::cout << " ################## " << std::endl;
      }
    }

    // scan to check
    {
      const unsigned int *ph=0, *pb=0;
      int nbuf = EventHeaderPointerContainer.size();
      for( int i=0; i<nbuf; i++ ){
	ph = EventHeaderPointerContainer[i];
	pb = EventBufferPointerContainer[i];

	const unsigned int *ph_scan = ph;
	if( *ph_scan!=SMP_Header ) std::cout << " !! " << std::endl;
	int num_footer = 0;
	while( *ph_scan != SMP_Footer ){
	  if( ((*ph_scan)&SD_ModuleDownBit)==SD_ModuleDownBit ) num_footer++;
	  ph_scan++;
	}
#if 0
	std::cout << " num_footer  : " << num_footer << std::endl;
#endif
	if( EventNum==0 && i==0 ){
	  NumEventInBlock = num_footer;
	}
	else{
	  if( NumEventInBlock != num_footer ){
#if 0
	    std::cout << " !! anomalous data !! This block(#=" << BlockEventNum << " ) is skipped! " << std::endl;
#endif
	    return false;
	  }
	}
      }
    }

    int nbuffer = EventBufferPointerContainer.size();
    while(1){
      retval = ProcessingNormal(confMan);
      EventNum++;

      int counter = 0;
      while( counter<nbuffer ){
	unsigned int *ptmp = EventBufferPointerContainer[0];
	int nh = SD_GetNumEv( (*ptmp) );
	ptmp += nh+1;
	EventBufferPointerContainer.erase( EventBufferPointerContainer.begin() );
	EventBufferPointerContainer.push_back(ptmp);
	counter++;
      }
      unsigned int *ptmp = EventBufferPointerContainer[0];
      if( *ptmp == SMP_Footer ) break;
    }

    EventHeaderPointerContainer.clear();
    EventBufferPointerContainer.clear();
  }
  else if( RecordType==UniRecTypeBegin ){
    retval = ProcessingBegin(confMan);
  }
  else if( RecordType==UniRecTypePause ){
    retval = ProcessingPause();
  }
  else if( RecordType==UniRecTypeResume ){
    retval = ProcessingResume();
  }
  else if( RecordType==UniRecTypeEnd ){
    retval = ProcessingEnd(confMan);
  }
  else{
    retval = ProcessingUnknown();
  }
  return retval;
}

bool VEvent::ProcessingNormal( ConfMan *confMan )
{
  header->SetRunNumber(0);
  header->SetEventNumber(EventNum);

  if( EventNum%10000==0 ){
    std::cout << " Event# : " << EventNum << std::endl;
  }

  const unsigned int *ph=0, *pb=0;
  int nbuf = EventHeaderPointerContainer.size();
  for( int i=0; i<nbuf; i++ ){
    ph = EventHeaderPointerContainer[i];
    pb = EventBufferPointerContainer[i];
    int address = (*(ph+1))&SMP_AddressMask;
    int ndat = SD_GetNumEv(*(pb)); pb++;
    int ndat_org = ndat;
    ndat--;
    int c = confMan->GetCounterMapManager()->GetCrateNum(address);
    while(0<ndat && SD_ModuleDownBitMask((*pb)) != SD_ModuleDownBit ){
      int n = SD_GetMA(*(pb));
      int a = SD_GetSA(*(pb));
      int data = SD_GetData(*(pb));
      int cid = confMan->GetCounterMapManager()->GetCID(c,n,a);
      if( cid<0 ){
	std::cout << " unknown cid is found. this event(#=" << EventNum << " ) is skipped." << std::endl;
	std::cout << " Block : " << BlockEventNum << "  EventNum : " << EventNum << std::endl;
	std::cout << " address : " << std::hex << address << std::dec << std::endl;
	std::cout << " ndat : " << ndat << "/" << ndat_org << std::endl;
	std::cout << " *pb : " << std::hex << *EventBufferPointerContainer[i] << std::dec << std::endl;
	std::cout << " DownBit : " << std::hex << *(pb-1) << "  " << *pb << "  " << *(pb+ndat-1) << std::dec << std::endl;
	tkoCol->Clear();
	header->Clear();
	cdsMan->Clear();
	return true;
      }

      if( cid == CID_CDC ){
	a = DrT_GetCH(data);
	data = DrT_GetData(data);
      }
#if 0
      std::cout << " *pb:" << std::hex << *pb << std::dec << " c:" << c << " n:" << n << " a:" << a << " data:" << data << std::endl;
#endif
      TKOHit *hit = new TKOHit( c, n, a, data );
      tkoCol->AddHit( *hit );
      delete hit;
      
      pb++; ndat--;
    }
  }

  header->Convert( tkoCol, confMan );
  cdsMan->Convert( tkoCol, confMan );
//   cdsMan->RemoveNoTDCData();
//   if( 50<cdsMan->nCDC() ){
//     cdsMan->Clear();
//   }

  cdstree->Fill();
  tkoCol->Clear();
  header->Clear();
  cdsMan->Clear();

#if 0
  std::cout << "  cdstree:" << cdstree
	    << "  tkoCol:" << tkoCol
	    << "  header:" << header
	    << "  cdsMan:" << cdsMan
	    << std::endl;
#endif
  return true;
}

bool VEvent::ProcessingBegin( ConfMan *confMan )
{
  std::cout << "=== Begin Event ===> RUN#"
            << std::setw(10) <<RunNum << std::endl;
  PrintComment();

  rtfile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  cdstree = new TTree( "CDSTree", "CDSTree" );
  
  tkoCol = new TKOHitCollection();
  if( tkoCol==NULL ){ std::cerr << "!!!!" << std::endl; return false; }
  cdstree->Branch( "TKOHitCol", &tkoCol );

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return false; }
  //cdstree->Branch( "EventHeader", &header );

  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return false; }
  //cdstree->Branch( "CDSHitMan", &cdsMan );
  
  return true;
}

bool VEvent::ProcessingPause( void )
{
  std::cout << "=== Pause Event ===> RUN#"
            << std::setw(10) << RunNum  << std::endl;
  PrintComment();
  return true;
}

bool VEvent::ProcessingResume( void )
{
  std::cout << "=== Resume Event ===> RUN#"
            << std::setw(10) << RunNum  << std::endl;
  PrintComment();
  return true;
}

bool VEvent::ProcessingEnd( ConfMan *confMan )
{
  std::cout << "=== End Event ===> RUN#"
            << std::setw(10) << RunNum << std::endl;
  PrintComment();

  gFile->Write();
  gFile->Close();

  delete cdsMan;
  delete tkoCol;
  delete header;
  
  return true;
}

bool VEvent::ProcessingUnknown( void )
{
  std::cerr << "=== Unknown Event Type ["
            << std::setw(8) << RecordType << "]" << std::endl;
  return false;
}

void VEvent::PrintComment()
{
  bool f = true;
  char c0,c1,c2,c3;

  for( int i=6; i<Size; i++ ){
    c0 = CHAR0(Buf[i]); c1 = CHAR1(Buf[i]);
    c2 = CHAR2(Buf[i]); c3 = CHAR3(Buf[i]);

    if( c0==0 ){
      if(f){ std::cout << std::endl; f = false; }
    }
    else if( isprint(c0) ){ std::cout << c0; f = true; }
    if( c1==0 ){
      if(f){ std::cout << std::endl; f = false; }
    }
    else if( isprint(c1) ){ std::cout << c1; f = true; }
    if( c2==0 ){
      if(f){ std::cout << std::endl; f = false; }
    }
    else if( isprint(c2) ){ std::cout << c2; f = true; }
    if( c3==0 ){
      if(f){ std::cout << std::endl; f = false; }
    }
    else if( isprint(c3) ){ std::cout << c3; f = true; }
  }
  std::cout << std::endl;
}
