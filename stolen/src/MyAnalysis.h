// MyAnalysis.h
// 2014.12.22
// T. Yamaga

#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

bool Analysis(CDSHitMan* cdsMan,BeamLineHitMan* blMan,EventHeader* header, ScalerMan* scaMan){

	// Multihit of CDH and IH //
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
	if(!header->IsTrig(Trig_Cosmic) || MulCDH!=2 || MulIH!=2){
		return true;
	}

	return true;
}
