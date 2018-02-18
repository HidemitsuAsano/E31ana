// CDHCluster.cpp
//
#include "CDHCluster.h"

void CDHCluster::Calc()
{
	EMean=0; CTMean=99999;
	for(int i=0; i<CDHContainer.size(); i++){
		EMean+=CDHContainer[i].emean();
		if(CTMean>CDHContainer[i].ctmean()){
			CTMean=CDHContainer[i].ctmean();
			HitPos=CDHContainer[i].pos();
			HitPosition=CDHContainer[i].hitpos();
		}
	}
}

void CDHCluster::Clear()
{
	CDHContainer.clear();
}
