#ifndef MYHISTFC_H
#define MYHISTFC_H 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistFC();
void fillFC(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo,
	    ForwardChargeInfo *fcInfo);

void fillFC_gamma(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo);

double calc_offset(const double &mm, const TLorentzVector tot_lmom, const TLorentzVector lmom, AnaInfo *anaInfo);

bool kd_pLpim(AnaInfo *anaInfo);
bool kd_pS0pim(AnaInfo *anaInfo);


#endif
