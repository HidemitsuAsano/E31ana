#include "IMPiSigmaUtil.h"
#include "Tools.h"

bool Util::EveSelectCDHMul(CDSHitMan *cdsman)
{
  //** # of CDH-hits cut **//
  int nCDH = 0;
  for( int i=0; i<cdsman->nCDH(); i++ ) {
    Tools::Fill2D(Form("CDHtime"),cdsman->CDH(i)->seg(),cdsman->CDH(i)->ctmean());
    //if( cdsman->CDH(i)->CheckRange() ) nCDH++; //** only requirement of TDC **//
    if( cdsman->CDH(i)->CheckRange() && cdsman->CDH(i)->ctmean()<cdscuts::tdc_cdh_max ) {
      nCDH++;
    }
  }
  Tools::Fill1D( Form("mul_CDH"), nCDH );
  if( nCDH != cdscuts::cdhmulti  ) {
    return false;
  }

  return true;
}
