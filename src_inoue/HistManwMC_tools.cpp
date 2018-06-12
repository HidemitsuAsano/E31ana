#include "HistManwMC.h"

void HistManwMC::printHit()
{
  std::cout<<"> T0 hit  : "<<blMan->nT0()<<" w/ TDC : "<<fT0_hit.size()<<std::endl;
  std::cout<<"> CVC hit : "<<blMan->nCVC()<<" w/ TDC : "<<fCVC_hit.size()<<std::endl;
  std::cout<<"> BVC hit : "<<blMan->nBVC()<<" w/ TDC : "<<fBVC_hit.size()<<std::endl;
  std::cout<<"> NC hit  : "<<blMan->nNC()<<" w/ TDC : "<<nNC()<<std::endl;

}
