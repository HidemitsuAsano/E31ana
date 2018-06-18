#include "MyHistTrig.h"

void initHistTrig()
{
  new TH1F("TrigMode0", "Trigger Mode",  20,  0.5, 20.5);
  new TH1F("TrigMode",  "Trigger Mode",  20,  0.5, 20.5);
  new TH1F("TrigMode2", "Trigger Mode2", 20,  0.5, 20.5);
  new TH1F("nTrig",     "nTrig",         21, -0.5, 20.5);

  for( int i=0; i<20; i++ ){
    new TH1F(Form("trig%d_TDC", i), Form("Trigger %d TDC", i), 4000, 0, 4000);
  }
}

void fillTrig(EventHeader *header)
{
  for( int i=0; i<20; i++ ){
    MyHistTools::fillTH(Form("trig%d_TDC", i), header->pattern(i));
  }

  int nTrig=0;
  for( int i=0; i<20; i++ ){
    MyHistTools::fillTH("TrigMode0", header->trigmode());
    if( header->trigmode(i)  ){
      MyHistTools::fillTH("TrigMode", i);
      nTrig++;
    }
    if( header->trigmode2(i) )  MyHistTools::fillTH("TrigMode2", i);
  }
  MyHistTools::fillTH("nTrig", nTrig);
}
