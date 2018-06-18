#include "ConfMan.h"

int main(int argc, char** argv)
{
  if( argc!=2 ) return 0;

  int runnum=atoi(argv[1]);

  ConfMan *conf = new ConfMan("conf/Run68/analyzer.conf", runnum);
  conf-> Initialize();
  std::cout<<"runnum : "<<runnum<<" USWK current : "<<conf->GetUshiwakaCurrent()<<std::endl;

  return 0;
}
