#include <iostream>

#include "FPAnaMan.h"

using namespace std;

int main(int argc, char** argv)
{
  if( argc!=2 ){
    cout<<"Please input conffile"<<endl;
    return 0;
  }

  FPAnaMan *anaMan = new FPAnaMan();
  anaMan->init(argv[1]);
  cout<<"Initialization finish"<<endl;
  anaMan-> finit();

  return 0;
}
