#include "FNAnaMan.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
  cout<<"===== anaFN START ====="<<endl;
  if( argc!=2 ){
    cout<<"Please input $(conffile)"<<endl;
    return 0;
  }
  FNAnaMan *anaMan=new FNAnaMan();
  anaMan-> init(argv[1]);
  anaMan-> fitKNpi_all();
  anaMan-> write("ite_0");

  for( int ite=1; ite<=5; ite++ ){
    anaMan-> fillHistMC();
    anaMan-> fitAll();
    anaMan-> write(Form("ite_%d", ite));
  }

  // anaMan-> fillHistMC();
  // anaMan-> fitNpipiIM_K0(5);
  // anaMan-> write("ite_fitK0");
  // anaMan-> fillHistMC();
  anaMan-> postAna();
  anaMan-> finit();
  cout<<"===== anaFN FINISH ====="<<endl;
  return 0;
}
