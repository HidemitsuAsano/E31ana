#include "FitConf.h"
#include <vector>

using namespace std;

int main(int argc, char** argv)
{
  cout<<"===== d(K-, n) Fitting start ====="<<endl;
  FitConf *conf=new  FitConf();
  if( argc==2 ) conf-> init(argv[1]);

  conf-> finit();
  cout<<"===== d(K-, n) Fitting finish ====="<<endl;
}
