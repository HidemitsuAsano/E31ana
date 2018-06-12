#include <iostream>

#include "ScalerMapMan.h"

ClassImp(ScalerMapMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ScalerMapMan::ScalerMapMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ScalerMapMan::ScalerMapMan( const std::string &name )
{
  FileName = name;
  NumModules = 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ScalerMapMan::SetFileName( const std::string &name )
{
  FileName = name;
  NumModules = 0;
}

const int MAXCHAR = 144;
bool ScalerMapMan::Initialize()
{
  static const std::string funcname = "ScalerMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";

  int i;
  char name[MAXCHAR];
  char str[MAXCHAR];
  FILE *fp;

  if( (fp=fopen(FileName.c_str(), "r"))==0 ){
    std::cerr << " File open fail. [" << FileName << "]" << std::endl;
    exit(-1);
  }
  
  Clear();

  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;

    if( sscanf(str,"NSCA: %d", &i)==1 ){
      NumModules = i;
    }
    if( sscanf(str,"%d %s", &i, name)==2 ){
      ScalerContainer.push_back(name);
    }
    else{
      std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
      std::cerr << std::string(str) << std::endl;
    }
  }

  fclose(fp);
  
  //PrintMap();

  //  std::cout << "[" << funcname << "] Initialization finish." << std::endl;
 std::cout << " finish." << std::endl;
  return true;
}

std::string ScalerMapMan::GetName( const int &i )
{
  if( 0<=i && i<(int)ScalerContainer.size() ){
    return ScalerContainer[i];
  }
  else{
    return "";
  }
}

void ScalerMapMan::PrintMap()
{
  for( int i=0; i<(int)ScalerContainer.size(); i++ ){
    std::cout << i << "  " << ScalerContainer[i] << std::endl;
  }
}

void ScalerMapMan::Clear()
{
  ScalerContainer.clear();
}
