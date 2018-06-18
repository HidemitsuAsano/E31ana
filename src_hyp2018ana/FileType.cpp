// FileType.cpp

#include <cstdio>
#include <cstdlib>
#include <string>

#include "FileType.h"

FileCompressedType FileType( const char *filename )
{
  FILE *fp;

  if( (fp=fopen(filename,"r"))==0 ){
    return NotExist;
  }

  FileCompressedType type = NotCompressed;

  std::string com1 = "file ";
  com1 += filename; com1 += "| grep bzip2"; com1 += " > /dev/null";
  std::string com2 = "file ";
  com2 += filename; com2 += "| grep gzip";  com2 += " > /dev/null";

  if( system(com1.c_str())==0 ) type = Bzip2Compressed;
  else if( system(com2.c_str())==0 ) type = GzipCompressed;

  return type;
}
