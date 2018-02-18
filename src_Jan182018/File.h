// File.h

#ifndef File_h
#define File_h 1

#include <string>
#include <cstdio>

#include "ConfMan.h"
//class ConfMan;

class File
{
 private:
  std::string InputFileName;
  //ConfMan *confManager;
 public:
  File( const std::string &InputFileName );
  ~File();

 private:
  File();
  File( const File & );
  File & operator = ( const File & );

 public:
  void SetInputFileName( const std::string infile );

  //bool Processing( void );
  bool Processing( ConfMan *confMan );

 private:
  int RunNum, EventNum;
  int ReadOneBlockEvent( unsigned int *buf, FILE *fp );
  int ReadOneEvent( unsigned int *buf, FILE *fp );
  unsigned int *EventBuffer;
  void PrintBuffer( unsigned int *buf );
};

#endif
