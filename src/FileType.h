// FileType.h

#ifndef FileType_h
#define FileType_h 1

enum FileCompressedType { NotExist=0, NotCompressed, Bzip2Compressed, GzipCompressed };

FileCompressedType FileType( const char *filename );

#endif
