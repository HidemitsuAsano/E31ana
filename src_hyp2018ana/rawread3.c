static float rawread_version=3.0;
/*
 * V1.0 simple and short program to read E471, E549, E570 data structure.
 *      Please modify e471_tko.h to change number of scaler channels.
 * V1.1 Policy modified in handling wrong data structure from loose to strict.
 * V1.2 small modifications on variable names etc..
 * V1.3 Modified rawread.c to suit to e471 data
 *      Modified analysis.h to add file_prefix
 * V1.4 e471, e549, e570 - auto data type recognition
 * V1.5 bug fix
 * V1.6 E570 Run 61 -> 62 scaler added as seen in E570 Logbook 2 p41.
 * V1.7 Added feature of setting number of events to analyze.
 * V1.8 Error handling loosened. Number of scalers for runs 82-84 of E570.
 * V1.9 uevent_id function added.
 * V2.0 Event numbering scheme changed to fit to TS's routine, hopefully 
 * V2.1 Origin of the event number changed from 0 to 1
 * V2.2 header file added to suit to gcc4
 * V3.0 change file open and read function open() ,read() => fopen(), fread();and readable bzip2 and gzip file
 */
//#define DEBUG
/* UNIDAQ READOUT ROUTINE */
//
// Please prepare 
//   uinit, ufini, uana, and uopt functions.
//   Samples can be found in user.c
//    
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "e471_tko.h"
#include "analysis.h"

#define DOWNBIT 0
#define BLEN (65536*4)
#define HLEN 6
//#define HLEN 42 // 2010/09/27 by k.tsukada

extern char *optarg;
extern int optind, opterr, optopt;

char file[128];
FILE *fd;
int unihead[HLEN];
int head[MAXSMP][TKO_HEADER_LEN];
unsigned int* buf[MAXSMP];
unsigned int scaler[NSCALER];
int bs[MAXSMP]; // buffer size

static int crate_type[MAXSMP][23];

static int number_of_smp = MAXSMP;
static int number_of_scaler_channels = NSCALER;
static int itype;

int ana_size=-1;
int ana_event=-1;
unsigned long event_number=1;
int end_of_analysis=0;
int nw;
s_raw raw[BLEN];
s_header subheader[MAXSMP]; // <-- kore tsukatte inai dou siyou .

static int number;
static s_run run;

//////////////////////////////////////////////////////////////////

// External User Functions
int uana(int, s_raw*); // Routine called for every event.
int uinit();   // Routine called at the beginning of analysis.
int ufini();   // Routine called at the end of analysis.
int uopt(int); // Routine called at beginning/end of a buffer. 
int usca(int nsca, unsigned int *sca); // added by tsukada '09/11/26
int utime(int time); // added by hashimoto '11/09/12
int set_num_modules( int nscach, int nsmp ); // added by tsukada '10/10/04
int set_crate_type( int ncr, int nsl, int type);
// Interface Functions
int uscaler(int); // returns i-th scaler value
int uhead(int,int); // return i-th crate's j-th header value
s_run info_run();   // return run information
unsigned long uevent_id();

//////////////////////////////////////////////////////////////////

unsigned long uevent_id()
{
  return(event_number);
}

s_run info_run()
{
  static s_run run;

  run.number = number;
  run.type = itype;
  //  strncpy(run.comment, file_prefix, 4);
  strncpy(run.fname, file, 80 * sizeof(char));
  //  fprintf(stderr,"$I rawread: run number -> %d\n",number);
  return (run);
}

int uscaler(int i)
{
  if (i<number_of_scaler_channels){
    return(scaler[i]);
  }else{
    return(-1);
  }
}

int how_many_smps()
{
  return(number_of_smp);
}

int how_many_scaler_channels()
{
  return(number_of_scaler_channels);
}

int uhead(int cr, int i)
{
  if (cr<MAXSMP && i<TKO_HEADER_LEN){
    return(head[cr][i]);
  }else{
    return(-1);
  }
}

static int compose_an_event()
{
  int i,j;
  unsigned int cr,ma,sa,val;
  unsigned int bp[MAXSMP]; // buffer pointer
  static unsigned int prev_serial[MAXSMP];
  static unsigned int serial[MAXSMP];
  unsigned int nword[MAXSMP];
  int cond1, cond2;
  int module_bit_status;
  int status;

  for(cr=0;cr<number_of_smp;cr++){
    bp[cr]=0;
    prev_serial[cr]=-1;
  }

  for(;;){
    nw = 0;
    module_bit_status = 0;
    for(cr=0;cr<number_of_smp;cr++){
      serial[cr]= 0xffff & (buf[cr][ bp[cr] ]>>16);
      nword[cr] = 0xffff & buf[cr][ bp[cr] ];

      if (prev_serial[cr]!=-1 && prev_serial[cr] < serial[cr]){
	if (prev_serial[cr]+1 != serial[cr]){
	  fprintf(stderr,"$E Serial number mismatch. Skipping. Data may be broken.\n");
	  return(1);
	}
      }
      prev_serial[cr] = serial[cr];

      for(i=1;i<nword[cr];i++){ // i=0 -> serial<<16 | nword 
	raw[nw].c = cr;
	raw[nw].m  = (buf[cr][bp[cr]+i] >> 27) & 0x1f;
	if (crate_type[cr][raw[nw].m-1]==NORMAL) {
	  raw[nw].s  = (buf[cr][bp[cr]+i] >> 16) & 0xff;
	  raw[nw].v =   buf[cr][bp[cr]+i]        & 0xffff;
	  //fprintf(stdout,"Normal buf[%d][%2d]:%8x  m:%2d s:%2d v:%12d\n",cr,i,buf[cr][bp[cr]+i], raw[nw].m, raw[nw].s, raw[nw].v);
	}else{
	  raw[nw].s  = (buf[cr][bp[cr]+i] >> 11) & 0x1f;
	  raw[nw].v =   buf[cr][bp[cr]+i]        & 0x7ff;
	  //fprintf(stdout,"DRT    buf[%d][%2d]:%8x  m:%2d s:%2d v:%12d\n",cr,i,buf[cr][bp[cr]+i], raw[nw].m, raw[nw].s, raw[nw].v);
	}
	nw++;
	if (nw > BLEN){
	  fprintf(stderr,"$E Size limit exceeded. Skipping.\n");
	  return(2);
	}
	// printf("%d:%d->%x\n",cr,i, buf[cr][bp[cr]+i] );
      }

      if (buf[cr][bp[cr]+nword[cr]]!=0xff000000){
	fprintf(stderr,"$E rawread: DATA may be broken 1. Maybe module down bit in cr %d. Skipping.\n",cr);
	//	return(1);
#if DOWNBIT
	for(i=1;i<nword[cr];i++){ // i=0 -> serial<<16 | nword 
	  raw[nw].c = cr;
	  raw[nw].m  = (buf[cr][bp[cr]+i] >> 27) & 0x1f;
	  if (crate_type[cr][raw[nw].m-1]==NORMAL) {
	    raw[nw].s  = (buf[cr][bp[cr]+i] >> 16) & 0xff;
	    raw[nw].v =   buf[cr][bp[cr]+i]        & 0xffff;
	    //fprintf(stdout,"Normal buf[%d][%2d]:%8x  m:%2d s:%2d v:%12d\n",cr,i,buf[cr][bp[cr]+i], raw[nw].m, raw[nw].s, raw[nw].v);
	  }else{
	    raw[nw].m  = (buf[cr][bp[cr]+i] >> 27) & 0x1f;
	    raw[nw].s  = (buf[cr][bp[cr]+i] >> 11) & 0x1f;
	    raw[nw].v =   buf[cr][bp[cr]+i]        & 0x7ff;
	    fprintf(stdout,"DRT    buf[%d][%2d]:%8x  m:%2d s:%2d v:%12d\n",cr,i,buf[cr][bp[cr]+i], raw[nw].m, raw[nw].s, raw[nw].v);
	  }
	  nw++;
	  if (nw > BLEN){
	    fprintf(stderr,"$E Size limit exceeded. Skipping.\n");
	    return(2);
	  }       
	}
	raw[nw].m  = (buf[cr][bp[cr]+nword[cr]] >> 27) & 0x1f;
	raw[nw].s  = (buf[cr][bp[cr]+nword[cr]] >> 11) & 0x1f;
	raw[nw].v =   buf[cr][bp[cr]+nword[cr]]        & 0x7ff;
	fprintf(stdout,"DRT    buf[%d][%2d]:%8x  m:%2d s:%2d v:%12d\n",cr,i,buf[cr][bp[cr]+nword[cr]], raw[nw].m, raw[nw].s, raw[nw].v);
	//	printf("%d:%d->%x\n",cr,i, buf[cr][bp[cr]+nword[cr]] );
#endif
      
	module_bit_status = 1;
      }
      //      printf("--> %x\n",buf[cr][bp[cr]+nword[cr]]);

      bp[cr] += nword[cr]+1;
      //      fprintf(stdout,"%d: %d <-> %d\n",cr,bp[cr],bs[cr]);

    } // for(cr....)

    // call user routin

    if ( 0 == module_bit_status){
      status = uana(nw,raw);
      if( status<0 ){
	end_of_analysis=1;
	return(0);
      }
      if (event_number>=ana_event && ana_event>0){
	end_of_analysis=1;
	return(0);
      }
    }
    event_number++; // <--- 本当に外に出しますか?

    // Check end of buffer condition
    // bs[cr]: buffer size
    // bp[cr]: current position
    //    for(cr=0;cr<number_of_smp;cr++)fprintf(stderr,"%d %d %d\n",cr,bp[cr],bs[cr]);
    cond1=0; cond2=0;
    for(cr=0;cr<number_of_smp;cr++){
      if (bp[cr] == bs[cr])cond1++; // Correct number of SDS for each SMP
      if (bp[cr]  > bs[cr])cond2++; // Too many data for this SMP
    }
    if (cond1 == number_of_smp)return(0); // Normal completion
    if (cond1 >= 1){
      fprintf(stderr,"$I rawread: DATA be broken 3. Skipping. %d,%d\n",cond1,cond2);
      return(3);
    }
    if (cond2 >= 1){
      fprintf(stderr,"$I rawread: DATA be broken 2. Skipping.\n");
      return(1);
    }
  } // for(;;)
}

int set_num_modules( int nscach, int nsmp )
{
  number_of_smp = nsmp;
  number_of_scaler_channels = nscach;
}

int set_crate_type( int ncr, int nsl, int type )
{
  // sl 1 origin
  crate_type[ncr][nsl-1]=type;
}

//int main(int argc, char* argv[])
int cmain(int argc, char* argv[]) // change by tsukada '09/11/24
{
  int i,j,k;
  int status;
  int val;
  unsigned int cr;
  static int isds=0;

  fprintf(stderr,"$I rawread: Version %3.1f\n",rawread_version);

  for(cr=0;cr<MAXSMP;cr++){
    buf[cr] = (int*) malloc(sizeof(int)*BLEN);
    if (buf[cr]==NULL){
      fprintf(stderr,"$E Failed allocating memory\n");
      exit(1);
    }
  }

  strcpy(file,argv[4]);

  //check file exit
  FILE *fp;
  if( (fp=fopen(file,"r"))==0 ){
      fprintf(stderr,"$E rawread: fileopen error\n"); 
    return 1;
  }
  fclose(fp);

  int type =0; // NotCompressed=0,Bzip2Com=1,GzipCom=2;
  char comtmp1[64] = "file ";
  char comtmp2[64] =" | grep bzip2 > /dev/null";
  char comtmp3[64] =" | grep gzip > /dev/null";
  char com1[256];  char com2[256];
  sprintf(com1,"%s%s%s",comtmp1,file,comtmp2);
  sprintf(com2,"%s%s%s",comtmp1,file,comtmp3);

  if( system(com1)==0 ) type = 1;//Bzip2Compressed;
  else if( system(com2)==0 ) type =2;// GzipCompressed;


  //  FileCompressedType type = FileType( file );
  //  if( type == NotExist ) return false;
  fprintf(stderr,"$I rawread: %s filetype %d\n",file,type); 
  if( type == 1/*Bzip2Compressed*/ ){
     char com[32] = "bzcat "; strcat(com,file);
     fprintf(stdout,"rawread: bzcat fileopen\n"); 
    if( (fd=popen(com,"r"))==0 ){
      fprintf(stderr,"$E rawread: fileopen error\n"); 
      return 1;
    }

  }
  else if( type ==2/* GzipCompressed */){
     char com[32] = "zcat ";strcat(com,file);
    if( (fd=popen(com,"r"))==0 ){
      fprintf(stderr,"$E rawread: fileopen error\n"); 
      return 1;
    }

  }
  else{
    if( (fd=fopen(file,"r"))==0 ){
      fprintf(stderr,"$E rawread: fileopen error\n"); 
      return 1;
    }
  }


/*   if (fd=fopen(file,"r")==0){ */
/*     fprintf(stderr,"$E rawread: fileopen error\n"); */
/*     exit(1); */
/*   } */

  fprintf(stderr,"$I rawread: Analysis started with%2d SMPs and %2dch scalers\n",
	  number_of_smp, number_of_scaler_channels);

  uinit();

  for(;;){
    if (isds == ana_size && ana_size >0){
      fprintf(stderr,"$I rawread: End of Analysis reached_1\n");
      break;
    }

    // read UNIDAQ header
    /*     if (read(fd, unihead, HLEN*4)<HLEN*4){ */
    /*       fprintf(stderr,"$I rawread: End of File reached\n"); */
    /*       break; */
    /*     } */

    // modified by k.tsukada 2010/9/27
    // if the unidaq header is not 0 (meaning begin,resume,pause,end), the data is skipped.
    // from here
    status=0;
    //    if (read(fd, unihead, 2*4)<2*4){
    if (fread( unihead, 2*4,1,fd)<1){
      fprintf(stderr,"$I rawread: End of File reached_2\n");
      break;
    }
    //fprintf(stdout,"unihead[0]:%x,unihead[1]:%x\n",unihead[0],unihead[1]);
    while( unihead[1]!=0 ){
      for( i=2; i<unihead[0]; i++ ){
	//	if (read(fd, &val, 1*4)<1*4){
	if (fread( &val, 1*4,1,fd)<1){
	  fprintf(stderr,"$I rawread: End of File reached_3_test\n");
	  status=1;
	  break;
	}
      }
      if (fread( unihead, 2*4,1,fd)<1){
	fprintf(stderr,"$I rawread: End of File reached_4\n");
	status=1;
      }
      if( status==1 ) break;
    }
    if( status==1 ) break;
    if( fread( &unihead[2], (HLEN-2)*4,1,fd )<1 ){
      fprintf(stderr,"$I rawread: End of File reached_5\n");
      break;
    }
    // to here

#ifdef DEBUG
    fprintf(stdout,
	    "Num of Words:     %8d\n"
	    "Record Number:    %8d\n"
	    "Run Number:       %8d\n"
	    "Event Set Number: %8d\n"
	    "Mode:             %8d\n",
	    unihead[0],unihead[1],unihead[2],unihead[3],unihead[4]);
#endif

    for(cr=0;cr<number_of_smp;cr++){
      // read TKO header
      if (fread( head[cr], TKO_HEADER_LEN*4,1,fd)<0){
	fprintf(stderr,"$E rawread: Unexpected EOF 0\n");
	ufini();
	exit(3); // eof
      }
#ifdef DEBUG
      fprintf(stdout,
	      "  ## CRATE %d ##\n"
	      "    STAR (status)         : 0x%x\n"
	      "    CMR  (command)        : 0x%x\n"
	      "    SNCR (serial number)  : 0x%x\n"
	      "    DSR  (data size)      : 0x%x\n"
	      "    ENR  (event number)   : 0x%x\n"
	      "    BCR  (buffer control) : 0x%x\n"
	      "    SMP buffer length     : 0x%x\n"
	      "    SMP Crate ID          : %d\n"
	      "    SMP Marker            : 0x%x\n"
	      "    tv sec                : %d\n"
	      "    tv usec               : %d\n"
	      ,cr,head[cr][0],head[cr][1],head[cr][2],head[cr][3],head[cr][4],
	      head[cr][5],head[cr][6],head[cr][7],head[cr][8],head[cr][9],
	      head[cr][10]);
#endif
      subheader[cr].star = head[cr][0];
      subheader[cr].cmr  = head[cr][1];
      subheader[cr].sncr = head[cr][2];
      subheader[cr].dsr  = head[cr][3];
      subheader[cr].enr  = head[cr][4];
      subheader[cr].bcr  = head[cr][5];
      subheader[cr].isize= head[cr][6];
      subheader[cr].cid  = head[cr][7];
      subheader[cr].sec  = head[cr][9];
      subheader[cr].usec = head[cr][10];

      /*for( i=0; i<11; i++ ){ // for debug by k.tsukada 2010/09/27
	fprintf( stdout, " head[%d][%d]=%x\n",cr,i,head[cr][i] );
	}*/
      
      if (head[cr][8]!=0x1234){
	fprintf(stderr,"$E rawread: Marker mismatch: %x check num. of SMP and scaler chs\n",head[cr][8]);
	ufini();
	exit(2);
      }
      if (head[cr][6]>BLEN){
	fprintf(stderr,"$E rawread: Buffer length too long %d > %d\n",
		head[cr][6],BLEN);
	ufini();
	exit(2);
      }
      // read data into buffer
      bs[cr] = head[cr][6];
      if (fread( buf[cr], (head[cr][6])*4,1,fd) < 0){
	fprintf(stderr,"$E rawread: Unexpected EOF 1\n");
	ufini();
	exit(3); // eof
      }

      /*for( i=bs[cr]-5; i<bs[cr]; i++ ){ // for debug by k.tsukada 2010/09/27
	fprintf( stdout, "buf[%d]:%x\n",i,buf[cr][i] );
	}*/
      
    } // loop for cr 
    utime(subheader[0].sec); // added by hashimoto on '11/09/12

    if (fread(scaler, number_of_scaler_channels*4,1,fd)<1){
	fprintf(stderr,"$E rawread: Unexpected EOF 2\n");
	ufini();
	exit(3); // eof
    }
    usca(number_of_scaler_channels, scaler); // added by tsukada '09/11/26
    /*    {
      // Debugging by KI on 7/Feb/2006
      int i;
      for (i=0;i<number_of_scaler_channels;i++)
	fprintf(stderr,"%x\n",scaler[i]);
	}*/

    uopt(0); // Optional Hook

    // After reading one buffer switch, we will compose an SDS event
    compose_an_event();

    uopt(1); // Optional Hook
    
    if (end_of_analysis == 1) break;

    isds++;
  }
  ufini();
  fclose(fd);
  for(cr=0;cr<MAXSMP;cr++){
    free(buf[cr]);
  }
}

