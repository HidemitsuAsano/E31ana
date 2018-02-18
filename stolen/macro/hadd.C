#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );

void hadd() {
  // in an interactive ROOT session, edit the file names
  // Target and FileList, then
  // root > .L hadd.C
  // root > hadd()
  
  Target = TFile::Open( "root_out/all.root", "RECREATE" );
  FileList = new TList();

  int run[] = {3057,3058,3059,
               3060,3061,     3063,3064,3065,3066,3067,3068,3069,
	       3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,
	       3080,3081,3082,    ,3084,3085,3086,3087,3088,3089.
	       3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,
	       3100,3101,3102,3103,3104,3105,3106};
  int nrun = sizeof(run)/sizeof(int);
  for( int i=0; i<nrun; i++ ){
    FileList->Add(  TFile::Open( Form("root_out/trout_new_%d.root",run[i])) );
  }

  MergeRootfile( Target, FileList );

}   

void MergeRootfile( TDirectory *target, TList *sourcelist ) {
  
  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );
  
  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;
  cout << "Source path: " << current_sourcedir->GetPath() << endl;
  
  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;
  while ( (key = (TKey*)nextkey())) {
    
    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();
    
    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> merge it
      
      cout << "Merging histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;
      
      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {
	
        // make sure we are at the correct directory level by cd'ing to path
        nextsource->cd( path );
        TH1 *h2 = (TH1*)gDirectory->Get( h1->GetName() );
        if ( h2 ) {
          h1->Add( h2 );
          delete h2; // don't know if this is necessary, i.e. if 
                     // h2 is created by the call to gDirectory above.
        }
	
        nextsource = (TFile*)sourcelist->After( nextsource );
      }
    }
    else if ( obj->IsA()->InheritsFrom( "TTree" ) ) {
      
      // loop over all source files create a chain of Trees "globChain"
      const char* obj_name= obj->GetName();
      
      globChain = new TChain(obj_name);
      globChain->Add(first_source->GetName());
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
//       const char* file_name = nextsource->GetName();
//       cout << "file name  " << file_name << endl;
      while ( nextsource ) {
	
	globChain->Add(nextsource->GetName());
	nextsource = (TFile*)sourcelist->After( nextsource );
      }
      
    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory
      
      cout << "Found subdirectory " << obj->GetName() << endl;
      
      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
      
      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      MergeRootfile( newdir, sourcelist );
      
    } else {
      
      // object is of no type that we know or can handle
      cout << "Unknown object type, name: " 
           << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }
    
    // now write the merged histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();
      
      //!!if the object is a tree, it is stored in globChain...
      if(obj->IsA()->InheritsFrom( "TTree" ))
	globChain->Write( key->GetName() );
      else
	obj->Write( key->GetName() );
    }
    
  } // while ( ( TKey *key = (TKey*)nextkey() ) )
  
    // save modifications to target file
  target->Write();

}
  
