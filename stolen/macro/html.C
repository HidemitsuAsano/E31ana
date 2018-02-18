void html()
{
  gROOT->Reset();

  gSystem->Load("./lib/libAll.so");

  THtml html;

  Char_t *cmd = "rm -rf  ./htmldoc/*";
  gSystem->Exec( cmd );

  html.SetSourceDir( ".:src:src" );
  html.SetInputDir(  "./src" );
  //html.SetInputDir(  "./src:/w/e15/common/root/include" );
  html.SetRootURL( "http://root.cern.ch/root/html" );
  html.SetViewCVS( "http://ag.riken.jp/viewvc/viewvc.cgi/k18ana/src/" );
  //html.SetViewCVS( "http://ag.riken.jp/viewvc/viewvc.cgi" );
  html.SetCharset( "EUC-JP" );
  html.SetProductName( "K18BR analyzer" );

  //html.SetOutputDir( "./html/" );
  //html.SetInputDir( ".:src/:src/" );
  //html.SetSourceDir( "./src/" );
  //html.SetInputPath("./src/");
  //html.SetIncludePath("./src/" );

  // main
  html.MakeClass( "GlobalVariables" );
  html.MakeClass( "HodoscopeLikeHit" );
  html.MakeClass( "CherenkovLikeHit" );
  html.MakeClass( "ChamberLikeHit" );
  html.MakeClass( "CDCHit" );
  html.MakeClass( "BeamLineHitMan" );
  html.MakeClass( "CDSHitMan" );
  html.MakeClass( "EventHeader" );
  html.MakeClass( "Scaler" );
  html.MakeClass( "ScalerMan" );
  html.MakeClass( "TKOHit" );
  html.MakeClass( "TKOHitCollection" );
  html.MakeClass( "HelixFit" );
  html.MakeClass( "CircleFit" );
  html.MakeClass( "HitCluster" );
  html.MakeClass( "CDSTrack" );
  html.MakeClass( "CDSTrackingMan" );
  html.MakeClass( "SDDHit" );
  html.MakeClass( "SDDHitMan" );

  // parameters
  html.MakeClass( "ConfMan" );
  html.MakeClass( "CounterMapMan" );
  html.MakeClass( "GainMap" );
  html.MakeClass( "GainMapMan" );
  html.MakeClass( "GeomMap" );
  html.MakeClass( "GeomMapMan" );
  html.MakeClass( "ScalerMapMan" );
  html.MakeClass( "SlewingMap" );
  html.MakeClass( "SlewingMapMan" );
  html.MakeClass( "CDCWireMap" );
  html.MakeClass( "CDCWireMapMan" );
  html.MakeClass( "BLDCWireMap" );
  html.MakeClass( "BLDCWireMapMan" );
  html.MakeClass( "XTMap" );
  html.MakeClass( "XTMapMan" );
  html.MakeClass( "ReslMap" );
  html.MakeClass( "ReslMapMan" );

  // others
  html.MakeClass( "SimDataMan" );
  html.MakeClass( "MathTools" );
  html.MakeClass( "Display" );
  html.MakeClass( "Display3D" );
  html.MakeClass( "Shape" );
  html.MakeClass( "Sphere" );
  html.MakeClass( "Box" );
  html.MakeClass( "Tube" );
  html.MakeClass( "Line" );
  html.MakeClass( "Mark" );

  html.MakeIndex();

  //Char_t *cmd = "cd ./html; mv USER__Index.html K1.8BRClass.html";
  //Char_t *cmd = "cd ./html; mv MAIN_Index.html K1.8BRClass.html";
  //gSystem->Exec( cmd );
}
