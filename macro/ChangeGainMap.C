//H. Asano
//a macro to change GainMap 

void ChangeGainMap(){

  gSystem->Load("lib/libAll.so");
  
  //FDC
  const int crfdc=9;
  const int slfdc[12]={9,12,13,14,16,17,18,19,20,21,22,23};
  const int cr[7]   ={ 1, 0, 0, 1, 3, 4, 5};
  const int slini[7]={18, 1, 7, 1, 1, 1, 1};
  const int slend[7]={22, 4,22,16,20,20,20};

  ConfMan *conf = new ConfMan("conf/Run78/analyzer.conf",27);
  conf->Initialize();

  
  
  ofstream ofsbl("GainMapBL_0001_0028.param");
  ofstream ofscds("GainMapCDS_0001_0164.param");
  
  //FDC slot
  for(int isl=0;isl<12;isl++){
    double offs,temp1;
    conf->GetGainMapManager()->GetParam( cr,sl,ch, temp1, offs );
    conf->GetGainMapManager()->SetParam( cr,sl,ch, offs, p1, p2 ); // set a new parameter into map
  }
}
