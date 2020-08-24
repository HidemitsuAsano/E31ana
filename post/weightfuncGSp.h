//Missing mass



Double_t func_MMnmiss_mod(Double_t *x,Double_t *par)
{
   //connection using woods-saxon
   if(0.0<x[0] && x[0]<1.5){
     double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0))*(1./(1.0+exp((x[0]-1.075)/par[10])))
       +(1.0-1./(1.0+exp((x[0]-1.075)/par[10])))*(par[11]*exp(-0.5*pow(((x[0]-par[12])/par[13]),2.0)))*(1./(1.0+exp((x[0]-1.14)/par[14])))
       +(1.0-1./(1.0+exp((x[0]-1.14)/par[14])))*
        (par[15]+par[16]*x[0]+par[17]*pow(x[0],2.0)+par[18]*pow(x[0],3.0)+par[19]*pow(x[0],4.0));
     if(ret<0) return 0;
     else      return ret;
   }else{
     return 1.0;
   }
}



Double_t param_MMnmiss_mod[20]={
-0.00176863,
3.01634,
-4.15562,
8.69902,
-14.1173,
147.881,
-611.02,
1027.05,
-759.249,
204.731,
0.018,
11.7534,
0.632306,
0.257789,
0.01,
303.948,
-866.895,
946.229,
-465.747,
86.7758
};



Double_t func_MMnmiss_wK0_mod(Double_t *x,Double_t *par)
{
  double ret= (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0))
       *(1./(1.0+exp((x[0]-1.08)/par[6])))
       +(1. - 1./(1.0+exp((x[0]-1.08)/par[6])))
       *(par[7]+par[8]*x[0]+par[9]*pow(x[0],3.0)+par[10]*pow(x[0],4.0))
       ;
   if(ret<0) return 0;
   else return ret;
}

Double_t param_MMnmiss_wK0_mod[11]={
-0.00714053,
0.991907,
4.4545,
-20.0309,
35.2696,
-19.0939,
0.01,
-7.59543,
8.25753,
4.54261,
-4.65572
};


Double_t func_IMnpip_mod(Double_t *x,Double_t *par)
{
  if(1.06<=x[0] && x[0]<2.00){
    //gaus (1.06-1.23) 
    double ret= // (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))
                (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
                *(1./(1.0+exp((x[0]-1.23)/par[3])))
                +(1.0-1./(1.0+exp((x[0]-1.23)/par[3])))
                *(exp(par[4]+par[5]*x[0]));
    if(ret<0) return 0;
    else return ret;
  }else{
    return 1.;
  }
}

Double_t param_IMnpip_mod[6]={
2.17496,
1.16759,
0.0445005,
0.01,
14.4465,
-11.5632
};


Double_t func_Mompippim(Double_t *x,Double_t *par)
{
  return 
    par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)
    +par[5]*pow(x[0],5.0)    
    +par[6]*pow(x[0],6.0);    
}

Double_t param_Mompippim[7]={
0.826769,
-5.3402,
54.5078,
-149.974,
144.084,
-31.3555,
-12.6739
};



Double_t func_IMnpip_wK0_mod(Double_t *x,Double_t *par)
{
  if(1.08<x[0] && x[0]<1.63){
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))
    *(1./(1.0+exp((x[0]-1.19)/par[4])))
    +(1.0-1./(1.0+exp((x[0]-1.19)/par[4])))
    //*(par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)+par[9]*pow(x[0],4.0))
    *(par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0))
    //*(1./(1.0+exp((x[0]-1.25)/par[10])))
    *(1./(1.0+exp((x[0]-1.25)/par[8])))
    //+(1.0-1./(1.0+exp((x[0]-1.25)/par[10])))
    +(1.0-1./(1.0+exp((x[0]-1.25)/par[8])))
    //*(par[11]+par[12]*x[0]+par[13]*pow(x[0],2.0)+par[14]*pow(x[0],3.0)+par[15]*pow(x[0],4.0));
    *(par[9]+par[10]*x[0]+par[11]*pow(x[0],2.0)+par[12]*pow(x[0],3.0)+par[13]*pow(x[0],4.0));
  
    if(ret<0) return 0;
    else return ret;
  }else{
    return 0.0;
  }
};

Double_t param_IMnpip_wK0_mod[14]={
-232.137,
179.349,
244.028,
-194.096,
0.01,
-14.2499,
24.9294,
-9.69115,
0.01,
37.9034,
-47.0971,
14.3774,
0.154323,
-0.00105421,
};




Double_t func_IMnpim_mod(Double_t *x,Double_t *par)
{
  if(x[0]<1.07){
    return 0.0;
  }else if(1.07 <= x[0]){
    double ret = (par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0)))
          //*(1./(1.0+exp((x[0]-1.10)/par[3])))
          *(1./(1.0+exp((x[0]-1.12)/par[3])))
          //+(1.0 - 1./(1.0+exp((x[0]-1.10)/par[3])))
          +(1.0 - 1./(1.0+exp((x[0]-1.12)/par[3])))
          *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)
           +par[9]*pow(x[0],5.0)); //+par[10]*pow(x[0],6.0)+par[11]*pow(x[0],7.0));
          //*exp(par[4]+par[5]*x[0]);
    if(ret>0) return ret;
    else      return 0;
 }else{
   return 0;
  }
}


Double_t param_IMnpim_mod[10]={
1.08686,
1.10987,
0.0169923,
0.015,
1433.46,
-5025.37,
6969.36,
-4775.93,
1617.75,
-216.861
};


Double_t func_IMnpim_wK0_mod(Double_t *x,Double_t *par)
{
  if(1.00 <= x[0] && x[0]<2.0) {
   //1.07-1.11  gaus    1.11-1.24 pol5            1.24-2.0 exp
    return (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
         *(1./(1.0+exp((x[0]-1.092)/par[3])))
         +(1.0 - (1./(1.0+exp((x[0]-1.092)/par[3]))))
         //*(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0))
         //*(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0))
         *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0))
         //*(1./(1.0+exp((x[0]-1.22)/par[8])))
         *(1./(1.0+exp((x[0]-1.22)/par[7])))
         //+(1.0 - (1./(1.0+exp((x[0]-1.22)/par[8]))))
         +(1.0 - (1./(1.0+exp((x[0]-1.22)/par[7]))))
         //*exp(par[9]+par[10]*x[0]);
         *exp(par[8]+par[9]*x[0]);
  } else {
    return 0.;
  }
}


Double_t param_IMnpim_wK0_mod[10]={
0.0858578,
1.10463,
0.00904288,
0.005,
424.378,
-705.599,
293.465,
0.005,
8.48006,
-6.98101
};




Double_t func_nmom_mod(Double_t *x,Double_t *par)
{
  if(x[0]<0.140){
    return 1.0;
  }else if(0.14<= x[0] && x[0]<1.0){
 //   return exp(par[0]+par[1]*x[0])*((1./(1.0+exp((x[0]-0.2)/par[2]))))
 //   + (par[3]+par[4]*x[0]+par[5]*pow(x[0],2.0)+par[6]*pow(x[0],3.0)+par[7]*pow(x[0],4.0)+par[8]*pow(x[0],5.0)
 //   +par[9]*pow(x[0],6.0))*(1.-1./(1.0+exp((x[0]-0.2)/par[2])));
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))*((1./(1.0+exp((x[0]-0.2)/par[4]))))
    + (par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)+par[9]*pow(x[0],4.0)+par[10]*pow(x[0],5.0)
    +par[11]*pow(x[0],6.0))*(1.-1./(1.0+exp((x[0]-0.2)/par[4])));
    if(ret<0) return 0;
    else return ret;
  }else{
    return 1.0;
  }
}

Double_t param_nmom_mod[12]={
128.949,
-1835.61,
9281.78,
-16181.3,
0.01,
20.1085,
-160.068,
544.765,
-989.2,
1001.7,
-533.68,
116.749
};


Double_t param_nmom_wK0_mod[12]={
2946.78,
-50380.9,
290158,
-557585,
0.01,
106.424,
-612.966,
1283.35,
-651.591,
-1466.61,
2230.77,
-891.24
};




Double_t func_IMpippim_mod(Double_t *x,Double_t *par)
{
  if(0.27<x[0]){
  double ret = (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
        *(1./(1.0+exp((x[0]-0.30)/par[3])))            
      +(1.0 - 1./(1.0+exp((x[0]-0.30)/par[3])))
       *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))
       *(1./(1.0+exp((x[0]-0.73)/par[11])))
      +(1.0 - 1./(1.0+exp((x[0]-0.73)/par[11])))
       *(par[12]*exp(-0.5*pow(((x[0]-par[13])/par[14]),2.0)))
       ;
  if(ret<0) return 0;
  else      return ret;
  }else{
    return 0.0;
  }
};


Double_t param_IMpippim_mod[15]={
0.360188,
0.296949,
0.00982654,
0.015,
6.39021,
-32.8895,
35.8731,
88.5842,
-70.0974,
-283.483,
296.511,
0.03,
0.562649,
0.786209,
0.0525424
};




Double_t func_q_mod(Double_t *x,Double_t *par)
{
  if(0<x[0]) {  //&& x[0]<1.17){
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0));
    if(ret>0) return ret;
    else      return 0;
  }else{
   return 0;
  }
}

Double_t param_q_mod[8]={
3.0229,
-2.85293,
3.70457,
-2.08965,
-2.97093,
1.50674,
2.97921,
-1.88548
};

Double_t func_q_wK0_mod(Double_t *x,Double_t *par)
{
  if(0<x[0] && x[0]<1.3) {  //&& x[0]<1.17){
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0));
    
    
    if(ret <0 ) return 0;
    else return ret;
  }else{
    return 0;
  }
}



Double_t param_q_wK0_mod[8]={
0.287913,
6.54975,
-34.3679,
64.8546,
-38.1852,
-29.9088,
47.7638,
-16.515
};


