//Missing mass



Double_t func_MMnmiss_mod(Double_t *x,Double_t *par)
{
   //connection using woods-saxon
   if(0.0<x[0] && x[0]<1.5){
     return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0))*(1./(1.0+exp((x[0]-1.075)/par[10])))
       +(1.0-1./(1.0+exp((x[0]-1.075)/par[10])))*(par[11]*exp(-0.5*pow(((x[0]-par[12])/par[13]),2.0)))*(1./(1.0+exp((x[0]-1.14)/par[14])))
       +(1.0-1./(1.0+exp((x[0]-1.14)/par[14])))*
        (par[15]+par[16]*x[0]+par[17]*pow(x[0],2.0)+par[18]*pow(x[0],3.0)+par[19]*pow(x[0],4.0));
   }else{
     return 1.0;
   }
}



Double_t param_MMnmiss_mod[20]={
-0.000828238,
1.05546,
-0.141575,
-3.80951,
-2.82955,
159.784,
-616.575,
1010.3,
-768.422,
222.31,
0.018,
2.25934,
1.11012,
0.0390585,
0.01,
300.939,
-865.414,
947.955,
-465.283,
85.8418
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
    double ret=  (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))
                *(1./(1.0+exp((x[0]-1.25)/par[7])))
                //+(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))
                +(par[4]*exp(-0.5*pow(((x[0]-par[5])/par[6]),2.0)))
                *(1.0-1./(1.0+exp((x[0]-1.25)/par[7])))
                *(1./(1.0+exp((x[0]-1.52)/par[8])))
                +(1.0-1./(1.0+exp((x[0]-1.52)/par[8])))
                *(par[9]*exp(-0.5*pow(((x[0]-par[10])/par[11]),2.0)))
                ;
    if(ret<0) return 0;
    else return ret;
  }else{
    return 1.;
  }
}

Double_t param_IMnpip_mod[12]={
-334.885,
756.849,
-559.954,
135.605,
1.22469,
1.23201,
0.166981,
0.01,
0.01,
0.686025,
1.32384,
0.190184
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
  //pol3(1.08-1.14)+pol4(1.14-1.25)+pol5(1.25-1.93)
  if(1.08<x[0] && x[0]<1.63){
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))
    *(1./(1.0+exp((x[0]-1.13)/par[4])))
    +(1.0-1./(1.0+exp((x[0]-1.13)/par[4])))
    *(par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)+par[9]*pow(x[0],4.0))
    *(1./(1.0+exp((x[0]-1.25)/par[10])))
    +(1.0-1./(1.0+exp((x[0]-1.25)/par[10])))
    *(par[11]+par[12]*x[0]+par[13]*pow(x[0],2.0)+par[14]*pow(x[0],3.0)+par[15]*pow(x[0],4.0));
  
    if(ret<0) return 0;
    else return ret;
  }else{
    return 0.0;
  }
};

Double_t param_IMnpip_wK0_mod[16]={
-88509,
242125,
-220810,
67133.3,
0.01,
-3776.85,
6211.08,
27.1978,
-4245.76,
1746.82,
0.01,
32.9473,
-35.4828,
-0.449256,
10.4611,
-2.73038
};




Double_t func_IMnpim_mod(Double_t *x,Double_t *par)
{
  if(x[0]<1.07){
    return 0.0;
  }else if(1.07 <= x[0]){
    return (par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0)))
          //*(1./(1.0+exp((x[0]-1.10)/par[3])))
          *(1./(1.0+exp((x[0]-1.140)/par[3])))
          //+(1.0 - 1./(1.0+exp((x[0]-1.10)/par[3])))
          +(1.0 - 1./(1.0+exp((x[0]-1.140)/par[3])))
          //*(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)
          // +par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0)+par[11]*pow(x[0],7.0));
          *exp(par[4]+par[5]*x[0]);
 }else{
   return 0;
  }
}


Double_t param_IMnpim_mod[6]={
1.93283,
1.11567,
0.0295707,
0.005,
8.47662,
-7.04357
};


Double_t func_IMnpim_wK0_mod(Double_t *x,Double_t *par)
{
  if(1.00 <= x[0] && x[0]<2.0) {
   //1.07-1.11  gaus    1.11-1.24 pol5            1.24-2.0 exp
    return (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
         *(1./(1.0+exp((x[0]-1.11)/par[3])))
         +(1.0 - (1./(1.0+exp((x[0]-1.11)/par[3]))))
         //*(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0))
         *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0))
         *(1./(1.0+exp((x[0]-1.24)/par[8])))
         +(1.0 - (1./(1.0+exp((x[0]-1.24)/par[8]))))
         *exp(par[9]+par[10]*x[0]);
  } else {
    return 0.;
  }
}


Double_t param_IMnpim_wK0_mod[11]={
3.80053,
1.10253,
0.00869531,
0.02,
4266.42,
-10484.8,
8583.89,
-2340.8,
0.01,
8.57184,
-7.00872
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
122.936,
-1815.6,
9404.75,
-16663.4,
0.01,
10.8499,
-61.4766,
123.645,
-59.6072,
-121.095,
170.816,
-62.9352
};


Double_t param_nmom_wK0_mod[12]={
2674.21,
-45316.8,
259551,
-496799,
0.01,
107.933,
-617.738,
1281.31,
-646.596,
-1456.37,
2235.93,
-909.037
};




Double_t func_IMpippim_mod(Double_t *x,Double_t *par)
{
  if(0.27<x[0]){
  return (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
        *(1./(1.0+exp((x[0]-0.30)/par[3])))            
      +(1.0 - 1./(1.0+exp((x[0]-0.30)/par[3])))
       *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))
       *(1./(1.0+exp((x[0]-0.73)/par[11])))
      +(1.0 - 1./(1.0+exp((x[0]-0.73)/par[11])))
       *(par[12]*exp(-0.5*pow(((x[0]-par[13])/par[14]),2.0)))
       ;
   }else{
    return 0.0;
   }
};


Double_t param_IMpippim_mod[15]={
0.55615,
0.297704,
0.0146312,
0.015,
-0.313228,
4.23276,
-5.95692,
-5.32242,
41.5461,
61.3149,
-153.941,
0.03,
3.00187,
0.697503,
0.0787875
};




Double_t func_q_mod(Double_t *x,Double_t *par)
{
  if(0<x[0]) {  //&& x[0]<1.17){
    return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0));
  }else{
   return 0;
  }
}

Double_t param_q_mod[8]={
4.09016,
1.3497,
-50.8339,
210.801,
-416.419,
422.52,
-206.684,
37.9385
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
0.190352,
6.57234,
-33.3275,
63.9037,
-40.3182,
-26.5057,
46.8204,
-16.7796
};


