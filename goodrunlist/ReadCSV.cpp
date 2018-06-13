//H. Asano
//This is a code to read a CSV file "RunSummary_run78.csv" 
//and make a good run list for RUN78

#include <iostream>
#include <fstream>
//#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;


//split a line to words using delimiter
vector <string> split(const string& src, const char* delimiter = ","){
  vector <string> wordvec;
  string::size_type len = src.length();

  for(string::size_type i=0,n; i< len;i=n+1){
    n =src.find_first_of(delimiter,i);
    if( n== string::npos){
      n = len;
    }
    wordvec.push_back(src.substr(i,n-i));
  }

  return wordvec;
}

int main(){
  
  ifstream ifs("RunSummary_run78.csv");
  std::string str;
  
  const int maxline=1024;
  int runnum[maxline];
  std::string duration[maxline];
  int duration_sec[maxline];
  std::string comment[maxline];
  
  int iline=0;
  while(getline(ifs,str)){
    vector<string> arr = split(str);
    for(int i=0;i<arr.size();++i){
      //cout << arr[i] << "\t";
    }
    std::cout << std::endl;
    runnum[iline] = atoi(arr[0].c_str());
    duration[iline] = arr[7];
    //std::cout << arr.size() << std::endl;
    if(arr.size()==70) comment[iline] = arr[37];
    else if(arr.size()== 71) comment[iline] = arr[38];
    else if(arr.size()== 72) comment[iline] = arr[39];
    else if(arr.size()== 79) comment[iline] = arr[46];
    else if(arr.size()== 80) comment[iline] = arr[47];
    else if(arr.size()== 81) comment[iline] = arr[48];
    else if(arr.size()== 82) comment[iline] = arr[49];
    else { 
      std::cout << "invalid format !! " << arr.size() << std::endl;
    }
    if(duration[iline]=="") duration_sec[iline]=0;
    else{
      vector<string> stime = split(duration[iline],":");
      if(stime.size()==1){
        duration_sec[iline] = atoi(stime[0].c_str());
      }else if(stime.size()==2){
        duration_sec[iline] = atoi(stime[1].c_str())+ 60.*atoi(stime[0].c_str());
      }else if(stime.size()==3){
        duration_sec[iline] = atoi(stime[2].c_str())+ 60.*atoi(stime[1].c_str())
        + 3600.*atoi(stime[0].c_str());
      }
    }
    cout << runnum[iline] << "\t" << duration[iline] << "\t" <<  duration_sec[iline] << "\t" << comment[iline];
    
    iline++;
  }

  
  ofstream ofs_pro("production_list");
  ofstream ofs_emp("empty_list");
  ofstream ofs_d5("D5scan_list");
  ofstream ofs_uswk("USWKscan_list");
  ofstream ofs_short("shortrun_list");
  int pro_sumtime=0;
  int pro_nruns=0;
  int emp_sumtime=0;
  int d5_sumtime=0;
  int uswk_sumtime=0;
  int junk_sumtime=0;

  //generate shell script for hadd
  ofstream ofs_pro_haddcsh("hadd.csh");
  ofs_pro_haddcsh << "#!/bin/tcsh -f" << endl;

  for(int i=0; i<iline;i++){
    
    int pos_pro = comment[i].find("production run");
    int pos_emp = comment[i].find("empty");
    int pos_d5 = comment[i].find("D5");
    int pos_uswk = comment[i].find("USWK");
    //good run & more than 20 minutus
    if(pos_pro!=-1){
     if(duration_sec[i]>20*60){
       if(pro_nruns%50==0){
         ofs_pro_haddcsh << endl;
         ofs_pro_haddcsh << "hadd -To " ;
         char sumname[256];
         sprintf(sumname,"evanaRead_M%d.root",(int)pro_nruns/50);
         ofs_pro_haddcsh << sumname << " "  ;
       }
       //cout << runnum[i] <<  duration_sec[i] << "\t" << comment[i];
       ofs_pro << runnum[i] << endl; 
       pro_sumtime+=duration_sec[i];
       char rootname[256];
       sprintf(rootname,"evanaReadAna_0%03d.root",runnum[i]);
       //cout << rootname << endl;
       ofs_pro_haddcsh << rootname << " "; 
       pro_nruns++;
     }else{
       ofs_short << runnum[i] << endl; 
       junk_sumtime +=duration_sec[i];
     }
    }else if(pos_emp!=-1){
      ofs_emp << runnum[i] << endl;
      emp_sumtime+=duration_sec[i];
    }else if(pos_d5!=-1){
      ofs_d5 << runnum[i] << endl;
      d5_sumtime+=duration_sec[i];
    }else if(pos_uswk!=-1){
      ofs_uswk << runnum[i] << endl;
      uswk_sumtime+=duration_sec[i] ;
    }else{
      cout << runnum[i] << "\t" << duration_sec[i] << "\t" << comment[i] << endl;
      junk_sumtime+=duration_sec[i]; 
    }
  }

  cout << endl;
  cout << "production run: " << pro_sumtime << endl;
  cout << "empty run " << emp_sumtime << endl;
  cout << "D5 Scan: " << d5_sumtime << endl;
  cout << "USWK Scan: " << uswk_sumtime << endl;
  cout << "discard: " << junk_sumtime << endl;

}


