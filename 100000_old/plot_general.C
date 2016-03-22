#include "libPerso.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void plot_general (){
   gSystem->Exec("find -name '*.root' -print > data.dat");
   ifstream datafile ("./data.dat");
   TObjArray RootFile(0);
   int numLines=0;
   TFile *rootfile;

   string str;
   TString TStr;
   TFile *myfiles;

   while(getline(datafile,str)){
       if(!datafile){cout<<"No root file found"<<endl;}
       TStr=str;
       myfiles=new TFile( TStr.Data(),"RECREATE" );
       RootFile.Add(myfiles);
       numLines++;
   }
   datafile.close();
   cout<<numLines<<" root file(s) found"<<endl;
   gSystem->Exec("rm data.dat");

   rootfile=(TFile)RootFile[0].At(0);

}
