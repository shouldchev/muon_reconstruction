#include "LiRate.h"
#include "libPerso.h"
#include "plot.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#define C 299792458

using namespace std;

void LiRate_raw_2 (){

    //simulation datafile
    ifstream junk ("./raw_2.dat");

    //Getting the number of lines in the simulation datafile
    Int_t numLines0=0;
    if(numLines0==0){string unused1;
        while (getline(junk,unused1)){
                if(!junk){cout<<"error opening file"<<endl;}//if junk
                ++numLines0;
            }//while
    }//if numlines
    junk.close();
    numLines0=(numLines0+1)/2.;
    cout<<"number of generated muons : "<<numLines0<<endl;
    const int numLines=numLines0;

    //simulation datafile
    ifstream file;

    Double_t numInDet=0;
    Int_t evtID=0;
    Double_t mass=0.105658;
    Int_t check=0;

    Double_t R=35.4/2.;
    Double_t Rate=3.5;
    Double_t factor=0.0215;

    //Opening the datafile and creating rootfile and tree
    file.open("./raw_2.dat", ios::in);
    rootfile=new TFile ("./rootfile_raw_2.root","RECREATE");
    myTree=new TTree ("myTree","");

    //Creating branches for the tree
    TBranch *b_evtID=myTree->Branch("evtID",&evtID,"evtID/i");
    TBranch *b_check=myTree->Branch("check",&check,"check/i");
    TBranch *b_posX=myTree->Branch("posX",&Xm,"posX/d");
    TBranch *b_posY=myTree->Branch("posY",&Ym,"posY/d");
    TBranch *b_posZ=myTree->Branch("posZ",&Zm,"posZ/d");
    TBranch *b_momX=myTree->Branch("momX",&Xp,"momX/d");
    TBranch *b_momY=myTree->Branch("momY",&Yp,"momY/d");
    TBranch *b_momZ=myTree->Branch("momZ",&Zp,"momZ/d");
    TBranch *b_E=myTree->Branch("Energy",&E,"Energy/d");
    TBranch *b_Length=myTree->Branch("Length",&Length,"Length[check]/d");
    TBranch *b_LiRate=myTree->Branch("LiRate",&LiRate,"LiRate[check]/d");
    TBranch *b_LiRateInt=myTree->Branch("LiRateInt",&LiRateInt,"LiRateInt/d");
    TBranch *b_numInDet=myTree->Branch("numInDet",&numInDet,"numInDet/d");

    for(int i=0 ; i<numLines ; i++){
        LiRate[0]=0;
        check=0;
        evtID=i;
        file>>junk2>>signe>>junk2>>junk2>>Xp>>Yp>>Zp>>junk2>>junk2>>Xm>>Ym>>Zm>>junk2;
        Xm=Xm/1000.;
        Ym=Ym/1000.;
        Zm=Zm/1000.;
        a=Xp*Xp+Yp*Yp+Zp*Zp;
        b=2.*(Xp*Xm+Yp*Ym+Zp*Zm);
        c=Xm*Xm+Ym*Ym+Zm*Zm-R*R;
        d=b*b-4*a*c;
        E=TMath::Sqrt(mass*mass+(Xp*Xp+Yp*Yp+Zp*Zp));
        if(d>0){
            check=1;
            K[0]=(-b+TMath::Sqrt(d))/(2*a);
            K[1]=(-b-TMath::Sqrt(d))/(2*a);
            numInDet++;
            Length[0]=TMath::Sqrt( Xp*Xp*(K[1]-K[0])*(K[1]-K[0]) +  Yp*Yp*(K[1]-K[0])*(K[1]-K[0]) +Zp*Zp*(K[1]-K[0])*(K[1]-K[0]) );
            LiRate[0]=factor*pow(E,0.74)*Length[0]*Rate;            
        }//if d
        LiRateInt+=LiRate[0];
        if(i==numLines-1){LiRateInt=LiRateInt/numInDet;}
        myTree->Fill();
    }// for i
    cout<<"Number of muons passing in detector : "<<numInDet<<endl;
    cout<<"Lithium rate : "<<LiRateInt<<endl;
    cout<<"Time : "<<numInDet/Rate<<endl;

    file.close();
    myTree->Write();
    rootfile->Close();

    //plot("./rootfile_raw_2.root",0);

    return;

}//main

