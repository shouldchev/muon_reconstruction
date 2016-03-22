#include "LiRate.h"
#include "libPerso.h"
#include "plot.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#define C 299792458

using namespace std;

void LiRate_raw (){

    ifstream junk ("./raw.dat");
    Int_t numLines0=0;

    if(numLines0==0){string unused1;
        while (getline(junk,unused1)){
                if(!junk){cout<<"error opening file"<<endl;}
                ++numLines0;
            }
    }    
    junk.close();
    numLines0=(numLines0+1)/2.;
    cout<<"number of generated muons : "<<numLines0<<endl;
    const int numLines=numLines0;

    ifstream file;

    Int_t numInDet=0;
    Int_t TrackID=0;
    Double_t mass=0.105658;

    Double_t R=35.4/2.;
    Double_t Rate=3.5;
    Double_t factor=0.0215;
    Double_t LiRate=0;

    file.open("./raw.dat", ios::in);
    rootfile=new TFile ("./rootfile_raw.root","RECREATE");
    myTree=new TTree ("myTree","");
    TBranch *b_TrackID=myTree->Branch("TrackID",&TrackID,"TrackID/i");
    TBranch *b_posX=myTree->Branch("posX",&Xm,"posX/d");
    TBranch *b_posY=myTree->Branch("posY",&Ym,"posY/d");
    TBranch *b_posZ=myTree->Branch("posZ",&Zm,"posZ/d");
    TBranch *b_momX=myTree->Branch("momX",&Xp,"momX/d");
    TBranch *b_momY=myTree->Branch("momY",&Yp,"momY/d");
    TBranch *b_momZ=myTree->Branch("momZ",&Zp,"momZ/d");
    TBranch *b_E=myTree->Branch("Energy",&E,"Energy/d");
    TBranch *b_Length=myTree->Branch("Length",&Length,"Length/d");
    TBranch *b_LiRate=myTree->Branch("LiRate",&LiRate,"LiRate/d");

    for(int i=0 ; i<numLines ; i++){
        TrackID=i;
        file>>junk2>>signe>>junk2>>junk2>>Xp>>Yp>>Zp>>junk2>>junk2>>Xm>>Ym>>Zm>>junk2;
        Xm=Xm/1000.;
        Ym=Ym/1000.;
        Zm=Zm/1000.;
        a=Xp*Xp+Yp*Yp+Zp*Zp;
        b=2.*(Xp*Xm+Yp*Ym+Zp*Zm);
        c=Xm*Xm+Ym*Ym+Zm*Zm-R*R;
        d=b*b-4*a*c;
        if(d>0){
            K[0]=(-b+TMath::Sqrt(d))/(2*a);
            K[1]=(-b-TMath::Sqrt(d))/(2*a);
            numInDet++;
            Length=TMath::Sqrt( Xp*Xp*(K[1]-K[0])*(K[1]-K[0]) +  Yp*Yp*(K[1]-K[0])*(K[1]-K[0]) +Zp*Zp*(K[1]-K[0])*(K[1]-K[0]) );                        
            E=TMath::Sqrt(mass*mass+(Xp*Xp+Yp*Yp+Zp*Zp));
            LiRate=factor*pow(E,0.74)*Length*Rate;
            myTree->Fill();
        }
    }
    cout<<"Number of muons passing in detector : "<<numInDet<<endl;
    file.close();
    myTree->Write();
    rootfile->Close();

    plot("./rootfile_raw.root");

    return;

}

