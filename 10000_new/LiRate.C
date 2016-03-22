#include "LiRate.h"
#include "libPerso.h"
#include "plot.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#define C 299792458

using namespace std;

void LiRate (){

    //simulation datafile
    ifstream junk ("./muon_data_Mars2016.dat");

    //Getting the number of lines in the simulation datafile
    Int_t numLines0=0;
    if(numLines0==0){string unused1; int kt[300000]; int k2=0; int k=0;
        while (getline(junk,unused1)){
            kt[numLines0]=k;
            k2=k/2;
            if(!junk.is_open()){cout<<"Error opening file"<<endl; return;}
            if(k==2*k2 && unused1=="1"){cout<<"Problem in file at line "<<k<<endl; return;}
            ++numLines0;
            k++;
       }//while
    }//if numlines
    junk.close();
    numLines0=(numLines0+1)/2.;
    cout<<"number of generated muons : "<<numLines0<<endl;
    const int numLines=numLines0;

    //simulation datafile
    ifstream file;

    Double_t numInDet_O=0;
    Double_t numInDet_WP=0;
    Int_t evtID=0;
    Double_t mass=0.105658;

    Double_t R_O=35.4/2.;
    Double_t H=42.5;
    Double_t R_WP=H/2.;
    Double_t Rate=13;
    Double_t factor=0.0215;

    //Opening the datafile and creating rootfile and tree
    file.open("./muon_data_Mars2016.dat", ios::in);
    if(!file.is_open()){cout<<"Error opening file"<<endl; return;}
    rootfile=new TFile ("./rootfile_raw_2.root","RECREATE");
    myTree=new TTree ("myTree","");

    //Creating branches for the tree
    TBranch *b_evtID=myTree->Branch("evtID",&evtID,"evtID/i");
    TBranch *b_check_O=myTree->Branch("check_O",&check_O,"check_O/i");
    TBranch *b_check_WP=myTree->Branch("check_WP",&check_WP,"check_WP/i");
    TBranch *b_posX=myTree->Branch("posX",&Xm,"posX/d");
    TBranch *b_posY=myTree->Branch("posY",&Ym,"posY/d");
    TBranch *b_posZ=myTree->Branch("posZ",&Zm,"posZ/d");
    TBranch *b_momX=myTree->Branch("momX",&Xp,"momX/d");
    TBranch *b_momY=myTree->Branch("momY",&Yp,"momY/d");
    TBranch *b_momZ=myTree->Branch("momZ",&Zp,"momZ/d");
    TBranch *b_E=myTree->Branch("Energy",&E,"Energy/d");
    TBranch *b_Length_O=myTree->Branch("Length_O",&Length_O,"Length_O/d");
    TBranch *b_Length_WP=myTree->Branch("Length_WP",&Length_WP,"Length_WP/d");
    TBranch *b_LiRate=myTree->Branch("LiRate",&LiRate,"LiRate/d");
    TBranch *b_LiRateInt=myTree->Branch("LiRateInt",&LiRateInt,"LiRateInt/d");
    TBranch *b_numInDet_O=myTree->Branch("numInDet_O",&numInDet_O,"numInDet_O/d");
    TBranch *b_numInDet_WP=myTree->Branch("numInDet_WP",&numInDet_WP,"numInDet_WP/d");

    for(int i=0 ; i<numLines ; i++){

        //Reinitialisation of variables
        LiRate[0]=0;
        Length_O=0;
        Length_WP=0;
        check_WP=0;
        check_O=0;
        a=0;
        b=0;
        c=0;
        d=0;
        E=0;
        Xm=0;
        Ym=0;
        Zm=0;
        Xp=0;
        Yp=0;
        Zp=0;
        evtID=i;

        //filling variables
        if(i!=numLines-1){file>>junk2>>signe>>junk2>>junk2>>Xp>>Yp>>Zp>>junk2>>junk2>>Xm>>Ym>>Zm>>junk2;}
        if(i==numLines-1){file>>junk2>>signe>>junk2>>junk2>>Xp>>Yp>>Zp>>junk2>>junk2>>Xm>>Ym>>Zm;}
        Xm=Xm/1000.;
        Ym=Ym/1000.;
        Zm=Zm/1000.;

        //checking if muon is in central detector
        a=Xp*Xp+Yp*Yp+Zp*Zp;
        b=2.*(Xp*Xm+Yp*Ym+Zp*Zm);
        c=Xm*Xm+Ym*Ym+Zm*Zm-R_O*R_O;
        d=b*b-4*a*c;
        E=TMath::Sqrt(mass*mass+(Xp*Xp+Yp*Yp+Zp*Zp));
        if(d>0){
            K[0]=(-b+TMath::Sqrt(d))/(2*a);
            K[1]=(-b-TMath::Sqrt(d))/(2*a);
            if(K[0]*K[1]<0){check_0=3; cout<<i<<" muon generated in the central detector : not taken into account"<<endl;}
            if(check_O!=3){
                numInDet_O++;
                Length_O=TMath::Sqrt( Xp*Xp*(K[1]-K[0])*(K[1]-K[0]) +  Yp*Yp*(K[1]-K[0])*(K[1]-K[0]) +Zp*Zp*(K[1]-K[0])*(K[1]-K[0]) );
                LiRate=factor*pow(E,0.74)*Length_O*Rate;
                check_O=1;
            }
        }//if d
        if(d<=0){Length_O=0;LiRate=0;}

        //Reinitialisation of variables
        K[0]=0;
        K[1]=0;
        a=0;
        b=0;
        c=0;
        d=0;

        //Checking if muon is in WP
        a=Xp*Xp+Yp*Yp;
        b=2.*(Xm*Xp+Ym*Yp);
        c=Xm*Xm+Ym*Ym-R_WP*R_WP;
        d=b*b-4*a*c;
        if(d>0){

            //Reinitialisation of variables
            Double_t A=0;
            Double_t B=0;

            K[0]=(-b-TMath::Sqrt(d))/(2*a);
            K[1]=(-b+TMath::Sqrt(d))/(2*a);
            if(Zm+K[0]*Zp>Zm+K[1]*Zp){Double_t chien=0; chien=K[0]; K[0]=K[1]; K[1]=chien;}     //We impose Z1 < Z2

            //Muon passing through both sides of WP
            if(Zm+K[0]*Zp>-H/2. && Zm+K[1]*Zp<H/2.){
                check_WP=1;
                if(K[0]*K[1]<0){check_WP=3;}
            }// if side

            //Muon passing in theoretical infinite cylinder, but above or under the real WP
            if(check_WP==0 && (Zm+K[0]*Zp>H/2. || Zm+K[1]*Zp<-H/2.)){check_WP=3;}

            //Muon passing through the top or the bottom of the WP
            if(check_WP==0){

                if(Zm+K[0]*Zp<-H/2.){
                    A=((-H/2.)-Zm)/Zp;
                    check_WP=2;
                }// if < H/2
                if(Zm+K[1]*Zp>H/2.){
                    B=((H/2.)-Zm)/Zp;
                    check_WP+=4;
                }// if > H/2
                if( (A*B<0 && check_WP==6) || (A*K[1]<0 && check_WP==2) || (B*K[0]<0 && check_WP==4) ){check_WP=3;}
            }// if check_WP 0

            //Muon passing in the WP
            if(check_WP!=3){
                numInDet_WP++;
                if(check_WP==1){
                    Length_WP=TMath::Sqrt( Xp*Xp*(K[1]-K[0])*(K[1]-K[0]) +  Yp*Yp*(K[1]-K[0])*(K[1]-K[0]) +Zp*Zp*(K[1]-K[0])*(K[1]-K[0]) ) - Length_O;
                }
                if(check_WP==2){
                    Length_WP=TMath::Sqrt( Xp*Xp*(K[1]-A)*(K[1]-A) +  Yp*Yp*(K[1]-A)*(K[1]-A) +Zp*Zp*(K[1]-A)*(K[1]-A) ) - Length_O;
                }
                if(check_WP==4){
                    Length_WP=TMath::Sqrt( Xp*Xp*(B-K[0])*(B-K[0]) +  Yp*Yp*(B-K[0])*(B-K[0]) +Zp*Zp*(B-K[0])*(B-K[0]) ) - Length_O;
                }
                if(check_WP==6){
                    Length_WP=TMath::Sqrt( Xp*Xp*(B-A)*(B-A) +  Yp*Yp*(B-A)*(B-A) +Zp*Zp*(B-A)*(B-A) ) - Length_O;
                }
                A=0; B=0;
                cout<<Length_WP<<endl;
                //LiRate=factor*pow(E,0.74)*Length_O*Rate;
            }// if check_WP!=3

        } //if d
        if(d<=0){Length_WP=0;}

        LiRateInt+=LiRate;
        if(i==numLines-1){LiRateInt=LiRateInt/numLines;}
        myTree->Fill();
    }// for i
    cout<<"Number of muons passing in central detector : "<<numInDet_O<<endl;
    cout<<"Number of muons passing in water pool : "<<numInDet_WP<<endl;
    cout<<"Lithium rate : "<<LiRateInt<<endl;
    cout<<"Time : "<<numLines/Rate<<endl;

    file.close();
    myTree->Write();
    rootfile->Close();

    //plot("./rootfile_raw_2.root",0);

    return;

}//main

