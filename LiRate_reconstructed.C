#include "LiRate.h"
#include "libPerso.h"
#include "plot.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#define C 299792458

using namespace std;

void LiRate_reconstructed (){
    //reconstruction rootfile and tree
    TFile *file=new TFile("./Tracking_detsim_3wall.root","");
    TTree *myChain=(TTree*)file->Get("TT");
    //simulation rootfile and tree
    TFile *file_raw=new TFile("./rootfile_raw_2.root","");    
    TTree *tree_raw=(TTree*)file_raw->Get("myTree");
    const Int_t nentries=myChain->GetEntries();
    const Int_t nentries_raw=tree_raw->GetEntries();
    cout<<"Entries in simulation : "<<nentries_raw<<" Entries in reconstruction : "<<nentries<<endl;

    //Declaring variables from simulation
    Double_t LiRateInt_raw;
    Double_t numInDet_raw;

    //Declaring the variables from reconstruction
    Int_t evtID, NTracks;
    Double_t Coeff0[50];
    Double_t Coeff1[50];
    Double_t Coeff2[50];
    Double_t Coeff3[50];
    Double_t Coeff4[50];
    Double_t Coeff5[50];

    //Creating branches from reconstruction
    TBranch *b_evtID, *b_NTracks, *b_Coeff0, *b_Coeff1, *b_Coeff2, *b_Coeff3, *b_Coeff4, *b_Coeff5;
    //Creating branches from simulation
    TBranch *b_E, *b_check, *b_LiRateInt_raw, *b_numInDet_raw;

    //Getting the branches from the reconstruction rootfile
    myChain->SetBranchAddress ("evtID",&evtID,&b_evtID);
    myChain->SetBranchAddress ("NTracks",&NTracks, &b_NTracks);
    myChain->SetBranchAddress("Coeff0", &Coeff0,&b_Coeff0);
    myChain->SetBranchAddress("Coeff1", &Coeff1,&b_Coeff1);
    myChain->SetBranchAddress("Coeff2", &Coeff2,&b_Coeff2);
    myChain->SetBranchAddress("Coeff3", &Coeff3,&b_Coeff3);
    myChain->SetBranchAddress("Coeff4", &Coeff4,&b_Coeff4);
    myChain->SetBranchAddress("Coeff5", &Coeff5,&b_Coeff5);

    //Getting the branches from the simulation rootfile
    tree_raw->SetBranchAddress("Energy",&E,&b_E);
    tree_raw->SetBranchAddress("check",&check,&b_check);
    tree_raw->SetBranchAddress("LiRateInt",&LiRateInt_raw,&b_LiRateInt_raw);
    tree_raw->SetBranchAddress("numInDet",&numInDet_raw,&b_numInDet_raw);

    Double_t R=35.4/2.;
    Double_t factor=0.0215;
    Int_t numInDet=0;

    //new rootfile
    rootfile=new TFile("./rootfile_reconstructed.root","RECREATE");
    myTree=new TTree ("myTree","");

    //Creating the branches for the new rootfile
    TBranch *b_evtID_2=myTree->Branch("evtID",&evtID,"evtID/i");
    TBranch *b_posX=myTree->Branch("posX",&Xm,"posX/d");
    TBranch *b_posY=myTree->Branch("posY",&Ym,"posY/d");
    TBranch *b_posZ=myTree->Branch("posZ",&Zm,"posZ/d");
    TBranch *b_momX=myTree->Branch("momX",&Xp,"momX/d");
    TBranch *b_momY=myTree->Branch("momY",&Yp,"momY/d");
    TBranch *b_momZ=myTree->Branch("momZ",&Zp,"momZ/d");
    TBranch *b_e=myTree->Branch("Energy",&E,"Energy/d");
    TBranch *b_check_2=myTree->Branch("check",&check,"check/i");
    TBranch *b_Length=myTree->Branch("Length",&Length,"Length/d");
    TBranch *b_LiRate=myTree->Branch("LiRate",&LiRate,"LiRate/d");

    for (int i = 0 ; i<nentries ; i++){
        myChain->GetEntry(i);
        if (NTracks!=0){
            tree_raw->GetEntry(i);
            Xm=Coeff0[0];
            Ym=Coeff1[0];
            Zm=Coeff2[0];
            Xp=Coeff3[0];
            Yp=Coeff4[0];
            Zp=Coeff5[0];
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
                LiRate=factor*pow(E,0.74)*Length;
                LiRateInt+=LiRate;
                myTree->Fill();
            }//if d
        }//if NTracks
    }//for i
    tree_raw->GetEntry(tree_raw->GetEntries()-1);
    cout<<"Number of muons in detector : "<<numInDet_raw<<endl;
    cout<<"Lithium rate from simulation : "<<LiRateInt_raw<<endl;
    LiRateInt=(LiRateInt*13)/nentries;
    cout<<"Number of reconstructed muons : "<<numInDet<<endl;
    cout<<"Lithium rate from TT : "<<LiRateInt<<endl;
    cout<<"ratio : "<<LiRateInt/LiRateInt_raw<<endl;
    myTree->Write();
    rootfile->Close();
    file_raw->Close();
    file->Close();

    //plot("./rootfile_reconstructed.root",1);

    return;
}//main
