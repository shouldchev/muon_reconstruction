#include "libPerso.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void plot (const char *filepath, int ID){

    Double_t Rate=3.5;
    Double_t factor=0.0215;
    Int_t count=0;

    //Name of the PNG file
    TString output;
    output=filepath;
    output.Append(".png");

    //Tree from which to plot
    TChain *myChain=new TChain ("myTree");
    myChain->Add(filepath);

    //Declaration of the variables from the tree
    Double_t e;
    Int_t check;
    Double_t length[10];
    //Double_t liRate[10];

    //Declaration of their branches
    TBranch *b_e;
    TBranch *b_length;
    //TBranch *b_liRate;
    TBranch *b_check;

    //Assigning the branches
    myChain->SetBranchAddress ("check",&check,&b_check);
    myChain->SetBranchAddress ("Energy",&e,&b_e);
    //myChain->SetBranchAddress ("LiRate",&liRate,&b_liRate);
    myChain->SetBranchAddress ("Length",&length,&b_length);

    //Declaring histograms
    TH1D *h_energy = new TH1D("h_energy","",50,0.,1500);
    TH1D *h_length = new TH1D("h_length","",50.,0.,40);
    //TH1D *h_LiRate = new TH1D("h_LiRate","",100.,0.,400);

    for (Int_t i=0 ; i<myChain->GetEntries() ; i++){
        myChain->GetEntry(i);
        h_energy->Fill(e);
        if(check==1){
            count++;
            h_length->Fill(length[0]);
            //h_LiRate->Fill(liRate[0]);
        }//if check
    }//for i
    cout<<"Mean energy : "<<h_energy->GetMean(1)<<endl;
    if(ID==0){
        cout<<"Mean simulated lithium rate : "<<factor*pow(h_energy->GetMean(1),0.74)*h_length->GetMean(1)*Rate<<endl;
        ofstream lithium_rate("./lithiumrate.dat");
        lithium_rate<<factor*pow(h_energy->GetMean(1),0.74)*h_length->GetMean(1)*Rate<<" "<<count<<endl;
        lithium_rate.close();
    }//if ID
    if(ID==1){
        ifstream lithium_rate;
        lithium_rate.open("./lithiumrate.dat",ios::in);
        Double_t Li;
        Int_t nevents=myChain->GetEntries();
        Double_t Nevents;
        lithium_rate>>Li>>Nevents;
        cout<<nevents<<" "<<Nevents<<endl;
        lithium_rate.close();
        Double_t li=(nevents/Nevents);
        li=factor*pow(h_energy->GetMean(1),0.74)*h_length->GetMean(1)*Rate*li;
        cout<<"Mean reconstructed lithium rate : "<<li<<endl;
        cout<<"Ratio reconstructed/simulated : "<<li/Li<<endl;
    }//if ID
    else{cout<<"Wrong ID"<<endl; return;}//else

    //Creating Canvas
    /*TCanvas *mycan=new TCanvas ("mycan","mycan",1500,500);
    mycan->Divide(2,1);
    mycan->cd(1);
    h_energy->Draw();
    h_energy->SetTitle("Energy;Energy (GeV);counts");
    mycan->cd(2);
    h_length->Draw();
    h_length->SetTitle("Length;Length (m);counts");
    mycan->cd(3);
    h_LiRate->Draw();
    h_LiRate->SetTitle("Lithium production rate in Juno;Li rate (events/day);counts");

    mycan->Print(output.Data());*/

    return;
}//main
