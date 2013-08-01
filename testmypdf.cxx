
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

#include "MyPDF.h"

main()
{
    string gridName="grid--TTbar_yttatlas.root";
    string steeringFileName = "steering_mypdf.txt";
    //TFile *fout;
    //string filename="testmypdf_output";
  	//filename+=".root";
  	//fout=new TFile(filename.c_str(),"recreate");
  	
  	TGraphAsymmErrors *h_PDFBand_results;
  	TGraphAsymmErrors *h_AlphaS_results;
    
    
    /////////////////////////////////////////////    
    std::cout<<"****TEST 1 START"<<std::endl;

    MyPDF *mypdf1 = new MyPDF(false);

    mypdf1->Print();
    mypdf1->SetGridName(gridName);
    mypdf1->SetSteeringFilePath(steeringFileName);
    mypdf1->ReadSteering(steeringFileName);
    mypdf1->Initialize();
    
    mypdf1->Print();
    
    mypdf1->InitializeErrorGraphs();
    mypdf1->CalcSystErrors();
        
    h_PDFBand_results=mypdf1->h_PDFBand_results;
        
    h_PDFBand_results->SetTitle("mypdf1_h_pdfband");
    h_PDFBand_results->SetName("mypdf1_h_pdfband");
    h_PDFBand_results->Print("all");
    //h_PDFBand_results->Write();
        

    std::cout<<"****TEST 1 END"<<std::endl;
    /////////////////////////////////////////////    
    std::cout<<"****TEST 2 START"<<std::endl;
    
    MyPDF *mypdf2 = new MyPDF(gridName, 1.0, string("steering_mypdf_pdfband.txt"), false);
    //mypdf2->setDebug(true);
    
    mypdf2->Print();
    mypdf2->InitializeErrorGraphs();
    mypdf2->CalcSystErrors();
    
    h_PDFBand_results=mypdf2->h_PDFBand_results;
    
    h_PDFBand_results->SetTitle("mypdf2_h_pdfband");
    h_PDFBand_results->SetName("mypdf2_h_pdfband");
    h_PDFBand_results->Print("all");
    //h_PDFBand_results->Write();


    std::cout<<"****TEST 2 END"<<std::endl;
    /////////////////////////////////////////////    
    std::cout<<"****TEST 3 START"<<std::endl;

    MyPDF *mypdf3 = new MyPDF(gridName, 1.0, string("steering/steering_mypdf.txt"), true);
    //mypdf3->setDebug(true);
    
    mypdf3->Print();
    mypdf3->InitializeErrorGraphs();
    mypdf3->CalcSystErrors();
    
    h_AlphaS_results=mypdf3->h_AlphaS_results;
    
    h_AlphaS_results->SetTitle("mypdf3_h_alphaS");
    h_AlphaS_results->SetName("mypdf3_h_alphaS");
    h_AlphaS_results->Print("all");
    //h_AlphaS_results->Write();
    
    delete mypdf3;
    
    
    std::cout<<"****TEST 3 END"<<std::endl;
    /////////////////////////////////////////////    
    
    
    //fout->Write();
  	//fout->Close();
  	
}
