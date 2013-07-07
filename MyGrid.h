#ifndef __MY_GRID_H
#define __MY_GRID_H

#include <stdlib.h> // exit()
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <map>

using std::string;
#include <sys/time.h> 

#include "root.h"
//#include "MySubProcess.h"
#include "MyEvent.h"
#include "MyData.h"
#include "appl_grid/appl_grid.h"
#include "appl_grid/generic_pdf.h"

//#include "appl_grid/appl_pdf.h"

//static const int __nf__   = 5;
//static const int __nf2__  = 11;


class MyGrid {
 private:

  bool debug; 
  bool bookrefsubprocess;

  long unsigned int alluncorrevents;
  long unsigned int allevents;
  std::vector<long unsigned int> uncorrevents;
  std::vector<long unsigned int> events;

  bool newevent;
  int eventid;
  string glabel;

  string ntupdiroutput;
  string ntupdirinput;
  string ntupname;
  //  string filename;
  string subprocesssteername;

  string gridnamedir;
  string steername;
  string versionname; //<--added as optional argument
  
  

  int processnumber; // process number foresee setting from reaction
  int nsub; // number of subprocesses
 
  double *weight;
  // grid architecture -> could be in steering file
  //
  int    nXbins;
  double xLow;
  double xUp;
  int xorder;
  int nQ2bins;
  double q2Low; 
  double q2Up;
  int    qorder;
  int    iorder;

  bool incljets;
  
  // lowest order in alphas	
  int lowest_order;
  // how many loops
  int nloops;

  string pdf_function;

  std::vector<int> istylemarkerdata; // marker of theory prediction
  std::vector<int> icolormarkerdata; // marker of theory prediction
  std::vector<appl::grid*> mygrid;     // grid vector
  std::vector<bool> isBooked;          // flag grid already booked 

  std::vector<string> gridname;

  typedef std::vector<TH1D*> vhref;

  std::vector<vhref> hrefLOsubprocesseshistos;      //saved in outputHistos
  std::vector<vhref> hrefsubprocesseshistos;        //saved in outputHistos
  std::vector<vhref> hrefLOsubprocesseshistossum2;
  std::vector<vhref> hrefsubprocesseshistossum2;
  vhref hreference;                                 //saved in outputHistos
  vhref hreferencesum2;
  vhref hreferencefine;                             //saved in outputHistos
  vhref hreferencefinesum2;

  std::vector<vhref> hrefLOsubprocesseshistostmp;
  std::vector<vhref> hrefsubprocesseshistostmp;
  std::vector<vhref> hrefLOsubprocesseshistossum2tmp;
  std::vector<vhref> hrefsubprocesseshistossum2tmp;
  vhref hreferencetmp;
  vhref hreferencesum2tmp;
  vhref hreferencefinetmp;
  vhref hreferencefinesum2tmp;
  
  //holds vector of vector of reference histograms for retrievel later if two grids are added.
  std::vector< std::vector< vhref > > outputHistos; 
 
  std::vector<MyData*> mydata;        // information about data from steering file

  generic_pdf *mypdf;
  
  int decideSubProcess(int iflav1,int iflav2);

 public:
 //MyGrid(char name[100]);
  MyGrid(string name, string version="");
  
  void AddGrid(MyGrid *myOtherGrid); //add the reference histograms of another grid to this one.

  int GetSubProcess() { return mypdf->GetCurrentSubProcess();};
  string GetSubProcessName(int isub) { return mypdf->GetSubProcessName(isub);};
  void SetSubProcess(generic_pdf *subpro) { mypdf=subpro; return;};

  void Initialize();
  //void ReadSteering(char fname[100]);
  void ReadSteering(string);
  bool file_exists(const string& s);
  void book_grid(int igrid);  // inital grid booking  
  void fill(MyEvent *myevent);
  bool eventcuts(MyEvent *myevent,int igrid);
  void write_grid();   // writes out grid after some events
  void DivideByBinWidth(TH1D* h);
  void Normalise(TH1D* h1, double yscale, double xscale, bool normtot); 
  void getObservable(MyEvent *myevent);
  int GetNGrid(){return gridname.size();};
  string GetNtupDirOutput(){ return ntupdiroutput;};
  string GetInputNtupDir() { return ntupdirinput;};
  void SetInputNtupDir(TString dirinput){ ntupdirinput=dirinput; return;};
  string GetNtupName(){ return ntupname;};

  string GetGridName(int i){
   TString name=TString(gridname[i]); 
   name.ReplaceAll(".txt","");
   return string(name.Data());
  };
  string GetGridFileName(int i){
    //return string(this->GetGridName(i).Data())+".root";
    return this->GetGridName(i)+".root";
  }
  string GetGridVersionName(){ //<--added accessor method to versionname
    return versionname;
  }

  string GetGridFullFileName(int i){
    return this->GetNtupDirOutput()+"/"+this->GetGridFileName(i);
  }

  MyData *GetMyData(int i){ return mydata[i];};

  int GetProcessNumber() {return processnumber;};
  void SetProcessNumber(int nproc) {processnumber=nproc; return;};
  //int GetNumberSubProcess() {return nsub;};
  //void SetNumberSubProcess(int nproc) {nsub=nproc; return;};


  string GetGridDir(){ return gridnamedir;};
  string GetDataName(int i){ return gridnamedir+"/"+gridname[i];};

  int GetStyleMarker(int i) {return istylemarkerdata[i];};
  void SetStyleMarker (int i) {istylemarkerdata.push_back(i); return;};
  int GetColorMarker(int i) {return icolormarkerdata[i];};
  void SetColorMarker (int i) {icolormarkerdata.push_back(i); return;};

  double *GetWeight(){return weight;};
  double GetWeight(int i){return weight[i];};
  void SetWeight(int i,double wgt){weight[i]=wgt; return;};
  void PrintWeight(){ 
   cout<<" MyGrid::PrintWeight: nsub= "<<nsub<<" (zero suppressed print out)"<<endl;
   for (int i=0; i<nsub; i++) 
    if (weight[i]!=0.)
     cout<< " weight["<<i<<"]="<<weight[i]<<endl; 
   return;
  }  

  void ResetWeight(){ for (int i=0; i<nsub; i++) weight[i]=0.; return;};

  appl::grid *GetGrid(int igrid){
   if (!mygrid[igrid]) cout<<" MyGrid::GetGrid mygrid["<<igrid<<"] not filled ! "<<endl;
   return mygrid[igrid];
  }

  TH1D *GetReference(int igrid){
   if (!mygrid[igrid]) cout<<" MyGrid::GetReference mygrid["<<igrid<<"] not filled ! "<<endl;
   TH1D *href=(TH1D*)mygrid[igrid]->getReference();
   if (!href) {
    cout<<" MyGrid::GetReference: reference histo not found igrid="<<igrid<<endl;
   }
   return href;
  }
  
  
  //Accessor method for getting vector<vector<vector<TH1>>> that contains histograms needed for adding
  std::vector<std::vector<vhref> > GetOutputHistograms(){
   if (outputHistos.size()<=0) {
    cout<<" MyGrid::GetOutputHistograms: no output histograms found!"<<endl;
   }
   return outputHistos; 
  }
  

 double GetObsMax(int igrid){
  if (!mydata[igrid])  cout<<" MyGrid::GetObsMax: grid not found= mydata["<<igrid<<"]"<<endl;
  //int nobs=mydata[igrid]->GetNBins(); 
  //cout<<" GetObsMax nobs= "<<nobs<<endl;
  //double *xbins=mydata[igrid]->GetBins();
  double obsmax=mydata[igrid]->GetMaxX(); 
  //if (debug) cout<<" obsmax= "<<obsmax<<endl;
  return obsmax;
 }

  TH1D *GetHisto(int igrid,string name="htest"){
   if (!mydata[igrid]) cout<<" MyGrid::GetHisto mydata not filled ! "<<endl;
   int nobs=mydata[igrid]->GetNBins(); 
   if (debug) 
    cout<<" MyGrid::GetHisto: nobs= "<<nobs<<", for igrid: "<<igrid<<endl;
   double *xbins=mydata[igrid]->GetBins();

   string hname=name+this->GetGridName(igrid)+this->GetGridVersionName(); //<--added +this->GetGridVersionName() to make unique name
   if(debug)
    cout<<" MyGrid::GetHisto: name set to: "<<name<<endl;
   //cout<<" book histogram hname= "<<hname<<endl;
   TH1D *htest = new TH1D(hname.c_str(),hname.c_str(),nobs,xbins);
   if (!htest) cout<<" htest["<<igrid<<"] not found "<<endl;
   //htest->Print("all");
   return htest;
  }; 

  void SetEventId(int evid){
   newevent=false;
   //cout<<" MyGrid::NewEvent: present eventid= "<<eventid<<" evid= "<<evid<<endl;
   if (eventid!=evid){
    eventid=evid;
    newevent=true;
   }
   return;
  };
  bool NewEvent(){return newevent;};

  double GetUncorrellatedEventNumber(){return alluncorrevents;};
  double GetUncorrellatedEventNumber(int igrid){return uncorrevents[igrid];};
  double GetTotalEventNumber(){return allevents;};
  double GetTotalEventNumber(int igrid){return events[igrid];};

  void Print();
  void PrintSubprocessRefHistos(int igrid);
  void PrintRefHistos(int igrid);
  void NormRefHistos(int igrid, double norm);
  TH1D *TH1NormStatError(TH1D *hsum, TH1D *hsum2, double norm);
};


#endif
