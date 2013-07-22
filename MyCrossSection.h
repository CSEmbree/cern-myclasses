#ifndef __MY_CROSS_SECTION_H
#define __MY_CROSS_SECTION_H

#include <stdlib.h> // exit()
#include <iostream> // needed for io
#include <sstream>  // needed for internal io
#include <vector>
#include <string>
//#include <map>

using std::string;
#include <sys/time.h> 

#include "root.h"
//#include "MySubProcess.h"
#include "MyData.h"
#include "MyFrame.h"
#include "MyPDF.h"

#include "appl_grid/appl_grid.h"
#include "appl_grid/mcfmw_pdf.h"
#include "appl_grid/generic_pdf.h"
//#include "theory_error_info.h"

#include "OptionHandler.h"
#include "MyPDF.h"

//#include "appl_grid/appl_pdf.h"



class MyCrossSection {
 private:

  bool debug; 

  string crosssectionname;
  string glabel;
  string subprocesssteername;


            
  double renScaleValUp;
  double renScaleValDefault;
  double renScaleValDown;
  double facScaleValUp;
  double facScaleValDefault;
  double facScaleValDown;
  //string ntupdirinput;
  //string ntupdiroutput;

  string ntupname;
  string steername;

  string datanamedir;
  string gridnamedir;


  double xlabel;
  double ylabel;

  std::vector<string> gridname;
  std::vector<string> dataname;
  std::vector<string> corrname;
  std::vector<string> vardesc;
  std::vector<int> markerstyle;
  std::vector<int> markercolor;
  std::vector<double> mcscalex; // factor to scale x-values
  std::vector<double> datascalex; // factor to scale x-values
  std::vector<double> scaley; // factor to scale y-values
  std::vector<int> frameid;   // grids/data with same frameid will be overlayed
  std::vector<int> divideid;  // grids/data with same divideid will be divided
  std::vector<int> refhistlinestyle; // line style of reference histogram
  std::vector<int> refhistlinecolor; // line color of reference histogram
  std::vector<long unsigned int> events; // number of events from grid
  std::vector<MyFrame*> framepointer; // pointer to MyFrame
  
  //std::vector<string> pdfdata;
  std::vector<std::vector<string> > pdfdata;
  /*
  std::vector<std::vector<string> > pdfdata;
  std::vector<string> pdfsteering;
  */
  std::vector<std::vector<MyPDF*> > t_mypdf; //first part of vec is grid id, second is pdf id


  int processnumber; // process number foresee setting from reaction

  string pdf_function;

  std::vector<appl::grid*> mygrid;   // grid vector  
  //std::vector<appl::grid> mygrid;  // grid vector
  std::vector<bool> isBooked;        // flag grid already booked 
  std::vector<MyData*> mydata;       // information about data from steering file
  generic_pdf *mypdf;
  
    std::vector<string>* ParseString(std::string rawData, char delimeter); //For parsing steering input


 public:
  //std::vector<int> PDFSetCodes_vec;
  bool do_PDFBand;
  bool do_AlphaS;
  bool do_RenormalizationScale;
  bool do_FactorizationScale;
  int ErrorSize;

  MyCrossSection(char name[100]);

  generic_pdf * GetSubProcess() { return mypdf;};
  void SetSubProcess(generic_pdf *subpro) { mypdf=subpro; return;};
  void split_string(std::string str, std::vector<std::string>& split_results, std::string delimiters);

  void Initialize();
  void ReadSteering(char fname[100]);
  //void ReadSteeringOptions(char fname[100]);
  bool file_exists(const string& s);
  
    //return number pdf pdfs for a particular grid
  int GetNPDF(int igrid) {
    return t_mypdf.at(igrid).size();
  };
  
  //return the total number of pdfs (can have multiple per grid name)
  int GetNPDF() {
    int pdfCount=0;
    for(int x=0; x<t_mypdf.size(); x++) {
        for(int y=0; y<t_mypdf.at(x).size(); y++) {
            pdfCount++;
        }
    }
    
    return pdfCount;
  };
  
  /*
  int GetNumPDF() {
    return pdfdata.size();
  };
  
  int GetNumPDFForGrid(int igrid) {
    if(igrid>pdfdata.size()) {
        std::cout<<" MyCrossSection:: GetNumPDFForGrid(int): ERROR: mypdf not found for igrid: "<<igrid<<std::endl;
        return -1;
    }
    return pdfdata.at(igrid).size();
  };
  */

  int reflinestyle;
  void SetLineColor(int igrid, int rl) {refhistlinecolor[igrid]=rl; return;};
  void SetLineStyle(int igrid, int rl) {refhistlinestyle[igrid]=rl; return;};

  void Normalise(TH1D* h, double yscale, double xscale, bool normtot);
  int GetNGrid(){return gridname.size();};
  
  void mypdfInitializeErrorGraphs(int igrid, int ipdf=-1);
  void mypdfCalcSystErrors(int igrid, int ipdf=-1);
  void mypdfGetRatioToTH1(TH1D* href, int igrid, int ipdf=-1);

  
  //string GetNtupDirInput(){ return ntupdirinput;};
  //string GetNtupDirOutput(){ return ntupdiroutput;};
  string GetNtupName(){ return ntupname;};

  string GetGridName(int igrid){ 
    if(gridnamedir.compare("")==0)
        return gridname[igrid];
    else
        return gridnamedir+"/"+gridname[igrid];
  };
  
  
  std::vector<string> *GetPDFData(int igrid) {
    if(igrid>pdfdata.size()) {
        std::cout<<" MyCrossSection:: GetPDFData: ERROR: pdfdata not found for igrid: "<<igrid<<std::endl;
        return NULL;
    }
    else
        return &(pdfdata.at(igrid));
  };
  
  vector<std::vector<std::string> > *GetPDFData() {
    if(pdfdata.size()==0)
        return NULL;
    else
        return &pdfdata;
  };
  
  
  std::vector<MyPDF*> *GetMyPDF(int igrid) {
    if(igrid>t_mypdf.size()) {
        std::cout<<" MyCrossSection:: GetMyPDF(int): ERROR: mypdf not found for igrid: "<<igrid<<std::endl;
        return NULL;
    }
    else
        return &(t_mypdf.at(igrid));
  };
  
  MyPDF *GetMyPDF(int igrid, int ipdf) {
    if(igrid>t_mypdf.size() || ipdf>t_mypdf.at(igrid).size()) {
        std::cout<<" MyCrossSection:: GetMyPDF(int, int): ERROR: mypdf not found for igrid: "<<igrid<<" and ipdf: "<<ipdf<<std::endl;
        return NULL;
    }
    else
        return t_mypdf.at(igrid).at(ipdf);
  };
  
  std::vector< std::vector<MyPDF*> > *GetMyPDF() {
    if(t_mypdf.size()==0) {
        std::cout<<" MyCrossSection:: GetMyPDF: ERROR: mypdf not found for igrid: "<<std::endl;
        return NULL;
    }
    else
        return &t_mypdf;
  };
  
  
  TString GetTStringGridName(int igrid){ return ((TString) (gridnamedir+"/"+gridname[igrid]));  };
  TString GetVarDesc(int igrid){ return ((TString) (vardesc.at(igrid))); };
 /*
  string GetGridName(int i){
   TString name=TString(gridname[i]); 
   name.ReplaceAll(".txt","");
   return string(name.Data());
  };

  string GetGridFileName(int i){
    //return string(this->GetGridName(i).Data())+".root";
    return this->GetGridName(i)+".root";
  }


  string GetGridFullFileName(int i){
    //return ntupdirinput+"/"+this->GetGridFileName(i);
    return this->GetGridDir()+"/"+this->GetGridFileName(i);
  }
 */
  MyData *GetMyData(int igrid){ 
   if (igrid>mydata.size()){
    cout<<" MyCrossSection::GetMyData something is wrong mydata too short igrid= "<<igrid<<endl;
    return 0;
   }
   return mydata[igrid];
  };

  int GetProcessNumber() {return processnumber;};
  void SetProcessNumber(int nproc) {processnumber=nproc; return;};

  string GetGridDir(){ return gridnamedir;};
  bool GetDataOk(int igrid){
   bool flag=true;
   //if (debug) cout<<" MyCrossSection::GetDataOk i= "<<i<<" dataname.size()= "<<dataname.size()<<endl;
   if (dataname.size()<=igrid) flag=false;
   //if (debug) cout<<" MyCrossSection::GetDataOk flag= "<<flag<<endl;
   return flag;
  };
  
  string GetDataName(int igrid){
   string bad="";
   if (!this->GetDataOk(igrid)) return bad;
   else return datanamedir+"/"+dataname[igrid];
  };
  
  string GetCorrelationsName(int igrid){ 
   if (corrname.size()<=igrid) {
    cout<<" MyCrossSection::GetCorrelationsName corrname vector too short or correlation matrix not found ! igrid= "<<igrid<<endl;
    return "";
   }
   return datanamedir+"/"+corrname[igrid]; 
  };

  int GetMarkerStyle(int igrid) {
   if (markerstyle.size()<=igrid) {
    cout<<" MyCrossSection::GetMarkerStyle something is wrong markserstyle too short igrid= "<<igrid<<endl;
    return 0;
   }
   return markerstyle[igrid];
  };
  
  int GetMarkerColor(int igrid) {
   if (markercolor.size()<=igrid) {
    cout<<" MyCrossSection::GetMarkerColor something is wrong markercolor too short igrid= "<<igrid<<endl;
    return 0;
   }
   return markercolor[igrid];
  };
  
  void SetStyleMarker (int igrid, int ms) {
   if (igrid>markerstyle.size()){
     cout<<" MyCrossSection::SetStyleMarker something is wrong markerstyle too short igrid= "<<igrid<<endl;
   } else markerstyle[igrid]=ms; 
   return;
  };
  
  void SetColorMarker (int igrid, int mc) {
   if (igrid>markerstyle.size()){
     cout<<" MyCrossSection::SetColorMarker something is wrong markercolor too short igrid= "<<igrid<<endl;
   } else markercolor[igrid]=mc; 
   return;
  };

  double GetDataScaleX(int igrid){ 
   if (datascalex.size()<=igrid) {
     cout<<" MyCrossSection::GetScaleX something is wrong markercolor too short igrid= "<<igrid<<endl;
     return 1.;
   }
   return datascalex[igrid];
  };
  
  double GetMCScaleX(int igrid){ 
   if (datascalex.size()<=igrid) {
     cout<<" MyCrossSection::GetScaleX something is wrong markercolor too short igrid= "<<igrid<<endl;
     return 1.;
   }
   return mcscalex[igrid];
  };
  
  double GetScaleY(int igrid){ 
   if (scaley.size()<=igrid) {
     cout<<" MyCrossSection::GetScaleY something is wrong scaley too short igrid= "<<igrid<<endl;
     return 1.;
   }
   return scaley[igrid];
  };
  
  int GetFrameID(int igrid){ 
   if (frameid.size()<=igrid) {
     cout<<" MyCrossSection::GetFrameID something is wrong frameid too short igrid= "<<igrid<<endl;
     return -1;
   }
   return frameid[igrid];
  };
  
  MyFrame *GetMyFrame(int igrid ) {
   if (framepointer.size()<=igrid) {
     cout<<" MyCrossSection::GetMyFrame something is wrong divideid too short igrid= "<<igrid<<endl;
     return 0;
   }
   return framepointer[igrid];
  };
  
  int GetDivideID(int igrid){ 
   if (divideid.size()<=igrid) {
     cout<<" MyCrossSection::GetDivideID something is wrong divideid too short igrid= "<<igrid<<endl;
     return 1.;
   }
   return divideid[igrid];
  };

  int GetFrameNumber();

  appl::grid *GetGrid(int igrid){
    //  if (!mygrid[igrid]) cout<<" MyCrossSection::GetReference mygrid["<<igrid<<"] not filled ! "<<endl;
   return mygrid[igrid];
  };

  TH1D *GetReference(int igrid);
  void Draw(int igrid);

  
  void DrawErrors(TString x_title, float x_min, float x_max, bool first_of_canv, int error_code);
  void DrawError(int igrid, int ipdf, TString x_title, float x_min, float x_max, bool first_of_canv, int error_code);
  
  void DrawData(int igrid){
   mydata[igrid]->DrawData();
   return;
  };
 void SetLabelX(double x) {xlabel=x; return;};
 void SetLabelY(double y) {ylabel=y; return;};

 void DrawReference(int igrid);
 TGraphAsymmErrors* GetReferenceRatio(int igrid);
 TH1D* GetNormalisedReference(int igrid);

 void DrawinFrame(int iframe);

 double GetObsMax(int igrid){
  if (!mydata[igrid])  cout<<" GetObsMax grid not found= mydata["<<igrid<<"]"<<endl;
  double obsmax=mydata[igrid]->GetMaxX(); 
  //if (debug) cout<<" obsmax= "<<obsmax<<endl;
  return obsmax;
 }

  TH1D *GetHisto(int igrid,string name="htest"){
   if (!mydata[igrid]) cout<<" MyCrossSection::GetHisto mydata not filled ! "<<endl;
   int nobs=mydata[igrid]->GetNBins(); 
   if (debug) 
    cout<<" MyCrossSection::GetHisto nobs= "<<nobs<<endl;
   double *xbins=mydata[igrid]->GetBins();

   string hname=name+this->GetGridName(igrid);
   //cout<<" book histogram hname= "<<hname<<endl;
   TH1D *htest = new TH1D(hname.c_str(),hname.c_str(),nobs,xbins);
   if (!htest) cout<<" htest["<<igrid<<"] not found "<<endl;
   //htest->Print("all");
   return htest;
  }; 

  double GetTotalEventNumber(int igrid){
   if (events.size()<=igrid) return 0; 
   return events[igrid];
  };

  TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2, Int_t noerr=1);
  TGraphAsymmErrors* TH1TOTGraphAsymm(TH1 *h1);

  void Print();
};


#endif
