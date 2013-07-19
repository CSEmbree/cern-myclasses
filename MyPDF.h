#ifndef __MY_PDF_H
#define __MY_PDF_H

#include <iostream>
#include <iomanip>
#include <string>
#include "OptionHandler.h"

//root
#include <TH1D.h>
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TMatrixT.h"

//appl_grid
#include "appl_grid/appl_grid.h"
#include "appl_grid/generic_pdf.h"

//LHAPDF
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF.h"

//TEMP
//#include "VariableDefinitions.h"

#define DEFAULT -1

using namespace std;


//TEMP delcared in Variables.h file needed for stand alone myPDF use but not when using older MyCrossSection file
enum enum_RenScales {e_RenScale0p5, e_RenScale1p0, e_RenScale2p0, e_n_RenScaleVals}; //e_n_* dictates future size
enum enum_FacScales {e_FacScale0p5, e_FacScale1p0, e_FacScale2p0, e_n_FacScaleVals};
enum scales{UP=0, DEF, DOWN, n_SCALES};


class MyPDF {

    public:
        //VARIABLES
        string optionsFileName;
        const string defaultOptionsFileName; //default is to look in current directory for this file
        TH1D* h_qqbar_prenorm;
        TH1D* h_gg_prenorm;
        TH1D* h_tot_prenorm;
        TH1D* h_qqbar;
        TH1D* h_gg;
        TH1D* h_tot;
        TH1D* h_gg_frac;
        TH1D* h_qqbar_frac;
        TGraphAsymmErrors *h_PDFBand_results;
        TGraphAsymmErrors *h_PDFBand_results_ratio_to_ref;
        TGraphAsymmErrors *h_AlphaS_results;
        TGraphAsymmErrors *h_AlphaS_results_ratio_to_ref;
        TGraphAsymmErrors *h_RenormalizationScale_results;
        TGraphAsymmErrors *h_RenormalizationScale_results_ratio_to_ref;
        TGraphAsymmErrors *h_FactorizationScale_results;
        TGraphAsymmErrors *h_FactorizationScale_results_ratio_to_ref;
        TGraphAsymmErrors *h_TotError_results;
        TGraphAsymmErrors *h_TotError_results_ratio_to_ref;
        string calc_desc;
        
        //METHODS
        MyPDF(bool _debug=false); //default constructor
        MyPDF(string _gridName, double _xscale=1.0, string _steeringFileName="steering_mypdf.txt", bool _debug=false);
        virtual ~MyPDF() { CleanUpMyPDF(); }; //destructor
        void CleanUpMyPDF();
        void Initialize();
        void Print();
        void PrintKnownOptionsForSteering();
        void PrintFoundOptionsFromSteering();
        void ReadSteering(const string _fileName);
        TH1 *TH1NormToTot(TH1 *h1, double _yscale=1.0, double _xscale=1.);
        void InitializeErrorGraphs();
        void CalcSystErrors();
        void CalcPDFBandErrors();
        void CalcAlphaSErrors();
        void CalcRenormalizationScaleErrors();
        void CalcFactorizationScaleErrors();
        void CalcTotErrors();
        void GetRatioToTH1(TH1D* href);
        TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2, Int_t noerr=1);
        TGraphAsymmErrors* TH1TOTGraphAsymm(TH1 *h1);
        void CalcChi2(TGraphAsymmErrors *g_theory, TGraphAsymmErrors *g_data, TMatrixT<double> data_cov_matrix, double &chi2);
        
        
        //accessor methods
        bool isDebugOn() const{return debug;};
        string getSteeringFilePath() const{return steeringFilePath;};
        string getSteeringFileDir() const{return steeringFileDir;};
        string getSteeringFileName() const{return steeringFileName;};
        string getPDFtype() const{return PDFtype;};
        string getPDFname() const{return PDFname;};
        int getNumPDFMembers() const{return n_PDFMembers;};
        int getFillStyleCode() const{return fillStyleCode;};
        int getFillColorCode() const{return fillColorCode;};
        string getPDFBandType() const{return PDFBandType;};
        string getPDFErrorType() const{return PDFErrorType;};
        string getPDFErrorSize() const{return PDFErrorSize;};
        
        string getRenScaleNameUp() const{return renScaleNameUp;};
        string getRenScaleNameDefault() const{return renScaleNameDefault;};
        string getRenScaleNameDown() const{return renScaleNameDown;};
        double getRenScaleValUp() const{return renScaleValUp;};
        double getRenScaleValDefault() const{return renScaleValDefault;};
        double getRenScaleValDown() const{return renScaleValDown;};
        string getFacScaleNameUp() const{return facScaleNameUp;};
        string getFacScaleNameDefault() const{return facScaleNameDefault;};
        string getFacScaleNameDown() const{return facScaleNameDown;};
        double getFacScaleValUp() const{return facScaleValUp;};
        double getFacScaleValDefault() const{return facScaleValDefault;};
        double getFacScaleValDown() const{return facScaleValDown;};
        
        bool getDoPDFBand() const{return do_PDFBand;};
        bool getDoAlphaS() const{return do_AlphaS;};
        bool getDoRenormalizationScale() const{return do_RenormalizationScale;};
        bool getDoFactorizationScale() const{return do_FactorizationScale;};
        bool getDoTotError() const{return do_TotError;};
        
        //mutator methods
        void setDebug(bool _debug);
        void setGridName(string _gridName);
        void setSteeringFilePath(string _steeringFilePath);
        void setSteeringFileDir(string _steeringFileDir);
        void setSteeringFileName(string _steeringFileName);
        void setPDFtype(string _PDFtype);
        void setPDFname(string _PDFname);
        void setNumPDFMembers(int _n_PDFMembers);
        void setFillStyleCode(int _fillStyleCode);
        void setFillColorCode(int _fillColorCode);
        void setPDFBandType(string _PDFBandType);
        void setPDFErrorType(string _PDFErrorType);
        void setPDFErrorSize(string _PDFErrorSize);
        void setRenScaleValUp(double _renScaleVal);
        void setRenScaleValDefault(double _renScaleVal);
        void setRenScaleValDown(double _renScaleVal);
        void setFacScaleValUp(double _facScaleVal);
        void setFacScaleValDefault(double _facScaleVal);
        void setFacScaleValDown(double _facScaleVal);
        void setOptionsFileName(string _optionsFileName);

    private:
        //VARIABLES
        bool debug;
        string steeringFilePath;
        string steeringFileDir;
        string steeringFileName;    //name of steering file
        string PDFtype;             //general name for PDF EX: "MSTW2008nlo"
        string PDFname;             //specific name for PDF EX: "MSTW2008nlo68cl"
        int n_PDFMembers;
        int fillStyleCode;
        int fillColorCode;
        string PDFBandType;
        string PDFErrorType;
        string PDFErrorSize;
        string renScaleNameUp;
        string renScaleNameDefault;
        string renScaleNameDown;
        string renScaleNames[n_SCALES];
        double renScaleValUp;
        double renScaleValDefault;
        double renScaleValDown;
        double renScaleVals[n_SCALES];
        string facScaleNameUp;
        string facScaleNameDefault;
        string facScaleNameDown;
        string facScaleNames[n_SCALES];
        double facScaleValUp;
        double facScaleValDefault;
        double facScaleValDown;
        double facScaleVals[n_SCALES];
        OptionHandler *myOptions;
        bool f_PDFBandType;
        bool f_PDFErrorSize;
        string pdfSetPath;
        
        //start old from therory_error_info/calc
        appl::grid *my_grid;
    
        std::vector<TH1D*> h_errors_RenormalizationScale;
        std::vector<TH1D*> h_errors_FactorizationScale;
        std::vector<TH1D*> h_errors_PDFBand;
        std::vector<TH1D*> h_errors_prenorm;
        std::vector<TH1D*> h_errors_AlphaS;
        std::vector<TH1D*> h_errors_AlphaS_prenorm;
        string gridName;
        
        double xscale;
        bool do_PDFBand;
        bool do_AlphaS;
        bool do_RenormalizationScale;
        bool do_FactorizationScale;
        bool do_TotError;
        //end old from therory_error_info/calc
        
        //METHODS
        bool fileExists(const string _fileName);
        void setVariablesDefault();
        void setSteeringFileNameAndDir(const string _path);
        string GetEnv( const string & var);
};

const string defaultOptionsFileName="options_mypdf.txt";
const string deafultPDFSetPath="PDFsets";

/*
//currently decalred in local directory LHAPDF.h
extern "C" void evolvepdf_(const double& , const double& , double* );
extern "C" double alphaspdf_(const double& Q);
*/

/*
void getPDF(const double& x, const double& Q, double* xf) {
    evolvepdf_(x, Q, xf);
    //evolvePDF( x, Q, xf);        //// calls LHAPDF
}

double alphasPDF(const double& Q) {
    return alphaspdf_(Q);
}
*/

#endif
