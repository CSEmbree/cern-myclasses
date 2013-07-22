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
//enum enum_RenScales {e_RenScale0p5, e_RenScale1p0, e_RenScale2p0, e_n_RenScaleVals}; //e_n_* dictates future size
//enum enum_FacScales {e_FacScale0p5, e_FacScale1p0, e_FacScale2p0, e_n_FacScaleVals};
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
        TGraphAsymmErrors* MyTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2, Int_t noerr=1);
        TGraphAsymmErrors* TH1TOTGraphAsymm(TH1 *h1);
        void CalcChi2(TGraphAsymmErrors *g_theory, TGraphAsymmErrors *g_data, TMatrixT<double> data_cov_matrix, double &chi2);
        
        
        //accessor methods
        bool IsDebugOn() const{return debug;};
        string GetSteeringFilePath() const{return steeringFilePath;};
        string GetSteeringFileDir() const{return steeringFileDir;};
        string GetSteeringFileName() const{return steeringFileName;};
        string GetPDFtype() const{return PDFtype;};
        string GetPDFname() const{return PDFname;};
        int GetNumPDFMembers() const{return n_PDFMembers;};
        int GetFillStyleCode() const{return fillStyleCode;};
        int GetFillColorCode() const{return fillColorCode;};
        string GetPDFBandType() const{return PDFBandType;};
        string GetPDFErrorType() const{return PDFErrorType;};
        string GetPDFErrorSize() const{return PDFErrorSize;};
        
        string GetRenScaleNameUp() const{return renScaleNameUp;};
        string GetRenScaleNameDefault() const{return renScaleNameDefault;};
        string GetRenScaleNameDown() const{return renScaleNameDown;};
        double GetRenScaleValUp() const{return renScaleValUp;};
        double GetRenScaleValDefault() const{return renScaleValDefault;};
        double GetRenScaleValDown() const{return renScaleValDown;};
        string GetFacScaleNameUp() const{return facScaleNameUp;};
        string GetFacScaleNameDefault() const{return facScaleNameDefault;};
        string GetFacScaleNameDown() const{return facScaleNameDown;};
        double GetFacScaleValUp() const{return facScaleValUp;};
        double GetFacScaleValDefault() const{return facScaleValDefault;};
        double GetFacScaleValDown() const{return facScaleValDown;};
        
        bool GetDoPDFBand() const{return do_PDFBand;};
        bool GetDoAlphaS() const{return do_AlphaS;};
        bool GetDoRenormalizationScale() const{return do_RenormalizationScale;};
        bool GetDoFactorizationScale() const{return do_FactorizationScale;};
        bool GetDoTotError() const{return do_TotError;};
        
        //mutator methods
        void SetDebug(bool _debug);
        void SetGridName(string _gridName);
        void SetSteeringFilePath(string _steeringFilePath);
        void SetSteeringFileDir(string _steeringFileDir);
        void SetSteeringFileName(string _steeringFileName);
        void SetPDFtype(string _PDFtype);
        void SetPDFname(string _PDFname);
        void SetNumPDFMembers(int _n_PDFMembers);
        void SetFillStyleCode(int _fillStyleCode);
        void SetFillColorCode(int _fillColorCode);
        void SetPDFBandType(string _PDFBandType);
        void SetPDFErrorType(string _PDFErrorType);
        void SetPDFErrorSize(string _PDFErrorSize);
        void SetRenScaleValUp(double _renScaleVal);
        void SetRenScaleValDefault(double _renScaleVal);
        void SetRenScaleValDown(double _renScaleVal);
        void SetFacScaleValUp(double _facScaleVal);
        void SetFacScaleValDefault(double _facScaleVal);
        void SetFacScaleValDown(double _facScaleVal);
        void SetOptionsFileName(string _optionsFileName);
        void SetDoPDFBand(bool _doit);
        void SetDoAplphaS(bool _doit);
        void SetDoRenormalizationScale(bool _doit);
        void SetDoFactorizationScale(bool _doit);
        void SetDoTotError(bool _doit);

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
        bool FileExists(const string _fileName);
        void SetVariablesDefault();
        void SetSteeringFileNameAndDir(const string _path);
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
