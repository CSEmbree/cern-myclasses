/*
 * Title:    MyPDF
 * Author:   Cameron S. Embree
 * Contact:  CSEmbree@gmail.com
 * Created:  01-Jun-2013
 * Edited:   22-Jun-2013
 * Notes:    Class implimentation suggested by Dr. Carli based on the "theory_error_info.{cxx,h}" class
 */

/*
    TEST AND RUN WITH:

    make testmypdf
    ./testmypdf

    OR

    make plot_new
    ./plot_new
*/

#include "MyPDF.h"


//clean wrapper for evolvepdf_ fortran call
void getPDF(const double& x, const double& Q, double* xf) {
    evolvepdf_(&x, &Q, xf);
    //evolvePDF( x, Q, xf);        //// calls LHAPDF
}

//find and returns envrionment variable contents if it exists
string MyPDF::GetEnv( const string & var ) {
    const char* res= getenv( var.c_str() );
    std::string s = res!=NULL? res:"";
    return s;
}

/*
//currently declared in local directories LHAPDF.h file
double alphasPDF(const double& Q) {
    return alphaspdf_(Q);
}
*/



/******************************************************************
 ** Method Implementations
 ******************************************************************/

//default constructor
MyPDF::MyPDF(bool _debug)
{
    debug=_debug;
    if(debug)std::cout<<" MyPDF::MyPDF: start"<<std::endl;

    //no values are being set for default, so assign default values to avoid crashes
    SetVariablesDefault();

    if(debug)std::cout<<" MyPDF::MyPDF: End"<<std::endl;
}


//overloaded constructor - Use if PDFErrorType is declared in MyPDF steering file.
MyPDF::MyPDF(string _gridName, double _xscale, string _steeringFile, bool _debug)
{
    debug=_debug;
    if(debug)std::cout<<" MyPDF::MyPDF: Start overloaded constructor"<<std::endl;

    SetVariablesDefault();   
    gridName=_gridName;
    xscale=_xscale;
    steeringFilePath=_steeringFile;

    if(debug)
        std::cout<<"Setting up MyPDF with:<<<<<<<<<<<<<<<<<<<<<<<<<<"
                 <<"\n\tDebug: "<<(debug? "ON":"OFF")
                 <<"\n\tGridName: "<<gridName
                 <<"\n\tSteeringFile: "<<steeringFilePath
                 <<"\n\tOptionsFile: "<<optionsFileName
                 <<"\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<std::endl;

    if(FileExists(_steeringFile)==true) {
        ReadSteering(_steeringFile);
        Initialize();
    } else {
        std::cout<<"MyPDF::MyPDF: WARNING: Couldn't find file names: "<<_steeringFile<<std::endl;
        SetVariablesDefault();
    }

    if(debug)std::cout<<" MyPDF::MyPDF: End overloaded constructor"<<std::endl;
}


//overloaded constructor - Use if PDFErrorType is not declared in MyPDF steering file
MyPDF::MyPDF(string _gridName, double _xscale, bool _do_PDFBand, bool _do_AlphaS, bool _do_RenScale, bool _do_FacScale, bool _do_TotError, string _steeringFile, bool _debug)
{
    debug=_debug;
    if(debug)std::cout<<" MyPDF::MyPDF: Start overloaded constructor"<<std::endl;
    
    SetVariablesDefault();
    do_PDFBand = _do_PDFBand;
    do_AlphaS = _do_AlphaS;
    do_RenormalizationScale = _do_RenScale;
    do_FactorizationScale = _do_FacScale;
    do_TotError = _do_TotError;
    gridName = _gridName;
    xscale = _xscale;
    steeringFilePath = _steeringFile;

    if(debug)
        std::cout<<"Setting up MyPDF with:<<<<<<<<<<<<<<<<<<<<<<<<<<"
           <<"\n\tDebug: "<<(debug? "ON":"OFF")
           <<"\n\tGridName: "<<gridName
           <<"\n\tSteeringFile: "<<steeringFilePath
           <<"\n\tOptionsFile: "<<optionsFileName
           <<"\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<std::endl;

    if(FileExists(_steeringFile)==true) {
        ReadSteering(_steeringFile);
        //Initialize(); //User needs to call Initialize() explicitly after they have set all variables they want explictly
    } else {
        std::cout<<"MyPDF::MyPDF: WARNING: Couldn't find file names: "<<_steeringFile<<std::endl;
        SetVariablesDefault();
    }
    
    if(debug)std::cout<<" MyPDF::MyPDF: End overloaded constructor"<<std::endl;
}



//perform additional setup after all variables have ben set from steering and constructors
void MyPDF::Initialize()
{
    if (debug) std::cout<<" MyPDF::Initialize: Performing Initialization"<<std::endl;

    if(steeringFilePath.size()>0) SetSteeringFileNameAndDir(steeringFilePath);
    else std::cout<<" MyPDF::Initialize: ERROR: No steering file was provided!"<<std::endl;


    calc_desc = "mypdf-theory_errors";
    if(do_PDFBand)
        calc_desc+="_PDFBand";
    if(do_AlphaS)
        calc_desc+="_AlphaS";
    if(do_RenormalizationScale)
        calc_desc+="_RS";
    if(do_FactorizationScale)
        calc_desc+="_FS";
    if( !do_PDFBand && !do_AlphaS && !do_RenormalizationScale && !do_FactorizationScale ) {
        std::cout<<" MyPDF::Initialize: ERORR: All theory uncertainties disabled. Possible steering file error? Check file: "<<steeringFilePath<<std::endl;
        exit(0); //TEST
    }

    my_grid = new appl::grid(gridName.c_str());
    my_grid->trim();

    static const int nLoops    = 1;
    static const int nFlavours = 5;
    h_errors_PDFBand.clear();
    h_errors_AlphaS.clear();
    h_errors_RenormalizationScale.clear();
    h_errors_FactorizationScale.clear();

    double stev=1000.;
    TH1D* temp_hist;
    TH1D* temp_hist_prenorm;



    std::cout<<" MyPDF::Initialize: Fill PDF errors for PDFType: "<<PDFtype<<", PDFName: "<<PDFname<<", from PDFPath: "<<pdfSetPath<<std::endl;
    string default_pdf_set_name = (std::string) (pdfSetPath+"/"+PDFname+".LHgrid");
    if (PDFtype.compare("HERAPDF15NLO")==0) default_pdf_set_name =pdfSetPath+"/"+PDFtype+"_EIG.LHgrid"; //neededs the extra "_EIG"???
    std::cout<<" MyPDF::Initialize: init PDF set called: "<<default_pdf_set_name.c_str()<<std::endl;



    LHAPDF::initPDFSet(default_pdf_set_name.c_str(), 0);

    h_qqbar_prenorm = (TH1D*)my_grid->convolute_subproc(6, getPDF, alphasPDF, nLoops);      ///// maybe also subprocess 5?
    h_qqbar_prenorm->SetName((TString) ("h_qqbar_prenorm_" + calc_desc));

    TH1D* h_qqbar_prenorm2 = (TH1D*) my_grid->convolute_subproc(5, getPDF, alphasPDF, nLoops);
    h_qqbar_prenorm2->SetName((TString) ("h_qqbar_prenorm_" + calc_desc));
    h_qqbar_prenorm->Add(h_qqbar_prenorm2);
    h_qqbar_prenorm->SetLineColor(fillColorCode);
    h_qqbar_prenorm->SetMarkerColor(fillColorCode);

    h_qqbar = (TH1D*) TH1NormToTot(h_qqbar_prenorm, 1. / 1000., 1000.*xscale/stev);
    h_qqbar->SetName((TString) ("h_qqbar_" + calc_desc));
    h_qqbar->SetLineColor(fillColorCode);
    h_qqbar->SetMarkerColor(fillColorCode);

    h_gg_prenorm = (TH1D*) my_grid->convolute_subproc(0, getPDF, alphasPDF, nLoops);
    h_gg_prenorm->SetName((TString) ("h_gg_prenorm_" + calc_desc));
    h_gg_prenorm->SetLineColor(fillColorCode);
    h_gg_prenorm->SetMarkerColor(fillColorCode);

    h_gg = (TH1D*) TH1NormToTot(h_gg_prenorm, 1. / 1000., 1000.*xscale/stev);
    h_gg->SetName((TString) ("h_gg_" + calc_desc));
    h_gg->SetLineColor(fillColorCode);
    h_gg->SetMarkerColor(fillColorCode);

    h_tot_prenorm = (TH1D*) my_grid->convolute(getPDF, alphasPDF, nLoops);
    h_tot_prenorm->SetName((TString) ("h_tot_prenorm_" + calc_desc));
    h_tot_prenorm->SetLineColor(fillColorCode);
    h_tot_prenorm->SetMarkerColor(fillColorCode);

    h_tot = (TH1D*) TH1NormToTot(h_tot_prenorm, 1. / 1000., 1000.*xscale/stev);
    h_tot->SetName((TString) ("h_tot_" + calc_desc));
    h_tot->SetLineColor(fillColorCode);
    h_tot->SetMarkerColor(fillColorCode);

    h_gg_frac = (TH1D*) h_gg_prenorm->Clone((TString) ("h_gg_frac_" + calc_desc));
    h_gg_frac->Divide(h_tot_prenorm);

    h_qqbar_frac = (TH1D*) h_qqbar_prenorm->Clone((TString) ("h_qqbar_frac_" + calc_desc));
    h_qqbar_frac->Divide(h_tot_prenorm);


    if( do_RenormalizationScale ) {
        //check for used vars before continuing
        for(int icheck=0; icheck<n_SCALES; icheck++) {
            if(renScaleVals[icheck]==DEFAULT) {
                std::cout<<" MyPDF::Initialize: ERROR: A factorisation scale was not provided in steer file: "<<steeringFilePath<<std::endl;
                exit(0);
            }
        }
        
        for(int i_renScale = 0; i_renScale < n_SCALES; i_renScale++) {
            temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops, renScaleVals[i_renScale], 1.);
            temp_hist_prenorm->SetName((TString) ("h_xsec_rscale_" + renScaleNames[i_renScale]));
            temp_hist = (TH1D*) TH1NormToTot(temp_hist_prenorm, 1. / 1000., 1000.*xscale/stev);
            temp_hist->SetName((TString) ("h_xsec_rscale_" + renScaleNames[i_renScale] + "_norm"));
            h_errors_RenormalizationScale.push_back(temp_hist);
        }
    }

    if( do_FactorizationScale ) {
        //check for used vars before continuing
        for(int icheck=0; icheck<n_SCALES; icheck++) {
            if(facScaleVals[icheck]==DEFAULT) {
                std::cout<<" MyPDF::Initialize: ERROR: A factorisation scale was not provided in steer file: "<<steeringFilePath<<std::endl;
                exit(0);
            }
        }
        
        for(int i_facScale = 0; i_facScale < n_SCALES; i_facScale++) {
            temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops, 1., facScaleVals[i_facScale]);
            temp_hist_prenorm->SetName((TString) ("h_xsec_rscale_" + facScaleNames[i_facScale]));
            temp_hist = (TH1D*) TH1NormToTot(temp_hist_prenorm, 1. / 1000., 1000.*xscale/stev);
            temp_hist->SetName((TString) ("h_xsec_rscale_" + facScaleNames[i_facScale] + "_norm"));
            h_errors_FactorizationScale.push_back(temp_hist);
        }
    }

    /*
        if( do_AlphaS ) {
            LHAPDF::initPDFSet(default_pdf_set_name.c_str(), 0);
            temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
            temp_hist_prenorm->SetName((TString) ("h_xsec_default"));
            h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
            if(PDFtype.compare("CT10")==0) {
                LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/CT10as.LHgrid")).c_str(), 3);  //// alphaS down
                temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
                temp_hist_prenorm->SetName((TString) ("h_xsec_CT10as116_prenorm"));
                h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
                LHAPDF::initPDF(7);    /// alphaS up
                temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
                temp_hist_prenorm->SetName((TString) ("h_xsec_CT10as120_prenorm"));
                h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
            }
            else if(PDFtype.compare("MSTW2008nlo")==0) {
                LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/MSTW2008nlo68cl_asmz-68cl.LHgrid")).c_str(), 0);   //// alphaS down
                temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
                temp_hist_prenorm->SetName((TString) ("h_xsec_MSTW2008nloAsDown_prenorm"));
                h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
                LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/MSTW2008nlo68cl_asmz+68cl.LHgrid")).c_str(), 0);   //// alphaS up
                temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
                temp_hist_prenorm->SetName((TString) ("h_xsec_MSTW2008nloAsUp_prenorm"));
                h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
            }
            else if(PDFtype.compare("NNPDF23nlo")==0) {
                LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/NNPDF23_nlo_as_0116.LHgrid")).c_str(), 0);  //// alphaS down
                temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
                temp_hist_prenorm->SetName((TString) ("h_xsec_NNPDF23as117_prenorm"));
                h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
                LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/NNPDF23_nlo_as_0120.LHgrid")).c_str(), 0);  //// alphaS up
                temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
                temp_hist_prenorm->SetName((TString) ("h_xsec_NNPDF23as123_prenorm"));
                h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
            }
            else if(PDFtype.compare("HERAPDF15NLO")==0) {
                LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/HERAPDF15NLO_ALPHAS.LHgrid")).c_str(), 9);  //// alphaS down
                temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
                temp_hist_prenorm->SetName((TString) ("h_xsec_HERAPDF15NLOas1156_prenorm"));
                h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
                LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/HERAPDF15NLO_ALPHAS.LHgrid")).c_str(), 11);  //// alphaS up
                temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
                temp_hist_prenorm->SetName((TString) ("h_xsec_NNPDF23as1196_prenorm")); //looks wrong? copied from previous incorrectly??
                h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);
            }
            else {
                std::cout<<" MyPDF::Initialize: unsupported pdfCode '"<<PDFtype<<"' encountered."<<std::endl;
                exit(0);
            }
    */

    if( do_AlphaS ) {
        //check for necessary names before continuing
        if(AlphaSmemberNumDown==DEFAULT) {
            std::cout<<" MyPDF::Initialize: ERROR: 'AlphaSmemberNumDown' not found. Check steering: "<<steeringFilePath<<std::endl; exit(0);
        }
        if(AlphaSmemberNumUp==DEFAULT) {
            std::cout<<" MyPDF::Initialize: ERROR: 'AlphaSmemberNumUp' not found. Check steering: "<<steeringFilePath<<std::endl; exit(0);
        }
        if(AlphaSPDFSetNameUp.compare("")==0) {
            std::cout<<" MyPDF::Initialize: ERROR: 'AlphaSPDFSetNameUp' not found. Check steering: "<<steeringFilePath<<std::endl; exit(0);
        }
        if(AlphaSPDFSetNameDown.compare("")==0) {
            std::cout<<" MyPDF::Initialize: ERROR: 'AlphaSPDFSetNameDown' not found. Check steering: "<<steeringFilePath<<std::endl; exit(0);
        }
        if(AlphaSPDFSetHistNameUp.compare("")==0) {
            std::cout<<" MyPDF::Initialize: ERROR: 'AlphaSPDFSetHistNameUp' not found. Check steering: "<<steeringFilePath<<std::endl; exit(0);
        }
        if(AlphaSPDFSetHistNameDown.compare("")==0) {
            std::cout<<" MyPDF::Initialize: ERROR: 'AlphaSPDFSetHistNameDown' not found. Check steering: "<<steeringFilePath<<std::endl; exit(0);
        }

        //new alphaS set-up, check previously commented out section above for old implimentation
        LHAPDF::initPDFSet(default_pdf_set_name.c_str(), 0);
        temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
        temp_hist_prenorm->SetName((TString) ("h_xsec_default"));
        h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);

        LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/"+AlphaSPDFSetNameDown+".LHgrid")).c_str(), AlphaSmemberNumDown);  //// alphaS down
        temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
        temp_hist_prenorm->SetName((TString) ("h_xsec_"+AlphaSPDFSetHistNameDown));
        h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);

        if(PDFtype.compare("CT10")==0)
            LHAPDF::initPDF(7);    /// alphaS up
        else
            LHAPDF::initPDFSet(((std::string) (pdfSetPath+"/"+AlphaSPDFSetNameUp+".LHgrid")).c_str(), AlphaSmemberNumUp);  //// alphaS up
        temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
        temp_hist_prenorm->SetName((TString) ("h_xsec_"+AlphaSPDFSetHistNameUp));
        h_errors_AlphaS_prenorm.push_back(temp_hist_prenorm);


        for(int alphai = 0; alphai < h_errors_AlphaS_prenorm.size(); alphai++) {
            double stev=1000.;
            if(debug) std::cout<<" MyPDF::Initialize: use an xscale of "<<xscale<<" #evs assumed: "<<stev<<"\n";
            temp_hist = (TH1D*) TH1NormToTot(h_errors_AlphaS_prenorm.at(alphai), 1. / 1000., 1000.*xscale/stev);
            TString the_name = h_errors_AlphaS_prenorm.at(alphai)->GetName();
            temp_hist->SetName((TString) (the_name + "_normalized"));
            h_errors_AlphaS.push_back(temp_hist);
        }
    }

    if( do_PDFBand ) {
        LHAPDF::initPDFSet(default_pdf_set_name.c_str(), 0);
        //// Calculate PDF errors using standard PDF error band
        double stev;
        if(debug) std::cout<<" MyPDF::Initialize: Calc PDF errors"<<std::endl;
        for(int pdferri = 0; pdferri < n_PDFMembers; pdferri++) {
            if(debug) std::cout<<" MyPDF::Initialize: pdferri: "<<pdferri<<" of "<<n_PDFMembers<<std::endl;
            if(PDFtype.compare("HERAPDF15NLO")!=0) LHAPDF::initPDF(pdferri); //change "pdferri" to be fore name
            else {   //// e_HERAPDF15NLO
                if( pdferri <= 20 ) {
                    LHAPDF::initPDF(pdferri);
                } else if( pdferri == 21 ) {
                    LHAPDF::initPDFSet(pdfSetPath+"/HERAPDF15NLO_VAR.LHgrid", 0);
                } else if( pdferri > 21 ) {
                    LHAPDF::initPDF(pdferri - 21);
                }
            }
            TH1D* temp_hist_prenorm = (TH1D*) my_grid->convolute( getPDF, alphasPDF, nLoops);
            TH1D* temp_hist;
            TString this_pdf_err_code = "";
            this_pdf_err_code += pdferri;
            temp_hist_prenorm->SetName((TString) ("h_xsec_" + PDFname + "_" + this_pdf_err_code + "_prenorm"));
            h_errors_prenorm.push_back(temp_hist_prenorm);
            stev=1000.;
            TString the_name = temp_hist_prenorm->GetName();
            if(debug) std::cout<<" MyPDF::Initialize: use an xscale of "<<xscale<<" #evs assumed: "<<stev<<"\n";
            temp_hist = (TH1D*) TH1NormToTot(temp_hist_prenorm, 1. / 1000., 1000.*xscale/stev);
            temp_hist->SetName((TString) (the_name + "_normalized"));

            h_errors_PDFBand.push_back(temp_hist);
        }   /// pdf errors loop
    }  /// do_PDFBand
    if(debug) std::cout<<" MyPDF::Initialize: End of PDF errors loop"<<std::endl;
}



// divide histogram h1 scale by total number and scale bin width x
// assumes that histogram is divided by bin-width
TH1 *MyPDF::TH1NormToTot (TH1 *h1, double _yscale, double _xscale)
{
    if(debug) {
        std::cout<<" MyPDF::TH1NormToTot: Start"<<std::endl;
        std::cout<<" MyPDF::TH1NormToTot: renormalize with xscale = "<<_xscale<<", yscale = "<<_yscale<<std::endl;
    }

    if (!h1) std::cout<<" MyPDF::TH1NormToTot: TH1NormTotot h1 not found "<<std::endl;
    TH1D* h1new=(TH1D*) h1->Clone(h1->GetName());

    Double_t x, y, ey;
    Double_t sigtot=0.;
    for (Int_t i=1; i<=h1->GetNbinsX(); i++) {
        y=h1->GetBinContent(i)*_yscale;
        x=h1->GetBinWidth(i);
        sigtot+=y*x;
    }

    if (debug) std::cout<<" MyPDF::TH1NormToTot: sigtot= "<<sigtot<<std::endl;

    for (Int_t i=1; i<=h1->GetNbinsX(); i++) {
        x =h1->GetBinWidth(i);
        y =h1->GetBinContent(i)*_yscale*x;
        ey=h1->GetBinError(i)  *_yscale*x;
        x =h1->GetBinWidth(i)  *_xscale;
        if (x!=0) h1new->SetBinContent(i,y/x);
        else      h1new->SetBinContent(i,0.);
        if(debug)std::cout<<" MyPDF::TH1NormToTot: bin "<<i<<", center = "<<x<<", y-val = "<<y<<", fill with y/x = "<<h1new->GetBinContent(i)<<std::endl;

        if (x!=0) h1new->SetBinError(i,ey/x);
        else      h1new->SetBinError(i,0.);
    }

    if (sigtot!=0.)
        h1new->Scale(1./sigtot);
    if(debug) std::cout<<" MyPDF::TH1NormToTot: Scale by 1 / sigtot to get the final result "<<std::endl;

    float integral_so_far = 0.;
    for (Int_t i=1; i<=h1->GetNbinsX(); i++) {
        y =h1new->GetBinContent(i);
        x =h1new->GetBinWidth(i);
        integral_so_far += y*x;
        if(debug) cout<<" MyPDF::TH1NormToTot: i:"<<i<<" bincenter= "<<h1new->GetBinCenter(i)<<" Binw = "<<x<<" y= "<<y<<", integral so far: "<<integral_so_far<<std::endl;
    }

    if(debug) std::cout<<" MyPDF::TH1NormToTot: End"<<std::endl;
    return h1new;
}


void MyPDF::InitializeErrorGraphs()
{
    if(debug) std::cout<<" MyPDF::InitializeErrorGraphs: Start of InitializePDFErrorGraphs"<<std::endl;
    string this_name = "syst_band_" + PDFtype;

    if(debug) std::cout<<" MyPDF::InitializeErrorGraphs: PDFtype: "<<PDFtype<<", do_AlphaS: "<<do_AlphaS<<", do_PDFBand: "<<do_PDFBand<<std::endl;
    int n_bins = 0;
    Print();
    if(debug)  std::cout<<"MyPDF::InitializeErrorGraphs: Init PDFBand_results TGraphs. "
                            <<" AlphaS size: "<<h_errors_AlphaS.size()
                            <<", RenScale size: "<<h_errors_RenormalizationScale.size()
                            <<", FactScale size: "<<h_errors_FactorizationScale.size()
                            <<", PDFBand size: "<<h_errors_PDFBand.size()<<std::endl;
    if(debug)  std::cout<<"MyPDF::InitializeErrorGraphs: errors on: "
                            <<" AlphaS: "<<do_AlphaS
                            <<", RenScale: "<<do_RenormalizationScale
                            <<", FactScale: "<<do_FactorizationScale
                            <<", PDFBand: "<<do_PDFBand<<std::endl;

    if( do_AlphaS ) n_bins = h_errors_AlphaS.at(0)->GetNbinsX();
    else if( do_RenormalizationScale ) n_bins = h_errors_RenormalizationScale.at(0)->GetNbinsX();
    else if( do_FactorizationScale ) n_bins = h_errors_FactorizationScale.at(0)->GetNbinsX();
    else if( do_PDFBand ) n_bins = h_errors_PDFBand.at(0)->GetNbinsX();
    else {
        std::cout<<" MyPDF::InitializeErrorGraphs: ERORR: no theory uncertainties set. Configuration file error?"<<std::endl;
        exit(0); //TEST
    }

    double x_vals       [n_bins];
    double x_errs_low   [n_bins];
    double x_errs_high  [n_bins];
    double y_vals       [n_bins];
    double y_errs_low   [n_bins];
    double y_errs_high  [n_bins];


    if(debug) std::cout<<" MyPDF::InitializeErrorGraphs: Starting initialization"<<std::endl;
    float x_width;
    for(int pi = 0; pi < n_bins; pi++) {
        x_vals[pi] = 0;
        x_width = 0;
        if( do_AlphaS ) {
            x_vals[pi] = h_errors_AlphaS.at(0)->GetBinCenter(pi+1);
            x_width = h_errors_AlphaS.at(0)->GetBinWidth(pi+1);
            y_vals[pi] =  h_errors_AlphaS.at(0)->GetBinContent(pi+1);
        } else if( do_PDFBand ) {
            x_vals[pi] = h_errors_PDFBand.at(0)->GetBinCenter(pi+1);
            x_width = h_errors_PDFBand.at(0)->GetBinWidth(pi+1);
            y_vals[pi] = h_errors_PDFBand.at(0)->GetBinContent(pi+1);
        } else if( do_RenormalizationScale ) {
            x_vals[pi] = h_errors_RenormalizationScale.at(0)->GetBinCenter(pi+1);
            x_width = h_errors_RenormalizationScale.at(0)->GetBinWidth(pi+1);
            y_vals[pi] = h_errors_RenormalizationScale.at(0)->GetBinContent(pi+1);
        } else if( do_FactorizationScale ) {
            x_vals[pi] = h_errors_FactorizationScale.at(0)->GetBinCenter(pi+1);
            x_width = h_errors_FactorizationScale.at(0)->GetBinWidth(pi+1);
            y_vals[pi] = h_errors_FactorizationScale.at(0)->GetBinContent(pi+1);
        }

        x_errs_low[pi] = x_width / 2.;
        x_errs_high[pi] = x_width / 2.;
        y_errs_low[pi] = 0.;
        y_errs_high[pi] = 0.;
    }

    if(debug) std::cout<<" MyPDF::InitializeErrorGraphs: Making results graphs"<<std::endl;
    h_PDFBand_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
    h_PDFBand_results->SetName((TString) (this_name + "_PDFBand_results"));
    h_PDFBand_results->SetFillColor(fillColorCode);
    h_PDFBand_results->SetFillStyle(fillStyleCode);

    h_AlphaS_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
    h_AlphaS_results->SetName((TString) (this_name + "_AlphaS_results"));
    h_AlphaS_results->SetFillColor(fillColorCode);
    h_AlphaS_results->SetFillStyle(fillStyleCode);

    h_RenormalizationScale_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
    h_RenormalizationScale_results->SetName((TString) (this_name + "_RenScale_results"));
    h_RenormalizationScale_results->SetFillColor(fillColorCode);
    h_RenormalizationScale_results->SetFillStyle(fillStyleCode);

    h_FactorizationScale_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
    h_FactorizationScale_results->SetName((TString) (this_name + "_FacScale_results"));
    h_FactorizationScale_results->SetFillColor(fillColorCode);
    h_FactorizationScale_results->SetFillStyle(fillStyleCode);

    h_TotError_results = new TGraphAsymmErrors(n_bins, x_vals, y_vals, x_errs_low, x_errs_high, y_errs_low, y_errs_high);
    h_TotError_results->SetName((TString) (this_name + "_TotError_results"));
    h_TotError_results->SetFillColor(fillColorCode);
    h_TotError_results->SetFillStyle(fillStyleCode);

    if(debug) std::cout<<" MyPDF::InitializeErrorGraphs: End of InitializePDFErrorGraphs"<<std::endl;
}



void MyPDF::CalcSystErrors()
{
    std::cout<<" MyPDF::CalcSystErrors: Start syst error calc for: "<<PDFtype
             <<"\n\tPDFBand: "<<do_PDFBand
             <<", do_AlphaS: "<<do_AlphaS
             <<", do_RenScale: "<<do_RenormalizationScale
             <<", FacScale: "<<do_FactorizationScale<<std::endl;

    if( do_PDFBand ) CalcPDFBandErrors();
    if( do_AlphaS ) CalcAlphaSErrors();
    if( do_RenormalizationScale )  CalcRenormalizationScaleErrors();
    if( do_FactorizationScale )  CalcFactorizationScaleErrors();
    CalcTotErrors();

    if(debug) std::cout<<" MyPDF::CalcSystErrors: End syst error calc for: "<<PDFtype<<std::endl;
}


void MyPDF::CalcPDFBandErrors()
{
    if(debug) std::cout<<" MyPDF::CalcPDFBandErrors: Start calc of PDFBandErrors for: "<<PDFtype<<std::endl;

    double this_err_up         = 0.;
    double this_err_down       = 0.;
    double central_val         = 0.;    // needed for MSTW2008nlo and HERAPDF15NLO
    double average             = 0.;    // needed for NNPDF
    double extreme_pos_diff    = 0.;    // needed for HERAPDF15NLO
    double extreme_neg_diff    = 0.;    // needed for HERAPDF15NLO
    double diff_central        = 0.;    // needed for HERAPDF15NLO
    double mod_val             = 0.;    // needed for MSTW2008nlo


    for(int bi = 1; bi <= h_errors_PDFBand.at(0)->GetNbinsX(); bi++) {
        if(PDFtype.compare("NNPDF23nlo")==0) {
            for(int pdferri = 1; pdferri < (int) h_errors_PDFBand.size(); pdferri++) {
                average += h_errors_PDFBand.at(pdferri)->GetBinContent(bi);
            }
            average /= h_errors_PDFBand.size()-1;
            for(int pdferri = 1; pdferri < (int) h_errors_PDFBand.size(); pdferri++)  {
                this_err_up += pow(h_errors_PDFBand.at(pdferri)->GetBinContent(bi)-average, 2.);
            }
            this_err_up = TMath::Sqrt(this_err_up / (h_errors_PDFBand.size()-1));
            this_err_down = this_err_up;
            if( PDFErrorSize.compare("90Percent")==0 ) {
                this_err_up *= 1.645;
                this_err_down *= 1.645;
            }
        }
        else if(PDFtype.compare("CT10")==0 ) {
            for(int pdferri = 1; pdferri < (int) h_errors_PDFBand.size()-1; pdferri += 2) {
                this_err_up += pow( h_errors_PDFBand.at(pdferri)->GetBinContent(bi) - h_errors_PDFBand.at(pdferri+1)->GetBinContent(bi), 2.);
            }
            this_err_up = 0.5*TMath::Sqrt(this_err_up);
            if( PDFErrorSize.compare("OneSigma") ) this_err_up /= 1.645;
            this_err_down = this_err_up;
        }
        else if(PDFtype.compare("MSTW2008nlo")==0 ) {
            central_val = h_errors_PDFBand.at(0)->GetBinContent(bi);
            for(int pdferri = 1; pdferri < (int) h_errors_PDFBand.size(); pdferri ++) {
                mod_val = h_errors_PDFBand.at(pdferri)->GetBinContent(bi);
                if( mod_val > central_val ) this_err_up += pow(mod_val-central_val, 2.);
                else this_err_down += pow(central_val - mod_val, 2.);
            }
            this_err_down = TMath::Sqrt(this_err_down);
            this_err_up = TMath::Sqrt(this_err_up);
            if( PDFErrorSize.compare("90Percent")==0 ) {
                this_err_up *= 1.645;
                this_err_down *= 1.645;
            }
        }
        else if(PDFtype.compare("HERAPDF15NLO")==0) {
            central_val = h_errors_PDFBand.at(0)->GetBinContent(bi);
            for(int pdferri = 1; pdferri < 20; pdferri += 2) {   //// experimental errors
                this_err_up += pow( 0.5*(h_errors_PDFBand.at(pdferri)->GetBinContent(bi) - h_errors_PDFBand.at(pdferri+1)->GetBinContent(bi)), 2.);
            }
            this_err_down = this_err_up;
            for(int pdferri = 21; pdferri < 29; pdferri++) {   /// model errors
                if( h_errors_PDFBand.at(pdferri)->GetBinContent(bi) > central_val ) {
                    this_err_up  += pow( h_errors_PDFBand.at(pdferri)->GetBinContent(bi) - central_val, 2.);
                    this_err_down += 0.;
                } else {
                    this_err_up  += 0.;
                    this_err_down += pow( central_val - h_errors_PDFBand.at(pdferri)->GetBinContent(bi), 2.);
                }
            }

            for(int pdferri = 29; pdferri < 33; pdferri++) {     //// parameterization errors
                diff_central = h_errors_PDFBand.at(pdferri)->GetBinContent(bi) - central_val;
                if( diff_central > 0 && diff_central > extreme_pos_diff ) extreme_pos_diff = diff_central;
                if( diff_central < 0 && diff_central < extreme_neg_diff ) extreme_neg_diff = diff_central;
            }
            if( extreme_pos_diff > 0. ) this_err_up += pow(extreme_pos_diff, 2.);
            if( extreme_neg_diff < 0. ) this_err_down += pow(extreme_neg_diff, 2.);
            this_err_up = TMath::Sqrt(this_err_up);
            this_err_down = TMath::Sqrt(this_err_down);
            if( PDFErrorSize.compare("90Percent")==0 ) {
                this_err_up *= 1.645;
                this_err_down *= 1.645;
            }
        } else {
            std::cout<<" MyPDF::CalcPDFBandErrors: unsupported pdfCode encountered."<<std::endl;
            exit(0); //TEST
        }

        h_PDFBand_results->SetPointEYhigh(bi-1, this_err_up);
        h_PDFBand_results->SetPointEYlow(bi-1, this_err_down);
        Double_t x_val;
        Double_t y_val;
        h_PDFBand_results->GetPoint(bi-1, x_val, y_val);
    }  /// loop over bins

    if(debug) std::cout<<" MyPDF::CalcPDFBandErrors: End cal of PDFBandErrors for: "<<PDFtype<<std::endl;
}



void MyPDF::CalcAlphaSErrors()
{
    if(debug) std::cout<<" MyPDF::CalcAlphaSErrors: Starting calc of PDFAlphaSErrors for: "<<PDFtype<<std::endl;

    //assert(h_errors_AlphaS.size() == 3);
    double this_default_val = 0.;
    double this_err_down = 0.;
    double this_err_up = 0.;
    double error = 0.;

    for(int bi = 1; bi <= h_errors_AlphaS.at(0)->GetNbinsX(); bi++) {
        this_default_val = h_errors_AlphaS.at(0)->GetBinContent(bi);
        this_err_down = h_errors_AlphaS.at(1)->GetBinContent(bi);
        this_err_up = h_errors_AlphaS.at(2)->GetBinContent(bi);
        if(debug) std::cout<<" MyPDF::CalcAlphaSErrors: bi = "<<bi<<", default val = "<<this_default_val<<" +"<<this_err_up<<" -"<<this_err_down<<std::endl;

        error = 0.5*fabs(this_err_up-this_err_down);
        if( PDFErrorSize.compare("90Percent")==0 ) error *= 1.645;
        Double_t init_x_val;
        Double_t init_y_val;
        h_AlphaS_results->GetPoint(bi-1, init_x_val, init_y_val);
        h_AlphaS_results->SetPoint(bi-1, init_x_val, this_default_val);
        h_AlphaS_results->SetPointEYhigh(bi-1, error);
        h_AlphaS_results->SetPointEYlow(bi-1, error);

    } /// bi

    if(debug) std::cout<<" MyPDF::CalcAlphaSErrors: End calc of PDFAlphaSErrors for: "<<PDFtype<<std::endl;
}


void MyPDF::CalcRenormalizationScaleErrors()
{
    if(debug) std::cout<<" MyPDF::CalcRenormalizationScaleErrors: Start of calc of PDFRenormalizationScaleErrors for :"<<PDFtype<<std::endl;
    //assert(h_errors_RenormalizationScale.size() >= 3);

    for(int bi = 1; bi <= h_errors_RenormalizationScale.at(0)->GetNbinsX(); bi++) {
        double this_default_val = h_errors_RenormalizationScale.at(DEF)->GetBinContent(bi);
        double this_err_down = h_errors_RenormalizationScale.at(DOWN)->GetBinContent(bi);
        double this_err_up = h_errors_RenormalizationScale.at(UP)->GetBinContent(bi);
        if(debug)std::cout<<" MyPDF::CalcRenormalizationScaleErrors: bi = "<<bi<<", default val = "<<this_default_val<<" +"<<this_err_up  <<" -"<<this_err_down<<"\n";

        double error = 0.5*fabs(this_err_up-this_err_down);
        if( PDFErrorSize.compare("90Percent")==0 ) error *= 1.645;
        Double_t init_x_val;
        Double_t init_y_val;

        h_RenormalizationScale_results->GetPoint(bi-1, init_x_val, init_y_val);
        h_RenormalizationScale_results->SetPoint(bi-1, init_x_val, this_default_val);
        h_RenormalizationScale_results->SetPointEYhigh(bi-1, error);
        h_RenormalizationScale_results->SetPointEYlow(bi-1, error);

    } /// bi

    if(debug) std::cout<<" MyPDF::CalcRenormalizationScaleErrors: End of calc of PDFRenormalizationScaleErrors for: "<<PDFtype<<std::endl;
}


void MyPDF::CalcFactorizationScaleErrors()
{
    if(debug) std::cout<<" MyPDF::CalcFactorizationScaleErrors: Start of calc of PDFFactorizationScaleErrors for: "<<PDFtype<<std::endl;
    //assert(h_errors_FactorizationScale.size() >= 3);

    for(int bi = 1; bi <= h_errors_FactorizationScale.at(0)->GetNbinsX(); bi++) {
        double this_default_val = h_errors_FactorizationScale.at(DEF)->GetBinContent(bi);
        double this_err_down = h_errors_FactorizationScale.at(DOWN)->GetBinContent(bi);
        double this_err_up = h_errors_FactorizationScale.at(UP)->GetBinContent(bi);
        if(debug)std::cout<<" MyPDF::CalcFactorizationScaleErrors: bi = "<<bi<<", default val = "<<this_default_val<<" +"<<this_err_up  <<" -"<<this_err_down<<"\n";

        double error = 0.5*fabs(this_err_up-this_err_down);
        if( PDFErrorSize.compare("90Percent")==0 ) error *= 1.645;
        Double_t init_x_val;
        Double_t init_y_val;

        h_FactorizationScale_results->GetPoint(bi-1, init_x_val, init_y_val);
        h_FactorizationScale_results->SetPoint(bi-1, init_x_val, this_default_val);
        h_FactorizationScale_results->SetPointEYhigh(bi-1, error);
        h_FactorizationScale_results->SetPointEYlow(bi-1, error);

    } /// bi

    if(debug) std::cout<<" MyPDF::CalcFactorizationScaleErrors: End of calc for ThisFactorizationScalePDFSyst for: "<<PDFtype<<std::endl;
}


void MyPDF::CalcTotErrors()
{
    if(debug) std::cout<<" MyPDF::CalcTotErrors: Start of calc for total errors for: "<<PDFtype<<std::endl;

    double PDFBand_err_high = 0.;
    double PDFBand_err_low  = 0.;
    double AlphaS_err_high  = 0.;
    double AlphaS_err_low   = 0.;
    double RenormalizationScale_err_high    = 0.;
    double RenormalizationScale_err_low     = 0.;
    double FactorizationScale_err_high      = 0.;
    double FactorizationScale_err_low       = 0.;

    for(int pi = 0; pi < h_TotError_results->GetN(); pi++) {
        Double_t x_val=-999;
        Double_t y_val=-999;
        if( do_PDFBand ) {
            PDFBand_err_high = h_PDFBand_results->GetErrorYhigh(pi);
            PDFBand_err_low = h_PDFBand_results->GetErrorYlow(pi);
            h_PDFBand_results->GetPoint(pi, x_val, y_val);
        }
        if( do_AlphaS ) {
            AlphaS_err_high = h_AlphaS_results->GetErrorYhigh(pi);
            AlphaS_err_low = h_AlphaS_results->GetErrorYlow(pi);
            h_AlphaS_results->GetPoint(pi, x_val, y_val);
        }
        if( do_RenormalizationScale ) {
            RenormalizationScale_err_high = h_RenormalizationScale_results->GetErrorYhigh(pi);
            RenormalizationScale_err_low = h_RenormalizationScale_results->GetErrorYlow(pi);
            h_RenormalizationScale_results->GetPoint(pi, x_val, y_val);
        }
        if( do_FactorizationScale ) {
            FactorizationScale_err_high = h_FactorizationScale_results->GetErrorYhigh(pi);
            FactorizationScale_err_low = h_FactorizationScale_results->GetErrorYlow(pi);
            h_FactorizationScale_results->GetPoint(pi, x_val, y_val);
        }
        Double_t Tot_err_high = TMath::Sqrt(pow(PDFBand_err_high, 2.) + pow(AlphaS_err_high, 2.) + pow(RenormalizationScale_err_high, 2.) + pow(FactorizationScale_err_high, 2.));
        Double_t Tot_err_low = TMath::Sqrt(pow(PDFBand_err_low, 2.) + pow(AlphaS_err_low, 2.) + pow(RenormalizationScale_err_low, 2.) + pow(FactorizationScale_err_low, 2.));

        h_TotError_results->SetPoint(pi, x_val, y_val);
        h_TotError_results->SetPointEYhigh(pi, Tot_err_high);
        h_TotError_results->SetPointEYlow(pi, Tot_err_low);
    }  /// loop over points

    if(debug) std::cout<<" MyPDF::CalcTotErrors: End of calc for total errors for: "<<PDFtype<<std::endl;
}


void MyPDF::GetRatioToTH1(TH1D* href)
{
    if(debug) std::cout<<" MyPDF::GetRatioToTH1: start"<<std::endl;

    TGraphAsymmErrors* tgraph_href = TH1TOTGraphAsymm(href);
    TString ratio_to_ref_name = (TString) h_PDFBand_results->GetName() + "_ratio_to_ref";
    h_PDFBand_results_ratio_to_ref = MyTGraphErrorsDivide(h_PDFBand_results, tgraph_href);
    h_PDFBand_results_ratio_to_ref->SetName(ratio_to_ref_name);
    h_PDFBand_results_ratio_to_ref->SetFillColor(fillColorCode);
    h_PDFBand_results_ratio_to_ref->SetFillStyle(fillStyleCode);

    ratio_to_ref_name = (TString) h_AlphaS_results->GetName() + "_ratio_to_ref";
    h_AlphaS_results_ratio_to_ref = MyTGraphErrorsDivide(h_AlphaS_results, tgraph_href);
    h_AlphaS_results_ratio_to_ref->SetName(ratio_to_ref_name);
    h_AlphaS_results_ratio_to_ref->SetFillColor(fillColorCode);
    h_AlphaS_results_ratio_to_ref->SetFillStyle(fillStyleCode);

    ratio_to_ref_name = (TString) h_RenormalizationScale_results->GetName() + "_ratio_to_ref";
    h_RenormalizationScale_results_ratio_to_ref = MyTGraphErrorsDivide(h_RenormalizationScale_results, tgraph_href);
    h_RenormalizationScale_results_ratio_to_ref->SetName(ratio_to_ref_name);
    h_RenormalizationScale_results_ratio_to_ref->SetFillColor(fillColorCode);
    h_RenormalizationScale_results_ratio_to_ref->SetFillStyle(fillStyleCode);

    ratio_to_ref_name = (TString) h_FactorizationScale_results->GetName() + "_ratio_to_ref";
    h_FactorizationScale_results_ratio_to_ref = MyTGraphErrorsDivide(h_FactorizationScale_results, tgraph_href);
    h_FactorizationScale_results_ratio_to_ref->SetName(ratio_to_ref_name);
    h_FactorizationScale_results_ratio_to_ref->SetFillColor(fillColorCode);
    h_FactorizationScale_results_ratio_to_ref->SetFillStyle(fillStyleCode);

    ratio_to_ref_name = (TString) h_TotError_results->GetName() + "_ratio_to_ref";
    h_TotError_results_ratio_to_ref = MyTGraphErrorsDivide(h_TotError_results, tgraph_href);
    h_TotError_results_ratio_to_ref->SetName(ratio_to_ref_name);
    h_TotError_results_ratio_to_ref->SetFillColor(fillColorCode);
    h_TotError_results_ratio_to_ref->SetFillStyle(fillStyleCode);

    if(debug) std::cout<<" MyPDF::GetRatioToTH1: End"<<std::endl;
}


//convert the histogram h1 into a graph
TGraphAsymmErrors* MyPDF::TH1TOTGraphAsymm(TH1 *h1)
{
    if(debug) std::cout<<" MyPDF::TH1TOTGraphAsymm: start"<<std::endl;

    if (!h1) {
        std::cout<<" MyPDF::TH1TOTGraphAsymm: histogram not found !"<<std::endl;
        exit(0); //TEST
    }
    TGraphAsymmErrors* g1= new TGraphAsymmErrors();

    Double_t x, y, ex, ey;
    for (Int_t i=0; i<h1->GetNbinsX(); i++) {
        y=h1->GetBinContent(i+1);
        ey=h1->GetBinError(i+1);
        x=h1->GetBinCenter(i+1);
        ex=h1->GetBinWidth(i+1)/2.;
        //if(debug) std::cout<<i<<" x,y = "<<x<<" "<<y<<" ex,ey = "<<ex<<" "<<ey<<std::endl;

        g1->SetPoint(i,x,y);
        g1->SetPointError(i,ex,ex,ey,ey);
    }

    if(debug) std::cout<<" MyPDF::TH1TOTGraphAsymm: End"<<std::endl;
    return g1;
}


// Divide two TGraphAsymmErrors: new=g1/g2
//
// noerr=0: put all errors to zero
//       1: add errors from two graph quadrically
//       2: set errors from graph 2 to zero
TGraphAsymmErrors* MyPDF::MyTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2, Int_t noerr) {
    if(debug) std::cout<<" MyPDF::MyTGraphErrorsDivide: started"<<std::endl;
    if (!g1) std::cout<<" MyPDF::MyTGraphErrorsDivide: g1 does not exist ! "<<std::endl;
    if (!g2) std::cout<<" MyPDF::MyTGraphErrorsDivide: g2 does not exist ! "<<std::endl;

    Int_t n1=g1->GetN();
    Int_t n2=g2->GetN();

    if (n1!=n2) {
        std::cout<<" MyPDF::MyTGraphErrorsDivide: vector do not have the same number of entries!"
                 <<"\n\tg1: "<<g1->GetName()<<" n1= "<<n1
                 <<"\n\tg2: "<<g2->GetName()<<" n2= "<<n2<<std::endl;
    }

    TGraphAsymmErrors* g3= new TGraphAsymmErrors();
    if (!g3) std::cout<<" MyPDF::MyTGraphErrorsDivide: problem to make new vector ! "<<std::endl;
    g3->SetName       (g1->GetName());
    g3->SetMarkerStyle(g1->GetMarkerStyle());
    g3->SetMarkerColor(g1->GetMarkerColor());
    g3->SetLineColor  (g1->GetLineColor());

    Double_t   x1=0.,   y1=0., x2=0., y2=0.;
    Double_t dx1h=0., dx1l=0.;
    Double_t dy1h=0., dy1l=0.;
    Double_t dy2h=0., dy2l=0.;

    Double_t* EXhigh1 = g1->GetEXhigh();
    Double_t* EXlow1 =  g1->GetEXlow();
    Double_t* EYhigh1 = g1->GetEYhigh();
    Double_t* EYlow1 =  g1->GetEYlow();

    Double_t* EXhigh2 = g2->GetEXhigh();
    Double_t* EXlow2 =  g2->GetEXlow();
    Double_t* EYhigh2 = g2->GetEYhigh();
    Double_t* EYlow2 =  g2->GetEYlow();

    Int_t iv=0;
    for (Int_t i1=0; i1<n1; i1++) {  //loop over point of graph1
        Int_t matchcount=0;
        for (Int_t i2=0; i2<n2; i2++) {//loop over point of graph2
            g1->GetPoint(i1,x1,y1);
            g2->GetPoint(i2,x2,y2);
            Double_t emean=(EXhigh1[i1]+EXhigh2[i2]+EXlow1[i1]+EXlow2[i2])/4.;
            if (fabs(x1-x2)>emean) {
                //std::cout<<" MyPDF::MyTGraphErrorsDivide: x1 and x2 not the same x1= "<<x1<<" x2= "<<x2<<std::endl;
            } else { // do something only if x1=x2
                matchcount++;
                //std::cout<<" MyPDF::MyTGraphErrorsDivide: x1 and x2 match x1= "<<x1<<" x2= "<<x2<<std::endl;
                dx1h  = EXhigh1[i1];
                dx1l  = EXlow1[i1];
                if (y1!=0.) dy1h  = EYhigh1[i1]/y1;
                else        dy1h  = 0.;
                if (y2!=0.) dy2h  = EYhigh2[i2]/y2;
                else        dy2h  = 0.;
                if (y1!=0.) dy1l  = EYlow1 [i1]/y1;
                else        dy1l  = 0.;
                if (y2!=0.) dy2l  = EYlow2 [i2]/y2;
                else        dy2l  = 0.;

                if (debug) {
                    std::cout<<"MyPDF::MyTGraphErrorsDivide: "
                             <<"\n\ti1: "    <<i1<<", i2: "<<i2
                             <<"\n\tdy1l: "  <<dy1l<<", dy1h: "<<dy1h
                             <<"\n\tdy2l: "  <<dy2l<<", dy2h "<<dy2h
                             <<"\n\tsqrt: " <<sqrt(dy1l*dy1l+dy2l*dy2l)<<", "<<sqrt(dy1h*dy1h+dy2h*dy2h)<<std::endl;
                }

                if (y2!=0.) g3->SetPoint(iv, x1,y1/y2);
                else        g3->SetPoint(iv, x1,y2);
                Double_t el=0.;
                Double_t eh=0.;

                if (noerr==2) {
                    dy2l=0.;
                    dy2h=0.;
                }
                if (y1!=0. && y2!=0.) el=sqrt(dy1l*dy1l+dy2l*dy2l)*(y1/y2);
                if (y1!=0. && y2!=0.) eh=sqrt(dy1h*dy1h+dy2h*dy2h)*(y1/y2);


                if (debug) std::cout<<"\tdx1h: "<<dx1h<<", dx1l: "<<dx1l<<", el: "<<el<<", eh: "<<eh<<std::endl;

                if (noerr==0)   g3->SetPointError(iv,dx1l,dx1h,0,0);
                else            g3->SetPointError(iv,dx1l,dx1h,el,eh);

                iv++;
            }
        }
        if (matchcount>1) {
            std::cout<<" MyPDF::MyTGraphErrorsDivide: too many x-points matched ! "<<std::endl;
            exit (1);
        }
    }

    if(debug) std::cout<<" MyPDF::MyTGraphErrorsDivide: End"<<std::endl;
    return g3;
}


void MyPDF::CalcChi2(TGraphAsymmErrors *g_theory, TGraphAsymmErrors *g_data, TMatrixT<double> data_cov_matrix, double &chi2)
{
    // https://cds.cern.ch/record/1470588/files/ATL-COM-PHYS-2012-1137.pdf
    if(debug) std::cout<<" MyPDF::CalcChi2: CalcChi2"<<std::endl;
    //assert( g_theory->GetN() == g_data->GetN() );
    //assert( g_theory->GetN() == data_cov_matrix.GetNcols() );
    //assert( g_theory->GetN() == data_cov_matrix.GetNrows() );

    /// Fill in the theory covariance matrix
    TMatrixT<double> theory_cov_matrix(g_theory->GetN(), g_theory->GetN());
    for(int pi = 0; pi < g_theory->GetN(); pi++) {
        for(int pi2 = 0; pi2 < g_theory->GetN(); pi2++) {
            if( pi != pi2 ) theory_cov_matrix(pi, pi2) = 0;
            if( pi == pi2 ) {
                double theory_uncertainty = 0.5*(g_theory->GetErrorYhigh(pi) + g_theory->GetErrorYlow(pi));
                theory_cov_matrix(pi, pi2) = theory_uncertainty*theory_uncertainty;
            }
        }
    }

    if(debug)std::cout<<" MyPDF::CalcChi2: At Start, dump contents of data cov matrix: "<<std::endl;
    for(int pi = 0; pi < g_theory->GetN(); pi++) {
        for(int pi2 = 0; pi2 < g_theory->GetN(); pi2++) {
            std::cout<<data_cov_matrix(pi,pi2)<<"\t";
        }
        std::cout<<"\n"; //formatting
    }

    TMatrixT<double> tot_cov_matrix = theory_cov_matrix + data_cov_matrix;
    if(debug)std::cout<<" MyPDF::CalcChi2: After adding theory, dump contents of cov matrix: "<<std::endl;
    for(int pi = 0; pi < g_theory->GetN(); pi++) {
        for(int pi2 = 0; pi2 < g_theory->GetN(); pi2++) {
            std::cout<<tot_cov_matrix(pi,pi2)<<"\t";
        }
        std::cout<<"\n"; //formatting
    }

    TMatrixT<double> invertex_cov_matrix = tot_cov_matrix.Invert();    //// Now it includes the theory errors in the diagonal elements ...
    if(debug)std::cout<<" MyPDF::CalcChi2: After inversion, dump contents of cov matrix: "<<std::endl;
    for(int pi = 0; pi < g_theory->GetN(); pi++) {
        for(int pi2 = 0; pi2 < g_theory->GetN(); pi2++) {
            std::cout<<invertex_cov_matrix(pi,pi2)<<"\t";
        }
        std::cout<<"\n"; //formatting
    }

    // Loop over bins and determine data-theory matrices
    TMatrixT<double> row_data_minus_theory(1, g_theory->GetN());
    TMatrixT<double> col_data_minus_theory(g_theory->GetN(), 1);
    for(int pi = 0; pi < g_theory->GetN(); pi++) {
        Double_t data_val;
        Double_t theory_val;
        Double_t x_val;

        g_theory->GetPoint(pi, x_val, theory_val);
        g_data->GetPoint(pi, x_val, data_val);
        row_data_minus_theory(0,pi) = data_val - theory_val;
        col_data_minus_theory(pi,0) = data_val - theory_val;
        if(debug)std::cout<<" MyPDF::CalcChi2: At "<<x_val<<", data = "<<data_val<<", theory = "<<theory_val<<", content = "<<row_data_minus_theory(0,pi)<<std::endl;
    }  // pi

    TMatrixT<double> cov_times_col = invertex_cov_matrix*col_data_minus_theory;
    //assert( cov_times_col.GetNrows() == g_theory->GetN());
    //assert( cov_times_col.GetNcols() == 1);
    if(debug)std::cout<<" MyPDF::CalcChi2: After first multiplication matrix is:"<<std::endl;
    for(int pi = 0; pi < g_theory->GetN(); pi++) {
        std::cout<<cov_times_col(pi, 0)<<"\t";
    }
    std::cout<<"\n";
    TMatrixT<double> result = row_data_minus_theory*cov_times_col;
    //assert( result.GetNrows() == 1);
    //assert( result.GetNcols() == 1);

    chi2 = result(0,0);
    if(debug)std::cout<<" MyPDF::CalcChi2: End chi2 = "<<chi2<<std::endl;
}



//read the provided steering file and set internal variables depending on what is read
void MyPDF::ReadSteering(const string _fileName)
{
    string fName="";
    if(_fileName.size()>0)  fName=_fileName;
    else                    fName=steeringFilePath;

    if (debug) std::cout<<" MyPDF::ReadSteering: reading steering file named: "<<fName<<std::endl;

    //Open the file for reading if it can be read/found
    ifstream infile(fName.c_str(), ios::in);
    if(!infile) {
        cerr<<" MyPDF::ReadSteering: WARNING: Can't open "<<fName<<std::endl;
        infile.close();
        exit (1); //TEMP TEST
    } else {
        if (debug) std::cout<<" MyPDF::ReadSteering: Steering file named successfuly opened."<<std::endl;
    }


    pdfSetPath = deafultPDFSetPath;
    string pdfSetDefaultPath = GetEnv("LHAPATH");

    if(pdfSetDefaultPath.size()>0) {
        if(pdfSetDefaultPath.find_last_of("/") == pdfSetDefaultPath.size()-1)
            pdfSetDefaultPath = pdfSetDefaultPath.substr(0,pdfSetDefaultPath.size()-1); //remove trailing slashes if there are any

        pdfSetPath=pdfSetDefaultPath;

        std::cout<<" makegridfromsherpa::main: LHAPATH environment varaible found, using path: "<<pdfSetPath<<std::endl;
    } else {
        std::cout<<" makegridfromsherpa::main: LHAPATH environment varaible not set, using default: "<<pdfSetPath<<std::endl;
    }


    string line;        //read a whole line from steering file and save in this buffer
    string optionName;  //holds the name of the option read from the line buffer
    string text;        //holds remaining text after option name

    //load in all valid options that could show up in the steering file
    myOptions= new OptionHandler(optionsFileName);
    //mycsoptions->generateResultFile(); //generates *.txt file showing results of OptionHandler::isKnownOption execution
    int w=20; //arbitrary width number for printing nicely formatted output

    //read and set all options and data from the steering file
    while (infile.good()) {
        getline(infile, line); //read a line from steering to the buffer

        optionName=line.substr(0, line.find(' ')); //extract option
        text=line.substr(line.find(' ')+1,line.size()); //'text' after the steering option could be broken up further if desired

        if(debug) {
            std::cout<<"\n MyPDF::ReadSteering: Read in:<<<<<<<<<<<<<<<<<<<"
                     <<"\n"<<setw(w)<<"line:"<<line
                     <<"\n"<<setw(w)<<"optionName:"<<optionName
                     <<"\n"<<setw(w)<<"text:"<<text
                     <<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
        }
        if(line[0] != '%') { //ignore comments
            if(myOptions->isKnownOption(optionName)==false) {
                if(debug) std::cout<<" MyPDF::ReadSteering: Found unsupported option: '"<<optionName<<"'"<<std::endl;
            } else if (optionName.compare("debug")==0) {
                debug=true;
            } else if (optionName.compare("PDFtype")==0) {
                PDFtype=text;
            } else if (optionName.compare("PDFname")==0) {
                PDFname=text;
            } else if (optionName.compare("numPDFMembers")==0) {
                sscanf(text.c_str(), "%d", &n_PDFMembers);
            } else if (optionName.compare("fillStyleCode")==0) {
                sscanf(text.c_str(), "%d", &fillStyleCode);
            } else if (optionName.compare("fillColorCode")==0) {
                sscanf(text.c_str(), "%d", &fillColorCode);
            } else if (optionName.compare("PDFBandType")==0) {
                PDFBandType=text;
            } else if (optionName.compare("PDFErrorType")==0) {
                PDFErrorType=text;
                //possible error types: PDFBand, AlphaS, RenormalizationScale, FactorizationScale, TotError
                if(PDFErrorType.compare("PDFBand")==0) {
                    do_PDFBand=true;
                } else if(PDFErrorType.compare("AlphaS")==0) {
                    do_AlphaS=true;
                } else if(PDFErrorType.compare("RenormalizationScale")==0) {
                    do_RenormalizationScale=true;
                } else if(PDFErrorType.compare("FactorizationScale")==0) {
                    do_FactorizationScale=true;
                } else if(PDFErrorType.compare("TotError")==0) {
                    do_TotError=true;
                }
            } else if (optionName.compare("PDFErrorSize")==0) {
                PDFErrorSize=text;
            } else if (optionName.compare("renScaleValUp")==0) {
                sscanf(text.c_str(), "%lf", &renScaleValUp);
                SetRenScaleValUp(renScaleValUp);
            } else if (optionName.compare("renScaleValDefault")==0) {
                sscanf(text.c_str(), "%lf", &renScaleValDefault);
                SetRenScaleValDefault(renScaleValDefault);
            } else if (optionName.compare("renScaleValDown")==0) {
                sscanf(text.c_str(), "%lf", &renScaleValDown);
                SetRenScaleValDown(renScaleValDown);
            } else if (optionName.compare("facScaleValUp")==0) {
                sscanf(text.c_str(), "%lf", &facScaleValUp);
                SetFacScaleValUp(facScaleValUp);
            } else if (optionName.compare("facScaleValDefault")==0) {
                sscanf(text.c_str(), "%lf", &facScaleValDefault);
                SetFacScaleValDefault(facScaleValDefault);
            } else if (optionName.compare("facScaleValDown")==0) {
                sscanf(text.c_str(), "%lf", &facScaleValDown);
                SetFacScaleValUp(facScaleValDown);
            } else if (optionName.compare("pdfSetPath")==0) {
                pdfSetPath=text;
                if(pdfSetPath.find_last_of("/") == pdfSetPath.size()-1)
                    pdfSetPath = pdfSetPath.substr(0,pdfSetPath.size()-1); //remove trailing slashes if there are any
                std::cout<<" MyPDF::ReadSteering: new PDF set path found! Now path is: "<<pdfSetPath<<std::endl;
            } else if (optionName.compare("AlphaSmemberNumDown")==0) {
                sscanf(text.c_str(), "%d", &AlphaSmemberNumDown);
            } else if (optionName.compare("AlphaSmemberNumUp")==0) {
                sscanf(text.c_str(), "%d", &AlphaSmemberNumUp);
            } else if (optionName.compare("AlphaSPDFSetNameDown")==0) {
                AlphaSPDFSetNameDown = text;
            } else if (optionName.compare("AlphaSPDFSetNameUp")==0) {
                AlphaSPDFSetNameUp = text;
            } else if (optionName.compare("AlphaSPDFSetHistNameDown")==0) {
                AlphaSPDFSetHistNameDown = text;
            } else if (optionName.compare("AlphaSPDFSetHistNameUp")==0) {
                AlphaSPDFSetHistNameUp = text;
            }
        }
    }
}


//Print all relevant internal variable values
void MyPDF::Print(string level)
{
    int w=30;               //arbitrary size that makes the formatting look pretty
    string empty="<empty>"; //print this if no input has been provided for that variable
    string ON="ON";         //bool true
    string OFF="OFF";       //bool false

    std::cout<<" MyPDF::Print: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
             <<"\n"<<setw(w)<<"debug:"                <<setw(w)<<(debug? "ON":"OFF")
             <<"\n"<<setw(w)<<"steeringFilePath:"     <<setw(w)<<(steeringFilePath.size()>0? steeringFilePath:empty)
             <<"\n"<<setw(w)<<"steeringFileDir:"      <<setw(w)<<(steeringFileDir.size()>0? steeringFileDir:empty)
             <<"\n"<<setw(w)<<"steeringFileName:"     <<setw(w)<<(steeringFileName.size()>0? steeringFileName:empty)
             <<"\n"<<setw(w)<<"optionsFile:"          <<setw(w)<<(optionsFileName.size()>0? optionsFileName:empty)
             <<"\n"<<setw(w)<<"gridName:"             <<setw(w)<<(gridName.size()>0? gridName:empty)
             <<std::endl;
             
    if(level.compare("all")==0){ //print all internal variables if requested
        std::cout<<"\n"<<setw(w)<<"PDFtype:"              <<setw(w)<<(PDFtype.size()>0? PDFtype:empty)
             <<"\n"<<setw(w)<<"PDFname:"              <<setw(w)<<(PDFname.size()>0? PDFname:empty)
             <<"\n"<<setw(w)<<"numPDFMembers:"        <<setw(w)<<(n_PDFMembers!=DEFAULT? to_string(n_PDFMembers):empty)
             <<"\n"<<setw(w)<<"fillStyleCode:"        <<setw(w)<<(fillStyleCode!=DEFAULT? to_string(fillStyleCode):empty)
             <<"\n"<<setw(w)<<"fillColorCode:"        <<setw(w)<<(fillColorCode!=DEFAULT? to_string(fillColorCode):empty)
             <<"\n"<<setw(w)<<"PDFBandType:"          <<setw(w)<<(PDFBandType.size()>0? PDFBandType:empty)
             <<"\n"<<setw(w)<<"**PDF ERROR TYPE(s) ACTIVE**"
             <<"\n"<<setw(w)<<"PDFBand:"              <<setw(w)<<(do_PDFBand? ON:OFF)
             <<"\n"<<setw(w)<<"AlphaS:"               <<setw(w)<<(do_AlphaS? ON:OFF)
             <<"\n"<<setw(w)<<"RenScale:"             <<setw(w)<<(do_RenormalizationScale? ON:OFF)
             <<"\n"<<setw(w)<<"FactScale:"            <<setw(w)<<(do_FactorizationScale? ON:OFF)
             <<"\n"<<setw(w)<<"TotError:"             <<setw(w)<<(do_TotError? ON:OFF)
             <<"\n"<<setw(w)<<"PDFErrorSize:"         <<setw(w)<<(PDFErrorSize.size()>0? PDFErrorSize:empty)
             <<"\n"<<setw(w)<<"**REN SCALE VALS**"
             <<"\n"<<setw(w)<<"renScaleVal-Up:"       <<setw(w)<<(renScaleValUp!=DEFAULT? to_string(renScaleValUp):empty)
             <<"\n"<<setw(w)<<"renScaleVal-Default:"  <<setw(w)<<(renScaleValDefault!=DEFAULT? to_string(renScaleValDefault):empty)
             <<"\n"<<setw(w)<<"renScaleVal-Down:"     <<setw(w)<<(renScaleValDown!=DEFAULT? to_string(renScaleValDown):empty)
             <<"\n"<<setw(w)<<"**FAC SCALE VALS**"
             <<"\n"<<setw(w)<<"facScaleVals-Up:"      <<setw(w)<<(renScaleValUp!=DEFAULT? to_string(renScaleValUp):empty)
             <<"\n"<<setw(w)<<"facScaleVals-Default:" <<setw(w)<<(renScaleValDefault!=DEFAULT? to_string(renScaleValDefault):empty)
             <<"\n"<<setw(w)<<"facScaleVals-Down:"    <<setw(w)<<(renScaleValDown!=DEFAULT? to_string(renScaleValDown):empty)<<std::endl;
    }
    
    std::cout<<"MyPDF::Print:<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<std::endl;
}


//always check for file existence before usage
bool MyPDF::FileExists(const string _fileName)
{
    bool exists;

    if ( FILE* file=fopen(_fileName.c_str(),"r") ) {
        fclose(file); 
        return true;
    }
    else return false;
}


//give default vals for ALL variables to provide predictable behavior even if accidental use before initialisation.
//allows for proper check for expected setup deafults before doing anything
void MyPDF::SetVariablesDefault()
{
    if(debug) std::cout<<" MyPDF::setVariablesDefault: Start default values being set."<<std::endl;
    string defaultString="";
    //int defaultInt=-1;

    //debug=false;
    //optionsFileName=defaultOptionsFileName;
    optionsFileName="options_mypdf.txt"; // hardcoded, should use some default string
    myOptions=NULL;
    steeringFileDir=defaultString;
    steeringFileName=defaultString;
    PDFtype=defaultString;
    PDFname=defaultString;
    n_PDFMembers=DEFAULT;
    fillStyleCode=DEFAULT;
    fillColorCode=DEFAULT;
    PDFBandType=defaultString;
    PDFErrorType=defaultString;
    PDFErrorSize=defaultString;

    renScaleNameUp=defaultString;
    renScaleNameDefault=defaultString;
    renScaleNameDown=defaultString;
    
    renScaleValUp=DEFAULT;
    renScaleValDefault=DEFAULT;
    renScaleValDown=DEFAULT;
    
    facScaleNameUp=defaultString;
    facScaleNameDefault=defaultString;
    facScaleNameDown=defaultString;
    
    facScaleValUp=DEFAULT;
    facScaleValDefault=DEFAULT;
    facScaleValDown=DEFAULT;
    
    for(int i=0; i<n_SCALES; i++) {
        renScaleNames[i] = defaultString;
        facScaleNames[i] = defaultString;
        renScaleVals[i] = DEFAULT;
        renScaleVals[i] = DEFAULT;
    }

    do_PDFBand=false;
    do_AlphaS=false;
    do_RenormalizationScale=false;
    do_FactorizationScale=false;
    do_TotError=false;

    AlphaSmemberNumDown=DEFAULT;
    AlphaSmemberNumUp=DEFAULT;
    AlphaSPDFSetNameDown=defaultString;
    
    AlphaSPDFSetNameUp=defaultString;
    AlphaSPDFSetHistNameDown=defaultString;
    AlphaSPDFSetHistNameUp=defaultString;

    if(debug) std::cout<<" MyPDF::setVariablesDefault: End default values are set."<<std::endl;
}


//determine the steering file name and directory from the provided path
void MyPDF::SetSteeringFileNameAndDir(const string _path)
{
    if(debug) std::cout<<" MyPDF::setSteeringFileNameAndDir: Start extracting file name and directory name from path: "<<_path<<std::endl;
    int pathLength=_path.length();

    if(pathLength>0) {
        int found=_path.find_last_of("/");
        if(found==-1) {
            steeringFileDir="<current dir>"; //if there are no slashes, then the file is in the current directory
            steeringFileName=_path;
        } else {
            steeringFileDir =_path.substr(0,found);
            steeringFileName=_path.substr(found+1);
        }

        if(debug) {
            std::cout<<" MyPDF::setSteeringFileNameAndPath:"
                     <<" \n\tSplitting: '"<<_path<<"'"
                     <<" \n\tpath: '"<<_path.substr(0,found)<<"'"
                     <<" \n\tfile: '"<<_path.substr(found+1)<<"'"<<std::endl;
        }
    } else {
        std::cout<<" MyPDF::setSteeringFileNameAndDir: Steering file path '"<<_path<<"' is invalid."<<std::endl;
    }
}


//Print out all known supported options that COULD be read from the steeting file
void MyPDF::PrintKnownOptionsForSteering()
{
    vector<string> *knownOptions;
    int w=25;               //arbitrary size that makes the formatting look pretty

    if(myOptions!=NULL) {
        knownOptions = myOptions->getKnownOptions();
        std::cout<<" MyPDF::PrintKnownOptionsForSteering: Start printing "<<knownOptions->size()<<" known options..."<<std::endl;
        for(int i=0; i<knownOptions->size(); i++)
            std::cout<<"\t("<<i<<"):"<<setw(w)<<knownOptions->at(i)<<std::endl;
        std::cout<<" MyPDF::PrintKnownOptionsForSteering: End printing known options."<<std::endl;
    } else {
        std::cout<<" MyPDF::PrintKnownOptionsForSteering: ERROR: Steering file has not yet been read."<<std::endl;
    }
}

//Print out all found options that HAVE been read from the steeting file, if they have been read
void MyPDF::PrintFoundOptionsFromSteering()
{
    vector<string> *foundOptions;
    int w=25;               //arbitrary size that makes the formatting look pretty

    if(myOptions!=NULL) {
        foundOptions = myOptions->getFoundOptions();
        std::cout<<" MyPDF::PrintFoundOptionsFromSteering: Start printing "<<foundOptions->size()<<" found options..."<<std::endl;
        for(int i=0; i<foundOptions->size(); i++)
            std::cout<<"\t("<<i<<"):"<<setw(w)<<foundOptions->at(i)<<std::endl;
        std::cout<<" MyPDF::PrintFoundOptionsFromSteering: End printing found options."<<std::endl;
    } else {
        std::cout<<" MyPDF::PrintFoundOptionsFromSteering: ERROR: Steering file has not yet been read."<<std::endl;
    }
}


//mutator methods
void MyPDF::SetDebug(bool _debug) {
    debug=_debug;
}
void MyPDF::SetGridName(string _gridName) {
    gridName=_gridName;
}
void MyPDF::SetSteeringFilePath(string _steeringFilePath) {
    steeringFilePath=_steeringFilePath;
    //udate the Dir location and file name of the steering file if the path is changed
    SetSteeringFileNameAndDir(steeringFilePath);
}
void MyPDF::SetSteeringFileDir(string _steeringFileDir) {
    steeringFileDir=_steeringFileDir;
}
void MyPDF::SetSteeringFileName(string _steeringFileName) {
    steeringFileName=_steeringFileName;
}
void MyPDF::SetPDFtype(string _PDFtype) {
    PDFtype=_PDFtype;
}
void MyPDF::SetPDFname(string _PDFname) {
    PDFname=_PDFname;
}
void MyPDF::SetNumPDFMembers(int _n_PDFMembers) {
    n_PDFMembers=_n_PDFMembers;
}
void MyPDF::SetFillStyleCode(int _fillStyleCode) {
    fillStyleCode=_fillStyleCode;
}
void MyPDF::SetFillColorCode(int _fillColorCode) {
    fillColorCode=_fillColorCode;
}
void MyPDF::SetPDFBandType(string _PDFBandType) {
    PDFBandType=_PDFBandType;
}
void MyPDF::SetPDFErrorType(string _PDFErrorType) {
    PDFErrorType=_PDFErrorType;
}
void MyPDF::SetPDFErrorSize(string _PDFErrorSize) {
    PDFErrorSize=_PDFErrorSize;
}
void MyPDF::SetRenScaleValUp(double _renScaleVal) {
    renScaleValUp=_renScaleVal;
    renScaleVals[UP]=_renScaleVal;
    renScaleNameUp=to_string(_renScaleVal);
    renScaleNames[UP]="RenScale("+to_string(_renScaleVal)+")";
}
void MyPDF::SetRenScaleValDefault(double _renScaleVal) {
    renScaleValDefault=_renScaleVal;
    renScaleVals[DEF]=_renScaleVal;
    renScaleNameDefault=to_string(_renScaleVal);
    renScaleNames[DEF]="RenScale("+to_string(_renScaleVal)+")";
}
void MyPDF::SetRenScaleValDown(double _renScaleVal) {
    renScaleValDown=_renScaleVal;
    renScaleVals[DOWN]=_renScaleVal;
    renScaleNameDown=to_string(_renScaleVal);
    renScaleNames[DOWN]="RenScale("+to_string(_renScaleVal)+")";
}
void MyPDF::SetFacScaleValUp(double _facScaleVal) {
    facScaleValUp=_facScaleVal;
    facScaleVals[UP]=_facScaleVal;
    facScaleNameUp=to_string(_facScaleVal);
    facScaleNames[UP]="FacScale("+to_string(_facScaleVal)+")";
}
void MyPDF::SetFacScaleValDefault(double _facScaleVal) {
    facScaleValDefault=_facScaleVal;
    facScaleVals[DEF]=_facScaleVal;
    facScaleNameDefault=to_string(_facScaleVal);
    facScaleNames[DEF]="FacScale("+to_string(_facScaleVal)+")";
}
void MyPDF::SetFacScaleValDown(double _facScaleVal) {
    facScaleValDown=_facScaleVal;
    facScaleVals[DOWN]=_facScaleVal;
    facScaleNameDown=to_string(_facScaleVal);
    facScaleNames[DOWN]="FacScale("+to_string(_facScaleVal)+")";
}
void MyPDF::SetOptionsFileName(string _optionsFileName) {
    optionsFileName=_optionsFileName;
}
void MyPDF::SetDoPDFBand(bool _doit) {
    do_PDFBand = _doit;
}
void MyPDF::SetDoAplphaS(bool _doit) {
    do_AlphaS = _doit;
}
void MyPDF::SetDoRenormalizationScale(bool _doit) {
    do_RenormalizationScale = _doit;
}
void MyPDF::SetDoFactorizationScale(bool _doit) {
    do_FactorizationScale = _doit;
}
void MyPDF::SetDoTotError(bool _doit) {
    do_TotError = _doit;
}
void MyPDF::SetAlphaSmemberNumDown(int _memberNum) {
    AlphaSmemberNumDown=_memberNum;
}
void MyPDF::SetAlphaSmemberNumUp(int _memberNum) {
    AlphaSmemberNumUp=_memberNum;
}
void MyPDF::SetAlphaSPDFSetNameDown(string _name) {
    AlphaSPDFSetNameDown= _name;
}
void MyPDF::SetAlphaSPDFSetNameUp(string _name) {
    AlphaSPDFSetNameUp= _name;
}
void MyPDF::SetAlphaSPDFSetHistNameDown(string _name) {
    AlphaSPDFSetHistNameDown= _name;
}
void MyPDF::SetAlphaSPDFSetHistNameUp(string _name) {
    AlphaSPDFSetHistNameUp= _name;
}

void MyPDF::CleanUpMyPDF() {
    if(debug) std::cout<<"MyPDF::CleanUpMyPDF: Starting to clean up..."<<std::endl;

    if(h_errors_RenormalizationScale.size()>0) {
        for(int i=0; i<h_errors_RenormalizationScale.size(); ++i) {
            delete h_errors_RenormalizationScale.at(i);
        }
    }

    if(h_errors_FactorizationScale.size()>0) {
        for(int i=0; i<h_errors_FactorizationScale.size(); ++i) {
            delete h_errors_FactorizationScale.at(i);
        }
    }

    if(h_errors_PDFBand.size()>0) {
        for(int i=0; i<h_errors_PDFBand.size(); ++i) {
            delete h_errors_PDFBand.at(i);
        }
    }

    if(h_errors_prenorm.size()>0) {
        for(int i=0; i<h_errors_prenorm.size(); ++i) {
            delete h_errors_prenorm.at(i);
        }
    }

    if(h_errors_AlphaS.size()>0) {
        for(int i=0; i<h_errors_AlphaS.size(); ++i) {
            delete h_errors_AlphaS.at(i);
        }
    }

    if(h_errors_AlphaS_prenorm.size()>0) {
        for(int i=0; i<h_errors_AlphaS_prenorm.size(); ++i) {
            delete h_errors_AlphaS_prenorm.at(i);
        }
    }

    if(debug) std::cout<<"MyPDF::CleanUpMyPDF: Finished clean up!"<<std::endl;
}

