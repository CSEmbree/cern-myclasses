
#include "root.h"
#include <vector>
#include "utils.h"

#include "MyFrame.h"
#include "MyFrameData.h"
#include "MyData.h"
#include "MyCrossSection.h"
#include "MyPDF.h"

#include "AtlasStyle.h"
#include "TLegend.h"


void DumpTGraphAsymmErrors(TGraphAsymmErrors* my_graph ) ;
void DumpTH1D(TH1D* my_hist);
//void DrawRatioPlot(MyFrameData *myframe, TGraphAsymmErrors* reference_ratio, MyCrossSection *mycross, float x_min, float x_max, int error_code);
void DrawRatioPlot(MyFrameData *myframe, TGraphAsymmErrors* reference_ratio, MyCrossSection *mycross, float x_min, float x_max, int error_code, int igrid);
void DrawAndSaveResults(MyCrossSection *mycross, int igrid, int error_code);
void SaveSubprocessResults(MyCrossSection *mycross, int igrid, TString var_desc, TString x_title, bool logy, bool logx);
void SaveThisSubprocess(std::vector<TH1D*> hists_to_plot, MyCrossSection *mycross, TString plot_desc, TString x_title, bool logy, bool logx);
std::vector<bool> first_time_pdf; //vector indicating if a paritcular index has been drawn for the first time yet



int main(int argc, char** argv)
{
    SetAtlasStyle();
    first_time_pdf.clear();
    
    

    bool debug=false;
    string inputname="atlas2012_top.txt";
    if ( argc>1 ) inputname = TString(argv[1]);

    if (inputname.compare("atlas2012_top.txt")!=0) {

        TString psfile=inputname.c_str();
        psfile+=".eps";
        TString psfile2;

        MyCrossSection *mycross2= new MyCrossSection( (char*) inputname.c_str());
        mycross2->Print();
        int nframe=mycross2->GetFrameNumber();
        cout<<" Number of frames= "<<nframe<<endl;

        for (int iframe=0; iframe<nframe; iframe++) {
            std::cout << "DrawinFrame for " << iframe << std::endl;
            mycross2->DrawinFrame(iframe);
            std::cout << "After DrawinFrame for " << iframe << std::endl;
            if (iframe==0) psfile2=psfile+"(";
            else           psfile2=psfile;
            if (iframe==nframe-1) psfile2=psfile+")";
            std::cout << "Print mycross2" << std::endl;

            mycross2->GetMyFrame(iframe)->GetCanvas()->Print(psfile2);
            std::cout << "Printed" << std::endl;

        }
        std::cout << "Done with frame loop " << std::endl;
    }


    inputname="atlas2012_top.txt";

    MyCrossSection *mycross= new MyCrossSection( (char*) inputname.c_str());
    std::cout<<"mycross has this many pdfs: "<<mycross->getNPDF()<<std::endl;
    std::cout<<"Printing Crosssection..."<<std::endl;
    mycross->Print();
    std::cout << "Printed Cross section" << std::endl;
    
            int numPDFtypes = mycross->GetNumPDF();

    for (int igrid=0; igrid<numPDFtypes; igrid++) {
        std::cout << "MyFrameData for igrid " << igrid <<", of NGrids: "<<mycross->GetNGrid()<<", and numPDFtypes: "<<numPDFtypes<<std::endl;

        first_time_pdf.push_back(true);
        
        TH1D* h_reference = mycross->GetReference(igrid);        

        double xscale_factor = mycross->GetMyData(igrid)->GetUnitGeVFactor();
        std::cout << "for grid: " << mycross->GetGridName(igrid) << ", xscale_factor is " << xscale_factor << std::endl;
        mycross->GetMyPDF(igrid)->Print();

        mycross->mypdfInitializeErrorGraphs(igrid);
        mycross->mypdfCalcSystErrors(igrid);
        mycross->mypdfGetRatioToTH1(igrid, h_reference);

       
        DrawAndSaveResults(mycross, igrid, 0);

        SaveSubprocessResults(mycross, igrid, mycross->GetVarDesc(igrid), mycross->GetMyData(igrid)->GetTitleX(), mycross->GetMyData(igrid)->GetLogY1(), mycross->GetMyData(igrid)->GetLogX());
    }  /// igrid

    return 0;
}




/////////////////////////////////////////////////////////////////////////////////////////////
void SaveSubprocessResults(MyCrossSection *mycross, int igrid, TString var_desc, TString x_title, bool logy, bool logx)
{
    std::cout << "Save sub1" << std::endl;
    std::vector<TH1D*> hists_to_plot_gg_prenorm;
    hists_to_plot_gg_prenorm.clear();
    std::vector<TH1D*> hists_to_plot_qqbar_prenorm;
    hists_to_plot_qqbar_prenorm.clear();
    std::vector<TH1D*> hists_to_plot_tot_prenorm;
    hists_to_plot_tot_prenorm.clear();
    std::vector<TH1D*> hists_to_plot_gg;
    hists_to_plot_gg.clear();
    std::vector<TH1D*> hists_to_plot_qqbar;
    hists_to_plot_qqbar.clear();
    std::vector<TH1D*> hists_to_plot_tot;
    hists_to_plot_tot.clear();
    std::vector<TH1D*> hists_to_plot_gg_frac;
    hists_to_plot_gg_frac.clear();
    std::vector<TH1D*> hists_to_plot_qqbar_frac;
    hists_to_plot_qqbar_frac.clear();
    std::vector<int> PDF_codes;
    std::cout << "Save sub2" << std::endl;

    //int numPDFtypes = mypdf[0]->getNumPDFtypes();
    int numPDFtypes = mycross->GetNumPDF();
    
    for(int pdfi = 0; pdfi < numPDFtypes; pdfi++)
    {
        MyPDF *current_pdf = mycross->GetMyPDF(pdfi);
        std::cout << "Save sub2, pdfi: " << pdfi << std::endl;
        hists_to_plot_gg_prenorm.push_back      (current_pdf->h_gg_prenorm);
        hists_to_plot_qqbar_prenorm.push_back   (current_pdf->h_qqbar_prenorm);
        hists_to_plot_tot_prenorm.push_back     (current_pdf->h_tot_prenorm);
        hists_to_plot_gg.push_back              (current_pdf->h_gg);
        hists_to_plot_gg_frac.push_back         (current_pdf->h_gg_frac);
        hists_to_plot_qqbar_frac.push_back      (current_pdf->h_qqbar_frac);
        hists_to_plot_qqbar.push_back           (current_pdf->h_qqbar);
        hists_to_plot_tot.push_back             (current_pdf->h_tot);
        PDF_codes.push_back(pdfi);
    }

    SaveThisSubprocess(hists_to_plot_gg_prenorm,    mycross, TString(var_desc + "_gg_prenorm"),    x_title, logy, logx);
    SaveThisSubprocess(hists_to_plot_qqbar_prenorm, mycross, TString(var_desc + "_qqbar_prenorm"), x_title, logy, logx);
    SaveThisSubprocess(hists_to_plot_tot_prenorm,   mycross, TString(var_desc + "_tot_prenorm"),   x_title, logy, logx);
    SaveThisSubprocess(hists_to_plot_gg,            mycross, TString(var_desc + "_gg"),            x_title, logy, logx);
    SaveThisSubprocess(hists_to_plot_gg_frac,       mycross, TString(var_desc + "_gg_frac"),       x_title, false, logx);
    SaveThisSubprocess(hists_to_plot_qqbar_frac,    mycross, TString(var_desc + "_qqbar_frac"),    x_title, false, logx);
    SaveThisSubprocess(hists_to_plot_qqbar,         mycross, TString(var_desc + "_qqbar"),         x_title, logy, logx);
    SaveThisSubprocess(hists_to_plot_tot,           mycross, TString(var_desc + "_tot"),           x_title, logy, logx);


    for(int pdfi = 0; pdfi < numPDFtypes; pdfi++) 
    {
        std::cout << "Get inclusive cross-section for " << var_desc << ", PDF: " << mycross->GetMyPDF(pdfi)->getPDFBandType() << ": " << std::endl;
        float inclusive_cross_section = 0.;
        for(int bi = 1; bi <= hists_to_plot_tot_prenorm.at(pdfi)->GetNbinsX(); bi++) {
            float bin_width = hists_to_plot_tot_prenorm.at(pdfi)->GetBinWidth(bi) / 1000.;
            float bin_center = hists_to_plot_tot_prenorm.at(pdfi)->GetBinCenter(bi);
            float bin_content = hists_to_plot_tot_prenorm.at(pdfi)->GetBinContent(bi);
            inclusive_cross_section += (bin_width*bin_content);
            std::cout << "center: " << bin_center << ", width: " << bin_width << ", content: " << bin_content << ", Add " << (bin_width*bin_content) << ", cross-section now is " << inclusive_cross_section << "\n";

        }
        std::cout << "For " << mycross->GetMyPDF(pdfi)->getPDFBandType() << ", inclusive cross-section  is " << inclusive_cross_section << "\n";
    }
}

//void SaveThisSubprocess(std::vector<TH1D*> hists_to_plot, std::vector<int> PDF_codes, TString plot_desc, TString x_title, bool logy, bool logx)
void SaveThisSubprocess(std::vector<TH1D*> hists_to_plot, MyCrossSection *mycross, TString plot_desc, TString x_title, bool logy, bool logx)
{
    //int numPDFtypes = mypdf[0]->getNumPDFtypes();
    int numPDFtypes = mycross->GetNumPDF();

    TString plot_dir = "images/subprocesses/";
    TCanvas *print_canv = new TCanvas("print_canv", "", 600, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    TLegend *my_leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    my_leg->SetFillColor(0);
    my_leg->SetBorderSize(0);

    float max_val = 0.;
    float min_val = 9999999.;
    std::cout << "Spew hist contents for " << plot_desc << "\n";

    //for(int pdfi = 0; pdfi < PDF_codes.size(); pdfi++) {
    for(int pdfi = 0; pdfi < numPDFtypes; pdfi++) {
        //std::cout << "PDF: " << mypdf[pdfi]->getPDFname() << ": "<<endl;
        std::cout << "PDF: " << mycross->GetMyPDF(pdfi)->getPDFname() << ": "<<endl;
        
        
        std::cout << "TEST: hists_to_plot.size(): "<<hists_to_plot.size()<<std::endl;
        //exit(0); //TEST

        for(int bi = 1; bi <= hists_to_plot.at(pdfi)->GetNbinsX(); bi++) {
            
            float center = hists_to_plot.at(pdfi)->GetBinCenter(bi);
            float content = hists_to_plot.at(pdfi)->GetBinContent(bi);

            if( content > max_val ) max_val = content;
            if( content < min_val && content > 0.0000000001 ) min_val = content;
            std::cout << "center: " << center << ", content: " << content << ", max so far: " << max_val << ", min: " << min_val << "\n";
        }
    }

    max_val *= 1.4;
    if( logy ) {
        max_val *= 10.;
        min_val *= 0.1;
    }

    if( !logy ) min_val = 0.;
    std::cout << "logx: " << logx << ", logy: " << logy << ", max val: " << max_val << ", min_val: " << min_val << "\n";
    hists_to_plot.at(0)->GetYaxis()->SetRangeUser(min_val, max_val);

    if( logy ) gPad->SetLogy();
    else  gPad->SetLogy(0);
    if( logx ) gPad->SetLogx();
    else gPad->SetLogx(0);
    hists_to_plot.at(0)->SetXTitle(x_title);
    hists_to_plot.at(0)->Draw();

    for(int pdfi = 1; pdfi < numPDFtypes; pdfi++) {
        hists_to_plot.at(pdfi)->Draw("same");
    }
    for(int pdfi = numPDFtypes-1; pdfi >= 0; pdfi--) {
        my_leg->AddEntry(hists_to_plot.at(pdfi), mycross->GetMyPDF(pdfi)->getPDFname().c_str());
    }
    my_leg->Draw();

    std::cout << "cd to pad2" << std::endl;
    gPad->RedrawAxis();
    print_canv->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
    pad2->SetTopMargin(0);
    pad2->Draw();
    pad2->cd();
    std::vector<TH1D*> ratio_hists;

    float ratio_miny = 999.;
    float ratio_maxy = -999.;
    //for(int pdfi = 1; pdfi < PDF_codes.size(); pdfi++) {
    for(int pdfi = 1; pdfi < numPDFtypes; pdfi++) {
        TString this_ratio_name = hists_to_plot.at(pdfi)->GetName();
        this_ratio_name += "_ratio";
        TH1D* this_ratio_hist = (TH1D*) hists_to_plot.at(pdfi)->Clone(this_ratio_name);
        this_ratio_hist->Divide(hists_to_plot.at(0));

        for(int bi = 1; bi <= this_ratio_hist->GetNbinsX(); bi++) {
            float content = this_ratio_hist->GetBinContent(bi);
            if( content < ratio_miny ) ratio_miny = content;
            if( content > ratio_maxy ) ratio_maxy = content;
        }
        ratio_hists.push_back(this_ratio_hist);
    }  /// pdfi

    std::cout << "collected ratio hists" << std::endl;
    float range = ratio_maxy - ratio_miny;
    if( range > 0. ) {
        ratio_miny -= 0.2 * range;
        ratio_maxy += 0.2 * range;
    }

    std::cout << "Draw ratio hists" << std::endl;
    if( ratio_miny < 0. || ratio_miny > 2. ) ratio_miny = 0;
    if( ratio_maxy < 0.5 || ratio_maxy > 3. ) ratio_maxy = 2.;
    if( logx ) gPad->SetLogx();
    else gPad->SetLogx(0);
    ratio_hists.at(0)->GetYaxis()->SetRangeUser(ratio_miny, ratio_maxy);
    ratio_hists.at(0)->Draw();
    ratio_hists.at(0)->SetXTitle(x_title);
    ratio_hists.at(0)->SetYTitle("Ratio to CT10");

    for(int pdfi = 1; pdfi < ratio_hists.size(); pdfi++) {
        ratio_hists.at(pdfi)->Draw("same");
    }
    print_canv->cd();

    print_canv->Print((TString) (plot_dir + plot_desc + ".eps"));
}


void DrawAndSaveResults(MyCrossSection *mycross, int igrid, int error_code)
{
    std::cout << "Draw and save result" << std::endl;
    TString experiment_and_year = (TString) (mycross->GetMyData(igrid)->GetExperiment() + "_");
    experiment_and_year += mycross->GetMyData(igrid)->GetYear();
    experiment_and_year += " Data";

    TLegend *this_leg = new TLegend(0.4, 0.55, 0.9, 0.9);
    this_leg->SetBorderSize(0);
    this_leg->SetFillColor(0);
    this_leg->AddEntry(mycross->GetMyData(igrid)->GetTGraphTotErr(), experiment_and_year, "pl");

    //int numPDFtypes = mypdf[0]->getNumPDFtypes();
    int numPDFtypes = mycross->GetNumPDF();

    TString the_desc = "";
    for(int pdfi = 0; pdfi < numPDFtypes; pdfi++) {
        the_desc = mycross->GetMyPDF(pdfi)->calc_desc;
        the_desc += (TString) ("_" + mycross->GetMyPDF(pdfi)->getPDFErrorType() + "_");
        the_desc += mycross->GetVarDesc(igrid);
    }
    std::cout << "Draw and save results for " << the_desc << std::endl;

    double y=0.7;
    
    //potential issue - createing a new TLegend each time when we only want to update the "igrid" 
    MyFrameData *myframe= new MyFrameData(600,600,mycross->GetMyData(igrid));
    if (!myframe) cout<<" myframe not found "<<endl;
    else          cout<<" myframe created "<<endl;
    std::cout << "Now will draw igrid " << igrid << std::endl;

    mycross->SetLabelY(y);
    bool first_of_canv = true;
    float x_min = mycross->GetMyData(igrid)->GetTGraphTotErr()->GetXaxis()->GetBinLowEdge(1);
    int n_bins = mycross->GetMyData(igrid)->GetTGraphTotErr()->GetXaxis()->GetNbins();
    float x_max = mycross->GetMyData(igrid)->GetTGraphTotErr()->GetXaxis()->GetBinUpEdge(n_bins);
    std::cout << "Draw errors for desc: " << the_desc << std::endl;

    //mycross->DrawErrors(mycross->GetMyData(igrid)->GetTitleX(), x_min, x_max, first_of_canv, error_code);
    mycross->DrawError(igrid, mycross->GetMyData(igrid)->GetTitleX(), x_min, x_max, first_of_canv, error_code);
    std::cout << "Draw individual" << std::endl;

    mycross->Draw(igrid);
    std::cout << "Draw legend" << std::endl;
    std::cout << "Now will get frame " << std::endl;

    TGraphAsymmErrors* reference_ratio=mycross->GetReferenceRatio(igrid);
    reference_ratio->SetName((TString) (the_desc + "_reference_ratio"));
    TGraphAsymmErrors* reference_graph = TH1TOTGraphAsymm(mycross->GetReference(igrid));
    reference_graph->SetName((TString) (the_desc + "_reference"));

    double chi2[numPDFtypes];
    double chi2perdof[numPDFtypes];
    double chi2prob[numPDFtypes];
    TString chi2prob_str[numPDFtypes];
    /*
    for(int pdfi = 0; pdfi < numPDFtypes; pdfi++) {
        MyPDF *current_pdf = mycross->GetMyPDF(pdfi);
        chi2[pdfi] = -999;
        //std::cout << "For " << the_desc << ", " << mypdf[pdfi]->getPDFBandType() << ", determine chi2 to the data:" << std::endl;
        std::cout << "For " << the_desc << ", " << current_pdf->getPDFBandType() << ", at: "<<pdfi<<", determine chi2 to the data:" << std::endl;
        std::cout << "PDFBand results has name: "<< current_pdf->h_PDFBand_results << std::endl;
        std::cout << "data has name: "<< (mycross->GetMyData(igrid)->GetTGraphTotErr())->GetName() << std::endl;
        std::cout << "Covariance matrix ncols: " << ( mycross->GetMyData(igrid)->GetCovarianceMatrix()).GetNcols() << std::endl;
        std::cout << "chi2: " << chi2[pdfi] << std::endl;
        
        if( current_pdf->getDoPDFBand() ) 
            current_pdf->CalcChi2(current_pdf->h_PDFBand_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[pdfi]);
        if( current_pdf->getDoAlphaS() ) 
            current_pdf->CalcChi2(current_pdf->h_AlphaS_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[pdfi]);
        if( current_pdf->getDoRenormalizationScale() ) 
            current_pdf->CalcChi2(current_pdf->h_RenormalizationScale_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[pdfi]);
        if( current_pdf->getDoFactorizationScale() ) 
            current_pdf->CalcChi2(current_pdf->h_FactorizationScale_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[pdfi]);
        if( current_pdf->getDoTotError() ) 
            current_pdf->CalcChi2(current_pdf->h_TotError_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[pdfi]);

        chi2perdof[pdfi] = chi2[pdfi] / current_pdf->h_TotError_results->GetN();
        chi2prob[pdfi] = TMath::Prob(chi2[pdfi], current_pdf->h_TotError_results->GetN());
        std::stringstream strstream;
        strstream << setprecision(2) << chi2prob[pdfi];
        std::string cpp_format = strstream.str();
        chi2prob_str[pdfi] = "";
        chi2prob_str[pdfi] += cpp_format;
        this_leg->AddEntry(current_pdf->h_TotError_results, (TString) (current_pdf->getPDFBandType()+", #chi^{2} prob = "+chi2prob_str[pdfi]), "f");
        std::cout << "For " << the_desc << ", " << current_pdf->getPDFBandType() << ", chi2 = " << chi2[pdfi] << ", / dof = " << chi2perdof[pdfi] << ", prob = " << chi2prob[pdfi] << ", string format: " << chi2prob_str[pdfi] << std::endl;
    }
    */
    

        MyPDF *current_pdf = mycross->GetMyPDF(igrid);
        chi2[igrid] = -999;
        //std::cout << "For " << the_desc << ", " << mypdf[igrid]->getPDFBandType() << ", determine chi2 to the data:" << std::endl;
        std::cout << "For " << the_desc << ", " << current_pdf->getPDFBandType() << ", at: "<<igrid<<", determine chi2 to the data:" << std::endl;
        std::cout << "PDFBand results has name: "<< current_pdf->h_PDFBand_results << std::endl;
        std::cout << "data has name: "<< (mycross->GetMyData(igrid)->GetTGraphTotErr())->GetName() << std::endl;
        std::cout << "Covariance matrix ncols: " << ( mycross->GetMyData(igrid)->GetCovarianceMatrix()).GetNcols() << std::endl;
        std::cout << "chi2: " << chi2[igrid] << std::endl;
        
        if( current_pdf->getDoPDFBand() ) 
            current_pdf->CalcChi2(current_pdf->h_PDFBand_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[igrid]);
        if( current_pdf->getDoAlphaS() ) 
            current_pdf->CalcChi2(current_pdf->h_AlphaS_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[igrid]);
        if( current_pdf->getDoRenormalizationScale() ) 
            current_pdf->CalcChi2(current_pdf->h_RenormalizationScale_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[igrid]);
        if( current_pdf->getDoFactorizationScale() ) 
            current_pdf->CalcChi2(current_pdf->h_FactorizationScale_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[igrid]);
        if( current_pdf->getDoTotError() ) 
            current_pdf->CalcChi2(current_pdf->h_TotError_results, mycross->GetMyData(igrid)->GetTGraphTotErr(), mycross->GetMyData(igrid)->GetCovarianceMatrix(), chi2[igrid]);

        chi2perdof[igrid] = chi2[igrid] / current_pdf->h_TotError_results->GetN();
        chi2prob[igrid] = TMath::Prob(chi2[igrid], current_pdf->h_TotError_results->GetN());
        std::stringstream strstream;
        strstream << setprecision(2) << chi2prob[igrid];
        std::string cpp_format = strstream.str();
        chi2prob_str[igrid] = "";
        chi2prob_str[igrid] += cpp_format;
        this_leg->AddEntry(current_pdf->h_TotError_results, (TString) (current_pdf->getPDFBandType()+", #chi^{2} prob = "+chi2prob_str[igrid]), "f");
        std::cout << "For " << the_desc << ", " << current_pdf->getPDFBandType() << ", chi2 = " << chi2[igrid] << ", / dof = " << chi2perdof[igrid] << ", prob = " << chi2prob[igrid] << ", string format: " << chi2prob_str[igrid] << std::endl;
   
    
    
    this_leg->Draw();

    myframe->GetMyFrame()->GetSubPad2()->cd();
    std::cout << "Now get ref ratio" << std::endl;
    std::cout << "Got ref ratio" << std::endl;
    //DrawRatioPlot(myframe, reference_ratio, mycross, x_min, x_max, error_code);
        DrawRatioPlot(myframe, reference_ratio, mycross, x_min, x_max, error_code, igrid);
    gPad->Update();
    myframe->SaveFile((TString) ("images/ErrorResults/Results_" + the_desc  + "_overlay.eps"));
}

void DumpTH1D(TH1D* my_hist)
{
    std::cout << "hist has " << my_hist->GetNbinsX() << " bins" << std::endl;
    for(int bi = 0; bi <= my_hist->GetNbinsX()+1; bi++) {
        float center = my_hist->GetBinCenter(bi);
        std::cout << "bin " << bi << ", center: " << center << ", content = " << my_hist->GetBinContent(bi) << " +/- " << my_hist->GetBinError(bi) << "\n";
    }
}

void DumpTGraphAsymmErrors(TGraphAsymmErrors* my_graph )
{
    float sum = 0.;
    std::cout << "Dump TGraphAsymm\n";
    for(int pi = 0; pi < my_graph->GetN(); pi++) {

        Double_t x_val;
        Double_t y_val;
        my_graph->GetPoint(pi, x_val, y_val);
        Double_t exh = my_graph->GetErrorXhigh(pi);
        Double_t exl = my_graph->GetErrorXlow(pi);
        Double_t eyh = my_graph->GetErrorYhigh(pi);
        Double_t eyl = my_graph->GetErrorYlow(pi);

        double width = 2.*exh;
        double product = width*y_val;
        sum += product;
        std::cout << "bin center: " << x_val << " +" << exh << " -" << exl << ", y: " << y_val << " +" << eyh << " -" << eyl << "\n";
        std::cout << "product: " << product << ", sum so far: " << sum << "\n";

    }
    std::cout << "Total sum for this TGraph: " << sum << "\n";
}

void DrawRatioPlot(MyFrameData *myframe, TGraphAsymmErrors* reference_ratio, MyCrossSection *mycross, float x_min, float x_max, int error_code, int igrid)
{
    //int numPDFtypes = mypdf[0]->getNumPDFtypes();
    int numPDFtypes = mycross->GetNumPDF();

    std::cout << "Now compute range" << std::endl;
    double xmin=0., xmax=0., ymin=0., ymax=0.;

    ymin = 999.;
    ymax = -999.;
    for(int pi = 0; pi < reference_ratio->GetN(); pi++) {
        Double_t y_val;
        Double_t x_val;
        reference_ratio->GetPoint(pi, x_val, y_val);
        Double_t eyl = reference_ratio->GetErrorYlow(pi);
        Double_t eyh = reference_ratio->GetErrorYhigh(pi);
        if( ymin > y_val-eyl ) ymin = y_val-eyl;
        if( ymax < y_val+eyh ) ymax = y_val+eyh;
    }

    double ydiff = ymax - ymin;
    ymax += 0.2*ydiff;
    ymin -= 0.2*ydiff;

    if( x_min < x_max ) {
        myframe->GetMyFrame()->GetXAxis2()->SetRangeUser(x_min,x_max);
    }
    myframe->GetMyFrame()->GetYAxis2()->SetRangeUser(ymin,ymax);
    myframe->GetMyFrame()->GetYAxis2()->SetTitle("Hist / NLO Pred.");

/*
    bool first_time_pdf = true;
    
    for(int pdfi = 0; pdfi < numPDFtypes; pdfi++) 
    {
        MyPDF *current_pdf = mycross->GetMyPDF(pdfi);
        
        if( current_pdf->getDoPDFBand() ) {
            if( first_time_pdf )    current_pdf->h_PDFBand_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_PDFBand_results_ratio_to_ref->Draw("e2 same");
         }
         else if( current_pdf->getDoAlphaS() ) {
            if( first_time_pdf )    current_pdf->h_AlphaS_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_AlphaS_results_ratio_to_ref->Draw("e2 same");
        }
        else if( current_pdf->getDoRenormalizationScale() ) {
            if( first_time_pdf )    current_pdf->h_RenormalizationScale_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_RenormalizationScale_results_ratio_to_ref->Draw("e2 same");
        }
        else if( current_pdf->getDoFactorizationScale() ) {
            if( first_time_pdf )    current_pdf->h_FactorizationScale_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_FactorizationScale_results_ratio_to_ref->Draw("e2 same");
        }
        else if( current_pdf->getDoTotError() ) {
            if( first_time_pdf )    current_pdf->h_TotError_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_TotError_results_ratio_to_ref->Draw("e2 same");
         }
         
       if(first_time_pdf) first_time_pdf = false;
    }
    */
    
    
        MyPDF *current_pdf = mycross->GetMyPDF(igrid);
        
        if( current_pdf->getDoPDFBand() ) {
            if( first_time_pdf.at(igrid) )      current_pdf->h_PDFBand_results_ratio_to_ref->Draw("e2");
            else                                current_pdf->h_PDFBand_results_ratio_to_ref->Draw("e2 same");
         }
         else if( current_pdf->getDoAlphaS() ) {
            if( first_time_pdf.at(igrid) )    current_pdf->h_AlphaS_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_AlphaS_results_ratio_to_ref->Draw("e2 same");
        }
        else if( current_pdf->getDoRenormalizationScale() ) {
            if( first_time_pdf.at(igrid) )    current_pdf->h_RenormalizationScale_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_RenormalizationScale_results_ratio_to_ref->Draw("e2 same");
        }
        else if( current_pdf->getDoFactorizationScale() ) {
            if( first_time_pdf.at(igrid) )    current_pdf->h_FactorizationScale_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_FactorizationScale_results_ratio_to_ref->Draw("e2 same");
        }
        else if( current_pdf->getDoTotError() ) {
            if( first_time_pdf.at(igrid) )    current_pdf->h_TotError_results_ratio_to_ref->Draw("e2");
            else                    current_pdf->h_TotError_results_ratio_to_ref->Draw("e2 same");
         }
         
         //v1.insert(v1.begin()+i, v2[i])
         first_time_pdf.insert(first_time_pdf.begin()+igrid,false);
       //if(first_time_pdf) first_time_pdf = false;
    
    

    reference_ratio->Draw("ep same");
}
