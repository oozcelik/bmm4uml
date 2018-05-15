#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TFileCollection.h>
#include <TFileInfo.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TNtupleD.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TTree.h>

#include <Math/Util.h>

#include <RooAcceptReject.h>
#include <RooAddition.h>
#include <RooAddModel.h>
#include <RooAddPdf.h>
#include <RooBernstein.h>
#include <RooBifurGauss.h>
#include <RooCategory.h>
#include <RooCBShape.h>
#include <RooCFunction1Binding.h>
#include <RooChebychev.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooDecay.h>
#include <RooEffProd.h>
#include <RooExponential.h>
#include <RooFormulaVar.h>
#include <RooFitResult.h>
#include <RooGamma.h>
#include <RooGaussian.h>
#include <RooGaussModel.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooHistPdf.h>
#include <RooKeysPdf.h>
#include <RooLinkedList.h>
#include <RooLognormal.h>
#include <RooMinuit.h>
#include <RooPlot.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooRandom.h>
#include <RooRealConstant.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooTruthModel.h>
#include <RooWorkspace.h>

using namespace RooFit;
using namespace std;

// global configuration (enable one of these below)
#define CONFIG_BMM3     0
#define CONFIG_BMM4     1

#define NCPU            1
#define MC_EVENT_LIMIT  200000
#define USING_MINUIT2   0

// variables that are essential, included in the reduced minos fit & ploting
const vector<TString> POI_list = {"BF_bs", "BF_bd", "dblmu_corr_scale", "EffTau_bs"};

// boundaries for observables
const double Mass_bound[2] = {4.9, 5.9};
const double ReducedMassRes_bound[2] = {0.0009, 0.045};
const double BDT_bound[2] = {0.1, 1.0};
const double Tau_bound[2] = {1., 12.}; // in unit of ps
const double TauRes_bound[2] = {0.01, 0.25}; // in unit of ps

// bin boundries for decay time / decay time resolution sPlot
const vector<double> Tau_bins = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
const vector<double> TauRes_bins = {
     0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
     0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20,
     0.21, 0.22, 0.23, 0.24, 0.25};

// available BMM3 samples
namespace bmm3 {
    const int ndecays = 20;
    const TString eras[] = {"2011", "2012"};
    const TString decays[ndecays] = {"SgMc", "BdMc",
        "bgBs2KK", "bgBs2KPi", "bgBs2PiPi", "bgBd2KK", "bgBd2KPi", "bgBd2PiPi", "bgLb2PiP", "bgLb2KP",
        "bgBs2KMuNu", "bgBd2PiMuNu", "bgLb2PMuNu",
        "bgBu2PiMuMu", "bgBu2KMuMu", "bgBd2Pi0MuMu", "bgBd2K0MuMu", "bgBd2MuMuGamma", "bgBs2MuMuGamma",
        "SgData"};
    const TString fullprocess[ndecays] = {"B_{s}#rightarrow #mu^{+}#mu^{-}", "B^{0}#rightarrow #mu^{+}#mu^{-}",
        "B_{s}#rightarrow K^{+}K^{-}", "B_{s}#rightarrow K^{+}#pi^{-}", "B_{s}#rightarrow #pi^{+}#pi^{-}",
        "B^{0}#rightarrow K^{+}K^{-}", "B^{0}#rightarrow K^{+}#pi^{-}", "B^{0}#rightarrow #pi^{+}#pi^{-}",
        "#Lambda_{b}#rightarrow p^{+}#pi^{-}", "#Lambda_{b}#rightarrow p^{+}K^{-}",
        "B_{s}#rightarrow K^{+}#mu^{-}#nu", "B^{0}#rightarrow #pi^{+}#mu^{-}#nu", "#Lambda_{b}#rightarrow p^{+}#mu^{-}#nu",
        "B^{+}#rightarrow #pi^{+}#mu^{+}#mu^{-}", "B^{+}#rightarrow K^{+}#mu^{+}#mu^{-}", "B^{0}#rightarrow #pi^{0}#mu^{+}#mu^{-}",
        "B^{0}#rightarrow #pi^{0}#mu^{+}#mu^{-}", "B^{0}#rightarrow #mu^{+}#mu^{-}#gamma", "B_{s}#rightarrow #mu^{+}#mu^{-}#gamma",
        "Data"};
    enum { _SgMc, _BdMc,
        _bgBs2KK, _bgBs2KPi, _bgBs2PiPi, _bgBd2KK, _bgBd2KPi, _bgBd2PiPi, _bgLb2PiP, _bgLb2KP,
        _bgBs2KMuNu, _bgBd2PiMuNu, _bgLb2PMuNu,
        _bgBu2PiMuMu, _bgBu2KMuMu, _bgBd2Pi0MuMu, _bgBd2K0MuMu, _bgBd2MuMuGamma, _bgBs2MuMuGamma,
        _SgData};
}

// available BMM4 samples
namespace bmm4 {
    const int ndecays = 20;
    const TString eras[] = {"2011s01","2012s01","2016BFs01","2016GHs01"};
    const TString decays[ndecays] = {"bsmmMc", "bdmmMc",
        "bskkMcBg", "bskpiMcBg", "bspipiMcBg", "bdkkMcBg", "bdkpiMcBg", "bdpipiMcBg", "lbppiMcBg", "lbpkMcBg",
        "bskmunuMcBg", "bdpimunuMcBg", "lbpmunuMcBg",
        "bdpimumuMcBg", "bupimumuMcBg",
        "bmmData",
        "bspsiphiMc", "bspsiphiData",
        "bupsikMc", "bupsikData"};
    const TString texlabels[ndecays] = {"bsmm", "bdmm",
        "bskkBg", "bskpiBg", "bspipiBg", "bdkkBg", "bdkpiBg", "bdpipiBg", "lbppiBg", "lbpkBg",
        "bskmunuBg", "bdpimunuBg", "lbpmunuBg",
        "bdpimumuBg", "bupimumuBg",
        "bmm",
        "bspsiphi", "bspsiphi",
        "bupsik", "bupsik"};
    const TString fullprocess[ndecays] = {"B_{s}#rightarrow #mu^{+}#mu^{-}", "B^{0}#rightarrow #mu^{+}#mu^{-}",
        "B_{s}#rightarrow K^{+}K^{-}", "B_{s}#rightarrow K^{+}#pi^{-}", "B_{s}#rightarrow #pi^{+}#pi^{-}",
        "B^{0}#rightarrow K^{+}K^{-}", "B^{0}#rightarrow K^{+}#pi^{-}", "B^{0}#rightarrow #pi^{+}#pi^{-}",
        "#Lambda_{b}#rightarrow p^{+}#pi^{-}", "#Lambda_{b}#rightarrow p^{+}K^{-}",
        "B_{s}#rightarrow K^{+}#mu^{-}#nu", "B^{0}#rightarrow #pi^{+}#mu^{-}#nu", "#Lambda_{b}#rightarrow p^{+}#mu^{-}#nu",
        "B^{0}#rightarrow #pi^{0}#mu^{+}#mu^{-}","B^{+}#rightarrow #pi^{+}#mu^{+}#mu^{-}",
        "Data",
        "B_{s}#rightarrow J/#psi#phi", "B_{s}#rightarrow J/#psi#phi",
        "B^{+}#rightarrow J/#psi K^{+}", "B^{+}#rightarrow J/#psi K^{+}"};
    enum { _bsmmMc, _bdmmMc,
        _bskkMcBg, _bskpiMcBg, _bspipiMcBg, _bdkkMcBg, _bdkpiMcBg, _bdpipiMcBg, _lbppiMcBg, _lbpkMcBg,
        _bskmunuMcBg, _bdpimunuMcBg, _lbpmunuMcBg,
        _bdpimumuMcBg, _bupimumuMcBg,
        _bmmData,
        _bspsiphiMc, _bspsiphiData,
        _bupsikMc, _bupsikData};
}

// CategoryMananger as the binning control
class CategoryMananger {
public:
    
    struct CateDef{
        TString id;    // unique ID string
        int index;     // global binning index
        TString era;   // era = 2011, 2012, 2016BF, 2016GH
        int region;    // channel separation based on detector region
        int bdt_bin;   // binning index in BDT
        double bdt_min, bdt_max;
    };
    
    vector<CateDef> cats;
    
    CategoryMananger() {}
    
    void RegisterCategory(TString era, int region, int bdt_bin, double bdt_min, double bdt_max) {
        CateDef cat;
        cat.id      = Form("%s_%d_%d",era.Data(), region, bdt_bin);
        cat.index   = cats.size();
        cat.era     = era;
        cat.region  = region;
        cat.bdt_bin = bdt_bin;
        cat.bdt_min = bdt_min;
        cat.bdt_max = bdt_max;
        cats.push_back(cat);
    }
    
    void RegisterBMM3Categories() {
        const int run1_channels = 4; // BDT-categories strategy
        const int run1_index_max[run1_channels] = {2,2,4,4};
        const double run1_bdt_bins[run1_channels][5] = {
            {0.10,0.31,1.00,0.,0.},      // 2011 - barrel
            {0.10,0.26,1.00,0.,0.},      // 2011 - endcap
            {0.10,0.23,0.33,0.44,1.00},  // 2012 - barrel
            {0.10,0.22,0.33,0.45,1.00}}; // 2012 - endcap
        
         /*const int run1_channels = 4; // BDT-cut strategy
         const int run1_index_max[run1_channels] = {1,1,1,1};
         const double run1_bdt_bins[run1_channels][2] = {
         {0.29,1.00},  // 2011 - barrel
         {0.29,1.00},  // 2011 - endcap
         {0.36,1.00},  // 2012 - barrel
         {0.38,1.00}}; // 2012 - endcap*/
        
        for (int i = 0; i < run1_channels; i++)
            for (int j = 0; j < run1_index_max[i]; j++)
                RegisterCategory(i<=1?"2011":"2012", i%2, j, run1_bdt_bins[i][j], run1_bdt_bins[i][j+1]);
    }
    void RegisterBMM4Categories() {
        
        /*
        // BDT-cut ref
        RegisterCategory("2011s01", 0, 0, 0.290, 1.00);
        RegisterCategory("2011s01", 1, 0, 0.290, 1.00);
        RegisterCategory("2012s01", 0, 0, 0.360, 1.00);
        RegisterCategory("2012s01", 1, 0, 0.380, 1.00);
        RegisterCategory("2016BFs01", 0, 0, 0.300, 1.00);
        RegisterCategory("2016BFs01", 1, 0, 0.440, 1.00);
        RegisterCategory("2016GHs01", 0, 0, 0.320, 1.00);
        RegisterCategory("2016GHs01", 1, 0, 0.380, 1.00);*/
        
        // BDT-cut
        RegisterCategory("2011s01", 0, 0, 0.34, 1.00);
        RegisterCategory("2011s01", 1, 0, 0.23, 1.00);
        RegisterCategory("2012s01", 0, 0, 0.34, 1.00);
        RegisterCategory("2012s01", 1, 0, 0.35, 1.00);
        RegisterCategory("2016BFs01", 0, 0, 0.28, 1.00);
        RegisterCategory("2016BFs01", 1, 0, 0.28, 1.00);
        RegisterCategory("2016GHs01", 0, 0, 0.26, 1.00);
        RegisterCategory("2016GHs01", 1, 0, 0.29, 1.00);
        
        /*
        // BDT-cat 1
        RegisterCategory("2011s01",   0, 0, 0.25, 0.38);
        RegisterCategory("2011s01",   1, 0, 0.21, 0.39);
        RegisterCategory("2012s01",   0, 0, 0.28, 0.37);
        RegisterCategory("2012s01",   1, 0, 0.32, 0.47);
        RegisterCategory("2016BFs01", 0, 0, 0.19, 0.33);
        RegisterCategory("2016BFs01", 1, 0, 0.19, 0.32);
        RegisterCategory("2016GHs01", 0, 0, 0.23, 0.34);
        RegisterCategory("2016GHs01", 1, 0, 0.21, 0.36);
        RegisterCategory("2011s01",   0, 1, 0.38, 1.00);
        RegisterCategory("2011s01",   1, 1, 0.39, 1.00);
        RegisterCategory("2012s01",   0, 1, 0.37, 1.00);
        RegisterCategory("2012s01",   1, 1, 0.47, 1.00);
        RegisterCategory("2016BFs01", 0, 1, 0.33, 1.00);
        RegisterCategory("2016BFs01", 1, 1, 0.32, 1.00);
        RegisterCategory("2016GHs01", 0, 1, 0.34, 1.00);
        RegisterCategory("2016GHs01", 1, 1, 0.36, 1.00);*/
        
        /*
        // BDT-cat 2
        RegisterCategory("2011s01",   0, 0, 0.25, 0.38);
        RegisterCategory("2011s01",   1, 0, 0.21, 0.39);
        RegisterCategory("2011s01",   0, 1, 0.38, 1.00);
        RegisterCategory("2011s01",   1, 1, 0.39, 1.00);
        RegisterCategory("2012s01",   0, 0, 0.25, 0.32);
        RegisterCategory("2012s01",   1, 0, 0.25, 0.35);
        RegisterCategory("2016BFs01", 0, 0, 0.14, 0.23);
        RegisterCategory("2016BFs01", 1, 0, 0.18, 0.31);
        RegisterCategory("2016GHs01", 0, 0, 0.17, 0.26);
        RegisterCategory("2016GHs01", 1, 0, 0.19, 0.30);
        RegisterCategory("2012s01",   0, 1, 0.32, 0.38);
        RegisterCategory("2012s01",   1, 1, 0.35, 0.48);
        RegisterCategory("2016BFs01", 0, 1, 0.23, 0.34);
        RegisterCategory("2016BFs01", 1, 1, 0.31, 0.44);
        RegisterCategory("2016GHs01", 0, 1, 0.26, 0.35);
        RegisterCategory("2016GHs01", 1, 1, 0.30, 0.40);
        RegisterCategory("2012s01",   0, 2, 0.38, 1.00);
        RegisterCategory("2012s01",   1, 2, 0.48, 1.00);
        RegisterCategory("2016BFs01", 0, 2, 0.34, 1.00);
        RegisterCategory("2016BFs01", 1, 2, 0.44, 1.00);
        RegisterCategory("2016GHs01", 0, 2, 0.35, 1.00);
        RegisterCategory("2016GHs01", 1, 2, 0.40, 1.00);*/

    }
    
    void Print() {
        cout << ">>> CategoryMananger: Summary" << endl;
        for (auto& cat : cats) {
            cout << ">>> CATEGORY " << cat.index << "(id: " << cat.id << ") - ";
            cout << cat.era << "/" << cat.region << "/" << cat.bdt_bin;
            cout << " - BDT range [" << cat.bdt_min << ", " << cat.bdt_max << "]" << endl;
        }
    }
    
    int index(TString era, int region, int bdt_bin) {
        for (int i=0; i<(int)cats.size(); i++) {
            if (cats[i].era == era && cats[i].region == region && cats[i].bdt_bin == bdt_bin)
                return cats[i].index;
        }
        return -1;
    }
    int index(TString era, int region, double bdt_val) {
        for (int i=0; i<(int)cats.size(); i++) {
            if (cats[i].era == era && cats[i].region == region &&
                bdt_val>=cats[i].bdt_min && bdt_val<=cats[i].bdt_max)
                return cats[i].index;
        }
        return -1;
    }
}CatMan;

// protecting non-converged fits
void converge_protection()
{
    cout << ">>> WARNING: fit does not fully converge!" << endl;
    cout << ">>> continue? (y/n)" << endl;
    string buffer;
    cin >> buffer;
    if (buffer[0]!='y' && buffer[0]!='Y') exit(1);
}

// protecting non-existing input file
void exist_protection(TString filename)
{
    if (FILE *file = fopen(filename.Data(), "r")) {
        fclose(file);
        return;
    }else {
        cout << ">>> ERROR: " << filename << " does not exist!" << endl;
        exit(1);
    }
}

// repeating the fit until converge or limited trial times
RooFitResult* fit_with_retry(RooAbsPdf& pdf, RooAbsData& data, RooLinkedList& list, int retry = 3)
{
    RooFitResult *res = pdf.fitTo(data,list);
    while (res->status()!=0 && retry>0) {
        delete res;
        res = pdf.fitTo(data,list);
        retry--;
    }
    return res;
}

// utilities to read BMM3-type fitter data cards (plain text)
void ReadValuesFromPlainText(TString filename, vector<TString> keys, vector<double>& values, vector<double>& errors)
{
    values.assign(keys.size(),0.);
    errors.assign(keys.size(),0.);
    vector<bool> flag(keys.size(), false);
    
    TString buffer;
    
    ifstream fin(filename);
    while(buffer.ReadLine(fin)) {
        if (buffer[0]=='#') continue;
        
        TString tok, keystring = "";
        Ssiz_t from = 0;
        double val = 0.,err = 0.;
        
        if (buffer.Tokenize(tok, from, "[ \t]+")) keystring = tok;
        if (buffer.Tokenize(tok, from, "[ \t]+")) val = tok.Atof();
        if (buffer.Tokenize(tok, from, "[ \t]+")) err = tok.Atof();
        
        for (unsigned int i=0; i<keys.size(); i++) {
            if (keystring.Contains(keys[i])) {
                values[i] = val;
                errors[i] = err;
                flag[i] = true;
            }
        }
    }
    for (unsigned int i=0; i<flag.size(); i++) { // force break if any key is missing
        if (!flag[i]) {
            cout << ">>> ERROR: key '" << keys[i] << "' not found in '" << filename << "'." << endl;
            exit(1);
        }
    }
}

// utilities to read BMM3/4-type fitter data cards (tex files)
void ReadValuesFromTex(TString filename, vector<TString> keys, vector<double>& values)
{
    values.assign(keys.size(), 0.);
    vector<bool> flag(keys.size(), false);
    
    TString buffer;
    ifstream fin(filename);
    while(buffer.ReadLine(fin)) {
        double val = 0.;
        char keystring[1024] = {0};
        sscanf(buffer.Data(), "%s   {\\ensuremath{{%lf } } }", keystring, &val);
        for (unsigned int i=0; i<keys.size(); i++) {
            if (strstr(keystring,keys[i].Data())) {
                values[i] = val;
                flag[i] = true;
            }
        }
    }
    for (unsigned int i=0; i<flag.size(); i++) { // force break if any key is missing
        if (!flag[i]) {
            cout << ">>> ERROR: key '" << keys[i] << "' not found in '" << filename << "'." << endl;
            exit(1);
        }
    }
}

class TexVar {
public:
    double val, estat, esyst, etot;
    
    TexVar() {val = estat = esyst = etot = 0.;}
    TexVar(TString filename, TString key) {
        ReadValues(filename,key);
    }
    
    TexVar& operator=(const TexVar& in)
    {
        val   = in.val;
        estat = in.estat;
        esyst = in.esyst;
        etot  = in.etot;
        return *this;
    }
    
    void ReadValues(TString filename, TString key) {
        vector<TString> keys;
        vector<double> values;

        keys.push_back(key+TString(":val"));
        keys.push_back(key+TString(":estat"));
        keys.push_back(key+TString(":esyst"));
        keys.push_back(key+TString(":etot"));
        
        ReadValuesFromTex(filename, keys, values);
        val   = values[0];
        estat = values[1];
        esyst = values[2];
        etot  = values[3];
    }
    
    void AddVar(TexVar in) {
        val += in.val;
        estat = sqrt(pow(estat,2) + pow(in.estat,2));
        esyst += in.esyst;
        etot = sqrt(pow(estat,2) + pow(esyst,2));
    }
    void SubVar(TexVar in, bool nonnegative = true) {
        val -= in.val;
        if (nonnegative && val<0.) val = 0.;
        estat = sqrt(max(0.,pow(estat,2) - pow(in.estat,2)));
        esyst -= in.esyst;
        if (nonnegative && esyst<0.) esyst = 0.;
        etot = sqrt(pow(estat,2) + pow(esyst,2));
    }
};

// sPlot fitting tool

TH1D *Fit_hist = NULL;
RooAbsPdf *Fit_pdf = NULL;
RooRealVar *Fit_obs = NULL;
bool Fit_w_weight = false;

double Fit_NLLVar(double)
{
    double f = 0.;
    
    Fit_obs->setRange("full", Fit_obs->getMin(), Fit_obs->getMax());
    RooAbsReal* area_full = Fit_pdf->createIntegral(*Fit_obs,Range("full"));
    double norm = area_full->getVal();
    delete area_full;
    
    double sum = 0.;
    for(int bin=1;bin<=Fit_hist->GetNbinsX();bin++)
        sum += Fit_hist->GetBinContent(bin);
    
    for(int bin=1;bin<=Fit_hist->GetNbinsX();bin++) {
        
        double x_min = Fit_hist->GetBinLowEdge(bin);
        double x_max = x_min + Fit_hist->GetBinWidth(bin);
        double y     = Fit_hist->GetBinContent(bin);
        double error = Fit_hist->GetBinError(bin);
        
        Fit_obs->setRange("bin", x_min, x_max);
        RooAbsReal* area = Fit_pdf->createIntegral(*Fit_obs,Range("bin"));
        
        double val   = max(area->getVal()/norm,0.);
        delete area;
        
        // match the implementation as in ROOT::Fit::FitUtil
        double nll = 0.;
        bool extended = true;
        if (Fit_w_weight) {
            if (y!=0.) {
                double weight = error*error/y;
                if (extended) nll = val * weight;
                nll -= weight * y * ROOT::Math::Util::EvalLog(val);
            }
        }else {
            if (extended) nll = val - y;
            if (y>0.) nll += y * (ROOT::Math::Util::EvalLog(y)-ROOT::Math::Util::EvalLog(val));
        }
        f+= nll;
    }
    return f;
}

// h_tau   - input histgram
// Tau_pdf - fitting PDF
// Tau     - decay time observable
// EffTau  - fitting effective lifetime
// res1    - results for step 1 fit
// res2    - results for step 2 fit
//
void Fit_sPlot(TH1D *h_tau, RooAbsPdf *Tau_pdf, RooRealVar *Tau, RooRealVar *EffTau, RooFitResult **res1, RooFitResult **res2)
{
    // Binned weighted likelihood fit w/ PDF bin integration
    Fit_hist = h_tau;
    Fit_pdf = Tau_pdf;
    Fit_obs = Tau;
    Fit_w_weight = false;
    
    RooAbsReal* NLLVar = bindFunction("NLLVar",Fit_NLLVar,*EffTau);
    RooMinuit *RooMin = new RooMinuit(*NLLVar);
    RooMin->setErrorLevel(0.5);
    RooMin->migrad();
    RooMin->hesse();
    RooMin->minos();
    (*res1) = RooMin->save();
    TMatrixDSym Cov1 = (*res1)->covarianceMatrix();
    
    Fit_w_weight = true;
    RooMin->hesse();
    (*res2) = RooMin->save();
    TMatrixDSym Cov2 = (*res2)->covarianceMatrix();
    TMatrixDSym Cov2inv(TMatrixDSym::kInverted,Cov2);
    TMatrixD Cov = Cov1*Cov2inv*Cov1;
    
    double scale = sqrt(Cov(0,0))/sqrt(Cov1(0,0)); // ~ original error / new weighted error
    RooRealVar* final_EffTau = (RooRealVar*)(*res1)->floatParsFinal().find(EffTau->GetName());
    final_EffTau->setError(final_EffTau->getError()*scale);
    final_EffTau->setAsymError(final_EffTau->getErrorLo()*scale, final_EffTau->getErrorHi()*scale);
    
    // Show fit results
    cout << ">>> Fit results:" << endl;
    cout << ">>> EffTau: " << final_EffTau->getVal() << " +- " <<
    final_EffTau->getError() << " (+" << final_EffTau->getErrorHi() << "/" << final_EffTau->getErrorLo() << ")" << endl;
    
    delete RooMin;
    delete NLLVar;
}
