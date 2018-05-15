#include "bmm4common.h"

// ----------------------------------------------------
// Prepare the resolution functions for lifetime fit
// current model: double-Gaussian, with <taue> as the main scaling factor
// available options:
// bsmm      - build B->mumu resolution model
// bupsik    - build B->J/psi K+ resolution model
// --
// triple    - using a triple Guassian model
// double    - reduced to a double Guassian model
// --
// even      - only use the even events
// odd       - only use the odd events
//
void PrepareLifetimeResolutionModel(RooWorkspace *wspace, TString opt = "bsmm:triple", TString target_cat = "2016BF", double alter_BDT = -2.)
{
    cout << ">>> PrepareLifetimeResolutionModel() start" << endl;
    
    TString tag, treename;
    
    if (opt.Contains("bsmm")) {
        tag = "bsmm";
        treename = "bsmmMc";
    }else if (opt.Contains("bupsik")) {
        tag = "bupsik";
        treename = "bupsikMc";
    }else {
        cout << ">>> Undefined option: " << opt << endl;
        return;
    }
    tag += Form("_%s", target_cat.Data());
    
    TString target_era = target_cat;
    for (auto& cat: CatMan.cats)
        if (cat.id == target_cat) target_era = cat.era;
    
    RooRealVar delTau("delTau","",-0.6,+0.6);
    TH1D *h_global_deltau = new TH1D("h_global_deltau","",240,-0.6,+0.6);
    TH1D *h_global_taue = new TH1D("h_global_taue","",192,TauRes_bound[0],TauRes_bound[1]);
    h_global_deltau->Sumw2();
    h_global_taue->Sumw2();
    
    // load the corresponding MC sample
    TString filename = Form("input/bmm4/small%s-%s.root",target_era.Data(),treename.Data());
    TChain *events = new TChain(treename);
    exist_protection(filename);
    events->Add(filename);
    cout << ">>> Loading from " << filename << ", with " << events->GetEntries() << " entries." << endl;
        
    double m_t, me_t, pt_t, tau_t, taue_t, gtau_t, bdt_t;
    int chan_t;
    bool muid_t, cnc_t;
        
    events->SetBranchAddress("m", &m_t);
    events->SetBranchAddress("me", &me_t);
    events->SetBranchAddress("pt", &pt_t);
    events->SetBranchAddress("tau", &tau_t);
    events->SetBranchAddress("taue", &taue_t);
    events->SetBranchAddress("gtau", &gtau_t);
    events->SetBranchAddress("bdt", &bdt_t);
    events->SetBranchAddress("chan", &chan_t);
    events->SetBranchAddress("muid", &muid_t);
    events->SetBranchAddress("cnc", &cnc_t);
    
    for (int evt=0; evt<events->GetEntries();evt++) {
        events->GetEntry(evt);
        
        if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
        if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
        if (bdt_t < -1.) continue;
        
        int index = CatMan.index(target_era, chan_t, bdt_t);
        
        if (alter_BDT>-1.) {
            if (bdt_t<alter_BDT) continue;
        }else if (index<0) continue;
        
        if (target_cat!=target_era && index>=0 && target_cat!=CatMan.cats[index].id) continue;
        
        if (tau_t*1E12 < Tau_bound[0] || tau_t*1E12 > Tau_bound[1]) continue;
        if (taue_t*1E12 < TauRes_bound[0] || taue_t*1E12 > TauRes_bound[1]) continue;
        
        if (opt.Contains("even") && (evt%2==1)) continue;
        if (opt.Contains("odd") && (evt%2==0)) continue;
        
        h_global_deltau->Fill((tau_t-gtau_t)*1E12);
        h_global_taue->Fill(taue_t*1E12);
    }
    
    delete events;
    
    double mean_taue = h_global_taue->GetMean();
    cout << ">>> Average tau_err = " << mean_taue << ", to be used in the resolution function." << endl;
    
    RooDataHist *h_deltau_data = new RooDataHist("h_deltau_data", "", RooArgList(delTau), h_global_deltau);
    
    RooRealVar par_taue(Form("respar_%s_taue",tag.Data()),"",mean_taue);
    RooRealVar par_mean(Form("respar_%s_mean",tag.Data()),"",0.0,-0.5,0.5);
    RooRealVar par_sig1(Form("respar_%s_sig1",tag.Data()),"",0.7,0.3,1.0);
    RooRealVar par_sig2(Form("respar_%s_sig2",tag.Data()),"",1.2,0.7,3.5);
    RooRealVar par_sig3(Form("respar_%s_sig3",tag.Data()),"",2.6,1.5,5.5);
    RooRealVar par_frac1(Form("respar_%s_frac1",tag.Data()),"",0.4,0.0,0.6);
    RooRealVar par_frac2(Form("respar_%s_frac2",tag.Data()),"",0.4,0.0,0.6);
    
    RooProduct par_smean(Form("respar_%s_smean",tag.Data()),"",RooArgList(par_taue,par_mean));
    RooProduct par_ssig1(Form("respar_%s_ssig1",tag.Data()),"",RooArgList(par_taue,par_sig1));
    RooProduct par_ssig2(Form("respar_%s_ssig2",tag.Data()),"",RooArgList(par_taue,par_sig2));
    RooProduct par_ssig3(Form("respar_%s_ssig3",tag.Data()),"",RooArgList(par_taue,par_sig3));
    
    RooGaussian model_g1("model_g1","",delTau,par_smean,par_ssig1);
    RooGaussian model_g2("model_g2","",delTau,par_smean,par_ssig2);
    RooGaussian model_g3("model_g3","",delTau,par_smean,par_ssig3);
    
    RooAddPdf *model = 0;
    if (opt.Contains("triple"))
        model = new RooAddPdf("model","",RooArgList(model_g1,model_g2,model_g3),RooArgList(par_frac1,par_frac2));
    if (opt.Contains("double")) {
        par_frac1.setMax(1.0);
        model = new RooAddPdf("model","",RooArgList(model_g1,model_g2),RooArgList(par_frac1));
    }
    
    RooFitResult *res = model->fitTo(*h_deltau_data, Extended(false), NumCPU(NCPU), Hesse(false), Save(true), SumW2Error(true));
    if (res->status()!=0) converge_protection();
    delete res;
    
    RooPlot* frame = delTau.frame(Title(" "));
    
    h_deltau_data->plotOn(frame);
    model->plotOn(frame, LineColor(kBlue), LineWidth(3), Name("model"));
    
    frame->SetMinimum(0.);
    frame->SetMaximum(frame->GetMaximum()*1.4);
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.14,0.06,0.13,0.07);
    
    frame->GetYaxis()->SetTitleOffset(1.48);
    frame->GetYaxis()->SetTitle("Arbitary Unit");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("tau(rec)-tau(gen) [ps]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->Draw();
    
    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.035);
    tex.SetTextAlign(11);
    tex.SetNDC();
    tex.DrawLatex(0.14,0.94,"CMS simulation");
    
    TLegend *leg1 = new TLegend(0.45,0.65,0.93,0.91);
    leg1->SetFillColor(kWhite);
    leg1->SetFillStyle(0);
    leg1->SetLineColor(kWhite);
    leg1->SetLineWidth(0);
    leg1->AddEntry((TObject*)0,Form("<#sigma_{#tau}> = %.4g ps",mean_taue),"");
    leg1->AddEntry((TObject*)0,Form("mean = %.4g #pm %.4g",par_mean.getVal(),par_mean.getError()),"");
    leg1->AddEntry((TObject*)0,Form("sig1 = %.4g #pm %.4g",par_sig1.getVal(),par_sig1.getError()),"");
    leg1->AddEntry((TObject*)0,Form("sig2 = %.4g #pm %.4g",par_sig2.getVal(),par_sig2.getError()),"");
    if (opt.Contains("triple")) leg1->AddEntry((TObject*)0,Form("sig3 = %.4g #pm %.4g",par_sig3.getVal(),par_sig3.getError()),"");
    leg1->AddEntry((TObject*)0,Form("frac1 = %.4g #pm %.4g",par_frac1.getVal(),par_frac1.getError()),"");
    if (opt.Contains("triple")) leg1->AddEntry((TObject*)0,Form("frac2 = %.4g #pm %.4g",par_frac2.getVal(),par_frac2.getError()),"");
    leg1->Draw();
    
    canvas->Print(Form("fig/model_taures_%s.pdf",tag.Data()));
    
    RooRealVar *Tau = wspace->var("Tau");
    
    // fixed the parameters and import the model
    par_mean.setConstant(true);
    par_sig1.setConstant(true);
    par_sig2.setConstant(true);
    par_sig3.setConstant(true);
    par_frac1.setConstant(true);
    par_frac2.setConstant(true);
    RooGaussModel TauRes_g1(Form("TauRes_g1_%s",tag.Data()),"",*Tau,par_smean,par_ssig1);
    RooGaussModel TauRes_g2(Form("TauRes_g2_%s",tag.Data()),"",*Tau,par_smean,par_ssig2);
    RooGaussModel TauRes_g3(Form("TauRes_g3_%s",tag.Data()),"",*Tau,par_smean,par_ssig3);
    
    if (opt.Contains("triple")) {
        RooAddModel TauRes_Model(Form("TauRes_Model_%s",tag.Data()),"",RooArgList(TauRes_g1,TauRes_g2,TauRes_g3),RooArgList(par_frac1,par_frac2));
        wspace->import(TauRes_Model);
    }
    if (opt.Contains("double")) {
        RooAddModel TauRes_Model(Form("TauRes_Model_%s",tag.Data()),"",RooArgList(TauRes_g1,TauRes_g2),RooArgList(par_frac1));
        wspace->import(TauRes_Model);
    }
    
    // Also produce the projection for taue
    h_global_taue->GetYaxis()->SetTitleOffset(1.48);
    h_global_taue->GetYaxis()->SetTitle("Arbitary Unit");
    h_global_taue->GetXaxis()->SetTitleOffset(1.15);
    h_global_taue->GetXaxis()->SetLabelOffset(0.01);
    h_global_taue->GetXaxis()->SetTitle("decay time uncertainty #sigma_{#tau} [ps]");
    h_global_taue->GetXaxis()->SetTitleSize(0.043);
    h_global_taue->GetYaxis()->SetTitleSize(0.043);
    h_global_taue->SetFillColor(50);
    h_global_taue->SetStats(false);
    h_global_taue->Draw("hist");
    
    tex.DrawLatex(0.14,0.94,"CMS simulation");
    
    TLegend *leg2 = new TLegend(0.43,0.80,0.93,0.91);
    leg2->SetFillColor(kWhite);
    leg2->SetFillStyle(0);
    leg2->SetLineColor(kWhite);
    leg2->SetLineWidth(0);
    leg2->AddEntry((TObject*)0,Form("<#sigma_{#tau}> = %.4g ps",mean_taue),"");
    leg2->Draw();
    
    canvas->Print(Form("fig/taures_%s_mc_reco.pdf",tag.Data()));
    
    delete model;
    delete h_global_taue;
    delete h_global_deltau;
    delete h_deltau_data;
    delete leg1;
    delete leg2;
    delete canvas;
}

// ----------------------------------------------------
// Prepare efficiency functions used for lifetime fit
// also produce a histogram PDF for systematic study
//
// available options:
// bsmm      - build B->mumu efficiency model
// bupsik    - build B->J/psi K+ efficiency model
// --
// one_over_exp - modeling by [0]+[1]*x+[2]*x*x+[3]/(1.+exp(-[4]*x))
// threshold    - modeling by threshold function
// --
// taucorr      - apply tau-dependent correction
// --
// even      - only use the even events
// odd       - only use the odd events
//
void PrepareLifetimeEfficiencyModel(RooWorkspace *wspace, TString opt = "bsmm:threshold", TString target_cat = "2016BF", double alter_BDT = -2.)
{
    cout << ">>> PrepareLifetimeEfficiencyModel() start" << endl;
    
    TString tag, treename;
    
    if (opt.Contains("bsmm")) {
        tag = "bsmm";
        treename = "bsmmMc";
    }else if (opt.Contains("bupsik")) {
        tag = "bupsik";
        treename = "bupsikMc";
    }else {
        cout << ">>> Undefined option: " << opt << endl;
        return;
    }
    tag += Form("_%s", target_cat.Data());
    
    TH1D *h_taucorr = 0;
    if (opt.Contains("taucorr")) {
        h_taucorr = (TH1D *)wspace->obj(Form("h_taucorr_%s",tag.Data()));
        tag += "_taucorr";
    }
    
    TString target_era = target_cat;
    for (auto& cat: CatMan.cats)
        if (cat.id == target_cat) target_era = cat.era;

    vector<double> xbins = {
        0.5,0.625,0.75,0.875,
        1.,1.125,1.25,1.375,1.5,1.625,1.75,1.875,
        2.,2.125,2.25,2.375,2.5,2.75,
        3.,3.25,3.5,3.75,
        4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.};
    
    TH1D *h_global_taueff = new TH1D("h_global_taueff","",xbins.size()-1,xbins.data());
    h_global_taueff->Sumw2();
    TH1D *h_taureco_mc = new TH1D(Form("h_taureco_mc_%s",tag.Data()),"",xbins.size()-1,xbins.data());
    h_taureco_mc->Sumw2();
    
    // load the corresponding MC sample
    TString filename = Form("input/bmm4/small%s-%s.root",target_era.Data(),treename.Data());
    TChain *events = new TChain(treename);
    exist_protection(filename);
    events->Add(filename);
    cout << ">>> Loading from " << filename << ", with " << events->GetEntries() << " entries." << endl;
    
    double m_t, me_t, pt_t, tau_t, taue_t, gtau_t, bdt_t;
    int chan_t;
    bool muid_t, cnc_t;
    
    events->SetBranchAddress("m", &m_t);
    events->SetBranchAddress("me", &me_t);
    events->SetBranchAddress("pt", &pt_t);
    events->SetBranchAddress("tau", &tau_t);
    events->SetBranchAddress("taue", &taue_t);
    events->SetBranchAddress("gtau", &gtau_t);
    events->SetBranchAddress("bdt", &bdt_t);
    events->SetBranchAddress("chan", &chan_t);
    events->SetBranchAddress("muid", &muid_t);
    events->SetBranchAddress("cnc", &cnc_t);
    
    for (int evt=0; evt<events->GetEntries();evt++) {
        events->GetEntry(evt);
        
        if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
        if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
        if (bdt_t < -1.) continue;
        
        int index = CatMan.index(target_era, chan_t, bdt_t);
        
        if (alter_BDT>-1.) {
            if (bdt_t<alter_BDT) continue;
        }else if (index<0) continue;
        
        if (target_cat!=target_era && index>=0 && target_cat!=CatMan.cats[index].id) continue;
        
        //if (tau_t*1E12 < Tau_bound[0] || tau_t*1E12 > Tau_bound[1]) continue;
        if (tau_t*1E12 < 0.5 || tau_t*1E12 > Tau_bound[1]) continue; // allow a litte bit more to the lower side
        if (taue_t*1E12 < TauRes_bound[0] || taue_t*1E12 > TauRes_bound[1]) continue;
        
        if (opt.Contains("even") && (evt%2==1)) continue;
        if (opt.Contains("odd") && (evt%2==0)) continue;
        
        h_global_taueff->Fill(gtau_t*1E12);
        h_taureco_mc->Fill(tau_t*1E12);
    }
    
    delete events;
    
    TH1D *h_taueff_norm = new TH1D("h_taueff_norm","",xbins.size()-1,xbins.data());
    h_taueff_norm->Sumw2();
    
    // extract lifetime & resolution
    RooRealVar *Tau = wspace->var("Tau");
    //RooResolutionModel *TauRes_Model = (RooResolutionModel*)wspace->obj(Form("TauRes_Model_%s",tag.Data()));
    RooTruthModel *TauRes_Model = new RooTruthModel("TauRes_Model","",*Tau); // ideal model
    RooRealVar EffTau("EffTau","",1.6);
    RooDecay RawDecay("RawDecay","",*Tau,EffTau,*TauRes_Model,RooDecay::SingleSided);

    if (tag.Contains("bsmm")) EffTau.setVal(1.472);
    if (tag.Contains("bupsik") && (target_era.Contains("2011") || target_era.Contains("2012"))) EffTau.setVal(1.671);
    if (tag.Contains("bupsik") && target_era.Contains("2016"))  EffTau.setVal(1.638);
    
    double scale_factor = 0.;
    for (int bin=1; bin<=h_taueff_norm->GetNbinsX(); bin++) {
        
        Tau->setRange("bin",h_taueff_norm->GetBinLowEdge(bin),h_taueff_norm->GetBinLowEdge(bin)+h_taueff_norm->GetBinWidth(bin));
        RooAbsReal* area = RawDecay.createIntegral(*Tau,Range("bin"));
        h_taueff_norm->SetBinContent(bin,area->getVal());
        delete area;
        
        double ratio = h_global_taueff->GetBinContent(bin)/h_taueff_norm->GetBinContent(bin);
        if (ratio>scale_factor) scale_factor = ratio;
    }
    h_taueff_norm->Scale(scale_factor);
    
    delete TauRes_Model;
    
    for (int bin=1; bin<=h_global_taueff->GetNbinsX(); bin++) {
        double val = h_global_taueff->GetBinContent(bin)/h_taueff_norm->GetBinContent(bin);
        double err = h_global_taueff->GetBinError(bin)/h_taueff_norm->GetBinContent(bin);
        h_global_taueff->SetBinContent(bin,val);
        h_global_taueff->SetBinError(bin,err);
    }
    
    if (h_taucorr!=0) h_global_taueff->Multiply(h_taucorr);
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.14,0.06,0.13,0.07);
    
    h_global_taueff->GetYaxis()->SetTitleOffset(1.48);
    h_global_taueff->GetYaxis()->SetTitle("Arbitary Unit");
    h_global_taueff->GetXaxis()->SetTitleOffset(1.15);
    h_global_taueff->GetXaxis()->SetLabelOffset(0.01);
    h_global_taueff->GetXaxis()->SetTitle("#tau [ps]");
    h_global_taueff->GetXaxis()->SetTitleSize(0.043);
    h_global_taueff->GetYaxis()->SetTitleSize(0.043);
    h_global_taueff->SetStats(false);
    h_global_taueff->SetLineColor(kBlack);
    h_global_taueff->SetMarkerStyle(20);
    h_global_taueff->SetMaximum(h_global_taueff->GetMaximum()*1.2);
    
    TF1 *effmodel = 0;
    // Chandi's model
    if (opt.Contains("one_over_exp")) {
        effmodel = new TF1("effmodel","[0]+[1]*x+[2]*x*x+[3]/(1.+exp(-[4]*x))",xbins.front(),xbins.back());
        effmodel->SetParameters(0.12, 0.006, 0., 0.2, 1.2);
    }
    // Threshold func.
    if (opt.Contains("threshold")) {
        effmodel = new TF1("effmodel","[0]+[1]*pow(x,[2])*exp([3]*x+[4]*x*x)",xbins.front(),xbins.back());
        effmodel->SetParameters(-4., 4., 0.3, 0., 0.);
    }
    
    h_global_taueff->Fit("effmodel","VI","", 0.875, Tau_bound[1]); // 1st fit
    h_global_taueff->Fit("effmodel","VI","", 0.875, Tau_bound[1]); // 2nd fit
    int fitstat = h_global_taueff->Fit("effmodel","VI","", 0.875, Tau_bound[1]); // 3rd fit
    cout << ">>> Fit status = " << fitstat << endl;
    if (fitstat!=0 && fitstat!=4000) converge_protection();
    
    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.035);
    tex.SetTextAlign(11);
    tex.SetNDC();
    tex.DrawLatex(0.14,0.94,"CMS simulation");
    
    TLegend *leg1 = new TLegend(0.52,0.60,0.93,0.91);
    if (opt.Contains("one_over_exp")) leg1->SetHeader("Fit to f(t) = k_{1} + k_{2}t + k_{3}t^{2} + #frac{k_{4}}{1+exp(-k_{5}t)}","C");
    if (opt.Contains("threshold")) leg1->SetHeader("Fit to f(t) = k_{1} + k_{2}*t^{k_{3}}exp(k_{4}t+k_{5}t^{2})]","C");
    //leg1->AddEntry((TObject*)0,"","");
    leg1->SetFillColor(kWhite);
    leg1->SetFillStyle(0);
    leg1->SetLineColor(kWhite);
    leg1->SetLineWidth(0);
    for (int i=0; i<5; i++)
        leg1->AddEntry((TObject*)0,Form("k_{%d} = %.4g #pm %.4g",i+1,effmodel->GetParameter(i),effmodel->GetParError(i)),"");
    leg1->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3g/%d",effmodel->GetChisquare(),effmodel->GetNDF()),"");
    leg1->Draw();
    
    canvas->Print(Form("fig/model_globaleff_%s.pdf",tag.Data()));
    
    RooRealVar *par[5];
    for (int i=0; i<5; i++) {
        par[i] = new RooRealVar(Form("effpar_%s_k%d",tag.Data(),i+1),"",effmodel->GetParameter(i));
        par[i]->setError(effmodel->GetParError(i));
    }
    TString formula = "";
    if (opt.Contains("one_over_exp")) {
        formula += Form("max(effpar_%s_k1",tag.Data());
        formula += Form("+effpar_%s_k2*Tau",tag.Data());
        formula += Form("+effpar_%s_k3*Tau*Tau",tag.Data());
        formula += Form("+effpar_%s_k4/(1.+exp(-Tau*effpar_%s_k5)),1E-5)",tag.Data(),tag.Data());
    }
    
    // the factor 0.5 is to protect the "efficiency" should be always below 1
    if (opt.Contains("threshold")) {
        formula += Form("max(effpar_%s_k1+",tag.Data());
        formula += Form("effpar_%s_k2*",tag.Data());
        formula += Form("pow(Tau,effpar_%s_k3)*",tag.Data());
        formula += Form("exp(effpar_%s_k4*Tau+effpar_%s_k5*Tau*Tau),1E-5)*0.5",tag.Data(),tag.Data());
    }
    
    RooArgList varlist;
    varlist.add(*Tau);
    for (int i=0; i<5; i++) varlist.add(*par[i]);
    RooFormulaVar TauEff_Model(Form("TauEff_Model_%s",tag.Data()),formula,varlist);
    wspace->import(TauEff_Model);
    
    // alternative histogram model
    TH1D *h_tau_tmp = new TH1D("h_tau_tmp","",(int)((Tau_bound[1]-Tau_bound[0])*8.),Tau_bound[0],Tau_bound[1]); // convert to fixed bin width histogram
    for (int i=1; i<=h_tau_tmp->GetNbinsX(); i++) {
        for (int bin=1; bin<=h_global_taueff->GetNbinsX(); bin++) {
            double x = h_tau_tmp->GetBinCenter(i);
            double min = h_global_taueff->GetBinLowEdge(bin);
            double max = min+h_global_taueff->GetBinWidth(bin);
            if (x>=min && x<max) h_tau_tmp->SetBinContent(i,h_global_taueff->GetBinContent(bin));
        }
    }
    RooDataHist h_global_taueff_data(Form("h_global_taueff_data_%s",tag.Data()), "", *Tau, h_tau_tmp);
    RooHistPdf TauEff_Model_Hist(Form("TauEff_Model_Hist_%s",tag.Data()),"",*Tau,h_global_taueff_data,0);
    wspace->import(TauEff_Model_Hist);
    
    RooPlot* frame = Tau->frame(Title(" "));
    TauEff_Model_Hist.plotOn(frame, DrawOption("L"), LineColor(kBlue), LineWidth(2), LineStyle(1), NumCPU(NCPU));
    
    frame->GetYaxis()->SetTitleOffset(1.48);
    frame->GetYaxis()->SetTitle("Arbitary Unit");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("#tau [ps]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->SetStats(false);
    frame->Draw();
    
    canvas->Print(Form("fig/model_globaleff_histpdf_%s.pdf",tag.Data()));
    
    if (!opt.Contains("taucorr")) wspace->import(*h_taureco_mc);
    delete h_taureco_mc;
    
    for (int i=0; i<5; i++) delete par[i];
    delete h_global_taueff;
    delete effmodel;
    delete leg1;
    delete canvas;
    delete frame;
}

// ----------------------------------------------------
// available options:
// bsmm      - test with B->mumu MC
// bupsik    - test with B->J/psi K+ MC
// --
// genonly   - fit to pure generated decay time
// geneff    - fit to generated decay time, with efficiency correction
// reco      - fit to reco decay time
// toybkg    - mixing with toy combinatorial background
//
void PerformMCTauSPlotStudy(RooWorkspace *wspace, TString opt = "bsmm:reco", TString target_cat = "2016BF", double alter_BDT = -2.)
{
    cout << ">>> PerformMCTauSPlotStudy() start" << endl;
    
    TString tag, filename, treename, filename_eff;
    
    if (opt.Contains("bsmm")) {
        tag = "bsmm";
        treename = "bsmmMc";
    }else if (opt.Contains("bupsik")) {
        tag = "bupsik";
        treename = "bupsikMc";
    }else {
        cout << ">>> Undefined option: " << opt << endl;
        return;
    }
    tag += Form("_%s", target_cat.Data());
    
    TString target_era = target_cat;
    for (auto& cat: CatMan.cats)
        if (cat.id == target_cat) target_era = cat.era;
    
    filename = Form("input/bmm4/small%s-%s.root",target_era.Data(),treename.Data());
    filename_eff = Form("input/bmm4/eff%s-%s.root",target_era.Data(),treename.Data());
    
    TChain *events = new TChain(treename);
    exist_protection(filename);
    events->Add(filename);
    
    TChain *effTree = new TChain("effTree");
    //exist_protection(filename_eff);
    effTree->Add(filename_eff);
    
    double m_t, me_t, pt_t, tau_t, taue_t, gtau_t, bdt_t;
    int chan_t;
    bool muid_t, cnc_t;
    
    events->SetBranchAddress("m", &m_t);
    events->SetBranchAddress("me", &me_t);
    events->SetBranchAddress("pt", &pt_t);
    events->SetBranchAddress("tau", &tau_t);
    events->SetBranchAddress("taue", &taue_t);
    events->SetBranchAddress("gtau", &gtau_t);
    events->SetBranchAddress("bdt", &bdt_t);
    events->SetBranchAddress("chan", &chan_t);
    events->SetBranchAddress("muid", &muid_t);
    events->SetBranchAddress("cnc", &cnc_t);
    
    float effgtau_t;
    
    effTree->SetBranchAddress("gtau", &effgtau_t);
    
    // EffTau model
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *Tau = wspace->var("Tau");
    
    RooResolutionModel *TauRes_Model = (RooResolutionModel*)wspace->obj(Form("TauRes_Model_%s",tag.Data()));
    RooFormulaVar *TauEff_Model = (RooFormulaVar*)wspace->obj(Form("TauEff_Model_%s",tag.Data()));
    
    RooRealVar EffTau("EffTau","",1.6,0.1,5.0);
    RooDecay RawDecay_reco("RawDecay_reco","",*Tau,EffTau,*TauRes_Model,RooDecay::SingleSided);
    RooEffProd Tau_pdf_reco("Tau_pdf_reco","",RawDecay_reco,*TauEff_Model);

    RooTruthModel TauRes_TruthModel("TauRes_TruthModel","",*Tau); // ideal model
    RooDecay Tau_pdf_genonly("Tau_pdf_genonly","",*Tau,EffTau,TauRes_TruthModel,RooDecay::SingleSided);
    
    RooEffProd Tau_pdf_geneff("Tau_pdf_geneff","",Tau_pdf_genonly,*TauEff_Model);
    
    RooAbsPdf *Tau_pdf = &Tau_pdf_reco;
    
    if (opt.Contains("geneff")) Tau_pdf = &Tau_pdf_geneff;
    
    // place holder for sPlot
    //TH1D *h_tau = new TH1D("h_tau", "", Tau_bins.size()-1, Tau_bins.data());
    TH1D *h_tau = new TH1D("h_tau", "", (int)((Tau_bound[1]-Tau_bound[0])*10.), Tau_bound[0], Tau_bound[1]);
    h_tau->Sumw2();
    
    if (opt.Contains("genonly")) { // generater tau
        for (int evt=0; evt<effTree->GetEntries();evt++) {
            effTree->GetEntry(evt);
            h_tau->Fill(effgtau_t*1E12);
        }
        Tau_pdf = &Tau_pdf_genonly;
    }else {
        RooDataSet *rds = new RooDataSet("rds","",RooArgSet(*Tau));
    
        for (int evt=0; evt<events->GetEntries();evt++) {
            events->GetEntry(evt);
            
            if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
            if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
            if (bdt_t < -1.) continue;
            
            int index = CatMan.index(target_era, chan_t, bdt_t);
            
            if (alter_BDT>-1.) {
                if (bdt_t<alter_BDT) continue;
            }else if (index<0) continue;
            
            if (target_cat!=target_era && index>=0 && target_cat!=CatMan.cats[index].id) continue;
            
            // tau cuts
            if (tau_t*1E12 < Tau_bound[0] || tau_t*1E12 > Tau_bound[1]) continue;
            if (taue_t*1E12 < TauRes_bound[0] || taue_t*1E12 > TauRes_bound[1]) continue;
            
            if (opt.Contains("geneff")) Tau->setVal(gtau_t*1E12);
            else Tau->setVal(tau_t*1E12);
            rds->add(RooArgSet(*Tau));
        }
        cout << ">>> " << rds->numEntries() << " MC events to be included in the fit" << endl;
        
        // Add toy background
        if (opt.Contains("toybkg")) {
        
            int ns = rds->numEntries();
            int nb = rds->numEntries()*50; // assume a 50x background
            double mass_mean_sig = 5.35;
            if (opt.Contains("bupsik")) {
                nb = rds->numEntries()*5; // assume a 5x background for J/psi K+
                mass_mean_sig = 5.28;
            }
        
            RooGaussian Mass_pdf_sig("Mass_pdf_sig", "", *Mass,RooConst(mass_mean_sig),RooConst(0.035));
            RooDataSet *rds_toymass = Mass_pdf_sig.generate(RooArgSet(*Mass),ns);
            rds->merge(rds_toymass);
            delete rds_toymass;
        
            RooChebychev Mass_pdf_comb("Mass_pdf_comb", "", *Mass, RooArgList(RooConst(0.)));
            RooRealVar TauComb("TauComb","",1.2); // assume an 1.2 ps lifetime for combintorial PDF
            RooDecay RawDecay_comb("RawDecay_comb","",*Tau,TauComb,*TauRes_Model,RooDecay::SingleSided);
            RooEffProd Tau_pdf_comb("Tau_pdf_comb","",RawDecay_comb,*TauEff_Model);
        
            RooProdPdf pdf_comb("pdf_comb","",Mass_pdf_comb,Tau_pdf_comb);
            RooDataSet *rds_toybkg = pdf_comb.generate(RooArgSet(*Mass,*Tau),nb);
            rds->append(*rds_toybkg);
            delete rds_toybkg;
            
            // Make a projection for toy mass distribution
            RooAddPdf Mass_pdf("Mass_pdf","",
                               RooArgList(Mass_pdf_sig, Mass_pdf_comb),
                               RooArgList(RooConst(ns), RooConst(nb)));
            
            RooPlot* frame_m = Mass->frame();
            
            rds->plotOn(frame_m,Binning(160));
            Mass_pdf.plotOn(frame_m);
            Mass_pdf.plotOn(frame_m,Components("Mass_pdf_sig"),LineColor(kRed-10),LineWidth(2),LineStyle(kSolid),FillColor(kRed-10), VLines(), DrawOption("F"));
            Mass_pdf.plotOn(frame_m,Components("Mass_pdf_comb"),LineColor(kCyan+2),LineWidth(3),LineStyle(7));

            TCanvas* canvas_m = new TCanvas("canvas_m", "", 600, 600);
            canvas_m->SetMargin(0.14,0.06,0.13,0.07);
            
            frame_m->GetYaxis()->SetTitleOffset(1.48);
            frame_m->GetYaxis()->SetTitle("Entries");
            frame_m->GetXaxis()->SetTitleOffset(1.15);
            frame_m->GetXaxis()->SetLabelOffset(0.01);
            frame_m->GetXaxis()->SetTitle("Mass [GeV]");
            frame_m->GetXaxis()->SetTitleSize(0.043);
            frame_m->GetYaxis()->SetTitleSize(0.043);
            frame_m->SetStats(false);
            frame_m->SetTitle("");
            
            frame_m->Draw();
            canvas_m->Print(Form("fig/proj_mass_%s_mcstudy_toybkg.pdf",tag.Data()));
            
            delete frame_m;
            delete canvas_m;
            
            // Now calculate sWeights
            enum {_sig, _comb, _nspec};
            
            double yield[_nspec];
            yield[_sig] = ns;
            yield[_comb] = nb;
            
            TMatrixD covInv(_nspec, _nspec);
            covInv = 0.;
            RooArgSet norm(*Mass);
            
            for (int evt=0; evt<rds->numEntries(); evt++) {
                const RooArgSet* arg = rds->get(evt);
                Mass->setVal(arg->getRealValue("Mass"));
                Tau->setVal(arg->getRealValue("Tau"));
            
                double pdf[_nspec];
                pdf[_sig] = Mass_pdf_sig.getVal(&norm);
                pdf[_comb] = Mass_pdf_comb.getVal(&norm);
            
                double pdf_total = 0.;
                for (int idx = 0; idx<_nspec; idx++) pdf_total += yield[idx]*pdf[idx];
            
                for (int row = 0; row<_nspec; row++)
                    for (int col = 0; col<_nspec; col++)
                        covInv(row,col) += pdf[row]*pdf[col]/(pdf_total*pdf_total);
            }
        
            TMatrixD covMatrix(TMatrixD::kInverted,covInv);
        
            for (int evt=0; evt<rds->numEntries(); evt++) {
                const RooArgSet* arg = rds->get(evt);
                Mass->setVal(arg->getRealValue("Mass"));
                Tau->setVal(arg->getRealValue("Tau"));
            
                double pdf[_nspec];
                pdf[_sig] = Mass_pdf_sig.getVal(&norm);
                pdf[_comb] = Mass_pdf_comb.getVal(&norm);
            
                double denominator = 0.;
                for (int idx = 0; idx<_nspec; idx++) denominator += yield[idx]*pdf[idx];
            
                double numerator = 0.;
                for (int idx = 0; idx<_nspec; idx++) numerator += covMatrix(_sig,idx)*pdf[idx];
            
                double weight = numerator/denominator;
                h_tau->Fill(arg->getRealValue("Tau"),weight);
            }
        }else { // signal RECO only
            for (int evt=0; evt<rds->numEntries(); evt++) {
                const RooArgSet* arg = rds->get(evt);
                h_tau->Fill(arg->getRealValue("Tau"));
            }
        }
        
        delete rds;
    }
    
    // Binned weighted likelihood fit w/ PDF bin integration
    RooFitResult *res1 = NULL, *res2 = NULL;
    Fit_sPlot(h_tau, Tau_pdf, Tau, &EffTau, &res1, &res2);
    
    // Projection
    RooPlot* frame = Tau->frame(Title(" "));
    h_tau->SetMarkerStyle(1);
    h_tau->SetLineColor(kBlack);
    RooDataHist *h_tau_data = new RooDataHist("h_tau_data", "", RooArgList(*Tau), h_tau);
    h_tau_data->plotOn(frame,MarkerStyle(1),LineColor(kBlack));
    
    Tau_pdf->plotOn(frame, DrawOption("L"), LineColor(kBlue), LineWidth(2), LineStyle(1), NumCPU(NCPU));
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.14,0.06,0.13,0.07);
    
    frame->GetYaxis()->SetTitleOffset(1.48);
    frame->GetYaxis()->SetTitle("Entries");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("Decay time [ps]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->SetStats(false);
    frame->Draw("E");
    
    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.035);
    tex.SetTextAlign(11);
    tex.SetNDC();
    tex.DrawLatex(0.14,0.94,"CMS Preliminary");
    
    TLine lin;
    lin.SetLineColor(kGray+1);
    lin.SetLineWidth(2);
    lin.SetLineStyle(7);
    lin.DrawLine(1.,0.,12.,0.);
    
    canvas->Update();
    
    TLegend *leg1 = new TLegend(0.50,0.86,0.91,0.91);
    leg1->SetNColumns(1);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(kWhite);
    if (opt.Contains("bsmm")) leg1->AddEntry(h_tau, "B_{s} #rightarrow #mu^{+}#mu^{-} MC", "lep");
    if (opt.Contains("bupsik")) leg1->AddEntry(h_tau, "B^{+} #rightarrow J/#psi K^{+} MC", "lep");
    leg1->Draw();
    
    TString pdf_filename = Form("fig/tau_splot_%s_mc_reco.pdf",tag.Data());
    if (opt.Contains("genonly")) pdf_filename = Form("fig/tau_splot_%s_mc_genonly.pdf",tag.Data());
    if (opt.Contains("geneff"))  pdf_filename = Form("fig/tau_splot_%s_mc_geneff.pdf",tag.Data());
    if (opt.Contains("toybkg"))  pdf_filename = Form("fig/tau_splot_%s_mc_toybkg.pdf",tag.Data());
    
    canvas->Print(pdf_filename);
    
    delete h_tau_data;
    delete h_tau;
    delete leg1;
    delete frame;
    delete canvas;
    delete res1;
    delete res2;
    delete events;
    delete effTree;
}

// ----------------------------------------------------
// available options:
// bupsik    - test with B->J/psi K+ data
// --
// taucorr   - apply tau-dependent correction
// --
// even      - only use the even events
// odd       - only use the odd events
//
void PerformBuDataTauSPlotStudy(RooWorkspace *wspace, TString opt = "bupsik", TString target_cat = "2016BF", double alter_BDT = -2.)
{
    cout << ">>> PerformBuDataTauSPlotStudy() start" << endl;
    
    TString tag, treename;
    
    if (opt.Contains("bupsik")) {
        tag = "bupsik";
        treename = "bupsikData";
    }else {
        cout << ">>> Undefined option: " << opt << endl;
        return;
    }
    tag += Form("_%s", target_cat.Data());
    
    RooResolutionModel *TauRes_Model = (RooResolutionModel*)wspace->obj(Form("TauRes_Model_%s",tag.Data()));
    RooFormulaVar *TauEff_Model = (RooFormulaVar*)wspace->obj(Form("TauEff_Model_%s",tag.Data()));
    
    if (opt.Contains("taucorr")) {
        tag += "_taucorr";
        TauEff_Model = (RooFormulaVar*)wspace->obj(Form("TauEff_Model_%s",tag.Data()));
    }
    
    TString target_era = target_cat;
    for (auto& cat: CatMan.cats)
        if (cat.id == target_cat) target_era = cat.era;
    
    // extract observables
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *Tau = wspace->var("Tau");
    RooRealVar *TauRes = wspace->var("TauRes");
    RooRealVar *Weight = wspace->var("Weight");
    
    // Tigher boundary for Mass and extend the boundary for Tau
    Mass->setMin(5.0);
    Mass->setMax(5.8);
    Tau->setMin(0.5);
    
    RooRealVar EffTau("EffTau","",1.6,0.1,5.0);
    RooDecay RawDecay_reco("RawDecay_reco","",*Tau,EffTau,*TauRes_Model,RooDecay::SingleSided);
    RooEffProd Tau_pdf_reco("Tau_pdf_reco","",RawDecay_reco,*TauEff_Model);
    
    // place holder for sPlot
    //TH1D *h_tau = new TH1D("h_tau", "", Tau_bins.size()-1, Tau_bins.data());
    TH1D *h_tau = new TH1D("h_tau", "", (int)((Tau_bound[1]-Tau_bound[0])*10.), Tau_bound[0], Tau_bound[1]);
    h_tau->Sumw2();
    TH1D *h_taures = new TH1D("h_taures", "", TauRes_bins.size()-1, TauRes_bins.data());
    h_taures->Sumw2();
    
    vector<double> xbins = {
        0.5,0.625,0.75,0.875,
        1.,1.125,1.25,1.375,1.5,1.625,1.75,1.875,
        2.,2.125,2.25,2.375,2.5,2.75,
        3.,3.25,3.5,3.75,
        4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.};
    TH1D *h_taureco = new TH1D(Form("h_taureco_%s",tag.Data()),"",xbins.size()-1,xbins.data());
    h_taureco->Sumw2();
    
    RooDataSet *rds = new RooDataSet("rds","",RooArgSet(*Mass,*Tau,*TauRes,*Weight));

    TString filename;
    if (target_era.Contains("2016")) filename = Form("input/bmm4/small2016*-%s.root",treename.Data());
    else filename = Form("input/bmm4/small%s-%s.root",target_era.Data(),treename.Data());

    TChain *events = new TChain(treename);
    //exist_protection(filename);
    events->Add(filename);
    
    double m_t, me_t, pt_t, tau_t, taue_t, gtau_t, bdt_t, cw8_t;
    int chan_t, run_t;
    bool muid_t, cnc_t;
    
    events->SetBranchAddress("run", &run_t);
    events->SetBranchAddress("m", &m_t);
    events->SetBranchAddress("me", &me_t);
    events->SetBranchAddress("pt", &pt_t);
    events->SetBranchAddress("tau", &tau_t);
    events->SetBranchAddress("taue", &taue_t);
    events->SetBranchAddress("gtau", &gtau_t);
    events->SetBranchAddress("bdt", &bdt_t);
    events->SetBranchAddress("cw8", &cw8_t);
    events->SetBranchAddress("chan", &chan_t);
    events->SetBranchAddress("muid", &muid_t);
    events->SetBranchAddress("cnc", &cnc_t);
    
    cout << ">>> Loading from " << filename << ", with " << events->GetEntries() << " entries." << endl;
    
    for (int evt=0; evt<events->GetEntries();evt++) {
        events->GetEntry(evt);
        
        if (m_t < 5.0 || m_t > 5.8 || !muid_t) continue;
        if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
        if (bdt_t < -1.) continue;
        
        int index = CatMan.index(target_era, chan_t, bdt_t);
        
        if (alter_BDT>-1.) {
            if (bdt_t<alter_BDT) continue;
        }else if (index<0) continue;
        
        if (target_cat!=target_era && index>=0 && target_cat!=CatMan.cats[index].id) continue;
        
        // tau cuts
        //if (tau_t*1E12 < Tau_bound[0] || tau_t*1E12 > Tau_bound[1]) continue;
        if (tau_t*1E12 < 0.5 || tau_t*1E12 > Tau_bound[1]) continue; // allow a litte bit more to the lower side
        if (taue_t*1E12 < TauRes_bound[0] || taue_t*1E12 > TauRes_bound[1]) continue;
        
        if (opt.Contains("even") && (evt%2==1)) continue;
        if (opt.Contains("odd") && (evt%2==0)) continue;
        
        if (target_era=="2016B" && (run_t<273150 || run_t>275376)) continue;
        if (target_era=="2016C" && (run_t<275657 || run_t>276283)) continue;
        if (target_era=="2016D" && (run_t<276315 || run_t>276811)) continue;
        if (target_era=="2016E" && (run_t<276831 || run_t>277420)) continue;
        if (target_era=="2016F" && (run_t<277981 || run_t>278808)) continue;
        if (target_era=="2016G" && (run_t<278820 || run_t>280385)) continue;
        if (target_era=="2016H" && (run_t<281613 || run_t>284044)) continue;
        if (target_era=="2016BF" && (run_t<273150 || run_t>278808)) continue;
        if (target_era=="2016GH" && (run_t<278820 || run_t>284044)) continue;
        
        Mass->setVal(m_t); // unconstrained mass
        Tau->setVal(tau_t*1E12);
        TauRes->setVal(taue_t*1E12);
        Weight->setVal(cw8_t);
        rds->add(RooArgSet(*Mass,*Tau,*TauRes,*Weight));
    }
    
    delete events;
    cout << ">>> " << rds->numEntries() << " events to be included in the fit" << endl;
    
    // B+ mass signal model
    RooRealVar Mass_mean("Mass_mean","",5.27926,5.0,5.8);
    RooRealVar Mass_sigma("Mass_sigma","",0.012,0.001,0.050);
    
    RooRealVar Mass_g2scale("Mass_g2scale","",1.8,1.0,8.0);
    RooRealVar Mass_g2fraction("Mass_g2fraction","",0.30,0.0,0.6);
    RooRealVar Mass_g2shift("Mass_g2shift","",-0.07,-2.5,+2.5);

    RooProduct Mass_g2sigma("Mass_g2sigma","",RooArgList(Mass_sigma,Mass_g2scale));
    RooProduct Mass_g2deltam("Mass_g2deltam","",RooArgList(Mass_sigma,Mass_g2shift));
    RooAddition Mass_g2mean("Mass_g2mean","",RooArgList(Mass_mean,Mass_g2deltam));
    
    RooGaussian Mass_g1("Mass_g1","",*Mass,Mass_mean,Mass_sigma);
    RooGaussian Mass_g2("Mass_g2","",*Mass,Mass_g2mean,Mass_g2sigma);
    
    RooAddPdf Mass_pdf_sig("Mass_pdf_sig","",RooArgList(Mass_g2,Mass_g1),RooArgList(Mass_g2fraction));
    
    // Combintorial model
    RooRealVar Mass_exp("Mass_exp","",-4.0,-20.,+20.);
    RooExponential Mass_pdf_comb("Mass_pdf_comb","",*Mass,Mass_exp);
    
    // J/psi pi+ model (from BPH-15-004 setup)
    // w/o mass constraint
    // 1  m_jpsipi_fraction2   2.54323e-01   2.22693e-02  -2.38567e-02   2.14143e-02
    // 2  m_jpsipi_fraction3   1.26181e-01   1.07897e-02  -1.10564e-02   1.05980e-02
    // 3  m_jpsipi_mean1   5.35555e+00   1.34882e-03  -1.34331e-03   1.35031e-03
    // 4  m_jpsipi_mean2   5.46283e+00   6.46652e-03  -5.95385e-03   7.31768e-03
    // 5  m_jpsipi_mean3   5.47169e+00   7.98208e-03  -8.51978e-03   7.77263e-03
    // 6  m_jpsipi_sigma1l   4.30865e-02   1.13666e-03  -1.12339e-03   1.14618e-03
    // 7  m_jpsipi_sigma1r   6.36127e-02   2.61743e-03  -2.60383e-03   2.63277e-03
    // 8  m_jpsipi_sigma2   9.75406e-02   5.12110e-03  -5.12012e-03   5.25654e-03
    // 9  m_jpsipi_sigma3   3.16814e-01   2.45739e-02  -2.13032e-02   2.95821e-02*/
    
    RooRealVar JpsiPi_mean1("JpsiPi_mean1","",5.35555e+00);
    RooRealVar JpsiPi_mean2("JpsiPi_mean2","",5.46283e+00);
    RooRealVar JpsiPi_mean3("JpsiPi_mean3","",5.47169e+00);
    RooRealVar JpsiPi_sigma1l("JpsiPi_sigma1l","",4.30865e-02);
    RooRealVar JpsiPi_sigma1r("JpsiPi_sigma1r","",6.36127e-02);
    RooRealVar JpsiPi_sigma2("JpsiPi_sigma2","",9.75406e-02);
    RooRealVar JpsiPi_sigma3("JpsiPi_sigma3","",3.16814e-01);
    RooRealVar JpsiPi_fraction2("JpsiPi_fraction2","",2.54323e-01);
    RooRealVar JpsiPi_fraction3("JpsiPi_fraction3","",1.26181e-01);
    
    RooBifurGauss JpsiPi_g1("JpsiPi_g1","",*Mass,JpsiPi_mean1,JpsiPi_sigma1l,JpsiPi_sigma1r);
    RooGaussian JpsiPi_g2("JpsiPi_g2","",*Mass,JpsiPi_mean2,JpsiPi_sigma2);
    RooGaussian JpsiPi_g3("JpsiPi_g3","",*Mass,JpsiPi_mean3,JpsiPi_sigma3);
    
    RooAddPdf Mass_pdf_jpsipi("Mass_pdf_jpsipi","",RooArgList(JpsiPi_g3,JpsiPi_g2,JpsiPi_g1),RooArgList(JpsiPi_fraction3,JpsiPi_fraction2));

    // Non-propmpt model
    RooRealVar NonPrompt_scale("NonPrompt_scale","",2.2e-02,0.001,0.08);
    RooRealVar NonPrompt_shift("NonPrompt_shift","",5.14,5.12,5.16);
    RooGenericPdf Mass_pdf_nonprompt("Mass_pdf_nonprompt","","TMath::Erfc((Mass-NonPrompt_shift)/NonPrompt_scale)",
                                     RooArgList(*Mass,NonPrompt_scale,NonPrompt_shift));
    
    RooRealVar n_sig("n_sig","",3.0E5,0.,1E6);
    RooRealVar n_comb("n_comb","",5E4,0.,1E6);
    RooRealVar n_nonprompt("n_nonprompt","",4E4,0.,1E6);
    
    RooRealVar f_jpsipi("f_jpsipi","",4.1E-5/1.026E-3);  // BF(Jpsi pi) = (4.1±0.4)×10−5 / BF(Jpsi K) = (1.026±0.031)×10−3
    RooProduct n_jpsipi("n_jpsipi","",RooArgList(n_sig,f_jpsipi));
    
    // Full mass model
    RooAddPdf Mass_pdf("Mass_pdf","",
                    RooArgList(Mass_pdf_sig, Mass_pdf_comb, Mass_pdf_jpsipi, Mass_pdf_nonprompt),
                    RooArgList(n_sig, n_comb, n_jpsipi, n_nonprompt));
    
    TH1D* histo_data = (TH1D*)rds->createHistogram("histo_data", *Mass, Binning(320,5.0,5.8));
    RooDataHist *rds_hist = new RooDataHist("rds_hist", "", RooArgList(*Mass), histo_data);
    
    Mass_g2shift.setConstant(true);
    Mass_g2scale.setConstant(true);
    Mass_g2fraction.setConstant(true);
    Mass_pdf.fitTo(*rds_hist);
    Mass_g2shift.setConstant(false);
    Mass_g2scale.setConstant(false);
    Mass_g2fraction.setConstant(false);
    Mass_pdf.fitTo(*rds_hist);
    //Mass_pdf.fitTo(*rds); // bin fit is already good enough
    
    RooPlot* frame_m = Mass->frame();

    histo_data->Sumw2(false);
    histo_data->SetBinErrorOption(TH1::kPoisson);
    histo_data->SetMarkerStyle(20);
    histo_data->SetMarkerSize(0.4);
    histo_data->SetLineColor(kBlack);

    rds->plotOn(frame_m,Name("data"),Binning(320),Invisible());
    Mass_pdf.plotOn(frame_m,Name("model"),Precision(2E-4));
    
    Mass_pdf.plotOn(frame_m,Name("sig"),Precision(2E-4),Components("Mass_pdf_sig"),LineColor(kRed-10),LineWidth(2),LineStyle(kSolid),FillColor(kRed-10), VLines(), DrawOption("F"));
    
    Mass_pdf.plotOn(frame_m,Name("comb"),Precision(2E-4),Components("Mass_pdf_comb"),LineColor(kCyan+2),LineWidth(3),LineStyle(7));
    
    gStyle->SetHatchesSpacing(2.);
    gStyle->SetHatchesLineWidth(2);
    
    Mass_pdf.plotOn(frame_m,Name("jpsipi"),Precision(2E-4),Components("Mass_pdf_jpsipi"),LineColor(kViolet),LineWidth(1),LineStyle(kSolid),FillStyle(3145),FillColor(kViolet), VLines(), DrawOption("F"));
    Mass_pdf.plotOn(frame_m,Name("jpsipi_L"),Precision(2E-4),Components("Mass_pdf_jpsipi"),LineColor(kViolet),LineWidth(1),LineStyle(kSolid));
    
    Mass_pdf.plotOn(frame_m,Name("nonprompt"),Precision(2E-4),Components("Mass_pdf_nonprompt"),LineColor(kOrange-6),LineWidth(3),LineStyle(3));
    
    TCanvas* canvas_m = new TCanvas("canvas_m", "", 600, 600);
    canvas_m->SetMargin(0.14,0.06,0.13,0.07);
    
    frame_m->GetYaxis()->SetTitleOffset(1.48);
    frame_m->GetYaxis()->SetTitle("Entries");
    frame_m->GetXaxis()->SetTitleOffset(1.15);
    frame_m->GetXaxis()->SetLabelOffset(0.01);
    frame_m->GetXaxis()->SetTitle("Mass [GeV]");
    frame_m->GetXaxis()->SetTitleSize(0.043);
    frame_m->GetYaxis()->SetTitleSize(0.043);
    frame_m->SetStats(false);
    frame_m->SetTitle("");
    
    frame_m->Draw();
    histo_data->Draw("Esame");
    
    TLegend *leg1 = new TLegend(0.55,0.64,0.88,0.91);
    leg1->SetTextSize(0.04);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(kWhite);
    leg1->AddEntry(histo_data,"Data", "EP");
    leg1->AddEntry("model","Fit to all entries", "L");
    leg1->AddEntry("sig","B^{+} #rightarrow J/#psi K^{+} signal", "F");
    leg1->AddEntry("comb","Combinatorial", "L");
    leg1->AddEntry("jpsipi","B^{+} #rightarrow J/#psi #pi^{+}", "F");
    leg1->AddEntry("nonprompt","B #rightarrow J/#psi + hadrons", "L");
    leg1->Draw();
    
    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.035);
    tex.SetTextAlign(11);
    tex.SetNDC();
    tex.DrawLatex(0.14,0.94,"CMS Preliminary");
    
    canvas_m->Print(Form("fig/proj_mass_%s_data_normal.pdf",tag.Data()));
    
    canvas_m->SetLogy();
    canvas_m->Update();
    
    canvas_m->Print(Form("fig/proj_mass_%s_data_log.pdf",tag.Data()));
    
    // Now –– produce the sPlots
    
    enum {_sig, _comb, _jpsipi, _nonprompt, _nspec};
    
    double yield[_nspec];
    yield[_sig] = n_sig.getVal();
    yield[_comb] = n_comb.getVal();
    yield[_jpsipi] = n_jpsipi.getVal();
    yield[_nonprompt] = n_nonprompt.getVal();
    
    TMatrixD covInv(_nspec, _nspec);
    covInv = 0.;
    RooArgSet norm(*Mass);
    
    for (int evt=0; evt<rds->numEntries(); evt++) {
        const RooArgSet* arg = rds->get(evt);
        
        Mass->setVal(arg->getRealValue("Mass"));
        Tau->setVal(arg->getRealValue("Tau"));
        TauRes->setVal(arg->getRealValue("TauRes"));
        
        double pdf[_nspec];
        pdf[_sig] = Mass_pdf_sig.getVal(&norm);
        pdf[_comb] = Mass_pdf_comb.getVal(&norm);
        pdf[_jpsipi] = Mass_pdf_jpsipi.getVal(&norm);
        pdf[_nonprompt] = Mass_pdf_nonprompt.getVal(&norm);
        
        double pdf_total = 0.;
        for (int idx = 0; idx<_nspec; idx++) pdf_total += yield[idx]*pdf[idx];
        
        for (int row = 0; row<_nspec; row++)
            for (int col = 0; col<_nspec; col++)
                covInv(row,col) += pdf[row]*pdf[col]/(pdf_total*pdf_total);
    }
    
    TMatrixD covMatrix(TMatrixD::kInverted,covInv);
    
    double mean_taue = 0.;
    double count_taue = 0;
    for (int evt=0; evt<rds->numEntries(); evt++) {
        const RooArgSet* arg = rds->get(evt);
        
        Mass->setVal(arg->getRealValue("Mass"));
        Tau->setVal(arg->getRealValue("Tau"));
        TauRes->setVal(arg->getRealValue("TauRes"));
        
        double pdf[_nspec];
        pdf[_sig] = Mass_pdf_sig.getVal(&norm);
        pdf[_comb] = Mass_pdf_comb.getVal(&norm);
        pdf[_jpsipi] = Mass_pdf_jpsipi.getVal(&norm);
        pdf[_nonprompt] = Mass_pdf_nonprompt.getVal(&norm);
        
        double denominator = 0.;
        for (int idx = 0; idx<_nspec; idx++) denominator += yield[idx]*pdf[idx];
        
        double numerator = 0.;
        for (int idx = 0; idx<_nspec; idx++) numerator += covMatrix(_sig,idx)*pdf[idx];
        
        double weight = numerator/denominator*arg->getRealValue("Weight");
        h_tau->Fill(arg->getRealValue("Tau"),weight);
        h_taures->Fill(arg->getRealValue("TauRes"),weight);
        h_taureco->Fill(arg->getRealValue("Tau"),weight);
        
        mean_taue += arg->getRealValue("TauRes")*weight;
        count_taue += weight;
    }
    mean_taue /= count_taue;
    cout << ">>> Average tau_err = " << mean_taue << endl;
  
    // ISSUE: Consider to set the <taue> to the data value?
    //wspace->var("respar_bupsik_taue")->setVal(mean_taue);
    
    // Set back the usual boundary for Tau
    Tau->setMin(Tau_bound[0]);

    // Binned weighted likelihood fit w/ PDF bin integration
    RooFitResult *res1 = NULL, *res2 = NULL;
    Fit_sPlot(h_tau, &Tau_pdf_reco, Tau, &EffTau, &res1, &res2);
    
    // Projection
    RooPlot* frame_t = Tau->frame(Title(" "));
    h_tau->SetMarkerStyle(1);
    h_tau->SetLineColor(kBlack);
    RooDataHist *h_tau_data = new RooDataHist("h_tau_data", "", RooArgList(*Tau), h_tau);
    h_tau_data->plotOn(frame_t,MarkerStyle(1),LineColor(kBlack));
    
    Tau_pdf_reco.plotOn(frame_t, DrawOption("L"), LineColor(kBlue), LineWidth(2), LineStyle(1), NumCPU(NCPU));
    
    TCanvas* canvas_t = new TCanvas("canvas_t", "", 600, 600);
    canvas_t->SetMargin(0.14,0.06,0.13,0.07);
    
    frame_t->GetYaxis()->SetTitleOffset(1.48);
    frame_t->GetYaxis()->SetTitle("Entries");
    frame_t->GetXaxis()->SetTitleOffset(1.15);
    frame_t->GetXaxis()->SetLabelOffset(0.01);
    frame_t->GetXaxis()->SetTitle("#tau [ps]");
    frame_t->GetXaxis()->SetTitleSize(0.043);
    frame_t->GetYaxis()->SetTitleSize(0.043);
    frame_t->SetStats(false);
    frame_t->Draw("E");
    
    tex.DrawLatex(0.14,0.94,"CMS Preliminary");
    
    TLine lin;
    lin.SetLineColor(kGray+1);
    lin.SetLineWidth(2);
    lin.SetLineStyle(7);
    lin.DrawLine(1.,0.,12.,0.);
    
    canvas_t->Update();
    
    TLegend *leg2 = new TLegend(0.50,0.86,0.91,0.91);
    leg2->SetNColumns(1);
    leg2->SetFillColor(kWhite);
    leg2->SetLineColor(kWhite);
    leg2->AddEntry(h_tau, "B^{+} #rightarrow J/#psi K^{+} Data", "lep");
    leg2->Draw();
    
    canvas_t->Print(Form("fig/tau_splot_%s_data_reco.pdf",tag.Data()));
    
    h_taures->GetYaxis()->SetTitleOffset(1.48);
    h_taures->GetYaxis()->SetTitle("Entries");
    h_taures->GetXaxis()->SetTitleOffset(1.15);
    h_taures->GetXaxis()->SetLabelOffset(0.01);
    h_taures->GetXaxis()->SetTitle("#sigma_{#tau} [ps]");
    h_taures->GetXaxis()->SetTitleSize(0.043);
    h_taures->GetYaxis()->SetTitleSize(0.043);
    h_taures->SetStats(false);
    h_taures->SetMarkerStyle(20);
    h_taures->SetLineColor(kBlack);
    h_taures->Draw("E");
    
    tex.DrawLatex(0.14,0.94,"CMS Preliminary");
    
    canvas_t->Print(Form("fig/taures_splot_%s_data_reco.pdf",tag.Data()));
    
    if (!opt.Contains("taucorr")) wspace->import(*h_taureco);
    delete h_taureco;
    
    delete h_tau;
    delete h_taures;
    delete rds;
    delete canvas_m;
    delete canvas_t;
    delete leg1;
    delete leg2;
    delete frame_m;
    delete frame_t;
    delete h_tau_data;
}

// ----------------------------------------------------
// available options:
// bupsik    - produce the corrections based on B->J/psi K+ data and MC
//
void PrepareLifetimeEfficiencyCorrection(RooWorkspace *wspace, TString opt = "bupsik", TString target_cat = "2016BF")
{
    cout << ">>> PrepareLifetimeEfficiencyCorrection() start" << endl;
    
    TString tag;
    if (opt.Contains("bupsik")) tag = "bupsik";
    else {
        cout << ">>> Undefined option: " << opt << endl;
        return;
    }
    tag += Form("_%s", target_cat.Data());
    
    TString target_era = target_cat;
    for (auto& cat: CatMan.cats)
        if (cat.id == target_cat) target_era = cat.era;
    
    // extract h_taureco and h_taureco_mc
    TH1D *h_taureco = (TH1D*)wspace->obj(Form("h_taureco_%s",tag.Data()));
    TH1D *h_taureco_mc = (TH1D*)wspace->obj(Form("h_taureco_mc_%s",tag.Data()));
    
    // extract lifetime & resolution
    RooRealVar *Tau = wspace->var("Tau");
    RooResolutionModel *TauRes_Model = (RooResolutionModel*)wspace->obj(Form("TauRes_Model_%s",tag.Data()));
    //RooTruthModel *TauRes_Model = new RooTruthModel("TauRes_Model","",*Tau); // ideal model
    RooRealVar EffTau("EffTau","",1.6);
    RooDecay RawDecay("RawDecay","",*Tau,EffTau,*TauRes_Model,RooDecay::SingleSided);
    
    TH1D *h_norm = (TH1D*)h_taureco->Clone("h_norm");
    
    for (TH1D* h_reco: {h_taureco,h_taureco_mc}) {
        if (h_reco==h_taureco) EffTau.setVal(1.638); // PDG = 1.638 +/- 0.004 ps
        else {
            if (target_era.Contains("2011") || target_era.Contains("2012")) EffTau.setVal(1.671);
            if (target_era.Contains("2016")) EffTau.setVal(1.638);
        }
        
        h_norm->Reset();
        double scale_factor = 0.;
        for (int bin=1; bin<=h_norm->GetNbinsX(); bin++) {
            
            Tau->setRange("bin",h_norm->GetBinLowEdge(bin),h_norm->GetBinLowEdge(bin)+h_norm->GetBinWidth(bin));
            RooAbsReal* area = RawDecay.createIntegral(*Tau,Range("bin"));
            h_norm->SetBinContent(bin,area->getVal());
            delete area;
            
            double ratio = h_reco->GetBinContent(bin)/h_norm->GetBinContent(bin);
            if (ratio>scale_factor) scale_factor = ratio;
        }
        h_norm->Scale(scale_factor);
     
        for (int bin=1; bin<=h_reco->GetNbinsX(); bin++) {
            double val = h_reco->GetBinContent(bin)/h_norm->GetBinContent(bin);
            double err = h_reco->GetBinError(bin)/h_norm->GetBinContent(bin);
            h_reco->SetBinContent(bin,val);
            h_reco->SetBinError(bin,err);
        }
    }
    
    TH1D *h_taucorr = (TH1D*)h_taureco->Clone(Form("h_taucorr_%s",tag.Data()));
    h_taucorr->Divide(h_taureco_mc);
    wspace->import(*h_taucorr);
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.14,0.06,0.13,0.07);
    canvas->SetGrid();
    
    for (TH1D* h_reco: {h_taureco,h_taureco_mc}) {
        h_reco->GetYaxis()->SetTitleOffset(1.48);
        h_reco->GetYaxis()->SetTitle("Efficiency Function [Arbitary Unit]");
        h_reco->GetXaxis()->SetTitleOffset(1.15);
        h_reco->GetXaxis()->SetLabelOffset(0.01);
        h_reco->GetXaxis()->SetTitle("#tau [ps]");
        h_reco->GetXaxis()->SetTitleSize(0.043);
        h_reco->GetYaxis()->SetTitleSize(0.043);
        h_reco->SetStats(0);
        h_reco->SetMaximum(1.3);
        h_reco->SetMinimum(0.0);
        h_reco->SetMarkerSize(0.5);
    }
    
    h_taureco_mc->SetMarkerStyle(20);
    h_taureco_mc->SetLineColor(kBlack);
    h_taureco_mc->SetMarkerColor(kBlack);
    h_taureco_mc->Draw();
    
    h_taureco->SetMarkerStyle(21);
    h_taureco->SetLineColor(kRed);
    h_taureco->SetMarkerColor(kRed);
    h_taureco->Draw("same");
    
    TLatex tex;
    tex.SetTextSize(0.05);
    tex.DrawLatex(1.5,1.2,target_cat);
    
    TLegend *leg1 = new TLegend(0.5,0.80,0.93,0.90);
    leg1->SetFillColor(kWhite);
    leg1->SetFillStyle(0);
    leg1->SetLineColor(kWhite);
    leg1->SetLineWidth(0);
    leg1->AddEntry(h_taureco,"J/#psi K^{+} data","PL");
    leg1->AddEntry(h_taureco_mc,"J/#psi K^{+} MC","PL");
    leg1->Draw();
    
    canvas->Print(Form("fig/taueff_comp_%s.pdf",tag.Data()));
    
    h_taucorr->GetYaxis()->SetTitleOffset(1.48);
    h_taucorr->GetYaxis()->SetTitle("Ratio of Efficiency Functions [Arbitary Unit]");
    h_taucorr->GetXaxis()->SetTitleOffset(1.15);
    h_taucorr->GetXaxis()->SetLabelOffset(0.01);
    h_taucorr->GetXaxis()->SetTitle("#tau [ps]");
    h_taucorr->GetXaxis()->SetTitleSize(0.043);
    h_taucorr->GetYaxis()->SetTitleSize(0.043);
    h_taucorr->SetStats(0);
    h_taucorr->SetMaximum(1.5);
    h_taucorr->SetMinimum(0.5);
    h_taucorr->SetMarkerStyle(20);
    h_taucorr->SetMarkerSize(0.5);
    h_taucorr->SetLineColor(kBlack);
    h_taucorr->SetMarkerColor(kBlack);
    h_taucorr->Draw();
    
    tex.SetTextSize(0.05);
    tex.DrawLatex(1.5,1.4,target_cat);
    
    canvas->Print(Form("fig/taucorr_%s.pdf",tag.Data()));
    
    delete leg1;
    delete canvas;
    delete h_norm;
    delete h_taucorr;
}

void bmm4validation(TString commands = "")
{
    // -----------------------------------------------------------
    // parse the commands
    bool do_tausplot_mcstudy = false, do_tausplot_datastudy = false, do_tausplot_correction = false;
    double alter_BDT = -2.;
    TString select_era = "2016BF";
    TString select_proc = "bsmm";
    TString select_type = "reco";
    TString file_valid = "wspace_valid.root";
    
    cout << ">>> -------------------------" << endl;
    cout << ">>> BMM4 validation start" << endl;
    cout << ">>> -------------------------" << endl;
    cout << ">>> commands:" << endl;
    cout << ">>> - tausplot_mcstudy                   : perform MC study for lifetime fit" << endl;
    cout << ">>> - tausplot_datastudy                 : perform data study for lifetime fit (J/psi K+ only)" << endl;
    cout << ">>> - tausplot_correction                : compute correction factors for each year/era" << endl;
    cout << ">>> - era=[2016BF]                       : set the target sample (from which year/era)" << endl;
    cout << ">>> - proc=[bsmm]                        : set the target process, Bs->mumu (bsmm) or B+->J/psiK+ (bupsik)" << endl;
    cout << ">>> - bdt=[-2.]                          : alternative BDT threshold, -2: using the category definition." << endl;
    cout << ">>> - type=[reco]                        : set the target type " << endl;
    cout << ">>>         genonly   - fit to pure generated decay time" << endl;
    cout << ">>>         geneff    - fit to generated decay time, with efficiency correction" << endl;
    cout << ">>>         reco      - fit to reco decay time" << endl;
    cout << ">>>         toybkg    - mixing with toy combinatorial background" << endl;
    cout << ">>> - ws_valid=[wspace_valid.root]       : set output workspace" << endl;
    cout << ">>> parsing commands: [" << commands << "]" << endl;
    if (commands=="") return;
    
    TString tok;
    Ssiz_t from = 0;
    while(commands.Tokenize(tok, from, "[ \t;=:]")) {
             if (tok=="tausplot_mcstudy") do_tausplot_mcstudy = true;       // perform MC study for lifetime fit
        else if (tok=="tausplot_datastudy") do_tausplot_datastudy = true;   // perform data study for lifetime fit, J/psi K+ only
        else if (tok=="tausplot_correction") do_tausplot_correction = true; // compute correction factors for each year/era
        else if (tok=="era") {                                              // set the target era
            commands.Tokenize(tok, from, "[ \t;=:]");
            select_era = tok;
        }else if (tok=="bdt") {                                             // alternative BDT threshold
            commands.Tokenize(tok, from, "[ \t;=:]");
            alter_BDT = tok.Atof();
        }else if (tok=="ws_valid") {                                        // set source model workspace
            commands.Tokenize(tok, from, "[ \t;=:]");
            file_valid = tok;
        }else if (tok=="proc") {                                            // set the target process
            commands.Tokenize(tok, from, "[ \t;=:]");
            select_proc = tok;
        }else if (tok=="type") {                                            // set the target type
            commands.Tokenize(tok, from, "[ \t;=:]");
            select_type = tok;
        }else {
            cout << ">>> unknown command '" << tok << "'" << endl;
            return;
        }
    }
    
    if (do_tausplot_datastudy && select_proc=="bsmm") {
        cout << ">>> tausplot_datastudy only works for B+->J/psiK+ (bupsik)." << endl;
        cout << ">>> proc=bupsik set." << endl;
        select_proc = "bupsik";
    }
    if (do_tausplot_correction && select_proc=="bsmm") {
        cout << ">>> tausplot_correction only works for B+->J/psiK+ (bupsik)." << endl;
        cout << ">>> proc=bupsik set." << endl;
        select_proc = "bupsik";
    }
    
    // -----------------------------------------------------------
    // Initializing categories definition
    
    if (CONFIG_BMM3) CatMan.RegisterBMM3Categories();
    if (CONFIG_BMM4) CatMan.RegisterBMM4Categories();
    CatMan.Print();
    
    if (do_tausplot_mcstudy || do_tausplot_datastudy || do_tausplot_correction) {
    
        // -----------------------------------------------------------
        // Create the shared workspace
        RooWorkspace *wspace = new RooWorkspace("wspace");
    
        // Observables
        RooRealVar Mass("Mass", "", Mass_bound[0], Mass_bound[1]);
        RooRealVar Tau("Tau", "", Tau_bound[0], Tau_bound[1]);
        RooRealVar TauRes("TauRes", "", TauRes_bound[0], TauRes_bound[1]); // in unit of ps
        RooRealVar Weight("Weight", "", 0., 1E10); // event weighting for internal use
        
        wspace->import(Mass);
        wspace->import(Tau);
        wspace->import(TauRes);
        wspace->import(Weight);
    
        if (do_tausplot_mcstudy) {
            // rebuild the efficiency & resolution models
            PrepareLifetimeResolutionModel(wspace, select_proc+":triple", select_era, alter_BDT);
            PrepareLifetimeEfficiencyModel(wspace, select_proc+":threshold", select_era, alter_BDT);
            
            // run MC study
            PerformMCTauSPlotStudy(wspace,select_proc+":"+select_type, select_era, alter_BDT);
        }
    
        if (do_tausplot_datastudy) {
            // rebuild the efficiency & resolution models
            PrepareLifetimeResolutionModel(wspace, select_proc+":triple", select_era, alter_BDT);
            PrepareLifetimeEfficiencyModel(wspace, select_proc+":threshold", select_era, alter_BDT);

            // run data study
            PerformBuDataTauSPlotStudy(wspace, select_proc, select_era, alter_BDT);
            
            // rebuild the efficiency model w/ correction
            //PrepareLifetimeEfficiencyCorrection(wspace, select_proc, select_era);
            //PrepareLifetimeEfficiencyModel(wspace, select_proc+":threshold"+":taucorr", select_era, alter_BDT);
            
            // run data study again
            //PerformBuDataTauSPlotStudy(wspace, select_proc+":taucorr", select_era, alter_BDT);
        }
        
        if (do_tausplot_correction) {
            //for (TString era : bmm4::eras) {
            for (TString era : {"2016BF","2016GH"}) {
                // rebuild the efficiency & resolution models
                PrepareLifetimeResolutionModel(wspace, select_proc+":triple", era, alter_BDT);
                PrepareLifetimeEfficiencyModel(wspace, select_proc+":threshold", era, alter_BDT);
                
                // run data study
                PerformBuDataTauSPlotStudy(wspace, select_proc, era, alter_BDT);
                
                // calculate the correction histograms
                PrepareLifetimeEfficiencyCorrection(wspace, select_proc, era);
            }
        }
        
        wspace->writeToFile(file_valid);
    }
    cout << ">>> BMM4 validation end." << endl;
}

int main(int argc, char *argv[]) {
    TString cmd = "";
    for (int i=1; i<argc; i++) {
        cmd += TString(argv[i])+";";
    }
    bmm4validation(cmd);
}
