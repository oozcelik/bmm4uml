#include "bmm4common.h"

void ProduceToyTauSPlot(RooWorkspace *wspace_res, RooWorkspace *wspace, RooDataSet *toy_data, int toy_idx, bool display = false)
{
    cout << ">>> ProduceToyTauSPlot() start" << endl;
    
    enum {_bs, _bd, _peak, _semi, _h2mu, _comb, _nspec};
    vector<TString> specs= {"bs", "bd", "peak", "semi", "h2mu", "comb"};
    
    TH1D *h_tau = new TH1D("h_tau_bs", "", Tau_bins.size()-1, Tau_bins.data());
    h_tau->Sumw2();
    
    for (auto& cat: CatMan.cats) {
        
        RooRealVar *Mass = wspace->var("Mass");
        RooRealVar *ReducedMassRes = wspace->var("ReducedMassRes");
        RooCategory *PairCat = wspace->cat("PairCat");
        RooArgSet norm(*Mass,*ReducedMassRes,*PairCat);
        
        double yield[_nspec];
        yield[_bs]   = wspace->function(Form("N_bs_formula_%s",cat.id.Data()))->getVal();
        yield[_bd]   = wspace->function(Form("N_bd_formula_%s",cat.id.Data()))->getVal();
        yield[_peak] = wspace->function(Form("N_peak_formula_%s",cat.id.Data()))->getVal();
        yield[_semi] = wspace->function(Form("N_semi_formula_%s",cat.id.Data()))->getVal();
        yield[_h2mu] = wspace->var(Form("N_h2mu_%s",cat.id.Data()))->getVal();
        yield[_comb] = wspace->var(Form("N_comb_%s",cat.id.Data()))->getVal();
        
        TMatrixD covInv(_nspec, _nspec);
        covInv = 0.;
        
        for (int evt=0; evt<toy_data->numEntries(); evt++) {
            const RooArgSet* arg = toy_data->get(evt);
            if (arg->getCatIndex("GlobalCat")!=cat.index) continue;
            
            Mass->setVal(arg->getRealValue("Mass"));
            ReducedMassRes->setVal(arg->getRealValue("ReducedMassRes"));
            PairCat->setIndex(arg->getCatIndex("PairCat"));
            
            double pdf[_nspec];
            for (int idx = 0; idx<_nspec; idx++)
                pdf[idx] = wspace->pdf(Form("pdf_%s_%s",specs[idx].Data(),cat.id.Data()))->getVal(&norm);
            
            double pdf_total = 0.;
            for (int idx = 0; idx<_nspec; idx++) pdf_total += yield[idx]*pdf[idx];
            
            for (int row = 0; row<_nspec; row++)
                for (int col = 0; col<_nspec; col++)
                    covInv(row,col) += pdf[row]*pdf[col]/(pdf_total*pdf_total);
        }
        
        TMatrixD covMatrix(TMatrixD::kInverted,covInv);
        
        for (int evt=0; evt<toy_data->numEntries(); evt++) {
            const RooArgSet* arg = toy_data->get(evt);
            if (arg->getCatIndex("GlobalCat")!=cat.index) continue;
            if (arg->getCatIndex("SelCat")!=1) continue;
            
            Mass->setVal(arg->getRealValue("Mass"));
            ReducedMassRes->setVal(arg->getRealValue("ReducedMassRes"));
            PairCat->setIndex(arg->getCatIndex("PairCat"));
            
            double pdf[_nspec];
            for (int idx = 0; idx<_nspec; idx++)
                pdf[idx] = wspace->pdf(Form("pdf_%s_%s",specs[idx].Data(),cat.id.Data()))->getVal(&norm);
            
            double denominator = 0.;
            for (int idx = 0; idx<_nspec; idx++) denominator += yield[idx]*pdf[idx];

            double numerator = 0.;
            for (int idx = 0; idx<_nspec; idx++) numerator += covMatrix(_bs,idx)*pdf[idx];
            
            double weight = numerator/denominator;
            
            h_tau->Fill(arg->getRealValue("Tau"),weight);
        }
    }
    
    RooRealVar *Tau = wspace->var("Tau");
    RooRealVar *EffTau_bs = wspace->var("EffTau_bs");
    RooAbsPdf *Tau_pdf_bs = wspace->pdf("Tau_pdf_bs_mix");
    
    EffTau_bs->setConstant(false);
    EffTau_bs->setMax(8.0);
    
    // Binned weighted likelihood fit w/ PDF bin integration
    RooFitResult *res1 = NULL, *res2 = NULL;
    Fit_sPlot(h_tau, Tau_pdf_bs, Tau, EffTau_bs, &res1, &res2);
    
    wspace_res->import(*res1,Form("fitresult_taubs_toy_%d",toy_idx)); // results with corrected uncertainties
    wspace_res->import(*res2,Form("fitresult_taubs_wl_toy_%d",toy_idx)); // results from the 2nd weighted likelihood fit
    
    if (display) {
        RooDataHist *h_tau_data = new RooDataHist("h_tau_data", "", RooArgList(*Tau), h_tau);
        
        TString title = "CMS Preliminary";
        
        RooPlot* frame = Tau->frame(Title(" "));
        h_tau->SetMarkerStyle(20);
        h_tau->SetLineColor(kBlack);
        h_tau_data->plotOn(frame,MarkerStyle(20),LineColor(kBlack));
        
        Tau_pdf_bs->plotOn(frame, DrawOption("L"), LineColor(kBlue), LineWidth(2), LineStyle(1), NumCPU(NCPU));
        
        TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetMargin(0.14,0.06,0.13,0.07);
        
        frame->GetYaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitle("Entries");
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetXaxis()->SetLabelOffset(0.01);
        frame->GetXaxis()->SetTitle("#tau [ps]");
        frame->GetXaxis()->SetTitleSize(0.043);
        frame->GetYaxis()->SetTitleSize(0.043);
        frame->SetStats(false);
        frame->Draw("E");
        gStyle->SetErrorX(0);
        
        TLatex tex;
        tex.SetTextFont(42);
        tex.SetTextSize(0.035);
        tex.SetTextAlign(11);
        tex.SetNDC();
        tex.DrawLatex(0.14,0.94,title);
        
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
        leg1->AddEntry(h_tau, "weighted toy (B_{s})", "lep");
        leg1->Draw();
        
        canvas->Print("fig/tau_splot_fit_toy_bs.pdf");
        
        delete leg1;
        delete h_tau_data;
        delete canvas;
    }
    
    delete h_tau;
    delete res1;
    delete res2;
}

// Available options:
// ---------
// prefit   - adopt pre-fit nuisances, randomize all of the nuisances according to the constraints
// freq     - generate frequentist toy, adopt post-fit nuisances, randomize all of the constrained means
// ---------
// signif   - do significance estimation
// bdstat   - Calculate profile Likelihood test statistics (F&C study or upper limit calculation for Bd)
// ---------
// tausplot - add decay time component, produce splot, and fit
//
void PerformToyStudy(RooWorkspace *wspace_res, RooWorkspace* wspace_gen, RooWorkspace* wspace_fit, int iterations = 10, bool do_minos = false, TString opt = "prefit;signif")
{
    cout << ">>> PerformToyStudy() start" << endl;
    
    RooRealVar *Mass = wspace_gen->var("Mass");
    RooRealVar *ReducedMassRes = wspace_gen->var("ReducedMassRes");
    RooRealVar *Tau = wspace_gen->var("Tau");
    RooRealVar *TauRes = wspace_gen->var("TauRes");
    RooCategory *SelCat = wspace_gen->cat("SelCat");
    RooCategory *PairCat = wspace_gen->cat("PairCat");
    RooCategory *GlobalCat = wspace_gen->cat("GlobalCat");
    RooArgSet varlist(*Mass,*ReducedMassRes,*PairCat,*Tau,*TauRes,*SelCat,*GlobalCat);
    
    for(int idx=0; idx<iterations; idx++) {
        
        cout << ">>> # iterations: " << idx+1 << "/" << iterations << endl;
        
        RooWorkspace* wspace = new RooWorkspace(*wspace_fit); // clone wspace_fit
        
        if (opt.Contains("prefit")) { // adopt pre-fit nuisances, randomize all of the nuisances according to the constraints
            RooDataSet *tmp_evt;
            
            // change the wspace_gen
            RooRealVar *fs_over_fu = wspace_gen->var("fs_over_fu");
            RooRealVar *one_over_BRBR = wspace_gen->var("one_over_BRBR");
            RooAbsPdf *fs_over_fu_gau = wspace_gen->pdf("fs_over_fu_gau");
            RooAbsPdf *one_over_BRBR_gau = wspace_gen->pdf("one_over_BRBR_gau");
            
            tmp_evt = fs_over_fu_gau->generate(RooArgSet(*fs_over_fu),1);
            fs_over_fu->setVal(tmp_evt->get(0)->getRealValue(fs_over_fu->GetName())); delete tmp_evt;
            
            tmp_evt = one_over_BRBR_gau->generate(RooArgSet(*one_over_BRBR),1);
            one_over_BRBR->setVal(tmp_evt->get(0)->getRealValue(one_over_BRBR->GetName())); delete tmp_evt;
            
            for (auto& cat: CatMan.cats) {
                
                // change the wspace_gen
                RooRealVar *N_bu = wspace_gen->var(Form("N_bu_%s", cat.id.Data()));
                RooRealVar *N_peak = wspace_gen->var(Form("N_peak_%s", cat.id.Data()));
                RooRealVar *N_semi = wspace_gen->var(Form("N_semi_%s", cat.id.Data()));
                RooRealVar *N_h2mu = wspace_gen->var(Form("N_h2mu_%s", cat.id.Data()));
                RooRealVar *effratio_bs = wspace_gen->var(Form("effratio_bs_%s", cat.id.Data()));
                RooRealVar *effratio_bd = wspace_gen->var(Form("effratio_bd_%s", cat.id.Data()));
                
                RooAbsPdf *N_bu_gau = wspace_gen->pdf(Form("N_bu_%s_gau", cat.id.Data()));
                RooAbsPdf *N_peak_lnn = wspace_gen->pdf(Form("N_peak_%s_lnn", cat.id.Data()));
                RooAbsPdf *N_semi_gau = wspace_gen->pdf(Form("N_semi_%s_gau", cat.id.Data()));
                RooAbsPdf *N_h2mu_lnn = wspace_gen->pdf(Form("N_h2mu_%s_lnn", cat.id.Data()));
                RooAbsPdf *effratio_bs_gau = wspace_gen->pdf(Form("effratio_bs_%s_gau", cat.id.Data()));
                RooAbsPdf *effratio_bd_gau = wspace_gen->pdf(Form("effratio_bd_%s_gau", cat.id.Data()));
                
                tmp_evt = N_bu_gau->generate(RooArgSet(*N_bu),1);
                N_bu->setVal(tmp_evt->get(0)->getRealValue(N_bu->GetName())); delete tmp_evt;
                
                tmp_evt = N_peak_lnn->generate(RooArgSet(*N_peak),1);
                N_peak->setVal(tmp_evt->get(0)->getRealValue(N_peak->GetName())); delete tmp_evt;
                
                tmp_evt = N_semi_gau->generate(RooArgSet(*N_semi),1);
                N_semi->setVal(tmp_evt->get(0)->getRealValue(N_semi->GetName())); delete tmp_evt;
                
                tmp_evt = N_h2mu_lnn->generate(RooArgSet(*N_h2mu),1);
                N_h2mu->setVal(tmp_evt->get(0)->getRealValue(N_h2mu->GetName())); delete tmp_evt;
                
                tmp_evt = effratio_bs_gau->generate(RooArgSet(*effratio_bs),1);
                effratio_bs->setVal(tmp_evt->get(0)->getRealValue(effratio_bs->GetName())); delete tmp_evt;
                
                tmp_evt = effratio_bd_gau->generate(RooArgSet(*effratio_bd),1);
                effratio_bd->setVal(tmp_evt->get(0)->getRealValue(effratio_bd->GetName())); delete tmp_evt;
            }
        }else if (opt.Contains("freq")) { // generate frequentist toy, adopt post-fit nuisances, randomize all of the constrained means
            RooDataSet *tmp_evt;
            
            // change the clone of wspace_fit
            RooRealVar *fs_over_fu = wspace->var("fs_over_fu");
            RooRealVar *one_over_BRBR = wspace->var("one_over_BRBR");
            RooRealVar *fs_over_fu_mean = wspace->var("fs_over_fu_mean");
            RooRealVar *one_over_BRBR_mean = wspace->var("one_over_BRBR_mean");
            RooAbsPdf *fs_over_fu_gau = wspace->pdf("fs_over_fu_gau");
            RooAbsPdf *one_over_BRBR_gau = wspace->pdf("one_over_BRBR_gau");
            
            fs_over_fu_mean->setVal(fs_over_fu->getVal());
            tmp_evt = fs_over_fu_gau->generate(RooArgSet(*fs_over_fu),1);
            fs_over_fu_mean->setVal(tmp_evt->get(0)->getRealValue(fs_over_fu->GetName())); delete tmp_evt;

            one_over_BRBR_mean->setVal(one_over_BRBR->getVal());
            tmp_evt = one_over_BRBR_gau->generate(RooArgSet(*one_over_BRBR),1);
            one_over_BRBR_mean->setVal(tmp_evt->get(0)->getRealValue(one_over_BRBR->GetName())); delete tmp_evt;
            
            for (auto& cat: CatMan.cats) {
                
                // change the clone of wspace_fit
                RooRealVar *N_bu = wspace->var(Form("N_bu_%s", cat.id.Data()));
                RooRealVar *N_peak = wspace->var(Form("N_peak_%s", cat.id.Data()));
                RooRealVar *N_semi = wspace->var(Form("N_semi_%s", cat.id.Data()));
                RooRealVar *N_h2mu = wspace->var(Form("N_h2mu_%s", cat.id.Data()));
                RooRealVar *effratio_bs = wspace->var(Form("effratio_bs_%s", cat.id.Data()));
                RooRealVar *effratio_bd = wspace->var(Form("effratio_bd_%s", cat.id.Data()));

                RooRealVar *N_bu_mean = wspace->var(Form("N_bu_%s_mean", cat.id.Data()));
                RooRealVar *N_peak_mean = wspace->var(Form("N_peak_%s_mean", cat.id.Data()));
                RooRealVar *N_semi_mean = wspace->var(Form("N_semi_%s_mean", cat.id.Data()));
                RooRealVar *N_h2mu_mean = wspace->var(Form("N_h2mu_%s_mean", cat.id.Data()));
                RooRealVar *effratio_bs_mean = wspace->var(Form("effratio_bs_%s_mean", cat.id.Data()));
                RooRealVar *effratio_bd_mean = wspace->var(Form("effratio_bd_%s_mean", cat.id.Data()));
                
                RooAbsPdf *N_bu_gau = wspace->pdf(Form("N_bu_%s_gau", cat.id.Data()));
                RooAbsPdf *N_peak_lnn = wspace->pdf(Form("N_peak_%s_lnn", cat.id.Data()));
                RooAbsPdf *N_semi_gau = wspace->pdf(Form("N_semi_%s_gau", cat.id.Data()));
                RooAbsPdf *N_h2mu_lnn = wspace->pdf(Form("N_h2mu_%s_lnn", cat.id.Data()));
                RooAbsPdf *effratio_bs_gau = wspace->pdf(Form("effratio_bs_%s_gau", cat.id.Data()));
                RooAbsPdf *effratio_bd_gau = wspace->pdf(Form("effratio_bd_%s_gau", cat.id.Data()));
                
                N_bu_mean->setVal(N_bu->getVal());
                tmp_evt = N_bu_gau->generate(RooArgSet(*N_bu),1);
                N_bu_mean->setVal(tmp_evt->get(0)->getRealValue(N_bu->GetName())); delete tmp_evt;
                
                N_peak_mean->setVal(N_peak->getVal());
                tmp_evt = N_peak_lnn->generate(RooArgSet(*N_peak),1);
                N_peak_mean->setVal(tmp_evt->get(0)->getRealValue(N_peak->GetName())); delete tmp_evt;
                
                N_semi_mean->setVal(N_semi->getVal());
                tmp_evt = N_semi_gau->generate(RooArgSet(*N_semi),1);
                N_semi_mean->setVal(tmp_evt->get(0)->getRealValue(N_semi->GetName())); delete tmp_evt;
                
                N_h2mu_mean->setVal(N_h2mu->getVal());
                tmp_evt = N_h2mu_lnn->generate(RooArgSet(*N_h2mu),1);
                N_h2mu_mean->setVal(tmp_evt->get(0)->getRealValue(N_h2mu->GetName())); delete tmp_evt;
                
                effratio_bs_mean->setVal(effratio_bs->getVal());
                tmp_evt = effratio_bs_gau->generate(RooArgSet(*effratio_bs),1);
                effratio_bs_mean->setVal(tmp_evt->get(0)->getRealValue(effratio_bs->GetName())); delete tmp_evt;
                
                effratio_bd_mean->setVal(effratio_bd->getVal());
                tmp_evt = effratio_bd_gau->generate(RooArgSet(*effratio_bd),1);
                effratio_bd_mean->setVal(tmp_evt->get(0)->getRealValue(effratio_bd->GetName())); delete tmp_evt;
            }
        }
        
        RooDataSet *toy_data = new RooDataSet("toy_data","",varlist);
        
        // do category-by-category, component-by-component generation, replace the buggy CBShape
        for (auto& cat: CatMan.cats) {
            GlobalCat->setIndex(cat.index);
            
            RooProdPdf *TauSelCat_pdf_bs = 0, *TauSelCat_pdf_bd = 0, *TauSelCat_pdf_peak = 0, *TauSelCat_pdf_semi = 0, *TauSelCat_pdf_h2mu = 0, *TauSelCat_pdf_comb = 0;
            if (opt.Contains("tausplot")) {
                TauSelCat_pdf_bs   = new RooProdPdf("TauSelCat_pdf_bs",  "", *wspace_gen->pdf(Form("Tau_pdf_bs_%s",cat.id.Data())),
                                                    *wspace_gen->pdf(Form("SelCat_pdf_bs_%s",cat.id.Data())));
                TauSelCat_pdf_bd   = new RooProdPdf("TauSelCat_pdf_bd",  "", *wspace_gen->pdf(Form("Tau_pdf_bd_%s",cat.id.Data())),
                                                    *wspace_gen->pdf(Form("SelCat_pdf_bd_%s",cat.id.Data())));
                TauSelCat_pdf_peak = new RooProdPdf("TauSelCat_pdf_peak","", *wspace_gen->pdf(Form("Tau_pdf_peak_%s",cat.id.Data())),
                                                    *wspace_gen->pdf(Form("SelCat_pdf_peak_%s",cat.id.Data())));
                TauSelCat_pdf_semi = new RooProdPdf("TauSelCat_pdf_semi","", *wspace_gen->pdf(Form("Tau_pdf_semi_%s",cat.id.Data())),
                                                    *wspace_gen->pdf(Form("SelCat_pdf_semi_%s",cat.id.Data())));
                TauSelCat_pdf_h2mu = new RooProdPdf("TauSelCat_pdf_h2mu","", *wspace_gen->pdf(Form("Tau_pdf_h2mu_%s",cat.id.Data())),
                                                    *wspace_gen->pdf(Form("SelCat_pdf_h2mu_%s",cat.id.Data())));
                TauSelCat_pdf_comb = new RooProdPdf("TauSelCat_pdf_comb","", *wspace_gen->pdf(Form("Tau_pdf_comb_%s",cat.id.Data())),
                                                    *wspace_gen->pdf(Form("SelCat_pdf_comb_%s",cat.id.Data())));
            }
            
            RooAbsPdf *pdf_bs = wspace_gen->pdf(Form("pdf_bs_%s", cat.id.Data()));
            RooAbsPdf *pdf_bd = wspace_gen->pdf(Form("pdf_bd_%s", cat.id.Data()));
            RooAbsPdf *pdf_comb = wspace_gen->pdf(Form("pdf_comb_%s", cat.id.Data()));
            RooAbsPdf *pdf_semi = wspace_gen->pdf(Form("pdf_semi_%s", cat.id.Data()));
            RooAbsPdf *pdf_h2mu = wspace_gen->pdf(Form("pdf_h2mu_%s", cat.id.Data()));
            RooAbsPdf *pdf_peak = wspace_gen->pdf(Form("pdf_peak_%s", cat.id.Data()));
            
            int N_bs = RooRandom::randomGenerator()->Poisson(wspace_gen->function(Form("N_bs_formula_%s", cat.id.Data()))->getVal());
            int N_bd = RooRandom::randomGenerator()->Poisson(wspace_gen->function(Form("N_bd_formula_%s", cat.id.Data()))->getVal());
            int N_comb = RooRandom::randomGenerator()->Poisson(wspace_gen->var(Form("N_comb_%s", cat.id.Data()))->getVal());
            int N_semi = RooRandom::randomGenerator()->Poisson(wspace_gen->function(Form("N_semi_formula_%s", cat.id.Data()))->getVal());
            int N_h2mu = RooRandom::randomGenerator()->Poisson(wspace_gen->var(Form("N_h2mu_%s", cat.id.Data()))->getVal());
            int N_peak = RooRandom::randomGenerator()->Poisson(wspace_gen->function(Form("N_peak_formula_%s", cat.id.Data()))->getVal());

            RooDataSet *toy_gen, *toy_tau;
            
            toy_gen = pdf_comb->generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),N_comb);
            if (opt.Contains("tausplot")) {
                toy_tau = TauSelCat_pdf_comb->generate(RooArgSet(*Tau,*SelCat),N_comb);
                toy_gen->merge(toy_tau);
                delete toy_tau;
            }
            toy_gen->addColumn(*GlobalCat);
            toy_data->append(*toy_gen);
            delete toy_gen;

            toy_gen = pdf_semi->generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),N_semi);
            if (opt.Contains("tausplot")) {
                toy_tau = TauSelCat_pdf_semi->generate(RooArgSet(*Tau,*SelCat),N_semi);
                toy_gen->merge(toy_tau);
                delete toy_tau;
            }
            toy_gen->addColumn(*GlobalCat);
            toy_data->append(*toy_gen);
            delete toy_gen;
            
            toy_gen = pdf_h2mu->generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),N_h2mu);
            if (opt.Contains("tausplot")) {
                toy_tau = TauSelCat_pdf_h2mu->generate(RooArgSet(*Tau,*SelCat),N_h2mu);
                toy_gen->merge(toy_tau);
                delete toy_tau;
            }
            toy_gen->addColumn(*GlobalCat);
            toy_data->append(*toy_gen);
            delete toy_gen;
            
            toy_gen = pdf_peak->generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),N_peak);
            if (opt.Contains("tausplot")) {
                toy_tau = TauSelCat_pdf_peak->generate(RooArgSet(*Tau,*SelCat),N_peak);
                toy_gen->merge(toy_tau);
                delete toy_tau;
            }
            toy_gen->addColumn(*GlobalCat);
            toy_data->append(*toy_gen);
            delete toy_gen;
            
            double resample;
            RooAcceptReject AR_bs(*pdf_bs, RooArgSet(*Mass,*ReducedMassRes,*PairCat), *pdf_bs->getGeneratorConfig(), false, 0);
            RooAcceptReject AR_bd(*pdf_bd, RooArgSet(*Mass,*ReducedMassRes,*PairCat), *pdf_bd->getGeneratorConfig(), false, 0);
            
            toy_gen = new RooDataSet("toy_gen","",RooArgSet(*Mass,*ReducedMassRes,*PairCat));
            resample = 0.;
            for (int count=0; count<N_bs; count++)
                toy_gen->add(*AR_bs.generateEvent(N_bs-count,resample));
            if (opt.Contains("tausplot")) {
                toy_tau = TauSelCat_pdf_bs->generate(RooArgSet(*Tau,*SelCat),N_bs);
                toy_gen->merge(toy_tau);
                delete toy_tau;
            }
            toy_gen->addColumn(*GlobalCat);
            toy_data->append(*toy_gen);
            delete toy_gen;
            
            toy_gen = new RooDataSet("toy_gen","",RooArgSet(*Mass,*ReducedMassRes,*PairCat));
            resample = 0.;
            for (int count=0; count<N_bd; count++)
                toy_gen->add(*AR_bd.generateEvent(N_bd-count,resample));
            if (opt.Contains("tausplot")) {
                toy_tau = TauSelCat_pdf_bd->generate(RooArgSet(*Tau,*SelCat),N_bd);
                toy_gen->merge(toy_tau);
                delete toy_tau;
            }
            toy_gen->addColumn(*GlobalCat);
            toy_data->append(*toy_gen);
            delete toy_gen;
            
            if (opt.Contains("tausplot")) {
                delete TauSelCat_pdf_bs;
                delete TauSelCat_pdf_bd;
                delete TauSelCat_pdf_peak;
                delete TauSelCat_pdf_semi;
                delete TauSelCat_pdf_h2mu;
                delete TauSelCat_pdf_comb;
            }
        }
        
        RooArgSet global_ext_constr(*wspace->pdf("fs_over_fu_gau"),*wspace->pdf("one_over_BRBR_gau"));
        RooArgSet minos_list;
        for (auto& POI: POI_list) minos_list.add(*wspace->var(POI));
        
        vector<RooCmdArg> fit_cmd;
        fit_cmd.push_back(ExternalConstraints(global_ext_constr));
        fit_cmd.push_back(Strategy(2));
        fit_cmd.push_back(Extended(true));
        fit_cmd.push_back(NumCPU(NCPU));
        fit_cmd.push_back(Hesse(true));
        fit_cmd.push_back(Optimize(true));
        if (USING_MINUIT2) fit_cmd.push_back(Minimizer("Minuit2","migrad"));
        //fit_cmd.push_back(Offset(true));
        fit_cmd.push_back(Save(true));
        if (do_minos) fit_cmd.push_back(Minos(minos_list));
        
        RooLinkedList fit_linkl;
        for (auto& cmd : fit_cmd) fit_linkl.Add(&cmd);
        RooFitResult *res = fit_with_retry(*wspace->pdf("global_pdf"),*toy_data,fit_linkl,1);
        
        res->Print("v");
        wspace_res->import(*res,Form("fitresult_toy_%d",idx));
        
        const RooArgSet *all_floats = wspace->set("all_floats");
        wspace->saveSnapshot("best_fit",*all_floats);
        
        if (opt.Contains("signif")) {
            // Significance estimation for Bd
            wspace->var("BF_bd")->setVal(0.);
            wspace->var("BF_bd")->setConstant(true);
            RooFitResult *res_zerobd = fit_with_retry(*wspace->pdf("global_pdf"),*toy_data,fit_linkl,1);
            wspace_res->import(*res_zerobd,Form("fitresult_zerobd_toy_%d",idx));
            wspace->loadSnapshot("best_fit");

            // Significance estimation for Bs
            wspace->var("BF_bs")->setVal(0.);
            wspace->var("BF_bs")->setConstant(true);
            RooFitResult *res_zerobs = fit_with_retry(*wspace->pdf("global_pdf"),*toy_data,fit_linkl,1);
            wspace_res->import(*res_zerobs,Form("fitresult_zerobs_toy_%d",idx));
            wspace->loadSnapshot("best_fit");

            delete res_zerobd;
            delete res_zerobs;
            
        }else if (opt.Contains("bdstat")) { // assume the default fit is carried out with BF_bd fixed to the target branching fraction
            // release Bd BF
            wspace->var("BF_bd")->setConstant(false);
            RooFitResult *res_floatbd = fit_with_retry(*wspace->pdf("global_pdf"),*toy_data,fit_linkl,1);
            wspace_res->import(*res_floatbd,Form("fitresult_floatbd_toy_%d",idx));
            wspace->loadSnapshot("best_fit");
            
            delete res_floatbd;
        }
        
        if (opt.Contains("tausplot")) ProduceToyTauSPlot(wspace_res, wspace, toy_data, idx);
        
        delete toy_data;
        delete wspace;
        delete res;
    }
}

void ProducePullPlots(TString path)
{
    cout << ">>> ProducePullPlots() start" << endl;
    
    TFileCollection fcoll("fcoll");
    fcoll.Add(path);
    
    for (auto& POI: POI_list) {
        
        TNtupleD *nt = new TNtupleD("nt","","pull");
        
        for (int file = 0; file < fcoll.GetNFiles(); file++) {
            
            TFileInfo* finfo = (TFileInfo*)fcoll.GetList()->At(file);
            TFile fin_res(finfo->GetCurrentUrl()->GetFile());
            
            RooWorkspace *wspace_res = (RooWorkspace *)fin_res.Get("wspace");
        
            int idx = 0;
            while (true) {
                RooFitResult *res = (RooFitResult*)wspace_res->obj(Form("fitresult_toy_%d",idx));
                if (POI == "EffTau_bs") {
                    res = (RooFitResult*)wspace_res->obj(Form("fitresult_taubs_toy_%d",idx));
                    RooFitResult *res_wl = (RooFitResult*)wspace_res->obj(Form("fitresult_taubs_wl_toy_%d",idx));
                    if (res_wl==NULL || res_wl->status()!=0) break;
                }
                if (res==NULL || res->status()!=0) break;
                idx++;
                
                if (res->status()==0) {
                    RooRealVar* init = (RooRealVar*)res->floatParsInit().find(POI);
                    RooRealVar* final = (RooRealVar*)res->floatParsFinal().find(POI);
                    if (init==NULL || final==NULL) break;
                    
                    if (POI == "EffTau_bs") {
                        double var = (final->getVal()-init->getVal())/final->getError();
                        if (final->getError()<4. && final->getError()>0.1)
                            nt->Fill(&var);
                    }else {
                        double err_hi = fabs(final->getErrorHi());
                        double err_lo = fabs(final->getErrorLo());
                        if (err_hi==0.) err_hi = final->getError();
                        if (err_lo==0.) err_lo = final->getError();
                    
                        double var = (final->getVal()-init->getVal());
                        if (var>0.) var /= err_lo;
                        if (var<0.) var /= err_hi;
                        if (err_hi>0. && err_lo>0.) nt->Fill(&var);
                    }
                }
            }
            
            wspace_res->Delete();
        }
        if (nt->GetEntries()==0) continue;
        
        RooRealVar pull("pull","",-6.,6.);
        
        TString title = Form("CMS simluation");
        
        RooPlot* frame = pull.frame(Bins(30), Title(" "));
        
        RooDataSet *rds_pull = new RooDataSet("rds_pull","",nt,RooArgSet(pull));
        
        RooRealVar Mean("Mean", "", 0., -3.0, 3.0);
        RooRealVar Sigma("Sigma","", 1., 0.5, 2.0);
        RooGaussian model("model","", pull, Mean, Sigma);
        
        model.fitTo(*rds_pull,Minos(true));
        
        rds_pull->plotOn(frame);
        model.plotOn(frame, LineColor(kBlue), LineWidth(3), NumCPU(NCPU), Name("model"));
        
        frame->SetMinimum(0.);
        frame->SetMaximum(frame->GetMaximum()*1.3);
        
        TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetMargin(0.14,0.06,0.13,0.07);
        
        frame->GetYaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitle("Entries");
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetXaxis()->SetLabelOffset(0.01);
        frame->GetXaxis()->SetTitle(Form("Pull [%s]",POI.Data()));
        frame->GetXaxis()->SetTitleSize(0.043);
        frame->GetYaxis()->SetTitleSize(0.043);
        frame->Draw();
        gStyle->SetErrorX(0);
        
        TLatex tex;
        tex.SetTextFont(42);
        tex.SetTextSize(0.035);
        tex.SetTextAlign(11);
        tex.SetNDC();
        tex.DrawLatex(0.14,0.94,title);
        
        canvas->Update();
        
        TLegend *leg1 = new TLegend(0.20,0.82,0.91,0.91);
        leg1->SetTextSize(0.042);
        leg1->SetHeader(Form("Pull fit: #mu = %.2f #pm%.2f, #sigma = %.2f #pm%.2f",
                             Mean.getVal(),Mean.getError(),Sigma.getVal(),Sigma.getError()));
        leg1->SetFillColor(kWhite);
        leg1->SetLineColor(kWhite);
        leg1->Draw();
        
        canvas->Print(Form("fig/pull_%s.pdf",POI.Data()));
        
        delete frame;
        delete canvas;
        
        pull.setRange("center",-3.,3.);
        
        frame = pull.frame(Bins(30), Title(" "));
        
        model.fitTo(*rds_pull,Minos(true),Range("center"));
        
        rds_pull->plotOn(frame);
        model.plotOn(frame, LineColor(kBlue), LineWidth(3), NumCPU(NCPU), Name("model"));
        
        frame->SetMinimum(0.);
        frame->SetMaximum(frame->GetMaximum()*1.3);
        
        canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetMargin(0.14,0.06,0.13,0.07);
        
        frame->GetYaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitle("Entries");
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetXaxis()->SetLabelOffset(0.01);
        frame->GetXaxis()->SetTitle(Form("Pull [%s]",POI.Data()));
        frame->GetXaxis()->SetTitleSize(0.043);
        frame->GetYaxis()->SetTitleSize(0.043);
        frame->Draw();
        gStyle->SetErrorX(0);
        
        tex.DrawLatex(0.14,0.94,title);
        
        canvas->Update();
        
        leg1->SetHeader(Form("Pull fit [-3,+3]: #mu = %.2f #pm%.2f, #sigma = %.2f #pm%.2f",
                             Mean.getVal(),Mean.getError(),Sigma.getVal(),Sigma.getError()));
        leg1->Draw();
        
        canvas->Print(Form("fig/pull_%s_narrow.pdf",POI.Data()));
        
        delete nt;
        delete rds_pull;
        delete leg1;
        delete frame;
        delete canvas;
    }
}

void ProduceMeanErrorPlots(TString path)
{
    cout << ">>> ProduceMeanErrorPlots() start" << endl;
    
    TFileCollection fcoll("fcoll");
    fcoll.Add(path);
    
    for (auto& POI: POI_list) {
        
        vector<double> x_val[4];
        
        for (int file = 0; file < fcoll.GetNFiles(); file++) {
            
            TFileInfo* finfo = (TFileInfo*)fcoll.GetList()->At(file);
            TFile fin_res(finfo->GetCurrentUrl()->GetFile());
            
            RooWorkspace *wspace_res = (RooWorkspace *)fin_res.Get("wspace");
            
            int idx = 0;
            while (true) {
                RooFitResult *res = (RooFitResult*)wspace_res->obj(Form("fitresult_toy_%d",idx));
                if (POI == "EffTau_bs") {
                    res = (RooFitResult*)wspace_res->obj(Form("fitresult_taubs_toy_%d",idx));
                    RooFitResult *res_wl = (RooFitResult*)wspace_res->obj(Form("fitresult_taubs_wl_toy_%d",idx));
                    if (res_wl==NULL || res_wl->status()!=0) break;
                }
                if (res==NULL || res->status()!=0) break;
                idx++;
                
                if (res->status()==0) {
                    RooRealVar* final = (RooRealVar*)res->floatParsFinal().find(POI);
                    if (final==NULL) break;
                    
                    if (POI == "EffTau_bs")
                        if (final->getError()>8. || final->getError()<0.05) continue;
                    
                    x_val[0].push_back(final->getVal());
                    x_val[1].push_back(final->getError());
                    x_val[2].push_back(final->getErrorHi());
                    x_val[3].push_back(final->getErrorLo());
                }
            }
            wspace_res->Delete();
        }
        if (x_val[0].size()==0) continue;
        
        TString title = Form("CMS simluation");
        
        vector<TString> types = {"mean","error","errorhi","errorlo"};
        vector<TString> xtitles = {"Fitted","Error","Error(+)","Error(-)"};
        
        for(int type=0; type<4; type++) {
            
            double x_min    = TMath::MinElement(x_val[type].size(),x_val[type].data());
            double x_max    = TMath::MaxElement(x_val[type].size(),x_val[type].data());
            double x_mean   = TMath::Mean(x_val[type].size(),x_val[type].data());
            double x_rms    = TMath::RMS(x_val[type].size(),x_val[type].data());
            double x_median = TMath::Median(x_val[type].size(),x_val[type].data());
            
            if (POI == "EffTau_bs") {
                if (x_min<-8.) x_min = -8.;
                if (x_max>+8.) x_max = +8.;
            }
            
            TH1D* frame = new TH1D("frame","", 24, x_min-(x_max-x_min)*0.1, x_max+(x_max-x_min)*0.1);
            for(auto& val : x_val[type]) frame->Fill(val);
            
            frame->SetMinimum(0.);
            frame->SetStats(false);
            frame->SetFillColor(50);
            frame->SetMaximum(frame->GetMaximum()*1.3);
            
            TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
            canvas->SetMargin(0.14,0.06,0.13,0.07);
            
            frame->GetYaxis()->SetTitleOffset(1.15);
            frame->GetYaxis()->SetTitle("Entries");
            frame->GetXaxis()->SetTitleOffset(1.15);
            frame->GetXaxis()->SetLabelOffset(0.01);
            frame->GetXaxis()->SetTitle(Form("%s [%s]",xtitles[type].Data(),POI.Data()));
            frame->GetXaxis()->SetTitleSize(0.043);
            frame->GetYaxis()->SetTitleSize(0.043);
            frame->Draw();
            gStyle->SetErrorX(0);
            
            TLatex tex;
            tex.SetTextFont(42);
            tex.SetTextSize(0.035);
            tex.SetTextAlign(11);
            tex.SetNDC();
            tex.DrawLatex(0.14,0.94,title);
            
            canvas->Update();
            
            TLegend *leg1 = new TLegend(0.20,0.82,0.91,0.91);
            leg1->SetTextSize(0.042);
            leg1->SetHeader(Form("median: %.3g, mean: %.3g, rms: %.3g",x_median,x_mean,x_rms));
            leg1->SetFillColor(kWhite);
            leg1->SetLineColor(kWhite);
            leg1->Draw();
            
            canvas->Print(Form("fig/toy_%s_%s.pdf",types[type].Data(),POI.Data()));
            
            delete leg1;
            delete frame;
            delete canvas;
        }
    }
}

void ProduceSignificancePlots(TString path)
{
    cout << ">>> ProduceSignificancePlots() start" << endl;
    
    TFileCollection fcoll("fcoll");
    fcoll.Add(path);
    
    vector<double> x_val[2];
    
    for (int file = 0; file < fcoll.GetNFiles(); file++) {
        
        TFileInfo* finfo = (TFileInfo*)fcoll.GetList()->At(file);
        TFile fin_res(finfo->GetCurrentUrl()->GetFile());
        
        RooWorkspace *wspace_res = (RooWorkspace *)fin_res.Get("wspace");
        
        int idx = 0;
        while (true) {
            RooFitResult *res = (RooFitResult*)wspace_res->obj(Form("fitresult_toy_%d",idx));
            RooFitResult *res_zerobd = (RooFitResult*)wspace_res->obj(Form("fitresult_zerobd_toy_%d",idx));
            RooFitResult *res_zerobs = (RooFitResult*)wspace_res->obj(Form("fitresult_zerobs_toy_%d",idx));
            if (res==NULL || res_zerobd==NULL || res_zerobs==NULL) break;
            
            if (res->status()==0 && res_zerobd->status()==0)
                x_val[0].push_back(sqrt(max(0.,res_zerobd->minNll() - res->minNll())*2.));
        
            if (res->status()==0 && res_zerobs->status()==0)
                x_val[1].push_back(sqrt(max(0.,res_zerobs->minNll() - res->minNll())*2.));
            
            idx++;
        }
        wspace_res->Delete();
    }
    
    TString title = Form("CMS simluation");
    
    vector<TString> types = {"sigbd","sigbs"};
    vector<TString> xtitles = {"Significance (B_{d})","Significance (B_{s})"};
    
    for(int type=0; type<2; type++) {
        if (x_val[type].size()==0) continue;
        
        double x_max    = TMath::MaxElement(x_val[type].size(),x_val[type].data());
        double x_mean   = TMath::Mean(x_val[type].size(),x_val[type].data());
        double x_rms    = TMath::RMS(x_val[type].size(),x_val[type].data());
        double x_median = TMath::Median(x_val[type].size(),x_val[type].data());
        
        TH1D* frame;
        if (x_max<6.) frame = new TH1D("frame","", 32, -0.2, 6.2);
        else frame = new TH1D("frame","", 52, -0.2, 10.2);
        for(auto& val : x_val[type]) frame->Fill(val);
        
        frame->SetMinimum(0.);
        frame->SetStats(false);
        frame->SetFillColor(50);
        frame->SetMaximum(frame->GetMaximum()*1.3);
        
        TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetMargin(0.14,0.06,0.13,0.07);
        
        frame->GetYaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitle("Entries");
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetXaxis()->SetLabelOffset(0.01);
        frame->GetXaxis()->SetTitle(xtitles[type]);
        frame->GetXaxis()->SetTitleSize(0.043);
        frame->GetYaxis()->SetTitleSize(0.043);
        frame->Draw();
        gStyle->SetErrorX(0);
        
        TLatex tex;
        tex.SetTextFont(42);
        tex.SetTextSize(0.035);
        tex.SetTextAlign(11);
        tex.SetNDC();
        tex.DrawLatex(0.14,0.94,title);
        
        canvas->Update();
        
        TLegend *leg1 = new TLegend(0.20,0.82,0.91,0.91);
        leg1->SetTextSize(0.042);
        leg1->SetHeader(Form("median: %.3g, mean: %.3g, rms: %.3g",x_median,x_mean,x_rms));
        leg1->SetFillColor(kWhite);
        leg1->SetLineColor(kWhite);
        leg1->Draw();
        
        canvas->Print(Form("fig/toy_significance_%s.pdf",types[type].Data()));
        
        delete leg1;
        delete frame;
        delete canvas;
    }
}

void ProduceDemoSubPlots(RooWorkspace *wspace)
{
    cout << ">>> ProduceDemoSubPlots() start" << endl;
    
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *ReducedMassRes = wspace->var("ReducedMassRes");
    RooCategory *PairCat = wspace->cat("PairCat");
    RooCategory *GlobalCat = wspace->cat("GlobalCat");
    RooArgSet varlist(*Mass,*ReducedMassRes,*PairCat,*GlobalCat);
    
    RooDataSet *toy_data = new RooDataSet("toy_data","",varlist);
    
    // do category-by-category generation, replace the buggy CBShape
    for (auto& cat: CatMan.cats) {
        GlobalCat->setIndex(cat.index);
        
        RooArgList pdf_list;
        pdf_list.add(*wspace->pdf(Form("pdf_comb_%s", cat.id.Data())));
        pdf_list.add(*wspace->pdf(Form("pdf_semi_%s", cat.id.Data())));
        pdf_list.add(*wspace->pdf(Form("pdf_h2mu_%s", cat.id.Data())));
        pdf_list.add(*wspace->pdf(Form("pdf_peak_%s", cat.id.Data())));
        
        RooArgList N_list;
        N_list.add(*wspace->var(Form("N_comb_%s", cat.id.Data())));
        N_list.add(*wspace->function(Form("N_semi_formula_%s", cat.id.Data())));
        N_list.add(*wspace->var(Form("N_h2mu_%s", cat.id.Data())));
        N_list.add(*wspace->function(Form("N_peak_formula_%s", cat.id.Data())));
        
        RooAddPdf pdf_bkg_sum(Form("pdf_bkg_sum_%s", cat.id.Data()), "", pdf_list, N_list);
        
        RooDataSet *toy_gen = pdf_bkg_sum.generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),Extended(false));
        
        RooAbsReal *N_bs_formula = wspace->function(Form("N_bs_formula_%s", cat.id.Data()));
        RooAbsReal *N_bd_formula = wspace->function(Form("N_bd_formula_%s", cat.id.Data()));
        
        int N_bs = round(N_bs_formula->getVal());
        int N_bd = round(N_bd_formula->getVal());
        
        RooAbsPdf *pdf_bs = wspace->pdf(Form("pdf_bs_%s", cat.id.Data()));
        RooAbsPdf *pdf_bd = wspace->pdf(Form("pdf_bd_%s", cat.id.Data()));
        
        RooAcceptReject AR_bs(*pdf_bs, RooArgSet(*Mass,*ReducedMassRes,*PairCat), *pdf_bs->getGeneratorConfig(), false, 0);
        RooAcceptReject AR_bd(*pdf_bd, RooArgSet(*Mass,*ReducedMassRes,*PairCat), *pdf_bd->getGeneratorConfig(), false, 0);
        
        double resample = 0.;
        for (int count=0; count<N_bs; count++)
            toy_gen->add(*AR_bs.generateEvent(N_bs-count,resample));
        for (int count=0; count<N_bd; count++)
            toy_gen->add(*AR_bd.generateEvent(N_bd-count,resample));
        
        toy_gen->addColumn(*GlobalCat);
        toy_data->append(*toy_gen);
        delete toy_gen;
    }
    
    for (auto& cat: CatMan.cats) {
        
        RooRealVar *Mass = wspace->var("Mass");
        
        TString cut = Form("GlobalCat==%d",cat.index);
        TString title = Form("Category: %s", cat.id.Data());
        
        RooPlot* frame = Mass->frame(Bins(25), Title(" "));
        toy_data->plotOn(frame, Cut(cut), Invisible());
        
        RooAbsPdf *pdf = wspace->pdf(Form("pdf_ext_sum_%s",cat.id.Data()));
        
        double norm = pdf->expectedEvents(*Mass);
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), LineColor(kBlue), LineWidth(3), NumCPU(NCPU), Name("fullpdf"));
        
        norm = wspace->function(Form("N_bs_formula_%s",cat.id.Data()))->getVal();
        pdf = wspace->pdf(Form("pdf_bs_%s",cat.id.Data()));
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kRed), FillStyle(3365), NumCPU(NCPU), Name("bs"));
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kRed), LineWidth(2), LineStyle(1), NumCPU(NCPU));
        
        norm = wspace->function(Form("N_bd_formula_%s",cat.id.Data()))->getVal();
        pdf = wspace->pdf(Form("pdf_bd_%s",cat.id.Data()));
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3344), NumCPU(NCPU), Name("bd"));
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kViolet - 4), LineWidth(2), LineStyle(1), NumCPU(NCPU));
        
        norm = wspace->function(Form("N_peak_formula_%s",cat.id.Data()))->getVal();
        pdf = wspace->pdf(Form("pdf_peak_%s",cat.id.Data()));
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kBlack), LineWidth(3), LineStyle(5), NumCPU(NCPU), Name("peak"));
        
        norm = wspace->function(Form("N_semi_formula_%s",cat.id.Data()))->getVal();
        pdf = wspace->pdf(Form("pdf_semi_%s",cat.id.Data()));
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kGreen -3), LineWidth(4), LineStyle(2), NumCPU(NCPU), Name("semi"));
        
        norm = wspace->var(Form("N_h2mu_%s",cat.id.Data()))->getVal();
        pdf = wspace->pdf(Form("pdf_h2mu_%s",cat.id.Data()));
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kOrange -3), LineWidth(4), LineStyle(3), NumCPU(NCPU), Name("h2mu"));
        
        norm = wspace->var(Form("N_comb_%s",cat.id.Data()))->getVal();
        pdf = wspace->pdf(Form("pdf_comb_%s",cat.id.Data()));
        pdf->plotOn(frame, Normalization(norm, RooAbsReal::NumEvent), DrawOption("L"), LineColor(kBlue - 1), LineWidth(3),  LineStyle(2), NumCPU(NCPU), Name("comb"));
        
        TH1D* hist_data = (TH1D*)toy_data->createHistogram(Form("hist_data_%s",cat.id.Data()), *Mass, Cut(cut), Binning(25, Mass_bound[0], Mass_bound[1]));
        hist_data->Sumw2(false);
        hist_data->SetBinErrorOption(TH1::kPoisson);
        hist_data->SetMarkerStyle(20);
        hist_data->SetLineColor(kBlack);
        
        frame->SetMinimum(0.);
        frame->SetMaximum((hist_data->GetMaximum()+sqrt(hist_data->GetMaximum()*3.))*1.3);
        
        TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetMargin(0.14,0.06,0.13,0.07);
        
        frame->GetYaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetTitle("Entries / 0.04 GeV");
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetXaxis()->SetLabelOffset(0.01);
        frame->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
        frame->GetXaxis()->SetTitleSize(0.043);
        frame->GetYaxis()->SetTitleSize(0.043);
        frame->Draw();
        gStyle->SetErrorX(0);
        hist_data->Draw("Esame");
        
        TLatex tex;
        tex.SetTextFont(42);
        tex.SetTextSize(0.035);
        tex.SetTextAlign(11);
        tex.SetNDC();
        tex.DrawLatex(0.14,0.94,title);
        
        canvas->Update();
        TLegend *leg1 = new TLegend(0.20,0.77,0.91,0.91);
        leg1->SetNColumns(2);
        leg1->SetFillColor(kWhite);
        leg1->SetLineColor(kWhite);
        TLegendEntry *data_entry = new TLegendEntry(hist_data, "toy events", "lep");
        data_entry->SetMarkerStyle(20);
        leg1->AddEntry(data_entry, "toy events", "ep");
        leg1->AddEntry(frame->findObject("fullpdf"),"full PDF","L");
        TLegendEntry *bs_entry = new TLegendEntry(frame->findObject("bs"), "B^{0}_{s}#rightarrow#mu^{+}#mu^{-}", "f");
        bs_entry->SetLineColor(kRed);
        bs_entry->SetFillColor(kRed);
        bs_entry->SetFillStyle(3365);
        leg1->AddEntry(bs_entry,"B^{0}_{s}#rightarrow#mu^{+}#mu^{-}","f");
        TLegendEntry *bd_entry = new TLegendEntry(frame->findObject("bd"), "B^{0}#rightarrow#mu^{+}#mu^{-}", "f");
        bd_entry->SetLineColor(kViolet - 4);
        bd_entry->SetFillColor(kViolet - 4);
        bd_entry->SetFillStyle(3344);
        leg1->AddEntry(bd_entry,"B^{0}#rightarrow#mu^{+}#mu^{-}","f");
        leg1->AddEntry(frame->findObject("comb"),"combinatorial bkg","L");
        leg1->AddEntry(frame->findObject("semi"),"semileptonic bkg","L");
        leg1->AddEntry(frame->findObject("h2mu"),"B#rightarrow h#mu^{+}#mu^{-} bkg","L");
        leg1->AddEntry(frame->findObject("peak"),"peaking bkg","L");
        leg1->Draw();
        
        canvas->Print(Form("fig/proj_toy_bin_%s.pdf",cat.id.Data()));
        
        delete leg1;
        delete hist_data;
        delete frame;
        delete canvas;
    }
    
    delete toy_data;
}

void ProduceDemoTauSPlot(RooWorkspace *wspace)
{
    cout << ">>> ProduceDemoTauSPlot() start" << endl;
    
    enum {_bs, _bd, _peak, _semi, _h2mu, _comb, _nspec};
    vector<TString> specs= {"bs", "bd", "peak", "semi", "h2mu", "comb"};
    
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *ReducedMassRes = wspace->var("ReducedMassRes");
    RooCategory *PairCat = wspace->cat("PairCat");
    RooRealVar *Tau = wspace->var("Tau");
    RooCategory *SelCat = wspace->cat("SelCat");
    RooCategory *GlobalCat = wspace->cat("GlobalCat");
    RooArgSet varlist(*Mass,*ReducedMassRes,*PairCat,*Tau,*SelCat,*GlobalCat);
    
    RooDataSet *toy_data = new RooDataSet("toy_data","",varlist);
    
    // do category-by-category, component-by-component generation, replace the buggy CBShape
    for (auto& cat: CatMan.cats) {
        GlobalCat->setIndex(cat.index);
        
        RooProdPdf *TauSelCat_pdf_bs   = new RooProdPdf("TauSelCat_pdf_bs",  "", *wspace->pdf(Form("Tau_pdf_bs_%s",cat.id.Data())),
                                            *wspace->pdf(Form("SelCat_pdf_bs_%s",cat.id.Data())));
        RooProdPdf *TauSelCat_pdf_bd   = new RooProdPdf("TauSelCat_pdf_bd",  "", *wspace->pdf(Form("Tau_pdf_bd_%s",cat.id.Data())),
                                            *wspace->pdf(Form("SelCat_pdf_bd_%s",cat.id.Data())));
        RooProdPdf *TauSelCat_pdf_peak = new RooProdPdf("TauSelCat_pdf_peak","", *wspace->pdf(Form("Tau_pdf_peak_%s",cat.id.Data())),
                                            *wspace->pdf(Form("SelCat_pdf_peak_%s",cat.id.Data())));
        RooProdPdf *TauSelCat_pdf_semi = new RooProdPdf("TauSelCat_pdf_semi","", *wspace->pdf(Form("Tau_pdf_semi_%s",cat.id.Data())),
                                            *wspace->pdf(Form("SelCat_pdf_semi_%s",cat.id.Data())));
        RooProdPdf *TauSelCat_pdf_h2mu = new RooProdPdf("TauSelCat_pdf_h2mu","", *wspace->pdf(Form("Tau_pdf_h2mu_%s",cat.id.Data())),
                                            *wspace->pdf(Form("SelCat_pdf_h2mu_%s",cat.id.Data())));
        RooProdPdf *TauSelCat_pdf_comb = new RooProdPdf("TauSelCat_pdf_comb","", *wspace->pdf(Form("Tau_pdf_comb_%s",cat.id.Data())),
                                            *wspace->pdf(Form("SelCat_pdf_comb_%s",cat.id.Data())));
        
        RooAbsPdf *pdf_bs = wspace->pdf(Form("pdf_bs_%s", cat.id.Data()));
        RooAbsPdf *pdf_bd = wspace->pdf(Form("pdf_bd_%s", cat.id.Data()));
        RooAbsPdf *pdf_comb = wspace->pdf(Form("pdf_comb_%s", cat.id.Data()));
        RooAbsPdf *pdf_semi = wspace->pdf(Form("pdf_semi_%s", cat.id.Data()));
        RooAbsPdf *pdf_h2mu = wspace->pdf(Form("pdf_h2mu_%s", cat.id.Data()));
        RooAbsPdf *pdf_peak = wspace->pdf(Form("pdf_peak_%s", cat.id.Data()));
        
        int N_bs = (int)(wspace->function(Form("N_bs_formula_%s", cat.id.Data()))->getVal());
        int N_bd = (int)(wspace->function(Form("N_bd_formula_%s", cat.id.Data()))->getVal());
        int N_comb = (int)(wspace->var(Form("N_comb_%s", cat.id.Data()))->getVal());
        int N_semi = (int)(wspace->function(Form("N_semi_formula_%s", cat.id.Data()))->getVal());
        int N_h2mu = (int)(wspace->var(Form("N_h2mu_%s", cat.id.Data()))->getVal());
        int N_peak = (int)(wspace->function(Form("N_peak_formula_%s", cat.id.Data()))->getVal());
        
        RooDataSet *toy_gen, *toy_tau;
        
        toy_gen = pdf_comb->generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),N_comb);
        toy_tau = TauSelCat_pdf_comb->generate(RooArgSet(*Tau,*SelCat),N_comb);
        toy_gen->merge(toy_tau);
        toy_gen->addColumn(*GlobalCat);
        toy_data->append(*toy_gen);
        delete toy_tau; delete toy_gen;
        
        toy_gen = pdf_semi->generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),N_semi);
        toy_tau = TauSelCat_pdf_semi->generate(RooArgSet(*Tau,*SelCat),N_semi);
        toy_gen->merge(toy_tau);
        toy_gen->addColumn(*GlobalCat);
        toy_data->append(*toy_gen);
        delete toy_tau; delete toy_gen;
        
        toy_gen = pdf_h2mu->generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),N_h2mu);
        toy_tau = TauSelCat_pdf_h2mu->generate(RooArgSet(*Tau,*SelCat),N_h2mu);
        toy_gen->merge(toy_tau);
        toy_gen->addColumn(*GlobalCat);
        toy_data->append(*toy_gen);
        delete toy_tau; delete toy_gen;
        
        toy_gen = pdf_peak->generate(RooArgSet(*Mass,*ReducedMassRes,*PairCat),N_peak);
        toy_tau = TauSelCat_pdf_peak->generate(RooArgSet(*Tau,*SelCat),N_peak);
        toy_gen->merge(toy_tau);
        toy_gen->addColumn(*GlobalCat);
        toy_data->append(*toy_gen);
        delete toy_tau; delete toy_gen;
        
        double resample;
        RooAcceptReject AR_bs(*pdf_bs, RooArgSet(*Mass,*ReducedMassRes,*PairCat), *pdf_bs->getGeneratorConfig(), false, 0);
        RooAcceptReject AR_bd(*pdf_bd, RooArgSet(*Mass,*ReducedMassRes,*PairCat), *pdf_bd->getGeneratorConfig(), false, 0);
        
        toy_gen = new RooDataSet("toy_gen","",RooArgSet(*Mass,*ReducedMassRes,*PairCat));
        resample = 0.;
        for (int count=0; count<N_bs; count++)
            toy_gen->add(*AR_bs.generateEvent(N_bs-count,resample));
        toy_tau = TauSelCat_pdf_bs->generate(RooArgSet(*Tau,*SelCat),N_bs);
        toy_gen->merge(toy_tau);
        toy_gen->addColumn(*GlobalCat);
        toy_data->append(*toy_gen);
        delete toy_tau; delete toy_gen;
        
        toy_gen = new RooDataSet("toy_gen","",RooArgSet(*Mass,*ReducedMassRes,*PairCat));
        resample = 0.;
        for (int count=0; count<N_bd; count++)
            toy_gen->add(*AR_bd.generateEvent(N_bd-count,resample));
        toy_tau = TauSelCat_pdf_bd->generate(RooArgSet(*Tau,*SelCat),N_bd);
        toy_gen->merge(toy_tau);
        toy_gen->addColumn(*GlobalCat);
        toy_data->append(*toy_gen);
        delete toy_tau; delete toy_gen;
        
        delete TauSelCat_pdf_bs;
        delete TauSelCat_pdf_bd;
        delete TauSelCat_pdf_peak;
        delete TauSelCat_pdf_semi;
        delete TauSelCat_pdf_h2mu;
        delete TauSelCat_pdf_comb;
    }
    
    // produce sPlot
    TH1D *h_tau = new TH1D("h_tau_bs", "", Tau_bins.size()-1, Tau_bins.data());
    h_tau->Sumw2();
    
    for (auto& cat: CatMan.cats) {
        
        RooArgSet norm(*Mass,*ReducedMassRes,*PairCat);
        
        double yield[_nspec];
        yield[_bs]   = wspace->function(Form("N_bs_formula_%s",cat.id.Data()))->getVal();
        yield[_bd]   = wspace->function(Form("N_bd_formula_%s",cat.id.Data()))->getVal();
        yield[_peak] = wspace->function(Form("N_peak_formula_%s",cat.id.Data()))->getVal();
        yield[_semi] = wspace->function(Form("N_semi_formula_%s",cat.id.Data()))->getVal();
        yield[_h2mu] = wspace->var(Form("N_h2mu_%s",cat.id.Data()))->getVal();
        yield[_comb] = wspace->var(Form("N_comb_%s",cat.id.Data()))->getVal();
        
        TMatrixD covInv(_nspec, _nspec);
        covInv = 0.;
        
        for (int evt=0; evt<toy_data->numEntries(); evt++) {
            const RooArgSet* arg = toy_data->get(evt);
            if (arg->getCatIndex("GlobalCat")!=cat.index) continue;
            
            Mass->setVal(arg->getRealValue("Mass"));
            ReducedMassRes->setVal(arg->getRealValue("ReducedMassRes"));
            PairCat->setIndex(arg->getCatIndex("PairCat"));
            
            double pdf[_nspec];
            for (int idx = 0; idx<_nspec; idx++)
                pdf[idx] = wspace->pdf(Form("pdf_%s_%s",specs[idx].Data(),cat.id.Data()))->getVal(&norm);
            
            double pdf_total = 0.;
            for (int idx = 0; idx<_nspec; idx++) pdf_total += yield[idx]*pdf[idx];
            
            for (int row = 0; row<_nspec; row++)
                for (int col = 0; col<_nspec; col++)
                    covInv(row,col) += pdf[row]*pdf[col]/(pdf_total*pdf_total);
        }
        
        TMatrixD covMatrix(TMatrixD::kInverted,covInv);
        
        for (int evt=0; evt<toy_data->numEntries(); evt++) {
            const RooArgSet* arg = toy_data->get(evt);
            if (arg->getCatIndex("GlobalCat")!=cat.index) continue;
            if (arg->getCatIndex("SelCat")!=1) continue;
            
            Mass->setVal(arg->getRealValue("Mass"));
            ReducedMassRes->setVal(arg->getRealValue("ReducedMassRes"));
            PairCat->setIndex(arg->getCatIndex("PairCat"));
            
            double pdf[_nspec];
            for (int idx = 0; idx<_nspec; idx++)
                pdf[idx] = wspace->pdf(Form("pdf_%s_%s",specs[idx].Data(),cat.id.Data()))->getVal(&norm);
            
            double denominator = 0.;
            for (int idx = 0; idx<_nspec; idx++) denominator += yield[idx]*pdf[idx];
            
            double numerator = 0.;
            for (int idx = 0; idx<_nspec; idx++) numerator += covMatrix(_bs,idx)*pdf[idx];
            
            double weight = numerator/denominator;
            h_tau->Fill(arg->getRealValue("Tau"),weight);
        }
    }
    
    RooDataHist *h_tau_data = new RooDataHist("h_tau_data", "", RooArgList(*Tau), h_tau);
    
    TString title = "CMS Preliminary";
    
    RooPlot* frame = Tau->frame(Title(" "));
    h_tau->SetMarkerStyle(20);
    h_tau->SetLineColor(kBlack);
    h_tau_data->plotOn(frame,MarkerStyle(20),LineColor(kBlack));
    
    RooAbsPdf *Tau_pdf_bs = wspace->pdf("Tau_pdf_bs_mix");
    
    Tau_pdf_bs->plotOn(frame, DrawOption("L"), LineColor(kBlue), LineWidth(2), LineStyle(1), NumCPU(NCPU));
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.14,0.06,0.13,0.07);
    
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetYaxis()->SetTitle("Entries");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("#tau [ps]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->SetStats(false);
    frame->Draw("E");
    gStyle->SetErrorX(0);
    
    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.035);
    tex.SetTextAlign(11);
    tex.SetNDC();
    tex.DrawLatex(0.14,0.94,title);
    
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
    leg1->AddEntry(h_tau, "weighted toy (B_{s})", "lep");
    leg1->Draw();
    
    canvas->Print("fig/tau_splot_toy_bs.pdf");
    
    delete leg1;
    delete h_tau;
    delete h_tau_data;
    delete toy_data;
    delete frame;
    delete canvas;
}

void ApplyYieldScale(RooWorkspace *wspace, double yield_scale)
{
    cout << ">>> ApplyYieldScale() start" << endl;
    
    for (auto& cat: CatMan.cats) {
        
        RooRealVar *N_bu = wspace->var(Form("N_bu_%s", cat.id.Data()));
        RooRealVar *N_bu_mean = wspace->var(Form("N_bu_%s_mean", cat.id.Data()));
        RooRealVar *N_bu_sigma = wspace->var(Form("N_bu_%s_sigma", cat.id.Data()));
        N_bu_mean->setVal(N_bu_mean->getVal()*yield_scale);
        N_bu_sigma->setVal(N_bu_sigma->getVal()*yield_scale);
        N_bu->setMax(N_bu->getMax()*yield_scale);
        N_bu->setVal(N_bu->getVal()*yield_scale);
        N_bu->setMin(N_bu->getMin()*yield_scale);
        
        RooRealVar *N_semi = wspace->var(Form("N_semi_%s", cat.id.Data()));
        RooRealVar *N_semi_mean = wspace->var(Form("N_semi_%s_mean", cat.id.Data()));
        RooRealVar *N_semi_sigma = wspace->var(Form("N_semi_%s_sigma", cat.id.Data()));
        N_semi_mean->setVal(N_semi_mean->getVal()*yield_scale);
        N_semi_sigma->setVal(N_semi_sigma->getVal()*yield_scale);
        N_semi->setMax(N_semi->getMax()*yield_scale);
        N_semi->setVal(N_semi->getVal()*yield_scale);
        N_semi->setMin(N_semi->getMin()*yield_scale);
        
        RooRealVar *N_h2mu = wspace->var(Form("N_h2mu_%s", cat.id.Data()));
        RooRealVar *N_h2mu_mean = wspace->var(Form("N_h2mu_%s_mean", cat.id.Data()));
        N_h2mu_mean->setVal(N_h2mu_mean->getVal()*yield_scale);
        N_h2mu->setMax(N_h2mu->getMax()*yield_scale);
        N_h2mu->setVal(N_h2mu->getVal()*yield_scale);
        N_h2mu->setMin(N_h2mu->getMin()*yield_scale);
        
        RooRealVar *N_peak = wspace->var(Form("N_peak_%s", cat.id.Data()));
        RooRealVar *N_peak_mean = wspace->var(Form("N_peak_%s_mean", cat.id.Data()));
        N_peak_mean->setVal(N_peak_mean->getVal()*yield_scale);
        N_peak->setMax(N_peak->getMax()*yield_scale);
        N_peak->setVal(N_peak->getVal()*yield_scale);
        N_peak->setMin(N_peak->getMin()*yield_scale);
        
        RooRealVar *N_comb = wspace->var(Form("N_comb_%s", cat.id.Data()));
        N_comb->setMax(N_comb->getMax()*yield_scale);
        N_comb->setVal(N_comb->getVal()*yield_scale);
        N_comb->setMin(N_comb->getMin()*yield_scale);
    }
}

void bmm4toystudy(TString commands = "")
{
    // -----------------------------------------------------------
    // parse the commands
    bool do_genfit = false, do_minos = false, do_tausplot = false;
    bool do_make_plots = false, do_make_demo = false;
    bool do_freq = false, do_bdstat = false;
    int seed = 1234;
    int iterations = 10;
    double yield_scale = 1.0;
    TString file_gen = "wspace_prefit.root";
    TString file_fit = "n/a";
    TString file_res = "wspace_toyresult.root";
    vector<TString> set_genpar_names, set_fitpar_names;
    vector<double> set_genpar_values, set_fitpar_values;
    vector<int> set_genpar_states, set_fitpar_states;
    
    cout << ">>> -------------------------" << endl;
    cout << ">>> BMM4 toy study start" << endl;
    cout << ">>> -------------------------" << endl;
    cout << ">>> commands:" << endl;
    cout << ">>> - genfit                                : gen & fit to toy" << endl;
    cout << ">>> - minos                                 : call MINOS for BFs" << endl;
    cout << ">>> - tausplot                              : enable lifetime sPlot study" << endl;
    cout << ">>> - freq                                  : produce frequentist toy (w/ postfit nuisances)" << endl;
    cout << ">>> - bdstat                                : calculate profile likelihood test statistics (F&C study or upper limit for Bd)" << endl;
    cout << ">>> - make_plots                            : produce resulting plots from toy" << endl;
    cout << ">>> - make_demo                             : produce demo projection plots" << endl;
    cout << ">>> - seed=[1234]                           : set the random seed" << endl;
    cout << ">>> - iterations=[10]                       : set # of iterations" << endl;
    cout << ">>> - yield_scale=[1.0]                     : set scaling factor for the yields" << endl;
    cout << ">>> - ws_gen=[wspace_prefit.root]           : set generation workspace" << endl;
    cout << ">>> - ws_fit=[n/a]                          : fit model workspace" << endl;
    cout << ">>> - ws_res=[wspace_toyresult.root]        : resulting workspace file" << endl;
    cout << ">>> - set_genpar [name]=[value]=[float/fix] : set parameter value & state (for gen workspace)" << endl;
    cout << ">>> - set_fitpar [name]=[value]=[float/fix] : set parameter value & state (for fit workspace)" << endl;
    cout << ">>> parsing commands: [" << commands << "]" << endl;
    if (commands=="") return;
    
    TString tok;
    Ssiz_t from = 0;
    while(commands.Tokenize(tok, from, "[ \t;=:]")) {
        if (tok=="genfit")              do_genfit = true;       // gen & fit to toy
        else if (tok=="minos")          do_minos = true;        // call MINOS for BFs
        else if (tok=="tausplot")       do_tausplot = true;     // enable lifetime sPlot study
        else if (tok=="freq")           do_freq = true;         // switch to postfit nuisances
        else if (tok=="bdstat")         do_bdstat = true;       // calculate profile likelihood test statistics for Bd
        else if (tok=="make_plots")     do_make_plots = true;   // produce resulting plots from toy
        else if (tok=="make_demo")      do_make_demo = true;    // produce demo projection plots
        else if (tok=="seed") {                                 // set the random seed
            commands.Tokenize(tok, from, "[ \t;=:]");
            seed = tok.Atoi();
        }else if (tok=="iterations") {                          // set # of iterations
            commands.Tokenize(tok, from, "[ \t;=:]");
            iterations = tok.Atoi();
        }else if (tok=="yield_scale") {                         // set scaling factor for the yields
            commands.Tokenize(tok, from, "[ \t;=:]");
            yield_scale = tok.Atof();
        }else if (tok=="ws_gen") {                              // generation base workspace file
            commands.Tokenize(tok, from, "[ \t;=:]");
            file_gen = tok;
        }else if (tok=="ws_fit") {                              // fit model base workspace file
            commands.Tokenize(tok, from, "[ \t;=:]");
            file_fit = tok;
        }else if (tok=="ws_res") {                              // result workspace file
            commands.Tokenize(tok, from, "[ \t;=:]");
            file_res = tok;
        }else if (tok=="set_genpar") {                          // set parameter value & fix/float for gen workspace
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_genpar_names.push_back(tok);
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_genpar_values.push_back(tok.Atof());
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_genpar_states.push_back(tok=="fix"?1:0); // 1 - fix, 0 - float
        }else if (tok=="set_fitpar") {                          // set parameter value & fix/float for fit workspace
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_fitpar_names.push_back(tok);
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_fitpar_values.push_back(tok.Atof());
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_fitpar_states.push_back(tok=="fix"?1:0); // 1 - fix, 0 - float
        }else {
            cout << ">>> unknown command '" << tok << "'" << endl;
            return;
        }
    }
    
    // set the RooFit seed
    cout << ">>> set the seed to " << seed << endl;
    RooRandom::randomGenerator()->SetSeed(seed);
    
    // -----------------------------------------------------------
    // Initializing categories definition
    
    if (CONFIG_BMM3) CatMan.RegisterBMM3Categories();
    if (CONFIG_BMM4) CatMan.RegisterBMM4Categories();
    CatMan.Print();
    
    if (do_genfit) {
        // workspace for store results
        RooWorkspace *wspace_res = new RooWorkspace("wspace");
    
        // loading the generation model from the given workspace file
        cout << ">>> read the generation model from '" << file_gen << "'" << endl;
        exist_protection(file_gen);
        TFile *fin_gen = new TFile(file_gen);
        RooWorkspace *wspace_gen = (RooWorkspace *)fin_gen->Get("wspace");
        
        wspace_gen->var("PeeK_bs_2016GHs01_0_0")->setVal(wspace_gen->var("PeeK_bs_2016GHs01_0_0")->getVal()*0.794);
        wspace_gen->var("PeeK_bs_2016GHs01_1_0")->setVal(wspace_gen->var("PeeK_bs_2016GHs01_1_0")->getVal()*0.692);
        wspace_gen->var("PeeK_bd_2016GHs01_0_0")->setVal(wspace_gen->var("PeeK_bd_2016GHs01_0_0")->getVal()*0.806);
        wspace_gen->var("PeeK_bd_2016GHs01_1_0")->setVal(wspace_gen->var("PeeK_bd_2016GHs01_1_0")->getVal()*0.696);
        
        // parameter modifications (gen workspace)
        for (int idx = 0; idx<(int)set_genpar_names.size(); idx++) {
            cout << ">>> set gen parameter " << set_genpar_names[idx] << " to " << set_genpar_values[idx] << (set_genpar_states[idx]?" (fixed)":" (floated)") << endl;
            wspace_gen->var(set_genpar_names[idx])->setVal(set_genpar_values[idx]);
            wspace_gen->var(set_genpar_names[idx])->setConstant(set_genpar_states[idx]);
        }
        
        RooWorkspace *wspace_fit = NULL;
        if (file_fit=="n/a" || file_fit==file_gen) {
            cout << ">>> clone the fitting model from the generation model" << endl;
            wspace_fit = new RooWorkspace(*wspace_gen);
        } else {
            cout << ">>> read the fitting model from '" << file_fit << "'" << endl;
            exist_protection(file_fit);
            TFile *fin_fit = new TFile(file_fit);
            wspace_fit = (RooWorkspace *)fin_fit->Get("wspace");
        }
        
        // parameter modifications (fit workspace)
        for (int idx = 0; idx<(int)set_fitpar_names.size(); idx++) {
            cout << ">>> set fit parameter " << set_fitpar_names[idx] << " to " << set_fitpar_values[idx] << (set_fitpar_states[idx]?" (fixed)":" (floated)") << endl;
            wspace_fit->var(set_fitpar_names[idx])->setVal(set_fitpar_values[idx]);
            wspace_fit->var(set_fitpar_names[idx])->setConstant(set_fitpar_states[idx]);
        }
        
        if (yield_scale!=1.0) {
            ApplyYieldScale(wspace_gen, yield_scale);
            ApplyYieldScale(wspace_fit, yield_scale);
        }
        
        TString toy_opt = "";
        if (do_freq) toy_opt += "freq;";
        else toy_opt += "prefit;";
        
        if (do_bdstat) toy_opt += "bdstat;";
        else toy_opt += "signif;";
        
        if (do_tausplot) toy_opt += "tausplot;";
        
        PerformToyStudy(wspace_res, wspace_gen, wspace_fit, iterations, do_minos, toy_opt);
    
        cout << ">>> save the results to '" << file_res << "'" << endl;
        wspace_res->writeToFile(file_res);
        
        delete fin_gen;
    }
    
    if (do_make_plots) {
        cout << ">>> read the results from '" << file_res << "'" << endl;
        ProducePullPlots(file_res);
        ProduceMeanErrorPlots(file_res);
        ProduceSignificancePlots(file_res);
    }
    
    if (do_make_demo && file_gen!="n/a") { // toy projection plots
        cout << ">>> read the generation model from '" << file_gen << "'" << endl;
        exist_protection(file_gen);
        TFile *fin_gen = new TFile(file_gen);
        RooWorkspace *wspace_gen = (RooWorkspace *)fin_gen->Get("wspace");
        
        if (yield_scale!=1.0)
            ApplyYieldScale(wspace_gen, yield_scale);
        
        ProduceDemoSubPlots(wspace_gen);
        ProduceDemoTauSPlot(wspace_gen);
        delete fin_gen;
    }
    
    cout << ">>> BMM4 toy study end." << endl;
}

int main(int argc, char *argv[]) {
    TString cmd = "";
    for (int i=1; i<argc; i++) {
        cmd += TString(argv[i])+";";
    }
    bmm4toystudy(cmd);
}
