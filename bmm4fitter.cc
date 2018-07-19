#include "bmm4common.h"

void BuildGaussianConstaint(RooWorkspace* wspace, TString key, RooAbsReal& mean, RooAbsReal& sigma,
                              double min = 0., double max = 0., bool always_positive = true)
{
    if (min==max) { // auto boundaries, set to +-10 sigma
        min = mean.getVal() - sigma.getVal()*10.;
        max = mean.getVal() + sigma.getVal()*10.;
        if (always_positive && min<0.) min = 0.;
    }
    RooRealVar main_var(key, "", mean.getVal(), min, max);
    main_var.setError(sigma.getVal());
    RooGaussian gauss(Form("%s_gau",key.Data()), "", main_var, mean, sigma);
    
    cout << ">>> BuildGaussianConstaint: " << key << ": mu = " << mean.getVal() << ", sigma = " << sigma.getVal() << " [" << min << ", " << max << "]" << endl;
    wspace->import(gauss);
}

void BuildGaussianConstaint(RooWorkspace* wspace, TString key, double mean, double sigma,
                            double min = 0., double max = 0., bool always_positive = true)
{
    RooRealVar mean_const(Form("%s_mean",key.Data()), "", mean);
    RooRealVar sigma_const(Form("%s_sigma",key.Data()), "", sigma);
    BuildGaussianConstaint(wspace, key, mean_const, sigma_const, min, max, always_positive);
}

void BuildLognormalConstaint(RooWorkspace* wspace, TString key, RooAbsReal& mean, RooAbsReal& kappa,
                              double min = 0., double max = 0., bool always_positive = true)
{
    if (min==max) { // auto boundaries, set to +-10 sigma
        min = mean.getVal() - (kappa.getVal()-1.)*mean.getVal()*10.;
        max = mean.getVal() + (kappa.getVal()-1.)*mean.getVal()*10.;
        if (always_positive && min<0.) min = 0.;
    }
    RooRealVar main_var(key, "", mean.getVal(), min, max);
    main_var.setError((kappa.getVal()-1.)*mean.getVal());
    RooLognormal Lognormal(Form("%s_lnn",key.Data()), "", main_var, mean, kappa);
    
    cout << ">>> BuildLognormalConstaint: " << key << ": mu = " << mean.getVal() << ", kappa = " << kappa.getVal() << " [" << min << ", " << max << "]" << endl;
    wspace->import(Lognormal);
}

void BuildLognormalConstaint(RooWorkspace* wspace, TString key, double mean, double kappa,
                            double min = 0., double max = 0., bool always_positive = true)
{
    RooRealVar mean_const(Form("%s_mean",key.Data()), "", mean);
    RooRealVar kappa_const(Form("%s_kappa",key.Data()), "", kappa);
    BuildLognormalConstaint(wspace, key, mean_const, kappa_const, min, max, always_positive);
}

void PrepareGlobalVariables(RooWorkspace* wspace, RooWorkspace* wspace_base = 0)
{
    cout << ">>> PrepareGlobalVariables() start" << endl;
    
    if (wspace_base!=0) {
        cout << ">>> import from existing workspace." << endl;
        wspace->import(*wspace_base->var("BF_bs"));
        wspace->import(*wspace_base->var("BF_bd"));
        wspace->import(*wspace_base->var("Mass"));
        wspace->import(*wspace_base->var("ReducedMassRes"));
        wspace->import(*wspace_base->var("BDT"));
        wspace->import(*wspace_base->var("Tau"));
        wspace->import(*wspace_base->var("TauRes"));
        wspace->import(*wspace_base->cat("GlobalCat"));
        wspace->import(*wspace_base->cat("SpCat"));
        wspace->import(*wspace_base->cat("SelCat"));
        wspace->import(*wspace_base->cat("PairCat"));
        wspace->import(*wspace_base->var("Weight"));
        wspace->import(*wspace_base->var("dblmu_corr_scale"));
        wspace->import(*wspace_base->pdf("one_over_BRBR_gau"));
        wspace->import(*wspace_base->pdf("fs_over_fu_gau"));
        wspace->import(*wspace_base->pdf("fs_over_fu_S13_gau"));
        wspace->import(*wspace_base->var("EffTau_bs"));
        wspace->import(*wspace_base->var("EffTau_bd"));
        return;
    }
    
    cout << ">>> prepare new variables & constraints." << endl;
    
    vector<TString> keys;
    vector<double> values;
    
    keys.push_back("fsfu:val");
    keys.push_back("fsfu:err");
    keys.push_back("fsfu_S13:val");
    keys.push_back("fsfu_S13:err");
    keys.push_back("BF_BuToKpsiK:val");
    keys.push_back("BF_BuToKpsiK:err");
    keys.push_back("BF_JpsiToMuMu:val");
    keys.push_back("BF_JpsiToMuMu:err");
    keys.push_back("BF_BsToMuMu:val");
    keys.push_back("BF_BdToMuMu:val");
    keys.push_back("EffTau_BsToMuMu:val");
    keys.push_back("EffTau_BdToMuMu:val");
    
    enum {
        _fsfu, _fsfu_err,
        _fsfu_S13, _fsfu_S13_err,
        _BF_bu2jpsik, _BF_bu2jpsik_err,
        _BF_jpsi2mumu, _BF_jpsi2mumu_err,
        _BF_bs2mumu, _BF_bd2mumu,
        _EffTau_bs2mumu, _EffTau_bd2mumu
    };
    
    ReadValuesFromTex("input/external_parameters.tex",keys,values);
    
    RooRealVar BF_bs("BF_bs", "", values[_BF_bs2mumu], 0., 1e-8);
    RooRealVar BF_bd("BF_bd", "", values[_BF_bd2mumu], 0., 1e-8);
    
    wspace->import(BF_bs);
    wspace->import(BF_bd);
    
    double one_over_BRBR_val = 1./ (values[_BF_bu2jpsik] * values[_BF_jpsi2mumu]);
    double one_over_BRBR_err = one_over_BRBR_val *
        sqrt(pow(values[_BF_bu2jpsik_err]/values[_BF_bu2jpsik],2) + pow(values[_BF_jpsi2mumu_err]/values[_BF_jpsi2mumu],2));

    BuildGaussianConstaint(wspace, "one_over_BRBR", one_over_BRBR_val, one_over_BRBR_err);
    BuildGaussianConstaint(wspace, "fs_over_fu", values[_fsfu], values[_fsfu_err]);
    BuildGaussianConstaint(wspace, "fs_over_fu_S13", values[_fsfu_S13], values[_fsfu_S13_err]);
    
    // Observables
    RooRealVar Mass("Mass", "", Mass_bound[0], Mass_bound[1]);
    RooRealVar ReducedMassRes("ReducedMassRes", "", ReducedMassRes_bound[0], ReducedMassRes_bound[1]);
    RooRealVar BDT("BDT", "", BDT_bound[0], BDT_bound[1]);
    RooRealVar Tau("Tau", "", Tau_bound[0], Tau_bound[1]); // in unit of ps
    RooRealVar TauRes("TauRes", "", TauRes_bound[0], TauRes_bound[1]); // in unit of ps

    // Event category
    RooCategory GlobalCat("GlobalCat", "");
    GlobalCat.defineType("undefined",-1); // uncategorized events
    for (auto& cat: CatMan.cats)
        GlobalCat.defineType(cat.id,cat.index);
    
    // Internal observables
    RooCategory SpCat("SpCat", ""); // sample index for internal use
    for (int i=0;i<max(bmm3::ndecays,bmm4::ndecays);i++)
        SpCat.defineType(Form("sp%d",i),i);
    RooRealVar Weight("Weight", "", 0., 1E10); // event weighting for internal use
    RooCategory SelCat("SelCat", ""); // internal selection (applied for tau & taures boundaries so far)
    SelCat.defineType("failed",0);
    SelCat.defineType("passed",1);
    
    wspace->import(Mass);
    wspace->import(ReducedMassRes);
    wspace->import(BDT);
    wspace->import(Tau);
    wspace->import(TauRes);
    wspace->import(GlobalCat);
    wspace->import(SpCat);
    wspace->import(SelCat);
    wspace->import(Weight);
    
    RooRealVar dblmu_corr_scale("dblmu_corr_scale", "", 1.0, 0.1, 3.0); // in-situ "double-mu" correction
    //dblmu_corr_scale.setConstant(true);
    wspace->import(dblmu_corr_scale);
    
    RooCategory PairCat("PairCat",""); // seagull vs. cowboy
    PairCat.defineType("seagull",0);
    PairCat.defineType("cowboy",1);
    wspace->import(PairCat);
    
    RooRealVar EffTau_bs("EffTau_bs", "", values[_EffTau_bs2mumu], 0.1, 5.0); // in unit of ps, PDG: BH: 1.610 +- 0.012 ps / BL: 1.422 +- 0.008 ps
    RooRealVar EffTau_bd("EffTau_bd", "", values[_EffTau_bs2mumu], 0.1, 5.0); // in unit of ps, PDG: 1.520 +- 0.004 ps
    EffTau_bs.setConstant(true);
    EffTau_bd.setConstant(true);
    wspace->import(EffTau_bs);
    wspace->import(EffTau_bd);
}

void PrepareBMM3SubVariables(RooWorkspace* wspace, RooWorkspace* wspace_base = 0)
{
    cout << ">>> PrepareBMM3SubVariables() start" << endl;
    
    if (wspace_base!=0) {
        cout << ">>> import from existing workspace." << endl;
        for (auto& cat: CatMan.cats) {
            if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
            wspace->import(*wspace_base->pdf(Form("N_bu_%s_gau",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("effratio_bs_%s_gau",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("effratio_bd_%s_gau",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("N_peak_%s_lnn",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("N_semi_%s_gau",cat.id.Data()))); // ISSUE: bmm4 takes lnn for semi background
            wspace->import(*wspace_base->pdf(Form("N_h2mu_%s_lnn",cat.id.Data())));
            wspace->import(*wspace_base->var(Form("N_comb_%s",cat.id.Data())));
        }
        return;
    }
    
    cout << ">>> build variables & constraints." << endl;
    
    for (auto& cat: CatMan.cats) {
        if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
        
        vector<TString> keys;
        vector<double> values;
        vector<double> errors;
        
        keys.push_back(Form("N-OBS-BPLUS%d:val",cat.region));
        keys.push_back(Form("N-OBS-BPLUS%d:tot",cat.region));
        keys.push_back(Form("N-EFF-TOT-BSMM%d:val",cat.region));
        keys.push_back(Form("N-EFF-TOT-BSMM%d:tot",cat.region));
        keys.push_back(Form("N-EFF-TOT-BDMM%d:val",cat.region));
        keys.push_back(Form("N-EFF-TOT-BDMM%d:tot",cat.region));
        keys.push_back(Form("N-EFF-TOT-BPLUS%d:val",cat.region));
        keys.push_back(Form("N-EFF-TOT-BPLUS%d:tot",cat.region));
        
        enum {
            _N_bu, _N_bu_err,
            _eff_bs, _eff_bs_err,
            _eff_bd, _eff_bd_err,
            _eff_bu, _eff_bu_err
        };
        
        ReadValuesFromTex(Form("input/%s/anaBmm.plotResults.%s-cat2%d.tex",cat.era.Data(),cat.era.Data(),cat.bdt_bin),keys,values);
        
        if (values[_N_bu_err]>values[_N_bu]) values[_N_bu_err] = values[_N_bu]*0.05; // hotfix for N_bu uncertainties
        BuildGaussianConstaint(wspace, Form("N_bu_%s",cat.id.Data()),values[_N_bu],values[_N_bu_err]);
        
        double eff_rel_err = 0.;
        
        vector<TString> sys_keys;
        vector<double> sys_values;
        vector<double> sys_errors;
        
        sys_keys.push_back(Form("Acceptance_%d",cat.region));
        sys_keys.push_back(Form("MassScale_%d",cat.region));
        sys_keys.push_back(Form("KaonTrack_%d",cat.region));
        sys_keys.push_back(Form("Trigger_%d",cat.region));
        sys_keys.push_back(Form("MuonID_%d",cat.region));
        
        enum { _Acceptance, _MassScale, _KaonTrack, _Trigger, _MuonID };
        
        ReadValuesFromPlainText("input/external_numbers.txt",sys_keys,sys_values,sys_errors);
        
        eff_rel_err += pow(sys_values[_Acceptance],2);
        eff_rel_err += pow(sys_values[_MassScale],2);
        eff_rel_err += pow(sys_values[_KaonTrack],2);
        eff_rel_err += pow(sys_values[_Trigger],2);
        eff_rel_err += pow(sys_values[_MuonID],2);
        
        sys_keys.clear();
        sys_keys.push_back(Form("relDeltaEpsNoDataNoMcchan%d:val",cat.region));
        sys_keys.push_back(Form("relDeltaEpsCsDataCsMcchan%d:val",cat.region));
        
        enum { _NO_err, _CS_err };
        
        ReadValuesFromTex(Form("input/%s/anaBmm.plotReducedOverlays.%s.tex",cat.era.Data(),cat.era.Data()),sys_keys,sys_values);
        
        eff_rel_err += pow(sys_values[_NO_err],2);
        eff_rel_err += pow(sys_values[_CS_err],2);
        
        eff_rel_err = sqrt(eff_rel_err);
        double effratio_bs = values[_eff_bs]/values[_eff_bu];
        double effratio_bd = values[_eff_bd]/values[_eff_bu];
        double effratio_bs_err = eff_rel_err * effratio_bs;
        double effratio_bd_err = eff_rel_err * effratio_bd;
        
        BuildGaussianConstaint(wspace, Form("effratio_bs_%s",cat.id.Data()),effratio_bs,effratio_bs_err);
        BuildGaussianConstaint(wspace, Form("effratio_bd_%s",cat.id.Data()),effratio_bd,effratio_bd_err);
        
        vector<TString> types = {"Peak", "Rsl", "Comb"};
        vector<TString> fields = {"val", "e1", "e2"};

        enum { _N_peak, _N_semi, _N_comb };
        double N_val[types.size()], N_err[types.size()];
        
        for (unsigned int i=0; i<types.size(); i++) {
            TString& type = types[i];
            
            keys.clear();
            for (auto& field : fields) {
                keys.push_back(Form("Bg%sLo%d:%s",type.Data(),cat.region,field.Data()));
                keys.push_back(Form("Bg%sBd%d:%s",type.Data(),cat.region,field.Data()));
                keys.push_back(Form("Bg%sBs%d:%s",type.Data(),cat.region,field.Data()));
                keys.push_back(Form("Bg%sHi%d:%s",type.Data(),cat.region,field.Data()));
            }
            ReadValuesFromTex(Form("input/%s/anaBmm.plotResults.%s-cat2%d.tex",cat.era.Data(),cat.era.Data(),cat.bdt_bin),keys,values);
        
            double val = 0., e1 = 0., e2 = 0.;
            for (int j=0; j< 4; j++) {
                val += values[j];
                e1  += values[j+4];
                e2  += values[j+8];
            }
            N_val[i] = val;
            N_err[i] = sqrt(e1*e1 + e2*e2);
        }
        
        BuildLognormalConstaint(wspace, Form("N_peak_%s",cat.id.Data()),N_val[_N_peak],1.+N_err[_N_peak]/N_val[_N_peak]);
        
        RooRealVar N_comb(Form("N_comb_%s",cat.id.Data()), "", N_val[_N_comb], 0., N_val[_N_comb]+N_err[_N_comb]*10.);
        wspace->import(N_comb);
        
        double N_semi = 0., N_semi_err = 0.;
        // read the semileptonic background yields (process with 1 real muon)
        for (int idx_sample=bmm3::_bgBs2KMuNu; idx_sample<=bmm3::_bgLb2PMuNu; idx_sample++) {
            keys.clear();
            keys.push_back(Form("%s:loSideband%d:val",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bdRare%d}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bsRare%d}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:hiSideband%d:val",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:loSideband%d:err",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bdRare%dE}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bsRare%dE}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:hiSideband%d:err",bmm3::decays[idx_sample].Data(),cat.region));
            
            ReadValuesFromTex(Form("input/%s/anaBmm.plotResults.%s-cat2%d.tex",cat.era.Data(),cat.era.Data(),cat.bdt_bin),keys,values);
            
            N_semi += values[0]+values[1]+values[2]+values[3]; // yields
            N_semi_err += pow(values[4]+values[5]+values[6]+values[7],2); // error
        }
        N_semi_err = sqrt(N_semi_err);
        
        double N_h2mu = 0., N_h2mu_err = 0.;
        // read the hmumu background yields (process with 2 real muons)
        for (int idx_sample=bmm3::_bgBu2PiMuMu; idx_sample<=bmm3::_bgBs2MuMuGamma; idx_sample++) {
            keys.clear();
            keys.push_back(Form("%s:loSideband%d:val",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bdRare%d}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bsRare%d}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:hiSideband%d:val",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:loSideband%d:err",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bdRare%dE}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bsRare%dE}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:hiSideband%d:err",bmm3::decays[idx_sample].Data(),cat.region));
            
            ReadValuesFromTex(Form("input/%s/anaBmm.plotResults.%s-cat2%d.tex",cat.era.Data(),cat.era.Data(),cat.bdt_bin),keys,values);
            
            N_h2mu += values[0]+values[1]+values[2]+values[3]; // yields
            N_h2mu_err += pow(values[4]+values[5]+values[6]+values[7],2); // error
        }
        N_h2mu_err = sqrt(N_h2mu_err);
        
        BuildGaussianConstaint(wspace, Form("N_semi_%s",cat.id.Data()), N_semi, N_semi_err);
        BuildLognormalConstaint(wspace, Form("N_h2mu_%s",cat.id.Data()), N_h2mu, 1.+N_h2mu_err/N_h2mu);
    }
}

void PrepareBMM4SubVariables(RooWorkspace* wspace, RooWorkspace* wspace_base = 0)
{
    cout << ">>> PrepareBMM4SubVariables() start" << endl;
    
    if (wspace_base!=0) {
        cout << ">>> import from existing workspace." << endl;
        for (auto& cat: CatMan.cats) {
            wspace->import(*wspace_base->pdf(Form("N_bu_%s_gau",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("effratio_bs_%s_gau",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("effratio_bd_%s_gau",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("N_peak_%s_lnn",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("N_semi_%s_lnn",cat.id.Data())));
            wspace->import(*wspace_base->pdf(Form("N_h2mu_%s_lnn",cat.id.Data())));
            wspace->import(*wspace_base->var(Form("N_comb_%s",cat.id.Data())));
            wspace->import(*wspace_base->var(Form("DeltaMass_%s",cat.id.Data())));
            wspace->import(*wspace_base->function(Form("ScaledMass_%s",cat.id.Data())));
        }
        return;
    }
    
    cout << ">>> build variables & constraints." << endl;
    
    for (auto& cat: CatMan.cats) {
        
        cout << ">>> Prepare for category: " << cat.id << endl;
        
        vector<TString> keys;
        vector<double> values;
        
        keys.push_back(Form("%s:DeltaMass-chan%d:val",cat.era.Data(),cat.region));
        ReadValuesFromTex("input/external_parameters.tex",keys,values);

        RooRealVar *Mass = wspace->var("Mass");
        RooRealVar DeltaMass(Form("DeltaMass_%s", cat.id.Data()), "", values[0]); // data-MC mass difference
        RooFormulaVar ScaledMass(Form("ScaledMass_%s", cat.id.Data()), "", "@0-@1", RooArgList(*Mass, DeltaMass));
        wspace->import(ScaledMass);
        
        TString TexSource = Form("input/bmm4/scanBDT-%s.tex",cat.era.Data());
        int bdt_min = (int)(cat.bdt_min*100.);
        int bdt_max = (int)(cat.bdt_max*100.);
        
        TexVar N_bu(TexSource, Form("bdt_%d_%s:N-W8OBS-bupsik-chan%d",bdt_min,cat.era.Data(),cat.region));
        TexVar eff_bs(TexSource, Form("bdt_%d_%s:EFF-TOT-bsmm-chan%d",bdt_min,cat.era.Data(),cat.region));
        TexVar eff_bd(TexSource, Form("bdt_%d_%s:EFF-TOT-bdmm-chan%d",bdt_min,cat.era.Data(),cat.region));
        TexVar eff_bu(TexSource, Form("bdt_%d_%s:EFF-TOT-bupsik-chan%d",bdt_min,cat.era.Data(),cat.region));
        
        if (bdt_max<100) {
            N_bu.SubVar(TexVar(TexSource, Form("bdt_%d_%s:N-W8OBS-bupsik-chan%d",bdt_max,cat.era.Data(),cat.region)));
            eff_bs.SubVar(TexVar(TexSource, Form("bdt_%d_%s:EFF-TOT-bsmm-chan%d",bdt_max,cat.era.Data(),cat.region)));
            eff_bd.SubVar(TexVar(TexSource, Form("bdt_%d_%s:EFF-TOT-bdmm-chan%d",bdt_max,cat.era.Data(),cat.region)));
            eff_bu.SubVar(TexVar(TexSource, Form("bdt_%d_%s:EFF-TOT-bupsik-chan%d",bdt_max,cat.era.Data(),cat.region)));
        }
        
        BuildGaussianConstaint(wspace, Form("N_bu_%s",cat.id.Data()),N_bu.val,N_bu.etot);
        
        double effratio_bs = eff_bs.val/eff_bu.val;
        double effratio_bd = eff_bd.val/eff_bu.val;
        double effratio_bs_err = pow(eff_bs.estat/eff_bs.val,2);
        double effratio_bd_err = pow(eff_bd.estat/eff_bd.val,2);
        
        double eff_rel_err = 0.; // common relative uncertainties
        
        keys.clear();
        keys.push_back(Form("effmuid4r_systematics:sys")); // muon ID uncertainty for efficiency ratios
        keys.push_back(Form("efftrig4r_systematics:sys")); // trigger uncertainty for efficiency ratios
        keys.push_back(Form("efftrack_systematics:sys")); // kaon tracking efficiency
        keys.push_back(Form("%s:yieldcorr_systematics:sys",cat.era.Data())); // yield instability correction
        ReadValuesFromTex("input/external_parameters.tex",keys,values);
        
        eff_rel_err += pow(values[0],2); // muon ID
        eff_rel_err += pow(values[1],2); // trigger
        eff_rel_err += pow(values[2],2); // kaon
        eff_rel_err += pow(values[3],2); // yield instability
        
        keys.clear();
        keys.push_back(Form("%s:acceptance_systematics_chan%d:sys",cat.era.Data(),cat.region)); // acceptance
        keys.push_back(Form("%s:effana4r_systematics_chan%d:sys",cat.era.Data(),cat.region)); // analysis efficiency
        keys.push_back(Form("%s:effcandbupsik_systematics_chan%d:sys",cat.era.Data(),cat.region)); // candidate selection: B+->J/psiK+
        keys.push_back(Form("%s:effcandbsmm_systematics_chan%d:sys",cat.era.Data(),cat.region)); // candidate selection: Bs->mumu
        keys.push_back(Form("%s:effcandbdmm_systematics_chan%d:sys",cat.era.Data(),cat.region)); // candidate selection: Bd->mumu
        ReadValuesFromTex(Form("input/bmm4/plotSystematics.%s.tex",cat.era.Data()),keys,values);
        
        eff_rel_err += pow(values[0],2); // acceptance
        eff_rel_err += pow(values[1],2); // analysis efficiency
        eff_rel_err += pow(values[2],2); // B+->J/psiK+ candidate
        
        effratio_bs_err += eff_rel_err;
        effratio_bd_err += eff_rel_err;
        
        effratio_bs_err += pow(values[3],2); // Bs->mumu candidate
        effratio_bd_err += pow(values[4],2); // Bd->mumu candidate
        
        // Adding effective lifetime correction on efficiencies & errors
        keys.clear();
        keys.push_back(Form("%s:effTau_efficienecy_correction:val",cat.id.Data())); // correction factor
        keys.push_back(Form("%s:effTau_efficienecy_correction:sys",cat.id.Data())); // systematic uncertainty
        ReadValuesFromTex(binsetup_parameter,keys,values);
        effratio_bs *= values[0];
        effratio_bs_err += pow(values[1],2);
        
        // convert relative error to absolute error
        effratio_bs_err = sqrt(effratio_bs_err)*effratio_bs;
        effratio_bd_err = sqrt(effratio_bd_err)*effratio_bd;
        
        BuildGaussianConstaint(wspace, Form("effratio_bs_%s",cat.id.Data()),effratio_bs,effratio_bs_err);
        BuildGaussianConstaint(wspace, Form("effratio_bd_%s",cat.id.Data()),effratio_bd,effratio_bd_err);
        
        // peaking background
        TexVar N_peak;
        for (int idx_sample=bmm4::_bskkMcBg; idx_sample<=bmm4::_lbpkMcBg; idx_sample++) {
            TexVar var;
            for (int bin=0; bin<=3; bin++) {
                TexVar v(TexSource, Form("bdt_%d_%s:N-SCALEDYIELD-MBIN%d-%s-chan%d",bdt_min,cat.era.Data(),bin,bmm4::texlabels[idx_sample].Data(),cat.region));
                if (bdt_max<100) v.SubVar(TexVar(TexSource, Form("bdt_%d_%s:N-SCALEDYIELD-MBIN%d-%s-chan%d",
                                                                 bdt_max,cat.era.Data(),bin,bmm4::texlabels[idx_sample].Data(),cat.region)));
                var.AddVar(v);
            }
            cout << ">>> Expected yield for " << cat.id << ": " << bmm4::texlabels[idx_sample] << ": " << var.val << " +- " << var.estat << " +- " << var.esyst << endl;
            N_peak.AddVar(var);
        }
        
        BuildLognormalConstaint(wspace, Form("N_peak_%s",cat.id.Data()),N_peak.val,1.+N_peak.etot/N_peak.val);
        
        // semileptonic background
        TexVar N_semi;
        for (int idx_sample=bmm4::_bskmunuMcBg; idx_sample<=bmm4::_lbpmunuMcBg; idx_sample++) {
            TexVar var;
            for (int bin=0; bin<=3; bin++) {
                TexVar v(TexSource, Form("bdt_%d_%s:N-SCALEDYIELD-MBIN%d-%s-chan%d",bdt_min,cat.era.Data(),bin,bmm4::texlabels[idx_sample].Data(),cat.region));
                if (bdt_max<100) v.SubVar(TexVar(TexSource, Form("bdt_%d_%s:N-SCALEDYIELD-MBIN%d-%s-chan%d",
                                                                 bdt_max,cat.era.Data(),bin,bmm4::texlabels[idx_sample].Data(),cat.region)));
                var.AddVar(v);
            }
            cout << ">>> Expected yield for " << cat.id << ": " << bmm4::texlabels[idx_sample] << ": " << var.val << " +- " << var.estat << " +- " << var.esyst << endl;
            N_semi.AddVar(var);
        }
        
        BuildLognormalConstaint(wspace, Form("N_semi_%s",cat.id.Data()),N_semi.val,1.+N_semi.etot/N_semi.val);
        
        // h2mu background
        TexVar N_h2mu;
        for (int idx_sample=bmm4::_bdpimumuMcBg; idx_sample<=bmm4::_bupimumuMcBg; idx_sample++) {
            TexVar var;
            for (int bin=0; bin<=3; bin++) {
                TexVar v(TexSource, Form("bdt_%d_%s:N-SCALEDYIELD-MBIN%d-%s-chan%d",bdt_min,cat.era.Data(),bin,bmm4::texlabels[idx_sample].Data(),cat.region));
                if (bdt_max<100) v.SubVar(TexVar(TexSource, Form("bdt_%d_%s:N-SCALEDYIELD-MBIN%d-%s-chan%d",
                                                                 bdt_max,cat.era.Data(),bin,bmm4::texlabels[idx_sample].Data(),cat.region)));
                var.AddVar(v);
            }
            cout << ">>> Expected yield for " << cat.id << ": " << bmm4::texlabels[idx_sample] << ": " << var.val << " +- " << var.estat << " +- " << var.esyst << endl;
            N_h2mu.AddVar(var);
        }
        
        BuildLognormalConstaint(wspace, Form("N_h2mu_%s",cat.id.Data()),N_h2mu.val,1.+N_h2mu.etot/N_h2mu.val);
        
        // combinatorial background
        TexVar N_comb;
        for (int bin=0; bin<=3; bin++) {
            TexVar v(TexSource, Form("bdt_%d_%s:N-FIT-MBIN%d-CB-chan%d",bdt_min,cat.era.Data(),bin,cat.region));
            if (bdt_max<100) v.SubVar(TexVar(TexSource, Form("bdt_%d_%s:N-FIT-MBIN%d-CB-chan%d",bdt_max,cat.era.Data(),bin,cat.region)));
            N_comb.AddVar(v);
        }
        if (N_comb.val<1.) N_comb.val = 1.; // ISSUE: hot fix for too small background
        wspace->import(RooRealVar(Form("N_comb_%s",cat.id.Data()), "", N_comb.val, 0., max(10.,N_comb.val+sqrt(N_comb.val)*10.)));
    }
}

void LoadDataFromBMM3Tree(RooWorkspace* wspace, RooAbsData *dataset, TString filename, TString era)
{
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *ReducedMassRes = wspace->var("ReducedMassRes");
    RooRealVar *BDT = wspace->var("BDT");
    RooRealVar *Tau = wspace->var("Tau");
    RooRealVar *TauRes = wspace->var("TauRes");
    RooCategory *SelCat = wspace->cat("SelCat");
    RooCategory *PairCat = wspace->cat("PairCat");
    RooCategory *GlobalCat = wspace->cat("GlobalCat");
    RooArgSet varlist(*Mass,*ReducedMassRes,*BDT,*Tau,*TauRes,*SelCat,*PairCat,*GlobalCat);
    
    cout << ">>> Load data from tree: " << filename << endl;
    
    exist_protection(filename);
    TFile *fin = new TFile(filename);
    TTree *tin = (TTree*)fin->Get("SgData_bdt");
    
    double m1eta_t, m2eta_t, m_t, bdt_t, me_t, tau_t;
    bool muid_t;
    tin->SetBranchAddress("m1eta", &m1eta_t);
    tin->SetBranchAddress("m2eta", &m2eta_t);
    tin->SetBranchAddress("m", &m_t);
    tin->SetBranchAddress("me", &me_t);
    tin->SetBranchAddress("bdt", &bdt_t);
    tin->SetBranchAddress("tau", &tau_t);
    tin->SetBranchAddress("muid",  &muid_t);
    
    cout << ">>> parsing " << tin->GetEntries() << " entries." << endl;
    
    for (int evt=0; evt<tin->GetEntries();evt++) {
        tin->GetEntry(evt);
        
        bool isBarrel = fabs(m1eta_t)<1.4 && fabs(m2eta_t)<1.4;
        if (isBarrel) m_t += 0.006; // shift mass to MC
        else m_t += 0.007;
        
        if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
        if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
        if (bdt_t < -1.) continue;
            
        int index = CatMan.index(era, isBarrel?0:1, bdt_t);
        if (index<0) continue;
        
        if (tau_t*1E12 > Tau_bound[0] && tau_t*1E12 < Tau_bound[1]) SelCat->setIndex(1); // ISSUE: no tau resolution information
        else SelCat->setIndex(0);
        
        Mass->setVal(m_t);
        ReducedMassRes->setVal(me_t/m_t);
        BDT->setVal(bdt_t);
        Tau->setVal(tau_t*1E12);
        TauRes->setVal(0.); // ISSUE: no tau resolution information 
        PairCat->setIndex(evt%2); // ISSUE: no seagull/cowboy information
        GlobalCat->setIndex(index);

        dataset->add(varlist);
    }
    fin->Close();
}

void LoadDataFromBMM4Tree(RooWorkspace* wspace, RooAbsData *dataset, TString filename, TString era)
{
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *ReducedMassRes = wspace->var("ReducedMassRes");
    RooRealVar *BDT = wspace->var("BDT");
    RooRealVar *Tau = wspace->var("Tau");
    RooRealVar *TauRes = wspace->var("TauRes");
    RooCategory *SelCat = wspace->cat("SelCat");
    RooCategory *PairCat = wspace->cat("PairCat");
    RooCategory *GlobalCat = wspace->cat("GlobalCat");
    RooArgSet varlist(*Mass,*ReducedMassRes,*BDT,*Tau,*TauRes,*SelCat,*PairCat,*GlobalCat);
    
    cout << ">>> Load data from tree: " << filename << endl;
    
    exist_protection(filename);
    TFile *fin = new TFile(filename);
    TTree *tin = (TTree*)fin->Get("bmmData");
    
    double m_t, me_t, bdt_t, m1phi_t, m2phi_t, tau_t, taue_t;
    int chan_t, m1q_t;
    bool muid_t;
    tin->SetBranchAddress("m", &m_t);
    tin->SetBranchAddress("me", &me_t);
    tin->SetBranchAddress("bdt", &bdt_t);
    tin->SetBranchAddress("tau",  &tau_t);
    tin->SetBranchAddress("taue",  &taue_t);
    tin->SetBranchAddress("chan",  &chan_t);
    tin->SetBranchAddress("muid",  &muid_t);
    tin->SetBranchAddress("m1phi",  &m1phi_t);
    tin->SetBranchAddress("m2phi",  &m2phi_t);
    tin->SetBranchAddress("m1q",  &m1q_t);
    
    cout << ">>> parsing " << tin->GetEntries() << " entries." << endl;
    
    for (int evt=0; evt<tin->GetEntries();evt++) {
        tin->GetEntry(evt);
        
        if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
        if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
        if (bdt_t < -1.) continue;
        
        int index = CatMan.index(era, chan_t, bdt_t);
        if (index<0) continue;
        
        if (tau_t*1E12 > Tau_bound[0] && tau_t*1E12 < Tau_bound[1] &&
            taue_t*1E12 > TauRes_bound[0] && taue_t*1E12 < TauRes_bound[1]) SelCat->setIndex(1);
        else SelCat->setIndex(0);
        
        double dPhi = m1phi_t-m2phi_t;
        while (dPhi >= M_PI) dPhi -= M_PI*2;
        while (dPhi < -M_PI) dPhi += M_PI*2;
        bool isCowboy = (m1q_t*dPhi > 0);
        
        Mass->setVal(m_t);
        ReducedMassRes->setVal(me_t/m_t);
        BDT->setVal(bdt_t);
        Tau->setVal(tau_t*1E12);
        TauRes->setVal(taue_t*1E12);
        PairCat->setIndex(isCowboy);
        GlobalCat->setIndex(index);
        
        dataset->add(varlist);
    }
    fin->Close();
}

void PrepareData(RooWorkspace* wspace, RooWorkspace* wspace_base = 0)
{
    cout << ">>> PrepareData() start" << endl;
    if (wspace_base!=0) {
        cout << ">>> import from existing workspace." << endl;
        wspace->import(*wspace_base->data("global_data"));
        return;
    }

    cout << ">>> load events from trees." << endl;
    
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *ReducedMassRes = wspace->var("ReducedMassRes");
    RooRealVar *BDT = wspace->var("BDT");
    RooRealVar *Tau = wspace->var("Tau");
    RooRealVar *TauRes = wspace->var("TauRes");
    RooCategory *SelCat = wspace->cat("SelCat");
    RooCategory *PairCat = wspace->cat("PairCat");
    RooCategory *GlobalCat = wspace->cat("GlobalCat");
    RooArgSet varlist(*Mass,*ReducedMassRes,*BDT,*Tau,*TauRes,*SelCat,*PairCat,*GlobalCat);
    
    RooAbsData *global_data = new RooDataSet("global_data", "", varlist);
/*    if (CONFIG_BMM3) {
        LoadDataFromBMM3Tree(wspace, global_data, "input/2011/small-SgData-unblinded.root", "2011");
        LoadDataFromBMM3Tree(wspace, global_data, "input/2012/small-SgData-unblinded.root", "2012");
    }*/
    if (CONFIG_BMM4) {
        LoadDataFromBMM4Tree(wspace, global_data, "input/bmm4/small2016BFs01-bmmData.root", "2016BFs01");
        LoadDataFromBMM4Tree(wspace, global_data, "input/bmm4/small2016GHs01-bmmData.root", "2016GHs01");
    }

    wspace->import(*global_data);
}

// ----------------------------------------------------
// Prepare the resolution functions for lifetime fit
// current model: double-Gaussian, with <taue> as the main scaling factor; error on <taue> has been set to 12%
// available options:
// bsmm      - build B->mumu resolution model
// bupsik    - build B->J/psi K+ resolution model
// --
// triple    - using a triple Guassian model
// double    - reduced to a double Guassian model
// --
// available target_cat (other than nominal category IDs):
// mix       - using the expected yield to combine MC samples
//
void PrepareLifetimeResolutionModel(RooWorkspace *wspace, TString opt = "bsmm:triple", TString target_cat = "mix", RooWorkspace* wspace_base = 0)
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
    tag += TString("_")+target_cat;
    
    if (wspace_base!=0) {
        cout << ">>> import from existing workspace." << endl;
        RooResolutionModel *TauRes_Model = (RooResolutionModel*)wspace_base->obj(Form("TauRes_Model_%s",tag.Data()));
        wspace->import(*TauRes_Model);
        return;
    }
    
    RooRealVar delTau("delTau","",-0.6,+0.6);
    TH1D *h_global_deltau = new TH1D("h_global_deltau","",240,-0.6,+0.6);
    TH1D *h_global_taue = new TH1D("h_global_taue","",192,TauRes_bound[0],TauRes_bound[1]);
    h_global_deltau->Sumw2();
    h_global_taue->Sumw2();
    
    for (auto& cat: CatMan.cats) {
        if (target_cat!="mix" && cat.id!=target_cat) continue;
        
        // obtain the expected yield
        double yield = 0.;
        RooRealVar *N_bu = wspace->var(Form("N_bu_%s", cat.id.Data()));
        RooRealVar *effratio_bs = wspace->var(Form("effratio_bs_%s", cat.id.Data()));
        if (tag.Contains("bsmm")) yield = N_bu->getVal()*effratio_bs->getVal();
        if (tag.Contains("bupsik")) yield = N_bu->getVal();
        
        TH1D *h_deltau = new TH1D("h_deltau","",240,-0.6,+0.6);
        TH1D *h_taue = new TH1D("h_taue","",192,TauRes_bound[0],TauRes_bound[1]);
        h_deltau->Sumw2();
        h_taue->Sumw2();
        
        // load the corresponding MC sample
        TString filename = Form("input/bmm4/small%s-%s.root",cat.era.Data(),treename.Data());
        TChain *events = new TChain(treename);
        exist_protection(filename);
        events->Add(filename);
        cout << ">>> Loading from " << filename << ", with " << events->GetEntries() << " entries." << endl;
    
        double m_t, me_t, tau_t, taue_t, gtau_t, bdt_t;
        int chan_t;
        bool muid_t;
    
        events->SetBranchAddress("m", &m_t);
        events->SetBranchAddress("me", &me_t);
        events->SetBranchAddress("tau", &tau_t);
        events->SetBranchAddress("taue", &taue_t);
        events->SetBranchAddress("gtau", &gtau_t);
        events->SetBranchAddress("bdt", &bdt_t);
        events->SetBranchAddress("chan", &chan_t);
        events->SetBranchAddress("muid", &muid_t);
        
        for (int evt=0; evt<events->GetEntries();evt++) {
            events->GetEntry(evt);
            
            if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
            if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
            if (bdt_t < -1.) continue;
            
            int index = CatMan.index(cat.era, chan_t, bdt_t);
            if (index!=cat.index) continue;
            
            if (tau_t*1E12 < Tau_bound[0] || tau_t*1E12 > Tau_bound[1]) continue;
            if (taue_t*1E12 < TauRes_bound[0] || taue_t*1E12 > TauRes_bound[1]) continue;
            
            h_deltau->Fill((tau_t-gtau_t)*1E12);
            h_taue->Fill(taue_t*1E12);
        }
        
        cout << ">>> Category ID: " << cat.id << endl;
        
        if (target_cat=="mix") {
            cout << ">>> " << h_deltau->GetEntries() << " entries filled, to be scaled to " << yield << endl;
            double weight = yield/(double)h_deltau->GetEntries();
            h_global_deltau->Add(h_deltau,weight);
            h_global_taue->Add(h_taue,weight);
        }else {
            cout << ">>> " << h_deltau->GetEntries() << " entries filled/appended." << endl;
            h_global_deltau->Add(h_deltau);
            h_global_taue->Add(h_taue);
        }
        
        delete h_deltau;
        delete h_taue;
        delete events;
    }
    
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
    
    // set the taue error to 12%, which is the maximum difference between J/psi K+ data and MC in difference era
    par_taue.setError(mean_taue*0.12);
    
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
// bsmm         - build B->mumu efficiency model
// bupsik       - build B->J/psi K+ efficiency model
// --
// one_over_exp - modeling by [0]+[1]*x+[2]*x*x+[3]/(1.+exp(-[4]*x))
// threshold    - modeling by threshold function
// --
// available target_cat (other than nominal category IDs):
// mix       - using the expected yield to combine MC samples
//
void PrepareLifetimeEfficiencyModel(RooWorkspace *wspace, TString opt = "bsmm:threshold", TString target_cat = "mix", RooWorkspace* wspace_base = 0)
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
    tag += TString("_")+target_cat;
    
    if (wspace_base!=0) {
        cout << ">>> import from existing workspace." << endl;
        
        RooFormulaVar *TauEff_Model = (RooFormulaVar*)wspace_base->obj(Form("TauEff_Model_%s",tag.Data()));
        wspace->import(*TauEff_Model);

        RooHistPdf *TauEff_Model_Hist = (RooHistPdf*)wspace_base->obj(Form("TauEff_Model_Hist_%s",tag.Data()));
        wspace->import(*TauEff_Model_Hist);
        return;
    }

    vector<double> xbins = {
        0.5,0.625,0.75,0.875,
        1.,1.125,1.25,1.375,1.5,1.625,1.75,1.875,
        2.,2.125,2.25,2.375,2.5,2.75,
        3.,3.25,3.5,3.75,
        4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.};
    
    TH1D *h_global_taueff = new TH1D("h_global_taueff","",xbins.size()-1,xbins.data());
    h_global_taueff->Sumw2();
    
    for (auto& cat: CatMan.cats) {
        if (target_cat!="mix" && cat.id!=target_cat) continue;
        
        // obtain the expected yield
        double yield = 0.;
        RooRealVar *N_bu = wspace->var(Form("N_bu_%s", cat.id.Data()));
        RooRealVar *effratio_bs = wspace->var(Form("effratio_bs_%s", cat.id.Data()));
        if (tag.Contains("bsmm")) yield = N_bu->getVal()*effratio_bs->getVal();
        if (tag.Contains("bupsik")) yield = N_bu->getVal();
        
        TH1D *h_taueff = new TH1D("h_taueff","",xbins.size()-1,xbins.data());
        h_taueff->Sumw2();
    
        // load the corresponding MC sample
        TString filename = Form("input/bmm4/small%s-%s.root",cat.era.Data(),treename.Data());
        TChain *events = new TChain(treename);
        exist_protection(filename);
        events->Add(filename);
        cout << ">>> Loading from " << filename << ", with " << events->GetEntries() << " entries." << endl;
        
        double m_t, me_t, tau_t, taue_t, gtau_t, bdt_t;
        int chan_t;
        bool muid_t;
    
        events->SetBranchAddress("m", &m_t);
        events->SetBranchAddress("me", &me_t);
        events->SetBranchAddress("tau", &tau_t);
        events->SetBranchAddress("taue", &taue_t);
        events->SetBranchAddress("gtau", &gtau_t);
        events->SetBranchAddress("bdt", &bdt_t);
        events->SetBranchAddress("chan", &chan_t);
        events->SetBranchAddress("muid", &muid_t);
    
        for (int evt=0; evt<events->GetEntries();evt++) {
            events->GetEntry(evt);
        
            if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
            if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
            if (bdt_t < -1.) continue;
        
            int index = CatMan.index(cat.era, chan_t, bdt_t);
            if (index!=cat.index) continue;
        
            //if (tau_t*1E12 < Tau_bound[0] || tau_t*1E12 > Tau_bound[1]) continue;
            if (tau_t*1E12 < 0.5 || tau_t*1E12 > Tau_bound[1]) continue; // allow a litte bit more to the lower side
            if (taue_t*1E12 < TauRes_bound[0] || taue_t*1E12 > TauRes_bound[1]) continue;
        
            h_taueff->Fill(gtau_t*1E12);
        }

        cout << ">>> Category ID: " << cat.id << endl;
        
        if (target_cat=="mix") {
            cout << ">>> " << h_taueff->GetEntries() << " entries filled, to be scaled to " << yield << endl;
            double weight = yield/(double)h_taueff->GetEntries();
            h_global_taueff->Add(h_taueff, weight);
        }else {
            cout << ">>> " << h_taueff->GetEntries() << " entries filled/appended." << endl;
            h_global_taueff->Add(h_taueff);
        }
        
        delete events;
        delete h_taueff;
    }
    
    TH1D *h_taueff_norm = new TH1D("h_taueff_norm","",xbins.size()-1,xbins.data());
    h_taueff_norm->Sumw2();
    
    // extract lifetime & resolution
    RooRealVar *Tau = wspace->var("Tau");
    //RooResolutionModel *TauRes_Model = (RooResolutionModel*)wspace->obj(Form("TauRes_Model_%s",tag.Data()));
    RooTruthModel *TauRes_Model = new RooTruthModel("TauRes_Model","",*Tau); // ideal model
    RooRealVar EffTau("EffTau","",1.6);
    RooDecay RawDecay("RawDecay","",*Tau,EffTau,*TauRes_Model,RooDecay::SingleSided);
    
    if (tag.Contains("bsmm")) EffTau.setVal(1.472);
    if (tag.Contains("bupsik") && (target_cat.Contains("2011") || target_cat.Contains("2012"))) EffTau.setVal(1.671);
    if (tag.Contains("bupsik") && target_cat.Contains("2016"))  EffTau.setVal(1.638);
    if (tag.Contains("bupsik") && target_cat=="mix") {
        // B+ gen lifetime was 1.671 ps in 2011/2012 MC, was 1.638 ps in 2016 MC
        // Have to estimate the expected gen distribution (roughly):
        // 2011: Lumi = 5/fb x 7 TeV 1.1730e+07 pb
        // 2012: Lumi = 20/fb x 8 TeV 1.3660e+07 pb
        // 2016: Lumi = 36/fb x 13 TeV 2.3120e+07 pb
        // 2011+2012 vs 2016: 1:2.51
        cout << ">>> 'mix' option for J/psi K+ is not available; exit." << endl;
        exit(1);
    }
    
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
    TH1D *h_tau_tmp = new TH1D("h_tau_tmp","",(int)((Tau_bound[1]-Tau_bound[0])*8),Tau_bound[0],Tau_bound[1]); // convert to fixed bin width histogram
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
    
    for (int i=0; i<5; i++) delete par[i];
    delete h_global_taueff;
    delete effmodel;
    delete leg1;
    delete canvas;
    delete frame;
}

void PrepareLifetimeModel(RooWorkspace *wspace, RooWorkspace* wspace_base = 0)
{
    PrepareLifetimeResolutionModel(wspace, "bsmm:triple", "mix", wspace_base);
    PrepareLifetimeEfficiencyModel(wspace, "bsmm:threshold", "mix", wspace_base);
    
    for (auto& cat: CatMan.cats) {
        PrepareLifetimeResolutionModel(wspace, "bsmm:double", cat.id, wspace_base);
        PrepareLifetimeEfficiencyModel(wspace, "bsmm:threshold", cat.id, wspace_base);
    }
}

void BuildReducedMassResPDF(RooWorkspace* wspace, TString key, RooDataSet *rds)
{
    RooRealVar* ReducedMassRes = wspace->var("ReducedMassRes");
    
    RooDataSet *rds_keys = NULL;
    if (rds->numEntries() > 20000) {
        TH1* hist = rds->createHistogram("ReducedMassRes",1000);
        RooRealVar* Weight = wspace->var("Weight");
        RooArgSet varlist(*ReducedMassRes,*Weight);
        rds_keys = new RooDataSet(Form("ReducedMassRes_rds_%s",key.Data()), "", varlist, "Weight");

        for (int idx=1; idx<=hist->GetNbinsX(); idx++) {
            int count = (int)ceil(20000.*hist->GetBinContent(idx)/hist->GetSumOfWeights());
            double delta = hist->GetBinWidth(idx)/count;
            
            for (int n=0; n<count; n++) {
                ReducedMassRes->setVal(hist->GetBinLowEdge(idx)+delta*n+0.5*delta);
                Weight->setVal(hist->GetBinContent(idx)/count);
                rds_keys->add(varlist, Weight->getVal());
            }
        }
        delete hist;
    }else rds_keys = (RooDataSet*)rds->Clone(Form("ReducedMassRes_rds_%s",key.Data()));
    RooKeysPdf model(Form("ReducedMassRes_pdf_%s",key.Data()), "", *ReducedMassRes, *rds_keys);
    
    wspace->import(model);
}

void BuildBDTPDF(RooWorkspace* wspace, TString key, RooDataSet *rds)
{
    RooRealVar* BDT = wspace->var("BDT");
    
    RooDataSet *rds_keys = NULL;
    if (rds->numEntries() > 20000) {
        TH1* hist = rds->createHistogram("BDT",1000);
        RooRealVar* Weight = wspace->var("Weight");
        RooArgSet varlist(*BDT,*Weight);
        rds_keys = new RooDataSet(Form("BDT_rds_%s",key.Data()), "", varlist, "Weight");
        
        for (int idx=1; idx<=hist->GetNbinsX(); idx++) {
            int count = (int)ceil(20000.*hist->GetBinContent(idx)/hist->GetSumOfWeights());
            double delta = hist->GetBinWidth(idx)/count;
            
            for (int n=0; n<count; n++) {
                BDT->setVal(hist->GetBinLowEdge(idx)+delta*n+0.5*delta);
                Weight->setVal(hist->GetBinContent(idx)/count);
                rds_keys->add(varlist, Weight->getVal());
            }
        }
        delete hist;
    }else rds_keys = (RooDataSet*)rds->Clone(Form("BDT_rds_%s",key.Data()));
    RooKeysPdf model(Form("BDT_pdf_%s",key.Data()), "", *BDT, *rds_keys, RooKeysPdf::MirrorBoth);
    
    wspace->import(model);
}

// ----------------------------------------------------
// available options:
// single       - signal exponetial decay
// double       - double exponetial decay
// --
// hist         - take the histogram efficiency function
//
void BuildDecayTimePDF(RooWorkspace* wspace, TString key, TString target_cat, RooDataSet *rds, TString opt = "single", RooRealVar *EffTauCommon = 0)
{
    RooRealVar* Tau = wspace->var("Tau");
    
    RooResolutionModel *TauRes_Model = (RooResolutionModel*)wspace->obj(Form("TauRes_Model_bsmm_%s",target_cat.Data()));
    RooFormulaVar *TauEff_Model = (RooFormulaVar*)wspace->obj(Form("TauEff_Model_bsmm_%s",target_cat.Data()));
    
    if (opt.Contains("hist")) { // histogram efficiency
        TauEff_Model = (RooFormulaVar*)wspace->obj(Form("TauEff_Model_Hist_bsmm_%s",target_cat.Data()));
    }
    
    if (opt.Contains("single")) { // single Exp
        RooRealVar EffTau(Form("EffTau_%s",key.Data()),"",1.6,0.1,5.0);
        RooDecay RawDecay(Form("RawDecay_%s",key.Data()),"",*Tau,(EffTauCommon!=0?*EffTauCommon:EffTau),*TauRes_Model,RooDecay::SingleSided);
        RooEffProd Tau_pdf(Form("Tau_pdf_%s",key.Data()),"",RawDecay,*TauEff_Model);

        if (EffTauCommon==0) {
            RooFitResult *res = Tau_pdf.fitTo(*rds, Extended(false), SumW2Error(true), NumCPU(NCPU), Hesse(false), Save(true));
            if (res->status()!=0) converge_protection();
            delete res;
        }
    
        EffTau.setConstant(true);
        wspace->import(Tau_pdf);
    }else if (opt.Contains("double")) { // double Exp
        RooRealVar EffTau1(Form("EffTau1_%s",key.Data()),"",1.0,0.1,5.0);
        RooRealVar EffTau2(Form("EffTau2_%s",key.Data()),"",1.6,0.1,15.0);
        RooRealVar EffTau_frac(Form("EffTau_frac_%s",key.Data()),"",0.5,0.,1.);
        RooDecay RawDecay1(Form("RawDecay1_%s",key.Data()),"",*Tau,EffTau1,*TauRes_Model,RooDecay::SingleSided);
        RooDecay RawDecay2(Form("RawDecay2_%s",key.Data()),"",*Tau,EffTau2,*TauRes_Model,RooDecay::SingleSided);
        RooAddPdf RawDecay(Form("RawDecay_%s",key.Data()), "", RooArgList(RawDecay1, RawDecay2),  EffTau_frac);
        RooEffProd Tau_pdf(Form("Tau_pdf_%s",key.Data()),"",RawDecay,*TauEff_Model);
        
        RooFitResult *res = Tau_pdf.fitTo(*rds, Extended(false), SumW2Error(true), NumCPU(NCPU), Hesse(false), Save(true));
        if (res->status()!=0) converge_protection();
        delete res;
        
        EffTau1.setConstant(true);
        EffTau2.setConstant(true);
        EffTau_frac.setConstant(true);
        wspace->import(Tau_pdf);
    }
}

// ----------------------------------------------------
// available options:
// flat   - uniform along seagull/cowboy
// mu     - seagull = 1 / cowboy = f
// musq   - seagull = 1 / cowboy = f^2
//
void BuildPairCategoryPDF(RooWorkspace* wspace, TString key, TString opt = "flat")
{
    RooCategory* PairCat = wspace->cat("PairCat");
    RooRealVar* dblmu_corr_scale = wspace->var("dblmu_corr_scale");
    
    if (opt=="flat") {
        RooGenericPdf PairCat_pdf(Form("PairCat_pdf_%s",key.Data()), "", "1.", RooArgList(*PairCat));
        wspace->import(PairCat_pdf);
    }else if (opt=="mu") {
        RooGenericPdf PairCat_pdf(Form("PairCat_pdf_%s",key.Data()), "", "@0?@1:1.", RooArgList(*PairCat,*dblmu_corr_scale));
        wspace->import(PairCat_pdf);
    }else if (opt=="musq") {
        RooGenericPdf PairCat_pdf(Form("PairCat_pdf_%s",key.Data()), "", "@0?@1*@1:1.", RooArgList(*PairCat,*dblmu_corr_scale));
        wspace->import(PairCat_pdf);
    }
}

// choice of flavour
enum{_bs_signal, _bd_signal};

void BuildSignalMassPDF(RooWorkspace* wspace, TString key, int flavour, RooDataSet *rds)
{
    //RooRealVar* Mass = wspace->var("Mass");
    RooRealVar* ReducedMassRes = wspace->var("ReducedMassRes");
    RooAbsPdf* ReducedMassRes_pdf = wspace->pdf(Form("ReducedMassRes_pdf_%s",key.Data()));
    
    TString id = key(key.First('_')+1,key.Length());
    RooRealVar* DeltaMass = wspace->var(Form("DeltaMass_%s",id.Data()));
    RooFormulaVar* ScaledMass = (RooFormulaVar*)wspace->function(Form("ScaledMass_%s",id.Data()));
    double DeltaMass_reserve = DeltaMass->getVal();
    DeltaMass->setVal(0.);
    
    double mean_init = 0.;
    if (flavour==_bs_signal) mean_init = 5.35;
    if (flavour==_bd_signal) mean_init = 5.25;
    
    RooRealVar Mean(Form("Mean_%s",key.Data()), "", mean_init, mean_init-0.05, mean_init+0.05);
    RooRealVar Alpha(Form("Alpha_%s",key.Data()), "", 2.8, 0.1, 3.0);
    RooRealVar Enne(Form("Enne_%s",key.Data()), "", 1., 0., 10.);
    RooRealVar PeeK(Form("PeeK_%s",key.Data()), "", 1., 0.1, 10.);
    
    RooFormulaVar SigmaRes(Form("SigmaRes_%s",key.Data()), "@0*@1*@2", RooArgList(*ReducedMassRes, *ScaledMass, PeeK));
    
    RooCBShape CB(Form("CB_%s",key.Data()), "", *ScaledMass, Mean, SigmaRes, Alpha, Enne);
    
    RooProdPdf Mass_pdf(Form("Mass_pdf_%s",key.Data()), "", *ReducedMassRes_pdf, Conditional(CB, *ScaledMass));
    
    RooFitResult *res = Mass_pdf.fitTo(*rds, ConditionalObservables(*ReducedMassRes), Extended(false), SumW2Error(true), NumCPU(NCPU), Hesse(false), Save(true));
    if (res->status()!=0) converge_protection();
    delete res;
    
    Mean.setConstant(true);
    Alpha.setConstant(true);
    Enne.setConstant(true);
    PeeK.setConstant(true);
    
    DeltaMass->setVal(DeltaMass_reserve);
    
    wspace->import(Mass_pdf);
}

void BuildSemiMassPDF(RooWorkspace* wspace, TString key, RooDataSet *rds)
{
    RooRealVar* Mass = wspace->var("Mass");
    RooAbsPdf* ReducedMassRes_pdf = wspace->pdf(Form("ReducedMassRes_pdf_%s",key.Data()));
    
    RooDataSet *rds_keys = NULL;
    if (rds->numEntries() > 20000) {
        TH1* hist = rds->createHistogram("Mass",1000);
        RooRealVar* Weight = wspace->var("Weight");
        RooArgSet varlist(*Mass,*Weight);
        rds_keys = new RooDataSet(Form("Mass_rds_%s",key.Data()), "", varlist, "Weight");
        
        for (int idx=1; idx<=hist->GetNbinsX(); idx++) {
            int count = (int)ceil(20000.*hist->GetBinContent(idx)/hist->GetSumOfWeights());
            double delta = hist->GetBinWidth(idx)/count;
            
            for (int n=0; n<count; n++) {
                Mass->setVal(hist->GetBinLowEdge(idx)+delta*n+0.5*delta);
                Weight->setVal(hist->GetBinContent(idx)/count);
                rds_keys->add(varlist, Weight->getVal());
            }
        }
        delete hist;
    }else rds_keys = (RooDataSet*)rds->Clone(Form("Mass_rds_%s",key.Data()));
    RooKeysPdf Keys(Form("Keys_%s",key.Data()), "", *Mass, *rds_keys, RooKeysPdf::MirrorBoth, 2.5);
    
    RooProdPdf Mass_pdf(Form("Mass_pdf_%s",key.Data()), "", *ReducedMassRes_pdf, Conditional(Keys, *Mass));
    
    wspace->import(Mass_pdf);
}

void BuildH2muMassPDF(RooWorkspace* wspace, TString key, RooDataSet *rds)
{
    RooRealVar* Mass = wspace->var("Mass");
    RooAbsPdf* ReducedMassRes_pdf = wspace->pdf(Form("ReducedMassRes_pdf_%s",key.Data()));
    
    RooDataSet *rds_keys = NULL;
    if (rds->numEntries() > 20000) {
        TH1* hist = rds->createHistogram("Mass",1000);
        RooRealVar* Weight = wspace->var("Weight");
        RooArgSet varlist(*Mass,*Weight);
        rds_keys = new RooDataSet(Form("Mass_rds_%s",key.Data()), "", varlist, "Weight");
        
        for (int idx=1; idx<=hist->GetNbinsX(); idx++) {
            int count = (int)ceil(20000.*hist->GetBinContent(idx)/hist->GetSumOfWeights());
            double delta = hist->GetBinWidth(idx)/count;
            
            for (int n=0; n<count; n++) {
                Mass->setVal(hist->GetBinLowEdge(idx)+delta*n+0.5*delta);
                Weight->setVal(hist->GetBinContent(idx)/count);
                rds_keys->add(varlist, Weight->getVal());
            }
        }
        delete hist;
    }else rds_keys = (RooDataSet*)rds->Clone(Form("Mass_rds_%s",key.Data()));
    RooKeysPdf Keys(Form("Keys_%s",key.Data()), "", *Mass, *rds_keys, RooKeysPdf::MirrorBoth, 2.5);
    
    RooProdPdf Mass_pdf(Form("Mass_pdf_%s",key.Data()), "", *ReducedMassRes_pdf, Conditional(Keys, *Mass));
    
    wspace->import(Mass_pdf);
}

void BuildPeakMassPDF(RooWorkspace* wspace, TString key, RooDataSet *rds)
{
    //RooRealVar* Mass = wspace->var("Mass");
    RooRealVar* ReducedMassRes = wspace->var("ReducedMassRes");
    RooAbsPdf* ReducedMassRes_pdf = wspace->pdf(Form("ReducedMassRes_pdf_%s",key.Data()));
    
    TString id = key(key.First('_')+1,key.Length());
    RooRealVar* DeltaMass = wspace->var(Form("DeltaMass_%s",id.Data()));
    RooFormulaVar* ScaledMass = (RooFormulaVar*)wspace->function(Form("ScaledMass_%s",id.Data()));
    double DeltaMass_reserve = DeltaMass->getVal();
    DeltaMass->setVal(0.);
    
    RooRealVar Mean(Form("Mean_%s",key.Data()), "", 5.1, 4.9, 5.4);
    RooRealVar Sigma(Form("Sigma_%s",key.Data()), "", 0.02, 0.005, 0.2);
    RooRealVar Sigma2(Form("Sigma2_%s",key.Data()), "", 0.04, 0.005, 0.2);
    RooRealVar Alpha(Form("Alpha_%s",key.Data()), "", 2.8, 0., 100.0);
    RooRealVar Enne(Form("Enne_%s",key.Data()), "", 1., 0., 10.);
    RooRealVar CoeffGauss(Form("CoeffGauss_%s",key.Data()), "", 0.5, 0., 1.);
    
    RooGaussian Gauss(Form("Gauss_%s",key.Data()), "", *ScaledMass, Mean, Sigma);
    
    // add a prefit to stablize the fit
    Gauss.fitTo(*rds, Extended(false), SumW2Error(true), NumCPU(NCPU), Hesse(false));
    
    RooCBShape CB(Form("CB_%s",key.Data()), "", *ScaledMass, Mean, Sigma2, Alpha, Enne);
    
    RooAddPdf CBPlusGau(Form("CBPlusGau_%s",key.Data()), "", RooArgList(Gauss, CB),  CoeffGauss);

    RooProdPdf Mass_pdf(Form("Mass_pdf_%s",key.Data()), "", *ReducedMassRes_pdf, Conditional(CBPlusGau, *ScaledMass));
    
    RooFitResult *res = Mass_pdf.fitTo(*rds, ConditionalObservables(*ReducedMassRes), Extended(false), SumW2Error(true), NumCPU(NCPU), Hesse(false), Save(true));
    if (res->status()!=0) converge_protection();
    delete res;
    
    Mean.setConstant(true);
    Sigma.setConstant(true);
    Sigma2.setConstant(true);
    Alpha.setConstant(true);
    Enne.setConstant(true);
    CoeffGauss.setConstant(true);
    
    DeltaMass->setVal(DeltaMass_reserve);
    
    wspace->import(Mass_pdf);
}

void BuildCombMassPDF(RooWorkspace* wspace, TString key)
{
    RooRealVar* Mass = wspace->var("Mass");
    RooAbsPdf* ReducedMassRes_pdf = wspace->pdf(Form("ReducedMassRes_pdf_%s",key.Data()));
    
    RooRealVar Br1(Form("Br1_%s",key.Data()), "", 0.5, 0. , 1);
    
    RooFormulaVar Br2(Form("Br2_%s",key.Data()), "", "1.-@0", RooArgList(Br1));
    RooBernstein Bern(Form("Bern_%s",key.Data()), "", *Mass, RooArgList(Br1, Br2));
    
    RooProdPdf Mass_pdf(Form("Mass_pdf_%s",key.Data()), "", *ReducedMassRes_pdf, Conditional(Bern, *Mass));
    
    wspace->import(Mass_pdf);
}

void ProducePDFProjectionPlot(RooWorkspace* wspace, TString key, RooDataSet *rds, TString opt = "me:bmm4:dblbins")
{
    RooRealVar* POI = NULL;
    RooAbsPdf *pdf = NULL;
    RooCategory *SpCat = wspace->cat("SpCat");
    TString var, title, filename;
    int nbins = 40;
    
    if (opt.Contains("dblbins")) nbins *= 2;
    if (opt.Contains("me")) {
        var = "ReducedMassRes";
        POI = wspace->var(var);
        pdf = wspace->pdf(Form("ReducedMassRes_pdf_%s",key.Data()));
        title = Form("ReducedMassRes_pdf_%s",key.Data());
        filename = Form("fig/proj_ReducedMassRes_pdf_%s.pdf",key.Data());
    } else if (opt.Contains("bdt")) {
        var = "BDT";
        POI = wspace->var(var);
        pdf = wspace->pdf(Form("BDT_pdf_%s",key.Data()));
        title = Form("BDT_pdf_%s",key.Data());
        filename = Form("fig/proj_BDT_pdf_%s.pdf",key.Data());
    } else if (opt.Contains("tau")) {
        var = "Tau";
        POI = wspace->var(var);
        pdf = wspace->pdf(Form("Tau_pdf_%s",key.Data()));
        title = Form("Tau_pdf_%s",key.Data());
        filename = Form("fig/proj_Tau_pdf_%s.pdf",key.Data());
    } else if (opt.Contains("m")) {
        var = "Mass";
        POI = wspace->var(var);
        pdf = wspace->pdf(Form("Mass_pdf_%s",key.Data()));
        title = Form("Mass_pdf_%s",key.Data());
        filename = Form("fig/proj_Mass_pdf_%s.pdf",key.Data());
    }

    RooPlot* frame = POI->frame(Bins(nbins), Title(" "));
    rds->plotOn(frame, Invisible());
    pdf->plotOn(frame, DrawOption("L"), LineColor(kBlack), LineWidth(3), LineStyle(1), NumCPU(NCPU), Name("pdf"));
    
    THStack *hs = new THStack("histstack","");
    TH1D* hist_sub[SpCat->numTypes()];
    
    for (int i=0; i<SpCat->numTypes(); i++)
        hist_sub[i] = new TH1D(Form("hist_sub%d",i),"",nbins,POI->getMin(),POI->getMax());
    for (int i=0; i<rds->numEntries(); i++) {
        const RooArgSet* arg = rds->get(i);
        hist_sub[arg->getCatIndex("SpCat")]->Fill(arg->getRealValue(var),rds->weight());
    }
    for (int i=0; i<SpCat->numTypes(); i++) {
        hist_sub[i]->SetFillColor(i+20);
        hist_sub[i]->SetLineColor(kBlack);
        hist_sub[i]->SetLineWidth(0);
        hs->Add(hist_sub[i],"hist");
    }
    
    frame->SetMinimum(0.);
    frame->SetMaximum(hs->GetMaximum()*1.35);
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.14,0.06,0.13,0.07);
    
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetYaxis()->SetTitle("Entries");
    
    if (var == "ReducedMassRes") frame->GetXaxis()->SetTitle("#sigma_{M(#mu#mu)}/M(#mu#mu)");
    else if (var == "Mass") frame->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    else if (var == "BDT") frame->GetXaxis()->SetTitle("BDT");
    else if (var == "Tau") frame->GetXaxis()->SetTitle("Decay time [ps]");

    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->Draw("");
    hs->Draw("same");
    frame->Draw("same");
    
    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.035);
    tex.SetTextAlign(11);
    tex.SetNDC();
    tex.DrawLatex(0.14,0.94,title);

    int n_entry = 0;
    for (int i=0; i<SpCat->numTypes(); i++)
        if (hist_sub[i]->GetEntries()>0) n_entry++;
    TLegend *leg1 = new TLegend(0.20,(n_entry<3?0.86:(n_entry<6?0.81:0.76)),0.91,0.91);
    leg1->SetNColumns(3);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(kWhite);
    leg1->AddEntry(frame->findObject("pdf"),"PDF","L");
    if (opt.Contains("bmm4")) {
        for (int i=0; i<bmm4::ndecays; i++) {
            if (hist_sub[i]->GetEntries()>0)
                leg1->AddEntry(hist_sub[i],bmm4::fullprocess[i],"f");
        }
    }else if (opt.Contains("bmm3")) {
        for (int i=0; i<bmm3::ndecays; i++) {
            if (hist_sub[i]->GetEntries()>0)
                leg1->AddEntry(hist_sub[i],bmm3::fullprocess[i],"f");
        }
    }
    leg1->Draw();
    
    canvas->Update();
    canvas->Print(filename);
    
    delete hs;
    for (int i=0; i<SpCat->numTypes(); i++) delete hist_sub[i];
    delete frame;
    delete canvas;
}

void PrepareBMM3SubPDF(RooWorkspace* wspace, RooWorkspace* wspace_base = 0)
{
    cout << ">>> PrepareBMM3SubPDF() start" << endl;
    
    if (wspace_base!=0) {
        cout << ">>> import from existing workspace." << endl;
        for (auto& cat: CatMan.cats) {
            if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
            for (TString type : {"Mass","BDT","PairCat","Tau","SelCat"}) {
                wspace->import(*wspace_base->pdf(Form("%s_pdf_bs_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_bd_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_peak_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_semi_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_h2mu_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_comb_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
            }
        }
        return;
    }
    cout << ">>> build PDFs from trees & data cards." << endl;
    
    vector<RooDataSet*> samples;
    samples.assign(bmm3::ndecays,NULL);
    
    // create the RooDataSet
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *ReducedMassRes = wspace->var("ReducedMassRes");
    RooRealVar *BDT = wspace->var("BDT");
    RooRealVar *Tau = wspace->var("Tau");
    RooRealVar *TauRes = wspace->var("TauRes");
    RooCategory *GlobalCat = wspace->cat("GlobalCat");
    RooCategory *SpCat = wspace->cat("SpCat");
    RooCategory *SelCat = wspace->cat("SelCat");
    RooRealVar *Weight = wspace->var("Weight");
    RooArgSet varlist(*Mass,*ReducedMassRes,*BDT,*Tau,*TauRes,*GlobalCat,*SpCat,*SelCat,*Weight);
    
    for (int i=0; i<bmm3::ndecays; i++)
        samples[i] = new RooDataSet(Form("rds_%s",bmm3::decays[i].Data()), "", varlist, "Weight");
    RooDataSet *rds_data_lowbdt = new RooDataSet("rds_data_lowbdt", "", varlist, "Weight");
    RooDataSet *rds_data_himass = new RooDataSet("rds_data_himass", "", varlist, "Weight");
    
    // allocate counters for processing each categories
    vector<double> yields[bmm3::ndecays], counts[bmm3::ndecays];
    for (int idx_sample=0; idx_sample<bmm3::ndecays; idx_sample++) {
        yields[idx_sample].assign(CatMan.cats.size(),0.);
        counts[idx_sample].assign(CatMan.cats.size(),0.);
    }
    
    // read the expected yields from data cards
    for (int idx_sample=bmm3::_bgBs2KK; idx_sample<=bmm3::_bgBs2MuMuGamma; idx_sample++) {
        for (auto& cat: CatMan.cats) {
            if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
            
            vector<TString> keys;
            vector<double> values;
            
            keys.push_back(Form("%s:loSideband%d:val",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bdRare%d}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:bsRare%d}",bmm3::decays[idx_sample].Data(),cat.region));
            keys.push_back(Form("%s:hiSideband%d:val",bmm3::decays[idx_sample].Data(),cat.region));
            
            ReadValuesFromTex(Form("input/%s/anaBmm.plotResults.%s-cat2%d.tex",cat.era.Data(),cat.era.Data(),cat.bdt_bin),keys,values);
            
            yields[idx_sample][cat.index] = values[0]+values[1]+values[2]+values[3]; // yields
        }
    }
    
    // loop over the samples
    for (int idx_sample=0; idx_sample<bmm3::ndecays; idx_sample++) {
        for (TString era : bmm3::eras) {
            
            bool applyWeights = true;
            if (idx_sample==bmm3::_SgMc || idx_sample==bmm3::_BdMc || idx_sample==bmm3::_SgData) applyWeights = false;
            
            TString filename = Form("input/%s/small-%s.root",era.Data(),bmm3::decays[idx_sample].Data());
            TString treename = bmm3::decays[idx_sample]+TString("_bdt");
            
            exist_protection(filename);
            TFile *fin = new TFile(filename);
            TTree *tin = (TTree*)fin->Get(treename);
            
            cout << ">>> Loading from " << filename << ", with " << tin->GetEntries() << " entries." << endl;
            
            double m1eta_t, m2eta_t, m_t, bdt_t, me_t, tau_t;
            bool muid_t;
            tin->SetBranchAddress("m1eta", &m1eta_t);
            tin->SetBranchAddress("m2eta", &m2eta_t);
            tin->SetBranchAddress("m", &m_t);
            tin->SetBranchAddress("me", &me_t);
            tin->SetBranchAddress("bdt", &bdt_t);
            tin->SetBranchAddress("tau", &tau_t);
            tin->SetBranchAddress("muid",  &muid_t);
            
            // first loop: event counting, peaking/semi MC only
            int nevt_limit = 0;
            if (applyWeights)
                for (int evt=0; evt<tin->GetEntries();evt++) {
                    tin->GetEntry(evt);
                    
                    if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
                    if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
                    if (bdt_t < -1.) continue;
                    
                    bool isBarrel = fabs(m1eta_t)<1.4 && fabs(m2eta_t)<1.4;
                    int index = CatMan.index(era, isBarrel?0:1, bdt_t);
                    if (index<0) continue;
                    
                    counts[idx_sample][index] += 1.; // adding counts
                    
                    nevt_limit++;
                    if (nevt_limit>=MC_EVENT_LIMIT) break; // read up to MC_EVENT_LIMIT per sample
                }
            
            // second loop: event filling
            nevt_limit = 0;
            for (int evt=0; evt<tin->GetEntries();evt++) {
                tin->GetEntry(evt);
                
                if (m_t < Mass_bound[0] || m_t > Mass_bound[1]) continue;
                if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
                if (bdt_t < -1.) continue;
                
                bool isBarrel = fabs(m1eta_t)<1.4 && fabs(m2eta_t)<1.4;
                int index = CatMan.index(era, isBarrel?0:1, bdt_t);
                
                if (tau_t*1E12 > Tau_bound[0] && tau_t*1E12 < Tau_bound[1]) SelCat->setIndex(1);
                else SelCat->setIndex(0);
                
                Mass->setVal(m_t);
                ReducedMassRes->setVal(me_t/m_t);
                BDT->setVal(bdt_t);
                Tau->setVal(tau_t*1E12);
                TauRes->setVal(0.07); // ISSUE: TauRes is not available in BMM3
                GlobalCat->setIndex(index);
                Weight->setVal(1.);
                if (applyWeights)
                    Weight->setVal(yields[idx_sample][index]/counts[idx_sample][index]);
                SpCat->setIndex(idx_sample);
                
                if (idx_sample == bmm3::_SgData) {
                    
                    for (auto& cat: CatMan.cats) { // fill event regardless of BDT binning for low BDT sideband
                        if (cat.era == era && cat.region == (isBarrel?0:1)) {
                            if (!muid_t) continue; // reject failed muon ID events
                            if (m_t > 5.20 && m_t < 5.45) continue; // exclude the signal for data sideband
                            if (bdt_t < -0.2 || bdt_t > 0.1) continue; // exclude good BDT region
                            GlobalCat->setIndex(cat.index);
                            rds_data_lowbdt->add(varlist, Weight->getVal());
                        }
                    }
                    
                    for (auto& cat: CatMan.cats) { // fill event regardless of BDT binning for high mass sideband
                        if (cat.era == era && cat.region == (isBarrel?0:1)) {
                            if (m_t < 5.45) continue; // exclude the signal & low mass regions
                            if (bdt_t < 0.1) continue; // exclude low BDT region
                            GlobalCat->setIndex(cat.index);
                            rds_data_himass->add(varlist, Weight->getVal());
                        }
                    }
                    
                    nevt_limit++;
                }else {
                    if (!muid_t) continue; // reject failed muon ID events
                    if (index<0) continue; // reject non-categorized events
                    
                    samples[idx_sample]->add(varlist, Weight->getVal());
                    
                    nevt_limit++;
                    if (idx_sample!=bmm3::_SgData && nevt_limit>=MC_EVENT_LIMIT) break; // read up to MC_EVENT_LIMIT per sample
                }
            } // end of evt loop
            
            cout << ">>> " << nevt_limit << " events loaded." << endl;
            
            fin->Close();
            
        } // end of "era" loop
    } // end of "idx_sample" loop
    
    // alias & merge
    RooDataSet *rds_bs = samples[bmm3::_SgMc];
    RooDataSet *rds_bd = samples[bmm3::_BdMc];
    RooDataSet* rds_peak = (RooDataSet*)samples[bmm3::_bgBs2KK]->Clone("rds_peak");
    for (int i=bmm3::_bgBs2KK+1; i<=bmm3::_bgLb2KP; i++) rds_peak->append(*samples[i]);
    RooDataSet* rds_semi = (RooDataSet*)samples[bmm3::_bgBs2KMuNu]->Clone("rds_semi");
    for (int i=bmm3::_bgBs2KMuNu+1; i<=bmm3::_bgLb2PMuNu; i++) rds_semi->append(*samples[i]);
    RooDataSet* rds_h2mu = (RooDataSet*)samples[bmm3::_bgBu2PiMuMu]->Clone("rds_h2mu");
    for (int i=bmm3::_bgBu2PiMuMu+1; i<=bmm3::_bgBs2MuMuGamma; i++) rds_h2mu->append(*samples[i]);
    
    // build ReducedMassRes PDF
    for (auto& cat: CatMan.cats) {
        if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
        
        RooDataSet *rds_bs_res = (RooDataSet*)rds_bs->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_bd_res = (RooDataSet*)rds_bd->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_peak_res = (RooDataSet*)rds_peak->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_semi_res = (RooDataSet*)rds_semi->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_h2mu_res = (RooDataSet*)rds_h2mu->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_comb_res = (RooDataSet*)rds_data_lowbdt->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        
        BuildReducedMassResPDF(wspace, Form("bs_%s",cat.id.Data()), rds_bs_res);
        BuildReducedMassResPDF(wspace, Form("bd_%s",cat.id.Data()), rds_bd_res);
        BuildReducedMassResPDF(wspace, Form("peak_%s",cat.id.Data()), rds_peak_res);
        BuildReducedMassResPDF(wspace, Form("semi_%s",cat.id.Data()), rds_semi_res);
        BuildReducedMassResPDF(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_res);
        BuildReducedMassResPDF(wspace, Form("comb_%s",cat.id.Data()), rds_comb_res);
        
        ProducePDFProjectionPlot(wspace, Form("bs_%s",cat.id.Data()), rds_bs_res, "me:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("bd_%s",cat.id.Data()), rds_bd_res, "me:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("peak_%s",cat.id.Data()), rds_peak_res, "me:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("semi_%s",cat.id.Data()), rds_semi_res, "me:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_res, "me:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("comb_%s",cat.id.Data()), rds_comb_res, "me:bmm3:dblbins");
        
        delete rds_bs_res;
        delete rds_bd_res;
        delete rds_peak_res;
        delete rds_semi_res;
        delete rds_h2mu_res;
        delete rds_comb_res;
    }
    
    // build BDT PDF
    for (auto& cat: CatMan.cats) {
        if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
        
        RooDataSet *rds_bs_bdt = (RooDataSet*)rds_bs->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_bd_bdt = (RooDataSet*)rds_bd->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_peak_bdt = (RooDataSet*)rds_peak->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_semi_bdt = (RooDataSet*)rds_semi->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_h2mu_bdt = (RooDataSet*)rds_h2mu->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_comb_bdt = (RooDataSet*)rds_data_himass->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        
        BuildBDTPDF(wspace, Form("bs_%s",cat.id.Data()), rds_bs_bdt);
        BuildBDTPDF(wspace, Form("bd_%s",cat.id.Data()), rds_bd_bdt);
        BuildBDTPDF(wspace, Form("peak_%s",cat.id.Data()), rds_peak_bdt);
        BuildBDTPDF(wspace, Form("semi_%s",cat.id.Data()), rds_semi_bdt);
        BuildBDTPDF(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_bdt);
        BuildBDTPDF(wspace, Form("comb_%s",cat.id.Data()), rds_comb_bdt);
        
        ProducePDFProjectionPlot(wspace, Form("bs_%s",cat.id.Data()), rds_bs_bdt, "bdt:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("bd_%s",cat.id.Data()), rds_bd_bdt, "bdt:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("peak_%s",cat.id.Data()), rds_peak_bdt, "bdt:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("semi_%s",cat.id.Data()), rds_semi_bdt, "bdt:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_bdt, "bdt:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("comb_%s",cat.id.Data()), rds_comb_bdt, "bdt:bmm3:dblbins");
        
        delete rds_bs_bdt;
        delete rds_bd_bdt;
        delete rds_peak_bdt;
        delete rds_semi_bdt;
        delete rds_h2mu_bdt;
        delete rds_comb_bdt;
    }
    
    // build decay time PDFs
    for (auto& cat: CatMan.cats) {
        if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
        
        RooDataSet *rds_bs_tau = (RooDataSet*)rds_bs->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_bd_tau = (RooDataSet*)rds_bd->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_peak_tau = (RooDataSet*)rds_peak->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_semi_tau = (RooDataSet*)rds_semi->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_h2mu_tau = (RooDataSet*)rds_h2mu->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_comb_tau = (RooDataSet*)rds_data_himass->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        
        RooRealVar *EffTau_bs = wspace->var("EffTau_bs");
        RooRealVar *EffTau_bd = wspace->var("EffTau_bd");
        
        BuildDecayTimePDF(wspace, Form("bs_%s",cat.id.Data()), cat.id, rds_bs_tau, "single", EffTau_bs);
        BuildDecayTimePDF(wspace, Form("bd_%s",cat.id.Data()), cat.id, rds_bd_tau, "single", EffTau_bd);
        BuildDecayTimePDF(wspace, Form("peak_%s",cat.id.Data()), cat.id, rds_peak_tau, "single");
        BuildDecayTimePDF(wspace, Form("semi_%s",cat.id.Data()), cat.id, rds_semi_tau, "single");
        BuildDecayTimePDF(wspace, Form("h2mu_%s",cat.id.Data()), cat.id, rds_h2mu_tau, "single");
        BuildDecayTimePDF(wspace, Form("comb_%s",cat.id.Data()), cat.id, rds_comb_tau, "double");
        
        ProducePDFProjectionPlot(wspace, Form("bs_%s",cat.id.Data()), rds_bs_tau, "tau:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("bd_%s",cat.id.Data()), rds_bd_tau, "tau:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("peak_%s",cat.id.Data()), rds_peak_tau, "tau:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("semi_%s",cat.id.Data()), rds_semi_tau, "tau:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_tau, "tau:bmm3:dblbins");
        ProducePDFProjectionPlot(wspace, Form("comb_%s",cat.id.Data()), rds_comb_tau, "tau:bmm3:dblbins");
        
        // PDF for passing/not-passing the boundaries for tau
        wspace->import(RooGenericPdf(Form("SelCat_pdf_bs_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_bs->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_bs->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_bd_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_bd->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_bd->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_peak_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_peak->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_peak->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_semi_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_semi->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_semi->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_h2mu_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_h2mu->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_h2mu->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_comb_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_data_himass->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_data_himass->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        
        delete rds_bs_tau;
        delete rds_bd_tau;
        delete rds_peak_tau;
        delete rds_semi_tau;
        delete rds_h2mu_tau;
        delete rds_comb_tau;
    }
    
    {// also build the "mix" PDF
        RooDataSet *rds_bs_tau = (RooDataSet*)rds_bs->reduce(RooArgSet(*Tau,*SpCat),"SelCat==1");
        RooRealVar *EffTau_bs = wspace->var("EffTau_bs");
        BuildDecayTimePDF(wspace, Form("bs_mix"), "mix", rds_bs_tau, "single", EffTau_bs);

        ProducePDFProjectionPlot(wspace, Form("bs_mix"), rds_bs_tau, "tau:bmm3:dblbins");
    
        delete rds_bs_tau;
    }
    
    // build pair category PDFs
    for (auto& cat: CatMan.cats) {
        if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
        
        BuildPairCategoryPDF(wspace, Form("bs_%s",cat.id.Data()));
        BuildPairCategoryPDF(wspace, Form("bd_%s",cat.id.Data()));
        BuildPairCategoryPDF(wspace, Form("peak_%s",cat.id.Data()),"musq");
        BuildPairCategoryPDF(wspace, Form("semi_%s",cat.id.Data()),"mu");
        BuildPairCategoryPDF(wspace, Form("h2mu_%s",cat.id.Data()));
        BuildPairCategoryPDF(wspace, Form("comb_%s",cat.id.Data()));
    }
    
    // build mass PDFs
    for (auto& cat: CatMan.cats) {
        if (!cat.era.Contains("2011") && !cat.era.Contains("2012")) continue;
        
        // Bs->mumu
        RooDataSet *rds_bs_reduced = (RooDataSet*)rds_bs->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildSignalMassPDF(wspace, Form("bs_%s",cat.id.Data()), _bs_signal, rds_bs_reduced);
        ProducePDFProjectionPlot(wspace, Form("bs_%s",cat.id.Data()), rds_bs_reduced, "m:bmm3:dblbins");
        delete rds_bs_reduced;
        
        // Bd->mumu
        RooDataSet *rds_bd_reduced = (RooDataSet*)rds_bd->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildSignalMassPDF(wspace, Form("bd_%s",cat.id.Data()), _bd_signal, rds_bd_reduced);
        ProducePDFProjectionPlot(wspace, Form("bd_%s",cat.id.Data()), rds_bd_reduced, "m:bmm3:dblbins");
        delete rds_bd_reduced;
        
        // semi backgroud PDF
        RooDataSet *rds_semi_reduced = (RooDataSet*)rds_semi->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildSemiMassPDF(wspace, Form("semi_%s",cat.id.Data()), rds_semi_reduced);
        ProducePDFProjectionPlot(wspace, Form("semi_%s",cat.id.Data()), rds_semi_reduced, "m:bmm3");
        delete rds_semi_reduced;
        
        // h2mu backgroud PDF
        RooDataSet *rds_h2mu_reduced = (RooDataSet*)rds_h2mu->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildH2muMassPDF(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_reduced);
        ProducePDFProjectionPlot(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_reduced, "m:bmm3");
        delete rds_h2mu_reduced;
        
        // peak backgroud PDF
        RooDataSet *rds_peak_reduced = (RooDataSet*)rds_peak->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildPeakMassPDF(wspace, Form("peak_%s",cat.id.Data()), rds_peak_reduced);
        ProducePDFProjectionPlot(wspace, Form("peak_%s",cat.id.Data()), rds_peak_reduced, "m:bmm3");
        delete rds_peak_reduced;
        
        // comb background PDF
        BuildCombMassPDF(wspace, Form("comb_%s",cat.id.Data()));
    }
    
    for (int i=0; i<bmm3::ndecays; i++) // clean up
        delete samples[i];
    delete rds_data_lowbdt;
    delete rds_data_himass;
    delete rds_peak;
    delete rds_semi;
    delete rds_h2mu;
}

void PrepareBMM4SubPDF(RooWorkspace* wspace, RooWorkspace* wspace_base = 0)
{
    cout << ">>> PrepareBMM4SubPDF() start" << endl;
    
    if (wspace_base!=0) {
        cout << ">>> import from existing workspace." << endl;
        for (auto& cat: CatMan.cats) {
            for (TString type : {"Mass","BDT","PairCat","Tau","SelCat"}) {
                wspace->import(*wspace_base->pdf(Form("%s_pdf_bs_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_bd_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_peak_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_semi_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_h2mu_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
                wspace->import(*wspace_base->pdf(Form("%s_pdf_comb_%s",type.Data(),cat.id.Data())),RecycleConflictNodes());
            }
        }
        return;
    }
    cout << ">>> build PDFs from trees & data cards." << endl;
    
    vector<RooDataSet*> samples;
    samples.assign(bmm4::ndecays,NULL);
    
    // create the RooDataSet
    RooRealVar *Mass = wspace->var("Mass");
    RooRealVar *ReducedMassRes = wspace->var("ReducedMassRes");
    RooRealVar *BDT = wspace->var("BDT");
    RooRealVar *Tau = wspace->var("Tau");
    RooRealVar *TauRes = wspace->var("TauRes");
    RooCategory *GlobalCat = wspace->cat("GlobalCat");
    RooCategory *SpCat = wspace->cat("SpCat");
    RooCategory *SelCat = wspace->cat("SelCat");
    RooRealVar *Weight = wspace->var("Weight");
    
    RooArgSet varlist(*Mass,*ReducedMassRes,*BDT,*Tau,*TauRes,*GlobalCat,*SpCat,*SelCat,*Weight);

    for (int i=0; i<bmm4::ndecays; i++)
        samples[i] = new RooDataSet(Form("rds_%s",bmm4::decays[i].Data()), "", varlist, "Weight");
    RooDataSet *rds_data_lowbdt = new RooDataSet("rds_data_lowbdt", "", varlist, "Weight");
    RooDataSet *rds_data_himass = new RooDataSet("rds_data_himass", "", varlist, "Weight");

    // allocate counters for processing each categories
    vector<double> yields[bmm4::ndecays], counts[bmm4::ndecays];
    for (int idx_sample=0; idx_sample<bmm4::ndecays; idx_sample++) {
        yields[idx_sample].assign(CatMan.cats.size(),0.);
        counts[idx_sample].assign(CatMan.cats.size(),0.);
    }
    
    // read the expected yields from data cards
    for (int idx_sample=bmm4::_bskkMcBg; idx_sample<=bmm4::_bupimumuMcBg; idx_sample++) {
        for (auto& cat: CatMan.cats) {
            TString TexSource = Form("input/bmm4/scanBDT-%s.tex",cat.era.Data());
            int bdt_min = (int)(cat.bdt_min*100.);
            int bdt_max = (int)(cat.bdt_max*100.);
            
            TexVar var;
            for (int bin=0; bin<=3; bin++) {
                TexVar v(TexSource, Form("bdt_%d_%s:N-SCALEDYIELD-MBIN%d-%s-chan%d",bdt_min,cat.era.Data(),bin,bmm4::texlabels[idx_sample].Data(),cat.region));
                if (bdt_max<100) v.SubVar(TexVar(TexSource, Form("bdt_%d_%s:N-SCALEDYIELD-MBIN%d-%s-chan%d",
                                                                 bdt_max,cat.era.Data(),bin,bmm4::texlabels[idx_sample].Data(),cat.region)));
                var.AddVar(v);
            }
            cout << ">>> Expected yield for " << cat.id << ": " << bmm4::texlabels[idx_sample] << ": " << var.val << " +- " << var.estat << " +- " << var.esyst << endl;
            yields[idx_sample][cat.index] = var.val; // yields
        }
    }
    
    // loop over the samples
    for (int idx_sample=bmm4::_bsmmMc; idx_sample<=bmm4::_bmmData; idx_sample++) {
        for (TString era : bmm4::eras) {
            
            bool applyWeights = true;
            if (idx_sample==bmm4::_bsmmMc || idx_sample==bmm4::_bdmmMc || idx_sample==bmm4::_bmmData) applyWeights = false;
            
            TString filename = Form("input/bmm4/small%s-%s.root",era.Data(),bmm4::decays[idx_sample].Data());
            TString treename = bmm4::decays[idx_sample];
            
            exist_protection(filename);
            TFile *fin = new TFile(filename);
            TTree *tin = (TTree*)fin->Get(treename);
            
            cout << ">>> Loading from " << filename << ", with " << tin->GetEntries() << " entries." << endl;
            
            double m_t, me_t, bdt_t, tau_t, taue_t;
            int chan_t;
            bool muid_t;
            tin->SetBranchAddress("m", &m_t);
            tin->SetBranchAddress("me", &me_t);
            tin->SetBranchAddress("bdt", &bdt_t);
            tin->SetBranchAddress("tau", &tau_t);
            tin->SetBranchAddress("taue", &taue_t);
            tin->SetBranchAddress("chan",  &chan_t);
            tin->SetBranchAddress("muid",  &muid_t);
            
            // first loop: event counting, peaking/semi MC only
            int nevt_limit = 0;
            if (applyWeights)
                for (int evt=0; evt<tin->GetEntries();evt++) {
                    tin->GetEntry(evt);
                    
                    if (m_t < Mass_bound[0] || m_t > Mass_bound[1] || !muid_t) continue;
                    if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
                    if (bdt_t < -1.) continue;
                    
                    int index = CatMan.index(era, chan_t, bdt_t);
                    if (index<0) continue;
                    
                    counts[idx_sample][index] += 1.; // adding counts
                    
                    nevt_limit++;
                    if (nevt_limit>=MC_EVENT_LIMIT) break; // read up to MC_EVENT_LIMIT per sample
                }
            
            // second loop: event filling
            nevt_limit = 0;
            for (int evt=0; evt<tin->GetEntries();evt++) {
                tin->GetEntry(evt);

                if (m_t < Mass_bound[0] || m_t > Mass_bound[1]) continue;
                if (me_t/m_t < ReducedMassRes_bound[0] || me_t/m_t > ReducedMassRes_bound[1]) continue;
                if (bdt_t < -1.) continue;
                
                int index = CatMan.index(era, chan_t, bdt_t);

                if (tau_t*1E12 > Tau_bound[0] && tau_t*1E12 < Tau_bound[1] &&
                    taue_t*1E12 > TauRes_bound[0] && taue_t*1E12 < TauRes_bound[1]) SelCat->setIndex(1);
                else SelCat->setIndex(0);
                
                Mass->setVal(m_t);
                ReducedMassRes->setVal(me_t/m_t);
                BDT->setVal(bdt_t);
                Tau->setVal(tau_t*1E12);
                TauRes->setVal(taue_t*1E12);
                GlobalCat->setIndex(index);
                Weight->setVal(1.);
                if (applyWeights)
                    Weight->setVal(yields[idx_sample][index]/counts[idx_sample][index]);
                SpCat->setIndex(idx_sample);
                
                if (idx_sample == bmm4::_bmmData) {

                    for (auto& cat: CatMan.cats) { // fill event regardless of BDT binning for low BDT sideband
                        if (cat.era == era && cat.region == chan_t) {
                            if (!muid_t) continue; // reject failed muon ID events
                            if (m_t > 5.20 && m_t < 5.45) continue; // exclude the signal for data sideband
                            if (bdt_t < BDT_bound[0]-0.3 || bdt_t > BDT_bound[0]) continue; // exclude good BDT region
                            GlobalCat->setIndex(cat.index);
                            rds_data_lowbdt->add(varlist, Weight->getVal());
                        }
                    }
                    
                    for (auto& cat: CatMan.cats) { // fill event regardless of BDT binning for high mass sideband
                        if (cat.era == era && cat.region == chan_t) {
                            if (m_t < 5.45) continue; // exclude the signal & low mass regions
                            if (bdt_t < BDT_bound[0]) continue; // exclude low BDT region
                            GlobalCat->setIndex(cat.index);
                            rds_data_himass->add(varlist, Weight->getVal());
                        }
                    }
                    
                    nevt_limit++;
                }else {
                    if (!muid_t) continue; // reject failed muon ID events
                    if (index<0) continue; // reject non-categorized events
                    
                    samples[idx_sample]->add(varlist, Weight->getVal());
                    
                    nevt_limit++;
                    if (idx_sample!=bmm4::_bmmData && nevt_limit>=MC_EVENT_LIMIT) break; // read up to MC_EVENT_LIMIT per sample
                }
            } // end of evt loop
            
            cout << ">>> " << nevt_limit << " events loaded." << endl;
            
            fin->Close();
            
        } // end of "era" loop
    } // end of "idx_sample" loop
    
    // alias & merge
    RooDataSet *rds_bs = samples[bmm4::_bsmmMc];
    RooDataSet *rds_bd = samples[bmm4::_bdmmMc];
    RooDataSet* rds_peak = (RooDataSet*)samples[bmm4::_bskkMcBg]->Clone("rds_peak");
    for (int i=bmm4::_bskkMcBg+1; i<=bmm4::_lbpkMcBg; i++) rds_peak->append(*samples[i]);
    RooDataSet* rds_semi = (RooDataSet*)samples[bmm4::_bskmunuMcBg]->Clone("rds_semi");
    for (int i=bmm4::_bskmunuMcBg+1; i<=bmm4::_lbpmunuMcBg; i++) rds_semi->append(*samples[i]);
    RooDataSet* rds_h2mu = (RooDataSet*)samples[bmm4::_bdpimumuMcBg]->Clone("rds_h2mu");
    for (int i=bmm4::_bdpimumuMcBg+1; i<=bmm4::_bupimumuMcBg; i++) rds_h2mu->append(*samples[i]);
    
    // build ReducedMassRes PDF
    for (auto& cat: CatMan.cats) {
        
        RooDataSet *rds_bs_res = (RooDataSet*)rds_bs->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_bd_res = (RooDataSet*)rds_bd->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_peak_res = (RooDataSet*)rds_peak->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_semi_res = (RooDataSet*)rds_semi->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_h2mu_res = (RooDataSet*)rds_h2mu->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_comb_res = (RooDataSet*)rds_data_lowbdt->reduce(RooArgSet(*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        
        BuildReducedMassResPDF(wspace, Form("bs_%s",cat.id.Data()), rds_bs_res);
        BuildReducedMassResPDF(wspace, Form("bd_%s",cat.id.Data()), rds_bd_res);
        BuildReducedMassResPDF(wspace, Form("peak_%s",cat.id.Data()), rds_peak_res);
        BuildReducedMassResPDF(wspace, Form("semi_%s",cat.id.Data()), rds_semi_res);
        BuildReducedMassResPDF(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_res);
        BuildReducedMassResPDF(wspace, Form("comb_%s",cat.id.Data()), rds_comb_res);
        
        ProducePDFProjectionPlot(wspace, Form("bs_%s",cat.id.Data()), rds_bs_res, "me:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("bd_%s",cat.id.Data()), rds_bd_res, "me:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("peak_%s",cat.id.Data()), rds_peak_res, "me:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("semi_%s",cat.id.Data()), rds_semi_res, "me:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_res, "me:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("comb_%s",cat.id.Data()), rds_comb_res, "me:bmm4:dblbins");
        
        delete rds_bs_res;
        delete rds_bd_res;
        delete rds_peak_res;
        delete rds_semi_res;
        delete rds_h2mu_res;
        delete rds_comb_res;
    }
    
    // build BDT PDF
    for (auto& cat: CatMan.cats) {
        
        RooDataSet *rds_bs_bdt = (RooDataSet*)rds_bs->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_bd_bdt = (RooDataSet*)rds_bd->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_peak_bdt = (RooDataSet*)rds_peak->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_semi_bdt = (RooDataSet*)rds_semi->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_h2mu_bdt = (RooDataSet*)rds_h2mu->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        RooDataSet *rds_comb_bdt = (RooDataSet*)rds_data_himass->reduce(RooArgSet(*BDT,*SpCat),Form("GlobalCat==%d",cat.index));
        
        BuildBDTPDF(wspace, Form("bs_%s",cat.id.Data()), rds_bs_bdt);
        BuildBDTPDF(wspace, Form("bd_%s",cat.id.Data()), rds_bd_bdt);
        BuildBDTPDF(wspace, Form("peak_%s",cat.id.Data()), rds_peak_bdt);
        BuildBDTPDF(wspace, Form("semi_%s",cat.id.Data()), rds_semi_bdt);
        BuildBDTPDF(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_bdt);
        BuildBDTPDF(wspace, Form("comb_%s",cat.id.Data()), rds_comb_bdt);
        
        ProducePDFProjectionPlot(wspace, Form("bs_%s",cat.id.Data()), rds_bs_bdt, "bdt:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("bd_%s",cat.id.Data()), rds_bd_bdt, "bdt:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("peak_%s",cat.id.Data()), rds_peak_bdt, "bdt:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("semi_%s",cat.id.Data()), rds_semi_bdt, "bdt:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_bdt, "bdt:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("comb_%s",cat.id.Data()), rds_comb_bdt, "bdt:bmm4:dblbins");
        
        delete rds_bs_bdt;
        delete rds_bd_bdt;
        delete rds_peak_bdt;
        delete rds_semi_bdt;
        delete rds_h2mu_bdt;
        delete rds_comb_bdt;
    }
    
    // build decay time PDFs
    for (auto& cat: CatMan.cats) {
        
        RooDataSet *rds_bs_tau = (RooDataSet*)rds_bs->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_bd_tau = (RooDataSet*)rds_bd->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_peak_tau = (RooDataSet*)rds_peak->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_semi_tau = (RooDataSet*)rds_semi->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_h2mu_tau = (RooDataSet*)rds_h2mu->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        RooDataSet *rds_comb_tau = (RooDataSet*)rds_data_himass->reduce(RooArgSet(*Tau,*SpCat),Form("SelCat==1&&GlobalCat==%d",cat.index));
        
        RooRealVar *EffTau_bs = wspace->var("EffTau_bs");
        RooRealVar *EffTau_bd = wspace->var("EffTau_bd");
        
        BuildDecayTimePDF(wspace, Form("bs_%s",cat.id.Data()), cat.id, rds_bs_tau, "single", EffTau_bs);
        BuildDecayTimePDF(wspace, Form("bd_%s",cat.id.Data()), cat.id, rds_bd_tau, "single", EffTau_bd);
        BuildDecayTimePDF(wspace, Form("peak_%s",cat.id.Data()), cat.id, rds_peak_tau, "single");
        BuildDecayTimePDF(wspace, Form("semi_%s",cat.id.Data()), cat.id, rds_semi_tau, "single");
        BuildDecayTimePDF(wspace, Form("h2mu_%s",cat.id.Data()), cat.id, rds_h2mu_tau, "single");
        BuildDecayTimePDF(wspace, Form("comb_%s",cat.id.Data()), cat.id, rds_comb_tau, "double");
        
        ProducePDFProjectionPlot(wspace, Form("bs_%s",cat.id.Data()), rds_bs_tau, "tau:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("bd_%s",cat.id.Data()), rds_bd_tau, "tau:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("peak_%s",cat.id.Data()), rds_peak_tau, "tau:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("semi_%s",cat.id.Data()), rds_semi_tau, "tau:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_tau, "tau:bmm4:dblbins");
        ProducePDFProjectionPlot(wspace, Form("comb_%s",cat.id.Data()), rds_comb_tau, "tau:bmm4:dblbins");
        
        // PDF for passing/not-passing the boundaries for tau
        wspace->import(RooGenericPdf(Form("SelCat_pdf_bs_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_bs->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_bs->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_bd_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_bd->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_bd->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_peak_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_peak->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_peak->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_semi_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_semi->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_semi->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_h2mu_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_h2mu->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_h2mu->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        wspace->import(RooGenericPdf(Form("SelCat_pdf_comb_%s",cat.id.Data()),"","@0?@1:@2",
                                     RooArgList(*SelCat,RooConst(rds_data_himass->sumEntries(Form("SelCat==1&&GlobalCat==%d",cat.index))),
                                                RooConst(rds_data_himass->sumEntries(Form("SelCat==0&&GlobalCat==%d",cat.index))))));
        
        delete rds_bs_tau;
        delete rds_bd_tau;
        delete rds_peak_tau;
        delete rds_semi_tau;
        delete rds_h2mu_tau;
        delete rds_comb_tau;
    }
    
    {// also build the "mix" PDF
        RooDataSet *rds_bs_tau = (RooDataSet*)rds_bs->reduce(RooArgSet(*Tau,*SpCat),"SelCat==1");
        RooRealVar *EffTau_bs = wspace->var("EffTau_bs");
        BuildDecayTimePDF(wspace, Form("bs_mix"), "mix", rds_bs_tau, "single", EffTau_bs);
        
        ProducePDFProjectionPlot(wspace, Form("bs_mix"), rds_bs_tau, "tau:bmm4:dblbins");
        
        delete rds_bs_tau;
    }
    
    // build pair category PDFs
    for (auto& cat: CatMan.cats) {
        
        BuildPairCategoryPDF(wspace, Form("bs_%s",cat.id.Data()));
        BuildPairCategoryPDF(wspace, Form("bd_%s",cat.id.Data()));
        BuildPairCategoryPDF(wspace, Form("peak_%s",cat.id.Data()),"musq");
        BuildPairCategoryPDF(wspace, Form("semi_%s",cat.id.Data()),"mu");
        BuildPairCategoryPDF(wspace, Form("h2mu_%s",cat.id.Data()));
        BuildPairCategoryPDF(wspace, Form("comb_%s",cat.id.Data()));
    }
    
    // build mass PDFs
    for (auto& cat: CatMan.cats) {
        
        // Bs->mumu
        RooDataSet *rds_bs_reduced = (RooDataSet*)rds_bs->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildSignalMassPDF(wspace, Form("bs_%s",cat.id.Data()), _bs_signal, rds_bs_reduced);
        ProducePDFProjectionPlot(wspace, Form("bs_%s",cat.id.Data()), rds_bs_reduced, "m:bmm4:dblbins");
        delete rds_bs_reduced;
        
        // Bd->mumu
        RooDataSet *rds_bd_reduced = (RooDataSet*)rds_bd->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildSignalMassPDF(wspace, Form("bd_%s",cat.id.Data()), _bd_signal, rds_bd_reduced);
        ProducePDFProjectionPlot(wspace, Form("bd_%s",cat.id.Data()), rds_bd_reduced, "m:bmm4:dblbins");
        delete rds_bd_reduced;
        
        // semi backgroud PDF
        RooDataSet *rds_semi_reduced = (RooDataSet*)rds_semi->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildSemiMassPDF(wspace, Form("semi_%s",cat.id.Data()), rds_semi_reduced);
        ProducePDFProjectionPlot(wspace, Form("semi_%s",cat.id.Data()), rds_semi_reduced, "m:bmm4");
        delete rds_semi_reduced;
        
        // h2mu backgroud PDF
        RooDataSet *rds_h2mu_reduced = (RooDataSet*)rds_h2mu->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildH2muMassPDF(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_reduced);
        ProducePDFProjectionPlot(wspace, Form("h2mu_%s",cat.id.Data()), rds_h2mu_reduced, "m:bmm4");
        delete rds_h2mu_reduced;
        
        // peak backgroud PDF
        RooDataSet *rds_peak_reduced = (RooDataSet*)rds_peak->reduce(RooArgSet(*Mass,*ReducedMassRes,*SpCat),Form("GlobalCat==%d",cat.index));
        BuildPeakMassPDF(wspace, Form("peak_%s",cat.id.Data()), rds_peak_reduced);
        ProducePDFProjectionPlot(wspace, Form("peak_%s",cat.id.Data()), rds_peak_reduced, "m:bmm4");
        delete rds_peak_reduced;
        
        // comb background PDF
        BuildCombMassPDF(wspace, Form("comb_%s",cat.id.Data()));
    }
    
    for (int i=0; i<bmm4::ndecays; i++) // clean up
        delete samples[i];
    delete rds_data_lowbdt;
    delete rds_data_himass;
    delete rds_peak;
    delete rds_semi;
    delete rds_h2mu;
}

// bdtcat  - BDT-Cut or BDT-Category strategy
// bdtfit  - BDT-Fit strategy
void PrepareGlobalPDF(RooWorkspace* wspace, TString opt = "bdtcat")
{
    cout << ">>> PrepareGlobalPDF() start" << endl;
    
    RooRealVar *BF_bs = wspace->var("BF_bs");
    RooRealVar *BF_bd = wspace->var("BF_bd");
    RooRealVar *one_over_BRBR = wspace->var("one_over_BRBR");
    RooRealVar *fs_over_fu = wspace->var("fs_over_fu");
    RooRealVar *fs_over_fu_S13 = wspace->var("fs_over_fu_S13");
    RooRealVar *dblmu_corr_scale = wspace->var("dblmu_corr_scale");
    RooCategory *GlobalCat = wspace->cat("GlobalCat");
    
    RooSimultaneous *global_pdf = new RooSimultaneous("global_pdf", "", *GlobalCat);
    
    RooFormulaVar *N_bs_formula[CatMan.cats.size()],*N_bd_formula[CatMan.cats.size()];
    RooFormulaVar *N_semi_formula[CatMan.cats.size()],*N_peak_formula[CatMan.cats.size()];
    RooAddPdf  *pdf_ext_sum[CatMan.cats.size()];
    RooProdPdf *constraints_pdfs[CatMan.cats.size()];
    RooProdPdf *pdf_ext_total[CatMan.cats.size()];
    
    // Now build individual PDF
    for (auto& cat: CatMan.cats) {
        for (TString spec : {"bs","bd","peak","semi","h2mu","comb"}) {
            RooArgList subpdf_list;
            subpdf_list.add(*wspace->pdf(Form("Mass_pdf_%s_%s", spec.Data(), cat.id.Data())));
            subpdf_list.add(*wspace->pdf(Form("PairCat_pdf_%s_%s", spec.Data(), cat.id.Data())));
            if (opt=="bdtfit")
                subpdf_list.add(*wspace->pdf(Form("BDT_pdf_%s_%s", spec.Data(), cat.id.Data())));
            
            RooProdPdf subpdf(Form("pdf_%s_%s", spec.Data(), cat.id.Data()), "", subpdf_list);
            wspace->import(subpdf);
        }
    }
    
    for (auto& cat: CatMan.cats) {
        
        RooRealVar *N_bu = wspace->var(Form("N_bu_%s", cat.id.Data()));
        RooRealVar *effratio_bs = wspace->var(Form("effratio_bs_%s", cat.id.Data()));
        RooRealVar *effratio_bd = wspace->var(Form("effratio_bd_%s", cat.id.Data()));
        
        if (cat.era.Contains("2016"))
            N_bs_formula[cat.index] = new RooFormulaVar(Form("N_bs_formula_%s", cat.id.Data()), "", "@0*@1*@2*@3*@4*@5",
                                                        RooArgList(*BF_bs, *N_bu, *fs_over_fu, *fs_over_fu_S13, *effratio_bs, *one_over_BRBR));
        else
            N_bs_formula[cat.index] = new RooFormulaVar(Form("N_bs_formula_%s", cat.id.Data()), "", "@0*@1*@2*@3*@4",
                                                        RooArgList(*BF_bs, *N_bu, *fs_over_fu, *effratio_bs, *one_over_BRBR));
        N_bd_formula[cat.index] = new RooFormulaVar(Form("N_bd_formula_%s", cat.id.Data()), "", "@0*@1*@2*@3",
                                                  RooArgList(*BF_bd, *N_bu, *effratio_bd, *one_over_BRBR));
        
        RooRealVar *N_semi = wspace->var(Form("N_semi_%s", cat.id.Data()));
        RooRealVar *N_peak = wspace->var(Form("N_peak_%s", cat.id.Data()));
        
        N_semi_formula[cat.index] = new RooFormulaVar(Form("N_semi_formula_%s", cat.id.Data()), "", "@0*(1.+@1)*0.5",
                                                    RooArgList(*N_semi, *dblmu_corr_scale));
        N_peak_formula[cat.index] = new RooFormulaVar(Form("N_peak_formula_%s", cat.id.Data()), "", "@0*(1.+@1*@1)*0.5",
                                                    RooArgList(*N_peak, *dblmu_corr_scale));
        
        RooArgList pdf_list;
        pdf_list.add(*wspace->pdf(Form("pdf_bs_%s", cat.id.Data())));
        pdf_list.add(*wspace->pdf(Form("pdf_bd_%s", cat.id.Data())));
        pdf_list.add(*wspace->pdf(Form("pdf_comb_%s", cat.id.Data())));
        pdf_list.add(*wspace->pdf(Form("pdf_semi_%s", cat.id.Data())));
        pdf_list.add(*wspace->pdf(Form("pdf_h2mu_%s", cat.id.Data())));
        pdf_list.add(*wspace->pdf(Form("pdf_peak_%s", cat.id.Data())));
        
        RooArgList N_list;
        N_list.add(*N_bs_formula[cat.index]);
        N_list.add(*N_bd_formula[cat.index]);
        N_list.add(*wspace->var(Form("N_comb_%s", cat.id.Data())));
        N_list.add(*N_semi_formula[cat.index]);
        N_list.add(*wspace->var(Form("N_h2mu_%s", cat.id.Data())));
        N_list.add(*N_peak_formula[cat.index]);
        
        RooArgList constraints_list;
        constraints_list.add(*wspace->pdf(Form("N_bu_%s_gau", cat.id.Data())));
        constraints_list.add(*wspace->pdf(Form("N_peak_%s_lnn", cat.id.Data())));
        constraints_list.add(*wspace->pdf(Form("N_semi_%s_lnn", cat.id.Data())));
        constraints_list.add(*wspace->pdf(Form("N_h2mu_%s_lnn", cat.id.Data())));
        constraints_list.add(*wspace->pdf(Form("effratio_bs_%s_gau", cat.id.Data())));
        constraints_list.add(*wspace->pdf(Form("effratio_bd_%s_gau", cat.id.Data())));
        
        pdf_ext_sum[cat.index] = new RooAddPdf(Form("pdf_ext_sum_%s", cat.id.Data()), "", pdf_list, N_list);
        
        constraints_pdfs[cat.index] = new RooProdPdf(Form("pdf_constraints_%s", cat.id.Data()), "", constraints_list);
        
        pdf_ext_total[cat.index] = new RooProdPdf(Form("pdf_ext_total_%s", cat.id.Data()), "", RooArgList(*pdf_ext_sum[cat.index], *constraints_pdfs[cat.index]));
        
        global_pdf->addPdf(*pdf_ext_total[cat.index], cat.id);
    }
    
    wspace->import(*global_pdf);
}

void PrepareNamedArgSets(RooWorkspace* wspace)
{
    cout << ">>> PrepareNamedArgSets() start" << endl;
    
    RooArgSet all_nuisances;
    
    all_nuisances.add(*wspace->var("fs_over_fu"));
    all_nuisances.add(*wspace->var("fs_over_fu_S13"));
    all_nuisances.add(*wspace->var("one_over_BRBR"));
    
    for (auto& cat: CatMan.cats) {
        all_nuisances.add(*wspace->var(Form("N_bu_%s", cat.id.Data())));
        all_nuisances.add(*wspace->var(Form("effratio_bs_%s", cat.id.Data())));
        all_nuisances.add(*wspace->var(Form("effratio_bd_%s", cat.id.Data())));
        all_nuisances.add(*wspace->var(Form("N_semi_%s", cat.id.Data())));
        all_nuisances.add(*wspace->var(Form("N_h2mu_%s", cat.id.Data())));
        all_nuisances.add(*wspace->var(Form("N_peak_%s", cat.id.Data())));
    }
    
    RooArgSet all_floats(all_nuisances);
    
    all_floats.add(*wspace->var("BF_bs"));
    all_floats.add(*wspace->var("BF_bd"));
    all_floats.add(*wspace->var("dblmu_corr_scale"));
    
    for (auto& cat: CatMan.cats) {
        all_floats.add(*wspace->var(Form("N_comb_%s", cat.id.Data())));
        all_floats.add(*wspace->var(Form("Br1_comb_%s", cat.id.Data())));
    }
    
    wspace->defineSet("all_nuisances",all_nuisances); // all nuisances
    wspace->defineSet("all_floats",all_floats); // all floated parameters
    
    // Summary table for the yields
    for (auto& cat: CatMan.cats) {
        cout << ">>> " << string(32,'-') << endl;
        cout << ">>> CATEGORY ID " << cat.id << endl;
        cout << ">>> " << string(32,'-') << endl;
        
        RooAbsReal *N_bs_formula = wspace->function(Form("N_bs_formula_%s", cat.id.Data()));
        RooAbsReal *N_bd_formula = wspace->function(Form("N_bd_formula_%s", cat.id.Data()));
        
        RooRealVar *N_bu_mean = wspace->var(Form("N_bu_%s_mean", cat.id.Data()));
        RooRealVar *N_peak_mean = wspace->var(Form("N_peak_%s_mean", cat.id.Data()));
        RooRealVar *N_semi_mean = wspace->var(Form("N_semi_%s_mean", cat.id.Data()));
        RooRealVar *N_h2mu_mean = wspace->var(Form("N_h2mu_%s_mean", cat.id.Data()));

        RooRealVar *N_bu_sigma = wspace->var(Form("N_bu_%s_sigma", cat.id.Data()));
        RooRealVar *N_peak_kappa = wspace->var(Form("N_peak_%s_kappa", cat.id.Data()));
        RooRealVar *N_semi_kappa = wspace->var(Form("N_semi_%s_kappa", cat.id.Data()));
        RooRealVar *N_h2mu_kappa = wspace->var(Form("N_h2mu_%s_kappa", cat.id.Data()));
        
        RooRealVar *N_comb = wspace->var(Form("N_comb_%s", cat.id.Data()));
        
        cout << ">>> N_bu:   mu = " << N_bu_mean->getVal() << ", sigma = " << N_bu_sigma->getVal() << endl;
        cout << ">>> N_peak: mu = " << N_peak_mean->getVal() << ", kappa = " << N_peak_kappa->getVal() << endl;
        cout << ">>> N_semi: mu = " << N_semi_mean->getVal() << ", kappa = " << N_semi_kappa->getVal() << endl;
        cout << ">>> N_h2mu: mu = " << N_h2mu_mean->getVal() << ", kappa = " << N_h2mu_kappa->getVal() << endl;
        cout << ">>> N_bs:   mu = " << N_bs_formula->getVal() << endl;
        cout << ">>> N_bd:   mu = " << N_bd_formula->getVal() << endl;
        cout << ">>> N_comb: mu = " << N_comb->getVal() << endl;
    }
}

void ProduceSubPlots(RooWorkspace *wspace)
{
    cout << ">>> ProduceSubPlots() start" << endl;
    
    for (auto& cat: CatMan.cats) {
        
        TString title_base;
        if (cat.era.Contains("2011") && cat.region == 0) title_base = "CMS - L = 5 fb^{-1} #sqrt{s} = 7 TeV - 2011 Barrel";
        if (cat.era.Contains("2011") && cat.region == 1) title_base = "CMS - L = 5 fb^{-1} #sqrt{s} = 7 TeV - 2011 Endcap";
        if (cat.era.Contains("2012") && cat.region == 0) title_base = "CMS - L = 20 fb^{-1} #sqrt{s} = 8 TeV - 2012 Barrel";
        if (cat.era.Contains("2012") && cat.region == 1) title_base = "CMS - L = 20 fb^{-1} #sqrt{s} = 8 TeV - 2012 Endcap";
        if (cat.era.Contains("2016BF")) title_base = Form("CMS - L = 36 fb^{-1} #sqrt{s} = 13 TeV - 2016 B-F Channel %d",cat.region);
        if (cat.era.Contains("2016GH")) title_base = Form("CMS - L = 36 fb^{-1} #sqrt{s} = 13 TeV - 2016 G-H Channel %d",cat.region);
        
        RooRealVar *Mass = wspace->var("Mass");
        
        TString cut = Form("GlobalCat==%d",cat.index);
        TString title = Form("%s - BDT bin %d", title_base.Data(),cat.bdt_bin);
        
        RooPlot* frame = Mass->frame(Bins(25), Title(" "));
        wspace->data("global_data")->plotOn(frame, Cut(cut), Invisible());
        
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
        
        TH1D* hist_data = (TH1D*)wspace->data("global_data")->createHistogram(Form("hist_data_%s",cat.id.Data()), *Mass, Cut(cut), Binning(25, Mass_bound[0], Mass_bound[1]));
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
        TLegendEntry *data_entry = new TLegendEntry(hist_data, "data", "lep");
        data_entry->SetMarkerStyle(20);
        leg1->AddEntry(data_entry, "data", "ep");
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
        
        canvas->Print(Form("fig/proj_bin_%s.pdf",cat.id.Data()));
        
        delete frame;
        delete canvas;
    }
}

void ProduceTauSPlot(RooWorkspace *wspace, bool sys_study = false)
{
    cout << ">>> ProduceTauSPlot() start" << endl;
    
    enum {_bs, _bd, _peak, _semi, _h2mu, _comb, _nspec};
    vector<TString> specs= {"bs", "bd", "peak", "semi", "h2mu", "comb"};
    
    TH1D *h_tau[_nspec];
    for (int idx = 0; idx < _nspec; idx++) {
        h_tau[idx] = new TH1D(Form("h_tau_%s",specs[idx].Data()), "", Tau_bins.size()-1, Tau_bins.data());
        h_tau[idx]->Sumw2();
    }
    
    RooFitResult *res1 = NULL, *res2 = NULL;
    
    vector<TH1D*> h_tau_sys;
    vector<double> diff_tau_sys;
    vector<TString> cmd_tau_sys;
    
    if (sys_study) {
        cmd_tau_sys.push_back(Form("0HistEffModel"));
        cmd_tau_sys.push_back(Form("+respar_bsmm_mix_taue"));
        cmd_tau_sys.push_back(Form("-respar_bsmm_mix_taue"));
        cmd_tau_sys.push_back(Form("+fs_over_fu_S13"));
        cmd_tau_sys.push_back(Form("-fs_over_fu_S13"));
        for (auto& cat: CatMan.cats) {
            cmd_tau_sys.push_back(Form("+N_bu_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("-N_bu_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("+effratio_bs_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("-effratio_bs_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("+effratio_bd_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("-effratio_bd_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("+N_semi_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("-N_semi_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("+N_h2mu_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("-N_h2mu_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("+N_peak_%s", cat.id.Data()));
            cmd_tau_sys.push_back(Form("-N_peak_%s", cat.id.Data()));
        }
        
        h_tau_sys.assign(cmd_tau_sys.size(),NULL);
        diff_tau_sys.assign(cmd_tau_sys.size(),0.);
        for (int sysidx = 0; sysidx < (int)cmd_tau_sys.size(); sysidx++) {
            h_tau_sys[sysidx] = new TH1D(Form("h_tau_sys_%03d",sysidx), "", Tau_bins.size()-1, Tau_bins.data());
            h_tau_sys[sysidx]->Sumw2();
        }
    }
    
    for (int sysidx = -1; sysidx < (int)cmd_tau_sys.size(); sysidx++) {
        
        if (!sys_study && sysidx>=0) break; // skip nuisances variations
        
        double nuisance_preserve = 0.;
        if (sysidx>=0) {
            TString var = cmd_tau_sys[sysidx](0,1);
            TString key = cmd_tau_sys[sysidx](1,cmd_tau_sys[sysidx].Length()-1);
            
            // Take post-fit nuisance variation
            if (var!="0") {
                RooRealVar *nuisance = wspace->var(key);
                nuisance_preserve = nuisance->getVal();
                if (var=="-") nuisance->setVal(nuisance_preserve-nuisance->getError());
                if (var=="+") nuisance->setVal(nuisance_preserve+nuisance->getError());
            }
        }
        
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
            
            for (int evt=0; evt<wspace->data("global_data")->numEntries(); evt++) {
                const RooArgSet* arg = wspace->data("global_data")->get(evt);
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
            
            double weight_total[_nspec];
            for (int idx = 0; idx < _nspec; idx++) weight_total[idx] = 0.;
            
            for (int evt=0; evt<wspace->data("global_data")->numEntries(); evt++) {
                const RooArgSet* arg = wspace->data("global_data")->get(evt);
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
                
                if (sysidx<0) { // central sPlot
                    for (int target_spec = 0; target_spec<_nspec; target_spec++) {
                        double numerator = 0.;
                        for (int idx = 0; idx<_nspec; idx++) numerator += covMatrix(target_spec,idx)*pdf[idx];
                    
                        double weight = numerator/denominator;
                        h_tau[target_spec]->Fill(arg->getRealValue("Tau"),weight);
                        weight_total[target_spec] += weight;
                    }
                }else {
                    double numerator = 0.;
                    for (int idx = 0; idx<_nspec; idx++) numerator += covMatrix(_bs,idx)*pdf[idx];

                    double weight = numerator/denominator;
                    h_tau_sys[sysidx]->Fill(arg->getRealValue("Tau"),weight);
                }
            }
            if (sysidx<0) { // central sPlot
                cout << ">>> sPlot check sum for bin: " << cat.id << endl;
                for (int idx = 0; idx<_nspec; idx++)
                    cout << ">>> Spec " << specs[idx] << ": " << weight_total[idx] << " / " << yield[idx] << endl;
            }
        }
        
        RooRealVar *Tau = wspace->var("Tau");
        RooRealVar *EffTau_bs = wspace->var("EffTau_bs");
        RooAbsPdf *Tau_pdf_bs = wspace->pdf("Tau_pdf_bs_mix");
        
        EffTau_bs->setConstant(false);
        EffTau_bs->setVal(1.61); // always start from SM value
        
        if (sysidx<0) {
            Fit_sPlot(h_tau[_bs], Tau_pdf_bs, Tau, EffTau_bs, &res1, &res2);
            wspace->import(*res1,"fitresult_taubs"); // results with corrected uncertainties
            wspace->import(*res2,"fitresult_taubs_wl"); // results from the 2nd weighted likelihood fit
            
        }else {
            TString key = cmd_tau_sys[sysidx](1,cmd_tau_sys[sysidx].Length()-1);
            
            RooFitResult *res1sys = NULL, *res2sys = NULL;
            
            if (key=="HistEffModel") {
                RooResolutionModel *TauRes_Model = (RooResolutionModel*)wspace->obj(Form("TauRes_Model_bsmm_mix"));
                RooFormulaVar *TauEff_Model = (RooFormulaVar*)wspace->obj(Form("TauEff_Model_Hist_bsmm_mix"));
                
                RooDecay RawDecay("RawDecay_alter","",*Tau,*EffTau_bs,*TauRes_Model,RooDecay::SingleSided);
                RooEffProd Tau_pdf("Tau_pdf_alter","",RawDecay,*TauEff_Model);
                Fit_sPlot(h_tau_sys[sysidx], &Tau_pdf, Tau, EffTau_bs, &res1sys, &res2sys);
            }else
                Fit_sPlot(h_tau_sys[sysidx], Tau_pdf_bs, Tau, EffTau_bs, &res1sys, &res2sys);
            
            RooRealVar* tau0 = (RooRealVar*)res1->floatParsFinal().find("EffTau_bs");
            RooRealVar* tau1 = (RooRealVar*)res1sys->floatParsFinal().find("EffTau_bs");
            diff_tau_sys[sysidx] = tau1->getVal()-tau0->getVal();
            
            delete res1sys;
            delete res2sys;
        }
        
        if (sysidx<0) // central sPlot
        for (int target_spec = 0; target_spec<_nspec; target_spec++) {
            
            TString title = "CMS Preliminary";
            
            RooRealVar ps("ps","",0.,14.);
            RooDataHist *h_tau_data = new RooDataHist("h_tau_data", "", RooArgList(ps), h_tau[target_spec]);
            RooPlot* frame = ps.frame(Title(" "));
            h_tau[target_spec]->SetMarkerStyle(20);
            h_tau[target_spec]->SetLineColor(kBlack);
            h_tau_data->plotOn(frame,MarkerStyle(20),LineColor(kBlack));
            
            if (target_spec==_bs)
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
            lin.DrawLine(0.,0.,14.,0.);
            
            canvas->Update();
            
            TLegend *leg1 = new TLegend(0.50,0.86,0.91,0.91);
            leg1->SetNColumns(1);
            leg1->SetFillColor(kWhite);
            leg1->SetLineColor(kWhite);
            leg1->AddEntry(h_tau[target_spec], Form("weighted data (%s)",specs[target_spec].Data()), "ep");
            leg1->Draw();
            
            canvas->Print(Form("fig/tau_splot_%s.pdf",specs[target_spec].Data()));
            
            delete leg1;
            delete h_tau_data;
            delete frame;
            delete canvas;
        }
        
        if (sysidx>=0) { // restore the nuisance value
            TString var = cmd_tau_sys[sysidx](0,1);
            TString key = cmd_tau_sys[sysidx](1,cmd_tau_sys[sysidx].Length()-1);
            if (var!="0") {
                RooRealVar *nuisance = wspace->var(key);
                nuisance->setVal(nuisance_preserve);
            }
        }
    } // ending of systematics loop
    
    if (sys_study) {
        
        double err_p = 0., err_m = 0.;
        for (int sysidx = 0; sysidx < (int)cmd_tau_sys.size(); sysidx++) {
            TString var = cmd_tau_sys[sysidx](0,1);
            TString key = cmd_tau_sys[sysidx](1,cmd_tau_sys[sysidx].Length()-1);
            cout << ">>> SYS(" << var << ") " << key << ": " << diff_tau_sys[sysidx] << endl;
            
            if (diff_tau_sys[sysidx]>0.) err_p += diff_tau_sys[sysidx]*diff_tau_sys[sysidx];
            if (diff_tau_sys[sysidx]<0.) err_m += diff_tau_sys[sysidx]*diff_tau_sys[sysidx];
            
            delete h_tau_sys[sysidx];
        }
        err_p = sqrt(err_p);
        err_m = sqrt(err_m);
        
        RooFitResult res3(*res1);
        RooRealVar* tau1 = (RooRealVar*)res3.floatParsFinal().find("EffTau_bs");
        tau1->setError(max(err_p, err_m));
        tau1->setAsymError(-err_m, err_p);
        wspace->import(res3,"fitresult_taubs_syst"); // systematic results
        
        RooRealVar* tau0 = (RooRealVar*)res1->floatParsFinal().find("EffTau_bs");
        cout << ">>> EffTau central: " << tau0->getVal() << " +- " <<
        tau0->getError() << " (+" << tau0->getErrorHi() << "/" << tau0->getErrorLo() << ")" <<
        " (+" << err_p << "/-" << err_m << ")" << endl;
    }
    for (int idx = 0; idx < _nspec; idx++) delete h_tau[idx];
    delete res1;
    delete res2;
}

void bmm4fitter(TString commands = "")
{
    // -----------------------------------------------------------
    // parse the commands
    bool do_create_base = false;
    bool do_update_pdf = false, do_update_config = false, do_update_data = false;
    bool do_fit = false, do_minos = false, do_fullminos = false, do_bdtfit = false;
    bool do_signif = false, do_bdstat = false;
    bool do_make_plots = false;
    TString file_base = "wspace_base.root";
    TString file_prefit = "wspace_prefit.root";
    TString file_postfit = "wspace_postfit.root";
    vector<TString> set_par_names;
    vector<double> set_par_values;
    vector<int> set_par_states;
    
    cout << ">>> -------------------------" << endl;
    cout << ">>> BMM4 fitter start" << endl;
    cout << ">>> -------------------------" << endl;
    cout << ">>> commands:" << endl;
    cout << ">>> - create_base                        : reset & create the base workspace" << endl;
    cout << ">>> - update_pdf                         : update PDF (slow!)" << endl;
    cout << ">>> - update_config                      : update variables/constraints" << endl;
    cout << ">>> - update_data                        : update data set" << endl;
    cout << ">>> - fit                                : fit to data" << endl;
    cout << ">>> - minos                              : call MINOS for BFs" << endl;
    cout << ">>> - fullminos                          : call full MINOS" << endl;
    cout << ">>> - bdtfit                             : enable BDT-Fit strategy" << endl;
    cout << ">>> - signif                             : calculate significances after the fit" << endl;
    cout << ">>> - bdstat                             : calculate profile likelihood test statistics (F&C study or upper limit for Bd)" << endl;
    cout << ">>> - make_plots                         : produce projection plots" << endl;
    cout << ">>> - ws_base=[wspace_base.root]         : set base workspace" << endl;
    cout << ">>> - ws_prefit=[wspace_prefit.root]     : set prefit workspace" << endl;
    cout << ">>> - ws_postfit=[wspace_postfit.root]   : set postfit workspace" << endl;
    cout << ">>> - set_par [name]=[value]=[float/fix] : set parameter value & state" << endl;
    cout << ">>> parsing commands: [" << commands << "]" << endl;
    if (commands=="") return;
    
    TString tok;
    Ssiz_t from = 0;
    while(commands.Tokenize(tok, from, "[ \t;=:]")) {
        if (tok=="create_base")        do_create_base = true;   // reset & create the base workspace
        else if (tok=="update_pdf")    do_update_pdf = true;    // update PDF (slow!)
        else if (tok=="update_config") do_update_config = true; // update variables/constraints
        else if (tok=="update_data")   do_update_data = true;   // update data set
        
        else if (tok=="fit")           do_fit = true;           // fit to data
        else if (tok=="minos")         do_minos = true;         // call MINOS for BFs
        else if (tok=="fullminos")     do_fullminos = true;     // call full MINOS
        else if (tok=="bdtfit")        do_bdtfit = true;        // enable BDT-Fit strategy
        else if (tok=="signif")        do_signif = true;        // calculate significances after the fit
        else if (tok=="bdstat")        do_bdstat = true;        // calculate profile likelihood test statistics for Bd
        else if (tok=="make_plots")    do_make_plots = true;    // produce projection plots
        else if (tok=="ws_base") {                              // set base workspace
            commands.Tokenize(tok, from, "[ \t;=:]");
            file_base = tok;
        }else if (tok=="ws_prefit") {                           // set prefit workspace file
            commands.Tokenize(tok, from, "[ \t;=:]");
            file_prefit = tok;
        }else if (tok=="ws_postfit") {                          // set postfit workspace file
            commands.Tokenize(tok, from, "[ \t;=:]");
            file_postfit = tok;
        }else if (tok=="set_par") {                             // set parameter value & fix/float
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_par_names.push_back(tok);
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_par_values.push_back(tok.Atof());
            commands.Tokenize(tok, from, "[ \t;=:]");
            set_par_states.push_back(tok=="fix"?1:0); // 1 - fix, 0 - float
        }else {
            cout << ">>> unknown command '" << tok << "'" << endl;
            return;
        }
    }
    
    // -----------------------------------------------------------
    // Initializing categories definition
    
//    if (CONFIG_BMM3) CatMan.RegisterBMM3Categories();
    if (CONFIG_BMM4) CatMan.RegisterBMM4Categories();
    CatMan.Print();
    
    RooWorkspace *wspace_base = NULL;
    
    // -----------------------------------------------------------
    // read the base workspace
    if (file_base!="n/a" && !do_create_base) {
        cout << ">>> read the base model from '" << file_base << "'" << endl;
        exist_protection(file_base);
        TFile *fin_base = new TFile(file_base);
        wspace_base = (RooWorkspace *)fin_base->Get("wspace");
    }
    
    // global workspace
    RooWorkspace *wspace = new RooWorkspace("wspace");
    
    // define global variables & constraints
    PrepareGlobalVariables(wspace,do_update_config?0:wspace_base);
    
    // define subchannel variables & constraints
//    if (CONFIG_BMM3) PrepareBMM3SubVariables(wspace,do_update_config?0:wspace_base);
    if (CONFIG_BMM4) PrepareBMM4SubVariables(wspace,do_update_config?0:wspace_base);
    
    // define decay time related models
    PrepareLifetimeModel(wspace, do_update_pdf?0:wspace_base);
    
    // define subchannel PDF
//    if (CONFIG_BMM3) PrepareBMM3SubPDF(wspace,do_update_pdf?0:wspace_base);
    if (CONFIG_BMM4) PrepareBMM4SubPDF(wspace,do_update_pdf?0:wspace_base);
    
    // define global PDF
    PrepareGlobalPDF(wspace, do_bdtfit?"bdtfit":"bdtcat");
    
    // define list of nuisances/floats
    PrepareNamedArgSets(wspace);
    
    // prepare data
    PrepareData(wspace,do_update_data?0:wspace_base);
    
    // if base workspace reset
    if (do_create_base) {
        cout << ">>> save the base workspace to '" << file_base << "'" << endl;
        wspace->writeToFile(file_base);
    }
    
    // parameter modifications
    for (int idx = 0; idx<(int)set_par_names.size(); idx++) {
        cout << ">>> set parameter " << set_par_names[idx] << " to " << set_par_values[idx] << (set_par_states[idx]?" (fixed)":" (floated)") << endl;
        wspace->var(set_par_names[idx])->setVal(set_par_values[idx]);
        wspace->var(set_par_names[idx])->setConstant(set_par_states[idx]);
    }
    
    // full workspace prepared, save as the prefit model
    if (file_prefit!="n/a" && (do_fit || do_create_base)) {
        cout << ">>> save the pre-fit model to '" << file_prefit << "'" << endl;
        wspace->writeToFile(file_prefit);
    }
    
    // -----------------------------------------------------------
    // fit to data
    if (do_fit) {
        RooAbsData *global_data = wspace->data("global_data");
        RooAbsPdf *global_pdf = wspace->pdf("global_pdf");
        
        RooArgSet global_ext_constr(*wspace->pdf("fs_over_fu_gau"),*wspace->pdf("fs_over_fu_S13_gau"),*wspace->pdf("one_over_BRBR_gau")); // external Gaussian constraints
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
        if (do_fullminos) fit_cmd.push_back(Minos(true));
        
        RooLinkedList fit_linkl;
        for (auto& cmd : fit_cmd) fit_linkl.Add(&cmd);
        RooFitResult *res = fit_with_retry(*global_pdf,*global_data,fit_linkl);
        
        res->Print("v");
        wspace->import(*res,"fitresult");
        
        const RooArgSet *all_floats = wspace->set("all_floats");
        wspace->saveSnapshot("best_fit",*all_floats);
        
        if (do_signif) {
            // Significance estimation for Bd
            wspace->var("BF_bd")->setVal(0.);
            wspace->var("BF_bd")->setConstant(true);
            RooFitResult *res_zerobd = fit_with_retry(*global_pdf,*global_data,fit_linkl);
            wspace->import(*res_zerobd,"fitresult_zerobd");
            wspace->loadSnapshot("best_fit");
            
            // Significance estimation for Bs
            wspace->var("BF_bs")->setVal(0.);
            wspace->var("BF_bs")->setConstant(true);
            RooFitResult *res_zerobs = fit_with_retry(*global_pdf,*global_data,fit_linkl);
            wspace->import(*res_zerobs,"fitresult_zerobs");
            wspace->loadSnapshot("best_fit");
            
            double sigbs = sqrt(max(0.,res_zerobs->minNll() - res->minNll())*2.);
            double sigbd = sqrt(max(0.,res_zerobd->minNll() - res->minNll())*2.);
            
            cout << ">>> Significance for Bs: " << setprecision(3) << sigbs << " sigma." << endl;
            cout << ">>> Significance for Bd: " << setprecision(3) << sigbd << " sigma." << endl;
            
            delete res_zerobd;
            delete res_zerobs;
            
        }else if (do_bdstat) { // assume the default fit is carried out with BF_bd fixed to the target branching fraction
            // release Bd BF
            wspace->var("BF_bd")->setConstant(false);
            RooFitResult *res_floatbd = fit_with_retry(*global_pdf,*global_data,fit_linkl);
            wspace->import(*res_floatbd,"fitresult_floatbd");
            wspace->loadSnapshot("best_fit");
            
            delete res_floatbd;
        }
        
        if (file_postfit!="n/a") {
            cout << ">>> save the post-fit model and results to '" << file_postfit << "'" << endl;
            wspace->writeToFile(file_postfit);
        }
        
        delete res;
    }
    
    // -----------------------------------------------------------
    // produce the projection plots
    if (do_make_plots) {
        ProduceSubPlots(wspace);
        ProduceTauSPlot(wspace);
    }
    
    delete wspace;
    cout << ">>> BMM4 fitter end." << endl;
}

int main(int argc, char *argv[]) {
    TString cmd = "";
    for (int i=1; i<argc; i++) {
        cmd += TString(argv[i])+";";
    }
    bmm4fitter(cmd);
}
