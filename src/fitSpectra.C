#include "commonFunctions.h"

void fitSpectra() {
  globalSetup();
  ROOT::EnableImplicitMT();

  TFile *fIn = new TFile("output/AxEcorrected.root", "READ");

  TFile *fFits = new TFile("output/fitResults.root", "RECREATE");

  // load omega phase parameter from the previous global fit
  TFile *fGlobalFit = new TFile("output/omegaPhase.root", "READ");
  TH1D *hOmegaPhase = (TH1D*) fGlobalFit->Get("hOmegaPhase");
  double omegaPhase = hOmegaPhase->GetBinContent(1);

  // fit function for muons
  TF1 *fitFuncMuons = new TF1("fitFuncMuons", "[0]*x^(-[1])", fitMin, fitMax);
  fitFuncMuons->SetParNames("#it{a}_{#mu}", "#it{b}_{#mu}");
  TH1D *hFitParamMuons_a = new TH1D("hFitParamMuons_a", Form("hFitParamMuons_a;%s (rad);#it{a}_{#mu}", deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());
  TH1D *hFitParamMuons_b = new TH1D("hFitParamMuons_b", Form("hFitParamMuons_b;%s (rad);#it{b}_{#mu}", deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());
  TH1D *hFitParamMuons_chi2 = new TH1D("hFitParamMuons_chi2", Form("hFitParamMuons_chi2;%s (rad);#it{#chi^{2}}/ndf", deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());

  // fit function for soeding
  TF1 *fitFuncSoeding = new TF1("fitFuncSoeding", Soeding, fitMin, fitMax, 10);
  fitFuncSoeding->SetParNames("#it{A}_{#rho}", "#it{m}_{#rho}", "#it{#Gamma}_{#rho}", "#it{C}_{#omega}", "#it{m}_{#omega}", "#it{#Gamma}_{#omega}", "#it{#phi}_{#omega}", "#it{B}_{#pi#pi}", "#it{a}_{#mu}", "#it{b}_{#mu}");
  TString paramNames[11] = {"#it{A}_{#rho}", "#it{m}_{#rho}", "#it{#Gamma}_{#rho}", "#it{C}_{#omega}", "#it{m}_{#omega}", "#it{#Gamma}_{#omega}", "#it{#phi}_{#omega}", "#it{B}_{#pi#pi}", "#it{a}_{#mu}", "#it{b}_{#mu}", "#it{#chi}^{2}/ndf"};
  TH1D *hFitParamSoeding[11];
  for (int i = 0; i < 10; i++) hFitParamSoeding[i] = new TH1D(Form("hFitParamSoeding_%d", i), Form("hFitParamSoeding_%d;%s (rad);%s", i, deltaPhiLabel.c_str(), paramNames[i].Data()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());
  hFitParamSoeding[10] = new TH1D("hFitParamSoeding_chi2", Form("hFitParamSoeding_chi2;%s (rad);#it{#chi^{2}}/ndf", deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());

  // functions for individual fit function contributions
  TF1 *fitFuncRho = new TF1("fitFuncRho", rhoContrib, fitMin, fitMax, 10);
  TF1 *fitFuncOmega = new TF1("fitFuncOmega", omegaContrib, fitMin, fitMax, 10);
  TF1 *fitFuncDiPion = new TF1("fitFuncDiPion", diPionContrib, fitMin, fitMax, 10);
  TF1 *fitFuncRhoDiPionInt = new TF1("fitFuncRhoDiPionInt", rhoDiPionInterference, fitMin, fitMax, 10);
  TF1 *fitFuncRhoOmegaInt = new TF1("fitFuncRhoOmegaInt", rhoOmegaInterference, fitMin, fitMax, 10);
  TF1 *fitFuncOmegaDiPionInt = new TF1("fitFuncOmegaDiPionInt", omegaDiPionInterference, fitMin, fitMax, 10);
  TF1 *fitFuncBackground = new TF1("fitFuncBackground", powerBackground, fitMin, fitMax, 10);

  // histogram to store integrated rho yield
  TH1D *hRhoYield = new TH1D("hRhoYield", Form("hRhoYield;%s (rad);corrected #rho^{0} yield", deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());
  // histos to store yield fit parameters
  TH1D *hRhoYieldFitParam[5];
  TString yieldParamNames[5] = {"c", "#it{a}_{1}", "#it{a}_{2}", "#it{a}_{3}", "#it{a}_{4}"};
  for (int i = 0; i < 5; i++) hRhoYieldFitParam[i] = new TH1D(Form("hRhoYieldFitParam_%d", i), Form("hRhoYieldFitParam_%d;%s (rad);%s", i, deltaPhiLabel.c_str(), yieldParamNames[i].Data()), 1, -TMath::Pi(), TMath::Pi());

  // loop over neutron classes, pT bins and phi bins
  for (int i = 0; i < neutronClasses.size(); i++) {
    for (int j = 0; j < nBinsPt; j++) {
      for (int k = 1; k <= nBinsDeltaPhi; ++k) {
        // calculate Delta phi range corresponding to the phi bin
        const double phiBinWidth = (2.0 * TMath::Pi()) / nBinsDeltaPhi;
        const double phiBinMin = -TMath::Pi() + (k - 1) * phiBinWidth;
        const double phiBinMax = phiBinMin + phiBinWidth;
        // i would like to print this into the plot in multiples of pi, to make obvious that the range is really -pi to pi
        const TString phiBinLabel = Form("%s #in [%.2f #pi, %.2f #pi) rad", deltaPhiLabel.c_str(), phiBinMin/TMath::Pi(), phiBinMax/TMath::Pi());

        // get muon histogram and fit it with the muon fit funtion
        TH1D *hMuon = (TH1D*) fIn->Get(Form("%s/pTbin_%d/phibin_%d/correctedMuonSpectrum", neutronClasses[i].Data(), j, k));
        fitFuncMuons->SetParameters(1.0, -1.0);
        TFitResultPtr muonFit = hMuon->Fit(fitFuncMuons, "EMR0S", "", fitMin, fitMax);
        if (!muonFit) std::cerr << "Error: Muon fit failed for " << neutronClasses[i].Data() << " pT bin " << j << " phi bin " << k << std::endl;
        double *muonFitParameters = fitFuncMuons->GetParameters();
        const TMatrixDSym& muonCovMatrix = muonFit->GetCovarianceMatrix();
        if (muonCovMatrix.GetNrows() == 0) std::cerr << "Covariance matrix is empty for muon fit " << neutronClasses[i].Data() << " pT bin " << j << " phi bin " << k << std::endl;
        // store fit parameters
        hFitParamMuons_a->SetBinContent(k, fitFuncMuons->GetParameter(0));
        hFitParamMuons_a->SetBinError(k, fitFuncMuons->GetParError(0));
        hFitParamMuons_b->SetBinContent(k, fitFuncMuons->GetParameter(1));
        hFitParamMuons_b->SetBinError(k, fitFuncMuons->GetParError(1));
        hFitParamMuons_chi2->SetBinContent(k, fitFuncMuons->GetChisquare()/fitFuncMuons->GetNDF());
        // plot results
        TCanvas *cMu = new TCanvas(Form("cMu_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k), Form("cMu_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k), 1000, 600);
        cMu->SetLogx();
        cMu->SetLogy();
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);

        hMuon->SetLineColorAlpha(tabBlue, 1.0);
        hMuon->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
        hMuon->GetXaxis()->SetTitleSize(0.04);
        hMuon->GetYaxis()->SetTitleSize(0.04);
        hMuon->GetXaxis()->SetLabelSize(0.04);
        hMuon->GetYaxis()->SetLabelSize(0.04);
        hMuon->GetYaxis()->SetTitle(Form("corrected counts per %.3f GeV/ #it{c}^{2}", hMuon->GetBinWidth(1)));
        hMuon->SetMarkerColor(hMuon->GetLineColor());
        hMuon->SetLineWidth(2);
        hMuon->Draw("e");

        fitFuncMuons->SetLineColorAlpha(tabOrange, 1.0);
        fitFuncMuons->SetNpx(1000);
        fitFuncMuons->SetLineWidth(2);
        fitFuncMuons->Draw("same");
        hMuon->Draw("e same");

        TLegend* muonStats = new TLegend(0.5, 0.6, 0.9, 0.89);
        muonStats->SetTextAlign(32);
        muonStats->SetBorderSize(0);
        muonStats->SetFillStyle(0);
        muonStats->AddEntry((TObject*)0, "LHC25g15 kTwoGammaToMuLow", "");
        muonStats->AddEntry((TObject*)0, Form("misidentified #gamma#gamma #rightarrow #mu^{+}#mu^{-}, %s", neutronClasses[i].Data()), "");
        muonStats->AddEntry((TObject*)0, Form("%.2f #leq #it{p}_{T} #leq %.2f GeV/#it{c}, | #it{y}| #leq %.1f", pTbinEdges[j], pTbinEdges[j+1], maxY), "");
        muonStats->AddEntry((TObject*)0, phiBinLabel, "");
        // muonStats->AddEntry((TObject*)0, Form("%s #in (%.2f, %.2f) rad", deltaPhiLabel.c_str(), phiBinMin, phiBinMax), "");
        muonStats->AddEntry((TObject*)0, Form("#it{a}_{#mu} = %.0f #pm %.0f", fitFuncMuons->GetParameter(0), fitFuncMuons->GetParError(0)), "");
        muonStats->AddEntry((TObject*)0, Form("#it{b}_{#mu} = %.1f #pm %.1f", fitFuncMuons->GetParameter(1), fitFuncMuons->GetParError(1)), "");
        muonStats->AddEntry((TObject*)0, Form("#it{#chi}^{2}/ndf = %.1f", fitFuncMuons->GetChisquare()/fitFuncMuons->GetNDF()), "");
        muonStats->Draw();

        gSystem->Exec(Form("mkdir -p output/fitResults/%s/pTbin_%d/phiBin_%d", neutronClasses[i].Data(), j, k));
        cMu->SaveAs(Form("output/fitResults/%s/pTbin_%d/phiBin_%d/hMuon.pdf", neutronClasses[i].Data(), j, k));
        delete cMu;
        delete muonStats;

        // plot correlation matrix as a 2D histogram
        TCanvas *cMuCorr = new TCanvas(Form("cMuCorr_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k), Form("cMuCorr_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k), 800, 600);
        TMatrixDSym muonCorrMatrix = muonFit->GetCorrelationMatrix();
        int nMuPar = muonCorrMatrix.GetNrows();
        TH2D *hMuCorr = new TH2D("hMuCorr", ";;", nMuPar, 0, nMuPar, nMuPar, 0, nMuPar);
        for (int m = 0; m < nMuPar; ++m) {
          hMuCorr->GetXaxis()->SetBinLabel(m + 1, fitFuncMuons->GetParName(m));
          hMuCorr->GetXaxis()->SetLabelSize(0.05);
          hMuCorr->GetYaxis()->SetBinLabel(m + 1, fitFuncMuons->GetParName(m));
          hMuCorr->GetYaxis()->SetLabelSize(0.05);
          for (int n = 0; n < nMuPar; ++n) {
            hMuCorr->SetBinContent(m + 1, n + 1, muonCorrMatrix(m, n));
          }
        }
        hMuCorr->SetMarkerColor(kWhite);
        gStyle->SetPaintTextFormat(".2f");
        hMuCorr->GetZaxis()->SetRangeUser(-1, 1);
        hMuCorr->Draw("colz1 text");
        cMuCorr->SaveAs(Form("output/fitResults/%s/pTbin_%d/phiBin_%d/hMuonCorrelation.pdf", neutronClasses[i].Data(), j, k));
        delete cMuCorr;
        delete hMuon;
        delete hMuCorr;

        // get soeding histogram and fit it with the soeding fit function
        TH1D *hData = (TH1D*) fIn->Get(Form("%s/pTbin_%d/phibin_%d/correctedSpectrum", neutronClasses[i].Data(), j, k));
        // set initial parameters for soeding fit depending on the neutron class
        if (neutronClasses[i] == "AnAn") fitFuncSoeding->SetParameters(500.0, kMrho, kWrho, -10.0, kMomega, kWomega, omegaPhase, -315.0, 10e3);
        if (neutronClasses[i] == "0n0n") fitFuncSoeding->SetParameters(450.0, kMrho, kWrho, -10.0, kMomega, kWomega, omegaPhase, -310.0, 8e3);
        if (neutronClasses[i] == "Xn0n" || neutronClasses[i] == "0nXn") fitFuncSoeding->SetParameters(110.0, kMrho, kWrho, -2.0, kMomega, kWomega, omegaPhase, -55.0, 500.0);
        if (neutronClasses[i] == "XnXn") fitFuncSoeding->SetParameters(55.0, kMrho, kWrho, 0.0, kMomega, kWomega, omegaPhase, -17.0, 150.0);
        fitFuncSoeding->SetParLimits(3, -50.0, 0.0);
        // if (neutronClasses[i] == "XnXn") fitFuncSoeding->FixParameter(3, 0.0); // no omega contribution in XnXn class
        fitFuncSoeding->SetParLimits(6, -TMath::Pi(), 0.0); 
        fitFuncSoeding->SetParLimits(7, -1000.0, 0.0);
        fitFuncSoeding->FixParameter(9, fitFuncMuons->GetParameter(1)); // muon b
        if (fixPoles) {
          fitFuncSoeding->FixParameter(1, kMrho);
          fitFuncSoeding->FixParameter(2, kWrho);
          fitFuncSoeding->FixParameter(4, kMomega);
          fitFuncSoeding->FixParameter(5, kWomega);
        }
        fitFuncSoeding->FixParameter(6, omegaPhase); // fix omega phase to the value obtained from the global fit
        TFitResultPtr fitResult = hData->Fit(fitFuncSoeding, "EMR0S", "", fitMin, fitMax);
        double *fitParameters = fitFuncSoeding->GetParameters();
        const TMatrixDSym& covMatrix = fitResult->GetCovarianceMatrix();
        // store fit parameters
        for (int p = 0; p < fitFuncSoeding->GetNpar(); p++) {
          hFitParamSoeding[p]->SetBinContent(k, fitFuncSoeding->GetParameter(p));
          hFitParamSoeding[p]->SetBinError(k, fitFuncSoeding->GetParError(p));
        }
        hFitParamSoeding[10]->SetBinContent(k, fitFuncSoeding->GetChisquare()/fitFuncSoeding->GetNDF());

        // plot results
        TCanvas *cSoeding = new TCanvas(Form("cSoeding_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k), Form("cSoeding_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k), 1000, 600);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
        
        hData->SetLineColorAlpha(tabBlue, 1.0);
        hData->SetMarkerColor(hData->GetLineColor());
        hData->SetLineWidth(2);
        hData->SetMinimum(-hData->GetMaximum()*0.25);
        hData->GetXaxis()->SetTitle("#it{m}_{#pi#pi} (GeV/#it{c}^{2})");
        hData->GetXaxis()->SetTitleSize(0.04);
        hData->GetYaxis()->SetTitleSize(0.04);
        hData->GetXaxis()->SetLabelSize(0.04);
        hData->GetYaxis()->SetLabelSize(0.04);
        hData->GetYaxis()->SetTitle(Form("corrected counts per %.2f GeV/ #it{c}^{2}", hData->GetBinWidth(1)));
        // hData->GetYaxis()->SetTitle("corrected counts");
        hData->Draw("e");
        
        fitFuncSoeding->SetLineColorAlpha(tabOrange, 1.0);
        fitFuncSoeding->SetNpx(1000);
        fitFuncSoeding->SetLineWidth(2);
        fitFuncSoeding->Draw("same");
        
        fitFuncRho->SetParameters(fitFuncSoeding->GetParameters());
        fitFuncRho->SetLineColorAlpha(tabGreen, 1.0);
        fitFuncRho->SetNpx(1000);
        fitFuncRho->SetLineWidth(2);
        fitFuncRho->Draw("same");
        // get integral of rho contribution to extract the yield
        double rhoYield = fitFuncRho->Integral(fitMin, fitMax);
        double rhoYieldError = fitFuncRho->IntegralError(fitMin, fitMax, fitParameters, covMatrix.GetMatrixArray());
        hRhoYield->SetBinContent(k, rhoYield);
        hRhoYield->SetBinError(k, rhoYieldError);
        
        hData->Draw("same");

        fitFuncOmega->SetParameters(fitFuncSoeding->GetParameters());
        fitFuncOmega->SetLineColorAlpha(tabRed, 1.0);
        fitFuncOmega->SetNpx(1000);
        fitFuncOmega->SetLineWidth(2);
        fitFuncOmega->Draw("same");

        fitFuncDiPion->SetParameters(fitFuncSoeding->GetParameters());
        fitFuncDiPion->SetLineColorAlpha(tabPurple, 1.0);
        fitFuncDiPion->SetNpx(1000);
        fitFuncDiPion->SetLineWidth(2);
        fitFuncDiPion->Draw("same");

        fitFuncRhoDiPionInt->SetParameters(fitFuncSoeding->GetParameters());
        fitFuncRhoDiPionInt->SetLineColorAlpha(tabCyan, 1.0);
        fitFuncRhoDiPionInt->SetNpx(1000);
        fitFuncRhoDiPionInt->SetLineWidth(2);
        fitFuncRhoDiPionInt->Draw("same");

        fitFuncRhoOmegaInt->SetParameters(fitFuncSoeding->GetParameters());
        fitFuncRhoOmegaInt->SetLineColorAlpha(tabPink, 1.0);
        fitFuncRhoOmegaInt->SetNpx(1000);
        fitFuncRhoOmegaInt->SetLineWidth(2);
        fitFuncRhoOmegaInt->Draw("same");

        fitFuncOmegaDiPionInt->SetParameters(fitFuncSoeding->GetParameters());
        fitFuncOmegaDiPionInt->SetLineColorAlpha(tabOlive, 0.5);
        fitFuncOmegaDiPionInt->SetNpx(1000);
        fitFuncOmegaDiPionInt->SetLineWidth(2);
        fitFuncOmegaDiPionInt->Draw("same");

        fitFuncBackground->SetParameters(fitFuncSoeding->GetParameters());
        fitFuncBackground->SetLineColorAlpha(tabBrown, 1.0);
        fitFuncBackground->SetNpx(1000);
        fitFuncBackground->SetLineWidth(2);
        fitFuncBackground->Draw("same");

        hData->Draw("e same");

        TLegend* soedingStats = new TLegend(0.5, 0.45, 0.9, 0.89);
        soedingStats->SetTextAlign(32);
        soedingStats->SetBorderSize(0);
        soedingStats->SetFillStyle(0);
        soedingStats->AddEntry((TObject*)0, "ALICE Pb#font[122]{-}Pb UPC #sqrt{s_{NN}} = 5.36 TeV", "");
        // soedingStats->AddEntry((TObject*)0, "UD_LHC23_pass5_SingleGap, this thesis", "");
        soedingStats->AddEntry((TObject*)0, Form("#rho^{0} #rightarrow #pi^{+}#pi^{-}, %s", neutronClasses[i].Data()), "");
        soedingStats->AddEntry((TObject*)0, Form("%.2f #leq #it{p}_{T} #leq %.2f GeV/#it{c}, | #it{y}| #leq %.1f", pTbinEdges[j], pTbinEdges[j+1], maxY), "");
        soedingStats->AddEntry((TObject*)0, phiBinLabel, "");
        for (int p = 0; p < fitFuncSoeding->GetNpar(); p++) {
          if (fixPoles && (p == 1 || p == 2 || p == 4 || p == 5)) continue; // skip fixed parameters
          if (p == 9) continue; // skip muon b parameter
          soedingStats->AddEntry((TObject*)0, Form("%s = %.2f #pm %.2f", fitFuncSoeding->GetParName(p), fitFuncSoeding->GetParameter(p), fitFuncSoeding->GetParError(p)), "");
        }
        soedingStats->AddEntry((TObject*)0, Form("#it{#chi}^{2}/ndf = %.1f", fitFuncSoeding->GetChisquare()/fitFuncSoeding->GetNDF()), "");

        soedingStats->Draw();

        TLegend *soedingLegend = new TLegend(0.125, 0.55, 0.35, 0.89);
        soedingLegend->SetBorderSize(0);
        soedingLegend->SetFillStyle(0);
        soedingLegend->AddEntry(fitFuncSoeding, "Soeding fit", "l");
        soedingLegend->AddEntry(fitFuncRho, "#rho^{0}", "l");
        soedingLegend->AddEntry(fitFuncOmega, "#omega^{0}", "l");
        soedingLegend->AddEntry(fitFuncDiPion, "direct #pi^{+}#pi^{-}", "l");
        soedingLegend->AddEntry(fitFuncRhoDiPionInt, "#rho^{0}-#pi^{+}#pi^{-}", "l");
        soedingLegend->AddEntry(fitFuncRhoOmegaInt, "#rho^{0}-#omega^{0}", "l");
        soedingLegend->AddEntry(fitFuncOmegaDiPionInt, "#omega^{0}-#pi^{+}#pi^{-}", "l");
        soedingLegend->AddEntry(fitFuncBackground, "#mu^{+}#mu^{-} background", "l");
        soedingLegend->Draw();

        gSystem->Exec(Form("mkdir -p output/fitResults/%s/pTbin_%d/phiBin_%d", neutronClasses[i].Data(), j, k));
        cSoeding->SaveAs(Form("output/fitResults/%s/pTbin_%d/phiBin_%d/hSoeding.pdf", neutronClasses[i].Data(), j, k));
        delete cSoeding;
        delete soedingStats;
        delete soedingLegend;

        // plot correlation matrix as a 2D histogram
        TCanvas *cSoedingCorr = new TCanvas(Form("cSoedingCorr_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k), Form("cSoedingCorr_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k), 800, 600);
        TMatrixDSym soedingCorrMatrix = fitResult->GetCorrelationMatrix();
        int nSoedingPar = soedingCorrMatrix.GetNrows();
        TH2D *hSoedingCorr = new TH2D("hSoedingCorr", ";;", nSoedingPar, 0, nSoedingPar, nSoedingPar, 0, nSoedingPar);
        for (int m = 0; m < nSoedingPar; ++m) {
          hSoedingCorr->GetXaxis()->SetBinLabel(m + 1, fitFuncSoeding->GetParName(m));
          hSoedingCorr->GetXaxis()->SetLabelSize(0.05);
          hSoedingCorr->GetYaxis()->SetBinLabel(m + 1, fitFuncSoeding->GetParName(m));
          hSoedingCorr->GetYaxis()->SetLabelSize(0.05);
          for (int n = 0; n < nSoedingPar; ++n) {
            hSoedingCorr->SetBinContent(m + 1, n + 1, soedingCorrMatrix(m, n));
          }
        }
        hSoedingCorr->SetMarkerColor(kWhite);
        gStyle->SetPaintTextFormat(".2f");
        hSoedingCorr->GetZaxis()->SetRangeUser(-1, 1);
        hSoedingCorr->Draw("colz1 text");
        cSoedingCorr->SaveAs(Form("output/fitResults/%s/pTbin_%d/phiBin_%d/hSoedingCorrelation.pdf", neutronClasses[i].Data(), j, k));
        delete cSoedingCorr;
        delete hData;
        delete hSoedingCorr;
      }

      // fit yield histogram with 4th order Fourier decomposition
      TF1 *fitFuncYield = new TF1("fitFuncYield", fourierDecomp, -TMath::Pi(), TMath::Pi(), 5);
      fitFuncYield->SetParNames("c", "#it{a}_{1}", "#it{a}_{2}", "#it{a}_{3}", "#it{a}_{4}");
      fitFuncYield->SetParameters(fitFuncYield->GetMaximum(), 0.05, 0.03, 0.01, 0.0);
      fitFuncYield->SetParLimits(1, 0.0, 1.0);
      fitFuncYield->SetParLimits(2, 0.0, 1.0);
      fitFuncYield->SetParLimits(3, 0.0, 0.2);
      fitFuncYield->SetParLimits(4, 0.0, 0.1);
      TFitResultPtr yieldFit = hRhoYield->Fit(fitFuncYield, "EMR0S", "", -TMath::Pi(), TMath::Pi());
      if (!yieldFit) std::cerr << "Error: Yield fit failed for " << neutronClasses[i].Data() << " pT bin " << j << std::endl;

      // plot yield histogram and fit result
      TCanvas *cRhoYield = new TCanvas(Form("cYield_%s_pTbin_%d", neutronClasses[i].Data(), j), Form("cYield_%s_pTbin_%d", neutronClasses[i].Data(), j), 1000, 600);
      hRhoYield->SetLineColorAlpha(tabBlue, 1.0);
      hRhoYield->SetMarkerColor(hRhoYield->GetLineColor());
      hRhoYield->SetLineWidth(2);
      hRhoYield->GetXaxis()->SetTitle(Form("%s (rad)", deltaPhiLabel.c_str()));
      hRhoYield->GetXaxis()->SetTitleSize(0.04);
      hRhoYield->GetYaxis()->SetTitleSize(0.04);
      hRhoYield->GetXaxis()->SetLabelSize(0.04);
      hRhoYield->GetYaxis()->SetLabelSize(0.04);
      hRhoYield->GetYaxis()->SetTitle("corrected #rho^{0} yield");
      // hRhoYield->SetMaximum(hRhoYield->GetMaximum()*1.05);
      hRhoYield->Draw("e");
      fitFuncYield->SetLineColorAlpha(tabOrange, 1.0);
      fitFuncYield->SetNpx(1000);
      fitFuncYield->SetLineWidth(2);
      fitFuncYield->Draw("same");
      hRhoYield->Draw("e same");
      TLegend* yieldStats = new TLegend(0.5, 0.6, 0.9, 0.9);
      yieldStats->SetTextAlign(32);
      yieldStats->SetBorderSize(0);
      yieldStats->SetFillStyle(0);
      yieldStats->AddEntry((TObject*)0, "ALICE Pb#font[122]{-}Pb UPC #sqrt{s_{NN}} = 5.36 TeV", "");
      // yieldStats->AddEntry((TObject*)0, "UD_LHC23_pass5_SingleGap, this thesis", "");
      yieldStats->AddEntry((TObject*)0, Form("#rho^{0} #rightarrow #pi^{+}#pi^{-}, %s", neutronClasses[i].Data()), "");
      yieldStats->AddEntry((TObject*)0, Form("%.2f #leq #it{p}_{T} #leq %.2f GeV/#it{c}, | #it{y}| #leq %.1f", pTbinEdges[j], pTbinEdges[j+1], maxY), "");
      for (int p = 0; p < 5; p++) {
        if (p == 0) continue; // skip constant parameter
        yieldStats->AddEntry((TObject*)0, Form("%s = %.3f #pm %.3f", fitFuncYield->GetParName(p), fitFuncYield->GetParameter(p), fitFuncYield->GetParError(p)), "");
      }
      yieldStats->AddEntry((TObject*)0, Form("#it{#chi}^{2}/ndf = %.1f", fitFuncYield->GetChisquare()/fitFuncYield->GetNDF()), "");
      yieldStats->Draw();
      gSystem->Exec(Form("mkdir -p output/fitResults/%s/pTbin_%d", neutronClasses[i].Data(), j));
      cRhoYield->SaveAs(Form("output/fitResults/%s/pTbin_%d/hRhoYield.pdf", neutronClasses[i].Data(), j));
      // store yield fit parameters
      for (int p = 0; p < 5; p++) {
        hRhoYieldFitParam[p]->SetBinContent(1, fitFuncYield->GetParameter(p));
        hRhoYieldFitParam[p]->SetBinError(1, fitFuncYield->GetParError(p));
      }
      delete cRhoYield;
      delete yieldStats;

      // plot the correlation matrix of the yield fit as a 2D histogram
      TCanvas *cYieldCorr = new TCanvas(Form("cYieldCorr_%s_pTbin_%d", neutronClasses[i].Data(), j), Form("cYieldCorr_%s_pTbin_%d", neutronClasses[i].Data(), j), 800, 600);
      TMatrixDSym yieldCorrMatrix = yieldFit->GetCorrelationMatrix();
      int nYieldPar = yieldCorrMatrix.GetNrows();
      TH2D *hYieldCorr = new TH2D("hYieldCorr", ";;", nYieldPar, 0, nYieldPar, nYieldPar, 0, nYieldPar);
      for (int m = 0; m < nYieldPar; ++m) {
        hYieldCorr->GetXaxis()->SetBinLabel(m + 1, fitFuncYield->GetParName(m));
        hYieldCorr->GetXaxis()->SetLabelSize(0.05);
        hYieldCorr->GetYaxis()->SetBinLabel(m + 1, fitFuncYield->GetParName(m));
        hYieldCorr->GetYaxis()->SetLabelSize(0.05);
        for (int n = 0; n < nYieldPar; ++n) {
          hYieldCorr->SetBinContent(m + 1, n + 1, yieldCorrMatrix(m, n));
        }
      }
      hYieldCorr->SetMarkerColor(kWhite);
      gStyle->SetPaintTextFormat(".2f");
      hYieldCorr->GetZaxis()->SetRangeUser(-1, 1);
      hYieldCorr->Draw("colz1 text");
      cYieldCorr->SaveAs(Form("output/fitResults/%s/pTbin_%d/hRhoYieldCorrelation.pdf", neutronClasses[i].Data(), j));
      delete cYieldCorr;

      // create B/A and C/A ratio histograms
      TH1D *hBoverA = new TH1D("hBoverA", Form("hBoverA;%s (rad);|#it{B}/#it{A}|", deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());
      TH1D *hCoverA = new TH1D("hCoverA", Form("hCoverA;%s (rad);|#it{C}/#it{A}|", deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());
      for (int k = 1; k <= nBinsDeltaPhi; ++k) {
        double A = hFitParamSoeding[0]->GetBinContent(k);
        double B = hFitParamSoeding[7]->GetBinContent(k);
        double C = hFitParamSoeding[3]->GetBinContent(k);
        double Aerr = hFitParamSoeding[0]->GetBinError(k);
        double Berr = hFitParamSoeding[7]->GetBinError(k);
        double Cerr = hFitParamSoeding[3]->GetBinError(k);
        hBoverA->SetBinContent(k, std::abs(B/A));
        hCoverA->SetBinContent(k, std::abs(C/A));
        // error propagation for ratio
        if (A != 0) {
          hBoverA->SetBinError(k, std::abs(B/A)*std::sqrt(std::pow(Berr/B, 2) + std::pow(Aerr/A, 2)));
          hCoverA->SetBinError(k, std::abs(C/A)*std::sqrt(std::pow(Cerr/C, 2) + std::pow(Aerr/A, 2)));
        } else {
          hBoverA->SetBinError(k, 0.0);
          hCoverA->SetBinError(k, 0.0);
        }
      }

      TLegend *ratioLegend = new TLegend(0.5, 0.7, 0.9, 0.9);
      ratioLegend->SetTextAlign(32);
      ratioLegend->SetBorderSize(0);
      ratioLegend->SetFillStyle(0);
      ratioLegend->AddEntry((TObject*)0, "ALICE Pb#font[122]{-}Pb UPC #sqrt{s_{NN}} = 5.36 TeV", "");
      // ratioLegend->AddEntry((TObject*)0, "UD_LHC23_pass5_SingleGap, this thesis", "");
      ratioLegend->AddEntry((TObject*)0, Form("#rho^{0} #rightarrow #pi^{+}#pi^{-}, %s", neutronClasses[i].Data()), "");
      ratioLegend->AddEntry((TObject*)0, Form("%.2f #leq #it{p}_{T} #leq %.2f GeV/#it{c}, | #it{y}| #leq %.1f", pTbinEdges[j], pTbinEdges[j+1], maxY), "");
      
      TCanvas *cBoverA = new TCanvas(Form("cBoverA_%s_pTbin_%d", neutronClasses[i].Data(), j), Form("cBoverA_%s_pTbin_%d", neutronClasses[i].Data(), j), 1000, 600);
      hBoverA->SetLineColorAlpha(tabBlue, 1.0);
      hBoverA->SetMarkerColor(hBoverA->GetLineColor());
      hBoverA->SetLineWidth(2);
      // hBoverA->GetXaxis()->SetTitleSize(0.04);
      // hBoverA->GetYaxis()->SetTitleSize(0.04);
      // hBoverA->GetXaxis()->SetLabelSize(0.04);
      // hBoverA->GetYaxis()->SetLabelSize(0.04);
      hBoverA->SetMaximum(hBoverA->GetMaximum()*1.05);
      hBoverA->Draw("e");
      ratioLegend->Draw();
      cBoverA->SaveAs(Form("output/fitResults/%s/pTbin_%d/hBoverA.pdf", neutronClasses[i].Data(), j));
      delete cBoverA;

      TCanvas *cCoverA = new TCanvas(Form("cCoverA_%s_pTbin_%d", neutronClasses[i].Data(), j), Form("cCoverA_%s_pTbin_%d", neutronClasses[i].Data(), j), 1000, 600);
      hCoverA->SetLineColorAlpha(tabBlue, 1.0);
      hCoverA->SetMarkerColor(hCoverA->GetLineColor());
      hCoverA->SetLineWidth(2);
      // hCoverA->GetXaxis()->SetTitleSize(0.04);
      // hCoverA->GetYaxis()->SetTitleSize(0.04);
      // hCoverA->GetXaxis()->SetLabelSize(0.04);
      // hCoverA->GetYaxis()->SetLabelSize(0.04);
      hCoverA->SetMaximum(hCoverA->GetMaximum()*1.25);
      hCoverA->Draw("e");
      ratioLegend->Draw();
      cCoverA->SaveAs(Form("output/fitResults/%s/pTbin_%d/hCoverA.pdf", neutronClasses[i].Data(), j));
      delete cCoverA;
      delete ratioLegend;

      TDirectory *dir = fFits->mkdir(Form("%s/pTbin_%d", neutronClasses[i].Data(), j));
      dir->cd();
      hFitParamMuons_a->Write();
      hFitParamMuons_b->Write();
      hFitParamMuons_chi2->Write();
      for (int p = 0; p < 11; p++) hFitParamSoeding[p]->Write();
      hRhoYield->Write();
      for (int p = 0; p < 5; p++) hRhoYieldFitParam[p]->Write();
      hBoverA->Write();
      hCoverA->Write();
    }
  }
}