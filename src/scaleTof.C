#include "commonFunctions.h"

TH2D* getScalingCoefficients(TH2D* hNoTof, TH2D* hTof) {
  if (!hNoTof || !hTof) {
    std::cerr << "ERROR: One or both histograms are null!" << std::endl;
    return nullptr;
  }

  if (hNoTof->GetNbinsX() != hTof->GetNbinsX() || hNoTof->GetNbinsY() != hTof->GetNbinsY()) {
    std::cerr << "ERROR: Histograms have different number of bins!" << std::endl;
    return nullptr;
  }

  // divide the two histograms (assuming binomial statistics for error propagation)
  TH2D* hCoefficients = dynamic_cast<TH2D*>(hNoTof->Clone("hCoefficients"));
  hCoefficients->Reset();
  hCoefficients->SetTitle("TOF scaling coefficients; #it{m} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c}); scaling coefficient");
  hCoefficients->Divide(hNoTof, hTof, 1, 1, "B"); // "B" option for binomial errors
  return hCoefficients;
}

TH1D* correctSpectra(TH2D* hTof, TH2D* hCoefficients) {
  if (!hTof || !hCoefficients) {
    std::cerr << "ERROR: One or both histograms are null!" << std::endl;
    return nullptr;
  }

  if (hTof->GetNbinsX() != hCoefficients->GetNbinsX() || hTof->GetNbinsY() != hCoefficients->GetNbinsY()) {
    std::cerr << "ERROR: Histograms have different number of bins!" << std::endl;
    return nullptr;
  }

  TH2D* hCorrected = dynamic_cast<TH2D*>(hTof->Clone("hCorrected"));
  hCorrected->Reset();
  hCorrected->SetTitle("TOF-corrected spectrum; #it{m} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c}); counts");
  hCorrected->Multiply(hTof, hCoefficients, 1, 1);

  // project the corrected 2D histogram into a 1D histogram along the x-axis (m) for each pT bin
  TH1D* hCorrected1D = hCorrected->ProjectionX("hCorrected1D");
  hCorrected1D->SetTitle("TOF-corrected spectrum (1D projection); #it{m} (GeV/#it{c}^{2}); counts");

  return hCorrected1D;
}

void regulariseHistogram(TH2D* hist, double threshold) {
  if (!hist) {
    std::cerr << "ERROR: Histogram is null!" << std::endl;
    return;
  }

  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    for (int j = 1; j <= hist->GetNbinsY(); ++j) {
      if (hist->GetBinContent(i, j) < threshold) {
        hist->SetBinContent(i, j, 0);
        hist->SetBinError(i, j, 0);
      }
    }
  }
}

TH1D* applyAxE(TH1D *hist, TH1D *eff) {
  TH1D *hDiv = (TH1D*)hist->Clone("hDiv");
  hDiv->Reset();
  hDiv->Divide(hist, eff, 1, 1);

  // return the divided histogram
  return hDiv;
}

void scaleTof() {
  globalSetup();
  
  TFile* corrections = TFile::Open("output/AxEcorrected.root", "READ");
  TFile* data = TFile::Open("input/data.root", "READ");
  TFile* muons = TFile::Open("input/muonReco.root", "READ");

  if (!corrections || !data || !muons) {
    std::cerr << "ERROR: Cannot open one or more input files!" << std::endl;
    return;
  }

  TFile* scaledTof = TFile::Open("output/scaledTof.root", "RECREATE");

  for (int i = 0; i < neutronClasses.size(); i++) {
    TDirectory* neutronDir = scaledTof->mkdir(neutronClasses[i].Data());
    neutronDir->cd();
    for (int j = 0; j < nBinsPt; j++) {
      TDirectory* pTDir = neutronDir->mkdir(Form("pTbin_%d", j));
      pTDir->cd();

      TH2D* hDataTof = dynamic_cast<TH2D*>(data->Get(Form("AnAn/unlike/pTbin_%d/hPtVsMwithTof", j)));
      TH2D* hDataNoTof = dynamic_cast<TH2D*>(data->Get(Form("AnAn/unlike/pTbin_%d/hPtVsM", j)));
      TH2D* hDataCoefficients = getScalingCoefficients(hDataNoTof, hDataTof);
      hDataCoefficients->Write("hDataCoefficients");

      // load 3D histograms storing events in (m, deltaPhi, pT) space
      TH3D* hData3D = dynamic_cast<TH3D*>(data->Get(Form("%s/unlike/pTbin_%d/hDeltaPhiVsPtVsMwithTof", neutronClasses[i].Data(), j)));
      if (!hData3D) {
        std::cerr << "ERROR: Cannot find 3D histogram for neutron class " << neutronClasses[i].Data() << " and pT bin " << j << std::endl;
        continue;
      }

      for (int k = 1; k <= nBinsDeltaPhi; ++k) {
        TDirectory* phiDir = pTDir->mkdir(Form("phibin_%d", k));
        phiDir->cd();

        // project the 3D histogram into the (m, pT) plane for the current deltaPhi bin
        hData3D->GetZaxis()->SetRange(k, k);
        TH2D* hDataTofSlice = dynamic_cast<TH2D*>(hData3D->Project3D("yxe"));
        regulariseHistogram(hDataTofSlice, 25.0); // regularise the histogram to avoid low-statistics bins
        hDataTofSlice->Write("hDataTofSlice");

        // clone the projection and correct the spectrum using the coefficients
        TH1D* hDataTofCorrected = correctSpectra(hDataTofSlice, hDataCoefficients);
        hDataTofCorrected->Write("hDataTofCorrected");

        TH1D* hAxE = dynamic_cast<TH1D*>(corrections->Get(Form("%s/pTbin_%d/phibin_%d/efficiency", neutronClasses[i].Data(), j, k)));
        if (!hAxE) {
          std::cerr << "ERROR: Cannot find AxE histogram for neutron class " << neutronClasses[i].Data() << ", pT bin " << j << ", deltaPhi bin " << k << std::endl;
          continue;
        }
        hAxE->Write("hAxE");

        TH1D* hDataTofFinal = applyAxE(hDataTofCorrected, hAxE);
        hDataTofFinal->Write("correctedSpectrum");
      }
    }
  }
  scaledTof->Close();
}