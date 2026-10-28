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
      
      // TH2D* hMuonTof = dynamic_cast<TH2D*>(muons->Get(Form("AnAn/unlike/pTbin_%d/hPtVsMwithTof", j)));
      // TH2D* hMuonNoTof = dynamic_cast<TH2D*>(muons->Get(Form("AnAn/unlike/pTbin_%d/hPtVsM", j)));
      // TH2D* hMuonCoefficients = getScalingCoefficients(hMuonNoTof, hMuonTof);
      // hMuonCoefficients->Write("hMuonCoefficients");

      // load 3D histograms storing events in (m, deltaPhi, pT) space
      TH3D* hData3D = dynamic_cast<TH3D*>(data->Get(Form("%s/unlike/pTbin_%d/hDeltaPhiVsPtVsMwithTof", neutronClasses[i].Data(), j)));
      if (!hData3D) {
        std::cerr << "ERROR: Cannot find 3D histogram for neutron class " << neutronClasses[i].Data() << " and pT bin " << j << std::endl;
        continue;
      }

      // same for muons
      // TH3D* hMuon3D = dynamic_cast<TH3D*>(muons->Get(Form("AnAn/unlike/pTbin_%d/hDeltaPhiVsPtVsMwithTof", j)));
      // if (!hMuon3D) {
      //   std::cerr << "ERROR: Cannot find 3D histogram for muons and pT bin " << j << std::endl;
      //   continue;
      // }

      for (int k = 1; k <= nBinsDeltaPhi; ++k) {
        TDirectory* phiDir = pTDir->mkdir(Form("phibin_%d", k));
        phiDir->cd();

        // project the 3D histogram into the (m, pT) plane for the current deltaPhi bin
        hData3D->GetYaxis()->SetRange(k, k);
        TH2D* hDataTofSlice = dynamic_cast<TH2D*>(hData3D->Project3D("zxe"));
        hDataTofSlice->Write("hDataTofSlice");

        // hMuon3D->GetYaxis()->SetRange(k, k);
        // TH2D* hMuonTofSlice = dynamic_cast<TH2D*>(hMuon3D->Project3D("zxe"));
        // hMuonTofSlice->Write("hMuonTofSlice");

        // clone the projection and correct the spectrum using the coefficients
        TH1D* hDataTofCorrected = correctSpectra(hDataTofSlice, hDataCoefficients);
        hDataTofCorrected->Write("hDataTofCorrected");

        // TH1D* hMuonTofCorrected = correctSpectra(hMuonTofSlice, hMuonCoefficients);
        // hMuonTofCorrected->Write("hMuonTofCorrected");

        TH1D* hAxE = dynamic_cast<TH1D*>(corrections->Get(Form("%s/pTbin_%d/phibin_%d/efficiency", neutronClasses[i].Data(), j, k)));
        if (!hAxE) {
          std::cerr << "ERROR: Cannot find AxE histogram for neutron class " << neutronClasses[i].Data() << ", pT bin " << j << ", deltaPhi bin " << k << std::endl;
          continue;
        }
        hAxE->Write("hAxE");

        TH1D* hDataTofFinal = applyAxE(hDataTofCorrected, hAxE);
        hDataTofFinal->Write("correctedSpectrum");

        // TH1D* hMuonTofFinal = applyAxE(hMuonTofCorrected, hAxE);
        // hMuonTofFinal->Write("correctedMuonSpectrum");
      }
    }
  }
  scaledTof->Close();
}