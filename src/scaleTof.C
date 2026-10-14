#include "commonFunctions.h"

TH1D* getCoefficients(TH1D* hNoTof, TH1D* hTof) {
  if (!hNoTof || !hTof) {
    std::cerr << "ERROR: One or both histograms are null!" << std::endl;
    return nullptr;
  }

  if (hNoTof->GetNbinsX() != hTof->GetNbinsX()) {
    std::cerr << "ERROR: Histograms have different number of bins!" << std::endl;
    return nullptr;
  }

  TH1D* hCoefficients = dynamic_cast<TH1D*>(hNoTof->Clone("hCoefficients"));
  hCoefficients->SetTitle("TOF Scaling Coefficients;#it{m} (GeV/#it{c}^{2});scaling coefficient");
  hCoefficients->Divide(hTof);

  return hCoefficients;
}

TH1D* correctSpectra(TH1D* hTof, TH1D* hCoefficients) {
  if (!hTof || !hCoefficients) {
    std::cerr << "ERROR: One or both histograms are null!" << std::endl;
    return nullptr;
  }

  if (hTof->GetNbinsX() != hCoefficients->GetNbinsX()) {
    std::cerr << "ERROR: Histograms have different number of bins!" << std::endl;
    return nullptr;
  }

  TH1D* hCorrected = dynamic_cast<TH1D*>(hTof->Clone("hCorrected"));
  hCorrected->SetTitle("TOF-corrected spectrum;#it{m} (GeV/#it{c}^{2});counts");
  hCorrected->Multiply(hCoefficients);

  return hCorrected;
}

void scaleTof() {
  globalSetup();
  
  TFile* dataAxE = TFile::Open("output/AxEcorrected.root", "READ");
  TFile* dataBase = TFile::Open("input/data.root", "READ");
  TFile* muonBase = TFile::Open("input/muonReco.root", "READ");

  TFile* scaledTof = TFile::Open("output/scaledTof.root", "RECREATE");

  if (!dataAxE || !dataBase || !muonBase) {
    std::cerr << "ERROR: Cannot open one or more input files!" << std::endl;
    return;
  }

  for (int i = 0; i < neutronClasses.size(); i++) {
    TDirectory* neutronDir = scaledTof->mkdir(neutronClasses[i].Data());
    neutronDir->cd();
    for (int j = 0; j < nBinsPt; j++) {
      TDirectory* pTDir = neutronDir->mkdir(Form("pTbin_%d", j));
      pTDir->cd();

      TH2D* hMuonTof = dynamic_cast<TH2D*>(muonBase->Get(Form("AnAn/unlike/pTbin_%d/hDeltaPhiVsMwithTof", j)));
      TH2D* hDataTof = dynamic_cast<TH2D*>(dataBase->Get(Form("%s/unlike/pTbin_%d/hDeltaPhiVsMwithTof", neutronClasses[i].Data(), j)));
      TH2D* hDataNoTof = dynamic_cast<TH2D*>(dataBase->Get(Form("%s/unlike/pTbin_%d/hDeltaPhiVsM", neutronClasses[i].Data(), j)));

      for (int k = 1; k <= nBinsDeltaPhi; ++k) {
        TDirectory* phiDir = pTDir->mkdir(Form("phibin_%d", k));
        phiDir->cd();

        TH1D* hDataTofSlice = hDataTof->ProjectionX(Form("hDataTofSlice_pTbin_%d_deltaPhiBin_%d", j, k), k, k);
        TH1D* hMuonTofSlice = hMuonTof->ProjectionX(Form("hMuonTofSlice_pTbin_%d_deltaPhiBin_%d", j, k), k, k);
        TH1D* hDataNoTofSlice = dynamic_cast<TH1D*>(dataAxE->Get(Form("%s/pTbin_%d/phibin_%d/correctedSpectrum", neutronClasses[i].Data(), j, k)));
        TH1D* hMuonNoTofSlice = dynamic_cast<TH1D*>(dataAxE->Get(Form("AnAn/pTbin_%d/phibin_%d/correctedMuonSpectrum", j, k)));

        TH1D* hDataRawNoTofSlice = hDataNoTof->ProjectionX(Form("hDataRawNoTofSlice_pTbin_%d_deltaPhiBin_%d", j, k), k, k);
        TH1D* hTofFraction = getCoefficients(hDataTofSlice, hDataRawNoTofSlice);
        hTofFraction->Write("hTofFraction");

        TH1D* hDataCoefficients = getCoefficients(hDataNoTofSlice, hDataTofSlice);
        TH1D* hMuonCoefficients = getCoefficients(hMuonNoTofSlice, hMuonTofSlice);

        hDataCoefficients->Write("hDataCoefficients");
        hMuonCoefficients->Write("hMuonCoefficients");

        TH1D* hDataCorrected = correctSpectra(hDataTofSlice, hDataCoefficients);
        TH1D* hMuonCorrected = correctSpectra(hMuonTofSlice, hMuonCoefficients);

        hDataCorrected->Write("hCorrectedSpectrum");
        hMuonCorrected->Write("hCorrectedMuonSpectrum");
      }
    }
  }
  scaledTof->Close();
}