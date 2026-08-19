#include "commonFunctions.h"

void calculateWeights() {
  globalSetup();

  TFile *fUpc = TFile::Open("input/MCgen.root"); // central UPC production
  TH2* hUpc = dynamic_cast<TH2*>(fUpc->Get("AnAn/unlike/pTbin_0/hDeltaPhiVsPtWide")); // the "wide" histogram has finer binning + no cut in pT
  if (!hUpc) {
    std::cerr << "Histogram hDeltaPhiVsPtWide not found in MCgen.root!" << std::endl;
    return;
  }
  
  TFile *fStarLight = TFile::Open("input/starlight.root"); // STARlight productions with different settings

  TFile *fWeights = TFile::Open("output/weights.root", "RECREATE"); // output file to store the reweighting factors
  
  // loop through all histograms in the StarLight file and calculate the reweighting factors
  for (const auto& r : rValues) {
    std::string histoName = Form("hDeltaPhiVsPt_R%.2f", r);
    TH2* hStarLight = dynamic_cast<TH2*>(fStarLight->Get(histoName.c_str()));

    TH2* hWeights = dynamic_cast<TH2*>(hUpc->Clone(Form("hWeights_R%.2f", r)));
    
    // check for bin consistency between the two histograms
    if (!hStarLight) {
      std::cerr << "Histogram " << histoName << " not found in one of the files!" << std::endl;
      continue;
    }
    if (hUpc->GetNbinsX() != hStarLight->GetNbinsX() || hUpc->GetNbinsY() != hStarLight->GetNbinsY()) {
      std::cerr << "Bin mismatch between UPC and STARlight histograms for R = " << r << std::endl;
      continue;
    }
    
    for (int yBin = 1; yBin <= hWeights->GetNbinsY(); yBin++) {
      // take a projection of the 2D histos to properly calculate the weights per DeltaPhi bin
      TH1D *hUpcProj = hUpc->ProjectionX(Form("hUpcProj_R%.2f", r), yBin, yBin);
      TH1D *hStarLightProj = hStarLight->ProjectionX(Form("hStarLightProj_R%.2f", r), yBin, yBin);
      // scaling by the integral of the primary diffractive peak (pT < 0.1 GeV/c)
      hUpcProj->Scale(1.0 / hUpcProj->Integral(1, hUpcProj->FindBin(0.1)));
      hStarLightProj->Scale(1.0 / hStarLightProj->Integral(1, hStarLightProj->FindBin(0.1)));
      
      // calculate the weights for this DeltaPhi bin
      TH1* hWeightsRow = dynamic_cast<TH1*>(hStarLightProj->Clone(Form("hWeightsRow_R%.2f", r)));
      hWeightsRow->Divide(hUpcProj);
      
      // fill the weights into the 2D histogram
      for (int xBin = 1; xBin <= hWeights->GetNbinsX(); ++xBin) {
        hWeights->SetBinContent(xBin, yBin, hWeightsRow->GetBinContent(xBin));
        hWeights->SetBinError(xBin, yBin, hWeightsRow->GetBinError(xBin));
      }      
    }
    hWeights->Write();
  }
  fWeights->Close();
}