#include "commonFunctions.h"

using WeightHistogramPtr = std::shared_ptr<TH2>;

// helper function to store copies of the weight histograms for inexpensive access
std::vector<WeightHistogramPtr> loadWeightHistograms() {
  std::vector<WeightHistogramPtr> weightHistograms;
  weightHistograms.reserve(rValues.size());

  TFile *weightsFile = TFile::Open("output/weights.root", "READ");
  if (!weightsFile || weightsFile->IsZombie()) {
    std::cerr << "ERROR: Cannot open weights file!" << std::endl;
    return weightHistograms;
  }

  for (const auto& r : rValues) {
    const std::string histoName = Form("hWeights_R%.2f", r);
    TH2* hWeights = dynamic_cast<TH2*>(weightsFile->Get(histoName.c_str()));
    if (!hWeights) {
      std::cerr << "Histogram " << histoName << " not found in weights file!" << std::endl;
      weightHistograms.emplace_back();
      continue;
    }

    TH2* cloned = dynamic_cast<TH2*>(hWeights->Clone());
    cloned->SetDirectory(nullptr);
    weightHistograms.emplace_back(WeightHistogramPtr(cloned));
  }

  return weightHistograms;
}

void reweightReco() {
  const TString treePath = findTreePath("input/resolutionTree.root", "O2resolutiontree");
  if (treePath.IsNull()) {
    return;
  }
  std::cout << "Found tree: " << treePath << std::endl;

  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df(treePath, "input/resolutionTree.root");
  ROOT::RDF::Experimental::AddProgressBar(df);

  const auto weightHistograms = loadWeightHistograms();

  // base DataFrame with calcuated gen and reco pT
  auto dfBase = df.Define("genSystemLV", getRho, {"fLeadingGenPt", "fSubleadingGenPt", "fLeadingGenEta", "fSubleadingGenEta", "fLeadingGenPhi", "fSubleadingGenPhi"})
                  .Define("genPt", [](const TLorentzVector& lv) { return lv.Pt(); }, {"genSystemLV"})
                  .Define("genDeltaPhi", getDeltaPhi, {"fLeadingSign", "fSubleadingSign", "fLeadingGenPt", "fSubleadingGenPt", "fLeadingGenEta", "fSubleadingGenEta", "fLeadingGenPhi", "fSubleadingGenPhi"})
                  .Define("recoSystemLV", getRho, {"fLeadingRecoPt", "fSubleadingRecoPt", "fLeadingRecoEta", "fSubleadingRecoEta", "fLeadingRecoPhi", "fSubleadingRecoPhi"})
                  .Define("recoPt", [](const TLorentzVector& lv) { return lv.Pt(); }, {"recoSystemLV"})
                  .Filter("recoPt <= 0.1")
                  .Define("recoY", [](const TLorentzVector& lv) { return lv.Rapidity(); }, {"recoSystemLV"})
                  .Filter("std::abs(recoY) <= 0.9")
                  .Define("recoDeltaPhi", getDeltaPhi, {"fLeadingSign", "fSubleadingSign", "fLeadingRecoPt", "fSubleadingRecoPt", "fLeadingRecoEta", "fSubleadingRecoEta", "fLeadingRecoPhi", "fSubleadingRecoPhi"});

  std::vector<ROOT::RDF::RResultPtr<TH2D>> hRecoPtWeighted; // storage for the weighted histograms
  for (int i = 0; i < rValues.size(); ++i) {
    // get the corresponding weight histogram
    const auto& weightHist = weightHistograms[i];
    if (!weightHist) continue;
    const TString weightColumnName = Form("weight_R%d", i);

    // define a new column in the DataFrame that contains the weights based on genPt and genDeltaPhi
    auto dfWithWeight = dfBase.Define(weightColumnName, [weightHist](double genPt, double genDeltaPhi) { return weightHist->GetBinContent(weightHist->FindBin(genPt, genDeltaPhi)); }, {"genPt", "genDeltaPhi"});
    // create a weighted histogram of recoPt using the weights
    auto hRecoPt = dfWithWeight.Histo2D<double, double>({Form("hRecoPtWeighted_R%.2f", rValues[i]), Form("hRecoPtWeighted for #it{R} = %.2f; #it{p}_{T} (GeV/#it{c}); #Delta#phi (rad); weighted counts", rValues[i]), 100, 0, 0.1, 12, -TMath::Pi(), TMath::Pi()}, "recoPt" , "recoDeltaPhi", weightColumnName);
    hRecoPtWeighted.push_back(hRecoPt);
  }

  // save the histograms to a new ROOT file
  TFile* outFile = TFile::Open("output/reweightedReco.root", "RECREATE");
  for (int i = 0; i < hRecoPtWeighted.size(); ++i) hRecoPtWeighted[i]->Write();
  outFile->Close();
}