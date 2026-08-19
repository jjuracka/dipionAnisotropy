#include "commonFunctions.h"

using WeightHistogramPtr = std::shared_ptr<TH2>;

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

double getBestR(TFile* fOptimalR, const double genDeltaPhi) {
  TH1* hAveragedR = dynamic_cast<TH1*>(fOptimalR->Get("0n0n/hAveragedR"));
  if (!hAveragedR) {
    std::cerr << "Error: Could not find histogram hAveragedR in file " << fOptimalR->GetName() << std::endl;
    return 0.0;
  }
  // get the bin corresponding to the genDeltaPhi
  int bin = hAveragedR->GetXaxis()->FindBin(genDeltaPhi);
  double bestR = hAveragedR->GetBinContent(bin);
  // round it to nearest .1 (granularity of the available R values)
  bestR = std::round(bestR * 10.0) / 10.0;
  return bestR;
}

double getWeight(const std::vector<WeightHistogramPtr>& weightHistograms, const double genPt, const double genDeltaPhi, const double bestR) {
  // find the index of the bestR in rValues
  auto it = std::find(rValues.begin(), rValues.end(), bestR);
  if (it == rValues.end()) {
    std::cerr << "Error: bestR " << bestR << " not found in rValues!" << std::endl;
    return 1.0;
  }
  size_t index = std::distance(rValues.begin(), it);
  if (index >= weightHistograms.size()) {
    std::cerr << "Error: index " << index << " out of range for weightHistograms!" << std::endl;
    return 1.0;
  }
  const auto& hWeights = weightHistograms[index];
  if (!hWeights) {
    std::cerr << "Error: weight histogram for R=" << bestR << " is null!" << std::endl;
    return 1.0;
  }
  // get the bin corresponding to genPt and genDeltaPhi
  int binX = hWeights->GetXaxis()->FindBin(genPt);
  int binY = hWeights->GetYaxis()->FindBin(genDeltaPhi);
  double weight = hWeights->GetBinContent(binX, binY);
  return weight;
}

TH2D* getNumerator(int pTbin) {
  TFile* fResolution = new TFile("input/resolutionTree.root", "READ");
  TFile* fOptimalR = new TFile("output/optimalR.root", "READ");
  // the file includes a tree with generated and reconstructed pion pairs
  const TString treePath = findTreePath("input/resolutionTree.root", "O2resolutiontree");
  if (treePath.IsNull()) {
    return nullptr;
  }
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df(treePath, "input/resolutionTree.root");
  ROOT::RDF::Experimental::AddProgressBar(df);

  const auto weightHistograms = loadWeightHistograms();
  // build an in-memory lookup of averaged R per delta-phi bin (thread-safe)
  TH1* hAveragedR = dynamic_cast<TH1*>(fOptimalR->Get("0n0n/hAveragedR"));
  if (!hAveragedR) {
    std::cerr << "Error: hAveragedR not found in optimalR.root" << std::endl;
    return nullptr;
  }

  double minPt = pTbinEdges[pTbin];
  double maxPt = pTbinEdges[pTbin + 1];

  // base DataFrame with calculated gen and reco pT
  auto dfBase = df.Define("genSystemLV", getRho, {"fLeadingGenPt", "fSubleadingGenPt", "fLeadingGenEta", "fSubleadingGenEta", "fLeadingGenPhi", "fSubleadingGenPhi"})
                  .Define("genPt", [](const TLorentzVector& lv) { return lv.Pt(); }, {"genSystemLV"})
                  .Define("genDeltaPhi", getDeltaPhi, {"fLeadingSign", "fSubleadingSign", "fLeadingGenPt", "fSubleadingGenPt", "fLeadingGenEta", "fSubleadingGenEta", "fLeadingGenPhi", "fSubleadingGenPhi"})
                  .Define("recoSystemLV", getRho, {"fLeadingRecoPt", "fSubleadingRecoPt", "fLeadingRecoEta", "fSubleadingRecoEta", "fLeadingRecoPhi", "fSubleadingRecoPhi"})
                  .Define("recoPt", [](const TLorentzVector& lv) { return lv.Pt(); }, {"recoSystemLV"})
                  .Filter([minPt, maxPt](double recoPt) { return recoPt >= minPt && recoPt < maxPt; }, {"recoPt"})
                  .Define("recoDeltaPhi", getDeltaPhi, {"fLeadingSign", "fSubleadingSign", "fLeadingRecoPt", "fSubleadingRecoPt", "fLeadingRecoEta", "fSubleadingRecoEta", "fLeadingRecoPhi", "fSubleadingRecoPhi"})
                  .Define("recoY", [](const TLorentzVector& lv) { return lv.Rapidity(); }, {"recoSystemLV"})
                  .Filter("std::abs(recoY) <= 0.9")
                  .Define("recoM", [](const TLorentzVector& lv) { return lv.M(); }, {"recoSystemLV"})
                  // define bestR depending on the generated deltaPhi, using the averagedR histogram
                  .Define("bestR", [hAveragedR](const double genDeltaPhi) {
                    int bin = hAveragedR->GetXaxis()->FindBin(genDeltaPhi);
                    double bestR = hAveragedR->GetBinContent(bin);
                    bestR = std::round(bestR * 10.0) / 10.0; // round to nearest .1
                    return bestR;
                  }, {"genDeltaPhi"})
                  .Define("weight", [weightHistograms](const double genPt, const double genDeltaPhi, const double bestR) { return getWeight(weightHistograms, genPt, genDeltaPhi, bestR); }, {"genPt", "genDeltaPhi", "bestR"});

  // we need to apply a weight to the reco events based on the generated pT and deltaPhi, depending on the deltaPhi, a different weight histogram is used
  // create weighted numerator histogram via RDataFrame
  auto hNumRes = dfBase.Histo2D({"hNumerator", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi} (rad)", nBinsMass, minMass, maxMass, nBinsDeltaPhi, -TMath::Pi(), TMath::Pi()}, "recoM", "recoDeltaPhi", "weight");
  // clone the result to own it outside the RDF lifetime
  TH2D* hNumerator = dynamic_cast<TH2D*>(hNumRes->Clone("hNumerator_cloned"));
  hNumerator->SetDirectory(nullptr);
  return hNumerator;
}

TH2D* getDenominator(int pTbin) {
  TFile *fOptimalR = new TFile("output/optimalR.root", "READ");
  TH2D *hDenominator = new TH2D("hDenominator", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi} (rad)", nBinsMass, minMass, maxMass, nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());
  hDenominator->Sumw2();
  TH1* hAveragedR = dynamic_cast<TH1*>(fOptimalR->Get("0n0n/hAveragedR"));
  if (!hAveragedR) {
    std::cerr << "Error: hAveragedR not found in optimalR.root" << std::endl;
    return nullptr;
  }
  if (hDenominator->GetNbinsY() != hAveragedR->GetNbinsX()) {
    std::cerr << "Error: Number of deltaPhi bins in denominator histogram does not match number of bins in hAveragedR!" << std::endl;
    return nullptr;
  }
  // the different deltaPhi bins will basically store 1D histograms of the invariant mass according to the averaged R value for that deltaPhi bin
  for (int i = 1; i <= hDenominator->GetNbinsY(); ++i) {
    double bestR = hAveragedR->GetBinContent(i);
    // round to nearest .1
    bestR = std::round(bestR * 10.0) / 10.0;
    TFile *fStarlight = new TFile(Form("input/starlight_R%.2f.root", bestR), "READ");
    const std::string histoName = Form("AnAn/unlike/pTbin_%d/hDeltaPhiVsM", pTbin);
    TH2D* hMass = dynamic_cast<TH2D*>(fStarlight->Get(histoName.c_str()));
    if (!hMass) {
      std::cerr << "Error: Histogram " << histoName << " not found in starlight.root!" << std::endl;
      continue;
    }
    if (hMass->GetNbinsX() != hDenominator->GetNbinsX() || hMass->GetNbinsY() != hDenominator->GetNbinsY()) {
      std::cerr << "Error: Histogram " << histoName << " has different number of bins than denominator histogram!" << std::endl;
      continue;
    }
    // fill the corresponding deltaPhi bin in the 2D histogram with the mass histogram
    for (int j = 1; j <= hMass->GetNbinsX(); ++j) {
      hDenominator->SetBinContent(j, i, hMass->GetBinContent(j, i));
      hDenominator->SetBinError(j, i, hMass->GetBinError(j, i));
    }
  }
  return hDenominator;
}

void reweightAxE() {
  globalSetup();
  // output file
  TFile *fOut = new TFile("output/reweightedAxE.root", "RECREATE");
  for (int pTbin = 0; pTbin < nBinsPt; ++pTbin) {
    TH2D* hNumerator = getNumerator(pTbin);
    TH2D* hDenominator = getDenominator(pTbin);
    if (!hNumerator || !hDenominator) {
      std::cerr << "Error: Could not get numerator or denominator for pT bin " << pTbin << std::endl;
      continue;
    }
    fOut->cd();
    TDirectory* dir = fOut->mkdir(Form("AnAn/unlike/pTbin_%d", pTbin));
    dir->cd();
    hNumerator->Write("hNumerator");
    hDenominator->Write("hDenominator");
  }
  fOut->Close();
}