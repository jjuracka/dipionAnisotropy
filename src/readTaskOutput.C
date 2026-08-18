// ROOT macro to read trees from file using RDataFrame

#include "commonFunctions.h"

// PID selection
bool tracksPassPID(const float leadingTpcPi, const float subleadingTpcPi) {
  const double radius = leadingTpcPi * leadingTpcPi + subleadingTpcPi * subleadingTpcPi;
  return radius < nSigmaCut * nSigmaCut;
}

// selections at system level
bool systemPassesCuts(const TLorentzVector& system, double minPt, double maxPt) {
  if (system.M() < minMass || system.M() > maxMass) return false;
  if (system.Pt() < minPt || system.Pt() >= maxPt) return false;
  if (std::abs(system.Rapidity()) > maxY) return false;
  return true;
}

// process reconstructed-level information (MC/data)
void processRecoTree(TString path, TString tree, TString outFileName, bool isMC = false) {
  // find tree path without having to look inside
  const TString treePath = findTreePath(path, tree);
  if (treePath.IsNull()) return;
  std::cout << "Found tree: " << treePath << std::endl;

  // RDataFrame set-up
  ROOT::RDataFrame df(treePath, path);
  ROOT::RDF::Experimental::AddProgressBar(df);

  // start applying filters and define variables
  auto cut = df.Filter(tracksPassPID, {"fLeadingTrackPiPID", "fSubleadingTrackPiPID"})
               .Define("systemLV", getRho, {"fLeadingTrackPt", "fSubleadingTrackPt", "fLeadingTrackEta", "fSubleadingTrackEta", "fLeadingTrackPhi", "fSubleadingTrackPhi"})
               .Define("mass", [](const TLorentzVector& lv) { return lv.M(); }, {"systemLV"})
               .Define("pT", [](const TLorentzVector& lv) { return lv.Pt(); }, {"systemLV"})
               .Define("rapidity", [](const TLorentzVector& lv) { return lv.Rapidity(); }, {"systemLV"})
               .Define("deltaPhi", getDeltaPhi, {"fLeadingTrackSign", "fSubleadingTrackSign", "fLeadingTrackPt", "fSubleadingTrackPt", "fLeadingTrackEta", "fSubleadingTrackEta", "fLeadingTrackPhi", "fSubleadingTrackPhi"});

  const auto neutronFilters = getNeutronFilters(isMC);
  const auto chargeFilters = getChargeFilters(isMC);

  // vectors for histo initialisation
  std::vector<ROOT::RDF::RResultPtr<TH1D>> vecM, vecPt, vecY, vecDeltaPhi;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> vecTotalFT0AmplitudeA, vecTotalFT0AmplitudeC, vecTotalFV0AmplitudeA, vecTotalFDDAmplitudeA, vecTotalFDDAmplitudeC, vecTimeFT0A, vecTimeFT0C, vecTimeFV0A, vecTimeFDDA, vecTimeFDDC;
  std::vector<ROOT::RDF::RResultPtr<TH2D>> vecDeltaPhiVsM, vecPtVsM, vecDeltaPhiVsPt, vecDeltaPhiVsPtWide;
  std::vector<std::tuple<std::string, std::string, int>> hNames;

  for (const auto& neutronFilter : neutronFilters) {
    for (const auto& chargeFilter : chargeFilters) {
      for (int i = 0; i < nBinsPt; ++i) {
        // select pT range
        const double minPt = pTbinEdges[i];
        const double maxPt = pTbinEdges[i + 1];
        auto fullCutWide = cut.Filter(neutronFilter.expression)
                              .Filter(chargeFilter.expression);
        auto fullCut = fullCutWide.Filter([minPt, maxPt](const TLorentzVector& lv) { return systemPassesCuts(lv, minPt, maxPt); }, {"systemLV"});

        // initialise histograms for storage in prepared vectors
        auto hM = fullCut.Histo1D<double>({"hM", Form("hM for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{m} (GeV/#it{c}^{2});counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), nBinsMass, minMass, maxMass}, "mass");
        vecM.push_back(hM);

        auto hPt = fullCut.Histo1D<double>({"hPt", Form("hPt for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{p}_{T} (GeV/#it{c});counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 100, 0.0, 0.1}, "pT");
        vecPt.push_back(hPt);

        auto hY = fullCut.Histo1D<double>({"hY", Form("hY for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{y};counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), nBinsY, -maxY, maxY}, "rapidity");
        vecY.push_back(hY);

        auto hDeltaPhi = fullCut.Histo1D<double>({"hDeltaPhi", Form("hDeltaPhi for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); %s (rad);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i, deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi()}, "deltaPhi");
        vecDeltaPhi.push_back(hDeltaPhi);

        auto hTotalFT0AmplitudeA = fullCut.Histo1D<float>({"hTotalFT0AmplitudeA", Form("hTotalFT0AmplitudeA for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); Total FT0 Amplitude A;counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 300, -100.0, 200.0}, "fTotalFT0AmplitudeA");
        vecTotalFT0AmplitudeA.push_back(hTotalFT0AmplitudeA);

        auto hTotalFT0AmplitudeC = fullCut.Histo1D<float>({"hTotalFT0AmplitudeC", Form("hTotalFT0AmplitudeC for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); Total FT0 Amplitude C;counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 300, -100.0, 200.0}, "fTotalFT0AmplitudeC");
        vecTotalFT0AmplitudeC.push_back(hTotalFT0AmplitudeC);

        auto hTotalFV0AmplitudeA = fullCut.Histo1D<float>({"hTotalFV0AmplitudeA", Form("hTotalFV0AmplitudeA for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); Total FV0 Amplitude A;counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 300, -100.0, 200.0}, "fTotalFV0AmplitudeA");
        vecTotalFV0AmplitudeA.push_back(hTotalFV0AmplitudeA);

        auto hTotalFDDAmplitudeA = fullCut.Histo1D<float>({"hTotalFDDAmplitudeA", Form("hTotalFDDAmplitudeA for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); Total FDD Amplitude A;counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 300, -100.0, 200.0}, "fTotalFDDAmplitudeA");
        vecTotalFDDAmplitudeA.push_back(hTotalFDDAmplitudeA);

        auto hTotalFDDAmplitudeC = fullCut.Histo1D<float>({"hTotalFDDAmplitudeC", Form("hTotalFDDAmplitudeC for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); Total FDD Amplitude C;counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 300, -100.0, 200.0}, "fTotalFDDAmplitudeC");
        vecTotalFDDAmplitudeC.push_back(hTotalFDDAmplitudeC);

        auto hTimeFT0A = fullCut.Histo1D<float>({"hTimeFT0A", Form("hTimeFT0A for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); FT0 A time (ns);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 400, -20.0, 20.0}, "fTimeFT0A");
        vecTimeFT0A.push_back(hTimeFT0A);

        auto hTimeFT0C = fullCut.Histo1D<float>({"hTimeFT0C", Form("hTimeFT0C for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); FT0 C time (ns);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 400, -20.0, 20.0}, "fTimeFT0C");
        vecTimeFT0C.push_back(hTimeFT0C);

        auto hTimeFV0A = fullCut.Histo1D<float>({"hTimeFV0A", Form("hTimeFV0A for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); FV0 A time (ns);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 400, -20.0, 20.0}, "fTimeFV0A");
        vecTimeFV0A.push_back(hTimeFV0A);

        auto hTimeFDDA = fullCut.Histo1D<float>({"hTimeFDDA", Form("hTimeFDDA for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); FDD A time (ns);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 400, -20.0, 20.0}, "fTimeFDDA");
        vecTimeFDDA.push_back(hTimeFDDA);

        auto hTimeFDDC = fullCut.Histo1D<float>({"hTimeFDDC", Form("hTimeFDDC for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); FDD C time (ns);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 400, -20.0, 20.0}, "fTimeFDDC");
        vecTimeFDDC.push_back(hTimeFDDC);

        auto hDeltaPhiVsM = fullCut.Histo2D<double, double>({"hDeltaPhiVsM", Form("hDeltaPhiVsM for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{m} (GeV/#it{c}^{2}); %s (rad);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i, deltaPhiLabel.c_str()), nBinsMass, minMass, maxMass, nBinsDeltaPhi, -TMath::Pi(), TMath::Pi()}, "mass", "deltaPhi");
        vecDeltaPhiVsM.push_back(hDeltaPhiVsM);

        auto hPtVsM = fullCut.Histo2D<double, double>({"hPtVsM", Form("hPtVsM for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{m} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c});counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), nBinsMass, minMass, maxMass, 100, 0.0, 0.1}, "mass", "pT");
        vecPtVsM.push_back(hPtVsM);

        auto hDeltaPhiVsPt = fullCut.Histo2D<double, double>({"hDeltaPhiVsPt", Form("hDeltaPhiVsPt for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{p}_{T} (GeV/#it{c}); %s (rad);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i, deltaPhiLabel.c_str()), 100, 0.0, 0.1, nBinsDeltaPhi, -TMath::Pi(), TMath::Pi()}, "pT", "deltaPhi");
        vecDeltaPhiVsPt.push_back(hDeltaPhiVsPt);

        // for reweighting
        auto hDeltaPhiVsPtWide = fullCutWide.Histo2D<double, double>({"hDeltaPhiVsPtWide", Form("hDeltaPhiVsPtWide for %s %s-signed events for #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{p}_{T} (GeV/#it{c}); %s (rad);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i, deltaPhiLabel.c_str()), 1000, 0.0, 1.0, 120, -TMath::Pi(), TMath::Pi()}, "pT", "deltaPhi");
        vecDeltaPhiVsPtWide.push_back(hDeltaPhiVsPtWide);

        // takes care of naming convention + directory creation
        hNames.emplace_back(neutronFilter.name, chargeFilter.name, i);
      }
    }
  }

  // make output file
  TFile* outFile = TFile::Open(Form("%s.root", outFileName.Data()), "RECREATE");
  if (!outFile) {
    std::cerr << "ERROR: Cannot create output file " << outFileName << std::endl;
    return;
  }

  // save histos in appropriate directory
  for (size_t i = 0; i < vecM.size(); ++i) {
    const auto& name = hNames[i];
    const TString neutronDirName = std::get<0>(name).c_str();
    const TString chargeDirName = std::get<1>(name).c_str();
    const int pTbin = std::get<2>(name);

    TDirectory* targetDir = createDirectoryHierarchy(outFile, {neutronDirName.Data(), chargeDirName.Data(), Form("pTbin_%d", pTbin)});
    targetDir->cd();

    vecM[i]->Write("hM");
    vecPt[i]->Write("hPt");
    vecY[i]->Write("hY");
    vecDeltaPhi[i]->Write("hDeltaPhi");
    vecDeltaPhiVsM[i]->Write("hDeltaPhiVsM");
    vecPtVsM[i]->Write("hPtVsM");
    vecDeltaPhiVsPt[i]->Write("hDeltaPhiVsPt");
    vecDeltaPhiVsPtWide[i]->Write("hDeltaPhiVsPtWide");
    vecTotalFT0AmplitudeA[i]->Write("hTotalFT0AmplitudeA");
    vecTotalFT0AmplitudeC[i]->Write("hTotalFT0AmplitudeC");
    vecTotalFV0AmplitudeA[i]->Write("hTotalFV0AmplitudeA");
    vecTotalFDDAmplitudeA[i]->Write("hTotalFDDAmplitudeA");
    vecTotalFDDAmplitudeC[i]->Write("hTotalFDDAmplitudeC");
    vecTimeFT0A[i]->Write("hTimeFT0A");
    vecTimeFT0C[i]->Write("hTimeFT0C");
    vecTimeFV0A[i]->Write("hTimeFV0A");
    vecTimeFDDA[i]->Write("hTimeFDDA");
    vecTimeFDDC[i]->Write("hTimeFDDC");
  }

  outFile->Close();
}

// process generated-level information (MC only)
void processGenTree(TString path, TString tree, TString outFileName) {
  const TString treePath = findTreePath(path, tree);
  if (treePath.IsNull()) return;
  std::cout << "Found tree: " << treePath << std::endl;

  ROOT::RDataFrame df(treePath, path);
  ROOT::RDF::Experimental::AddProgressBar(df);

  auto cut = df.Define("systemLV", getRho, {"fLeadingTrackPt", "fSubleadingTrackPt", "fLeadingTrackEta", "fSubleadingTrackEta", "fLeadingTrackPhi", "fSubleadingTrackPhi"})
               .Define("mass", [](const TLorentzVector& lv) { return lv.M(); }, {"systemLV"})
               .Define("pT", [](const TLorentzVector& lv) { return lv.Pt(); }, {"systemLV"})
               .Define("rapidity", [](const TLorentzVector& lv) { return lv.Rapidity(); }, {"systemLV"})
               .Define("deltaPhi", getDeltaPhi, {"fLeadingTrackSign", "fSubleadingTrackSign", "fLeadingTrackPt", "fSubleadingTrackPt", "fLeadingTrackEta", "fSubleadingTrackEta", "fLeadingTrackPhi", "fSubleadingTrackPhi"});

  const auto neutronFilters = getNeutronFilters(true);
  const auto chargeFilters = getChargeFilters(true);

  std::vector<ROOT::RDF::RResultPtr<TH1D>> vecM, vecPt, vecY, vecDeltaPhi;
  std::vector<ROOT::RDF::RResultPtr<TH2D>> vecDeltaPhiVsM, vecPtVsM, vecDeltaPhiVsPt, vecDeltaPhiVsPtWide;
  std::vector<std::tuple<std::string, std::string, int>> hNames;

  for (const auto& neutronFilter : neutronFilters) {
    for (const auto& chargeFilter : chargeFilters) {
      for (int i = 0; i < nBinsPt; ++i) {
        const double minPt = pTbinEdges[i];
        const double maxPt = pTbinEdges[i + 1];
        auto fullCut = cut.Filter([minPt, maxPt](const TLorentzVector& lv) { return systemPassesCuts(lv, minPt, maxPt); }, {"systemLV"});

        auto hM = fullCut.Histo1D<double>({"hM", Form("hM for %s %s-signed events in #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{m} (GeV/#it{c}^{2});counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), nBinsMass, minMass, maxMass}, "mass");
        vecM.push_back(hM);

        auto hPt = fullCut.Histo1D<double>({"hPt", Form("hPt for %s %s-signed events in #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{p}_{T} (GeV/#it{c});counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), 100, 0.0, 0.1}, "pT");
        vecPt.push_back(hPt);

        auto hY = fullCut.Histo1D<double>({"hY", Form("hY for %s %s-signed events in #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{y};counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), nBinsY, -maxY, maxY}, "rapidity");
        vecY.push_back(hY);

        auto hDeltaPhi = fullCut.Histo1D<double>({"hDeltaPhi", Form("hDeltaPhi for %s %s-signed events in #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); %s (rad);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i, deltaPhiLabel.c_str()), nBinsDeltaPhi, -TMath::Pi(), TMath::Pi()}, "deltaPhi");
        vecDeltaPhi.push_back(hDeltaPhi);

        auto hDeltaPhiVsM = fullCut.Histo2D<double, double>({"hDeltaPhiVsM", Form("hDeltaPhiVsM for %s %s-signed events in #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{m} (GeV/#it{c}^{2}); %s (rad);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i, deltaPhiLabel.c_str()), nBinsMass, minMass, maxMass, nBinsDeltaPhi, -TMath::Pi(), TMath::Pi()}, "mass", "deltaPhi");
        vecDeltaPhiVsM.push_back(hDeltaPhiVsM);

        auto hPtVsM = fullCut.Histo2D<double, double>({"hPtVsM", Form("hPtVsM for %s %s-signed events in #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{m} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c});counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i), nBinsMass, minMass, maxMass, 100, 0.0, 0.1}, "mass", "pT");
        vecPtVsM.push_back(hPtVsM);

        auto hDeltaPhiVsPt = fullCut.Histo2D<double, double>({"hDeltaPhiVsPt", Form("hDeltaPhiVsPt for %s %s-signed events in #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{p}_{T} (GeV/#it{c}); %s (rad);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i, deltaPhiLabel.c_str()), 100, 0.0, 0.1, nBinsDeltaPhi, -TMath::Pi(), TMath::Pi()}, "pT", "deltaPhi");
        vecDeltaPhiVsPt.push_back(hDeltaPhiVsPt);

        auto hDeltaPhiVsPtWide = cut.Histo2D<double, double>({"hDeltaPhiVsPtWide", Form("hDeltaPhiVsPtWide for %s %s-signed events in #it{p}_{T} #in (%.2f, %.2f) GeV/#it{c} (bin %d); #it{p}_{T} (GeV/#it{c}); %s (rad);counts", neutronFilter.name.c_str(), chargeFilter.name.c_str(), minPt, maxPt, i, deltaPhiLabel.c_str()), 1000, 0.0, 1.0, 120, -TMath::Pi(), TMath::Pi()}, "pT", "deltaPhi");
        vecDeltaPhiVsPtWide.push_back(hDeltaPhiVsPtWide);

        hNames.emplace_back(neutronFilter.name, chargeFilter.name, i);
      }
    }
  }

  TFile* outFile = TFile::Open(Form("%s.root", outFileName.Data()), "RECREATE");
  if (!outFile) {
    std::cerr << "ERROR: Cannot create output file " << outFileName << std::endl;
    return;
  }

  for (size_t i = 0; i < vecM.size(); ++i) {
    const auto& name = hNames[i];
    const TString neutronDirName = std::get<0>(name).c_str();
    const TString chargeDirName = std::get<1>(name).c_str();
    const int pTbin = std::get<2>(name);

    TDirectory* targetDir = createDirectoryHierarchy(outFile, {neutronDirName.Data(), chargeDirName.Data(), Form("pTbin_%d", pTbin)});
    targetDir->cd();

    vecM[i]->Write();
    vecPt[i]->Write();
    vecY[i]->Write();
    vecDeltaPhi[i]->Write();
    vecDeltaPhiVsM[i]->Write();
    vecPtVsM[i]->Write();
    vecDeltaPhiVsPt[i]->Write();
    vecDeltaPhiVsPtWide[i]->Write();
  }

  outFile->Close();
}

void readTaskOutput() {
  // global set-up
  globalSetup();

  if (pTbinEdges.size() < 2) {
    std::cerr << "ERROR: pTbinEdges must contain at least two elements to define a bin." << std::endl;
    return;
  }

  ROOT::EnableImplicitMT();
  gSystem->Exec("mkdir -p output");
  gSystem->ChangeDirectory("input");

  // data
  processRecoTree("MergedTree.root", "O2recotree", "data");

  // pion MC
  processRecoTree("MergedTreeMc.root", "O2recotree", "MCreco", true);
  processGenTree("MergedTreeMc.root", "O2mctree", "MCgen");

  // muon MC
  processRecoTree("MergedTreeMuon.root", "O2recotree", "muonReco", true);
  // no point in running processGenTree for muons -- since there are no pions in the generated MC, it will be empty
  
  // STARlight for reweighting
  for (const auto& r : rValues) {
    processGenTree("starlight.root", Form("tree_R%.2f", r), Form("starlight_R%.2f", r));
  }
}