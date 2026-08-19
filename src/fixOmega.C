#include "commonFunctions.h"

const double localFitMin = 0.73;
const double localFitMax = 0.83;

void fixOmega() {
  globalSetup();
  ROOT::EnableImplicitMT();

  TFile *fIn = new TFile("output/AxEcorrected.root", "READ");

  // fit the mass spectrum of the whole sample to fix the omega phase
  // first add up the muon and pion mass spectra for all pT and phi bins
  TH1D *hMuonSum = new TH1D("hMuonSum", ";#it{m}_{#mu#mu} (GeV/#it{c}^{2});corrected counts", nBinsMass, minMass, maxMass);
  TH1D *hPionSum = new TH1D("hPionSum", ";#it{m}_{#pi#pi} (GeV/#it{c}^{2});corrected counts", nBinsMass, minMass, maxMass);
  // no need to separate into neutron classes, just take all the AnAn events
  for (int j = 0; j < nBinsPt; j++) {
    for (int k = 1; k <= nBinsDeltaPhi; ++k) {
      TH1D *hMuon = (TH1D*) fIn->Get(Form("AnAn/pTbin_%d/phibin_%d/correctedMuonSpectrum", j, k));
      TH1D *hPion = (TH1D*) fIn->Get(Form("AnAn/pTbin_%d/phibin_%d/correctedSpectrum", j, k));
      hMuonSum->Add(hMuon);
      hPionSum->Add(hPion);
    }
  }

  // fit function for muons
  TF1 *fitFuncMuons = new TF1("fitFuncMuons", "[0]*x^(-[1])", localFitMin, localFitMax);
  fitFuncMuons->SetParNames("#it{a}_{#mu}", "#it{b}_{#mu}");
  fitFuncMuons->SetParameters(1000.0, -10.0);

  // fit the muons to get the slope
  TFitResultPtr muonFit = hMuonSum->Fit(fitFuncMuons, "EMR0S", "", localFitMin, localFitMax);
  if (!muonFit) std::cerr << "Error: Muon fit failed for the whole sample" << std::endl;
  double muonSlope = fitFuncMuons->GetParameter(1);
  // plot the muon fit result
  TCanvas *cMuon = new TCanvas("cMuon", "cMuon", 1000, 600);
  cMuon->SetLogy();
  hMuonSum->SetLineColorAlpha(tabBlue, 1.0);
  hMuonSum->SetMarkerColor(hMuonSum->GetLineColor());
  hMuonSum->SetLineWidth(2);
  hMuonSum->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
  hMuonSum->GetYaxis()->SetTitle("corrected counts");
  hMuonSum->Draw("e");
  fitFuncMuons->SetLineColorAlpha(tabOrange, 1.0);
  fitFuncMuons->SetNpx(1000);
  fitFuncMuons->SetLineWidth(2);
  fitFuncMuons->Draw("same");
  // cMuon->SaveAs("output/fitResults/muonFitGlobal.pdf");

  // fit function for soeding
  TF1 *fitFuncSoeding = new TF1("fitFuncSoeding", Soeding, localFitMin, localFitMax, 10);
  fitFuncSoeding->SetParNames("#it{A}_{#rho}", "#it{m}_{#rho}", "#it{#Gamma}_{#rho}", "#it{C}_{#omega}", "#it{m}_{#omega}", "#it{#Gamma}_{#omega}", "#it{#phi}_{#omega}", "#it{B}_{#pi#pi}", "#it{a}_{#mu}", "#it{b}_{#mu}");
  fitFuncSoeding->SetParameters(500.0, kMrho, kWrho, -10.0, kMomega, kWomega, -1.2, -315.0, muonSlope);
  fitFuncSoeding->FixParameter(1, kMrho);
  fitFuncSoeding->FixParameter(2, kWrho);
  fitFuncSoeding->FixParameter(4, kMomega);
  fitFuncSoeding->FixParameter(5, kWomega);
  fitFuncSoeding->SetParLimits(6, -TMath::Pi(), TMath::Pi());
  fitFuncSoeding->FixParameter(9, muonSlope);
  TFitResultPtr soedingFit = hPionSum->Fit(fitFuncSoeding, "EMR0S", "", localFitMin, localFitMax);
  if (!soedingFit) std::cerr << "Error: Soeding fit failed for the whole sample" << std::endl;
  double omegaPhase = fitFuncSoeding->GetParameter(6);
  double omegaPhaseError = fitFuncSoeding->GetParError(6);
  std::cout << "Omega phase: " << omegaPhase << " +/- " << omegaPhaseError << " rad" << std::endl;
  // plot the soeding fit result
  TCanvas *cSoeding = new TCanvas("cSoeding", "cSoeding", 1000, 600);
  // cSoeding->SetLogy();
  hPionSum->SetLineColorAlpha(tabBlue, 1.0);
  hPionSum->SetMarkerColor(hPionSum->GetLineColor());
  hPionSum->SetLineWidth(2);
  hPionSum->GetXaxis()->SetTitle("#it{m}_{#pi#pi} (GeV/#it{c}^{2})");
  hPionSum->GetYaxis()->SetTitle("corrected counts");
  hPionSum->Draw("e");
  fitFuncSoeding->SetLineColorAlpha(tabOrange, 1.0);
  fitFuncSoeding->SetNpx(1000);
  fitFuncSoeding->SetLineWidth(2);
  fitFuncSoeding->Draw("same");
  // cSoeding->SaveAs("output/fitResults/soedingFitGlobal.pdf");

  // save result into a root file
  TFile *fOut = new TFile("output/omegaPhase.root", "RECREATE");
  fOut->cd();
  TH1D *hOmegaPhase = new TH1D("hOmegaPhase", "hOmegaPhase;#Delta#it{#phi} (rad);#it{#phi}_{#omega} (rad)", 1, -TMath::Pi(), TMath::Pi());
  hOmegaPhase->SetBinContent(1, omegaPhase);
  hOmegaPhase->SetBinError(1, fitFuncSoeding->GetParError(6));
  hOmegaPhase->Write();
  cSoeding->Write();
  cMuon->Write();
  fOut->Close();
}