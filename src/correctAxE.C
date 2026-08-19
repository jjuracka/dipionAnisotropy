#include "commonFunctions.h"

// safeguard function
void reconcileHistos(TH1D *hGen, TH1D *hPassed) {
  // ensure the histograms have the same number of bins
  if (hGen->GetNbinsX() != hPassed->GetNbinsX()) {
    std::cerr << "Error: Histograms have different number of bins!" << std::endl;
    return;
  }

  // loop over all bins and make sure the passed histogram does not have more entries than the generated histogram
  for (int bin = 1; bin <= hGen->GetNbinsX(); ++bin) {
    double genContent = hGen->GetBinContent(bin);
    double passedContent = hPassed->GetBinContent(bin);
    if (passedContent > genContent) { 
      hPassed->SetBinContent(bin, genContent);
      std::cout << "Warning: Bin " << bin << " in passed histogram has more entries than generated histogram. Setting " << passedContent << " to " << genContent << "." << std::endl;
    }
  }
}

// TEfficiency cannot handle weighted histograms, so we need to build the efficiency histogram manually
TH1D* buildEfficiencyHistogram(TH1D *hPassed, TH1D *hTotal, const TString &name) {
  TH1D *hEff = (TH1D*)hPassed->Clone(name);
  hEff->Reset();

  for (int bin = 1; bin <= hEff->GetNbinsX(); ++bin) {
    const double total = hTotal->GetBinContent(bin);
    const double passed = hPassed->GetBinContent(bin);
    const double totalErr = hTotal->GetBinError(bin);
    const double passedErr = hPassed->GetBinError(bin);

    if (total > 0) {
      const double efficiency = passed / total;
      const double efficiencyErr = TMath::Sqrt(
        TMath::Power(passedErr / total, 2) +
        TMath::Power(passed * totalErr / (total * total), 2)
      );
      hEff->SetBinContent(bin, efficiency);
      hEff->SetBinError(bin, efficiencyErr);
    } else {
      hEff->SetBinContent(bin, 0);
      hEff->SetBinError(bin, 0);
    }
  }

  return hEff;
}

TH1D* correctSpectrum(TH1D *hist, TH1D *eff) {
  TH1D *hDiv = (TH1D*)hist->Clone("hDiv");
  hDiv->Reset();

  for (int bin = 1; bin <= hDiv->GetNbinsX(); ++bin) {
    // get the efficiency for the current bin
    double efficiency = eff->GetBinContent(bin);
    if (efficiency > 0) {
      // get the bin content and error
      double content = hist->GetBinContent(bin);
      double error = hist->GetBinError(bin);
      // calculate the corrected content and error
      double newContent = content / efficiency;
      double effError = eff->GetBinError(bin);
      double newError = TMath::Sqrt(
        TMath::Power(error / efficiency, 2) +
        TMath::Power(content * effError / (efficiency * efficiency), 2)
      );
      // set the corrected content and error
      hDiv->SetBinContent(bin, newContent);
      hDiv->SetBinError(bin, newError);
    } else { // if efficiency is zero, set content and error to zero
      hDiv->SetBinContent(bin, 0);
      hDiv->SetBinError(bin, 0);
    }
  }
  // return the divided histogram
  return hDiv;
}

void correctAxE() {
  globalSetup();
  // load necessary files
  // data files
  TFile *fData = new TFile("output/subtracted.root", "READ");
  TFile *fMuon = new TFile("input/muonReco.root", "READ");
  // for AxE calculation
  TFile *fReweighted = new TFile("output/reweightedAxE.root", "READ");
  if (!fData || !fMuon || !fReweighted) {
    std::cerr << "Error: One or more input files could not be opened!" << std::endl;
    return;
  }
  // output
  TFile *fOut = new TFile("output/AxEcorrected.root", "RECREATE");

  for (int i = 0; i < neutronClasses.size(); i++) {
    for (int j = 0; j < nBinsPt; j++) {
      TH2D *hRecMc = (TH2D*)fReweighted->Get(Form("AnAn/unlike/pTbin_%d/hNumerator", j));
      TH2D *hGenMc = (TH2D*)fReweighted->Get(Form("AnAn/unlike/pTbin_%d/hDenominator", j));
      TH2D *hData = (TH2D*)fData->Get(Form("%s/pTbin_%d/hDeltaPhiVsM", neutronClasses[i].Data(), j));
      TH2D *hMuon = (TH2D*)fMuon->Get(Form("AnAn/unlike/pTbin_%d/hDeltaPhiVsM", j));
      for (int k = 1; k <= hGenMc->GetNbinsY(); k++) {
        gSystem->Exec(Form("mkdir -p output/fitResults/%s/pTbin_%d/phiBin_%d", neutronClasses[i].Data(), j, k));
        TDirectory *dir = fOut->mkdir(Form("%s/pTbin_%d/phibin_%d", neutronClasses[i].Data(), j, k));
        dir->cd();
        // make efficiency
        TH1D *hGenMc_proj = hGenMc->ProjectionX("hGenMc_proj", k, k);
        TH1D *hRecMc_proj = hRecMc->ProjectionX("hRecMc_proj", k, k);
        reconcileHistos(hGenMc_proj, hRecMc_proj);
        TH1D *eff = buildEfficiencyHistogram(hRecMc_proj, hGenMc_proj, Form("hEfficiency_%s_pTbin_%d_phiBin_%d", neutronClasses[i].Data(), j, k));
        eff->Write("efficiency");
        // project data histograms and correct for efficiency
        TH1D *hData_proj = hData->ProjectionX("hData_proj", k, k);
        TH1D *hData_corr = correctSpectrum(hData_proj, eff);
        TH1D *hMuon_proj = hMuon->ProjectionX("hMuon_proj", k, k);
        TH1D *hMuon_corr = correctSpectrum(hMuon_proj, eff);
        // save histograms
        hData_corr->Write("correctedSpectrum");
        hMuon_corr->Write("correctedMuonSpectrum");
      }
    }
  }
}