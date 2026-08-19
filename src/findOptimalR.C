// root macro to copare MC with data to find optimal setting of the R parameter

#include "commonFunctions.h"

// range for the calculation of chi2
double rangeMin = 0.03;
double rangeMax = 0.10;
// TGaxis::SetMaxDigits(4);

double getDeltaPhiBinCenter(int bin) {
  double deltaPhiMin = -TMath::Pi();
  double deltaPhiMax = TMath::Pi();
  double binWidth = (deltaPhiMax - deltaPhiMin) / nBinsDeltaPhi;
  return deltaPhiMin + (bin + 0.5) * binWidth;
}

double deltaPhiBinWidth = (2.0 * TMath::Pi()) / nBinsDeltaPhi;

void findOptimalR() {
  globalSetup();
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1e6); // increase possible number of iterations for fitting
  
  TFile *fOptimalR = TFile::Open("output/optimalR.root", "RECREATE"); // output file to store the optimal R values
  // load data histogram
  TFile *fData = TFile::Open("input/data.root");

  for (const auto& neutronClass : neutronClasses) {
    if (neutronClass != "0n0n") continue; // only do 0n0n
    // create a directory in the optimalR.root file for each neutron class
    fOptimalR->mkdir(neutronClass);
    fOptimalR->cd(neutronClass);

    TH2* hData = dynamic_cast<TH2*>(fData->Get(Form("%s/unlike/pTbin_0/hDeltaPhiVsPt", neutronClass.Data())));
  
    // load weighted MC histograms for different R values
    TFile *fMC = TFile::Open("output/reweightedReco.root");
    std::vector<TH2*> mcHistos;
    for (const auto& r : rValues) {
      TH2* hMc = dynamic_cast<TH2*>(fMC->Get(Form("hRecoPtWeighted_R%.2f", r)));
      mcHistos.push_back(hMc);
    }
    // "output" histogram to store the optimal R values for each delta phi bin
    TH1D *hOptimalR = new TH1D("hOptimalR", ";#Delta#it{#phi} (rad);optimal #it{R} (fm)", nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());

    for (int i = 0; i < nBinsDeltaPhi; ++i) {
      TH1D* hDataProj = hData->ProjectionX(Form("hDataProj_bin%d", i + 1), i + 1, i + 1);
      hDataProj->Scale(1.0 / hDataProj->Integral(hDataProj->FindBin(rangeMin), hDataProj->FindBin(rangeMax)));

      std::vector<double> chi2Values;
      for (size_t j= 0; j < mcHistos.size(); ++j) {
        TH1D* hMcProj = mcHistos[j]->ProjectionX(Form("hMcProj_R%.2f_bin%d", rValues[j], i + 1), i + 1, i + 1);
        hMcProj->Scale(1.0 / hMcProj->Integral(hMcProj->FindBin(rangeMin), hMcProj->FindBin(rangeMax)));
        double chi2 = 0.0;
        int ndf = 0;
        for (int bin = 1; bin <= hDataProj->GetNbinsX(); ++bin) {
          // only consider bins in which the pT is in the range of interest
          double binCenter = hDataProj->GetBinCenter(bin);
          if (binCenter < rangeMin || binCenter > rangeMax) continue;
          // in this case the data histogram is the expected value and the MC histogram is the measured value
          double dataValue = hDataProj->GetBinContent(bin);
          double dataError = hDataProj->GetBinError(bin);
          double mcValue = hMcProj->GetBinContent(bin);
          double mcError = hMcProj->GetBinError(bin);
          // calculate chi2
          chi2 += std::pow(mcValue - dataValue, 2) / (std::pow(dataError, 2));
          ndf++;
        }
        chi2 /= ndf;
        chi2Values.push_back(chi2);
        // std::cout << "R = " << rValues[j] << ", chi2/ndf = " << chi2 << std::endl;
      }
      // make a graph of chi2 values vs R values
      TGraph *gChi2 = new TGraph(rValues.size());
      for (size_t j = 0; j < rValues.size(); ++j) {
        gChi2->SetPoint(j, rValues[j], chi2Values[j]);
      }

      // fit the graph with a "parabola" to find the minimum
      TF1 *fFit = new TF1("fFit", "pol6", 7.0, 9.0);
      fFit->SetNpx(10000);
      gChi2->Fit(fFit, "EMR0QS+");
      double optimalR = fFit->GetMinimumX();
      
      TCanvas *c = new TCanvas(Form("c_bin%d", i), "c", 800, 500);
      c->SetLogy();
      gChi2->Draw("AP");
      gChi2->SetTitle(";#it{R} (fm);#it{#chi}^{2}/ndf");
      gChi2->SetMarkerStyle(20);
      gChi2->SetMarkerSize(0.75);
      gChi2->SetMarkerColor(tabBlue);
      gChi2->SetMinimum(fFit->GetMinimum()-0.5);
      TLine *line = new TLine(optimalR, fFit->GetMinimum() - 10, optimalR, fFit->GetMinimum() + 10);
      line->SetLineColorAlpha(tabGreen, 0.7);
      fFit->SetLineColorAlpha(tabOrange, 1.0);
      fFit->SetLineWidth(2);
      fFit->Draw("same");
      gChi2->Draw("P same");
      // draw a rectangle around the minimum in the inset to indicate the +-1 sigma range
      double minimumY = fFit->GetMinimum();
      double minimumX = optimalR;
      double sigmaY = minimumY + 1.0; // for chi2, the 1-sigma range corresponds to chi2 + 1
      double sigmaXLeft = fFit->GetX(sigmaY, minimumX - 0.4, minimumX);
      // if it is NaN, set it to zero
      if (std::isnan(sigmaXLeft)) sigmaXLeft = 0.0;
      double sigmaXRight = fFit->GetX(sigmaY, minimumX, minimumX + 0.4);
      if (std::isnan(sigmaXRight)) sigmaXRight = 0.0;
      TBox *box = new TBox(sigmaXLeft, minimumY, sigmaXRight, sigmaY);
      box->SetFillStyle(0);
      box->SetLineColorAlpha(tabGreen, 0.7);
      box->SetLineWidth(2);
      box->Draw("same");
      c->cd();
      TLegend *legendChi2 = new TLegend(0.55, 0.7, 0.8, 0.89);
      legendChi2->SetBorderSize(0);
      legendChi2->SetFillStyle(0);
      legendChi2->AddEntry(line, Form("optimal #it{R} = %.3f ^{+%.3f}_{-%.3f} fm", optimalR, sigmaXRight - optimalR, optimalR - sigmaXLeft), "l");
      legendChi2->Draw();

      fOptimalR->cd(neutronClass);
      c->Write(Form("cChi2_bin%d", i));

      hOptimalR->SetBinContent(i + 1, optimalR);
      hOptimalR->SetBinError(i + 1, std::max(optimalR - sigmaXLeft, sigmaXRight - optimalR));

    }
    fOptimalR->cd(neutronClass);
    hOptimalR->Write("hOptimalR");
    
    // now we fit the optimal R values as a function of delta phi with the fourier decomp function
    TF1 *fFourier = new TF1("fFourier", fourierDecomp, -TMath::Pi(), TMath::Pi(), 5);
    fFourier->SetParNames("#it{c}", "#it{a}_{1}", "#it{a}_{2}", "#it{a}_{3}", "#it{a}_{4}");
    fFourier->SetParameters(0.1, 0.01, 0.01, 0.01, 0.01);
    fFourier->SetNpx(10000);

    hOptimalR->Fit(fFourier, "EMR0QS+");
    TCanvas *cFit = new TCanvas("cFit", "cFit", 800, 500);
    hOptimalR->Draw();
    hOptimalR->SetTitle(";#Delta#it{#phi} (rad);optimal #it{R} (fm)");
    hOptimalR->SetMarkerStyle(20);
    hOptimalR->SetMarkerSize(0.75);
    fFourier->SetLineColorAlpha(tabOrange, 1.0);
    fFourier->SetLineWidth(2);
    fFourier->Draw("same");
    cFit->Write("cFit");

    // what is to be saved is the bin average of the fitted function
    TH1D *hAveragedR = new TH1D("hAveragedR", ";#Delta#it{#phi} (rad); averaged optimal #it{R} (fm)", nBinsDeltaPhi, -TMath::Pi(), TMath::Pi());
    for (int i = 1; i <= hAveragedR->GetNbinsX(); ++i) {
      double binLow = hAveragedR->GetXaxis()->GetBinLowEdge(i);
      double binHigh = hAveragedR->GetXaxis()->GetBinUpEdge(i);
      double binWidth = binHigh - binLow;
      double averagedR = fFourier->Integral(binLow, binHigh) / binWidth;
      hAveragedR->SetBinContent(i, averagedR);
    }
    hAveragedR->Write();
  }
  fOptimalR->Close();
}