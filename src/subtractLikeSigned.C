#include "commonFunctions.h"

TH2D* subtractBackground(TH2D *hUnlike, TH2D *hPos, TH2D *hNeg) {
  // clone the unlike-sign histogram
  TH2D *hSub = (TH2D*)hUnlike->Clone("hSub");

  // loop over all bins in the 2D histogram
  for (int xBin = 1; xBin <= hSub->GetNbinsX(); ++xBin) {
    for (int yBin = 1; yBin <= hSub->GetNbinsY(); ++yBin) {
      // get bin contents and errors
      double nUnlike = hUnlike->GetBinContent(xBin, yBin);
      double nPos = hPos->GetBinContent(xBin, yBin);
      double nNeg = hNeg->GetBinContent(xBin, yBin);
      double errUnlike = hUnlike->GetBinError(xBin, yBin);
      double errPos = hPos->GetBinError(xBin, yBin);
      double errNeg = hNeg->GetBinError(xBin, yBin);

      // calculate the subtracted value
      double nSub = nUnlike - 4 * TMath::Sqrt(nPos * nNeg);
      
      // propagate errors if all values are positive
      double errSub = 0.0;
      if (nUnlike > 0 && nPos > 0 && nNeg > 0) {
        errSub = TMath::Sqrt(
          errUnlike * errUnlike +
          4 * nNeg / nPos * errPos * errPos +
          4 * nPos / nNeg * errNeg * errNeg
        );
      }

      // handle cases where nPos or nNeg is zero
      if (nPos == 0 || nNeg == 0) {
        nSub = nUnlike;
        errSub = errUnlike;
      }

      // set the bin content and error in the subtracted histogram
      hSub->SetBinContent(xBin, yBin, nSub);
      hSub->SetBinError(xBin, yBin, errSub);
    }
  }

  return hSub;
}

void subtractLikeSigned() {
  globalSetup();
  // open all necessary files
  TFile *fData = new TFile("input/data.root", "READ");
  TFile *fOut = new TFile("output/subtracted.root", "RECREATE");

  for (int i = 0; i < neutronClasses.size(); i++) {
    for (int j = 0; j < nBinsPt; j++) {
      TH2D *hUnlike = (TH2D*)fData->Get(Form("%s/unlike/pTbin_%d/hDeltaPhiVsM", neutronClasses[i].Data(), j));
      if (hUnlike == nullptr) {
        std::cout << "Histogram " << Form("%s/unlike/pTbin_%d/hDeltaPhiVsM", neutronClasses[i].Data(), j) << " not found!" << std::endl;
        continue;
      }
      TH2D *hPos = (TH2D*)fData->Get(Form("%s/positive/pTbin_%d/hDeltaPhiVsM", neutronClasses[i].Data(), j));
      if (hPos == nullptr) {
        std::cout << "Histogram " << Form("%s/positive/pTbin_%d/hDeltaPhiVsM", neutronClasses[i].Data(), j) << " not found!" << std::endl;
        continue;
      }
      TH2D *hNeg = (TH2D*)fData->Get(Form("%s/negative/pTbin_%d/hDeltaPhiVsM", neutronClasses[i].Data(), j));
      if (hNeg == nullptr) {
        std::cout << "Histogram " << Form("%s/negative/pTbin_%d/hDeltaPhiVsM", neutronClasses[i].Data(), j) << " not found!" << std::endl;
        continue;
      }
      // subtract combinatorial background from the unlike-sign histogram
      TH2D *hSub = subtractBackground(hUnlike, hPos, hNeg);
      // write the subtracted histograms to the output file
      TDirectory *dir = fOut->mkdir(Form("%s/pTbin_%d", neutronClasses[i].Data(), j));
      dir->cd();
      hSub->Write("hDeltaPhiVsM");
    }
  }
  fOut->Close();

  // // plot one example
  // globalSetup();
  // gSystem->Exec("mkdir -p plots");
  // TH1D *hMunlike = (TH1D*)fData->Get("unlike/AnAn/hM");
  // TH1D *hMpositive = (TH1D*)fData->Get("positive/AnAn/hM");
  // TH1D *hMnegative = (TH1D*)fData->Get("negative/AnAn/hM");
  // // do the subtraction for 1D histogram
  // TH1D *hM_subtracted = (TH1D*)hMunlike->Clone("hM_subtracted");
  // hM_subtracted->SetName("hM_subtracted");
  // for (int i = 1; i <= hM_subtracted->GetNbinsX(); i++) {
  //   double nUnlike = hMunlike->GetBinContent(i);
  //   double nPos = hMpositive->GetBinContent(i);
  //   double nNeg = hMnegative->GetBinContent(i);
  //   double errUnlike = hMunlike->GetBinError(i);
  //   double errPos = hMpositive->GetBinError(i);
  //   double errNeg = hMnegative->GetBinError(i);

  //   // Calculate the subtracted value
  //   double nSub = nUnlike - 4 * TMath::Sqrt(nPos * nNeg);
  //   double errSub = 0.0;

  //   // Propagate errors if all values are positive
  //   if (nUnlike > 0 && nPos > 0 && nNeg > 0) {
  //     errSub = TMath::Sqrt(
  //       errUnlike * errUnlike +
  //       4 * nNeg / nPos * errPos * errPos +
  //       4 * nPos / nNeg * errNeg * errNeg
  //     );
  //   }

  //   // Handle cases where nPos or nNeg is zero
  //   if (nPos == 0 || nNeg == 0) {
  //     nSub = nUnlike;
  //     errSub = errUnlike;
  //   }

  //   // Set the bin content and error in the subtracted histogram
  //   hM_subtracted->SetBinContent(i, nSub);
  //   hM_subtracted->SetBinError(i, errSub);
  // }
  // // plot the subtracted histogram
  // TCanvas *c = new TCanvas("c", "c", 1000, 600);
  // c->SetLogy();
  // hMunlike->SetMinimum(hMpositive->GetMinimum());
  // hMunlike->SetMaximum(2e8);
  // hMunlike->SetLineColor(tabOrange);
  // // hMunlike->SetLineWidth(2);
  // hMunlike->SetMarkerColor(tabOrange);
  // hMpositive->SetLineColor(tabGreen);
  // // hMpositive->SetLineWidth(2);
  // hMpositive->SetMarkerColor(tabGreen);
  // hMnegative->SetLineColor(tabRed);
  // // hMnegative->SetLineWidth(2);
  // hMnegative->SetMarkerColor(tabRed);
  // hM_subtracted->SetLineColor(tabBlue);
  // // hM_subtracted->SetLineWidth(2);
  // hM_subtracted->SetMarkerColor(tabBlue);

  // hMunlike->GetXaxis()->SetTitle("#it{m}_{#pi#pi} (GeV/#it{c}^{2})");
  // hMunlike->Draw("E");
  // hMpositive->Draw("E SAME");
  // hMnegative->Draw("E SAME");
  // hM_subtracted->Draw("E SAME");

  // // make header
  // TLegend* header = new TLegend(0.75, 0.9, 0.9, 0.95);
  // header->SetTextAlign(32);
  // header->AddEntry((TObject*)0, "this work", "");
  // header->SetBorderSize(0);
  // header->SetFillStyle(0);
  // // header->Draw();

  // // make legend
  // TLegend* stats = new TLegend(0.5, 0.59, 0.9, 0.89);
  // stats->SetTextAlign(32);
  // stats->SetBorderSize(0);
  // stats->SetFillStyle(0);
  // stats->AddEntry((TObject*)0, "ALICE Pb#font[122]{-}Pb UPC #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
  // stats->AddEntry((TObject*)0, "UD_LHC23_pass4_SingleGap", "");
  // stats->AddEntry((TObject*)0, "#rho^{0} #rightarrow #pi^{+}#pi^{-}", "");
  // // stats->AddEntry((TObject*)0, "m #in (0.5, 1.2) GeV/#it{c}^{2}", "");
  // stats->AddEntry((TObject*)0, "#it{p}_{T} #leq 0.1 GeV/#it{c}", "");
  // stats->AddEntry((TObject*)0, "| #it{y}| #leq 0.9", "");
  // stats->AddEntry((TObject*)0, "this thesis", "");
  // stats->Draw();
  // TLegend* legend = new TLegend(0.125, 0.7, 0.4, 0.89);
  // // legend->SetTextAlign(32);
  // legend->SetTextSize(stats->GetTextSize());
  // legend->AddEntry(hMunlike, "unlike-sign pairs", "le");
  // legend->AddEntry(hMpositive, "positive pairs", "le");
  // legend->AddEntry(hMnegative, "negative pairs", "le");
  // legend->AddEntry(hM_subtracted, "background-subtracted", "le");
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);
  // legend->Draw();

  // c->SaveAs("plots/hMsubtracted.pdf");

  // // make plot of ratios
  // TGaxis::SetMaxDigits(2);

  // TCanvas *c2 = new TCanvas("c2", "c2", 1000, 600);
  // c2->Divide(2,1);
  // // adjust pad margins
  // c2->cd(1);
  // c2->GetPad(1)->SetLeftMargin(0.13);
  // c2->GetPad(1)->SetRightMargin(0.05);
  // // c2->GetPad(1)->SetBottomMargin(0.15);
  // c2->cd(2);
  // c2->GetPad(2)->SetLeftMargin(0.13);
  // c2->GetPad(2)->SetRightMargin(0.05);
  // // c2->GetPad(2)->SetBottomMargin(0.15);

  // TH1D *hMpositive_ratio = (TH1D*)hMpositive->Clone("hMpositive_ratio");
  // hMpositive_ratio->Divide(hMpositive, hMunlike, 1, 1, "B");
  // hMpositive_ratio->SetLineColor(tabGreen);
  // hMpositive_ratio->SetMarkerColor(tabGreen);
  // TH1D *hMnegative_ratio = (TH1D*)hMnegative->Clone("hMnegative_ratio");
  // hMnegative_ratio->Divide(hMnegative, hMunlike, 1, 1, "B");
  // hMnegative_ratio->SetLineColor(tabRed);
  // hMnegative_ratio->SetMarkerColor(tabRed);
  // TH1D *hM_subtracted_ratio = (TH1D*)hM_subtracted->Clone("hM_subtracted_ratio");
  // hM_subtracted_ratio->Divide(hM_subtracted, hMunlike, 1, 1, "B");
  // hM_subtracted_ratio->SetLineColor(tabBlue);
  // hM_subtracted_ratio->SetMarkerColor(tabBlue);
  // hMpositive_ratio->GetXaxis()->SetTitle("#it{m}_{#pi#pi} (GeV/#it{c}^{2})");
  // hMnegative_ratio->GetXaxis()->SetTitle("#it{m}_{#pi#pi} (GeV/#it{c}^{2})");
  // hM_subtracted_ratio->GetXaxis()->SetTitle("#it{m}_{#pi#pi} (GeV/#it{c}^{2})");
  // hMpositive_ratio->GetYaxis()->SetTitle("ratio to reconstructed unlike-sign pairs");
  // hMnegative_ratio->GetYaxis()->SetTitle("ratio to reconstructed unlike-sign pairs");
  // hM_subtracted_ratio->GetYaxis()->SetTitle("ratio to reconstructed unlike-sign pairs");
  
  // // plot all
  // c2->cd(1);
  // hMnegative_ratio->Draw("E");
  // hMpositive_ratio->Draw("E SAME");
  // // hMpositive_ratio->GetYaxis()->SetRangeUser(2e-4, 1);

  // // make legend
  // TLegend* legend2 = new TLegend(0.15, 0.79, 0.5, 0.89);
  // // legend2->SetTextAlign(32);
  // legend2->SetTextSize(stats->GetTextSize());
  // legend2->AddEntry(hMpositive_ratio, "positive pairs", "le");
  // legend2->AddEntry(hMnegative_ratio, "negative pairs", "le");
  // // legend2->AddEntry(hM_subtracted_ratio, "background-subtracted", "le");
  // legend2->SetBorderSize(0);
  // legend2->SetFillStyle(0);
  // legend2->Draw();

  // c2->cd(2);
  // hM_subtracted_ratio->Draw("E");
  // TLegend* legend3 = new TLegend(0.15, 0.83, 0.5, 0.89);
  // // legend3->SetTextAlign(32);
  // legend3->SetTextSize(stats->GetTextSize());
  // legend3->AddEntry(hM_subtracted_ratio, "subtracted", "le");
  // legend3->SetBorderSize(0);
  // legend3->SetFillStyle(0);
  // legend3->Draw();

  // c2->SaveAs("plots/hMsubtracted_ratio.pdf");
}