// header file with common variables and functions for the analysis macros

#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

// C++ includes
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <iostream>
#include <regex>
#include <chrono>
#include <memory>
#include <algorithm>

// standard ROOT headers
#include <TCanvas.h>
#include <TColor.h>
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Math/MinimizerOptions.h>
#include <TDatabasePDG.h>
#include <TGaxis.h>
#include <TComplex.h>
#include <THStack.h>
#include <TRandom.h>
#include <TTree.h>
#include <TKey.h>

// headers for RDataFrame stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/TThreadExecutor.hxx>
#include <ROOT/RResultPtr.hxx>

// settable constants for analysis
const double nSigmaCut = 3.0;

const double minMass = 0.5;
const double maxMass = 1.0;
const int nBinsMass = 50;

const double fitMin = 0.60;
const double fitMax = 1.00;

const double maxY = 0.9;
const int nBinsY = 180;

const int nBinsDeltaPhi = 12;

const bool useDeltaPhiRandom = false;
const std::string deltaPhiLabel = useDeltaPhiRandom ? "#Delta#it{#phi}_{random}" : "#Delta#it{#phi}_{charge}";

const std::vector<double> pTbinEdges = {0.0, 0.1};
const int nBinsPt = pTbinEdges.size() - 1;

const bool separateXn0n = true;

const bool fixPoles = true;

std::vector<TString> neutronClasses;
if (separateXn0n) neutronClasses = {"AnAn", "0n0n", "Xn0n", "0nXn", "XnXn"};
else neutronClasses = {"AnAn", "0n0n", "Xn0n", "XnXn"};
 
// mimic matplotlib tableau colours
const auto tabBlue = TColor::GetColor("#1f77b4");
const auto tabOrange = TColor::GetColor("#ff7f0e");
const auto tabGreen = TColor::GetColor("#2ca02c");
const auto tabRed = TColor::GetColor("#d62728");
const auto tabPurple = TColor::GetColor("#9467bd");
const auto tabBrown = TColor::GetColor("#8c564b");
const auto tabPink = TColor::GetColor("#e377c2");
const auto tabGray = TColor::GetColor("#7f7f7f");
const auto tabOlive = TColor::GetColor("#bcbd22");
const auto tabCyan = TColor::GetColor("#17becf");
// gather these into a vector for easy access
const std::vector<int> tabColours = {tabBlue, tabOrange, tabGreen, tabRed, tabPurple, tabBrown, tabPink, tabGray, tabOlive, tabCyan};

// for reweighting
std::vector<double> rValues = {7.00, 7.10, 7.20, 7.30, 7.40, 7.50, 7.60, 7.70, 7.80, 7.90, 8.00, 8.10, 8.20, 8.30, 8.40, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00};

// constants for fitting
// auto PDG = TDatabasePDG::Instance();
const double kMpi = 0.13957039; //PDG->GetParticle(211)->Mass();
const double kM2pi = kMpi * kMpi;
const double kMrho = 0.7692; //PDG->GetParticle(113)->Mass();
const double kWrho = 0.1515; //PDG->GetParticle(113)->Width();
const double kMomega = 0.78266; //PDG->GetParticle(223)->Mass();
const double kWomega = 0.00868; //PDG->GetParticle(223)->Width();

// helper struct to hold filter specifications
struct FilterSpec {
  std::string name;
  std::string expression;
};

// miscellaneous set-up to unify plot graphics and such
void globalSetup() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetStripDecimals(0); // equivalent to SetDecimals
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(999); // smoothen the colour gradient for Z axes
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1e6); // increase possible number of iterations for fitting
  TGaxis::SetMaxDigits(4);
  TH1::SetDefaultSumw2(kTRUE); // enable proper error calculation for histograms by default
}

// reconstruct rho LV
TLorentzVector getRho(const float leadingPt, const float subleadingPt,
                      const float leadingEta, const float subleadingEta,
                      const float leadingPhi, const float subleadingPhi) {
  TLorentzVector lv1, lv2;
  lv1.SetPtEtaPhiM(leadingPt, leadingEta, leadingPhi, kMpi);
  lv2.SetPtEtaPhiM(subleadingPt, subleadingEta, subleadingPhi, kMpi);

  return lv1 + lv2;
}

// calculate Delta phi
double getDeltaPhi(const int leadingSign, const int subleadingSign,
                  const float leadingPt, const float subleadingPt,
                  const float leadingEta, const float subleadingEta,
                  const float leadingPhi, const float subleadingPhi) {
  const bool useLeadingFirst = useDeltaPhiRandom ? (gRandom->Uniform() < 0.5) : (leadingSign > 0);
  TLorentzVector lv1, lv2;

  if (useLeadingFirst) {
    lv1.SetPtEtaPhiM(leadingPt, leadingEta, leadingPhi, kMpi);
    lv2.SetPtEtaPhiM(subleadingPt, subleadingEta, subleadingPhi, kMpi);
  } else {
    lv1.SetPtEtaPhiM(subleadingPt, subleadingEta, subleadingPhi, kMpi);
    lv2.SetPtEtaPhiM(leadingPt, leadingEta, leadingPhi, kMpi);
  }

  const TLorentzVector lvPlus = lv1 + lv2;
  const TLorentzVector lvMinus = lv1 - lv2;
  return lvPlus.DeltaPhi(lvMinus);
}

// dark magic to find the correct path to the tree without knowing the DF number in advance
TString findTreePath(const TString& path, const TString& tree) {
  // find the ROOT file
  TFile file(path, "READ");
  if (!file.IsOpen()) {
    std::cerr << "ERROR: Cannot open file " << path << std::endl;
    return "";
  }
  std::cout << "Opened file: " << path << std::endl;

  // recursive lambda function to search for the tree in the directory structure
  auto searchDirectory = [&](auto&& self, TDirectory* dir, const TString& prefix, TString& foundPath) -> bool {
    TIter next(dir->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = static_cast<TKey*>(next()))) {
      TString name = key->GetName();
      const std::string className = key->GetClassName();
      const TString currentPath = prefix.IsNull() ? name : prefix + "/" + name;

      if (name == tree && className == "TTree") {
        foundPath = currentPath;
        return true;
      }

      if (className == "TDirectory" || className == "TDirectoryFile") {
        TObject* obj = key->ReadObj();
        if (auto subDir = dynamic_cast<TDirectory*>(obj)) {
          if (self(self, subDir, currentPath, foundPath)) {
            return true;
          }
        }
      }
    }
    return false;
  };

  TString foundPath;
  if (searchDirectory(searchDirectory, &file, "", foundPath)) {
    return foundPath;
  }

  std::cerr << "ERROR: Cannot find tree matching " << tree << " in file " << path << std::endl;
  return "";
}

// decide which neutron classes are to be used
std::vector<FilterSpec> getNeutronFilters(bool isMC) {
  // no ZDC in MC, so only AnAn
  if (isMC) {
    return {{"AnAn", "(fNeutronClass == -1) || (fNeutronClass == 0) || (fNeutronClass == 1) || (fNeutronClass == 2) || (fNeutronClass == 3)"}};
  }

  if (separateXn0n) {
    return {
      {"AnAn", "(fNeutronClass == -1) || (fNeutronClass == 0) || (fNeutronClass == 1) || (fNeutronClass == 2) || (fNeutronClass == 3)"},
      {"0n0n", "fNeutronClass == 0"},
      {"Xn0n", "fNeutronClass == 1"},
      {"0nXn", "fNeutronClass == 2"},
      {"XnXn", "fNeutronClass == 3"}
    };
  }

  // if not separating Xn0n and 0nXn, combine them into one class (Xn0n)
  return {
    {"AnAn", "(fNeutronClass == -1) || (fNeutronClass == 0) || (fNeutronClass == 1) || (fNeutronClass == 2) || (fNeutronClass == 3)"},
    {"0n0n", "fNeutronClass == 0"},
    {"Xn0n", "(fNeutronClass == 1) || (fNeutronClass == 2)"},
    {"XnXn", "fNeutronClass == 3"}
  };
}

// similar thing but for filtering on charge
std::vector<FilterSpec> getChargeFilters(bool isMC) {
  // no like-signed events in MC
  if (isMC) {
    return {{"unlike", "fLeadingTrackSign + fSubleadingTrackSign == 0"}};
  }

  // positive and negative used for corrections
  return {
    {"unlike", "fLeadingTrackSign + fSubleadingTrackSign == 0"},
    {"positive", "fLeadingTrackSign + fSubleadingTrackSign == 2"},
    {"negative", "fLeadingTrackSign + fSubleadingTrackSign == -2"}
  };
}

// helper function to create a directory hierarchy in an output ROOT file
TDirectory* createDirectoryHierarchy(TDirectory* parent, const std::vector<std::string>& names) {
  TDirectory* current = parent;
  for (const auto& name : names) {
    TDirectory* child = current->GetDirectory(name.c_str());
    if (!child) child = current->mkdir(name.c_str());
    current = child;
  }
  return current;
}

// functions for fitting
double Gamma_rho(double pipi, double mass, double width) {
  return width*(mass/pipi)*TMath::Power((pipi*pipi-4.0*kM2pi)/(mass*mass-4.0*kM2pi), 3./2.);
}

double Gamma_omega(double pipi, double mass, double width) {
  return width*(mass/pipi)*TMath::Power((pipi*pipi-9.0*kM2pi)/(mass*mass-9.0*kM2pi), 3./2.);
}

double Gamma_omega_pion(double pipi, double mass, double width) {
  return 0.0153 * Gamma_rho(pipi, mass, width); // neglecting errors for the BR
}

double x1(double *x, double *par) {
  double m_pipi = x[0];
  double A_rho = par[0];
  double m_rho = par[1];
  double gamma0_rho = par[2];
  double a = std::sqrt(m_pipi * m_rho * Gamma_rho(m_pipi, m_rho, gamma0_rho));
  double b = m_pipi * m_pipi - m_rho * m_rho;
  double c = m_rho * Gamma_rho(m_pipi, m_rho, gamma0_rho);
  return A_rho * a * b / (b * b + c * c);
}

double x2(double *x, double *par) {
  double B_pipi = par[7];
  return B_pipi;
}

double x3(double *x, double *par) {
  double m_pipi = x[0];
  double C_omega = par[3];
  double m_omega = par[4];
  double gamma0_omega = par[5];
  double phi_omega = par[6];
  double d = std::sqrt(m_pipi * m_omega * Gamma_omega(m_pipi, m_omega, gamma0_omega));
  double e = m_pipi * m_pipi - m_omega * m_omega;
  double f = m_omega * Gamma_omega(m_pipi, m_omega, gamma0_omega);
  return C_omega * d * e * std::cos(phi_omega) / (e * e + f * f);
}

double x4(double *x, double *par) {
  double m_pipi = x[0];
  double C_omega = par[3];
  double m_omega = par[4];
  double gamma0_omega = par[5];
  double phi_omega = par[6];
  double d = std::sqrt(m_pipi * m_omega * Gamma_omega(m_pipi, m_omega, gamma0_omega));
  double e = m_pipi * m_pipi - m_omega * m_omega;
  double f = m_omega * Gamma_omega(m_pipi, m_omega, gamma0_omega);
  return C_omega * d * f * std::sin(phi_omega) / (e * e + f * f);
}
  
double y1(double *x, double *par) {
  double m_pipi = x[0];
  double A_rho = par[0];
  double m_rho = par[1];
  double gamma0_rho = par[2];
  double a = std::sqrt(m_pipi * m_rho * Gamma_rho(m_pipi, m_rho, gamma0_rho));
  double b = m_pipi * m_pipi - m_rho * m_rho;
  double c = m_rho * Gamma_rho(m_pipi, m_rho, gamma0_rho);
  return -1.0 * A_rho * a * c / (b * b + c * c);
}

double y2(double *x, double *par) {
  double m_pipi = x[0];
  double C_omega = par[3];
  double m_omega = par[4];
  double gamma0_omega = par[5];
  double phi_omega = par[6];
  double d = std::sqrt(m_pipi * m_omega * Gamma_omega(m_pipi, m_omega, gamma0_omega));
  double e = m_pipi * m_pipi - m_omega * m_omega;
  double f = m_omega * Gamma_omega(m_pipi, m_omega, gamma0_omega);
  return C_omega * d * e * std::sin(phi_omega) / (e * e + f * f);
}

double y3(double *x, double *par) {
  double m_pipi = x[0];
  double C_omega = par[3];
  double m_omega = par[4];
  double gamma0_omega = par[5];
  double phi_omega = par[6];
  double d = std::sqrt(m_pipi * m_omega * Gamma_omega(m_pipi, m_omega, gamma0_omega));
  double e = m_pipi * m_pipi - m_omega * m_omega;
  double f = m_omega * Gamma_omega(m_pipi, m_omega, gamma0_omega);
  return -1.0 * C_omega * d * f * std::cos(phi_omega) / (e * e + f * f);
}

double rhoContrib(double *x, double *par) {
  return std::pow(x1(x, par), 2) + std::pow(y1(x, par), 2);
}

double omegaContrib(double *x, double *par) {
  return std::pow(x3(x, par), 2) + std::pow(x4(x, par), 2) + std::pow(y2(x, par), 2) + std::pow(y3(x, par), 2) + 2 * (x3(x, par) * x4(x, par) + y2(x, par) * y3(x, par));
}

double diPionContrib(double *x, double *par) {
  return std::pow(x2(x, par), 2);
}

double rhoDiPionInterference(double *x, double *par) {
  return 2 * x1(x, par) * x2(x, par);
}

double rhoOmegaInterference(double *x, double *par) {
  return 2 * x1(x, par) * x3(x, par) + 2 * x1(x, par) * x4(x, par) + 2 * y1(x, par) * y2(x, par) + 2 * y1(x, par) * y3(x, par);
}

double omegaDiPionInterference(double *x, double *par) {
  return 2 * x2(x, par) * x3(x, par) + 2 * x2(x, par) * x4(x, par);
}

double powerBackground(double *x, double *par) {
  return par[8] * std::pow(x[0], -1.0 * par[9]);
}

// build soeding function from the individual contributions
double Soeding(double *x, double *par) {
  return rhoContrib(x, par) + omegaContrib(x, par) + diPionContrib(x, par) + rhoDiPionInterference(x, par) + rhoOmegaInterference(x, par) + omegaDiPionInterference(x, par) + powerBackground(x, par);
}

double SoedingNoBackground(double *x, double *par) {
  return rhoContrib(x, par) + omegaContrib(x, par) + diPionContrib(x, par) + rhoDiPionInterference(x, par) + rhoOmegaInterference(x, par) + omegaDiPionInterference(x, par);
}

// relativistic Breit-Wigner function for rho0
TComplex BreitWigner_rho(double *x, double *par) {
  double m_pipi = x[0];
  double m_rho = par[0];
  double gamma0_rho = par[1];
  TComplex numerator(std::sqrt(m_pipi * m_rho * Gamma_rho(m_pipi, m_rho, gamma0_rho)), 0.0);
  TComplex denominator(m_rho * m_rho - m_pipi * m_pipi, -1.0 * m_rho * Gamma_rho(m_pipi, m_rho, gamma0_rho));
  return numerator / denominator;
}

// Ross-Stodolsky function for rho0
double RossStodolsky(double *x, double *par) {
  double m_pipi = x[0];
  double m_rho = par[0];
  double gamma0_rho = par[1];
  double f = par[2];
  double k = par[3];
  return f * BreitWigner_rho(x, par).Rho2() * std::pow(m_rho / m_pipi, k);
}

// 4th order Fourier decomp of yield histogram
double fourierDecomp(double *x, double *par) {
  double deltaPhi = x[0];
  double c = par[0];
  double a1 = par[1];
  double a2 = par[2];
  double a3 = par[3];
  double a4 = par[4];
  return c * (1 + a1 * TMath::Cos(deltaPhi) + a2 * TMath::Cos(2 * deltaPhi) + a3 * TMath::Cos(3 * deltaPhi) + a4 * TMath::Cos(4 * deltaPhi));
}

#endif
