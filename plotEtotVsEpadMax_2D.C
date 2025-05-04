#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>

TGraph* GetGraph(const TString inFileName);

void plotEtotVsEpadMax(const TString inFileName)
{
  TGraph* EtotVsEpadMax = GetGraph(inFileName);

  if(!EtotVsEpadMax) {
    std::cerr << "Error retrieving TGraph." << std::endl;
    return;
  }

  // Define the ranges
  Double_t xMin = 1e-2;
  Double_t xMax = 2.0;
  Double_t yMin = 0.2e-2;
  Double_t yMax = 0.3;

  // Number of bins
  Int_t nBinsX = 150;
  Int_t nBinsY = 150;

  // Create logarithmically spaced bin edges for X
  // We want uniform spacing in log10 space
  Double_t logXmin = TMath::Log10(xMin);
  Double_t logXmax = TMath::Log10(xMax);
  Double_t dlogX = (logXmax - logXmin) / nBinsX;
  std::vector<Double_t> xEdges(nBinsX+1);
  for (int i = 0; i <= nBinsX; ++i) {
    xEdges[i] = TMath::Power(10, logXmin + i * dlogX);
  }

  // Similarly for Y if you also want exponential/log spacing for Y:
  // (If you only want uniform linear spacing in Y, skip this and do a linear array)
  Double_t logYmin = TMath::Log10(yMin);
  Double_t logYmax = TMath::Log10(yMax);
  Double_t dlogY = (logYmax - logYmin) / nBinsY;
  std::vector<Double_t> yEdges(nBinsY+1);
  for (int j = 0; j <= nBinsY; ++j) {
    yEdges[j] = TMath::Power(10, logYmin + j * dlogY);
  }

  // Create a TH2D with variable bin widths based on these edges
  TH2D* h2 = new TH2D("h2","Density Plot;E_{summed}^{pad} (MeV);E_{max}^{pad} (MeV)",
                      nBinsX, &xEdges[0], nBinsY, &yEdges[0]);

  // Fill the TH2D with points from the TGraph
  Int_t nPoints = EtotVsEpadMax->GetN();
  Double_t *x = EtotVsEpadMax->GetX();
  Double_t *y = EtotVsEpadMax->GetY();
  for (Int_t i = 0; i < nPoints; i++) {
    h2->Fill(x[i], y[i]);
  }

  gStyle->SetPalette(kBird);
  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas("c1","Density Plot",800,800);
  c1->SetLogx();
  c1->SetLogy();
  
  if (gPad)
  {
    gPad->SetLeftMargin(.18);
    gPad->SetRightMargin(.18);
    gPad->SetBottomMargin(.18);
    gPad->SetTopMargin(.09);  
  }

  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->GetXaxis()->SetTitleOffset(1.2);
  h2->GetYaxis()->SetTitleOffset(1.5);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetXaxis()->SetTitleSize(.055);
  h2->GetYaxis()->SetTitleSize(.055);

  // Draw the 2D histogram as a colz plot
  h2->Draw("colz");

  c1->Print("dEvsE_map_2D_logBins.png");
}

TGraph* GetGraph(const TString inFileName)
{
  TFile *file = TFile::Open(inFileName);
  
  if (!file || file->IsZombie()) 
  {
    std::cerr << "Error: File could not be opened!" << std::endl;
    return nullptr;
  }
  
  file->cd();
  TGraph* graph = (TGraph*)file->Get("EtotvsEpadMax");
  file->Close();

  return graph;
}

