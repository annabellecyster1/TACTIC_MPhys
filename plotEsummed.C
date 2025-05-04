#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>

#include <iostream>

TH1D* GetHistogram(const char* fileName, const char* histName);

void plotEsummed(const TString inFileName)
{
  // Get the histogram
  TH1D* hist = GetHistogram(inFileName, "EsummedHist");

  if (!hist)
  {
    return;
  }

  // Create a canvas
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 900);

  if (gPad)
  {
    gPad->SetLeftMargin(.145);
    gPad->SetBottomMargin(.15);
    gPad->SetTopMargin(.1);
  }
  
  hist->SetStats(false);

  // Apply axis decorations
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlue);
  hist->GetXaxis()->SetRangeUser(0., 7.);
  // hist->GetYaxis()->SetRangeUser(1, 200000);
  hist->GetXaxis()->SetTitle("E_{summed} (keV)");
  hist->GetYaxis()->SetTitle("Count");
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetLabelSize(.04);
  hist->GetYaxis()->SetLabelSize(.04);
  hist->GetXaxis()->SetTitleSize(.045);
  hist->GetYaxis()->SetTitleSize(.045);

  // Draw the histogram
  hist->Draw();

  // Save the canvas as an image (optional)
  c1->SetLogy();
  c1->SaveAs("Esummed_55Fe.png");

  // Clean up
  delete hist;
  delete c1;
}

TH1D* GetHistogram(const char* fileName, const char* histName)
{
  // Open the ROOT file
  TFile* file = TFile::Open(fileName, "READ");

  if (!file || file->IsZombie())
  {
    std::cerr << "Error: Unable to open file " << fileName << "!"
              << std::endl;
    return nullptr;
  }

  // Retrieve the histogram
  TH1D* hist = dynamic_cast<TH1D*>(file->Get(histName));

  if (!hist)
  {
    std::cerr << "Error: Histogram '" << histName 
              << "' not found in file " << fileName << "!" 
              << std::endl;
    file->Close();
    delete file;
    return nullptr;
  }

  // Clone the histogram to decouple it from the file
  TH1D* histClone = dynamic_cast<TH1D*>(hist->Clone());
  histClone->SetDirectory(0); // Detach from the file's memory

  // Close the file
  file->Close();

  return histClone;
}