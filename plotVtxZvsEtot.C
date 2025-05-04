#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>

void plotContour(TString fileName = "output.root") {
    // Open the ROOT file containing the histogram
    TFile *inFile = new TFile(fileName, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open the ROOT file " << fileName << "!" << std::endl;
        return;
    }

    // Retrieve the 2D histogram
    TH2D *VtxZvsEtot = (TH2D*)inFile->Get("VtxZvsEtot");
    if (!VtxZvsEtot) {
        std::cerr << "Error: Histogram 'VtxZvsEtot' not found in the file!" << std::endl;
        return;
    }

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Vertex Z position vs total energy Contour", 800, 600);
    //c1->SetLogy();

    
    // Improve visualization
    gStyle->SetOptStat(0);  // Disable stats box
    gStyle->SetPalette(kBird); // Use a nice color palette
    VtxZvsEtot->GetYaxis()->SetRangeUser(0.0, 1);
    VtxZvsEtot->GetXaxis()->SetTitle("Vertex Z (mm)");
    VtxZvsEtot->GetYaxis()->SetTitle("Total energy (MeV)");

	
    // Draw the histogram as a contour plot
    VtxZvsEtot->Draw("CONTZ");  // "CONTZ" gives colored contour lines

 

    // Close the ROOT file
   // inFile->Close();
}

