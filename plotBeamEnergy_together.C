#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <iostream>

// Usage example (in a shell):
//   root -l 'plotBeamEnergyVsZ.C("file1.root","file2.root","file3.root","file4.root","Label1","Label2","Label3","Label4")'

void plotBeamEnergyVsZ(const char* inFileName1 = "", const char* inFileName2 = "",
                       const char* inFileName3 = "", const char* inFileName4 = "",
                       const char* label1      = "", const char* label2      = "",
                       const char* label3      = "", const char* label4      = "")
{
    // Check if the file names are provided
    if (std::string(inFileName1).empty() ||
        std::string(inFileName2).empty() ||
        std::string(inFileName3).empty() ||
        std::string(inFileName4).empty())
    {
        std::cerr << "Usage: plotBeamEnergyVsZ(\"file1.root\", \"file2.root\", "
                     "\"file3.root\", \"file4.root\", \"label1\", \"label2\", "
                     "\"label3\", \"label4\")"
                  << std::endl;
        return;
    }

    // Array of file names and labels
    const char* fileNames[4] = { inFileName1, inFileName2, inFileName3, inFileName4 };
    const char* labels[4]    = { label1, label2, label3, label4 };

    // Colors for each graph (adjust as you wish)
    Color_t colors[4] = { kBlue, kRed, kGreen+2, kMagenta+1 };

    // Array of TGraph pointers
    TGraphErrors* graphs[4] = {nullptr};

    // Open each file and retrieve the TGraph
    for (int i = 0; i < 4; ++i) {
        TFile* inFile = TFile::Open(fileNames[i], "READ");
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << fileNames[i] << std::endl;
            if (inFile) inFile->Close();
            return;
        }

        graphs[i] = dynamic_cast<TGraphErrors*>(inFile->Get("beamEnergyVsZ"));
        if (!graphs[i]) {
            std::cerr << "Error: Could not find TGraphErrors 'beamEnergyVsZ' in " << fileNames[i] << std::endl;
            inFile->Close();
            return;
        }

        // Set style for each graph
        graphs[i]->SetMarkerStyle(20 + i);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);

        inFile->Close();
    }

    // Create a canvas
    TCanvas* c = new TCanvas("c", "Beam Energy vs Z", 800, 600);

    // Combine graphs into a TMultiGraph so axes auto-scale to *all* points
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle("Mean Beam Energy vs Z;Z (mm);Mean Beam Energy (MeV)");
    for (int i = 0; i < 4; ++i) {
        mg->Add(graphs[i], "PL");
    }
   
    // Draw the multigraph (the "A" option forces drawing of axes)
    mg->Draw("A");
    
    
   // TLine *line1 = new TLine(-46, 8.1, 127, 8.1);
   // line1->SetLineStyle(2);  // Dashed line
    //line1->Draw();
    //line1->SetLineColor(kRed);
    
    TLine *line84 = new TLine(84,0,84,35.3);
    line84->SetLineColor(kBlue); // Line color
    line84->SetLineStyle(1); // Dashed line style
    line84->Draw();  // Draw the line at y = 84 mm

   // TLine *line2 = new TLine(-46, 22.95, 102, 22.95);
   // line2->SetLineColor(kRed);
    //line2->SetLineStyle(2);  
    //line2->Draw();
    
    c->Update();
    gPad->Update();
    
    // Create a legend
    TLegend* legend = new TLegend(0.65, 0.7, 0.85, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(kWhite, 0);
    for (int i = 0; i < 4; ++i) {
        legend->AddEntry(graphs[i], labels[i], "lp");
    }
    legend->Draw("same");

    // Optionally save as PNG
    // c->SaveAs("beamEnergyVsZ_comparison.png");
}
