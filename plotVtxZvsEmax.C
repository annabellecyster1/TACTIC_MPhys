void plotEpadMaxVsZ(const char* filename = "outputFile.root")
{
  // Set style options if desired
  gStyle->SetOptStat(0); // No stats box
  gStyle->SetPalette(kBird); // Nice color palette

  // Open the ROOT file specified by the user
  TFile* f = TFile::Open(filename, "READ");
  if (!f || !f->IsOpen()) {
    std::cerr << "Could not open file: " << filename << std::endl;
    return;
  }

  // Retrieve the TH2D histogram
  TH2D* h2 = dynamic_cast<TH2D*>(f->Get("VtxZvsEmax")); 
  if (!h2) {
    std::cerr << "Histogram 'VtxZvsEmax' not found in file: " << filename << std::endl;
    f->Close();
    return;
  }

  // Remove white lines between bins
  h2->SetLineColor(kNone);
  h2->SetLineWidth(0);

  // No log scale on Y-axis this time
  // We can set a linear range as needed
  // For example, from 0.0 to 0.2 MeV:
  //h2->GetYaxis()->SetRangeUser(0.0, 0.2);

  // Create a canvas to draw the histogram
  TCanvas* c = new TCanvas("c", "Z vs EpadMax", 800, 600);
  c->SetRightMargin(0.15); // give space for z-axis color scale
  c->SetLogy();
  c->SetGrid(0,0);  // No grid lines

  // Draw the histogram with a color scale
  h2->Draw("COLZ");
  h2->GetXaxis()->SetRangeUser(-40, 60.);
  h2->GetYaxis()->SetRangeUser(0, 0.15);


  // Axis titles
  h2->GetXaxis()->SetTitle("Vertex Z [mm]");
  h2->GetYaxis()->SetTitle("Max Pad Energy [MeV]");
  h2->GetZaxis()->SetTitle("Counts");

  // Update the canvas
  c->Update();

}

