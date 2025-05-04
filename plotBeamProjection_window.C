void plot_beam_projection_y2(const char* filename = "processedData.root") {
    // Open the specified ROOT file
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: couldn't open file " << filename << std::endl;
        return;
    }

    // Load the second version of the histogram
    TH1D *h = (TH1D*)f->Get("BeamProjectionY;2");
    if (!h) {
        std::cerr << "Error: Histogram BeamProjectionY;2 not found in file " << filename << std::endl;
        return;
    }

    // Set axis labels and title
    h->SetTitle("Beam Energy Profile at 20 MeV;Beam Energy Profile (MeV);Counts");

    // Style options (optional)
    h->SetLineColor(kBlue + 1);
    h->SetLineWidth(2);
    h->GetXaxis()->SetRangeUser(12, 20);

    gStyle->SetOptStat(1110); // Show stats box with Name, Entries, Mean, RMS

    // Define Gaussian fit function
    TF1 *gausFit = new TF1("gausFit", "gaus", 14, 18); // Fit range 18-22
    gausFit->SetLineColor(kRed);  // Red color for the Gaussian fit

    // Perform the fit in the range 18 to 22
    h->Fit(gausFit, "R");

    // Extract the standard deviation (sigma) from the fit
    double sigma = gausFit->GetParameter(2); // Standard deviation (sigma)

    // Calculate the Full Width at Half Maximum (FWHM)
    double fwhm = 2 * TMath::Sqrt(2 * TMath::Log(2)) * sigma;

    // Add FWHM to the statistics box manually
    TString stats = TString::Format("FWHM = %.2f MeV", fwhm);
    
    // Find the statistics box and add the FWHM
    TPaveStats* statsBox = (TPaveStats*)h->FindObject("stats");
    if (statsBox) {
        statsBox->AddText(stats);  // Add FWHM information
        statsBox->SetX1NDC(0.7);   // Position stats box (optional, adjust as needed)
        statsBox->SetY1NDC(0.7);   // Position stats box (optional, adjust as needed)
    }

    // Create and draw on canvas
    TCanvas *c = new TCanvas("c", "Beam Projection Y", 800, 600);
    h->Draw("HIST");

    // Draw the Gaussian fit on top of the histogram
    gausFit->Draw("SAME");

   
}

