void plot_multiple_beam_profiles_with_gauss_fit_and_fwhm(
    const char* filename1 = "file1.root", const char* filename2 = "file2.root", const char* filename3 = "file3.root",
    const char* label1 = "Beam Profile 1", const char* label2 = "Beam Profile 2", const char* label3 = "Beam Profile 3",
    double offset1 = 0.0, double offset2 = 0.0, double offset3 = 0.0) {

    // Open the ROOT files
    TFile *f1 = TFile::Open(filename1);
    TFile *f2 = TFile::Open(filename2);
    TFile *f3 = TFile::Open(filename3);

    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie() || !f3 || f3->IsZombie()) {
        std::cerr << "Error: Couldn't open one of the files." << std::endl;
        return;
    }

    // Load histograms from each file
    TH1D *h1 = (TH1D*)f1->Get("BeamProjectionY;2");
    TH1D *h2 = (TH1D*)f2->Get("BeamProjectionY;2");
    TH1D *h3 = (TH1D*)f3->Get("BeamProjectionY;2");

    if (!h1 || !h2 || !h3) {
        std::cerr << "Error: Histogram BeamProjectionY;2 not found in one of the files." << std::endl;
        return;
    }

    // Normalize histograms so the peaks are at 1 (scale individually)
    h1->Scale(1.0 / h1->GetMaximum());  
    h2->Scale(1.0 / h2->GetMaximum());  
    h3->Scale(1.0 / h3->GetMaximum());  
    // Apply manual offsets to histograms
    h1->Add(h1, offset1); 
    h2->Add(h2, offset2); 
    h3->Add(h3, offset3); 


    h1->SetLineColor(kRed);  
    h2->SetLineColor(kBlue);
    h3->SetLineColor(kGreen); 

    // Set axis labels and titles
    h1->SetTitle("Beam Energy Profile;Beam Energy (MeV);Normalised Counts");
    h2->SetTitle("Beam Energy Profile;Beam Energy (MeV);Normalised Counts");
    h3->SetTitle("Beam Energy Profile;Beam Energy (MeV);Normalised Counts");

  
    TCanvas *c = new TCanvas("c", "Multiple Beam Profiles with Gaussian Fit and FWHM", 800, 600);

    // Set x-axis range to limit between 15 and 25 MeV
    h1->GetXaxis()->SetRangeUser(15, 25);
    h2->GetXaxis()->SetRangeUser(15, 25);
    h3->GetXaxis()->SetRangeUser(15, 25);

    double maxY = std::max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum()});
    h1->GetYaxis()->SetRangeUser(0, maxY * 1.1); // Slightly extend y-axis 

    // Draw the first histogram (beam profile)
    h1->SetLineWidth(2);
    h1->Draw("HIST");

    // Fit a Gaussian to the first histogram
    TF1 *gausFit1 = new TF1("gausFit1", "gaus", 15, 25);  
    h1->Fit(gausFit1, "R");

    // Draw the Gaussian fit over the histogram
    gausFit1->SetLineColor(kRed);  
    gausFit1->Draw("SAME");

    // Extract FWHM from the Gaussian fit
    double sigma1 = gausFit1->GetParameter(2);  // Standard deviation (sigma)
    double fwhm1 = 2.355 * sigma1;  // FWHM = 2.355 * sigma

    // Draw the second histogram with an offset
    h2->SetLineWidth(2);
    h2->Draw("HIST SAME");

    // Fit a Gaussian to the second histogram
    TF1 *gausFit2 = new TF1("gausFit2", "gaus", 15, 25);  // Fit Gaussian in range 15-25 MeV
    h2->Fit(gausFit2, "R");

    // Draw the Gaussian fit over the second histogram
    gausFit2->SetLineColor(kBlue);  // Color for the Gaussian fit
    gausFit2->Draw("SAME");

    // Extract FWHM from the Gaussian fit
    double sigma2 = gausFit2->GetParameter(2);  // Standard deviation (sigma)
    double fwhm2 = 2.355 * sigma2;  // FWHM = 2.355 * sigma

    // Draw the third histogram with an offset
    h3->SetLineWidth(2);
    h3->Draw("HIST SAME");

    // Fit a Gaussian to the third histogram
    TF1 *gausFit3 = new TF1("gausFit3", "gaus", 15, 25);  // Fit Gaussian in range 15-25 MeV
    h3->Fit(gausFit3, "R");

    // Draw the Gaussian fit over the third histogram
    gausFit3->SetLineColor(kGreen);  // Color for the Gaussian fit
    gausFit3->Draw("SAME");

    // Extract FWHM from the Gaussian fit
    double sigma3 = gausFit3->GetParameter(2);  // Standard deviation (sigma)
    double fwhm3 = 2.355 * sigma3;  // FWHM = 2.355 * sigma

    // Add a legend to the plot
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h1, label1, "l");
    legend->AddEntry(h2, label2, "l");
    legend->AddEntry(h3, label3, "l");
    legend->Draw();

    // Display FWHM values on the plot
    TLatex *fwhmLabel1 = new TLatex(15.5, maxY * 0.9, Form("FWHM1 = %.2f MeV", fwhm1));
    fwhmLabel1->SetTextColor(kRed);
    fwhmLabel1->Draw();

    TLatex *fwhmLabel2 = new TLatex(15.5, maxY * 0.8, Form("FWHM2 = %.2f MeV", fwhm2));
    fwhmLabel2->SetTextColor(kBlue);
    fwhmLabel2->Draw();

    TLatex *fwhmLabel3 = new TLatex(15.5, maxY * 0.7, Form("FWHM3 = %.2f MeV", fwhm3));
    fwhmLabel3->SetTextColor(kGreen);
    fwhmLabel3->Draw();

    // Update the canvas to show all histograms
    c->Update();
}

