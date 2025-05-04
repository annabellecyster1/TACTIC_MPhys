TGraph* GetGraph(const TString);

void plotBeamEprofile(const TString f1, const TString f2, 
                      const TString label1, const TString label2) 
{
  TGraph* graph1 = GetGraph(f1);
  TGraph* graph2 = GetGraph(f2);

  graph1->SetMarkerColor(kBlue);
  graph1->SetMarkerStyle(21);
  graph1->SetMarkerSize(1.8);
  graph1->SetLineWidth(3);
  graph1->SetLineColor(kBlue);

  graph2->SetMarkerColor(kRed);
  graph2->SetMarkerStyle(21);
  graph2->SetMarkerSize(1.8);
  graph2->SetLineWidth(3);
  graph2->SetLineColor(kRed);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 700);
  if (gPad)
  {
    gPad->SetLeftMargin(.145);
    gPad->SetBottomMargin(.15);
    gPad->SetTopMargin(.1);
  }

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(graph1);
  mg->Add(graph2);
  mg->GetXaxis()->SetRangeUser(-40., 40.);
  mg->GetYaxis()->SetRangeUser(0., 1.2);
  mg->GetXaxis()->SetTitle("Z (mm)");
  mg->GetYaxis()->SetTitle("Beam energy profile FWHM (MeV)");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.3);
  mg->GetYaxis()->SetTitleOffset(1.5);
  mg->GetXaxis()->SetLabelSize(.04);
  mg->GetYaxis()->SetLabelSize(.04);
  mg->GetXaxis()->SetTitleSize(.045);
  mg->GetYaxis()->SetTitleSize(.04);

  mg->Draw("al");

  TLegend *legend = new TLegend(.7, .7, .87, .85);
  legend->SetTextSize(.04);
  legend->AddEntry(graph1, "^{14}N", "lf");
  legend->AddEntry(graph2, "^{23}Na", "lf");
  legend->Draw();

  c1->Print("beamEProfileVsZ.png");
}

TGraph* GetGraph(const TString inputFileName)
{
  TFile *file = TFile::Open(inputFileName);
  file->cd();
  TGraph* graph = (TGraph*)file->Get("beamEProfileVsZ");
  file->Close();

  return graph;
}
