#include "processSim.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>
#include <cmath>

static constexpr int nRings = 61;
static constexpr int nSectors = 8;

void PrintProgressBar(const int progress, const int total, const int barWidth = 70);
bool IsValidIndex(const int ringIndex, const int sectorIndex);
int GetMultiplicity(const std::vector<std::vector<double>>& Epads, const double& Ethr_pad);
double FindEpadMax(const std::vector<std::vector<double>>& Epads, const double& Ethr_pad);
TGraphErrors* GetEffVsMCut(const std::vector<int>& multiplicities, const int nEvts);
void PlotEpadHists(const std::vector<std::vector<TH1D*>>& EpadHists);
std::vector<std::vector<TH1D*>> GetEpadSpectra(const std::vector<std::vector<TH1D*>>& EpadHists);

void processSim(
  const TString& inFileName, 
  const TString& outFileName,
  const int mCut,
  const double Ethr_pad,
  const double protonVtxZCut, 
  const bool checkProtonPadHit)
{
  // Open the ROOT file
  TFile* inFile = new TFile(inFileName, "READ");

  if (!inFile->IsOpen() || inFile->IsZombie()) 
  {
    std::cerr << "Error: Could not open file " << inFileName << std::endl;
    return;
  }

  // Reading the data (in TTree format)
  TTree* tree = dynamic_cast<TTree*>(inFile->Get("TACTIC;1"));

  if (!tree) 
  {
    std::cerr << "Error: Could not find TTree in the file " << inFileName << std::endl;
    inFile->Close();
    return;
  }

  // Variables to read from the data
  double eventID, trackID, atomicNb, time, E_kin, E_ionization, sectorID, 
         ringID, hitX, hitY, hitZ, radius, postRadius, vtxZ, vtxE_kin, 
         vtxTheta, vtxPhi, trackL;

  // Connect variables to TTree data
  tree->SetBranchAddress("EventID", &eventID);
  tree->SetBranchAddress("TrackID", &trackID);
  tree->SetBranchAddress("ParticleID", &atomicNb);
  tree->SetBranchAddress("Time", &time);
  tree->SetBranchAddress("E_kin", &E_kin);
  tree->SetBranchAddress("E_ionization", &E_ionization);
  tree->SetBranchAddress("Sector", &sectorID);
  tree->SetBranchAddress("Pad", &ringID);
  tree->SetBranchAddress("X", &hitX);
  tree->SetBranchAddress("Y", &hitY);
  tree->SetBranchAddress("Z", &hitZ);
  tree->SetBranchAddress("R", &radius);
  tree->SetBranchAddress("PostR", &postRadius);
  tree->SetBranchAddress("VertexZ", &vtxZ);
  tree->SetBranchAddress("VertexE_kin", &vtxE_kin);
  tree->SetBranchAddress("VertexTheta", &vtxTheta);
  tree->SetBranchAddress("VertexPhi", &vtxPhi);
  tree->SetBranchAddress("TrackLength", &trackL);

  // Output file and histograms
  TFile* outFile = new TFile(outFileName, "RECREATE");
  
  TH1D* vtxZHist = new TH1D("vtxZ", "", 27, protonVtxZCut, 84.);
  TH1D* vtxZ_det = new TH1D("vtxZ_det", "", 27, protonVtxZCut, 84.);
  TH2D* beamProfile = new TH2D("beamProfile", "", 60, 0., 60., 200, -10., 10.);
  TH2D* beamEprofile = new TH2D("beamEprofile", "", 60, 0., 60., 2000, 0., 40.);
  TH1D* ZHist = new TH1D("ZHist", "", 500, -100., 45.);

  // For EpadMax vs Etot
  std::vector<double> Etot;
  std::vector<double> EpadMax;
  std::vector<double> VertexZ;

  std::vector<std::vector<double>> Epads(nRings, std::vector<double>(nSectors, 0.));
  std::vector<std::vector<TH1D*>> EpadHists(nRings, std::vector<TH1D*>(nSectors, nullptr));

  for (int i = 0; i < nRings; ++i) 
  {
    for (int j = 0; j < nSectors; ++j) 
    {
      const TString padHistName = Form("EpadHist_%i_%i", i, j);
      EpadHists[i][j] = new TH1D(padHistName, "", 100, 0., .2);
    }
  }

  bool padHit = false;
  double pVtxZ = -9999.;
  int thisEventID = -1;
  int nP_det = 0;
  std::vector<int> multiplicities;

  const int nHits = tree->GetEntries();
  std::cout << nHits << " hits found." << std::endl;

  for (int i = 0; i < nHits; ++i)
  {
    tree->GetEntry(i);

    if (i % 10000 == 0 || i == nHits - 1) 
    {
      PrintProgressBar(i, nHits);
    }

    // Beginning of a new event
    if (eventID != thisEventID)
    {
      // Vertex of the created proton
      vtxZHist->Fill(pVtxZ);
      const int multiplicity = GetMultiplicity(Epads, Ethr_pad);

      if (pVtxZ > protonVtxZCut && (!checkProtonPadHit || padHit))
      {
        multiplicities.push_back(multiplicity);

        // Check the multiplicity cut
        if (multiplicity >= mCut) 
        {
          nP_det++;
          vtxZ_det->Fill(pVtxZ);
          VertexZ.push_back(pVtxZ);

          // Record the energy deposited in the individual pads
          for (int ii = 0; ii < nRings; ++ii) 
          {
            for (int jj = 0; jj < nSectors; ++jj) 
            {
              if (Epads[ii][jj] > Ethr_pad) 
              {
                EpadHists[ii][jj]->Fill(Epads[ii][jj]);
              }
            }
          }

          double Etot_evt = 0.;
          for (const auto& ring : Epads) 
          {
            Etot_evt += std::accumulate(ring.begin(), ring.end(), 0.);
          }
          Etot.push_back(Etot_evt);

          const double EpadMax_evt = FindEpadMax(Epads, Ethr_pad);
          EpadMax.push_back(EpadMax_evt);
        }
      }

      padHit = false;
      pVtxZ = -9999.;
      thisEventID = (int)eventID;

      // Reset the Epad to 0
      for (auto& ring : Epads)
      {
        std::fill(ring.begin(), ring.end(), 0.0);
      }
    }

    // Store beam energy and position
    if (trackID == 1)
    {
      beamProfile->Fill(ringID, hitX);
      beamEprofile->Fill(ringID, E_kin);
      ZHist->Fill(ringID);
    }

    // Store proton vertex Z
    if (atomicNb == 1 && pVtxZ < -9000.) 
    {
      pVtxZ = vtxZ;
    }

    const int ringIndex = static_cast<int>(ringID);
    const int sectorIndex = static_cast<int>(sectorID);

    if (IsValidIndex(ringIndex, sectorIndex)) 
    {
      Epads[ringIndex][sectorIndex] += E_ionization;

      // Check that the proton hit a pad
      if (atomicNb == 1 && postRadius >= 50.) 
      {
        padHit = true;
      }
    } 
    else if (ringIndex < 0 && ringIndex >= nRings) 
    {
      std::cerr << "Out-of-bounds access: ringID=" << ringIndex 
                << ", sectorID=" << sectorIndex << std::endl;
    }
  }
  // End loop over the events

  // To display 100%
  PrintProgressBar(nHits, nHits);

  // Calculate TACTIC detection efficiency
  const double nEvts = tree->GetMaximum("EventID");
  std::cout << "\n\nDetection efficiency = " << std::setprecision(2)
    << (double)nP_det * 100. / nEvts << "% \n";

  // Instead of using GetBeamFWHM, we now directly compute mean beam energy vs Z
  std::vector<double> Zpositions;
  std::vector<double> meanEnergies;
  std::vector<double>stdDevs;

  // Calculate mean energy for each ring bin
  for (int i = 1; i <= beamEprofile->GetNbinsX(); ++i) 
  {
    TH1D* projY = beamEprofile->ProjectionY(Form("projY_%d", i), i, i);
    if (projY->GetEntries() > 0) 
    {
      double meanE = projY->GetMean();
      double stdDev= projY->GetStdDev();
      double Zposition = (i + 0.5) * (251.9 / 60.) - 251.9 / 2.;
      Zpositions.push_back(Zposition);
      meanEnergies.push_back(meanE);
      stdDevs.push_back(stdDev);
    }
    delete projY;
  }

  TGraphErrors *beamEnergyVsZ = new TGraphErrors((int)Zpositions.size(), Zpositions.data(), meanEnergies.data(),nullptr, stdDevs.data());
  beamEnergyVsZ->SetName("beamEnergyVsZ");
  beamEnergyVsZ->SetTitle("Mean Beam Energy vs Z;Z Position (mm);Mean Beam Energy (MeV)");

  TGraphErrors* effVsMCut = GetEffVsMCut(multiplicities, nEvts);
  TGraph *EtotvsEpadMax = new TGraph((int)Etot.size(), Etot.data(), EpadMax.data());
  TGraph *VertexZvsEpadMax = new TGraph((int)VertexZ.size(), VertexZ.data(), EpadMax.data());

  effVsMCut->SetName("effVsMCut");
  EtotvsEpadMax->SetName("EtotvsEpadMax");
  VertexZvsEpadMax->SetName("VertexZvsEpadMax");
  
  outFile->cd();
  beamEnergyVsZ->Write();
  effVsMCut->Write();
  EtotvsEpadMax->Write();
  VertexZvsEpadMax->Write();

  PlotEpadHists(EpadHists);
  
  std::vector<std::vector<TH1D*>> EpadSpectra = GetEpadSpectra(EpadHists);

  for (int i = 0; i < nRings; ++i) 
  {
    for (int j = 0; j < nSectors; ++j) 
    {
      TH1D* EPadSpectrum = EpadSpectra[i][j];
      if (EPadSpectrum) 
      {
        TString histName = TString::Format("EPadSpectrum_%d_%d", i, j);
        EPadSpectrum->SetName(histName);
        EPadSpectrum->Write();
      }
    }
  }
  
  delete beamEnergyVsZ;
  delete effVsMCut;
  delete EtotvsEpadMax;
  delete VertexZvsEpadMax;

  for (int i = 0; i < nRings; ++i) 
  {
    for (int j = 0; j < nSectors; ++j) 
    {
      delete EpadSpectra[i][j];
    }
  }

  outFile->Close();
  inFile->Close();
}

void PrintProgressBar(const int progress, const int total, const int barWidth) 
{
  const float percentage = (float)progress / total;
  const int pos = (int)(barWidth * percentage);

  std::cout << "[";

  for (int i = 0; i < barWidth; ++i) 
  {
    if (i < pos)
    {
      std::cout << "=";
    }
    else if (i == pos)
    {
      std::cout << ">";
    }
    else 
    {
      std::cout << " ";
    }
  }
  
  std::cout << "] " << (int)(percentage * 100.) << "%\r";
  std::cout.flush();
}

bool IsValidIndex(const int ringIndex, const int sectorIndex) 
{
  return ringIndex >= 0 && 
         ringIndex < nRings && 
         sectorIndex >= 0 && 
         sectorIndex < nSectors;
}

int GetMultiplicity(
  const std::vector<std::vector<double>>& Epads, 
  const double& Ethr_pad)
{
  int multiplicity = 0;
  for (const auto& ring : Epads) 
  {
    for (const auto& energy : ring) 
    {
      if (energy > Ethr_pad) ++multiplicity;
    }
  }
  return multiplicity;
}

double FindEpadMax(
  const std::vector<std::vector<double>>& Epads, 
  const double& /*Ethr_pad*/) 
{
  double EpadMax = std::numeric_limits<double>::lowest();
  for (const auto& row : Epads) 
  {
    for (double EPad : row) 
    {
      if (EPad > EpadMax) 
      {
        EpadMax = EPad;
      }
    }
  }
  return EpadMax;
}


TGraphErrors* GetEffVsMCut(
  const std::vector<int>& multiplicities, 
  const int nEvts)
{
  std::vector<int> counts(11, 0);
  for (int multiplicity : multiplicities) 
  {
    if (multiplicity <= 10) 
    {
      for (int i = multiplicity; i >= 0; i--) counts[i]++;
    }
  }

  std::vector<double> mCuts;
  std::vector<double> efficiencies;
  std::vector<double> efficiency_errors;

  for (size_t i = 0; i < counts.size(); ++i) 
  {
    mCuts.push_back((double)i);    
    const double efficiency = (double) counts[i] * 100. / nEvts;
    efficiencies.push_back(efficiency);
    const double efficiency_error =
     100. / nEvts * std::sqrt((double)counts[i] * (nEvts - counts[i]) / nEvts);
    efficiency_errors.push_back(efficiency_error);
  }

  TGraphErrors* effVsM = new TGraphErrors((int)mCuts.size(), 
                                          mCuts.data(), 
                                          efficiencies.data(), 
                                          nullptr,
                                          efficiency_errors.data());
  return effVsM;
}

void PlotEpadHists(const std::vector<std::vector<TH1D*>>& EpadHists)
{  
  const int totalPads = nRings * (nSectors + 1); 
  TH2D* hist2D = new TH2D("hist2D", "Energy Distribution; Energy; Pad Number", 
                          20, 0, .1, totalPads, 0, totalPads);

  for (int ring = 0; ring < nRings; ++ring) 
  {
    for (int sector = 0; sector < nSectors; ++sector) 
    {
      TH1D* hist = EpadHists[ring][sector];
      if (hist) 
      {
        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) 
        {
          const double energy = hist->GetBinCenter(bin);
          const double count = hist->GetBinContent(bin);
          int padNumber = ring * (nSectors + 1) + sector;
          hist2D->Fill(energy, padNumber, count);
        }
      }
    }
  }
}

std::vector<std::vector<TH1D*>> GetEpadSpectra(
  const std::vector<std::vector<TH1D*>>& EpadHists) 
{
  std::vector<std::vector<TH1D*>> EpadSpectra(nRings, std::vector<TH1D*>(nSectors, nullptr));
  for (int i = 0; i < nRings; ++i) 
  {
    for (int j = 0; j < nSectors; ++j) 
    {
      TH1D* EpadHist = EpadHists[i][j];
      int nBins = EpadHist->GetNbinsX();
      TString histName = TString::Format("hist_%d_%d", i, j);
      TH1D* newEpadHist = new TH1D(histName, "", 
                                   EpadHist->GetNbinsX(), 
                                   EpadHist->GetXaxis()->GetXmin(), 
                                   EpadHist->GetXaxis()->GetXmax());
      for (int bin = 1; bin <= nBins; ++bin) 
      {
        newEpadHist->SetBinContent(bin, EpadHist->GetBinContent(bin));
      }
      EpadSpectra[i][j] = newEpadHist;
    }
  }
  return EpadSpectra;
}

