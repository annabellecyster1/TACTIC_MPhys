#include "processSim.hh"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>     // Added for std::unique_ptr
#include <numeric>
#include <vector>
#include <cmath>

void processSim(
  const TString& inFileName, 
  const TString& outFileName,
  const int mCut,
  const double Ethr_pad,
  const double protonVtxZCut, 
  const bool checkProtonPadHit)
{
  TFile* inFile = nullptr;
  std::unique_ptr<TTreeReader> reader;
  std::unique_ptr<TreeDataReader> dataReader;

  // Open the file and set up the TTreeReader.
  if (!OpenAndSetupTreeReader(inFileName, inFile, reader, dataReader)) {
      std::cerr << "Failed to open and set up the tree reader." << "\n";
      return;
  }

  // Output file and histograms
  std::unique_ptr<TFile> outFile(new TFile(outFileName, "RECREATE"));

  static constexpr int multiplicity_min = 0;  
  static constexpr int multiplicity_max = 10;

  static constexpr int Esummed_nBins = 100;  
  static constexpr double Esummed_min = 0.;  
  static constexpr double Esummed_max = 1.5;
  
  static constexpr int EpadMax_nBins = 100;  
  static constexpr double EpadMax_min = 0.;  
  static constexpr double EpadMax_max = 0.2;
  
  static constexpr int Theta_nBins = 100;
  static constexpr double Theta_min = 0.;
  static constexpr double Theta_max = 180.;

  static constexpr int vtxZ_nBins = 100;  
  static const double vtxZ_min = protonVtxZCut;  
  static constexpr double vtxZ_max = 84.;

  static constexpr int beamProfile_nBinsX = 60;  
  static constexpr double beamProfile_minX = 0.;  
  static constexpr double beamProfile_maxX = 60.;
  static constexpr int beamProfile_nBinsY = 200;  
  static constexpr double beamProfile_minY = -10.;  
  static constexpr double beamProfile_maxY = -beamProfile_minY;

  static constexpr int beamEProfile_nBinsX = 60;  
  static constexpr double beamEProfile_minX = 0.;  
  static constexpr double beamEProfile_maxX = 60.;
  static constexpr int beamEProfile_nBinsY = 2000;  
  static constexpr double beamEProfile_minY = 0.;  
  static constexpr double beamEProfile_maxY = 40.;

  std::unique_ptr<TH1I> multiplicityHist(new TH1I("multiplicityHist", "", 
                                                   kMultiplicity_nBins, 
                                                   multiplicity_min, 
                                                   multiplicity_max));

  std::unique_ptr<TH1D> EsummedHist(new TH1D("EsummedHist", "", 
                                             Esummed_nBins, 
                                             Esummed_min, 
                                             Esummed_max));

  std::unique_ptr<TH1D> vtxZHist(new TH1D("vtxZ", "", vtxZ_nBins, vtxZ_min, vtxZ_max));
  
  std::unique_ptr<TH1D> ThetaHist(new TH1D("ThetaHist", "", Theta_nBins,Theta_min, Theta_max));
  
  std::unique_ptr<TH1D> RadiusHist(new TH1D("RadiusHist", "", 100,0, 80));

  std::unique_ptr<TH2D> beamProfile(new TH2D("beamProfile", "", 
                                             beamProfile_nBinsX, beamProfile_minX, beamProfile_maxX,
                                             beamProfile_nBinsY, beamProfile_minY, beamProfile_maxY));

  std::unique_ptr<TH2D> beamEprofile(new TH2D("beamEprofile", "", 
                                              beamEProfile_nBinsX, beamEProfile_minX, beamEProfile_maxX, 
                                              beamEProfile_nBinsY, beamEProfile_minY, beamEProfile_maxY));
                                              
  std::unique_ptr<TH2D> beamEnergy = std::make_unique<TH2D>("beamEnergy", "", 100, 0, 80, 100, 0, 39);                                         
                                           
  std::unique_ptr<TH2D> ThetavsEmax(new TH2D("ThetavsEmax", "", 
                                              20,Theta_min, Theta_max, 
                                              20 , EpadMax_min , EpadMax_max));
   
  std::unique_ptr<TH2D> ThetavsEsummed(new TH2D("ThetavsEsummed", "", 
                                              Theta_nBins,Theta_min, Theta_max, 
                                              Esummed_nBins , 0 , 1.5)); 
                                              
  std::unique_ptr<TH2D> ThetavsVtxZ(new TH2D("ThetavsVtxZ", "", 
                                              20,Theta_min, Theta_max, 
                                              20, vtxZ_min, vtxZ_max));  
                                                                                         
  std::unique_ptr<TH2D> VtxZvsEtot(new TH2D("VtxZvsEtot", "", 
                                              20,vtxZ_min, vtxZ_max, 
                                              20, 0, 1.5)); 

  std::unique_ptr<TH2D> VtxZvsEmax(new TH2D("VtxZvsEmax", "", 
                                              100,vtxZ_min, vtxZ_max, 
                                              100, 0, 0.2)); 
                                              
  std::unique_ptr<TH2D> EtotEmaxVtx(new TH2D("EtotEmaxVtx", "", 
                                              100,0, 1, 
                                              100, 0, 0.2)); 
                                             
                                              
  std::unique_ptr<TH2D> VtxZreconvsEmax(new TH2D("VtxZreconvsEmax", "", 
                                              100,vtxZ_min, vtxZ_max, 
                                              100, 0, 0.2));                                             
                                              
  std::unique_ptr<TH1D> ReconstructedvtxZ(new TH1D("ReconstructedvtxZ", "", 50, -50, vtxZ_max));            
  
  std::unique_ptr<TH1D> Radius100keV(new TH1D("Radius100keV", "", 100,0, 80));                                                                             
                                            
  //std::unique_ptr<TH1D> VertexError(new TH1D("VertexError", "", 50, -10, 10));
                                             
  // For EpadMax vs Etot
  std::vector<double> Etot;
  std::vector<double> EpadMax;

  std::vector<std::vector<double>> Epads(kNRings, std::vector<double>(kNSectors, 0.));

  std::vector<std::vector<TH1D*>> padSpectra = InitializePadSpectra();

  bool padHit = false;
  double pVtxZ = kUndefinedVtxZ;
  double Theta = kUndefinedTheta;
  double trvertex = kUndefinedtrvertex;
  //double vertex_error = kUndefinedvertexErr;
  double RadiusSim = kUndefinedradius;
  int thisEventID = -1;
  int nP_det = 0;
  std::vector<int> multiplicities;

  TTree* tree = reader->GetTree();
  const Long64_t nHits = tree->GetEntries();
  std::cout << nHits << " hits found.\n";

  Long64_t progressCount = 0;

  // Loop over the TTree entries using TTreeReader.
  while (reader->Next()) 
  {
    const int eventID     = static_cast<int>(*(dataReader->eventID));
    const int trackID     = static_cast<int>(*(dataReader->trackID));
    const int atomicNb    = static_cast<int>(*(dataReader->atomicNb));
    const int ringID      = static_cast<int>(*(dataReader->ringID));
    const int sectorID    = static_cast<int>(*(dataReader->sectorID));
    const double time     = *(dataReader->time);
    const double E_kin    = *(dataReader->E_kin);
    const double E_ion    = *(dataReader->E_ionization);
    const double hitX     = *(dataReader->hitX);
    const double hitY     = *(dataReader->hitY);
    const double hitZ     = *(dataReader->hitZ);
    const double radius   = *(dataReader->radius);
    const double postR    = *(dataReader->postRadius);
    const double vtxZ     = *(dataReader->vtxZ);
    const double vtxE_kin = *(dataReader->vtxE_kin);
    const double vtxTheta = *(dataReader->vtxTheta);
    const double vtxPhi   = *(dataReader->vtxPhi);
    const double trackL   = *(dataReader->trackL);

    if (progressCount % kProgressUpdateInterval == 0 || 
        progressCount == nHits - 1) 
    {
      PrintProgressBar(progressCount, nHits);
    }

    progressCount++;

    // Beginning of a new event - record data from the previous event and reset 
    // the variables and vectors
    if (eventID != thisEventID)
    {
      
      const int multiplicity = GetMultiplicity(Epads, Ethr_pad);
      multiplicityHist->Fill(multiplicity);
      vtxZHist->Fill(pVtxZ); //pVtxZ

      // Check if the vertex Z satisfies the cut 
      // Then check the proton pad hit condition:
      // - If proton pad hit check is disabled, continue
      // - If proton pad hit check is enabled, ensure padHit is true
      if (pVtxZ > protonVtxZCut && (!checkProtonPadHit || padHit))
      {
        multiplicities.push_back(multiplicity);

        // Check the multiplicity cut
        if (multiplicity >= mCut) 
        {
          // Count the number of protons detected
          nP_det++;
          
          ThetaHist->Fill(Theta);
          
          RadiusHist->Fill(RadiusSim);
          
         
          if (12 <= radius && radius <= 50){
          trvertex = ReconstructVertex(radius, multiplicity, ringID);
          //vertex_error = ComputeVertexError(radius,multiplicity);
          ReconstructedvtxZ->Fill(trvertex);
          //VertexError->Fill(vertex_error);
          }
          
          // Vertex of the created proton
          
          
          
          ThetavsVtxZ->Fill(Theta,pVtxZ);

          // Record the energy deposited in the individual pads during the event
          FillPadSpectra(Epads, Ethr_pad, padSpectra);

          RecordSummedEnergy(Epads, Etot, EsummedHist.get());
          
          for (const auto& energy : Etot) {
              VtxZvsEtot->Fill(pVtxZ, energy); //pVtxZ original
		}


          // Get the energy max deposited within a pad
          const double EpadMax_evt = FindEpadMax(Epads, Ethr_pad);
          EpadMax.push_back(EpadMax_evt);
          
          
          for (int pVtxZ = 18; pVtxZ <= 20; ++pVtxZ) {
          	for (const auto& energy : Etot) {
             	EtotEmaxVtx->Fill(energy,EpadMax_evt);
             	}
          }
          
          if (EpadMax_evt > 0.12){
          Radius100keV->Fill(RadiusSim);
          }
          
          ThetavsEmax->Fill(Theta,EpadMax_evt);
          VtxZvsEmax->Fill(pVtxZ, EpadMax_evt);
          VtxZreconvsEmax->Fill(trvertex, EpadMax_evt);
          
          
        }
      }

      padHit = false;
      pVtxZ = kUndefinedVtxZ;
      trvertex = kUndefinedtrvertex;
      //vertex_error = kUndefinedvertexErr;
      Theta = kUndefinedTheta;
      RadiusSim = kUndefinedradius;
      thisEventID = eventID;

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
      beamEnergy->Fill(ringID, E_kin);
      beamEprofile->Fill(ringID, E_kin);
    }

    // Store proton vertex Z
    if (atomicNb == 1 && pVtxZ == kUndefinedVtxZ) 
    {
      pVtxZ = vtxZ;
    }
    
    if (RadiusSim == kUndefinedradius)
    {
    RadiusSim = radius;
    }
    
    if (Theta == kUndefinedTheta) 
    {
      Theta = vtxTheta;
    }

    const int ringIndex = ringID;
    const int sectorIndex = sectorID;

    if (IsValidIndex(ringID, sectorID)) 
    {
      Epads[ringIndex][sectorIndex] += E_ion;

      // Check that the proton hit a pad
      if (atomicNb == 1 && postR >= kDriftChamberR) 
      {
        padHit = true;
      }
    } 
    else if (ringIndex < 0 || ringIndex >= kNRings) 
    {
      std::cerr << "Out-of-bounds access: ringID=" << ringIndex 
                << ", sectorID=" << sectorIndex << "\n";
    }
  }
  // End loop over the events
  
  // To display 100%
  PrintProgressBar(nHits, nHits);

  //  TACTIC detection efficiency
  const double nEvts = tree->GetMaximum("EventID");
  const double efficiency = static_cast<double>(nP_det) * 100. / nEvts;
  
  std::cout << "\n\nDetection efficiency = " << std::setprecision(2) 
            << efficiency << "% \n";

  TGraphErrors* beamProfileVsZ = GetBeamFWHM(beamProfile.get());
  TGraphErrors* beamEprofileVsZ = GetBeamFWHM(beamEprofile.get());
  TGraphErrors* effVsMCut = GetEffVsMCut(multiplicities, nEvts);
  std::unique_ptr<TGraph> EtotvsEpadMax(new TGraph(Etot.size(), Etot.data(), EpadMax.data()));


  beamProfileVsZ->SetName("beamProfileVsZ");
  beamEprofileVsZ->SetName("beamEProfileVsZ");
  effVsMCut->SetName("effVsMCut");
  EtotvsEpadMax->SetName("EtotvsEpadMax");
  
  outFile->cd();
  multiplicityHist->Write();
  EsummedHist->Write();
  beamProfileVsZ->Write();
  beamEprofileVsZ->Write();
  effVsMCut->Write();
  beamEnergy->Write();
  ThetaHist->Write();
  vtxZHist->Write();
  EtotvsEpadMax->Write();
  ThetavsEmax->Write();
  ThetavsVtxZ->Write();
  VtxZvsEtot->Write();
  ReconstructedvtxZ->Write();
  RadiusHist->Write();
  VtxZvsEmax->Write();
  VtxZreconvsEmax->Write();
  EtotEmaxVtx->Write();
  Radius100keV->Write();
 // VertexError->Write();
  //ThetavsEsummed->Write();

  WriteEpadSpectra(padSpectra);
  
  outFile->Close();
  inFile->Close();
}

// This function opens the ROOT file, retrieves the TTree named "TACTIC;1",
// creates a TTreeReader, and initializes our branch readers.
static bool OpenAndSetupTreeReader(
  const TString &inFileName,
  TFile*& inFile,
  std::unique_ptr<TTreeReader>& reader,
  std::unique_ptr<TreeDataReader>& dataReader)
{
   // Open the ROOT file in READ mode.
   inFile = TFile::Open(inFileName, "READ");
   if (inFile == nullptr || !inFile->IsOpen() || inFile->IsZombie()) {
      std::cerr << "Error: Could not open file " << inFileName << "\n";
      return false;
   }

   // Retrieve the TTree named "TACTIC;1" from the file.
   auto* tree = dynamic_cast<TTree*>(inFile->Get("TACTIC;1"));
   if (!tree) {
      std::cerr << "Error: Could not find TTree 'TACTIC;1' in file " 
                << inFileName << "\n";
      inFile->Close();
      return false;
   }

   // Create a TTreeReader from the TTree.
   reader = std::unique_ptr<TTreeReader>(new TTreeReader(tree));

   // Initialize the branch readers.
   dataReader = std::unique_ptr<TreeDataReader>(new TreeDataReader(*reader));

   return true;
}

// Function to initialize pad spectra (histograms) for a given number of rings and sectors.
std::vector<std::vector<TH1D*>> InitializePadSpectra() 
{
  // Create a 2D vector of TH1D* with dimensions kNRings x kNSectors.
  std::vector<std::vector<TH1D*>>
   padSpectra(kNRings, std::vector<TH1D*>(kNSectors, nullptr));

  // Loop over each ring and sector to create a histogram.
  for (int i = 0; i < kNRings; ++i) 
  {
    for (int j = 0; j < kNSectors; ++j)
    {
      // Create a unique name for each histogram using the ring and sector indices.
      const TString padHistName = Form("EpadHist_%i_%i", i, j);
            
      // Create a new TH1D histogram.
      padSpectra[i][j] =
       new TH1D(padHistName, "", kPadSpectra_nBins, kPadSpectra_min, kPadSpectra_max);
    }
  }
    
  return padSpectra;
}

void FillPadSpectra(const std::vector<std::vector<double>> &Epads,
                    double Ethr_pad,
                    const std::vector<std::vector<TH1D*>> &padSpectra)
{
    for (int i = 0; i < kNRings; ++i) {
        for (int j = 0; j < kNSectors; ++j) {
            // Apply the energy threshold before filling the histogram.
            if (Epads[i][j] > Ethr_pad) {
                padSpectra[i][j]->Fill(Epads[i][j]);
            }
        }
    }
}

double ReconstructVertex(double radius, double multiplicity, int ringID)
{
    double trbase, trangle, trangle_rad, trbasetotal, trvertex, posdiff, fpadpos;
    double ds,trangle_error;
    
    
    trbase = 4.2 * multiplicity;
    
    trangle = atan(radius / trbase) * 57.325;
    
    trangle_rad = trangle * 3.14159/180;
    
    double radius2 = radius + 12; // double radius2 = radius + 12;
    
    trbasetotal = radius2 / tan(trangle_rad);
    
    //int maxpad = *std::max_element(ringID.begin(), ringID.end());

    fpadpos = ringID * 4.2;
    
    posdiff = fpadpos - trbasetotal;
    trvertex = posdiff - 126/2;  

    return trvertex;  // Return the reconstructed vertex
  
}

//double ComputeVertexError(double radius, double multiplicity) {

    //double trbase = 4.2 * multiplicity;
    //double trangle_rad = atan(radius / trbase);
// Error calculations//
    //const double perr = 2.1 * 2.1;  // Beam spot error (mm^2)
   // const double rcerr = (2.4 * 2.4) + (2 * 2);  // Beam staggering error (mm^2)
    
   // double rderr = 0.02 * radius2;
    
   // ds=(rderr/trbase)*(cos(trangle_rad)*cos(trangle_rad));
    
    //trangle_error=ds*57.32;
    
    //double rterr = (rcerr +(rderr*rderr))/(tan(trangle_rad)*tan(trangle_rad));
   // double rzerr=(ds*ds*radius2*radius2)/(cos(trangle_rad)*cos(trangle_rad));
   // vertex_error=sqrt(perr+rterr+rzerr);
    
    //return vertex_error;
   // }

void RecordSummedEnergy(const std::vector<std::vector<double>>& Epads,
                        std::vector<double>& Etot, 
                        TH1D* EsummedHist)
{
    double Esummed = 0.0;

    // Sum the energy deposited in each pad
    for (const auto& ring : Epads) {
        Esummed += std::accumulate(ring.begin(), ring.end(), 0.0);
    }

    // Record the total energy in the vector
    Etot.push_back(Esummed);

    // If some energy was deposited, fill the histogram (scaling by 1000)
    if (Esummed > 0.0) {
        EsummedHist->Fill(Esummed); //take away scaling by 1000
    }
   
}

void PrintProgressBar(const Long64_t progress, const Long64_t nHits) 
{
  const float percentage = static_cast<float>(progress) / nHits;
  const int pos = kBarWidth * percentage;

  std::cout << "[";

  for (int i = 0; i < kBarWidth; ++i) 
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
  
  std::cout << "] " << static_cast<int>(percentage * 100.) << "%\r";
  std::cout.flush();
}

bool IsValidIndex(const int ringIndex, const int sectorIndex) 
{
  return ringIndex >= 0 && ringIndex < kNRings &&
         sectorIndex >= 0 && sectorIndex < kNSectors;
}

static int GetMultiplicity(
  const std::vector<std::vector<double>>& Epads, 
  const double& Ethr_pad)
{
  int multiplicity = 0;

  // Iterate through each ring in Epads
  for (const auto& ring : Epads) 
  {
    // Iterate through each sector in the current ring
    for (const auto& energy : ring) 
    {
      // Apply the energy threshold before counting
      if (energy > Ethr_pad)
        ++multiplicity;
    }
  }

  return multiplicity;
}

static double FindEpadMax(
  const std::vector<std::vector<double>>& Epads, 
  const double& Ethr_pad) 
{
  double EpadMax = std::numeric_limits<double>::lowest();

  for (const auto& row : Epads) 
  {
    for (double EPad : row) 
    {
      if (EPad > EpadMax)
        EpadMax = EPad;
    }
  }

  return EpadMax;
}

TGraphErrors* GetBeamFWHM(const TH2D* beamProfile) 
{
  std::vector<double> beamFWHMs;
  std::vector<double> beamFWHM_errors;
  std::vector<double> Zpositions;

  // Gaussian fitting
  TF1* gausFit = new TF1("gausFit", "gaus");

  // Define constants outside the loop
  constexpr double padLength = 251.9 / 60.;
  constexpr double TACTIC_halfLength = 251.9 / 2.;
  constexpr double sigmaToFWHM = 2.35482;

  // Loop over all x-bins
  for (int i = 1; i <= beamProfile->GetNbinsX(); ++i) 
  {
    // Project the y-axis slice at the given bin
    TH1D* projY = beamProfile->ProjectionY(Form("_py_%d", i), i, i);

    if (projY == nullptr) 
    {
      std::cerr << "Could not project in the y-axis for bin " << i << std::endl;
      continue;
    }

    // Fit the projection if it has entries
    if (projY->GetEntries() > 0) 
    {
      projY->Fit(gausFit, "Q", "", projY->GetXaxis()->GetXmin(),
                 projY->GetXaxis()->GetXmax());

      // Check fit results
      if (gausFit->GetParameter(2) > 0.) 
      {
        // Z-position for the center of bin i
        const double Zposition = ((i + 0.5) * padLength) - TACTIC_halfLength;
        const double sigma = gausFit->GetParameter(2);
        const double fwhm = sigmaToFWHM * sigma;
        const double sigma_error = gausFit->GetParError(2);
        const double fwhm_error = sigmaToFWHM * sigma_error;

        if (fwhm > 0.1) //0.1
        {
          Zpositions.push_back(Zposition); 
          beamFWHMs.push_back(fwhm);
          beamFWHM_errors.push_back(fwhm_error);    
          
          const int desiredBinIndex = 31; // Set this to whatever bin you need
          if (i==desiredBinIndex){
	  const int binLower = desiredBinIndex ; // Lower bound for projection (1 bin before)
   	  const int binUpper = desiredBinIndex; // Upper bound for projection (1 bin after)

         // Project the Y-axis slice for the desired bin range (binLower, binUpper)
          TH1D* projY = beamProfile->ProjectionY("BeamProjectionY", binLower, binUpper);
            // Clone to avoid memory management issues (if necessary)
          //TH1D* projY_clone = (TH1D*)projY->Clone("BeamProjectionY");
          projY->Fit(gausFit, "Q", "");
          //projY->GetXaxis()->SetRangeUser(20,25);

            // Write the projection to the ROOT file
          projY->Write();
          gausFit->Write();
          
          }  
          
          
        }
      }
    } 

   // delete projY;
  }

  //delete gausFit;

  TGraphErrors* beamFWHMVsZ = new TGraphErrors(Zpositions.size(), 
                                               Zpositions.data(), 
                                               beamFWHMs.data(), 
                                               nullptr, 
                                               beamFWHM_errors.data());
  return beamFWHMVsZ;
}


TGraphErrors* GetEffVsMCut(const std::vector<int>& multiplicities, 
                           const int nEvts)
{
  // Vector to store counts for multiplicities between 0 and 10
  std::vector<int> counts(kMultiplicity_nBins, 0);

  for (int multiplicity : multiplicities) 
  {
    if (multiplicity < kMultiplicity_nBins) 
    {
      for (int i = multiplicity; i >= 0; i--) 
        counts.at(i)++;
    }
  }

  std::vector<double> mCuts;
  std::vector<double> efficiencies;
  std::vector<double> efficiency_errors;

  for (size_t i = 0; i < counts.size(); ++i) 
  {
    mCuts.push_back(i);    
    const double efficiency = static_cast<double>(counts[i]) * 100. / nEvts;
    efficiencies.push_back(efficiency);
    const double efficiency_error =
      100. / nEvts * std::sqrt(counts[i] * (nEvts - counts[i]) / static_cast<double>(nEvts));
    efficiency_errors.push_back(efficiency_error);
  }

  TGraphErrors* effVsM = new TGraphErrors(mCuts.size(), 
                                          mCuts.data(), 
                                          efficiencies.data(), 
                                          nullptr,
                                          efficiency_errors.data());
  return effVsM;
}

static std::vector<std::vector<TH1D*>>
 GetEpadSpectra(const std::vector<std::vector<TH1D*>>& padSpectra) 
{
  std::vector<std::vector<TH1D*>>
   EpadSpectra(kNRings, std::vector<TH1D*>(kNSectors, nullptr));

  for (int i = 0; i < kNRings; ++i) 
  {
    for (int j = 0; j < kNSectors; ++j) 
    {
      TH1D* EpadHist = padSpectra[i][j];

      // Create a new histogram with the same binning
      const TString histName = TString::Format("hist_%d_%d", i, j);
      
      TH1D* newEpadHist = new TH1D(histName, "", 
                                   kPadSpectra_nBins, 
                                   kPadSpectra_min, 
                                   kPadSpectra_max);

      // Copy contents from the original histogram
      for (int bin = 1; bin <= kPadSpectra_nBins; ++bin) 
      {
        newEpadHist->SetBinContent(bin, EpadHist->GetBinContent(bin));
      }
      EpadSpectra[i][j] = newEpadHist;
    }
  }

  return EpadSpectra;
}

// Updated to match the declaration (no extra parameters)
static void WriteEpadSpectra(const std::vector<std::vector<TH1D*>>& padSpectra)
{
    // Obtain the pad spectra from the existing padSpectra.
    std::vector<std::vector<TH1D*>> EpadSpectra = GetEpadSpectra(padSpectra);

    // Loop over rings and sectors.
    for (int i = 0; i < kNRings; ++i) {
        for (int j = 0; j < kNSectors; ++j) {
            TH1D* EPadSpectrum = EpadSpectra[i][j];
            if (EPadSpectrum) {
                // Create a name for the histogram in the format "EPadSpectrum_i_j".
                const TString histName = TString::Format("EPadSpectrum_%d_%d", i, j);
                EPadSpectrum->SetName(histName);
                EPadSpectrum->Write();
            }
        }
    }
}
