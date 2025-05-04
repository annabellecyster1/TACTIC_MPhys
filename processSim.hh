#ifndef PROCESSSIM_H
#define PROCESSSIM_H

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <vector>

// Structure to hold TTreeReaderValues for each branch.
struct TreeDataReader {
   TTreeReaderValue<double> eventID;
   TTreeReaderValue<double> trackID;
   TTreeReaderValue<double> atomicNb;
   TTreeReaderValue<double> time;
   TTreeReaderValue<double> E_kin;
   TTreeReaderValue<double> E_ionization;
   TTreeReaderValue<double> sectorID;
   TTreeReaderValue<double> ringID;
   TTreeReaderValue<double> hitX;
   TTreeReaderValue<double> hitY;
   TTreeReaderValue<double> hitZ;
   TTreeReaderValue<double> radius;
   TTreeReaderValue<double> postRadius;
   TTreeReaderValue<double> vtxZ;
   TTreeReaderValue<double> vtxE_kin;
   TTreeReaderValue<double> vtxTheta;
   TTreeReaderValue<double> vtxPhi;
   TTreeReaderValue<double> trackL;

   // Constructor: each TTreeReaderValue is initialized with the provided reader and branch name.
   explicit TreeDataReader(TTreeReader &reader)
      : eventID(reader, "EventID"),
        trackID(reader, "TrackID"),
        atomicNb(reader, "ParticleID"),
        time(reader, "Time"),
        E_kin(reader, "E_kin"),
        E_ionization(reader, "E_ionization"),
        sectorID(reader, "Sector"),
        ringID(reader, "Pad"),
        hitX(reader, "X"),
        hitY(reader, "Y"),
        hitZ(reader, "Z"),
        radius(reader, "R"),
        postRadius(reader, "PostR"),
        vtxZ(reader, "VertexZ"),
        vtxE_kin(reader, "VertexE_kin"),
        vtxTheta(reader, "VertexTheta"),
        vtxPhi(reader, "VertexPhi"),
        trackL(reader, "TrackLength")
   {}
};

// Pass smart pointers by reference so ownership can be modified.
static bool OpenAndSetupTreeReader(
  const TString &inFileName,
  TFile*& inFile,
  std::unique_ptr<TTreeReader>& reader,
  std::unique_ptr<TreeDataReader>& dataReader);

std::vector<std::vector<TH1D*>> InitializePadSpectra();

void FillPadSpectra(const std::vector<std::vector<double>> &Epads,
                    double Ethr_pad,
                    const std::vector<std::vector<TH1D*>> &padSpectra);

void RecordSummedEnergy(const std::vector<std::vector<double>>& Epads,
                        std::vector<double>& Etot,
                        TH1D* EsummedHist);

// Function to show the progress bar
void PrintProgressBar(
  const Long64_t progress, 
  const Long64_t nHits);

bool IsValidIndex(
  const int ringIndex, 
  const int sectorIndex);

// Function to calculate multiplicity based on energy deposits.
static int GetMultiplicity(
  const std::vector<std::vector<double>>& Epads, 
  const double& Ethr_pad);

// Function to find the maximum energy deposited in any pad.
static double FindEpadMax(
  const std::vector<std::vector<double>>& Epads, 
  const double& Ethr_pad);

double ReconstructVertex(double radius, double multiplicity, int ringID);

// Function to compute the beam profile Full Width at Half Maximum (FWHM).
TGraphErrors* GetBeamFWHM(const TH2D* beamProfile);

// Function to calculate efficiency versus multiplicity cut.
TGraphErrors* GetEffVsMCut(
  const std::vector<int>& multiplicities, 
  const int nEvts);

// Function to get energy deposit spectra from histograms.
static std::vector<std::vector<TH1D*>>
 GetEpadSpectra(const std::vector<std::vector<TH1D*>>& padSpectra);

static void WriteEpadSpectra(const std::vector<std::vector<TH1D*>>& padSpectra);

// Main function to process simulation data.
void processSim(
  const TString& inFileName = "TACTIC.root",
  const TString& outFileName = "processedSimData.root",
  const int mCut = 5,
  const double Ethr_pad = 0.02,
  const double protonVtxZCut = -38.,
  const bool checkProtonPadHit = true);

static constexpr int kNRings = 61;
static constexpr int kNSectors = 8;
static constexpr double kUndefinedVtxZ = -9999.;
static constexpr double kUndefinedTheta = -9999.;
static constexpr double kUndefinedvertexErr = -9999.;
static constexpr double kUndefinedtrvertex = -9999;
static constexpr double kUndefinedradius = -9999;
static constexpr int kProgressUpdateInterval = 10000;
static constexpr double kDriftChamberR = -9999.;
static constexpr int kBarWidth = 50;

static constexpr int kMultiplicity_nBins = 11;
static constexpr int kPadSpectra_nBins = 100;  
static constexpr double kPadSpectra_min = 0.;  
static constexpr double kPadSpectra_max = 10.;

#endif // PROCESSSIM_H
