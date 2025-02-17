#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include <iostream>
#include <array>

void fillHistogramsFromTree()
{
    // Open the input file containing the tree
    TFile *inFile = TFile::Open("HFTagging.root", "READ");
    if (!inFile || inFile->IsZombie())
    {
        std::cerr << "Error: Could not open input file" << std::endl;
        return;
    }

    // Get the tree
    TTree *tree = (TTree *)inFile->Get("jetTree");
    if (!tree)
    {
        std::cerr << "Error: Could not find tree" << std::endl;
        return;
    }

    // Create output file
    TFile *outFile = new TFile("AnalysisResultsFromTree_AFter.root", "RECREATE");

    // Create histograms (same binning as in bjet_analysis.C)
    // Basic jet kinematics
    TH1F *hJetPt = new TH1F("hJetPt", "Jet p_{T};p_{T} (GeV/c);Counts", 100, 0, 100);
    TH1F *hJetEta = new TH1F("hJetEta", "Jet #eta;#eta;Counts", 100, -5, 5);

    // Flavor-specific jet pT
    TH1F *hJetPt_b = new TH1F("hJetPt_b", "b-Jet p_{T};p_{T} (GeV/c);Counts", 100, 0, 100);
    TH1F *hJetPt_c = new TH1F("hJetPt_c", "c-Jet p_{T};p_{T} (GeV/c);Counts", 100, 0, 100);

    // SV properties
    TH1F *hSVMass = new TH1F("hSVMass", "SV Mass (inclusive);Mass (GeV/c^{2});Counts", 100, 0, 10);
    TH1F *hSVMass_b = new TH1F("hSVMass_b", "SV Mass (b-jets);Mass (GeV/c^{2});Counts", 100, 0, 10);
    TH1F *hSVMass_c = new TH1F("hSVMass_c", "SV Mass (c-jets);Mass (GeV/c^{2});Counts", 100, 0, 10);

    // Decay lengths
    TH1F *hDecayLength2D = new TH1F("hDecayLength2D", "2D Decay Length;L_{xy} (cm);Counts", 100, 0, 30);
    TH1F *hDecayLength2D_b = new TH1F("hDecayLength2D_b", "2D Decay Length (b-jets);L_{xy} (cm);Counts", 100, 0, 30);
    TH1F *hDecayLength2D_c = new TH1F("hDecayLength2D_c", "2D Decay Length (c-jets);L_{xy} (cm);Counts", 100, 0, 30);

    TH1F *hDecayLength3D = new TH1F("hDecayLength3D", "3D Decay Length;L_{xyz} (cm);Counts", 100, 0, 30);
    TH1F *hDecayLength3D_b = new TH1F("hDecayLength3D_b", "3D Decay Length (b-jets);L_{xyz} (cm);Counts", 100, 0, 30);
    TH1F *hDecayLength3D_c = new TH1F("hDecayLength3D_c", "3D Decay Length (c-jets);L_{xyz} (cm);Counts", 100, 0, 30);

    // Number of constituents and SVs
    TH1F *hNConstituents = new TH1F("hNConstituents", "Number of Jet Constituents (inclusive);N_{constituents};Counts", 50, 0, 50);
    TH1F *hNConstituents_b = new TH1F("hNConstituents_b", "Number of Jet Constituents (b-jets);N_{constituents};Counts", 50, 0, 50);
    TH1F *hNConstituents_c = new TH1F("hNConstituents_c", "Number of Jet Constituents (c-jets);N_{constituents};Counts", 50, 0, 50);

    TH1F *hNSVperJet = new TH1F("hNSVperJet", "Number of SV per Jet (inclusive);N_{SV};Counts", 10, 0, 10);
    TH1F *hNSVperJet_b = new TH1F("hNSVperJet_b", "Number of SV per Jet (b-jets);N_{SV};Counts", 10, 0, 10);
    TH1F *hNSVperJet_c = new TH1F("hNSVperJet_c", "Number of SV per Jet (c-jets);N_{SV};Counts", 10, 0, 10);

    // Set up branch variables
    float jetPt, jetEta, jetPhi, jetMass;
    int nTracks, nSV, jetFlavor;
    std::array<float, 13> svMass;
    std::array<float, 13> decayLength2D;
    std::array<float, 13> decayLength3D;

    // Set up branches
    tree->SetBranchAddress("jetPt", &jetPt);
    tree->SetBranchAddress("jetEta", &jetEta);
    tree->SetBranchAddress("jetPhi", &jetPhi);
    tree->SetBranchAddress("jetMass", &jetMass);
    tree->SetBranchAddress("nTracks", &nTracks);
    tree->SetBranchAddress("nSV", &nSV);
    tree->SetBranchAddress("jetFlavor", &jetFlavor);
    tree->SetBranchAddress("svMass", svMass.data());
    tree->SetBranchAddress("svDecayLen2D", decayLength2D.data());
    tree->SetBranchAddress("svDecayLen3D", decayLength3D.data());

    // Process all entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++)
    {
        tree->GetEntry(i);

        // Fill basic jet kinematics
        hJetPt->Fill(jetPt);
        hJetEta->Fill(jetEta);

        // Fill flavor-specific histograms
        if (jetFlavor == 5)
        { // b-jet
            hJetPt_b->Fill(jetPt);
            hNConstituents_b->Fill(nTracks);
            hNSVperJet_b->Fill(nSV);
        }
        else if (jetFlavor == 4)
        { // c-jet
            hJetPt_c->Fill(jetPt);
            hNConstituents_c->Fill(nTracks);
            hNSVperJet_c->Fill(nSV);
        }

        // Fill constituent multiplicity
        hNConstituents->Fill(nTracks);
        hNSVperJet->Fill(nSV);

        // Process secondary vertices
        for (int iSV = 0; iSV < nSV && iSV < 13; iSV++)
        {
            hSVMass->Fill(svMass[iSV]);
            hDecayLength2D->Fill(decayLength2D[iSV]);
            hDecayLength3D->Fill(decayLength3D[iSV]);

            if (jetFlavor == 5)
            {
                hSVMass_b->Fill(svMass[iSV]);
                hDecayLength2D_b->Fill(decayLength2D[iSV]);
                hDecayLength3D_b->Fill(decayLength3D[iSV]);
            }
            else if (jetFlavor == 4)
            {
                hSVMass_c->Fill(svMass[iSV]);
                hDecayLength2D_c->Fill(decayLength2D[iSV]);
                hDecayLength3D_c->Fill(decayLength3D[iSV]);
            }
        }
    }

    // Write and close files
    outFile->Write();
    outFile->Close();
    inFile->Close();

    std::cout << "Processed " << nEntries << " jets" << std::endl;
}