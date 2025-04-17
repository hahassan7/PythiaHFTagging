#include "TaggingUtilities.h"
#include "JTreeHFFeatures.h"
#include "JetTreeGenerator.h"
#include "THnSparse.h"
#include "TSystem.h"

class BjetAnalysis
{
private:
    std::unique_ptr<JetTreeGenerator<10>> mTreeGen;
    DetectorSimulation mDetector;
    SecondaryVertexFinder mSvFinder;
    DebugLevel mDebugLevel;

public:
    BjetAnalysis(std::vector<std::string> features, DebugLevel debugLevel = DebugLevel::INFO)
        : mTreeGen(std::make_unique<JetTreeGenerator<10>>("HFTagging.root", features, debugLevel)),
          mDetector(), // Default constructor
          mSvFinder(debugLevel),
          mDebugLevel(debugLevel)
    {
    }
    BjetAnalysis(std::string modelpath, std::vector<std::string> features, DebugLevel debugLevel = DebugLevel::INFO)
        : mTreeGen(std::make_unique<JetTreeGenerator<10>>(features, modelpath.c_str(), debugLevel)),
          mDetector(), // Default constructor
          mSvFinder(debugLevel),
          mDebugLevel(debugLevel)
    {
    }

    void analyze(int nEvents = 10000, double pTHatmin = 20, double pTHatmax = -1, bool writeTHnSpase = false, bool doNormalization = false)
    {
        auto logLocal = [this](DebugLevel level, const std::string &message)
        {
            log(mDebugLevel, level, message);
        };

        logLocal(DebugLevel::INFO, "Starting analysis with " + std::to_string(nEvents) + " events...");

        // Add progress tracking variables
        const int progressStep = std::max(1, nEvents / 100); // Show progress every 1%
        int lastProgress = -1;
        auto showProgress = [&](int current)
        {
            int progress = (current * 100) / nEvents;
            if (progress > lastProgress)
            {
                std::string progressBar = "[";
                for (int i = 0; i < 50; i++)
                {
                    if (i < progress / 2)
                        progressBar += "=";
                    else if (i == progress / 2)
                        progressBar += ">";
                    else
                        progressBar += " ";
                }
                progressBar += "]";

                logLocal(DebugLevel::INFO, "\rProgress: " + progressBar + " " +
                                               std::to_string(progress) + "% (" +
                                               std::to_string(current) + "/" +
                                               std::to_string(nEvents) + " events)");
                lastProgress = progress;
            }
        };

        // Create output file
        TFile *outFile = new TFile("AnalysisResults.root", "RECREATE");

        // Create histograms
        // Basic jet kinematics
        TH1F *hJetPt = new TH1F("hJetPt", "Jet p_{T};p_{T} (GeV/c);Counts", 200, 0, 200);
        TH1F *hJetEta = new TH1F("hJetEta", "Jet #eta;#eta;Counts", 100, -5, 5);

        // Flavor-specific jet pT
        TH1F *hJetPt_b = new TH1F("hJetPt_b", "b-Jet p_{T};p_{T} (GeV/c);Counts", 200, 0, 200);
        TH1F *hJetPt_c = new TH1F("hJetPt_c", "c-Jet p_{T};p_{T} (GeV/c);Counts", 200, 0, 200);

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

        // Number of constituents
        TH1F *hNConstituents = new TH1F("hNConstituents", "Number of Jet Constituents (inclusive);N_{constituents};Counts", 50, 0, 50);
        TH1F *hNConstituents_b = new TH1F("hNConstituents_b", "Number of Jet Constituents (b-jets);N_{constituents};Counts", 50, 0, 50);
        TH1F *hNConstituents_c = new TH1F("hNConstituents_c", "Number of Jet Constituents (c-jets);N_{constituents};Counts", 50, 0, 50);

        // Add these histogram declarations
        TH1F *hNSVperJet = new TH1F("hNSVperJet", "Number of SV per Jet (inclusive);N_{SV};Counts", 10, 0, 10);
        TH1F *hNSVperJet_b = new TH1F("hNSVperJet_b", "Number of SV per Jet (b-jets);N_{SV};Counts", 10, 0, 10);
        TH1F *hNSVperJet_c = new TH1F("hNSVperJet_c", "Number of SV per Jet (c-jets);N_{SV};Counts", 10, 0, 10);

        // Add this with other histogram declarations
        TH1F *hPrimaryVertexZ = new TH1F("hPrimaryVertexZ", "Primary Vertex Z Position;z (cm);Counts", 100, -10, 10);
        TH1F *hPrimaryVertexZDiff = new TH1F("hPrimaryVertexZDiff", "Primary Vertex Z Position;z (cm);Counts", 1000, -5, 5);

        TH2D *hIPJetpTN2_ljet;
        TH2D *hIPxyJetpTN3_ljet;
        TH2D *hIPzJetpTN3_ljet;
        TH2D *hIPJetpTN3_ljet;

        TH2D *hIPJetpTN2_bjet;
        TH2D *hIPxyJetpTN3_bjet;
        TH2D *hIPzJetpTN3_bjet;
        TH2D *hIPJetpTN3_bjet;

        TH2D *hIPJetpTN2_cjet;
        TH2D *hIPxyJetpTN3_cjet;
        TH2D *hIPzJetpTN3_cjet;
        TH2D *hIPJetpTN3_cjet;

        TH3D *hDecayLengthDispersionJetpT_ljet;
        TH3D *hDecayLengthDispersionJetpT_bjet;
        TH3D *hDecayLengthDispersionJetpT_cjet;

        TH2D *hCPAvsJetpT_ljet;
        TH2D *hCPAvsJetpT_bjet;
        TH2D *hCPAvsJetpT_cjet;

        // Define the number of dimensions and the binning for each dimension
        const int nDims = 7;
        int bins[nDims] = {200, 3, 500, 500, 100, 100, 100}; // Adjust bin numbers as needed
        double xmin[nDims] = {0, 0, 0, 0, 0, 0, 0};
        double xmax[nDims] = {200, 3, 1, 30, 100, 10, 1}; // Adjust ranges as needed

        THnSparseF *hJetTaggingInfo = new THnSparseF("hJetTaggingInfo", "Jet information;jetpT;jetFlavor;Score;-log(1-score);jetMass;svMass;svfE", nDims, bins, xmin, xmax);

        TH2F *hScoreVsJetPt_incl;
        TH2F *hLogScoreVsJetPt_incl;

        TH2F *hScoreVsJetPt_b;
        TH2F *hLogScoreVsJetPt_b;

        TH2F *hScoreVsJetPt_c;
        TH2F *hLogScoreVsJetPt_c;

        TH2F *hScoreVsJetPt_lf;
        TH2F *hLogScoreVsJetPt_lf;

        // 2D histograms for Score vs jet pT and -log(1-score) vs jet pT
        if (mTreeGen->getDoPrediction())
        {
            hScoreVsJetPt_incl = new TH2F("h2_score_jetpT", "Score vs Jet pT (inclusive);Jet pT (GeV/c);Score", 200, 0, 200, 240, -0.1, 1.1);
            hLogScoreVsJetPt_incl = new TH2F("h2_logscore_jetpT", "-log(1-Score) vs Jet pT (inclusive);Jet pT (GeV/c);-log(1-Score)", 200, 0, 200, 240, 0, 30);

            hScoreVsJetPt_b = new TH2F("h2_score_jetpT_bjet", "Score vs Jet pT (b-jets);Jet pT (GeV/c);Score", 200, 0, 200, 240, -0.1, 1.1);
            hLogScoreVsJetPt_b = new TH2F("h2_logscore_jetpT_bjet", "-log(1-Score) vs Jet pT (b-jets);Jet pT (GeV/c);-log(1-Score)", 200, 0, 200, 240, 0, 30);

            hScoreVsJetPt_c = new TH2F("h2_score_jetpT_cjet", "Score vs Jet pT (c-jets);Jet pT (GeV/c);Score", 200, 0, 200, 240, -0.1, 1.1);
            hLogScoreVsJetPt_c = new TH2F("h2_logscore_jetpT_cjet", "-log(1-Score) vs Jet pT (c-jets);Jet pT (GeV/c);-log(1-Score)", 200, 0, 200, 240, 0, 30);

            hScoreVsJetPt_lf = new TH2F("h2_score_jetpT_lfjet", "Score vs Jet pT (light-flavor jets);Jet pT (GeV/c);Score", 200, 0, 200, 240, -0.1, 1.1);
            hLogScoreVsJetPt_lf = new TH2F("h2_logscore_jetpT_lfjet", "-log(1-Score) vs Jet pT (light-flavor jets);Jet pT (GeV/c);-log(1-Score)", 200, 0, 200, 240, 0, 30);

            hIPJetpTN2_ljet = new TH2D("h2IPJetpTN2_lfjet", "IP vs. Jet pT (light-flavor jets)N2;Jet pT (GeV/c);IP (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPxyJetpTN3_ljet = new TH2D("h2IPxyJetpTN3_lfjet", "IP_{xy} vs. Jet pT (light-flavor jets)N3;Jet pT (GeV/c);IP_{xy} (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPzJetpTN3_ljet = new TH2D("h2IPzJetpTN3_lfjet", "IP_{z} vs. Jet pT (light-flavor jets)N3;Jet pT (GeV/c);IP_{z} (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPJetpTN3_ljet = new TH2D("h2IPJetpTN3_lfjet", "IP vs. Jet pT (light-flavor jets)N3;Jet pT (GeV/c);IP (cm)", 200, 0, 200, 600, -1.2, 1.2);

            hIPJetpTN2_bjet = new TH2D("h2IPJetpTN2_bjet", "IP vs. Jet pT (b-jets)N2;Jet pT (GeV/c);IP (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPxyJetpTN3_bjet = new TH2D("h2IPxyJetpTN3_bjet", "IP_{xy} vs. Jet pT (b-jets)N3;Jet pT (GeV/c);IP_{xy} (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPzJetpTN3_bjet = new TH2D("h2IPzJetpTN3_bjet", "IP_{z} vs. Jet pT (b-jets)N3;Jet pT (GeV/c);IP_{z} (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPJetpTN3_bjet = new TH2D("h2IPJetpTN3_bjet", "IP vs. Jet pT (b-jets)N3;Jet pT (GeV/c);IP (cm)", 200, 0, 200, 600, -1.2, 1.2);

            hIPJetpTN2_cjet = new TH2D("h2IPJetpTN2_cjet", "IP vs. Jet pT (c-jets)N2;Jet pT (GeV/c);IP (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPxyJetpTN3_cjet = new TH2D("h2IPxyJetpTN3_cjet", "IP_{xy} vs. Jet pT (c-jets)N3;Jet pT (GeV/c);IP_{xy} (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPzJetpTN3_cjet = new TH2D("h2IPzJetpTN3_cjet", "IP_{z} vs. Jet pT (c-jets)N3;Jet pT (GeV/c);IP_{z} (cm)", 200, 0, 200, 600, -1.2, 1.2);
            hIPJetpTN3_cjet = new TH2D("h2IPJetpTN3_cjet", "IP vs. Jet pT (c-jets)N3;Jet pT (GeV/c);IP (cm)", 200, 0, 200, 600, -1.2, 1.2);

            hDecayLengthDispersionJetpT_ljet = new TH3D("h3DecayLengthDispersionJetpT_lfjet", "Decay Length vs. Dispersion vs. Jet pT (light-flavor jets);Jet pT (GeV/c);Decay Length (cm);Dispersion", 200, 0, 200, 600, 0, 30, 200, 0, 2);
            hDecayLengthDispersionJetpT_bjet = new TH3D("h3DecayLengthDispersionJetpT_bjet", "Decay Length vs. Dispersion vs. Jet pT (b-jets);Jet pT (GeV/c);Decay Length (cm);Dispersion", 200, 0, 200, 600, 0, 30, 200, 0, 2);
            hDecayLengthDispersionJetpT_cjet = new TH3D("h3DecayLengthDispersionJetpT_cjet", "Decay Length vs. Dispersion vs. Jet pT (c-jets);Jet pT (GeV/c);Decay Length (cm);Dispersion", 200, 0, 200, 600, 0, 30, 200, 0, 2);

            hCPAvsJetpT_ljet = new TH2D("hCPAvsJetpT_lfjet", "CPA vs. Jet pT (light-flavor jets);Jet pT (GeV/c);CPA", 200, 0, 200, 210, -1.05, 1.05);
            hCPAvsJetpT_bjet = new TH2D("hCPAvsJetpT_bjet", "CPA vs. Jet pT (b-jets);Jet pT (GeV/c);CPA", 200, 0, 200, 210, -1.05, 1.05);
            hCPAvsJetpT_cjet = new TH2D("hCPAvsJetpT_cjet", "CPA vs. Jet pT (c-jets);Jet pT (GeV/c);CPA", 200, 0, 200, 210, -1.05, 1.05);
        }

        // Initialize Pythia with more detailed error checking
        logLocal(DebugLevel::INFO, "Initializing Pythia...");
        Pythia8::Pythia pythia;

        // Set random seed to current time modulo max allowed value
        // std::time_t currentTime = std::time(nullptr);
        // int pythiaSeed = currentTime % 900000000;  // Ensure seed is in valid range
        // pythia.readString("Random:setSeed = on");
        // pythia.readString("Random:seed = " + std::to_string(pythiaSeed));
        // std::cout << "Set Pythia random seed to: " << pythiaSeed << std::endl;

        if (gSystem->Getenv("CONFIG_SEED"))
        {
            int64_t pythiaSeed = atoi(gSystem->Getenv("CONFIG_SEED"));
            pythia.readString("Random:setSeed = on");
            pythia.readString("Random:seed = " + std::to_string(pythiaSeed % 900000000));
            logLocal(DebugLevel::INFO, "Set Pythia random seed to: " + std::to_string(pythiaSeed));
        }

        // Basic settings
        pythia.readString("Beams:eCM = 13000.");
        logLocal(DebugLevel::DEBUG, "Set beam energy");

        pythia.readString("Beams:allowVertexSpread = on");
        logLocal(DebugLevel::DEBUG, "Enabled vertex spread");

        // More reasonable vertex spread values (in mm)
        pythia.readString("Beams:sigmaVertexX = 10"); // Back to mm (= 1 cm)
        pythia.readString("Beams:sigmaVertexY = 10"); // Back to mm (= 1 cm)
        pythia.readString("Beams:sigmaVertexZ = 50"); // Back to mm (= 5 cm)
        logLocal(DebugLevel::DEBUG, "Set vertex spread parameters");

        pythia.readString("HardQCD:all = on");                           // Turn on hard QCD processes
        pythia.readString(Form("PhaseSpace:pTHatMin = %.1f", pTHatmin)); // Min pT for hard interaction
        pythia.readString(Form("PhaseSpace:pTHatMax = %.1f", pTHatmax)); // Max pT for hard interaction
        logLocal(DebugLevel::DEBUG, "Set physics process parameters");

        if (!pythia.init())
        {
            logLocal(DebugLevel::ERROR, "Pythia initialization failed!");
            return;
        }
        logLocal(DebugLevel::INFO, "Pythia initialized successfully");

        // Initialize FastJet
        double R = 0.4; // jet radius parameter
        fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);

        // Initialize detector simulation, vertex finder, and DCA smearing
        mDetector = DetectorSimulation(mDebugLevel);
        DCASmearing dcaSmearer(mDebugLevel); // Create instance of DCASmearing

        // Initialize histograms for DCA smearing
        dcaSmearer.initializeHistograms();

        // Event loop with progress
        long totalJets = 0;
        long totalEvents = 0;

        auto startTime = std::chrono::high_resolution_clock::now();

        for (int iEvent = 0; iEvent < nEvents; ++iEvent)
        {
            // Show progress
            if (iEvent % progressStep == 0)
            {
                showProgress(iEvent);

                // Show timing information
                if (iEvent > 0)
                {
                    auto currentTime = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime);
                    double eventsPerSecond = static_cast<double>(iEvent) / duration.count();
                    int remainingSeconds = static_cast<int>((nEvents - iEvent) / eventsPerSecond);

                    logLocal(DebugLevel::DEBUG,
                             "Processing speed: " + std::to_string(eventsPerSecond) + " events/s, " +
                                 "ETA: " + std::to_string(remainingSeconds) + "s");
                }
            }

            logLocal(DebugLevel::DEBUG, "=== Processing event " + std::to_string(iEvent) + " ===");

            if (!pythia.next())
            {
                logLocal(DebugLevel::ERROR, "Pythia event generation failed");
                continue;
            }

            // Collect tracks
            logLocal(DebugLevel::DEBUG, "Collecting tracks...");
            std::vector<Track> primaryTracks;
            std::vector<JParticle> mcparticles;

            // Process particles
            logLocal(DebugLevel::DEBUG, "Processing particles...");
            std::vector<fastjet::PseudoJet> particles;
            std::vector<Track> allTracks;
            std::vector<Track> truthTracks;

            // First, get true primary vertex from Pythia
            TVector3 truePrimaryVertex(pythia.event[5].xProd() * 0.1, // Convert mm to cm
                                       pythia.event[5].yProd() * 0.1,
                                       pythia.event[5].zProd() * 0.1);

            // Smear primary vertex
            TVector3 smearedPrimaryVertex = mDetector.smearPrimaryVertex(truePrimaryVertex);

            // Now process particles with respect to this vertex
            for (int i = 0; i < pythia.event.size(); ++i)
            {
                const Pythia8::Particle &part = pythia.event[i];

                if (mDebugLevel >= DebugLevel::VERBOSE)
                {
                    logLocal(DebugLevel::VERBOSE, "Processing particle " + std::to_string(i) +
                                                      " (ID: " + std::to_string(part.id()) +
                                                      ", pT: " + std::to_string(part.pT()) +
                                                      ", eta: " + std::to_string(part.eta()) +
                                                      ", status: " + std::to_string(part.status()) +
                                                      ", charge: " + std::to_string(part.charge()) +
                                                      ", zProd: " + std::to_string(part.zProd()) + " mm)");
                }

                // Parton Definition
                if (std::abs(part.id()) == 5 || std::abs(part.id()) == 4)
                {
                    mcparticles.push_back(convertToJParticle(part));
                }

                // Select only final state charged particles
                if (!part.isFinal() || part.charge() == 0)
                {
                    logLocal(DebugLevel::DEBUG, "Skipping particle: not final or not charged");
                    continue;
                }

                // After detector simulation
                logLocal(DebugLevel::DEBUG, "Applying detector simulation...");
                Track trackParams = Track();
                mDetector.smearTrack(part, trackParams, truePrimaryVertex);
                logLocal(DebugLevel::DEBUG, "Smeared z position: " + std::to_string(trackParams.getPosition().Z()) + " cm");

                // Check mother
                int mother1 = part.mother1();
                if (mother1 >= 0 && mother1 < pythia.event.size())
                {
                    const Pythia8::Particle &mother = pythia.event[mother1];
                    logLocal(DebugLevel::DEBUG, "Mother particle ID: " + std::to_string(mother.id()) +
                                                    ", isHadron: " + std::to_string(mother.isHadron()) +
                                                    ", tau: " + std::to_string(mother.tau()));

                    // Check if particle comes from hadron decay
                    if (!(mother.isHadron() && mother.tau() > 0))
                    {
                        primaryTracks.push_back(trackParams);
                        logLocal(DebugLevel::DEBUG, "Added to primary tracks");
                    }
                }

                // Basic acceptance cuts
                if (abs(part.eta()) > 0.9 || part.pT() < 0.15)
                {
                    logLocal(DebugLevel::DEBUG, "Skipping: outside acceptance");
                    continue;
                }

                allTracks.push_back(trackParams);
                allTracks.back().setIndex(allTracks.size() - 1);

                double xPos = part.xProd() * 0.1;
                double yPos = part.yProd() * 0.1;
                double zPos = part.zProd() * 0.1;
                double trueDCAz = zPos - truePrimaryVertex.Z();
                double trueXYdca = std::copysign(std::sqrt(std::pow(xPos - truePrimaryVertex.X(), 2) + std::pow(yPos - truePrimaryVertex.Y(), 2)), yPos - truePrimaryVertex.Y());

                // Create truth track
                Track truthTrack(TVector3(xPos, yPos, zPos),
                                 TVector3(part.px(), part.py(), part.pz()),
                                 trueXYdca,
                                 trueDCAz,
                                 part.charge());
                truthTracks.push_back(truthTrack);
                truthTracks.back().setIndex(truthTracks.size() - 1);

                // Add to jet finding
                particles.push_back(fastjet::PseudoJet(
                    trackParams.getMomentum().X(),
                    trackParams.getMomentum().Y(),
                    trackParams.getMomentum().Z(),
                    trackParams.getMomentum().Mag()));
                particles.back().set_user_index(allTracks.size() - 1);
            }

            // Find primary vertex with additional checks
            logLocal(DebugLevel::DEBUG, "Finding primary vertex...");
            logLocal(DebugLevel::DEBUG, "Number of primary tracks: " + std::to_string(primaryTracks.size()));

            if (primaryTracks.empty())
            {
                logLocal(DebugLevel::WARNING, "No primary tracks found, skipping vertex finding");
                continue;
            }

            if (primaryTracks.size() < 3)
            {
                logLocal(DebugLevel::WARNING, "Too few tracks for vertex finding: " + std::to_string(primaryTracks.size()));
                continue;
            }

            try
            {
                logLocal(DebugLevel::DEBUG, "\nStarting vertex finding process...");

                if (primaryTracks.empty())
                {
                    logLocal(DebugLevel::WARNING, "No primary tracks available for vertex finding");
                    continue;
                }

                TVector3 primaryVertex = VertexFinder::findPrimaryVertex(primaryTracks,
                                                                         0.001,        // 10 micron convergence
                                                                         50,           // max iterations
                                                                         1.0,          // temperature
                                                                         6.0,          // tighter chi2 cut
                                                                         mDebugLevel); // Pass debug level

                logLocal(DebugLevel::DEBUG, "Primary vertex found at: (" + std::to_string(primaryVertex.X()) + ", " + std::to_string(primaryVertex.Y()) + ", " + std::to_string(primaryVertex.Z()) + ")");

                float vertexDifference = smearedPrimaryVertex.Z() - primaryVertex.Z();
                hPrimaryVertexZDiff->Fill(vertexDifference);

                // Safety check for vertex position
                if (std::isnan(primaryVertex.X()) || std::isnan(primaryVertex.Y()) || std::isnan(primaryVertex.Z()))
                {
                    logLocal(DebugLevel::WARNING, "Warning: Invalid primary vertex position detected");
                    continue;
                }

                // Fill primary vertex Z position
                hPrimaryVertexZ->Fill(primaryVertex.Z());

                // Vertex position check
                if (std::abs(primaryVertex.Z()) > 10.0 ||
                    std::abs(primaryVertex.X()) > 2.0 ||
                    std::abs(primaryVertex.Y()) > 2.0)
                {
                    logLocal(DebugLevel::WARNING, "WARNING: Vertex position outside range!");
                    continue;
                }

                // Jet finding with safety checks
                logLocal(DebugLevel::DEBUG, "\nStarting jet finding...");

                if (particles.empty())
                {
                    logLocal(DebugLevel::WARNING, "No particles available for jet finding");
                    continue;
                }

                // Declare variables before try block
                std::vector<fastjet::PseudoJet> jets;
                std::unique_ptr<fastjet::ClusterSequence> cs = nullptr;

                // Run clustering with try-catch
                try
                {
                    cs = std::make_unique<fastjet::ClusterSequence>(particles, jetDef);
                    jets = fastjet::sorted_by_pt(cs->inclusive_jets(5.0));

                    logLocal(DebugLevel::DEBUG, "\nFound " + std::to_string(jets.size()) + " jets");

                    // Process jets with safety checks
                    for (const auto &jet : jets)
                    {
                        // Check if jet is valid using FastJet's recommended method
                        if (jet.E() <= 0 || std::isnan(jet.E()))
                        {
                            logLocal(DebugLevel::WARNING, "Invalid jet detected (E = " + std::to_string(jet.E()) + ")");
                            continue;
                        }

                        if (jet.eta() > 0.5 || jet.eta() < -0.5)
                        {
                            logLocal(DebugLevel::DEBUG, "Jet outside acceptance (eta = " + std::to_string(jet.eta()) + ")");
                            continue;
                        }

                        // Get jet flavor with safety checks
                        int jetFlavor = getJetFlavor(jet, mcparticles, R, mDebugLevel);

                        logLocal(DebugLevel::DEBUG, "Jet: pT=" + std::to_string(jet.pt()) +
                                                        " eta=" + std::to_string(jet.eta()) +
                                                        " phi=" + std::to_string(jet.phi()) +
                                                        " flavor=" + std::to_string(jetFlavor));

                        // Get constituents with safety checks
                        std::vector<fastjet::PseudoJet> constituents = jet.constituents();

                        // Fill constituent multiplicity histograms
                        if (hNConstituents)
                        {
                            hNConstituents->Fill(constituents.size());
                            if (jetFlavor == JetTaggingSpecies::beauty && hNConstituents_b)
                            {
                                hNConstituents_b->Fill(constituents.size());
                            }
                            else if (jetFlavor == JetTaggingSpecies::charm && hNConstituents_c)
                            {
                                hNConstituents_c->Fill(constituents.size());
                            }
                        }

                        // Collect tracks associated with this jet
                        std::vector<Track> jetTracks;
                        for (const auto &constituent : constituents)
                        {
                            size_t trackIndex = constituent.user_index();
                            if (trackIndex < allTracks.size())
                            {
                                jetTracks.push_back(allTracks[trackIndex]);
                            }
                        }

                        // Find all possible secondary vertices
                        std::vector<SecondaryVertexFinder::SecVtxInfo> secondaryVertices = mSvFinder.findSecondaryVertices(jetTracks);

                        mSvFinder.fillMonitoringHistograms(secondaryVertices);

                        // Fill SV-related histograms
                        if (!secondaryVertices.empty())
                        {
                            // Number of SVs per jet
                            if (hNSVperJet)
                                hNSVperJet->Fill(secondaryVertices.size());
                            if (jetFlavor == JetTaggingSpecies::beauty && hNSVperJet_b)
                            {
                                hNSVperJet_b->Fill(secondaryVertices.size());
                            }
                            else if (jetFlavor == JetTaggingSpecies::charm && hNSVperJet_c)
                            {
                                hNSVperJet_c->Fill(secondaryVertices.size());
                            }

                            // Process each secondary vertex
                            for (const auto &sv : secondaryVertices)
                            {
                                double mass = sv.momentum.M();
                                double decay2D = SVCalculations::calculateDecayLength2D(primaryVertex, sv.position);
                                double decay3D = SVCalculations::calculateDecayLength3D(primaryVertex, sv.position);

                                hSVMass->Fill(mass);
                                hDecayLength2D->Fill(decay2D);
                                hDecayLength3D->Fill(decay3D);

                                // Fill flavor-specific histograms
                                if (jetFlavor == JetTaggingSpecies::beauty)
                                {
                                    hSVMass_b->Fill(mass);
                                    hDecayLength2D_b->Fill(decay2D);
                                    hDecayLength3D_b->Fill(decay3D);
                                }
                                else if (jetFlavor == JetTaggingSpecies::charm)
                                {
                                    hSVMass_c->Fill(mass);
                                    hDecayLength2D_c->Fill(decay2D);
                                    hDecayLength3D_c->Fill(decay3D);
                                }

                                logLocal(DebugLevel::VERBOSE, "SV in jet: mass=" + std::to_string(mass) +
                                                                  " decay2D=" + std::to_string(decay2D) +
                                                                  " decay3D=" + std::to_string(decay3D));
                            }
                        }

                        // Fill other jet-related histograms
                        if (hJetPt)
                            hJetPt->Fill(jet.pt());
                        if (hJetEta)
                            hJetEta->Fill(jet.eta());

                        if (jetFlavor == JetTaggingSpecies::beauty)
                        {
                            if (hJetPt_b)
                                hJetPt_b->Fill(jet.pt());
                        }
                        else if (jetFlavor == JetTaggingSpecies::charm)
                        {
                            if (hJetPt_c)
                                hJetPt_c->Fill(jet.pt());
                        }

                        // Process jets
                        float score = 0.0;
                        if (mTreeGen->getDoPrediction())
                        {
                            auto [jetparams, trackparams, svparams] = mTreeGen->processJetFeatures(jet, allTracks, secondaryVertices, primaryVertex);
                            auto output = mTreeGen->predict(jetparams, trackparams, svparams);
                            score = output[0];

                            // Fill the THnSparse histogram
                            double logScore = -log(1 - score);
                            double svMass = secondaryVertices.empty() ? 0 : secondaryVertices[0].momentum.M();
                            double svfE = secondaryVertices.empty() ? 0 : secondaryVertices[0].momentum.E() / jet.e();
                            double values[nDims] = {jet.pt(), static_cast<double>(jetFlavor), score, logScore, jet.m(), svMass, svfE};

                            // Fill the THnSparse histogram
                            if (writeTHnSpase)
                            {
                                hJetTaggingInfo->Fill(values);
                            }

                            // Fill 2D histograms
                            hScoreVsJetPt_incl->Fill(jet.pt(), score);
                            hLogScoreVsJetPt_incl->Fill(jet.pt(), logScore);

                            if (jetFlavor == JetTaggingSpecies::beauty)
                            {
                                hScoreVsJetPt_b->Fill(jet.pt(), score);
                                hLogScoreVsJetPt_b->Fill(jet.pt(), logScore);

                                if (trackparams[1].trackpT > 0)
                                {
                                    hIPJetpTN2_bjet->Fill(jet.pt(), std::copysign(std::sqrt(trackparams[1].signedIP3D * trackparams[1].signedIP3D + trackparams[1].signedIP2D * trackparams[1].signedIP2D), trackparams[1].signedIP2D));
                                }
                                if (trackparams[2].trackpT > 0)
                                {
                                    hIPJetpTN3_bjet->Fill(jet.pt(), std::copysign(std::sqrt(trackparams[2].signedIP3D * trackparams[2].signedIP3D + trackparams[2].signedIP2D * trackparams[2].signedIP2D),trackparams[2].signedIP2D));
                                    hIPzJetpTN3_bjet->Fill(jet.pt(), trackparams[2].signedIP3D);
                                    hIPxyJetpTN3_bjet->Fill(jet.pt(), trackparams[2].signedIP2D);
                                }
                                if (svparams[0].svPt > 0)
                                {
                                    hDecayLengthDispersionJetpT_bjet->Fill(jet.pt(), svparams[0].svDecayLength3D, svparams[0].svDispersion);
                                    hCPAvsJetpT_bjet->Fill(jet.pt(), svparams[0].svCPA);
                                }
                            }
                            else if (jetFlavor == JetTaggingSpecies::charm)
                            {
                                hScoreVsJetPt_c->Fill(jet.pt(), score);
                                hLogScoreVsJetPt_c->Fill(jet.pt(), logScore);

                                if (trackparams[1].trackpT > 0)
                                {
                                    hIPJetpTN2_cjet->Fill(jet.pt(), std::copysign(std::sqrt(trackparams[1].signedIP3D * trackparams[1].signedIP3D + trackparams[1].signedIP2D * trackparams[1].signedIP2D), trackparams[1].signedIP2D));
                                }
                                if (trackparams[2].trackpT > 0)
                                {
                                    hIPJetpTN3_cjet->Fill(jet.pt(), std::copysign(std::sqrt(trackparams[2].signedIP3D * trackparams[2].signedIP3D + trackparams[2].signedIP2D * trackparams[2].signedIP2D), trackparams[2].signedIP2D));
                                    hIPzJetpTN3_cjet->Fill(jet.pt(), trackparams[2].signedIP3D);
                                    hIPxyJetpTN3_cjet->Fill(jet.pt(), trackparams[2].signedIP2D);
                                }
                                if (svparams[0].svPt > 0)
                                {
                                    hDecayLengthDispersionJetpT_cjet->Fill(jet.pt(), svparams[0].svDecayLength3D, svparams[0].svDispersion);
                                    hCPAvsJetpT_cjet->Fill(jet.pt(), svparams[0].svCPA);
                                }
                            }
                            else
                            {
                                hScoreVsJetPt_lf->Fill(jet.pt(), score);
                                hLogScoreVsJetPt_lf->Fill(jet.pt(), logScore);

                                if (trackparams[1].trackpT > 0)
                                {
                                    hIPJetpTN2_ljet->Fill(jet.pt(), std::copysign(std::sqrt(trackparams[1].signedIP3D * trackparams[1].signedIP3D + trackparams[1].signedIP2D * trackparams[1].signedIP2D), trackparams[1].signedIP2D));
                                }
                                if (trackparams[2].trackpT > 0)
                                {
                                    hIPJetpTN3_ljet->Fill(jet.pt(), std::copysign(std::sqrt(trackparams[2].signedIP3D * trackparams[2].signedIP3D + trackparams[2].signedIP2D * trackparams[2].signedIP2D), trackparams[2].signedIP2D));
                                    hIPzJetpTN3_ljet->Fill(jet.pt(), trackparams[2].signedIP3D);
                                    hIPxyJetpTN3_ljet->Fill(jet.pt(), trackparams[2].signedIP2D);
                                }
                                if (svparams[0].svPt > 0)
                                {
                                    hDecayLengthDispersionJetpT_ljet->Fill(jet.pt(), svparams[0].svDecayLength3D, svparams[0].svDispersion);
                                    hCPAvsJetpT_ljet->Fill(jet.pt(), svparams[0].svCPA);
                                }
                            }
                        }
                        else
                        {
                            mTreeGen->processJet(jet, allTracks, mcparticles, secondaryVertices, primaryVertex, R);
                        }
                    }

                    // Update counters safely
                    if (totalEvents < std::numeric_limits<long>::max())
                    {
                        totalEvents++;
                        totalJets += jets.size();
                    }
                }
                catch (const fastjet::Error &e)
                {
                    logLocal(DebugLevel::ERROR, "FastJet error: " + std::string(e.message()));
                }
                catch (const std::exception &e)
                {
                    logLocal(DebugLevel::ERROR, "Standard exception in jet finding: " + std::string(e.what()));
                }
                catch (...)
                {
                    logLocal(DebugLevel::ERROR, "Unknown exception in jet finding");
                }

                // Clean up
                if (cs)
                {
                    cs.reset();
                }

                // Remove duplicate jet finding code
                if (iEvent % 100 == 0)
                {
                    logLocal(DebugLevel::INFO, "Processed " + std::to_string(iEvent) + " events");
                }

                // Apply smearing to all tracks
                for (size_t i = 0; i < allTracks.size(); i++)
                {
                    if (i < truthTracks.size())
                    { // Add safety check
                        dcaSmearer.fillHistograms(truthTracks[i], allTracks[i]);
                    }
                    else
                    {
                        logLocal(DebugLevel::WARNING, "Truth track index out of range");
                    }
                }

                logLocal(DebugLevel::DEBUG, "Event " + std::to_string(iEvent) + " processed successfully");
            }
            catch (const std::exception &e)
            {
                logLocal(DebugLevel::ERROR, "Exception caught: " + std::string(e.what()));
            }
            catch (...)
            {
                logLocal(DebugLevel::ERROR, "Unknown exception caught");
            }
        }

        // Show final progress
        showProgress(nEvents);

        // Show final statistics
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);

        logLocal(DebugLevel::INFO, "\nProcessing completed:");
        logLocal(DebugLevel::INFO, "Total time: " + std::to_string(duration.count()) + " seconds");
        logLocal(DebugLevel::INFO, "Average processing speed: " +
                                       std::to_string(static_cast<double>(nEvents) / duration.count()) + " events/s");

        if (totalEvents > 0)
        {
            double avgJets = static_cast<double>(totalJets) / totalEvents;
            logLocal(DebugLevel::INFO, "Average number of jets per event: " + std::to_string(avgJets));
        }

        // Write histograms with safety checks
        logLocal(DebugLevel::INFO, "\nWriting histograms...");
        if (outFile && outFile->IsOpen())
        {
            if ((pTHatmax > 0 || pTHatmin > 20) && doNormalization)
            {
                hJetPt->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hJetPt_b->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hJetPt_c->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hNSVperJet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hNSVperJet_b->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hNSVperJet_c->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hSVMass->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hSVMass_b->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hSVMass_c->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hDecayLength2D->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hDecayLength2D_b->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hDecayLength2D_c->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hDecayLength3D->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hDecayLength3D_b->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hDecayLength3D_c->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hNConstituents->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hNConstituents_b->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                hNConstituents_c->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                if (mTreeGen->getDoPrediction())
                {
                    if (writeTHnSpase)
                    {
                        hJetTaggingInfo->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    }
                    hScoreVsJetPt_incl->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hLogScoreVsJetPt_incl->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hScoreVsJetPt_b->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hLogScoreVsJetPt_b->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hScoreVsJetPt_c->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hLogScoreVsJetPt_c->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hScoreVsJetPt_lf->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hLogScoreVsJetPt_lf->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPJetpTN2_ljet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPJetpTN2_bjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPJetpTN2_cjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPxyJetpTN3_ljet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPxyJetpTN3_bjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPxyJetpTN3_cjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPzJetpTN3_ljet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPzJetpTN3_bjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPzJetpTN3_cjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPJetpTN3_ljet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPJetpTN3_bjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hIPJetpTN3_cjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hDecayLengthDispersionJetpT_ljet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hDecayLengthDispersionJetpT_bjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hDecayLengthDispersionJetpT_cjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hCPAvsJetpT_ljet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hCPAvsJetpT_bjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                    hCPAvsJetpT_cjet->Scale(1.0 * pythia.info.sigmaGen() / pythia.info.nTried());
                }
            }

            outFile->cd();
            if (mTreeGen->getDoPrediction() && writeTHnSpase)
            {
                hJetTaggingInfo->Write(); // Write the THnSparse histogram
            }
            outFile->Write();
            outFile->Close();
        }

        TFile *xsecFile = new TFile("xsection.root", "RECREATE");
        xsecFile->cd();
        TProfile *pCrossSection = new TProfile("xsection", "Cross Section Profile; pT (GeV); Cross Section (pb)", 1, 0, 1);
        TH1F *hTrials = new TH1F("hTrials", "Number of Trials; Trials; Counts", 1, 0, 1);

        pCrossSection->Fill(0.5, pythia.info.sigmaGen());
        hTrials->Fill(0.5, pythia.info.nTried());

        xsecFile->Write();
        xsecFile->Close();

        // Write tree
        if (!mTreeGen->getDoPrediction())
        {
            mTreeGen->write();
        }

        logLocal(DebugLevel::INFO, "Analysis completed successfully!");
    }
};

void bjet_analysis(int nEvents = 10000, double pTHatMin = 20, double pTHatMax = -1, std::string modelpath = "", DebugLevel debugLevel = DebugLevel::INFO)
{
    // Initialize features
    std::vector<std::string> features = {"jetpT", "jetEta", "jetPhi", "jetMass", "nTracks", "nSV", "jetFlavor",
                                         "trackpT", "trackEta", "dotProdTrackJet", "dotProdTrackJetOverJet",
                                         "deltaRJetTrack", "signedIP2D", "signedIP3D", "momFraction", "deltaRTrackVertex",
                                         "svPt", "deltaRSVJet", "svMass", "svfE", "svIPxy", "svCPA", "svChi2PCA",
                                         "svDispersion", "svDecayLength2D", "svDecayLength3D"};

    if (!modelpath.empty())
    {
        features.erase(features.begin() + 6); // Remove jetFlavor
        BjetAnalysis bjetAnalysis(modelpath, features, debugLevel);
        bjetAnalysis.analyze(nEvents, pTHatMin, pTHatMax);
    }
    else
    {
        BjetAnalysis bjetAnalysis(features, debugLevel);
        bjetAnalysis.analyze(nEvents, pTHatMin, pTHatMax);
    }
}

// Update main to accept number of events
int main(int argc, char *argv[])
{
    int nEvents = 10000; // default value
    double pTHatMin = 20;
    double pTHatMax = -1;
    std::string modelpath = "";
    DebugLevel debugLevel = DebugLevel::INFO;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-n" || arg == "--nevents")
        {
            if (i + 1 < argc)
            {
                nEvents = std::atoi(argv[++i]);
                if (nEvents <= 0)
                {
                    std::cerr << "Error: Number of events must be positive\n";
                    return 1;
                }
            }
        }
        else if (arg == "-d" || arg == "--debug")
        {
            if (i + 1 < argc)
            {
                int level = std::atoi(argv[++i]);
                if (level >= 0 && level <= static_cast<int>(DebugLevel::VERBOSE))
                {
                    debugLevel = static_cast<DebugLevel>(level);
                }
            }
        }
        else if (arg == "-pmin" || arg == "--pTHatMin")
        {
            if (i + 1 < argc)
            {
                pTHatMin = std::atoi(argv[++i]);
            }
        }
        else if (arg == "-pmax" || arg == "--pTHatMax")
        {
            if (i + 1 < argc)
            {
                pTHatMax = std::atoi(argv[++i]);
            }
        }
        else if (arg == "-ml" || arg == "--modelPath")
        {
            if (i + 1 < argc)
            {
                modelpath = std::string(argv[++i]);
            }
        }
        else if (arg == "-h" || arg == "--help")
        {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  -n, --nevents N           Number of events to process (default: 10000)\n"
                      << "  -d, --debug N             Debug level (0-5, default: 3)\n"
                      << "  -pmin, --pTHatMin pt      Min pT of the hard scattering\n"
                      << "  -pmax, --pTHatMax pt      Max pT of the hard scattering\n"
                      << "  -ml, --modelPath path     ML model path\n"
                      << "  -h, --help                Show this help message\n";
            return 0;
        }
    }

    bjet_analysis(nEvents, pTHatMin, pTHatMax, modelpath, debugLevel);
    return 0;
}
