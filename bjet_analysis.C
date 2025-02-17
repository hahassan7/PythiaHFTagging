#include "TaggingUtilities.h"
#include "JTreeHFFeatures.h"
#include "JetTreeGenerator.h"

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

    void analyze(int nEvents = 10000)
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

        // Initialize Pythia with more detailed error checking
        logLocal(DebugLevel::INFO, "Initializing Pythia...");
        Pythia8::Pythia pythia;

        // Set random seed to current time modulo max allowed value
        // std::time_t currentTime = std::time(nullptr);
        // int pythiaSeed = currentTime % 900000000;  // Ensure seed is in valid range
        // pythia.readString("Random:setSeed = on");
        // pythia.readString("Random:seed = " + std::to_string(pythiaSeed));
        // std::cout << "Set Pythia random seed to: " << pythiaSeed << std::endl;

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

        pythia.readString("HardQCD:all = on");          // Turn on hard QCD processes
        pythia.readString("PhaseSpace:pTHatMin = 20."); // Min pT for hard interaction
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
                double trueXYdca = std::copysign(std::sqrt(std::pow(xPos - truePrimaryVertex.X(), 2) + std::pow(yPos - truePrimaryVertex.Y(), 2)), yPos);

                // Create truth track
                Track truthTrack(TVector3(xPos, yPos, zPos),
                                 TVector3(part.px(), part.py(), part.pz()),
                                 trueXYdca,
                                 trueDCAz,
                                 part.charge());
                truthTracks.push_back(truthTrack);

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
                            int trackIndex = constituent.user_index();
                            if (trackIndex >= 0 && trackIndex < allTracks.size())
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
                        mTreeGen->processJet(jet, allTracks, mcparticles, secondaryVertices,
                                             primaryVertex, R);
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
            outFile->Write();
            outFile->Close();
        }

        // Write tree
        mTreeGen->write();

        logLocal(DebugLevel::INFO, "Analysis completed successfully!");
    }
};

// Update main to accept number of events
int main(int argc, char *argv[])
{
    int nEvents = 10000; // default value
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
        else if (arg == "-h" || arg == "--help")
        {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  -n, --nevents N    Number of events to process (default: 10000)\n"
                      << "  -d, --debug N      Debug level (0-5, default: 3)\n"
                      << "  -h, --help         Show this help message\n";
            return 0;
        }
    }

    // Initialize features
    std::vector<std::string> features = {"JetpT", "JetEta", "JetPhi", "JetMass", "NTracks", "NSV", "JetFlavor",
                                         "TrackpT", "TrackEta", "DotProdTrackJet", "DotProdTrackJetOverJet", "DeltaRJetTrack", "SignedIP2D", "SignedIP3D",
                                         "MomFraction", "DeltaRTrackVertex",
                                         "SVpT", "DeltaRSVJet", "SVMass", "SVfE", "IPxy", "CPA", "Chi2PCA",
                                         "Dispersion", "DecayLength2D", "DecayLength3D"};

    BjetAnalysis bjetAnalysis(features, debugLevel);
    bjetAnalysis.analyze(nEvents);
    return 0;
}
