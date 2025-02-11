#include "TaggingUtilities.h"
#include "JTreeHFFeatures.h"

void bjet_analysis()
{
    std::cout << "Starting analysis..." << std::endl;

    // Initialize features
    std::cout << "Initializing features..." << std::endl;
    std::vector<std::string> features = {
        "JetpT", "JetEta", "JetPhi", "JetMass", "JetFlavor",
        "TrackpT", "TrackEta", "SignedIP2D", "SignedIP3D",
        "SVMass", "DecayLength2D", "DecayLength3D"};

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
    std::cout << "Initializing Pythia..." << std::endl;
    Pythia8::Pythia pythia;
    
    // Basic settings
    pythia.readString("Beams:eCM = 13000.");
    std::cout << "Set beam energy" << std::endl;
    
    pythia.readString("Beams:allowVertexSpread = on");
    std::cout << "Enabled vertex spread" << std::endl;
    
    // More reasonable vertex spread values (in mm)
    pythia.readString("Beams:sigmaVertexX = 0.015");  // 15 μm
    pythia.readString("Beams:sigmaVertexY = 0.015");  // 15 μm
    pythia.readString("Beams:sigmaVertexZ = 0.050");  // 50 μm
    std::cout << "Set vertex spread parameters" << std::endl;
    
    pythia.readString("HardQCD:all = on");          // Turn on hard QCD processes
    pythia.readString("PhaseSpace:pTHatMin = 20."); // Min pT for hard interaction
    std::cout << "Set physics process parameters" << std::endl;

    if (!pythia.init()) {
        std::cout << "ERROR: Pythia initialization failed!" << std::endl;
        return;
    }
    std::cout << "Pythia initialized successfully" << std::endl;

    // Initialize FastJet
    double R = 0.4; // jet radius parameter
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);

    // Initialize detector simulation and vertex finder
    DetectorSimulation detector;
    SecondaryVertexFinder svFinder;

    // Event loop
    int nEvents = 10000;
    long totalJets = 0;
    long totalEvents = 0;
    constexpr int maxConst = 13;

    JTreeHFFeatures<maxConst> treeWriter("bjet_features.root", features);

    // Create random number generator
    TRandom3 random(0); // seed with 0 for reproducibility

    // Define detector resolution (in cm)
    const double spatialResolution = 0.005; // 50 microns

    // Initialize histograms
    DCASmearing dcaSmearer;
    dcaSmearer.initializeHistograms();

    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        std::cout << "\n=== Processing event " << iEvent << " ===" << std::endl;
        
        if (!pythia.next())
        {
            std::cout << "Pythia event generation failed" << std::endl;
            continue;
        }

        // Collect tracks
        std::cout << "Collecting tracks..." << std::endl;
        std::vector<Track> primaryTracks;
        std::vector<TParticle> mcparticles;
        
        // Process particles
        std::cout << "Processing particles..." << std::endl;
        std::vector<fastjet::PseudoJet> particles;
        std::vector<Track> allTracks;
        std::vector<Track> truthTracks;

        for (int i = 0; i < pythia.event.size(); ++i)
        {
            const Pythia8::Particle &part = pythia.event[i];
            
            std::cout << "Processing particle " << i 
                      << " (ID: " << part.id() 
                      << ", pT: " << part.pT()
                      << ", eta: " << part.eta()
                      << ", zProd: " << part.zProd() << " mm)" << std::endl;

            // Parton Definition
            if ((part.id() == 5 || part.id() == -5))
            {
                mcparticles.push_back(convertToTParticle(part));
            }
            else if ((part.id() == 4 || part.id() == -4))
            {
                mcparticles.push_back(convertToTParticle(part));
            }

            // Select only final state charged particles
            if (!part.isFinal() || part.charge() == 0)
                continue;

            // After detector simulation
            std::cout << "Applying detector simulation..." << std::endl;
            auto trackParams = detector.smearTrack(part);
            std::cout << "Smeared z position: " << trackParams.pos.Z() << " cm" << std::endl;

            // Check mother
            const Pythia8::Particle &mother = pythia.event[part.mother1()];
            std::cout << "Mother particle ID: " << mother.id() 
                      << ", isHadron: " << mother.isHadron()
                      << ", tau: " << mother.tau() << std::endl;

            // Check if particle comes from hadron decay
            if (!(mother.isHadron() && mother.tau() > 0))
            {
                primaryTracks.push_back(trackParams);
            }

            if (abs(part.eta()) > 0.9)
                continue; // detector acceptance
            if (part.pT() < 0.15)
                continue; // pT cut

            allTracks.push_back(trackParams);
            truthTracks.push_back(Track(TVector3(part.xProd() * 0.1, part.yProd() * 0.1, part.zProd() * 0.1),
                                        TVector3(part.px(), part.py(), part.pz()),
                                        std::copysign(std::sqrt(part.xProd() * part.xProd() + part.yProd() * part.yProd()), part.yProd()) * 0.1,
                                        part.zProd() * 0.1,
                                        part.charge()));

            // Add to jet finding
            particles.push_back(fastjet::PseudoJet(trackParams.mom.X(),
                                                   trackParams.mom.Y(),
                                                   trackParams.mom.Z(),
                                                   trackParams.mom.Mag()));
            particles.back().set_user_index(allTracks.size() - 1);
        }

        // Find primary vertex
        std::cout << "Finding primary vertex..." << std::endl;
        if (primaryTracks.size() < 3)
        {
            std::cout << "Too few tracks for vertex finding: " << primaryTracks.size() << std::endl;
            continue;
        }

        TVector3 primaryVertex = VertexFinder::findPrimaryVertex(primaryTracks,
                                                                 0.001, // 10 micron convergence
                                                                 50,    // max iterations
                                                                 1.0,   // temperature
                                                                 6.0);  // tighter chi2 cut

        // Fill primary vertex Z position
        hPrimaryVertexZ->Fill(primaryVertex.Z());

        // Vertex position check
        if (std::abs(primaryVertex.Z()) > 10.0 ||
            std::abs(primaryVertex.X()) > 2.0 ||
            std::abs(primaryVertex.Y()) > 2.0)
        {
            std::cout << "WARNING: Vertex position outside range!" << std::endl;
            continue;
        }

        // Jet finding
        fastjet::ClusterSequence cs(particles, jetDef);
        std::vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(5.0)); // jets with pT > 20 GeV

        // Analyze each jet
        for (const auto &jet : jets)
        {
            hJetPt->Fill(jet.pt());
            hJetEta->Fill(jet.eta());

            // Get number of constituents
            std::vector<fastjet::PseudoJet> constituents = jet.constituents();
            int nConstituents = constituents.size();

            // Fill inclusive histogram
            hNConstituents->Fill(nConstituents);

            // Get jet flavor
            int jetFlavor = getJetFlavor(jet, mcparticles, R);

            if (jetFlavor == JetTaggingSpecies::beauty)
            {
                hNConstituents_b->Fill(nConstituents);
                hJetPt_b->Fill(jet.pt());
            }
            else if (jetFlavor == JetTaggingSpecies::charm)
            {
                hNConstituents_c->Fill(nConstituents);
                hJetPt_c->Fill(jet.pt());
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
            SecondaryVertexFinder svFinder;
            std::vector<SecondaryVertexFinder::SecVtxInfo> secondaryVertices =
                svFinder.findSecondaryVertices(jetTracks);

            SecondaryVertexFinder::fillMonitoringHistograms(secondaryVertices);
            // Fill SV multiplicity histograms
            hNSVperJet->Fill(secondaryVertices.size());

            if (jetFlavor == JetTaggingSpecies::beauty)
            {
                hNSVperJet_b->Fill(secondaryVertices.size());
            }
            else if (jetFlavor == JetTaggingSpecies::charm)
            {
                hNSVperJet_c->Fill(secondaryVertices.size());
            }

            // Process each secondary vertex
            for (const auto &sv : secondaryVertices)
            {
                hSVMass->Fill(sv.mass);
                hDecayLength2D->Fill(SVCalculations::calculateDecayLength2D(primaryVertex, sv.position));
                hDecayLength3D->Fill(SVCalculations::calculateDecayLength3D(primaryVertex, sv.position));

                // Fill flavor-specific histograms
                if (jetFlavor == JetTaggingSpecies::beauty)
                {
                    hSVMass_b->Fill(sv.mass);
                    hDecayLength2D_b->Fill(SVCalculations::calculateDecayLength2D(primaryVertex, sv.position));
                    hDecayLength3D_b->Fill(SVCalculations::calculateDecayLength3D(primaryVertex, sv.position));
                }
                else if (jetFlavor == JetTaggingSpecies::charm)
                {
                    hSVMass_c->Fill(sv.mass);
                    hDecayLength2D_c->Fill(SVCalculations::calculateDecayLength2D(primaryVertex, sv.position));
                    hDecayLength3D_c->Fill(SVCalculations::calculateDecayLength3D(primaryVertex, sv.position));
                }
            }

            // Fill jet features
            treeWriter.getData().mJetpT = jet.pt();
            treeWriter.getData().mJetEta = jet.eta();
            treeWriter.getData().mJetPhi = jet.phi();
            treeWriter.getData().mJetMass = jet.m();
            treeWriter.getData().mJetFlavor = jetFlavor;
            treeWriter.getData().mNTracks = jetTracks.size();
            treeWriter.getData().mNSV = secondaryVertices.size();

            // Fill track features
            for (int i = 0; i < treeWriter.getData().mNTracks; i++)
            {
                treeWriter.getData().mTrackpT[i] = jetTracks[i].mom.Pt();
                treeWriter.getData().mTrackEta[i] = jetTracks[i].mom.Eta();
                // treeWriter.getData().mSignedIP2D[i] = calculateSignedIP2D(jetTracks[i], primaryVertex);
                // treeWriter.getData().mSignedIP3D[i] = calculateSignedIP3D(jetTracks[i], primaryVertex);
            }

            // Fill SV features
            for (int i = 0; i < treeWriter.getData().mNSV; i++)
            {
                treeWriter.getData().mSVMass[i] = secondaryVertices[i].mass;
                treeWriter.getData().mDecayLength2D[i] = SVCalculations::calculateDecayLength2D(primaryVertex, secondaryVertices[i].position);
                treeWriter.getData().mDecayLength3D[i] = SVCalculations::calculateDecayLength3D(primaryVertex, secondaryVertices[i].position);
            }

            // Fill tree
            treeWriter.fill();
        }

        // After clustering jets, add to total
        jets = cs.inclusive_jets(5.0);
        totalJets += jets.size();
        totalEvents++;

        // Progress indicator
        if (iEvent % 100 == 0)
        {
            std::cout << "Processed " << iEvent << " events" << std::endl;
        }

        // Apply smearing to all tracks
        for (int i = 0; i < allTracks.size(); i++)
        {
            // Fill histograms
            dcaSmearer.fillHistograms(truthTracks[i], allTracks[i]);
        }

        std::cout << "Event " << iEvent << " processed successfully" << std::endl;
    }

    // Write and close
    outFile->cd();
    outFile->Write();
    outFile->Close();

    std::cout << "Analysis completed!" << std::endl;

    // Add safety check before calculating average
    if (totalEvents > 0)
    {
        double avgJets = static_cast<double>(totalJets) / totalEvents;
        std::cout << "Average number of jets with pT > 5 GeV per event: " << avgJets << std::endl;
    }
    else
    {
        std::cout << "No events processed!" << std::endl;
    }
}
