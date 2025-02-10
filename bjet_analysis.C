#include "TaggingUtilities.h"
#include "JTreeHFFeatures.h"

void bjet_analysis()
{
    // Initialize features we want to store
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

    // Initialize Pythia
    Pythia8::Pythia pythia;
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("Beams:allowVertexSpread = on");
    pythia.readString("Beams:sigmaVertexX = 10");   // Back to mm (= 1 cm)
    pythia.readString("Beams:sigmaVertexY = 10");   // Back to mm (= 1 cm)
    pythia.readString("Beams:sigmaVertexZ = 50");   // Back to mm (= 5 cm)
    pythia.readString("HardQCD:all = on");          // Turn on hard QCD processes
    pythia.readString("PhaseSpace:pTHatMin = 20."); // Min pT for hard interaction

    pythia.init();

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
    DCASmearing::initializeHistograms();

    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        if (!pythia.next())
            continue;

        // Collect tracks for primary vertex finding
        std::vector<TrackWithWeight> primaryTracks;

        std::vector<fastjet::PseudoJet> mcparticles;

        // Loop over all tracks in the event
        for (int i = 0; i < pythia.event.size(); ++i)
        {
            const Pythia8::Particle &part = pythia.event[i];

            // Parton Definition
            if ((part.id() == 5 || part.id() == -5))
            {
                mcparticles.push_back(fastjet::PseudoJet(part.px(), part.py(), part.pz(), part.e()));
            }
            else if ((part.id() == 4 || part.id() == -4))
            {
                mcparticles.push_back(fastjet::PseudoJet(part.px(), part.py(), part.pz(), part.e()));
            }

            // Select only charged final state particles
            if (!part.isFinal() || part.charge() == 0)
                continue;

            // Get production vertex coordinates
            double xProd = part.xProd() * 0.1; // Convert mm to cm
            double yProd = part.yProd() * 0.1;
            double zProd = part.zProd() * 0.1;

            // Add detector resolution smearing
            double xSmeared = random.Gaus(xProd, spatialResolution);
            double ySmeared = random.Gaus(yProd, spatialResolution);
            double zSmeared = random.Gaus(zProd, spatialResolution);

            // Check if particle comes from hadron decay
            bool isFromHadronDecay = false;
            const Pythia8::Particle &mother = pythia.event[part.mother1()];
            if (mother.isHadron() && mother.tau() > 1e-12)
            {
                isFromHadronDecay = true;
            }
            if (isFromHadronDecay)
                continue;

            // Create track parameters with smeared positions
            TVector3 pos(xSmeared, ySmeared, zSmeared);
            TVector3 truePos(xProd, yProd, zProd);
            TVector3 mom(part.px(), part.py(), part.pz());

            // Add to collection with appropriate uncertainty
            primaryTracks.push_back(VertexFinder::createTrackWithWeight(pos, mom, spatialResolution));
        }

        // Add check for minimum number of tracks
        if (primaryTracks.size() < 3)
        {
            std::cout << "Warning: Too few tracks for vertex finding" << std::endl;
            continue;
        }

        // Find primary vertex
        TVector3 primaryVertex = VertexFinder::findPrimaryVertex(primaryTracks,
                                                                 0.001, // 10 micron convergence
                                                                 50,    // max iterations
                                                                 1.0,   // temperature
                                                                 6.0);  // tighter chi2 cut

        // Add sanity check on found vertex position
        if (std::abs(primaryVertex.Z()) > 10.0 ||
            std::abs(primaryVertex.X()) > 2.0 ||
            std::abs(primaryVertex.Y()) > 2.0)
        {
            std::cout << "Warning: Found vertex position outside expected range" << std::endl;
            continue;
        }

        // Collect final state particles and apply detector simulation
        std::vector<fastjet::PseudoJet> particles;
        std::vector<DetectorSimulation::SmearTrackParams> allTracks;

        for (int i = 0; i < pythia.event.size(); ++i)
        {
            const Pythia8::Particle &part = pythia.event[i];

            // Select only final state charged particles
            if (!part.isFinal() || part.charge() == 0)
                continue;
            if (abs(part.eta()) > 0.9)
                continue; // detector acceptance
            if (part.pT() < 0.15)
                continue; // pT cut

            // Apply detector simulation
            auto trackParams = detector.smearTrack(part, primaryVertex);
            allTracks.push_back(trackParams);

            // Add to jet finding
            particles.push_back(fastjet::PseudoJet(trackParams.mom.X(),
                                                   trackParams.mom.Y(),
                                                   trackParams.mom.Z(),
                                                   trackParams.mom.Mag()));
            particles.back().set_user_index(allTracks.size() - 1);
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
            std::vector<DetectorSimulation::SmearTrackParams> jetTracks;
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
                treeWriter.getData().mSignedIP2D[i] = calculateSignedIP2D(jetTracks[i], primaryVertex);
                treeWriter.getData().mSignedIP3D[i] = calculateSignedIP3D(jetTracks[i], primaryVertex);
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
        for (const auto &track : allTracks)
        {
            // Create truth DCA parameters
            DCASmearing::DCAParams truthDCA;
            // Calculate signed DCA values (can be positive or negative)
            truthDCA.dcaXY = std::copysign(std::sqrt(track.pos.X() * track.pos.X() + track.pos.Y() * track.pos.Y()), track.pos.Y());
            truthDCA.dcaZ = track.pos.Z();

            // Apply smearing
            DCASmearing::DCAParams smearedDCA = DCASmearing::smearDCA(track.mom, truthDCA, DCASmearing::SmearingMethod::PT_DEPENDENT);

            // Fill histograms
            DCASmearing::fillHistograms(track.mom, truthDCA, smearedDCA);
        }
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
