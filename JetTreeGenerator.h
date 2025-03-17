#ifndef JETTREEGENERATOR_H
#define JETTREEGENERATOR_H

#include "JTreeHFFeatures.h"
#include "TaggingUtilities.h"
#include "MLForHFTagging.h"
#include "TVector2.h"
#include <memory>
#include <algorithm>
#include <random> // Include for random number generation

template <std::size_t N = 10>
class JetTreeGenerator
{
private:
    std::unique_ptr<JTreeHFFeatures<N>> mJetTree;
    MLForHFTagging mFeatureMaker;
    std::unique_ptr<JetModel<N>> mMLModel;
    DebugLevel mDebugLevel;
    bool mDoPrediction = false;
    const double minTrack = 0.5; // Minimum track pT for track selection
    const double maxDCAxy = 1.0; // Maximum DCAxy for track selection
    const double maxDCAz = 2.0;  // Maximum DCAz for track selection

    bool passesCuts(const Track &track)
    {
        return track.pt() > minTrack && std::abs(track.getDCAxy()) < maxDCAxy && std::abs(track.getDCAz()) < maxDCAz;
    }

public:
    JetTreeGenerator(const char *filename, std::vector<std::string> features, DebugLevel debugLevel = DebugLevel::INFO)
        : mDebugLevel(debugLevel), mDoPrediction(false)
    {
        mJetTree = std::make_unique<JTreeHFFeatures<N>>(filename, features);

        // Seed the random number generator using current time
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
    }

    JetTreeGenerator(std::vector<std::string> features, const char *modelpath, DebugLevel debugLevel = DebugLevel::INFO)
        : mDebugLevel(debugLevel), mDoPrediction(true)
    {
        mFeatureMaker.cacheInputFeaturesIndices(features);
        mMLModel = std::make_unique<JetModel<N>>(modelpath);
    }

    void setDoPrediction(bool doPrediction) { mDoPrediction = doPrediction; }
    bool getDoPrediction() const { return mDoPrediction; }

    // Comparison function for sorting secondary vertices based on 2D decay length
    bool compareDecayLength2D(const SecondaryVertexFinder::SecVtxInfo &a,
                              const SecondaryVertexFinder::SecVtxInfo &b,
                              const TVector3 &pVertex)
    {
        return SVCalculations::calculateDecayLength2D(pVertex, a.position) > SVCalculations::calculateDecayLength2D(pVertex, b.position);
    }

    // Comparison function for sorting tracks based on IP2D
    bool compareIP2D(const Track &a, const Track &b, const fastjet::PseudoJet &jet, const TVector3 &pVertex)
    {
        double dcaTrack1 = getGeoSign(jet, a, pVertex) * std::abs(a.getDCAxy());
        double dcaTrack2 = getGeoSign(jet, b, pVertex) * std::abs(b.getDCAxy());
        return dcaTrack1 > dcaTrack2; // Sort in descending order
    }

    // Function to determine rejection based on pT
    bool rejectLightFlavorJet(float jetPt)
    {
        // Define rejection thresholds based on pT bins
        if (jetPt < 10.0)
        {
            return (static_cast<double>(std::rand()) / RAND_MAX) < 0.85;
        }
        else if (jetPt < 20.0)
        {
            return (static_cast<double>(std::rand()) / RAND_MAX) < 0.82;
        }
        else if (jetPt < 100.0)
        {
            return (static_cast<double>(std::rand()) / RAND_MAX) < 0.80;
        }
        else if (jetPt < 200.0)
        {
            return (static_cast<double>(std::rand()) / RAND_MAX) < 0.78;
        }
        else
        {
            return (static_cast<double>(std::rand()) / RAND_MAX) < 0.70;
        }
    }

    void processJet(const fastjet::PseudoJet &jet,
                    const std::vector<Track> &allTracks,
                    const std::vector<JParticle> &mcparticles,
                    std::vector<SecondaryVertexFinder::SecVtxInfo> &secondaryVertices,
                    const TVector3 &primaryVertex,
                    double R)
    {
        auto &treeData = mJetTree->getData();

        // Fill jet features
        treeData.mJetpT = jet.pt();
        treeData.mJetEta = jet.eta();
        treeData.mJetPhi = jet.phi();
        treeData.mJetMass = jet.m();
        treeData.mJetFlavor = getJetFlavor(jet, mcparticles, R, mDebugLevel);

        // Reject light flavor jets based on pT
        if (treeData.mJetFlavor == 0 && rejectLightFlavorJet(treeData.mJetpT))
        {
            mJetTree->clear();
            return; // Skip processing for this jet if rejected
        }

        // Process constituents
        std::vector<fastjet::PseudoJet> constituents = jet.constituents();
        treeData.mNTracks = constituents.size();
        treeData.mNSV = secondaryVertices.size();

        // Process tracks
        std::vector<Track> sortedTracks; // Vector to hold sorted tracks
        for (size_t iTrack = 0; iTrack < constituents.size(); ++iTrack)
        {
            size_t trackIndex = constituents[iTrack].user_index();
            if (trackIndex < allTracks.size())
            {
                if (passesCuts(allTracks[trackIndex]))
                {
                    sortedTracks.push_back(allTracks[trackIndex]); // Add track to the sorted list
                }
            }
        }

        // Sort secondary vertices based on 2D decay length (descending order)
        std::sort(sortedTracks.begin(), sortedTracks.end(),
                  [&](const Track &a, const Track &b)
                  {
                      return compareIP2D(a, b, jet, primaryVertex);
                  });

        // Process tracks
        TVector3 jetVec(jet.px(), jet.py(), jet.pz());
        for (size_t iTrack = 0; iTrack < std::min(sortedTracks.size(), size_t(N)); ++iTrack)
        {
            const auto &track = sortedTracks[iTrack];

            // Basic track kinematics
            treeData.mTrackpT[iTrack] = track.pt();
            treeData.mTrackEta[iTrack] = track.eta();

            // Track-jet correlation features
            treeData.mDotProdTrackJet[iTrack] = track.getMomentum().Dot(jetVec);
            treeData.mDotProdTrackJetOverJet[iTrack] = track.getMomentum().Dot(jetVec) / jetVec.Mag();
            treeData.mDeltaRJetTrack[iTrack] = deltaR(track, jetVec);

            // Impact parameter features
            int sign = getGeoSign(jet, track, primaryVertex);
            treeData.mSignedIP2D[iTrack] = sign * std::abs(track.getDCAxy());
            treeData.mSignedIP3D[iTrack] = sign * std::abs(track.getDCAz());

            // Momentum fraction
            treeData.mMomFraction[iTrack] = track.pt() / jet.pt();

            // Track to primary vertex distance
            double minDistance = std::numeric_limits<double>::max();
            for (auto &vertex : secondaryVertices)
            {
                double distance = deltaR(track, vertex.momentum);
                if (distance < minDistance)
                {
                    minDistance = distance;
                }
            }
            treeData.mDeltaRTrackVertex[iTrack] = minDistance;
        }

        // Sort secondary vertices based on 2D decay length (descending order)
        std::sort(secondaryVertices.begin(), secondaryVertices.end(),
                  [&](const SecondaryVertexFinder::SecVtxInfo &a, const SecondaryVertexFinder::SecVtxInfo &b)
                  {
                      return compareDecayLength2D(a, b, primaryVertex);
                  });

        // Process secondary vertices
        for (size_t iSV = 0; iSV < std::min(secondaryVertices.size(), size_t(N)); ++iSV)
        {
            const auto &sv = secondaryVertices[iSV];
            TVector3 svMomentum(sv.momentum.Px(), sv.momentum.Py(), sv.momentum.Pz());
            TVector3 flightLine = sv.position - primaryVertex;

            // Basic SV kinematics
            treeData.mSVpT[iSV] = sv.momentum.Pt();
            treeData.mSVMass[iSV] = sv.momentum.M();

            // SV-jet correlation
            treeData.mDeltaRSVJet[iSV] = deltaR(jet, sv.momentum);

            // Energy fraction
            treeData.mSVfE[iSV] = sv.momentum.E() / jet.e();

            // Impact parameter
            treeData.mIPxy[iSV] = SVCalculations::calculateIPxy(primaryVertex, sv.position, sv.momentum.Vect());

            // Pointing angle
            treeData.mCPA[iSV] = SVCalculations::calculateCPA(primaryVertex, sv.position, sv.momentum.Vect());

            // Chi2 of PCA
            treeData.mChi2PCA[iSV] = sv.chi2;

            // Vertex dispersion
            treeData.mDispersion[iSV] = SVCalculations::calculateDispersion(sv, allTracks);

            // Decay lengths
            treeData.mDecayLength2D[iSV] = SVCalculations::calculateDecayLength2D(primaryVertex, sv.position);
            treeData.mDecayLength3D[iSV] = SVCalculations::calculateDecayLength3D(primaryVertex, sv.position);
        }

        mJetTree->fill();
    }

    std::tuple<BJetParams, std::vector<BJetTrackParams>, std::vector<BJetSVParams>> processJetFeatures(const fastjet::PseudoJet &jet,
                                                                                                       const std::vector<Track> &allTracks,
                                                                                                       std::vector<SecondaryVertexFinder::SecVtxInfo> &secondaryVertices,
                                                                                                       const TVector3 &primaryVertex)
    {

        BJetParams jetparams;
        std::vector<BJetTrackParams> trackparams;
        std::vector<BJetSVParams> svparams;

        // Fill jet features
        jetparams.jetpT = jet.pt();
        jetparams.jetEta = jet.eta();
        jetparams.jetPhi = jet.phi();
        jetparams.jetMass = jet.m();

        // Process constituents
        std::vector<fastjet::PseudoJet> constituents = jet.constituents();
        jetparams.nTracks = constituents.size();
        jetparams.nSV = secondaryVertices.size();

        // Process tracks
        std::vector<Track> sortedTracks; // Vector to hold sorted tracks
        for (size_t iTrack = 0; iTrack < constituents.size(); ++iTrack)
        {
            size_t trackIndex = constituents[iTrack].user_index();
            if (trackIndex < allTracks.size())
            {
                if (passesCuts(allTracks[trackIndex]))
                {
                    sortedTracks.push_back(allTracks[trackIndex]); // Add track to the sorted list
                }
            }
        }

        // Sort secondary vertices based on 2D decay length (descending order)
        std::sort(sortedTracks.begin(), sortedTracks.end(),
                  [&](const Track &a, const Track &b)
                  {
                      return compareIP2D(a, b, jet, primaryVertex);
                  });

        // Process tracks
        TVector3 jetVec(jet.px(), jet.py(), jet.pz());
        for (size_t iTrack = 0; iTrack < std::min(sortedTracks.size(), size_t(N)); ++iTrack)
        {
            const auto &track = sortedTracks[iTrack];

            // Basic track kinematics
            BJetTrackParams trackparam;
            trackparam.trackpT = track.pt();
            trackparam.trackEta = track.eta();

            // Track-jet correlation features
            trackparam.dotProdTrackJet = track.getMomentum().Dot(jetVec);
            trackparam.dotProdTrackJetOverJet = track.getMomentum().Dot(jetVec) / jetVec.Mag();
            trackparam.deltaRJetTrack = deltaR(track, jetVec);

            // Impact parameter features
            int sign = getGeoSign(jet, track, primaryVertex);
            trackparam.signedIP2D = sign * std::abs(track.getDCAxy());
            trackparam.signedIP3D = sign * std::abs(track.getDCAz());

            // Momentum fraction
            trackparam.momFraction = track.pt() / jet.pt();

            // Track to primary vertex distance
            double minDistance = std::numeric_limits<double>::max();
            for (auto &vertex : secondaryVertices)
            {
                double distance = deltaR(track, vertex.momentum);
                if (distance < minDistance)
                {
                    minDistance = distance;
                }
            }
            trackparam.deltaRTrackVertex = minDistance;
            trackparams.push_back(trackparam);
        }

        // Sort secondary vertices based on 2D decay length (descending order)
        std::sort(secondaryVertices.begin(), secondaryVertices.end(),
                  [&](const SecondaryVertexFinder::SecVtxInfo &a, const SecondaryVertexFinder::SecVtxInfo &b)
                  {
                      return compareDecayLength2D(a, b, primaryVertex);
                  });

        // Process secondary vertices
        for (size_t iSV = 0; iSV < std::min(secondaryVertices.size(), size_t(N)); ++iSV)
        {
            BJetSVParams svparam;
            const auto &sv = secondaryVertices[iSV];
            TVector3 svMomentum(sv.momentum.Px(), sv.momentum.Py(), sv.momentum.Pz());
            TVector3 flightLine = sv.position - primaryVertex;

            // Basic SV kinematics
            svparam.svPt = sv.momentum.Pt();
            svparam.svMass = sv.momentum.M();

            // SV-jet correlation
            svparam.deltaRSVJet = deltaR(jet, sv.momentum);

            // Energy fraction
            svparam.svfE = sv.momentum.E() / jet.e();

            // Impact parameter
            svparam.svIPxy = SVCalculations::calculateIPxy(primaryVertex, sv.position, sv.momentum.Vect());

            // Pointing angle
            svparam.svCPA = SVCalculations::calculateCPA(primaryVertex, sv.position, sv.momentum.Vect());

            // Chi2 of PCA
            svparam.svChi2PCA = sv.chi2;

            // Vertex dispersion
            svparam.svDispersion = SVCalculations::calculateDispersion(sv, allTracks);

            // Decay lengths
            svparam.svDecayLength2D = SVCalculations::calculateDecayLength2D(primaryVertex, sv.position);
            svparam.svDecayLength3D = SVCalculations::calculateDecayLength3D(primaryVertex, sv.position);
            svparams.push_back(svparam);
        }

        trackparams.resize(N);
        svparams.resize(N);
        return std::make_tuple(jetparams, trackparams, svparams);
    }

    std::vector<float> predict(const BJetParams &jetparams, const std::vector<BJetTrackParams> &trackparams, const std::vector<BJetSVParams> &svparams)
    {
        auto inputFeatures = mFeatureMaker.getInputFeatures2D(jetparams, trackparams, svparams);
        return mMLModel->predict(inputFeatures[0], inputFeatures[1], inputFeatures[2]);
    }

    void write()
    {
        mJetTree->write();
    }
};

#endif // JETTREEGENERATOR_H