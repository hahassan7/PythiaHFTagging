#ifndef JTREEHFFEATURES_H
#define JTREEHFFEATURES_H

#include "TFile.h"
#include "TTree.h"
#include <map>
#include <string>
#include <vector>

template <std::size_t N = 10>
struct HFjetTree
{
    float mJetpT = 0.0;
    float mJetEta = -99;
    float mJetPhi = -99;
    int mNTracks = -1;
    int mNSV = -1;
    float mJetMass = 0.0;
    int mJetFlavor = 0;

    std::array<float, N> mTrackpT;
    std::array<float, N> mTrackEta;
    std::array<float, N> mDotProdTrackJet;
    std::array<float, N> mDotProdTrackJetOverJet;
    std::array<float, N> mDeltaRJetTrack;
    std::array<float, N> mSignedIP2D;
    std::array<float, N> mSignedIP2DSign;
    std::array<float, N> mSignedIP3D;
    std::array<float, N> mSignedIP3DSign;
    std::array<float, N> mMomFraction;
    std::array<float, N> mDeltaRTrackVertex;

    std::array<float, N> mSVpT;
    std::array<float, N> mDeltaRSVJet;
    std::array<float, N> mSVMass;
    std::array<float, N> mSVfE;
    std::array<float, N> mIPxy;
    std::array<float, N> mCPA;
    std::array<float, N> mChi2PCA;
    std::array<float, N> mDispersion;
    std::array<float, N> mDecayLength2D;
    std::array<float, N> mDecayLength2DError;
    std::array<float, N> mDecayLength3D;
    std::array<float, N> mDecayLength3DError;

    // Constructor to initialize arrays with -999
    HFjetTree()
    {
        mJetpT = 0.0;
        mJetEta = -99;
        mJetPhi = -99;
        mNTracks = -1;
        mNSV = -1;
        mJetMass = 0.0;
        mJetFlavor = 0;

        mTrackpT.fill(0);
        mTrackEta.fill(-99);
        mDotProdTrackJet.fill(0);
        mDotProdTrackJetOverJet.fill(0);
        mDeltaRJetTrack.fill(-1);
        mSignedIP2D.fill(-999);
        mSignedIP2DSign.fill(-999);
        mSignedIP3D.fill(-999);
        mSignedIP3DSign.fill(-999);
        mMomFraction.fill(0);
        mDeltaRTrackVertex.fill(-1);

        mSVpT.fill(0);
        mDeltaRSVJet.fill(-1);
        mSVMass.fill(0);
        mSVfE.fill(0);
        mIPxy.fill(-999);
        mCPA.fill(-999);
        mChi2PCA.fill(999);
        mDispersion.fill(999);
        mDecayLength2D.fill(-999);
        mDecayLength2DError.fill(-999);
        mDecayLength3D.fill(-999);
        mDecayLength3DError.fill(-999);
    }
};

enum class HFFeatures : uint8_t
{
    // Jet features
    jetpT = 0,
    jetEta,
    jetPhi,
    nTracks,
    nSV,
    jetMass,
    jetFlavor,

    // Track features
    trackpT,
    trackEta,
    dotProdTrackJet,
    dotProdTrackJetOverJet,
    deltaRJetTrack,
    signedIP2D,
    signedIP2DSign,
    signedIP3D,
    signedIP3DSign,
    momFraction,
    deltaRTrackVertex,
    trackPhi,
    trackCharge,
    trackITSChi2NCl,
    trackTPCChi2NCl,
    trackITSNCls,
    trackTPCNCls,
    trackTPCNCrossedRows,
    trackOrigin,
    trackVtxIndex,

    // Secondary vertex features
    svPt,
    deltaRSVJet,
    svMass,
    svfE,
    svIPxy,
    svCPA,
    svChi2PCA,
    svDispersion,
    svDecayLength2D,
    svDecayLength2DError,
    svDecayLength3D,
    svDecayLength3DError,
};

// Define the macro for filling the feature map
#define FILL_MAP_HF(FEATURE)                                 \
  {                                                          \
    #FEATURE, static_cast<uint8_t>(HFFeatures::FEATURE)      \
  }

// Macro for scalar branches
#define CREATE_BRANCH_HF_SCALAR(FEATURE, FEATURE1, BRANCHNAME, TYPE)                                 \
    case static_cast<uint8_t>(HFFeatures::FEATURE):                                                  \
    {                                                                                                \
        tree->Branch(BRANCHNAME, &data.m##FEATURE1, (std::string(BRANCHNAME) + "/" + TYPE).c_str()); \
        break;                                                                                       \
    }

// Macro for array branches
#define CREATE_BRANCH_HF_ARRAY(FEATURE, FEATURE1, BRANCHNAME, TYPE)                                                             \
    case static_cast<uint8_t>(HFFeatures::FEATURE):                                                                             \
    {                                                                                                                           \
        tree->Branch(BRANCHNAME, &data.m##FEATURE1, (std::string(BRANCHNAME) + "[" + std::to_string(N) + "]/" + TYPE).c_str()); \
        break;                                                                                                                  \
    }

template <std::size_t N = 10>
class JTreeHFFeatures
{
private:
    TFile *outFile;
    TTree *tree;
    HFjetTree<N> data; // Using the existing struct
    std::map<std::string, uint8_t> mAvailableFeatures;
    std::vector<uint8_t> mSelectedFeatures;

    void initializeFeatureMap()
    {
        mAvailableFeatures = {
            // Jet features
            FILL_MAP_HF(jetpT),
            FILL_MAP_HF(jetEta),
            FILL_MAP_HF(jetPhi),
            FILL_MAP_HF(jetMass),
            FILL_MAP_HF(nTracks),
            FILL_MAP_HF(nSV),
            FILL_MAP_HF(jetFlavor),

            // Track features
            FILL_MAP_HF(trackpT),
            FILL_MAP_HF(trackEta),
            FILL_MAP_HF(dotProdTrackJet),
            FILL_MAP_HF(dotProdTrackJetOverJet),
            FILL_MAP_HF(deltaRJetTrack),
            FILL_MAP_HF(signedIP2D),
            FILL_MAP_HF(signedIP2DSign),
            FILL_MAP_HF(signedIP3D),
            FILL_MAP_HF(signedIP3DSign),
            FILL_MAP_HF(momFraction),
            FILL_MAP_HF(deltaRTrackVertex),

            // SV features
            FILL_MAP_HF(svPt),
            FILL_MAP_HF(deltaRSVJet),
            FILL_MAP_HF(svMass),
            FILL_MAP_HF(svfE),
            FILL_MAP_HF(svIPxy),
            FILL_MAP_HF(svCPA),
            FILL_MAP_HF(svChi2PCA),
            FILL_MAP_HF(svDispersion),
            FILL_MAP_HF(svDecayLength2D),
            FILL_MAP_HF(svDecayLength2DError),
            FILL_MAP_HF(svDecayLength3D),
            FILL_MAP_HF(svDecayLength3DError)};
    }

    void createSelectedBranches()
    {
        for (const auto &featureIdx : mSelectedFeatures)
        {
            switch (featureIdx)
            {
                // Jet features
                CREATE_BRANCH_HF_SCALAR(jetpT, JetpT, "jetPt", "F")
                CREATE_BRANCH_HF_SCALAR(jetEta, JetEta, "jetEta", "F")
                CREATE_BRANCH_HF_SCALAR(jetPhi, JetPhi, "jetPhi", "F")
                CREATE_BRANCH_HF_SCALAR(jetMass, JetMass, "jetMass", "F")
                CREATE_BRANCH_HF_SCALAR(nTracks, NTracks, "nTracks", "I")
                CREATE_BRANCH_HF_SCALAR(nSV, NSV, "nSV", "I")
                CREATE_BRANCH_HF_SCALAR(jetFlavor, JetFlavor, "jetFlavor", "I")

                // Track features
                CREATE_BRANCH_HF_ARRAY(trackpT, TrackpT, "trackPt", "F")
                CREATE_BRANCH_HF_ARRAY(trackEta, TrackEta, "trackEta", "F")
                CREATE_BRANCH_HF_ARRAY(dotProdTrackJet, DotProdTrackJet, "dotProdTrackJet", "F")
                CREATE_BRANCH_HF_ARRAY(dotProdTrackJetOverJet, DotProdTrackJetOverJet, "dotProdTrackJetOverJet", "F")
                CREATE_BRANCH_HF_ARRAY(deltaRJetTrack, DeltaRJetTrack, "deltaRJetTrack", "F")
                CREATE_BRANCH_HF_ARRAY(signedIP2D, SignedIP2D, "signedIP2D", "F")
                CREATE_BRANCH_HF_ARRAY(signedIP2DSign, SignedIP2DSign, "signedIP2DSign", "F")
                CREATE_BRANCH_HF_ARRAY(signedIP3D, SignedIP3D, "signedIP3D", "F")
                CREATE_BRANCH_HF_ARRAY(signedIP3DSign, SignedIP3DSign, "signedIP3DSign", "F")
                CREATE_BRANCH_HF_ARRAY(momFraction, MomFraction, "momFraction", "F")
                CREATE_BRANCH_HF_ARRAY(deltaRTrackVertex, DeltaRTrackVertex, "deltaRTrackVertex", "F")

                // SV features
                CREATE_BRANCH_HF_ARRAY(svPt, SVpT, "svPt", "F")
                CREATE_BRANCH_HF_ARRAY(deltaRSVJet, DeltaRSVJet, "deltaRSVJet", "F")
                CREATE_BRANCH_HF_ARRAY(svMass, SVMass, "svMass", "F")
                CREATE_BRANCH_HF_ARRAY(svfE, SVfE, "svfE", "F")
                CREATE_BRANCH_HF_ARRAY(svIPxy, IPxy, "svIPxy", "F")
                CREATE_BRANCH_HF_ARRAY(svCPA, CPA, "svCPA", "F")
                CREATE_BRANCH_HF_ARRAY(svChi2PCA, Chi2PCA, "svChi2PCA", "F")
                CREATE_BRANCH_HF_ARRAY(svDispersion, Dispersion, "svDispersion", "F")
                CREATE_BRANCH_HF_ARRAY(svDecayLength2D, DecayLength2D, "svDecayLength2D", "F")
                CREATE_BRANCH_HF_ARRAY(svDecayLength2DError, DecayLength2DError, "svDecayLength2DError", "F")
                CREATE_BRANCH_HF_ARRAY(svDecayLength3D, DecayLength3D, "svDecayLength3D", "F")
                CREATE_BRANCH_HF_ARRAY(svDecayLength3DError, DecayLength3DError, "svDecayLength3DError", "F")
            }
        }
    }

public:
    JTreeHFFeatures(const char *filename, const std::vector<std::string> &features)
    {
        outFile = new TFile(filename, "RECREATE");
        tree = new TTree("jetTree", "Jet features tree");

        initializeFeatureMap();

        // Convert feature names to indices
        for (const auto &feature : features)
        {
            auto it = mAvailableFeatures.find(feature);
            if (it != mAvailableFeatures.end())
            {
                mSelectedFeatures.push_back(it->second);
            } else {
                throw std::runtime_error("Feature " + feature + " not found in the available features list");
            }
        }

        createSelectedBranches();
    }

    void fill()
    {
        tree->Fill();
        // Reset data structure for next jet
        data = HFjetTree<N>();
    }

    void clear()
    {
        data = HFjetTree<N>();
    }

    void write()
    {
        outFile->cd();
        tree->Write();
        outFile->Close();
    }

    ~JTreeHFFeatures()
    {
        if (outFile)
        {
            delete outFile;
        }
    }

    // Add these methods to access the data structure
    HFjetTree<N> &getData() { return data; }
    const HFjetTree<N> &getData() const { return data; }
};

#undef FILL_MAP_HF
#undef CREATE_BRANCH_HF_SCALAR
#undef CREATE_BRANCH_HF_ARRAY

#endif // JTREEHFFEATURES_H