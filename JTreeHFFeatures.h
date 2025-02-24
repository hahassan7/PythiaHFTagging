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
        mSignedIP2D.fill(-999);
        mSignedIP2DSign.fill(-999);
        mSignedIP3D.fill(-999);
        mSignedIP3DSign.fill(-999);

        mIPxy.fill(-999);
        mCPA.fill(-999);
        mChi2PCA.fill(999);
        mDispersion.fill(999);
        mDecayLength2D.fill(-999);
        mDecayLength3D.fill(-999);
    }
};

enum class HFFeatures : uint8_t
{
    // Jet features
    JetpT = 0,
    JetEta,
    JetPhi,
    JetMass,
    NTracks,
    NSV,
    JetFlavor,

    // Track features
    TrackpT,
    TrackEta,
    DotProdTrackJet,
    DotProdTrackJetOverJet,
    DeltaRJetTrack,
    SignedIP2D,
    SignedIP2DSign,
    SignedIP3D,
    SignedIP3DSign,
    MomFraction,
    DeltaRTrackVertex,

    // SV features
    SVpT,
    DeltaRSVJet,
    SVMass,
    SVfE,
    IPxy,
    CPA,
    Chi2PCA,
    Dispersion,
    DecayLength2D,
    DecayLength2DError,
    DecayLength3D,
    DecayLength3DError
};

// Define the macro for filling the feature map
#define FILL_MAP_HF(FEATURE)                                 \
  {                                                          \
    #FEATURE, static_cast<uint8_t>(HFFeatures::FEATURE)      \
  }

// Macro for scalar branches
#define CREATE_BRANCH_HF_SCALAR(FEATURE, BRANCHNAME, TYPE)                                          \
    case static_cast<uint8_t>(HFFeatures::FEATURE):                                                 \
    {                                                                                               \
        tree->Branch(BRANCHNAME, &data.m##FEATURE, (std::string(BRANCHNAME) + "/" + TYPE).c_str()); \
        break;                                                                                      \
    }

// Macro for array branches
#define CREATE_BRANCH_HF_ARRAY(FEATURE, BRANCHNAME, TYPE)                                                                      \
    case static_cast<uint8_t>(HFFeatures::FEATURE):                                                                            \
    {                                                                                                                          \
        tree->Branch(BRANCHNAME, &data.m##FEATURE, (std::string(BRANCHNAME) + "[" + std::to_string(N) + "]/" + TYPE).c_str()); \
        break;                                                                                                                 \
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
            FILL_MAP_HF(JetpT),
            FILL_MAP_HF(JetEta),
            FILL_MAP_HF(JetPhi),
            FILL_MAP_HF(JetMass),
            FILL_MAP_HF(NTracks),
            FILL_MAP_HF(NSV),
            FILL_MAP_HF(JetFlavor),

            // Track features
            FILL_MAP_HF(TrackpT),
            FILL_MAP_HF(TrackEta),
            FILL_MAP_HF(DotProdTrackJet),
            FILL_MAP_HF(DotProdTrackJetOverJet),
            FILL_MAP_HF(DeltaRJetTrack),
            FILL_MAP_HF(SignedIP2D),
            FILL_MAP_HF(SignedIP2DSign),
            FILL_MAP_HF(SignedIP3D),
            FILL_MAP_HF(SignedIP3DSign),
            FILL_MAP_HF(MomFraction),
            FILL_MAP_HF(DeltaRTrackVertex),

            // SV features
            FILL_MAP_HF(SVpT),
            FILL_MAP_HF(DeltaRSVJet),
            FILL_MAP_HF(SVMass),
            FILL_MAP_HF(SVfE),
            FILL_MAP_HF(IPxy),
            FILL_MAP_HF(CPA),
            FILL_MAP_HF(Chi2PCA),
            FILL_MAP_HF(Dispersion),
            FILL_MAP_HF(DecayLength2D),
            FILL_MAP_HF(DecayLength2DError),
            FILL_MAP_HF(DecayLength3D),
            FILL_MAP_HF(DecayLength3DError)};
    }

    void createSelectedBranches()
    {
        for (const auto &featureIdx : mSelectedFeatures)
        {
            switch (featureIdx)
            {
                // Jet features
                CREATE_BRANCH_HF_SCALAR(JetpT, "jetPt", "F")
                CREATE_BRANCH_HF_SCALAR(JetEta, "jetEta", "F")
                CREATE_BRANCH_HF_SCALAR(JetPhi, "jetPhi", "F")
                CREATE_BRANCH_HF_SCALAR(JetMass, "jetMass", "F")
                CREATE_BRANCH_HF_SCALAR(NTracks, "nTracks", "I")
                CREATE_BRANCH_HF_SCALAR(NSV, "nSV", "I")
                CREATE_BRANCH_HF_SCALAR(JetFlavor, "jetFlavor", "I")

                // Track features
                CREATE_BRANCH_HF_ARRAY(TrackpT, "trackPt", "F")
                CREATE_BRANCH_HF_ARRAY(TrackEta, "trackEta", "F")
                CREATE_BRANCH_HF_ARRAY(DotProdTrackJet, "dotProdTrackJet", "F")
                CREATE_BRANCH_HF_ARRAY(DotProdTrackJetOverJet, "dotProdTrackJetOverJet", "F")
                CREATE_BRANCH_HF_ARRAY(DeltaRJetTrack, "deltaRJetTrack", "F")
                CREATE_BRANCH_HF_ARRAY(SignedIP2D, "signedIP2D", "F")
                CREATE_BRANCH_HF_ARRAY(SignedIP2DSign, "signedIP2DSign", "F")
                CREATE_BRANCH_HF_ARRAY(SignedIP3D, "signedIP3D", "F")
                CREATE_BRANCH_HF_ARRAY(SignedIP3DSign, "signedIP3DSign", "F")
                CREATE_BRANCH_HF_ARRAY(MomFraction, "momFraction", "F")
                CREATE_BRANCH_HF_ARRAY(DeltaRTrackVertex, "deltaRTrackVertex", "F")

                // SV features
                CREATE_BRANCH_HF_ARRAY(SVpT, "svPt", "F")
                CREATE_BRANCH_HF_ARRAY(DeltaRSVJet, "deltaRSVJet", "F")
                CREATE_BRANCH_HF_ARRAY(SVMass, "svMass", "F")
                CREATE_BRANCH_HF_ARRAY(SVfE, "svfE", "F")
                CREATE_BRANCH_HF_ARRAY(IPxy, "svIPxy", "F")
                CREATE_BRANCH_HF_ARRAY(CPA, "svCPA", "F")
                CREATE_BRANCH_HF_ARRAY(Chi2PCA, "svChi2", "F")
                CREATE_BRANCH_HF_ARRAY(Dispersion, "svDispersion", "F")
                CREATE_BRANCH_HF_ARRAY(DecayLength2D, "svDecayLength2D", "F")
                CREATE_BRANCH_HF_ARRAY(DecayLength2DError, "svDecayLength2DError", "F")
                CREATE_BRANCH_HF_ARRAY(DecayLength3D, "svDecayLength3D", "F")
                CREATE_BRANCH_HF_ARRAY(DecayLength3DError, "svDecayLength3DError", "F")
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