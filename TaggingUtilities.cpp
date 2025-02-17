#include "TaggingUtilities.h"

// Static member initialization
TH1F *DCASmearing::hDCAxy_Truth = nullptr;
TH1F *DCASmearing::hDCAxy_Smeared = nullptr;
TH1F *DCASmearing::hDCAxy_Diff = nullptr;
TH1F *DCASmearing::hDCAz_Truth = nullptr;
TH1F *DCASmearing::hDCAz_Smeared = nullptr;
TH1F *DCASmearing::hDCAz_Diff = nullptr;
TH2F *DCASmearing::hDCAxy_vs_pT = nullptr;
TH2F *DCASmearing::hDCAz_vs_pT = nullptr;
TProfile *DCASmearing::hResolutionXY_vs_pT = nullptr;
TProfile *DCASmearing::hResolutionZ_vs_pT = nullptr;
TH2F *DCASmearing::hDCAxy_vs_eta = nullptr;
TH2F *DCASmearing::hDCAz_vs_eta = nullptr;

// Implementation of isBMeson
bool isBMeson(int pc)
{
    static const int bPdG[] = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523,
                               515, 525, 531, 10531, 533, 10533, 20533, 535, 541, 10541, 543,
                               10543, 20543, 545, 551, 10551, 100551, 110551, 200551, 210551,
                               553, 10553, 20553, 30553, 100553, 110553, 120553, 130553, 200553,
                               210553, 220553, 300553, 9000533, 9010553, 555, 10555, 20555, 100555,
                               110555, 120555, 200555, 557, 100557};
    for (int i = 0; i < (int)(sizeof(bPdG) / sizeof(int)); ++i)
        if (std::abs(pc) == bPdG[i])
            return true;
    return false;
}

// Implementation of isDMeson
bool isDMeson(int pc)
{
    static const int dPdG[] = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20431,
                               20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, 441, 10441,
                               100441, 443, 10443, 20443, 100443, 30443, 9000443, 9010443, 9020443,
                               445, 100445};
    for (int i = 0; i < (int)(sizeof(dPdG) / sizeof(int)); ++i)
        if (std::abs(pc) == dPdG[i])
            return true;
    return false;
}

// Function to smear pT based on ALICE ITS resolution
double DetectorSimulation::smear_pT(double pT)
{
    // Compute momentum resolution
    double sigma_pT_over_pT = std::sqrt(kPtResA * kPtResA + (kPtResB * kPtResB) / (pT * pT) + (kPtResC * pT) * (kPtResC * pT));

    double smeared_pT = pT * (1.0 + random.Gaus(0.0, sigma_pT_over_pT));

    return smeared_pT;
}

void DetectorSimulation::smearDCA(const Pythia8::Particle &parttobesmeared, Track &smearedTrack, const TVector3 &primaryVertex, SmearingMethod method)
{
    double xProd = parttobesmeared.xProd() * 0.1; // Convert to cm
    double yProd = parttobesmeared.yProd() * 0.1;
    double zProd = parttobesmeared.zProd() * 0.1;
    double pT = parttobesmeared.pT();
    double eta = parttobesmeared.eta();

    double etaFactor = 1.0 + 0.1 * std::abs(eta);

    double trueXYdca = std::copysign(std::sqrt(std::pow(xProd - primaryVertex.X(), 2) + std::pow(yProd - primaryVertex.Y(), 2)), yProd);
    double trueZdca = zProd - primaryVertex.Z();

    double sigmaXY, sigmaZ;

    if (pT < 0.1)
        pT = 0.1;

    switch (method)
    {
    case SmearingMethod::PT_DEPENDENT:
    {
        // Calculate resolutions
        sigmaXY = std::sqrt(std::pow(kDCAxyOffset, 2) + std::pow(kDCAxySlope, 2) / (pT * pT));
        sigmaZ = std::sqrt(std::pow(kDCAzOffset, 2) + std::pow(kDCAzSlope, 2) / (pT * pT));

        sigmaXY *= etaFactor;
        sigmaZ *= etaFactor;

        // Smear DCA
        smearedTrack.setDCAxy(random.Gaus(trueXYdca, sigmaXY));
        smearedTrack.setDCAz(random.Gaus(trueZdca, sigmaZ));

        smearedTrack.setPosition(smearPrimaryVertex(TVector3(xProd, yProd, zProd)));
        break;
    }

    case SmearingMethod::GAUSSIAN:
    {
        // Fixed resolution based on average ITS2 performance
        sigmaXY = 0.0010; // 10 μm -> cm
        sigmaZ = 0.0010;  // 10 μm -> cm

        smearedTrack.setDCAxy(random.Gaus(trueXYdca, sigmaXY));
        smearedTrack.setDCAz(random.Gaus(trueZdca, sigmaZ));

        smearedTrack.setPosition(smearPrimaryVertex(TVector3(xProd, yProd, zProd)));
        break;
    }

    default:
        smearedTrack.setDCAxy(trueXYdca);
        smearedTrack.setDCAz(trueZdca);
        smearedTrack.setPosition(TVector3(xProd, yProd, zProd));
        break;
    }
}

void DCASmearing::initializeHistograms()
{
    // DCA XY (finer binning near zero)
    hDCAxy_Truth = new TH1F("hDCAxy_Truth", "Truth DCA_{xy};DCA_{xy} (cm);Counts",
                            200, -1.0, 1.0); // ±1 mm range
    hDCAxy_Smeared = new TH1F("hDCAxy_Smeared", "Smeared DCA_{xy};DCA_{xy} (cm);Counts",
                              200, -1.0, 1.0);
    hDCAxy_Diff = new TH1F("hDCAxy_Diff", "Smeared - Truth DCA_{xy};#Delta DCA_{xy} (cm);Counts",
                           200, -0.01, 0.01); // ±100 μm range

    // DCA Z
    hDCAz_Truth = new TH1F("hDCAz_Truth", "Truth DCA_{z};DCA_{z} (cm);Counts",
                           200, -1.0, 1.0);
    hDCAz_Smeared = new TH1F("hDCAz_Smeared", "Smeared DCA_{z};DCA_{z} (cm);Counts",
                             200, -1.0, 1.0);
    hDCAz_Diff = new TH1F("hDCAz_Diff", "Smeared - Truth DCA_{z};#Delta DCA_{z} (cm);Counts",
                          200, -0.01, 0.01);

    // pT-dependent resolution (log scale in pT)
    hDCAxy_vs_pT = new TH2F("hDCAxy_vs_pT", "DCA_{xy} vs p_{T};p_{T} (GeV/c);DCA_{xy} (cm)",
                            100, 0.1, 20, 200, -0.1, 0.1);
    hDCAz_vs_pT = new TH2F("hDCAz_vs_pT", "DCA_{z} vs p_{T};p_{T} (GeV/c);DCA_{z} (cm)",
                           100, 0.1, 20, 200, -0.1, 0.1);

    // Fill resolution vs pT plots
    hResolutionXY_vs_pT = new TProfile("hResolutionXY_vs_pT",
                                       "DCA_{xy} Resolution vs p_{T};p_{T} (GeV/c);#sigma_{xy} (cm)",
                                       100, 0.1, 20);
    hResolutionZ_vs_pT = new TProfile("hResolutionZ_vs_pT",
                                      "DCA_{z} Resolution vs p_{T};p_{T} (GeV/c);#sigma_{z} (cm)",
                                      100, 0.1, 20);

    // Eta-dependent histograms
    hDCAxy_vs_eta = new TH2F("hDCAxy_vs_eta", "DCA_{xy} vs #eta;#eta;DCA_{xy} (cm)",
                             100, -2.5, 2.5, 200, -1.0, 1.0);
    hDCAz_vs_eta = new TH2F("hDCAz_vs_eta", "DCA_{z} vs #eta;#eta;DCA_{z} (cm)",
                            100, -2.5, 2.5, 200, -1.0, 1.0);
}

void DCASmearing::fillHistograms(const Track &truthPart, const Track &smearedPart)
{
    double pT = truthPart.pt();
    double eta = truthPart.eta();

    // Fill DCA distributions
    hDCAxy_Truth->Fill(truthPart.getDCAxy());
    hDCAxy_Smeared->Fill(smearedPart.getDCAxy());
    hDCAxy_Diff->Fill(smearedPart.getDCAxy() - truthPart.getDCAxy());

    hDCAz_Truth->Fill(truthPart.getDCAz());
    hDCAz_Smeared->Fill(smearedPart.getDCAz());
    hDCAz_Diff->Fill(smearedPart.getDCAz() - truthPart.getDCAz());

    // Fill pT-dependent histograms
    hDCAxy_vs_pT->Fill(pT, smearedPart.getDCAxy());
    hDCAz_vs_pT->Fill(pT, smearedPart.getDCAz());

    // Fill resolution profiles
    hResolutionXY_vs_pT->Fill(pT, std::abs(smearedPart.getDCAxy() - truthPart.getDCAxy()));
    hResolutionZ_vs_pT->Fill(pT, std::abs(smearedPart.getDCAz() - truthPart.getDCAz()));

    // Fill eta-dependent histograms
    hDCAxy_vs_eta->Fill(eta, smearedPart.getDCAxy());
    hDCAz_vs_eta->Fill(eta, smearedPart.getDCAz());
}
