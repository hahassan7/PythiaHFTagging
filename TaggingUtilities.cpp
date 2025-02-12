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
bool isBMeson(int pc) {
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
bool isDMeson(int pc) {
    static const int dPdG[] = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20431, 
                  20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, 441, 10441, 
                  100441, 443, 10443, 20443, 100443, 30443, 9000443, 9010443, 9020443, 
                  445, 100445};
    for (int i = 0; i < (int)(sizeof(dPdG) / sizeof(int)); ++i)
        if (std::abs(pc) == dPdG[i])
            return true;
    return false;
}

void DetectorSimulation::smearDCA(const Pythia8::Particle &parttobesmeared, Track &smearedTrack, SmearingMethod method)
{
    double xProd = parttobesmeared.xProd() * 0.1; // Convert to cm
    double yProd = parttobesmeared.yProd() * 0.1;
    double zProd = parttobesmeared.zProd() * 0.1;
    double pT = parttobesmeared.pT();

    double sigmaXY, sigmaZ;

    switch (method)
    {
    case SmearingMethod::GAUSSIAN:
        sigmaXY = dcaXYParam[0];
        sigmaZ = dcaZParam[0];
        smearedTrack.dxy = random.Gaus(std::copysign(std::sqrt(xProd * xProd +
                                                              yProd * yProd),
                                                    yProd),
                                     sigmaXY);
        smearedTrack.dz = random.Gaus(zProd, sigmaZ);
        break;

    case SmearingMethod::PT_DEPENDENT:
        sigmaXY = std::sqrt(dcaXYParam[0] * dcaXYParam[0] + (dcaXYParam[1] / pT) * (dcaXYParam[1] / pT));
        sigmaZ = std::sqrt(dcaZParam[0] * dcaZParam[0] + (dcaZParam[1] / pT) * (dcaZParam[1] / pT));
        smearedTrack.dxy = random.Gaus(std::copysign(std::sqrt(xProd * xProd +
                                                              yProd * yProd),
                                                    yProd),
                                     sigmaXY);
        smearedTrack.dz = random.Gaus(zProd, sigmaZ);
        break;

    default:
        smearedTrack.dxy = std::copysign(std::sqrt(xProd * xProd + yProd * yProd), yProd);
        smearedTrack.dz = zProd;
        break;
    }
}

void DCASmearing::initializeHistograms() {
    // ... implementation ...
}

void DCASmearing::fillHistograms(const Track &truthPart, const Track &smearedPart) {
    // ... implementation ...
} 