#include "TaggingUtilities.h"

// Initialize static histogram pointers for DCASmearing
TH1F* DCASmearing::hDCAxy_Truth = nullptr;
TH1F* DCASmearing::hDCAxy_Smeared = nullptr;
TH1F* DCASmearing::hDCAxy_Diff = nullptr;
TH1F* DCASmearing::hDCAz_Truth = nullptr;
TH1F* DCASmearing::hDCAz_Smeared = nullptr;
TH1F* DCASmearing::hDCAz_Diff = nullptr;
TH2F* DCASmearing::hDCAxy_vs_pT = nullptr;
TH2F* DCASmearing::hDCAz_vs_pT = nullptr;
TProfile* DCASmearing::hResolutionXY_vs_pT = nullptr;
TProfile* DCASmearing::hResolutionZ_vs_pT = nullptr;
TH2F* DCASmearing::hDCAxy_vs_eta = nullptr;
TH2F* DCASmearing::hDCAz_vs_eta = nullptr;

void DCASmearing::initializeHistograms() {
    hDCAxy_Truth = new TH1F("hDCAxy_Truth", "True DCA_{xy};DCA_{xy} [cm];Entries", 200, -1, 1);
    hDCAxy_Smeared = new TH1F("hDCAxy_Smeared", "Smeared DCA_{xy};DCA_{xy} [cm];Entries", 200, -1, 1);
    hDCAxy_Diff = new TH1F("hDCAxy_Diff", "DCA_{xy} Difference;#Delta DCA_{xy} [cm];Entries", 200, -0.1, 0.1);
    
    hDCAz_Truth = new TH1F("hDCAz_Truth", "True DCA_{z};DCA_{z} [cm];Entries", 200, -1, 1);
    hDCAz_Smeared = new TH1F("hDCAz_Smeared", "Smeared DCA_{z};DCA_{z} [cm];Entries", 200, -1, 1);
    hDCAz_Diff = new TH1F("hDCAz_Diff", "DCA_{z} Difference;#Delta DCA_{z} [cm];Entries", 200, -0.1, 0.1);
    
    hDCAxy_vs_pT = new TH2F("hDCAxy_vs_pT", "DCA_{xy} vs p_{T};p_{T} [GeV/c];DCA_{xy} [cm]", 100, 0, 50, 200, -1, 1);
    hDCAz_vs_pT = new TH2F("hDCAz_vs_pT", "DCA_{z} vs p_{T};p_{T} [GeV/c];DCA_{z} [cm]", 100, 0, 50, 200, -1, 1);
    
    hResolutionXY_vs_pT = new TProfile("hResolutionXY_vs_pT", "DCA_{xy} Resolution vs p_{T};p_{T} [GeV/c];#sigma(DCA_{xy}) [cm]", 100, 0, 50);
    hResolutionZ_vs_pT = new TProfile("hResolutionZ_vs_pT", "DCA_{z} Resolution vs p_{T};p_{T} [GeV/c];#sigma(DCA_{z}) [cm]", 100, 0, 50);
    
    hDCAxy_vs_eta = new TH2F("hDCAxy_vs_eta", "DCA_{xy} vs #eta;#eta;DCA_{xy} [cm]", 100, -4, 4, 200, -1, 1);
    hDCAz_vs_eta = new TH2F("hDCAz_vs_eta", "DCA_{z} vs #eta;#eta;DCA_{z} [cm]", 100, -4, 4, 200, -1, 1);
}

void DCASmearing::fillHistograms(const Track& truthTrack, const Track& smearedTrack) {
    hDCAxy_Truth->Fill(truthTrack.getDCAxy());
    hDCAxy_Smeared->Fill(smearedTrack.getDCAxy());
    hDCAxy_Diff->Fill(smearedTrack.getDCAxy() - truthTrack.getDCAxy());
    
    hDCAz_Truth->Fill(truthTrack.getDCAz());
    hDCAz_Smeared->Fill(smearedTrack.getDCAz());
    hDCAz_Diff->Fill(smearedTrack.getDCAz() - truthTrack.getDCAz());
    
    double pT = truthTrack.pt();
    hDCAxy_vs_pT->Fill(pT, smearedTrack.getDCAxy());
    hDCAz_vs_pT->Fill(pT, smearedTrack.getDCAz());
    
    hResolutionXY_vs_pT->Fill(pT, std::abs(smearedTrack.getDCAxy() - truthTrack.getDCAxy()));
    hResolutionZ_vs_pT->Fill(pT, std::abs(smearedTrack.getDCAz() - truthTrack.getDCAz()));
    
    double eta = truthTrack.eta();
    hDCAxy_vs_eta->Fill(eta, smearedTrack.getDCAxy());
    hDCAz_vs_eta->Fill(eta, smearedTrack.getDCAz());
}

double DetectorSimulation::smear_pT(double pT) {
    double resolution = std::sqrt(kPtResA * kPtResA + 
                                kPtResB * kPtResB / (pT * pT) + 
                                kPtResC * kPtResC * pT * pT);
    return random.Gaus(pT, pT * resolution);
}

void DetectorSimulation::smearDCA(const Pythia8::Particle& part, Track& track, 
                                 const TVector3& primaryVertex, SmearingMethod method) {
    double pT = part.pT();
    double dca_resolution_xy, dca_resolution_z;

    switch (method) {
        case SmearingMethod::GAUSSIAN:
            dca_resolution_xy = 0.001;
            dca_resolution_z = 0.001;
            break;

        case SmearingMethod::PT_DEPENDENT:
            dca_resolution_xy = std::sqrt(kDCAxyOffset * kDCAxyOffset + 
                                        kDCAxySlope * kDCAxySlope / (pT * pT));
            dca_resolution_z = std::sqrt(kDCAzOffset * kDCAzOffset + 
                                       kDCAzSlope * kDCAzSlope / (pT * pT));
            break;

        case SmearingMethod::NONE:
        default:
            return;
    }

    double posX = part.xProd() * 0.1 - primaryVertex.X();
    double posY = part.yProd() * 0.1 - primaryVertex.Y();
    double posZ = part.zProd() * 0.1 - primaryVertex.Z();

    double smeared_dcaxy = std::copysign(std::sqrt(posX * posX + posY * posY), posY);
    smeared_dcaxy = random.Gaus(smeared_dcaxy, dca_resolution_xy);
    double smeared_dcaz = random.Gaus(posZ, dca_resolution_z);

    track.setDCAxy(smeared_dcaxy);
    track.setDCAz(smeared_dcaz);

    TVector3 pos(part.xProd() * 0.1, part.yProd() * 0.1, part.zProd() * 0.1);
    track.setPosition(pos);
} 