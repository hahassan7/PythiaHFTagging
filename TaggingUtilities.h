#ifdef __CLING__
R__LOAD_LIBRARY(libpythia8.dylib)
R__ADD_INCLUDE_PATH($PYTHIA_ROOT / include)
#endif

#include "JTreeHFFeatures.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/ParticleData.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TProfile.h"
#include "TParticle.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <random>

#ifdef __CLING__
R__LOAD_LIBRARY(libfastjet)
R__LOAD_LIBRARY(libsiscone)
R__LOAD_LIBRARY(libsiscone_spherical)
R__LOAD_LIBRARY(libfastjetplugins)
R__ADD_INCLUDE_PATH($FASTJET_ROOT / include)
#endif
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#ifndef __FJCORE__
#include "fastjet/GhostedAreaSpec.hh" // for area support
#endif                                // __FJCORE__

template <typename T, typename U>
float deltaR(T const &A, U const &B)
{
    double dPhi, dEta;
    if constexpr (std::is_same<U, TVector3>::value || std::is_same<U, TParticle>::value)
    {
        dPhi = TVector2::Phi_mpi_pi(A.phi() - B.Phi());
        dEta = A.eta() - B.Eta();
    }
    else
    {
        dPhi = TVector2::Phi_mpi_pi(A.phi() - B.phi());
        dEta = A.eta() - B.eta();
    }

    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

enum JetTaggingSpecies
{
    none = 0,
    charm = 1,
    beauty = 2,
    lightflavour = 3,
    lightquark = 4,
    gluon = 5
};

struct Track
{
    TVector3 pos;
    TVector3 mom;
    double dcaXY;
    double dcaZ;
    double charge;

    Track(const TVector3 &p = TVector3(), const TVector3 &m = TVector3(),
          double dxy = 0, double dz = 0, double q = 0)
        : pos(p), mom(m), dcaXY(dxy), dcaZ(dz), charge(q) {}
};

template <typename AnyJet, typename AllMCParticles>
int16_t getJetFlavor(AnyJet const &jet, AllMCParticles const &mcparticles, float maxDistance = 0.4)
{
    bool charmQuark = false;
    for (auto const &mcpart : mcparticles)
    {
        int pdgcode = mcpart.GetPdgCode();
        if (std::abs(pdgcode) == 21 || (std::abs(pdgcode) >= 1 && std::abs(pdgcode) <= 5))
        {
            double dR = deltaR(jet, mcpart);

            if (dR < maxDistance)
            {
                if (std::abs(pdgcode) == 5)
                {
                    return JetTaggingSpecies::beauty; // Beauty jet
                }
                else if (std::abs(pdgcode) == 4)
                {
                    charmQuark = true;
                }
            }
        }
    }

    if (charmQuark)
    {
        return JetTaggingSpecies::charm; // Charm jet
    }

    return JetTaggingSpecies::lightflavour; // Light flavor jet
}

//________________________________________________________________________
bool isBMeson(int pc)
{
    int bPdG[] = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 20513, 20523, 515, 525, 531,
                  10531, 533, 10533, 20533, 535, 541, 10541, 543, 10543, 20543, 545, 551, 10551, 100551,
                  110551, 200551, 210551, 553, 10553, 20553, 30553, 100553, 110553, 120553, 130553, 200553, 210553, 220553,
                  300553, 9000533, 9010553, 555, 10555, 20555, 100555, 110555, 120555, 200555, 557, 100557, 5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 5312, 5322, 5314, 5324, 5332, 5334, 5142, 5242, 5412, 5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 5522, 5514, 5524, 5532,
                  5534, 5542, 5544, 5554};
    for (int i = 0; i < (int)(sizeof(bPdG) / sizeof(int)); ++i)
        if (std::abs(pc) == bPdG[i])
            return true;
    return false;
}
//________________________________________________________________________
bool isDMeson(int pc)
{
    int bPdG[] = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20431, 20423, 415,
                  425, 431, 10431, 433, 10433, 20433, 435, 441, 10441, 100441, 443, 10443, 20443,
                  100443, 30443, 9000443, 9010443, 9020443, 445, 100445, 4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};
    for (int i = 0; i < (int)(sizeof(bPdG) / sizeof(int)); ++i)
        if (std::abs(pc) == bPdG[i])
            return true;
    return false;
}

// Detector resolution simulation
class DetectorSimulation
{
private:
    TRandom3 *rand;

    // Track parameter resolutions
    const double dcaXYParam[2] = {4.0e-4, 20.0e-4};  // [µm] (offset, pt_slope)
    const double dcaZParam[2] = {4.0e-4, 20.0e-4};   // [µm]
    const double ptParam[3] = {0.005, 0.01, 0.0003}; // pt resolution
    const double phiRes = 0.001;                     // [rad]
    const double thetaRes = 0.001;                   // [rad]

    // Vertex resolutions
    const double vtxXYRes = 12.0; // [µm] in pp
    const double vtxZRes = 15.0;  // [µm] in pp

public:
    DetectorSimulation() : rand(new TRandom3(0)) {}

    // Different smearing methods
    enum class SmearingMethod
    {
        GAUSSIAN,
        PT_DEPENDENT,
        NONE // For truth values
    };

    void smearDCA(const Pythia8::Particle &parttobesmeared, Track &smearedTrack, SmearingMethod method = SmearingMethod::PT_DEPENDENT)
    {
        double xProd = parttobesmeared.xProd() * 0.1;
        double yProd = parttobesmeared.yProd() * 0.1;
        double zProd = parttobesmeared.zProd() * 0.1;
        double pT = parttobesmeared.pT();
        double eta = parttobesmeared.eta();

        switch (method)
        {
        case SmearingMethod::PT_DEPENDENT:
        {

            if (pT < 0.1)
            {
                pT = 0.1;
            }

            std::cout << "pT: " << std::fixed << std::setprecision(10) << pT << " eta: " << std::fixed << std::setprecision(10) << eta << std::endl;

            // Calculate resolutions
            double sigmaXY = std::sqrt(dcaXYParam[0] * dcaXYParam[0] + (dcaXYParam[1] * dcaXYParam[1]) / (pT * pT));
            double sigmaZ = std::sqrt(dcaZParam[0] * dcaZParam[0] + (dcaZParam[1] * dcaZParam[1]) / (pT * pT));

            std::cout << "sigmaXY: " << std::fixed << std::setprecision(10) << sigmaXY << " sigmaZ: " << std::fixed << std::setprecision(10) << sigmaZ << std::endl;

            // Add eta dependence (resolution degrades at higher eta)
            double etaFactor = 1.0 + 0.1 * std::abs(eta); // Approximate effect
            sigmaXY *= etaFactor;
            sigmaZ *= etaFactor;

            // Apply smearing
            smearedTrack.dcaXY = rand->Gaus(std::copysign(std::sqrt(xProd * xProd + yProd * yProd), yProd), sigmaXY);
            smearedTrack.dcaZ = rand->Gaus(zProd, sigmaZ);
            break;
        }

        case SmearingMethod::GAUSSIAN:
        {
            // Fixed resolution based on average ITS2 performance
            const double sigmaXY = 0.0010; // 10 μm -> cm
            const double sigmaZ = 0.0010;  // 10 μm -> cm

            smearedTrack.dcaXY = rand->Gaus(std::copysign(std::sqrt(xProd * xProd + yProd * yProd), yProd), sigmaXY);
            smearedTrack.dcaZ = rand->Gaus(zProd, sigmaZ);
            break;
        }

        case SmearingMethod::NONE:
        default:
            break;
        }
    }

    // Function to smear pT based on ALICE ITS resolution
    double smear_pT(double pT)
    {
        // Compute momentum resolution
        double sigma_pT_over_pT = std::sqrt(ptParam[0] * ptParam[0] + (ptParam[1] / pT) * (ptParam[1] / pT) + (ptParam[2] * pT) * (ptParam[2] * pT));

        // Apply smearing
        double smeared_pT = rand->Gaus(pT, sigma_pT_over_pT);
        return smeared_pT;
    }

    Track smearTrack(const Pythia8::Particle &part)
    {
        Track params = Track();

        smearDCA(part, params);

        // Smear momentum
        double pt = smear_pT(part.pT());
        double phi = part.phi() + rand->Gaus(0, phiRes);
        double theta = part.theta() + rand->Gaus(0, thetaRes);

        params.mom.SetPtThetaPhi(pt, theta, phi);
        params.charge = part.charge();

        // Calculate position at DCA
        // TODO: Check if this is correct
        params.pos = TVector3(std::abs(params.dcaXY) * sin(phi),
                              std::abs(params.dcaXY) * cos(phi),
                              params.dcaZ);

        return params;
    }

    void smearPrimaryVertex(TVector3 &primaryVtx)
    {
        primaryVtx.SetX(rand->Gaus(primaryVtx.X(), vtxXYRes * 1e-4));
        primaryVtx.SetY(rand->Gaus(primaryVtx.Y(), vtxXYRes * 1e-4));
        primaryVtx.SetZ(rand->Gaus(primaryVtx.Z(), vtxZRes * 1e-4));
    }

    ~DetectorSimulation()
    {
        delete rand;
    }
};

// Secondary vertex finder
class SecondaryVertexFinder
{
private:
    const double maxDCA = 1.0; // cm
    const double maxChi2 = 5.0;
    const double minVtxSignif = 0.3;
    const double maxRadialDist = 10.0; // cm
    const double maxZDist = 12.0;      // cm
    const double minTrackPt = 0.3;     // GeV/c

    struct TrackAtVertex
    {
        TVector3 pos;
        TVector3 mom;
        double weight;
    };

public:
    struct SecVtxInfo
    {
        TVector3 position;
        double mass;
        double significance;
        int nTracks;
        double chi2;
    };

    // Modified to return vector of SVs
    std::vector<SecVtxInfo> findSecondaryVertices(const std::vector<Track> &inputTracks)
    {
        std::vector<SecVtxInfo> vertices;

        // First apply track cuts
        std::vector<Track> goodTracks;
        for (const auto &track : inputTracks)
        {
            if (passesTrackCuts(track))
            {
                goodTracks.push_back(track);
            }
        }

        // Need at least 3 tracks after cuts
        if (goodTracks.size() < 3)
            return vertices;

        // Try all possible combinations of 3 tracks
        for (size_t i = 0; i < goodTracks.size() - 2; ++i)
        {
            for (size_t j = i + 1; j < goodTracks.size() - 1; ++j)
            {
                for (size_t k = j + 1; k < goodTracks.size(); ++k)
                {
                    std::vector<Track> trackTriplet = {
                        goodTracks[i], goodTracks[j], goodTracks[k]};

                    SecVtxInfo vtxInfo = findSingleSecondaryVertex(trackTriplet);

                    // Apply all quality cuts
                    if (vtxInfo.chi2 < maxChi2 &&
                        vtxInfo.significance > minVtxSignif &&
                        vtxInfo.nTracks == 3 &&
                        passesVertexCuts(vtxInfo))
                    {
                        vertices.push_back(vtxInfo);
                    }
                }
            }
        }

        return vertices;
    }

    // Add monitoring histograms
    static void fillMonitoringHistograms(const std::vector<SecVtxInfo> &vertices)
    {
        static TH1F *hSVRadialDist = new TH1F("hSVRadialDist", "SV Radial Distance;R (cm);Counts", 100, 0, 20);
        static TH1F *hSVZDist = new TH1F("hSVZDist", "SV Z Distance;|Z| (cm);Counts", 120, 0, 20);
        static TH1F *hTrackPt = new TH1F("hTrackPt", "Track p_{T} at SV;p_{T} (GeV/c);Counts", 100, 0, 100);

        for (const auto &vtx : vertices)
        {
            double r = std::sqrt(vtx.position.X() * vtx.position.X() +
                                 vtx.position.Y() * vtx.position.Y());
            hSVRadialDist->Fill(r);
            hSVZDist->Fill(std::abs(vtx.position.Z()));

            // Track pT histogram would need access to the tracks
            // You might need to modify the SecVtxInfo struct to store track information
        }
    }

private:
    // Renamed original function to handle single vertex finding
    SecVtxInfo findSingleSecondaryVertex(const std::vector<Track> &tracks)
    {
        SecVtxInfo vtxInfo;
        if (tracks.size() < 2)
            return vtxInfo;

        // Initial vertex position (weighted average of track DCAs)
        TVector3 vtxPos(0, 0, 0);
        double totalWeight = 0;

        for (const auto &track : tracks)
        {
            double weight = track.mom.Pt(); // Weight by pT
            vtxPos += track.pos * weight;
            totalWeight += weight;
        }
        vtxPos *= 1.0 / totalWeight;

        // Iterative vertex finding
        const int maxIter = 100;
        double prevChi2 = 1e9;

        for (int iter = 0; iter < maxIter; ++iter)
        {
            std::vector<TrackAtVertex> tracksAtVtx;

            // Extrapolate tracks to current vertex position
            for (const auto &track : tracks)
            {
                TVector3 diff = vtxPos - track.pos;
                if (diff.Mag() > maxDCA)
                    continue;

                TrackAtVertex trkVtx;
                trkVtx.pos = track.pos;
                trkVtx.mom = track.mom;
                trkVtx.weight = 1.0 / (diff.Mag2() + 0.01); // Weight by distance

                tracksAtVtx.push_back(trkVtx);
            }

            if (tracksAtVtx.size() < 2)
                break;

            // Update vertex position
            TVector3 newPos(0, 0, 0);
            double newWeight = 0;

            for (const auto &trk : tracksAtVtx)
            {
                newPos += trk.pos * trk.weight;
                newWeight += trk.weight;
            }
            newPos *= 1.0 / newWeight;

            // Check convergence
            double chi2 = calcVertexChi2(tracks, newPos);
            if (fabs(chi2 - prevChi2) < 0.01)
                break;

            vtxPos = newPos;
            prevChi2 = chi2;
        }

        // Calculate vertex properties
        vtxInfo.position = vtxPos;
        vtxInfo.chi2 = prevChi2;
        vtxInfo.nTracks = 0;

        TLorentzVector vtxMom;
        for (const auto &track : tracks)
        {
            if ((track.pos - vtxPos).Mag() < maxDCA)
            {
                vtxInfo.nTracks++;
                TLorentzVector p;
                p.SetXYZM(track.mom.X(), track.mom.Y(), track.mom.Z(), 0.139); // pion mass
                vtxMom += p;
            }
        }

        vtxInfo.mass = vtxMom.M();
        vtxInfo.significance = vtxPos.Perp() / sqrt(prevChi2);

        return vtxInfo;
    }

    double calcVertexChi2(const std::vector<Track> &tracks,
                          const TVector3 &vtxPos)
    {
        double chi2 = 0;
        for (const auto &track : tracks)
        {
            TVector3 diff = track.pos - vtxPos;
            // Simplified chi2 calculation
            chi2 += diff.Mag2();
        }
        return chi2;
    }

    bool passesVertexCuts(const SecVtxInfo &vtx)
    {
        // Check radial distance
        double r = std::sqrt(vtx.position.X() * vtx.position.X() +
                             vtx.position.Y() * vtx.position.Y());
        if (r >= maxRadialDist)
            return false;

        // Check Z distance
        if (std::abs(vtx.position.Z()) >= maxZDist)
            return false;

        return true;
    }

    bool passesTrackCuts(const Track &track)
    {
        // Check track pT
        if (track.mom.Pt() <= minTrackPt)
            return false;

        return true;
    }
};

struct SVCalculations
{
    static double calculateDecayLength2D(const TVector3 &primaryVtx, const TVector3 &secondaryVtx)
    {
        double dx = secondaryVtx.X() - primaryVtx.X();
        double dy = secondaryVtx.Y() - primaryVtx.Y();
        return std::sqrt(dx * dx + dy * dy);
    }

    static double calculateDecayLength3D(const TVector3 &primaryVtx, const TVector3 &secondaryVtx)
    {
        return (secondaryVtx - primaryVtx).Mag();
    }

    static double calculateSVPt(const std::vector<TVector3> &daughterMomenta)
    {
        TVector3 sumMom;
        for (const auto &mom : daughterMomenta)
        {
            sumMom += mom;
        }
        return sumMom.Pt();
    }

    static double calculateDeltaRToJet(const TVector3 &svPosition, const fastjet::PseudoJet &jet)
    {
        return deltaR(jet, svPosition);
    }

    static double calculateSVfE(double svEnergy, double jetEnergy)
    {
        return svEnergy / jetEnergy;
    }

    static double calculateIP2D(const TVector3 &primaryVtx, const TVector3 &secondaryVtx, const fastjet::PseudoJet &jet)
    {
        TVector2 sv2D(secondaryVtx.X() - primaryVtx.X(), secondaryVtx.Y() - primaryVtx.Y());
        TVector2 jet2D(jet.px(), jet.py());
        jet2D = jet2D.Unit();

        // IP2D is the magnitude of the projection perpendicular to the jet axis
        return std::abs(sv2D.X() * jet2D.Y() - sv2D.Y() * jet2D.X());
    }

    static double calculateCPA(const TVector3 &primaryVtx, const TVector3 &secondaryVtx, const TVector3 &svMomentum)
    {
        TVector3 flightLine = secondaryVtx - primaryVtx;
        return flightLine.Unit().Dot(svMomentum.Unit());
    }

    static double calculateDispersion(const std::vector<TVector3> &trackPositions, const TVector3 &svPosition)
    {
        if (trackPositions.empty())
            return 0.0;

        double sumDistanceSquared = 0.0;
        for (const auto &trackPos : trackPositions)
        {
            TVector3 diff = trackPos - svPosition;
            sumDistanceSquared += diff.Mag2();
        }
        return std::sqrt(sumDistanceSquared / trackPositions.size());
    }
};

class VertexFinder
{
public:
    static TVector3 findPrimaryVertex(const std::vector<Track> &tracks,
                                      double convergenceDist = 0.01, // in cm
                                      int maxIterations = 100,
                                      double temperatureScale = 2.0, // Controls how fast weights change
                                      double chi2Cut = 9.0)
    { // 3 sigma cut in 3D
        if (tracks.empty())
        {
            return TVector3(0., 0., 0.);
        }

        // Initial vertex position: mean of track DCA points
        TVector3 vtxPos(0., 0., 0.);
        for (const auto &track : tracks)
        {
            vtxPos += track.pos;
        }
        vtxPos *= (1.0 / tracks.size());

        // Iterative fitting
        TVector3 lastPos = vtxPos;
        int iter = 0;
        bool converged = false;

        while (!converged && iter < maxIterations)
        {
            // Reset sums for weighted average
            TVector3 newPos(0., 0., 0.);
            double weightSum = 0.;

            // Update weights and calculate new vertex position
            for (auto &track : tracks)
            {
                // Calculate chi2 distance to current vertex estimate
                TVector3 diff = track.pos - vtxPos;

                // Simplified chi2 calculation (you might want to use full covariance matrix)
                double chi2 = diff.Mag2(); // Simplified - assumes spherical errors

                // Calculate weight using temperature
                double weight = 1.0;
                if (chi2 < chi2Cut)
                {
                    weight = exp(-chi2 / (2.0 * temperatureScale));
                }
                else
                {
                    weight = 0.0; // Remove outliers
                }

                // Update sums
                newPos += track.pos * weight;
                weightSum += weight;
            }

            // Calculate new vertex position
            if (weightSum > 0)
            {
                newPos *= (1.0 / weightSum);
            }
            else
            {
                // If all tracks were rejected, keep the old position
                newPos = vtxPos;
            }

            // Check convergence
            TVector3 posDiff = newPos - lastPos;
            if (posDiff.Mag() < convergenceDist)
            {
                converged = true;
            }

            lastPos = vtxPos;
            vtxPos = newPos;
            iter++;
        }

        return vtxPos;
    }
};

class DCASmearing
{
public:
    // Update histogram ranges to match ALICE scales
    static void initializeHistograms()
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

    static void fillHistograms(const Track &truthPart, const Track &smearedPart)
    {
        double pT = truthPart.mom.Pt();
        double eta = truthPart.mom.Eta();

        // Fill DCA distributions
        hDCAxy_Truth->Fill(truthPart.dcaXY);
        hDCAxy_Smeared->Fill(smearedPart.dcaXY);
        hDCAxy_Diff->Fill(smearedPart.dcaXY - truthPart.dcaXY);

        hDCAz_Truth->Fill(truthPart.dcaZ);
        hDCAz_Smeared->Fill(smearedPart.dcaZ);
        hDCAz_Diff->Fill(smearedPart.dcaZ - truthPart.dcaZ);

        // Fill pT-dependent histograms
        hDCAxy_vs_pT->Fill(pT, smearedPart.dcaXY);
        hDCAz_vs_pT->Fill(pT, smearedPart.dcaZ);

        // Fill resolution profiles
        hResolutionXY_vs_pT->Fill(pT, std::abs(smearedPart.dcaXY - truthPart.dcaXY));
        hResolutionZ_vs_pT->Fill(pT, std::abs(smearedPart.dcaZ - truthPart.dcaZ));

        // Fill eta-dependent histograms
        hDCAxy_vs_eta->Fill(eta, smearedPart.dcaXY);
        hDCAz_vs_eta->Fill(eta, smearedPart.dcaZ);
    }

private:
    // All histogram declarations
    static TH1F *hDCAxy_Truth;
    static TH1F *hDCAxy_Smeared;
    static TH1F *hDCAxy_Diff;
    static TH1F *hDCAz_Truth;
    static TH1F *hDCAz_Smeared;
    static TH1F *hDCAz_Diff;
    static TH2F *hDCAxy_vs_pT;
    static TH2F *hDCAz_vs_pT;
    static TProfile *hResolutionXY_vs_pT;
    static TProfile *hResolutionZ_vs_pT;
    // Add eta-dependent histograms
    static TH2F *hDCAxy_vs_eta;
    static TH2F *hDCAz_vs_eta;
};

// Initialize all static histogram pointers
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

TParticle convertToTParticle(const Pythia8::Particle &part)
{
    return TParticle(
        part.id(),                                             // pdg code
        part.status(),                                         // status
        part.mother1(),                                        // mother1
        part.mother2(),                                        // mother2
        part.daughter1(),                                      // daughter1
        part.daughter2(),                                      // daughter2
        part.px(), part.py(), part.pz(), part.e(),             // momentum
        part.xProd(), part.yProd(), part.zProd(), part.tProd() // production vertex
    );
}
