#pragma once

// Include FastJet compatibility layer first
// #include "FastJetCompat.h"

// First include system headers
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <memory>
#include <array>

// Custom headers
#include "JParticle.h"

// ROOT headers
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TParticle.h"

// Pythia headers
#include "Pythia8/Pythia.h"
#include "Pythia8/ParticleData.h"

// FastJet headers
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#ifndef __FJCORE__
#include "fastjet/GhostedAreaSpec.hh"
#endif

#ifdef __CLING__
R__LOAD_LIBRARY(libpythia8)
R__ADD_INCLUDE_PATH($PYTHIA_ROOT / include)
#endif

#ifdef __CLING__
R__LOAD_LIBRARY(libfastjet)
R__LOAD_LIBRARY(libsiscone)
R__LOAD_LIBRARY(libsiscone_spherical)
R__LOAD_LIBRARY(libfastjetplugins)
R__ADD_INCLUDE_PATH($FASTJET_ROOT / include)
#endif

// Forward declarations
bool isBMeson(int pc);
bool isDMeson(int pc);

// Define DebugLevel at the very top
enum class DebugLevel
{
    SILENT = 0,  // No output
    ERROR = 1,   // Only errors
    WARNING = 2, // Errors and warnings
    INFO = 3,    // Basic run information
    DEBUG = 4,   // Detailed debugging information
    VERBOSE = 5  // Very detailed debugging information
};

// Add logging helper function right after enum
inline void log(DebugLevel currentLevel, DebugLevel messageLevel, const std::string &message)
{
    if (static_cast<int>(messageLevel) <= static_cast<int>(currentLevel))
    {
        switch (messageLevel)
        {
        case DebugLevel::ERROR:
            std::cout << "\033[1;31m[ERROR]\033[0m " << message << std::endl;
            break;
        case DebugLevel::WARNING:
            std::cout << "\033[1;33m[WARNING]\033[0m " << message << std::endl;
            break;
        case DebugLevel::INFO:
            std::cout << "\033[1;32m[INFO]\033[0m " << message << std::endl;
            break;
        case DebugLevel::DEBUG:
            std::cout << "\033[1;34m[DEBUG]\033[0m " << message << std::endl;
            break;
        case DebugLevel::VERBOSE:
            std::cout << "\033[1;35m[VERBOSE]\033[0m " << message << std::endl;
            break;
        default:
            break;
        }
    }
}

class Track
{
private:
    TVector3 pos;
    TVector3 mom;
    double dcaXY;
    double dcaZ;
    double charge;
    int index;

public:
    // Default constructor - initialize everything
    Track() : pos(TVector3()), mom(TVector3()), dcaXY(0), dcaZ(0), charge(0), index(-1) {}

    // Main constructor
    Track(const TVector3 &p, const TVector3 &m, double dxy = 0, double dz = 0, double q = 0)
        : pos(p), mom(m), dcaXY(dxy), dcaZ(dz), charge(q), index(-1) // Initialize index
    {
        // Validate inputs
        if (std::isnan(pos.X()) || std::isnan(pos.Y()) || std::isnan(pos.Z()))
        {
            throw std::runtime_error("Invalid position in Track constructor");
        }
        if (std::isnan(mom.X()) || std::isnan(mom.Y()) || std::isnan(mom.Z()))
        {
            throw std::runtime_error("Invalid momentum in Track constructor");
        }
        if (std::isnan(dxy) || std::isnan(dz) || std::isnan(charge))
        {
            throw std::runtime_error("Invalid impact parameters or charge in Track constructor");
        }
    }

    // Copy constructor
    Track(const Track &other)
        : pos(other.pos),
          mom(other.mom),
          dcaXY(other.dcaXY),
          dcaZ(other.dcaZ),
          charge(other.charge),
          index(other.index) // Copy the index
    {
        // Validate copied data
        if (std::isnan(pos.X()) || std::isnan(pos.Y()) || std::isnan(pos.Z()))
        {
            throw std::runtime_error("Invalid position in Track copy constructor");
        }
        if (std::isnan(mom.X()) || std::isnan(mom.Y()) || std::isnan(mom.Z()))
        {
            throw std::runtime_error("Invalid momentum in Track copy constructor");
        }
        if (std::isnan(dcaXY) || std::isnan(dcaZ) || std::isnan(charge))
        {
            throw std::runtime_error("Invalid impact parameters or charge in Track copy constructor");
        }
    }

    // Assignment operator
    Track &operator=(const Track &other)
    {
        if (this != &other)
        {
            pos = other.pos;
            mom = other.mom;
            dcaXY = other.dcaXY;
            dcaZ = other.dcaZ;
            charge = other.charge;
            index = other.index; // Copy the index
        }
        return *this;
    }

    // Getters
    const TVector3 &getPosition() const { return pos; }
    const TVector3 &getMomentum() const { return mom; }
    double getDCAxy() const { return dcaXY; }
    double getDCAz() const { return dcaZ; }
    double getCharge() const { return charge; }
    int getIndex() const { return index; }

    // Setters
    void setPosition(const TVector3 &p) { pos = p; }
    void setMomentum(const TVector3 &m) { mom = m; }
    void setDCAxy(double dxy) { dcaXY = dxy; }
    void setDCAz(double dz) { dcaZ = dz; }
    void setCharge(double q) { charge = q; }
    void setIndex(int idx)
    {
        if (idx < 0)
        {
            throw std::invalid_argument("Track index cannot be negative");
        }
        index = idx;
    }

    // Utility functions
    double pt() const { return mom.Pt(); }
    double eta() const { return mom.Eta(); }
    double phi() const { return mom.Phi(); }
    double p() const { return mom.Mag(); }

    // Friend class declaration for DetectorSimulation
    friend class DetectorSimulation;
};

// Change enum class to enum for implicit conversion to int16_t
enum JetTaggingSpecies
{
    lightflavour = 0,
    charm = 1,
    beauty = 2
};

// Update deltaR template function to handle different types
template <typename T, typename U>
float deltaR(T const &A, U const &B)
{
    double dPhi, dEta;
    if constexpr (std::is_same<U, TLorentzVector>::value)
    {
        dPhi = TVector2::Phi_mpi_pi(A.phi() - B.Phi());
        dEta = A.eta() - B.Eta();
    }
    else if constexpr (std::is_same<U, TVector3>::value)
    {
        // For TVector3, calculate phi and eta
        dPhi = TVector2::Phi_mpi_pi(A.phi() - B.Phi());
        dEta = A.eta() - B.PseudoRapidity();
    }
    else
    {
        dPhi = TVector2::Phi_mpi_pi(A.phi() - B.phi());
        dEta = A.eta() - B.eta();
    }
    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

/**
 * return geometric sign which is calculated scalar product between jet axis with DCA (track propagated to PV )
 * positive and negative value are expected from primary vertex
 * positive value is expected from secondary vertex
 *
 * @param jet
 * @param track which is needed aod::JTrackExtras
 */
template <typename T>
int getGeoSign(T const &jet, Track const &track, TVector3 const &pVertex,
               DebugLevel debugLevel = DebugLevel::INFO)
{
    double posX = track.getPosition().X() - pVertex.X();
    double posY = track.getPosition().Y() - pVertex.Y();
    double posZ = track.getDCAz();

    auto sign = TMath::Sign(1, posX * jet.px() + posY * jet.py() + posZ * jet.pz());
    if (sign < -1 || sign > 1)
        log(debugLevel, DebugLevel::INFO, Form("Sign is %d", sign));

    return sign;
}

// Update getJetFlavor template function to include debug level
template <typename AnyJet, typename AllMCParticles>
int16_t getJetFlavor(AnyJet const &jet,
                     AllMCParticles const &mcparticles,
                     float maxDistance = 0.4,
                     DebugLevel debugLevel = DebugLevel::INFO)
{
    bool charmQuark = false;
    for (auto const &mcpart : mcparticles)
    {
        int pdgcode = mcpart.pdgId();
        if (std::abs(pdgcode) == 21 || (std::abs(pdgcode) >= 1 && std::abs(pdgcode) <= 5))
        {
            double dR = deltaR(jet, mcpart);
            log(debugLevel, DebugLevel::VERBOSE,
                "Checking particle PDG=" + std::to_string(pdgcode) +
                    " deltaR=" + std::to_string(dR));

            if (dR < maxDistance)
            {
                if (std::abs(pdgcode) == 5)
                {
                    log(debugLevel, DebugLevel::DEBUG, "Found b-jet match");
                    return beauty; // Beauty jet
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
        log(debugLevel, DebugLevel::DEBUG, "Found c-jet match");
        return charm; // Charm jet
    }

    log(debugLevel, DebugLevel::DEBUG, "No heavy flavor match found");
    return lightflavour; // Light flavor jet
}

enum class SmearingMethod
{
    NONE,
    GAUSSIAN,
    PT_DEPENDENT
};

class DetectorSimulation
{
private:
    TRandom3 random;

    // Resolution parameters
    static constexpr double kDCAxyOffset = 7.0e-4; // 4 μm -> 0.004 cm
    static constexpr double kDCAxySlope = 23.0e-4; // 20 μm -> 0.02 cm
    static constexpr double kDCAzOffset = 6.0e-4;  // 4 μm -> 0.004 cm
    static constexpr double kDCAzSlope = 25.0e-4;  // 20 μm -> 0.02 cm
    static constexpr double kPtResA = 0.005;       // Constant term
    static constexpr double kPtResB = 0.0008;      // 1/pT term
    static constexpr double kPtResC = 0.001;       // pT term
    static constexpr double kPhiRes = 0.001;       // [rad]
    static constexpr double kThetaRes = 0.001;     // [rad]
    static constexpr double kVtxXYRes = 12e-4;     // 12 μm -> 0.0012 cm
    static constexpr double kVtxZRes = 15e-4;      // 15 μm -> 0.0015 cm

    DebugLevel debugLevel;

public:
    DetectorSimulation(DebugLevel level = DebugLevel::INFO)
        : random(0), debugLevel(level) {} // Initialize with seed 0 for reproducibility

    // New method to smear primary vertex
    TVector3 smearPrimaryVertex(const TVector3 &trueVertex)
    {
        return TVector3(
            random.Gaus(trueVertex.X(), kVtxXYRes),
            random.Gaus(trueVertex.Y(), kVtxXYRes),
            random.Gaus(trueVertex.Z(), kVtxZRes));
    }

    // Function to smear pT based on ALICE ITS resolution
    double smear_pT(double pT);

    void smearDCA(const Pythia8::Particle &parttobesmeared, Track &smearedTrack, const TVector3 &primaryVertex, SmearingMethod method);

    void smearTrack(const Pythia8::Particle &part, Track &params, const TVector3 &primaryVertex)
    {
        log(debugLevel, DebugLevel::DEBUG, "Smearing track with parameters:");
        log(debugLevel, DebugLevel::DEBUG, "Original position: (" + std::to_string(part.xProd() * 0.1) + ", " + std::to_string(part.yProd() * 0.1) + ", " + std::to_string(part.zProd() * 0.1) + ") cm");
        log(debugLevel, DebugLevel::DEBUG, "Original momentum: (" + std::to_string(part.px()) + ", " + std::to_string(part.py()) + ", " + std::to_string(part.pz()) + ") GeV");

        double smeared_pT = smear_pT(part.pT());
        double smeared_phi = part.phi() + random.Gaus(0.0, kPhiRes);
        double smeared_theta = part.theta() + random.Gaus(0.0, kThetaRes);

        TVector3 smearedMom;
        smearedMom.SetPtThetaPhi(smeared_pT, smeared_theta, smeared_phi);
        params.setMomentum(smearedMom);
        params.setCharge(part.charge());

        // Apply DCA smearing
        smearDCA(part, params, primaryVertex, SmearingMethod::PT_DEPENDENT);

        log(debugLevel, DebugLevel::DEBUG, "Created track with parameters:");
        log(debugLevel, DebugLevel::DEBUG, "Position: (" + std::to_string(params.pos.X()) + ", " + std::to_string(params.pos.Y()) + ", " + std::to_string(params.pos.Z()) + ") cm");
        log(debugLevel, DebugLevel::DEBUG, "Momentum: (" + std::to_string(params.mom.X()) + ", " + std::to_string(params.mom.Y()) + ", " + std::to_string(params.mom.Z()) + ") GeV");
        log(debugLevel, DebugLevel::DEBUG, "Impact parameters: dxy = " + std::to_string(params.getDCAxy()) + " cm, dz = " + std::to_string(params.getDCAz()) + " cm");
    }

    ~DetectorSimulation()
    {
        // Destructor implementation if needed
    }
};

// Secondary vertex finder
class SecondaryVertexFinder
{
private:
    const double maxDCA = 100.0; // cm
    const double maxChi2 = 20.0;
    const double minVtxSignif = 0.3;
    const double maxRadialDist = 10.0; // cm
    const double maxZDist = 12.0;      // cm
    const double minTrackPt = 0.3;     // GeV/c
    const double maxDCAxy = 100.0;       // cm
    const double maxDCAz = 100.0;        // cm
    const double minSVPt = 1.0;        // GeV
    DebugLevel debugLevel;

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
        TLorentzVector momentum;
        double significance = 0.0;
        int nTracks = 0;
        double chi2 = 0.0;
        std::vector<size_t> indices;

        SecVtxInfo() : position(0, 0, 0),
                       momentum(0, 0, 0, 0),
                       significance(0),
                       nTracks(0),
                       chi2(0),
                       indices() {}
    };

    SecondaryVertexFinder(DebugLevel level = DebugLevel::INFO) : debugLevel(level) {}

    // Modified to return vector of SVs
    std::vector<SecVtxInfo> findSecondaryVertices(const std::vector<Track> &inputTracks)
    {
        log(debugLevel, DebugLevel::DEBUG, "Starting secondary vertex finding...");
        log(debugLevel, DebugLevel::DEBUG, "Input tracks: " + std::to_string(inputTracks.size()));

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
                    std::vector<Track> trackTriplet = {goodTracks[i], goodTracks[j], goodTracks[k]};

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

        for (const auto &vtx : vertices)
        {
            double r = std::sqrt(vtx.position.X() * vtx.position.X() +
                                 vtx.position.Y() * vtx.position.Y());
            hSVRadialDist->Fill(r);
            hSVZDist->Fill(std::abs(vtx.position.Z()));
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
            double weight = track.pt() * track.pt(); // Weight by pT
            vtxPos += track.getPosition() * weight;
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
                TVector3 diff = vtxPos - track.getPosition();
                if (diff.Mag() > maxDCA)
                    continue;

                TrackAtVertex trkVtx;
                trkVtx.pos = track.getPosition();
                trkVtx.mom = track.getMomentum();
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
            if ((track.getPosition() - vtxPos).Mag() < maxDCA)
            {
                vtxInfo.nTracks++;
                TLorentzVector p;
                p.SetXYZM(track.getMomentum().X(), track.getMomentum().Y(), track.getMomentum().Z(), 0.139); // pion mass
                vtxMom += p;
            }
        }

        vtxInfo.momentum = vtxMom;
        vtxInfo.significance = vtxPos.Perp() / sqrt(prevChi2);
        vtxInfo.indices = std::vector<size_t>{};

        for (const auto &track : tracks)
        {
            vtxInfo.indices.push_back(track.getIndex());
        }

        return vtxInfo;
    }

    double calcVertexChi2(const std::vector<Track> &tracks,
                          const TVector3 &vtxPos)
    {
        double chi2 = 0;
        for (const auto &track : tracks)
        {
            TVector3 diff = track.getPosition() - vtxPos;
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

        // Check SV pT
        if (vtx.momentum.Pt() < minSVPt)
            return false;

        return true;
    }

    bool passesTrackCuts(const Track &track)
    {
        // Check track pT
        if (track.pt() <= minTrackPt)
            return false;

        if (std::abs(track.getDCAxy()) > maxDCAxy || std::abs(track.getDCAz()) > maxDCAz) {
            return false;
        }

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

    /// Calculates impact parameter in the bending plane of the particle w.r.t. a point
    /// \param point  position of the point
    /// \param posSV  position of the secondary vertex
    /// \param mom  particle momentum array
    /// \return impact parameter in {x, y}
    static double calculateIPxy(const TVector3 &point, const TVector3 &posSV, const TVector3 &mom)
    {
        TVector2 flightLineXY = posSV.XYvector() - point.XYvector();
        auto k = (flightLineXY.X() * mom.X() + flightLineXY.Y() * mom.Y()) / (mom.X() * mom.X() + mom.Y() * mom.Y());
        auto dx = flightLineXY.X() - k * static_cast<double>(mom.X());
        auto dy = flightLineXY.Y() - k * static_cast<double>(mom.Y());
        auto absImpPar = std::sqrt(dx * dx + dy * dy);
        auto flightLine = posSV - point;
        auto cross = mom.Cross(flightLine);
        return (cross.Z() > 0. ? absImpPar : -1. * absImpPar);
    }

    static double calculateDispersion(const SecondaryVertexFinder::SecVtxInfo &vtxInfo,
                                      const std::vector<Track> &allTracks,
                                      DebugLevel debugLevel = DebugLevel::INFO)
    {
        if (vtxInfo.indices.empty())
        {
            log(debugLevel, DebugLevel::DEBUG, "No indices in vtxInfo");
            return 0.0;
        }

        // Validate all indices before processing
        for (const auto idx : vtxInfo.indices)
        {
            if (idx >= allTracks.size())
            {
                log(debugLevel, DebugLevel::ERROR,
                    "Invalid track index found: " + std::to_string(idx) +
                        " (valid range: 0 to " + std::to_string(allTracks.size() - 1) + ")");
                return 0.0;
            }
        }

        std::vector<Track> prongs;
        for (const auto &idx : vtxInfo.indices)
        {
            try
            {
                prongs.push_back(allTracks[idx]);
                log(debugLevel, DebugLevel::VERBOSE,
                    "Added track " + std::to_string(idx));
            }
            catch (const std::exception &e)
            {
                log(debugLevel, DebugLevel::ERROR,
                    "Error copying track: " + std::string(e.what()));
                return 0.0;
            }
        }

        if (prongs.empty())
        {
            log(debugLevel, DebugLevel::DEBUG, "No valid prongs found");
            return 0.0;
        }

        double sumDistanceSquared = 0.0;
        for (const auto &prong : prongs)
        {
            TVector3 diff = prong.getPosition() - vtxInfo.position;
            if (std::isnan(diff.Mag2()))
            {
                log(debugLevel, DebugLevel::WARNING,
                    "Invalid distance calculation in dispersion");
                continue;
            }
            sumDistanceSquared += diff.Mag2();
        }

        return std::sqrt(sumDistanceSquared / prongs.size());
    }
};

class VertexFinder
{
public:
    static TVector3 findPrimaryVertex(const std::vector<Track> &tracks,
                                      double convergenceDist = 0.01, // in cm
                                      int maxIterations = 100,
                                      double temperatureScale = 2.0, // Controls how fast weights change
                                      double chi2Cut = 9.0,
                                      DebugLevel debugLevel = DebugLevel::INFO)
    { // 3 sigma cut in 3D
        log(debugLevel, DebugLevel::DEBUG, "Starting primary vertex finding...");
        log(debugLevel, DebugLevel::DEBUG, "Number of tracks: " + std::to_string(tracks.size()));

        if (tracks.empty())
        {
            log(debugLevel, DebugLevel::WARNING, "No tracks provided for vertex finding");
            return TVector3(0., 0., 0.);
        }

        // Initial vertex position: mean of track DCA points
        TVector3 vtxPos(0., 0., 0.);
        for (const auto &track : tracks)
        {
            vtxPos += track.getPosition();
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
                TVector3 diff = track.getPosition() - vtxPos;

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
                newPos += track.getPosition() * weight;
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
    DCASmearing(DebugLevel level = DebugLevel::INFO) : debugLevel(level) {}
    static void initializeHistograms();
    static void fillHistograms(const Track &truthPart, const Track &smearedPart);
    void setDebugLevel(DebugLevel level)
    {
        debugLevel = level;
    }

private:
    DebugLevel debugLevel;
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
    static TH2F *hDCAxy_vs_eta;
    static TH2F *hDCAz_vs_eta;
};

// Add the declaration for convertToTParticle
TParticle convertToTParticle(const Pythia8::Particle &part);
