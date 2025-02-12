#pragma once

#include "TLorentzVector.h"
#include "TVector3.h"
#include "Pythia8/Pythia.h"

class JParticle {
private:
    int pdgId_;
    int status_;
    int mother1_;
    int mother2_;
    int daughter1_;
    int daughter2_;
    TLorentzVector momentum_;
    TVector3 vertex_;
    double prodTime_;

public:
    // Default constructor
    JParticle() : 
        pdgId_(0), status_(0), mother1_(-1), mother2_(-1),
        daughter1_(-1), daughter2_(-1), prodTime_(0) {}

    // Constructor from Pythia particle
    JParticle(const Pythia8::Particle& part) :
        pdgId_(part.id()),
        status_(part.status()),
        mother1_(part.mother1()),
        mother2_(part.mother2()),
        daughter1_(part.daughter1()),
        daughter2_(part.daughter2()),
        momentum_(part.px(), part.py(), part.pz(), part.e()),
        vertex_(part.xProd(), part.yProd(), part.zProd()),
        prodTime_(part.tProd()) {}

    // Getters
    int pdgId() const { return pdgId_; }
    int status() const { return status_; }
    int mother1() const { return mother1_; }
    int mother2() const { return mother2_; }
    int daughter1() const { return daughter1_; }
    int daughter2() const { return daughter2_; }
    
    // Momentum getters
    double px() const { return momentum_.Px(); }
    double py() const { return momentum_.Py(); }
    double pz() const { return momentum_.Pz(); }
    double e() const { return momentum_.E(); }
    double p() const { return momentum_.P(); }
    double pt() const { return momentum_.Pt(); }
    double eta() const { return momentum_.Eta(); }
    double phi() const { return momentum_.Phi(); }
    double m() const { return momentum_.M(); }
    const TLorentzVector& p4() const { return momentum_; }

    // Vertex getters
    double vx() const { return vertex_.X(); }
    double vy() const { return vertex_.Y(); }
    double vz() const { return vertex_.Z(); }
    const TVector3& vertex() const { return vertex_; }
    double prodTime() const { return prodTime_; }
};

// Helper function to convert Pythia8::Particle to JParticle
inline JParticle convertToJParticle(const Pythia8::Particle& part) {
    return JParticle(part);
} 