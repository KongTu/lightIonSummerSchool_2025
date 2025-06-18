#include "RiceStyle.h"
using namespace std;
#define MASS_KAON     0.493667
#define MASS_PION     0.13957
#define MASS_ELECTRON 0.00051

struct MCp {
    float status, px, py, pz, mass;
};

struct Particle {
    float px, py, pz, charge;
};

struct Cluster_EEMC {
    float x, y, energy;
};

struct Cluster_ZDC {
    float x, y, z, energy;
};

struct Hit_RP {
    float x, y, z;
};

struct Hit_OMD {
    float x, y, z;
};

struct Event {
    Int_t nParticles, nMCParticles;
    std::vector<MCp> mcp;
    std::vector<Particle> particles;
    std::vector<Cluster_EEMC> clusters_eemc;
    std::vector<Cluster_ZDC> clusters_zdc;
    std::vector<Hit_RP> hit_rp;
    std::vector<Hit_OMD> hit_omd;
};