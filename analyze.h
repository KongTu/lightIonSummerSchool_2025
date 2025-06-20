#include "RiceStyle.h"
using namespace std;
#define MASS_KAON     0.493667
#define MASS_PION     0.13957
#define MASS_ELECTRON 0.00051
#define MASS_PROTON   0.93826

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

double giveme_t_method_E(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn, TLorentzVector vmOut)
{
    double method_E = -(eIn-eOut-vmOut).Mag2();
    return method_E;
}

double giveme_t_method_L(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn, TLorentzVector vmOut)
{
    TLorentzVector aInVec(pIn.Px(),pIn.Py(),pIn.Pz(),sqrt(pIn.Px()*pIn.Px() + pIn.Py()*pIn.Py() + pIn.Pz()*pIn.Pz() + MASS_PROTON*MASS_PROTON) );
    double method_L = 0;
    TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
    double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
    double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
    double p_Aminus = (MASS_PROTON*MASS_PROTON + p_TAsquared) / p_Aplus;
    TLorentzVector a_beam_scattered_corr; 
    a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
    method_L = -(a_beam_scattered_corr-aInVec).Mag2();
    return method_L;
}