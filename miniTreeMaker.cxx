#include "pleaseIncludeMe.h"
#include <TVector3.h>
#include <TVector.h>
#include <cmath>
#include "TF1.h"
#include "TF2.h"
#include "TSystem.h"
#include "TTreeReaderValue.h"

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

int miniTreeMaker(TString rec_file, TString outputfile, int doMC_=0)
{	

    
TString name_of_input = (TString) rec_file;
std::cout << "Input file = " << name_of_input << endl;
auto tree = new TChain("events");
tree->Add(name_of_input);
TTreeReader tree_reader(tree);       // !the tree reader

//input tree:    
TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
// MC particle pz array for each MC particle
TTreeReaderArray<double> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
TTreeReaderArray<double> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
TTreeReaderArray<double> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};

//Reconstructed EcalEndcapNClusters
TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};

//ReconstructedFarForwardZDCNeutrals 
TTreeReaderArray<float> zdc_x_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.x"};
TTreeReaderArray<float> zdc_y_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.y"};
TTreeReaderArray<float> zdc_z_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.z"};
TTreeReaderArray<float> zdc_energy_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.energy"};

//ForwardRomanPotRecHits
TTreeReaderArray<float> rp_x_array = {tree_reader, "ForwardRomanPotRecHits.position.x"};
TTreeReaderArray<float> rp_y_array = {tree_reader, "ForwardRomanPotRecHits.position.y"};
TTreeReaderArray<float> rp_z_array = {tree_reader, "ForwardRomanPotRecHits.position.z"};

//ForwardOffMTrackerRecHits
TTreeReaderArray<float> omd_x_array = {tree_reader, "ForwardOffMTrackerRecHits.position.x"};
TTreeReaderArray<float> omd_y_array = {tree_reader, "ForwardOffMTrackerRecHits.position.y"};
TTreeReaderArray<float> omd_z_array = {tree_reader, "ForwardOffMTrackerRecHits.position.z"};

//Reconstructed all central + B0 tracks
TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.charge"};

//define output file
TString output_name_dir = outputfile+"_output.root";
cout << "Output file = " << output_name_dir << endl;
TFile* output = new TFile(output_name_dir,"RECREATE");

//define output tree
TTree* outputTree = new TTree("miniTree", "Tree with structured event data");
Event event;
outputTree->Branch("event", &event);

//events
tree_reader.SetEntriesRange(0, tree->GetEntries());
//chain->GetEntries();
while (tree_reader.Next()) {

	int numberOfChargedParticles=0;
	int numberOfMCParticles=0;
	event.mcp.clear();
  event.particles.clear();
  event.clusters_eemc.clear();
  event.clusters_zdc.clear();
  event.hit_rp.clear();
  event.hit_omd.clear();

  if(doMC_){
	  for(int imc=0;imc<mc_px_array.GetSize();imc++){
	  	MCp mc;
	  	mc.status = mc_genStatus_array[imc];
	  	mc.px = mc_px_array[imc];
	  	mc.py = mc_py_array[imc];
	  	mc.pz = mc_pz_array[imc];
	  	mc.mass = mc_mass_array[imc];
	  	event.mcp.push_back(mc);
	  	numberOfMCParticles++;
	  }
	}

  for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
  	Particle p;
  	p.px = reco_px_array[itrk];
  	p.py = reco_py_array[itrk];
  	p.pz = reco_pz_array[itrk];
  	p.charge = reco_charge_array[itrk];
    event.particles.push_back(p);
    numberOfChargedParticles++;
  }

  //EEMC
  for(int iclus=0;iclus<em_energy_array.GetSize();iclus++){
  		Cluster_EEMC cluster;
			cluster.energy=em_energy_array[iclus];
			cluster.x=em_x_array[iclus];
			cluster.y=em_y_array[iclus];
	    event.clusters_eemc.push_back(cluster);
	}

	//RP
	for(int ihit=0;ihit<rp_x_array.GetSize();ihit++){
  		Hit_RP hit;
			hit.x=rp_x_array[ihit];
			hit.y=rp_y_array[ihit];
			hit.z=rp_z_array[ihit];
	    event.hit_rp.push_back(hit);
	}

	//OMD
	for(int ihit=0;ihit<omd_x_array.GetSize();ihit++){
  		Hit_OMD hit;
			hit.x=omd_x_array[ihit];
			hit.y=omd_y_array[ihit];
			hit.z=omd_z_array[ihit];
	    event.hit_omd.push_back(hit);
	}

  event.nParticles=numberOfChargedParticles;
  event.nMCParticles=numberOfMCParticles;
	outputTree->Fill();
}

output->cd();
outputTree->Write();
output->Close();

return 0;
}
