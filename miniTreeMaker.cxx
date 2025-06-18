#include "pleaseIncludeMe.h"
#include <TVector3.h>
#include <TVector.h>
#include <cmath>
#include "TF1.h"
#include "TF2.h"
#include "TSystem.h"

struct Particle {
    Int_t pid;
    Float_t px, py, pz, energy;
};

struct Cluster {
    Float_t eta, phi, energy;
};

struct Event {
    Int_t event_id;
    Float_t Q2_e, x_e, y_e;
    std::vector<Particle> particles;
    std::vector<Cluster> clusters;
};

int miniTreeMaker(TString rec_file, TString outputfile)
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

//InclusiveKinematicsElectron
TTreeReaderValue<float_t> inclusive_Q2_value = {tree_reader, "InclusiveKinematicsElectron.Q2"};
TTreeReaderValue<float_t> inclusive_x_value = {tree_reader, "InclusiveKinematicsElectron.x"};
TTreeReaderValue<float_t> inclusive_y_value = {tree_reader, "InclusiveKinematicsElectron.y"};

//Reconstructed EcalEndcapNClusters
TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};
TTreeReaderArray<float> emhits_x_array = {tree_reader, "EcalEndcapNRecHits.position.x"};
TTreeReaderArray<float> emhits_y_array = {tree_reader, "EcalEndcapNRecHits.position.y"};
TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};

// Reconstructed Truth Seeded Charged particles has all central + B0 tracks
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

	event.Q2_e = inclusive_Q2_value;
  event.x_e = inclusive_x_value;
  event.y_e = inclusive_y_value;
  event.particles.clear();
  event.clusters.clear();
	
	outputTree->Fill();
}

output->cd();
outputTree->Write();
output->Close();

return 0;
}
