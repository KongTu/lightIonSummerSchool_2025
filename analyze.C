
#include "analyze.h"
double giveme_t_method_L(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn, TLorentzVector vmOut)
{
    //MASS_Ca40 = 37.2 GeV
	TLorentzVector aInVec(pIn.Px()*40,pIn.Py()*40,pIn.Pz()*40,sqrt(pIn.Px()*40*pIn.Px()*40 + pIn.Py()*40*pIn.Py()*40 + pIn.Pz()*40*pIn.Pz()*40 + 37.2*37.2) );
	double method_L = 0;
	TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
	double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
	double p_Aminus = (37.2*37.2 + p_TAsquared) / p_Aplus;
	TLorentzVector a_beam_scattered_corr; 
	a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = -(a_beam_scattered_corr-aInVec).Mag2();
	return method_L;
}
void analyze(){

	//input file
	TFile* file = new TFile("tree_all.root");

	//read tree
	auto tree = (TTree*) file->Get("miniTree");
	Event* event = nullptr;
	tree->SetBranchAddress("event", &event);
	
	//number of events
	int nEvents=tree->GetEntries();
	
	//beam
	TLorentzVector eIn(0,0,-18,18);
	TLorentzVector aIn(-3.438,0,137.457,137.503);

	TFile* output = new TFile("results.root","RECREATE");
	//example hist
	TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{e,MC}",100,0,20);
	TH1D* h_y_e = new TH1D("h_y_e",";y_{e,MC}",100,0,1);
	TH1D* h_Q2REC_e = new TH1D("h_Q2REC_e",";Q^{2}_{e,REC}",100,0,20);
	TH1D* h_yREC_e = new TH1D("h_yREC_e",";y_{e,REC}",100,0,1);
	
	TH1D* h_EoverP_REC = new TH1D("h_EoverP_REC",";E/p",200,0,2);
	TH1D* h_Epz_REC = new TH1D("h_Epz_REC",";E - p_{z} (GeV)",200,0,50);

	TH1D* h_massMC = new TH1D("h_massMC",";mass",100,0.9,1.2);
	TH1D* h_massREC = new TH1D("h_massREC",";mass",100,0.9,1.2);

	TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} (GeV)",100,0,20);	
	TH1D* h_energy_REC = new TH1D("h_energy_REC",";E_{REC} (GeV)",100,0,20);

	TH1D* h_tMC = new TH1D("h_tMC",";t_{MC} (GeV^{2})",100,0,0.3);
	TH1D* h_tREC = new TH1D("h_tREC",";t_{REC} (GeV^{2})",100,0,0.3);

	//event loop
	for(int ievt=0;ievt<nEvents;ievt++){
		tree->GetEntry(ievt);

		if(event->clusters_zdc.size()>0) continue;
		if(event->hit_rp.size()>0) continue;
		if(event->hit_omd.size()>0) continue;
		
		//mc particle loop
		TLorentzVector escat(0,0,0,0);
		TLorentzVector Kaon(0,0,0,0);
		vector<TLorentzVector> kdaug;
		int incoherent=0;
		for(int i=0;i<event->mcp.size();i++){
			const MCp& p = event->mcp[i];
			int status=p.status;
			float mass=p.mass;
			TVector3 part(p.px,p.py,p.pz);
			double energy=sqrt(p.px*p.px+p.py*p.py+p.pz*p.pz+mass*mass);
			if(status==1){
				//select electrons
				if(mass<0.0006 && mass>0.0005){
					escat.SetPxPyPzE(p.px,p.py,p.pz,energy);
				}
				//select kaons
				if(mass>0.49&&mass<0.5){
					Kaon.SetPxPyPzE(p.px,p.py,p.pz,energy);
					kdaug.push_back(Kaon);
				}

				if(part.Eta()>3.5&&(p.mass>0.9383||p.mass<0.938)){
					incoherent++;
				}
			}
		}

		if(incoherent) continue;

	   	/*phase space*/	
		TLorentzVector qbeam=eIn-escat;
		double Q2=-(qbeam).Mag2();  
		double pq=aIn.Dot(qbeam);
		double y= pq/aIn.Dot(eIn);

	    //MC level phase space cut
		if(Q2<2.||Q2>10.) continue;
		if(y<0.01||y>0.85) continue;
		h_Q2_e->Fill(Q2);
		h_y_e->Fill(y);
		h_energy_MC->Fill( escat.E() );
		/*end phase space*/

		if(kdaug.size()>=2){
			if(fabs(kdaug[0].Eta())<1.5&&fabs(kdaug[1].Eta())<1.5){
				double mass = (kdaug[0]+kdaug[1]).M();
				h_massMC->Fill(mass);
				double t_MC = -(eIn-escat-(kdaug[0]+kdaug[1])).Mag2();
				if( fabs(mass-1.02)<0.02&& fabs((kdaug[0]+kdaug[1]).Rapidity())<1.5){
					double t=giveme_t_method_L(eIn,escat,aIn,(kdaug[0]+kdaug[1]));
					h_tMC->Fill(t);
				}
			}
		}

		//reco events:
		TLorentzVector escat_REC(0,0,0,0);
		TLorentzVector K1_REC(0,0,0,0);
		TLorentzVector K2_REC(0,0,0,0);

		double e_energy=0;
		//reconstructed EEMC loop
		for(int i=0;i<event->clusters_eemc.size();i++){
			const Cluster_EEMC& clus = event->clusters_eemc[i];
			double radius=sqrt(clus.x*clus.x+clus.y*clus.y);
			if( radius>550. ) continue;
			e_energy=clus.energy;
		}
		h_energy_REC->Fill( e_energy );

		double maxpz=0;
		TLorentzVector elec_trk(0,0,0,0);
		TLorentzVector hfs_e(0,0,0,0);
		TLorentzVector part(0,0,0,0);
		int incoherent_rec=0;
		//reconstructed charged particle loop
		for(int i=0;i<event->particles.size();i++){
			const Particle& p = event->particles[i];
			TVector3 trk(p.px,p.py,p.pz);
			if(fabs(p.pz)>maxpz && trk.Eta()<-1.5) {
				maxpz=p.pz;
				double pmom = sqrt(e_energy*e_energy- MASS_ELECTRON*MASS_ELECTRON );
				double eta=trk.Eta();
				double phi=trk.Phi();
				double pt = TMath::Sin(trk.Theta())*pmom;
				escat_REC.SetPtEtaPhiM(pt,eta,phi,MASS_ELECTRON);
				elec_trk.SetVectM(trk,MASS_ELECTRON);
			}

			// daughter of vm.
			if(fabs(trk.Eta())<1.5){
				if(p.charge==-1) K1_REC.SetVectM(trk,MASS_KAON);
				if(p.charge==+1) K2_REC.SetVectM(trk,MASS_KAON);
			}
			//hfs + scattered electron
			if(fabs(trk.Eta())<5) {
				part.SetVectM(trk,MASS_PION);
				hfs_e += part;
			}
			//incoherent
			if(trk.Eta()>3.5) incoherent_rec++;
		}

		//selections
		if(incoherent_rec) continue;

		//E over p
		double EoverP=escat_REC.E()/elec_trk.P();
		h_EoverP_REC->Fill( EoverP );

		//E - pz
		double EpzREC= (hfs_e-elec_trk+escat_REC).E() - (hfs_e-elec_trk+escat_REC).Pz();
		h_Epz_REC->Fill( EpzREC );

		TLorentzVector pbeam(0,0,137.5,137.5);
		TLorentzVector qbeamREC=eIn-escat_REC;
		double Q2REC=-(qbeamREC).Mag2();  
		double pqREC=pbeam.Dot(qbeamREC);
		double yREC= pqREC/pbeam.Dot(eIn);
		// 	double sigma=(hfs_e-elec_trk).E()-(hfs_e-elec_trk).Pz();
		// 	double cosTheta=TMath::Cos(escat_REC.Theta());
		// double yREC_sigma = sigma/(sigma+escat_REC.E()*(1-cosTheta));
		// yREC=yREC_sigma;

		//Event selection:
		if( EpzREC<27||EpzREC>40 ) continue;
		if( EoverP<0.8||EoverP>1.16 ) continue;		

		h_Q2REC_e->Fill(Q2REC);
		h_yREC_e->Fill(yREC);

		//MC level phase space cut
		if(Q2REC<2.||Q2REC>10.) continue;
		if(yREC<0.01||yREC>0.85) continue;
		
		double mass_REC=(K1_REC+K2_REC).M();
		h_massREC->Fill((K1_REC+K2_REC).M());
		double t_REC = -(eIn-escat_REC-(K1_REC+K2_REC)).Mag2();
		if( fabs(mass_REC-1.02)<0.02&& fabs((K1_REC+K2_REC).Rapidity())<1.5)
		{
			double t_L=giveme_t_method_L(eIn,escat_REC,aIn,(K1_REC+K2_REC));
			h_tREC->Fill(t_L);
		}

		// //reconstructed ZDC loop
		// for(int i=0;i<event->clusters_zdc.size();i++){
		// 	const Cluster_ZDC& clus = event->clusters_zdc[i];
		// }

		// //reconstructed RP hit loop
		// for(int i=0;i<event->hit_rp.size();i++){
		// 	const Hit_RP& hit = event->hit_rp[i];
		// }

		// //reconstructed OMD hit loop
		// for(int i=0;i<event->hit_omd.size();i++){
		// 	const Hit_OMD& hit = event->hit_omd[i];
		// }

	}

	output->Write();
	output->Close();

	

}