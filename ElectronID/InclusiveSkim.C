// Data model headers
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/utils/vector_utils.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "podio/Frame.h"
// #include "podio/ROOTFrameReader.h" -- for simulation campaign before March 2023 with eic-shell 25.03.2
#include "podio/ROOTReader.h"

// ROOT headers
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// Analysis headers
#include "InclusiveSkim.h"
#include "drawKinematics.C"

void InclusiveSkim() {

	// Set beam energy
	double Ee = 10.;
	double Eh = 100.;

	// Set what electron ID to use
	int eID_type = truthID;

	// access local file
	// vector<std::string> inFiles = {"pythia8NCDIS_10x100_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0001.eicrecon.tree.edm4eic.root"};

	// access remote file
	vector<std::string> inFiles = {"root://dtn-rucio.jlab.org:1094//volatile/eic/EPIC/RECO/25.05.0/epic_craterlake/DIS/NC/10x100/minQ2=1/pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0000.eicrecon.edm4eic.root"};

	// auto reader = podio::ROOTFrameReader();
	auto reader = podio::ROOTReader();
	reader.openFiles(inFiles);

	// setup eid
	ElectronID* eFinder = new ElectronID(Ee, Eh);

	// setup kinematic algorithms
	algorithm.push_back(Kinematics("E_Method", 2)); // electron method
    algorithm.push_back(Kinematics("JB_Method", 4)); // Jacquet-Blondel method
    algorithm.push_back(Kinematics("DA_Method", 1)); // Double-Angle method
    algorithm.push_back(Kinematics("Sig_Method", kOrange+1)); // Sigma method
    algorithm.push_back(Kinematics("ESig_Method", kGreen+3)); // E-Sigma method
	algorithm.push_back(Kinematics("MC", kGray)); // MC info

	// setup output file and histograms
	TString outFileName = Form("inclusive_skim_%.0fx%.0fGeV.root", Ee, Eh);
	CreateOutputTree(outFileName);
	DefineHistograms();

	// loop events
	for(size_t ev = 0; ev < reader.getEntries("events"); ev++) {

		if(ev%1000==0) 
			cout << "Event " << ev << "/" << reader.getEntries("events") << endl;

		ResetVariables();

		const auto event = podio::Frame(reader.readNextEntry("events"));
		eFinder->SetEvent(&event);

		edm4hep::MCParticleCollection e_mc = eFinder->GetMCElectron();

		// Skip events without scattered MC electron
		if(e_mc.size() == 0) 
			continue;

		// Set boost and get incoming particles
		TLorentzVector in_e;
		TLorentzVector in_n;
		LorentzRotation boost;
        eFinder->GetBeam(boost, in_e, in_n);
		in_e = boost(in_e);
		in_n =boost(in_n);

		// Calculate kinematic variables using MC electron
		TLorentzVector kprime;
    	kprime.SetXYZM(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, Me);
    	TLorentzVector true_e_lab = kprime;
		kprime = boost(kprime); // boost back to head-on frame

		// Get momentum vector elements from MC electron
		mc_p = edm4hep::utils::magnitude(e_mc[0].getMomentum());
		mc_eta = edm4hep::utils::eta(e_mc[0].getMomentum());
		mc_phi = edm4hep::utils::angleAzimuthal(e_mc[0].getMomentum());

		// Don't have to boost if we are just calulcating x, Q2, y using e method because they are Lorentz invariant
		// But since we boosted the outgoing ones earlier for later analysis, we need to boost the incoming ones as well
    	CalculateRealElectronKinematics(in_e, in_n, kprime, algorithm[MC]);

		if ( !pass_ke_cut(algorithm[MC]) ) 
            continue;

		// Set Range for plots
	    int Qrange = 0;
		if ( algorithm[MC].Q2 > 1000 )
			Qrange = 3;
		else if ( algorithm[MC].Q2 > 100 )
			Qrange = 2;
		else if ( algorithm[MC].Q2 > 10 )
			Qrange = 1;

		// Find scattered electrons using mc association
		auto e_truth = eFinder->GetTruthReconElectron();

		// Find scattered electrons using detector information
		auto e_candidates = eFinder->FindScatteredElectron();
		edm4eic::ReconstructedParticle e_rec;	

		// If there are multiple candidates, select one with highest pT
		if(e_candidates.size() > 0) {
			positive_eID = 1;			
			e_rec = eFinder->SelectHighestPT(e_candidates);
		}

		// select which eID we are using
		PxPyPzEVector scat_ele(0, 0, 0, -1);
		if ( eID_type == truthID ) {
			if(e_truth.size() < 1) 
				continue;
			
			scat_ele.SetPxPyPzE(e_truth[0].getMomentum().x, e_truth[0].getMomentum().y, e_truth[0].getMomentum().z, e_truth[0].getEnergy());
		}
		else if ( eID_type == mcID ) {
			scat_ele = true_e_lab;
		}
		else {
			if(!positive_eID)
				continue;

			scat_ele.SetPxPyPzE(e_rec.getMomentum().x, e_rec.getMomentum().y, e_rec.getMomentum().z, e_rec.getEnergy());
		}

		PxPyPzEVector rec_e_lab = scat_ele;
		if( eID_type!=2 ) 
		{
			// Calculate kinematics using track info only
			kprime.SetXYZM(rec_e_lab.Px(), rec_e_lab.Py(), rec_e_lab.Pz(), Me);
			TLorentzVector track_kprime = kprime;
			h_dEt[Qrange]->Fill((track_kprime.E()-true_e_lab.E())/true_e_lab.E());
			h_dt[Qrange]->Fill((track_kprime.Theta()-true_e_lab.Theta())/true_e_lab.Theta());
			h_dp[Qrange]->Fill((track_kprime.Phi()-true_e_lab.Phi())/true_e_lab.Phi());
			h_tde[1]->Fill(true_e_lab.Theta()*(180./M_PI), (track_kprime.E()-true_e_lab.E())/true_e_lab.E());
			h_tdp->Fill(true_e_lab.Theta()*(180./M_PI), (track_kprime.Phi()-true_e_lab.Phi())/true_e_lab.Phi());
			h_pde->Fill(true_e_lab.Pt(), (track_kprime.E()-true_e_lab.E())/true_e_lab.E());
			h_pt->Fill(true_e_lab.Theta()*(180./M_PI),true_e_lab.Pt());

			// Calculate kinematics using track + cluster energy
			TVector3 e3v = kprime.Vect();

			double track_clust_E;
			if ( eID_type == truthID )
				track_clust_E = eFinder->GetCalorimeterEnergy(e_truth[0]);
			else if ( eID_type == mcID )  
				track_clust_E = e_mc[0].getEnergy();
			else
				track_clust_E = eFinder->GetCalorimeterEnergy(e_rec);

			e3v.SetMag(track_clust_E);
			kprime.SetVectM(e3v, Me);
			h_dEc[Qrange]->Fill((track_clust_E-true_e_lab.E())/true_e_lab.E());
			h_tde[0]->Fill(true_e_lab.Theta()*(180./M_PI), (track_clust_E-true_e_lab.E())/true_e_lab.E());
			h_ede->Fill(true_e_lab.E(), (track_clust_E-true_e_lab.E())/true_e_lab.E());
			h_et->Fill(true_e_lab.Theta()*(180./M_PI),true_e_lab.E());
	
			TLorentzVector calo_kprime = kprime;

			// // cut example
			// if ( algorithm[MC].Q2 > 100 )
			// 	scat_ele = calo_kprime;
		}

		scat_ele = boost(scat_ele);
		float theta = scat_ele.Theta();
		float E = scat_ele.E();	

		// Find hadronic final state
		double pxsum = 0;
		double pysum = 0;
		double pzsum = 0;
		double Esum  = 0;

		// MC and Rec data structure are different.. not sure how to not repeat code here
		auto mc_hfsCollection = eFinder->GetMCHadronicFinalState();
		for (const auto p : mc_hfsCollection) {
			// Lorentz vector in lab frame
			PxPyPzEVector hf(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
			// Boost to colinear frame
			hf = boost(hf);

			pxsum += hf.Px();
			pysum += hf.Py();
			pzsum += hf.Pz();
			Esum += hf.E();
		}
		double true_sigma_h = Esum - pzsum;
		double true_pt_had = sqrt(pxsum*pxsum + pysum*pysum);

		pxsum = 0;
		pysum = 0;
		pzsum = 0;
		Esum  = 0;

		auto hfsCollection = eFinder->FindHadronicFinalState(eID_type, e_rec.getObjectID().index, boost);
		for (const auto p : hfsCollection) {
			// Lorentz vector in lab frame
			// double track_clust_E = eFinder->GetCalorimeterEnergy(p); 
			// cout << "energy " << track_clust_E << " " << p.getEnergy() << " " << track_clust_E - p.getEnergy() << endl;
			PxPyPzEVector hf(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
			// PxPyPzEVector hf(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, track_clust_E);
			// Boost to colinear frame
			hf = boost(hf);

			pxsum += hf.Px();
			pysum += hf.Py();
			pzsum += hf.Pz();
			Esum += hf.E();
		}

		double sigma_h = Esum - pzsum;
		double pt_had = sqrt(pxsum*pxsum + pysum*pysum);

		for ( int i = 0; i < eFinder->hfs_dpt.size(); i++ )
		{
			h_hfs_dpt->Fill(eFinder->hfs_dpt[i]);
			h_hfs_dpt_t->Fill(eFinder->hfs_theta[i], eFinder->hfs_dpt[i]);
		}

		for ( int i = 0; i < eFinder->hfs_dpz.size(); i++ )
		{
			h_hfs_dpz->Fill(eFinder->hfs_dpz[i]);
			h_hfs_dpt_t->Fill(eFinder->hfs_theta[i], eFinder->hfs_dpz[i]);
		}

		for ( int i = 0; i < eFinder->hfs_de.size(); i++ )
		{
			h_hfs_de->Fill(eFinder->hfs_de[i]);
			h_hfs_de_t->Fill(eFinder->hfs_theta[i], eFinder->hfs_de[i]);
		}

		if ( eID_type == mcID )
		{
			sigma_h = true_sigma_h;
			pt_had = true_pt_had;
		}

		// Calculate kinematics
		
		// example cut
		// if ( abs((scat_ele.E() - true_e_lab.E())/true_e_lab.E()) > 0.2 )
		//    return;

		if ( E > 0 )
			algorithm[EL].save_variables(calc_elec_method(E, theta, pt_had, sigma_h, Ee, Eh), Mp);

		if ( hfsCollection.size() != 0 ) 
		{
			h_dhf[Qrange]->Fill((sigma_h-true_sigma_h)/true_sigma_h);
			h_dpt[Qrange]->Fill((pt_had-true_pt_had)/true_pt_had);

			algorithm[JB].save_variables(calc_jb_method(E, theta, pt_had, sigma_h, Ee, Eh), Mp);

			if ( E > 0 )
			{
				algorithm[DA].save_variables(calc_da_method(E, theta, pt_had, sigma_h, Ee, Eh), Mp);
				algorithm[SIG].save_variables(calc_sig_method(E, theta, pt_had, sigma_h, Ee, Eh), Mp);
				algorithm[ESIG].save_variables(calc_esig_method(E, theta, pt_had, sigma_h, Ee, Eh), Mp);
			}
		}

		algorithm[MC].Process_and_QA(algorithm[MC].xB, algorithm[MC].y, algorithm[MC].Q2);
		for ( int i = 0; i < algorithm.size()-1; i ++ )
			algorithm[i].Process_and_QA(algorithm[MC].xB, algorithm[MC].y, algorithm[MC].Q2);
		

		outTree->Fill();
	}

	outFile->cd();
	outTree->Write(outTree->GetName(), 2);

	plot_energy_and_track();
	plot_kinematics_qa();

	// outFile->Close();

	return;
}

void CreateOutputTree(TString outFileName) {

	outFile = new TFile(outFileName, "RECREATE");
	outTree = new TTree("T", "Reconstructed and generated information from EICRecon");

	outTree->Branch("positive_eID",		&positive_eID);

	// Create branches for each algorithm
	for ( int i = 0; i < algorithm.size(); i ++ )
	{
		outTree->Branch(Form("%s_xB", algorithm[i].method_name.c_str()), &algorithm[i].xB);
		outTree->Branch(Form("%s_Q2", algorithm[i].method_name.c_str()), &algorithm[i].Q2);
		outTree->Branch(Form("%s_y", algorithm[i].method_name.c_str()), &algorithm[i].y);
		outTree->Branch(Form("%s_nu", algorithm[i].method_name.c_str()), &algorithm[i].nu);
		outTree->Branch(Form("%s_W2", algorithm[i].method_name.c_str()), &algorithm[i].W2);
	}

	outTree->Branch("mc_p",			&mc_p);
	outTree->Branch("mc_eta",		&mc_eta);
	outTree->Branch("mc_phi",		&mc_phi);

	outTree->Branch("track_p",		&track_p);	
	outTree->Branch("track_eta",	&track_eta);	
	outTree->Branch("track_phi",	&track_phi);
}

void ResetVariables() {

	positive_eID = 0;

	for ( auto &v : algorithm )
        v.reset();

	mc_p = -999;
	mc_eta = -999;
	mc_phi = -999;
	
	track_p = -999;
	track_eta = -999;
	track_phi = -999;
	      
}

// mc kinematics for comparison - use real incident particle kinematics
void CalculateRealElectronKinematics(TLorentzVector ki, TLorentzVector P, TLorentzVector kf, Kinematics &v) {
  // input is boosted
  TLorentzVector q = ki - kf;
  v.Q2 = -(q.Dot(q));
  v.nu = (q.Dot(P))/Mp;
  v.xB = v.Q2/(2.*Mp*v.nu);
  v.y  = (q.Dot(P))/(ki.Dot(P));
  v.W2  = (Mp*Mp + (2.*Mp*v.nu) - v.Q2);		
}

TLorentzVector GetHadronBeam(double fEh) {
 
	TLorentzVector hadron_beam;
	hadron_beam.SetX(fEh*sin(crossing_angle));
	hadron_beam.SetY(0.);
	hadron_beam.SetZ(fEh*cos(crossing_angle));
	hadron_beam.SetE(std::hypot(fEh, Mp));
	return hadron_beam;

}

bool pass_ke_cut(Kinematics v)
{
    if ( v.y <= 0.01 || v.y >= 0.95 )
        return false;

    if ( v.Q2 < 2 )
        return false;

    if ( v.W2 < 4 )
        return false;

    return true;
}