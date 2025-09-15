#include "ElectronID.hh"

#include "edm4hep/utils/vector_utils.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4eic/ClusterCollection.h"

#include <iostream>

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

ElectronID::ElectronID() {

	mEe = 10.;
	mEh = 100.;
	std::cout << "!!! ElectronID: You have not specified beam energies...defaulting to 10x100 GeV !!!" << std::endl;

	mEoP_min = 0.9;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.4;
	mIsoE = 0.9;


}

ElectronID::ElectronID(double Ee, double Eh) {

	mEe = Ee;
	mEh = Eh;

	mEoP_min = 0.9;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.4;
	mIsoE = 0.9;

}

ElectronID::~ElectronID() {
}


void ElectronID::SetEvent(const podio::Frame* event) {

	mEvent = event;
	hfs_dpt.clear();
	hfs_dpz.clear();
	hfs_de.clear();
	hfs_theta.clear();
}

edm4eic::ReconstructedParticleCollection ElectronID::FindHadronicFinalState(bool use_mc, int object_id, LorentzRotation boost) {

	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon->setSubsetCollection();

	edm4hep::MCParticleCollection meMiss;
	meMiss->setSubsetCollection();

	auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");

	if ( use_mc )
	{
		edm4hep::MCParticleCollection meMC = GetMCElectron();
		auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");
	
		for(const auto& assoc : RecoMC) 
		{
			if(assoc.getSim() != meMC[0] && assoc.getSim().getGeneratorStatus() == 1) 
			{
				PxPyPzEVector v(assoc.getSim().getMomentum().x, assoc.getSim().getMomentum().y, assoc.getSim().getMomentum().z, assoc.getSim().getEnergy());
				PxPyPzEVector u(assoc.getRec().getMomentum().x, assoc.getRec().getMomentum().y, assoc.getRec().getMomentum().z, assoc.getRec().getEnergy());
				PxPyPzEVector c(assoc.getRec().getMomentum().x, assoc.getRec().getMomentum().y, assoc.getRec().getMomentum().z, GetCalorimeterEnergy(assoc.getRec()));
				
				hfs_dpt.push_back((u.Pt()-v.Pt())/v.Pt());
				hfs_dpz.push_back((u.Z()-v.Z())/v.Z());
				hfs_de.push_back((u.E()-v.E())/v.E());
				hfs_theta.push_back(v.Theta()*(180./M_PI));				

				v = boost(v);
				u = boost(u);
				c = boost(c);

				meRecon.push_back(assoc.getRec());
			}
		}
	}
	else
	{ 
		for(const auto& mcp : rcparts) {
			if ( mcp.getObjectID().index != object_id )
				meRecon.push_back(mcp);
		}
	}

	return meRecon;
}

edm4hep::MCParticleCollection ElectronID::GetMCHadronicFinalState() {

	edm4hep::MCParticleCollection mhMC;
	mhMC->setSubsetCollection();

	auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");

	std::vector<edm4hep::MCParticle> mc_hadronic;
	edm4hep::MCParticleCollection meMC = GetMCElectron();
	
	bool found_scattered_e = false; 
	for(const auto& mcp : mcparts) {
		if (mcp.getGeneratorStatus() == 1 && mcp.getObjectID().index != meMC[0].getObjectID().index ) 
			mhMC.push_back(mcp);	
	}

	return mhMC;
}

edm4eic::ReconstructedParticleCollection ElectronID::FindScatteredElectron() {

	// Get all the edm4eic objects needed for electron ID
	auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
	
	// Create collection for storing scattered electron candidates
	// (subset collection of ReconstructedParticleCollection)
	edm4eic::ReconstructedParticleCollection scatteredElectronCandidates;
	scatteredElectronCandidates->setSubsetCollection();

	// Loop over all ReconstructedParticles for this event
	for (const auto& reconPart : rcparts) {

		// Require negative particle
		if(reconPart.getCharge() >= 0) continue;

		// Require at least one track and one cluster
		if(reconPart.getClusters().size() == 0 || reconPart.getTracks().size() == 0) continue;

		// Calculate rcpart_ member variables for this event
		CalculateParticleValues(reconPart, rcparts);

		// Calculate E/p and isolation fraction for this event
		// Note that the rcpart_ variables are set in CalculateParticleValues
		double recon_EoP = rcpart_sum_cluster_E / edm4hep::utils::magnitude(reconPart.getMomentum());
		double recon_isoE = rcpart_sum_cluster_E / rcpart_isolation_E;
		
		// Apply scattered electron ID cuts
		if(recon_EoP < mEoP_min || recon_EoP > mEoP_max) continue;
		if(recon_isoE < mIsoE) continue;

		// If particle passes cuts, add to output collection
		scatteredElectronCandidates.push_back(reconPart);

	}	

	return scatteredElectronCandidates;

}

edm4hep::MCParticleCollection ElectronID::GetMCElectron() {

	edm4hep::MCParticleCollection meMC;
	meMC->setSubsetCollection();

	auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");

	std::vector<edm4hep::MCParticle> mc_electrons;
	
	for(const auto& mcp : mcparts) {
		if(mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 1) {
			mc_electrons.push_back(mcp);
		}
	}

	// For now, just take first electron
	if(mc_electrons.size() > 0) {
		meMC.push_back(mc_electrons[0]);
	}

	return meMC;
	
}

edm4eic::ReconstructedParticleCollection ElectronID::GetTruthReconElectron() {

	edm4hep::MCParticleCollection meMC = GetMCElectron();
	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon->setSubsetCollection();

	auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

	for(const auto& assoc : RecoMC) {
		if(assoc.getSim() == meMC[0]) {
			meRecon.push_back(assoc.getRec());
			break;
		}
	}

	return meRecon;

}
	


void ElectronID::CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
		const edm4eic::ReconstructedParticleCollection& rcparts) {

	rcpart_sum_cluster_E = 0.;
	rcpart_lead_cluster_E = 0.;
	rcpart_isolation_E = 0.;
	rcpart_deltaH = 0.;

	const edm4eic::Cluster* lead_cluster = nullptr;

	for (const auto& cluster : rcp.getClusters()) {
		rcpart_sum_cluster_E += cluster.getEnergy();
		if(cluster.getEnergy() > rcpart_lead_cluster_E) {
			lead_cluster = &cluster;
			rcpart_lead_cluster_E = cluster.getEnergy();
		}
	}

	if(!lead_cluster) return;

	const auto& lead_pos = lead_cluster->getPosition();
	double lead_eta = edm4hep::utils::eta(lead_pos);
	double lead_phi = edm4hep::utils::angleAzimuthal(lead_pos);

	for (const auto& other_rcp : rcparts) {
		if (&other_rcp == &rcp) continue;  // Skip the same particle

		for (const auto& other_cluster : other_rcp.getClusters()) {

			const auto& other_pos = other_cluster.getPosition();
			double other_eta = edm4hep::utils::eta(other_pos);
			double other_phi = edm4hep::utils::angleAzimuthal(other_pos);

			double d_eta = other_eta - lead_eta;
			double d_phi = other_phi - lead_phi;

			// Adjust d_phi to be in the range (-pi, pi)
			if (d_phi > M_PI) d_phi-=2*M_PI;
			if (d_phi < -M_PI) d_phi+=2*M_PI;

			double dR = std::sqrt(std::pow(d_eta, 2) + std::pow(d_phi, 2));

			// Check if the cluster is within the isolation cone
			if (dR < mIsoR) {
				rcpart_isolation_E += other_cluster.getEnergy();
			}
		}
	}

	return;
}

edm4eic::ReconstructedParticle ElectronID::SelectHighestPT(const edm4eic::ReconstructedParticleCollection& ecandidates) {

	edm4eic::ReconstructedParticle erec;
	double max_pT = 0.;
	
	for(auto ecand : ecandidates) {
		double e_pT = edm4hep::utils::magnitudeTransverse(ecand.getMomentum());
		if(e_pT > max_pT) {
			erec = ecand;
			max_pT = e_pT;
		}
	}

	return erec;

}

double ElectronID::GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp) {

	double sum_cluster_E = 0.;
	for (const auto& cluster : rcp.getClusters()) {
		sum_cluster_E += cluster.getEnergy();
	}
	return sum_cluster_E;

}

void ElectronID::GetBeam(LorentzRotation &boost, TLorentzVector &in_e, TLorentzVector &in_n) 
{ 
    edm4hep::MCParticle mc_electron;
    edm4hep::MCParticle mc_nucleon;

    auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");
    vector<edm4hep::MCParticle> spec_protons;

	for(const auto& mcp : mcparts)
    {
        if ( mcp.getGeneratorStatus() == 4 )
        {
            if ( mcp.getPDG() == ID_ELECTRON )
                mc_electron = mcp;
            else 
                mc_nucleon = mcp;
        }
    }

    if ( !mc_electron.isAvailable() || !mc_nucleon.isAvailable() )
        return;

	in_e.SetPxPyPzE(mc_electron.getMomentum().x, mc_electron.getMomentum().y, mc_electron.getMomentum().z, mc_electron.getEnergy());
    in_n.SetPxPyPzE(mc_nucleon.getMomentum().x, mc_nucleon.getMomentum().y, mc_nucleon.getMomentum().z, mc_nucleon.getEnergy());

    // get boost matrix -- redo every run because the proton / neutron has different mass .. but really this should not change per event .. to be changed
    const PxPyPzEVector ei(
        eicrecon::round_beam_four_momentum(
            mc_electron.getMomentum(),
            mc_electron.getMass(),
            {-1*mEe},
            0.0)
        );

    const PxPyPzEVector pi(
        eicrecon::round_beam_four_momentum(
            mc_nucleon.getMomentum(),
            mc_nucleon.getMass(),
            {mEh},
            -0.025)
        );
    
    boost = eicrecon::determine_boost(ei, pi); // Get boost to colinear frame

	return;
}
