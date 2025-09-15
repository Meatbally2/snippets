// Analysis headers
#include "ElectronID.cc"
#include "kinematics.cc"
#include "draw_helper.h"
#include "reconMethod.hh"
#include "constants.h"

#pragma once

TFile* outFile = NULL;
TTree* outTree = NULL;

const double Me = 0.511e-3; // GeV
const double Mp = 0.938272; // GeV
const double crossing_angle = -0.025; // rad

// QA histograms
TH1F* h_dEc[4];
TH1F* h_dEt[4];
TH1F* h_dt[4];
TH1F* h_dp[4];
TH2F* h_tde[2];
TH2F* h_tdp;
TH2F* h_pde;
TH2F* h_ede;
TH2F* h_pt;
TH2F* h_et;

TH1F* h_dhf[4];
TH1F* h_dpt[4];

TH1F* h_hfs_dpz;
TH1F* h_hfs_dpt;
TH1F* h_hfs_de;
TH2F* h_hfs_dpz_t;
TH2F* h_hfs_dpt_t;
TH2F* h_hfs_de_t;

enum {EL, JB, DA, SIG, ESIG, MC};
enum {reconID, truthID, mcID};

void SetInputBranchAddresses();
void CreateOutputTree(TString outFileName);
void ResetVariables();
void CalculateRealElectronKinematics(TLorentzVector ki, TLorentzVector P, TLorentzVector kf, Kinematics &v);
TLorentzVector GetHadronBeam(double fEh);
void DefineHistograms();
void plot_kinematics_qa();
void plot_energy_and_track();
bool pass_ke_cut(Kinematics v);

using namespace std;

int positive_eID;
vector<Kinematics> algorithm;

float mc_p;
float mc_eta;
float mc_phi;

float track_p;
float track_eta;
float track_phi;