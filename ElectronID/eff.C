
#include "../GlobalUtil/ePICStyle.c"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"

const std::string group[3] = {"1to10", "10to100", "100to1000"};
const double total_lumi = 8.65*3; // fb^-1

double cs[3][2] = {{0.198440424611563, 0.205327493968226}, {4.04371412707044E-02, 4.41976212963417E-02}, {1.36416909784756E-03, 1.69583242740138E-03}};
double ev[3][2] = {{333675, 666325}, {333694, 666306}, {333365, 666640}};

enum { NO_FOUND, FOUND_E, FOUND_PI, FOUND_OTHERS };

void eff()
{
    SetePICStyle();

    TH2F* h_xq2_all = BookTH2(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_acp = BookTH2(Form("h_xq2_acp"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_eID = BookTH2(Form("h_xq2_eID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_piID = BookTH2(Form("h_xq2_piID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    for ( int i = 0; i < 3; i ++ )
    {
        double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));
        double n_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43));
        double p_lumi = ev[i][1]/(cs[i][1]*(1e-34/1e-43));

        TH2F* h_tmp_all = BookTH2(Form("h_tmp_all_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_acp = BookTH2(Form("h_tmp_acp_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_eID = BookTH2(Form("h_tmp_eID_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_piID = BookTH2(Form("h_tmp_piID_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

        TFile* file = new TFile(Form("../data/eID/lowerEoPcut/10x166_%s_eIDrecon.root", group[i].c_str()));

        TTreeReader reader("T_eID", file);
        TTreeReaderValue<int> eID_status(reader, "eID_status");
        TTreeReaderValue<double> mc_xB(reader, "mc_xB");
        TTreeReaderValue<double> mc_Q2(reader, "mc_Q2");
        TTreeReaderValue<double> mc_y(reader, "mc_y");
        TTreeReaderValue<double> mc_W2(reader, "mc_W2");
        TTreeReaderValue<double> mc_nu(reader, "mc_nu");

        Long64_t nEntries = reader.GetEntries();

        for( size_t ev = 0; ev < nEntries; ev++ ) 
        {
            reader.Next();
            // if(ev%100==0) 
            // cout << "Analysing file " << i << " event " << ev << "/" << nEntries << "\t\r" << std::flush;

            int status = *eID_status;
            float xB = *mc_xB;
            float Q2 = *mc_Q2;
            float y = *mc_y;
            float W2 = *mc_W2;
            float nu = *mc_nu;

            h_tmp_all->Fill(xB, Q2);
            if ( status != NO_FOUND )
                h_tmp_acp->Fill(xB, Q2);
            if ( status == FOUND_E )
                h_tmp_eID->Fill(xB, Q2);
            if ( status == FOUND_PI )
                h_tmp_piID->Fill(xB, Q2);
        }

        h_tmp_all->Scale(total_lumi/gen_lumi);
        h_tmp_acp->Scale(total_lumi/gen_lumi);
        h_tmp_eID->Scale(total_lumi/gen_lumi);
        h_tmp_piID->Scale(total_lumi/gen_lumi);

        h_xq2_all->Add(h_tmp_all);
        h_xq2_acp->Add(h_tmp_acp);
        h_xq2_eID->Add(h_tmp_eID);
        h_xq2_piID->Add(h_tmp_piID);

        file->Close();
    }

    set_2d_scale(h_xq2_all);
    TCanvas* c_xq2_all = draw_2d_standard(h_xq2_all, "c_xq2_all", "all events", 700, 600, true, true);
    TCanvas* c_xq2_acp = draw_2d_standard(h_xq2_acp, "c_xq2_acp", "acp events", 700, 600, true, true);
    TCanvas* c_xq2_eID = draw_2d_standard(h_xq2_eID, "c_xq2_eID", "eID events", 700, 600, true, true);
    TCanvas* c_xq2_piID = draw_2d_standard(h_xq2_piID, "c_xq2_piID", "piID events", 700, 600, true, true);

    TH2F* h_xq2_acp_copy = (TH2F*)h_xq2_acp->Clone();
    process_eff_hist(h_xq2_acp_copy, h_xq2_all);
    TCanvas* c_xq2_acp_eff = draw_2d_efficiency(h_xq2_acp_copy, "c_xq2_acp_eff", "xq2 acp eff", 1400, 600, false, true);

    TH2F* h_xq2_eID_copy = (TH2F*)h_xq2_eID->Clone();
    process_eff_hist(h_xq2_eID_copy, h_xq2_acp);
    TCanvas* c_xq2_eID_eff = draw_2d_efficiency(h_xq2_eID_copy, "c_xq2_eID_eff", "xq2 eID eff", 1400, 600, false, true);

    TH2F* h_xq2_piID_copy = (TH2F*)h_xq2_piID->Clone();
    process_eff_hist(h_xq2_piID_copy, h_xq2_acp);
    TCanvas* c_xq2_piID_eff = draw_2d_efficiency(h_xq2_piID_copy, "c_xq2_piID_eff", "xq2 piID eff", 1400, 600, false, true);

    return;
}