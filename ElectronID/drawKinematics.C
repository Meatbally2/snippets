#include "InclusiveSkim.h"

void DefineHistograms() {

	// assume to split to 4 Q2 regions, should make this more flexible
	std::string QrangeName[4] = {"lowQ", "midQ", "highQ", "ultraQ"};
    std::string Qrange[4] = {"Q^{2} #leq 10 (GeV/c^{2})^{2}", "10 (GeV/c^{2})^{2} < Q^{2} #leq 100 (GeV/c^{2})^{2}", "100 (GeV/c^{2})^{2} < Q^{2} #leq 1000 (GeV/c^{2})^{2}", "Q^{2} > 1000 (GeV/c^{2})^{2}"};
    for ( int i = 0; i < 4; i ++ )
    {
        h_dEc[i] = new TH1F(Form("H_%s_dE_cal",  QrangeName[i].c_str()), Form("%s;(Ee_{reco}-Ee_{true})/Ee_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        h_dEt[i] = new TH1F(Form("H_%s_dE_track",  QrangeName[i].c_str()), Form("%s;(Ee_{reco}-Ee_{true})/Ee_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        h_dt[i] = new TH1F(Form("H_%s_dtheta",  QrangeName[i].c_str()), Form("%s;(e_{reco}-e_{true})/e_{true};Counts", Qrange[i].c_str()), 40, -0.01, 0.01);
        h_dp[i] = new TH1F(Form("H_%s_dphi",  QrangeName[i].c_str()), Form("%s;(e_{reco}-e_{true})/e_{true};Counts", Qrange[i].c_str()), 40, -0.01, 0.01);

        h_dhf[i] = new TH1F(Form("H_%s_dhf",  QrangeName[i].c_str()), Form("%s;(#delta_{h, reco}-#delta_{h, true})/#delta_{h, true};Counts", Qrange[i].c_str()), 100, -2, 2);
        h_dpt[i] = new TH1F(Form("H_%s_dpt",  QrangeName[i].c_str()), Form("%s;(p_{t, reco}-p_{t, true})/p_{t, true};Counts", Qrange[i].c_str()), 100, -2, 2);
    }

    std::string etype[2] = {"cal", "track"};
    for ( int i = 0; i < 2; i ++ )
    {
        h_tde[i] = new TH2F(Form("H_dEvT_%s", etype[i].c_str()), Form("dE_{%s} vs #theta_{mc};#theta_{e};dE/E", etype[i].c_str()), 180, 0, 180, 100,  -1, 1);
    }
    h_tdp = new TH2F(Form("H_dPvT"), Form("d#phi vs #theta_{mc};#theta_{e};d#phi/#phi"), 180, 0, 180, 40, -0.01, 0.01);
    h_pde = new TH2F(Form("H_dEvP_%s", etype[1].c_str()), Form("dE_{%s} vs p_{mc};p_{t,e};dE/E", etype[1].c_str()), 180, 0, 60, 100,  -0.2, 0.2);
    h_ede = new TH2F(Form("H_dEvE_%s", etype[0].c_str()), Form("dE_{%s} vs E_{mc};E_{e};dE/E", etype[0].c_str()), 180, 0, 60, 100,  -0.2, 0.2);
    h_pt = new TH2F(Form("H_PvT_%s", etype[1].c_str()), Form("p_{t,%s} vs #theta_{mc};#theta_{e};p_{t}", etype[1].c_str()), 180, 0, 180, 240, 0, 80);
    h_et = new TH2F(Form("H_EvT_%s", etype[0].c_str()), Form("E_{%s} vs #theta_{mc};#theta_{e};E", etype[0].c_str()), 180, 0, 180, 240, 0, 80);

    h_hfs_dpz = new TH1F("H_HFS_dpz", Form(";(p_{z, reco}-p_{z, true})/p_{z, true};Counts"), 100, -1, 1);
    h_hfs_dpt = new TH1F("H_HFS_dpt", Form(";(p_{t, reco}-p_{t, true})/p_{t, true};Counts"), 100, -1, 1);
    h_hfs_de = new TH1F("H_HFS_de", Form(";(E_{reco}-E_{true})/E_{true};Counts"), 100, -1, 1);
    h_hfs_dpz_t = new TH2F("H_HFS_dpz_vT", Form(";#theta_{mc};(p_{z, reco}-p_{z, true})/p_{z, true}"), 180, 0, 180, 100, -1, 1);
    h_hfs_dpt_t = new TH2F("H_HFS_dpt_vT", Form(";#theta_{mc};(p_{t, reco}-p_{t, true})/p_{t, true}"), 180, 0, 180, 100, -1, 1);
    h_hfs_de_t = new TH2F("H_HFS_de_vT", Form(";#theta_{mc};(E_{reco}-E_{true})/E_{true}"), 180, 0, 180, 100, -1, 1);
}

void plot_kinematics_qa()
{
    TCanvas* c_eRecon_eff[algorithm.size()-1];

    TCanvas* c_qa_plots = new TCanvas("c_qa_plots", "QA", 1600, 1400);
    c_qa_plots->Divide(4, 3);
    
    TLegend* leg = new TLegend(0.60, 0.5, 0.88, 0.85);
    leg->SetTextSize(0.04);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

    for ( int i = 0; i < algorithm.size()-1; i ++ )
    {
        process_eff_hist(algorithm[i].h_xq_eff, algorithm[MC].h_xq);
        c_eRecon_eff[i] = draw_2d_efficiency(algorithm[i].h_xq_eff, Form("c_%s_eff", algorithm[i].method_name.c_str()), Form("%s Efficiency", algorithm[i].method_name.c_str()), 1400, 600, false, false);

        for ( int j = 0; j < 4; j ++ )
        {
            c_qa_plots->cd(j+1);
            algorithm[i].h_dq[j]->Draw("HIST SAME");

            c_qa_plots->cd(j+1+4);
            algorithm[i].h_dx[j]->Draw("HIST SAME");

            c_qa_plots->cd(j+1+8);
            algorithm[i].h_dy[j]->Draw("HIST SAME");
        }  
        leg->AddEntry(algorithm[i].h_dq[0], algorithm[i].method_name.c_str(), "L");
    }

    for ( int i = 0; i < 12; i ++ )
    {
        c_qa_plots->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.18);
        leg->Draw();
    }

    //

    TCanvas* c_ycut_qa = new TCanvas("c_ycut_qa", "QA_YCUT", 1800, 950);
    c_ycut_qa->Divide(6, 3);
    
    double ycut_peak_max[3][6] = {0};
    for ( int i = 0; i < algorithm.size()-1; i ++ )
    {
        for ( int j = 0; j < 6; j ++ )
        {
            c_ycut_qa->cd(j+1);
            algorithm[i].h_dq_ycut[j]->Draw("HIST SAME");
            if ( algorithm[i].h_dq_ycut[j]->GetMaximum() > ycut_peak_max[0][j] )
                ycut_peak_max[0][j] = algorithm[i].h_dq_ycut[j]->GetMaximum();

            c_ycut_qa->cd(j+1+6);
            algorithm[i].h_dx_ycut[j]->Draw("HIST SAME");
            if ( algorithm[i].h_dx_ycut[j]->GetMaximum() > ycut_peak_max[1][j] )
                ycut_peak_max[1][j] = algorithm[i].h_dx_ycut[j]->GetMaximum();

            c_ycut_qa->cd(j+1+12);
            algorithm[i].h_dy_ycut[j]->Draw("HIST SAME");
            if ( algorithm[i].h_dy_ycut[j]->GetMaximum() > ycut_peak_max[2][j] )
                ycut_peak_max[2][j] = algorithm[i].h_dy_ycut[j]->GetMaximum();
        }  
    }

    for ( int j = 0; j < 6; j ++ )
    {
        algorithm[0].h_dq_ycut[j]->GetYaxis()->SetRangeUser(0, ycut_peak_max[0][j]*1.3);
        algorithm[0].h_dx_ycut[j]->GetYaxis()->SetRangeUser(0, ycut_peak_max[1][j]*1.3);
        algorithm[0].h_dy_ycut[j]->GetYaxis()->SetRangeUser(0, ycut_peak_max[2][j]*1.3);
    }

    for ( int i = 0; i < 18; i ++ )
    {
        c_ycut_qa->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        // gPad->SetRightMargin(0.05);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.2);
        leg->Draw();
    }

    //

    TCanvas* c_reco_v_true = new TCanvas("c_reco_v_true", "RECO v TRUE", 1800, 1000);
    c_reco_v_true->Divide(algorithm.size()-1, 3);

    for ( int i = 0; i < algorithm.size()-1; i ++ )
    {
        c_reco_v_true->cd(i+1);
        algorithm[i].h_qvq->Draw("COLZ");

        c_reco_v_true->cd(i+1+(algorithm.size()-1));
        algorithm[i].h_xvx->Draw("COLZ");

        c_reco_v_true->cd(i+1+2*(algorithm.size()-1));
        algorithm[i].h_yvy->Draw("COLZ");
    }

    for ( int i = 0; i < (algorithm.size()-1)*3; i ++ )
    {
        c_reco_v_true->cd(i+1);
        gPad->SetLogz();
        if ( i < (algorithm.size()-1)*2 )
        {
            gPad->SetLogx();
            gPad->SetLogy();
        } 
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }

    c_reco_v_true->cd(1);
    gPad->Update();
	TLine* l_q_diagonal = new TLine(1, 1, 1e4, 1e4);
	l_q_diagonal->SetLineColor(kBlack);
	l_q_diagonal->SetLineStyle(7);

    c_reco_v_true->cd(1+(algorithm.size()-1));
    gPad->Update();
	TLine* l_x_diagonal = new TLine(1e-4, 1e-4, 1, 1);
	l_x_diagonal->SetLineColor(kBlack);
	l_x_diagonal->SetLineStyle(7);

    c_reco_v_true->cd(1+2*(algorithm.size()-1));
    gPad->Update();
	TLine* l_y_diagonal = new TLine(0, 0, 1, 1);
	l_y_diagonal->SetLineColor(kBlack);
	l_y_diagonal->SetLineStyle(7);
    
    for ( int i = 0; i < algorithm.size()-1; i ++ )
    {
        c_reco_v_true->cd(i+1);
        l_q_diagonal->Draw("SAME");

        c_reco_v_true->cd(i+1+(algorithm.size()-1));
        l_x_diagonal->Draw("SAME");

        c_reco_v_true->cd(i+1+2*(algorithm.size()-1));
        l_y_diagonal->Draw("SAME");
    }

    // save to tree

    c_qa_plots->Write();
    c_reco_v_true->Write();
    c_ycut_qa->Write();

    for ( int i = 0; i < algorithm.size()-1; i ++ )
        c_eRecon_eff[i]->Write();
}

void plot_energy_and_track()
{
    TCanvas* c1 = new TCanvas("c_cal_and_track_1", "Energy and Theta QA", 1600, 800);
    c1->Divide(4, 2);

    for ( int i = 0; i < 8; i ++ )
    {
        c1->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }
    
    TLegend* leg = new TLegend(0.60, 0.70, 0.88, 0.85);
    leg->SetTextSize(0.04);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_dEt[0], "Track", "L");
    leg->AddEntry(h_dEc[0], "Calorimeter", "L");

    TLegend* leg_a = new TLegend(0.70, 0.5, 0.88, 0.85);
    leg_a->SetTextSize(0.04);
	leg_a->SetBorderSize(0);
	leg_a->SetFillColor(0);
    leg_a->SetFillStyle(0);

    leg_a->AddEntry(h_dp[0], "#phi_{e}", "L");
    leg_a->AddEntry(h_dt[0], "#theta_{e}", "L");
    

    for ( int i = 0; i < 4; i ++ )
    {
        c1->cd(i+1);
        h_dEt[i]->Draw("HIST SAME");
        h_dEt[i]->SetStats(0);
        h_dEt[i]->SetLineColor(kRed);
        format_eid_histogram(h_dEt[i]);

        h_dEc[i]->Draw("HIST SAME");
        h_dEc[i]->SetStats(0);
        h_dEc[i]->SetLineColor(kBlue);
        format_eid_histogram(h_dEc[i]);

        h_dEt[i]->GetYaxis()->SetRangeUser(0, fmax(h_dEt[i]->GetMaximum(), h_dEc[i]->GetMaximum())*1.1);
        h_dEt[i]->SetStats(0);
    
        leg->Draw();

        c1->cd(i+1+4);
        h_dt[i]->Draw("HIST SAME");
        h_dt[i]->SetStats(0);
        format_eid_histogram(h_dt[i]);

        h_dp[i]->Draw("HIST SAME");
        h_dp[i]->SetStats(0);
        h_dp[i]->SetLineColor(kRed);

        leg_a->Draw();
    }

    TCanvas* c2 = new TCanvas("c_cal_and_track_2", "Energy and Theta QA", 1200, 900);
    c2->Divide(3, 3);

    for ( int i = 0; i < 9; i ++ )
    {
        c2->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }

    c2->cd(1);
    h_tde[0]->Draw("COLZ");
    h_tde[0]->SetStats(0);
    format_eid_histogram(h_tde[0]);

    c2->cd(2);
    h_tde[1]->Draw("COLZ");
    h_tde[1]->SetStats(0);
    format_eid_histogram(h_tde[1]);

    c2->cd(3);
    h_tdp->Draw("COLZ");
    h_tdp->SetStats(0);
    format_eid_histogram(h_tdp);

    c2->cd(4);
    h_et->Draw("COLZ");
    h_et->SetStats(0);
    format_eid_histogram(h_et);

    c2->cd(5);
    h_pt->Draw("COLZ");
    h_pt->SetStats(0);
    format_eid_histogram(h_pt);

    c2->cd(7);
    h_ede->Draw("COLZ");
    h_ede->SetStats(0);
    format_eid_histogram(h_ede);

    c2->cd(8);
    h_pde->Draw("COLZ");
    h_pde->SetStats(0);
    format_eid_histogram(h_pde);

    TCanvas* c3 = new TCanvas("c_cal_and_track_3", "Hadronic final state QA", 1600, 800);
    c3->Divide(4, 2);

    for ( int i = 0; i < 8; i ++ )
    {
        c3->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }

    for ( int i = 0; i < 4; i ++ )
    {
        c3->cd(i+1);
        h_dhf[i]->Draw("HIST SAME");
        h_dhf[i]->SetStats(0);
        // h_dhf[i]->SetLineColor(kRed);
        format_eid_histogram(h_dhf[i]);

        c3->cd(i+1+4);
        h_dpt[i]->Draw("HIST SAME");
        h_dpt[i]->SetStats(0);
        format_eid_histogram(h_dpt[i]);
    }

    TCanvas* c4 = new TCanvas("c_cal_and_track_4", "Hadronic final state QA 2", 1600, 800);
    c4->Divide(3, 2);

    for ( int i = 0; i < 6; i ++ )
    {
        c4->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }

    c4->cd(1);
    h_hfs_dpt->Draw("HIST");
    h_hfs_dpt->SetStats(0);
    format_eid_histogram(h_hfs_dpt);

    c4->cd(2);
    h_hfs_dpz->Draw("HIST");
    h_hfs_dpz->SetStats(0);
    format_eid_histogram(h_hfs_dpz);

    c4->cd(3);
    h_hfs_de->Draw("HIST");
    h_hfs_de->SetStats(0);
    format_eid_histogram(h_hfs_de);

    c4->cd(4);
    h_hfs_dpt_t->Draw("COLZ");
    h_hfs_dpt_t->SetStats(0);
    format_eid_histogram(h_hfs_dpt_t);

    c4->cd(5);
    h_hfs_dpt_t->Draw("COLZ");
    h_hfs_dpt_t->SetStats(0);
    format_eid_histogram(h_hfs_dpt_t);

    c4->cd(6);
    h_hfs_de_t->Draw("COLZ");  
    h_hfs_de_t->SetStats(0); 
    format_eid_histogram(h_hfs_de_t); 

    // save to tree
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
            
    return;
}