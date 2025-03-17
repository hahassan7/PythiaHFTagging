#include "iostream"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

void ProcessTHnSparse(bool doProcessTHnSparse = false, bool doProcessIP = false, bool doProcessSV = false)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile *inputFile = TFile::Open("AnalysisResults.root");
    if (!inputFile || inputFile->IsZombie())
    {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    gDirectory->cd();

    std::array<float, 15> workingPoints = {0.0, 0.101, 0.201, 0.301, 0.401, 0.501, 0.601, 0.701, 0.801, 0.901, 0.92, 0.95, 0.97, 0.98, 0.99};
    std::array<double, 6> binning = {5, 10, 20, 40, 70, 200};

    EColor mycolors[] = {kBlack, kRed, kGreen, kMagenta, kBlue, kGray, kOrange, kCyan};

    // Define tagging TGraphs
    TGraph *bjet_taggingGraph[5];

    THnSparseF *hSparse;
    TH2D *incljetScore_pT;
    TH2D *bjetScore_pT;
    TH2D *cjetScore_pT;
    TH2D *lfjetScore_pT;

    if (doProcessTHnSparse)
    {
        // Load the THnSparse histogram
        hSparse = (THnSparseF *)gDirectory->Get("hJetTaggingInfo");
        if (!hSparse)
        {
            std::cerr << "Error loading THnSparse histogram!" << std::endl;
            return;
        }

        hSparse->GetAxis(1)->SetRange(3, 3);
        TH1D *bjets_Total = (TH1D *)hSparse->Projection(0);
        bjets_Total = (TH1D *)bjets_Total->Rebin(binning.size() - 1, "total_bjets_rebinned", binning.data());
        hSparse->GetAxis(1)->SetRange(0, -1);

        for (size_t iwork = 0; iwork < workingPoints.size(); iwork++)
        {
            // Define the range for the working point
            double wpMin = workingPoints[iwork];
            double wpMax = 1.0;

            // Apply the working point range
            int wpMinBin = hSparse->GetAxis(2)->FindBin(wpMin);
            int wpMaxBin = hSparse->GetAxis(2)->FindBin(wpMax);
            hSparse->GetAxis(2)->SetRange(wpMinBin, wpMaxBin);

            // Tagged jets
            TH1D *tagged_incjets = (TH1D *)hSparse->Projection(0)->Clone(Form("tagged_incjets_iwork%d", (int)iwork));
            tagged_incjets = (TH1D *)tagged_incjets->Rebin(binning.size() - 1, Form("tagged_incjets_rebinned_iwork%d", (int)iwork), binning.data());

            hSparse->GetAxis(1)->SetRange(3, 3);
            TH1D *tagged_bjets = (TH1D *)hSparse->Projection(0)->Clone(Form("tagged_bjets_iwork%d", (int)iwork));
            tagged_bjets = (TH1D *)tagged_bjets->Rebin(binning.size() - 1, Form("tagged_bjets_rebinned_iwork%d", (int)iwork), binning.data());

            // Tagging efficiency
            TH1D *bjets_TaggingEfficiency = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingEfficiency_iwork%d", (int)iwork));
            bjets_TaggingEfficiency->Divide(bjets_TaggingEfficiency, bjets_Total, 1.0, 1.0, "B");

            // Tagging Purity
            TH1D *bjets_TaggingPurity = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingPurity_iwork%d", (int)iwork));
            bjets_TaggingPurity->Divide(bjets_TaggingPurity, tagged_incjets, 1.0, 1.0, "B");

            for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
            {
                if (iwork == 0)
                {
                    bjet_taggingGraph[ibin] = new TGraph(workingPoints.size());
                    bjet_taggingGraph[ibin]->SetLineWidth(3);
                }

                std::cout << "Working point: " << workingPoints[iwork] << " bin: " << hSparse->GetAxis(2)->FindBin(workingPoints[iwork]) << ", pT bin: " << binning[ibin] << " - " << binning[ibin + 1] << ", Eff: " << bjets_TaggingEfficiency->GetBinContent(ibin + 1) << ", Pur: " << bjets_TaggingPurity->GetBinContent(ibin + 1) << std::endl;
                bjet_taggingGraph[ibin]->SetPoint(iwork, bjets_TaggingEfficiency->GetBinContent(ibin + 1), bjets_TaggingPurity->GetBinContent(ibin + 1));
            }

            hSparse->GetAxis(0)->SetRange(0, -1);
            hSparse->GetAxis(1)->SetRange(0, -1);
            hSparse->GetAxis(2)->SetRange(0, -1);
        }
    }
    else
    {

        incljetScore_pT = (TH2D *)gDirectory->Get("h2_score_jetpT");
        bjetScore_pT = (TH2D *)gDirectory->Get("h2_score_jetpT_bjet");
        cjetScore_pT = (TH2D *)gDirectory->Get("h2_score_jetpT_cjet");
        lfjetScore_pT = (TH2D *)gDirectory->Get("h2_score_jetpT_lfjet");

        // Total jets
        TH1D *bjets_Total = (TH1D *)bjetScore_pT->ProjectionX("h_jetpT_detector_bjet");
        bjets_Total = (TH1D *)bjets_Total->Rebin(binning.size() - 1, "total_bjets_rebinned", binning.data());

        TCanvas *cTest = new TCanvas("cTest", "Test", 800, 800);
        cTest->cd();
        TH1D *fraction_bjets = (TH1D *)bjetScore_pT->ProjectionX("fraction_bjets", bjetScore_pT->GetYaxis()->GetNbins(), bjetScore_pT->GetYaxis()->GetNbins() + 1);
        fraction_bjets = (TH1D *)fraction_bjets->Rebin(binning.size() - 1, "fraction_bjets_rebinned", binning.data());
        fraction_bjets->Divide(fraction_bjets, bjets_Total, 1.0, 1.0, "B");
        fraction_bjets->Draw();

        for (size_t iwork = 0; iwork < workingPoints.size(); iwork++)
        {
            // Tagged jetsiwork
            TH1D *tagged_incjets = (TH1D *)incljetScore_pT->ProjectionX(Form("tagged_incjets_iwork%d", (int)iwork), incljetScore_pT->GetYaxis()->FindBin(workingPoints[iwork]), bjetScore_pT->GetYaxis()->GetNbins(), "e");
            tagged_incjets = (TH1D *)tagged_incjets->Rebin(binning.size() - 1, Form("tagged_incjets_rebinned_iwork%d", (int)iwork), binning.data());

            TH1D *tagged_bjets = (TH1D *)bjetScore_pT->ProjectionX(Form("tagged_bjets_iwork%d", (int)iwork), bjetScore_pT->GetYaxis()->FindBin(workingPoints[iwork]), bjetScore_pT->GetYaxis()->GetNbins(), "e");
            tagged_bjets = (TH1D *)tagged_bjets->Rebin(binning.size() - 1, Form("tagged_bjets_rebinned_iwork%d", (int)iwork), binning.data());

            // Tagging efficiency
            TH1D *bjets_TaggingEfficiency = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingEfficiency_iwork%d", (int)iwork));
            bjets_TaggingEfficiency->Divide(bjets_TaggingEfficiency, bjets_Total, 1.0, 1.0, "B");

            // Tagging Purity
            TH1D *bjets_TaggingPurity = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingPurity_iwork%d", (int)iwork));
            bjets_TaggingPurity->Divide(bjets_TaggingPurity, tagged_incjets, 1.0, 1.0, "B");

            for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
            {
                if (iwork == 0)
                {
                    bjet_taggingGraph[ibin] = new TGraph(workingPoints.size());
                    bjet_taggingGraph[ibin]->SetLineWidth(3);
                }

                std::cout << "Working point: " << workingPoints[iwork] << " bin: " << bjetScore_pT->GetYaxis()->FindBin(workingPoints[iwork]) << ", pT bin: " << binning[ibin] << " - " << binning[ibin + 1] << ", Eff: " << bjets_TaggingEfficiency->GetBinContent(ibin + 1) << ", Pur: " << bjets_TaggingPurity->GetBinContent(ibin + 1) << std::endl;
                bjet_taggingGraph[ibin]->SetPoint(iwork, bjets_TaggingEfficiency->GetBinContent(ibin + 1), bjets_TaggingPurity->GetBinContent(ibin + 1));
            }
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Purity vs Efficiency", 800, 800);
    c1->cd();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTicks();
    gPad->SetMargin(0.1, 0.01, 0.1, 0.01);
    TH2D *hTemp1 = new TH2D("htemp1", "htemp;Efficiency;Purity", 1, 0, 1, 1, 0, 1);

    hTemp1->Draw();
    TLegend *leg = new TLegend(0.14, 0.14, 0.44, 0.42);
    leg->SetLineWidth(0);
    // leg->SetFillStyle(0);
    for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
    {
        bjet_taggingGraph[ibin]->SetLineColor(mycolors[ibin]);
        bjet_taggingGraph[ibin]->SetMarkerColor(mycolors[ibin]);
        bjet_taggingGraph[ibin]->SetMarkerStyle(21);
        bjet_taggingGraph[ibin]->Draw("same");
        leg->AddEntry(bjet_taggingGraph[ibin], Form("#it{p}_{T}^{jet} : %.f - %.f GeV/#it{c}", binning[ibin], binning[ibin + 1]));
    }
    leg->Draw();

    TLatex *textpp = new TLatex();
    textpp->SetTextSizePixels(25);
    textpp->DrawLatexNDC(0.6, 0.94, "PYTHIA8 pp #sqrt{#it{s}} = 13 TeV");
    textpp->DrawLatexNDC(0.57, 0.89, "Anti-#it{k}_{T} jet, #it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");

    // Plot score distributions for different jet flavors
    TCanvas *c2 = new TCanvas("c2", "Score Distributions", 1200, 800);
    c2->Divide(3, 2);

    TH2D *hTemp2 = new TH2D("htemp2", "htemp;Efficiency;Purity", 500, 0, 1, 500, 1e-1, 1e3);

    for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
    {
        c2->cd(ibin + 1);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetTicks();

        TH1D *hScore_incl;
        TH1D *hScore_bjet;
        TH1D *hScore_cjet;
        TH1D *hScore_lfjet;

        if (doProcessTHnSparse)
        {
            // Apply pT bin range
            hSparse->GetAxis(0)->SetRangeUser(binning[ibin], binning[ibin + 1]);

            // Project score for different jet flavors
            hScore_incl = (TH1D *)hSparse->Projection(2)->Clone(Form("hScore_incl_ibin%d", (int)ibin));
            hSparse->GetAxis(1)->SetRange(3, 3);
            hScore_bjet = (TH1D *)hSparse->Projection(2)->Clone(Form("hScore_bjet_ibin%d", (int)ibin));
            hSparse->GetAxis(1)->SetRange(2, 2);
            hScore_cjet = (TH1D *)hSparse->Projection(2)->Clone(Form("hScore_cjet_ibin%d", (int)ibin));
            hSparse->GetAxis(1)->SetRange(1, 1);
            hScore_lfjet = (TH1D *)hSparse->Projection(2)->Clone(Form("hScore_ljet_ibin%d", (int)ibin));

            hScore_incl->SetLineColor(kBlack);
            hScore_bjet->SetLineColor(kRed);
            hScore_cjet->SetLineColor(kGreen);
            hScore_lfjet->SetLineColor(kBlue);

            hScore_incl->SetLineWidth(2);
            hScore_bjet->SetLineWidth(2);
            hScore_cjet->SetLineWidth(2);
            hScore_lfjet->SetLineWidth(2);

            hTemp2->Draw();

            hScore_incl->Draw("HIST SAME");
            hScore_bjet->Draw("HIST SAME");
            hScore_cjet->Draw("HIST SAME");
            hScore_lfjet->Draw("HIST SAME");

            hSparse->GetAxis(0)->SetRange(0, -1);
            hSparse->GetAxis(1)->SetRange(0, -1);
            hSparse->GetAxis(2)->SetRange(0, -1);
        }
        else
        {
            hScore_incl = incljetScore_pT->ProjectionY(Form("hScore_incl_ibin%d", (int)ibin), incljetScore_pT->GetXaxis()->FindBin(binning[ibin]), incljetScore_pT->GetXaxis()->FindBin(binning[ibin + 1]));
            hScore_bjet = bjetScore_pT->ProjectionY(Form("hScore_bjet_ibin%d", (int)ibin), bjetScore_pT->GetXaxis()->FindBin(binning[ibin]), bjetScore_pT->GetXaxis()->FindBin(binning[ibin + 1]));
            hScore_cjet = cjetScore_pT->ProjectionY(Form("hScore_cjet_ibin%d", (int)ibin), cjetScore_pT->GetXaxis()->FindBin(binning[ibin]), cjetScore_pT->GetXaxis()->FindBin(binning[ibin + 1]));
            hScore_lfjet = lfjetScore_pT->ProjectionY(Form("hScore_ljet_ibin%d", (int)ibin), lfjetScore_pT->GetXaxis()->FindBin(binning[ibin]), lfjetScore_pT->GetXaxis()->FindBin(binning[ibin + 1]));

            hScore_incl->Rebin(2);
            hScore_bjet->Rebin(2);
            hScore_cjet->Rebin(2);
            hScore_lfjet->Rebin(2);

            hScore_incl->SetLineColor(kBlack);
            hScore_bjet->SetLineColor(kRed);
            hScore_cjet->SetLineColor(kGreen);
            hScore_lfjet->SetLineColor(kBlue);

            hScore_incl->SetLineWidth(2);
            hScore_bjet->SetLineWidth(2);
            hScore_cjet->SetLineWidth(2);
            hScore_lfjet->SetLineWidth(2);

            hScore_incl->GetYaxis()->SetRangeUser(1e-9, hScore_incl->GetMaximum() * 2);
            hScore_incl->Draw("HIST");
            hScore_bjet->Draw("HIST SAME");
            hScore_cjet->Draw("HIST SAME");
            hScore_lfjet->Draw("HIST SAME");
        }

        TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg2->AddEntry(hScore_incl, "Inclusive Jets", "l");
        leg2->AddEntry(hScore_bjet, "b-Jets", "l");
        leg2->AddEntry(hScore_cjet, "c-Jets", "l");
        leg2->AddEntry(hScore_lfjet, "Light Flavor Jets", "l");
        leg2->Draw();
    }

    if (doProcessIP)
    {
        std::array<float, 15> workingPointsIP = {0.0, 0.001, 0.002, 0.004, 0.008, 0.01, 0.02, 0.04, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35};

        // Define tagging TGraphs
        TGraph *bjet_taggingGraph_IPs[5];

        TH2D *h2IPxyJetpTN2 = (TH2D *)gDirectory->Get("h2IPxyJetpTN3");
        TH2D *h2IPxyJetpTN2_bjet = (TH2D *)gDirectory->Get("h2IPxyJetpTN3_bjet");

        TH1D *bjets_Total = (TH1D *)gDirectory->Get("hJetPt_b");
        bjets_Total = (TH1D *)bjets_Total->Rebin(binning.size() - 1, "total_bjets_rebinned", binning.data());

        for (size_t iwork = 0; iwork < workingPointsIP.size(); iwork++)
        {
            // Tagged jetsiwork
            TH1D *tagged_incjets = (TH1D *)h2IPxyJetpTN2->ProjectionX(Form("tagged_incjets_iwork%d", (int)iwork), h2IPxyJetpTN2->GetYaxis()->FindBin(workingPointsIP[iwork]), h2IPxyJetpTN2_bjet->GetYaxis()->GetNbins(), "e");
            tagged_incjets = (TH1D *)tagged_incjets->Rebin(binning.size() - 1, Form("tagged_incjets_rebinned_iwork%d", (int)iwork), binning.data());

            TH1D *tagged_bjets = (TH1D *)h2IPxyJetpTN2_bjet->ProjectionX(Form("tagged_bjets_iwork%d", (int)iwork), h2IPxyJetpTN2_bjet->GetYaxis()->FindBin(workingPointsIP[iwork]), h2IPxyJetpTN2_bjet->GetYaxis()->GetNbins(), "e");
            tagged_bjets = (TH1D *)tagged_bjets->Rebin(binning.size() - 1, Form("tagged_bjets_rebinned_iwork%d", (int)iwork), binning.data());

            // Tagging efficiency
            TH1D *bjets_TaggingEfficiency = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingEfficiency_iwork%d", (int)iwork));
            bjets_TaggingEfficiency->Divide(bjets_TaggingEfficiency, bjets_Total, 1.0, 1.0, "B");

            // Tagging Purity
            TH1D *bjets_TaggingPurity = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingPurity_iwork%d", (int)iwork));
            bjets_TaggingPurity->Divide(bjets_TaggingPurity, tagged_incjets, 1.0, 1.0, "B");

            for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
            {
                if (iwork == 0)
                {
                    bjet_taggingGraph_IPs[ibin] = new TGraph(workingPointsIP.size());
                    bjet_taggingGraph_IPs[ibin]->SetLineWidth(3);
                }

                std::cout << "Working point: " << workingPointsIP[iwork] << " bin: " << h2IPxyJetpTN2_bjet->GetYaxis()->FindBin(workingPointsIP[iwork]) << ", pT bin: " << binning[ibin] << " - " << binning[ibin + 1] << ", Eff: " << bjets_TaggingEfficiency->GetBinContent(ibin + 1) << ", Pur: " << bjets_TaggingPurity->GetBinContent(ibin + 1) << std::endl;
                bjet_taggingGraph_IPs[ibin]->SetPoint(iwork, bjets_TaggingEfficiency->GetBinContent(ibin + 1), bjets_TaggingPurity->GetBinContent(ibin + 1));
            }
        }

        TCanvas *cIPs = new TCanvas("cIPs", "Purity vs Efficiency for IP", 800, 800);
        cIPs->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetTicks();
        gPad->SetMargin(0.1, 0.01, 0.1, 0.01);
        TH2D *hTempIP = new TH2D("htempIP", "htemp;Efficiency;Purity", 1, 0, 1, 1, 0, 1);

        hTempIP->Draw();
        TLegend *leg1 = new TLegend(0.14, 0.14, 0.44, 0.42);
        leg1->SetLineWidth(0);
        // leg->SetFillStyle(0);
        for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
        {
            bjet_taggingGraph_IPs[ibin]->SetLineColor(mycolors[ibin]);
            bjet_taggingGraph_IPs[ibin]->SetMarkerColor(mycolors[ibin]);
            bjet_taggingGraph_IPs[ibin]->SetMarkerStyle(21);
            bjet_taggingGraph_IPs[ibin]->Draw("same");
            leg1->AddEntry(bjet_taggingGraph[ibin], Form("#it{p}_{T}^{jet} : %.f - %.f GeV/#it{c}", binning[ibin], binning[ibin + 1]));
        }
        leg1->Draw();
    }

    if (doProcessSV)
    {
        std::array<float, 15> workingPointsSV = {0.0, 0.05, 0.1, 0.15, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.8, 2, 2.5, 3, 3.5};
        double maxDispersion = 0.1;

        // Define tagging TGraphs
        TGraph *bjet_taggingGraph_SV[5];

        TH3D *h3DecayLengthDispersionJetpT = (TH3D *)gDirectory->Get("h3DecayLengthDispersionJetpT");
        TH3D *h3DecayLengthDispersionJetpT_bjet = (TH3D *)gDirectory->Get("h3DecayLengthDispersionJetpT_bjet");

        TH1D *bjets_Total = (TH1D *)gDirectory->Get("hJetPt_b");
        bjets_Total = (TH1D *)bjets_Total->Rebin(binning.size() - 1, "total_bjets_rebinned", binning.data());

        for (size_t iwork = 0; iwork < workingPointsSV.size(); iwork++)
        {
            // Tagged jetsiwork
            TH1D *tagged_incjets = (TH1D *)h3DecayLengthDispersionJetpT->ProjectionX(Form("tagged_incjets_iwork%d", (int)iwork), h3DecayLengthDispersionJetpT->GetYaxis()->FindBin(workingPointsSV[iwork]), h3DecayLengthDispersionJetpT->GetYaxis()->GetNbins(), 1, h3DecayLengthDispersionJetpT->GetZaxis()->FindBin(maxDispersion), "e");
            tagged_incjets = (TH1D *)tagged_incjets->Rebin(binning.size() - 1, Form("tagged_incjets_rebinned_iwork%d", (int)iwork), binning.data());

            TH1D *tagged_bjets = (TH1D *)h3DecayLengthDispersionJetpT_bjet->ProjectionX(Form("tagged_bjets_iwork%d", (int)iwork), h3DecayLengthDispersionJetpT_bjet->GetYaxis()->FindBin(workingPointsSV[iwork]), h3DecayLengthDispersionJetpT_bjet->GetYaxis()->GetNbins(), 1, h3DecayLengthDispersionJetpT_bjet->GetZaxis()->FindBin(maxDispersion), "e");
            tagged_bjets = (TH1D *)tagged_bjets->Rebin(binning.size() - 1, Form("tagged_bjets_rebinned_iwork%d", (int)iwork), binning.data());

            // Tagging efficiency
            TH1D *bjets_TaggingEfficiency = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingEfficiency_iwork%d", (int)iwork));
            bjets_TaggingEfficiency->Divide(bjets_TaggingEfficiency, bjets_Total, 1.0, 1.0, "B");

            // Tagging Purity
            TH1D *bjets_TaggingPurity = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingPurity_iwork%d", (int)iwork));
            bjets_TaggingPurity->Divide(bjets_TaggingPurity, tagged_incjets, 1.0, 1.0, "B");

            for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
            {
                if (iwork == 0)
                {
                    bjet_taggingGraph_SV[ibin] = new TGraph(workingPointsSV.size());
                    bjet_taggingGraph_SV[ibin]->SetLineWidth(3);
                }

                std::cout << "Working point: " << workingPointsSV[iwork] << " bin: " << h3DecayLengthDispersionJetpT_bjet->GetYaxis()->FindBin(workingPointsSV[iwork]) << ", pT bin: " << binning[ibin] << " - " << binning[ibin + 1] << ", Eff: " << bjets_TaggingEfficiency->GetBinContent(ibin + 1) << ", Pur: " << bjets_TaggingPurity->GetBinContent(ibin + 1) << std::endl;
                bjet_taggingGraph_SV[ibin]->SetPoint(iwork, bjets_TaggingEfficiency->GetBinContent(ibin + 1), bjets_TaggingPurity->GetBinContent(ibin + 1));
            }
        }

        TCanvas *cSV = new TCanvas("cSV", "Purity vs Efficiency for SV", 800, 800);
        cSV->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetTicks();
        gPad->SetMargin(0.1, 0.01, 0.1, 0.01);
        TH2D *hTempSV = new TH2D("hTempSV", "htemp;Efficiency;Purity", 1, 0, 1, 1, 0, 1);

        hTempSV->Draw();
        TLegend *leg1 = new TLegend(0.14, 0.14, 0.44, 0.42);
        leg1->SetLineWidth(0);
        // leg->SetFillStyle(0);
        for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
        {
            bjet_taggingGraph_SV[ibin]->SetLineColor(mycolors[ibin]);
            bjet_taggingGraph_SV[ibin]->SetMarkerColor(mycolors[ibin]);
            bjet_taggingGraph_SV[ibin]->SetMarkerStyle(21);
            bjet_taggingGraph_SV[ibin]->Draw("same");
            leg1->AddEntry(bjet_taggingGraph[ibin], Form("#it{p}_{T}^{jet} : %.f - %.f GeV/#it{c}", binning[ibin], binning[ibin + 1]));
        }
        leg1->Draw();
    }
}