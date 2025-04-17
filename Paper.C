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

void Paper(bool saveFigs = false)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile *inputFile = TFile::Open("OutputFiles/Results_pTHat_ML_New.root");
    if (!inputFile || inputFile->IsZombie())
    {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    gDirectory->cd();

    // std::array<float, 15> workingPoints = {0.0, 0.101, 0.201, 0.301, 0.401, 0.501, 0.601, 0.701, 0.801, 0.901, 0.92, 0.95, 0.97, 0.98, 0.99};
    std::array<float, 16> workingPoints = {0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 6, 7, 8};
    std::array<double, 6> binning = {5, 10, 20, 40, 70, 200};

    // EColor mycolors[] = {kBlack, kRed, kGreen, kMagenta, kBlue, kGray, kOrange, kCyan};
    Color_t mycolors[] = {kAzure + 4, kOrange + 7, kTeal + 3, kViolet - 5, kSpring - 8};

    // Define tagging TGraphs
    TGraph *bjet_taggingGraph[5];

    // TH2D *incljetScore_pT = (TH2D *)gDirectory->Get("h2_score_jetpT");
    // TH2D *bjetScore_pT = (TH2D *)gDirectory->Get("h2_score_jetpT_bjet");
    // TH2D *cjetScore_pT = (TH2D *)gDirectory->Get("h2_score_jetpT_cjet");
    // TH2D *lfjetScore_pT = (TH2D *)gDirectory->Get("h2_score_jetpT_lfjet");

    TH2D *incljetScore_pT = (TH2D *)gDirectory->Get("h2_logscore_jetpT");
    TH2D *bjetScore_pT = (TH2D *)gDirectory->Get("h2_logscore_jetpT_bjet");
    TH2D *cjetScore_pT = (TH2D *)gDirectory->Get("h2_logscore_jetpT_cjet");
    TH2D *lfjetScore_pT = (TH2D *)gDirectory->Get("h2_logscore_jetpT_lfjet");

    // Total jets
    TH1D *bjets_Total = (TH1D *)bjetScore_pT->ProjectionX("h_jetpT_detector_bjet");
    bjets_Total = (TH1D *)bjets_Total->Rebin(binning.size() - 1, "total_bjets_rebinned", binning.data());

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

    TCanvas *c1 = new TCanvas("c1", "Purity vs Efficiency", 800, 800);
    c1->cd();
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetTicks();
    gPad->SetMargin(0.1, 0.01, 0.1, 0.01);
    TH2D *hTemp1 = new TH2D("htemp1", "htemp;b-jet tagging efficiency;b-jet purity", 1, 0, 1, 1, 0, 1.5);
    hTemp1->GetYaxis()->SetTitleSize(0.045);
    hTemp1->GetYaxis()->SetTitleOffset(1.1);
    hTemp1->GetXaxis()->SetTitleSize(0.045);
    hTemp1->GetXaxis()->SetTitleOffset(0.9);

    hTemp1->Draw();
    TLegend *leg = new TLegend(0.54, 0.67, 0.86, 0.96);
    leg->SetTextSizePixels(26);
    leg->SetLineWidth(0);
    // leg->SetFillStyle(0);
    for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
    {
        bjet_taggingGraph[ibin]->SetLineColor(mycolors[ibin]);
        bjet_taggingGraph[ibin]->SetMarkerColor(mycolors[ibin]);
        bjet_taggingGraph[ibin]->SetMarkerStyle(0);
        bjet_taggingGraph[ibin]->Draw("same");
        // leg->AddEntry(bjet_taggingGraph[ibin], Form("#it{p}_{T}^{jet} : %.f - %.f GeV/#it{c}", binning[ibin], binning[ibin + 1]));
        leg->AddEntry(bjet_taggingGraph[ibin], Form("%.f#kern[-0.5]{ }#leq#kern[-0.5]{ }#it{p}_{T,ch jet} < %.f GeV/#it{c}", binning[ibin], binning[ibin + 1]));
    }
    leg->Draw();

    TLatex *textpp = new TLatex();
    textpp->SetTextFont(42);
    textpp->SetTextSizePixels(27);
    textpp->DrawLatexNDC(0.14, 0.93, "PYTHIA 8, pp,#kern[-0.5]{ }#sqrt{#it{s}} = 13 TeV");
    textpp->DrawLatexNDC(0.14, 0.885, "Charged-particle jets");
    textpp->DrawLatexNDC(0.14, 0.84, "Anti-#it{k}_{T},#kern[-0.5]{ }#it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");

    if (saveFigs)
    {
        c1->SaveAs("PaperFigures/PurityVsEfficiency.pdf");
    }

    // Plot score distributions for different jet flavors
    TCanvas *c2 = new TCanvas("c2", "Score Distributions", 800, 800);
    c2->cd();
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetTicks();
    gPad->SetLogy();
    gPad->SetMargin(0.11, 0.02, 0.1, 0.01);

    int mainBin = 3;

    TH1D *hScore_incl = incljetScore_pT->ProjectionY(Form("hScore_incl_ibin%d", (int)mainBin), incljetScore_pT->GetXaxis()->FindBin(binning[mainBin]), incljetScore_pT->GetXaxis()->FindBin(binning[mainBin + 1]));
    TH1D *hScore_bjet = bjetScore_pT->ProjectionY(Form("hScore_bjet_ibin%d", (int)mainBin), bjetScore_pT->GetXaxis()->FindBin(binning[mainBin]), bjetScore_pT->GetXaxis()->FindBin(binning[mainBin + 1]));
    TH1D *hScore_cjet = cjetScore_pT->ProjectionY(Form("hScore_cjet_ibin%d", (int)mainBin), cjetScore_pT->GetXaxis()->FindBin(binning[mainBin]), cjetScore_pT->GetXaxis()->FindBin(binning[mainBin + 1]));
    TH1D *hScore_lfjet = lfjetScore_pT->ProjectionY(Form("hScore_ljet_ibin%d", (int)mainBin), lfjetScore_pT->GetXaxis()->FindBin(binning[mainBin]), lfjetScore_pT->GetXaxis()->FindBin(binning[mainBin + 1]));

    // hScore_incl->Rebin(2);
    // hScore_bjet->Rebin(2);
    // hScore_cjet->Rebin(2);
    // hScore_lfjet->Rebin(2);

    hScore_incl->SetLineColor(kBlack);
    hScore_bjet->SetLineColor(kRed);
    hScore_cjet->SetLineColor(kGreen + 2);
    hScore_lfjet->SetLineColor(kBlue);

    hScore_incl->GetXaxis()->SetRangeUser(0, 20);
    hScore_incl->GetXaxis()->SetTitle("#minus^{}ln(1#kern[-0.5]{ }#minus Score)");
    hScore_incl->GetYaxis()->SetTitle("Probability density");
    hScore_incl->GetYaxis()->SetTitleSize(0.045);
    hScore_incl->GetYaxis()->SetTitleOffset(1.2);
    hScore_incl->GetXaxis()->SetTitleSize(0.045);
    hScore_incl->GetXaxis()->SetTitleOffset(0.9);
    hScore_incl->SetLineWidth(2);
    hScore_bjet->SetLineWidth(2);
    hScore_cjet->SetLineWidth(2);
    hScore_lfjet->SetLineWidth(2);

    hScore_bjet->Scale(1.0 / hScore_incl->Integral());
    hScore_cjet->Scale(1.0 / hScore_incl->Integral());
    hScore_lfjet->Scale(1.0 / hScore_incl->Integral());
    hScore_incl->Scale(1.0 / hScore_incl->Integral());

    hScore_incl->Draw("HIST");
    hScore_bjet->Draw("HIST SAME");
    hScore_cjet->Draw("HIST SAME");
    hScore_lfjet->Draw("HIST SAME");

    TLegend *leg2 = new TLegend(0.57, 0.46, 0.91, 0.66);
    leg2->SetLineWidth(0);
    leg2->AddEntry(hScore_incl, "inclusive jets", "l");
    leg2->AddEntry(hScore_bjet, "beauty jets", "l");
    leg2->AddEntry(hScore_cjet, "charm jets", "l");
    leg2->AddEntry(hScore_lfjet, "light-flavor jets", "l");
    leg2->Draw();

    textpp->DrawLatexNDC(0.55, 0.93, "PYTHIA 8, pp,#kern[-0.5]{ }#sqrt{#it{s}} = 13 TeV");
    textpp->DrawLatexNDC(0.55, 0.885, "Charged-particle jets");
    textpp->DrawLatexNDC(0.55, 0.84, "Anti-#it{k}_{T},#kern[-0.5]{ }#it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");
    textpp->DrawLatexNDC(0.55, 0.795, Form("%.f#kern[-0.5]{ }#leq#kern[-0.5]{ }#it{p}_{T,ch jet} < %.f GeV/#it{c}", binning[mainBin], binning[mainBin + 1]));

    gPad->RedrawAxis();

    if (saveFigs)
    {
        c2->SaveAs("PaperFigures/ScoreDistributions.pdf");
    }

    std::array<float, 15> workingPointsIP = {0.0, 0.001, 0.002, 0.004, 0.008, 0.01, 0.02, 0.04, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35};

    // Define tagging TGraphs
    TGraph *bjet_taggingGraph_IPs[5];

    TH2D *h2IPxyJetpTN2 = (TH2D *)gDirectory->Get("h2IPzJetpTN3");
    TH2D *h2IPxyJetpTN2_bjet = (TH2D *)gDirectory->Get("h2IPzJetpTN3_bjet");
    TH2D *h2IPxyJetpTN2_cjet = (TH2D *)gDirectory->Get("h2IPzJetpTN3_cjet");
    TH2D *h2IPxyJetpTN2_ljet = (TH2D *)h2IPxyJetpTN2->Clone("h2IPxyJetpTN3_ljet");
    h2IPxyJetpTN2_ljet->Add(h2IPxyJetpTN2_bjet, -1);
    h2IPxyJetpTN2_ljet->Add(h2IPxyJetpTN2_cjet, -1);

    // TH1D *bjets_Total_IPs = (TH1D *)gDirectory->Get("hJetPt_b");
    TH1D *bjets_Total_IPs = (TH1D *)h2IPxyJetpTN2_bjet->ProjectionX("bjets_Total_IPs");
    bjets_Total_IPs = (TH1D *)bjets_Total_IPs->Rebin(binning.size() - 1, "total_bjets_rebinned", binning.data());

    for (size_t iwork = 0; iwork < workingPointsIP.size(); iwork++)
    {
        // Tagged jetsiwork
        TH1D *tagged_incjets = (TH1D *)h2IPxyJetpTN2->ProjectionX(Form("tagged_incjets_iwork%d", (int)iwork), h2IPxyJetpTN2->GetYaxis()->FindBin(workingPointsIP[iwork]), h2IPxyJetpTN2_bjet->GetYaxis()->GetNbins(), "e");
        tagged_incjets = (TH1D *)tagged_incjets->Rebin(binning.size() - 1, Form("tagged_incjets_rebinned_iwork%d", (int)iwork), binning.data());

        TH1D *tagged_bjets = (TH1D *)h2IPxyJetpTN2_bjet->ProjectionX(Form("tagged_bjets_iwork%d", (int)iwork), h2IPxyJetpTN2_bjet->GetYaxis()->FindBin(workingPointsIP[iwork]), h2IPxyJetpTN2_bjet->GetYaxis()->GetNbins(), "e");
        tagged_bjets = (TH1D *)tagged_bjets->Rebin(binning.size() - 1, Form("tagged_bjets_rebinned_iwork%d", (int)iwork), binning.data());

        // Tagging efficiency
        TH1D *bjets_TaggingEfficiency = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingEfficiency_iwork%d", (int)iwork));
        bjets_TaggingEfficiency->Divide(bjets_TaggingEfficiency, bjets_Total_IPs, 1.0, 1.0, "B");

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

            bjet_taggingGraph_IPs[ibin]->SetPoint(iwork, bjets_TaggingEfficiency->GetBinContent(ibin + 1), bjets_TaggingPurity->GetBinContent(ibin + 1));
        }
    }

    for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
    {
        bjet_taggingGraph_IPs[ibin]->SetLineColor(mycolors[ibin]);
        bjet_taggingGraph_IPs[ibin]->SetMarkerColor(mycolors[ibin]);
        bjet_taggingGraph_IPs[ibin]->SetMarkerStyle(21);
    }

    // Plot IPs distributions for different jet flavors
    TCanvas *cDCA = new TCanvas("cDCA", "DCA Distributions", 800, 800);
    cDCA->cd();
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetTicks();
    gPad->SetLogy();
    gPad->SetMargin(0.11, 0.01, 0.1, 0.01);

    int mainBinIPs = 4;

    TH1D *hDCA_bjet = h2IPxyJetpTN2_bjet->ProjectionY(Form("hDCA_bjet_ibin%d", (int)mainBinIPs), h2IPxyJetpTN2_bjet->GetXaxis()->FindBin(binning[mainBinIPs]), h2IPxyJetpTN2_bjet->GetXaxis()->FindBin(binning[mainBinIPs + 1]));
    TH1D *hDCA_cjet = h2IPxyJetpTN2_cjet->ProjectionY(Form("hDCA_cjet_ibin%d", (int)mainBinIPs), h2IPxyJetpTN2_cjet->GetXaxis()->FindBin(binning[mainBinIPs]), h2IPxyJetpTN2_cjet->GetXaxis()->FindBin(binning[mainBinIPs + 1]));
    TH1D *hDCA_ljet = h2IPxyJetpTN2_ljet->ProjectionY(Form("hDCA_ljet_ibin%d", (int)mainBinIPs), h2IPxyJetpTN2_ljet->GetXaxis()->FindBin(binning[mainBinIPs]), h2IPxyJetpTN2_ljet->GetXaxis()->FindBin(binning[mainBinIPs + 1]));

    hDCA_bjet->Scale(1.0 / hDCA_bjet->Integral());
    hDCA_cjet->Scale(1.0 / hDCA_cjet->Integral());
    hDCA_ljet->Scale(1.0 / hDCA_ljet->Integral());

    // hDCA_bjet->Rebin(2);
    // hDCA_cjet->Rebin(2);
    // hDCA_ljet->Rebin(2);

    hDCA_bjet->SetLineColor(kRed);
    hDCA_bjet->SetMarkerColor(kRed);
    hDCA_cjet->SetLineColor(kGreen + 2);
    hDCA_cjet->SetMarkerColor(kGreen + 2);
    hDCA_ljet->SetLineColor(kBlue);
    hDCA_ljet->SetMarkerColor(kBlue);

    hDCA_bjet->SetLineWidth(2);
    hDCA_cjet->SetLineWidth(2);
    hDCA_ljet->SetLineWidth(2);

    hDCA_ljet->GetXaxis()->SetRangeUser(-1, 1);
    hDCA_ljet->GetXaxis()->SetTitle("signed DCA [cm]");
    hDCA_ljet->GetYaxis()->SetTitle("Probability Distribution");
    hDCA_ljet->GetYaxis()->SetTitleSize(0.045);
    hDCA_ljet->GetYaxis()->SetTitleOffset(1.2);
    hDCA_ljet->GetXaxis()->SetTitleSize(0.045);
    hDCA_ljet->GetXaxis()->SetTitleOffset(0.9);

    hDCA_ljet->Draw("");
    hDCA_bjet->Draw("SAME");
    hDCA_cjet->Draw("SAME");

    TLegend *legIP = new TLegend(0.15, 0.46, 0.49, 0.66);
    legIP->SetLineWidth(0);
    legIP->AddEntry(hDCA_bjet, "beauty jets", "l");
    legIP->AddEntry(hDCA_cjet, "charm jets", "l");
    legIP->AddEntry(hDCA_ljet, "light-flavor jets", "l");
    legIP->Draw();

    textpp->DrawLatexNDC(0.14, 0.93, "PYTHIA 8, pp,#kern[-0.5]{ }#sqrt{#it{s}} = 13 TeV");
    textpp->DrawLatexNDC(0.14, 0.885, "Charged-particle jets");
    textpp->DrawLatexNDC(0.14, 0.84, "Anti-#it{k}_{T},#kern[-0.5]{ }#it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");
    textpp->DrawLatexNDC(0.14, 0.795, Form("%.f#kern[-0.5]{ }#leq#kern[-0.5]{ }#it{p}_{T,ch jet} < %.f GeV/#it{c}", binning[mainBin], binning[mainBin + 1]));

    gPad->RedrawAxis();

    if (saveFigs)
    {
        cDCA->SaveAs("PaperFigures/DCA_Distributions.pdf");
    }

    std::array<float, 16> workingPointsSV = {0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.8};
    // std::array<float, 15> workingPointsSV = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.8};
    double maxDispersion = 0.03;

    // Define tagging TGraphs
    TGraph *bjet_taggingGraph_SV[5];

    TH3D *h3DecayLengthDispersionJetpT = (TH3D *)gDirectory->Get("h3DecayLengthDispersionJetpT");
    TH3D *h3DecayLengthDispersionJetpT_bjet = (TH3D *)gDirectory->Get("h3DecayLengthDispersionJetpT_bjet");
    TH3D *h3DecayLengthDispersionJetpT_cjet = (TH3D *)gDirectory->Get("h3DecayLengthDispersionJetpT_cjet");
    TH3D *h3DecayLengthDispersionJetpT_lfjet = (TH3D *)h3DecayLengthDispersionJetpT->Clone("h3DecayLengthDispersionJetpT_lfjet");
    h3DecayLengthDispersionJetpT_lfjet->Add(h3DecayLengthDispersionJetpT_bjet, -1);
    h3DecayLengthDispersionJetpT_lfjet->Add(h3DecayLengthDispersionJetpT_cjet, -1);

    // TH1D *bjets_Total_SV = (TH1D *)gDirectory->Get("hJetPt_b");
    TH1D *bjets_Total_SV = (TH1D *)h3DecayLengthDispersionJetpT_bjet->ProjectionX("bjets_Total_SV");
    bjets_Total_SV = (TH1D *)bjets_Total_SV->Rebin(binning.size() - 1, "total_bjets_rebinned", binning.data());

    for (size_t iwork = 0; iwork < workingPointsSV.size(); iwork++)
    {
        // Tagged jetsiwork
        TH1D *tagged_incjets = (TH1D *)h3DecayLengthDispersionJetpT->ProjectionX(Form("tagged_incjets_iwork%d", (int)iwork), h3DecayLengthDispersionJetpT->GetYaxis()->FindBin(workingPointsSV[iwork]), h3DecayLengthDispersionJetpT->GetYaxis()->GetNbins(), 1, h3DecayLengthDispersionJetpT->GetZaxis()->FindBin(maxDispersion), "e");
        tagged_incjets = (TH1D *)tagged_incjets->Rebin(binning.size() - 1, Form("tagged_incjets_rebinned_iwork%d", (int)iwork), binning.data());

        TH1D *tagged_bjets = (TH1D *)h3DecayLengthDispersionJetpT_bjet->ProjectionX(Form("tagged_bjets_iwork%d", (int)iwork), h3DecayLengthDispersionJetpT_bjet->GetYaxis()->FindBin(workingPointsSV[iwork]), h3DecayLengthDispersionJetpT_bjet->GetYaxis()->GetNbins(), 1, h3DecayLengthDispersionJetpT_bjet->GetZaxis()->FindBin(maxDispersion), "e");
        tagged_bjets = (TH1D *)tagged_bjets->Rebin(binning.size() - 1, Form("tagged_bjets_rebinned_iwork%d", (int)iwork), binning.data());

        // Tagging efficiency
        TH1D *bjets_TaggingEfficiency = (TH1D *)tagged_bjets->Clone(Form("bjets_TaggingEfficiency_iwork%d", (int)iwork));
        bjets_TaggingEfficiency->Divide(bjets_TaggingEfficiency, bjets_Total_SV, 1.0, 1.0, "B");

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

            bjet_taggingGraph_SV[ibin]->SetPoint(iwork, bjets_TaggingEfficiency->GetBinContent(ibin + 1), bjets_TaggingPurity->GetBinContent(ibin + 1));
        }
    }

    for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
    {
        bjet_taggingGraph_SV[ibin]->SetLineColor(mycolors[ibin]);
        bjet_taggingGraph_SV[ibin]->SetMarkerColor(mycolors[ibin]);
        bjet_taggingGraph_SV[ibin]->SetMarkerStyle(21);
    }

    TCanvas *cDL = new TCanvas("cDL", "Decay Length Distributions", 800, 800);

    cDL->cd();
    gPad->SetMargin(0.11, 0.02, 0.1, 0.01);
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetLogy();
    gPad->SetTicks();

    int mainBinDL = 3;

    TH1D *hDecayLength_bjet = h3DecayLengthDispersionJetpT_bjet->ProjectionY(Form("hDecayLength_bjet_ibin%d", (int)mainBinIPs), h3DecayLengthDispersionJetpT_bjet->GetXaxis()->FindBin(binning[mainBinIPs]), h3DecayLengthDispersionJetpT_bjet->GetXaxis()->FindBin(binning[mainBinIPs + 1]), 1, h3DecayLengthDispersionJetpT_bjet->GetZaxis()->FindBin(maxDispersion), "e");
    TH1D *hDecayLength_cjet = h3DecayLengthDispersionJetpT_cjet->ProjectionY(Form("hDecayLength_cjet_ibin%d", (int)mainBinIPs), h3DecayLengthDispersionJetpT_cjet->GetXaxis()->FindBin(binning[mainBinIPs]), h3DecayLengthDispersionJetpT_cjet->GetXaxis()->FindBin(binning[mainBinIPs + 1]), 1, h3DecayLengthDispersionJetpT_cjet->GetZaxis()->FindBin(maxDispersion), "e");
    TH1D *hDecayLength_ljet = h3DecayLengthDispersionJetpT_lfjet->ProjectionY(Form("hDecayLength_ljet_ibin%d", (int)mainBinIPs), h3DecayLengthDispersionJetpT_lfjet->GetXaxis()->FindBin(binning[mainBinIPs]), h3DecayLengthDispersionJetpT_lfjet->GetXaxis()->FindBin(binning[mainBinIPs + 1]), 1, h3DecayLengthDispersionJetpT_lfjet->GetZaxis()->FindBin(maxDispersion), "e");

    hDecayLength_bjet->Scale(1.0 / hDecayLength_bjet->Integral());
    hDecayLength_cjet->Scale(1.0 / hDecayLength_cjet->Integral());
    hDecayLength_ljet->Scale(1.0 / hDecayLength_ljet->Integral());

    hDecayLength_bjet->Rebin(2);
    hDecayLength_cjet->Rebin(2);
    hDecayLength_ljet->Rebin(2);

    hDecayLength_bjet->SetMarkerColor(kRed);
    hDecayLength_bjet->SetLineColor(kRed);
    hDecayLength_cjet->SetMarkerColor(kGreen + 2);
    hDecayLength_cjet->SetLineColor(kGreen + 2);
    hDecayLength_ljet->SetMarkerColor(kBlue);
    hDecayLength_ljet->SetLineColor(kBlue);

    hDecayLength_bjet->SetLineWidth(2);
    hDecayLength_cjet->SetLineWidth(2);
    hDecayLength_ljet->SetLineWidth(2);

    hDecayLength_ljet->GetXaxis()->SetRangeUser(0, 20);
    hDecayLength_ljet->GetXaxis()->SetTitle("Decay Length [cm]");
    hDecayLength_ljet->GetYaxis()->SetTitle("Probability Distribution");
    hDecayLength_ljet->GetYaxis()->SetTitleSize(0.045);
    hDecayLength_ljet->GetYaxis()->SetTitleOffset(1.2);
    hDecayLength_ljet->GetXaxis()->SetTitleSize(0.045);
    hDecayLength_ljet->GetXaxis()->SetTitleOffset(0.9);

    hDecayLength_ljet->Draw();
    hDecayLength_bjet->Draw("same");
    hDecayLength_cjet->Draw("SAME");

    TLegend *legDL = new TLegend(0.57, 0.51, 0.91, 0.71);
    legDL->SetLineWidth(0);
    legDL->AddEntry(hDecayLength_bjet, "beauty jets", "l");
    legDL->AddEntry(hDecayLength_cjet, "charm jets", "l");
    legDL->AddEntry(hDecayLength_ljet, "light-flavor jets", "l");
    legDL->Draw();

    textpp->DrawLatexNDC(0.55, 0.93, "PYTHIA 8, pp,#kern[-0.5]{ }#sqrt{#it{s}} = 13 TeV");
    textpp->DrawLatexNDC(0.55, 0.885, "Charged-particle jets");
    textpp->DrawLatexNDC(0.55, 0.84, "Anti-#it{k}_{T},#kern[-0.5]{ }#it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");
    textpp->DrawLatexNDC(0.55, 0.795, Form("%.f#kern[-0.5]{ }#leq#kern[-0.5]{ }#it{p}_{T,ch jet} < %.f GeV/#it{c}", binning[mainBinDL], binning[mainBinDL + 1]));

    gPad->RedrawAxis();

    if (saveFigs)
    {
        cDL->SaveAs("PaperFigures/DecayLengthDistributions.pdf");
    }

    // Plot CPA distributions for different jet flavors
    TH2D *hCPAvsJetpT = (TH2D *)gDirectory->Get("hCPAvsJetpT");
    TH2D *hCPAvsJetpT_bjet = (TH2D *)gDirectory->Get("hCPAvsJetpT_bjet");
    TH2D *hCPAvsJetpT_cjet = (TH2D *)gDirectory->Get("hCPAvsJetpT_cjet");
    TH2D *hCPAvsJetpT_ljet = (TH2D *)hCPAvsJetpT->Clone("hCPAvsJetpT_ljet");
    hCPAvsJetpT_ljet->Add(hCPAvsJetpT_bjet, -1);
    hCPAvsJetpT_ljet->Add(hCPAvsJetpT_cjet, -1);

    TCanvas *cCPA = new TCanvas("cCPA", "CPA Distributions", 800, 800);
    cCPA->cd();
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetTicks();
    gPad->SetLogy();
    gPad->SetMargin(0.11, 0.01, 0.1, 0.01);

    int mainBinCPA = 4;

    TH1D *hCPA_bjet = hCPAvsJetpT_bjet->ProjectionY(Form("hCPA_bjet_ibin%d", (int)mainBinCPA), hCPAvsJetpT_bjet->GetXaxis()->FindBin(binning[mainBinCPA]), hCPAvsJetpT_bjet->GetXaxis()->FindBin(binning[mainBinCPA + 1]));
    TH1D *hCPA_cjet = hCPAvsJetpT_cjet->ProjectionY(Form("hCPA_cjet_ibin%d", (int)mainBinCPA), hCPAvsJetpT_cjet->GetXaxis()->FindBin(binning[mainBinCPA]), hCPAvsJetpT_cjet->GetXaxis()->FindBin(binning[mainBinCPA + 1]));
    TH1D *hCPA_ljet = hCPAvsJetpT_ljet->ProjectionY(Form("hCPA_ljet_ibin%d", (int)mainBinCPA), hCPAvsJetpT_ljet->GetXaxis()->FindBin(binning[mainBinCPA]), hCPAvsJetpT_ljet->GetXaxis()->FindBin(binning[mainBinCPA + 1]));

    hCPA_bjet->Scale(1.0 / hCPA_bjet->Integral());
    hCPA_cjet->Scale(1.0 / hCPA_cjet->Integral());
    hCPA_ljet->Scale(1.0 / hCPA_ljet->Integral());

    hCPA_bjet->Rebin(2);
    hCPA_cjet->Rebin(2);
    hCPA_ljet->Rebin(2);

    hCPA_bjet->SetLineColor(kRed);
    hCPA_bjet->SetMarkerColor(kRed);
    hCPA_cjet->SetLineColor(kGreen + 2);
    hCPA_cjet->SetMarkerColor(kGreen + 2);
    hCPA_ljet->SetLineColor(kBlue);
    hCPA_ljet->SetMarkerColor(kBlue);

    hCPA_bjet->SetLineWidth(2);
    hCPA_cjet->SetLineWidth(2);
    hCPA_ljet->SetLineWidth(2);

    hCPA_ljet->GetXaxis()->SetTitle("CPA");
    hCPA_ljet->GetYaxis()->SetTitle("Probability Distribution");
    hDecayLength_ljet->GetYaxis()->SetTitleSize(0.045);
    hDecayLength_ljet->GetYaxis()->SetTitleOffset(1.2);
    hDecayLength_ljet->GetXaxis()->SetTitleSize(0.045);
    hDecayLength_ljet->GetXaxis()->SetTitleOffset(0.9);

    hCPA_bjet->Draw("");
    hCPA_ljet->Draw("SAME");
    hCPA_cjet->Draw("SAME");

    TLegend *legCPA = new TLegend(0.15, 0.51, 0.49, 0.71);
    legCPA->SetLineWidth(0);
    legCPA->AddEntry(hCPA_bjet, "beauty jets", "l");
    legCPA->AddEntry(hCPA_cjet, "charm jets", "l");
    legCPA->AddEntry(hCPA_ljet, "light-flavor jets", "l");
    legCPA->Draw();

    textpp->DrawLatexNDC(0.14, 0.93, "PYTHIA 8, pp,#kern[-0.5]{ }#sqrt{#it{s}} = 13 TeV");
    textpp->DrawLatexNDC(0.14, 0.885, "Charged-particle jets");
    textpp->DrawLatexNDC(0.14, 0.84, "Anti-#it{k}_{T},#kern[-0.5]{ }#it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");
    textpp->DrawLatexNDC(0.14, 0.795, Form("%.f#kern[-0.5]{ }#leq#kern[-0.5]{ }#it{p}_{T,ch jet} < %.f GeV/#it{c}", binning[mainBinCPA], binning[mainBinCPA + 1]));

    gPad->RedrawAxis();

    if (saveFigs)
    {
        cCPA->SaveAs("PaperFigures/CPA_Distributions.pdf");
    }

    TCanvas *cComp[binning.size() - 1];

    TLatex *latex = new TLatex();
    latex->SetTextFont(42);
    latex->SetTextSizePixels(25);

    for (size_t ibin = 0; ibin < binning.size() - 1; ibin++)
    {
        cComp[ibin] = new TCanvas(Form("cComp%d", (int)ibin), Form("Purity vs Efficiency for IP %.f - %.f", binning[ibin], binning[ibin + 1]), 800, 800);
        cComp[ibin]->cd();
        // gPad->SetGridx();
        // gPad->SetGridy();
        gPad->SetTicks();
        gPad->SetMargin(0.1, 0.01, 0.1, 0.01);

        hTemp1->Draw();
        TLegend *leg1 = new TLegend(0.7, 0.75, 0.91, 0.96);
        leg1->SetTextSizePixels(26);
        leg1->SetLineWidth(0);

        bjet_taggingGraph[ibin]->SetLineColor(kBlue);
        bjet_taggingGraph_IPs[ibin]->SetLineColor(kRed);
        bjet_taggingGraph_SV[ibin]->SetLineColor(kGreen + 2);

        bjet_taggingGraph[ibin]->Draw("same");
        bjet_taggingGraph_IPs[ibin]->Draw("same");
        bjet_taggingGraph_SV[ibin]->Draw("same");

        leg1->AddEntry(bjet_taggingGraph[ibin], "ML model", "l");
        leg1->AddEntry(bjet_taggingGraph_IPs[ibin], "IP method", "l");
        leg1->AddEntry(bjet_taggingGraph_SV[ibin], "SV method", "l");
        leg1->Draw();

        textpp->DrawLatexNDC(0.14, 0.93, "PYTHIA 8, pp,#kern[-0.5]{ }#sqrt{#it{s}} = 13 TeV");
        textpp->DrawLatexNDC(0.14, 0.885, "Charged-particle jets");
        textpp->DrawLatexNDC(0.14, 0.84, "Anti-#it{k}_{T},#kern[-0.5]{ }#it{R} = 0.4, |#it{#eta}_{jet}| < 0.5");
        textpp->DrawLatexNDC(0.14, 0.795, Form("%.f#kern[-0.5]{ }#leq#kern[-0.5]{ }#it{p}_{T,ch jet} < %.f GeV/#it{c}", binning[ibin], binning[ibin + 1]));

        gPad->RedrawAxis();

        if (saveFigs)
        {
            cComp[ibin]->SaveAs(Form("PaperFigures/PurityEfficiencyComparison%.f_%.f.pdf", binning[ibin], binning[ibin + 1]));
        }
    }
}