#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

void DrawTTreeVariables()
{
    bool drawFeatures = false;
    bool drawEff = true;
    gStyle->SetOptStat(0);

    // Open the ROOT file and retrieve the TTree
    TFile *file = TFile::Open("ML_bjets/Data/MergedTTree_pTHardBins.root");
    // TFile *file = TFile::Open("HFTagging.root");
    // TFile *file = TFile::Open("ML_bjets/Data/MergedTTree_MB.root");
    // TFile *file = TFile::Open("ML_bjets/Data/MergedTTree_MB_Large.root");
    TTree *tree = (TTree *)file->Get("jetTree");

    // Jet variables to plot
    const char *variables[] = {"jetPt", "jetEta", "jetPhi", "jetMass", "nTracks", "nSV"};

    std::vector<std::array<float, 3>> binningJet = {{200, 0, 700}, {100, -1.2, 1.2}, {100, 0, 3 * M_PI}, {200, 0, 200}, {100, 0, 100}, {200, 0, 5000}};
    std::vector<std::array<float, 3>> binningTrack = {{200, 0, 200}, {100, -1.2, 1.2}, {200, 0, 1000}, {200, 0, 1000}, {100, 0, 1}, {100, -2, 2}, {100, -2, 2}, {200, 0, 1.2}, {200, 0, 1}};
    std::vector<std::array<float, 3>> binningSV = {{200, 0, 200}, {200, 0, 1}, {200, 0, 20}, {200, 0, 1.2}, {200, -1, 1}, {200, -1.2, 1.2}, {200, 0, 5}, {200, 0, 5}, {200, 0, 30}, {200, 0, 30}};

    // Track variables to plot
    const char *trackVariables[] = {"trackPt", "trackEta", "dotProdTrackJet", "dotProdTrackJetOverJet", "deltaRJetTrack", "signedIP2D", "signedIP3D", "momFraction", "deltaRTrackVertex"};

    // SV variables to plot
    const char *svVariables[] = {"svPt", "deltaRSVJet", "svMass", "svfE", "svIPxy", "svCPA", "svChi2PCA", "svDispersion", "svDecayLength2D", "svDecayLength3D"};

    const int nVars = sizeof(variables) / sizeof(variables[0]);
    const int nTrackVars = sizeof(trackVariables) / sizeof(trackVariables[0]);
    const int nSVVars = sizeof(svVariables) / sizeof(svVariables[0]);

    // Jet flavors to compare (0 = light, 1 = c, 2 = b jets)
    int flavors[] = {0, 1, 2};
    std::string flavorNames[] = {"lf-jets", "c-jets", "b-jets"};
    const int nFlavors = sizeof(flavors) / sizeof(flavors[0]);

    // Define colors for each flavor
    int colors[] = {kBlue, kGreen, kRed};

    // Helper function to draw variables
    auto drawVariables = [&](const char *varList[], int nVars, std::vector<std::array<float, 3>> binning)
    {
        for (int i = 0; i < nVars; ++i)
        {
            TCanvas *c = new TCanvas(varList[i], varList[i], 800, 600);
            c->SetLogy(); // Set log scale on Y-axis
            TLegend *leg = new TLegend(0.6, 0.7, 0.8, 0.9);

            for (int j = 0; j < nFlavors; ++j)
            {
                TString hname = Form("h_%s_flavor%d", varList[i], flavors[j]);
                TH1F *h = new TH1F(hname, varList[i], binning[i][0], binning[i][1], binning[i][2]);

                TString drawCmd = Form("%s >> %s", varList[i], hname.Data());
                TString cut = Form("jetFlavor == %d && jetPt > 10  && jetPt < 200 && trackPt > 0 && svPt > 0", flavors[j]);

                tree->Draw(drawCmd, cut, "same");

                h->SetLineColor(colors[j]);
                h->SetLineWidth(2);

                h->Scale(1.0 / h->Integral());

                if (j == 0)
                    h->Draw();
                else
                    h->Draw("SAME");

                leg->AddEntry(h, Form("%s", flavorNames[j].c_str()), "l");
            }

            leg->Draw();
            c->SaveAs(Form("figures/TreeFigures/%s_comparison.png", varList[i]));
        }
    };

    // Draw all variables
    if (drawFeatures)
    {
        drawVariables(variables, nVars, binningJet);
        drawVariables(trackVariables, nTrackVars, binningTrack);
        drawVariables(svVariables, nSVVars, binningSV);
    }

    // Optional block to create histograms, apply cuts, and draw ratios
    auto drawRatios = [&](const char *cutCondition, const char *cutName)
    {
        TCanvas *c = new TCanvas(Form("c_%s", cutName), Form("Jet pT Ratios with %s Cut", cutName), 800, 600);
        TLegend *leg = new TLegend(0.6, 0.7, 0.8, 0.9);

        for (int j = 0; j < nFlavors; ++j)
        {
            TString hnameCut = Form("h_jetPt_%s_cut_flavor%d", cutName, flavors[j]);
            TString hnameNoCut = Form("h_jetPt_no_cut_flavor%d", flavors[j]);

            TH1F *hCut = new TH1F(hnameCut, "Jet pT", 200, 0, 700);
            TH1F *hNoCut = new TH1F(hnameNoCut, "Jet pT", 200, 0, 700);

            TString drawCmd = Form("jetPt >> %s", hnameCut.Data());
            TString cut = Form("jetFlavor == %d && %s", flavors[j], cutCondition);
            tree->Draw(drawCmd, cut, "same");

            drawCmd = Form("jetPt >> %s", hnameNoCut.Data());
            cut = Form("jetFlavor == %d", flavors[j]);
            tree->Draw(drawCmd, cut, "same");

            hCut->SetLineColor(colors[j]);
            hCut->SetLineWidth(2);
            hNoCut->SetLineColor(colors[j]);
            hNoCut->SetLineWidth(2);
            hNoCut->SetLineStyle(2);

            TH1F *hRatio = (TH1F *)hCut->Clone(Form("h_ratio_%s_flavor%d", cutName, flavors[j]));
            hRatio->Divide(hNoCut);

            hRatio->SetLineColor(colors[j]);
            hRatio->SetLineWidth(2);

            if (j == 0)
                hRatio->Draw();
            else
                hRatio->Draw("SAME");

            leg->AddEntry(hRatio, Form("%s", flavorNames[j].c_str()), "l");
        }

        leg->Draw();
        c->SaveAs(Form("figures/TreeFigures/jetPt_ratio_%s.png", cutName));
    };

    // Draw ratios with cuts
    if (drawEff)
    {
        drawRatios("svDecayLength3D[0] > 0.5 && svDispersion[0] < 0.03", "svDecayLength3D");
        drawRatios("TMath::Sqrt(signedIP2D[2]*signedIP2D[2] + signedIP3D[2]*signedIP3D[2]) > 0.008", "signedIP3D");
        drawRatios("signedIP2D[1] > 0.008", "signedIP2D");
    }

    // file->Close();
}
