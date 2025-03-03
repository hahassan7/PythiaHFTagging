#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void DrawTTreeVariables()
{
    // Open the ROOT file and retrieve the TTree
    // TFile *file = TFile::Open("ML_bjets/Data/MergedTTree_pTHardBins.root");
    // TFile *file = TFile::Open("ML_bjets/Data/MergedTTree_MB.root");
    TFile *file = TFile::Open("ML_bjets/Data/MergedTTree_MB_Large.root");
    TTree *tree = (TTree *)file->Get("jetTree");

    // Jet variables to plot
    const char *variables[] = {"jetPt", "jetEta", "jetPhi", "jetMass", "nTracks", "nSV"};

    std::vector<std::array<float, 3>> binningJet = {{200, 0, 700}, {100, -1.2, 1.2}, {100, 0, 3 * M_PI}, {200, 0, 200}, {100, 0, 100}, {200, 0, 5000}};
    std::vector<std::array<float, 3>> binningTrack = {{200, 0, 200}, {100, -1.2, 1.2}, {200, 0, 1000}, {200, 0, 1000}, {100, 0, 1}, {200, -5, 5}, {200, -5, 5}, {200, 0, 1.2}, {200, 0, 1}};
    std::vector<std::array<float, 3>> binningSV = {{200, 0, 200}, {200, 0, 1}, {200, 0, 20}, {200, 0, 1.2}, {200, -1, 1}, {200, -1.2, 1.2}, {200, 0, 10}, {200, 0, 10}, {200, 0, 50}, {200, 0, 50}};

    // Track variables to plot
    const char *trackVariables[] = {"trackPt", "trackEta", "dotProdTrackJet", "dotProdTrackJetOverJet", "deltaRJetTrack", "signedIP2D", "signedIP3D", "momFraction", "deltaRTrackVertex"};

    // SV variables to plot
    const char *svVariables[] = {"svPt", "deltaRSVJet", "svMass", "svfE", "svIPxy", "svCPA", "svChi2", "svDispersion", "svDecayLength2D", "svDecayLength3D"};

    const int nVars = sizeof(variables) / sizeof(variables[0]);
    const int nTrackVars = sizeof(trackVariables) / sizeof(trackVariables[0]);
    const int nSVVars = sizeof(svVariables) / sizeof(svVariables[0]);

    // Jet flavors to compare (0 = light, 1 = c, 2 = b jets)
    int flavors[] = {0, 1, 2};
    const int nFlavors = sizeof(flavors) / sizeof(flavors[0]);

    // Define colors for each flavor
    int colors[] = {kRed, kBlue, kGreen};

    // Helper function to draw variables
    auto drawVariables = [&](const char *varList[], int nVars, std::vector<std::array<float, 3>> binning)
    {
        for (int i = 0; i < nVars; ++i)
        {
            TCanvas *c = new TCanvas(varList[i], varList[i], 800, 600);
            c->SetLogy(); // Set log scale on Y-axis
            TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

            for (int j = 0; j < nFlavors; ++j)
            {
                TString hname = Form("h_%s_flavor%d", varList[i], flavors[j]);
                TH1F *h = new TH1F(hname, varList[i], binning[i][0], binning[i][1], binning[i][2]);

                TString drawCmd = Form("%s >> %s", varList[i], hname.Data());
                TString cut = Form("jetFlavor == %d", flavors[j]);

                tree->Draw(drawCmd, cut, "same");

                h->SetLineColor(colors[j]);
                h->SetLineWidth(2);

                h->Scale(1.0 / h->Integral());

                if (j == 0)
                    h->Draw();
                else
                    h->Draw("SAME");

                leg->AddEntry(h, Form("Flavor %d", flavors[j]), "l");
            }

            leg->Draw();
            c->SaveAs(Form("figures/TreeFigures/%s_comparison.png", varList[i]));
        }
    };

    // Draw all variables
    drawVariables(variables, nVars, binningJet);
    drawVariables(trackVariables, nTrackVars, binningTrack);
    drawVariables(svVariables, nSVVars, binningSV);

    file->Close();
}
