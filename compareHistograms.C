#include "TSystem.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TKey.h"
#include <iostream>
#include <iomanip>
#include <vector>

void compareHistograms()
{
    // Create figures directory if it doesn't exist
    gSystem->Exec("mkdir -p figures");

    // Open both files
    TFile *file1 = TFile::Open("AnalysisResults.root", "READ");
    TFile *file2 = TFile::Open("AnalysisResultsFromTree_AFter.root", "READ");

    if (!file1 || !file2)
    {
        std::cerr << "Error: Could not open one or both files" << std::endl;
        return;
    }

    // Create output file for comparison plots
    TFile *outFile = new TFile("figures/HistogramComparisons.root", "RECREATE");

    // Get list of keys (histograms) from first file
    TList *keys = file1->GetListOfKeys();

    // Create a canvas for each comparison
    TCanvas *canvas = new TCanvas("canvas", "Histogram Comparisons", 1200, 800);
    canvas->Divide(2, 1);

    // Create vectors to store comparison results
    std::vector<TString> histNames;
    std::vector<double> chi2Values;
    std::vector<double> ksValues;
    std::vector<double> meanDiffs;
    std::vector<double> rmsDiffs;

    // Loop over all histograms
    for (int i = 0; i < keys->GetSize(); i++)
    {
        TKey *key = (TKey *)keys->At(i);
        TString histName = key->GetName();

        // Get histograms from both files
        TH1 *hist1 = (TH1 *)file1->Get(histName);
        TH1 *hist2 = (TH1 *)file2->Get(histName);

        if (!hist1 || !hist2)
        {
            std::cout << "Warning: Histogram " << histName << " not found in both files" << std::endl;
            continue;
        }

        // Set different colors and styles
        hist1->SetLineColor(kBlue);
        hist2->SetLineColor(kRed);
        hist1->SetMarkerStyle(20);
        hist2->SetMarkerStyle(21);
        hist1->SetMarkerColor(kBlue);
        hist2->SetMarkerColor(kRed);

        // Calculate statistical comparison
        double chi2 = hist1->Chi2Test(hist2, "UW CHI2/NDF");
        double ks = hist1->KolmogorovTest(hist2);

        // Draw histograms with error bars
        canvas->cd(1);
        hist1->Draw("E");
        hist2->Draw("E SAME");

        // Add legend
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hist1, "Original", "lep");
        legend->AddEntry(hist2, "BTree", "lep");
        legend->Draw();

        // Add statistics
        TPaveText *stats = new TPaveText(0.7, 0.5, 0.9, 0.7, "NDC");
        stats->AddText(Form("Chi2/NDF: %.3f", chi2));
        stats->AddText(Form("KS test: %.3f", ks));
        stats->Draw();

        // Draw ratio with proper error propagation
        canvas->cd(2);
        TH1 *ratio = (TH1 *)hist1->Clone("ratio");
        ratio->Divide(hist1, hist2, 1.0, 1.0, "B"); // "B" option for binomial errors
        ratio->SetTitle(Form("%s Ratio", histName.Data()));
        ratio->GetYaxis()->SetTitle("Original/BTree");
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerColor(kBlack);
        ratio->SetLineColor(kBlack);

        // Set reasonable Y-axis range for ratio plot
        double maxRatio = ratio->GetMaximum();
        double minRatio = ratio->GetMinimum();
        double padding = 0.5; // 50% padding
        ratio->GetYaxis()->SetRangeUser(std::max(0.0, minRatio - padding), maxRatio + padding);

        ratio->Draw("E");

        // Draw reference line at 1
        TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 1,
                                ratio->GetXaxis()->GetXmax(), 1);
        line->SetLineStyle(2);
        line->Draw();

        // Save canvas to figures directory
        canvas->SaveAs(Form("figures/comparison_%s.png", histName.Data()));
        canvas->Write(Form("comparison_%s", histName.Data()));

        // Store results for summary
        histNames.push_back(histName);
        chi2Values.push_back(chi2);
        ksValues.push_back(ks);
        meanDiffs.push_back(TMath::Abs(hist1->GetMean() - hist2->GetMean()));
        rmsDiffs.push_back(TMath::Abs(hist1->GetRMS() - hist2->GetRMS()));
    }

    // Print summary table
    std::cout << "\n=== Summary Table ===\n"
              << std::endl;
    std::cout << std::setw(30) << "Histogram" << " | "
              << std::setw(12) << "Chi2/NDF" << " | "
              << std::setw(12) << "KS test" << " | "
              << std::setw(12) << "Mean Diff" << " | "
              << std::setw(12) << "RMS Diff" << std::endl;
    std::cout << std::string(85, '-') << std::endl;

    for (size_t i = 0; i < histNames.size(); i++)
    {
        std::cout << std::setw(30) << histNames[i].Data() << " | "
                  << std::setw(12) << std::fixed << std::setprecision(3) << chi2Values[i] << " | "
                  << std::setw(12) << ksValues[i] << " | "
                  << std::setw(12) << meanDiffs[i] << " | "
                  << std::setw(12) << rmsDiffs[i] << std::endl;
    }

    // Calculate and print statistics
    double avgChi2 = TMath::Mean(chi2Values.begin(), chi2Values.end());
    double stdChi2 = TMath::RMS(chi2Values.begin(), chi2Values.end());
    double avgKS = TMath::Mean(ksValues.begin(), ksValues.end());
    double avgMeanDiff = TMath::Mean(meanDiffs.begin(), meanDiffs.end());
    double avgRMSDiff = TMath::Mean(rmsDiffs.begin(), rmsDiffs.end());

    std::cout << std::string(85, '-') << std::endl;
    std::cout << std::setw(30) << "Average" << " | "
              << std::setw(12) << avgChi2 << " | "
              << std::setw(12) << avgKS << " | "
              << std::setw(12) << avgMeanDiff << " | "
              << std::setw(12) << avgRMSDiff << std::endl;
    std::cout << std::setw(30) << "Std Dev (Chi2)" << " | "
              << std::setw(12) << stdChi2 << std::endl;

    outFile->Close();
    file1->Close();
    file2->Close();

    delete canvas;
}