#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;

// This code is to calculate the Fprompt of Rn222 from it's pulse shape//

int main(int argc, char *argv[]) {
  // Making plots nicer
  gROOT->ForceStyle(1);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetLabelSize(0.06, "xyz");
  gStyle->SetTitleSize(0.06, "xyz");
  gStyle->SetTitleOffset(0.9, "x");
  gStyle->SetTitleOffset(1.1, "y");
  gStyle->SetTitleOffset(0.9, "z");
  gStyle->SetStatX(0.8);
  gStyle->SetStatW(0.2);
  gStyle->SetStatY(0.85);
  gStyle->SetStatH(0.1);
  gStyle->SetOptStat(0);
  gStyle->SetHistLineWidth(3);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);

  //////////////////////
  // Reading parameter///
  //////////////////////

  if (argc < 3) {
    cout << "Usage: " << argv[0] << " <LightSim Output 1> <LightSim Output 1>"
         << endl;
    return 1;
  }
  cout << endl;
  cout << endl;
  cout << "Running with: " << argv[1] << " and  " << argv[2] << endl;
  cout << endl;

  string OutputDir = "Output/";
  bool SavePlots = false;

  ////////////////////////
  // Define needed number//
  ////////////////////////

  const int NUM_BINS = 10000; // 10000[ns]/NUM_BINS is the timing resolution
  const double T_MAX = 4;     //[us]
  cout << "Number of bins: " << NUM_BINS << endl;
  cout << "Maximum time: " << T_MAX << " mu s" << endl;
  cout << "Time resolution: " << T_MAX / NUM_BINS * 1000 << " ns" << endl;

  ///////////////////
  // Open root files//
  ///////////////////

  TFile *file_Sample_1 = new TFile(argv[1]);
  TTree *Tree_Sample_1 = (TTree *)file_Sample_1->Get(
      "ScintSim_tree"); // If it works on total photon, it works on vuv

  vector<vector<double>> *total_time_vuv_sample_1 = nullptr;
  Tree_Sample_1->SetBranchAddress("total_time_vuv", &total_time_vuv_sample_1);

  TFile *file_Sample_2 = new TFile(argv[2]);
  TTree *Tree_Sample_2 = (TTree *)file_Sample_2->Get(
      "ScintSim_tree"); // If it works on total photon, it works on vuv

  vector<vector<double>> *total_time_vuv_sample_2 = nullptr;
  Tree_Sample_2->SetBranchAddress("total_time_vuv", &total_time_vuv_sample_2);

  cout << "Files opened" << endl;

  ////////////
  // Analysis//
  ////////////

  int number_events_sample_1 =
      Tree_Sample_1->GetEntries(); // Number of events for 5.6 MeV electrons
  int number_events_sample_2 =
      Tree_Sample_2->GetEntries(); // Number of events for Rn222,for now it's
                                   // all 1 million events

  cout << "Number of events in sample 1: " << number_events_sample_1 << endl;
  cout << "Number of events in sample 2: " << number_events_sample_2 << endl;

  // For each event we want to store the earliest time a photon gets detected
  vector<double> earliest_time_sample_1(number_events_sample_1, 100000000);

  std::cout << "Finding t0 for " << argv[1] << std::endl;
  // Finding t0 over all events!
  for (int i = 0; i < number_events_sample_1; i++) {
    Tree_Sample_1->GetEntry(i);

    for (int j = 0; j < total_time_vuv_sample_1->size();
         j++) { // Loop through all pmts

      for (int k = 0; k < total_time_vuv_sample_1->at(j).size();
           k++) { // Loop thorugh all photons in each pmt

        if (total_time_vuv_sample_1->at(j).at(k) < earliest_time_sample_1[i]) {

          earliest_time_sample_1[i] = total_time_vuv_sample_1->at(j).at(k);
        }
      }
    }
  }

  vector<double> earliest_time_sample_2(number_events_sample_2, 100000000);
  std::cout << "Finding t0 for " << argv[2] << std::endl;
  // Finding t0 over all events!
  for (int i = 0; i < number_events_sample_2; i++) {
    Tree_Sample_2->GetEntry(i);

    for (int j = 0; j < total_time_vuv_sample_2->size();
         j++) { // Loop through all pmts

      for (int k = 0; k < total_time_vuv_sample_2->at(j).size();
           k++) { // Loop thorugh all photons in each pmt

        if (total_time_vuv_sample_2->at(j).at(k) < earliest_time_sample_2[i]) {

          earliest_time_sample_2[i] = total_time_vuv_sample_2->at(j).at(k);
        }
      }
    }
  }

  int nPMT = total_time_vuv_sample_1->size(); // Number of PMTs
  assert(nPMT == total_time_vuv_sample_2->size());

  TH1D *pulse_shape_sample_1 =
      new TH1D("pulse_shape_sample_1", "test", NUM_BINS, 0, T_MAX);
  TH1D *pulse_shape_sample_2 =
      new TH1D("pulse_shape_sample_2", "test", NUM_BINS, 0, T_MAX);

  // This creates the pulseshape for each event, combining all PMTs together
  cout << " Fill t0 corrected Histograms from " << argv[1] << endl;

  TH1D *pulse_shape_sample_1_event =
      new TH1D("pulse_shape_sample_1_event", "pulse_shape_sample_1_event",
               NUM_BINS, 0, T_MAX);
  for (int EventIt = 0; EventIt < number_events_sample_1; EventIt++) {
    pulse_shape_sample_1_event->Reset();
    Tree_Sample_1->GetEntry(EventIt);
    for (int pmtIt = 0; pmtIt < nPMT; pmtIt++) { // Loop through all pmts
      for (int photonIt = 0;
           photonIt < total_time_vuv_sample_1->at(pmtIt).size();
           photonIt++) { // Loop thorugh all photons in each pmt
        pulse_shape_sample_1->Fill(
            total_time_vuv_sample_1->at(pmtIt).at(photonIt) -
            earliest_time_sample_1[EventIt]);
        if (SavePlots)
          pulse_shape_sample_1_event->Fill(
              total_time_vuv_sample_1->at(pmtIt).at(photonIt) -
              earliest_time_sample_1[EventIt]);
      }
    }
    if (SavePlots) {
      TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
      pulse_shape_sample_1_event->SetTitle(
          ("Reconstructed Pulse Shape of event " + to_string(EventIt)).c_str());
      pulse_shape_sample_1_event->GetXaxis()->SetTitle("Time (#mu s)");
      pulse_shape_sample_1_event->GetYaxis()->SetTitle("Counts");
      pulse_shape_sample_1_event->Draw("HIST");
      c1->SaveAs((OutputDir + "EventPulseShapes/pulse_shape_sample_1_event_" +
                  to_string(EventIt) + ".png")
                     .c_str());
      delete c1;
    }
  } // End for loop over first file

  cout << " Fill t0 corrected Histograms from " << argv[2] << endl;
  TH1D *pulse_shape_sample_2_event =
      new TH1D("pulse_shape_sample_2_event", "pulse_shape_sample_2_event",
               NUM_BINS, 0, T_MAX);
  for (int EventIt = 0; EventIt < number_events_sample_2; EventIt++) {
    pulse_shape_sample_2_event->Reset();
    Tree_Sample_2->GetEntry(EventIt);
    for (int pmtIt = 0; pmtIt < nPMT; pmtIt++) { // Loop through all pmts
      for (int photonIt = 0;
           photonIt < total_time_vuv_sample_2->at(pmtIt).size();
           photonIt++) { // Loop thorugh all photons in each pmt
        pulse_shape_sample_2->Fill(
            total_time_vuv_sample_2->at(pmtIt).at(photonIt) -
            earliest_time_sample_2[EventIt]);
        if (SavePlots)
          pulse_shape_sample_2_event->Fill(
              total_time_vuv_sample_1->at(pmtIt).at(photonIt) -
              earliest_time_sample_1[EventIt]);
      }
    }
    if (SavePlots) {
      TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
      pulse_shape_sample_2_event->SetTitle(
          ("Reconstructed Pulse Shape of event " + to_string(EventIt)).c_str());
      pulse_shape_sample_2_event->GetXaxis()->SetTitle("Time (#mu s)");
      pulse_shape_sample_2_event->GetYaxis()->SetTitle("Counts");
      pulse_shape_sample_2_event->Draw("HIST");
      c2->SaveAs((OutputDir + "EventPulseShapes/pulse_shape_sample_2_event_" +
                  to_string(EventIt) + ".png")
                     .c_str());
      delete c2;
    }

  } // End for loop over second file

  double MaxAmp1 = 0;
  double MaxAmp2 = 0;
  for (int i = 0; i < NUM_BINS; i++) {
    if (pulse_shape_sample_1->GetBinContent(i) > MaxAmp1) {
      MaxAmp1 = pulse_shape_sample_1->GetBinContent(i);
    }
    if (pulse_shape_sample_2->GetBinContent(i) > MaxAmp2) {
      MaxAmp2 = pulse_shape_sample_2->GetBinContent(i);
    }
  }

  // Running integral//
  cout << "Looking at running integral..." << endl;
  double diff = 0;
  int max_Bin = 0;
  double biggest_diff = 0;
  vector<double> Qft_1(NUM_BINS, 0);
  vector<double> Qft_2(NUM_BINS, 0);
  TH1D *Qft_hist = new TH1D("Qft_1_hist", "Qft_hist", NUM_BINS, 0, T_MAX);
  for (int BinIt = 0; BinIt < NUM_BINS; BinIt++) {
    Qft_1[BinIt] = ((double)pulse_shape_sample_1->Integral(0, BinIt)) /
                   ((double)pulse_shape_sample_1->Integral(0, NUM_BINS));
    Qft_2[BinIt] = ((double)pulse_shape_sample_2->Integral(0, BinIt)) /
                   ((double)pulse_shape_sample_2->Integral(0, NUM_BINS));
    diff = abs(Qft_1[BinIt] - Qft_2[BinIt]);
    Qft_hist->SetBinContent(BinIt, diff);
    if (diff > biggest_diff) {
      biggest_diff = diff;
      max_Bin = BinIt;
    }

  } // End for loop filling running integral differences

  double biggest_diff_time = (double)max_Bin / (double)NUM_BINS;
  cout << "The biggest difference is at time "
       << biggest_diff_time * T_MAX * 1000 << " ns with" << biggest_diff
       << endl;

  cout << "Calculating Fprompt..." << endl; // debug
  vector<int> prompt_sample_1(number_events_sample_1, 0);
  vector<int> prompt_sample_2(number_events_sample_2, 0);
  vector<int> full_sample_1(number_events_sample_1, 0);
  vector<int> full_sample_2(number_events_sample_2, 0);

  for (int EventIt = 0; EventIt < number_events_sample_1; EventIt++) {
    Tree_Sample_1->GetEntry(EventIt);
    for (int pmtIt = 0; pmtIt < nPMT; pmtIt++) { // Loop through all pmts
      for (int photonIt = 0;
           photonIt < total_time_vuv_sample_1->at(pmtIt).size();
           photonIt++) { // Loop thorugh all photons in each pmt
        if (total_time_vuv_sample_1->at(pmtIt).at(photonIt) -
                earliest_time_sample_1[EventIt] <
            biggest_diff_time) {
          prompt_sample_1[EventIt]++;
        }
        full_sample_1[EventIt]++;
      }
    }
  }

  for (int EventIt = 0; EventIt < number_events_sample_2; EventIt++) {
    Tree_Sample_2->GetEntry(EventIt);
    for (int pmtIt = 0; pmtIt < nPMT; pmtIt++) { // Loop through all pmts
      for (int photonIt = 0;
           photonIt < total_time_vuv_sample_2->at(pmtIt).size();
           photonIt++) { // Loop thorugh all photons in each pmt
        if (total_time_vuv_sample_2->at(pmtIt).at(photonIt) -
                earliest_time_sample_2[EventIt] <
            biggest_diff_time) {
          prompt_sample_2[EventIt]++;
        }
        full_sample_2[EventIt]++;
      }
    }
  }

  TH1D *fprompt_sample_1 = new TH1D("fprompt_sample_1", "test", 50, 0, 1);
  TH1D *fprompt_sample_2 = new TH1D("fprompt_sample_2", "test", 50, 0, 1);

  vector<double> prompt_fraction_sample_1(number_events_sample_1, 0);
  for (int EventIt = 0; EventIt < number_events_sample_1; EventIt++) {
    prompt_fraction_sample_1[EventIt] =
        (prompt_sample_1[EventIt] * 1.0) / (1. * full_sample_1[EventIt]);
    fprompt_sample_1->Fill(prompt_fraction_sample_1[EventIt]);
  }

  vector<double> prompt_fraction_sample_2(number_events_sample_2, 0);
  for (int EventIt = 0; EventIt < number_events_sample_2; EventIt++) {
    prompt_fraction_sample_2[EventIt] =
        (prompt_sample_2[EventIt] * 1.0) / (1. * full_sample_2[EventIt]);
    fprompt_sample_2->Fill(prompt_fraction_sample_2[EventIt]);
  }

  /* //////// */
  /* //Draw// */
  /* //////// */

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);
  TLegend *leg_sample = new TLegend(0.25, 0.7, 0.45, 0.93);
  pulse_shape_sample_1->SetTitle("Pulse Shape for all PMTs");
  pulse_shape_sample_1->GetXaxis()->SetTitle("time (#mus)");
  /* pulse_shape_sample_2->Scale(1./pulse_shape_sample_2->Integral()); */
  /* pulse_shape_sample_1->Scale(1./pulse_shape_sample_1->Integral()); */
  /* pulse_shape_sample_2->Scale(0.035/MaxAmp2); */
  /* pulse_shape_sample_1->Scale(0.035/MaxAmp1); */
  /* pulse_shape_sample_1->GetXaxis()->SetRangeUser(0,0.4);*/
  pulse_shape_sample_1->SetLineColor(1);
  pulse_shape_sample_2->SetLineColor(2);
  pulse_shape_sample_1->SetLineWidth(3);
  pulse_shape_sample_2->SetLineWidth(3);
  leg_sample->AddEntry(pulse_shape_sample_1, "Sample 1");
  leg_sample->AddEntry(pulse_shape_sample_2, "Sample 2");
  pulse_shape_sample_1->Draw("HIST");
  pulse_shape_sample_2->Draw("SAME HIST");
  leg_sample->Draw("SAME");
  c1->SetLogy();
  c1->SaveAs(Form("%spulse_shape.png", OutputDir.c_str()));
  delete c1;

  TCanvas *c2 = new TCanvas("c3", "c3", 1000, 800);
  Qft_hist->SetTitle("Qft");
  Qft_hist->GetXaxis()->SetTitle("time");
  Qft_hist->GetXaxis()->SetRangeUser(0, 0.4);
  Qft_hist->SetLineColor(1);
  Qft_hist->SetLineWidth(3);
  Qft_hist->Draw("HIST");
  c2->SaveAs(Form("%sQft_diff.png", OutputDir.c_str()));
  delete c2;

  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 800);
  TLegend *leg_fPrompt = new TLegend(0.25, 0.7, 0.45, 0.93);
  fprompt_sample_1->SetLineColor(1);
  fprompt_sample_1->SetLineWidth(3);
  fprompt_sample_1->SetTitle("fPrompt distro");
  fprompt_sample_1->GetXaxis()->SetTitle("fPrompt");
  fprompt_sample_1->GetYaxis()->SetTitle("fraction of events");
  fprompt_sample_2->SetLineColor(2);
  fprompt_sample_2->SetLineWidth(3);
  // fprompt_sample_1->Scale(1.0/(double)fprompt_sample_1->Integral());
  // fprompt_sample_2->Scale(1.0/(double)fprompt_sample_2->Integral());
  leg_fPrompt->AddEntry(fprompt_sample_1, "Sample 1");
  leg_fPrompt->AddEntry(fprompt_sample_2, "Sample 2");
  fprompt_sample_1->Draw("HIST");
  fprompt_sample_2->Draw("SAME HIST");
  leg_fPrompt->Draw("SAME");
  c3->SaveAs(Form("%sfrpompt.png", OutputDir.c_str()));
  delete c3;
  delete leg_fPrompt;

  return 0;
}
