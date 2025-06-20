#include "RiceStyle.h"
using namespace std;
void getImage()
{
    
    const char* file = "results.root";
    TFile* input = new TFile(file);
    
    TH1D* hdsigmadt_MC = (TH1D*)input->Get("h_tMC");
    TH1D* hdsigmadt_REC = (TH1D*)input->Get("h_tREC");

    int nbins = hdsigmadt_MC->GetNbinsX();
    double dsigmadt_MC,dsigmadt_REC, tBinWidth,t,b,delta, F_b_MC,F_b_REC, result1=0,result2=0;

    // ************************************
    double t_cut = 0.2;
    // ************************************
    
    double bmin= -12;
    double bmax= 12;
    double noOfBins = 300;
    double hbarc = 0.197;
    
    TH1D* hF_b_MC = new TH1D("hF_b_MC", "F_b_MC", noOfBins, bmin, bmax);
    TH1D* hF_b_REC = new TH1D("hF_b_REC", "F_b_MC", noOfBins, bmin, bmax);


   for (int j=1; j<=noOfBins; j++)
    {
        F_b_MC= 0;F_b_REC=0;
        
        b = hF_b_MC->GetBinCenter(j);
    
        double prefactor = 1/6.28;

        for (int i=1; i<=nbins; i++)
        {
            tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            t = hdsigmadt_MC->GetBinCenter(i); // GeV2
            delta =  sqrt(fabs(t)); // GeV
            
            dsigmadt_MC = hdsigmadt_MC->GetBinContent(i);  //nb/GeV^2
            dsigmadt_REC = hdsigmadt_REC->GetBinContent(i);  //nb/GeV^2
            
            dsigmadt_MC/=1e7; //in fm^2/GeV2
            dsigmadt_REC/=1e7; //in fm^2/GeV2
            
            double bessel=TMath::BesselJ0(b*delta/hbarc); //no units
            double amp_MC = sqrt(dsigmadt_MC);
            double amp_REC = sqrt(dsigmadt_REC);
            
            if(t>t_cut)
                continue;
            if(t>0.0414)  { //1st minima
                amp_MC*=-1;
                amp_REC*=-1;
            }
            if(t>0.135)  {   //2nd minima
                amp_MC*=-1;
                amp_REC*=-1;
            }

            result1 =  amp_MC * bessel * tBinWidth/2 ;
            result2 =  amp_REC * bessel * tBinWidth/2 ;
         
            F_b_MC += result1;
            F_b_REC += result2;
            
        }
        
        F_b_MC*=prefactor;F_b_MC/=hbarc;
        F_b_REC*=prefactor;F_b_REC/=hbarc;

        hF_b_MC->SetBinContent(j, F_b_MC);hF_b_MC->SetBinError(j, F_b_MC*0.001);
        hF_b_REC->SetBinContent(j, F_b_REC);hF_b_REC->SetBinError(j, F_b_MC*0.001);
    }

    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetLogx(0);

    TH1D* base1 = makeHist("base1", "", "b [fm]", "F(b)/#scale[0.6]{#int} F(b) db", 100,-8,8,kBlack);
    base1->GetYaxis()->SetRangeUser(0, 0.17);
    base1->GetXaxis()->SetTitleColor(kBlack);
    
    fixedFontHist1D(base1,1.07,1.3);
    base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
    base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
    base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
    base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
    base1->GetXaxis()->SetNdivisions(4,6,0);
    base1->GetYaxis()->SetNdivisions(4,6,0);
    base1->Draw("");

    hF_b_MC->Scale(1.0 / hF_b_MC->Integral("width"));
    hF_b_REC->Scale(1.0 / hF_b_REC->Integral("width"));
    
    gStyle->SetOptStat(0); // to remove the stat box
    gStyle->SetTitleFontSize(.043); //adjust the histo title size
    gStyle->SetTitleX(.53);
    gStyle->SetTitleY(0.96);
    
    hF_b_MC->SetTitle(" #gamma* + Ca #rightarrow #phi + Ca  ");
    
    hF_b_MC->SetMarkerStyle(21); // filled circle
    hF_b_MC->SetMarkerColor(kBlack);
    hF_b_MC->SetLineColor(kBlack);
    hF_b_MC->SetMarkerSize(0.8);
    hF_b_MC->SetLineWidth(2);
    
    hF_b_REC->SetMarkerStyle(24); // filled circle
    hF_b_REC->SetMarkerColor(kBlue);
    hF_b_REC->SetLineColor(kBlue);
    hF_b_REC->SetMarkerSize(0.8);
    
    hF_b_MC->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");

    hF_b_MC->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC->GetYaxis()->SetTitleOffset(1.5);

    hF_b_MC->Draw("same");
    hF_b_REC->Draw("Psame");

    double max_MC=hF_b_MC->GetMaximum();
    double max_REC=hF_b_REC->GetMaximum();
    double FWHM_MC = 0.;
    double FWHM_REC = 0.;
    for(int i=0;i<hF_b_MC->GetNbinsX();i++){
        if(hF_b_MC->GetBinContent(i+1)>=max_MC/2){
            FWHM_MC=fabs(hF_b_MC->GetBinCenter(i+1));
            break;
        }
    }
    for(int i=0;i<hF_b_REC->GetNbinsX();i++){
        if(hF_b_REC->GetBinContent(i+1)>=max_REC/2){
            FWHM_REC=fabs(hF_b_REC->GetBinCenter(i+1));
            break;
        }
    }

    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.85); // (x1, y1, x2, y2)
    leg->AddEntry(hF_b_MC, "MC t_{max} = 0.25", "P"); // "l" = line
    leg->AddEntry(hF_b_REC, "REC t_{max} = 0.25", "P"); // "l" = line
    leg->SetBorderSize(0);      // No border box
    leg->SetFillStyle(0);       // Transparent background
    leg->SetTextFont(42);       // Nice readable font
    leg->SetTextSize(0.03);     // Optional: adjust size
    leg->Draw("same");
    
    c1->Print("F_b.pdf");

    string rootfile= "spaitialProfile.root";
    TFile *hfile = 0;
    hfile = new TFile(rootfile.c_str(), "RECREATE");
    hF_b_MC->Write();
    hF_b_REC->Write();
    
    hfile->Close();

    cout << "True ion size = " << FWHM_MC << " fm " << endl;
    cout << "Reco ion size = " << FWHM_REC  << " fm " << endl;

    cout << rootfile.c_str() <<" written." << endl;
    cout << "All done. Bye." << endl;

}



