#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph2DErrors.h"
#include "TLegend.h"
#include <TStyle.h>
#include <TROOT.h>


void ImportRadcor(vector<double> &v1,vector<double> &v2,vector<double> &v3,vector<double> &v4,vector<double> &v5,vector<double> &v6,vector<double> &v7,vector<double> &v8,vector<double> &v9, vector<double> &v10,vector<double> &v11,vector<double> &v12,vector<double> &v13,vector<double> &v14,vector<double> &v15,vector<double> &v16, const char * filename){


   double i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16;
   ifstream infile(filename);

   	if(infile.fail()){
      		cout << "Cannot open the file: " << filename << endl;
      		exit(1);
   	}
	else{
         	while(!infile.eof()){
         		infile >> i1 >> i2 >> i3 >> i4>>i5>>i6>>i7>>i8>>i9>>i10>>i11>>i12>>i13>>i14>>i15>>i16;
         		v1.push_back (i1);
         		v2.push_back (i2);
         		v3.push_back (i3);
         		v4.push_back (i4);
         		v5.push_back (i5);
         		v6.push_back (i6);
         		v7.push_back (i7);
         		v8.push_back (i8);
         		v9.push_back (i9);
         		v10.push_back(i10);
         		v11.push_back(i11);
        		v12.push_back(i12);
         		v13.push_back(i13);
         		v14.push_back(i14);
         		v15.push_back(i15);
         		v16.push_back(i16);
         	}

      		infile.close();
      		v1.pop_back();
      		v2.pop_back();
      		v3.pop_back();
      		v4.pop_back();
      		v5.pop_back();
      		v6.pop_back();
      		v7.pop_back();
      		v8.pop_back();
      		v9.pop_back();
      		v10.pop_back();
      		v11.pop_back();
      		v12.pop_back();
      		v13.pop_back();
      		v14.pop_back();
      		v15.pop_back();
      		v16.pop_back();
         }
}


void SHMSacc(TString mcfile, TString datfile){

TString inputroot, replayroot;
inputroot = mcfile + ".root";
replayroot = datfile + ".root";
TFile *f1 = new TFile(inputroot); 
TTree *t1 = (TTree*) f1->Get("h1411");
TFile* f = TFile::Open(replayroot);


float psxptar; //reconstructed target dx/dz vertical slope
float psyptar; //reconstructed target dy/dz vertical slope
float psdelta; //reconstructed 100*(p - pc)/pc with pc = central SHMS momentum
float psxpfp; // focal plane dx/dz vertical slope
float psypfp; // focal plane dy/dz vertical slope
float psxfp; // fp x pos
float psyfp; // fp y pos
float normfac;


t1->SetBranchAddress("psdelta",&psdelta);
t1->SetBranchAddress("stop_id",&normfac); 
t1->SetBranchAddress("psxpfp",&psxpfp);
t1->SetBranchAddress("psypfp",&psypfp);
t1->SetBranchAddress("psxfp",&psxfp);
t1->SetBranchAddress("psyfp",&psyfp);
t1->SetBranchAddress("psxptar",&psxptar);
t1->SetBranchAddress("psyptar",&psyptar);



double Z = 6.0;
double A = 12.0;
double Ei = 2.218; //Beam energy //GeV
double Mp = 0.93825; //proton's mass //GeV
double N_A = 6.02*1e+23; // Avogadro's number
double Q_E = 1.60*1e-19; // Conversion from charge to electron
double p_spec = 1.6; //spectrometer momentum //GeV 
double deg2rad = 3.14159/180.;
double ts = 15.; //central spectrometer angle //deg
double cos_ts = cos(ts * deg2rad);
double sin_ts = sin(ts * deg2rad);
double car_density = 2.2; // density of carbon g/cm3
double mass_tar = 12. * 931.5; //mass of the target
double cur = 2.45; //current uA
double thick = 0.1749; // target thickness g/cm2
double dpp_up = 40.0; // momentum acceptance upper limit
double dpp_down = -15.0; // momentum acceptance lower limit
double phi_up = 100.0; //mrad
double phi_down = -100.0; //mrad
double theta_up = 100.0; //mrad
double theta_down = -100.0; //mrad
double domega = (phi_up - phi_down)*(theta_up-theta_down) / 1000. / 1000.; // diff solid angle in sr
double run_time = 481.905403;// sec;
double live_time = 52.45 / 100.0;
double lumin = thick*cur/12.0*N_A/Q_E*1e-36*run_time*live_time;// luminosity per ub at 5.5uA
double ep_min = p_spec*(1.+0.01*dpp_down); //GeV
double ep_max = p_spec*(1.+0.01*dpp_up); //GeV

vector<double>   V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16;

   //calling the function
   ImportRadcor(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,Form("radiated.out"));        

   const int size= V1.size();
   TGraph2DErrors *gr2D = new TGraph2DErrors(size); 

   	for(Int_t j =0;j<size;j++){

       	   gr2D->SetPoint(j,V2[j],V3[j],V7[j]); 
           gr2D->SetPointError(j,0.,0.,0.); 

        }


//define MC histograms
TH1F *hW = new TH1F("hW", "hW", 50, 1.8, 2.8);
hW->GetXaxis()->SetTitle("W [GeV]");
hW->GetYaxis()->SetTitle("counts");
TH1F *hQ2 = new TH1F("hQ2", "hQ2", 50, 0.8, 2.0);
hQ2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
hQ2->GetYaxis()->SetTitle("counts");
TH1F *hom = new TH1F("hom", "hom",50, 2.2, 4.0);
hom->GetXaxis()->SetTitle("#nu [GeV]");
hom->GetYaxis()->SetTitle("counts");
TH1F *hth = new TH1F("hth", "hth",50, 13.0, 17.0);
hth->GetXaxis()->SetTitle("theta [deg]");
hth->GetYaxis()->SetTitle("counts");
//1D focal plane histograms
TH1F *xpfp = new TH1F("xpfp", "xpfp",100, -0.1, 0.1);
xpfp->GetXaxis()->SetNdivisions(4);
xpfp->GetXaxis()->SetTitle("Xpfp(rad)");
xpfp->GetYaxis()->SetTitle("counts");
TH1F *ypfp = new TH1F("ypfp", "ypfp",100, -0.05, 0.05);
ypfp->GetXaxis()->SetNdivisions(4);
ypfp->GetXaxis()->SetTitle("Ypfp(rad)");
ypfp->GetYaxis()->SetTitle("counts");
TH1F *xfp = new TH1F("xfp", "xfp",100, -40, 40);
xfp->GetXaxis()->SetTitle("Xfp(cm)");
xfp->GetYaxis()->SetTitle("counts");
TH1F *yfp = new TH1F("yfp", "yfp",100, -40, 40);
yfp->GetXaxis()->SetTitle("Yfp(cm)");
yfp->GetYaxis()->SetTitle("counts");
TH1F *delta = new TH1F("delta", "delta",100, -20, 40);
delta->GetXaxis()->SetTitle("dp/p(%)");
delta->GetYaxis()->SetTitle("counts");
//2D focal plane histograms
TH2F *thvsX = new TH2F("thvsX", "thvsX",100, -0.1, 0.1, 100,-40, 40);
thvsX->GetXaxis()->SetNdivisions(4);
thvsX->GetXaxis()->SetTitle("Ypfp(rad)");
thvsX->GetYaxis()->SetTitle("Xfp(cm)");
TH2F *phivsX = new TH2F("phivsX", "phivsX",100, -0.1, 0.1, 100,-40, 40);
phivsX->GetXaxis()->SetNdivisions(4);
phivsX->GetXaxis()->SetTitle("Xpfp(rad)");
phivsX->GetYaxis()->SetTitle("Xfp(cm)");
TH2F *phivsY = new TH2F("phivsY", "phivsY",100, -0.1, 0.1, 100,-40, 40);
phivsY->GetXaxis()->SetNdivisions(4);
phivsY->GetXaxis()->SetTitle("Xpfp(rad)");
phivsY->GetYaxis()->SetTitle("Yfp(cm)");
TH2F *thvsY = new TH2F("thvsY", "thvsY",100, -0.1, 0.1, 100,-40, 40);
thvsY->GetXaxis()->SetNdivisions(4);
thvsY->GetXaxis()->SetTitle("Ypfp(rad)");
thvsY->GetYaxis()->SetTitle("Yfp(cm)");
TH2F *thvsphifp = new TH2F("thvsphifp", "thvsphifp", 100,-0.1, 0.1, 100, -0.1, 0.1);
thvsphifp->GetXaxis()->SetNdivisions(4);
thvsphifp->GetXaxis()->SetTitle("Ypfp(rad)");
thvsphifp->GetYaxis()->SetTitle("Xpfp(rad)");
TH2F *XvsY = new TH2F("XvsY", "XvsY", 100,-40, 40,100, -40, 40);
XvsY->GetXaxis()->SetTitle("Xfp(cm)");
XvsY->GetYaxis()->SetTitle("Yfp(cm)");

//same histograms with different names to use in ratio plots
TH1F *hWr = new TH1F("hWr", "hWr", 50, 1.8, 2.8);
hWr->GetXaxis()->SetTitle("W [GeV]");
hWr->GetYaxis()->SetTitle("counts");
TH1F *hQ2r = new TH1F("hQ2r", "hQ2r", 50, 0.8, 2.0);
hQ2r->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
hQ2r->GetYaxis()->SetTitle("counts");
TH1F *homr = new TH1F("homr", "homr",50, 2.2, 4.0);
homr->GetXaxis()->SetTitle("#nu [GeV]");
homr->GetYaxis()->SetTitle("counts");
TH1F *hthr = new TH1F("hthr", "hthr",50, 13.0, 17.0);
hthr->GetXaxis()->SetTitle("theta [deg]");
hthr->GetYaxis()->SetTitle("counts");
//1D focal plane histograms
TH1F *xpfpr = new TH1F("xpfpr", "xpfpr",100, -0.1, 0.1);
xpfpr->GetXaxis()->SetNdivisions(4);
xpfpr->GetXaxis()->SetTitle("Xpfp(rad)");
xpfpr->GetYaxis()->SetTitle("counts");
TH1F *ypfpr = new TH1F("ypfpr", "ypfpr",100, -0.05, 0.05);
ypfpr->GetXaxis()->SetNdivisions(4);
ypfpr->GetXaxis()->SetTitle("Ypfp(rad)");
ypfpr->GetYaxis()->SetTitle("counts");
TH1F *xfpr = new TH1F("xfpr", "xfpr",100, -40, 40);
xfpr->GetXaxis()->SetTitle("Xfp(cm)");
xfpr->GetYaxis()->SetTitle("counts");
TH1F *yfpr = new TH1F("yfpr", "yfpr",100, -40, 40);
yfpr->GetXaxis()->SetTitle("Yfp(cm)");
yfpr->GetYaxis()->SetTitle("counts");
TH1F *deltar = new TH1F("deltar", "deltar",100, -20, 40);
deltar->GetXaxis()->SetTitle("dp/p(%)");
deltar->GetYaxis()->SetTitle("counts");
//2D focal plane histograms
TH2F *thvsXr = new TH2F("thvsXr", "thvsXr",100, -0.1, 0.1, 100,-40, 40);
thvsXr->GetXaxis()->SetNdivisions(4);
thvsXr->GetXaxis()->SetTitle("Ypfp(rad)");
thvsXr->GetYaxis()->SetTitle("Xfp(cm)");
TH2F *phivsXr = new TH2F("phivsXr", "phivsXr",100, -0.1, 0.1, 100,-40, 40);
phivsXr->GetXaxis()->SetNdivisions(4);
phivsXr->GetXaxis()->SetTitle("Xpfp(rad)");
phivsXr->GetYaxis()->SetTitle("Xfp(cm)");
TH2F *phivsYr = new TH2F("phivsYr", "phivsYr",100, -0.1, 0.1, 100,-40, 40);
phivsYr->GetXaxis()->SetNdivisions(4);
phivsYr->GetXaxis()->SetTitle("Xpfp(rad)");
phivsYr->GetYaxis()->SetTitle("Yfp(cm)");
TH2F *thvsYr = new TH2F("thvsYr", "thvsYr",100, -0.1, 0.1, 100,-40, 40);
thvsYr->GetXaxis()->SetNdivisions(4);
thvsYr->GetXaxis()->SetTitle("Ypfp(rad)");
thvsYr->GetYaxis()->SetTitle("Yfp(cm)");
TH2F *thvsphifpr = new TH2F("thvsphifpr", "thvsphifpr", 100,-0.1, 0.1, 100, -0.1, 0.1);
thvsphifpr->GetXaxis()->SetNdivisions(4);
thvsphifpr->GetXaxis()->SetTitle("Ypfp(rad)");
thvsphifpr->GetYaxis()->SetTitle("Xpfp(rad)");
TH2F *XvsYr = new TH2F("XvsYr", "XvsYr", 100,-40, 40,100, -40, 40);
XvsYr->GetXaxis()->SetTitle("Xfp(cm)");
XvsYr->GetYaxis()->SetTitle("Yfp(cm)");

//define replayed data histograms
//kinematics
TH1F *hQ2_dat = (TH1F*)f->Get("pkin_q2_elec");
TH1F *hW_dat = (TH1F*)f->Get("pkin_w_elec");
TH1F *hom_dat = (TH1F*)f->Get("pkin_omega_elec");
TH1F *hth_dat = (TH1F*)f->Get("pkin_theta_elec");
TH1F *hw2_dat = (TH1F*)f->Get("pkin_w2_elec");
//focal plane 1D
TH1F *xfpD = (TH1F*)f->Get("pdc_xfp");
xfpD->GetXaxis()->SetTitle("Xfp(cm)");
TH1F *xpfpD = (TH1F*)f->Get("pdc_xpfp");
xpfpD->GetXaxis()->SetTitle("Xpfp(rad)");
TH1F *yfpD = (TH1F*)f->Get("pdc_yfp");
yfpD->GetXaxis()->SetTitle("Yfp(cm)");
TH1F *ypfpD = (TH1F*)f->Get("pdc_ypfp");
ypfpD->GetXaxis()->SetTitle("Ypfp(rad)");
//focal plane 2D
TH2F *XvsYd= (TH2F*)f->Get("pdc_xfp_vs_yfp");
XvsYd->GetXaxis()->SetTitle("Xfp(cm)");
XvsYd->GetYaxis()->SetTitle("Yfp(cm)");
TH2F *thvsphid= (TH2F*)f->Get("pdc_xpfp_vs_ypfp");
thvsphid->GetXaxis()->SetTitle("Ypfp(rad)");
thvsphid->GetYaxis()->SetTitle("Xpfp(rad)");
TH2F *thvsXd= (TH2F*)f->Get("pdc_xfp_vs_ypfp");
thvsXd->GetXaxis()->SetTitle("Ypfp(rad)");
thvsXd->GetYaxis()->SetTitle("Xfp(cm)");
TH2F *phivsXd= (TH2F*)f->Get("pdc_xfp_vs_xpfp");
phivsXd->GetXaxis()->SetTitle("Xpfp(rad)");
phivsXd->GetYaxis()->SetTitle("Xfp(cm)");
TH2F *thvsYd= (TH2F*)f->Get("pdc_yfp_vs_ypfp");
thvsYd->GetXaxis()->SetTitle("Ypfp(rad)");
thvsYd->GetYaxis()->SetTitle("Yfp(cm)");
TH2F *phivsYd= (TH2F*)f->Get("pdc_yfp_vs_xpfp");
phivsYd->GetXaxis()->SetTitle("Xpfp(rad)");
phivsYd->GetYaxis()->SetTitle("Yfp(cm)");
TH1F *deltaD = (TH1F*)f->Get("p_gtr_dp");
deltaD->GetXaxis()->SetTitle("dp/p(%)");


   Long64_t nentries = t1->GetEntries();

	for (int i = 0; i < nentries; i++) {
      		t1->GetEntry(i);
		
		if (normfac == 0){
		// Define kinematics
		double Ef = p_spec * (1.0 + 0.01*psdelta); //scattered electron energy //GeV
		double theta = TMath::ACos(cos_ts - psyptar * sin_ts) / TMath::Sqrt( 1. + psxptar * psxptar + psyptar * psyptar); // polar scattering angle relative to the beam line //rad
		double thetaDeg = theta / deg2rad;
		double Q2 = 4.0 * Ei * Ef * (TMath::Sin(theta / 2.0) * TMath::Sin(theta / 2.0)); //GeV^2
		double nu = Ei - Ef; //GeV
		double W2 = -Q2 + Mp * Mp + 2.0 * Mp * nu; // GeV^2
		double W=0.;
		if (W2 > 0) W = TMath::Sqrt(-Q2 + Mp * Mp + 2.0 * Mp * nu); //GeV
		if ( W < 1.075) cout << " W2 = " << W2 << endl;
		//weight with radiated CS
		double weight = (gr2D->Interpolate(Ef,thetaDeg))*(1e-6)*lumin*domega*(ep_max-ep_min)*1000.0/nentries;
		
		//Fill histograms with weight
		hW->Fill(W,weight);
		hQ2->Fill(Q2,weight);
		hom->Fill(nu,weight);
		hth->Fill(thetaDeg,weight);
		xfp->Fill(psxfp,weight);
		yfp->Fill(psyfp,weight);
		xpfp->Fill(psxpfp,weight);
		ypfp->Fill(psypfp,weight);
		delta->Fill(psdelta,weight);
		thvsX->Fill(psypfp,psxfp,weight);
		thvsY->Fill(psypfp,psyfp,weight);
		phivsX->Fill(psxpfp,psxfp, weight);
		phivsY->Fill(psxpfp, psyfp,weight);
		thvsphifp->Fill(psypfp,psxpfp,weight);
		XvsY->Fill(psxfp,psyfp,weight);
		//histograms below will be used in ratio plots
		hWr->Fill(W,weight);
		hQ2r->Fill(Q2,weight);
		homr->Fill(nu,weight);
		hthr->Fill(thetaDeg,weight);
		xfpr->Fill(psxfp,weight);
		yfpr->Fill(psyfp,weight);
		xpfpr->Fill(psxpfp,weight);
		ypfpr->Fill(psypfp,weight);
		deltar->Fill(psdelta,weight);
		thvsXr->Fill(psypfp,psxfp,weight);
		thvsYr->Fill(psypfp,psyfp,weight);
		phivsXr->Fill(psxpfp,psxfp, weight);
		phivsYr->Fill(psxpfp, psyfp,weight);
		thvsphifpr->Fill(psypfp,psxpfp,weight);
		XvsYr->Fill(psxfp,psyfp,weight);		
		}		
}

//2D FP plots
gROOT->Reset();
TCanvas *c1 = new TCanvas("c1", "c1", 800, 1200);
c1->Divide(4,3);
c1->cd(1);
gPad->SetLogz();
gPad->SetGrid();
gPad->SetGridx();
gPad->SetGridy();
thvsY->Draw("colz");
c1->cd(2);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
thvsYd->Draw("colz");
c1->cd(3);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
phivsY->Draw("colz");
c1->cd(4);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
phivsYd->Draw("colz");
c1->cd(5);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
thvsX->Draw("colz");
c1->cd(6);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
thvsXd->Draw("colz");
c1->cd(7);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
phivsX->Draw("colz");
c1->cd(8);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
phivsXd->Draw("colz");
c1->cd(9);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
thvsphifp->Draw("colz");
c1->cd(10);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
thvsphid->Draw("colz");
c1->cd(11);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
XvsY->Draw("colz");
c1->cd(12);
gPad->SetLogz();
gPad->SetGridx();
gPad->SetGridy();
XvsYd->Draw("colz");

// rates for MC and data
TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
c2->Divide(2,2);
c2->cd(1);
TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextSize(0.045);
leg->AddEntry(hW,"weighted MC" ,"l");
leg->AddEntry(hW_dat,"KPP 481" ,"l");
hW->SetFillColor(8);
hW->SetLineColor(8);
hW->Draw("hist");
hW_dat->SetMarkerStyle(kFullTriangleUp);
hW_dat->SetMarkerColor(kBlue);
hW_dat->Draw("P same");
leg->Draw();
c2->cd(2);
TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
leg1->SetFillColor(0);
leg1->SetBorderSize(0);
leg1->SetTextSize(0.045);
leg1->AddEntry(hQ2,"weighted MC" ,"l");
leg1->AddEntry(hQ2_dat,"KPP 481" ,"l");
hQ2->SetFillColor(8);
hQ2->SetLineColor(8);
hQ2->Draw("hist");
hQ2_dat->SetMarkerStyle(kFullTriangleUp);
hQ2_dat->SetMarkerColor(kBlue);
hQ2_dat->Draw("P same");
leg1->Draw();
c2->cd(3);
TLegend *leg2 = new TLegend(0.7,0.6,0.9,0.8);
leg2->SetFillColor(0);
leg2->SetBorderSize(0);
leg2->SetTextSize(0.045);
leg2->AddEntry(hth,"weighted MC" ,"l");
leg2->AddEntry(hth_dat,"KPP 481" ,"l");
hth->SetFillColor(8);
hth->SetLineColor(8);
hth->Draw("hist");
hth_dat->SetMarkerStyle(kFullTriangleUp);
hth_dat->SetMarkerColor(kBlue);
hth_dat->Draw("P same");
leg2->Draw();
c2->cd(4);
TLegend *leg3 = new TLegend(0.3,0.6,0.5,0.8);
leg3->SetFillColor(0);
leg3->SetBorderSize(0);
leg3->SetTextSize(0.045);
leg3->AddEntry(hom,"weighted MC" ,"l");
leg3->AddEntry(hom_dat,"KPP 481" ,"l");
hom->SetFillColor(8);
hom->SetLineColor(8);
hom->Draw("hist");
hom_dat->SetMarkerStyle(kFullTriangleUp);
hom_dat->SetMarkerColor(kBlue);
hom_dat->Draw("P same");
leg3->Draw();
//Mc delta and golden track delta
TCanvas *c3 = new TCanvas("c3", "c3", 800, 1200);
TLegend *leg4 = new TLegend(0.6,0.6,0.8,0.8);
leg4->SetFillColor(0);
leg4->SetBorderSize(0);
leg4->SetTextSize(0.045);
leg4->AddEntry(delta,"weighted MC" ,"l");
leg4->AddEntry(deltaD,"KPP 481" ,"l");
delta->SetFillColor(8);
delta->SetLineColor(8);
delta->Draw("hist");
deltaD->SetMarkerStyle(kFullTriangleUp);
deltaD->SetMarkerColor(kBlue);
deltaD->Draw("P same");
leg4->Draw();
//focal plane position/angle MC and data comparison
TCanvas *c4 = new TCanvas("c4", "c4", 800, 1200);
c4->Divide(2,2);
c4->cd(1);
TLegend *leg5 = new TLegend(0.6,0.6,0.8,0.8);
leg5->SetFillColor(0);
leg5->SetBorderSize(0);
leg5->SetTextSize(0.045);
leg5->AddEntry(xfp,"weighted MC" ,"l");
leg5->AddEntry(xfpD,"KPP 481" ,"l");
xfp->SetFillColor(8);
xfp->SetLineColor(8);
xfp->Draw("hist");
xfpD->SetMarkerStyle(kFullTriangleUp);
xfpD->SetMarkerColor(kBlue);
xfpD->Draw("P same");
leg5->Draw();
c4->cd(2);
TLegend *leg6 = new TLegend(0.6,0.6,0.8,0.8);
leg6->SetFillColor(0);
leg6->SetBorderSize(0);
leg6->SetTextSize(0.045);
leg6->AddEntry(yfp,"weighted MC" ,"l");
leg6->AddEntry(yfpD,"KPP 481" ,"l");
yfp->SetFillColor(8);
yfp->SetLineColor(8);
yfp->Draw("hist");
yfpD->SetMarkerStyle(kFullTriangleUp);
yfpD->SetMarkerColor(kBlue);
yfpD->Draw("P same");
leg6->Draw();
c4->cd(3);
TLegend *leg7 = new TLegend(0.6,0.6,0.8,0.8);
leg7->SetFillColor(0);
leg7->SetBorderSize(0);
leg7->SetTextSize(0.045);
leg7->AddEntry(xpfp,"weighted MC" ,"l");
leg7->AddEntry(xpfpD,"KPP 481" ,"l");
xpfp->SetFillColor(8);
xpfp->SetLineColor(8);
xpfp->Draw("hist");
xpfpD->SetMarkerStyle(kFullTriangleUp);
xpfpD->SetMarkerColor(kBlue);
xpfpD->Draw("P same");
leg7->Draw();
c4->cd(4);
TLegend *leg8 = new TLegend(0.6,0.6,0.8,0.8);
leg8->SetFillColor(0);
leg8->SetBorderSize(0);
leg8->SetTextSize(0.045);
leg8->AddEntry(ypfp,"weighted MC" ,"l");
leg8->AddEntry(ypfpD,"KPP 481" ,"l");
ypfp->SetFillColor(8);
ypfp->SetLineColor(8);
ypfp->Draw("hist");
ypfpD->SetMarkerStyle(kFullTriangleUp);
ypfpD->SetMarkerColor(kBlue);
ypfpD->Draw("P same");
leg8->Draw();
// MC/data Ratio plots for 2D fp variables
TCanvas *c5 = new TCanvas("c5", "c5", 800, 1200);
c5->Divide(3,2);
c5->cd(1);
gPad->SetLogz();
thvsYr->Divide(thvsYd);
thvsYr->SetTitle("Ratio of MC to KPP 481");
thvsYr->GetYaxis()->SetTitle("ratio");
thvsYr->Draw("colz");
c5->cd(2);
gPad->SetLogz();
phivsYr->Divide(phivsYd);
phivsYr->SetTitle("Ratio of MC to KPP 481");
phivsYr->GetYaxis()->SetTitle("ratio");
phivsYr->Draw("colz");
c5->cd(3);
gPad->SetLogz();
thvsXr->Divide(thvsXd);
thvsXr->SetTitle("Ratio of MC to KPP 481");
thvsXr->GetYaxis()->SetTitle("ratio");
thvsXr->Draw("colz");
c5->cd(4);
gPad->SetLogz();
phivsXr->Divide(phivsXd);
phivsXr->SetTitle("Ratio of MC to KPP 481");
phivsXr->GetYaxis()->SetTitle("ratio");
phivsXr->Draw("colz");
c5->cd(5);
gPad->SetLogz();
thvsphifpr->Divide(thvsphid);
thvsphifpr->SetTitle("Ratio of MC to KPP 481");
thvsphifpr->GetYaxis()->SetTitle("ratio");
thvsphifpr->Draw("colz");
c5->cd(6);
gPad->SetLogz();
XvsYr->Divide(XvsYd);
XvsYr->SetTitle("Ratio of MC to KPP 481");
XvsYr->GetYaxis()->SetTitle("ratio");
XvsYr->Draw("colz");
// MC/data Ratio plots for 1D fp variables
TCanvas *c6 = new TCanvas("c6", "c6", 800, 1200);
c6->Divide(2,2);
c6->cd(1);
xfpr->Divide(xfpD);
xfpr->SetFillColor(11);
xfpr->SetLineColor(11);
xfpr->GetYaxis()->SetRangeUser(0.0, 1.5);
xfpr->GetYaxis()->SetTitle("ratio");
xfpr->Draw("hist");
c6->cd(2);
yfpr->Divide(yfpD);
yfpr->SetFillColor(11);
yfpr->SetLineColor(11);
yfpr->GetYaxis()->SetRangeUser(0.0, 1.5);
yfpr->GetYaxis()->SetTitle("ratio");
yfpr->Draw("hist");
c6->cd(4);
ypfpr->Divide(ypfpD);
ypfpr->SetFillColor(11);
ypfpr->SetLineColor(11);
ypfpr->GetYaxis()->SetRangeUser(0.0, 1.5);
ypfpr->GetYaxis()->SetTitle("ratio");
ypfpr->Draw("hist");
c6->cd(3);
xpfpr->Divide(xpfpD);
xpfpr->SetFillColor(11);
xpfpr->SetLineColor(11);
xpfpr->GetYaxis()->SetRangeUser(0.0, 1.5);
xpfpr->GetYaxis()->SetTitle("ratio");
xpfpr->Draw("hist");
//MC/data ratio for target delta
TCanvas *c7 = new TCanvas("c7", "c7", 800, 1200);
deltar->Divide(deltaD);
deltar->SetFillColor(11);
deltar->SetLineColor(11);
deltar->GetYaxis()->SetRangeUser(0.0, 1.5);
deltar->GetYaxis()->SetTitle("ratio");
deltar->Draw("hist");
// MC/data Ratio plots for kin variables
TCanvas *c8 = new TCanvas("c8", "c8", 800, 1200);
c8->Divide(2,2);
c8->cd(2);
hQ2r->Divide(hQ2_dat);
hQ2r->SetFillColor(11);
hQ2r->SetLineColor(11);
hQ2r->GetYaxis()->SetRangeUser(0.0, 1.5);
hQ2r->GetYaxis()->SetTitle("ratio");
hQ2r->Draw("hist");
c8->cd(1);
hWr->Divide(hW_dat);
hWr->SetFillColor(11);
hWr->SetLineColor(11);
hWr->GetYaxis()->SetRangeUser(0.0, 1.5);
hWr->GetYaxis()->SetTitle("ratio");
hWr->Draw("hist");
c8->cd(3);
hthr->Divide(hth_dat);
hthr->SetFillColor(11);
hthr->SetLineColor(11);
hthr->GetYaxis()->SetRangeUser(0.0, 1.5);
hthr->GetYaxis()->SetTitle("ratio");
hthr->Draw("hist");
c8->cd(4);
homr->Divide(hom_dat);
homr->SetFillColor(11);
homr->SetLineColor(11);
homr->GetYaxis()->SetRangeUser(0.0, 1.5);
homr->GetYaxis()->SetTitle("ratio");
homr->Draw("hist");
}





