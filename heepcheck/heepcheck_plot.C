#include <iostream>
#include <math.h>
#include <cmath>
#include<TStyle.h>
using namespace std;



Float_t loss_proton(Float_t t_z, Float_t t_a, Float_t betap, Float_t thick);
Float_t loss_electron(Float_t t_z, Float_t t_a, Float_t t_dens, Float_t t_thick);


int heepcheck (TString target,TString sign,TString hpart)  // arguments are : target cell type (target1/target2/target3); sign of eloss (add/subs/no); type of particle in HMS (e/p)


{ 
  gROOT->Reset();  
  Int_t run;
  Double_t  En, theta_ee, p_Ecen, theta_pp, p_Pcen;
  Double_t E_off, theta_Eoff, p_eoff, theta_poff, p_poff;
  Float_t tar_a = 1.0;  // parameters for hydrogen target
  Float_t tar_z= 1.0;
  Float_t radlength = 63.05;
  Float_t tar_dens = 0.0708;
  Float_t M = 0.000511; //electron mass
  Float_t Mp = 0.93827;  //proton mass

  Float_t hthicks = 0.016*2.54;   // Scattering chamber window 16 mil Al
  Float_t hthicka = 15.0;         // 15 cm of air
  Float_t hthickw1 = 0.017*2.54;   // HMS entrance window 17 mil Kevlar
  Float_t hthickw2 = 0.005*2.54;  //  HMS entrance window 5 mil mylar
  Float_t hthick3 = 0.020*2.54;   // thickness of HMS exit window, 20 mil Ti 
  Float_t hthick4 = 0.002*2.54;   // thickness of mylar windows on HMS DC, 1 + 1 mil
  Float_t hthick5 = 0.080*2.54; //40 + 40 mil entrance+ exit window of HMS Cerenkov tank
  
  
  Float_t sthicks = 0.008*2.54;   // Scattering chamber window 16 mil Al
  Float_t sthicka = 15.0;         // 15 cm of air
  Float_t sthickw1 = 0.020*2.54;   // SHMS entrance window 20 mil Al
  Float_t sthick3 = 0.020*2.54;   // thickness of SHMS exit window, 20 mil Al
  Float_t sthick4 = 0.002*2.54;   // thickess of mylar window on SHMS DC, 1+ 1 mil
  Float_t sthick5 = 0.004*2.54;  // 2+2 mil Al on SHMS Cerenkov tank 

  //***********read in kinematics and offsets *****************
  //*****Beam energy, e-angle, e-momentum, p-angle, p-momentum, beam_offset, e-angle offset, e-momentum affset, p-angle offset, p-momentum offset
    ifstream myfile;
    myfile.open("heepcheck.input");
    while (myfile.eof() == false) {
     myfile >> run >> En >> theta_ee >> p_Ecen >> theta_pp >> p_Pcen >> E_off >> theta_Eoff >> p_eoff >> theta_poff >> p_poff;
    }
    myfile.close();

    // all central values and offsets
  const Float_t theta_e = theta_ee + theta_Eoff;
  const Float_t ep_cen = p_Ecen + p_eoff;
  const Float_t theta_p = theta_pp + theta_poff;
  const Float_t pp_cen = p_Pcen + p_poff;
  const Float_t E = En + E_off;  // beam energy with offset
  //cout << En << " " << theta_ee << " " << p_Pcen << " " << endl;

  // define all histograms
  TH1D *h33 = new TH1D("h33","W", 150, 0.8, 1.1);
 
  TH1D *h44= new TH1D("h44","Em", 100, -0.05, .05);

  TH1D *h55= new TH1D("h55","Pm_par", 100, -0.1, 0.1);

  TH1D *h66= new TH1D("h66","Pm_per",100, -0.05, 0.1);  
  
  TCanvas* c1 = new TCanvas("c1", "My histo8", 1000,1000);
  c1->Divide(2,2);


  
  TChain * chain = new TChain(Form("h%d",9500));
  chain->Add(Form("coin%d.root", run));
  // cout <<run<<endl;

  Float_t deltae, deltap, eyptar, exptar, pyptar, pxptar, esphi, psphi, ecer_npe, pcer_npe;

  if (hpart == "e") {  // assign root tree variables if electron in HMS
   chain->SetBranchAddress("hsdelta",&deltae);
   chain->SetBranchAddress("ssdelta",&deltap);
   chain->SetBranchAddress("hsyptar",&eyptar);
   chain->SetBranchAddress("hsxptar",&exptar);
   chain->SetBranchAddress("ssyptar",&pyptar);
   chain->SetBranchAddress("ssxptar",&pxptar);
   chain->SetBranchAddress("hcer_npe",&ecer_npe);
   chain->SetBranchAddress("scer_npe",&pcer_npe);
  }
  else if (hpart == "p") {  // root tree variables if electron in SHMS
   chain->SetBranchAddress("hsdelta",&deltap);
   chain->SetBranchAddress("ssdelta",&deltae);
   chain->SetBranchAddress("hsyptar",&pyptar);
   chain->SetBranchAddress("hsxptar",&pxptar);
   chain->SetBranchAddress("ssyptar",&eyptar);
   chain->SetBranchAddress("ssxptar",&exptar);
   chain->SetBranchAddress("hcer_npe",&pcer_npe);
   chain->SetBranchAddress("scer_npe",&ecer_npe);
  }
  //   chain->SetBranchAddress("Em", &Emm);
  // chain->SetBranchAddress("sinvmass", &W);
  // chain->SetBranchAddress("Pmpar",&Pmpar);
  // chain->SetBranchAddress("Pmper",&Pmper);

  h33->Reset();
  h44->Reset();
  h55->Reset();
  h66->Reset();


  if (hpart == "e") { // electron in HMS
 
   if (target == "target1"){   // 10 cm. 2.65" diameter target cell
    Float_t tar_len = 10.0; 
    Float_t tar_thick = tar_len*tar_dens;
    Float_t hthkfr = abs(tar_thick/2./cos(theta_ee));//(0.00106342)
    Float_t sthkfr = abs(tar_thick/2./cos(theta_pp));//(0.0004965)
    Float_t hthkside = abs(1.325*2.54*tar_dens/sin(theta_ee));//(0.248047)
    Float_t sthkside = abs(1.325*2.54*tar_dens/sin(theta_pp));//(0.52993)
   }

   else if (target  == "target2"){  // 10 cm, 1.6" diameter cell 
    Float_t tar_len = 10.0; 
    Float_t tar_thick = tar_len*tar_dens;
    Float_t  hthkfr = abs(tar_thick/2./cos(theta_ee));//(0.00106342)
    Float_t  sthkfr = abs(tar_thick/2./cos(theta_pp));//(0.0004965)
    Float_t  hthkside = abs(0.8*2.54*tar_dens/sin(theta_ee));//(0.248047)
    Float_t  sthkside = abs(0.8*2.54*tar_dens/sin(theta_pp));//(0.52993)
   }

   else if (target == "target3"){  // 4 cm, 1.6" diameter cell
    Float_t tar_len = 4.0;
    Float_t tar_thick = tar_len*tar_dens;
    Float_t  hthkfr = abs(tar_thick/2./cos(theta_ee));//(0.00106342) 
    Float_t  sthkfr = abs(tar_thick/2./cos(theta_pp));//(0.0004965)
    Float_t  hthkside = abs(0.8*2.54*tar_dens/sin(theta_ee));//(0.248047)
    Float_t  sthkside = abs(0.8*2.54*tar_dens/sin(theta_pp));//(0.52993)
   }
   
   if(hthkfr>hthkside){    // does the particle in HMS go through the endcap or the sidewall
  
    Float_t  hthick1 = hthkside;
    Float_t  hthick2 = 0.005*2.54*2.70/sin(theta_ee);
	 
   }

   else{
   
     Float_t  hthick1 = hthkfr;
     Float_t  hthick2 = 0.005*2.54*2.70/cos(theta_ee);//(0.0811083)
   }

   if(sthkfr>sthkside){  // does the particle in SHMS go through the endcap or the sidewall
     Float_t sthick1 = sthkside;
     Float_t sthick2 = 0.005*2.54*2.70/sin(theta_pp);
   }

   else{
     Float_t  sthick1 = sthkfr;
     Float_t  sthick2 = 0.005*2.54*2.70/cos(theta_pp);
   }

  
   Float_t loss_01 = loss_electron(1.0, 1.0, 0.0708, tar_len/2.*0.0708);//energy loss of beam in hyrdogen target
   Float_t loss_02 = loss_electron(13.0, 27., 2.7, 0.005*2.54*2.7);//energy loss of beam in Al endcap

   Float_t loss_11 = loss_electron(1.0, 1.0, 0.0708, hthick1) + loss_electron(13.0, 27.0, 2.7, hthick2); //energy loss of scattered electron in target
   Float_t loss_12 = loss_electron(13., 27., 2.7, hthicks); // energy loss of scattered electon in scattering chamber exit window;, Al
   Float_t loss_13 = loss_electron(7.32, 14.68, 0.00121, hthicka); // energy loss of scattered electron in air between scatt chamber and HMS entrance
   Float_t loss_14 = loss_electron(2.67, 4.67, 0.74, hthickw1) + loss_electron(0.52, 1.0, 1.39, hthickw2); //loss in HMS entrance window kevlar+mylar
   Float_t loss_15 = loss_electron(22., 47.9, 4.54, hthick3); // energy loss of scattered electron in HMS spectrometer exit Titanium
   Float_t loss_16 = loss_electron(0.52, 1.0, 1.39, hthick4); // energy loss of scattered electron in HMS DC mylar
   Float_t loss_17 = loss_electron(13., 27., 2.7, hthick5); // energy loss of scattered electron in HMS Cer Al


   Float_t sbeta = p_Pcen/sqrt(p_Pcen*p_Pcen + Mp*Mp);  //proton beta for central momentum 

   Float_t loss_21 = loss_proton(1.0, 1.0, 0.0708, sbeta, sthick1) + loss_proton(13.0, 26.98, 2.7, sbeta, sthick2); // energy loss of scattered proton in target
   Float_t loss_22 = loss_proton(13.0, 27.0, 2.7, sbeta, sthicks); //energy loss of scattered proton in scattering chamber exit window
   Float_t loss_23 = loss_proton(7.32, 14.68, 0.00121, sbeta, sthicka); //energy loss of scattered proton in air between scatt. chamber and SHMS entrance
   Float_t loss_24 = loss_proton(13.0, 27.0, 2.7, sbeta, sthickw1); //energy loss of scattered proton in entance window of SHMS
   Float_t loss_25 = loss_proton(13.0, 27.0, 2.7, sbeta, sthick3); // energy loss of scattered proton in SHMS exit
   Float_t loss_26 = loss_proton(0.52, 1.0, 1.39, sbeta, sthick4);  // energy loss of scattered proton in SHMS DC
   Float_t loss_27 = loss_proton(13., 27.0, 2.7, sbeta, sthick5);  // energy loss of scattered proton in SHMS Cer

  }
 
  else {   //  p in HMS
 
   if (target == "target1"){   // 10 cm. 2.65" diameter target cell
    Float_t tar_len = 10.0; 
    Float_t tar_thick = tar_len*tar_dens;
    Float_t sthkfr = abs(tar_thick/2./cos(theta_ee));//(0.00106342)
    Float_t hthkfr = abs(tar_thick/2./cos(theta_pp));//(0.0004965)
    Float_t sthkside = abs(1.325*2.54*tar_dens/sin(theta_ee));//(0.248047)
    Float_t hthkside = abs(1.325*2.54*tar_dens/sin(theta_pp));//(0.52993)
   }

   else if (target  == "target2"){  // 10 cm, 1.6" diameter cell 
    Float_t tar_len = 10.0; 
    Float_t tar_thick = tar_len*tar_dens;
    Float_t  sthkfr = abs(tar_thick/2./cos(theta_ee));//(0.00106342)
    Float_t  hthkfr = abs(tar_thick/2./cos(theta_pp));//(0.0004965)
    Float_t  sthkside = abs(0.8*2.54*tar_dens/sin(theta_ee));//(0.248047)
    Float_t  hthkside = abs(0.8*2.54*tar_dens/sin(theta_pp));//(0.52993)
   }

   else if (target == "target3"){  // 4 cm, 1.6" diameter cell
    Float_t tar_len = 4.0;
    Float_t tar_thick = tar_len*tar_dens;
    Float_t  sthkfr = abs(tar_thick/2./cos(theta_ee));//(0.00106342) 
    Float_t  hthkfr = abs(tar_thick/2./cos(theta_pp));//(0.0004965)
    Float_t  sthkside = abs(0.8*2.54*tar_dens/sin(theta_ee));//(0.248047)
    Float_t  hthkside = abs(0.8*2.54*tar_dens/sin(theta_pp));//(0.52993)
   }
   
   if(hthkfr>hthkside){    // does the particle in HMS go through the endcap or the sidewall
  
    Float_t  hthick1 = hthkside;
    Float_t  hthick2 = 0.005*2.54*2.70/sin(theta_pp);
	 
   }

   else{
   
     Float_t  hthick1 = hthkfr;
     Float_t  hthick2 = 0.005*2.54*2.70/cos(theta_pp);//(0.0811083)
   }

   if(sthkfr>sthkside){  // does the particle in SHMS go through the endcap or the sidewall
     Float_t sthick1 = sthkside;
     Float_t sthick2 = 0.005*2.54*2.70/sin(theta_ee);
   }

   else{
     Float_t  sthick1 = sthkfr;
     Float_t  sthick2 = 0.005*2.54*2.70/cos(theta_ee);
   }

  
   Float_t loss_01 = loss_electron(1.0, 1.0, 0.0708, tar_len/2.*0.0708);//energy loss of beam in hyrdogen target
   Float_t loss_02 = loss_electron(13.0, 27, 2.7, 0.005*2.54*2.7);//energy loss of beam in Al endcap

   Float_t loss_11 = loss_electron(1.0, 1.0, 0.0708, sthick1) + loss_electron(13.0, 27.0, 2.7, sthick2); //energy loss of scattered electron in target
   Float_t loss_12 = loss_electron(13.0, 27.0, 2.7, sthicks); // energy loss in scatt. chamber exit foil
   Float_t loss_13 = loss_electron(7.32, 14.68, 0.00121, sthicka); //energy loss in air between scatt. chamber and SHMS entrance
   Float_t loss_14 = loss_electron(13.0, 27.0, 2.7, sthickw1); //energy loss in entrance window of SHMS 
   Float_t loss_15 = loss_electron(13., 27., 2.7, sthick3); // energy loss of scattered electron in SHMS spectrometer exit
   Float_t loss_16 = loss_electron(0.52, 1.0, 1.39, sthick4); // energy loss of scattered electron in SHMS DC
   Float_t loss_17 = loss_electron(13., 27., 2.7, sthick5); // energy loss of scattered electron in SHMS Cer

   Float_t hbeta = p_Pcen/sqrt(p_Pcen*p_Pcen + Mp*Mp);  //proton beta for central momentum 

   Float_t loss_21 = loss_proton(1.0, 1.0, 0.0708, hbeta, hthick1) + loss_proton(13.0, 26.98, 2.7, hbeta, hthick2); // energy loss of scattered proton in target
   Float_t loss_22 = loss_proton(13.0, 27.0, 2.7, hbeta, hthicks); // loss due to scatt. chamber exit window
   Float_t loss_23 = loss_proton(7.32, 14.68, 0.00121, hbeta, hthicka); // loss due to air between scatt. chamber and HMS entrance
   Float_t loss_24 = loss_proton(2.67, 4.67, 0.74, hbeta, hthickw1) + loss_proton(0.52, 1.0, 1.39, hbeta, hthickw2); //loss in the HMS entrace window 
   Float_t loss_25 = loss_proton(22.0, 47.9, 4.54, hbeta, hthick3); // energy loss of scattered proton in HMS exit
   Float_t loss_26 = loss_proton(0.52, 1.0, 1.39, hbeta, hthick4);  // energy loss of scattered proton in HMS DC
   Float_t loss_27 = loss_proton(13., 27., 2.7, hbeta, hthick5);  // energy loss of scattered proton in HMS Cer
  }

  //*********************** Calculate the total energy loss

   Double_t totLoss_beam = loss_01+loss_02;  // total eloss of incoming beam
   Double_t totLoss_e = loss_11+loss_12+loss_13+loss_14+loss_15+loss_16+loss_17;  // total eloss of scattered electron
   Double_t totLoss_p = loss_21+loss_22+loss_23+loss_24+loss_25+loss_26+loss_27;  // total eloss of scattered proton

   //*********************************************

  
   Long64_t nentries = chain->GetEntries();  // read in number of entries in the root tree             
      
   for ( Long64_t i=0; i < nentries; i++)   // loop over all events
    
      {

	Double_t q_x, q_y, q_z, Ep, Epp, Em, Q_2, Q_22, W_2, p_beam, Pm_par, Pm_per, Pm_t, E_corr;

	Double_t  esp_z, psp_z, p_ex, p_ey, p_ez, p_px, p_py, p_pz, estheta, pstheta;

	Double_t p_e, p_p; //variable of momentum range for e and p
	
	chain->GetEntry(i);

	 
	p_e = (1. + deltae/100) * ep_cen ;//calculating the momentum
      
	p_p = (1. + deltap/100) * pp_cen;

	// Taking care of energy loss based on what was requested
	 if (sign == "subs"){
	   // cout << "It seems like you want to subs the energy loss" <<endl;
           E_corr = E - totLoss_beam; // beam energy after energy loss	
           Ep = sqrt(p_e * p_e + M*M)-totLoss_e;//scattered electron energy
           Epp = sqrt(p_p*p_p + Mp*Mp)-totLoss_p;// energy of scattered proton
	 }
	 else if(sign == "add"){
	   // cout << "It seems you want to add" <<endl;
           E_corr = E - totLoss_beam; 
	   Ep = sqrt(p_e * p_e + M*M)+totLoss_e;
           Epp = sqrt(p_p*p_p + Mp*Mp)+totLoss_p;
	 }
	 else if(sign == "no"){
           E_corr = E ;
	   Ep = sqrt(p_e * p_e + M*M);
	   Epp = sqrt(p_p*p_p + Mp*Mp);
	 }


	esp_z = p_e/sqrt(1.0 + exptar*exptar + eyptar*eyptar);

	psp_z = p_p/sqrt(1.0 + pxptar*pxptar + pyptar*pyptar);

	//***************Component of momentum of e and p******************
	 p_ex =  esp_z * exptar; 
	 p_ey =  esp_z * (eyptar * cos(theta_e) - sin(theta_e));
	 p_px =  psp_z * pxptar; 
	 p_py =  psp_z * (pyptar * cos(theta_p) - sin(theta_p));
      
        if (hpart == "e") {
	 p_ez =  esp_z * (eyptar * sin(theta_e) + cos(theta_e));
	 p_pz =  psp_z * (-pyptar * sin(theta_p) + cos(theta_p));
        }  
        else if (hpart == "p") {
	 p_ez =  esp_z * (-eyptar * sin(theta_e) + cos(theta_e));
	 p_pz =  psp_z * (pyptar * sin(theta_p) + cos(theta_p));
        }
       //*******************component of q******************
	
	p_beam = sqrt(E_corr*E_corr-M*M); //momentum of initial electron
	q_x = -p_ex;
	q_y = -p_ey;
	q_z =  p_beam - p_ez;
  
	//************** recalculate scattering and azimuthal angles

	if (abs(p_ez/p_e)<=1.0)
	  {
	    estheta = acos(p_ez/p_e); 
	  }
      
	else
	  {
	    estheta = -10.0;
	  }


	if (abs(p_pz/p_p)<=1.0)
	  {
	    pstheta = acos(p_pz/p_p);
	  }
      
	else
	  {
	    pstheta = -10.0;
	  }      
	  
	esphi = atan(p_ez/p_ey);
	
	psphi = atan(p_pz/p_py);
	   

	 // calculate the physics variables Q^2, W^2, Em, Pmpar, Pmper
 
	Q_2 = -((E_corr-Ep) * (E_corr-Ep) - (q_x*q_x + q_y*q_y + q_z*q_z));

	Q_22 = 4* E_corr* Ep * sin(estheta/2) * sin(estheta/2);

	W_2 = sqrt(Mp * Mp + 2 * Mp * (E_corr - Ep) - Q_22);    
       
	Em = E_corr + Mp - (Ep + Epp); // missing Energy.

	Pm_par = q_z - p_pz; //  missing parallel momentum

	Pm_per = q_x - p_px ;//sqrt(q_x*q_x + q_y*q_y) - sqrt(p_px*p_px + p_py*p_py);// missing perpendicular momentum.

	Pm_t = sqrt(Pm_par*Pm_par + Pm_per*Pm_per);//total missing momentum
      
	
	if (abs(W_2 -0.938)<0.15){                  // apply cut around W = proton mass             // && hcer_npe>0 && scer_npe ==0){
	  
    	  h33->Fill(W_2);
	   
    	  h44->Fill(Em);

	  h55->Fill(Pm_par);

	  h66->Fill(Pm_per);
	}
	
      }	
  
 cout <<E_off<<endl;
  c1->cd(1);
  h33->Draw();
  h33->SetLineColor(kRed);
	
  c1->cd(2);
  h44->Draw();
  h44->SetLineColor(kRed);
          
  c1->cd(3);
  h55->Draw();
  h55->SetLineColor(kRed);
      
  c1->cd(4);
  h66->Draw();
  h66->SetLineColor(kRed);
   
  c1->SaveAs("Plots_heepcheck.pdf");
    return 0;   

}

Float_t loss_electron(Float_t t_z, Float_t t_a, Float_t t_dens, Float_t t_thick){  // loss for electrons
 
  Float_t loss_e = (0.0001536*t_z)/(t_a*t_thick*(19.26 +log(t_thick/t_dens)));
  //cout << loss<<endl;
  return loss_e;
}

Float_t loss_proton (Float_t t_z, Float_t t_a, Float_t t_dens, Float_t betap, Float_t thick){         // loss for proton.

  Float_t me_gev=0.000510999;
  Float_t pmass=0.93827;
  Float_t gammap=1./sqrt(1.-betap*betap);
  Float_t tmax=abs(2.*me_gev*betap*betap*gammap*gammap/(1+2.*abs(gammap)*me_gev/pmass+(me_gev/pmass)*(me_gev/pmass))); //max energy transfer to orbitl electrons
  Float_t hnup=28.816e-09*sqrt(abs(t_dens*t_z/t_a));  //plasma frequency (h*nu_p)
  Float_t logbg=log(betap*gammap)/log(10.);
  Float_t icon_ev = 1.0 ;
  if (t_z < 1.5) {
    icon_ev = 21.8;
  } else if (t_z < 13) {
    icon_ev = 12.*t_z+7;
  } else {
    icon_ev = t_z*(9.76+58.8*(pow(t_z, -1.19)));
  }
  Float_t icon_gev = icon_ev*1.0e-09;
  Float_t C0 = -2*(log(icon_gev)-log(hnup)+.5);
  Float_t denscorr = 0.; 

  if (logbg < 0.) { 
    denscorr = 0.;
  } else if (logbg < 3.) {
    denscorr = C0+2.*log(10.)*logbg+abs(C0/27.)*(3.-logbg)*(3.-logbg)*(3.-logbg) ;
  } else if (logbg < 4.7) {
    denscorr = C0+2*log(10.)*logbg ;
  } else {
    denscorr=C0+2.*log(10.)*4.7 ;
  }  
  
  Float_t loss_p = .5*(log(2.*me_gev)+2*log(betap)+2*log(gammap)+log(tmax)-2*log(icon_gev)) - betap*betap -denscorr/2. ; 
  Float_t loss_p = 2.*0.1536e-03*t_z/t_a*thick/betap/betap * loss_p;
  return loss_p;
}

 
