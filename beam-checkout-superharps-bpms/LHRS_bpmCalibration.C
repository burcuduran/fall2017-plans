{
//
//  Steering script to extract the bpm calibration out of a bulls eye scan
//
// by Bodo Reitz, April 2003
//
// needed: ascii table, containing run numbers and results from
//         corresponding harp scans
//         always adjust pedestals before analyzing bulls eye scans
//         never change pedestals without reanalyzing it
//         do not worry, if you dont have new pedestals
//         since the bpm calibration here corrects for them
//         anyhow

// not perfect yet, feel free to improve ;)
// currently the fits put an equal weight on each data point
// one could use a weighted fit instead

#include "TVector.h"

  //char harpfname[255];
char expnr[25];
char exppath[200];
char codafname[255];
char buf[255];
char *filestatus;
int numofscans=0;

TVector hax(1);
TVector hay(1);
TVector hbx(1);
TVector hby(1);

TVector bax(1);
TVector bay(1);
TVector bbx(1);
TVector bby(1);

TCanvas* c1 = new TCanvas("c1","BPM Pedestals");
c1->Divide(2,4);
  



int dum1,dum2;
double dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11,dum12,dum13,dum14;
FILE *fi;
TFile *filein;
  
double xp, xm, yp, ym;
// double pedestal[8] = {3578, 3415, 3351, 3363, 3581, 3502, 3526, 3383};//LHRS
 double pedestal[8] = {3674, 3499, 3468, 3478, 3766, 3666, 3570, 3405};//LHRS

double calib = 0.01887;

cout<<"Filename for Harp results: ";
char* harpfname = "harp_resultsL.text";
//cin>>harpfname;

fi=fopen(harpfname,"r");
  
filestatus=fgets( buf, 255, fi);
 while (filestatus !=NULL) {
   sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf",
	  &dum1,&dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&dum8,&dum9,&dum10);

   hax->ResizeTo(numofscans+1);
   hbx->ResizeTo(numofscans+1);
   hay->ResizeTo(numofscans+1);
   hby->ResizeTo(numofscans+1);
   bax->ResizeTo(numofscans+1);
   bbx->ResizeTo(numofscans+1);
   bay->ResizeTo(numofscans+1);
   bby->ResizeTo(numofscans+1);
   
   hax(numofscans)=dum3*1e-3;
   hay(numofscans)=dum5*1e-3;
   hbx(numofscans)=dum7*1e-3;
   hby(numofscans)=dum9*1e-3;
   

   cout << Form("/adaqfs/home/a-onl/rastersize/LeftHRS/save_old/run%d.root",dum1) << endl;   
   filein = TFile::Open(Form("/adaqfs/home/a-onl/rastersize/LeftHRS/save_old/run%dL.root",dum1));

 

   TCanvas *c1 = new TCanvas("c1","c1");
         c1->Divide(2,4);

	 c1->cd(1);
        
	Double_t xpa1 = bpmaraw1->GetBinCenter( bpmaraw1->GetMaximumBin() );
	bpmaraw1->Fit("gaus","","",xpa1-60, xpa1+60);
	xpa1 = bpmaraw1->GetFunction("gaus")->GetParameter(1);
        Double_t xpa = bpmaraw1->GetFunction("gaus")->GetParameter(2);
	bpmaraw1->Fit("gaus","","",xpa1-2.5*xpa, xpa1+2.5*xpa);
       	xpa1 = bpmaraw1->GetFunction("gaus")->GetParameter(1);
        xp =xpa1 - pedestal[0];     
       
       
        c1->cd(2);
        Double_t xma1 = bpmaraw2->GetBinCenter( bpmaraw2->GetMaximumBin() );
	bpmaraw2->Fit("gaus","","",xma1 -60, xma1+60);
	xma1 = bpmaraw2->GetFunction("gaus")->GetParameter(1);
        Double_t xma = bpmaraw2->GetFunction("gaus")->GetParameter(2);
        bpmaraw2->Fit("gaus","","",xma1-2.5*xma,xma1+2.5*xma);
        xma1 = bpmaraw2->GetFunction("gaus")->GetParameter(1);
         xm = xma1 - pedestal[1];
       
        c1->cd(3);       
        Double_t ypa1 = bpmaraw3->GetBinCenter( bpmaraw3->GetMaximumBin() );
	bpmaraw3->Fit("gaus","","",ypa1 -60, ypa1+60);
	ypa1 = bpmaraw3->GetFunction("gaus")->GetParameter(1);
        Double_t ypa = bpmaraw3->GetFunction("gaus")->GetParameter(2);
        bpmaraw3->Fit("gaus","","",ypa1-2.5*ypa,ypa1+2.5*ypa);
        ypa1 = bpmaraw3->GetFunction("gaus")->GetParameter(1);
        yp = ypa1 - pedestal[2];      
       
         c1->cd(4);
        Double_t yma1 = bpmaraw4->GetBinCenter( bpmaraw4->GetMaximumBin() );
	bpmaraw4->Fit("gaus","","",yma1 -60, yma1+60);
	yma1 = bpmaraw4->GetFunction("gaus")->GetParameter(1);
        Double_t yma = bpmaraw4->GetFunction("gaus")->GetParameter(2);
        bpmaraw4->Fit("gaus","","",yma1-2.5*yma,yma1+2.5*yma);
        yma1 = bpmaraw4->GetFunction("gaus")->GetParameter(1);
         ym = yma1 - pedestal[3];

   bax(numofscans)=calib*(xp-xm)/(xp+xm);
   bay(numofscans)=calib*(yp-ym)/(yp+ym);
   
  
       c1->cd(5);
       Double_t xpb1 = bpmbraw1->GetBinCenter( bpmbraw1->GetMaximumBin() );
	bpmbraw1->Fit("gaus","","",xpb1-60, xpb1+60);
	xpb1 = bpmbraw1->GetFunction("gaus")->GetParameter(1);
        Double_t xpb = bpmbraw1->GetFunction("gaus")->GetParameter(2);
	bpmbraw1->Fit("gaus","","",xpb1-2.5*xpb, xpb1+2.5*xpb);
       xpb1 = bpmbraw1->GetFunction("gaus")->GetParameter(1);
        xp =xpb1 - pedestal[4];     
        c1->cd(6);
        Double_t xmb1 = bpmbraw2->GetBinCenter( bpmbraw2->GetMaximumBin() );
	bpmbraw2->Fit("gaus","","",xmb1 -60, xmb1+60);
	xmb1 = bpmbraw2->GetFunction("gaus")->GetParameter(1);
        Double_t xmb = bpmbraw2->GetFunction("gaus")->GetParameter(2);
        bpmbraw2->Fit("gaus","","",xmb1-2.5*xmb,xmb1+2.5*xmb);
        xmb1 = bpmbraw2->GetFunction("gaus")->GetParameter(1);
         xm = xmb1 - pedestal[5];
        c1->cd(7);       
        Double_t ypb1 = bpmbraw3->GetBinCenter( bpmbraw3->GetMaximumBin() );
	bpmbraw3->Fit("gaus","","",ypb1 -60, ypb1+60);
	ypb1 = bpmbraw3->GetFunction("gaus")->GetParameter(1);
        Double_t ypb = bpmbraw3->GetFunction("gaus")->GetParameter(2);
        bpmbraw3->Fit("gaus","","",ypb1-2.5*ypb,ypb1+2.5*ypb);
        ypb1 = bpmbraw3->GetFunction("gaus")->GetParameter(1);
        yp = ypb1 - pedestal[6];      
        c1->cd(8);
        Double_t ymb1 = bpmbraw4->GetBinCenter( bpmbraw4->GetMaximumBin() );
	bpmbraw4->Fit("gaus","","",ymb1 -60, ymb1+60);
	ymb1 = bpmbraw4->GetFunction("gaus")->GetParameter(1);
        Double_t ymb = bpmbraw4->GetFunction("gaus")->GetParameter(2);
        bpmbraw4->Fit("gaus","","",ymb1-2.5*ymb,ymb1+2.5*ymb);
         ymb1 = bpmbraw4->GetFunction("gaus")->GetParameter(1);
         ym = ymb1 - pedestal[7];
   bbx(numofscans)=calib*(xp-xm)/(xp+xm);
   bby(numofscans)=calib*(yp-ym)/(yp+ym);
   



   numofscans++;
   filestatus=fgets( buf, 255, fi);
    }
 
 fclose(fi);

Double_t x11=0.;
Double_t x12=0.; 
Double_t x1yx=0.; 
Double_t x1yy=0.;
Double_t x21=0.; 
Double_t x22=0.; 
Double_t x2yx=0.; 
Double_t x2yy=0.;
Double_t x1=0.; 
Double_t x2=0.; 
Double_t yx=0.; 
Double_t yy=0.;

for(Int_t j=0;j<numofscans;j++) {
  x11=x11+bax(j)*bax(j);
  x12=x12+bax(j)*bay(j);
  x1yx=x1yx+bax(j)*hax(j);
  x1yy=x1yy+bax(j)*hay(j);
  x1=x1+bax(j);
  x2=x2+bay(j);
  yx=yx+hax(j);
  yy=yy+hay(j);
  x21=x21+bay(j)*bax(j);
  x22=x22+bay(j)*bay(j);
  x2yx=x2yx+bay(j)*hax(j);
  x2yy=x2yy+bay(j)*hay(j);
}

TMatrix trafo(3,3);
TVector inhomo(3);

trafo(0,0)=x11;
trafo(0,1)=x12;
trafo(0,2)=x1;
inhomo(0)=x1yx;

trafo(1,0)=x21;
trafo(1,1)=x22;
trafo(1,2)=x2;
inhomo(1)=x2yx;

trafo(2,0)=x1;
trafo(2,1)=x2;
trafo(2,2)=numofscans;
inhomo(2)=yx;

TMatrix itrafo(trafo);
itrafo.Invert();

TVector solu1(inhomo);
solu1*=itrafo;
  
inhomo(0)=x1yy;
inhomo(1)=x2yy;
inhomo(2)=yy;
TVector solu2(inhomo);
solu2*=itrafo;


x11=0.;
x12=0.; 
x1yx=0.; 
x1yy=0.;
x21=0.; 
x22=0.; 
x2yx=0.; 
x2yy=0.;
x1=0.; 
x2=0.; 
yx=0.; 
yy=0.;

for(Int_t j=0;j<numofscans;j++) {
  x11=x11+bbx(j)*bbx(j);
  x12=x12+bbx(j)*bby(j);
  x1yx=x1yx+bbx(j)*hbx(j);
  x1yy=x1yy+bbx(j)*hby(j);
  x1=x1+bbx(j);
  x2=x2+bby(j);
  yx=yx+hbx(j);
  yy=yy+hby(j);
  x21=x21+bby(j)*bbx(j);
  x22=x22+bby(j)*bby(j);
  x2yx=x2yx+bby(j)*hbx(j);
  x2yy=x2yy+bby(j)*hby(j);
}

trafo(0,0)=x11;
trafo(0,1)=x12;
trafo(0,2)=x1;
inhomo(0)=x1yx;

trafo(1,0)=x21;
trafo(1,1)=x22;
trafo(1,2)=x2;
inhomo(1)=x2yx;

trafo(2,0)=x1;
trafo(2,1)=x2;
trafo(2,2)=numofscans;
inhomo(2)=yx;

TMatrix itrafo(trafo);
itrafo.Invert();


TVector solu3(inhomo);
solu3*=itrafo;
  
inhomo(0)=x1yy;
inhomo(1)=x2yy;
inhomo(2)=yy;
TVector solu4(inhomo);
solu4*=itrafo;


cout<<"Please change the BPMA constants to:"<<endl;

cout<<solu1(0)<<" "<<solu1(1)<<" "<<solu2(0)<<" "<<solu2(1)<<" "<<solu1(2)<<" "<<solu2(2)<<endl;

cout<<"Please change the BPMB constants to:"<<endl;

cout<<solu3(0)<<" "<<solu3(1)<<" "<<solu4(0)<<" "<<solu4(1)<<" "<<solu3(2)<<" "<<solu4(2)<<endl;


//0.750796 -0.723848 0.77656 0.72107 -0.00143766 -6.81262e-06
  //0.687445 0.361161 0.00380886 0.943285 -0.00111306 -0.000418619


}

