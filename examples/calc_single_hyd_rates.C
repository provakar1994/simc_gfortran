#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;
void calc_single_hyd_rates(TString basename="",  Double_t normfac=0.166844E+08,Double_t cur=2., Double_t ngen=100000) {
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="worksim/"+basename+".root";
   TString outputhist;
   outputhist="worksim/"+basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf="worksim/"+basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("h10");
Float_t         psxfp; // position at focal plane ,+X is pointing down
 Float_t         psyfp; // X x Y = Z so +Y pointing central ray left
 Float_t         psxpfp; // dx/dz at focal plane
 Float_t         psypfp; //  dy/dz at focal plane
 Float_t         psztari; // thrown position along the beam direction
 Float_t         psytari;  //thrown  horizontal position X x Y = Z so +Y pointing central ray left at plane perpendicular to SHMS at z=0
 Float_t         psdeltai; // thrown  100*(p - pc)/pc with pc = central SHMS momentum 
 Float_t         psyptari; // thrown target dy/dz horizontal slope
 Float_t         psxptari; // thrown target dx/dz vertical slope
 Float_t         psytar; //reconstructed horizontal position
 Float_t         psdelta;//reconstructed
   Float_t         psyptar;//reconstructed
   Float_t         psxptar;//reconstructed
   //
  tsimc->SetBranchAddress("ssxfp",&psxfp);
   tsimc->SetBranchAddress("ssyfp",&psyfp);
   tsimc->SetBranchAddress("ssxpfp",&psxpfp);
   tsimc->SetBranchAddress("ssypfp",&psypfp);
   tsimc->SetBranchAddress("ssytari",&psytari);
   tsimc->SetBranchAddress("ssdeltai",&psdeltai);
   tsimc->SetBranchAddress("ssyptari",&psyptari);
   tsimc->SetBranchAddress("ssxptari",&psxptari);
   tsimc->SetBranchAddress("ssytar",&psytar);
   tsimc->SetBranchAddress("ssdelta",&psdelta);
   tsimc->SetBranchAddress("ssyptar",&psyptar);
   tsimc->SetBranchAddress("ssxptar",&psxptar);
   tsimc->SetBranchAddress("ssxptar",&psxptar);
  //
   Float_t Q2;
    Float_t W;
   Float_t Weight;
   tsimc->SetBranchAddress("Q2",&Q2);
   tsimc->SetBranchAddress("W",&W);
  tsimc->SetBranchAddress("Weight",&Weight);
   //
   //   Double_t normfac = 00.896997E+07;
   //      Double_t normfac ;
   //   Double_t ngen = 100000.;
   Double_t run_charge =1.;
   Double_t simc_charge= 1; // mC
   Double_t rat_charge;
   Double_t wfac;
  // Run 1929 0.136314E+08
   run_charge=cur/1000.*3600; // mC charge in one hour
   //
   rat_charge = run_charge/simc_charge;
  wfac=normfac/ngen*rat_charge;
  cout << " wfac = " << wfac << endl;
   //
   TH1F *hW = new TH1F("hW", "; W [Gev}; counts", 300, 0.9,1.5);
   TH1F *hQ2 = new TH1F("hQ2", "; Q2 {GeV2]", 200, 0, 5);
     TH1F *hpsdelta = new TH1F("hpsdelta", ";SHMS Delta [%]", 100, -10., 10.);
   TH1F *hpsyptar = new TH1F("hpsyptar", ";SHMS Yptar", 100, -0.05, .05);
   TH1F *hpsxptar = new TH1F("hpsxptar", ";SHMS Xptar", 100, -0.1, .1);
   TH1F *h_sxfp = new TH1F("h_sxfp", ";x_fp ", 100, -10., 10.);
   TH1F *h_syfp = new TH1F("h_syfp", ";y_fp ", 100, -10., 10.);
   TH1F *h_sxpfp = new TH1F("h_sxpfp", ";xp_fp ", 100, -.06, .06);
   TH1F *h_sypfp = new TH1F("h_sypfp", ";yp_fp ", 100, -.03, .03);
   //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		//		if (hsdelta>-8 && hsdelta<8) {
		hQ2->Fill(Q2,Weight*wfac);
		hW->Fill(W,Weight*wfac);
		if (W>0.8&&W<1.075) {
		  //		hW->Fill(W,Weight*wfac);
		hpsdelta->Fill(psdelta,Weight*wfac);
		hpsyptar->Fill(psyptar,Weight*wfac);
		hpsxptar->Fill(psxptar,Weight*wfac);
		h_sxfp->Fill(psxfp,Weight*wfac);
		h_sxpfp->Fill(psxpfp,Weight*wfac);
		h_syfp->Fill(psyfp,Weight*wfac);
		h_sypfp->Fill(psypfp,Weight*wfac);
		}
		//}
	}
	//
 TFile hsimc(outputhist,"recreate");
 HList.Add(hW);
 HList.Add(hQ2);
 HList.Add(hpsdelta);
 HList.Add(hpsyptar);
 HList.Add(hpsxptar);
 HList.Add(h_sxfp);
 HList.Add(h_syfp);
 HList.Add(h_sxpfp);
 HList.Add(h_sypfp);
 HList.Write();
 cout << " integral = " << hQ2->Integral() << " cur = " << cur << " uA" << endl;
 cout << " time for 100000 events = " << 100000/(hQ2->Integral()/3600)/60. << "min for  rate = " << hQ2->Integral()/3600 << "/sec  charge = " << run_charge << " mC in hour"<< endl;
 cout << " Plotted histograms put in root file = " << outputhist << endl;
  //
}
