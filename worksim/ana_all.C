#include <iostream>
#include <regex>

#include "LookUpTableReader.h"

double const Mpip = 0.13957061; //GeV

int ExtractKinNumber(const std::string& input);
bool extractNgenAndNormfac(const std::string& filename, double& ngen, double& normfac);

void ana_all() {

  LookUpTableReader runplan;
  runplan.readCSV("R-SIDIS_run_plan.csv");

  std::string process = "sidis18";
  std::string meson = "pip";    //pip/pim
  std::string target = "lh2";   //lh2/ld2
  TChain *C = new TChain("h10");
  C->Add(Form("%s_kin*_%s_%s.root",process.c_str(),meson.c_str(),target.c_str()));

  // Create a new output file and clone the tree structure
  TFile *outFile = new TFile(Form("%s_allkin_%s_%s.root",
				  process.c_str(),meson.c_str(),target.c_str()), "RECREATE");
  TTree *Tout = C->CloneTree(0);  // Clone structure, not entries

  // Setting branch addresses
  float hsdelta;            C->SetBranchAddress("hsdelta", &hsdelta);
  float hsyptar;            C->SetBranchAddress("hsyptar", &hsyptar);    
  float hsxptar;            C->SetBranchAddress("hsxptar", &hsxptar);      
  float ssdelta;            C->SetBranchAddress("ssdelta", &ssdelta);
  float ssyptar;            C->SetBranchAddress("ssyptar", &ssyptar);    
  float ssxptar;            C->SetBranchAddress("ssxptar", &ssxptar);      

  float q_tree;             C->SetBranchAddress("q", &q_tree);      
  float Q2;                 C->SetBranchAddress("Q2", &Q2);    
  float nu;                 C->SetBranchAddress("nu", &nu);  
  float ppi;                C->SetBranchAddress("ppi", &ppi);
  float Weight;             C->SetBranchAddress("Weight", &Weight);    

  // Defining new branch
  bool T_deltacut;     Tout->Branch("deltacut", &T_deltacut, "deltacut/O");  
  int T_kin;           Tout->Branch("kin", &T_kin, "kin/I");
  double T_ebeam;      Tout->Branch("ebeam", &T_ebeam, "ebeam/D");  
  double T_ep;         Tout->Branch("ep", &T_ep, "ep/D");
  double T_thec;       Tout->Branch("thec", &T_thec, "thec/D"); // central value
  double T_the;        Tout->Branch("the", &T_the, "the/D");
  double T_phie;       Tout->Branch("phie", &T_phie, "phie/D");
  double T_thpic;      Tout->Branch("thpic", &T_thpic, "thpic/D"); // central value   
  double T_thpi;       Tout->Branch("thpi", &T_thpi, "thpi/D");    
  double T_phipi;      Tout->Branch("phipi", &T_phipi, "phipi/D");  
  double T_fweight;    Tout->Branch("fweight", &T_fweight, "fweight/D");  
  double T_y_lab; Tout->Branch("y_lab", &T_y_lab, "y_lab/D");
  double T_y_breit;    Tout->Branch("y_breit", &T_y_breit, "y_breit/D");    

  std::cout << std::endl;
  Long64_t Nevents = C->GetEntries(), nevent=0; 
  double timekeeper=0., timeremains=0.;
  int treenum=0, currenttreenum=0;
  int kin;  double normfac, ngen;
  double ebeam, ep, thec, thpic, the, thpi, phie, phipi;
  while(C->GetEntry(nevent++)) {
    // if(nevent % 100 == 0) std::cout << nevent << "/" << Nevents  << "\r";;
    // std::cout.flush();
    // ------

    // Access the current tree
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      //GlobalCut->UpdateFormulaLeaves();

      // finding the histfile and extracting info
      TString hftemp = C->GetFile()->GetName();
      hftemp = hftemp.ReplaceAll("worksim","outfiles").ReplaceAll(".root",".hist");
      kin = ExtractKinNumber(std::string(hftemp));
      extractNgenAndNormfac(std::string(hftemp), ngen, normfac);
      // extracting info from the runplan
      ebeam = runplan.GetValueByKey(kin,0); // GeV
      ep = runplan.GetValueByKey(kin,6);
      thec = runplan.GetValueByKey(kin,7)*TMath::DegToRad();      
      thpic = runplan.GetValueByKey(kin,9)*TMath::DegToRad();
    }

    // calculating physics angles 
    double phiec = 270.*TMath::DegToRad();
    double phipic = 90.*TMath::DegToRad(); // SHMS at beam left => phi=90 deg    

    the = acos((cos(thec)-hsyptar*sin(thec)*sin(phiec)/(1+hsxptar*hsxptar+hsyptar*hsyptar)));
    thpi = acos((cos(thpic)-ssyptar*sin(thpic)*sin(phipic)/(1+ssxptar*ssxptar+ssyptar*ssyptar)));

    phie = atan((sin(thec)*sin(phiec)+hsyptar*cos(thec))/(sin(thec)*cos(phiec)+hsxptar));
    phipi = atan((sin(thpic)*sin(phipic)+ssyptar*cos(thpic))/(sin(thpic)*cos(phipic)+ssxptar));    
    
    // calculating rapidity (y) in lab frame
    double Epi = sqrt(pow(ppi,2) + pow(Mpip,2));
    double pzpi = ppi*cos(thpi);
    double y_lab = 0.5 * log((Epi+pzpi)/(Epi-pzpi));

    // formulate the lab frame four vectors
    double Ee_prime = ebeam - nu;
    TLorentzVector k(0, 0, ebeam, ebeam);
    TLorentzVector k_prime(
			   Ee_prime*sin(the),
			   0,
			   Ee_prime*cos(the),
			   Ee_prime
			   );
    TLorentzVector pion_lab(
			    ppi*sin(thpi)*cos(phipi),
			    ppi*sin(thpi)*sin(phipi),
			    ppi*cos(thpi),
			    Epi
			    );    
    // Virtual photon q = k - k'
    TLorentzVector q = k - k_prime;

    // ---
    // Rapidity calculation in the Breit frame
    // Breit frame boost: make q = (0, 0, 0, Q)
    TVector3 beta_dir = q.Vect().Unit(); // along the q vector
    double beta_mag = q.E()/q.P();
    TVector3 boost_vec = beta_mag * beta_dir;
    //   
    // TLorentzVector q_breit = q;
    // q_breit.Boost(-boost_vec);
    // if (kin==1) {
    //   std::cout << "q_breit: (" << q_breit.Px() << ", " << q_breit.Py() << ", " << q_breit.Pz() << ", " << q_breit.E() << ")\n";
    //   std::cout << "q: (" << q.Px() << ", " << q.Py() << ", " << q.Pz() << ", " << q.E() << ")\n";
    // }
    //
    TLorentzVector pion_breit = pion_lab;
    pion_breit.Boost(-boost_vec);  // Boost to Breit frame         
 
    // Rapidity in Breit frame
    double Epi_breit = pion_breit.E();
    double pipz_breit = pion_breit.Pz();
    double y_breit = 0.5 * log((Epi_breit + pipz_breit) / (Epi_breit - pipz_breit));
    // ---
   
    // filling branches
    T_deltacut = abs(hsdelta)<8.&&abs(ssdelta-5.)<15.;
    T_kin = kin;
    T_ebeam = ebeam;
    T_ep = ep;
    T_thec = thec;
    T_the = the;    
    T_phie = phie;    
    T_thpic = thpic;
    T_thpi = thpi;
    T_phipi = phipi;        
    T_fweight = Weight*(normfac/ngen);
    T_y_lab = y_lab;
    T_y_breit = y_breit;
    
    Tout->Fill();
  }
  std::cout << "\n\n";  

  // Write and close
  outFile->cd();
  Tout->Write();
  outFile->Close();  
}

//////////////////////////////////////////////////////
int ExtractKinNumber(const std::string& input) {
    std::smatch match;
    std::regex pattern("kin(\\d+)");

    if (std::regex_search(input, match, pattern) && match.size() > 1) {
        return std::stoi(match[1]);
    } else {
        std::cerr << "No 'kin' number found in string.\n";
        return -1;
    }
}

bool extractNgenAndNormfac(const std::string& filename, double& ngen, double& normfac)
{
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Error opening file: " << filename << "\n";
    return false;
  }

  std::string line;
  ngen = -1.0;
  normfac = -1.0;

  while (std::getline(infile, line)) {
    std::istringstream iss(line);

    if (line.find("Ngen (request)") != std::string::npos) {
      std::string token;
      while (iss >> token) {
	if (token == "=") {
	  iss >> ngen;
	  break;
	}
      }
    }

    if (line.find("normfac") != std::string::npos) {
      std::string token;
      while (iss >> token) {
	if (token == "=") {
	  iss >> normfac;
	  break;
	}
      }
    }

    if (ngen != -1.0 && normfac != -1.0) break;  // Both values found
  }

  infile.close();

  return (ngen != -1.0 && normfac > 0);
}

  // C->SetBranchStatus("*", 0);
  // double hsdelta;          C->SetBranchAddress("hsdelta", &hsdelta);
  // double hsyptar;          C->SetBranchAddress("hsyptar", &hsyptar);    
  // double hsxptar;          C->SetBranchAddress("hsxptar", &hsxptar);      
  // double hsytar;           C->SetBranchAddress("hsytar", &hsytar);    
  // double hsxfp;            C->SetBranchAddress("hsxfp", &hsxfp);      
  // double hsxpfp;           C->SetBranchAddress("hsxpfp", &hsxpfp);
  // double hsyfp;            C->SetBranchAddress("hsyfp", &hsyfp);      
  // double hsypfp;           C->SetBranchAddress("hsypfp", &hsypfp);        
  // double ssdelta;          C->SetBranchAddress("ssdelta", &ssdelta);
  // double ssyptar;          C->SetBranchAddress("ssyptar", &ssyptar);    
  // double ssxptar;          C->SetBranchAddress("ssxptar", &ssxptar);      
  // double ssytar;           C->SetBranchAddress("ssytar", &ssytar);    
  // double ssxfp;            C->SetBranchAddress("ssxfp", &ssxfp);      
  // double ssxpfp;           C->SetBranchAddress("ssxpfp", &ssxpfp);
  // double ssyfp;            C->SetBranchAddress("ssyfp", &ssyfp);      
  // double ssypfp;           C->SetBranchAddress("ssypfp", &ssypfp);        
  // double q;                C->SetBranchAddress("q", &q);
  // double nu;               C->SetBranchAddress("nu", &nu);
  // double Q2;               C->SetBranchAddress("Q2", &Q2);
  // double W;                C->SetBranchAddress("W", &W);
  // double epsilon;          C->SetBranchAddress("epsilon", &epsilon);
  // double Em;               C->SetBranchAddress("Em", &Pm);              
