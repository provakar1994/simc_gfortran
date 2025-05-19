#include <iostream>
#include <set>

#include "LookUpTableReader.h"

std::string replace_dot_with_p(double num);
double GetRate(std::string process, std::string infilestr, std::string anacuts, double normfac, double ngen, TCanvas *canv);
bool extractNgenAndNormfac(const std::string& filename, double& ngen, double& normfac);

void cp_trig_rate() {

  LookUpTableReader runplan;
  //runplan.readCSV("R-SIDIS_run_plan.csv");
  runplan.readCSV("new_R-SIDIS_run_plan_050625.csv");  

  //std::string outfile = Form("allrates.csv");
  std::string outfile = Form("new_allrates_050625.csv");  
  ofstream outfile_data; outfile_data.open(outfile);
  outfile_data << "index,sidisR,sidisRnc,totalR,totalnc\n";

  // Analysis cuts = "abs(ssdelta-5.)<15.&&abs(hsdelta)<8."
  std::string anacuts = "abs(ssdelta-5.)<15.&&abs(hsdelta)<8.";

  
  double ngen;
  double normfac;
  std::string process;
  std::string histfile;
  TCanvas *canv = new TCanvas("canv","canv",800,600);

  std::vector<int> new_kin_050625{5,6,7,8,9,10,11,12,13,14,17,18,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37};

  //for (int kin=1; kin<38; kin++) {
  for (int kin : new_kin_050625) {
    
    std::cout << "kin = " << kin << "\n";
    
    outfile_data << kin << ",";
    
    std::string Q2 = replace_dot_with_p(runplan.GetValueByKey(kin,2));
    std::string x = replace_dot_with_p(runplan.GetValueByKey(kin,1));
    std::string z = replace_dot_with_p(runplan.GetValueByKey(kin,3));
    std::string thpi = replace_dot_with_p(runplan.GetValueByKey(kin,9));
    double I = runplan.GetValueByKey(kin,11); //uA

    std::string infilestr = Form("kin%d_q2%s_x%s_z%s_thpi%s_pip_lh2",
				 kin,Q2.c_str(),x.c_str(),z.c_str(),thpi.c_str());    
    process = "sidis18";
    histfile = Form("../outfiles/%s_%s.hist",process.c_str(),infilestr.c_str());
    extractNgenAndNormfac(histfile, ngen, normfac);
    double sidis18_cpmc = GetRate(process,infilestr,anacuts,normfac,ngen,canv);
    double sidis18_cpmc_nc = GetRate(process,infilestr,"1",normfac,ngen,canv);  
    outfile_data << sidis18_cpmc*1e-3*I*3.6 << "," << sidis18_cpmc_nc*1e-3*I*3.6 << ",";

    // Delta
    process = "delta";
    histfile = Form("../outfiles/%s_%s.hist",process.c_str(),infilestr.c_str());
    extractNgenAndNormfac(histfile, ngen, normfac);
    //TCanvas *canv = new TCanvas("canv","canv",800,600);
    double delta_cpmc = GetRate(process,infilestr,anacuts,normfac,ngen,canv);
    double delta_cpmc_nc = GetRate(process,infilestr,"1",normfac,ngen,canv);  
    //outfile_data << delta_cpmc*1e-3*I*3.6 << "," << delta_cpmc_nc*1e-3*I*3.6 << ",";

    // Rho
    process = "rho";
    histfile = Form("../outfiles/%s_%s.hist",process.c_str(),infilestr.c_str());
    extractNgenAndNormfac(histfile, ngen, normfac);
    //TCanvas *canv = new TCanvas("canv","canv",800,600);
    double rho_cpmc = GetRate(process,infilestr,anacuts,normfac,ngen,canv);
    double rho_cpmc_nc = GetRate(process,infilestr,"1",normfac,ngen,canv);  
    //outfile_data << rho_cpmc*1e-3*I*3.6 << "," << rho_cpmc_nc*1e-3*I*3.6 << ",";

    // Exclusive
    process = "exclusive";
    histfile = Form("../outfiles/%s_%s.hist",process.c_str(),infilestr.c_str());
    extractNgenAndNormfac(histfile, ngen, normfac);
    // TCanvas *canv = new TCanvas("canv","canv",800,600);
    double exclusive_cpmc = GetRate(process,infilestr,anacuts,normfac,ngen,canv);
    double exclusive_cpmc_nc = GetRate(process,infilestr,"1",normfac,ngen,canv);  
    //outfile_data << exclusive_cpmc*1e-3*I*3.6 << "," << exclusive_cpmc_nc*1e-3*I*3.6 << ",";

    // Total
    double total_cpmc = sidis18_cpmc + delta_cpmc + rho_cpmc + exclusive_cpmc;
    double total_cpmc_nc = sidis18_cpmc_nc + delta_cpmc_nc + rho_cpmc_nc + exclusive_cpmc_nc;   
    outfile_data << total_cpmc*1e-3*I*3.6 << "," << total_cpmc_nc*1e-3*I*3.6 << "\n";
  }
}


std::string replace_dot_with_p(double num) {
  std::ostringstream oss;

  // Check if it's effectively an integer (to 1e-10 tolerance)
  if (std::fabs(num - static_cast<int>(num)) < 1e-10) {
    oss << std::fixed << std::setprecision(1);  // Force one decimal
  }

  oss << num;
  std::string str = oss.str();

  size_t dot_pos = str.find('.');
  if (dot_pos != std::string::npos) {
    str[dot_pos] = 'p';  // Replace '.' with 'p'
  }

  return str;
}

double GetRate(std::string process, std::string infilestr, std::string anacuts, double normfac, double ngen, TCanvas *canv)
{
  std::string infile = Form("%s_%s.root",process.c_str(),infilestr.c_str());
  //std::cout << infile << "\n"; 

  // reading ROOT files as df
  ROOT::EnableImplicitMT();  
  ROOT::RDataFrame simu_rdf_raw("h10",infile);
  
  auto columnNames = simu_rdf_raw.GetColumnNames();
  std::set<std::string> colSet(columnNames.begin(), columnNames.end());
  ROOT::RDF::RNode simu_rdf = simu_rdf_raw;
  if (colSet.count("phad") && colSet.count("nu")) {
    simu_rdf = simu_rdf
      .Define("z_custom", "sqrt(phad*phad+.14*.14)/nu")
      .Filter(anacuts.c_str());
  } else {
    simu_rdf = simu_rdf.Filter(anacuts.c_str());  
  }

  TH1F *hz;
  if (process.compare("rho")==0 || process.compare("sidis18")==0)
    hz = (TH1F*)simu_rdf.Histo1D({"hz","",200,0,1},"z","Weight")->Clone();
  else
    hz = (TH1F*)simu_rdf.Histo1D({"hz","",200,0,1},"z_custom","Weight")->Clone();

  // Scale properly
  hz->Scale(normfac/ngen);  

  // calculate rate (count/mC)
  double rate_cpmc = hz->Integral();
  
  canv->cd();    
  hz->DrawClone();

  return rate_cpmc;
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
