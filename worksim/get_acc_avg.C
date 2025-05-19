#include "TLatex.h"
#include "LookUpTableReader.h"

void get_acc_avg() {
  ROOT::RDataFrame df("h10", "sidis18_allkin_pip_lh2.root");  // tree name and file

  std::vector<ROOT::RDF::RResultPtr<TH1D>> hzs,hxbjs,hQ2s,hys;

  std::string outfile = Form("output_get_acc_avg.csv");
  ofstream outfile_data; outfile_data.open(outfile);
  outfile_data << "index,avgx,rmsx,avgQ2,rmsQ2,avgz,rmsz,avgy_breit,rmsy_breit\n";  

  for (int kin=1; kin<38; kin++) {

    std::cout << "kin = " << kin << "\n";
    outfile_data << kin << ",";
    
    auto cutName = Form("deltacut&&kin==%d",kin);
    auto df_cut = df.Filter(cutName);
    
    // x
    auto hxbjtemp = df_cut.Histo1D({"hxbjtemp", Form("hxbj_kin%d", kin), 100, 0, 0.8}, "xbj", "fweight");
    double avgx = hxbjtemp->GetMean();
    double rmsx = hxbjtemp->GetRMS();    
    hxbjs.push_back(hxbjtemp);
    // Q2
    auto hQ2temp = df_cut.Histo1D({"hQ2temp", Form("hQ2_kin%d", kin), 100, 0.5, 8.5}, "Q2", "fweight");
    double avgQ2 = hQ2temp->GetMean();
    double rmsQ2 = hQ2temp->GetRMS();    
    hQ2s.push_back(hQ2temp);
    // z
    auto hztemp = df_cut.Histo1D({"hztemp", Form("hz_kin%d", kin), 100, 0, 1}, "z", "fweight");
    double avgz = hztemp->GetMean();
    double rmsz = hztemp->GetRMS();        
    hzs.push_back(hztemp);
    // // y_lab
    // auto hytemp = df_cut.Histo1D({"hytemp", Form("hy_kin%d", kin), 100, 1.4, 3.9}, "y_breit", "fweight");
    // double avgy = hytemp->GetMean();
    // double rmsy = hytemp->GetRMS();        
    // hys.push_back(hytemp);
    // y_breit
    auto hytemp = df_cut.Histo1D({"hytemp", Form("hy_kin%d", kin), 100, -2, 2}, "y_breit", "fweight");
    double avgy = hytemp->GetMean();
    double rmsy = hytemp->GetRMS();        
    hys.push_back(hytemp);    

    outfile_data << Form("%f,%f,%f,%f,%f,%f,%f,%f\n",avgx,rmsx,avgQ2,rmsQ2,avgz,rmsz,avgy,rmsy);
  }

  // Optional: Draw or write histos
  TCanvas *c = new TCanvas("c", "", 800, 600);
  //c->Divide(2, 2);
  //for (size_t i = 0; i < histos.size(); ++i) {
    //c->cd(i + 1);
    hys[0]->DrawClone("HIST");
    //}
}

void plot_acceptance() {
  LookUpTableReader runplan;
  runplan.readCSV("R-SIDIS_run_plan.csv");
  
  ROOT::RDataFrame df("h10", "sidis18_allkin_pip_lh2.root");  // tree name and file
  auto df_new = df
    .Define("xacc","sqrt(pt2)*cos(phipq)")
    .Define("yacc","sqrt(pt2)*sin(phipq)");
  
  std::vector<ROOT::RDF::RResultPtr<TH2D>> h2acc1, h2acc2, h2acc3;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> h1acc1, h1acc2, h1acc3;  

  std::vector<int> kin1 = {1,3,19,23,26,29}; // Ebeam = 6.5 GeV
  std::vector<int> kin2 = {2,4,5,7,9,11,13,15,20,24,27,30}; // Ebeam = 8.6 GeV
  std::vector<int> kin3 = {6,8,10,12,14,16,17,18,21,22,25,28,31,32,33,34,35,36,37}; // Ebeam = 10.7 GeV

  // 6.5 GeV
  for (int kin : kin1) {
    double ebeam = runplan.GetValueByKey(kin,0);    
    double x = runplan.GetValueByKey(kin,1);
    double z = runplan.GetValueByKey(kin,3);
    double thpi = runplan.GetValueByKey(kin,9);
    std::string htitle = Form("Kin=%d, E_{e}=%.1f, x=%.2f, z=%.2f, #theta_{#pi}=%.1f",kin,ebeam,x,z,thpi);
    
    auto cutName = Form("deltacut&&kin==%d",kin);    
    auto df_cut = df_new.Filter(cutName);
    auto h2temp = df_cut.Histo2D({Form("h2acc_6p5GeV_kin%d",kin),htitle.c_str(),
	100,-1,.6,100,-.6,.6},"xacc","yacc","fweight");
    h2temp->SetStats(0);
    h2acc1.push_back(h2temp);
    auto h1temp = df_cut.Histo1D({Form("h1acc_6p5GeV_kin%d",kin),htitle.c_str(),
	100,-2,2},"y_breit","fweight");
    h1temp->SetStats(0);
    h1acc1.push_back(h1temp);    
  }

  // 8.6 GeV
  for (int kin : kin2) {
    double ebeam = runplan.GetValueByKey(kin,0);    
    double x = runplan.GetValueByKey(kin,1);
    double z = runplan.GetValueByKey(kin,3);
    double thpi = runplan.GetValueByKey(kin,9);
    std::string htitle = Form("Kin=%d, E_{e}=%.1f, x=%.2f, z=%.2f, #theta_{#pi}=%.1f",kin,ebeam,x,z,thpi);
    
    auto cutName = Form("deltacut&&kin==%d",kin);    
    auto df_cut = df_new.Filter(cutName);
    auto h2temp = df_cut.Histo2D({Form("h2acc_8p6GeV_kin%d",kin),htitle.c_str(),
	100,-1,.6,100,-.6,.6},"xacc","yacc","fweight");
    h2temp->SetStats(0);
    h2acc2.push_back(h2temp);
    auto h1temp = df_cut.Histo1D({Form("h1acc_8p6GeV_kin%d",kin),htitle.c_str(),
	100,-2,2},"y_breit","fweight");
    h1temp->SetStats(0);
    h1acc2.push_back(h1temp);        
  }

  // 10.7 GeV
  for (int kin : kin3) {
    double ebeam = runplan.GetValueByKey(kin,0);    
    double x = runplan.GetValueByKey(kin,1);
    double z = runplan.GetValueByKey(kin,3);
    double thpi = runplan.GetValueByKey(kin,9);
    std::string htitle = Form("Kin=%d, E_{e}=%.1f, x=%.2f, z=%.2f, #theta_{#pi}=%.1f",kin,ebeam,x,z,thpi);
    
    auto cutName = Form("deltacut&&kin==%d",kin);    
    auto df_cut = df_new.Filter(cutName);    
    auto h2temp = df_cut.Histo2D({Form("h2acc_10p7GeV_kin%d",kin),htitle.c_str(),
	100,-1,.6,100,-.6,.6},"xacc","yacc","fweight");
    h2temp->SetStats(0);
    h2acc3.push_back(h2temp);
    auto h1temp = df_cut.Histo1D({Form("h1acc_10p7GeV_kin%d",kin),htitle.c_str(),
	100,-2,2},"y_breit","fweight");
    h1temp->SetStats(0);
    h1acc3.push_back(h1temp);        
  }

  // Draw histos: Ee = 6.5 GeV
  TCanvas *ckin1 = new TCanvas("ckin1","kin1",1200,900);
  ckin1->Divide(3,2);
  for (size_t i = 0; i < h2acc1.size(); ++i) {
    ckin1->cd(i + 1);
    gPad->SetGridx();
    gPad->SetGridy();    
    //h2acc1[i]->DrawClone("colz");
    h1acc1[i]->DrawClone("HIST");    
  }

  // Draw histos: Ee = 8.6 GeV
  TCanvas *ckin2_1 = new TCanvas("ckin2_1","kin2_1",1200,900);
  ckin2_1->Divide(3,2);
  for (size_t i = 0; i < 6; ++i) {
    ckin2_1->cd(i + 1);
    gPad->SetGridx();
    gPad->SetGridy();    
    //h2acc2[i]->DrawClone("colz");
    h1acc2[i]->DrawClone("HIST");    
  }
  //
  TCanvas *ckin2_2 = new TCanvas("ckin2_2","kin2_2",1200,900);
  ckin2_2->Divide(3,2);
  for (size_t i = 0; i < 6; ++i) {
    ckin2_2->cd(i + 1);
    gPad->SetGridx();
    gPad->SetGridy();    
    //h2acc2[i+6]->DrawClone("colz");
    h1acc2[i+6]->DrawClone("HIST");    
  }

  // Draw histos: Ee = 10.7 GeV
  TCanvas *ckin3_1 = new TCanvas("ckin3_1","kin3_1",1200,900);
  ckin3_1->Divide(3,2);
  for (size_t i = 0; i < 6; ++i) {
    ckin3_1->cd(i + 1);
    gPad->SetGridx();
    gPad->SetGridy();    
    //h2acc3[i]->DrawClone("colz");
    h1acc3[i]->DrawClone("HIST");    
  }
  //
  TCanvas *ckin3_2 = new TCanvas("ckin3_2","kin3_2",1200,900);
  ckin3_2->Divide(3,2);
  for (size_t i = 0; i < 6; ++i) {
    ckin3_2->cd(i + 1);
    gPad->SetGridx();
    gPad->SetGridy();    
    //h2acc3[i+6]->DrawClone("colz");
    h1acc3[i+6]->DrawClone("HIST");    
  }
  //
  TCanvas *ckin3_3 = new TCanvas("ckin3_3","kin3_3",1200,900);
  ckin3_3->Divide(3,2);
  for (size_t i = 0; i < 6; ++i) {
    ckin3_3->cd(i + 1);
    gPad->SetGridx();
    gPad->SetGridy();    
    //h2acc3[i+12]->DrawClone("colz");
    h1acc3[i+12]->DrawClone("HIST");    
  }
  //
  TCanvas *ckin3_4 = new TCanvas("ckin3_4","kin3_4",1200,900);
  ckin3_4->Divide(3,2);
  for (size_t i = 0; i < 1; ++i) {
    ckin3_4->cd(i + 1);
    gPad->SetGridx();
    gPad->SetGridy();    
    //h2acc3[i+18]->DrawClone("colz");
    h1acc3[i+18]->DrawClone("HIST");    
  }

  //char const * outfilebase = Form("output_plot_acceptance");
  char const * outfilebase = Form("output_plot_y_breit");  
  TFile *fout = new TFile(Form("%s.root",outfilebase), "RECREATE");
  fout->cd();
  ckin1->SaveAs(Form("%s.pdf[",outfilebase));
  ckin1->SaveAs(Form("%s.pdf",outfilebase)); ckin1->Write();
  ckin2_1->SaveAs(Form("%s.pdf",outfilebase)); ckin2_1->Write();
  ckin2_2->SaveAs(Form("%s.pdf",outfilebase)); ckin2_2->Write();
  ckin3_1->SaveAs(Form("%s.pdf",outfilebase)); ckin3_1->Write();
  ckin3_2->SaveAs(Form("%s.pdf",outfilebase)); ckin3_2->Write();
  ckin3_3->SaveAs(Form("%s.pdf",outfilebase)); ckin3_3->Write();
  ckin3_4->SaveAs(Form("%s.pdf",outfilebase)); ckin3_4->Write();    
  ckin3_4->SaveAs(Form("%s.pdf]",outfilebase));      
}
