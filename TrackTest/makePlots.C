void makePlots(char * file, char * suffix) {
  //Get rid of stat-box
  gROOT->GetStyle("Plain")->SetOptStat(0);
  // AXIS LABELS STYLE
  gROOT->GetStyle("Plain")->SetLabelColor(kBlack, "XYZ");
  gROOT->GetStyle("Plain")->SetLabelFont(43, "XYZ");
  gROOT->GetStyle("Plain")->SetLabelSize(20, "XYZ");
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle(kTRUE);

  const char * folder = "demo";
  const char * to_plot[] = {"TrackerClusterMap", "TrackerClusterOntrackMap", "TrackerSameClusterOntrackMap", "TrackerOccupancyOntrackMap", "TrackerOccupancySameClusterOntrackMap", 0};
  TFile * f = new TFile(file);
  TCanvas * c = new TCanvas("C", "C", 1024, 1024);
  for (const char **profile = to_plot; *profile; ++profile) {
    TString obj_name(folder);
    obj_name.Append("/").Append(*profile);
    TProfile * p = (TProfile*)f->Get(obj_name.Data());
    p->Draw("colz");
    TString filename(*profile);
    filename.Append("_").Append(suffix).Append(".png");
    c->SaveAs(filename.Data());
  }
}
