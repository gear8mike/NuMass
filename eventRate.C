{
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat("");

  gStyle->SetLabelFont(102,"");
  gStyle->SetLabelSize(0.06,"");
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.001,"x");
  gStyle->SetLabelOffset(0.01,"y");

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(0.9,"x");
  gStyle->SetTitleOffset(1.2,"y");

  gStyle->SetStripDecimals(kFALSE);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.25);

  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(99);

  gStyle->SetHistLineWidth(2);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFuncWidth(2);

  gStyle->SetStatFont(42);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  gStyle->SetOptStat(000000);


  // estimated neutrion cross section: https://www.semanticscholar.org/paper/From-eV-to-EeV%3A-Neutrino-Cross-Sections-Across-Formaggio-Zeller/7ff94e6c9f4359d11f63e674a27a81c72be5191a/figure/0
  double xsec[11] = {1.e-28, 1.e-27, 1.e-24, 1.e-22, 1.e-21, 1.e-19, 1.e-17, 1.e-16, 1.e-15, 1.e-14, 1.e-13 };
  double ixsec[11] = {1.e0, 1.e1, 1.e2, 1.e3, 1.e4,1.e5,1.e6,1.e7,1.e8,1.e9,1.e10};
  for (int i = 0; i<11; i++) { ixsec[i] *= 1e-9; }

  auto g1 = new TGraph(11,ixsec, xsec);

  cout<<1<<endl;
  double dist = 1.5e8; //km
  double eff = 1;
  double detMass;
  double conv = 1; 
  double dune_pot = 7.5e13; 
  double sunMass = 2.192e27;
  double dune_perTon_perPOT = 1.5e7 * (7.5e13/1.1e21) / 8./ dune_pot * 100.; // 100 is to assume the beam is very concentrated on the detector
  //https://indico.bnl.gov/event/6284/contributions/29058/attachments/23985/35170/diwan-nu-basics-interactions2.pdf
  long double back_int = 1. / (1e5*1e5) * 1.e-27;// event per bunch
  double targetN = 6.e29; // 1 ton has 6e29 protons and neutrons;
  
  double power;
  double backangle = 0.5 * (TMath::Pi() * 0.1 * 0.1) / (4*TMath::Pi() * dist * dist);

  cout<<2<<endl;
  ifstream in;
  in.open("oscillation.txt");
  double aa ; 
  double bb;
  int count = 0; 
  double osc[100] = {};
  while (!in.eof()){
    in>>aa>>bb; // energy in gev; ener = i*0.001
    cout<<aa<<" "<<bb<<endl;
    osc[(int)aa] = bb;
    count ++;
  }

  double upperPower=2.e27;
  TH2F* h1 = new TH2F("","Event rate per bunch; POT/bunch x detector mass in tonne; Neutrino energy (MeV)", 100, 0, upperPower, 100, 0, 500);
  TH2F* h2 = new TH2F("","Event rate per bunch w/o oscillation; POT/bunch x detector mass in tonne; Neutrino energy (MeV)", 100, 0, upperPower, 100, 0, 500);
  cout<<3<<endl;
  double N[100][100];
  double N2[100][100];
  for (int i =0;i<100;i++){
    for (int j = 0; j< 100; j++){
      power = upperPower* (i/100.);  // this number is POT/bunch x detector mass in tonne
      double hxsec = g1->Eval(j*0.005); // in unit of GeV
      // N is the detected neutrino per POT 
      N[i][j] = dune_perTon_perPOT * (g1->Eval(j*0.005)/g1->Eval(1)) * power * sunMass * backangle * osc[j]* back_int * hxsec * targetN;
      N2[i][j] = dune_perTon_perPOT * (g1->Eval(j*0.005)/g1->Eval(1)) * power * sunMass * backangle * back_int * hxsec * targetN;
      cout<<"power "<< 7.5e13*i<<" POT per bunch, at energy "<<j*0.001<<"  event number:  "<<N[i][j]<<endl;
      h1->SetBinContent(i+1,j+1, N[i][j]);
      h2->SetBinContent(i+1,j+1, N2[i][j]);
    }
  }

  new TCanvas();
  h1->Draw("colz");

  new TCanvas();
  h2->Draw("colz");

  double biny1[101] = {};
  double biny2[101] = {};
  for (int i=0;i<101;i++){
    biny1[i] = TMath::Power(10., -7 + (3* i/100.));
    biny2[i] = TMath::Power(10., -11. + (3*i/100.));
  }
  TH2F* h4 = new TH2F("","Mass limit (68%); Mass (MeV); Timing resolution (s)",100,biny1, 100, biny2);

  //TH1F* h3 = new TH1F("","",100,0.06e-6,0.12e-6);
  for (int mm=0;mm<100;mm++){ //mass test loop
  double c = 299792458.; // m/s
  double m1 = TMath::Power(10., -7 + (3* mm/100.)); // MeV
  double m2 = 0.000000000001e-6;
  double test_ke = 10; // mev
  double test_l = 3.e11; // in m

  double v3 = sqrt(1- pow(m1/(test_ke + m1) ,2) )*c;
  double ToF3 = test_l/v3;
  
  double v4 = sqrt(1- pow(m2/(test_ke + m2) ,2) )*c;
  double ToF4 = test_l/v4;
  double timeS;
  double test_lv;

  for (double j=0;j<100;j++){ // timeS resolution
    TH1F* h3 = new TH1F("","",100,0, 1.);
    double it = TMath::Power(10, -11. + (5*j/100.) );
    for (int i =0;i<30000; i++){ // bunch number
      int nip = gRandom->Poisson(N[50][2] );
      for (int ip = 0;ip<nip;ip++){
        timeS = gRandom->Gaus(ToF3 , it );
        
	test_lv = test_l + gRandom->Gaus(0, 696.34e3);
        //cout<<sqrt(1-TMath::Power(((test_l/ timeS ) / c ),2)  )<<"   "<<test_ke/ ( 1/(sqrt(1-TMath::Power(((test_l/ timeS ) / c ),2)  )) -1  )<<endl;
        double currM = test_ke/ ( 1/(sqrt(1-TMath::Power(((test_lv/ timeS ) / c ),2)  )) -1  );
        h3->Fill(currM);
      }
    }
    cout<<"m1, timeS, "<<m1<<" "<<it<<" "<<"mean and rms : "<<h3->GetMean()<<" "<<h3->GetRMS()<<"     ToF3 - ToF4 "<<ToF3-ToF4<<endl;
    //if (h3->GetMean()+h3->GetRMS()  )
    h4 ->SetBinContent(mm+1, j+1, h3->GetMean()+1*h3->GetRMS());
    h3->Clear();
    h3->Delete();
  }

  }
  new TCanvas();
  h4->Draw("colz");

}
