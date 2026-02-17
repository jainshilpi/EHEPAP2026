#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TObject.h"
#include "TPostScript.h"
#include <TRandom3.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>

//#define SMEARING
#define PIGAMMASEP
using namespace std;
using namespace CLHEP;

struct phoparam {
  Hep3Vector ppho;
  int nhits;
  double enr1;
  double enr9;
  double enr25;
  double dthe2;
  double dphi2;
  double diag1;
  double diag2;  
};

struct emalgoparam{

  double mu_theta;
  double mu_phi;
  double amp;

};

vector<phoparam> allphoton;
vector<emalgoparam> allclusters_EMalgo;

const double pival = acos(-1.0);
const double twopi = 2*pival;
const double pibytwo = pival/2.;
const double pibyfour = pival/4.;
const double degtorad= pival/180.;
const double radtodeg= 180./pival;
static unsigned int mypow_2[32];

const int nNoiseHit=-100; //Constant noise hit in the detector;




  const double highestPhoThr = 0.5; ///GeV
const double allOtherPhoThr = 0.05; // GeV



//const double seedthres=1.0; //Threshold energy of the centre crystal
const double seedthres=0.5; //Threshold energy of the centre crystal
double pedwidth = 0.040; //40 MeV threshold peak
//const double otherthrs = 3.0*pedwidth; //three sigma threshold
double otherthrs = 3.0*pedwidth; //three sigma threshold
//const double otherthrs = Nsigma*pedwidth; //Nsigma sigma threshold


void fill_digipos_array(int nhit, unsigned int* detid, float pedwidth, float thesh, float& simenr, float& digienr, vector<Hep3Vector>& hitpos);

void ecal_clustering(vector<Hep3Vector> digienr, double seedthres, double thresh, vector<phoparam>& recopho, vector< vector < Hep3Vector > >& cluster_digihitposMod); /// SJ modified

bool find_number(int nn, vector<int>array);



TRandom3* gRandom3 = new TRandom3();






int main(int argc, char* argv[]) {
  /*
g++ `root-config --cflags` `root-config --glibs` -o test1 test1.C 

            - C : a character string terminated by the 0 character
            - B : an 8 bit signed integer (Char_t)
            - b : an 8 bit unsigned integer (UChar_t)
            - S : a 16 bit signed integer (Short_t)
            - s : a 16 bit unsigned integer (UShort_t)
            - I : a 32 bit signed integer (Int_t)
            - i : a 32 bit unsigned integer (UInt_t)
            - F : a 32 bit floating point (Float_t)
            - D : a 64 bit floating point (Double_t)
            - L : a 64 bit signed integer (Long64_t)
            - l : a 64 bit unsigned integer (ULong64_t)
            - O : a boolean (Bool_t)


*/
  //Initialise all these numbers and store in this array, calculaton of power again and again
  //is a CPU consuming process. 

  
  for (int ij=0; ij<32; ij++) {
    mypow_2[ij] = (int)pow(float(2), ij);
  }


  if (argc < 4 ) {
    // Insufficient command-line arguments
    std::cerr << "Usage: " << argv[0] << " Nsigma noise inputfileName" << std::endl;
    return 1;  // Exit with an error code
  }
  
  // Convert command-line arguments to integer and float
  double Nsigma = std::atof(argv[1]); 
  pedwidth = std::atof(argv[2]);    
  otherthrs = Nsigma * pedwidth;
  char rootfiles[1000];
  std::strcpy(rootfiles, argv[3]);
  cout<<"##################################Running for Noise : Nsigma "<<pedwidth<<" "<<Nsigma<<" #######################"<<endl;
  //  cout<<"mean "<<endl;
  gStyle->SetOptStat(111111);


  unsigned   irun;                // Run number of these events
  unsigned   ievt;                //Event number
  unsigned   ngent;
  float	ievt_wt;		//

  const unsigned int ngenmx=50;
  int   pidin[ngenmx]; 	  //PID of incident particle
  float momin[ngenmx]; 	  //Energy of incident particle
  float thein[ngenmx];	  //Initial polar angle of incident particle
  float phiin[ngenmx];     //Initial azimuthal angle of incident particle 
  float posxin[ngenmx];	  //Initial X-position
  float posyin[ngenmx];     //Initial Y-position
  float poszin[ngenmx];     //Initial Z-position

  // For Simulation output of ECAL
  unsigned int nsimhtEC;
  static const unsigned int nsimhtmxEC=2000;
  unsigned int detidEC[nsimhtmxEC];
  float energyEC[nsimhtmxEC];
  float thetaEC[nsimhtmxEC];
  float phiEC[nsimhtmxEC]; 

  unsigned int nsimhtTk; // Looked only for early showering

  unsigned   nphoton;

  float genE = -99;
  
#ifdef PIGAMMASEP
  float enpho;   //Measured energy of reconstrued photon
  float thpho;	  //Measured polar angle of photon
  float phpho;   //Measured azimuthal angle of photon
  float fclus;
  float enr1ph;
  float enr9ph;
  float enr25ph;
  float dthe2ph;
  float dphi2ph;
  float diag1ph;
  float diag2ph;

  float mominx;
  float theinx;
  float phiinx;


  ////SJ
  int nclus, nclus_emalgo, nclus_seed;
  float enpho1, thpho1, phpho1, enpho2, thpho2, phpho2, enpho1_emalgo, thpho1_emalgo, phpho1_emalgo, enpho2_emalgo, thpho2_emalgo, phpho2_emalgo;
  float enpho1_seed, thpho1_seed, phpho1_seed, enpho2_seed, thpho2_seed, phpho2_seed;
  
#else

  const unsigned int nphotmx=20;
  float enpho[nphotmx];   //Measured energy of reconstrued photon
  float thpho[nphotmx];	  //Measured polar angle of photon
  float phpho[nphotmx];   //Measured azimuthal angle of photon
  float fclus[nphotmx];
  float enr1ph[nphotmx];
  float enr9ph[nphotmx];
  float enr25ph[nphotmx];
  float dthe2ph[nphotmx];
  float dphi2ph[nphotmx];
  float diag1ph[nphotmx];
  float diag2ph[nphotmx];  
#endif
  /*
  const double seedthreh=1.0; //Threshold energy of the centre crystal
  const double pedwidth = 0.040; //40 MeV threshold peak
  const double otherthrs = 3.0*pedwidth; //three sigma threshold
  */
  vector<Hep3Vector> digihitpos;


  ifstream file_db;
  
  vector<int>tmpv;
  tmpv.push_back(1);
  tmpv.push_back(2);
  cout<<"tmp "<<tmpv[0]<<"" <<tmpv[1]<<endl;
  swap (tmpv[0], tmpv[1]);
  cout<<"tmp "<<tmpv[0]<<"" <<tmpv[1]<<endl;

  file_db.open(rootfiles);  
  while(!(file_db.eof())) {
    int nentrymx=-1;
    


    char datafile_char[100];
    
    file_db >> datafile_char>>nentrymx >> genE;

    string datafile(datafile_char);
    
    if (strstr(datafile_char,"#")) continue;
    
    cout <<"datafile = "<<datafile<<endl;


    // Find the position of the last '/' and '.'
    size_t lastSlash = datafile.find_last_of('/');
    size_t lastDot = datafile.find_last_of('.');

    // Extract the substring between the last '/' and '.'
    std::string result = datafile.substr(lastSlash + 1, lastDot - lastSlash - 1);

    
    //// name the output file here
    //char rootfiles[1000];
    char outfile[100];
    char outfilx[100];
    
    int len = strlen(rootfiles);
    strncpy(outfilx, rootfiles, len-4);
    outfilx[len-4]='\0';
    cout<<"outfilx is "<<outfilx<<endl;
    
#ifdef PIGAMMASEP
    //sprintf (outfile,"%s_singlex.root",outfilx);
    sprintf (outfile,"%s_singlex_%0.1fsigma_noise%0.3f.root",result.c_str(),Nsigma,pedwidth);
#else
    //sprintf (outfile,"%s_multix.root",outfilx);
    sprintf (outfile,"%s_multix_%0.1fsigma_noise%0.3f.root",result.c_str(),Nsigma,pedwidth);
#endif
    
    char name[100];
    char title[100];
    
    //  cout<<"1mean "<<endl;
    TFile* fileOut = new TFile(outfile, "recreate");  
    TTree* Tout = new TTree("T2", "test");;
    
    
#ifdef PIGAMMASEP
    Tout->Branch("genE", &genE,"genE/F");
    Tout->Branch("enpho", &enpho,"enpho/F");
    Tout->Branch("thpho", &thpho,"thpho/F");
    Tout->Branch("phpho", &phpho,"phpho/F");  
    Tout->Branch("fclus", &fclus,"fclus/F");
    Tout->Branch("enr1ph", &enr1ph,"enr1ph/F");  
    Tout->Branch("enr9ph", &enr9ph,"enr9ph/F");  
    Tout->Branch("enr25ph", &enr25ph,"enr25ph/F");  
    Tout->Branch("dthe2ph", &dthe2ph,"dthe2ph/F"); 
    Tout->Branch("dphi2ph", &dphi2ph,"dphi2ph/F");  
    Tout->Branch("diag1ph", &diag1ph,"diag1ph/F");  
    Tout->Branch("diag2ph", &diag2ph,"diag2ph/F");  
    
    Tout->Branch("mominx",&mominx,"mominx/F");
    Tout->Branch("theinx",&theinx,"theinx/F");
    Tout->Branch("phiinx",&phiinx,"phiinx/F");

    ///SJ - info of the first two highest pT photons for pi0 mass reconstruction as an e.g.

    Tout->Branch("nclus",&nclus,"nclus/I");
    
    Tout->Branch("enpho1",&enpho1,"enpho1/F");
    Tout->Branch("thpho1",&thpho1,"thpho1/F");
    Tout->Branch("phpho1",&phpho1,"phpho1/F");

    Tout->Branch("enpho2",&enpho2,"enpho2/F");
    Tout->Branch("thpho2",&thpho2,"thpho2/F");
    Tout->Branch("phpho2",&phpho2,"phpho2/F");


    Tout->Branch("nclus_emalgo",&nclus_emalgo,"nclus_emalgo/I");
    Tout->Branch("enpho1_emalgo",&enpho1_emalgo,"enpho1_emalgo/F");
    Tout->Branch("thpho1_emalgo",&thpho1_emalgo,"thpho1_emalgo/F");
    Tout->Branch("phpho1_emalgo",&phpho1_emalgo,"phpho1_emalgo/F");

    Tout->Branch("enpho2_emalgo",&enpho2_emalgo,"enpho2_emalgo/F");
    Tout->Branch("thpho2_emalgo",&thpho2_emalgo,"thpho2_emalgo/F");
    Tout->Branch("phpho2_emalgo",&phpho2_emalgo,"phpho2_emalgo/F");


    Tout->Branch("nclus_seed",&nclus_seed,"nclus_seed/I");
    
    Tout->Branch("enpho1_seed",&enpho1_seed,"enpho1_seed/F");
    Tout->Branch("thpho1_seed",&thpho1_seed,"thpho1_seed/F");
    Tout->Branch("phpho1_seed",&phpho1_seed,"phpho1_seed/F");

    Tout->Branch("enpho2_seed",&enpho2_seed,"enpho2_seed/F");
    Tout->Branch("thpho2_seed",&thpho2_seed,"thpho2_seed/F");
    Tout->Branch("phpho2_seed",&phpho2_seed,"phpho2_seed/F");

    const int nvar=9;
    const char* varname[nvar]={"Energy (GeV)", "# of cluster", "e1/e9", "e9/e25", "e25/E?{clust}","1000*d#theta^2", "1000*d#phi^2", "1000*Major", "1000*Minor"};
    float alow[nvar] = { 5.0,  3.5, 0.1, 0.5, 0.8, 0.0, 0.0, 0.0, 0.0};
    float ahgh[nvar] = {50.0, 33.5, 1.001, 1.001, 1.001, 0.5, 1.0, 2.0, 0.2};
    
    
    TH1F* hist_sel[nvar][3]; // [3] : 0 for without any cut, 1 one by one cut, 2 for (n-1) crit
    for (int ij=0; ij <nvar; ij++) {
      int nbin = (ij==1) ? 30 : 120;
      for (int jk=0; jk<3; jk++) {
	sprintf(name, "hist_sel_%i_%i", ij, jk);
	sprintf(title, "%s %s", varname[ij], (jk==0) ? "All" : ((jk==1) ? "any-by-one" :"(n-1) crit"));
	hist_sel[ij][jk] = new TH1F(name, title, nbin, alow[ij], ahgh[ij]);
      }
    }
    /// SJ
    TH2F *h2d_orig = new TH2F("h2d_orig", "h2d_orig", 130,0,130, 100, -50, 50);

    
    TH2F *h2d_cluster = new TH2F("h2d_cluster", "h2d_cluster", 130,0,130, 100, -50, 50);
    TH2F *h2d_emalgo = new TH2F("h2d_emalgo", "h2d_emalgo", 130,0,130, 100, -50, 50);
    TH1F* hTheta_digi = new TH1F("hTheta_digi", "", 500,0,2);
    
#else
    Tout->Branch("irun",&irun,"irun/i");
    Tout->Branch("ievt",&ievt,"ievt/i");
    
    Tout->Branch("ngent",&ngent,"ngent/i");
    Tout->Branch("pidin",pidin,"pidin[ngent]/I");
    Tout->Branch("ievt_wt",&ievt_wt,"ievt_wt/F"); //Generated events with different weight
    Tout->Branch("momin",momin,"momin[ngent]/F");
    Tout->Branch("thein",thein,"thein[ngent]/F");
    Tout->Branch("phiin",phiin,"phiin[ngent]/F");
    Tout->Branch("posxin",posxin,"posxin[ngent]/F");
    Tout->Branch("posyin",posyin,"posyin[ngent]/F");
    Tout->Branch("poszin",poszin,"poszin[ngent]/F");
    
    Tout->Branch("nphoton",&nphoton,"nphoton/i");
    Tout->Branch("enpho",enpho,"enpho[nphoton]/F");
    Tout->Branch("thpho",thpho,"thpho[nphoton]/F");
    Tout->Branch("phpho",phpho,"phpho[nphoton]/F");  
    Tout->Branch("fclus",fclus,"fclus[nphoton]/F");
    Tout->Branch("enr1ph",enr1ph,"enr1ph[nphoton]/F");  
    Tout->Branch("enr9ph",enr9ph,"enr9ph[nphoton]/F");  
    Tout->Branch("enr25ph",enr25ph,"enr25ph[nphoton]/F");  
    Tout->Branch("dthe2ph",dthe2ph,"dthe2ph[nphoton]/F");  
    Tout->Branch("dphi2ph",dphi2ph,"dphi2ph[nphoton]/F");  
    Tout->Branch("diag1ph",diag1ph,"diag1ph[nphoton]/F");  
    Tout->Branch("diag2ph",diag2ph,"diag2ph[nphoton]/F");
    
    TH2F* egamma_cal[4];
    for (int ij=0; ij<4; ij++) {
      sprintf(name, "egamma_cal_%i", ij);
      sprintf(title, "egamma_cal, ntrk==%i", 4*ij);
      egamma_cal[ij] = new TH2F(name, title, 120, 2., 100., 100, -20.0, 0.);
    }
    
#endif
    
    
    
    ////////////////////////////////
    
    
    TFile* fileIn = new TFile(datafile.c_str(), "read");
    TTree *Tin = (TTree*)fileIn->Get("T1");

    Tin->SetBranchAddress("irun",&irun);
    Tin->SetBranchAddress("ievt",&ievt);
    
    Tin->SetBranchAddress("ngent",&ngent);
    Tin->SetBranchAddress("pidin",pidin);
    Tin->SetBranchAddress("ievt_wt",&ievt_wt);
    Tin->SetBranchAddress("momin",momin);
    Tin->SetBranchAddress("thein",thein);
    Tin->SetBranchAddress("phiin",phiin);
    Tin->SetBranchAddress("posxin",posxin);
    Tin->SetBranchAddress("posyin",posyin);
    Tin->SetBranchAddress("poszin",poszin);

    //SImulation output of ECAL
    Tin->SetBranchAddress("nsimhtEC", &nsimhtEC);
    Tin->SetBranchAddress("detidEC", detidEC);
    Tin->SetBranchAddress("energyEC", energyEC);
    Tin->SetBranchAddress("thetaEC", thetaEC);
    Tin->SetBranchAddress("phiEC", phiEC);
    
    Tin->SetBranchAddress("nsimhtTk", &nsimhtTk);


    //////////////////////////////////////

  
    int nentries = Tin->GetEntries();
    cout <<"nentr "<< datafile<<" "<<nentries<<endl;;

    ///SJ
    bool fill2D = false;
    bool fill2D_orig = false;
    int nEvFilled = -99; //if I put it to 0, it fills the first event as well in no clus and EM algo histogram
    
    for (int iev=0; iev<min(nentries, nentrymx); iev++) {
      //fileIn->cd();
      Tin->GetEntry(iev);
      cout<<" SimHit Size "<<nsimhtTk<<" "<< momin[0]<<" "<<pidin[0]<<" "<<digihitpos.size()<<endl;
      if (nsimhtEC==0) continue;

      float totsimenr =0;
      float totdigienr =0;
      digihitpos.clear();
      
      
      fill_digipos_array(nsimhtEC, detidEC, pedwidth, otherthrs, totsimenr, totdigienr, digihitpos);

      // Use clustering algorithm to find photons
      if (digihitpos.size()==0 || digihitpos[0].mag()<seedthres) continue;

      vector<vector<Hep3Vector>> cluster_digihitposMod;

      
      allphoton.clear();      
      ecal_clustering(digihitpos, seedthres, otherthrs, allphoton, cluster_digihitposMod);



    
  } //end of file loop

  ///set proper range of the TH2Ds
  
  
  


}
////////////////////////////////////////////////////////////////////////////////////
////
////                            Digitisation
////
////////////////////////////////////////////////////////////////////////////////////
void fill_digipos_array(int nhit, unsigned int* detid, float pedwidth, float thresh, float& simenr, float& digienr, vector<Hep3Vector>& hitpos) {
  //Decode detid to position and energy and efficiency to select hitpoint for fit 
  Hep3Vector tmp3vect;

   //Randomly generated hit position and fill array
  // In real case, it should be indigitised form and also depend on detector condition
#ifdef SMEARING  
  for (int ithe=0; ithe<128; ithe++) { //theta
    for (int iphi=0; iphi<128; iphi++) { //phi
      double enr = gRandom3->Gaus(0, thresh/3.);

      for (int ij =0; ij <nhit; ij++) {
	if ((detid[ij]>>18)==(128*ithe + iphi)) {
	  double tmpen = (detid[ij]&0x3ffff)/1000.0;
	  simenr += tmpen;
	  enr += tmpen;
	  break;
	}
      }
      
      if (enr<thresh) continue;
      digienr +=enr;
      double the = (ithe + 45.5)*degtorad;
      double phi = (iphi - 44.5)*degtorad;
      
      tmp3vect.setRThetaPhi(enr, the, phi);
      
      hitpos.push_back(tmp3vect);
    }
  }
#else 
  
  for (int ij = 0;  ij<nhit; ij++) {
    int ienr = (detid[ij]&0x3ffff);     //Digienergy in MeV
    //One  can add noise here to, but that will be biased due to already ped suppression
    double enr = ienr/1000.0; //convert to GeV
    
    simenr += enr;
    if (enr <thresh) continue;
    digienr +=enr;

    int ithe = ((detid[ij]>>25)&0x7f);  //itheta during coding
    int iphi = ((detid[ij]>>18)&0x7f);  //iphi during coding
    
    double the = (ithe + 45.5)*degtorad;
    double phi = (iphi - 44.5)*degtorad;
    
    tmp3vect.setRThetaPhi(enr, the, phi);
    
    hitpos.push_back(tmp3vect);
  }
  
#endif

  //order the hits according to energy
  //use bubble sort methor
  //Must be done here, otherwise need to change algorithm for shape variables
  for (int ij=0; ij<hitpos.size(); ij++) {
    for (int jk=ij+1; jk<hitpos.size(); jk++) {
      if (hitpos[ij].mag() < hitpos[jk].mag()) {
	//	hitpos.swap(hitpos.begin()+ij, hitpos.begin()+jk);
	swap(hitpos[ij], hitpos[jk]);
      }
    }
  }
  
  //  for (int ij=0; ij<hitpos.size(); ij++) {
  //    cout<<"energy tower "<< hitpos[ij]<<" "<<hitpos[ij].mag()<<endl;
  //  }



}
////////////////////////////////////////////////////////////////////////////////////
////
////       Simple clustering of photon with seed
////
////////////////////////////////////////////////////////////////////////////////////
void ecal_clustering(vector<Hep3Vector> digienr, double seedthres, double thresh, vector<phoparam>& recopho, vector< vector < Hep3Vector > >& cluster_digihitposMod) {

  int nsize = digienr.size();
  vector<Hep3Vector> digienr_org = digienr;

  double meanmm=0;
  double meancc=0;
  double meanmm2=0;
  double meancc2=0;
  
  //A simple logic, but not an optimised one in the context of CPU
  // Use map/boundary logic to reduce the loop

  vector < vector < Hep3Vector > > allclusters; 
  vector < Hep3Vector > cluster;

  

}


bool find_number(int nn, vector<int>array) {
  for (int ij=0; ij<array.size(); ij++) {
    if (nn==array[ij]) return true;
  }
  return false;
}


