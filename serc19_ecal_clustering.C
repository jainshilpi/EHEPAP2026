//  /usr/include/c++/4.8/iostream
// /usr/include/c++/4.8/bits/stl_vector.h

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
//#include "/Users/shilpi2015/geant4/geant4-v11.2.0/source/externals/clhep/include/CLHEP/Matrix/Matrix.h"
//#include <CLHEP/Matrix/Matrix.h>


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



///thresholds on the photon --> used for filling the tre

///original values
/*const double highestPhoThr = 5; ///GeV
const double allOtherPhoThr = 2.0; // GeV
*/


/// use this --> used for school except for 150 MeV pions
/*
  const double highestPhoThr = 0.5; ///GeV
const double allOtherPhoThr = 0.05; // GeV



//const double seedthres=1.0; //Threshold energy of the centre crystal
const double seedthres=0.5; //Threshold energy of the centre crystal
double pedwidth = 0.040; //40 MeV threshold peak
//const double otherthrs = 3.0*pedwidth; //three sigma threshold
double otherthrs = 3.0*pedwidth; //three sigma threshold
//const double otherthrs = Nsigma*pedwidth; //Nsigma sigma threshold
*/


////150 MeV pion - use above for the rest

const double highestPhoThr = 0.0; ///GeV
const double allOtherPhoThr = 0.0; // GeV


//const double seedthres=1.0; //Threshold energy of the centre crystal
const double seedthres=0.000005; //Threshold energy of the centre crystal
double pedwidth = 0.0; //40 MeV threshold peak
//const double otherthrs = 3.0*pedwidth; //three sigma threshold
double otherthrs = pedwidth; //three sigma threshold
//const double otherthrs = Nsigma*pedwidth; //Nsigma sigma threshold


void fill_digipos_array(int nhit, unsigned int* detid, float pedwidth, float thesh, float& simenr, float& digienr, vector<Hep3Vector>& hitpos);

void ecal_clustering(vector<Hep3Vector> digienr, double seedthres, double thresh, vector<phoparam>& recopho, vector< vector < Hep3Vector > >& cluster_digihitposMod); /// SJ modified

bool find_number(int nn, vector<int>array);

void EMalgo(vector<Hep3Vector> digienr, vector<emalgoparam> &clusters, vector<vector<double>> &hitsFrac_clus);

bool isDigiInAllClusters(Hep3Vector digi, std::vector<Hep3Vector> allclusters);

bool isHighEnergyHitAround(vector<Hep3Vector> digienr, Hep3Vector seed);

void ecal_clustering_seed(vector<Hep3Vector> digienr, vector<vector<Hep3Vector>> &cluster_digihitposMod_seed);
/*
void fcnsg(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t flag) {
  double fval=0;
  for (int ij=0; ij<fit_input3v.size(); ij++) {
    fval += abs((fit_input3v[ij].x()-par[0])*(fit_input3v[ij].x()-par[0]) 
	       +(fit_input3v[ij].y()-par[1])*(fit_input3v[ij].y()-par[1])
	       - par[2]*par[2]);
    //    cout<<" fval = "<<fit_input3v[ij]<<" "<<fval<<endl;
  }
  f = 10000*fval;
  //  cout<<" fvalxxxxxxxxxxxxx = "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<f<<endl;
}
*/
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
    
    
    /*cout <<"Give the input file name"<<endl;
      cin>> rootfiles;
    */
    
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
    
    
    
    /*
    if(file_db.eof()) {
      cout<<"Its end of file"<<endl;
      break; 
    }
    */
    
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
      
      for (int ij=0; ij<ngent; ij++) {
      	cout<<"genmom "<<iev<<" "<<ij<<" "<<momin[ij]<<" "<<radtodeg*thein[ij]<<" "<<radtodeg*phiin[ij]<<endl;
      }
    
      fill_digipos_array(nsimhtEC, detidEC, pedwidth, otherthrs, totsimenr, totdigienr, digihitpos);
      cout<<"Digihit size "<<nsimhtTk<<" "<< digihitpos.size()<<" "<<momin[0]<<" "<<totsimenr<<" "<<totdigienr<<endl;
      //      for (int ij=0; ij<min(18,(int)digihitpos.size()); ij++) {
      //      	cout<<"clusters "<<ij<<" "<<digihitpos[ij]<<" "<<digihitpos[ij].mag()<<" "<<radtodeg*digihitpos[ij].theta()<<" "<<radtodeg*digihitpos[ij].phi()<<endl;
      //      }

      // Use clustering algorithm to find photons
      if (digihitpos.size()==0 || digihitpos[0].mag()<seedthres) continue;

      //cout<<"SJ!!!Passed seed threshold "<<endl;
      //// SJ - its a vector of digis of all the clusters that are formed
      vector<vector<Hep3Vector>> cluster_digihitposMod;
      //// SJ
      
      allphoton.clear();      
      ecal_clustering(digihitpos, seedthres, otherthrs, allphoton, cluster_digihitposMod);

      ///use EMAlgo - SJ
      allclusters_EMalgo.clear();
      vector<vector<double>> hitsFrac_clus;
      EMalgo(digihitpos, allclusters_EMalgo, hitsFrac_clus);


      ///ecal_clustering only but with diff seeds
      vector<vector<Hep3Vector>> cluster_digihitposMod_seed;
      ecal_clustering_seed(digihitpos, cluster_digihitposMod_seed);
      ///SJ
      
      /*cout<<""<<endl;
      
      cout<<"#digis original "<<digihitpos.size()<<endl;
      cout<<"# clusters formed "<<cluster_digihitposMod.size()<<endl;
      */

      nclus = 0;
      enpho1 = thpho1 = phpho1 = enpho2 = thpho2 = phpho2 = 0;
      
      for(int icl=0; icl<cluster_digihitposMod.size(); icl++){
	//cout<<"SJ!!! original clus, iclus "<<icl<<endl;
	
	vector < Hep3Vector > digi = cluster_digihitposMod[icl];
	//cout<<"Number of digis in this "<<icl<<"th cluster is "<<digi.size()<<endl;
	for(int idigi=0; idigi<digi.size(); idigi++){


	  double theta = digi[idigi].theta()*radtodeg;
	  double phi = digi[idigi].phi()*radtodeg;
	  double digiEn = digi[idigi].mag();
	  cout<<"SJ!!! clustered hits theta : phi "<<theta<<" "<<phi<<" and energy is "<<digiEn<<endl;

	  hTheta_digi->Fill(digi[idigi].theta());
	  //fill when original digi size is diff to see the effect of clustering
	  //if(!fill2D && digi.size() > 0 && digi.size()!=digihitpos.size() && allclusters_EMalgo.size()!=cluster_digihitposMod.size() && icl==0) { ///used till Jan 20
	  if(!fill2D && digi.size() > 0 && allclusters_EMalgo.size()!=cluster_digihitposMod.size() && icl==0) {
	  //if(!fill2D && digi.size() > 0 && allclusters_EMalgo.size()!=cluster_digihitposMod.size() && icl==1) {  
	  //if(iev == 311){
	  //if(!fill2D && digi.size() > 0 && digi.size()!=digihitpos.size()) {
	    //cout<<"SJ!!!! filled the histo for event number "<<iev<<endl;
	    nEvFilled = iev;
	    h2d_cluster->Fill(theta, phi, digiEn);
	    //cout<<"Digi info that went inside original clustering and in that histo, theta phi en "<<theta<<" "<<phi<<" "<<digiEn<<endl;
	    if(idigi==digi.size() || idigi==digi.size()-1) fill2D = true;
	  }
	  
	  /////form 2 photons out of it to form the pi0 mass andalso store the nclus
	  if(icl==0){
	    enpho1 += digiEn ;
	    thpho1 += theta*digiEn;
	    phpho1 += phi*digiEn;
	    
	  }

	  if(icl==1){
	    enpho2 += digiEn ;
	    thpho2 += theta*digiEn;
	    phpho2 += phi*digiEn;
	    
	  }
	  
	}//for(int idigi=0; idigi<digi.size(); idigi++)

	
      }/// loop over formed clusters

      if(enpho1 > 0){
	thpho1 = thpho1/enpho1;
	phpho1 = phpho1/enpho1;
      }

      if(enpho2 > 0){
	thpho2 = thpho2/enpho2;
	phpho2 = phpho2/enpho2;
      }

      nclus = cluster_digihitposMod.size();


      
      double totDigE = 0;
      for(int idigi=0; idigi<digihitpos.size(); idigi++){
	double theta = digihitpos[idigi].theta()*radtodeg;
	double phi = digihitpos[idigi].phi()*radtodeg;
	double digiEn = digihitpos[idigi].mag();
	cout<<"SJ!!! All hits theta : phi "<<theta<<" "<<phi<<" and energy is "<<digiEn<<endl;
	
	totDigE += digiEn;
	if(iev==nEvFilled) { //fill for the same event as above
	//if(!fill2D_orig && iev>1){
	  cout<<"Digi info that went inside no clustering theta phi en "<<theta<<" "<<phi<<" "<<digiEn<<endl;

	  
	  h2d_orig->Fill(theta, phi, digiEn);
	  if(idigi==digihitpos.size() || idigi==digihitpos.size()-1) fill2D_orig = true;
	}

      }

      cout<<"Total energy in ECAL in this iev "<<iev << " is = "<<totDigE<<endl;
      ////SJ EM algo
      if(allclusters_EMalgo.size()!= hitsFrac_clus.size()){

	cout<<"SJ!!! Bug somewhere in EM algo code, the allclusters_EMalgo.size! = hitsFrac_clus.size, size of allclusters_EMalgo : size of hitsFrac_clus "<<allclusters_EMalgo.size()<<" "<<hitsFrac_clus.size()<<endl;
      }

      nclus_emalgo = 0;
      enpho1_emalgo = thpho1_emalgo = phpho1_emalgo = enpho2_emalgo = thpho2_emalgo = phpho2_emalgo = 0;

      
      cout<<"===========SJ!!! number of orig digis : number of clusters from std algo : clusters from EM algo "<<" "<<digihitpos.size()<<" "<<cluster_digihitposMod.size()<<" "<<allclusters_EMalgo.size()<<endl;
      for(int iclus=0; iclus<hitsFrac_clus.size(); iclus++){

	vector<double> hitsFrac = hitsFrac_clus[iclus];
	
	for(int idigi=0; idigi<hitsFrac.size(); idigi++){

	  double theta = digihitpos[idigi].theta()*radtodeg;
	  double phi = digihitpos[idigi].phi()*radtodeg;
	  double digiEn = digihitpos[idigi].mag();

	  if(digiEn < otherthrs) continue;
	  
	  cout<<"SJ!! after calling EM algo, fraction of hit associated with "<<iclus<<"th cluster is "<<hitsFrac[idigi]<<" and the energy is "<<digiEn<<endl;
	  //if(iev==nEvFilled && iclus==0) { //fill for the same event as above
	  if(iev==nEvFilled && iclus==0) { //fill for the same event as above
	  //if(iev==nEvFilled && iclus==1) { //fill for the same event as above
	    
	    if(hitsFrac[idigi]!=0) h2d_emalgo->Fill(theta, phi, hitsFrac[idigi]*digiEn);
	  }

	  /// form 2 photons out of it to form the pi0 mass and also store the nclus
	  if(iclus==0){
	    enpho1_emalgo += hitsFrac[idigi]*digiEn ;
	    thpho1_emalgo += theta*hitsFrac[idigi]*digiEn;
	    phpho1_emalgo += phi*hitsFrac[idigi]*digiEn;
	    
	  }

	  if(iclus==1){
	    enpho2_emalgo += hitsFrac[idigi]*digiEn ;
	    thpho2_emalgo += theta*hitsFrac[idigi]*digiEn;
	    phpho2_emalgo += phi*hitsFrac[idigi]*digiEn;
	    
	  }

	}
	
	
	
      }//for(int iclus=0; iclus<hitsFrac_clus.size(); iclus++)

      if(enpho1_emalgo > 0){
	thpho1_emalgo = thpho1_emalgo/enpho1_emalgo;
	phpho1_emalgo = phpho1_emalgo/enpho1_emalgo;
      }

      if(enpho2_emalgo > 0){
	thpho2_emalgo = thpho2_emalgo/enpho2_emalgo;
	phpho2_emalgo = phpho2_emalgo/enpho2_emalgo;
      }

      nclus_emalgo = hitsFrac_clus.size();

      //      cout<<iev<<" E "<<allphoton[0].ppho<<endl;


      nclus_seed = 0;
      enpho1_seed = thpho1_seed = phpho1_seed = enpho2_seed = thpho2_seed = phpho2_seed = 0;
      
      ////ecal clustering but with being careful abt the seeds
      for(int icl=0; icl<cluster_digihitposMod_seed.size(); icl++){
	vector < Hep3Vector > digi = cluster_digihitposMod_seed[icl];
	//cout<<"Number of digis in this "<<icl<<"th cluster is "<<digi.size()<<endl;
	for(int idigi=0; idigi<digi.size(); idigi++){


	  double theta = digi[idigi].theta()*radtodeg;
	  double phi = digi[idigi].phi()*radtodeg;
	  double digiEn = digi[idigi].mag();

	  if(digiEn < otherthrs) continue;
	  
	  cout<<"SJ!! after calling ecal clustering with seed constraint "<<icl<<"th cluster "<<" the energy is "<<digiEn<<endl;
	  
	  /////form 2 photons out of it to form the pi0 mass andalso store the nclus
	  if(icl==0){
	    enpho1_seed += digiEn ;
	    thpho1_seed += theta*digiEn;
	    phpho1_seed += phi*digiEn;
	    
	  }

	  if(icl==1){
	    enpho2_seed += digiEn ;
	    thpho2_seed += theta*digiEn;
	    phpho2_seed += phi*digiEn;
	    
	  }
	  
	}//for(int idigi=0; idigi<digi.size(); idigi++)

	
      }/// loop over formed clusters

      if(enpho1 > 0){
	thpho1_seed = thpho1_seed/enpho1_seed;
	phpho1_seed = phpho1_seed/enpho1_seed;
      }

      if(enpho2 > 0){
	thpho2_seed = thpho2_seed/enpho2_seed;
	phpho2_seed = phpho2_seed/enpho2_seed;
      }

      nclus_seed = cluster_digihitposMod_seed.size();

      
      fileOut->cd();
      nphoton = 0;
#ifdef PIGAMMASEP

      float vars[nvar];
      int ips[nvar]={0};
      int ipsall = 0;

      cout<<"Towards the end, allphoton.size = "<<allphoton.size()<<endl;
      if (allphoton.size()>0 && allphoton[0].ppho.mag() > highestPhoThr) { 
	//only events with energy more than 5 GeV, but keep in mint this is 
	//before any corrections, which is required for noise and also for threshold
	nphoton = 1;
	for (int ij=1; ij<allphoton.size(); ij++) {
	  if (allphoton[0].ppho.mag() > allOtherPhoThr) {
	    //Select events with only one reconstruted photon
	    nphoton++;
	  }
	}
	if (nphoton ==1) { 
	  vars[0] = enpho =  allphoton[0].ppho.mag();
	  thpho = allphoton[0].ppho.theta();	
	  phpho = allphoton[0].ppho.phi();
	  vars[1] = fclus = allphoton[0].nhits;
	  //	  enr1ph = allphoton[0].enr1;
	  //	  enr9ph = allphoton[0].enr9;
	  vars[2] =  enr1ph = allphoton[0].enr1/allphoton[0].enr9;
	  //	  enr25ph = allphoton[0].enr25;
	  vars[3] = enr9ph = allphoton[0].enr9/allphoton[0].enr25;
	  vars[4] = enr25ph = allphoton[0].enr25/enpho;
	  vars[5] = dthe2ph = 1000*allphoton[0].dthe2; //Avoid very small number
	  vars[6] = dphi2ph = 1000*allphoton[0].dphi2;
	  vars[7] = diag1ph = 1000*allphoton[0].diag1;
	  vars[8] = diag2ph = 1000*allphoton[0].diag2;

	  //These can be put is a loop also
	  if (vars[0] > 10 ) { ips[0] =  mypow_2[0]; ipsall += ips[0];}
	  if (vars[1] > 5  ) { ips[1] =  mypow_2[1]; ipsall += ips[1];}
	  if (vars[2] >0.1 ) { ips[2] =  mypow_2[2]; ipsall += ips[2];}
	  if (vars[3] >0.8 ) { ips[3] =  mypow_2[3]; ipsall += ips[3];}
	  if (vars[4] >0.9 ) { ips[4] =  mypow_2[4]; ipsall += ips[4];}

	  if (vars[5] <0.4 ) { ips[5] =  mypow_2[5]; ipsall += ips[5];}
	  if (vars[6] <0.6 ) { ips[6] =  mypow_2[6]; ipsall += ips[6];}
	  if (vars[7] <0.5 ) { ips[7] =  mypow_2[7]; ipsall += ips[7];}
	  if (vars[8] <0.1 ) { ips[8] =  mypow_2[8]; ipsall += ips[8];}
	  for (int ij=0; ij<nvar; ij++) { cout<<" "<<ips[ij];}
	  cout<<" ipasll "<<ipsall<<endl;
	  
	  for (int ij=0; ij<nvar; ij++) {
	    hist_sel[ij][1]->Fill(vars[ij]);  //Filling after previous criterion
	    if (ips[ij]==0) break;
	  }

	  for (int ij=0; ij<nvar; ij++) {
	    hist_sel[ij][0]->Fill(vars[ij]); //All events
	    if (ipsall-ips[ij]==mypow_2[nvar]-1 - mypow_2[ij]) { 
	      hist_sel[ij][2]->Fill(vars[ij]); //N-1 plots
	    }
	  }
	  nphoton = 1;
		
	  mominx = momin[0];
	  theinx = thein[0];
	  phiinx = phiin[0];


	} else {
	  nphoton = 0; //Do not store this event
	}

	cout<<"Towards the end, allphoton.size(2) = "<<allphoton.size()<<endl;
      }
#else
    
      for (int ij=0; ij<allphoton.size(); ij++) {
	if (allphoton[ij].ppho.mag()<1.0) continue; //store photon with minimum energy
	enpho[nphoton] =  allphoton[ij].ppho.mag();
	thpho[nphoton] = allphoton[ij].ppho.theta();	
	phpho[nphoton] = allphoton[ij].ppho.phi();
	fclus[nphoton] = allphoton[ij].nhits;
	enr1ph[nphoton] = allphoton[ij].enr1;
	enr9ph[nphoton] = allphoton[ij].enr9;
	enr25ph[nphoton] = allphoton[ij].enr25;
	dthe2ph[nphoton] = allphoton[ij].dthe2;
	dphi2ph[nphoton] = allphoton[ij].dphi2;
	diag1ph[nphoton] = allphoton[ij].diag1;
	diag2ph[nphoton] = allphoton[ij].diag2;  
	if (++nphoton >=nphotmx) break;
      }
      
      if (nsimhtTk<16) {egamma_cal[int(nsimhtTk/4)]->Fill(momin[0], 100*(enpho[0]-momin[0])/momin[0]);}

      cout<<"Towards the end, allphoton.size(2) = "<<allphoton.size()<<endl;
#endif      

      //if (nphoton>0) {cout<<"Filled this event"<<endl;
      ///changed from above on 17th Feb, 2026 since above is somehow not filling events even when cluster number is 2 or 3
      if (allphoton.size() > 0) {cout<<"Filled this event"<<endl; Tout->Fill();} 
      
    } //end of events loop
    
    fileIn->cd();

    fileOut->cd();
    fileOut->Write();
    
    
    delete Tin;
    
    delete fileIn;


    
  } //end of file loop

  ///set proper range of the TH2Ds
  
  
  

  TCanvas *c1 = new TCanvas("c1","c1",600,800);
  c1->Divide(2,2);

  c1->cd(1);
  //  histx->Draw();

  c1->cd(2);
  // histxy->Draw();
  
  c1->cd(3);
  // pdatax->Draw();

  c1->cd(4);
  // readx->Draw();

}
////////////////////////////////////////////////////////////////////////////////////
////
////                            Digitisation
////
////////////////////////////////////////////////////////////////////////////////////
void fill_digipos_array(int nhit, unsigned int* detid, float pedwidth, float thresh, float& simenr, float& digienr, vector<Hep3Vector>& hitpos) {
  //Decode detid to position and energy and efficiency to select hitpoint for fit 
  Hep3Vector tmp3vect;
  //  cout<<"nhit "<<nhit<<endl;
  //cout<<"SJ!!! Inside fill_digipos_array thresh = "<<thresh<<endl;
   //Randomly generated hit position and fill array
  // In real case, it should be indigitised form and also depend on detector condition
#ifdef SMEARING  
  for (int ithe=0; ithe<128; ithe++) { //theta
    for (int iphi=0; iphi<128; iphi++) { //phi
      //double enr = gRandom3->Gaus(0, thresh/3.);
      //double enr = gRandom3->Gaus(0, thresh/3.);
      //double enr = gRandom3->Gaus(0, 0.04);
      double enr = gRandom3->Gaus(0, thresh/3.);
      //cout<<"en after smearing is "<<enr<<endl;
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

  
  /*cout<<""<<endl;
  cout<<"SJ!!! inside ecal_clustering"<<endl;
  cout<<"# digist "<<(int)digienr.size()<<endl;
  */
  
  for (int ij=0; ij<digienr.size(); ij++) {
    cluster.clear();
    if (digienr[ij].mag()<seedthres) continue;
      
    Hep3Vector tmp3v1 = digienr[ij].unit();
    cluster.push_back(digienr[ij]);
    digienr.erase(digienr.begin()+ij);
    ij--;
    int ncount=1;

    //cout<<"SJ!! ij : cluster size now "<<ij<<" "<<(int)cluster.size()<<endl;
    while (ncount >0 && cluster.size()>0) {
      ncount=0;
      //cout<<"SJ!! inside while loop, cluster size "<<(int)cluster.size()<<endl;
      for (int jk=0; jk<cluster.size(); jk++) {
	//cout<<"SJ!! inside for loop, cluster size "<<(int)cluster.size()<<endl;	
	for (int kl=0; kl<digienr.size(); kl++) {
	  //cout<<"SJ!! inside digiloop, cluster number : digienr.size "<<kl<<" "<<(int)digienr.size()<<endl;
	  if (digienr[kl].mag()<thresh) continue;
	  double angle = cluster[jk].angle(digienr[kl].unit());
	  //	  cout<<"angle "<<" "<<jk<<" "<<kl<<" "<<angle<<" "<<angle*radtodeg<<endl;
	  if (angle*radtodeg <1.5) {
	  //if (angle*radtodeg <0.5) { //to check
	    cluster.push_back(digienr[kl]);
	    digienr.erase(digienr.begin()+kl);
	    kl--;
	    ncount++;
	  }
	}//for (int kl=0; kl<digienr.size(); kl++)
      }//for (int jk=0; jk<cluster.size(); jk++)
    }///while
    allclusters.push_back(cluster);
  }///loop over digis
  
  cluster_digihitposMod = allclusters; /// SJ
  
  //cout<<"allclusters.size() "<<allclusters.size()<<endl;

  if(allclusters.size()>1) cout<<"SJ!!! number of clusters in this event > 1 and is = "<<allclusters.size()<<endl;
  
  recopho.clear();

  for (int ij=0; ij<allclusters.size(); ij++) {
    phoparam tmppho;
    double dphi=0;
    double dphi2 = 0;
    double dthe=0;
    double dthe2=0;
    double ddiag=0;
    double etot=0;
    int ncell = allclusters[ij].size();
    double sdenr, sdthe, sdphi, sd1, sd9, sd25;
    for (int jk=0; jk<ncell; jk++) {
      double enr = allclusters[ij][jk].mag();
      double the = allclusters[ij][jk].theta();
      double phi = allclusters[ij][jk].phi();
      //      cout<<"jk "<<jk<<" "<<ncell<<endl;
      etot += enr;
      if (jk==0) {
	sdenr = sd1 = sd9 = sd25 = enr;
	sdthe = the;
	sdphi = phi;
      } else {
	double angle = allclusters[ij][0].angle(allclusters[ij][jk]);
	//	cout<<"df "<< allclusters[ij][0]<<" "<<allclusters[ij][jk]<<" "<<angle*radtodeg<<endl;
	if (angle*radtodeg < 2.9) { // larger than 2sqrt(2), but less than 3
	  sd25 +=enr;
	  if (angle*radtodeg < 1.5) { //more than sqrt(2) but less than 2
	    sd9 +=enr;
	  }
	}
	
	dthe += enr*(the - sdthe);
	dthe2 +=enr*(the - sdthe)*(the - sdthe);
	dphi +=enr*(phi - sdphi);
	dphi2 +=enr*(phi - sdphi)*(phi - sdphi);
	ddiag +=enr*(the - sdthe)*(phi - sdphi);
      }
    }
    if (etot<seedthres) continue;
    
    dphi /=etot;
    dthe /=etot;
    dthe2 /=etot;
    dphi2 /=etot;
    ddiag /=etot;

    // diagonised 2x2 matrix in (the,phi) plane
    double b2m4ac = sqrt((dthe2 - dphi2) * (dthe2 - dphi2) + 4* ddiag* ddiag);
    double val1 = ((dthe2 + dphi2) + b2m4ac)/2.;
    double val2 = ((dthe2 + dphi2) - b2m4ac)/2.;
    
    dthe2 = dthe2 - dthe*dthe; ///variance like quantity and hence telling the spread
    dphi2  = dphi2 -dphi*dphi;
    
    Hep3Vector tmp3v1; //central one
    tmp3v1.setRThetaPhi(sd1, sdthe, sdphi); ///direction and energy of the seed crystal
    Hep3Vector tmp3v2; //All others
    tmp3v2.setRThetaPhi(etot-sd1, sdthe+dthe, sdphi+dphi); //remaings of these

    tmppho.ppho = tmp3v1+tmp3v2;
    tmppho.nhits = ncell;
    tmppho.dthe2 = dthe2;
    tmppho.dphi2 = dphi2;
    tmppho.diag1 = val1;
    tmppho.diag2 = val2;
    tmppho.enr1 = sd1;
    tmppho.enr9 = sd9;
    tmppho.enr25 = sd25;
    //    cout<<"er1 "<<sd1<<" "<<sd9<<" "<<sd25<<" "<<etot<<endl;
    recopho.push_back(tmppho);
  }

  //Arranging them according to energy
  
  for (int ij=0; ij <recopho.size(); ij++) {
    for (int jk=ij+1; jk <recopho.size(); jk++) {
      if (recopho[ij].ppho.mag() < recopho[jk].ppho.mag()) {
	//	recopho.swap(recopho.begin()+ij, recopho.begin()+jk);
	swap(recopho[ij], recopho[jk]);
      }
    }
  }

  for (int ij=0; ij<recopho.size(); ij++) {
    cout<<"reco photon "<< recopho[ij].ppho<<" "<<recopho[ij].ppho.mag()<<" "<<radtodeg*recopho[ij].ppho.theta()<<" "<<radtodeg*recopho[ij].ppho.phi()<<endl;
  }


}


bool find_number(int nn, vector<int>array) {
  for (int ij=0; ij<array.size(); ij++) {
    if (nn==array[ij]) return true;
  }
  return false;
}


/// SJ
/// expectation maximization method
 void EMalgo(vector<Hep3Vector> digienr, vector<emalgoparam> &seeds_emalgo, vector<vector<double>> &hitsFrac_clus){

   cout<<"=============Inside EMalgo================="<<endl;
  const int niter = 30;
  double sigma_theta = 1; // ECAL crystals are also about a degree in size
  double sigma_phi = 1;
  //double seedthres = 1.;
  
  //amplitude of each particle
  const int NPART = 200;
  //const int NCELLS = 20000;
  //vector<double> amp(NPART,0.5); //assign 50%
  vector<double> amp(NPART,1); //assign 50%
  //double frac[NCELLS][NPART], contrib[NCELLS][NPART];
  double mu_theta[NPART], mu_phi[NPART];
  double mu_energy[NPART];

  vector<vector<double>> frac;
  vector<vector<double>> contrib;
  
  int nsize = digienr.size();
  vector<Hep3Vector> digienr_org = digienr;

  double meanmm=0;
  double meancc=0;
  double meanmm2=0;
  double meancc2=0;
  
  //A simple logic, but not an optimised one in the context of CPU
  // Use map/boundary logic to reduce the loop

  vector < Hep3Vector > seeds;


  /*cout<<"###############################"<<endl;
  cout<<"SJ!!! Inside EMalgo"<<endl;
  */
  
  for (int ij=0; ij<digienr.size(); ij++) {

    /*cout<<"Inside checking which hits pass seed threshold"<<endl;
    cout<<"This seed theta : phi: energy "<<digienr[ij].theta()*radtodeg<<" "<<digienr[ij].phi()*radtodeg<<" "<<digienr[ij].mag()<<endl;
    */
    
    //seeds.clear();
    if (digienr[ij].mag()<seedthres) continue; //already sorted collection

    
    seeds.push_back(digienr[ij]);
    
    //digienr.erase(digienr.begin()+ij);
  }

  /// check if the seeds are not neighbouring. If neighbouring, then take the ones with the highest energy

  cout<<"SJ!!! EMalgo : number of initial seeds before removal of neighbours "<<seeds.size()<<endl;
  
  for(int iclus=0; iclus<seeds.size(); iclus++){
    //cout<<"iclus : theta : phi : energy "<<iclus<<" "<<seeds[iclus].theta()*radtodeg<<" "<<seeds[iclus].phi()*radtodeg<<" "<<seeds[iclus].mag()<<endl;
    
    for(int jclus=iclus+1; jclus<seeds.size(); jclus++){
      //cout<<"jclus : theta : phi : energy "<<jclus<<" "<<seeds[jclus].theta()*radtodeg<<" "<<seeds[jclus].phi()*radtodeg<<" "<<seeds[jclus].mag()<<endl;
      
      double angle = seeds[iclus].angle(seeds[jclus].unit());
      //cout<<"angle is "<<angle*radtodeg<<endl;

      bool ishighEnDigiAround = isHighEnergyHitAround(digienr,seeds[jclus]);
      //if( angle*radtodeg<1.5 && ishighEnDigiAround ) ///sqrt(2) --> neighboring, touching the sides. Extent of ECAL xtal in degree is 1 degree
      if( angle*radtodeg<1.5 || ishighEnDigiAround ) ///sqrt(2) --> neighboring, touching the sides. Extent of ECAL xtal in degree is 1 degree

	{
	  //cout<<"seed with angle = "<<angle<< "and jclus = "<<jclus<<" removed"<<endl; 
	  seeds.erase(seeds.begin()+jclus); ///sorted in energy, so remove the next seeds
	  jclus--;
	}
    }

  }//for(int iclus=0; iclus<seeds.size(); iclus++)

  //cluster = seeds;
    
  int npart = seeds.size();
  int ncells = digienr.size();

  cout<<"SJ!!! #seeds after the nearest neighbour algo "<<npart<<endl;


   for(int ipart=0; ipart<npart; ipart++){
     amp[ipart] = 0.5;
     mu_theta[ipart] = seeds[ipart].theta()*radtodeg;
     mu_phi[ipart] = seeds[ipart].phi()*radtodeg;
     mu_energy[ipart] = seeds[ipart].mag();
   }

   ////before the EM algo, form a set of hits which are neibours with each others - as is done in the traditional algorithm given in ecal_clustering

   vector < vector < Hep3Vector > > allclusters;
   vector < Hep3Vector > cluster;
    

   
   for (int iseed=0; iseed<seeds.size(); iseed++) {
     vector < Hep3Vector > cluster;
     cluster.push_back(seeds[iseed]);
     
     vector< Hep3Vector > digienr_tmp = digienr; //so that I can still keep using digienr in the code later and always start from the original collection of seed for each crystal

     //cout<<"SJ!!! making simple cluster first for each seed, seed info: en : theta : phi : "<<seeds[iseed].mag()<<" "<<seeds[iseed].theta()<<" "<<seeds[iseed].phi()<<endl;
     
     for(int jk=0; jk<cluster.size(); jk++){
       //cout<<"SJ!!! making simple cluster first for each seed, cluster size of seed number "<<iseed<<" for jk="<<jk<<" is "<<cluster.size()<<endl;
       for (int kl=0; kl<digienr_tmp.size(); kl++) {
	 if (digienr_tmp[kl].mag()<otherthrs) continue;
	 double angle = cluster[jk].angle(digienr_tmp[kl].unit());

	 //cout<<"SJ!!! For digi kl="<<kl<<" and en : theta : phi : "<<digienr_tmp[kl].mag()<<" "<<digienr_tmp[kl].theta()<<" "<<digienr_tmp[kl].phi()<< " angle = "<<angle<<endl;
	 if (angle*radtodeg <1.5) { /// take only neiboring cells to what are already in the cluster
	   //cout<<"SJ!!! This digi included"<<endl;
	   cluster.push_back(digienr_tmp[kl]);
	   digienr_tmp.erase(digienr_tmp.begin()+kl); 
	   kl--;
	   
	 }
       }//for (int kl=0; kl<digienr.size(); kl++)
     }//for (int jk=0; jk<cluster.size(); jk++)

     allclusters.push_back(cluster);
   }//for (int iseed=0; iseed<seeds.size(); iseed++)

   

   
   //real EM algo
   for(int iiter=0; iiter<niter; iiter++){
     /*cout<<"=========================="<<endl;
     cout<<"SJ!!! Iter "<<iiter<<endl;
     */
     
     frac.clear();
     contrib.clear();
     
    //expectation
     for(int icell=0; icell<ncells; icell++){
      
       //cout<<""<<endl;
       vector<double> contrib_tmp_gaus;
       vector<double> frac_tmp_gaus;
       
       double en = digienr[icell].mag();

      //if(en < otherthrs)  continue;
      
      double cell_theta = digienr[icell].theta()*radtodeg;
      double cell_phi   = digienr[icell].phi()*radtodeg;

      //cout<<"cell energy : theta : phi : "<<en<<" "<<cell_theta<<" "<<cell_phi<<endl;
      double sum = 0;    
      for(int ipart=0; ipart<npart; ipart++){

	/// check if the cell under consideration is a part of the cluster for that seed
	bool isdigiincluster = isDigiInAllClusters(digienr[icell], allclusters[ipart]);
	
	//cout<<"ipart : en : mu_theta : mu_phi "<<ipart<<" "<<mu_energy[ipart]<<" "<<mu_theta[ipart]<<" "<<mu_phi[ipart]<<endl;
	double arg_the = (cell_theta - mu_theta[ipart])/sigma_theta;
	double arg_phi = (cell_phi - mu_phi[ipart])/sigma_phi;


	//contrib[icell][ipart] = amp[ipart]* TMath::Exp(-0.5*arg_the*arg_the - 0.5*arg_phi*arg_phi);
	
	double tmpcont = amp[ipart]* TMath::Exp(-0.5*arg_the*arg_the - 0.5*arg_phi*arg_phi);
	
	//cout<<"arg_the : arg_phi : amp[ipart] : contrib "<<arg_the<<" "<<arg_phi<<" "<<amp[ipart]<<" "<<contrib[icell][ipart]<<endl;
	//cout<<"arg_the : arg_phi : amp[ipart] : contrib "<<arg_the<<" "<<arg_phi<<" "<<amp[ipart]<<" "<<tmpcont<<endl;
	
	if(en < otherthrs)
	  //contrib[icell][ipart] = 0;
	  tmpcont = 0;
	
	if(!isdigiincluster)   
	  //contrib[icell][ipart] = 0;
	  tmpcont = 0;
	
	//cout<<"ipart : en : mu_theta : mu_phi contrib["<<icell<<"]["<<ipart<<"] "<<ipart<<" "<<mu_energy[ipart]<<" "<<mu_theta[ipart]<<" "<<mu_phi[ipart]<<" "<<contrib[icell][ipart]<<endl;
	//cout<<"ipart : en : mu_theta : mu_phi contrib["<<icell<<"]["<<ipart<<"] "<<ipart<<" "<<mu_energy[ipart]<<" "<<mu_theta[ipart]<<" "<<mu_phi[ipart]<<" "<<tmpcont<<endl;
	
	//sum += contrib[icell][ipart];
	sum += tmpcont;
	contrib_tmp_gaus.push_back(tmpcont);
      }

      contrib.push_back(contrib_tmp_gaus);
     
      for(int ipart=0; ipart<npart; ipart++){
	
	//if(sum!=0) frac[icell][ipart] = contrib[icell][ipart] /sum;
	double frac_tmp = 0;
	if(sum!=0) frac_tmp = contrib[icell][ipart] /sum;
	else frac_tmp = 0;
	frac_tmp_gaus.push_back(frac_tmp);
	
	//cout<<"After normalization, Sum : frac["<<icell<<"]["<<ipart<<"] : "<<sum<<" "<<frac_tmp<<endl;
      }
      
      frac.push_back(frac_tmp_gaus);
    }//for(int icell=0; icell<ncells; icell++)
    
    
    //maximizationn 
    for(int ipart=0; ipart<npart; ipart++){
      
      double a_tmp = 0, mu_tmp_theta = 0, mu_tmp_phi = 0;
      
      for(int icell=0; icell<ncells; icell++){
	
	/// check if the cell under consideration is a part of the cluster for that seed
	/// in the loop above, I already set the frac of such a cell to be 0
	//bool isdigiincluster = isDigiInAllClusters(digienr[icell], allclusters[ipart]);
	
	//if(!isdigiincluster) continue;

	double en = digienr[icell].mag();
	
	///no need for this since frac is already 0 for low energy particles
	//if(en < otherthrs) continue;
	
	double theta_cell = digienr[icell].theta()*radtodeg;
	double phi_cell = digienr[icell].phi()*radtodeg;
	
	a_tmp += frac[icell][ipart] * en;  ///dont need to divide by sum of f's because norm of part is the sum of contribution from each cell
	
	mu_tmp_theta += frac[icell][ipart] * en * theta_cell; ///need to divide by sum of weights since this is the average of all the values

	mu_tmp_phi += frac[icell][ipart] * en * phi_cell; ///need to divide by sum of weights since this is the average of all the values
	
	
	//cout<<"ctheta : cphi : f : en : a : mu "<<theta_cell<< " "<<phi_cell<<" "<<frac[icell][ipart]<<" "<< en <<" "<<a_tmp<<" "<<mu_tmp_theta<<endl;
      }//for(int icell=0; icell<ncells; icell++)

      //a_tmp = a_tmp;
      mu_tmp_theta = mu_tmp_theta/a_tmp;
      mu_tmp_phi = mu_tmp_phi/a_tmp;
      
      //cout<<"Max step number "<<iiter<<" with a_tmp : mu_tmp_theta : mu_tmp_phi "<<a_tmp<<" "<<mu_tmp_theta<<" "<<mu_tmp_phi<<endl;
      
      amp[ipart] = a_tmp;
      mu_theta[ipart] = mu_tmp_theta;
      mu_phi[ipart] = mu_tmp_phi;

      //f_est[ipart]->SetParameters(a[ipart],mu[ipart],sigma);

      //f_est[ipart]->Draw("same");
      //c1->Update();
    }//for(int ipart=0; ipart<npart; ipart++)

    //c1->Print(Form("plot_iter%d.gif",iiter));
    //printf("\nDouble click in the bottom right corner of the pad to continue\n");
    //c1->WaitPrimitive();
    
    
  }//for(int iiter=0; iiter<niter; iiter++)

  emalgoparam seed_tmp;
  ///store the final information of all the seedss here
  for(int ipart=0; ipart<npart; ipart++){
    
    seed_tmp.mu_theta = mu_theta[ipart];
    seed_tmp.mu_phi = mu_phi[ipart];
    seed_tmp.amp = amp[ipart];

    seeds_emalgo.push_back(seed_tmp);
  }
  
  ///also return the fraction of each cell associated with that seed which is stored in frac.

  for(int ipart=0; ipart<npart; ipart++){
    vector<double> hitsFrac; 
    for(int icell=0; icell<ncells; icell++){

      //cout<<"SJ!!! final fracs for igaus : icell and frac : "<<ipart<<" "<<icell<<" "<<frac[icell][ipart]<<endl;
      hitsFrac.push_back(frac[icell][ipart]);
    }
    hitsFrac_clus.push_back(hitsFrac);
  }
  
}



bool isDigiInAllClusters(Hep3Vector digi,  std::vector<Hep3Vector> cluster) {
    for (const auto& digiInClus : cluster) {
        if (digiInClus == digi) {
	  return true;  // Found a match
        }
    }
    return false;  // Digi is not in any cluster
}




bool isHighEnergyHitAround(vector<Hep3Vector> digienr, Hep3Vector seed){

  /*cout<<"Inside isHighEnergyHitAround "<<endl;
  cout<<"Seed energy is "<<seed.mag()<<endl;
  */
  
  bool rejectSeed = false;
    for (int ij=0; ij<digienr.size(); ij++) {

      double angle = seed.angle(digienr[ij].unit());

      double en_digi = digienr[ij].mag();
      double en_seed = seed.mag();
      //cout<<"digi En : angle "<<en_digi<<" "<<angle*radtodeg<<endl;
      if(en_digi > en_seed && angle*radtodeg < 1.5){
	rejectSeed = true;
	break;
      }
    }

    //cout<<"rejectSeed "<<rejectSeed<<endl;    
    return rejectSeed;
}


/// ================ ecal_clustering but with separate seeds
void ecal_clustering_seed(vector<Hep3Vector> digienr, vector<vector<Hep3Vector>> &cluster_digihitposMod_seed) {

  //cout<<"SJ!!! inside ecal_clustering_seed "<<endl;
  
  int nsize = digienr.size();
  vector<Hep3Vector> digienr_org = digienr;

  double meanmm=0;
  double meancc=0;
  double meanmm2=0;
  double meancc2=0;
  
  //A simple logic, but not an optimised one in the context of CPU
  // Use map/boundary logic to reduce the loop

  vector < Hep3Vector > seeds;


  /*cout<<"###############################"<<endl;
  cout<<"SJ!!! Inside EMalgo"<<endl;
  */
  
  for (int ij=0; ij<digienr.size(); ij++) {

    if (digienr[ij].mag()<seedthres) continue; //already sorted collection
    seeds.push_back(digienr[ij]);
    
    //digienr.erase(digienr.begin()+ij);
  }

  /// check if the seeds are not neighbouring. If neighbouring, then take the ones with the highest energy

  for(int iclus=0; iclus<seeds.size(); iclus++){
    for(int jclus=iclus+1; jclus<seeds.size(); jclus++){
      double angle = seeds[iclus].angle(seeds[jclus].unit());
      bool ishighEnDigiAround = isHighEnergyHitAround(digienr,seeds[jclus]);
      if( angle*radtodeg<1.5 || ishighEnDigiAround ) ///sqrt(2) --> neighboring, touching the sides. Extent of ECAL xtal in degree is 1 degree

	{
	  seeds.erase(seeds.begin()+jclus); ///sorted in energy, so remove the next seeds
	  jclus--;
	}
    }

  }//for(int iclus=0; iclus<seeds.size(); iclus++)

  //cout<<"Number of final seeds "<<seeds.size()<<endl;
  
  //cluster = seeds;
    
  int npart = seeds.size();
  int ncells = digienr.size();

   ////before the EM algo, form a set of hits which are neibours with each others - as is done in the traditional algorithm given in ecal_clustering

   vector < vector < Hep3Vector > > allclusters;
   vector < Hep3Vector > cluster;
   
   for (int iseed=0; iseed<seeds.size(); iseed++) {
     vector < Hep3Vector > cluster;
     cluster.push_back(seeds[iseed]);
     
     //vector< Hep3Vector > digienr_tmp = digienr; //so that I can still keep using digienr in the code later and always start from the original collection of seed for each crystal

     //cout<<"Size of digi at this point "<<digienr.size()<<endl;
     
     for(int jk=0; jk<cluster.size(); jk++){
       //cout<<"SJ!!! making simple cluster first for each seed, cluster size of seed number "<<iseed<<" for jk="<<jk<<" is "<<cluster.size()<<endl;
       for (int kl=0; kl<digienr.size(); kl++) {
	 if (digienr[kl].mag()<otherthrs) continue;
	 double angle = cluster[jk].angle(digienr[kl].unit());
	 
	 //cout<<"SJ!!! For digi kl="<<kl<<" and en : theta : phi : "<<digienr[kl].mag()<<" "<<digienr[kl].theta()<<" "<<digienr[kl].phi()<< " angle = "<<angle<<endl;
	 if (angle*radtodeg <1.5 && digienr[kl].mag()<cluster[jk].mag()) { /// take only neiboring cells to what are already in the cluster
	   //cout<<"SJ!!! This digi included"<<endl;
	   cluster.push_back(digienr[kl]);
	   digienr.erase(digienr.begin()+kl); 
	   kl--;
	   
	 }
       }//for (int kl=0; kl<digienr.size(); kl++)
     }//for (int jk=0; jk<cluster.size(); jk++)

     allclusters.push_back(cluster);

   }//for (int iseed=0; iseed<seeds.size(); iseed++)

   cluster_digihitposMod_seed = allclusters;
  


}


