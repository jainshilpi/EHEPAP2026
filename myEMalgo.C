
const int ncells = 60;
double xmax = 30;

void myEMalgo(){

  
  const int ntrials = 300;
  const int niter = 30;
  const int ngaus = 3;
  


  double sigma = 1.;
  
  //TF1 *f = new TF1("f","gaus",0,ncells);
  TF1 *f1 = new TF1("f1",myGaus,0,ncells,3);
  f1->SetParameter(0,1);
  f1->SetParameter(1,5);
  f1->SetParameter(2, sigma);

  TF1 *f2 = new TF1("f2",myGaus,0,ncells,3);
  f2->SetParameter(0,1);
  f2->SetParameter(1,8);
  f2->SetParameter(2, sigma);


  TF1 *f3 = new TF1("f3",myGaus,0,ncells,3);
  f3->SetParameter(0,1);
  f3->SetParameter(1,15);
  f3->SetParameter(2, sigma);
  

  TF1 *f_est[ngaus];
  int icol[] = {1,2,3,4,5,6,7,8,9,10};
  for(int igaus=0; igaus<ngaus; igaus++){
    f_est[igaus] = new TF1(Form("f%d_est",igaus),myGausWithBW,0,ncells,3);
    f_est[igaus]->SetLineColor(icol[igaus]);
  }


  //TH1F * h = new TH1F("h","",ncells,0,ncells);
  TH1F * h = new TH1F("h","",ncells,0,xmax);

  for(int itrial=0; itrial<ntrials; itrial++){
    
    
    h->Fill(f1->GetRandom());
    
    h->Fill(f2->GetRandom());
    h->Fill(f2->GetRandom());


    h->Fill(f3->GetRandom());
    h->Fill(f3->GetRandom());
    h->Fill(f3->GetRandom());
    
    
  }


  double c[ncells];
  double en[ncells];
  
  for(int ibin=1; ibin<=ncells; ibin++){

    c[ibin-1] = h->GetBinCenter(ibin);
    en[ibin-1] = h->GetBinContent(ibin);

    //cout<<"ibin : binc : ene : "<<ibin<<" "<<c[ibin-1]<<" "<<en[ibin-1]<<endl;
  }

  //double binw = h->GetBinWidth(1);
  //cout<<"binw "<<binw<<endl;
  //h->Scale(binw);
  

  h->Draw();
  

  double f[100][100], a[100], mu[100];

  for(int igaus=0; igaus<ngaus; igaus++){

    ///it seems that assigning diff values to a and mu is imp else the result is same as the beginning
    a[igaus] = 1*(igaus+10);
    mu[igaus] = 1+igaus;
  }

  for(int iiter=0; iiter<niter; iiter++){
    
    //expectation
    for(int icell=0; icell<ncells; icell++){
    
      double sum = 0;    
      for(int igaus=0; igaus<ngaus; igaus++){
	
	//cout<<"igaus : mu : a "<<igaus<<" "<<mu[igaus]<<" "<<a[igaus]<<endl;
	double arg = (c[icell] - mu[igaus])/sigma;

	f[icell][igaus] = a[igaus]* TMath::Exp(-0.5*arg*arg);
	
	sum += f[icell][igaus];
      }
      
      for(int igaus=0; igaus<ngaus; igaus++){
	f[icell][igaus] = f[icell][igaus] /sum;
	//cout<<"f["<<icell<<"]["<<igaus<<"] "<<f[icell][igaus]<<endl;
      }
      
      
    }//for(int icell=0; icell<ncells; icell++)
  

    //maximizationn 
    for(int igaus=0; igaus<ngaus; igaus++){
      
      double a_tmp = 0, mu_tmp = 0;

      for(int icell=0; icell<ncells; icell++){

	
	a_tmp += f[icell][igaus] * en[icell];  ///dont need to divide by sum of f's because norm of gaus is the sum of contribution from each cell
	
	mu_tmp += f[icell][igaus] * en[icell] * c[icell]; ///need to divide by sum of weights since this is the average of all the values
	//cout<<"c : f : en : a : mu "<<c[icell]<< " "<<f[icell][igaus]<<" "<< en[icell]<<" "<<a_tmp<<" "<<mu_tmp<<endl;
      }//for(int icell=0; icell<ncells; icell++)

      //a_tmp = a_tmp;
      mu_tmp = mu_tmp/a_tmp;
      
      cout<<"a_tmp : mu_tmp "<<a_tmp<<" "<<mu_tmp<<endl;
      
      a[igaus] = a_tmp;
      mu[igaus] = mu_tmp;

      f_est[igaus]->SetParameters(a[igaus],mu[igaus],sigma);

      f_est[igaus]->Draw("same");
      c1->Update();
    }//for(int igaus=0; igaus<ngaus; igaus++)

    c1->Print(Form("plot_iter%d.gif",iiter));
    printf("\nDouble click in the bottom right corner of the pad to continue\n");
    c1->WaitPrimitive();
    

  }//for(int iiter=0; iiter<niter; iiter++)

  ///now try to estimate the Gaussian parameters with expectation maximization algorithm
  

}




Double_t myGaus(double *x, double *par){

  double root2pi = sqrt(2.*TMath::Pi());
  Double_t arg = 0;
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  
  return fitval/root2pi;

}



Double_t myGausWithBW(double *x, double *par){

  double binwidth = xmax/ncells;
  double root2pi = sqrt(2.*TMath::Pi());
  Double_t arg = 0;
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  
  return fitval*binwidth/root2pi;

}
