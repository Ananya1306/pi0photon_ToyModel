# pi0photon_ToyModel
This is a toy model on how to calculate pi0 photon yields by scaling and background subtraction menthods
void final(){

//n = 20 random numbers are generated in an interval ( 2.0, 1.0 ) for 10 6 times to replicate
//a million events detected by the calorimeter in an actual heavy ion collision experiment
//and stored in a 1D histogram h2

TH1D *h2 = new TH1D("h2","distance between two non-correlated random numbers", 3000,0,3);
TH1D *h1 = new TH1D("h1","distance between two correlated random numbers",3000,0,3);

int n =20;
int m = 1000000;
double x[n], y[n];
double dist=0.;



for(int j=0; j<m; j++){



        for(int i=0; i<n; i++){
        x[i] = gRandom->Uniform(0,2);
        y[i] = gRandom->Uniform(0,1);
	}

	for(int p=0; p<n-1; p++){
        for(int q=p+1; q<n; q++){
                dist = sqrt(pow(x[p]-x[q],2)+pow(y[p]-y[q],2));
                h2->Fill(dist);
		}
	}




//n = 20 random points are generated with a constraint for the n th point to lie at a 
//fixed distance (= 0.2, 0.5, 0.8, 1.1 ) with a gaussian distribution around the ( n − 1 ) th point
//to replicate the detection of π 0 particles from a common vertex and stored in a 1D histogram h1


	for(int i=0; i<n; i++){
	x[i] = gRandom->Uniform(0,2);
 	y[i] = gRandom->Uniform(0,1);
	
	if (i == n-1)
		{
		double a = gRandom->Uniform(-1.0,1.0);
		double phi = a*TMath::Pi();
		double R = 0.45 + (gRandom->Gaus(0.,0.01)) ;
		double delx = R*TMath::Cos(phi);
		double dely = R*TMath::Sin(phi);
		x[i]=x[i-1]+delx;
		y[i]=y[i-1]+dely;
		//cout<<"\t"<<delx<<"\t"<<dely<<"\t"<<sqrt(pow(delx,2)+pow(dely,2))<<endl;				
		}
	}

//Points which lie outside the interval ( 2.0, 1.0 ) are rejected

for(int p=0; p<n-1; p++){

if(x[n-1]<0.0 || x[n-1]>2.0) break;
  if(y[n-1]<0.0 || y[n-1]>1.0) break;
  
	for(int q=p+1; q<n; q++){
		dist = sqrt(pow(x[p]-x[q],2)+pow(y[p]-y[q],2));
		h1->Fill(dist);
		}
	}


}

//h1 is now integrated around the Gaussian peak (without including the values from − 3σ to 3σ of the
//distribution). h2 is also integrated around the same values. h1 is then scaled up to h2 by a factor obtained
//after dividing the two integrated values and superposed over h2.


double a = h1->Integral(0,170);
double b = h1->Integral(230,300);
double c = h2->Integral(0,170);
double d = h2->Integral(230,300);
double e = (c+d)/(a+b);
h1->Scale(e);
h1->Draw();
h2->SetLineColor(kRed);
h2->Draw("same");



}
