#include<iostream>
#include<fstream>
using namespace std;


void yield_error_matrix(int r0, int s0){
	//r0 is the ring number: 1,2,3,4,5,6
	//s0 is the sector number: 0=closed, 1=transition, 2=open
	
	const int dimension_matrix=3;   // A_ee, A_epelastic, A_epinelastic
	
	TFile *file=new TFile("Counts_vs_asy.root");
	TTree *T=file->Get("T");
	
	int ring_number, sector_number;
	double rate_ee, rate_epelastic, rate_epinelastic;
	double A_ee, A_epelastic, A_epinelastic;

	T->SetBranchAddress("ring_number",&ring_number);
	T->SetBranchAddress("sector_number",&sector_number);
	T->SetBranchAddress("rate_ee",&rate_ee);
	T->SetBranchAddress("rate_epelastic",&rate_epelastic);
	T->SetBranchAddress("rate_epinelastic",&rate_epinelastic);
	T->SetBranchAddress("A_ee",&A_ee);
	T->SetBranchAddress("A_epelastic",&A_epelastic);
	T->SetBranchAddress("A_epinelastic",&A_epinelastic);

	int entries=T->GetEntries();

	//to get the bin in which people want to extract error matrix
	double A0_ee, A0_epelastic, A0_epinelastic;
	double rate0_ee, rate0_epelastic, rate0_epinelastic;
	for(int i=0;i<entries;i++){
		T->GetEntry(i);
		if(ring_number==r0&&sector_number==s0){
			cout<<"find the bin to extract error matrix...moving on..."<<endl;
			///////////////////////////////////////
			//ee epelastic epinelastic
			//////////////////////////////////////
			rate0_ee=rate_ee;
			rate0_epelastic=rate_epelastic;
			rate0_epinelastic=rate_epinelastic;
			A0_ee=A_ee;
			A0_epelastic=A_epelastic;
			A0_epinelastic=A_epinelastic;
		}
		
	}
	

	// now loop tree again to calculate F matrix
	TMatrixD F(dimension_matrix, dimension_matrix);
	//initialization...
	for(int ii=0;ii<dimension_matrix;ii++){
		for(int jj=0;jj<dimension_matrix;jj++){
			F(ii,jj)=0;
		}
	}
	
	
	for(int j=0;j<entries;j++){
		T->GetEntry(j);
		//////////////////////////////////////////////////////
		//ee epelasitc epinelastic
		//////////////////////////////////////////////////////
		double rate_tot=rate_ee + rate_epelastic + rate_epinelastic;
		double sigma_Am=1.0/sqrt(rate_tot);
		double f_ee=rate_ee/rate_tot *A_ee/A0_ee;
		double f_epelastic=rate_epelastic/rate_tot * A_epelastic/A0_epelastic;
		double f_epinelastic=rate_epinelastic/rate_tot * A_epinelastic/A0_epinelastic;
		
		// get the info for one entry, then sum it up for the matrix F
		double f[dimension_matrix]={f_ee, f_epelastic, f_epinelastic};
		for(int k1=0;k1<dimension_matrix;k1++){
			for(int k2=0;k2<dimension_matrix;k2++){
				F(k1,k2) += f[k1]*f[k2]/sigma_Am/sigma_Am;
			}
		}
	}

	F.Invert();

	//////////////////////////////////////////////
	//output matrix (235+95+14) days, pol=0.8
	//////////////////////////////////////////////
//	cout<<" ee:          "<<sqrt(F(0,0))<<endl;
//	cout<<" epelastic :  "<<sqrt(F(1,1))<<endl;
//	cout<<" epinelastic: "<<sqrt(F(2,2))<<endl;


	double sigma_ee=sqrt(F(0,0))/0.8/sqrt((235+95+14)*24.0*60*60);
	double sigma_epelastic=sqrt(F(1,1))/0.8/sqrt((235+95+14)*24.0*60*60);
	double sigma_epinelastic=sqrt(F(2,2))/0.8/sqrt((235+95+14)*24.0*60*60);

	//overall measurement containing the information from all bins:
	cout<<"--------------- overall ------------------"<<endl;
	cout<<" ee:          "<<A0_ee<<"	"<<sigma_ee/1e-9<<"	"<<sigma_ee/A0_ee/1e-9<<endl;
	cout<<" epelastic :  "<<A0_epelastic<<"	"<<sigma_epelastic/1e-9<<"	"<<sigma_epelastic/A0_epelastic/1e-9<<endl;
	cout<<" epinelastic: "<<A0_epinelastic<<"	"<<sigma_epinelastic/1e-9<<"	"<<sigma_epinelastic/A0_epinelastic/1e-9<<endl;
	
	//in the selected bin
	cout<<"============== in the selected bin ================="<<endl;
	double sys_epelastic=rate0_epelastic/rate0_ee * sigma_epelastic;
	double sys_epinelastic=rate0_epinelastic/rate0_ee * sigma_epinelastic;
	double rate0_tot=rate0_ee + rate0_epelastic + rate0_epinelastic;
	
	double stat_ee=rate0_tot/rate0_ee * 1.0/sqrt(rate0_tot)/0.8/sqrt((235+95+14)*24.0*60*60);
	
	cout<<" stat ee: "<<stat_ee/1e-9<<"	"<<stat_ee/1e-9/33<<endl;
	cout<<" sys_epelastic: "<<sys_epelastic/1e-9<<"	"<<sys_epelastic/33/1e-9<<endl;
	cout<<" sys_epinelastic: "<<sys_epinelastic/1e-9<<"	"<<sys_epinelastic/33/1e-9<<endl;

	//cout<<rate0_epelastic/rate0_ee<<endl;
	//cout<<rate0_epinelastic/rate0_ee<<endl;

	
}



		

