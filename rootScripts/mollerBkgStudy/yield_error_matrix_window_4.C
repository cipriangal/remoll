#include<iostream>
#include<fstream>
using namespace std;


void yield_error_matrix_window_4(int r0, int s0){
	//r0 is the ring number: 1,2,3,4,5,6
	//s0 is the sector number: 0=closed, 1=transition, 2=open
	
	const int dimension_matrix=4;   // A_ee, A_epelastic, A_epinelastic, A_eAlelastic
	
	TFile *file=new TFile("Counts_vs_asy.root");
	TTree *T=file->Get("T");
	
	int ring_number, sector_number;
	double rate_ee, rate_epelastic, rate_epinelastic, rate_eAlelastic, rate_eAlquasielastic, rate_eAlinelastic, rate_pim;
	double A_ee, A_epelastic, A_epinelastic, A_eAlelastic, A_eAlquasielastic, A_eAlinelastic, A_pim;

	T->SetBranchAddress("ring_number",&ring_number);
	T->SetBranchAddress("sector_number",&sector_number);
	T->SetBranchAddress("rate_ee",&rate_ee);
	T->SetBranchAddress("rate_epelastic",&rate_epelastic);
	T->SetBranchAddress("rate_epinelastic",&rate_epinelastic);
	T->SetBranchAddress("rate_eAlelastic",&rate_eAlelastic);
	T->SetBranchAddress("rate_eAlquasielastic",&rate_eAlquasielastic);
	T->SetBranchAddress("rate_eAlinelastic",&rate_eAlinelastic);
	T->SetBranchAddress("rate_pim",&rate_pim);
	
	T->SetBranchAddress("A_ee",&A_ee);
	T->SetBranchAddress("A_epelastic",&A_epelastic);
	T->SetBranchAddress("A_epinelastic",&A_epinelastic);
	T->SetBranchAddress("A_eAlelastic",&A_eAlelastic);
	T->SetBranchAddress("A_eAlquasielastic",&A_eAlquasielastic);
	T->SetBranchAddress("A_eAlinelastic",&A_eAlinelastic);
	T->SetBranchAddress("A_pim",&A_pim);
	
	int entries=T->GetEntries();

	//to get the bin in which people want to extract error matrix
	double A0_ee, A0_epelastic, A0_epinelastic, A0_eAlelastic, A0_eAlquasielastic, A0_eAlinelastic, A0_pim;
	double rate0_ee, rate0_epelastic, rate0_epinelastic, rate0_eAlelastic, rate0_eAlquasielastic, rate0_eAlinelastic, rate0_pim;
	for(int i=0;i<entries;i++){
		T->GetEntry(i);
		if(ring_number==r0&&sector_number==s0){
			cout<<"find the bin to extract error matrix...moving on..."<<endl;
			///////////////////////////////////////
			//ee epelastic epinelastic, eAlelastic, eAlquasielastic, eAlinelastic
			//////////////////////////////////////
			rate0_ee=rate_ee;
			rate0_epelastic=rate_epelastic;
			rate0_epinelastic=rate_epinelastic;
			rate0_eAlelastic=rate_eAlelastic;
			rate0_eAlquasielastic=rate_eAlquasielastic;
			rate0_eAlinelastic=rate_eAlinelastic;
			rate0_pim=rate_pim;
			
			A0_ee=A_ee;
			A0_epelastic=A_epelastic;
			A0_epinelastic=A_epinelastic;
			A0_eAlelastic=A_eAlelastic;
			A0_eAlquasielastic=A_eAlquasielastic;
			A0_eAlinelastic=A_eAlinelastic;
			A0_pim=A_pim;
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
	
	double B[dimension_matrix]={0};
	for(int j=0;j<entries;j++){
		T->GetEntry(j);
		//////////////////////////////////////////////////////
		//ee epelasitc epinelastic, eAlelastic, eAlquasielastic, eAlinelastic
		//////////////////////////////////////////////////////
		double rate_tot=rate_ee + rate_epelastic + rate_epinelastic + rate_eAlelastic;
		double sigma_Am=1.0/sqrt(rate_tot);
		double f_ee=rate_ee/rate_tot *A_ee/A0_ee;
		double f_epelastic=rate_epelastic/rate_tot * A_epelastic/A0_epelastic;
		double f_epinelastic=rate_epinelastic/rate_tot * A_epinelastic/A0_epinelastic;
		
		double f_eAlelastic=rate_eAlelastic/rate_tot * A_eAlelastic/A0_eAlelastic;
		
		// for unbiased check
		double A_m_i=f_ee*A0_ee + f_epelastic*A0_epelastic + f_epinelastic*A0_epinelastic + f_eAlelastic*A0_eAlelastic;

		// get the info for one entry, then sum it up for the matrix F
		double f[dimension_matrix]={f_ee, f_epelastic, f_epinelastic, f_eAlelastic};
		for(int k1=0;k1<dimension_matrix;k1++){
			for(int k2=0;k2<dimension_matrix;k2++){
				F(k1,k2) += f[k1]*f[k2]/sigma_Am/sigma_Am;
			}
		}
		for(int k3=0;k3<dimension_matrix;k3++){
			B[k3] += A_m_i*f[k3]/sigma_Am/sigma_Am;
		}
	} //loop over 18 measurements

	F.Invert();

	cout<<"__________ asymmetry estimator__________________________________________________"<<endl;
	cout<<" A_ee : "<<F(0,0)*B[0]+F(0,1)*B[1]+F(0,2)*B[2]+F(0,3)*B[3]<<endl;
	cout<<" A_epelastic: "<<F(1,0)*B[0]+F(1,1)*B[1]+F(1,2)*B[2]+F(1,3)*B[3]<<endl;
	cout<<" A_epinelastic: "<<F(2,0)*B[0]+F(2,1)*B[1]+F(2,2)*B[2]+F(2,3)*B[3]<<endl;
	cout<<" A_eAlelastic: "<<F(3,0)*B[0]+F(3,1)*B[1]+F(3,2)*B[2]+F(3,3)*B[3]<<endl;
	cout<<endl;
	cout<<endl;
	cout<<endl;
	cout<<endl;

	//////////////////////////////////////////////
	//output matrix (235+95+14) days, pol=0.8
	//////////////////////////////////////////////
//	cout<<" ee:          "<<sqrt(F(0,0))<<endl;
//	cout<<" epelastic :  "<<sqrt(F(1,1))<<endl;
//	cout<<" epinelastic: "<<sqrt(F(2,2))<<endl;


	double sigma_ee=sqrt(F(0,0))/0.8/sqrt((235+95+14)*24.0*60*60);
	double sigma_epelastic=sqrt(F(1,1))/0.8/sqrt((235+95+14)*24.0*60*60);
	double sigma_epinelastic=sqrt(F(2,2))/0.8/sqrt((235+95+14)*24.0*60*60);

	double sigma_eAlelastic=sqrt(F(3,3))/0.8/sqrt((235+95+14)*24.0*60*60);
	cout.precision(10);	
	//overall measurement containing the information from all bins:
	cout<<"--------------- overall ------------------"<<endl;
	cout<<" ee:          "<<A0_ee<<"	"<<sigma_ee/1e-9<<"	"<<sigma_ee/A0_ee/1e-9<<endl;
	cout<<" epelastic :  "<<A0_epelastic<<"	"<<sigma_epelastic/1e-9<<"	"<<sigma_epelastic/A0_epelastic/1e-9<<endl;
	cout<<" epinelastic: "<<A0_epinelastic<<"	"<<sigma_epinelastic/1e-9<<"	"<<sigma_epinelastic/A0_epinelastic/1e-9<<endl;
	
	cout<<" eAlelastic :  "<<A0_eAlelastic<<"	"<<sigma_eAlelastic/1e-9<<"	"<<sigma_eAlelastic/A0_eAlelastic/1e-9<<endl;
	
	//in the selected bin
	cout<<"============== in the selected bin ================="<<endl;
	double sys_epelastic=rate0_epelastic/rate0_ee * sigma_epelastic;
	double sys_epinelastic=rate0_epinelastic/rate0_ee * sigma_epinelastic;
	
	double sys_eAlelastic=rate0_eAlelastic/rate0_ee * sigma_eAlelastic;

	//subtracting terms...
	double sys_eAlquasielastic=rate0_eAlquasielastic/rate0_ee*fabs(A0_eAlquasielastic);
	double sys_eAlinelastic=rate0_eAlinelastic/rate0_ee*fabs(A0_eAlinelastic);
	double sys_pim=rate0_pim/rate0_ee*fabs(A0_pim)*0.1;  //10% uncertainty

	double rate0_tot=rate0_ee + rate0_epelastic + rate0_epinelastic + rate0_eAlelastic+ rate0_eAlquasielastic + rate0_eAlinelastic + rate0_pim;
	double stat_ee=rate0_tot/rate0_ee * 1.0/sqrt(rate0_tot)/0.8/sqrt((235+95+14)*24.0*60*60);
	
	cout<<" stat ee: "<<stat_ee/1e-9<<"	"<<stat_ee/1e-9/33<<endl;
	cout<<" sys_epelastic: "<<sys_epelastic/1e-9<<"	"<<sys_epelastic/33/1e-9<<endl;
	cout<<" sys_epinelastic: "<<sys_epinelastic/1e-9<<"	"<<sys_epinelastic/33/1e-9<<endl;
	
	cout<<" sys_eAlelastic: "<<sys_eAlelastic/1e-9<<"	"<<sys_eAlelastic/33/1e-9<<endl;
	
	//subtracting terms...
	cout<<" sys_eAlquasielastic: "<<sys_eAlquasielastic<<"	"<<sys_eAlquasielastic/33.0<<endl;
	cout<<" sys_eAlinelastic: "<<sys_eAlinelastic<<"	"<<sys_eAlinelastic/33.0<<endl;
	cout<<" sys_pim: "<<sys_pim<<"	"<<sys_pim/33.0<<endl;
	
	cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
	cout<<rate0_epelastic/rate0_ee<<endl;
	cout<<rate0_epinelastic/rate0_ee<<endl;
	cout<<rate0_eAlelastic/rate0_ee<<endl;

	
}



		

