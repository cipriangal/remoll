#include<iostream>
#include<fstream>
#include "math.h"
using namespace std;

void yield_counts_asymmetry(TString runtype){

	//TString name_rootfile="/home/yxzhao/Contamination/build/remoll_"+runtype+".root";
	//TString name_rootfile="/home/yxzhao/test/build/remoll_"+runtype+".root";
	//TString name_rootfile="/home/yxzhao/test_10.0/build/remoll_"+runtype+".root";
	
	TString name_rootfile="/home/yxzhao/Contamination/analysis/generate_plots/rootfiles/remoll_"+runtype+".root";
	TString name_txtfile=runtype+".txt";

	const double ring1_rmin=0.69;   // unit in m
	const double ring1_rmax=0.730;
	
	const double ring2_rmin=0.730;
	const double ring2_rmax=0.780;
	
	const double ring3_rmin=0.780;
	const double ring3_rmax=0.855;
	
	const double ring4_rmin=0.855;
	const double ring4_rmax=0.930;

	const double ring5_rmin=0.935;
	const double ring5_rmax=1.1;
	const double ring5_rmin_open=0.935;
	const double ring5_rmax_open=1.04;
	const double ring5_rmin_transition=0.96;
	const double ring5_rmax_transition=1.075;
	const double ring5_rmin_closed=0.96;
	const double ring5_rmax_closed=1.1;

	const double ring6_rmin=1.1;
	const double ring6_rmax=1.2;

	//////////////////////////////////////////////////////////
	// moller rootfile (ee processes)
	/////////////////////////////////////////////////////////
	TChain *T=new TChain("T","T");
	T->AddFile(name_rootfile);

	double ev_A=0;
	double ev_xs=0;     //ub
	double rate=0;  //weighting factor
	int hit_n=0;
	int hit_det[100]={0};
	int hit_pid[100]={0};
	int hit_trid[100]={0};
	int hit_mtrid[100]={0};
	double hit_x[100]={0};
	double hit_y[100]={0};
	double hit_z[100]={0};
	double hit_r[100]={0};
	double hit_p[100]={0};
	double hit_px[100]={0};
	double hit_py[100]={0};
	double hit_pz[100]={0};
	double hit_vx[100]={0};
	double hit_vy[100]={0};
	double hit_vz[100]={0};
	double hit_e[100]={0};
	double hit_m[100]={0};

	T->SetBranchAddress("ev.A",&ev_A);
	T->SetBranchAddress("ev.xs",&ev_xs);
	T->SetBranchAddress("hit.n",&hit_n);
	T->SetBranchAddress("rate",&rate);

	T->SetBranchAddress("hit.det",hit_det);
	T->SetBranchAddress("hit.pid",hit_pid);
	T->SetBranchAddress("hit.trid",hit_trid);
	T->SetBranchAddress("hit.mtrid",hit_mtrid);
	T->SetBranchAddress("hit.x",hit_x);
	T->SetBranchAddress("hit.y",hit_y);
	T->SetBranchAddress("hit.z",hit_z);
	T->SetBranchAddress("hit.r",hit_r);
	T->SetBranchAddress("hit.p",hit_p);
	T->SetBranchAddress("hit.px",hit_px);
	T->SetBranchAddress("hit.py",hit_py);
	T->SetBranchAddress("hit.pz",hit_pz);
	T->SetBranchAddress("hit.vx",hit_vx);
	T->SetBranchAddress("hit.vy",hit_vy);
	T->SetBranchAddress("hit.vz",hit_vz);
	T->SetBranchAddress("hit.e",hit_e);
	T->SetBranchAddress("hit.m",hit_m);

	double N_T_entries=T->GetEntries();
	cout<<N_T_entries<<endl;
	//TH1F* h_xs_r=new TH1F("h_xs_r","xs vs hit r",100,-5,20);
	TH1F* h_asy_ring1[3];
	h_asy_ring1[0]=new TH1F("h_asy_ring1_0","scaler to measure asymmetry",100,-1000000,1000000);   //ppb	
	h_asy_ring1[1]=new TH1F("h_asy_ring1_1","scaler to measure asymmetry",100,-1000000,1000000);   	
	h_asy_ring1[2]=new TH1F("h_asy_ring1_2","scaler to measure asymmetry",100,-1000000,1000000);   
	
	TH1F* h_asy_ring2[3];
	h_asy_ring2[0]=new TH1F("h_asy_ring2_0","scaler to measure asymmetry",100,-1000000,1000000);   //ppb	
	h_asy_ring2[1]=new TH1F("h_asy_ring2_1","scaler to measure asymmetry",100,-1000000,1000000);   	
	h_asy_ring2[2]=new TH1F("h_asy_ring2_2","scaler to measure asymmetry",100,-1000000,1000000);   
	
	TH1F* h_asy_ring3[3];
	h_asy_ring3[0]=new TH1F("h_asy_ring3_0","scaler to measure asymmetry",100,-1000000,1000000);   //ppb	
	h_asy_ring3[1]=new TH1F("h_asy_ring3_1","scaler to measure asymmetry",100,-1000000,1000000);   	
	h_asy_ring3[2]=new TH1F("h_asy_ring3_2","scaler to measure asymmetry",100,-1000000,1000000);   
	
	TH1F* h_asy_ring4[3];
	h_asy_ring4[0]=new TH1F("h_asy_ring4_0","scaler to measure asymmetry",100,-1000000,1000000);   //ppb???	
	h_asy_ring4[1]=new TH1F("h_asy_ring4_1","scaler to measure asymmetry",100,-1000000,1000000);   	
	h_asy_ring4[2]=new TH1F("h_asy_ring4_2","scaler to measure asymmetry",100,-1000000,1000000);   
	
	TH1F* h_asy_ring5[3];
	h_asy_ring5[0]=new TH1F("h_asy_ring5_0","scaler to measure asymmetry",100,-1000000,1000000);   //ppb???	
	h_asy_ring5[1]=new TH1F("h_asy_ring5_1","scaler to measure asymmetry",100,-1000000,1000000);   	
	h_asy_ring5[2]=new TH1F("h_asy_ring5_2","scaler to measure asymmetry",100,-1000000,1000000);   
	
	TH1F* h_asy_ring6[3];
	h_asy_ring6[0]=new TH1F("h_asy_ring6_0","scaler to measure asymmetry",100,-1000000,1000000);   //ppb???	
	h_asy_ring6[1]=new TH1F("h_asy_ring6_1","scaler to measure asymmetry",100,-1000000,1000000);   	
	h_asy_ring6[2]=new TH1F("h_asy_ring6_2","scaler to measure asymmetry",100,-1000000,1000000);   
	
	////////////////////////////////////////////////
	//currently only use rate to calculate counts
	//assuming it is right
	////////////////////////////////////////////////
	double sum_rate_ring1[3]={0};  //[0] for close, [1] for transition, [2] for open
	double sum_rate_ring2[3]={0};  
	double sum_rate_ring3[3]={0};  
	double sum_rate_ring4[3]={0};  
	double sum_rate_ring5[3]={0};  
	double sum_rate_ring6[3]={0};  
	
	double hit_phi=0;
	double sec_phi=0;
	//go through events
	for(int i=0;i<N_T_entries;i++){
		T->GetEntry(i);
		// go through hits
		for(int j=0;j<hit_n;j++){
			hit_phi=0;
			sec_phi=0;

			if(hit_z[j]>26){   // all possible cuts here...
				
				
				//---ring 1
				if(hit_r[j]>=ring1_rmin && hit_r[j]<ring1_rmax){
					hit_phi=atan2(hit_y[j],hit_x[j]);
					if( hit_phi < 0 ) {
						hit_phi += 2.0*3.14159;   //now hit_phi is in [0,2pi]
					}
					sec_phi=fmod(hit_phi, 2.0*3.14159/7 );
					
					if(sec_phi< 3.14159/28.){  //closed
						sum_rate_ring1[0] += rate;
						h_asy_ring1[0]->Fill(ev_A,rate);
					}else if(sec_phi < 3.0*3.14159/28.){ //transition
						sum_rate_ring1[1] += rate;
						h_asy_ring1[1]->Fill(ev_A,rate);
					}else if(sec_phi < 5.0*3.14159/28.){ //open
						sum_rate_ring1[2] +=rate;
						h_asy_ring1[2]->Fill(ev_A,rate);
					}else if(sec_phi < 7.0*3.14159/28.){ //transition
						sum_rate_ring1[1] +=rate;
						h_asy_ring1[1]->Fill(ev_A,rate);
					}else if(sec_phi < 8.0*3.14159/28.){ //closed
						sum_rate_ring1[0] +=rate;
						h_asy_ring1[0]->Fill(ev_A,rate);
					}
				}
				
				//---ring2
				if(hit_r[j]>=ring2_rmin && hit_r[j]<ring2_rmax){
					hit_phi=atan2(hit_y[j],hit_x[j]);
					if( hit_phi < 0 ) {
						hit_phi += 2.0*3.14159;   //now hit_phi is in [0,2pi]
					}
					sec_phi=fmod(hit_phi, 2.0*3.14159/7 );
					
					if(sec_phi< 3.14159/28.){  //closed
						sum_rate_ring2[0] += rate;
						h_asy_ring2[0]->Fill(ev_A,rate);
					}else if(sec_phi < 3.0*3.14159/28.){ //transition
						sum_rate_ring2[1] += rate;
						h_asy_ring2[1]->Fill(ev_A,rate);
					}else if(sec_phi < 5.0*3.14159/28.){ //open
						sum_rate_ring2[2] +=rate;
						h_asy_ring2[2]->Fill(ev_A,rate);
					}else if(sec_phi < 7.0*3.14159/28.){ //transition
						sum_rate_ring2[1] +=rate;
						h_asy_ring2[1]->Fill(ev_A,rate);
					}else if(sec_phi < 8.0*3.14159/28.){ //closed
						sum_rate_ring2[0] +=rate;
						h_asy_ring2[0]->Fill(ev_A,rate);
					}
				}
				
				//---ring3	
				if(hit_r[j]>=ring3_rmin && hit_r[j]<ring3_rmax){
					hit_phi=atan2(hit_y[j],hit_x[j]);
					if( hit_phi < 0 ) {
						hit_phi += 2.0*3.14159;   //now hit_phi is in [0,2pi]
					}
					sec_phi=fmod(hit_phi, 2.0*3.14159/7 );
					
					if(sec_phi< 3.14159/28.){  //closed
						sum_rate_ring3[0] += rate;
						h_asy_ring3[0]->Fill(ev_A,rate);
					}else if(sec_phi < 3.0*3.14159/28.){ //transition
						sum_rate_ring3[1] += rate;
						h_asy_ring3[1]->Fill(ev_A,rate);
					}else if(sec_phi < 5.0*3.14159/28.){ //open
						sum_rate_ring3[2] +=rate;
						h_asy_ring3[2]->Fill(ev_A,rate);
					}else if(sec_phi < 7.0*3.14159/28.){ //transition
						sum_rate_ring3[1] +=rate;
						h_asy_ring3[1]->Fill(ev_A,rate);
					}else if(sec_phi < 8.0*3.14159/28.){ //closed
						sum_rate_ring3[0] +=rate;
						h_asy_ring3[0]->Fill(ev_A,rate);
					}
				}
				
				//---ring4	
				if(hit_r[j]>=ring4_rmin && hit_r[j]<ring4_rmax){
					hit_phi=atan2(hit_y[j],hit_x[j]);
					if( hit_phi < 0 ) {
						hit_phi += 2.0*3.14159;   //now hit_phi is in [0,2pi]
					}
					sec_phi=fmod(hit_phi, 2.0*3.14159/7 );
					
					if(sec_phi< 3.14159/28.){  //closed
						sum_rate_ring4[0] += rate;
						h_asy_ring4[0]->Fill(ev_A,rate);
					}else if(sec_phi < 3.0*3.14159/28.){ //transition
						sum_rate_ring4[1] += rate;
						h_asy_ring4[1]->Fill(ev_A,rate);
					}else if(sec_phi < 5.0*3.14159/28.){ //open
						sum_rate_ring4[2] +=rate;
						h_asy_ring4[2]->Fill(ev_A,rate);
					}else if(sec_phi < 7.0*3.14159/28.){ //transition
						sum_rate_ring4[1] +=rate;
						h_asy_ring4[1]->Fill(ev_A,rate);
					}else if(sec_phi < 8.0*3.14159/28.){ //closed
						sum_rate_ring4[0] +=rate;
						h_asy_ring4[0]->Fill(ev_A,rate);
					}
				}
				
				//---ring5
				//---ring 5 open, transition and closed have different r range 
				if(hit_r[j]>=ring5_rmin && hit_r[j]<ring5_rmax){
					hit_phi=atan2(hit_y[j],hit_x[j]);
					if( hit_phi < 0 ) {
						hit_phi += 2.0*3.14159;   //now hit_phi is in [0,2pi]
					}
					sec_phi=fmod(hit_phi, 2.0*3.14159/7 );
					
					if(sec_phi< 3.14159/28.){  //closed
						if(hit_r[j]>=ring5_rmin_closed && hit_r[j]<ring5_rmax_closed){ //additional r constrain
							sum_rate_ring5[0] += rate;
							h_asy_ring5[0]->Fill(ev_A,rate);
						}
					
					}else if(sec_phi < 3.0*3.14159/28.){ //transition
						if(hit_r[j]>=ring5_rmin_transition && hit_r[j]<ring5_rmax_transition){	
							sum_rate_ring5[1] += rate;
							h_asy_ring5[1]->Fill(ev_A,rate);
						}
					
					}else if(sec_phi < 5.0*3.14159/28.){ //open
						if(hit_r[j]>=ring5_rmin_open && hit_r[j]<ring5_rmax_open){
							sum_rate_ring5[2] +=rate;
							h_asy_ring5[2]->Fill(ev_A,rate);
						}
					}else if(sec_phi < 7.0*3.14159/28.){ //transition
						if(hit_r[j]>=ring5_rmin_transition && hit_r[j]<ring5_rmax_transition){
							sum_rate_ring5[1] +=rate;
							h_asy_ring5[1]->Fill(ev_A,rate);
						}
					}else if(sec_phi < 8.0*3.14159/28.){ //closed
						if(hit_r[j]>=ring5_rmin_closed && hit_r[j]<ring5_rmax_closed){
							sum_rate_ring5[0] +=rate;
							h_asy_ring5[0]->Fill(ev_A,rate);
						}
					}
				}
			
				//---ring6	
				if(hit_r[j]>=ring6_rmin && hit_r[j]<ring6_rmax){
					hit_phi=atan2(hit_y[j],hit_x[j]);
					if( hit_phi < 0 ) {
						hit_phi += 2.0*3.14159;   //now hit_phi is in [0,2pi]
					}
					sec_phi=fmod(hit_phi, 2.0*3.14159/7 );
					
					if(sec_phi< 3.14159/28.){  //closed
						sum_rate_ring6[0] += rate;
						h_asy_ring6[0]->Fill(ev_A,rate);
					}else if(sec_phi < 3.0*3.14159/28.){ //transition
						sum_rate_ring6[1] += rate;
						h_asy_ring6[1]->Fill(ev_A,rate);
					}else if(sec_phi < 5.0*3.14159/28.){ //open
						sum_rate_ring6[2] +=rate;
						h_asy_ring6[2]->Fill(ev_A,rate);
					}else if(sec_phi < 7.0*3.14159/28.){ //transition
						sum_rate_ring6[1] +=rate;
						h_asy_ring6[1]->Fill(ev_A,rate);
					}else if(sec_phi < 8.0*3.14159/28.){ //closed
						sum_rate_ring6[0] +=rate;
						h_asy_ring6[0]->Fill(ev_A,rate);
					}
				}
			

				
			} //end cuts
		}//end hit loop
	}//end event loop


	//////////////////////////////////////////////////////
	//Draw plots and do outputs
	/////////////////////////////////////////////////////
	
//	TCanvas *can = new TCanvas("can","can");
//	can->cd();
//	h_xs_r->Draw();
	//h_xs_p->Draw();
	//h_xs_n->Draw();
//	cout<<sum<<endl;
//	cout<<h_xs_r->Integral()<<endl;

	ofstream output_counts_asy;
	output_counts_asy.open(name_txtfile);
	if(!output_counts_asy){
		cout<<"Error!!! Can't open output file"<<endl;
		return 1;
	}
	output_counts_asy.precision(10);
	
	//output ring 1, in each ring closed=0, transition=1, open=2
	output_counts_asy<<1<<"	"<<0<<"	"<<sum_rate_ring1[0]<<"	"<<h_asy_ring1[0]->GetMean()<<endl;
	output_counts_asy<<1<<"	"<<1<<"	"<<sum_rate_ring1[1]<<"	"<<h_asy_ring1[1]->GetMean()<<endl;
	output_counts_asy<<1<<"	"<<2<<"	"<<sum_rate_ring1[2]<<"	"<<h_asy_ring1[2]->GetMean()<<endl;

	//output ring 2, in each ring closed=0, transition=1, open=2
	output_counts_asy<<2<<"	"<<0<<"	"<<sum_rate_ring2[0]<<"	"<<h_asy_ring2[0]->GetMean()<<endl;
	output_counts_asy<<2<<"	"<<1<<"	"<<sum_rate_ring2[1]<<"	"<<h_asy_ring2[1]->GetMean()<<endl;
	output_counts_asy<<2<<"	"<<2<<"	"<<sum_rate_ring2[2]<<"	"<<h_asy_ring2[2]->GetMean()<<endl;

	//output ring 3, in each ring closed=0, transition=1, open=2
	output_counts_asy<<3<<"	"<<0<<"	"<<sum_rate_ring3[0]<<"	"<<h_asy_ring3[0]->GetMean()<<endl;
	output_counts_asy<<3<<"	"<<1<<"	"<<sum_rate_ring3[1]<<"	"<<h_asy_ring3[1]->GetMean()<<endl;
	output_counts_asy<<3<<"	"<<2<<"	"<<sum_rate_ring3[2]<<"	"<<h_asy_ring3[2]->GetMean()<<endl;

	//output ring 4, in each ring closed=0, transition=1, open=2
	output_counts_asy<<4<<"	"<<0<<"	"<<sum_rate_ring4[0]<<"	"<<h_asy_ring4[0]->GetMean()<<endl;
	output_counts_asy<<4<<"	"<<1<<"	"<<sum_rate_ring4[1]<<"	"<<h_asy_ring4[1]->GetMean()<<endl;
	output_counts_asy<<4<<"	"<<2<<"	"<<sum_rate_ring4[2]<<"	"<<h_asy_ring4[2]->GetMean()<<endl;

	//output ring 5, in each ring closed=0, transition=1, open=2
	output_counts_asy<<5<<"	"<<0<<"	"<<sum_rate_ring5[0]<<"	"<<h_asy_ring5[0]->GetMean()<<endl;
	output_counts_asy<<5<<"	"<<1<<"	"<<sum_rate_ring5[1]<<"	"<<h_asy_ring5[1]->GetMean()<<endl;
	output_counts_asy<<5<<"	"<<2<<"	"<<sum_rate_ring5[2]<<"	"<<h_asy_ring5[2]->GetMean()<<endl;

	//output ring 6, in each ring closed=0, transition=1, open=2
	output_counts_asy<<6<<"	"<<0<<"	"<<sum_rate_ring6[0]<<"	"<<h_asy_ring6[0]->GetMean()<<endl;
	output_counts_asy<<6<<"	"<<1<<"	"<<sum_rate_ring6[1]<<"	"<<h_asy_ring6[1]->GetMean()<<endl;
	output_counts_asy<<6<<"	"<<2<<"	"<<sum_rate_ring6[2]<<"	"<<h_asy_ring6[2]->GetMean()<<endl;

	
}




