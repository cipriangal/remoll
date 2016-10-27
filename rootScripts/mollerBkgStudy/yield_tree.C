#include<iostream>
#include<fstream>
using namespace std;


void yield_tree(){

	ifstream INPUT_moller;
	INPUT_moller.open("moller.txt");
	if(!INPUT_moller){
		cout<<"Error!!! Can't find input file"<<endl;
		return 1;
	}
	
	ifstream INPUT_epelastic;
	INPUT_epelastic.open("epelastic.txt");
	if(!INPUT_epelastic){
		cout<<"Error!!! Can't find input file"<<endl;
		return 1;
	}

	ifstream INPUT_epinelastic;
	INPUT_epinelastic.open("epinelastic.txt");
	if(!INPUT_epinelastic){
		cout<<"Error!!! Can't find input file"<<endl;
		return 1;
	}

	int ring_number=0;
	int sector_number=0;
	double rate_ee=0, rate_epelastic=0, rate_epinelastic=0;
	double A_ee=0, A_epelastic=0, A_epinelastic=0;

	TFile *myfile=new TFile("Counts_vs_asy.root","RECREATE");
	TTree *T=new TTree("T","T");
	
	T->Branch("ring_number",&ring_number,"ring_number/I");
	T->Branch("sector_number",&sector_number,"sector_number/I");
	T->Branch("rate_ee",&rate_ee,"rate_ee/D");
	T->Branch("rate_epelastic",&rate_epelastic,"rate_epelastic/D");
	T->Branch("rate_epinelastic",&rate_epinelastic,"rate_epinelastic/D");
	T->Branch("A_ee",&A_ee,"A_ee/D");
	T->Branch("A_epelastic",&A_epelastic,"A_epelastic/D");
	T->Branch("A_epinelastic",&A_epinelastic,"A_epinelastic/D");

	int tmp_ring_number_ee, tmp_ring_number_epelastic, tmp_ring_number_epinelastic;
	int tmp_sector_number_ee, tmp_sector_number_epelastic, tmp_sector_number_epinelastic;

	while (INPUT_moller>>tmp_ring_number_ee>>tmp_sector_number_ee>>rate_ee>>A_ee
		&& INPUT_epelastic>>tmp_ring_number_epelastic>>tmp_sector_number_epelastic>>rate_epelastic>>A_epelastic
		&& INPUT_epinelastic>>tmp_ring_number_epinelastic>>tmp_sector_number_epinelastic>>rate_epinelastic>>A_epinelastic){

		//double check ring number and sector number
		if(tmp_ring_number_ee==tmp_ring_number_epelastic && tmp_ring_number_ee==tmp_ring_number_epinelastic){
			ring_number=tmp_ring_number_ee;
		}else{
			cout<<" Something should be wrong !!! ring number doesn't match!!!"<<endl;
		}

		if(tmp_sector_number_ee==tmp_sector_number_epelastic && tmp_sector_number_ee==tmp_sector_number_epinelastic){
			sector_number=tmp_sector_number_ee;
		}else{
			cout<<" Something should be wrong !!! sector number doesn't match!!!"<<endl;
		}

		T->Fill();
	}

	myfile->Write();
	myfile->Close();
}



