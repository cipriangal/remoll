#include<iostream>
#include<fstream>
using namespace std;


void yield_tree_window(){

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
	
	ifstream INPUT_eAlelastic;
	INPUT_eAlelastic.open("eAlelastic.txt");
	if(!INPUT_eAlelastic){
		cout<<"Error!!! Can't find input file"<<endl;
		return 1;
	}
	
	ifstream INPUT_eAlquasielastic;
	INPUT_eAlquasielastic.open("eAlquasielastic.txt");
	if(!INPUT_eAlquasielastic){
		cout<<"Error!!! Can't find input file"<<endl;
		return 1;
	}
	
	ifstream INPUT_eAlinelastic;
	INPUT_eAlinelastic.open("eAlinelastic.txt");
	if(!INPUT_eAlinelastic){
		cout<<"Error!!! Can't find input file"<<endl;
		return 1;
	}

	ifstream INPUT_pim;
	INPUT_pim.open("pim.txt");
	if(!INPUT_pim){
		cout<<"Error!!! Can't find input file"<<endl;
		return 1;
	}




	int ring_number=0;
	int sector_number=0;
	double rate_ee=0, rate_epelastic=0, rate_epinelastic=0, rate_eAlelastic=0, rate_eAlquasielastic=0, rate_eAlinelastic=0,rate_pim=0;
	double A_ee=0, A_epelastic=0, A_epinelastic=0, A_eAlelastic=0, A_eAlquasielastic=0, A_eAlinelastic=0, A_pim;
	double rate_tot=0;
	TFile *myfile=new TFile("Counts_vs_asy.root","RECREATE");
	TTree *T=new TTree("T","T");
	
	T->Branch("ring_number",&ring_number,"ring_number/I");
	T->Branch("sector_number",&sector_number,"sector_number/I");
	T->Branch("rate_ee",&rate_ee,"rate_ee/D");
	T->Branch("rate_epelastic",&rate_epelastic,"rate_epelastic/D");
	T->Branch("rate_epinelastic",&rate_epinelastic,"rate_epinelastic/D");
	T->Branch("rate_eAlelastic",&rate_eAlelastic,"rate_eAlelastic/D");
	T->Branch("rate_eAlquasielastic",&rate_eAlquasielastic,"rate_eAlquasielastic/D");
	T->Branch("rate_eAlinelastic",&rate_eAlinelastic,"rate_eAlinelastic/D");
	T->Branch("rate_pim",&rate_pim,"rate_pim/D");
	T->Branch("rate_tot",&rate_tot,"rate_tot/D");
	
	T->Branch("A_ee",&A_ee,"A_ee/D");
	T->Branch("A_epelastic",&A_epelastic,"A_epelastic/D");
	T->Branch("A_epinelastic",&A_epinelastic,"A_epinelastic/D");
	T->Branch("A_eAlelastic",&A_eAlelastic,"A_eAlelastic/D");
	T->Branch("A_eAlquasielastic",&A_eAlquasielastic,"A_eAlquasielastic/D");
	T->Branch("A_eAlinelastic",&A_eAlinelastic,"A_eAlinelastic/D");
	T->Branch("A_pim",&A_pim,"A_pim/D");

	int tmp_ring_number_ee, tmp_ring_number_epelastic, tmp_ring_number_epinelastic, tmp_ring_number_pim;
	int tmp_ring_number_eAlelastic, tmp_ring_number_eAlquasielastic, tmp_ring_number_eAlinelastic;
	
	int tmp_sector_number_ee, tmp_sector_number_epelastic, tmp_sector_number_epinelastic, tmp_sector_number_pim;
	int tmp_sector_number_eAlelastic, tmp_sector_number_eAlquasielastic, tmp_sector_number_eAlinelastic;

	while (INPUT_moller>>tmp_ring_number_ee>>tmp_sector_number_ee>>rate_ee>>A_ee
		&& INPUT_epelastic>>tmp_ring_number_epelastic>>tmp_sector_number_epelastic>>rate_epelastic>>A_epelastic
		&& INPUT_epinelastic>>tmp_ring_number_epinelastic>>tmp_sector_number_epinelastic>>rate_epinelastic>>A_epinelastic
		&& INPUT_eAlelastic>>tmp_ring_number_eAlelastic>>tmp_sector_number_eAlelastic>>rate_eAlelastic>>A_eAlelastic
		&& INPUT_eAlquasielastic>>tmp_ring_number_eAlquasielastic>>tmp_sector_number_eAlquasielastic>>rate_eAlquasielastic>>A_eAlquasielastic

		&& INPUT_eAlinelastic>>tmp_ring_number_eAlinelastic>>tmp_sector_number_eAlinelastic>>rate_eAlinelastic>>A_eAlinelastic
		&& INPUT_pim>>tmp_ring_number_pim>>tmp_sector_number_pim>>rate_pim>>A_pim
		
		){

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
		
		rate_tot=rate_ee + rate_epelastic + rate_epinelastic + rate_eAlelastic + rate_eAlquasielastic + rate_eAlinelastic + rate_pim;
		T->Fill();
	}

	myfile->Write();
	myfile->Close();
}



