double Rate[3][6];

void rate(){

  gStyle->SetOptStat(0);
  string fnm;
  string hnm[3]={"open","transition","closed"};
  TH1D *rMol[3],*rPiM[3];
  for(int i=0;i<3;i++){
    cout<<i<<" "<<Form("Moller %s;r[m];rate[Hz]",hnm[i].c_str())<<" "<<endl;
    rMol[i]=new TH1D(Form("rMol_%d",i),Form("Moller %s;r[m];rate[Hz]",hnm[i].c_str()),400,0.6,1.3);
    rPiM[i]=new TH1D(Form("rPiM_%d",i),Form("#pi -  %s;r[m];rate[Hz]",hnm[i].c_str()),400,0.6,1.3);
  }


  
  cout<<rPiM[0]->GetTitle()<<endl;
  fnm="../output/remollout_Pion_2e6_g4963.root";
  //fnm="../output/remollout_1e6_Pion.root";
  rate1(rPiM,fnm,-211);//pi-
  integrate(rPiM);

  //fnm="../output/remollout_1e7_Moller.root";
  fnm="../output/remollout_Moller_2e6_g4963.root";
  rate1(rMol,fnm,11);//e-
  integrate(rMol);
  
  draw(rMol,rPiM);

}

void rate1(TH1 *rt[3], string fnm,int partID){

  for(int i=0;i<3;i++)
    for(int j=0;j<6;j++) Rate[i][j]=0;
  cout<<"running "<<rt[0]->GetTitle()<<" "<<fnm.c_str()<<" "<<partID<<endl;

  TFile *fin=TFile::Open(fnm.c_str(),"READ");
  TTree *t=(TTree*)fin->Get("T");

  double rate,hitR[10000],hitZ[10000],hitPh[10000];
  int pid[10000],det[10000];
  int nhit,colCut;
  t->SetBranchAddress("hit.n",&nhit);
  t->SetBranchAddress("rate",&rate);
  t->SetBranchAddress("hit.r",&hitR);
  t->SetBranchAddress("hit.z",&hitZ);
  t->SetBranchAddress("hit.pid",&pid);
  t->SetBranchAddress("hit.det",&det);
  t->SetBranchAddress("hit.ph",&hitPh);
  t->SetBranchAddress("hit.colCut",&colCut);

  int nev=t->GetEntries();
  for(int i=0;i<nev;i++){    
    t->GetEntry(i);
    if(nhit==0) continue;
    for(int j=0;j<nhit;j++){
      if(colCut!=1) continue;
      if(det[j]!=28) continue;
      if(pid[j]!=partID) continue;
      if(hitZ[j]<=26) continue;
      if(hitR[j]<0.69 || hitR[j]>1.2) continue;

      int phSect=phiSect(hitPh[j]);
      rt[phSect]->Fill(hitR[j],rate);
      sumUp(hitR[j],phSect,rate);
    }
  }
  fin->Close();
  for(int i=0;i<6;i++){
    cout<<i<<endl;
    for(int j=0;j<3;j++)
      cout<<" "<<Rate[j][i];
    cout<<endl;
  }
}

void draw(TH1* r1[3],TH1* r2[3]){
  TCanvas *c1=new TCanvas("c1","c1",1600,900);
  TCanvas *c2=new TCanvas("c2","c2",1600,900);

  c1->cd();
  c1->Divide(2);
  for(int i=0;i<3;i++){
    c1->cd(1);
    r1[i]->SetLineColor(i+1);
    r1[i]->SetLineStyle(1);
    if(i==0) {
      r1[i]->GetYaxis()->SetRangeUser(1,8000000000);
      r1[i]->DrawCopy();
    }
    else r1[i]->DrawCopy("same");

    c1->cd(2);
    r2[i]->SetLineColor(i+1);
    r2[i]->SetLineStyle(2);
    if(i==0) {
      r2[i]->GetYaxis()->SetRangeUser(1,2500000);
      r2[i]->DrawCopy();
    }
    else r2[i]->DrawCopy("same");
  }

  c2->cd();
  r1[0]->DrawCopy();
  for(int i=0;i<3;i++){
    r1[i]->DrawCopy("same");
    r2[i]->DrawCopy("same");
  }

  
}

void sumUp(double r,int phi, double val){
  double edge[3][7]={{0.690,0.730,0.780,0.855,0.935,1.04 ,1.2},//open
  		     {0.690,0.730,0.780,0.855,0.960,1.075,1.2},//transition
  		     {0.690,0.730,0.780,0.855,0.960,1.1  ,1.2}};//closed
  int n=-1;
  for(int i=0;i<6;i++)
    if(r>=edge[phi][i] && r<edge[phi][i+1]){
      //cout<<n<<" "<<r<<" "<<edge[phi][i]<<" "<<edge[phi][i+1]<<endl;
      n=i;
    }
  //cout<<endl<<endl<<n<<endl<<endl;
  if(n==-1) cout<<"problems "<<r<<endl;
  Rate[phi][n]+=val;
}

void integrate(TH1 *rt[3]){

  double edge[3][7]={{0.690,0.730,0.780,0.855,0.935,1.04 ,1.2},//open
  		     {0.690,0.730,0.780,0.855,0.960,1.075,1.2},//transition
  		     {0.690,0.730,0.780,0.855,0.960,1.1  ,1.2}};//closed

  // double edge[3][7]={{0.690,0.730,0.780,0.855,0.930,1.1 ,1.2},//open
  // 		     {0.690,0.730,0.780,0.855,0.930,1.1,1.2},//transition
  // 		     {0.690,0.730,0.780,0.855,0.930,1.1  ,1.2}};//closed
  cout<<rt[0]->GetTitle()<<endl;
  for(int i=0;i<6;i++){
    cout<<" "<<i;    
    for(int j=0;j<3;j++){
      int b1=rt[j]->GetXaxis()->FindBin(edge[j][i]);
      int b2=rt[j]->GetXaxis()->FindBin(edge[j][i+1]);
      cout<<" "<<rt[j]->Integral(b1,b2);
      //cout<<i<<" ~~ "<<rt->Integral(b1,b2,"width")<<endl; //this is bin dependent
    }
    cout<<endl;
  }
}

int phiSect(double ph){
  double phi=ph;
  if(phi<0) phi+=360;
  double sectphi=fmod(phi,360./7.);

  if( sectphi < 180./28 ) return 0;
  else if( sectphi < 3.*180./28. ) return 1;
  else if( sectphi < 5.*180./28. ) return 2;
  else if( sectphi < 7.*180./28. ) return 1;
  else return 0;
}
