#include "TFile.h"
#include "TH1F.h"
#include "TString.h"

using std::cout;
using std::endl;

void addHistos(){
	unsigned int nBack = 7;
	//Double_t lumi = 1.0;
	Double_t lumi = 27.27;// run BtoG in data	

	TString path = "ZprimeFV500deltabs0Files/";

	TFile *f[nBack];
	f[0] = TFile::Open(path+"mjj_QCD_HT2000toinftrigger.root");
	f[1] = TFile::Open(path+"mjj_QCD_HT200to300trigger.root");
	f[2] = TFile::Open(path+"mjj_QCD_HT300to500trigger.root");
	f[3] = TFile::Open(path+"mjj_QCD_HT500to700trigger.root");
	f[4] = TFile::Open(path+"mjj_QCD_HT700to1000trigger.root");
	f[5] = TFile::Open(path+"mjj_QCD_HT1000to1500trigger.root");
	f[6] = TFile::Open(path+"mjj_QCD_HT1500to2000trigger.root");

	TH1F*h1[nBack];
	TH1F*h2[nBack];
	TH1F*h3[nBack];
	TH1F*h4[nBack];
	TH1F*h5[nBack];
	TH1F*h6[nBack];					

	//std::cout<<"Path to file 2000toinf: " << f[0]->GetPath() << std::endl;
	h1[0] = (TH1F*) f[0]->Get("mjjsr1");
	h2[0] = (TH1F*) f[0]->Get("mjjsr2");
	h3[0] = (TH1F*) f[0]->Get("mjjsr3");
	h4[0] = (TH1F*) f[0]->Get("mjjcr1");
	h5[0] = (TH1F*) f[0]->Get("mjjcr2");
	h6[0] = (TH1F*) f[0]->Get("mjjcr3");					
	//std::cout << "n=0: " << " Entries in h1n: " << h1[0]->GetEntries()<<std::endl;
	h1[0]->Scale(lumi*25240./1990837.);
	h2[0]->Scale(lumi*25240./1990837.);
	h3[0]->Scale(lumi*25240./1990837.);
	h4[0]->Scale(lumi*25240./1990837.);
	h5[0]->Scale(lumi*25240./1990837.);
	h6[0]->Scale(lumi*25240./1990837.);					
	std::cout << "n=0 Entries after re-weight cr1: " << h4[0]->Integral(0,76)<<std::endl;	
	std::cout << "n=0 Entries after re-weight cr2: " << h5[0]->Integral(0,76)<<std::endl;	
	std::cout << "n=0 Entries after re-weight cr3: " << h6[0]->Integral(0,76)<<std::endl;			

	for(unsigned int n=1; n<nBack; n++){
		h1[n] = (TH1F*) f[n]->Get("mjjsr1");
		h2[n] = (TH1F*) f[n]->Get("mjjsr2");
		h3[n] = (TH1F*) f[n]->Get("mjjsr3");
		h4[n] = (TH1F*) f[n]->Get("mjjcr1");
		h5[n] = (TH1F*) f[n]->Get("mjjcr2");
		h6[n] = (TH1F*) f[n]->Get("mjjcr3");
		//std::cout << "n: " << n << " Entries in h1n: " << h1[n]->GetEntries()<<std::endl;
		//std::cout << "n="<<n<<" Entries cr1: " << h4[n]->GetEntries()<<std::endl;	
		//std::cout << "n="<<n<<" Entries cr2: " << h5[n]->GetEntries()<<std::endl;	
		//std::cout << "n="<<n<<" Entries cr3: " << h6[n]->GetEntries()<<std::endl;	

		Double_t weight = 0;
		if(n==1) weight = lumi*1712000000./16805485.;
		if(n==2) weight = lumi*347700000./16808178.;//QCD_HT300to500
		if(n==3) weight = lumi*32100000./10013260.;//QCD_HT500to700
		if(n==4) weight = lumi*6831000./15615291.;//QCD_HT700to1000
		if(n==5) weight = lumi*1207000./4664281.;//QCD_HT1000to1500
		if(n==6) weight = lumi*119900./3968837.;//QCD_HT1500to2000
		
		h1[n]->Scale(weight);
		h2[n]->Scale(weight);
		h3[n]->Scale(weight);
		h4[n]->Scale(weight);
		h5[n]->Scale(weight);
		h6[n]->Scale(weight);										

		h1[0]->Add(h1[n]);
		h2[0]->Add(h2[n]);
		h3[0]->Add(h3[n]);
		h4[0]->Add(h4[n]);
		h5[0]->Add(h5[n]);
		h6[0]->Add(h6[n]);

	std::cout << "n="<<n<<" Entries after re-weight cr1: " << h4[n]->Integral(0,76)<<std::endl;	
	std::cout << "n="<<n<<" Entries after re-weight cr2: " << h5[n]->Integral(0,76)<<std::endl;	
	std::cout << "n="<<n<<" Entries after re-weight cr3: " << h6[n]->Integral(0,76)<<std::endl;			

	}
	std::cout << "Entries in sr1 sum: " << h1[0]->GetEntries() << " Integral: " <<h1[0]->Integral(0,76) << std:: endl;
	std::cout << "Entries in sr2 sum: " << h2[0]->GetEntries() << " Integral: " <<h2[0]->Integral(0,76) << std:: endl;
	std::cout << "Entries in sr3 sum: " << h3[0]->GetEntries() << " Integral: " <<h3[0]->Integral(0,76) << std:: endl;
	std::cout << "Entries in cr1 sum: " << h4[0]->GetEntries() << " Integral: " <<h4[0]->Integral(0,76) << std:: endl;
	std::cout << "Entries in cr2 sum: " << h5[0]->GetEntries() << " Integral: " <<h5[0]->Integral(0,76) << std:: endl;
	std::cout << "Entries in cr3 sum: " << h6[0]->GetEntries() << " Integral: " <<h6[0]->Integral(0,76) << std:: endl;					

	TFile* outfile = new TFile(path+"mjj_QCDOldAllB2Gtrigger.root","RECREATE");
	h1[0]->Write();
	h2[0]->Write();
	h3[0]->Write();
	h4[0]->Write();
	h5[0]->Write();
	h6[0]->Write();					
	outfile->Close();
	delete outfile;

}
