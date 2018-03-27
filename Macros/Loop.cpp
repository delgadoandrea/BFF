#include <vector>
#include <iostream>
#include <map>
#include <iterator>

#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TChain.h"

using std::cout;
using std::endl;
using std::vector;

class dRComb{
public:

	dRComb(int a, int b, double c): i(a), j(b), dR(c){}
	int i;
	int j;
	double dR;

};

bool sortBySmallestDR (dRComb i, dRComb j){
	return (i.dR < j.dR);
}

bool sortByLargestM (dRComb i, dRComb j){
	return (i.dR > j.dR);
}

map <string, TH1*> histograms;
map <string, TH1*> histograms2D;

//void Loop (TChain *chain, TString name){
void Loop (TString file, TString name){	
	TFile *f = TFile::Open(file);
	TTree *t = (TTree*) f->Get("tree");

   	if (t == 0 ){
		std::cout << "No tree found!" <<std::endl;
		return;
   	} 
	Float_t GenJet_wNuPhi[50];
	Float_t GenJet_wNuEta[50];
	Float_t GenJet_wNuPt[50];
	Float_t GenJet_wNuM[50];
	Float_t GenJet_phi[50];
	Float_t GenJet_eta[50];
	Float_t GenJet_pt[50];
	Float_t GenJet_mass[50];	
	Float_t Jet_pt[50];		
	Float_t Jet_eta[50];		
	Float_t Jet_phi[50];		
	Float_t Jet_mass[50];	
	Float_t Jet_btagCMVA[50];				
	Float_t Jet_btagCSV[50];		
	Int_t Jet_hadronFlavour[50];
	Int_t HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v;

	Int_t nJet;
	Int_t nGenJet;	
	Int_t nGenStatus2bHad;

	t->SetBranchAddress("GenJet_wNuEta",GenJet_wNuEta);			
	t->SetBranchAddress("GenJet_wNuPhi",GenJet_wNuPhi);			
	t->SetBranchAddress("GenJet_wNuPt",GenJet_wNuPt);			
	t->SetBranchAddress("GenJet_wNuM",GenJet_wNuM);	
	t->SetBranchAddress("GenJet_eta",GenJet_eta);			
	t->SetBranchAddress("GenJet_phi",GenJet_phi);			
	t->SetBranchAddress("GenJet_pt",GenJet_pt);			
	t->SetBranchAddress("GenJet_mass",GenJet_mass);		
	t->SetBranchAddress("Jet_pt", Jet_pt);	
	t->SetBranchAddress("Jet_eta", Jet_eta);	
	t->SetBranchAddress("Jet_phi", Jet_phi);	
	t->SetBranchAddress("Jet_mass", Jet_mass);	
	t->SetBranchAddress("Jet_btagCSV", Jet_btagCSV);			
	t->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA);
	t->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);				
	t->SetBranchAddress("nJet",&nJet);			
	t->SetBranchAddress("nGenJet",&nGenJet);			
	t->SetBranchAddress("nGenStatus2bHad",&nGenStatus2bHad);
	t->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v", &HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v);
	t->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v", &HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v);
	t->SetBranchAddress("HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v", &HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v);
	t->SetBranchAddress("HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v", &HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v);	
			
	histograms["mjjsr1"] = new TH1F("mjjsr1", "mjj from two leading jets in event for signal region 1", 75,0, 1500);	
	histograms["mjjsr2"] = new TH1F("mjjsr2", "mjj from two leading jets in event for signal region 2", 75,0, 1500);	
	histograms["mjjsr3"] = new TH1F("mjjsr3", "mjj from two leading jets in event for signal region 3", 75,0, 1500);		
	histograms["mjjcr1"] = new TH1F("mjjcr1", "mjj from two leading jets in event for control region 1", 75,0, 1500);	
	histograms["mjjcr2"] = new TH1F("mjjcr2", "mjj from two leading jets in event for control region 2", 75,0, 1500);		
	histograms["mjjcr3"] = new TH1F("mjjcr3", "mjj from two leading jets in event for control region 3", 75,0, 1500);			

	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Sumw2();
	}

	//Int_t nPass = 0;
	cout<< "Analyzing sample: " << name << endl;
	cout<< "Number of events in t: " << t->GetEntries() << endl;						
	for (Int_t indx=0;indx<t->GetEntries();indx++){ 	

		t->GetEntry(indx);

		Int_t tagCL = 0;
		Int_t tagCM = 0;
		Int_t tagCT = 0;
		Int_t nLoose = 0;
		Int_t nMedium = 0;
		Int_t nTight = 0;
		//bool pass = false;
		Int_t jetT = 0;

		//Int_t reqTight = 2;
		//Int_t reqMed = 1;
		//Int_t reqLoose = 1;		

		for (int j=0; j<nJet; j++){
			if(Jet_pt[j]>30){
				if(fabs(Jet_eta[j])<2.4) jetT++;
				if(Jet_btagCSV[j]>0.9535) tagCT++;
				if(Jet_btagCSV[j]>0.8484) tagCM++;
				if(Jet_btagCSV[j]>0.5426) tagCL++;
				if(Jet_btagCSV[j]>0.9535)++nTight;
				else if(Jet_btagCSV[j]>0.8484)++nMedium;
				else if(Jet_btagCSV[j]>0.5426)++nLoose;		 
			}
		}

		/*while(nLoose+nTight+nMedium!=0){
			if(nTight){
				if(reqTight){--reqTight;--nTight;}
				else if(reqMed){--reqMed;--nTight;}
				else if(reqLoose){--reqLoose;--nTight;}
				else nTight = 0;
			}else if(nMedium){
				if(reqMed){--reqMed;--nMedium;}
				else if(reqLoose){--reqLoose;--nMedium;}
				else nMedium=0;
			}else if(nLoose){
				if(reqLoose){--reqLoose;--nLoose;}
				else nLoose=0;
			}

			if((reqTight+reqMed+reqLoose==0)){pass=true;break;}
//			if(!(reqTight+reqMed+reqLoose)){pass=true;break;}			
		}//while*/

		if(HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v==1 && nJet>3 && Jet_pt[3]>30 ){						
			TLorentzVector j1;
			TLorentzVector j2;
			j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);	
			TLorentzVector s = j1+j2;
			//histograms["mjjsr1"]->Fill(s.M(),1./50.);	
			if(tagCL>2){
				if(Jet_btagCSV[0]>0.9535 && Jet_btagCSV[1]>0.9535) histograms["mjjsr1"]->Fill(s.M());
				if(Jet_btagCSV[0]>0.8484 && Jet_btagCSV[0]< 0.9535 && Jet_btagCSV[1]> 0.8484 && Jet_btagCSV[1]<0.9535) histograms["mjjsr2"]->Fill(s.M());	
				if(Jet_btagCSV[0]>0.5426 && Jet_btagCSV[0]< 0.8484 && Jet_btagCSV[1]> 0.5426 && Jet_btagCSV[1]<0.8484) histograms["mjjsr3"]->Fill(s.M());																																				
			}

			//Loop over elements of A and B finding the dR between all possible combinations
			/*vector <dRComb> mCombs;
	
			for(int b=0; b<nJet; b++){
				for(int a=b+1; a<nJet; a++){

				TLorentzVector m1;
				TLorentzVector m2;
				m1.SetPtEtaPhiM(Jet_pt[b], Jet_eta[b], Jet_phi[b], Jet_mass[b]);
				m2.SetPtEtaPhiM(Jet_pt[a], Jet_eta[a], Jet_phi[a], Jet_mass[a]);				
				TLorentzVector sm = m1+m2;
				mCombs.push_back(dRComb(b,a,sm.M()));
				//cout << "genJ: " << b << " bhad: " << a << " dR: " << dR << endl;
				}
			}

			//sort the vector by smallest M.()
			std::sort(mCombs.begin(), mCombs.end(), sortByLargestM);
			histograms["mjjlm"]->Fill(mCombs.at(0).dR);
			*/					

		}//if pass
		if(nJet>3 && Jet_pt[3]>30 && tagCL<3){
			TLorentzVector j1;
			TLorentzVector j2;
			j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);	
			TLorentzVector s = j1+j2;
			if(HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v==1 || HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v==1){
				if(HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v!=1){
					if(Jet_btagCSV[0]>0.9535 && Jet_btagCSV[1]>0.9535) histograms["mjjcr1"]->Fill(s.M(),0.02);
					if(Jet_btagCSV[0]>0.8484 && Jet_btagCSV[0]< 0.9535 && Jet_btagCSV[1]> 0.8484 && Jet_btagCSV[1]<0.9535) histograms["mjjcr2"]->Fill(s.M(),0.02); 
					if(Jet_btagCSV[0]>0.5426 && Jet_btagCSV[0]< 0.8484 && Jet_btagCSV[1]> 0.5426 && Jet_btagCSV[1]<0.8484) histograms["mjjcr3"]->Fill(s.M(),0.02); 					
				}
				else 
					if(Jet_btagCSV[0]>0.9535 && Jet_btagCSV[1]>0.9535) histograms["mjjcr1"]->Fill(s.M());
					if(Jet_btagCSV[0]>0.8484 && Jet_btagCSV[0]< 0.9535 && Jet_btagCSV[1]> 0.8484 && Jet_btagCSV[1]<0.9535) histograms["mjjcr2"]->Fill(s.M()); 
					if(Jet_btagCSV[0]>0.5426 && Jet_btagCSV[0]< 0.8484 && Jet_btagCSV[1]> 0.5426 && Jet_btagCSV[1]<0.8484) histograms["mjjcr3"]->Fill(s.M()); 										
			} 
		}
	}//loop over entries
	//cout << "Events passing selection: " << nPass << endl;

	//Save invariant mass histogram to a file
	TFile* outfile = new TFile("ZprimeFV500deltabs0Files/mjj_"+name+".root","RECREATE");
//	TFile* outfile = new TFile("mjj_QCD_HT500to700.root","RECREATE");
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Write();	
	}
	outfile->Close();
	delete outfile;

	/*for(std::map<string,TH1*>::iterator it=histograms2D.begin(); it!=histograms2D.end(); ++it){
		histograms2D[it->first]->Sumw2();
		TString canname = Form("%s",it->first.c_str());		
		TCanvas *c = new TCanvas(canname,canname,700,700);		
		histograms2D[it->first]->Draw("colz text");
		TString canname2pdf = "06202017_"+it->first + ".pdf";
		//c->Print(canname2pdf);
	}
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Sumw2();
		Double_t nsw = histograms[it->first]->GetEntries();
		histograms[it->first]->Scale(1./nsw);
		TString canname = Form("%s",it->first.c_str());		
		TCanvas *c = new TCanvas(canname,canname,700,700);
		histograms[it->first]->Draw("hist E");
		TString canname2pdf = Form("06252017_%iL%iM%iT%s.pdf",1,1,1,it->first.c_str());
//		TString canname2pdf = Form("06252017_%iL%iM%iT%s_HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v.pdf",0,1,3,it->first.c_str());		
		//c->Print(canname2pdf);
	}*/

}
