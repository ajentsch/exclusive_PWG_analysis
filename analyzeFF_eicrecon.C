//-------------------------
//
// Simple analysis code to analyze EPIC simulation output for Roman Pots
//
// NOTE: This is for tracks reconstructed using the Roman Pots "online"
//
// Author: Alex Jentsch
//
// Date of last author update: 7/10/2023
//
//------------------------


using namespace std;


void analyzeFF_eicrecon(){

	TString fileList = "./inputFileList_particleGun.list";
	
	TString outputName = "ePIC_fullReco_LOCAL_COORD_RP_Output_";	

	TString date = "7_10_2023_";
	
	TString run  = "run_0";

	cout << "Input FileList: " << fileList << endl;
	TString fileType_ROOT = ".root";
	TString outputFileName = outputName + date + run + fileType_ROOT;
	string fileName;
	TFile * inputRootFile;
	TTree * rootTree;
	cout << "Output file: " << outputFileName << endl;


	ifstream fileListStream;
	fileListStream.open(fileList);
	if(!fileListStream) { cout << "NO_LIST_FILE " << fileList << endl; return;}

	
    //---------------------Roman Pots reconstruction constants------------------
	
	//N.B. this is all bullshit for right now while we solve the online reco problem
	
   

	//--------------------------------------------------------------------------

	//histograms -- only a few for now
	
	//MC information
	TH1D* h_eta_MC = new TH1D("h_eta",";Pseudorapidity, #eta",100,0.0,15.0);
	TH1D* h_px_MC = new TH1D("px_MC", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_py_MC = new TH1D("py_MC", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_pt_MC = new TH1D("pt_MC", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
	TH1D* h_pz_MC = new TH1D("pz_MC", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	
	//Roman pots
	TH1D* h_px_RomanPots = new TH1D("px_RomanPots", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_py_RomanPots = new TH1D("py_RomanPots", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_pt_RomanPots = new TH1D("pt_RomanPots", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
	TH1D* h_pz_RomanPots = new TH1D("pz_RomanPots", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	TH2D* h_rp_occupancy_map = new TH2D("Roman_pots_occupancy_map", "hit y [mm];hit x [mm]", 100, -150, 150, 100, -70, -70);
	
	//B0 tracker hits
	TH2D* h_B0_occupancy_map_layer_0 = new TH2D("B0_occupancy_map_0", "B0_occupancy_map_0", 100, -400, 0, 100, -170, 170);
	TH2D* h_B0_occupancy_map_layer_1 = new TH2D("B0_occupancy_map_1", "B0_occupancy_map_1", 100, -400, 0, 100, -170, 170);
	TH2D* h_B0_occupancy_map_layer_2 = new TH2D("B0_occupancy_map_2", "B0_occupancy_map_2", 100, -400, 0, 100, -170, 170);
	TH2D* h_B0_occupancy_map_layer_3 = new TH2D("B0_occupancy_map_3", "B0_occupancy_map_3", 100, -400, 0, 100, -170, 170);
	TH1D* h_B0_hit_energy_deposit = new TH1D("B0_tracker_hit_energy_deposit", ";Deposited Energy [keV]", 100, 0.0, 500.0);
	
	//B0 EMCAL clusters
	TH2D* h_B0_emcal_occupancy_map = new TH2D("B0_emcal_occupancy_map", "B0_emcal_occupancy_map", 100, -400, 0, 100, -170, 170);
	TH1D* h_B0_emcal_cluster_energy = new TH1D("B0_emcal_cluster_energy", ";Cluster Energy [GeV]", 100, 0.0, 100.0);
	
	//Reconstructed tracks (for usage with B0 too!!)
	TH1D* h_px_reco_track = new TH1D("px_reco_track", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_py_reco_track = new TH1D("py_reco_track", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_pt_reco_track = new TH1D("pt_reco_track", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
	TH1D* h_pz_reco_track = new TH1D("pz_reco_track", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	
	//ZDC EMCAL clusters
	TH2D* h_ZDC_emcal_occupancy_map = new TH2D("ZDC_emcal_occupancy_map", "ZDC_emcal_occupancy_map", 100, -1150, -1050, 100, -60, 60);
	TH1D* h_ZDC_emcal_cluster_energy = new TH1D("ZDC_emcal_cluster_energy", ";Cluster Energy [GeV]", 100, 0.0, 100.0);
	
	int fileCounter = 0;
	int iEvent = 0;

	while(getline(fileListStream, fileName)){

	    TString tmp = fileName;

	    cout << "Input file " << fileCounter << ": " << fileName << endl;

	    inputRootFile = new TFile(tmp);
	    if(!inputRootFile){ cout << "MISSING_ROOT_FILE"<< fileName << endl; continue;}
		
		fileCounter++;

		TTree * evtTree = (TTree*)inputRootFile->Get("events");

		int numEvents = evtTree->GetEntries();

    	TTreeReader tree_reader(evtTree);       // !the tree reader

		//MC particles
    
    	TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
    	TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
    	TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
    	TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
    	TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};
		
		// Roman Pots
	
		//momentum vector
   	 	TTreeReaderArray<float> reco_RP_px = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
    	TTreeReaderArray<float> reco_RP_py = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
    	TTreeReaderArray<float> reco_RP_pz = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
		
		//hit locations (for debugging)
   	 	TTreeReaderArray<float> global_hit_RP_x = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.x"};
    	TTreeReaderArray<float> global_hit_RP_y = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.y"};
    	TTreeReaderArray<float> global_hit_RP_z = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.z"};
		
		//b0 tracker hits
		TTreeReaderArray<float> b0_hits_x = {tree_reader, "B0TrackerRecHits.position.x"};
    	TTreeReaderArray<float> b0_hits_y = {tree_reader, "B0TrackerRecHits.position.y"};
    	TTreeReaderArray<float> b0_hits_z = {tree_reader, "B0TrackerRecHits.position.z"};
		TTreeReaderArray<float> b0_hits_eDep = {tree_reader, "B0TrackerRecHits.edep"}; //deposited energy per hit
		
		//b0 EMCAL
		TTreeReaderArray<float> b0_cluster_x = {tree_reader, "B0ECalClusters.position.x"};
    	TTreeReaderArray<float> b0_cluster_y = {tree_reader, "B0ECalClusters.position.y"};
    	TTreeReaderArray<float> b0_cluster_z = {tree_reader, "B0ECalClusters.position.z"};
		TTreeReaderArray<float>  b0_cluster_energy = {tree_reader, "B0ECalClusters.energy"}; //deposited energy in cluster
		
		//reco tracks (where b0 tracks live!!!)
		TTreeReaderArray<float> reco_track_x = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    	TTreeReaderArray<float> reco_track_y = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    	TTreeReaderArray<float> reco_track_z = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
		
		//ZDC EMCAL
		TTreeReaderArray<float> zdc_ecal_cluster_x = {tree_reader, "ZDCEcalClusters.position.x"};
    	TTreeReaderArray<float> zdc_ecal_cluster_y = {tree_reader, "ZDCEcalClusters.position.y"};
    	TTreeReaderArray<float> zdc_ecal_cluster_z = {tree_reader, "ZDCEcalClusters.position.z"};
		TTreeReaderArray<float>  zdc_ecal_cluster_energy = {tree_reader, "ZDCEcalClusters.energy"}; //deposited energy in cluster
		

		cout << "file has " << evtTree->GetEntries() <<  " events..." << endl;

		tree_reader.SetEntriesRange(0, evtTree->GetEntries());
		while (tree_reader.Next()) {

			cout << "Reading event: " << iEvent << endl;

	    	//MCParticles
	        //finding the far-forward proton;
	    	//TLorentzVector scatMC(0,0,0,0);
			TVector3 mctrk;
			TVector3 rptrk;
			
			double maxPt=-99.;
	    	for(int imc=0;imc<mc_px_array.GetSize();imc++){
	    		mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);	
				
	    		if(mc_pdg_array[imc] == 2212){ //only checking for protons here -- change as desired
	    			
					mctrk.RotateY(0.025);
					
					h_eta_MC->Fill(mctrk.Eta());
					
					h_px_MC->Fill(mctrk.Px());
					h_py_MC->Fill(mctrk.Py());
					h_pt_MC->Fill(mctrk.Perp());
					h_pz_MC->Fill(mctrk.Pz());
				}
				
	    	}			
	    	
			//roman pots reco tracks
			for(int iRPPart = 0; iRPPart < reco_RP_px.GetSize(); iRPPart++){
    		
				TVector3 prec_romanpots(reco_RP_px[iRPPart], reco_RP_py[iRPPart], reco_RP_pz[iRPPart]);	
    		
				//prec_romanpots.RotateY(-0.025);
			
				h_px_RomanPots->Fill(prec_romanpots.Px());
				h_py_RomanPots->Fill(prec_romanpots.Py());
				h_pt_RomanPots->Fill(prec_romanpots.Perp());
				h_pz_RomanPots->Fill(prec_romanpots.Pz());
				
				h_rp_occupancy_map->Fill(global_hit_RP_x[iRPPart], global_hit_RP_y[iRPPart]);
			}
				
			
			double hit_x = -9999.;
			double hit_y = -9999.;
			double hit_z = -9999.;
			double hit_deposited_energy = -9999.;
			
			for(int b0cluster = 0; b0cluster < b0_cluster_x.GetSize(); b0cluster++){
			
				hit_x = b0_cluster_x[b0cluster];
				hit_y = b0_cluster_y[b0cluster];
				hit_z = b0_cluster_z[b0cluster];
				hit_deposited_energy = b0_cluster_energy[b0cluster]*1.246; //poor man's calibration constant, for now
							
				h_B0_emcal_occupancy_map->Fill(hit_x, hit_y);
				h_B0_emcal_cluster_energy->Fill(hit_deposited_energy);
				
			}


			//b0 tracker hits -- for debugging or external tracking
			for(int b0hit = 0; b0hit < b0_hits_x.GetSize(); b0hit++){
    		
				hit_x = b0_hits_x[b0hit];
				hit_y = b0_hits_y[b0hit];
				hit_z = b0_hits_z[b0hit];
				hit_deposited_energy = b0_hits_eDep[b0hit]*1e6; //convert GeV --> keV
				
				h_B0_hit_energy_deposit->Fill(hit_deposited_energy);
				
				if(hit_deposited_energy < 10.0){ continue; } //threshold value -- 10 keV, arbitrary for now
				
    			if(hit_z > 5800 && hit_z < 6000){ h_B0_occupancy_map_layer_0->Fill(hit_x, hit_y); }
				if(hit_z > 6000 && hit_z < 6200){ h_B0_occupancy_map_layer_1->Fill(hit_x, hit_y); }
				if(hit_z > 6200 && hit_z < 6400){ h_B0_occupancy_map_layer_2->Fill(hit_x, hit_y); }
				if(hit_z > 6400 && hit_z < 6600){ h_B0_occupancy_map_layer_3->Fill(hit_x, hit_y); }
				
			}
		
			
			//reconstructed tracks with ACTS -- used for B0
			for(int iRecoTrk = 0; iRecoTrk < reco_track_x.GetSize(); iRecoTrk++){
    		
				TVector3 prec_reco_tracks(reco_track_x[iRecoTrk], reco_track_y[iRecoTrk], reco_track_z[iRecoTrk]);	
    					
				h_px_reco_track->Fill(prec_reco_tracks.Px());
				h_py_reco_track->Fill(prec_reco_tracks.Py());
				h_pt_reco_track->Fill(prec_reco_tracks.Perp());
				h_pz_reco_track->Fill(prec_reco_tracks.Pz());
			}
				
			for(int iZdcEMCALcluster = 0; iZdcEMCALcluster < zdc_ecal_cluster_x.GetSize(); iZdcEMCALcluster++){
			
				hit_x = zdc_ecal_cluster_x[iZdcEMCALcluster];
				hit_y = zdc_ecal_cluster_y[iZdcEMCALcluster];
				hit_z = zdc_ecal_cluster_z[iZdcEMCALcluster];
				
				hit_deposited_energy = zdc_ecal_cluster_energy[iZdcEMCALcluster]*1.246; //poor man's calibration constant, for now
							
				h_ZDC_emcal_occupancy_map->Fill(hit_x, hit_y);
				h_ZDC_emcal_cluster_energy->Fill(hit_deposited_energy);
				
			}
				
			iEvent++;
		}// event loop
		
		inputRootFile->Close();
		
	}// input file loop
		
	cout << "Check integrals: " << endl;
	cout << "pt_mc integral = " << h_pt_MC->Integral() << endl;
	cout << "pt_RP_reco integral = " << h_pt_RomanPots->Integral() << endl;

	TFile * outputFile = new TFile(outputFileName, "RECREATE");

	h_px_MC->Write();
	h_py_MC->Write();
	h_pt_MC->Write();
	h_pz_MC->Write();
	
	h_px_RomanPots->Write();
	h_py_RomanPots->Write();
	h_pt_RomanPots->Write();
	h_pz_RomanPots->Write();
	h_rp_occupancy_map->Write();
	
	h_B0_occupancy_map_layer_0->Write();
	h_B0_occupancy_map_layer_1->Write();
	h_B0_occupancy_map_layer_2->Write();
	h_B0_occupancy_map_layer_3->Write();
	h_B0_hit_energy_deposit->Write();
	
	h_B0_emcal_occupancy_map->Write();
	h_B0_emcal_cluster_energy->Write();

	h_px_reco_track->Write();
	h_py_reco_track->Write();
	h_pt_reco_track->Write();
	h_pz_reco_track->Write();

	h_ZDC_emcal_occupancy_map->Write();
	h_ZDC_emcal_cluster_energy->Write();

	outputFile->Close();

	

    return;

}

