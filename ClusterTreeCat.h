// -*- c++ -*-
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 26 15:59:48 2020 by ROOT version 6.19/01
// from TTree ClusterTree/ClusterTree
// found on file: results.root
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// File modified by: Peter Loch <loch@physics.arizona.edu>
// File modified on: July 15, 2021
//////////////////////////////////////////////////////////

#ifndef ClusterTree_h
#define ClusterTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>

#include <iostream>
#include <cstdio>
#include <string>

// Header file for the classes stored in the TTree if any.

class ClusterTree {
    public :
        TTree*          fChain   = { (TTree*)0 };   //! pointer to the analyzed TTree (input and modified output)
	TTree*          EM_tree  = { (TTree*)0 };   //! pointer to output EM tree  
	TTree*          Had_tree = { (TTree*)0 };   //! pointer to outout Had tree
        Int_t           fCurrent = { -1 }       ;   //!current Tree number in a TChain

        // -- DO NOT CHANGE THE VARIABLE TYPES!
	//
	//    Long64_t -> ROOT type definition for "long long int" (64 bits/8 bytes) [large range integer number] 
	//    Double_t -> ROOT type definition for "double"        (64 bits/8 bytes) [high precision floating point number]
        //
        // ----------------------------------------------------------
        // Data structure from results.root 07.15.2021 14:10 CEST
        //
	// Long64_t        entry;
	// Long64_t        runNumber;
	// Long64_t        eventNumber;
	// Double_t        truthE;
	// Double_t        truthPt;
	// Double_t        truthEta;
	// Double_t        truthPhi;
	// Long64_t        truthPDG;
	// Long64_t        nCluster;
	// Long64_t        clusterIndex;
	// Long64_t        cluster_nCells;
	// Long64_t        cluster_nCells_tot;
	// Double_t        clusterECalib;
	// Double_t        clusterPtCalib;
	// Double_t        clusterEtaCalib;
	// Double_t        clusterPhiCalib;
	// Double_t        cluster_sumCellECalib;
	// Double_t        cluster_fracECalib;
	// Double_t        cluster_fracECalib_ref;
	// Double_t        clusterE;
	// Double_t        clusterPt;
	// Double_t        clusterEta;
	// Double_t        clusterPhi;
	// Double_t        cluster_sumCellE;
	// Double_t        cluster_time;
	// Double_t        cluster_fracE;
	// Double_t        cluster_fracE_ref;
	// Double_t        cluster_EM_PROBABILITY;
	// Double_t        cluster_HAD_WEIGHT;
	// Double_t        cluster_OOC_WEIGHT;
	// Double_t        cluster_DM_WEIGHT;
	// Double_t        cluster_ENG_CALIB_TOT;
	// Double_t        cluster_ENG_CALIB_OUT_T;
	// Double_t        cluster_ENG_CALIB_DEAD_TOT;
	// Double_t        cluster_CENTER_MAG;
	// Double_t        cluster_FIRST_ENG_DENS;
	// Double_t        cluster_FIRST_PHI;
	// Double_t        cluster_FIRST_ETA;
	// Double_t        cluster_SECOND_R;
	// Double_t        cluster_SECOND_LAMBDA;
	// Double_t        cluster_DELTA_PHI;
	// Double_t        cluster_DELTA_THETA;
	// Double_t        cluster_DELTA_ALPHA;
	// Double_t        cluster_CENTER_X;
	// Double_t        cluster_CENTER_Y;
	// Double_t        cluster_CENTER_Z;
	// Double_t        cluster_CENTER_LAMBDA;
	// Double_t        cluster_LATERAL;
	// Double_t        cluster_LONGITUDINAL;
	// Double_t        cluster_ENG_FRAC_EM;
	// Double_t        cluster_ENG_FRAC_MAX;
	// Double_t        cluster_ENG_FRAC_CORE;
	// Double_t        cluster_SECOND_ENG_DENS;
	// Double_t        cluster_ISOLATION;
	// Double_t        cluster_ENG_BAD_CELLS;
	// Double_t        cluster_N_BAD_CELLS;
	// Double_t        cluster_N_BAD_CELLS_CORR;
	// Double_t        cluster_BAD_CELLS_CORR_E;
	// Double_t        cluster_BADLARQ_FRAC;
	// Double_t        cluster_ENG_POS;
	// Double_t        cluster_SIGNIFICANCE;
	// Double_t        cluster_CELL_SIGNIFICANCE;
	// Double_t        cluster_CELL_SIG_SAMPLING;
	// Double_t        cluster_AVG_LAR_Q;
	// Double_t        cluster_AVG_TILE_Q;
	// Double_t        cluster_ENG_BAD_HV_CELLS;
	// Double_t        cluster_N_BAD_HV_CELLS;
	// Double_t        cluster_PTD;
	// Double_t        cluster_MASS;
	// Double_t        cluster_SECOND_TIME;
	// Double_t        CalibratedE;                    
        // ---------------------------------------------------

	///////////////////////////////////////
        // Declaration of variables (leaves) //
	///////////////////////////////////////

	// -- book keeping
        Long64_t        runNumber;                   // run number (integral number)
        Long64_t        eventNumber;                 // event number (integral number)
	// -- truth particle kinematics and idenfication
        Double_t        truthE;                      // energy
        Double_t        truthPt;                     // transverse momentum
        Double_t        truthEta;                    // pseudo-rapidity
        Double_t        truthPhi;                    // azimuth
        Long64_t        truthPDG;                    // Particle Data Group indentifier for particle type
	// -- topo-cluster information
        Long64_t        nCluster;                    // number of clusters per particle 
        Long64_t        clusterIndex;                // index of cluster [0,...,nCluster-1]
        Long64_t        cluster_nCells;              // number of cells with E > 0 in cluster
        Long64_t        cluster_nCells_tot;          // number of all cells in cluster
	// -- topo-cluster signal
        Double_t        clusterECalib;               // energy calibrated at LCW scale
        Double_t        clusterPtCalib;              // transverse momentum calibrated at LCW scale
        Double_t        clusterEtaCalib;             // pseudo-rapidity reconstructed at LCW scale 
        Double_t        clusterPhiCalib;             // pseudo-rapidity reconstructed at LCW scale
        Double_t        cluster_sumCellECalib;       // sum of LCW-weighted cell energies  
	Double_t        cluster_fracECalib;          // NEW: LCW energy of cluster divided by the sum of LCW energies of all clusters associated with this particle
	Double_t        cluster_fracECalib_ref;      // NEW: LCW energy of cluster divided by the truth particle energy         
        Double_t        clusterE;                    // energy at EM scale
        Double_t        clusterPt;                   // transverse momentum calibrated at EM scale
        Double_t        clusterEta;                  // pseudo-rapidity reconstructed at EM scale
        Double_t        clusterPhi;		     // pseudo-rapidity reconstructed at EM scale
        Double_t        cluster_sumCellE;	     // sum of geometrically weighted cell energies
      	Double_t        cluster_fracE;               // NEW: EM energy of cluster divided by the sum of EM energies of all clusters associated with this particle
	Double_t        cluster_fracE_ref;           // NEW: EM energy of cluster divided by the truth particle energy         
	// -- topo-cluster properties                      
	Double_t        cluster_time;                // first moment of cell-energy-squared-weighted cell timing distribution
        Double_t        cluster_EM_PROBABILITY;      // EM probability from LCW calibration
	// -- topo-cluster calibration weights
        Double_t        cluster_HAD_WEIGHT;          // LCW hadronic calibration weight
        Double_t        cluster_OOC_WEIGHT;          // LCW out-of-cluster calibration weight
        Double_t        cluster_DM_WEIGHT;           // LCW dead material calibration weight
	// -- topo-cluster calibration hit energies
        Double_t        cluster_ENG_CALIB_TOT;       // energy deposited in cells in cluster
        Double_t        cluster_ENG_CALIB_OUT_T;     // energy deposited outside of cluster but assigned to it
        Double_t        cluster_ENG_CALIB_DEAD_TOT;  // energy in dead material around cluster
	// -- topo-cluster moments: shapes and location
        Double_t        cluster_CENTER_MAG;          // cluster distance to vertex
        Double_t        cluster_FIRST_ENG_DENS;      // first moment of cell energy densities in cluster
        Double_t        cluster_FIRST_PHI;           // first azimuthal moment
        Double_t        cluster_FIRST_ETA;           // first pseudo-rapidity moment
        Double_t        cluster_SECOND_R;            // second moment of radial cell distances from cluster center-of-gravity
        Double_t        cluster_SECOND_LAMBDA;       // second moment of longitudinal cell distances form cluster center-of-gravity 
        Double_t        cluster_DELTA_PHI;           // azimuthal distance between nominal cluster direction and principal cluster axis
        Double_t        cluster_DELTA_THETA;         // polar angle distance between nominal cluster direction and principal cluster axis
        Double_t        cluster_DELTA_ALPHA;         // angular distance between nominal cluster direction and principal cluster axis
        Double_t        cluster_CENTER_X;            // X coordinate of cluster center-of-gravity
        Double_t        cluster_CENTER_Y;            // Y coordinate of cluster center-of-gravity
        Double_t        cluster_CENTER_Z;            // Z coordinate of cluster center-of-gravity
        Double_t        cluster_CENTER_LAMBDA;       // distance of cluster center-of-gravity from front of calorimeter
        Double_t        cluster_LATERAL;             // normalized lateral energy dispersion
        Double_t        cluster_LONGITUDINAL;        // normalized longitudinal energy dispersion
        Double_t        cluster_ENG_FRAC_EM;         // fraction of cluster energy in electromagnetic calorimeter 
        Double_t        cluster_ENG_FRAC_MAX;        // fraction of cluster energy in highest signal cell
        Double_t        cluster_ENG_FRAC_CORE;       // fraction of cluster energy in four highest energy cells
        Double_t        cluster_SECOND_ENG_DENS;     // second moment of cell energy density distribution in cluster
        Double_t        cluster_ISOLATION;           // isolation moment
	// -- topo-cluster moments: signal quality and relevance
        Double_t        cluster_ENG_BAD_CELLS;       // energy in bad cells in cluster
        Long64_t        cluster_N_BAD_CELLS;         // number of bad cells in cluster
        Long64_t        cluster_N_BAD_CELLS_CORR;    // number of corrected bad cells in cluster
        Double_t        cluster_BAD_CELLS_CORR_E;    // corrected energy in bad cells in cluster
        Long64_t        cluster_BADLARQ_FRAC;        // fraction of cells with bad LAr signal quality
        Double_t        cluster_ENG_POS;             // sum of energies of cells with E > 0 
        Double_t        cluster_SIGNIFICANCE;        // cluster signal significance
        Double_t        cluster_CELL_SIGNIFICANCE;   // significance of cell signal with largest signal-over-noise
        Long64_t        cluster_CELL_SIG_SAMPLING;   // sampling of cell signal with largest signal-over-noise
        Double_t        cluster_AVG_LAR_Q;           // average LAr cell signal quality
        Double_t        cluster_AVG_TILE_Q;          // average Tile cell signal quality
        Double_t        cluster_ENG_BAD_HV_CELLS;    // energy in cells with bad high voltage
        Long64_t        cluster_N_BAD_HV_CELLS;      // number of cells with bad high voltage
	// -- topo-cluster moments: other shower and signal characteristics
        Double_t        cluster_PTD;                 // longitudinal fragmentation function 
        Double_t        cluster_MASS;                // cluster mass at EM scale calculated from cells with E > 0
	Double_t        cluster_SECOND_TIME;         // second moment of cell timing distribution in cluster 
        Long64_t        EM_Shower;                   // *** NOT IN INPUT FILE *** indicator for EM shower
        Double_t        EM_Pro;                      // *** NOT IN INPUT FILE *** EM probability from ML  
        Double_t        CalibratedE;                 // ML calibrated energy
        Double_t        Delta_E;                     // *** NOT IN INPUT FILE *** difference between true deposited energy and ML calibrated energy
        Double_t        Delta_Calib_E;               // *** NOT IN INPUT FILE *** difference between true deposited energy and LCW calibrated energy 

	///////////////////////////////////////////////
	// Branches for input tree and appended data //
	///////////////////////////////////////////////

	// -- book keeping
        TBranch* b_runNumber;   //!
        TBranch* b_eventNumber; //!
	// -- truth particle kinematics and idenfication
        TBranch* b_truthE;   //!
        TBranch* b_truthPt;  //!
        TBranch* b_truthEta; //!
        TBranch* b_truthPhi; //!
        TBranch* b_truthPDG; //!
	// -- topo-cluster information
        TBranch* b_nCluster;           //!
        TBranch* b_clusterIndex;       //!
        TBranch* b_cluster_nCells;     //!
        TBranch* b_cluster_nCells_tot; //!
	// -- topo-cluster signal
        TBranch* b_clusterECalib;          //!
        TBranch* b_clusterPtCalib;         //!
        TBranch* b_clusterEtaCalib;        //!
        TBranch* b_clusterPhiCalib;        //!
        TBranch* b_cluster_sumCellECalib;  //!
	TBranch* b_cluster_fracECalib;     //!
	TBranch* b_cluster_fracECalib_ref; //!
        TBranch* b_clusterE;               //!
        TBranch* b_clusterPt;              //!
        TBranch* b_clusterEta;             //!
        TBranch* b_clusterPhi;	           //!
        TBranch* b_cluster_sumCellE;	   //!
      	TBranch* b_cluster_fracE;          //!
	TBranch* b_cluster_fracE_ref;      //!
	// -- topo-cluster properties                      
	TBranch* b_cluster_time;           //!
        TBranch* b_cluster_EM_PROBABILITY; //!
	// -- topo-cluster calibration weights
        TBranch* b_cluster_HAD_WEIGHT; //!
        TBranch* b_cluster_OOC_WEIGHT; //!
        TBranch* b_cluster_DM_WEIGHT;  //!
	// -- topo-cluster calibration hit energies
        TBranch* b_cluster_ENG_CALIB_TOT;      //!
        TBranch* b_cluster_ENG_CALIB_OUT_T;    //!
        TBranch* b_cluster_ENG_CALIB_DEAD_TOT; //!
	// -- topo-cluster moments: shapes and location
        TBranch* b_cluster_CENTER_MAG;      //!
        TBranch* b_cluster_FIRST_ENG_DENS;  //!
        TBranch* b_cluster_FIRST_PHI;       //!
        TBranch* b_cluster_FIRST_ETA;       //!
        TBranch* b_cluster_SECOND_R;        //!
        TBranch* b_cluster_SECOND_LAMBDA;   //!
        TBranch* b_cluster_DELTA_PHI;       //!
        TBranch* b_cluster_DELTA_THETA;     //!
        TBranch* b_cluster_DELTA_ALPHA;     //!
        TBranch* b_cluster_CENTER_X;        //!
        TBranch* b_cluster_CENTER_Y;        //!
        TBranch* b_cluster_CENTER_Z;        //!
        TBranch* b_cluster_CENTER_LAMBDA;   //!
        TBranch* b_cluster_LATERAL;         //!
        TBranch* b_cluster_LONGITUDINAL;    //!
        TBranch* b_cluster_ENG_FRAC_EM;     //! 
        TBranch* b_cluster_ENG_FRAC_MAX;    //!
        TBranch* b_cluster_ENG_FRAC_CORE;   //!
        TBranch* b_cluster_SECOND_ENG_DENS; //!
        TBranch* b_cluster_ISOLATION;       //!
	// -- topo-cluster moments: signal quality and relevance
        TBranch* b_cluster_ENG_BAD_CELLS;     //!
        TBranch* b_cluster_N_BAD_CELLS;       //!
        TBranch* b_cluster_N_BAD_CELLS_CORR;  //!
        TBranch* b_cluster_BAD_CELLS_CORR_E;  //!
        TBranch* b_cluster_BADLARQ_FRAC;      //!
        TBranch* b_cluster_ENG_POS;           //!
        TBranch* b_cluster_SIGNIFICANCE;      //!
        TBranch* b_cluster_CELL_SIGNIFICANCE; //!
        TBranch* b_cluster_CELL_SIG_SAMPLING; //!
        TBranch* b_cluster_AVG_LAR_Q;         //!
        TBranch* b_cluster_AVG_TILE_Q;        //!
        TBranch* b_cluster_ENG_BAD_HV_CELLS;  //!
        TBranch* b_cluster_N_BAD_HV_CELLS;    //!
	// -- topo-cluster moments: other shower and signal characteristics
        TBranch* b_cluster_PTD;         //!
        TBranch* b_cluster_MASS;        //!
	TBranch* b_cluster_SECOND_TIME; //!
        TBranch* b_EM_Shower;           //!
        TBranch* b_EM_Pro;              //!
        TBranch* b_CalibratedE;         //!
        TBranch* b_Delta_E;             //!
        TBranch* b_Delta_Calib_E;       //! 

	//////////////////////////
	// Branches for EM_tree //
	//////////////////////////

	// -- book keeping
        TBranch* EM_b_runNumber;   //!
        TBranch* EM_b_eventNumber; //!
	// -- truth particle kinematics and idenfication
        TBranch* EM_b_truthE;   //!
        TBranch* EM_b_truthPt;  //!
        TBranch* EM_b_truthEta; //!
        TBranch* EM_b_truthPhi; //!
        TBranch* EM_b_truthPDG; //!
	// -- topo-cluster information
        TBranch* EM_b_nCluster;           //!
        TBranch* EM_b_clusterIndex;       //!
        TBranch* EM_b_cluster_nCells;     //!
        TBranch* EM_b_cluster_nCells_tot; //!
	// -- topo-cluster signal
        TBranch* EM_b_clusterECalib;          //!
        TBranch* EM_b_clusterPtCalib;         //!
        TBranch* EM_b_clusterEtaCalib;        //!
        TBranch* EM_b_clusterPhiCalib;        //!
        TBranch* EM_b_cluster_sumCellECalib;  //!
	TBranch* EM_b_cluster_fracECalib;     //!
	TBranch* EM_b_cluster_fracECalib_ref; //!
        TBranch* EM_b_clusterE;               //!
        TBranch* EM_b_clusterPt;              //!
        TBranch* EM_b_clusterEta;             //!
        TBranch* EM_b_clusterPhi;	           //!
        TBranch* EM_b_cluster_sumCellE;	   //!
      	TBranch* EM_b_cluster_fracE;          //!
	TBranch* EM_b_cluster_fracE_ref;      //!
	// -- topo-cluster properties                      
	TBranch* EM_b_cluster_time;           //!
        TBranch* EM_b_cluster_EM_PROBABILITY; //!
	// -- topo-cluster calibration weights
        TBranch* EM_b_cluster_HAD_WEIGHT; //!
        TBranch* EM_b_cluster_OOC_WEIGHT; //!
        TBranch* EM_b_cluster_DM_WEIGHT;  //!
	// -- topo-cluster calibration hit energies
        TBranch* EM_b_cluster_ENG_CALIB_TOT;      //!
        TBranch* EM_b_cluster_ENG_CALIB_OUT_T;    //!
        TBranch* EM_b_cluster_ENG_CALIB_DEAD_TOT; //!
	// -- topo-cluster moments: shapes and location
        TBranch* EM_b_cluster_CENTER_MAG;      //!
        TBranch* EM_b_cluster_FIRST_ENG_DENS;  //!
        TBranch* EM_b_cluster_FIRST_PHI;       //!
        TBranch* EM_b_cluster_FIRST_ETA;       //!
        TBranch* EM_b_cluster_SECOND_R;        //!
        TBranch* EM_b_cluster_SECOND_LAMBDA;   //!
        TBranch* EM_b_cluster_DELTA_PHI;       //!
        TBranch* EM_b_cluster_DELTA_THETA;     //!
        TBranch* EM_b_cluster_DELTA_ALPHA;     //!
        TBranch* EM_b_cluster_CENTER_X;        //!
        TBranch* EM_b_cluster_CENTER_Y;        //!
        TBranch* EM_b_cluster_CENTER_Z;        //!
        TBranch* EM_b_cluster_CENTER_LAMBDA;   //!
        TBranch* EM_b_cluster_LATERAL;         //!
        TBranch* EM_b_cluster_LONGITUDINAL;    //!
        TBranch* EM_b_cluster_ENG_FRAC_EM;     //! 
        TBranch* EM_b_cluster_ENG_FRAC_MAX;    //!
        TBranch* EM_b_cluster_ENG_FRAC_CORE;   //!
        TBranch* EM_b_cluster_SECOND_ENG_DENS; //!
        TBranch* EM_b_cluster_ISOLATION;       //!
	// -- topo-cluster moments: signal quality and relevance
        TBranch* EM_b_cluster_ENG_BAD_CELLS;     //!
        TBranch* EM_b_cluster_N_BAD_CELLS;       //!
        TBranch* EM_b_cluster_N_BAD_CELLS_CORR;  //!
        TBranch* EM_b_cluster_BAD_CELLS_CORR_E;  //!
        TBranch* EM_b_cluster_BADLARQ_FRAC;      //!
        TBranch* EM_b_cluster_ENG_POS;           //!
        TBranch* EM_b_cluster_SIGNIFICANCE;      //!
        TBranch* EM_b_cluster_CELL_SIGNIFICANCE; //!
        TBranch* EM_b_cluster_CELL_SIG_SAMPLING; //!
        TBranch* EM_b_cluster_AVG_LAR_Q;         //!
        TBranch* EM_b_cluster_AVG_TILE_Q;        //!
        TBranch* EM_b_cluster_ENG_BAD_HV_CELLS;  //!
        TBranch* EM_b_cluster_N_BAD_HV_CELLS;    //!
	// -- topo-cluster moments: other shower and signal characteristics
        TBranch* EM_b_cluster_PTD;         //!
        TBranch* EM_b_cluster_MASS;        //!
	TBranch* EM_b_cluster_SECOND_TIME; //!
        TBranch* EM_b_EM_Shower;           //!
        TBranch* EM_b_EM_Pro;              //!
        TBranch* EM_b_CalibratedE;         //!
        TBranch* EM_b_Delta_E;             //!
        TBranch* EM_b_Delta_Calib_E;       //! 

	///////////////////////////
	// Branches for Had_tree //
	///////////////////////////

	// -- book keeping
        TBranch* Had_b_runNumber;   //!
        TBranch* Had_b_eventNumber; //!
	// -- truth particle kinematics and idenfication
        TBranch* Had_b_truthE;   //!
        TBranch* Had_b_truthPt;  //!
        TBranch* Had_b_truthEta; //!
        TBranch* Had_b_truthPhi; //!
        TBranch* Had_b_truthPDG; //!
	// -- topo-cluster information
        TBranch* Had_b_nCluster;           //!
        TBranch* Had_b_clusterIndex;       //!
        TBranch* Had_b_cluster_nCells;     //!
        TBranch* Had_b_cluster_nCells_tot; //!
	// -- topo-cluster signal
        TBranch* Had_b_clusterECalib;          //!
        TBranch* Had_b_clusterPtCalib;         //!
        TBranch* Had_b_clusterEtaCalib;        //!
        TBranch* Had_b_clusterPhiCalib;        //!
        TBranch* Had_b_cluster_sumCellECalib;  //!
	TBranch* Had_b_cluster_fracECalib;     //!
	TBranch* Had_b_cluster_fracECalib_ref; //!
        TBranch* Had_b_clusterE;               //!
        TBranch* Had_b_clusterPt;              //!
        TBranch* Had_b_clusterEta;             //!
        TBranch* Had_b_clusterPhi;             //!
        TBranch* Had_b_cluster_sumCellE;       //!
      	TBranch* Had_b_cluster_fracE;          //!
	TBranch* Had_b_cluster_fracE_ref;      //!
	// -- topo-cluster properties                      
	TBranch* Had_b_cluster_time;           //!
        TBranch* Had_b_cluster_EM_PROBABILITY; //!
	// -- topo-cluster calibration weights
        TBranch* Had_b_cluster_HAD_WEIGHT; //!
        TBranch* Had_b_cluster_OOC_WEIGHT; //!
        TBranch* Had_b_cluster_DM_WEIGHT;  //!
	// -- topo-cluster calibration hit energies
        TBranch* Had_b_cluster_ENG_CALIB_TOT;      //!
        TBranch* Had_b_cluster_ENG_CALIB_OUT_T;    //!
        TBranch* Had_b_cluster_ENG_CALIB_DEAD_TOT; //!
	// -- topo-cluster moments: shapes and location
        TBranch* Had_b_cluster_CENTER_MAG;      //!
        TBranch* Had_b_cluster_FIRST_ENG_DENS;  //!
        TBranch* Had_b_cluster_FIRST_PHI;       //!
        TBranch* Had_b_cluster_FIRST_ETA;       //!
        TBranch* Had_b_cluster_SECOND_R;        //!
        TBranch* Had_b_cluster_SECOND_LAMBDA;   //!
        TBranch* Had_b_cluster_DELTA_PHI;       //!
        TBranch* Had_b_cluster_DELTA_THETA;     //!
        TBranch* Had_b_cluster_DELTA_ALPHA;     //!
        TBranch* Had_b_cluster_CENTER_X;        //!
        TBranch* Had_b_cluster_CENTER_Y;        //!
        TBranch* Had_b_cluster_CENTER_Z;        //!
        TBranch* Had_b_cluster_CENTER_LAMBDA;   //!
        TBranch* Had_b_cluster_LATERAL;         //!
        TBranch* Had_b_cluster_LONGITUDINAL;    //!
        TBranch* Had_b_cluster_ENG_FRAC_EM;     //! 
        TBranch* Had_b_cluster_ENG_FRAC_MAX;    //!
        TBranch* Had_b_cluster_ENG_FRAC_CORE;   //!
        TBranch* Had_b_cluster_SECOND_ENG_DENS; //!
        TBranch* Had_b_cluster_ISOLATION;       //!
	// -- topo-cluster moments: signal quality and relevance
        TBranch* Had_b_cluster_ENG_BAD_CELLS;     //!
        TBranch* Had_b_cluster_N_BAD_CELLS;       //!
        TBranch* Had_b_cluster_N_BAD_CELLS_CORR;  //!
        TBranch* Had_b_cluster_BAD_CELLS_CORR_E;  //!
        TBranch* Had_b_cluster_BADLARQ_FRAC;      //!
        TBranch* Had_b_cluster_ENG_POS;           //!
        TBranch* Had_b_cluster_SIGNIFICANCE;      //!
        TBranch* Had_b_cluster_CELL_SIGNIFICANCE; //!
        TBranch* Had_b_cluster_CELL_SIG_SAMPLING; //!
        TBranch* Had_b_cluster_AVG_LAR_Q;         //!
        TBranch* Had_b_cluster_AVG_TILE_Q;        //!
        TBranch* Had_b_cluster_ENG_BAD_HV_CELLS;  //!
        TBranch* Had_b_cluster_N_BAD_HV_CELLS;    //!
	// -- topo-cluster moments: other shower and signal characteristics
        TBranch* Had_b_cluster_PTD;         //!
        TBranch* Had_b_cluster_MASS;        //!
	TBranch* Had_b_cluster_SECOND_TIME; //!
        TBranch* Had_b_EM_Shower;           //!
        TBranch* Had_b_EM_Pro;              //!
        TBranch* Had_b_CalibratedE;         //!
        TBranch* Had_b_Delta_E;             //!
        TBranch* Had_b_Delta_Calib_E;       //! 

        ClusterTree(TTree *tree=0);
        virtual ~ClusterTree();
        virtual Int_t    Cut(Long64_t entry);
        virtual Int_t    GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Init(TTree *tree);
        virtual void     Init_EM(TTree *tree);
        virtual void     Init_Had(TTree *tree);
        virtual void     Loop(const std::string& outFileName="ml_results.root");
        virtual Bool_t   Notify();
        virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ClusterTree_cxx
ClusterTree::ClusterTree(TTree *tree) : fChain(0)
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("results.root");
        }
        f->GetObject("ClusterTree",tree);
    }
    // initialize input tree
    Init(tree);
    // initialize output trees -> now in Loop()!!!!!
    // Init_Had(new TTree("Had_tree", "Hadron shower Tree")); 
    // Init_EM (new TTree("EM_tree",   "EM shower Tree"   ));
}

ClusterTree::~ClusterTree()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t ClusterTree::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t ClusterTree::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void ClusterTree::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;

    fChain = tree;
    fChain->SetMakeClass(1);

    TObjArray* bPtrs = fChain->GetListOfBranches(); 
    if ( bPtrs == nullptr ) { std::cout << "Tree <" << fChain->GetName() << "> has no branches, initialization failed" << std::endl; return; } 
    // -- branch addresses are only set for branches that exist on input
    if ( bPtrs->FindObject("runNumber"                  ) != nullptr ) { fChain->SetBranchAddress("runNumber",                 &runNumber,                  &b_runNumber);                 }
    if ( bPtrs->FindObject("eventNumber"                ) != nullptr ) { fChain->SetBranchAddress("eventNumber",               &eventNumber,                &b_eventNumber);               }
    if ( bPtrs->FindObject("truthE"                     ) != nullptr ) { fChain->SetBranchAddress("truthE",                    &truthE,                     &b_truthE);                    }
    if ( bPtrs->FindObject("truthPt"                    ) != nullptr ) { fChain->SetBranchAddress("truthPt",                   &truthPt,                    &b_truthPt);                   }
    if ( bPtrs->FindObject("truthEta"                   ) != nullptr ) { fChain->SetBranchAddress("truthEta",                  &truthEta,                   &b_truthEta);                  }
    if ( bPtrs->FindObject("truthPhi"                   ) != nullptr ) { fChain->SetBranchAddress("truthPhi",                  &truthPhi,                   &b_truthPhi);                  }
    if ( bPtrs->FindObject("truthPDG"                   ) != nullptr ) { fChain->SetBranchAddress("truthPDG",                  &truthPDG,                   &b_truthPDG);                  }
    if ( bPtrs->FindObject("nCluster"                   ) != nullptr ) { fChain->SetBranchAddress("nCluster",                  &nCluster,                   &b_nCluster);                  }
    if ( bPtrs->FindObject("clusterIndex"               ) != nullptr ) { fChain->SetBranchAddress("clusterIndex",              &clusterIndex,               &b_clusterIndex);              }
    if ( bPtrs->FindObject("cluster_nCells"             ) != nullptr ) { fChain->SetBranchAddress("cluster_nCells",            &cluster_nCells,             &b_cluster_nCells);            }
    if ( bPtrs->FindObject("cluster_nCells_tot"         ) != nullptr ) { fChain->SetBranchAddress("cluster_nCells_tot",        &cluster_nCells_tot,         &b_cluster_nCells_tot);        }
    if ( bPtrs->FindObject("clusterECalib"              ) != nullptr ) { fChain->SetBranchAddress("clusterECalib",             &clusterECalib,              &b_clusterECalib);             }
    if ( bPtrs->FindObject("clusterPtCalib"             ) != nullptr ) { fChain->SetBranchAddress("clusterPtCalib",            &clusterPtCalib,             &b_clusterPtCalib);            }
    if ( bPtrs->FindObject("clusterEtaCalib"            ) != nullptr ) { fChain->SetBranchAddress("clusterEtaCalib",           &clusterEtaCalib,            &b_clusterEtaCalib);           }
    if ( bPtrs->FindObject("clusterPhiCalib"            ) != nullptr ) { fChain->SetBranchAddress("clusterPhiCalib",           &clusterPhiCalib,            &b_clusterPhiCalib);           }
    if ( bPtrs->FindObject("cluster_sumCellECalib"      ) != nullptr ) { fChain->SetBranchAddress("cluster_sumCellECalib",     &cluster_sumCellECalib,      &b_cluster_sumCellECalib);     }
    if ( bPtrs->FindObject("cluster_fracECalib"         ) != nullptr ) { fChain->SetBranchAddress("cluster_fracECalib",        &cluster_fracECalib,         &b_cluster_fracECalib);        }
    if ( bPtrs->FindObject("cluster_fracECalib_ref"     ) != nullptr ) { fChain->SetBranchAddress("cluster_fracECalib_ref",    &cluster_fracECalib_ref,     &b_cluster_fracECalib_ref);    }
    if ( bPtrs->FindObject("clusterE"                   ) != nullptr ) { fChain->SetBranchAddress("clusterE",                  &clusterE,                   &b_clusterE);                  }
    if ( bPtrs->FindObject("clusterPt"                  ) != nullptr ) { fChain->SetBranchAddress("clusterPt",                 &clusterPt,                  &b_clusterPt);                 }
    if ( bPtrs->FindObject("clusterEta"                 ) != nullptr ) { fChain->SetBranchAddress("clusterEta",                &clusterEta,                 &b_clusterEta);                }
    if ( bPtrs->FindObject("clusterPhi"                 ) != nullptr ) { fChain->SetBranchAddress("clusterPhi",                &clusterPhi,                 &b_clusterPhi);                }
    if ( bPtrs->FindObject("cluster_sumCellE"           ) != nullptr ) { fChain->SetBranchAddress("cluster_sumCellE",          &cluster_sumCellE,           &b_cluster_sumCellE);          }
    if ( bPtrs->FindObject("cluster_fracE"              ) != nullptr ) { fChain->SetBranchAddress("cluster_fracE",             &cluster_fracE,              &b_cluster_fracE);             }
    if ( bPtrs->FindObject("cluster_fracE_ref"          ) != nullptr ) { fChain->SetBranchAddress("cluster_fracE_ref",         &cluster_fracE_ref,          &b_cluster_fracE_ref);         }
    if ( bPtrs->FindObject("cluster_EM_PROBABILITY"     ) != nullptr ) { fChain->SetBranchAddress("cluster_EM_PROBABILITY",    &cluster_EM_PROBABILITY,     &b_cluster_EM_PROBABILITY);    }
    if ( bPtrs->FindObject("cluster_time"               ) != nullptr ) { fChain->SetBranchAddress("cluster_time",              &cluster_time,               &b_cluster_time);              }
    if ( bPtrs->FindObject("cluster_HAD_WEIGHT"         ) != nullptr ) { fChain->SetBranchAddress("cluster_HAD_WEIGHT",        &cluster_HAD_WEIGHT,         &b_cluster_HAD_WEIGHT);        }
    if ( bPtrs->FindObject("cluster_OOC_WEIGHT"         ) != nullptr ) { fChain->SetBranchAddress("cluster_OOC_WEIGHT",        &cluster_OOC_WEIGHT,         &b_cluster_OOC_WEIGHT);        }
    if ( bPtrs->FindObject("cluster_DM_WEIGHT"          ) != nullptr ) { fChain->SetBranchAddress("cluster_DM_WEIGHT",         &cluster_DM_WEIGHT,          &b_cluster_DM_WEIGHT);         }
    if ( bPtrs->FindObject("cluster_ENG_CALIB_TOT"      ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_CALIB_TOT",     &cluster_ENG_CALIB_TOT,      &b_cluster_ENG_CALIB_TOT);     }
    if ( bPtrs->FindObject("cluster_ENG_CALIB_OUT_T"    ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_CALIB_OUT_T",   &cluster_ENG_CALIB_OUT_T,    &b_cluster_ENG_CALIB_OUT_T);   }
    if ( bPtrs->FindObject("cluster_ENG_CALIB_DEAD_TOT" ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_CALIB_DEAD_TOT",&cluster_ENG_CALIB_DEAD_TOT, &b_cluster_ENG_CALIB_DEAD_TOT);}
    if ( bPtrs->FindObject("cluster_CENTER_MAG"         ) != nullptr ) { fChain->SetBranchAddress("cluster_CENTER_MAG",        &cluster_CENTER_MAG,         &b_cluster_CENTER_MAG);        }
    if ( bPtrs->FindObject("cluster_FIRST_ENG_DENS"     ) != nullptr ) { fChain->SetBranchAddress("cluster_FIRST_ENG_DENS",    &cluster_FIRST_ENG_DENS,     &b_cluster_FIRST_ENG_DENS);    }
    if ( bPtrs->FindObject("cluster_FIRST_PHI"          ) != nullptr ) { fChain->SetBranchAddress("cluster_FIRST_PHI",         &cluster_FIRST_PHI,          &b_cluster_FIRST_PHI);         }
    if ( bPtrs->FindObject("cluster_FIRST_ETA"          ) != nullptr ) { fChain->SetBranchAddress("cluster_FIRST_ETA",         &cluster_FIRST_ETA,          &b_cluster_FIRST_ETA);         }
    if ( bPtrs->FindObject("cluster_SECOND_R"           ) != nullptr ) { fChain->SetBranchAddress("cluster_SECOND_R",          &cluster_SECOND_R,           &b_cluster_SECOND_R);          }
    if ( bPtrs->FindObject("cluster_SECOND_LAMBDA"      ) != nullptr ) { fChain->SetBranchAddress("cluster_SECOND_LAMBDA",     &cluster_SECOND_LAMBDA,      &b_cluster_SECOND_LAMBDA);     }
    if ( bPtrs->FindObject("cluster_DELTA_PHI"          ) != nullptr ) { fChain->SetBranchAddress("cluster_DELTA_PHI",         &cluster_DELTA_PHI,          &b_cluster_DELTA_PHI);         }
    if ( bPtrs->FindObject("cluster_DELTA_THETA"        ) != nullptr ) { fChain->SetBranchAddress("cluster_DELTA_THETA",       &cluster_DELTA_THETA,        &b_cluster_DELTA_THETA);       }
    if ( bPtrs->FindObject("cluster_DELTA_ALPHA"        ) != nullptr ) { fChain->SetBranchAddress("cluster_DELTA_ALPHA",       &cluster_DELTA_ALPHA,        &b_cluster_DELTA_ALPHA);       }
    if ( bPtrs->FindObject("cluster_CENTER_X"           ) != nullptr ) { fChain->SetBranchAddress("cluster_CENTER_X",          &cluster_CENTER_X,           &b_cluster_CENTER_X);          }
    if ( bPtrs->FindObject("cluster_CENTER_Y"           ) != nullptr ) { fChain->SetBranchAddress("cluster_CENTER_Y",          &cluster_CENTER_Y,           &b_cluster_CENTER_Y);          }
    if ( bPtrs->FindObject("cluster_CENTER_Z"           ) != nullptr ) { fChain->SetBranchAddress("cluster_CENTER_Z",          &cluster_CENTER_Z,           &b_cluster_CENTER_Z);          }
    if ( bPtrs->FindObject("cluster_CENTER_LAMBDA"      ) != nullptr ) { fChain->SetBranchAddress("cluster_CENTER_LAMBDA",     &cluster_CENTER_LAMBDA,      &b_cluster_CENTER_LAMBDA);     }
    if ( bPtrs->FindObject("cluster_LATERAL"            ) != nullptr ) { fChain->SetBranchAddress("cluster_LATERAL",           &cluster_LATERAL,            &b_cluster_LATERAL);           }
    if ( bPtrs->FindObject("cluster_LONGITUDINAL"       ) != nullptr ) { fChain->SetBranchAddress("cluster_LONGITUDINAL",      &cluster_LONGITUDINAL,       &b_cluster_LONGITUDINAL);      }
    if ( bPtrs->FindObject("cluster_ENG_FRAC_EM"        ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_FRAC_EM",       &cluster_ENG_FRAC_EM,        &b_cluster_ENG_FRAC_EM);       }
    if ( bPtrs->FindObject("cluster_ENG_FRAC_MAX"       ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_FRAC_MAX",      &cluster_ENG_FRAC_MAX,       &b_cluster_ENG_FRAC_MAX);      }
    if ( bPtrs->FindObject("cluster_ENG_FRAC_CORE"      ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_FRAC_CORE",     &cluster_ENG_FRAC_CORE,      &b_cluster_ENG_FRAC_CORE);     }
    if ( bPtrs->FindObject("cluster_SECOND_ENG_DENS"    ) != nullptr ) { fChain->SetBranchAddress("cluster_SECOND_ENG_DENS",   &cluster_SECOND_ENG_DENS,    &b_cluster_SECOND_ENG_DENS);   }
    if ( bPtrs->FindObject("cluster_ISOLATION"          ) != nullptr ) { fChain->SetBranchAddress("cluster_ISOLATION",         &cluster_ISOLATION,          &b_cluster_ISOLATION);         }
    if ( bPtrs->FindObject("cluster_ENG_BAD_CELLS"      ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_BAD_CELLS",     &cluster_ENG_BAD_CELLS,      &b_cluster_ENG_BAD_CELLS);     }
    if ( bPtrs->FindObject("cluster_N_BAD_CELLS"        ) != nullptr ) { fChain->SetBranchAddress("cluster_N_BAD_CELLS",       &cluster_N_BAD_CELLS,        &b_cluster_N_BAD_CELLS);       }
    if ( bPtrs->FindObject("cluster_N_BAD_CELLS_CORR"   ) != nullptr ) { fChain->SetBranchAddress("cluster_N_BAD_CELLS_CORR",  &cluster_N_BAD_CELLS_CORR,   &b_cluster_N_BAD_CELLS_CORR);  }
    if ( bPtrs->FindObject("cluster_BAD_CELLS_CORR_E"   ) != nullptr ) { fChain->SetBranchAddress("cluster_BAD_CELLS_CORR_E",  &cluster_BAD_CELLS_CORR_E,   &b_cluster_BAD_CELLS_CORR_E);  }
    if ( bPtrs->FindObject("cluster_BADLARQ_FRAC"       ) != nullptr ) { fChain->SetBranchAddress("cluster_BADLARQ_FRAC",      &cluster_BADLARQ_FRAC,       &b_cluster_BADLARQ_FRAC);      }
    if ( bPtrs->FindObject("cluster_ENG_POS"            ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_POS",           &cluster_ENG_POS,            &b_cluster_ENG_POS);           }
    if ( bPtrs->FindObject("cluster_SIGNIFICANCE"       ) != nullptr ) { fChain->SetBranchAddress("cluster_SIGNIFICANCE",      &cluster_SIGNIFICANCE,       &b_cluster_SIGNIFICANCE);      }
    if ( bPtrs->FindObject("cluster_CELL_SIGNIFICANCE"  ) != nullptr ) { fChain->SetBranchAddress("cluster_CELL_SIGNIFICANCE", &cluster_CELL_SIGNIFICANCE,  &b_cluster_CELL_SIGNIFICANCE); }
    if ( bPtrs->FindObject("cluster_CELL_SIG_SAMPLING"  ) != nullptr ) { fChain->SetBranchAddress("cluster_CELL_SIG_SAMPLING", &cluster_CELL_SIG_SAMPLING,  &b_cluster_CELL_SIG_SAMPLING); }
    if ( bPtrs->FindObject("cluster_AVG_LAR_Q"          ) != nullptr ) { fChain->SetBranchAddress("cluster_AVG_LAR_Q",         &cluster_AVG_LAR_Q,          &b_cluster_AVG_LAR_Q);         }
    if ( bPtrs->FindObject("cluster_AVG_TILE_Q"         ) != nullptr ) { fChain->SetBranchAddress("cluster_AVG_TILE_Q",        &cluster_AVG_TILE_Q,         &b_cluster_AVG_TILE_Q);        }
    if ( bPtrs->FindObject("cluster_ENG_BAD_HV_CELLS"   ) != nullptr ) { fChain->SetBranchAddress("cluster_ENG_BAD_HV_CELLS",  &cluster_ENG_BAD_HV_CELLS,   &b_cluster_ENG_BAD_HV_CELLS);  }
    if ( bPtrs->FindObject("cluster_N_BAD_HV_CELLS"     ) != nullptr ) { fChain->SetBranchAddress("cluster_N_BAD_HV_CELLS",    &cluster_N_BAD_HV_CELLS,     &b_cluster_N_BAD_HV_CELLS);    }
    if ( bPtrs->FindObject("cluster_PTD"                ) != nullptr ) { fChain->SetBranchAddress("cluster_PTD",               &cluster_PTD,                &b_cluster_PTD);               }
    if ( bPtrs->FindObject("cluster_MASS"               ) != nullptr ) { fChain->SetBranchAddress("cluster_MASS",              &cluster_MASS,               &b_cluster_MASS);              }
    if ( bPtrs->FindObject("cluster_SECOND_TIME"        ) != nullptr ) { fChain->SetBranchAddress("cluster_SECOND_TIME",       &cluster_SECOND_TIME,        &b_cluster_SECOND_TIME);       }
    if ( bPtrs->FindObject("CalibratedE"                ) != nullptr ) { fChain->SetBranchAddress("CalibratedE",               &CalibratedE,                &b_CalibratedE);               }
    // add EM_Shower branch if not already existing
    if ( bPtrs->FindObject("EM_Shower") == nullptr ) { 
      printf("[ClusterTree::Init(...)] WARNING no EM_Shower branch found in input tree \042%s\042\n",fChain->GetName()); 
    } else {
      fChain->SetBranchAddress("EM_Shower", &EM_Shower, &b_EM_Shower);
    }
    // add EM_Pro branch if not already existing
    if ( bPtrs->FindObject("EM_Pro")    == nullptr ) { 
      printf("[ClusterTree::Init(...)] WARNING no EM_Pro branch found in input tree \042%s\042\n",fChain->GetName());
    } else {  
      fChain->SetBranchAddress("EM_Pro",    &EM_Pro,    &b_EM_Pro);
    }
    // new branches 
    b_Delta_E       = fChain->Branch("Delta_E",       &Delta_E,       "Delta_E/D")      ; /* FIXME */ fChain->SetBranchAddress("Delta_E",       &Delta_E,       &b_Delta_E      );
    b_Delta_Calib_E = fChain->Branch("Delta_Calib_E", &Delta_Calib_E, "Delta_Calib_E/D"); /* FIXME */ fChain->SetBranchAddress("Delta_Calib_E", &Delta_Calib_E, &b_Delta_Calib_E);

    Notify();
}

void ClusterTree::Init_EM(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;

    // -- all new branches linked to the same variables in memory as the input tree
    EM_b_runNumber                   = tree->Branch("runNumber",                 &runNumber,                  "runNumber/L"                  );
    EM_b_eventNumber                 = tree->Branch("eventNumber",               &eventNumber,                "eventNumber/L"                );
    EM_b_truthE                      = tree->Branch("truthE",                    &truthE,                     "truthE/D"                     );
    EM_b_truthPt                     = tree->Branch("truthPt",                   &truthPt,                    "truthPt/D"                    );
    EM_b_truthEta                    = tree->Branch("truthEta",                  &truthEta,                   "truthEta/D"                   );
    EM_b_truthPhi                    = tree->Branch("truthPhi",                  &truthPhi,                   "truthPhi/D"                   );
    EM_b_truthPDG                    = tree->Branch("truthPDG",                  &truthPDG,                   "truthPDG/L"                   );
    EM_b_nCluster                    = tree->Branch("nCluster",                  &nCluster,                   "nCluster/L"                   );
    EM_b_clusterIndex                = tree->Branch("clusterIndex",              &clusterIndex,               "clusterIndex/L"               );
    EM_b_cluster_nCells              = tree->Branch("cluster_nCells",            &cluster_nCells,             "cluster_nCells/L"             );
    EM_b_cluster_nCells_tot          = tree->Branch("cluster_nCells_tot",        &cluster_nCells_tot,         "cluster_nCells_tot/L"         );
    EM_b_clusterECalib               = tree->Branch("clusterECalib",             &clusterECalib,              "clusterECalib/D"              );
    EM_b_clusterPtCalib              = tree->Branch("clusterPtCalib",            &clusterPtCalib,             "clusterPtCalib/D"             );
    EM_b_clusterEtaCalib             = tree->Branch("clusterEtaCalib",           &clusterEtaCalib,            "clusterEtaCalib/D"            );
    EM_b_clusterPhiCalib             = tree->Branch("clusterPhiCalib",           &clusterPhiCalib,            "clusterPhiCalib/D"            );
    EM_b_cluster_sumCellECalib       = tree->Branch("cluster_sumCellECalib",     &cluster_sumCellECalib,      "cluster_sumCellECalib/D"      );
    EM_b_cluster_fracECalib          = tree->Branch("cluster_fracECalib",        &cluster_fracECalib,         "cluster_fracECalib/D"         );
    EM_b_cluster_fracECalib_ref      = tree->Branch("cluster_fracECalib_ref",    &cluster_fracECalib_ref,     "cluster_fracECalib_ref/D"     );
    EM_b_clusterE                    = tree->Branch("clusterE",                  &clusterE,                   "clusterE/D"                   );
    EM_b_clusterPt                   = tree->Branch("clusterPt",                 &clusterPt,                  "clusterPt/D"                  );
    EM_b_clusterEta                  = tree->Branch("clusterEta",                &clusterEta,                 "clusterEta/D"                 );
    EM_b_clusterPhi                  = tree->Branch("clusterPhi",                &clusterPhi,                 "clusterPhi/D"                 );
    EM_b_cluster_sumCellE            = tree->Branch("cluster_sumCellE",          &cluster_sumCellE,           "cluster_sumCellE/D"           );			    
    EM_b_cluster_EM_PROBABILITY      = tree->Branch("cluster_EM_PROBABILITY",    &cluster_EM_PROBABILITY,     "cluster_EM_PROBABILITY/D"     );
    EM_b_cluster_sumCellE            = tree->Branch("cluster_sumCellE",          &cluster_sumCellE,           "cluster_sumCellE/D"           );
    EM_b_cluster_fracE               = tree->Branch("cluster_fracE",             &cluster_fracE,              "cluster_fracE/D"              );    
    EM_b_cluster_fracE_ref           = tree->Branch("cluster_fracE_ref",         &cluster_fracE_ref,          "cluster_fracE_ref/D"          );
    EM_b_cluster_EM_PROBABILITY      = tree->Branch("cluster_EM_PROBABILITY",    &cluster_EM_PROBABILITY,     "cluster_EM_PROBABILITY/D"     );
    EM_b_cluster_time                = tree->Branch("cluster_time",              &cluster_time,               "cluster_time/D"               );
    EM_b_cluster_HAD_WEIGHT          = tree->Branch("cluster_HAD_WEIGHT",        &cluster_HAD_WEIGHT,         "cluster_HAD_WEIGHT/D"         );
    EM_b_cluster_OOC_WEIGHT          = tree->Branch("cluster_OOC_WEIGHT",        &cluster_OOC_WEIGHT,         "cluster_OOC_WEIGHT/D"         );
    EM_b_cluster_DM_WEIGHT           = tree->Branch("cluster_DM_WEIGHT",         &cluster_DM_WEIGHT,          "cluster_DM_WEIGHT/D"          );
    EM_b_cluster_ENG_CALIB_TOT       = tree->Branch("cluster_ENG_CALIB_TOT",     &cluster_ENG_CALIB_TOT,      "cluster_ENG_CALIB_TOT/D"      );
    EM_b_cluster_ENG_CALIB_OUT_T     = tree->Branch("cluster_ENG_CALIB_OUT_T",   &cluster_ENG_CALIB_OUT_T,    "cluster_ENG_CALIB_OUT_T/D"    );
    EM_b_cluster_ENG_CALIB_DEAD_TOT  = tree->Branch("cluster_ENG_CALIB_DEAD_TOT",&cluster_ENG_CALIB_DEAD_TOT, "cluster_ENG_CALIB_DEAD_TOT/D" );
    EM_b_cluster_CENTER_MAG          = tree->Branch("cluster_CENTER_MAG",        &cluster_CENTER_MAG,         "cluster_CENTER_MAG/D"         );
    EM_b_cluster_FIRST_ENG_DENS      = tree->Branch("cluster_FIRST_ENG_DENS",    &cluster_FIRST_ENG_DENS,     "cluster_FIRST_ENG_DENS/D"     );
    EM_b_cluster_FIRST_PHI           = tree->Branch("cluster_FIRST_PHI",         &cluster_FIRST_PHI,          "cluster_FIRST_PHI/D"          );
    EM_b_cluster_FIRST_ETA           = tree->Branch("cluster_FIRST_ETA",         &cluster_FIRST_ETA,          "cluster_FIRST_ETA/D"          );
    EM_b_cluster_SECOND_R            = tree->Branch("cluster_SECOND_R",          &cluster_SECOND_R,           "cluster_SECOND_R/D"           );
    EM_b_cluster_SECOND_LAMBDA       = tree->Branch("cluster_SECOND_LAMBDA",     &cluster_SECOND_LAMBDA,      "cluster_SECOND_LAMBDA/D"      );
    EM_b_cluster_DELTA_PHI           = tree->Branch("cluster_DELTA_PHI",         &cluster_DELTA_PHI,          "cluster_DELTA_PHI/D"          );
    EM_b_cluster_DELTA_THETA         = tree->Branch("cluster_DELTA_THETA",       &cluster_DELTA_THETA,        "cluster_DELTA_THETA/D"        );
    EM_b_cluster_DELTA_ALPHA         = tree->Branch("cluster_DELTA_ALPHA",       &cluster_DELTA_ALPHA,        "cluster_DELTA_ALPHA/D"        );
    EM_b_cluster_CENTER_X            = tree->Branch("cluster_CENTER_X",          &cluster_CENTER_X,           "cluster_CENTER_X/D"           );
    EM_b_cluster_CENTER_Y            = tree->Branch("cluster_CENTER_Y",          &cluster_CENTER_Y,           "cluster_CENTER_Y/D"           );
    EM_b_cluster_CENTER_Z            = tree->Branch("cluster_CENTER_Z",          &cluster_CENTER_Z,           "cluster_CENTER_Z/D"           );
    EM_b_cluster_CENTER_LAMBDA       = tree->Branch("cluster_CENTER_LAMBDA",     &cluster_CENTER_LAMBDA,      "cluster_CENTER_LAMBDA/D"      );
    EM_b_cluster_LATERAL             = tree->Branch("cluster_LATERAL",           &cluster_LATERAL,            "cluster_LATERAL/D"            );
    EM_b_cluster_LONGITUDINAL        = tree->Branch("cluster_LONGITUDINAL",      &cluster_LONGITUDINAL,       "cluster_LONGITUDINAL/D"       );
    EM_b_cluster_ENG_FRAC_EM         = tree->Branch("cluster_ENG_FRAC_EM",       &cluster_ENG_FRAC_EM,        "cluster_ENG_FRAC_EM/D"        );
    EM_b_cluster_ENG_FRAC_MAX        = tree->Branch("cluster_ENG_FRAC_MAX",      &cluster_ENG_FRAC_MAX,       "cluster_ENG_FRAC_MAX/D"       );
    EM_b_cluster_ENG_FRAC_CORE       = tree->Branch("cluster_ENG_FRAC_CORE",     &cluster_ENG_FRAC_CORE,      "cluster_ENG_FRAC_CORE/D"      );
    EM_b_cluster_SECOND_ENG_DENS     = tree->Branch("cluster_SECOND_ENG_DENS",   &cluster_SECOND_ENG_DENS,    "cluster_SECOND_ENG_DENS/D"    );
    EM_b_cluster_ISOLATION           = tree->Branch("cluster_ISOLATION",         &cluster_ISOLATION,          "cluster_ISOLATION/D"          );
    EM_b_cluster_ENG_BAD_CELLS       = tree->Branch("cluster_ENG_BAD_CELLS",     &cluster_ENG_BAD_CELLS,      "cluster_ENG_BAD_CELLS/D"      );
    EM_b_cluster_N_BAD_CELLS         = tree->Branch("cluster_N_BAD_CELLS",       &cluster_N_BAD_CELLS,        "cluster_N_BAD_CELLS/L"        );
    EM_b_cluster_N_BAD_CELLS_CORR    = tree->Branch("cluster_N_BAD_CELLS_CORR",  &cluster_N_BAD_CELLS_CORR,   "cluster_N_BAD_CELLS_CORR/L"   );
    EM_b_cluster_BAD_CELLS_CORR_E    = tree->Branch("cluster_BAD_CELLS_CORR_E",  &cluster_BAD_CELLS_CORR_E,   "cluster_BAD_CELLS_CORR_E/D"   );
    EM_b_cluster_BADLARQ_FRAC        = tree->Branch("cluster_BADLARQ_FRAC",      &cluster_BADLARQ_FRAC,       "cluster_BADLARQ_FRAC/L"       );
    EM_b_cluster_ENG_POS             = tree->Branch("cluster_ENG_POS",           &cluster_ENG_POS,            "cluster_ENG_POS/D"            );
    EM_b_cluster_SIGNIFICANCE        = tree->Branch("cluster_SIGNIFICANCE",      &cluster_SIGNIFICANCE,       "cluster_SIGNIFICANCE/D"       );
    EM_b_cluster_CELL_SIGNIFICANCE   = tree->Branch("cluster_CELL_SIGNIFICANCE", &cluster_CELL_SIGNIFICANCE,  "cluster_CELL_SIGNIFICANCE/D"  );
    EM_b_cluster_CELL_SIG_SAMPLING   = tree->Branch("cluster_CELL_SIG_SAMPLING", &cluster_CELL_SIG_SAMPLING,  "cluster_CELL_SIG_SAMPLING/D"  );
    EM_b_cluster_AVG_LAR_Q           = tree->Branch("cluster_AVG_LAR_Q",         &cluster_AVG_LAR_Q,          "cluster_AVG_LAR_Q/D"          );
    EM_b_cluster_AVG_TILE_Q          = tree->Branch("cluster_AVG_TILE_Q",        &cluster_AVG_TILE_Q,         "cluster_AVG_TILE_Q/D"         );
    EM_b_cluster_ENG_BAD_HV_CELLS    = tree->Branch("cluster_ENG_BAD_HV_CELLS",  &cluster_ENG_BAD_HV_CELLS,   "cluster_ENG_BAD_HV_CELLS/D"   );
    EM_b_cluster_N_BAD_HV_CELLS      = tree->Branch("cluster_N_BAD_HV_CELLS",    &cluster_N_BAD_HV_CELLS,     "cluster_N_BAD_HV_CELLS/L"     );
    EM_b_cluster_PTD                 = tree->Branch("cluster_PTD",               &cluster_PTD,                "cluster_PTD/D"                );
    EM_b_cluster_MASS                = tree->Branch("cluster_MASS",              &cluster_MASS,               "cluster_MASS/D"               );
    EM_b_cluster_SECOND_TIME         = tree->Branch("cluster_SECOND_TIME",       &cluster_SECOND_TIME,        "cluster_SECOND_TIME/D"        );      
    EM_b_CalibratedE                 = tree->Branch("CalibratedE",               &CalibratedE,                "CalibratedE/D"                );              
    EM_b_EM_Shower                   = tree->Branch("EM_Shower",                 &EM_Shower,                  "EM_Shower/L"                  );
    EM_b_EM_Pro                      = tree->Branch("EM_Pro",                    &EM_Pro,                     "EM_Pro/D"                     );
    EM_b_Delta_E                     = tree->Branch("Delta_E",                   &Delta_E,                    "Delta_E/D"                    );
    EM_b_Delta_Calib_E               = tree->Branch("Delta_Calib_E",             &Delta_Calib_E,              "Delta_Calib_E/D"              );

    Notify();
    return;
}

void ClusterTree::Init_Had(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // Init() will be called many times when running on PROOF
    // code, but the routine can be extended by the user if needed.
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;

    // -- all new branches linked to the same variables in memory as the input tree
    Had_b_runNumber                   = tree->Branch("runNumber",                 &runNumber,                  "runNumber/L"                  );
    Had_b_eventNumber                 = tree->Branch("eventNumber",               &eventNumber,                "eventNumber/L"                );
    Had_b_truthE                      = tree->Branch("truthE",                    &truthE,                     "truthE/D"                     );
    Had_b_truthPt                     = tree->Branch("truthPt",                   &truthPt,                    "truthPt/D"                    );
    Had_b_truthEta                    = tree->Branch("truthEta",                  &truthEta,                   "truthEta/D"                   );
    Had_b_truthPhi                    = tree->Branch("truthPhi",                  &truthPhi,                   "truthPhi/D"                   );
    Had_b_truthPDG                    = tree->Branch("truthPDG",                  &truthPDG,                   "truthPDG/L"                   );
    Had_b_nCluster                    = tree->Branch("nCluster",                  &nCluster,                   "nCluster/L"                   );
    Had_b_clusterIndex                = tree->Branch("clusterIndex",              &clusterIndex,               "clusterIndex/L"               );
    Had_b_cluster_nCells              = tree->Branch("cluster_nCells",            &cluster_nCells,             "cluster_nCells/L"             );
    Had_b_cluster_nCells_tot          = tree->Branch("cluster_nCells_tot",        &cluster_nCells_tot,         "cluster_nCells_tot/L"         );
    Had_b_clusterECalib               = tree->Branch("clusterECalib",             &clusterECalib,              "clusterECalib/D"              );
    Had_b_clusterPtCalib              = tree->Branch("clusterPtCalib",            &clusterPtCalib,             "clusterPtCalib/D"             );
    Had_b_clusterEtaCalib             = tree->Branch("clusterEtaCalib",           &clusterEtaCalib,            "clusterEtaCalib/D"            );
    Had_b_clusterPhiCalib             = tree->Branch("clusterPhiCalib",           &clusterPhiCalib,            "clusterPhiCalib/D"            );
    Had_b_cluster_sumCellECalib       = tree->Branch("cluster_sumCellECalib",     &cluster_sumCellECalib,      "cluster_sumCellECalib/D"      );
    Had_b_cluster_fracECalib          = tree->Branch("cluster_fracECalib",        &cluster_fracECalib,         "cluster_fracECalib/D"         );
    Had_b_cluster_fracECalib_ref      = tree->Branch("cluster_fracECalib_ref",    &cluster_fracECalib_ref,     "cluster_fracECalib_ref/D"     );
    Had_b_clusterE                    = tree->Branch("clusterE",                  &clusterE,                   "clusterE/D"                   );
    Had_b_clusterPt                   = tree->Branch("clusterPt",                 &clusterPt,                  "clusterPt/D"                  );
    Had_b_clusterEta                  = tree->Branch("clusterEta",                &clusterEta,                 "clusterEta/D"                 );
    Had_b_clusterPhi                  = tree->Branch("clusterPhi",                &clusterPhi,                 "clusterPhi/D"                 );
    Had_b_cluster_sumCellE            = tree->Branch("cluster_sumCellE",          &cluster_sumCellE,           "cluster_sumCellE/D"           );			    
    Had_b_cluster_EM_PROBABILITY      = tree->Branch("cluster_EM_PROBABILITY",    &cluster_EM_PROBABILITY,     "cluster_EM_PROBABILITY/D"     );
    Had_b_cluster_sumCellE            = tree->Branch("cluster_sumCellE",          &cluster_sumCellE,           "cluster_sumCellE/D"           );
    Had_b_cluster_fracE               = tree->Branch("cluster_fracE",             &cluster_fracE,              "cluster_fracE/D"              );    
    Had_b_cluster_fracE_ref           = tree->Branch("cluster_fracE_ref",         &cluster_fracE_ref,          "cluster_fracE_ref/D"          );
    Had_b_cluster_EM_PROBABILITY      = tree->Branch("cluster_EM_PROBABILITY",    &cluster_EM_PROBABILITY,     "cluster_EM_PROBABILITY/D"     );
    Had_b_cluster_time                = tree->Branch("cluster_time",              &cluster_time,               "cluster_time/D"               );
    Had_b_cluster_HAD_WEIGHT          = tree->Branch("cluster_HAD_WEIGHT",        &cluster_HAD_WEIGHT,         "cluster_HAD_WEIGHT/D"         );
    Had_b_cluster_OOC_WEIGHT          = tree->Branch("cluster_OOC_WEIGHT",        &cluster_OOC_WEIGHT,         "cluster_OOC_WEIGHT/D"         );
    Had_b_cluster_DM_WEIGHT           = tree->Branch("cluster_DM_WEIGHT",         &cluster_DM_WEIGHT,          "cluster_DM_WEIGHT/D"          );
    Had_b_cluster_ENG_CALIB_TOT       = tree->Branch("cluster_ENG_CALIB_TOT",     &cluster_ENG_CALIB_TOT,      "cluster_ENG_CALIB_TOT/D"      );
    Had_b_cluster_ENG_CALIB_OUT_T     = tree->Branch("cluster_ENG_CALIB_OUT_T",   &cluster_ENG_CALIB_OUT_T,    "cluster_ENG_CALIB_OUT_T/D"    );
    Had_b_cluster_ENG_CALIB_DEAD_TOT  = tree->Branch("cluster_ENG_CALIB_DEAD_TOT",&cluster_ENG_CALIB_DEAD_TOT, "cluster_ENG_CALIB_DEAD_TOT/D" );
    Had_b_cluster_CENTER_MAG          = tree->Branch("cluster_CENTER_MAG",        &cluster_CENTER_MAG,         "cluster_CENTER_MAG/D"         );
    Had_b_cluster_FIRST_ENG_DENS      = tree->Branch("cluster_FIRST_ENG_DENS",    &cluster_FIRST_ENG_DENS,     "cluster_FIRST_ENG_DENS/D"     );
    Had_b_cluster_FIRST_PHI           = tree->Branch("cluster_FIRST_PHI",         &cluster_FIRST_PHI,          "cluster_FIRST_PHI/D"          );
    Had_b_cluster_FIRST_ETA           = tree->Branch("cluster_FIRST_ETA",         &cluster_FIRST_ETA,          "cluster_FIRST_ETA/D"          );
    Had_b_cluster_SECOND_R            = tree->Branch("cluster_SECOND_R",          &cluster_SECOND_R,           "cluster_SECOND_R/D"           );
    Had_b_cluster_SECOND_LAMBDA       = tree->Branch("cluster_SECOND_LAMBDA",     &cluster_SECOND_LAMBDA,      "cluster_SECOND_LAMBDA/D"      );
    Had_b_cluster_DELTA_PHI           = tree->Branch("cluster_DELTA_PHI",         &cluster_DELTA_PHI,          "cluster_DELTA_PHI/D"          );
    Had_b_cluster_DELTA_THETA         = tree->Branch("cluster_DELTA_THETA",       &cluster_DELTA_THETA,        "cluster_DELTA_THETA/D"        );
    Had_b_cluster_DELTA_ALPHA         = tree->Branch("cluster_DELTA_ALPHA",       &cluster_DELTA_ALPHA,        "cluster_DELTA_ALPHA/D"        );
    Had_b_cluster_CENTER_X            = tree->Branch("cluster_CENTER_X",          &cluster_CENTER_X,           "cluster_CENTER_X/D"           );
    Had_b_cluster_CENTER_Y            = tree->Branch("cluster_CENTER_Y",          &cluster_CENTER_Y,           "cluster_CENTER_Y/D"           );
    Had_b_cluster_CENTER_Z            = tree->Branch("cluster_CENTER_Z",          &cluster_CENTER_Z,           "cluster_CENTER_Z/D"           );
    Had_b_cluster_CENTER_LAMBDA       = tree->Branch("cluster_CENTER_LAMBDA",     &cluster_CENTER_LAMBDA,      "cluster_CENTER_LAMBDA/D"      );
    Had_b_cluster_LATERAL             = tree->Branch("cluster_LATERAL",           &cluster_LATERAL,            "cluster_LATERAL/D"            );
    Had_b_cluster_LONGITUDINAL        = tree->Branch("cluster_LONGITUDINAL",      &cluster_LONGITUDINAL,       "cluster_LONGITUDINAL/D"       );
    Had_b_cluster_ENG_FRAC_EM         = tree->Branch("cluster_ENG_FRAC_EM",       &cluster_ENG_FRAC_EM,        "cluster_ENG_FRAC_EM/D"        );
    Had_b_cluster_ENG_FRAC_MAX        = tree->Branch("cluster_ENG_FRAC_MAX",      &cluster_ENG_FRAC_MAX,       "cluster_ENG_FRAC_MAX/D"       );
    Had_b_cluster_ENG_FRAC_CORE       = tree->Branch("cluster_ENG_FRAC_CORE",     &cluster_ENG_FRAC_CORE,      "cluster_ENG_FRAC_CORE/D"      );
    Had_b_cluster_SECOND_ENG_DENS     = tree->Branch("cluster_SECOND_ENG_DENS",   &cluster_SECOND_ENG_DENS,    "cluster_SECOND_ENG_DENS/D"    );
    Had_b_cluster_ISOLATION           = tree->Branch("cluster_ISOLATION",         &cluster_ISOLATION,          "cluster_ISOLATION/D"          );
    Had_b_cluster_ENG_BAD_CELLS       = tree->Branch("cluster_ENG_BAD_CELLS",     &cluster_ENG_BAD_CELLS,      "cluster_ENG_BAD_CELLS/D"      );
    Had_b_cluster_N_BAD_CELLS         = tree->Branch("cluster_N_BAD_CELLS",       &cluster_N_BAD_CELLS,        "cluster_N_BAD_CELLS/L"        );
    Had_b_cluster_N_BAD_CELLS_CORR    = tree->Branch("cluster_N_BAD_CELLS_CORR",  &cluster_N_BAD_CELLS_CORR,   "cluster_N_BAD_CELLS_CORR/L"   );
    Had_b_cluster_BAD_CELLS_CORR_E    = tree->Branch("cluster_BAD_CELLS_CORR_E",  &cluster_BAD_CELLS_CORR_E,   "cluster_BAD_CELLS_CORR_E/D"   );
    Had_b_cluster_BADLARQ_FRAC        = tree->Branch("cluster_BADLARQ_FRAC",      &cluster_BADLARQ_FRAC,       "cluster_BADLARQ_FRAC/L"       );
    Had_b_cluster_ENG_POS             = tree->Branch("cluster_ENG_POS",           &cluster_ENG_POS,            "cluster_ENG_POS/D"            );
    Had_b_cluster_SIGNIFICANCE        = tree->Branch("cluster_SIGNIFICANCE",      &cluster_SIGNIFICANCE,       "cluster_SIGNIFICANCE/D"       );
    Had_b_cluster_CELL_SIGNIFICANCE   = tree->Branch("cluster_CELL_SIGNIFICANCE", &cluster_CELL_SIGNIFICANCE,  "cluster_CELL_SIGNIFICANCE/D"  );
    Had_b_cluster_CELL_SIG_SAMPLING   = tree->Branch("cluster_CELL_SIG_SAMPLING", &cluster_CELL_SIG_SAMPLING,  "cluster_CELL_SIG_SAMPLING/D"  );
    Had_b_cluster_AVG_LAR_Q           = tree->Branch("cluster_AVG_LAR_Q",         &cluster_AVG_LAR_Q,          "cluster_AVG_LAR_Q/D"          );
    Had_b_cluster_AVG_TILE_Q          = tree->Branch("cluster_AVG_TILE_Q",        &cluster_AVG_TILE_Q,         "cluster_AVG_TILE_Q/D"         );
    Had_b_cluster_ENG_BAD_HV_CELLS    = tree->Branch("cluster_ENG_BAD_HV_CELLS",  &cluster_ENG_BAD_HV_CELLS,   "cluster_ENG_BAD_HV_CELLS/D"   );
    Had_b_cluster_N_BAD_HV_CELLS      = tree->Branch("cluster_N_BAD_HV_CELLS",    &cluster_N_BAD_HV_CELLS,     "cluster_N_BAD_HV_CELLS/L"     );
    Had_b_cluster_PTD                 = tree->Branch("cluster_PTD",               &cluster_PTD,                "cluster_PTD/D"                );
    Had_b_cluster_MASS                = tree->Branch("cluster_MASS",              &cluster_MASS,               "cluster_MASS/D"               );
    Had_b_cluster_SECOND_TIME         = tree->Branch("cluster_SECOND_TIME",       &cluster_SECOND_TIME,        "cluster_SECOND_TIME/D"        );      
    Had_b_CalibratedE                 = tree->Branch("CalibratedE",               &CalibratedE,                "CalibratedE/D"                );              
    Had_b_EM_Shower                   = tree->Branch("EM_Shower",                 &EM_Shower,                  "EM_Shower/L"                  );
    Had_b_EM_Pro                      = tree->Branch("EM_Pro",                    &EM_Pro,                     "EM_Pro/D"                     );
    Had_b_Delta_E                     = tree->Branch("Delta_E",                   &Delta_E,                    "Delta_E/D"                    );
    Had_b_Delta_Calib_E               = tree->Branch("Delta_Calib_E",             &Delta_Calib_E,              "Delta_Calib_E/D"              );

    Notify();

    return;
}

Bool_t ClusterTree::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void ClusterTree::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t ClusterTree::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef ClusterTree_cxx
