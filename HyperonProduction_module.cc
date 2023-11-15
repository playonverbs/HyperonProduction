////////////////////////////////////////////////////////////////////////
// Class:       HyperonProduction
// Plugin Type: analyzer (art v3_01_02)
// File:        HyperonProduction_module.cc
//
// Generated at Mon Nov 13 10:53:34 2023 by Niam Patel using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

#include <string>
#include <vector>

#define FNLOG(msg)	std::cout << "[" << __PRETTY_FUNCTION__ << "] " << msg << std::endl;

namespace hyperon {
	class HyperonProduction;
}

class hyperon::HyperonProduction : public art::EDAnalyzer {
	public:
		explicit HyperonProduction(fhicl::ParameterSet const& p);
		// The compiler-generated destructor is fine for non-base
		// classes without bare pointers or other resource use.

		// Plugins should not be copied or assigned.
		HyperonProduction(HyperonProduction const&) = delete;
		HyperonProduction(HyperonProduction&&) = delete;
		HyperonProduction& operator=(HyperonProduction const&) = delete;
		HyperonProduction& operator=(HyperonProduction&&) = delete;

		// Required functions.
		void analyze(art::Event const& e) override;

		// Selected optional functions.
		void beginJob() override;
		void endJob() override;

	private:

		// Declare member data here.

		TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
		TParticlePDG *muon   = TDatabasePDG::Instance()->GetParticle(13);
		TParticlePDG *pion   = TDatabasePDG::Instance()->GetParticle(211);
		TParticlePDG *photon = TDatabasePDG::Instance()->GetParticle(22);

		// FHICL Params
		std::string f_pfpProducer;
		std::string f_trackProducer;
		std::string f_showerProducer;
		std::string f_vertexProducer;
		std::string f_pidProducer;
		std::string f_hitProducer;
		std::string f_hitTruthAssnProducer;
		std::string f_trackHitAssnProducer;
		std::string f_showerHitAssnProducer;
		std::string f_metadataProducer;
		std::string f_generatorProducer;
		std::string f_g4Producer;

		std::string f_recoProducer;

		bool f_isData;
		bool f_debug = false;

		// output tree values here:

		unsigned int _run;
		unsigned int _subrun;
		unsigned int _event;

		unsigned int _n_primary_tracks;
		unsigned int _n_primary_showers;

		// RecoParticle fields
		std::vector<int>    _pdg;
		std::vector<double> _trk_shr_score;
		std::vector<double> _x;
		std::vector<double> _y;
		std::vector<double> _z;
		std::vector<double> _displacement;

		std::vector<double> _trk_length;
		std::vector<double> _trk_dir_x;
		std::vector<double> _trk_dir_y;
		std::vector<double> _trk_dir_z;
		std::vector<double> _trk_start_x;
		std::vector<double> _trk_start_y;
		std::vector<double> _trk_start_z;
		std::vector<double> _trk_end_x;
		std::vector<double> _trk_end_y;
		std::vector<double> _trk_end_z;
		std::vector<double> _trk_mean_dedx_plane0;
		std::vector<double> _trk_mean_dedx_plane1;
		std::vector<double> _trk_mean_dedx_plane2;
		std::vector<double> _trk_llrpid;

		std::vector<double> _shr_length;
		std::vector<double> _shr_dir_x;
		std::vector<double> _shr_dir_y;
		std::vector<double> _shr_dir_z;
		std::vector<double> _shr_start_x;
		std::vector<double> _shr_start_y;
		std::vector<double> _shr_start_z;
		std::vector<double> _shr_energy_plane0;
		std::vector<double> _shr_energy_plane1;
		std::vector<double> _shr_energy_plane2;
		std::vector<double> _shr_dedx_plane0;
		std::vector<double> _shr_dedx_plane1;
		std::vector<double> _shr_dedx_plane2;
		std::vector<double> _shr_open_angle;

		std::vector<bool>   _has_truth;
		std::vector<int>    _mc_truth_index;
		std::vector<int>    _trk_true_pdg;
		std::vector<double> _trk_true_energy;
		std::vector<double> _trk_true_ke;
		std::vector<double> _trk_true_px;
		std::vector<double> _trk_true_py;
		std::vector<double> _trk_true_pz;
		std::vector<double> _trk_true_length;
		std::vector<int>    _trk_true_origin;
		std::vector<double> _trk_true_purity;

		TTree* f_tree;

};


hyperon::HyperonProduction::HyperonProduction(fhicl::ParameterSet const& p)
	: EDAnalyzer{p}  // ,
	// More initializers here.
{
	// Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hyperon::HyperonProduction::analyze(art::Event const& e)
{
	// Implementation of required member function here.

	_run    = e.run();
	_subrun = e.subRun();
	_event  = e.id().event();

	/* art::ValidHandle<std::vector<simb::MCTruth>> mctruthListHandle; */
	/* std::vector<art::Ptr<simb::MCTruth>>    mcTrVect; */

	/* if (!e.getByLabel(f_generatorProducer, )) */

	art::ValidHandle<std::vector<recob::Slice>> sliceHandle =
		e.getValidHandle<std::vector<recob::Slice>>(f_recoProducer);
	std::vector<art::Ptr<recob::Slice>> sliceVector;

	if (sliceHandle.isValid())
		art::fill_ptr_vector(sliceVector, sliceHandle);

	f_tree->Fill();
}

void hyperon::HyperonProduction::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	f_tree = tfs->make<TTree>("tree", "Output TTree");

	f_tree->Branch("run",    &_run);
	f_tree->Branch("subrun", &_subrun);
	f_tree->Branch("event",  &_event);
}

void hyperon::HyperonProduction::endJob()
{

}

DEFINE_ART_MODULE(HyperonProduction)
