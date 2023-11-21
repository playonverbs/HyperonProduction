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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"

#include "ubana/HyperonProduction/Utils.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

#include <string>
#include <vector>

#define FNLOG(msg)	std::cout << "[" << __PRETTY_FUNCTION__ << "] " << msg << std::endl

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

		unsigned int fNPFParticles;
		unsigned int fNPrimaryChildren;

		// FHICL Params
		std::string fSliceLabel;
		std::string fPFParticleLabel;
		std::string fTrackLabel;
		std::string fShowerLabel;

		bool fIsData;
		bool fDebug = false;

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

		TTree* fTree;
};


hyperon::HyperonProduction::HyperonProduction(fhicl::ParameterSet const& p)
	: EDAnalyzer{p},
	fSliceLabel(p.get<std::string>("SliceLabel")),
	fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
	fTrackLabel(p.get<std::string>("TrackLabel")),
	fShowerLabel(p.get<std::string>("ShowerLabel"))
	// More initializers here.
{
	// Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hyperon::HyperonProduction::analyze(art::Event const& e)
{
	// Implementation of required member function here.

	_trk_length.clear();
	_shr_length.clear();
	
	_n_primary_tracks  = 0;
	_n_primary_showers = 0;

	_run    = e.run();
	_subrun = e.subRun();
	_event  = e.id().event();

	/* art::ValidHandle<std::vector<simb::MCTruth>> mctruthListHandle; */
	/* std::vector<art::Ptr<simb::MCTruth>>    mcTrVect; */

	/* if (!e.getByLabel(f_generatorProducer, )) */

	art::ValidHandle<std::vector<recob::Slice>> sliceHandle =
		e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
	std::vector<art::Ptr<recob::Slice>> sliceVector;

	if (sliceHandle.isValid())
		art::fill_ptr_vector(sliceVector, sliceHandle);

	art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);

	int nuID = -1;
	int nuSliceKey = -1;

	// iterate through all collected Slice objects.
	for (const art::Ptr<recob::Slice> &slice : sliceVector)
	{
		// collect all PFPs associated with the current slice key
		std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));

		for (const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs)
		{
			const bool isPrimary(slicePFP->IsPrimary());
			const bool isNeutrino(std::abs(slicePFP->PdgCode()) == 12 || std::abs(slicePFP->PdgCode()) == 14);

			if (!(isPrimary && isNeutrino))
				continue;

			nuSliceKey = slice.key();
			nuID = slicePFP->Self();
			fNPFParticles = slicePFPs.size();
			fNPrimaryChildren = slicePFP->NumDaughters();

			break;
		}

		// this assumes only one neutrino hierarchy among all slice.
		if (nuID >= 0)
			break;
	}

	if (nuSliceKey < 0)
	{
		if (fDebug)
			FNLOG("no nuSliceKey found");
		// TODO handle failure case properly
		return;
	}

	art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle =
		e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
	art::FindManyP<recob::Track>  pfpTrackAssoc(pfpHandle, e, fTrackLabel);
	art::FindManyP<recob::Shower> pfpShowerAssoc(pfpHandle, e, fShowerLabel);

	std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));

	for (const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
	{
		if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
			continue;

		// Handle Tracks

		std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());

		if (tracks.size() != 1)
			continue;

		_n_primary_tracks++;
		art::Ptr<recob::Track> track = tracks.at(0);
		_trk_length.push_back(track->Length());

		// Handle Showers

		std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(nuSlicePFP.key());

		if (showers.size() != 1)
			continue;

		_n_primary_showers++;
		art::Ptr<recob::Shower> shower = showers.at(0);
		
		if (shower->has_length()) {
			_shr_length.push_back(shower->Length());
		} else {
			_shr_length.push_back(-1.0);
		}
	}

	fTree->Fill();
}

void hyperon::HyperonProduction::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	fTree = tfs->make<TTree>("tree", "Output TTree");

	fTree->Branch("run",    &_run);
	fTree->Branch("subrun", &_subrun);
	fTree->Branch("event",  &_event);
	fTree->Branch("n_primary_tracks",  &_n_primary_tracks);
	fTree->Branch("trk_length", &_trk_length);

	fTree->Branch("n_primary_showers", &_n_primary_showers);
	fTree->Branch("shr_length", &_shr_length);
}

void hyperon::HyperonProduction::endJob()
{
	if (fDebug)
		FNLOG("ending job");
	return;
}

DEFINE_ART_MODULE(hyperon::HyperonProduction)
