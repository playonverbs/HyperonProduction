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
#include "lardataobj/RecoBase/PFParticleMetadata.h"
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

	struct Config;

	// define some default error values;
	namespace bogus {
		constexpr double LENGTH = -999.0;
		constexpr double POS    = -999.9;
	}
}

// Config is an art-compatible container for fhicl parameters: this allows for
// parameter validation and descriptions.
struct hyperon::Config {
	using Name    = fhicl::Name;
	using Comment = fhicl::Comment;

	template<typename T>
	using Atom = fhicl::Atom<T>;

	template<typename T, std::size_t SZ>
	using Sequence = fhicl::Sequence<T, SZ>;

	Atom<std::string> fSliceLabel      { Name("SliceLabel"),
										 Comment("Label for recob::Slice") };
	Atom<std::string> fPFParticleLabel { Name("PFParticleLabel"),
   										 Comment("Label for recob::PFParticle") };
	Atom<std::string> fTrackLabel      { Name("TrackLabel"),
										 Comment("Label for recob::Track") };
	Atom<std::string> fShowerLabel     { Name("ShowerLabel"),
										 Comment("Label for recob::Shower") };
	Atom<std::string> fCaloLabel       { Name("CaloLabel"),
										 Comment("Label for anab::Calorimetry") };
	Atom<std::string> fGeneratorLabel  { Name("GeneratorLabel"),
										 Comment("Label for simb::MCTruth") };
	Atom<std::string> fG4Label  	   { Name("G4Label"),
										 Comment("Label for simb::MCParticle") };
	Atom<bool>        fIsData          { Name("IsData"),
										 Comment("Flag to indicate if the input is Data") };
	Atom<bool>        fDebug           { Name("Debug"),
										 Comment("Flag to enable debug messages"),
										 false };
};

class hyperon::HyperonProduction : public art::EDAnalyzer {

	public:
		using Parameters = art::EDAnalyzer::Table<Config>;

		explicit HyperonProduction(Parameters const& config);
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
		void getTrackVariables(art::Ptr<recob::Track> &track);
		void getShowerVariables(art::Ptr<recob::Shower> &shower);
		void clearTreeVariables();
		void fillNull();

		// Declare member data here.

		unsigned int fNPFParticles;
		unsigned int fNPrimaryChildren;

		// FHICL Params
		std::string fSliceLabel;
		std::string fPFParticleLabel;
		std::string fTrackLabel;
		std::string fShowerLabel;
		std::string fCaloLabel;
		std::string fGeneratorLabel;
		std::string fG4Label;
		bool fIsData;
		bool fDebug;

		// output tree values here:

		unsigned int _run;
		unsigned int _subrun;
		unsigned int _event;

		int         _mc_nu_pdg;
		double      _mc_nu_q2;
		std::string _mc_ccnc;
		std::string _mc_mode;

		unsigned int _n_slices;
		unsigned int _n_primary_tracks;
		unsigned int _n_primary_showers;

		// RecoParticle fields
		std::vector<int>    _pdg;
		std::vector<double> _x;
		std::vector<double> _y;
		std::vector<double> _z;

		std::vector<double> _trk_length;
		std::vector<double> _trk_shr_score;
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


hyperon::HyperonProduction::HyperonProduction(Parameters const& config)
	: EDAnalyzer{config},
	fSliceLabel(config().fSliceLabel()),
	fPFParticleLabel(config().fPFParticleLabel()),
	fTrackLabel(config().fTrackLabel()),
	fShowerLabel(config().fShowerLabel()),
	fCaloLabel(config().fCaloLabel()),
	fGeneratorLabel(config().fGeneratorLabel()),
	fG4Label(config().fG4Label()),
	fIsData(config().fIsData()),
	fDebug(config().fDebug())
{
	// Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hyperon::HyperonProduction::analyze(art::Event const& e)
{
	// Implementation of required member function here.
	clearTreeVariables();

	_run    = e.run();
	_subrun = e.subRun();
	_event  = e.id().event();

	// Check if this is an MC file.
	if (!fIsData) {
		art::ValidHandle<std::vector<simb::MCTruth>> mcTruthHandle =
			e.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);
		std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;

		if (mcTruthHandle.isValid())
			art::fill_ptr_vector(mcTruthVector, mcTruthHandle);

		// Fill truth information
		
		for (const art::Ptr<simb::MCTruth> &truth : mcTruthVector)
		{
			_mc_nu_pdg = truth->GetNeutrino().Nu().PdgCode();
			_mc_nu_q2  = truth->GetNeutrino().QSqr();
			_mc_ccnc   = util::GetCCNC(truth->GetNeutrino().CCNC());
			_mc_mode   = util::GetEventType(truth->GetNeutrino().Mode());

			break; // escape after first MCNeutrino
		}
	}

	art::ValidHandle<std::vector<recob::Slice>> sliceHandle =
		e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
	std::vector<art::Ptr<recob::Slice>> sliceVector;

	if (sliceHandle.isValid())
		art::fill_ptr_vector(sliceVector, sliceHandle);

	art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);

	int nuID = -1;
	int nuSliceKey = -1;

	_n_slices = sliceVector.size();

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
		fillNull();
		return;
	}

	art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle =
		e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
	art::FindManyP<recob::Track>  pfpTrackAssoc(pfpHandle, e, fTrackLabel);
	art::FindManyP<recob::Shower> pfpShowerAssoc(pfpHandle, e, fShowerLabel);
	art::FindManyP<larpandoraobj::PFParticleMetadata>
		pfpMetaAssoc(pfpHandle, e, fPFParticleLabel);

	std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));

	for (const art::Ptr<recob::PFParticle> &nuSlicePFP : nuSlicePFPs)
	{
		if (nuSlicePFP->Parent() != static_cast<long unsigned int>(nuID))
			continue;

		_pdg.push_back(nuSlicePFP->PdgCode());

		std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetas =
			pfpMetaAssoc.at(nuSlicePFP.key());

		if (!pfpMetas.empty())
		{
			art::Ptr<larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetas.at(0);

			if (pfpMeta->GetPropertiesMap().find("TrackScore")
					!= pfpMeta->GetPropertiesMap().end())
				_trk_shr_score.push_back(pfpMeta->GetPropertiesMap().at("TrackScore"));
		}
		
		// Handle Tracks
		std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(nuSlicePFP.key());

		if (tracks.size() == 1)
		{
			_n_primary_tracks++;
			art::Ptr<recob::Track> track = tracks.at(0);

			getTrackVariables(track);
		}
		// Handle Showers

		std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(nuSlicePFP.key());

		if (showers.size() == 1)
		{
			_n_primary_showers++;
			art::Ptr<recob::Shower> shower = showers.at(0);

			getShowerVariables(shower);
		}
	}

	fTree->Fill();
}

void hyperon::HyperonProduction::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	fTree = tfs->make<TTree>("OutputTree", "Output TTree");

	fTree->Branch("run",    &_run);
	fTree->Branch("subrun", &_subrun);
	fTree->Branch("event",  &_event);

	fTree->Branch("mc_nu_pdg", &_mc_nu_pdg);
	fTree->Branch("mc_nu_q2",  &_mc_nu_q2);
	fTree->Branch("mc_ccnc",   &_mc_ccnc);
	fTree->Branch("mc_mode",   &_mc_mode);

	fTree->Branch("pdg",  &_pdg);

	fTree->Branch("n_slices",         &_n_slices);
	fTree->Branch("n_primary_tracks", &_n_primary_tracks);
	fTree->Branch("trk_length",       &_trk_length);
	fTree->Branch("trk_shr_score",    &_trk_shr_score);
	fTree->Branch("trk_start_x",      &_trk_start_x);
	fTree->Branch("trk_start_y",      &_trk_start_y);
	fTree->Branch("trk_start_z",      &_trk_start_z);
	fTree->Branch("trk_end_x",        &_trk_end_x);
	fTree->Branch("trk_end_y",        &_trk_end_y);
	fTree->Branch("trk_end_z",        &_trk_end_z);
	fTree->Branch("trk_dir_x",        &_trk_dir_x);
	fTree->Branch("trk_dir_y",        &_trk_dir_y);
	fTree->Branch("trk_dir_z",        &_trk_dir_z);

	fTree->Branch("n_primary_showers", &_n_primary_showers);
	fTree->Branch("shr_length",        &_shr_length);
	fTree->Branch("shr_start_x",       &_shr_start_x);
	fTree->Branch("shr_start_y",       &_shr_start_y);
	fTree->Branch("shr_start_z",       &_shr_start_z);
}

void hyperon::HyperonProduction::endJob()
{
	if (fDebug)
		FNLOG("ending job");
	return;
}

void hyperon::HyperonProduction::getTrackVariables(art::Ptr<recob::Track> &track)
{
	// TODO: apply SCE correction for points.
	_trk_length.push_back(track->Length());
	_trk_start_x.push_back(track->Start().X());
	_trk_start_y.push_back(track->Start().Y());
	_trk_start_z.push_back(track->Start().Z());
	_trk_end_x.push_back(track->End().X());
	_trk_end_y.push_back(track->End().Y());
	_trk_end_z.push_back(track->End().Z());
	_trk_dir_x.push_back(track->StartDirection().X());
	_trk_dir_y.push_back(track->StartDirection().Y());
	_trk_dir_z.push_back(track->StartDirection().Z());
}

void hyperon::HyperonProduction::getShowerVariables(art::Ptr<recob::Shower> &shower)
{	
	if (shower->has_length())
	{
		_shr_length.push_back(shower->Length());
	}
	else
	{
		_shr_length.push_back(bogus::LENGTH);
	}

	// TODO: apply SCE correction for start point.
	_shr_start_x.push_back(shower->ShowerStart().X());
	_shr_start_y.push_back(shower->ShowerStart().Y());
	_shr_start_z.push_back(shower->ShowerStart().Z());
}

// TODO: define and set NULL values for failure modes accessing slice, track, etc..
void hyperon::HyperonProduction::clearTreeVariables()
{
	_n_slices          = 0;
	_n_primary_tracks  = 0;
	_n_primary_showers = 0;

	_mc_nu_pdg = -999;
	_mc_nu_q2  = -1.0;
	_mc_ccnc   = "";
	_mc_mode   = "";

	_pdg.clear();
	_x.clear();
	_y.clear();
	_z.clear();

	_trk_length.clear();
	_trk_shr_score.clear();
	_trk_dir_x.clear();
	_trk_dir_y.clear();
	_trk_dir_z.clear();
	_trk_start_x.clear();
	_trk_start_y.clear();
	_trk_start_z.clear();
	_trk_end_x.clear();
	_trk_end_y.clear();
	_trk_end_z.clear();
	_trk_mean_dedx_plane0.clear();
	_trk_mean_dedx_plane1.clear();
	_trk_mean_dedx_plane2.clear();
	_trk_llrpid.clear();

	_shr_length.clear();
	_shr_dir_x.clear();
	_shr_dir_y.clear();
	_shr_dir_z.clear();
	_shr_start_x.clear();
	_shr_start_y.clear();
	_shr_start_z.clear();
	_shr_energy_plane0.clear();
	_shr_energy_plane1.clear();
	_shr_energy_plane2.clear();
	_shr_dedx_plane0.clear();
	_shr_dedx_plane1.clear();
	_shr_dedx_plane2.clear();
	_shr_open_angle.clear();

	_has_truth.clear();
	_mc_truth_index.clear();
	_trk_true_pdg.clear();
	_trk_true_energy.clear();
	_trk_true_ke.clear();
	_trk_true_px.clear();
	_trk_true_py.clear();
	_trk_true_pz.clear();
	_trk_true_length.clear();
	_trk_true_origin.clear();
	_trk_true_purity.clear();
}

// fillNull is fills the output ttree with null(ish) values when accessing data
// products fails.
void hyperon::HyperonProduction::fillNull()
{
	clearTreeVariables();

	fTree->Fill();
}

DEFINE_ART_MODULE(hyperon::HyperonProduction)
