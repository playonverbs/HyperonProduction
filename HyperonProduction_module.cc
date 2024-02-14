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

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"

#include "ubana/searchingfornues/Selection/CommonDefs/PIDFuncs.h"
#include "ubana/HyperonProduction/Utils.h"
#include "ubana/HyperonProduction/Alg.h"

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
		constexpr double ANGLE  = -999.0;
		constexpr int PDG       = -999;
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

	Atom<std::string> fSliceLabel         { Name("SliceLabel"),
									        Comment("Label for recob::Slice") };
	Atom<std::string> fPFParticleLabel    { Name("PFParticleLabel"),
									        Comment("Label for recob::PFParticle") };
	Atom<std::string> fTrackLabel         { Name("TrackLabel"),
									        Comment("Label for recob::Track") };
	Atom<std::string> fShowerLabel        { Name("ShowerLabel"),
									        Comment("Label for recob::Shower") };
	Atom<std::string> fCaloLabel          { Name("CaloLabel"),
									        Comment("Label for anab::Calorimetry") };
	Atom<std::string> fGeneratorLabel     { Name("GeneratorLabel"),
									        Comment("Label for simb::MCTruth") };
	Atom<std::string> fG4Label  	      { Name("G4Label"),
									        Comment("Label for simb::MCParticle") };
	Atom<std::string> fPIDLabel		      { Name("PIDLabel"),
									        Comment("Label for anab::ParticleID") };
	Atom<std::string> fHitLabel		      { Name("HitLabel"),
										    Comment("Label for recob::Hit") };
	Atom<std::string> fTrackHitAssnsLabel { Name("TrackHitAssnsLabel"),
										    Comment("Label for Assns between recob::Track and recob::Hit") };
	Atom<std::string> fHitTruthAssnsLabel { Name("HitTruthAssnsLabel"),
										    Comment("Label for Assns between MCParticle, Hit and BackTrackerHitMatchingData") };
	Atom<bool>        fIsData             { Name("IsData"),
									        Comment("Flag to indicate if the input is Data") };
	Atom<bool>        fDebug              { Name("Debug"),
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
		void getMCParticleVariables(art::Ptr<simb::MCParticle> &particle, float purity);
		void clearTreeVariables();
		void fillNull();

		void buildTruthHierarchy(std::vector<art::Ptr<simb::MCParticle>> particleVector);
		std::vector<int> getChildIds(const art::Ptr<simb::MCParticle> &p, bool IsNeutron = false);
		int getOrigin(int trackid);

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
		std::string fPIDLabel;
		std::string fHitLabel;
		std::string fTrackHitAssnsLabel;
		std::string fHitTruthAssnsLabel;
		bool fIsData;
		bool fDebug;

		std::map<int, art::Ptr<simb::MCParticle>> particleMap;

		std::vector<int> primary_ids;
		std::vector<int> sigmaZeroDaughter_ids;
		std::vector<int> lambdaDaughter_ids;

		// output tree values here:

		unsigned int _run;
		unsigned int _subrun;
		unsigned int _event;

		int    _mc_nu_pdg;
		double _mc_nu_pos_x;
		double _mc_nu_pos_y;
		double _mc_nu_pos_z;
		double _mc_nu_q2;
		int    _mc_lepton_pdg;
		double _mc_lepton_mom;

		std::string _mc_ccnc;
		std::string _mc_mode;

		unsigned int _n_mctruths;

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
		std::vector<double> _trk_three_plane_dedx;
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
	fPIDLabel(config().fPIDLabel()),
	fHitLabel(config().fHitLabel()),
	fTrackHitAssnsLabel(config().fTrackHitAssnsLabel()),
	fHitTruthAssnsLabel(config().fHitTruthAssnsLabel()),
	fIsData(config().fIsData()),
	fDebug(config().fDebug())
{
	// Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hyperon::HyperonProduction::analyze(art::Event const& e)
{
	// Implementation of required member function here.
	particleMap.clear();
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
			_n_mctruths++;

			_mc_nu_pdg   = truth->GetNeutrino().Nu().PdgCode();
			_mc_nu_q2    = truth->GetNeutrino().QSqr();
			_mc_nu_pos_x = truth->GetNeutrino().Nu().EndX();
			_mc_nu_pos_y = truth->GetNeutrino().Nu().EndY();
			_mc_nu_pos_z = truth->GetNeutrino().Nu().EndZ();

			_mc_lepton_pdg = truth->GetNeutrino().Lepton().PdgCode();
			_mc_lepton_mom = truth->GetNeutrino().Lepton().Momentum().P();

			_mc_ccnc = util::GetCCNC(truth->GetNeutrino().CCNC());
			_mc_mode = util::GetEventType(truth->GetNeutrino().Mode());

			// look at geant4 particles

			const std::vector<art::Ptr<simb::MCParticle>> g4particles =
				util::GetAssocProductVector<simb::MCParticle>(truth, e, fGeneratorLabel, fG4Label);

			for (const art::Ptr<simb::MCParticle> &g4p : g4particles) {
				particleMap.insert(std::make_pair(g4p->TrackId(), g4p));
			}

			buildTruthHierarchy(g4particles);

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

	art::ValidHandle<std::vector<recob::Hit>> hitHandle =
		e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
	art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
		partHitMatchAssoc(hitHandle, e, fHitTruthAssnsLabel);

	art::ValidHandle<std::vector<recob::Track>> trackHandle =
		e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);

	const std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs(slicePFPAssoc.at(nuSliceKey));

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
		art::FindManyP<recob::Hit> trackHitAssoc(trackHandle, e, fTrackHitAssnsLabel);

		if (tracks.size() == 1)
		{
			_n_primary_tracks++;
			art::Ptr<recob::Track> track = tracks.at(0);

			const std::vector<art::Ptr<anab::Calorimetry>> calos =
				util::GetAssocProductVector<anab::Calorimetry>(track, e, fTrackLabel, fCaloLabel);

			auto dedx = alg::ThreePlaneMeandEdX(track, calos);
			
			_trk_mean_dedx_plane0.push_back(dedx.plane0);
			_trk_mean_dedx_plane1.push_back(dedx.plane1);
			_trk_mean_dedx_plane2.push_back(dedx.plane2);
			_trk_three_plane_dedx.push_back(dedx.three_plane_average);

			getTrackVariables(track);

			std::vector<art::Ptr<recob::Hit>> hits = trackHitAssoc.at(track.key());
			
			float pur = -1.0;
			art::Ptr<simb::MCParticle> p = util::getAssocMCParticle(partHitMatchAssoc, hits, pur);
			getMCParticleVariables(p, pur);


		}

		// Handle Showers
		std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(nuSlicePFP.key());

		if (showers.size() == 1)
		{
			_n_primary_showers++;
			art::Ptr<recob::Shower> shower = showers.at(0);

			getShowerVariables(shower);
		}
	} // end nuSlicePFPs loop

	fTree->Fill();
}

void hyperon::HyperonProduction::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	fTree = tfs->make<TTree>("OutputTree", "Output TTree");

	fTree->Branch("run",    &_run);
	fTree->Branch("subrun", &_subrun);
	fTree->Branch("event",  &_event);

	fTree->Branch("mc_nu_pdg",     &_mc_nu_pdg);
	fTree->Branch("mc_nu_q2",      &_mc_nu_q2);
	fTree->Branch("mc_ccnc",       &_mc_ccnc);
	fTree->Branch("mc_mode",       &_mc_mode);
	fTree->Branch("mc_nu_pos_x",   &_mc_nu_pos_x);
	fTree->Branch("mc_nu_pos_y",   &_mc_nu_pos_y);
	fTree->Branch("mc_nu_pos_z",   &_mc_nu_pos_z);
	fTree->Branch("mc_lepton_pdg", &_mc_lepton_pdg);
	fTree->Branch("mc_lepton_mom", &_mc_lepton_mom);

	fTree->Branch("n_mctruths", &_n_mctruths);

	fTree->Branch("pdg",  &_pdg);

	fTree->Branch("n_slices",             & _n_slices);
	fTree->Branch("n_primary_tracks",     & _n_primary_tracks);
	fTree->Branch("trk_length",           & _trk_length);
	fTree->Branch("trk_shr_score",        & _trk_shr_score);
	fTree->Branch("trk_start_x",          & _trk_start_x);
	fTree->Branch("trk_start_y",          & _trk_start_y);
	fTree->Branch("trk_start_z",          & _trk_start_z);
	fTree->Branch("trk_end_x",            & _trk_end_x);
	fTree->Branch("trk_end_y",            & _trk_end_y);
	fTree->Branch("trk_end_z",            & _trk_end_z);
	fTree->Branch("trk_dir_x",            & _trk_dir_x);
	fTree->Branch("trk_dir_y",            & _trk_dir_y);
	fTree->Branch("trk_dir_z",            & _trk_dir_z);
	fTree->Branch("trk_mean_dedx_plane0", & _trk_mean_dedx_plane0);
	fTree->Branch("trk_mean_dedx_plane1", & _trk_mean_dedx_plane1);
	fTree->Branch("trk_mean_dedx_plane2", & _trk_mean_dedx_plane2);
	fTree->Branch("trk_three_plane_dedx", & _trk_three_plane_dedx);

	fTree->Branch("n_primary_showers", &_n_primary_showers);
	fTree->Branch("shr_length",        &_shr_length);
	fTree->Branch("shr_open_angle",    &_shr_open_angle);
	fTree->Branch("shr_start_x",       &_shr_start_x);
	fTree->Branch("shr_start_y",       &_shr_start_y);
	fTree->Branch("shr_start_z",       &_shr_start_z);
	fTree->Branch("shr_dir_x",         &_shr_dir_x);
	fTree->Branch("shr_dir_y",         &_shr_dir_y);
	fTree->Branch("shr_dir_z",         &_shr_dir_z);

	fTree->Branch("trk_true_pdg",    &_trk_true_pdg);
	fTree->Branch("trk_true_energy", &_trk_true_energy);
	fTree->Branch("trk_true_ke",     &_trk_true_ke);
	fTree->Branch("trk_true_px",     &_trk_true_px);
	fTree->Branch("trk_true_py",     &_trk_true_py);
	fTree->Branch("trk_true_pz",     &_trk_true_pz);
	fTree->Branch("trk_true_length", &_trk_true_length);
	fTree->Branch("trk_true_origin", &_trk_true_origin);
	fTree->Branch("trk_true_purity", &_trk_true_purity);
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

	if (shower->has_open_angle())
	{
		_shr_open_angle.push_back(shower->OpenAngle());
	}
	else
	{
		_shr_open_angle.push_back(bogus::ANGLE);
	}

	// TODO: apply SCE correction for start point.
	_shr_start_x.push_back(shower->ShowerStart().X());
	_shr_start_y.push_back(shower->ShowerStart().Y());
	_shr_start_z.push_back(shower->ShowerStart().Z());

	_shr_dir_x.push_back(shower->Direction().X());
	_shr_dir_y.push_back(shower->Direction().Y());
	_shr_dir_z.push_back(shower->Direction().Z());
}

void hyperon::HyperonProduction::getMCParticleVariables(art::Ptr<simb::MCParticle> &particle, float purity)
{
	if (particle.isNonnull()) {
		_trk_true_pdg.push_back(particle->PdgCode());
		_trk_true_energy.push_back(particle->E());
		_trk_true_ke.push_back(particle->T());
		_trk_true_px.push_back(particle->Px());
		_trk_true_py.push_back(particle->Py());
		_trk_true_pz.push_back(particle->Pz());
		// TODO: compute length and origin values.
		_trk_true_length.push_back(0.0);
		_trk_true_origin.push_back(-1);
		_trk_true_purity.push_back(purity);
	}
}

// TODO: define and set NULL values for failure modes accessing slice, track, etc..
void hyperon::HyperonProduction::clearTreeVariables()
{
	_n_mctruths        = 0;
	_n_slices          = 0;
	_n_primary_tracks  = 0;
	_n_primary_showers = 0;

	_mc_nu_pdg   = bogus::PDG;
	_mc_nu_q2    = -1.0;
	_mc_nu_pos_x = bogus::POS;
	_mc_nu_pos_y = bogus::POS;
	_mc_nu_pos_z = bogus::POS;

	_mc_lepton_pdg = bogus::PDG;
	_mc_lepton_mom = -1.0;

	_mc_ccnc   = "";
	_mc_mode   = "";

	primary_ids.clear();
	lambdaDaughter_ids.clear();
	sigmaZeroDaughter_ids.clear();

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
	_trk_three_plane_dedx.clear();
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

void hyperon::HyperonProduction::buildTruthHierarchy(const std::vector<art::Ptr<simb::MCParticle>> particleVector)
{
	for (const art::Ptr<simb::MCParticle> &p : particleVector)
	{
		if (p->Mother() != 0) continue;

		primary_ids.push_back(p->TrackId());

		std::vector<int> _ids = getChildIds(p);

		if (pdg::isHyperon(p)) {
			if (p->PdgCode() == pdg::SigmaZero && p->EndProcess() == "Decay")
				sigmaZeroDaughter_ids.insert(sigmaZeroDaughter_ids.begin(), _ids.begin(), _ids.end());
			else if (p->EndProcess() == "Decay")
				lambdaDaughter_ids.insert(lambdaDaughter_ids.begin(), _ids.begin(), _ids.end());
		}
	}

	for (size_t i_d = 0; i_d < sigmaZeroDaughter_ids.size(); i_d++)
	{
		if (particleMap.find(sigmaZeroDaughter_ids.at(i_d)) == particleMap.end()) continue;

		auto p = particleMap.at(sigmaZeroDaughter_ids.at(i_d));

		if (p->PdgCode() == pdg::Lambda) {
			std::vector<int> _ids = getChildIds(p);
			lambdaDaughter_ids.insert(lambdaDaughter_ids.begin(), _ids.begin(), _ids.end());
		}
	}
}

std::vector<int> hyperon::HyperonProduction::getChildIds(const art::Ptr<simb::MCParticle> &p, bool IsNeutron)
{
	std::vector<int> _decay_ids;

	if (p->EndProcess() != "Decay" && !IsNeutron && !pdg::isKaon(p)) return _decay_ids;

	for (int i_d = 0; i_d < p->NumberDaughters(); i_d++) {
		if (particleMap.find(p->Daughter(i_d)) == particleMap.end()) continue;

		art::Ptr<simb::MCParticle> daughter = particleMap.at(p->Daughter(i_d));

		if (daughter->PdgCode() > 10000) continue;

		if (!util::posMatch(
					TVector3(daughter->Position().X(),
						daughter->Position().Y(),
						daughter->Position().Z()),
					TVector3(p->EndPosition().X(),
						p->EndPosition().Y(),
						p->EndPosition().Z())))
			continue;

		_decay_ids.push_back(p->Daughter(i_d));
	}

	return _decay_ids;
}

int hyperon::HyperonProduction::getOrigin(int trackid)
{
	if (std::find(primary_ids.begin(), primary_ids.end(), trackid) != primary_ids.end()) return 1;
	else if (std::find(lambdaDaughter_ids.begin(), lambdaDaughter_ids.end(), trackid) != lambdaDaughter_ids.end()) return 2;
	else if (std::find(sigmaZeroDaughter_ids.begin(), sigmaZeroDaughter_ids.end(), trackid) != sigmaZeroDaughter_ids.end()) return 5;
	else return 3;
}

DEFINE_ART_MODULE(hyperon::HyperonProduction)
