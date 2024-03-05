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
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"

// larpandora
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "ubana/searchingfornues/Selection/CommonDefs/LLR_PID.h"
#include "ubana/searchingfornues/Selection/CommonDefs/LLRPID_correction_lookup.h"
#include "ubana/searchingfornues/Selection/CommonDefs/LLRPID_proton_muon_lookup.h"
#include "ubana/HyperonProduction/Utils.h"
#include "ubana/HyperonProduction/Alg.h"
#include "ubana/HyperonProduction/FV.h"

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
        constexpr double DOUBLE = -999.0;
        constexpr double LENGTH = -999.0;
        constexpr double POS    = -999.9;
        constexpr double ANGLE  = -999.0;
        constexpr int PDG       = -999;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Config is an art-compatible container for fhicl parameters: this allows for
// parameter validation and descriptions.
struct hyperon::Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    template<typename T>
    using Atom = fhicl::Atom<T>;

    template<typename T, std::size_t SZ>
    using Sequence = fhicl::Sequence<T, SZ>;

    Atom<std::string> fPandoraRecoLabel    { Name("PandoraRecoLabel"),
                                             Comment("Label for pandoraPatRec reconstruction products") };
    Atom<std::string> fFlashMatchRecoLabel { Name("FlashMatchRecoLabel"),
                                             Comment("Label for pandora reconstruction products") };
    Atom<std::string> fTrackLabel          { Name("TrackLabel"),
                                             Comment("Label for recob::Track") };
    Atom<std::string> fShowerLabel         { Name("ShowerLabel"),
                                             Comment("Label for recob::Shower") };
    Atom<std::string> fCaloLabel           { Name("CaloLabel"),
                                             Comment("Label for anab::Calorimetry") };
    Atom<std::string> fGeneratorLabel      { Name("GeneratorLabel"),
                                             Comment("Label for simb::MCTruth") };
    Atom<std::string> fG4Label             { Name("G4Label"),
                                             Comment("Label for simb::MCParticle") };
    Atom<std::string> fPIDLabel            { Name("PIDLabel"),
                                             Comment("Label for anab::ParticleID") };
    Atom<std::string> fHitLabel            { Name("HitLabel"),
                                             Comment("Label for recob::Hit") };
    Atom<std::string> fTrackHitAssnsLabel  { Name("TrackHitAssnsLabel"),
                                             Comment("Label for Assns between recob::Track and recob::Hit") };
    Atom<std::string> fHitTruthAssnsLabel  { Name("HitTruthAssnsLabel"),
                                             Comment("Label for Assns between MCParticle, Hit and BackTrackerHitMatchingData") };
    Atom<bool>        fIsData              { Name("IsData"),
                                             Comment("Flag to indicate if the input is Data") };
    Atom<bool>        fDebug               { Name("Debug"),
                                             Comment("Flag to enable debug messages"),
                                             false };
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        void analyze(art::Event const& evt) override;

        // Selected optional functions.
        void beginJob() override;
        void endJob() override;

    private:
        void fillPandoraMaps(art::Event const& evt);
        void fillMCParticleHitMaps(art::Event const& evt);
        bool isEM(const art::Ptr<simb::MCParticle> &mc_particle);
        int getLeadEMTrackID(const art::Ptr<simb::MCParticle> &mc_particle);
        void getEventMCInfo(art::Event const& evt);
        void getTrueNuSliceID(art::Event const& evt);
        std::vector<art::Ptr<recob::Hit>> collectHitsFromClusters(art::Event const& evt,
            const art::Ptr<recob::PFParticle> &pfparticle);
        void getFlashMatchNuSliceID(art::Event const& evt);
        void getTopologicalScoreNuSliceID(art::Event const& evt);
        void getEventRecoInfo(art::Event const& evt, const int nuSliceID);
        void getPFPRecoInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle);
        void getMCParticleVariables(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle);
        void getBogusMCParticleVariables();
        void getTrackVariables(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle);
        void getBogusTrackVariables();
        void getShowerVariables(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle);
        void getBogusShowerVariables();
        void clearTreeVariables();
        void fillNull();

        void buildTruthHierarchy(std::vector<art::Ptr<simb::MCParticle>> particleVector);
        std::vector<int> getChildIds(const art::Ptr<simb::MCParticle> &p, bool IsNeutron = false);
        int getOrigin(int trackid);

        // Declare member data here.
        searchingfornues::LLRPID llr_pid_calc;
        searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;
        searchingfornues::CorrectionLookUpParameters correction_parameters;

        // FHICL Params
        std::string fPandoraRecoLabel;
        std::string fFlashMatchRecoLabel;
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

        std::vector<int> primary_ids;
        std::vector<int> sigmaZeroDaughter_ids;
        std::vector<int> lambdaDaughter_ids;

        /////////////////////////////
        // Event ID
        /////////////////////////////
        unsigned int _run;
        unsigned int _subrun;
        unsigned int _event;

        /////////////////////////////
        // Event MC Info
        /////////////////////////////
        unsigned int _n_mctruths;
        std::string  _mc_ccnc;
        std::string  _mc_mode;
        int          _mc_nu_pdg;
        double       _mc_nu_pos_x;
        double       _mc_nu_pos_y;
        double       _mc_nu_pos_z;
        double       _mc_nu_q2;
        int          _mc_lepton_pdg;
        double       _mc_lepton_mom;
        int          _true_nu_slice_ID;
        double       _true_nu_slice_completeness;
        double       _true_nu_slice_purity;

        unsigned int _n_slices;
        //unsigned int _n_primary_tracks;          // do we need these if we have multiple slices?
        //unsigned int _n_primary_showers;         // do we need these if we have multiple slices?

        /////////////////////////////
        // FlashMatch Slice Info
        /////////////////////////////
        int _flash_match_nu_slice_ID;

        /////////////////////////////
        // Pandora Slice Info
        /////////////////////////////
        int _pandora_nu_slice_ID;

        /////////////////////////////
        // PFParticle Variables
        /////////////////////////////
        // True stuff
        std::vector<double> _pfp_purity;
        std::vector<double> _pfp_completeness;
        std::vector<bool>   _pfp_has_truth;
        std::vector<int>    _pfp_trackID;
        std::vector<int>    _pfp_true_pdg;
        std::vector<double> _pfp_true_energy;
        std::vector<double> _pfp_true_ke;
        std::vector<double> _pfp_true_px;
        std::vector<double> _pfp_true_py;
        std::vector<double> _pfp_true_pz;
        std::vector<double> _pfp_true_length;
        std::vector<int>    _pfp_true_origin;

        // Reco pfp stuff
        std::vector<int>    _pfp_pdg;
        std::vector<double> _pfp_trk_shr_score;
        std::vector<double> _pfp_x;
        std::vector<double> _pfp_y;
        std::vector<double> _pfp_z;

        /////////////////////////////
        // Track Variables
        /////////////////////////////
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
        std::vector<double> _trk_three_plane_dedx;
        std::vector<double> _trk_llrpid; // TODO: Set this!

        /////////////////////////////
        // Shower Variables
        /////////////////////////////
        std::vector<double> _shr_length;
        std::vector<double> _shr_dir_x;
        std::vector<double> _shr_dir_y;
        std::vector<double> _shr_dir_z;
        std::vector<double> _shr_start_x;
        std::vector<double> _shr_start_y;
        std::vector<double> _shr_start_z;
        std::vector<double> _shr_energy_plane0; // TODO: Set this!
        std::vector<double> _shr_energy_plane1; // TODO: Set this!
        std::vector<double> _shr_energy_plane2; // TODO: Set this!
        std::vector<double> _shr_dedx_plane0; // TODO: Set this!
        std::vector<double> _shr_dedx_plane1; // TODO: Set this!
        std::vector<double> _shr_dedx_plane2; // TODO: Set this!
        std::vector<double> _shr_open_angle;

        /////////////////////////////
        // Tree
        /////////////////////////////
        TTree* fTree;

        /////////////////////////////
        // Internal analyzer maps
        /////////////////////////////
        std::map<int, int> _hit_to_trackID;                 // Linking hit -> trackID of true owner
        std::map<int, art::Ptr<recob::Slice>> _slice_map;   // Linking sliceID -> slice
        std::map<int, std::vector<int>> _trackID_to_hits;   // Linking trackID -> nTrueHits
        lar_pandora::MCParticleMap _mc_particle_map;        // Linking TrackID -> MCParticle
        lar_pandora::PFParticleMap _pfp_map;                // Linking Self() -> PFParticle
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

hyperon::HyperonProduction::HyperonProduction(Parameters const& config)
    : EDAnalyzer{config},
    fPandoraRecoLabel(config().fPandoraRecoLabel()),
    fFlashMatchRecoLabel(config().fFlashMatchRecoLabel()),
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
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::analyze(art::Event const& evt)
{
    // Make sure our true->reco maps are clear
    _hit_to_trackID.clear();
    _trackID_to_hits.clear();
    _mc_particle_map.clear();
    _pfp_map.clear();

    // And tree variables...
    clearTreeVariables();

    // Need to fill some 'Pandora maps' to assist true->reco matching
    if (fDebug) std::cout << "Filling Pandora Maps..." << std::endl;
    fillPandoraMaps(evt);
    if (fDebug) std::cout << "Filling MCParticle Hit Maps..." << std::endl;
    fillMCParticleHitMaps(evt);

    // Fill those branches!
    _run    = evt.run();
    _subrun = evt.subRun();
    _event  = evt.id().event();

    // Fill the MC info for the event
    if (fDebug) std::cout << "Filling MC Event Info..." << std::endl;
    getEventMCInfo(evt);

    // Find the flash matched neutrino slice (if one exists)
    if (fDebug) std::cout << "Getting flash match slice ID..." << std::endl;
    getFlashMatchNuSliceID(evt);

    // Find the Pandora (highest topological score) neutrino slice (if one exists)
    if (fDebug) std::cout << "Getting the highest topological score slice ID..." << std::endl;
    getTopologicalScoreNuSliceID(evt);

    // Fill the reconstructed neutrino hierarchy variables (using flash match neutrino slice)
    if (fDebug) std::cout << "Filling Reconstructed Particle Variables..." << std::endl;
    getEventRecoInfo(evt, _flash_match_nu_slice_ID);

    fTree->Fill();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::beginJob()
{
    // Setup searchingfornues::LLRPID with parameters.
    llr_pid_calc.set_dedx_binning(0,      protonmuon_parameters.dedx_edges_pl_0);
    llr_pid_calc.set_par_binning(0,       protonmuon_parameters.parameters_edges_pl_0);
    llr_pid_calc.set_lookup_tables(0,     protonmuon_parameters.dedx_pdf_pl_0);
    llr_pid_calc.set_corr_par_binning(0,  correction_parameters.parameter_correction_edges_pl_0);
    llr_pid_calc.set_correction_tables(0, correction_parameters.correction_table_pl_0);

    llr_pid_calc.set_dedx_binning(1,      protonmuon_parameters.dedx_edges_pl_1);
    llr_pid_calc.set_par_binning(1,       protonmuon_parameters.parameters_edges_pl_1);
    llr_pid_calc.set_lookup_tables(1,     protonmuon_parameters.dedx_pdf_pl_1);
    llr_pid_calc.set_corr_par_binning(1,  correction_parameters.parameter_correction_edges_pl_1);
    llr_pid_calc.set_correction_tables(1, correction_parameters.correction_table_pl_1);

    llr_pid_calc.set_dedx_binning(2,      protonmuon_parameters.dedx_edges_pl_2);
    llr_pid_calc.set_par_binning(2,       protonmuon_parameters.parameters_edges_pl_2);
    llr_pid_calc.set_lookup_tables(2,     protonmuon_parameters.dedx_pdf_pl_2);
    llr_pid_calc.set_corr_par_binning(2,  correction_parameters.parameter_correction_edges_pl_2);
    llr_pid_calc.set_correction_tables(2, correction_parameters.correction_table_pl_2);

    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("OutputTree", "Output TTree");

    /////////////////////////////
    // Event ID
    /////////////////////////////
    fTree->Branch("run",    &_run);
    fTree->Branch("subrun", &_subrun);
    fTree->Branch("event",  &_event);

    /////////////////////////////
    // Event MC Info
    /////////////////////////////
    fTree->Branch("n_mctruths",                  &_n_mctruths);
    fTree->Branch("mc_nu_pdg",                   &_mc_nu_pdg);
    fTree->Branch("mc_nu_q2",                    &_mc_nu_q2);
    fTree->Branch("mc_ccnc",                     &_mc_ccnc);
    fTree->Branch("mc_mode",                     &_mc_mode);
    fTree->Branch("mc_nu_pos_x",                 &_mc_nu_pos_x);
    fTree->Branch("mc_nu_pos_y",                 &_mc_nu_pos_y);
    fTree->Branch("mc_nu_pos_z",                 &_mc_nu_pos_z);
    fTree->Branch("mc_lepton_pdg",               &_mc_lepton_pdg);
    fTree->Branch("mc_lepton_mom",               &_mc_lepton_mom);
    fTree->Branch("true_nu_slice_ID",            &_true_nu_slice_ID);
    fTree->Branch("true_nu_slice_completeness",  &_true_nu_slice_completeness);
    fTree->Branch("true_nu_slice_purity",        &_true_nu_slice_purity);
    fTree->Branch("n_slices",                    &_n_slices);

    /////////////////////////////
    // FlashMatch Slice Info
    /////////////////////////////
    fTree->Branch("flash_match_nu_slice_ID",     &_flash_match_nu_slice_ID);

    /////////////////////////////
    // Pandora Slice Info
    /////////////////////////////
    fTree->Branch("pandora_nu_slice_ID",         &_pandora_nu_slice_ID);

    /////////////////////////////
    // PFParticle Variables
    /////////////////////////////
    // True stuff
    fTree->Branch("pfp_purity",        &_pfp_purity);
    fTree->Branch("pfp_completeness",  &_pfp_completeness);
    fTree->Branch("pfp_has_truth",     &_pfp_has_truth);
    fTree->Branch("pfp_trackID",       &_pfp_trackID);
    fTree->Branch("pfp_true_pdg",      &_pfp_true_pdg);
    fTree->Branch("pfp_true_energy",   &_pfp_true_energy);
    fTree->Branch("pfp_true_ke",       &_pfp_true_ke);
    fTree->Branch("pfp_true_px",       &_pfp_true_px);
    fTree->Branch("pfp_true_py",       &_pfp_true_py);
    fTree->Branch("pfp_true_pz",       &_pfp_true_pz);
    fTree->Branch("pfp_true_length",   &_pfp_true_length);
    fTree->Branch("pfp_true_origin",   &_pfp_true_origin);

    // Reco pfp stuff
    fTree->Branch("pfp_pdg",            &_pfp_pdg);
    fTree->Branch("pfp_trk_shr_score",  &_pfp_trk_shr_score);
    fTree->Branch("pfp_x",              &_pfp_x);
    fTree->Branch("pfp_y",              &_pfp_y);
    fTree->Branch("pfp_z",              &_pfp_z);

    /////////////////////////////
    // Track Variables
    /////////////////////////////
    fTree->Branch("trk_length",              &_trk_length);
    fTree->Branch("trk_start_x",             &_trk_start_x);
    fTree->Branch("trk_start_y",             &_trk_start_y);
    fTree->Branch("trk_start_z",             &_trk_start_z);
    fTree->Branch("trk_end_x",               &_trk_end_x);
    fTree->Branch("trk_end_y",               &_trk_end_y);
    fTree->Branch("trk_end_z",               &_trk_end_z);
    fTree->Branch("trk_dir_x",               &_trk_dir_x);
    fTree->Branch("trk_dir_y",               &_trk_dir_y);
    fTree->Branch("trk_dir_z",               &_trk_dir_z);
    fTree->Branch("trk_mean_dedx_plane0",    &_trk_mean_dedx_plane0);
    fTree->Branch("trk_mean_dedx_plane1",    &_trk_mean_dedx_plane1);
    fTree->Branch("trk_mean_dedx_plane2",    &_trk_mean_dedx_plane2);
    fTree->Branch("trk_three_plane_dedx",    &_trk_three_plane_dedx);
    fTree->Branch("trk_llrpid",              &_trk_llrpid);

    /////////////////////////////
    // Shower Variables
    /////////////////////////////
    fTree->Branch("shr_length",         &_shr_length);
    fTree->Branch("shr_open_angle",     &_shr_open_angle);
    fTree->Branch("shr_start_x",        &_shr_start_x);
    fTree->Branch("shr_start_y",        &_shr_start_y);
    fTree->Branch("shr_start_z",        &_shr_start_z);
    fTree->Branch("shr_dir_x",          &_shr_dir_x);
    fTree->Branch("shr_dir_y",          &_shr_dir_y);
    fTree->Branch("shr_dir_z",          &_shr_dir_z);
    fTree->Branch("shr_energy_plane0",  &_shr_energy_plane0);
    fTree->Branch("shr_energy_plane1",  &_shr_energy_plane1);
    fTree->Branch("shr_energy_plane2",  &_shr_energy_plane2);
    fTree->Branch("shr_dedx_plane0",    &_shr_dedx_plane0);
    fTree->Branch("shr_dedx_plane1",    &_shr_dedx_plane1);
    fTree->Branch("shr_dedx_plane2",    &_shr_dedx_plane2);

    //fTree->Branch("n_primary_tracks",        & _n_primary_tracks);
    //fTree->Branch("n_primary_showers", &_n_primary_showers);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::endJob()
{
    if (fDebug)
        FNLOG("ending job");
    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::fillPandoraMaps(art::Event const& evt)
{
    // MCParticle map
    art::Handle<std::vector<simb::MCParticle>> mc_particle_handle;
    std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;

    if (!evt.getByLabel(fG4Label, mc_particle_handle))
        throw cet::exception("HyperonProduction::fillPandoraMaps") << "No MCParticle Data Products Found! :(" << std::endl;

    art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, _mc_particle_map);

    // PFParticle map
    art::Handle<std::vector<recob::PFParticle>> pfp_handle;
    std::vector<art::Ptr<recob::PFParticle>> pfp_vector;

    if (!evt.getByLabel(fPandoraRecoLabel, pfp_handle))
        throw cet::exception("HyperonProduction::fillPandoraMaps") << "No PFParticle Data Products Found! :(" << std::endl;

    art::fill_ptr_vector(pfp_vector, pfp_handle);

    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfp_vector, _pfp_map);

    // Slice map
    art::Handle<std::vector<recob::Slice>> slice_handle;
    std::vector<art::Ptr<recob::Slice>> slice_vector;

    if (!evt.getByLabel(fPandoraRecoLabel, slice_handle))
        throw cet::exception("HyperonProduction::fillPandoraMaps") << "No Slice Data Products Found! :(" << std::endl;

    art::fill_ptr_vector(slice_vector, slice_handle);

    for (const art::Ptr<recob::Slice> &slice : slice_vector)
        _slice_map[slice->ID()] = slice;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::fillMCParticleHitMaps(art::Event const& evt)
{
    // Get event hits
    art::Handle<std::vector<recob::Hit>> hit_handle;
    std::vector<art::Ptr<recob::Hit>> hit_vector;

    if (!evt.getByLabel(fHitLabel, hit_handle))
        throw cet::exception("HyperonProduction::fillMCParticleHitMaps") << "No Hit Data Products Found! :(" << std::endl;

    art::fill_ptr_vector(hit_vector, hit_handle);

    // Get backtracker info
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_particle =
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hit_handle, evt, fHitTruthAssnsLabel);

    // Truth match
    for (unsigned int hit_index = 0; hit_index < hit_vector.size(); hit_index++)
    {
        const art::Ptr<recob::Hit> &hit = hit_vector.at(hit_index);
        const std::vector<art::Ptr<simb::MCParticle>> &matched_mc_particles = assoc_mc_particle.at(hit.key());
        auto matched_datas = assoc_mc_particle.data(hit.key());

        for (unsigned int mc_particle_index = 0; mc_particle_index < matched_mc_particles.size(); ++mc_particle_index)
        {
            const art::Ptr<simb::MCParticle> &matched_mc_particle = matched_mc_particles.at(mc_particle_index);
            auto matched_data = matched_datas.at(mc_particle_index);

            if (matched_data->isMaxIDE != 1)
                continue;

            // Get trackID (EM showers need to be folded back)
            const int trackID = isEM(matched_mc_particle) ? getLeadEMTrackID(matched_mc_particle) : matched_mc_particle->TrackId();

            _hit_to_trackID[hit.key()] = trackID;
            _trackID_to_hits[trackID].push_back(hit.key());
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::HyperonProduction::isEM(const art::Ptr<simb::MCParticle> &mc_particle)
{
    return ((std::abs(mc_particle->PdgCode()) == 11) || (mc_particle->PdgCode() == 22));
}

///////////////////////////////////////////////////////////////////////////////////////////

// If its an EM particle, we have to move up the EM hierarchy
int hyperon::HyperonProduction::getLeadEMTrackID(const art::Ptr<simb::MCParticle> &mc_particle)
{
    int trackID = mc_particle->TrackId();
    art::Ptr<simb::MCParticle> parent_mc_particle = mc_particle;

    do
    {
        trackID = parent_mc_particle->TrackId();
        const int parentID = parent_mc_particle->Mother();

        if (_mc_particle_map.find(parentID) == _mc_particle_map.end())
            break;

        parent_mc_particle = _mc_particle_map.at(parentID);
    }
    while (isEM(parent_mc_particle));

    return trackID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getEventMCInfo(art::Event const& evt)
{
    // Check if this is an MC file.
    if (fIsData)
        return;

    art::ValidHandle<std::vector<simb::MCTruth>> mcTruthHandle =
        evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);
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
            util::GetAssocProductVector<simb::MCParticle>(truth, evt, fGeneratorLabel, fG4Label);

        buildTruthHierarchy(g4particles);

        break; // escape after first MCNeutrino
    }

    if (fDebug) std::cout << "Filling MC Slice Info..." << std::endl;
    getTrueNuSliceID(evt);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getTrueNuSliceID(art::Event const& evt)
{
    // Get slice information
    art::ValidHandle<std::vector<recob::Slice>> slice_handle =
        evt.getValidHandle<std::vector<recob::Slice>>(fPandoraRecoLabel);
    std::vector<art::Ptr<recob::Slice>> slice_vector;

    if (slice_handle.isValid())
        art::fill_ptr_vector(slice_vector, slice_handle);

    _n_slices = slice_vector.size();

    art::FindManyP<recob::Hit> hit_assoc = art::FindManyP<recob::Hit>(slice_handle, evt, fPandoraRecoLabel);

    // Now find true nu slice ID
    int true_slice_n_hits(-1);
    int true_slice_n_signal_hits(-1);
    std::map<int, int> slice_signal_hit_map;
    unsigned int total_signal_hits(0);

    for (art::Ptr<recob::Slice> &slice : slice_vector)
    {
        slice_signal_hit_map[slice->ID()] = 0;

        const std::vector<art::Ptr<recob::Hit>> &slice_hits(hit_assoc.at(slice.key()));

        for (const art::Ptr<recob::Hit> &slice_hit : slice_hits)
        {
            // Hits from cosmics won't have an associated trackID
            if (_hit_to_trackID.find(slice_hit.key()) == _hit_to_trackID.end())
                continue;

            ++slice_signal_hit_map[slice->ID()];
            ++total_signal_hits;
        }

        const int slice_n_signal_hits = slice_signal_hit_map[slice->ID()];

        if ((slice_n_signal_hits > true_slice_n_hits) && (slice_n_signal_hits > 0))
        {
            true_slice_n_signal_hits = slice_n_signal_hits;
            true_slice_n_hits = slice_hits.size();
            _true_nu_slice_ID = slice->ID();
        }
    }

    // Calculate slice completeness and purity
    if (true_slice_n_signal_hits > 0)
    {
        _true_nu_slice_completeness = static_cast<double>(true_slice_n_signal_hits) / static_cast<double>(total_signal_hits);
        _true_nu_slice_purity = static_cast<double>(true_slice_n_signal_hits) / static_cast<double>(true_slice_n_hits);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

std::vector<art::Ptr<recob::Hit>> hyperon::HyperonProduction::collectHitsFromClusters(art::Event const& evt,
    const art::Ptr<recob::PFParticle> &pfparticle)
{
    std::vector<art::Ptr<recob::Hit>> hits;

    art::Handle<std::vector<recob::PFParticle>> pfp_handle;

    if (!evt.getByLabel(fPandoraRecoLabel, pfp_handle))
        throw cet::exception("HyperonProduction::CollectHitsFromClusters") << "No PFParticle Data Products Found! :(" << std::endl;

    art::Handle<std::vector<recob::Cluster>> cluster_handle;

    if (!evt.getByLabel(fPandoraRecoLabel, cluster_handle))
        throw cet::exception("HyperonProduction::CollectHitsFromClusters") << "No Cluster Data Products Found! :(" << std::endl;

    art::FindManyP<recob::Cluster> pfp_clusters_assoc = art::FindManyP<recob::Cluster>(pfp_handle, evt, fPandoraRecoLabel);
    art::FindManyP<recob::Hit> cluster_hit_assoc = art::FindManyP<recob::Hit>(cluster_handle, evt, fPandoraRecoLabel);

    std::vector<art::Ptr<recob::Cluster>> clusters = pfp_clusters_assoc.at(pfparticle.key());

    for (const art::Ptr<recob::Cluster> cluster : clusters)
    {
        std::vector<art::Ptr<recob::Hit>> cluster_hits = cluster_hit_assoc.at(cluster.key());
        hits.insert(hits.end(), cluster_hits.begin(), cluster_hits.end());
    }

    return hits;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getFlashMatchNuSliceID(art::Event const& evt)
{
    art::Handle<std::vector<recob::PFParticle>> fm_pfp_handle;
    std::vector<art::Ptr<recob::PFParticle>> fm_pfp_vector;

    if (!evt.getByLabel(fFlashMatchRecoLabel, fm_pfp_handle))
        throw cet::exception("HyperonProduction::getFlashMatchNuSliceID") << "No PFParticle Data Products Found! :(" << std::endl;

    art::fill_ptr_vector(fm_pfp_vector, fm_pfp_handle);

    std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(fm_pfp_vector, neutrinoPFPs);

    if (neutrinoPFPs.size() > 1)
    {
        throw cet::exception("HyperonProduction::getFlashMatchNuSliceID") << "Too many neutrinos found!" << std::endl;
    }
    else if (neutrinoPFPs.size() == 1)
    {
        art::FindManyP<recob::Slice> fm_slice_assoc = art::FindManyP<recob::Slice>(fm_pfp_handle, evt, fFlashMatchRecoLabel);
        const std::vector<art::Ptr<recob::Slice>> &fm_slices = fm_slice_assoc.at(neutrinoPFPs[0].key());

        if (fm_slices.empty())
            return;

        // Get hits, and then work out if they also live in the pandoraPatRec reco?
        const art::Ptr<recob::Slice> &fm_slice(fm_slices.at(0));

        art::Handle<std::vector<recob::Slice>> fm_slice_handle;

        if (!evt.getByLabel(fFlashMatchRecoLabel, fm_slice_handle))
            throw cet::exception("HyperonProduction::getFlashMatchNuSliceID") << "No Flash Match Slice Data Products Found!" << std::endl;

        art::FindManyP<recob::Hit> fm_hit_assoc = art::FindManyP<recob::Hit>(fm_slice_handle, evt, fFlashMatchRecoLabel);
        const std::vector<art::Ptr<recob::Hit>> &fm_slice_hits(fm_hit_assoc.at(fm_slice.key()));

        if (fm_slice_hits.empty())
            return;

        // Get equivalent pandoraPatRec products
        art::ValidHandle<std::vector<recob::Slice>> slice_handle =
            evt.getValidHandle<std::vector<recob::Slice>>(fPandoraRecoLabel);

        std::vector<art::Ptr<recob::Slice>> slice_vector;

        if (slice_handle.isValid())
            art::fill_ptr_vector(slice_vector, slice_handle);

        art::FindManyP<recob::Hit> hit_assoc = art::FindManyP<recob::Hit>(slice_handle, evt, fPandoraRecoLabel);

        // Loop through pandoraPatRecSlices
        for (art::Ptr<recob::Slice> &slice : slice_vector)
        {
            const std::vector<art::Ptr<recob::Hit>> &slice_hits(hit_assoc.at(slice.key()));

            // could loop over whichever is smaller - to save time...
            const std::vector<art::Ptr<recob::Hit>> &hitsToLoop(slice_hits.size() < fm_slice_hits.size() ? slice_hits : fm_slice_hits);
            const art::Ptr<recob::Hit> &hitToCompare(slice_hits.size() < fm_slice_hits.size() ? fm_slice_hits.front() : slice_hits.front());

            for (const art::Ptr<recob::Hit> &hit : hitsToLoop)
            {
                if (hit.key() == hitToCompare.key())
                {
                    _flash_match_nu_slice_ID = slice->ID();
                    return;
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getTopologicalScoreNuSliceID(art::Event const& evt)
{
    art::ValidHandle<std::vector<recob::Slice>> slice_handle =
        evt.getValidHandle<std::vector<recob::Slice>>(fPandoraRecoLabel);

    std::vector<art::Ptr<recob::Slice>> slice_vector;

    if (slice_handle.isValid())
        art::fill_ptr_vector(slice_vector, slice_handle);

    art::ValidHandle<std::vector<recob::PFParticle>> pfp_handle =
        evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraRecoLabel);

    art::FindManyP<recob::PFParticle> pfp_assoc = art::FindManyP<recob::PFParticle>(slice_handle, evt, fPandoraRecoLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> metadata_assn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfp_handle, evt, fPandoraRecoLabel);

    double best_topological_score(-std::numeric_limits<double>::max());

    for (const art::Ptr<recob::Slice> &slice : slice_vector)
    {
        const std::vector<art::Ptr<recob::PFParticle>> slicePFPs = pfp_assoc.at(slice.key());

        for (const art::Ptr<recob::PFParticle> &pfp : slicePFPs)
        {
            // only topological score for the primary pfp
            if (!pfp->IsPrimary())
                continue;

            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfp_meta = metadata_assn.at(pfp.key());

            if (pfp_meta.empty())
                continue;

            const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties_map(pfp_meta.at(0)->GetPropertiesMap());

            if (pfp_properties_map.find("NuScore") != pfp_properties_map.end())
            {
                const double this_topological_score = pfp_properties_map.at("NuScore");

                if (this_topological_score > best_topological_score)
                {
                    best_topological_score = this_topological_score;
                    _pandora_nu_slice_ID = slice->ID();
                }
            }
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getEventRecoInfo(art::Event const& evt, const int nu_sliceID)
{
    if ((nu_sliceID < 0) || (_slice_map.find(nu_sliceID) == _slice_map.end()))
    {
        if (fDebug)
            FNLOG("flash match, or pandora slice not found");

        fillNull();
        return;
    }

    // Get our products
    art::ValidHandle<std::vector<recob::Slice>> slice_handle =
        evt.getValidHandle<std::vector<recob::Slice>>(fPandoraRecoLabel);
    art::FindManyP<recob::PFParticle> slice_pfp_assoc(slice_handle, evt, fPandoraRecoLabel);

    const std::vector<art::Ptr<recob::PFParticle>> nu_slice_pfps(slice_pfp_assoc.at(_slice_map.at(nu_sliceID).key()));

    for (const art::Ptr<recob::PFParticle> &nu_slice_pfp : nu_slice_pfps)
    {
        // Is the PFP in the neutrino hierarchy?
        // This is important, so we don't pick up the CR hypothesis output
        const art::Ptr<recob::PFParticle> parentPFP = lar_pandora::LArPandoraHelper::GetParentPFParticle(_pfp_map, nu_slice_pfp);

        if (!lar_pandora::LArPandoraHelper::IsNeutrino(parentPFP))
            continue;

        // Let's just looking at primaries
        const unsigned int generation = lar_pandora::LArPandoraHelper::GetGeneration(_pfp_map, nu_slice_pfp);

        if (generation != 2)
            continue;

        // Get PFP reco info
        getPFPRecoInfo(evt, nu_slice_pfp);

        // Get PFP truth matching variables
        if (!fIsData)
            getMCParticleVariables(evt, nu_slice_pfp);
        else
            getBogusMCParticleVariables();

        // Handle Tracks
        getTrackVariables(evt, nu_slice_pfp);

        // Handle Showers
        getShowerVariables(evt, nu_slice_pfp);

    } // end nu_slice_pfps loop
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getPFPRecoInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    // Pandora 'PDG'
    _pfp_pdg.push_back(pfparticle->PdgCode());

    // Track/shower score
    art::ValidHandle<std::vector<recob::PFParticle>> pfp_handle =
        evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraRecoLabel);

    art::FindManyP<larpandoraobj::PFParticleMetadata>
        pfp_meta_assoc(pfp_handle, evt, fPandoraRecoLabel);

    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfp_metas =
        pfp_meta_assoc.at(pfparticle.key());

    if ((!pfp_metas.empty()) && (pfp_metas.at(0)->GetPropertiesMap().find("TrackScore") !=
        pfp_metas.at(0)->GetPropertiesMap().end()))
    {
        _pfp_trk_shr_score.push_back(pfp_metas.at(0)->GetPropertiesMap().at("TrackScore"));
    }
    else
    {
        _pfp_trk_shr_score.push_back(bogus::DOUBLE);
    }

    // Vertex
    art::FindManyP<recob::Vertex>
        pfp_vertex_assoc(pfp_handle, evt, fPandoraRecoLabel);

    std::vector<art::Ptr<recob::Vertex>> pfp_vertices =
        pfp_vertex_assoc.at(pfparticle.key());

    if (!pfp_vertices.empty())
    {
        const art::Ptr<recob::Vertex> &vertex = pfp_vertices.at(0);

        _pfp_x.push_back(vertex->position().X());
        _pfp_y.push_back(vertex->position().Y());
        _pfp_z.push_back(vertex->position().Z());
    }
    else
    {
        _pfp_x.push_back(bogus::DOUBLE);
        _pfp_y.push_back(bogus::DOUBLE);
        _pfp_z.push_back(bogus::DOUBLE);
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getMCParticleVariables(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    std::vector<art::Ptr<recob::Hit>> pfp_hits = collectHitsFromClusters(evt, pfparticle);

    std::map<int, int> counting_map;
    for (const art::Ptr<recob::Hit> pfp_hit : pfp_hits)
     {
         if (_hit_to_trackID.find(pfp_hit.key()) == _hit_to_trackID.end())
             continue;

         const int trackID(_hit_to_trackID.at(pfp_hit.key()));

         if (counting_map.find(trackID) == counting_map.end())
             counting_map[trackID] = 1;
         else
             ++counting_map[trackID];
     }

    int max_hits = -1;
    int matched_trackID = -1;

    // pragmatic tie-breaker
    for (auto &entry : counting_map)
    {
        if ((entry.second > max_hits) || ((entry.second == max_hits) && (entry.first > matched_trackID)))
        {
            max_hits = entry.second;
            matched_trackID = entry.first;
        }
    }

    if (_mc_particle_map.find(matched_trackID) != _mc_particle_map.end())
    {
        const art::Ptr<simb::MCParticle> &matched_mc_particle(_mc_particle_map.at(matched_trackID));

        _pfp_has_truth.push_back(true);
        _pfp_trackID.push_back(matched_mc_particle->TrackId());
        _pfp_true_pdg.push_back(matched_mc_particle->PdgCode());
        _pfp_true_energy.push_back(matched_mc_particle->E());
        _pfp_true_ke.push_back(matched_mc_particle->T());
        _pfp_true_px.push_back(matched_mc_particle->Px());
        _pfp_true_py.push_back(matched_mc_particle->Py());
        _pfp_true_pz.push_back(matched_mc_particle->Pz());
        // TODO: compute length value.
        _pfp_true_length.push_back(0.0);
        _pfp_true_origin.push_back(getOrigin(matched_mc_particle->TrackId()));

        const int n_true_hits = _trackID_to_hits.at(matched_trackID).size();

        _pfp_completeness.push_back(static_cast<double>(max_hits) / static_cast<double>(n_true_hits));
        _pfp_purity.push_back(static_cast<double>(max_hits) / pfp_hits.size());
    }
    else
    {
        getBogusMCParticleVariables();
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getBogusMCParticleVariables()
{
    // fill with bogus values...
    _pfp_has_truth.push_back(false);
    _pfp_trackID.push_back(-1);
    _pfp_true_pdg.push_back(-1);
    _pfp_true_energy.push_back(bogus::DOUBLE);
    _pfp_true_ke.push_back(bogus::DOUBLE);
    _pfp_true_px.push_back(bogus::DOUBLE);
    _pfp_true_py.push_back(bogus::DOUBLE);
    _pfp_true_pz.push_back(bogus::DOUBLE);
    _pfp_true_length.push_back(bogus::DOUBLE);
    _pfp_true_origin.push_back(-1);
    _pfp_completeness.push_back(bogus::DOUBLE);
    _pfp_purity.push_back(bogus::DOUBLE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getTrackVariables(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    art::ValidHandle<std::vector<recob::PFParticle>> pfp_handle =
        evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraRecoLabel);

    art::FindManyP<recob::Track>  pfp_track_assoc(pfp_handle, evt, fTrackLabel);

    std::vector<art::Ptr<recob::Track>> tracks = pfp_track_assoc.at(pfparticle.key());

    if (tracks.size() == 0)
    {
        getBogusTrackVariables();
        return;
    }

    art::Ptr<recob::Track> track = tracks.at(0);

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

    // Calorimetry
    const std::vector<art::Ptr<anab::Calorimetry>> calos =
        util::GetAssocProductVector<anab::Calorimetry>(track, evt, fTrackLabel, fCaloLabel);

    auto dedx = alg::ThreePlaneMeandEdX(track, calos);

    _trk_mean_dedx_plane0.push_back(dedx.plane0);
    _trk_mean_dedx_plane1.push_back(dedx.plane1);
    _trk_mean_dedx_plane2.push_back(dedx.plane2);
    _trk_three_plane_dedx.push_back(dedx.three_plane_average);

    _trk_llrpid.push_back(alg::LLRPID(calos, llr_pid_calc));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getBogusTrackVariables()
{
    _trk_length.push_back(bogus::DOUBLE);
    _trk_start_x.push_back(bogus::DOUBLE);
    _trk_start_y.push_back(bogus::DOUBLE);
    _trk_start_z.push_back(bogus::DOUBLE);
    _trk_end_x.push_back(bogus::DOUBLE);
    _trk_end_y.push_back(bogus::DOUBLE);
    _trk_end_z.push_back(bogus::DOUBLE);
    _trk_dir_x.push_back(bogus::DOUBLE);
    _trk_dir_y.push_back(bogus::DOUBLE);
    _trk_dir_z.push_back(bogus::DOUBLE);
    _trk_mean_dedx_plane0.push_back(bogus::DOUBLE);
    _trk_mean_dedx_plane1.push_back(bogus::DOUBLE);
    _trk_mean_dedx_plane2.push_back(bogus::DOUBLE);
    _trk_three_plane_dedx.push_back(bogus::DOUBLE);
    _trk_llrpid.push_back(bogus::DOUBLE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getShowerVariables(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    art::ValidHandle<std::vector<recob::PFParticle>> pfp_handle =
        evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraRecoLabel);

    art::FindManyP<recob::Shower>  pfp_shower_assoc(pfp_handle, evt, fShowerLabel);

    std::vector<art::Ptr<recob::Shower>> showers = pfp_shower_assoc.at(pfparticle.key());

    if (showers.size() == 0)
    {
        getBogusShowerVariables();
        return;
    }

    art::Ptr<recob::Shower> shower = showers.at(0);

    if (shower->has_length())
        _shr_length.push_back(shower->Length());
    else
        _shr_length.push_back(bogus::LENGTH);

    if (shower->has_open_angle())
        _shr_open_angle.push_back(shower->OpenAngle());
    else
        _shr_open_angle.push_back(bogus::ANGLE);

    // TODO: apply SCE correction for start point.
    _shr_start_x.push_back(shower->ShowerStart().X());
    _shr_start_y.push_back(shower->ShowerStart().Y());
    _shr_start_z.push_back(shower->ShowerStart().Z());
    _shr_dir_x.push_back(shower->Direction().X());
    _shr_dir_y.push_back(shower->Direction().Y());
    _shr_dir_z.push_back(shower->Direction().Z());
    _shr_energy_plane0.push_back(bogus::DOUBLE);
    _shr_energy_plane1.push_back(bogus::DOUBLE);
    _shr_energy_plane2.push_back(bogus::DOUBLE);
    _shr_dedx_plane0.push_back(bogus::DOUBLE);
    _shr_dedx_plane1.push_back(bogus::DOUBLE);
    _shr_dedx_plane2.push_back(bogus::DOUBLE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::HyperonProduction::getBogusShowerVariables()
{
    _shr_length.push_back(bogus::LENGTH);
    _shr_open_angle.push_back(bogus::ANGLE);
    _shr_start_x.push_back(bogus::DOUBLE);
    _shr_start_y.push_back(bogus::DOUBLE);
    _shr_start_z.push_back(bogus::DOUBLE);
    _shr_dir_x.push_back(bogus::DOUBLE);
    _shr_dir_y.push_back(bogus::DOUBLE);
    _shr_dir_z.push_back(bogus::DOUBLE);
    _shr_energy_plane0.push_back(bogus::DOUBLE);
    _shr_energy_plane1.push_back(bogus::DOUBLE);
    _shr_energy_plane2.push_back(bogus::DOUBLE);
    _shr_dedx_plane0.push_back(bogus::DOUBLE);
    _shr_dedx_plane1.push_back(bogus::DOUBLE);
    _shr_dedx_plane2.push_back(bogus::DOUBLE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: define and set NULL values for failure modes accessing slice, track, etc..
void hyperon::HyperonProduction::clearTreeVariables()
{
    // ATTN: I DONT THINK THEY SHOULD BE HERE
    //i.e. this function is called when there is no reco info
    primary_ids.clear();
    lambdaDaughter_ids.clear();
    sigmaZeroDaughter_ids.clear();

    /////////////////////////////
    // Event ID
    /////////////////////////////
    _run = -999;
    _subrun = -999;
    _event = -999;

    /////////////////////////////
    // Event MC Info
    /////////////////////////////
    _n_mctruths  = 0;
    _mc_nu_pdg   = bogus::PDG;
    _mc_nu_q2    = bogus::DOUBLE;
    _mc_nu_pos_x = bogus::POS;
    _mc_nu_pos_y = bogus::POS;
    _mc_nu_pos_z = bogus::POS;
    _mc_lepton_pdg = bogus::PDG;
    _mc_lepton_mom = bogus::DOUBLE;
    _mc_ccnc   = "";
    _mc_mode   = "";
    _true_nu_slice_ID           = -1;
    _true_nu_slice_completeness = bogus::DOUBLE;
    _true_nu_slice_purity       = bogus::DOUBLE;

    _n_slices                   = 0;
    //_n_primary_tracks  = 0;
    //_n_primary_showers = 0;

    /////////////////////////////
    // FlashMatch Slice Info
    /////////////////////////////
    _flash_match_nu_slice_ID    = -1;

    /////////////////////////////
    // Pandora Slice Info
    /////////////////////////////
    _pandora_nu_slice_ID        = -1;

    /////////////////////////////
    // PFParticle Variables
    /////////////////////////////
    // True stuff
    _pfp_purity.clear();
    _pfp_completeness.clear();
    _pfp_has_truth.clear();
    _pfp_trackID.clear();
    _pfp_true_pdg.clear();
    _pfp_true_energy.clear();
    _pfp_true_ke.clear();
    _pfp_true_px.clear();
    _pfp_true_py.clear();
    _pfp_true_pz.clear();
    _pfp_true_length.clear();
    _pfp_true_origin.clear();

    // Reco pfp stuff
    _pfp_pdg.clear();
    _pfp_trk_shr_score.clear();
    _pfp_x.clear();
    _pfp_y.clear();
    _pfp_z.clear();

    /////////////////////////////
    // Track Variables
    /////////////////////////////
    _trk_length.clear();
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

    /////////////////////////////
    // Shower Variables
    /////////////////////////////
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
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// fillNull is fills the output ttree with null(ish) values when accessing data
// products fails.
void hyperon::HyperonProduction::fillNull()
{
    clearTreeVariables();

    fTree->Fill();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        if (_mc_particle_map.find(sigmaZeroDaughter_ids.at(i_d)) == _mc_particle_map.end()) continue;

        auto p = _mc_particle_map.at(sigmaZeroDaughter_ids.at(i_d));

        if (p->PdgCode() == pdg::Lambda) {
            std::vector<int> _ids = getChildIds(p);
            lambdaDaughter_ids.insert(lambdaDaughter_ids.begin(), _ids.begin(), _ids.end());
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> hyperon::HyperonProduction::getChildIds(const art::Ptr<simb::MCParticle> &p, bool IsNeutron)
{
    std::vector<int> _decay_ids;

    if (p->EndProcess() != "Decay" && !IsNeutron && !pdg::isKaon(p)) return _decay_ids;

    for (int i_d = 0; i_d < p->NumberDaughters(); i_d++) {
        if (_mc_particle_map.find(p->Daughter(i_d)) == _mc_particle_map.end()) continue;

        art::Ptr<simb::MCParticle> daughter = _mc_particle_map.at(p->Daughter(i_d));

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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

int hyperon::HyperonProduction::getOrigin(int trackid)
{
    if (std::find(primary_ids.begin(), primary_ids.end(), trackid) != primary_ids.end()) return 1;
    else if (std::find(lambdaDaughter_ids.begin(), lambdaDaughter_ids.end(), trackid) != lambdaDaughter_ids.end()) return 2;
    else if (std::find(sigmaZeroDaughter_ids.begin(), sigmaZeroDaughter_ids.end(), trackid) != sigmaZeroDaughter_ids.end()) return 5;
    else return 3;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(hyperon::HyperonProduction)
