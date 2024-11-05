#pragma once

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include <string>
#include <vector>

namespace hyperon {
    namespace util {

        // TODO: Add prefix when not using scoped enums
        enum struct GenEventType {
            QEL = 0,
            RES = 1,
            DIS = 2,
            COH = 3,
            ElectronScatter = 5,
            MEC = 10,
            Diffractive = 11,
            HYP = 1095,
            Other
        };

        // Convert the CCNC code from a neutrino generator to a string
        std::string GetCCNC(int ccnc) {
            switch (ccnc) {
                case 0:  return "CC";
                case 1:  return "NC";
                default: return "None";
            }

            return "None";
        }

        // Convert the Event Mode code from a neutrino generator to a string
        // ranges 0 - 11 account for the GENIE event generator.
        // 1095 accounts for the NuWro generator hyperon code.
        std::string GetEventType(int mode) {
            /* GenEventType type = static_cast<GenEventType>(mode); */

            switch (mode) {
                case 0:    return "QEL";
                case 1:    return "RES";
                case 2:    return "DIS";
                case 3:    return "COH";
                case 5:    return "ElectronScatter";
                case 10:   return "MEC";
                case 11:   return "Diffractive";
                case 1095: return "HYP";

                default:   return "Other";
            }

            return "Other";
        }

        template <typename T>
        std::vector<art::Ptr<T>> GetProductVector(const art::Event &e, const std::string &label)
        {
            auto productHandle = e.getValidHandle<std::vector<T>>(label);
            bool success = productHandle.isValid();

            if (!success)
            {
                return std::vector<art::Ptr<T>>();
            }

            std::vector<art::Ptr<T>> productVector;
            art::fill_ptr_vector(productVector, productHandle);

            return productVector;
        }

        template <typename T, typename U>
        std::vector<art::Ptr<T>> GetAssocProductVector(const art::Ptr<U> &pProd, const art::Event &e, const std::string &label, const std::string &assocLabel)
        {
            auto productHandle = e.getValidHandle<std::vector<U>>(label);
            bool success = productHandle.isValid();

            if (!success)
            {
                return std::vector<art::Ptr<T>>();
            }

            const art::FindManyP<T> findParticleAssocs(productHandle, e, assocLabel);

            return findParticleAssocs.at(pProd.key());
        }

        template <typename T, typename U>
        art::Ptr<T> GetAssocProduct(const art::Ptr<U> &pProd, const art::Event &e, const std::string &label, const std::string &assocLabel)
        {
            std::vector<art::Ptr<T>> assocProducts = GetAssocProductVector<T>(pProd, e, label, assocLabel);
            if (assocProducts.empty())
            {
                return art::Ptr<T>();
            }

            return assocProducts.at(0);
        }

        // Adapted from PeLEE code (searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h)
        art::Ptr<simb::MCParticle> getAssocMCParticle(art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &hittruth,
                                                        const std::vector<art::Ptr<recob::Hit>> &hits,
                                                        float &purity)
        {
            float pfpcharge = 0; // total hit charge from clusters
            float maxcharge = 0; // charge backtracked to best match

            std::unordered_map<int, float> trkide;
            std::unordered_map<int, float> trkq;
            double maxe = -1, tote = 0;
            art::Ptr<simb::MCParticle> maxp_me; // pointer for the matched particle;

            for (auto h : hits)
            {
                pfpcharge += h->Integral();

                std::vector<art::Ptr<simb::MCParticle>> particles = hittruth.at(h.key());
                std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth.data(h.key());

                for (size_t i_p = 0; i_p < particles.size(); ++i_p)
                {
                    trkide[particles[i_p]->TrackId()] += match_vec[i_p]->energy;
                    trkq[particles[i_p]->TrackId()] += h->Integral() * match_vec[i_p]->ideFraction;
                    tote += match_vec[i_p]->energy;

                    if (trkide[particles[i_p]->TrackId()] > maxe)
                    {
                        maxe = trkide[particles[i_p]->TrackId()];
                        maxp_me = particles[i_p];
                        maxcharge = trkq[particles[i_p]->TrackId()];
                    }
                }
            }

            purity = maxcharge / pfpcharge;
            /* completeness = 0; */

            // number of true hits in the reco'd particle vs true hit count.

            return maxp_me;
        }

        bool posMatch(TVector3 p1, TVector3 p2, const double _epsilon = 0.0001) {
            return (p1 - p2).Mag() < _epsilon;
        }

        enum struct OriginType {
            Neutrino = 1,
            Lambda = 2,
            Other = 3,
            SigmaZero = 5,
        };
    }

    namespace pdg {
        constexpr int Lambda      = 3122;
        constexpr int NeutralKaon = 311;
        constexpr int SigmaZero   = 3212;

        inline bool isHyperon(const int pdgCode) {
            const int absPDG = abs(pdgCode);
            return absPDG == 3122 || absPDG == 3212 || absPDG == 3222 || absPDG == 3112;
        }
        inline bool isHyperon(const art::Ptr<simb::MCParticle> p) { return isHyperon(p->PdgCode()); }

        inline bool isPion(const int pdgCode) {
            const int absPDG = abs(pdgCode);
            return absPDG == 111  || absPDG == 211;
        }
        inline bool isPion(const art::Ptr<simb::MCParticle> p) { return isPion(p->PdgCode()); }

        inline bool isNucleon(const int pdgCode) {
            const int absPDG = abs(pdgCode);
            return absPDG == 111  || absPDG == 211;
        }
        inline bool isNucleon(const art::Ptr<simb::MCParticle> p) { return isNucleon(p->PdgCode()); }

        inline bool isLepton(const int pdgCode) {
            const int absPDG = abs(pdgCode);
            return absPDG == 11 || absPDG == 13 || absPDG == 15;
        }
        inline bool isLepton(const art::Ptr<simb::MCParticle> p) { return isLepton(p->PdgCode()); }

        inline bool isNeutrino(const int pdgCode) {
            const int absPDG = abs(pdgCode);
            return absPDG == 12 || absPDG == 14 || absPDG == 16;
        }
        inline bool isNeutrino(const art::Ptr<simb::MCParticle> p) { return isNeutrino(p->PdgCode()); }

        inline bool isKaon(const int pdgCode) {
            const int absPDG = abs(pdgCode);
            return absPDG == 321 || absPDG == 311 || absPDG == 130 || absPDG == 310;
        }
        inline bool isKaon(const art::Ptr<simb::MCParticle> p) { return isKaon(p->PdgCode()); }

        inline bool isPhoton(const int pdgCode) {
            const int absPDG = abs(pdgCode);
            return absPDG == 22;
        }
        inline bool isPhoton(const art::Ptr<simb::MCParticle> p) { return isPhoton(p->PdgCode()); }
    }
}
