#pragma once

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "ubana/searchingfornues/Selection/CommonDefs/LLR_PID.h"

#include <vector>

namespace hyperon {
    namespace alg {
        struct Dedx {
            double weight_plane0 = 0.0f;
            double weight_plane1 = 0.0f;
            double weight_plane2 = 0.0f;

            double plane0 = -1.0f;
            double plane1 = -1.0f;
            double plane2 = -1.0f;

            double three_plane_average = -1.0f;
        };

        constexpr double p0_wireangle = 30 * 6.28 / 360.0;
        constexpr double p1_wireangle = -30 * 6.28 / 360.0;
        constexpr double p2_wireangle = 90 * 6.28 / 360.0;

        // TODO: abstract to class to allow for parameterisation. 
        constexpr double threshold = 0.175;

        double GetMeandEdX(art::Ptr<anab::Calorimetry> calo) {
            double totalE = 0.0f,
                   totalX = 0.0f;

            // fail if we only have 1 calorimetry point
            if (calo->XYZ().size() < 2) return -1.0f;

            for (size_t i_p = 0; i_p < calo->XYZ().size() - 1; i_p++) {
                anab::Point_t thisPos = calo->XYZ().at(i_p);
                anab::Point_t nextPos = calo->XYZ().at(i_p + 1);

                TVector3 step(thisPos.X() - nextPos.X(),
                              thisPos.X() - nextPos.X(),
                              thisPos.X() - nextPos.X());

                totalX += step.Mag();
                totalE += calo->dEdx().at(i_p) * step.Mag();
            }

            return totalE / totalX;
        }

        double dEdXPlaneWeight(art::Ptr<recob::Track> track, int i_plane) {
            double dy = track->End().y() - track->Start().y();
            double dz = track->End().z() - track->Start().z();

            TVector3 trackvec(0, dy, dz);
            trackvec = trackvec.Unit();
            TVector3 zaxis(0, 0, 1);
            double costheta_yz = trackvec.Dot(zaxis);
            double theta_yz = TMath::ACos(costheta_yz);
            if ((dy < 0) && (theta_yz > 0)) theta_yz *= -1;

            auto getThetaWires = [=](double wireangle) {
                return std::min(std::abs(wireangle - theta_yz), std::abs((-1*(6.28 - wireangle) - theta_yz)));
            };
            double theta_to_wires = 0.0f;

            switch (i_plane) {
                case 0u:
                    theta_to_wires = getThetaWires(p0_wireangle);
                    break;
                case 1u:
                    theta_to_wires = getThetaWires(p1_wireangle);
                    break;
                case 2u:
                    theta_to_wires = getThetaWires(p2_wireangle);
                    break;
                default: // uh oh
                    break;
            }

            double angle_planeweight = sin(theta_to_wires) * sin(theta_to_wires);
            if (angle_planeweight < threshold) angle_planeweight = 0;
            if (angle_planeweight != 0) angle_planeweight = 1;

            return angle_planeweight;
        }

        Dedx ThreePlaneMeandEdX(art::Ptr<recob::Track> track, const std::vector<art::Ptr<anab::Calorimetry>> calos) {
            Dedx data;

            double totaldEdX = 0.0f,
                   totalWeight = 0.0f;

            for (const art::Ptr<anab::Calorimetry> &calo : calos) {
                const int plane = calo->PlaneID().Plane;

                if (plane != 0 && plane != 1 && plane != 2) continue;

                double dEdx = GetMeandEdX(calo);

                if (dEdx < 0) continue;

                double plane_weight = dEdXPlaneWeight(track, plane);

                switch (plane) {
                    case 0u:
                        data.weight_plane0 = plane_weight;
                        data.plane0 = dEdx;
                        break;
                    case 1u:
                        data.weight_plane1 = plane_weight;
                        data.plane1 = dEdx;
                        break;
                    case 2u:
                        data.weight_plane2 = plane_weight;
                        data.plane2 = dEdx;
                        break;
                }

                totaldEdX += dEdx * plane_weight;
                totalWeight += plane_weight;
            }

            if (totalWeight > 0.0f) data.three_plane_average = totaldEdX / totalWeight;

            return data;
        }

        double LLRPID(
                const std::vector<art::Ptr<anab::Calorimetry>> calos,
                searchingfornues::LLRPID &pid_calc, // TODO: consider changing to LLRPID*
                double ResRangeCutoff = 5.0) {
            double this_llr_pid = 0;
            double this_llr_pid_score = 0;

            for (auto const &calo : calos) {
                auto const &plane = calo->PlaneID().Plane;
                auto const &dedx_values = calo->dEdx();
                auto const &rr = calo->ResidualRange();
                auto const &pitch = calo->TrkPitchVec();
                std::vector<std::vector<float>> par_values;
                par_values.push_back(rr);
                par_values.push_back(pitch);

                // Get partial length PIDs
                std::vector<std::vector<float>> par_values_partial;
                std::vector<float> dedx_values_partial, rr_partial, pitch_partial;

                if (dedx_values.size() != rr.size() || rr.size() != pitch.size())
                    throw cet::exception("HyperonProduction") << "Track calo point size mismatch\n";

                for (size_t i_p = 0; i_p < calo->dEdx().size(); i_p++) {
                    if (rr.at(i_p) > ResRangeCutoff) continue; // TODO: impl ResRangeCutoff

                    dedx_values_partial.push_back(calo->dEdx().at(i_p));
                    rr_partial.push_back(calo->ResidualRange().at(i_p));
                    pitch_partial.push_back(calo->TrkPitchVec().at(i_p));
                }
                par_values_partial.push_back(rr_partial);
                par_values_partial.push_back(pitch_partial);

                if (calo->ResidualRange().size() == 0) continue;

                float calo_energy = 0;
                for (size_t i = 0; i < dedx_values.size(); i++)
                    calo_energy += dedx_values[i] * pitch[i];

                // TODO: do searchingfornues::llr_pid_calculator setup
                float llr_pid = pid_calc.LLR_many_hits_one_plane(dedx_values,par_values,plane);
                this_llr_pid += llr_pid;

                float calo_energy_partial = 0;
                for (size_t i = 0; i < dedx_values_partial.size(); i++)
                    calo_energy_partial += dedx_values_partial[i] * pitch_partial[i];
            }

            this_llr_pid_score = atan(this_llr_pid/100.) * 2./3.14159266;

            return this_llr_pid_score;
        }

        // Adapted from searchingfornues
        void GetMoliereRadius(art::Event const& evt, const art::Ptr<recob::PFParticle> pfp, double& medangle, double& rmsangle, const std::string &productLabel, const std::string &showerLabel, const std::string &spcpntsLabel) {
            const art::Ptr<recob::Shower> shower =
                util::GetAssocProduct<recob::Shower>(pfp, evt, productLabel, showerLabel);

            const std::vector<art::Ptr<recob::SpacePoint>> spcpnts =
                util::GetAssocProductVector<recob::SpacePoint>(pfp, evt, productLabel, spcpntsLabel);

            if (shower.isNull() || spcpnts.size() == 0) {
                medangle = -999.0;
                rmsangle = -999.0;
                return;
            }

            TVector3 shrdir(shower->Direction().X(),
                            shower->Direction().Y(),
                            shower->Direction().Z());

            TVector3 shrvtx(shower->ShowerStart().X(),
                            shower->ShowerStart().Y(),
                            shower->ShowerStart().Z());

            std::vector<float> angle_v;

            medangle = 0;
            rmsangle = 0;

            for (auto &sp: spcpnts) {
                auto spxyz = sp->XYZ();

                TVector3 sppos(spxyz[0], spxyz[1], spxyz[2]);

                TVector3 sptovtx = sppos - shrvtx;

                if (sptovtx.Mag() == 0) continue;

                float angle = acos( sptovtx.Dot( shrdir ) / (sptovtx.Mag() * shrdir.Mag() ));
                angle *= (180./3.14);

                angle_v.push_back( angle );
            }

            if (angle_v.size() == 0) {
                medangle = std::numeric_limits<double>::max();
                rmsangle = std::numeric_limits<double>::max();
                return;
            }

            for (size_t d = 0; d < angle_v.size(); d++) {
                medangle += angle_v[d];
            }

            medangle /= angle_v.size();

            for (size_t d = 0; d < angle_v.size(); d++) {
                rmsangle += (angle_v[d] - medangle) * (angle_v[d] - medangle);
            }

            rmsangle /= sqrt( angle_v.size() );

            return;
        }
    }
}
