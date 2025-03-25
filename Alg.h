#pragma once

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "ubana/searchingfornues/Selection/CommonDefs/LLR_PID.h"

#include <vector>

namespace hyperon {
    namespace alg {
        // conversion from 3D detector position to wire-time coordinate
        // taken from https://github.com/cthorpe123/HyperonProduction/blob/master/Alg/Position_To_Wire.h
        namespace geometry {
            constexpr double A_w = 3.33328;
            constexpr double C_U = 338.140;
            constexpr double C_V = 2732.53;
            constexpr double C_Y = 4799.19;
            constexpr double A_t = 18.2148;
            constexpr double C_t = 818.351;

            constexpr double cos60 = 0.5;
            constexpr double sin60 = sqrt(3)/2.0;

            int U_wire(const TVector3 &pos) { return A_w*(-sin60*pos.Y()+cos60*pos.Z())+C_U; }
            int V_wire(const TVector3 &pos) { return A_w*(sin60*pos.Y()+cos60*pos.Z())+C_V; }
            int Y_wire(const TVector3 &pos) { return A_w*pos.Z() + C_Y; }
            int tick(const TVector3 &pos) { return A_t*pos.X() + C_t; }

            double dUdt(const TVector3 &dir) { return A_w/A_t*(-sin60*dir.Y()/dir.X()+cos60*dir.Z()/dir.X()); }
            double dVdt(const TVector3 &dir) { return A_w/A_t*(sin60*dir.Y()/dir.X()+cos60*dir.Z()/dir.X()); }
            double dYdt(const TVector3 &dir) { return A_w/A_t*dir.Z()/dir.X(); }

            double AngleU(const TVector3 &dir){
                bool invert = dir.X() < 0;
                double angle = (180/3.141)*atan(dUdt(dir));
                angle -= 2*(angle-45.0);
                if(invert && angle < 0) angle += 180;
                if(invert && angle > 0) angle -= 180;
                return angle;
            }
            double AngleV(const TVector3 &dir){
                bool invert = dir.X() < 0;
                double angle = (180/3.141)*atan(dVdt(dir));
                angle -= 2*(angle-45.0);
                if(invert && angle < 0) angle += 180;
                if(invert && angle > 0) angle -= 180;
                return angle;
            }
            double AngleY(const TVector3 &dir){
                bool invert = dir.X() < 0;
                double angle = (180/3.141)*atan(dYdt(dir));
                angle -= 2*(angle-45.0);
                if(invert && angle < 0) angle += 180;
                if(invert && angle > 0) angle -= 180;
                return angle;
            }
        }

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

                float llr_pid = pid_calc.LLR_many_hits_one_plane(dedx_values,par_values,plane);
                this_llr_pid += llr_pid;

                float calo_energy_partial = 0;
                for (size_t i = 0; i < dedx_values_partial.size(); i++)
                    calo_energy_partial += dedx_values_partial[i] * pitch_partial[i];
            }

            this_llr_pid_score = atan(this_llr_pid/100.) * 2./3.14159266;

            return this_llr_pid_score;
        }

        bool ADCThreshold(const art::Ptr<recob::Hit> hit, float min_adc_value) {
            return hit->PeakAmplitude() >= min_adc_value;
        }

        void BuildCTWindow(
                const std::vector<art::Ptr<recob::Hit>> hits,
                std::vector<std::vector<float>> &window_plane0,
                std::vector<std::vector<float>> &window_plane1,
                std::vector<std::vector<float>> &window_plane2,
                const TVector3 &nu_vtx,
                const unsigned int window_wires,
                const unsigned int window_ticks, // XXX: ignored for now
                float min_adc_value) {
            const unsigned int nu_wire_u = geometry::U_wire(nu_vtx);
            const unsigned int nu_wire_v = geometry::V_wire(nu_vtx);
            const unsigned int nu_wire_y = geometry::Y_wire(nu_vtx);
            /* const unsigned int nu_time   = geometry::tick(nu_vtx); */

            // Computes the start wire number given the number of wires we're
            // looking at and the located reconstructed neutrino vertex.
            const unsigned int start_wire_u = nu_wire_u - ((int)window_wires / 2);
            const unsigned int start_wire_v = nu_wire_v - ((int)window_wires / 2);
            const unsigned int start_wire_y = nu_wire_y - ((int)window_wires / 2);

            for (size_t i_hit = 0; i_hit < hits.size(); ++i_hit) {
                const art::Ptr<recob::Hit> hit = hits.at(i_hit);

                if (!ADCThreshold(hit, min_adc_value))
                    continue;

                const int wireID = hit->WireID().Wire;
                const geo::View_t view = hit->View(); // TODO: iterate through enum constants, don't cast to int
                const float start_tick = hit->PeakTime();

                const unsigned int wire_index_u = wireID - start_wire_u;
                const unsigned int wire_index_v = wireID - start_wire_v;
                const unsigned int wire_index_y = wireID - start_wire_y;

                switch (view) {
                    case geo::kU:
                        if ((wire_index_u) < window_wires)
                            window_plane0.at(wire_index_u).push_back(start_tick);
                        break;
                    case geo::kV:
                        if ((wire_index_v) < window_wires)
                            window_plane1.at(wire_index_v).push_back(start_tick);
                        break;
                    case geo::kW:
                        if ((wire_index_y) < window_wires)
                            window_plane2.at(wire_index_y).push_back(start_tick);
                        break;
                    default:
                        break;
                }
            }

            window_plane0.shrink_to_fit();
            window_plane1.shrink_to_fit();
            window_plane2.shrink_to_fit();

            return;
        }
    }
}
