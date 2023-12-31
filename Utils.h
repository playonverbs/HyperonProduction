#pragma once

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"

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
	}
}
