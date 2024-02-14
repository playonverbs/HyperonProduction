#pragma once

#include "TVector3.h"

#include <utility>
#include <vector>

namespace hyperon {
	namespace fv {
		const std::vector<double> TPCCenter = {
			126.625, 0.97, 518.5 // center of active TPC.
		};

		const std::vector<double> TPCSideLengths = {
			236.35, 233.0, 1036.8
		};

		using sides = std::pair<double, double>;

		const sides FVx = {0.0,     256.35};
		const sides FVy = {-115.53, 117.47};
		const sides FVz = {0.1,     1036.9};

		inline bool inActiveTPC(TVector3 pos) {
			if (pos.X() > FVx.second || pos.X() < FVx.first) return false;
			if (pos.Y() > FVy.second || pos.Y() < FVy.first) return false;
			if (pos.Z() > FVz.second || pos.Z() < FVz.first) return false;

			return true;
		}

		inline bool inActiveTPC(double x, double y, double z) {
			if (x > FVx.second || x < FVx.first) return false;
			if (y > FVy.second || y < FVy.first) return false;
			if (z > FVz.second || z < FVz.first) return false;

			return true;
		}
	}
}
