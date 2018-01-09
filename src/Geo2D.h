#pragma once

namespace Geo2D {
	struct AABB {
		double min_x, min_y;
		double max_x, max_y;
		double Width() const { return max_x - min_x; }
		double Height() const { return max_y - min_y; }
	};

	struct ClipResult {
		double x1, y1, x2, y2;
		bool success;
	};

	ClipResult ClipSegment(double x1, double y1, double x2, double y2, const AABB & bb);

}
