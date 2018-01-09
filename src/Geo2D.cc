#include "Geo2D.h"
#include <limits>
#include <algorithm>

Geo2D::ClipResult Geo2D::ClipSegment(double x1, double y1, double x2, double y2, const AABB & bb)
{
	ClipResult result {x1,y1, x2,y2, false};
	if ((x1 > bb.max_x && x2 > bb.max_x)
		|| (x1 < bb.min_x && x2 < bb.min_x)
		|| (y1 > bb.max_y && y2 > bb.max_y)
		|| (y1 < bb.min_y && y2 < bb.min_y))
	{
		result.success = false;
	}
	else
	{
		double nx = x2 - x1;
		double ny = y2 - y1;

		double tmin = std::numeric_limits<double>::lowest();
		double tmax = std::numeric_limits<double>::infinity();

		if (nx != 0.0)
		{
			double tx1 = (bb.min_x - x1) / nx;
			double tx2 = (bb.max_x - x1) / nx;

			tmin = std::max(tmin, std::min(tx1, tx2));
			tmax = std::min(tmax, std::max(tx1, tx2));
		}

		if (ny != 0.0)
		{
			double ty1 = (bb.min_y - y1) / ny;
			double ty2 = (bb.max_y - y1) / ny;

			tmin = std::max(tmin, std::min(ty1, ty2));
			tmax = std::min(tmax, std::max(ty1, ty2));
		}

		if (tmax >= tmin)
		{

			if (tmax < 1)
			{
				//Clip end point
				result.x2 = nx * tmax + x1;
				result.y2 = ny * tmax + y1;
			}

			if (tmin > 0)
			{
				//Clip start point
				result.x1 = nx * tmin + x1;
				result.y1 = ny * tmin + y1;
			}
			result.success = true;
		}
		else {
			//Something has gone very wrong if we hit this point
			result.success = false;
		}
	}
	return result;
}

