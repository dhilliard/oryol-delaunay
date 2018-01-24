#include "Geo2D.h"
#include <limits>
#include <algorithm>
using namespace glm;

Geo2D::ClipResult Geo2D::ClipSegment(const glm::dvec2 & a, const glm::dvec2 & b, const AABB & bb)
{
	ClipResult result;

	result.a = a;
	result.b = b;
	
	result.success = false;

	if ((a.x > bb.max.x && b.x > bb.max.x)
		|| (a.x < bb.min.x && b.x < bb.min.x)
		|| (a.y > bb.max.y && b.y > bb.max.y)
		|| (a.y < bb.min.y && b.y < bb.min.y))
	{
		result.success = false;
	}
	else
	{
		glm::dvec2 n = b - a;

		double tmin = std::numeric_limits<double>::lowest();
		double tmax = std::numeric_limits<double>::infinity();

		if (n.x != 0.0)
		{
			double tx1 = (bb.min.x - a.x) / n.x;
			double tx2 = (bb.max.x - a.x) / n.x;

			tmin = std::max(tmin, std::min(tx1, tx2));
			tmax = std::min(tmax, std::max(tx1, tx2));
		}

		if (n.y != 0.0)
		{
			double ty1 = (bb.min.y - a.y) / n.y;
			double ty2 = (bb.max.y - a.y) / n.y;

			tmin = std::max(tmin, std::min(ty1, ty2));
			tmax = std::min(tmax, std::max(ty1, ty2));
		}

		if (tmax >= tmin)
		{

			if (tmax < 1)
			{
				//Clip end point
				result.b = n * tmax + a;
			}

			if (tmin > 0)
			{
				//Clip start point
				result.a = n * tmin + a;
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


// Uses code sourced from: http://www.randygaul.net/2014/07/23/distance-point-to-line-segment/
double Geo2D::DistanceSquaredPointToLineSegment(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p) {
	glm::dvec2 n = b - a;
	glm::dvec2 pa = a - p;
	glm::dvec2 c = n * (glm::dot(pa, n) / glm::dot(n, n));
	glm::dvec2 d = pa - c;
	return glm::dot(d, d);
}

glm::dvec2 Geo2D::OrthogonallyProjectPointOnLineSegment(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p) {
	auto ap = p - a;
	auto ab = b - a;
	return a + glm::dot(ap, ab) / dot(ab, ab) * ab;
}

glm::dvec2 Geo2D::ComputeCircumcenter(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c) {
	
	dvec2 diffAB = b - a;
	dvec2 diffCA = c - a;
	double lengthSquaredAB = DistanceSquared(diffAB);
	double lengthSquaredCA = DistanceSquared(diffCA);
	double denominator = 0.5 / Sign(b, c, a);
	return a + denominator * dvec2(diffCA.y * lengthSquaredAB - diffAB.y * lengthSquaredCA, diffAB.x * lengthSquaredCA - diffCA.x * lengthSquaredAB);
}

