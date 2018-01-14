#pragma once
#include "glm/vec2.hpp"
#include "glm/geometric.hpp"
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

	//Computes the side of the line (a - b) and point p
	inline double Sign(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c) {
		return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
	}
	inline int DetermineSide(glm::dvec2 a, glm::dvec2 b, glm::dvec2 c) {
		double dot = Sign(a,b,c);
		return (dot == 0) ? 0 : (dot > 0 ? 1 : -1);
	}
	inline double DistanceSquared(glm::dvec2 v) {
		return (v.x)*(v.x) + (v.y)*(v.y);
	}
	// Uses code sourced from: http://www.randygaul.net/2014/07/23/distance-point-to-line-segment/
	inline double DistanceSquaredPointToLineSegment(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p) {
		auto n = b - a;
		auto pa = a - p;
		auto c = n * (glm::dot(pa, n) / glm::dot(n, n));
		auto d = pa - c;
		return glm::dot(d, d);
	}
	inline glm::dvec2 OrthogonallyProjectPointOnLineSegment(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p) {
		auto ap = p - a;
		auto ab = b - a;
		return a + glm::dot(ap,ab) / dot(ab,ab) * ab;
	}
}
