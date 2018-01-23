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

	//This function happens to compute the determinant of the matrix to solve the intersection point
	//of the vectors AB and AC
	//Otherwise;
	//		Returns 0; Point C is on the line AB
	//		Returns -ve : Point C is right of the line AB
	//		Returns +ve : Point C is left of the line AB
	inline double Sign(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c) {
		return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
	}

	inline double DistanceSquared(glm::dvec2 v) {
		return (v.x)*(v.x) + (v.y)*(v.y);
	}
	// Uses code sourced from: http://www.randygaul.net/2014/07/23/distance-point-to-line-segment/
	double DistanceSquaredPointToLineSegment(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p);

	glm::dvec2 OrthogonallyProjectPointOnLineSegment(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p);
	glm::dvec2 ComputeCircumcenter(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c);
    
    inline bool CounterClockwise(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c){
        return (c.y-a.y) * (b.x-a.x) > (b.y-a.y) * (c.x-a.x);
    }
    inline bool SegmentsIntersect(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c, const glm::dvec2 & d){
        return CounterClockwise(a, c, d) != CounterClockwise(b,c,d) and CounterClockwise(a,b,c) != CounterClockwise(a, b, d);
    }
    inline bool ComputeIntersection(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c, const glm::dvec2 & d, glm::dvec2 * intersection = nullptr){
        //Compute determinant for the matrix;
        double divisor = (a.x - b.x)*(c.y - d.y) + (b.y - a.y)*(c.x - d.x);

        if(divisor != 0.0){
            double t1 = (a.x*(c.y - d.y) + a.y*(d.x - c.x) + c.x*d.y - c.y*d.x) / divisor;
            double t2 = (a.x*(c.y - b.y) + a.y*(b.x - c.x) - b.x*c.y + b.y*c.x) / divisor;
            if(!(0 <= t1 && t1 <=1 && 0 <= t2 && t2 <=1)){
                if(intersection != nullptr){
                    *intersection = a + t1*(b - a);
                }
                return true;
            }
        }
        return false;
    }
}
