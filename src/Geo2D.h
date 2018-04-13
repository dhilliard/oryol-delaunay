/*
	Copyright 2013-2018 Denis Hilliard <denis.z.hilliard@gmail.com>

	Permission is hereby granted, free of charge, to any person obtaining a copy of this 
	software and associated documentation files(the "Software"), to deal in the Software 
	without restriction, including without limitation the rights to use, copy, modify, merge, 
	publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
	to whom the Software is furnished to do so, subject to the following conditions :

	The above copyright notice and this permission notice shall be included in all copies or 
	substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once
#include "glm/vec2.hpp"
#include "glm/geometric.hpp"
namespace Geo2D {
    
    template <typename T> inline bool IsInRange(T min, T value, T max) {
        return value >= min && value <= max;
    }
    
	struct AABB {
		glm::dvec2 min;
		glm::dvec2 max;
		double Width() const { return max.x - min.x; }
		double Height() const { return max.y - min.y; }
        inline bool IsPointInside(const glm::dvec2 & p) const {
            return IsInRange(min.x,p.x,max.x) && IsInRange(min.y,p.y,max.y);
        }
	};

	struct ClipResult {
		glm::dvec2 a, b;
		bool success;
	};

	ClipResult ClipSegment(const glm::dvec2 & a, const glm::dvec2 & b, const AABB & bb);
	
	//Computes the cross product of the vectors formed between AB and AC; 
	//	Returns 0; Point C is on the line AB
	//	Returns -ve : Point C is right of the line AB
	//	Returns +ve : Point C is left of the line AB
	//	Otherwise magnitude of cross product divided by 2 can be used as the area formed by the triangle
	inline double Sign(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c) {
		return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
	}

	inline double DistanceSquared(glm::dvec2 v) {
		return (v.x)*(v.x) + (v.y)*(v.y);
	}
	// Uses code sourced from: http://www.randygaul.net/2014/07/23/distance-point-to-line-segment/
    double DistanceSquaredPointToLine(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p);
	double DistanceSquaredPointToLineSegment(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p);

    glm::dvec2 OrthogonallyProjectPointOnLine(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p);
	glm::dvec2 OrthogonallyProjectPointOnLineSegment(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & p);
    
	glm::dvec2 ComputeCircumcenter(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c);
    
    inline bool CounterClockwise(const glm::dvec2 & p1, const glm::dvec2 & p2, const glm::dvec2 & p3){
        return ((p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y)) < 0;
    }
    inline bool SegmentsIntersect(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c, const glm::dvec2 & d){
        return CounterClockwise(a, c, d) != CounterClockwise(b,c,d) && CounterClockwise(a,b,c) != CounterClockwise(a, b, d);
    }
    bool ComputeIntersection(const glm::dvec2 & a, const glm::dvec2 & b, const glm::dvec2 & c, const glm::dvec2 & d, glm::dvec2 * intersection = nullptr);
}
