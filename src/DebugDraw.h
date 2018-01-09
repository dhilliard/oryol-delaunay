#pragma once

#include "Core/Containers/Array.h"
#include "Gfx/Gfx.h"
#include "glm/mat4x4.hpp"

struct Color {
	float r, g, b, a;
	Color(float r, float g, float b, float a = 1.0f)
		: r(r), g(g), b(b), a(a) {}
};
class DebugDraw {
public:
	void Setup(const Oryol::GfxSetup & setup);
	void Discard();
	void Triangle(float x1, float y1, float x2, float y2, float x3, float y3, const Color & color);
	void Line(float x1, float y1, float x2, float y2, const Color & color);
	void Point(float x, float y, float size, const Color & color);
	void Draw(glm::mat4x4 projectionMatrix);
	const int MaxNumTriangleVertices = 3 * 1024;
	const int MaxNumLineVertices = 2 * 1024;
	const int MaxNumPointVertices = 1 * 1024;
private:
	struct vertex_t {
		float x, y;
		uint32_t color;
	};
	struct point_t {
		float x, y, size;
		uint32_t color;
	};
	Oryol::Array<vertex_t> triangles;
	Oryol::Array<vertex_t> lines;
	Oryol::Array<point_t> points;
	Oryol::DrawState triangleDrawState;
	Oryol::DrawState lineDrawState;
	Oryol::DrawState pointDrawState;
	Oryol::ResourceLabel resourceLabel;
};