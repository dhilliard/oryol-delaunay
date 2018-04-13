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

#include "Core/Containers/Array.h"
#include "Gfx/Gfx.h"
#include "glm/mat4x4.hpp"

struct Color {
	float r, g, b, a;
	Color(float r, float g, float b, float a = 1.0f)
		: r(r), g(g), b(b), a(a) {}
};
class DebugBatch {
public:
	void Setup(const Oryol::GfxSetup & setup);
	void Discard();

	void Triangle(float x1, float y1, float x2, float y2, float x3, float y3, const Color & color);
	void Line(float x1, float y1, float x2, float y2, const Color & color);
	void Point(float x, float y, float size, const Color & color);

	void Draw(glm::mat4x4 projectionMatrix);
	const int MaxNumTriangleVertices = 3 * 1024;
	const int MaxNumLineVertices = 2 * 4 * 1024;
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
