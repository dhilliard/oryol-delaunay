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
#include "Core/Containers/Array.h"

namespace Delaunay {
    class Mesh;
    namespace Path {
        bool FindPath(Mesh & mesh, const glm::dvec2 & start, const glm::dvec2 & end, const double radius, Oryol::Array<uint32_t> & pathFaces, Oryol::Array<uint32_t> & pathEdges);
        void RefinePath(Mesh & mesh, const glm::dvec2 & start, const glm::dvec2 & end, const double radius, const Oryol::Array<uint32_t> & pathFaces, const Oryol::Array<uint32_t> & pathEdges, Oryol::Array<glm::vec2> & refinedPath);
    }
}
