//
//  Path.h
//  oryol-delaunay
//
//  Created by Denis Hilliard on 9/2/18.
//
//

#pragma once

#include "glm/vec2.hpp"
#include "Core/Containers/Array.h"

namespace Delaunay {
    class Mesh;
    namespace Path {
        bool FindPath(Mesh & mesh, const glm::dvec2 & start, const glm::dvec2 & end, const double radius, Oryol::Array<uint32_t> & pathFaces, Oryol::Array<uint32_t> & pathEdges);
    }
}
