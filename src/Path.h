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
        static bool FindPath(const glm::dvec2 & start, const glm::dvec2 & end, double radius, Oryol::Array<uint32_t> & path);
    }
}
