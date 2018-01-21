#pragma once

#include "Core/Containers/Array.h"
#include "Core/Containers/Map.h"

#include "glm/vec2.hpp"

//Uses concepts from 
// * https://infoscience.epfl.ch/record/100269/files/Kallmann_and_al_Geometric_Modeling_03.pdf -> For the overall implementation strategy
// * http://www.dtecta.com/files/GDC17_VanDenBergen_Gino_Brep_Triangle_Meshes.pdf -> For the low level Face/HalfEdge data structure
// * http://2.3jachtuches.pagesperso-orange.fr/dossiers/triangul/doc/fast.pdf -> For Mesh::Locate()
//Minimal delaunay implementation would implement Insert/DeleteVertex functions
//Minimal constrained implementation would implement Insert/DeleteConstraintSegment functions
//TODO: Write set of iterators for Mesh: Vertices/Edges/Faces
//TODO: Replace index type returned from modifier functions with a custom handle type with the 2 MSBs reserved for type.
//TODO: Replace ObjectPool with Array and implement function to use EraseSwap and patch half-edge indices to point to relocated face.
//TODO: Replace ObjectPool with Array and implement function using EraseSwap + patch destinationVertex indices to point to relocated vertex.

namespace Delaunay {
	
    class DebugDraw;
    
	class Mesh {
	public:
		
		struct HalfEdge {
            typedef uint32_t Index;
			Index destinationVertex; //End Vertex Index
			Index oppositeHalfEdge; //Opposite HalfEdge
			static const Index InvalidIndex = -1;
			struct IncomingHalfEdgeIterator {
				Index operator*();
				void operator++();
				bool operator!=(const IncomingHalfEdgeIterator & rhs);
				IncomingHalfEdgeIterator & begin();
				IncomingHalfEdgeIterator & end();
				IncomingHalfEdgeIterator(Mesh & mesh, Index first);
				const Mesh & mesh;
			private:
				Index current;
				const Index first;
			};
			struct OutgoingHalfEdgeIterator {
				Index operator*();
				void operator++();
				bool operator!=(const OutgoingHalfEdgeIterator & rhs);
				OutgoingHalfEdgeIterator & begin();
				OutgoingHalfEdgeIterator & end();
				OutgoingHalfEdgeIterator(Mesh & mesh, Index first);
				const Mesh & mesh;
			private:
				Index current;
				const Index first;
			};
		};
		struct Face {
			//The lowest 4 bits of flags are reserved for use by the mesh implementation
            HalfEdge::Index flags;
            HalfEdge::Index matId;
			HalfEdge edges[3];
			//Returns next (CCW) half-edge
            static inline HalfEdge::Index nextHalfEdge(HalfEdge::Index h);
			//Returns prev (CW) half-edge
            static inline HalfEdge::Index prevHalfEdge(HalfEdge::Index h);
		};
		static_assert(sizeof(Face) == sizeof(HalfEdge) * 4, "Face struct must be 4x the size of the HalfEdge");
		struct Vertex {
			glm::dvec2 position;
            HalfEdge::Index edge; //Can be either incoming or outgoing edge
			size_t constraintCount;
			size_t endPointCount;
			
			HalfEdge::IncomingHalfEdgeIterator IncomingEdges(Mesh & mesh);
			HalfEdge::OutgoingHalfEdgeIterator OutgoingEdges(Mesh & mesh);
		};
        struct ConstraintSegment {
            Oryol::Array<size_t> edgePairs;
            HalfEdge::Index startVertex;
            HalfEdge::Index endVertex;
        };
		struct LocateResult {
			size_t object;
			enum Code { None, Vertex, Edge, Face } type;
			bool operator!() const {
				return type == None;
			}
		};

		//Initialises the Delaunay Triangulation with a square mesh with specified width and height
		//Creates 5 vertices, and 6 faces. Vertex with index 0 is an infinite vertex
		void Setup(double width, double height);
		//Inserts a vertex by splitting an existing face/edge or returning an existing vertex 
		//if one exists at the specified point
		size_t InsertVertex(double x, double y);
		
		//Find which primitive the specified point is inside
		//Will only return primitives which are deemed to be "real"
		LocateResult Locate(double x, double y);
        
        void SetDebugDraw(DebugDraw * debug);
        void DrawDebugData();

	private:
		struct Impl;
		friend struct Impl;
		
		
        inline HalfEdge & edgeAt(HalfEdge::Index index);
        inline const HalfEdge & edgeAt(HalfEdge::Index index) const;
		
        inline bool isHalfEdgeConstrained(HalfEdge::Index h) {
			Face & f = this->faces[h / 4];
			return f.flags & (1 << (h & 3));
		}
        inline bool isHalfEdgeReal(HalfEdge::Index h) {
			Face & f = this->faces[h / 4];
			return (f.flags & 0x1) == 0;
		}


        Oryol::Array<Face> faces;
        Oryol::Array<Vertex> vertices;
        Oryol::Array<ConstraintSegment> constraints;
        
        Oryol::Map<HalfEdge::Index,size_t> halfEdgeToEdgePairMapping;
        Oryol::Array<Oryol::Array<size_t>> edgePairToConstraintMapping;
        
        DebugDraw * debugDraw;

	};
    
    class DebugDraw {
    public:
        virtual void DrawVertex(glm::vec2 position) = 0;
        virtual void DrawEdge(glm::vec2 origin, glm::vec2 destination, bool constrained) = 0;
    };


}
