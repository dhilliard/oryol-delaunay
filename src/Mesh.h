#pragma once

#include "Core/Containers/Array.h"
#include "Core/Containers/Map.h"
#include "Geo2D.h"
#include "glm/vec2.hpp"
#include "ObjectPool.h"

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
        struct Impl;

		struct HalfEdge {
            typedef uint32_t Index;
            enum IndexBits {
                Size = sizeof(Index) * 8
            };
			Index destinationVertex; //End Vertex Index
			Index oppositeHalfEdge; //Opposite HalfEdge
            Index constrained : 1; //Cache constraint state on the half edge
            Index edgePair : IndexBits::Size - 1; //Use this to references EdgeInfo struct
			static const Index InvalidIndex = -1;
		};
		struct Face {
            HalfEdge::Index userData[3];
			HalfEdge edges[3];
            bool isReal() const {
                return  edges[0].destinationVertex != 0 &&
                        edges[1].destinationVertex != 0 &&
                        edges[2].destinationVertex != 0;
            }
			//Returns next (CCW) half-edge
			static inline HalfEdge::Index nextHalfEdge(HalfEdge::Index h) {
				++h;
				return (h & 3) != 0 ? h : h - 3;
			}
			//Returns prev (CW) half-edge
			static inline HalfEdge::Index prevHalfEdge(HalfEdge::Index h) {
				--h;
				return (h & 3) != 0 ? h : h + 3;
			}
		};
		static_assert(sizeof(Face) == sizeof(HalfEdge) * 4, "Face struct must be 4x the size of the HalfEdge");
        struct EdgeInfo {
            HalfEdge::Index edge;
            Oryol::Array<HalfEdge::Index> constraints;
        };
		struct Vertex {
        public:
            Mesh * mesh;
			glm::dvec2 position;
            HalfEdge::Index edge; //Can be either incoming or outgoing edge
			size_t constraintCount;
			size_t endPointCount;
            inline HalfEdge::Index GetIncomingEdge() const {
                HalfEdge::Index vertexID = mesh->vertices.Distance(*this);
                return mesh->edgeAt(edge).destinationVertex == vertexID ? edge : mesh->edgeAt(edge).oppositeHalfEdge;
            }
            inline HalfEdge::Index GetOutgoingEdge() const {
                HalfEdge::Index vertexID = mesh->vertices.Distance(*this);
                return mesh->edgeAt(edge).destinationVertex != vertexID ? edge : mesh->edgeAt(edge).oppositeHalfEdge;
            }
            inline HalfEdge::Index GetNextIncomingEdge(HalfEdge::Index current) const {
                return mesh->edgeAt(Face::nextHalfEdge(current)).oppositeHalfEdge;
            }
            inline HalfEdge::Index GetNextOutgoingEdge(HalfEdge::Index current) const {
                return Face::nextHalfEdge(mesh->edgeAt(current).oppositeHalfEdge);
            }
            inline HalfEdge::Index GetPrevOutgoingEdge(HalfEdge::Index current) const {
                return mesh->edgeAt(Face::prevHalfEdge(current)).oppositeHalfEdge;
            }
            //TODO: Implement a get prev incoming edge method
		};
        struct ConstraintSegment {
            HalfEdge::Index startVertex;
            HalfEdge::Index endVertex;
            Oryol::Array<HalfEdge::Index> edgePairs;
        };
		struct ObjectRef {
            bool operator!() const {
                return type == None;
            }
            HalfEdge::Index object;
			enum Code { None, Vertex, Edge, Face } type;
            inline ObjectRef(size_t o, Code t): object(o), type(t) {}

		};

		//Initialises the Delaunay Triangulation with a square mesh with specified width and height
		//Creates 5 vertices, and 6 faces. Vertex with index 0 is an infinite vertex
		void Setup(double width, double height);
		//Inserts a vertex by splitting an existing face/edge or returning an existing vertex 
		//if one exists at the specified point
        size_t InsertVertex(const glm::dvec2 & p);

		bool RemoveVertex(const size_t vertexID);
		
        size_t InsertConstraintSegment(const glm::dvec2 & start, const glm::dvec2 & end);
        
		//Find which primitive the specified point is inside
		//Will only return primitives which are deemed to be "real"
        ObjectRef Locate(const glm::dvec2 & p);
        
        void SetDebugDraw(DebugDraw * debug);
        void DrawDebugData();

	private:
				
		
        HalfEdge & edgeAt(HalfEdge::Index index);
        const HalfEdge & edgeAt(HalfEdge::Index index) const;

		Geo2D::AABB boundingBox;
        ObjectPool<Face> faces;
        ObjectPool<Vertex> vertices;
        ObjectPool<ConstraintSegment> constraints;
        
        ObjectPool<EdgeInfo> edgeInfo;
        
        DebugDraw * debugDraw;

	};
    
    class DebugDraw {
    public:
        virtual void DrawVertex(glm::vec2 position) = 0;
        virtual void DrawEdge(glm::vec2 origin, glm::vec2 destination, bool constrained) = 0;
    };


}
