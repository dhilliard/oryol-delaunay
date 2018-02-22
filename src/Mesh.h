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
//TODO
// * Refactor Face/Half-Edge into MeshTypes.h



namespace Delaunay {
    
	class Mesh {
	public:
        
        
		struct HalfEdge {
            typedef uint32_t Index;
            enum {
                IndexBits = sizeof(Index) * 8
            };
			Index destinationVertex; //End Vertex Index
			Index oppositeHalfEdge; //Opposite HalfEdge
            Index constrained : 1; //Cache constraint state on the half edge
            Index edgePair : IndexBits - 1; //Use this to references EdgeInfo struct
			static const Index InvalidIndex = -1;
            //It may be somewhat hacky but you can obtain the address of the face that this edge belongs to by using;
            //reinterpret_cast<Face*>(reinterpret_cast<size_t>(this)&~((1<<4)-1))
		};
		struct Face {
            HalfEdge::Index flags;
            HalfEdge::Index matID;
            HalfEdge::Index userData;
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

		struct Vertex {
        public:
			glm::dvec2 position;
            HalfEdge::Index edge; //Can be either incoming or outgoing edge
			size_t constraintCount;
			size_t endPointCount;
		};
        struct ConstraintSegment {
            HalfEdge::Index startVertex;
            HalfEdge::Index endVertex;
            Oryol::Array<HalfEdge::Index> edgePairs;
        };
        //Always use this when providing references to internal objects
        struct LocateRef {
            bool operator!() const {
                return type == None;
            }
            uint32_t object;
            enum Code { None, Vertex, Edge, Face } type;
            inline LocateRef(uint32_t o, Code t): object(o), type(t) {}
            inline LocateRef() : object(-1), type(Code::None) {}
        };
		

		//Initialises the Delaunay Triangulation with a square mesh with specified width and height
		//Creates 5 vertices, and 6 faces. Vertex with index 0 is an infinite vertex
		void Setup(double width, double height);
        //Inserts a vertex by splitting an existing face/edge or returning an existing vertex
        //if one exists at the specified point
        uint32_t InsertVertex(const glm::dvec2 & p);
        bool RemoveVertex(const uint32_t vertexID);
        
        uint32_t InsertConstraintSegment(const glm::dvec2 & start, const glm::dvec2 & end);
        void RemoveConstraintSegment(const uint32_t constraintID);
        
        //Find which primitive the specified point is inside
        //Will only return primitives which are deemed to be "real"
        LocateRef Locate(const glm::dvec2 & p) const;

		bool CircleIntersectsConstraints(const glm::dvec2 & center, double radius) const;
        
        inline HalfEdge::Index GetIncomingEdgeFor(uint32_t vertexID) const {
            const Vertex & vertex = vertices[vertexID];
            return EdgeAt(vertex.edge).destinationVertex == vertexID ? vertex.edge : EdgeAt(vertex.edge).oppositeHalfEdge;
        }
        inline HalfEdge::Index GetOutgoingEdgeFor(uint32_t vertexID) const {
            const Vertex & vertex = vertices[vertexID];
            return EdgeAt(vertex.edge).destinationVertex != vertexID ? vertex.edge : EdgeAt(vertex.edge).oppositeHalfEdge;
        }
        inline HalfEdge::Index GetNextIncomingEdge(HalfEdge::Index current) const {
            return EdgeAt(Face::nextHalfEdge(current)).oppositeHalfEdge;
        }
        inline HalfEdge::Index GetNextOutgoingEdge(HalfEdge::Index current) const {
            return Face::nextHalfEdge(EdgeAt(current).oppositeHalfEdge);
        }
        inline HalfEdge::Index GetPrevOutgoingEdge(HalfEdge::Index current) const {
            return EdgeAt(Face::prevHalfEdge(current)).oppositeHalfEdge;
        }
        inline HalfEdge::Index GetPrevIncomingEdge(HalfEdge::Index current) const {
            return Face::prevHalfEdge(EdgeAt(current).oppositeHalfEdge);
        }
        inline const HalfEdge & EdgeAt(HalfEdge::Index index) const {
            //return faces[index / 4].edges[(index & 3) - 1];
            return faces.GetAs<const HalfEdge>(index);
        }
        inline const Vertex & VertexAt(uint32_t index) const {
            return vertices[index];
        }
        inline const Face & FaceAt(uint32_t index) const {
            return faces[index];
        }
        inline const ConstraintSegment & SegmentAt(uint32_t index) const {
            return segments[index];
        }
        const Oryol::Set<HalfEdge::Index> & ActiveVertexIndices() const {
            return vertices.ActiveIndices();
        }
        const Geo2D::AABB & GetBoundingBox() const {
            return boundingBox;
        }
		
	private:
        struct Impl;
        struct EdgeInfo {
            HalfEdge::Index edge;
            Oryol::Set<HalfEdge::Index> constraints;
        };
        HalfEdge & edgeAt(HalfEdge::Index index);
        
		Geo2D::AABB boundingBox;
        ObjectPool<Face> faces;
        ObjectPool<Vertex> vertices;
        ObjectPool<ConstraintSegment> segments;
        ObjectPool<EdgeInfo> edgeInfo;


	};

}
