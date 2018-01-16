#pragma once

#include "Core\Containers\Array.h"
#include "Core\Containers\Queue.h"
#include "Core\Containers\Set.h"
#include "glm/vec2.hpp"
//Uses concepts from 
// * https://infoscience.epfl.ch/record/100269/files/Kallmann_and_al_Geometric_Modeling_03.pdf -> For the overall implementation strategy
// * http://www.dtecta.com/files/GDC17_VanDenBergen_Gino_Brep_Triangle_Meshes.pdf -> For the low level Face/HalfEdge data structure
// * http://2.3jachtuches.pagesperso-orange.fr/dossiers/triangul/doc/fast.pdf -> For Mesh::Locate()
//Minimal delaunay implementation would implement Insert/DeleteVertex functions
//Minimal constrained implementation would implement Insert/DeleteConstraintSegment functions
//TODO: Write set of iterators for Mesh: Vertices/Edges/Faces
//TODO: Replace index type returned from modifier functions with a custom handle type with the 2 MSBs reserved for type.

namespace Delaunay {
	typedef uint32_t Index;
	class Mesh;


	struct IncomingHalfEdgeIterator {
		Index operator*();
		void operator++();
		bool operator!=(const IncomingHalfEdgeIterator & rhs);
		IncomingHalfEdgeIterator & begin();
		IncomingHalfEdgeIterator & end();
		IncomingHalfEdgeIterator(Mesh & mesh, Index vertex);
		const Mesh & mesh;
	private:
		Index current;
		Index first;
	};

	struct OutgoingHalfEdgeIterator {
		Index operator*();
		void operator++();
		bool operator!=(const OutgoingHalfEdgeIterator & rhs);
		OutgoingHalfEdgeIterator & begin();
		OutgoingHalfEdgeIterator & end();
		OutgoingHalfEdgeIterator(Mesh & mesh, Index vertex);
		const Mesh & mesh;
	private:
		Index current;
		Index first;
	};

	struct NeighboringVertexIterator {
		Index operator*();
		void operator++();
		bool operator!=(NeighboringVertexIterator & rhs);
		NeighboringVertexIterator & begin();
		NeighboringVertexIterator & end();
		NeighboringVertexIterator(Mesh & mesh, Index vertex);
	private:
		OutgoingHalfEdgeIterator impl;
	};

	class DebugDraw;

	class Mesh {
	public:
		
		struct LocateResult {
			Index object;
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
		Index InsertVertex(double x, double y);
		//Deletes the vertex with the specified index only if no constraints depend on said vertex
		bool DeleteVertex(Index index);

		Index InsertConstraintSegment(double x1, double y1, double x2, double y2);
		void DeleteConstraintSegment(Index index);
		
		//Find which primitive the specified point is inside
		//Will only return primitives which are deemed to be "real"
		LocateResult Locate(double x, double y);

		void SetDebugDraw(DebugDraw * debug);
		void DrawDebugData();
	private:

		friend struct IncomingHalfEdgeIterator;
		friend struct OutgoingHalfEdgeIterator;
		friend struct NeighboringVertexIterator;

		struct Vertex {
			glm::dvec2 position;
			Index edge; //Can be either incoming or outgoing edge
			static constexpr Index InvalidIndex = -1;
		};
		struct HalfEdge {
			Index destinationVertex; //End Vertex Index
			Index oppositeHalfEdge; //Opposite HalfEdge
			static constexpr Index InvalidIndex = -1;
			HalfEdge();
		};
		struct Face {
			//The highest 3 bits of flags is used to determine which edges are constrained
			Index flags;
			Index matId;
			HalfEdge edges[3];
			static constexpr Index InvalidIndex = -1;
			Face();
			//Returns next (CCW) half-edge
			static inline Index nextHalfEdge(Index h);
			//Returns prev (CW) half-edge
			static inline Index prevHalfEdge(Index h);
		};
		static_assert(sizeof(Face) == sizeof(HalfEdge) * 4, "Face struct must be 4x the size of the HalfEdge");
		struct ConstraintSegment {
			Oryol::Array<Index> vertices;
		};
		//Deletes 2 faces + creates 4 faces
		//Returns new vertex that is created as a result of splitting the edge
		Index splitEdge(Index h, double x, double y);
		//Returns halfedge from new face pair created as a result of the flip
		Index flipEdge(Index h);
		//Deletes one face + creates 3 new faces
		//Returns new vertex that is created as a result of splitting the face
		Index splitFace(Index f, double x, double y);
		//Must be called after Split/Flip functions are called to restore delaunay condition
		void RestoreAsDelaunay();

		LocateResult isInFace(double x, double y, Face & face);
		bool isDelaunay(Index h);

		inline HalfEdge & edgeAt(Index index);
		inline const HalfEdge & edgeAt(Index index) const;
		
		//Utilities for computing index based on an object's type.
		inline Index indexFor(Face & face);
		inline Index indexFor(Vertex & vertex);
		inline Index indexFor(HalfEdge & edge);
		inline Index indexFor(Face & face, Index edge);
		inline bool isHalfEdgeConstrained(Index h) {
			Face & f = this->faces[h / 4];
			return f.flags & (1 << (h & 3));
		}
		inline bool isHalfEdgeReal(Index h) {
			Face & f = this->faces[h / 4];
			return (f.flags & 0x1) == 0;
		}
		Face & requestFace();
		void recycleFace(Index f);
		Vertex & requestVertex(double x, double y,bool cache = true);

		Oryol::Array<Face> faces;
		Oryol::Queue<Index> freeFaces;
		//First vertex is a special vertex as it is a _infinite_ vertex
		Oryol::Array<Vertex> vertices;
		Oryol::Queue<Index> freeVertices;
		Oryol::Set<Index> cachedVertices;

		DebugDraw * debugDraw;

		Oryol::Array<Index> edgesToCheck;
		Index centerVertex;
	};
	class DebugDraw {
	public:
		virtual void DrawVertex(glm::vec2 position) = 0;
		virtual void DrawEdge(glm::vec2 origin, glm::vec2 destination, bool constrained) = 0;
	};

}