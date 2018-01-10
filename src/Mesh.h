#pragma once
#include "Core\Containers\Array.h"
#include "Core\Containers\Queue.h"
#include "Geo2D.h"
//Uses concepts from 
// * https://infoscience.epfl.ch/record/100269/files/Kallmann_and_al_Geometric_Modeling_03.pdf
// * http://www.dtecta.com/files/GDC17_VanDenBergen_Gino_Brep_Triangle_Meshes.pdf
//Minimal delaunay implementation would implement Insert/DeleteVertex functions
//Minimal constrained implementation would implement Insert/DeleteConstraintSegment functions
namespace Delaunay {
	

	class Mesh {
	public:
		typedef uint32_t Index;
		struct LocateResult {
			Index object;
			enum Code {None, Vertex, Edge, Face} type;
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
		//Deletes 2 faces + creates 4 faces
		//Returns new vertex that is created as a result of splitting the edge
		Index SplitEdge(Index h, double x, double y);
		//Returns halfedge from new face pair created as a result of the flip
		Index FlipEdge(Index h);
		//Deletes one face + creates 3 new faces
		//Returns new vertex that is created as a result of splitting the face
		Index SplitFace(Index f, double x, double y);
		//Must be called after Split/Flip functions are called to restore delaunay condition
		void RestoreAsDelaunay();
		//Find which primitive the specified point is inside
		//Will only return primitives which are deemed to be "real"
		LocateResult Locate(double x, double y);
	private:
		
		struct Vertex {
			double x, y;
			Index incomingEdge;
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
		};
		static_assert(sizeof(Face) == sizeof(HalfEdge) * 4, "Face struct must be 4x the size of the HalfEdge");

		inline HalfEdge & edgeAt(Index index);
		//Returns next (CCW) half-edge
		inline Index nextHalfEdge(Index h);
		//Returns prev (CW) half-edge
		inline Index prevHalfEdge(Index h);
		inline bool isFaceReal(const Face & f) const;
		inline bool isHalfEdgeConstrained(Index h);
		inline Index indexFor(Face & face);
		inline Index indexFor(Vertex & vertex);
		//Compute index for edge within a face.
		inline Index indexFor(Face & face, Index edge);
		Face & requestFace();
		Vertex & requestVertex(double x, double y);

		Oryol::Array<Face> faces;
		Oryol::Queue<Index> freeFaces;
		//First vertex is a special vertex as it is a _infinite_ vertex
		Oryol::Array<Vertex> vertices;
		Oryol::Queue<Index> freeVertices;
	};
	

}