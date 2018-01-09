#pragma once
#include "Core\Containers\Array.h"
#include "Core\Containers\Queue.h"
#include "Geo2D.h"
//Uses concepts from 
// * http://totologic.blogspot.com.au/2013/11/core-quad-edge-implementation-explained.html
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
		void Setup(double width, double height);
		
		Index InsertVertex(double x, double y);
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
		void RestoreAsDelaunay();
		LocateResult Locate(double x, double y);
	private:
		
		struct Vertex {
			Index edge;
			double x, y;
		};
		struct HalfEdge {
			Index destinationVertex; //End Vertex Index
			Index oppositeHalfEdge; //Opposite HalfEdge
		};
		//TODO: store constraint bitmask in flags inside each face
		struct Face {
			Index flags;
			Index matId;
			HalfEdge edges[3];
		};
		static_assert(sizeof(Face) == sizeof(HalfEdge) * 4, "Face struct must be 4x the size of the HalfEdge");

		inline HalfEdge & edgeAt(Index index) {
			return *reinterpret_cast<HalfEdge*>(faces.begin() + index);
		}
		//Returns next (CCW) half-edge
		inline Index nextHalfEdge(Index h) {
			++h;
			return (h & 3) != 0 ? h : h - 3;
		}
		//Returns prev (CW) half-edge
		inline Index prevHalfEdge(Index h) {
			--h;
			return (h & 3) != 0 ? h : h + 3;
		}
		
		inline bool isHalfEdgeConstrained(Index h) {
			return false;
		}
		Face & requestFace() {
			if (!freeFaces.Empty())
				return faces[freeFaces.Dequeue()];
			else
				return faces.Add();
		}
		Vertex & requestVertex() {
			if (!freeVertices.Empty())
				return vertices[freeVertices.Dequeue()];
			else
				return vertices.Add();
		}
		Oryol::Array<Face> faces;
		Oryol::Queue<Index> freeFaces;
		Oryol::Array<Vertex> vertices;
		Oryol::Queue<Index> freeVertices;
	};
	

}