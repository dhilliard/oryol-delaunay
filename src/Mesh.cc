#include "Mesh.h"
#include "Geo2D.h"
constexpr double EPSILON = 0.01;
constexpr double EPSILON_SQUARED = 0.0001;
using namespace Geo2D;

//TODO: Investigate a better approach for initialising everything.
void Delaunay::Mesh::Setup(double width, double height)
{
	//Clear everything
	vertices.Clear();
	freeVertices.Clear();
	faces.Clear();
	freeFaces.Clear();
	//Infinite vertex
	Vertex & vInf = requestVertex(width*0.5, height * 0.5);
	//Other vertices in CW order
	Vertex & vBL = requestVertex(0, 0);
	Vertex & vTL = requestVertex(0, height);
	Vertex & vTR = requestVertex(width, height);
	Vertex & vBR = requestVertex(width, 0);
	//Request faces
	Face & fTL_BL_BR = requestFace(); //eBR_TL & eTL_BL & eBL_BR
	Face & fBR_TR_TL = requestFace(); //eTL_BR & eBR_TR & eTR_TL
	Face & fTL_TR_vInf = requestFace(); //eInf_TL & eTL_TR & eTR_Inf
	Face & fBL_TL_vInf = requestFace(); //eInf_BL & eBL_TL & eTL_Inf
	Face & fBR_BL_vInf = requestFace(); //eInf_BR & eBR_BL & eBL_Inf
	Face & fTR_BR_vInf = requestFace(); //eInf_TR & eTR_BR & eBR_Inf
	//Setup vertex data
	fTL_BL_BR.edges[0].destinationVertex = indexFor(vTL);
	fTL_BL_BR.edges[1].destinationVertex = indexFor(vBL);
	fTL_BL_BR.edges[2].destinationVertex = indexFor(vBR);
	
	fBR_TR_TL.edges[0].destinationVertex = indexFor(vBR);
	fBR_TR_TL.edges[1].destinationVertex = indexFor(vTR);
	fBR_TR_TL.edges[2].destinationVertex = indexFor(vTL);
	
	fTL_TR_vInf.edges[0].destinationVertex = indexFor(vTL);
	fTL_TR_vInf.edges[1].destinationVertex = indexFor(vTR);
	fTL_TR_vInf.edges[2].destinationVertex = indexFor(vInf);

	fBL_TL_vInf.edges[0].destinationVertex = indexFor(vBL);
	fBL_TL_vInf.edges[1].destinationVertex = indexFor(vTL);
	fBL_TL_vInf.edges[2].destinationVertex = indexFor(vInf);

	fBR_BL_vInf.edges[0].destinationVertex = indexFor(vBR);
	fBR_BL_vInf.edges[1].destinationVertex = indexFor(vBL);
	fBR_BL_vInf.edges[2].destinationVertex = indexFor(vInf);

	fTR_BR_vInf.edges[0].destinationVertex = indexFor(vTR);
	fTR_BR_vInf.edges[1].destinationVertex = indexFor(vBR);
	fTR_BR_vInf.edges[2].destinationVertex = indexFor(vInf);
	//Setup shared edges
	//Weld bottom of pyramid together
	fTL_BL_BR.edges[0].oppositeHalfEdge = indexFor(fBR_TR_TL, 0); // eBR_TL.opposite = eTL_BR
	fBR_TR_TL.edges[0].oppositeHalfEdge = indexFor(fTL_BL_BR, 0); // eTL_BR.opposite = eBR_TL
	//Weld outer triangles to inner square
	fTL_TR_vInf.edges[1].oppositeHalfEdge = indexFor(fBR_TR_TL, 2); //eTL_TR.opposite = eTR_TL
	fBR_TR_TL.edges[2].oppositeHalfEdge = indexFor(fTL_TR_vInf, 1); //eTR_TL.opposite = eTL_TR

	fBL_TL_vInf.edges[1].oppositeHalfEdge = indexFor(fTL_BL_BR, 1); //eBL_TL.opposite = eTL_BL
	fTL_BL_BR.edges[1].oppositeHalfEdge = indexFor(fBL_TL_vInf, 1); //eTL_BL.opposite = eBL_TL

	fBR_BL_vInf.edges[1].oppositeHalfEdge = indexFor(fTL_BL_BR, 2); //eBR_BL.opposite = eBL_BR
	fTL_BL_BR.edges[2].oppositeHalfEdge = indexFor(fBR_BL_vInf, 1); //eBL_BR.opposite = eBR_BL
	
	fTR_BR_vInf.edges[1].oppositeHalfEdge = indexFor(fBR_TR_TL, 1); //eTR_BR.opposite = eBR_TR
	fBR_TR_TL.edges[1].oppositeHalfEdge = indexFor(fTR_BR_vInf, 1); //eBR_TR.opposite = eTR_BR

	//Weld neighbouring edges of pyramid together
	fTL_TR_vInf.edges[0].oppositeHalfEdge = indexFor(fBL_TL_vInf, 2); //eInf_TL.opposite = eTL_Inf
	fBL_TL_vInf.edges[2].oppositeHalfEdge = indexFor(fTL_TR_vInf, 0); //eTL_Inf.opposite = eInf_TL
	
	fBL_TL_vInf.edges[0].oppositeHalfEdge = indexFor(fBR_BL_vInf, 2); //eInf_BL.opposite = eBL_Inf
	fBR_BL_vInf.edges[2].oppositeHalfEdge = indexFor(fBL_TL_vInf, 0); //eBL_Inf.opposite = eInf_BL
	
	fBR_BL_vInf.edges[0].oppositeHalfEdge = indexFor(fTR_BR_vInf, 2); // eInf_BR.opposite = eBR_Inf
	fTR_BR_vInf.edges[2].oppositeHalfEdge = indexFor(fBR_BL_vInf, 0); // eBR_Inf.opposite = eInf_BR

	fTR_BR_vInf.edges[0].oppositeHalfEdge = indexFor(fTL_TR_vInf, 2); // eInf_TR.opposite = eTR_Inf
	fTL_TR_vInf.edges[2].oppositeHalfEdge = indexFor(fTR_BR_vInf, 0); // eTR_Inf.opposite = eInf_TR

	//Set edge on each vertex
	vInf.incomingEdge = indexFor(fTL_TR_vInf, 2);
	vBL.incomingEdge = indexFor(fTL_BL_BR, 1);
	vBR.incomingEdge = indexFor(fTL_BL_BR, 2);
	vTR.incomingEdge = indexFor(fBR_TR_TL, 1);
	vTL.incomingEdge = indexFor(fBR_TR_TL, 2);
}

//TODO: Ensure that constraint state is also copied over when the new faces are created
Delaunay::Mesh::Index Delaunay::Mesh::FlipEdge(Index h)
{
	o_error("Not a valid HalfEdge index", h % 4 == 0);
	//The triangles owning these half-edges will be replaced with a triangle pair with a common edge running left to right
	HalfEdge & eUp_Down = edgeAt(h);
	HalfEdge & eDown_Up = edgeAt(eUp_Down.oppositeHalfEdge);

	//These belong to the two faces which will be removed so their opposites must be patched.
	HalfEdge & eDown_Right = edgeAt(nextHalfEdge(indexFor(eUp_Down)));
	HalfEdge & eUp_Left = edgeAt(nextHalfEdge(indexFor(eDown_Up)));
	HalfEdge & eRight_Up = edgeAt(prevHalfEdge(indexFor(eUp_Down)));
	HalfEdge & eLeft_Down = edgeAt(prevHalfEdge(indexFor(eDown_Up)));

	//Make sure we store references to all these vertices because it will make things much easier.
	Vertex & vDown = this->vertices[eUp_Down.destinationVertex];
	Vertex & vUp = this->vertices[eDown_Up.destinationVertex];
	Vertex & vLeft = this->vertices[eUp_Left.destinationVertex];
	Vertex & vRight = this->vertices[eDown_Right.destinationVertex];

	//Request new faces; Ideally you should be able to request a pair/triplet of faces at the same time
	Face & fLeft_Right_Up = requestFace();
	Face & fRight_Left_Down = requestFace();

	//Set up vertex references for both new faces
	fLeft_Right_Up.edges[0].destinationVertex = indexFor(vLeft);
	fLeft_Right_Up.edges[1].destinationVertex = indexFor(vRight);
	fLeft_Right_Up.edges[2].destinationVertex = indexFor(vUp);

	fRight_Left_Down.edges[0].destinationVertex = indexFor(vRight);
	fRight_Left_Down.edges[1].destinationVertex = indexFor(vLeft);
	fRight_Left_Down.edges[2].destinationVertex = indexFor(vDown);

	//Weld the two new faces together along eLeft_Right and eRight_Left
	fLeft_Right_Up.edges[1].oppositeHalfEdge = indexFor(fRight_Left_Down, 1);
	fRight_Left_Down.edges[1].oppositeHalfEdge = indexFor(fLeft_Right_Up, 1);

	//Patch the new faces with the correct opposite references
	fLeft_Right_Up.edges[0].oppositeHalfEdge = eUp_Left.oppositeHalfEdge;
	edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = indexFor(fLeft_Right_Up, 0);

	fLeft_Right_Up.edges[2].oppositeHalfEdge = eRight_Up.oppositeHalfEdge;
	edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = indexFor(fLeft_Right_Up, 2);

	fRight_Left_Down.edges[0].oppositeHalfEdge = eDown_Right.oppositeHalfEdge;
	edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = indexFor(fRight_Left_Down, 0);

	fRight_Left_Down.edges[2].oppositeHalfEdge = eLeft_Down.oppositeHalfEdge;
	edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = indexFor(fRight_Left_Down, 2);

	//Patch edge references in each vertex to refer to the newly created faces (if needed)
	if (indexFor(eDown_Right) == vRight.incomingEdge)
		vRight.incomingEdge = indexFor(fRight_Left_Down, 0);
	if (indexFor(eLeft_Down) == vDown.incomingEdge)
		vDown.incomingEdge = indexFor(fRight_Left_Down, 2);
	if (indexFor(eUp_Left) == vLeft.incomingEdge)
		vLeft.incomingEdge = indexFor(fLeft_Right_Up, 0);
	if (indexFor(eRight_Up) == vUp.incomingEdge)
		vUp.incomingEdge = indexFor(fLeft_Right_Up, 2);

	//Release the old faces
	freeFaces.Enqueue(eDown_Up.oppositeHalfEdge / 4); //Left face
	freeFaces.Enqueue(eUp_Down.oppositeHalfEdge / 4); //Right face;

	return indexFor(fLeft_Right_Up, 1); //Returns the Half-Edge running left to right for the newly created triangle pair
}

inline Delaunay::Mesh::HalfEdge & Delaunay::Mesh::edgeAt(Index index) {
	return *reinterpret_cast<HalfEdge*>(faces.begin() + index);
}

inline Delaunay::Mesh::Index Delaunay::Mesh::nextHalfEdge(Index h) {
	++h;
	return (h & 3) != 0 ? h : h - 3;
}

inline Delaunay::Mesh::Index Delaunay::Mesh::prevHalfEdge(Index h) {
	--h;
	return (h & 3) != 0 ? h : h + 3;
}

inline bool Delaunay::Mesh::isFaceReal(const Face & f) const {
	return f.edges[0].destinationVertex != 0 && f.edges[1].destinationVertex != 0 && f.edges[2].destinationVertex != 0;
}

inline bool Delaunay::Mesh::isHalfEdgeConstrained(Index h) {
	Face & f = faces[h / 4];
	constexpr Index offset = 8 * sizeof(Index) - 1;
	Index mask = 1 << (offset - h & 3);
	return f.flags & mask;
}

inline Delaunay::Mesh::Index Delaunay::Mesh::indexFor(Face & face) {
	return &face - faces.begin();
}

inline Delaunay::Mesh::Index Delaunay::Mesh::indexFor(Vertex & vertex) {
	return &vertex - vertices.begin();
}

inline Delaunay::Mesh::Index Delaunay::Mesh::indexFor(HalfEdge & edge) {
	return &edge - reinterpret_cast<HalfEdge*>(this->faces.begin());
}

//Compute index for edge within a face.
inline Delaunay::Mesh::Index Delaunay::Mesh::indexFor(Face & face, Index edge) {
	return indexFor(face) * 4 + edge + 1;
}

inline Delaunay::Mesh::Face & Delaunay::Mesh::requestFace() {
	if (!freeFaces.Empty())
		return faces[freeFaces.Dequeue()] = {};
	else
		return faces.Add();
}

inline Delaunay::Mesh::Vertex & Delaunay::Mesh::requestVertex(double x, double y) {
	if (!freeVertices.Empty())
		return vertices[freeVertices.Dequeue()] = { x,y,HalfEdge::InvalidIndex };
	else
		return vertices.Add({ x,y,HalfEdge::InvalidIndex });
}

Delaunay::Mesh::HalfEdge::HalfEdge()
	: destinationVertex(Vertex::InvalidIndex), oppositeHalfEdge(HalfEdge::InvalidIndex){}

Delaunay::Mesh::Face::Face() 
	: flags(0), matId(-1) {}


