#include "Mesh.h"
#include "Geo2D.h"
#include <cmath>
#include <limits>
#include "Core/Containers/Set.h"
#include "Core/Containers/Queue.h"
#include "glm/geometric.hpp"

constexpr double EPSILON = 0.01;
constexpr double EPSILON_SQUARED = 0.0001;
using namespace Geo2D;
using namespace Delaunay;
using Index = Mesh::HalfEdge::Index;

int rand_range(int min, int max) {
	return min + (float)(rand() / RAND_MAX)*(max - min);
}

struct Mesh::Impl {
    static Mesh::LocateResult IsInFace(Mesh & mesh, double x, double y, Face & face)
    {
        LocateResult result{ HalfEdge::InvalidIndex, LocateResult::None };
        size_t faceIndex = &face - mesh.faces.begin();
        glm::dvec2 p = { x,y };
        
        glm::dvec2 v1 = mesh.vertices[face.edges[0].destinationVertex].position;
        glm::dvec2 v2 = mesh.vertices[face.edges[1].destinationVertex].position;
        glm::dvec2 v3 = mesh.vertices[face.edges[2].destinationVertex].position;
        
        if (Geo2D::Sign(v3, v1, p) >= 0.0 && Geo2D::Sign(v1, v2, p) >= 0.0 && Geo2D::Sign(v2, v3, p) >= 0.0) {
            //Cache the proximity info rather than calculate it more than necessary
            bool proximity[3] = {
                Geo2D::DistanceSquaredPointToLineSegment(v3,v1,p) <= EPSILON_SQUARED,
                Geo2D::DistanceSquaredPointToLineSegment(v1,v2,p) <= EPSILON_SQUARED,
                Geo2D::DistanceSquaredPointToLineSegment(v2,v3,p) <= EPSILON_SQUARED
            };
            
            if (proximity[0]) {
                if (proximity[1]) {
                    result.object = face.edges[0].destinationVertex;
                    result.type = LocateResult::Vertex;
                }
                else if (proximity[2]) {
                    result.object = face.edges[2].destinationVertex;
                    result.type = LocateResult::Vertex;
                }
                else {
                    result.object = faceIndex * 4 + 1; //eV3_V1
                    result.type = LocateResult::Edge;
                }
            }
            else if (proximity[1]) {
                if (proximity[2]) {
                    result.object = face.edges[1].destinationVertex;
                    result.type = LocateResult::Vertex;
                }
                else {
                    result.object = faceIndex * 4 + 2; //eV1_V2
                    result.type = LocateResult::Edge;
                }
            }
            else if (proximity[2]) {
                result.object = faceIndex * 4 + 3; //eV2_V3
                result.type = LocateResult::Edge;
            }
            else {
                result.object = faceIndex;
                result.type = LocateResult::Face;
            }
        }
        return result;
    }
/*
	//Deletes 2 faces + creates 4 faces
	//Returns new vertex that is created as a result of splitting the edge
	static Index SplitEdge(Mesh & mesh, Index h, double x, double y, Index & centerVertex, Oryol::Array<Index> & edgesToCheck)
	{

		edgesToCheck.Clear();
		//Have to specify this otherwise code that uses the indexFor() shortcut for generating indices will fail.
		//The ideal approach would be to use indices the entire way through and would avoid this particular "hack"
		mesh.faces.Reserve(4);
		HalfEdge & eUp_Down = mesh.edgeAt(h);
		HalfEdge & eDown_Up = mesh.edgeAt(eUp_Down.oppositeHalfEdge);

		Vertex & vDown = mesh.vertices[eUp_Down.destinationVertex];
		Vertex & vUp = mesh.vertices[eDown_Up.destinationVertex];
		//Check to see if the point is close enough to a vertex and return the corresponding halfedge
		if (Geo2D::DistanceSquared(vUp.position - glm::dvec2{ x, y }) <= EPSILON_SQUARED)
			return h;
		if (Geo2D::DistanceSquared(vDown.position - glm::dvec2{ x,y }) <= EPSILON_SQUARED)
			return eUp_Down.oppositeHalfEdge;
		//These are the edges which surround the pair of faces we will be removing
		HalfEdge & eDown_Right = mesh.edgeAt(Face::nextHalfEdge(mesh.indexFor(eUp_Down)));
		//Oryol::Log::Info("eDown_Right: %d Opposite: %d\n", indexFor(eDown_Right),eDown_Right.oppositeHalfEdge);
		HalfEdge & eUp_Left = mesh.edgeAt(Face::nextHalfEdge(mesh.indexFor(eDown_Up)));
		//Oryol::Log::Info("eUp_Left: %d Opposite: %d\n", indexFor(eUp_Left), eUp_Left.oppositeHalfEdge);
		HalfEdge & eRight_Up = mesh.edgeAt(Face::prevHalfEdge(mesh.indexFor(eUp_Down)));
		//Oryol::Log::Info("eRight_Up: %d Opposite: %d\n", indexFor(eRight_Up), eRight_Up.oppositeHalfEdge);
		HalfEdge & eLeft_Down = mesh.edgeAt(Face::prevHalfEdge(mesh.indexFor(eDown_Up)));

		//Make sure we store references to all these vertices because it will make things much easier.
		Vertex & vLeft = mesh.vertices[eUp_Left.destinationVertex];
		Vertex & vRight = mesh.vertices[eDown_Right.destinationVertex];

		//Create the 4 new faces that we need
		Face & fUp_Left_Center = mesh.requestFace();
		Face & fLeft_Down_Center = mesh.requestFace();
		Face & fDown_Right_Center = mesh.requestFace();
		Face & fRight_Up_Center = mesh.requestFace();
		//Create the new vertex we need
		//As the point may not be exactly on the line it must be projected on the line segment
		glm::dvec2 p = Geo2D::OrthogonallyProjectPointOnLineSegment(vDown.position, vUp.position, { x,y });
		Vertex & vCenter = mesh.requestVertex(p.x, p.y);
		vCenter.edge = mesh.indexFor(fUp_Left_Center, 2);
		//Set vertex data on the face edges
		fUp_Left_Center.edges[0].destinationVertex = mesh.indexFor(vUp);
		fUp_Left_Center.edges[1].destinationVertex = mesh.indexFor(vLeft);
		fUp_Left_Center.edges[2].destinationVertex = mesh.indexFor(vCenter);

		fLeft_Down_Center.edges[0].destinationVertex = mesh.indexFor(vLeft);
		fLeft_Down_Center.edges[1].destinationVertex = mesh.indexFor(vDown);
		fLeft_Down_Center.edges[2].destinationVertex = mesh.indexFor(vCenter);

		fDown_Right_Center.edges[0].destinationVertex = mesh.indexFor(vDown);
		fDown_Right_Center.edges[1].destinationVertex = mesh.indexFor(vRight);
		fDown_Right_Center.edges[2].destinationVertex = mesh.indexFor(vCenter);

		fRight_Up_Center.edges[0].destinationVertex = mesh.indexFor(vRight);
		fRight_Up_Center.edges[1].destinationVertex = mesh.indexFor(vUp);
		fRight_Up_Center.edges[2].destinationVertex = mesh.indexFor(vCenter);

		//Weld the new faces together
		fUp_Left_Center.edges[2].oppositeHalfEdge = mesh.indexFor(fLeft_Down_Center, 0); //eLeft_Center.opposite = eCenter_Left
		fLeft_Down_Center.edges[0].oppositeHalfEdge = mesh.indexFor(fUp_Left_Center, 2); //eCenter_Left.opposite = eLeft_Center

		fLeft_Down_Center.edges[2].oppositeHalfEdge = mesh.indexFor(fDown_Right_Center, 0); //eDown_Center.opposite = eCenter_Down
		fDown_Right_Center.edges[0].oppositeHalfEdge = mesh.indexFor(fLeft_Down_Center, 2); //eCenter_Down.opposite = eDown_Center

		fDown_Right_Center.edges[2].oppositeHalfEdge = mesh.indexFor(fRight_Up_Center, 0); //eRight_Center.opposite = eCenter_Right
		fRight_Up_Center.edges[0].oppositeHalfEdge = mesh.indexFor(fDown_Right_Center, 2); //eCenter_Right.opposite = eRight_Center

		fRight_Up_Center.edges[2].oppositeHalfEdge = mesh.indexFor(fUp_Left_Center, 0); //eUp_Center.opposite = eCenter_Up
		fUp_Left_Center.edges[0].oppositeHalfEdge = mesh.indexFor(fRight_Up_Center, 2); //eCenter_Up.opposite = eUp_Center

																				   //Weld the new faces to the existing faces
		fLeft_Down_Center.edges[1].oppositeHalfEdge = eLeft_Down.oppositeHalfEdge;
		mesh.edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fLeft_Down_Center, 1);

		fRight_Up_Center.edges[1].oppositeHalfEdge = eRight_Up.oppositeHalfEdge;
		mesh.edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fRight_Up_Center, 1);

		fUp_Left_Center.edges[1].oppositeHalfEdge = eUp_Left.oppositeHalfEdge;
		mesh.edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fUp_Left_Center, 1);

		fDown_Right_Center.edges[1].oppositeHalfEdge = eDown_Right.oppositeHalfEdge;
		mesh.edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fDown_Right_Center, 1);

		//Patch edge references in each vertex to refer to the newly created faces (if needed)
		if (mesh.indexFor(eDown_Right) == vRight.edge)
			vRight.edge = mesh.indexFor(fDown_Right_Center, 1);
		if (mesh.indexFor(eLeft_Down) == vDown.edge || mesh.indexFor(eUp_Down) == vDown.edge)
			vDown.edge = mesh.indexFor(fLeft_Down_Center, 1);
		if (mesh.indexFor(eUp_Left) == vLeft.edge)
			vLeft.edge = mesh.indexFor(fUp_Left_Center, 1);
		if (mesh.indexFor(eRight_Up) == vUp.edge || mesh.indexFor(eDown_Up) == vUp.edge)
			vUp.edge = mesh.indexFor(fRight_Up_Center, 1);

		//Copy constraint state to the newly constructed faces
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(eDown_Right)))
			fDown_Right_Center.flags |= (1 << 2);
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(eLeft_Down)))
			fLeft_Down_Center.flags |= (1 << 2);
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(eUp_Left)))
			fUp_Left_Center.flags |= (1 << 2);
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(eRight_Up)))
			fRight_Up_Center.flags |= (1 << 2);

		if (!mesh.isHalfEdgeReal(eDown_Up.oppositeHalfEdge)) {
			//eUp_Down is not real therefore tag the appropriate new faces as not real
			fRight_Up_Center.flags |= 0x1;
			fDown_Right_Center.flags |= 0x1;
		}
		else if (!mesh.isHalfEdgeReal(eUp_Down.oppositeHalfEdge)) {
			//eDown_Up is not real therefore tag the new faces as not real 
			fUp_Left_Center.flags |= 0x1;
			fLeft_Down_Center.flags |= 0x1;
		}
		//Handle the fact that the split edge is constrained
		if (mesh.isHalfEdgeConstrained(h)) {
			//Tag the half-edges as constrained
			fRight_Up_Center.flags |= (1 << 1);
			fDown_Right_Center.flags |= (1 << 1);
			fUp_Left_Center.flags |= (1 << 1);
			fDown_Right_Center.flags |= (1 << 1);

			//Find the constraint segments associated with the original half-edge and insert the new vertex between the two existing vertices
			//A typical approach would be to have a list on each half-edge but since that wouldn't work with the Face == 4xHalf-Edge constraint
			//We check the Map for vertices to constraint segments
		}
		//Recycle the old faces
		mesh.recycleFace(eUp_Down.oppositeHalfEdge / 4);
		mesh.recycleFace(h / 4);
		//Tag the edges for delaunay condition checking
		centerVertex = mesh.indexFor(vCenter);

		edgesToCheck.Add(mesh.indexFor(fRight_Up_Center, 1));
		edgesToCheck.Add(mesh.indexFor(fUp_Left_Center, 1));
		edgesToCheck.Add(mesh.indexFor(fDown_Right_Center, 1));
		edgesToCheck.Add(mesh.indexFor(fRight_Up_Center, 1));

		return centerVertex;
	}
	//Returns halfedge from new face pair created as a result of the flip
	static Index FlipEdge(Mesh & mesh, Index h) {
		if (h % 4 == 0)
			o_error("Not a valid half-edge index");
		mesh.faces.Reserve(2);
		//The triangles owning these half-edges will be replaced with a triangle pair with a common edge running left to right
		HalfEdge & eUp_Down = mesh.edgeAt(h);
		HalfEdge & eDown_Up = mesh.edgeAt(eUp_Down.oppositeHalfEdge);

		//These belong to the two faces which will be removed so their opposites must be patched.
		HalfEdge & eDown_Right = mesh.edgeAt(Face::nextHalfEdge(mesh.indexFor(eUp_Down)));
		HalfEdge & eUp_Left = mesh.edgeAt(Face::nextHalfEdge(mesh.indexFor(eDown_Up)));
		HalfEdge & eRight_Up = mesh.edgeAt(Face::prevHalfEdge(mesh.indexFor(eUp_Down)));
		HalfEdge & eLeft_Down = mesh.edgeAt(Face::prevHalfEdge(mesh.indexFor(eDown_Up)));

		//Make sure we store references to all these vertices because it will make things much easier.
		Vertex & vDown = mesh.vertices[eUp_Down.destinationVertex];
		Vertex & vUp = mesh.vertices[eDown_Up.destinationVertex];
		Vertex & vLeft = mesh.vertices[eUp_Left.destinationVertex];
		Vertex & vRight = mesh.vertices[eDown_Right.destinationVertex];

		//Request new faces; Ideally you should be able to request a pair/triplet of faces at the same time
		Face & fLeft_Right_Up = mesh.requestFace();
		Face & fRight_Left_Down = mesh.requestFace();

		//Set up vertex references for both new faces
		fLeft_Right_Up.edges[0].destinationVertex = mesh.indexFor(vLeft);
		fLeft_Right_Up.edges[1].destinationVertex = mesh.indexFor(vRight);
		fLeft_Right_Up.edges[2].destinationVertex = mesh.indexFor(vUp);

		fRight_Left_Down.edges[0].destinationVertex = mesh.indexFor(vRight);
		fRight_Left_Down.edges[1].destinationVertex = mesh.indexFor(vLeft);
		fRight_Left_Down.edges[2].destinationVertex = mesh.indexFor(vDown);

		//Weld the two new faces together along eLeft_Right and eRight_Left
		fLeft_Right_Up.edges[1].oppositeHalfEdge = mesh.indexFor(fRight_Left_Down, 1);
		fRight_Left_Down.edges[1].oppositeHalfEdge = mesh.indexFor(fLeft_Right_Up, 1);

		//Patch the new faces with the correct opposite references
		fLeft_Right_Up.edges[0].oppositeHalfEdge = eUp_Left.oppositeHalfEdge;
		mesh.edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fLeft_Right_Up, 0);

		fLeft_Right_Up.edges[2].oppositeHalfEdge = eRight_Up.oppositeHalfEdge;
		mesh.edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fLeft_Right_Up, 2);

		fRight_Left_Down.edges[0].oppositeHalfEdge = eDown_Right.oppositeHalfEdge;
		mesh.edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fRight_Left_Down, 0);

		fRight_Left_Down.edges[2].oppositeHalfEdge = eLeft_Down.oppositeHalfEdge;
		mesh.edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fRight_Left_Down, 2);

		//Patch edge references in each vertex to refer to the newly created faces (if needed)
		if (mesh.indexFor(eDown_Right) == vRight.edge)
			vRight.edge = mesh.indexFor(fRight_Left_Down, 0);
		if (mesh.indexFor(eLeft_Down) == vDown.edge || mesh.indexFor(eUp_Down) == vDown.edge)
			vDown.edge = mesh.indexFor(fRight_Left_Down, 2);
		if (mesh.indexFor(eUp_Left) == vLeft.edge)
			vLeft.edge = mesh.indexFor(fLeft_Right_Up, 0);
		if (mesh.indexFor(eRight_Up) == vUp.edge || mesh.indexFor(eDown_Up) == vUp.edge);
		vUp.edge = mesh.indexFor(fLeft_Right_Up, 2);
		//Copy constraint state to the newly constructed faces
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(eDown_Right)))
			fRight_Left_Down.flags |= (1 << 1);
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(eLeft_Down)))
			fRight_Left_Down.flags |= (1 << 3);
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(eUp_Left)))
			fLeft_Right_Up.flags |= (1 << 1);
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(eRight_Up)))
			fLeft_Right_Up.flags |= (1 << 3);

		//Release the old faces
		mesh.recycleFace(mesh.indexFor(eDown_Up) / 4); //Left face
		mesh.recycleFace(mesh.indexFor(eUp_Down) / 4); //Right face;

		return mesh.indexFor(fLeft_Right_Up, 1); //Returns the Half-Edge running left to right for the newly created triangle pair
	}

	//TODO: Replace any vertex references with their corresponding index to avoid unnecessary calculations
	//Deletes one face + creates 3 new faces
	//Returns new vertex that is created as a result of splitting the face
	static Index SplitFace(Mesh & mesh, Index f, double x, double y, Oryol::Array<Index> & edgesToCheck) {
		edgesToCheck.Clear();
		mesh.faces.Reserve(3);
		Face & fA_B_New = mesh.requestFace();
		Face & fB_C_New = mesh.requestFace();
		Face & fC_A_New = mesh.requestFace();
		Face & fA_B_C = mesh.faces[f];

		Vertex & vNew = mesh.requestVertex(x, y);
		Vertex & vA = mesh.vertices[fA_B_C.edges[0].destinationVertex];
		Vertex & vB = mesh.vertices[fA_B_C.edges[1].destinationVertex];
		Vertex & vC = mesh.vertices[fA_B_C.edges[2].destinationVertex];

		vNew.edge = mesh.indexFor(fA_B_New, 2);

		fA_B_New.edges[0].destinationVertex = mesh.indexFor(vA);
		fA_B_New.edges[1].destinationVertex = mesh.indexFor(vB);
		fA_B_New.edges[2].destinationVertex = mesh.indexFor(vNew);

		fB_C_New.edges[0].destinationVertex = mesh.indexFor(vB);
		fB_C_New.edges[1].destinationVertex = mesh.indexFor(vC);
		fB_C_New.edges[2].destinationVertex = mesh.indexFor(vNew);

		fC_A_New.edges[0].destinationVertex = mesh.indexFor(vC);
		fC_A_New.edges[1].destinationVertex = mesh.indexFor(vA);
		fC_A_New.edges[2].destinationVertex = mesh.indexFor(vNew);

		//Copy eC_A from old face
		fC_A_New.edges[1].oppositeHalfEdge = fA_B_C.edges[0].oppositeHalfEdge;
		mesh.edgeAt(fA_B_C.edges[0].oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fC_A_New, 1);

		//Copy eA_B from old face
		fA_B_New.edges[1].oppositeHalfEdge = fA_B_C.edges[1].oppositeHalfEdge;
		mesh.edgeAt(fA_B_C.edges[1].oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fA_B_New, 1);

		//Copy eB_C from old face
		fB_C_New.edges[1].oppositeHalfEdge = fA_B_C.edges[2].oppositeHalfEdge;
		mesh.edgeAt(fA_B_C.edges[2].oppositeHalfEdge).oppositeHalfEdge = mesh.indexFor(fB_C_New, 1);


		fA_B_New.edges[0].oppositeHalfEdge = mesh.indexFor(fC_A_New, 2); //eNew_A.opposite = eA_New
		fC_A_New.edges[2].oppositeHalfEdge = mesh.indexFor(fA_B_New, 0); //eA_New.opposite = eNew_A

		fB_C_New.edges[0].oppositeHalfEdge = mesh.indexFor(fA_B_New, 2); //eNew_B.opposite = eB_New
		fA_B_New.edges[2].oppositeHalfEdge = mesh.indexFor(fB_C_New, 0); //eB_New.opposite = eNew_B

		fC_A_New.edges[0].oppositeHalfEdge = mesh.indexFor(fB_C_New, 2); //eNew_C.opposite = eC_New
		fB_C_New.edges[2].oppositeHalfEdge = mesh.indexFor(fC_A_New, 0); //eC_New.opposite = eNew_C

		//Patch edge references in each vertex to refer to the newly created faces (if needed)
		if (mesh.indexFor(fA_B_C, 0) == vA.edge)
			vA.edge = mesh.indexFor(fC_A_New, 1);
		if (mesh.indexFor(fA_B_C, 1) == vB.edge)
			vB.edge = mesh.indexFor(fA_B_New, 1);
		if (mesh.indexFor(fA_B_C, 2) == vC.edge)
			vC.edge = mesh.indexFor(fB_C_New, 1);

		//Copy across constraint state to the new faces
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(fA_B_C, 0)))
			fC_A_New.flags |= (1 << 2);
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(fA_B_C, 1)))
			fA_B_New.flags |= (1 << 2);
		if (mesh.isHalfEdgeConstrained(mesh.indexFor(fA_B_C, 2)))
			fB_C_New.flags |= (1 << 2);

		//Free the old face
		mesh.recycleFace(mesh.indexFor(fA_B_C));

		//Tag the modified edges for delaunay condition checking
		edgesToCheck.Add(mesh.indexFor(fA_B_New, 1));
		edgesToCheck.Add(mesh.indexFor(fB_C_New, 1));
		edgesToCheck.Add(mesh.indexFor(fC_A_New, 1));
		return mesh.indexFor(vNew);
	}


	static bool IsDelaunay(Mesh & mesh,Index h)
	{
		Face & face = mesh.faces[h / 4];
		Vertex & vA = mesh.vertices[face.edges[0].destinationVertex];
		Vertex & vB = mesh.vertices[face.edges[1].destinationVertex];
		Vertex & vC = mesh.vertices[face.edges[2].destinationVertex];

		Vertex & vOpp = mesh.vertices[mesh.edgeAt(Face::nextHalfEdge(mesh.edgeAt(h).oppositeHalfEdge)).destinationVertex];


		glm::dvec2 circumcenter = Geo2D::ComputeCircumcenter(vA.position, vB.position, vC.position);
		double squaredRadius = Geo2D::DistanceSquared(vA.position - circumcenter);
		double squaredDistance = Geo2D::DistanceSquared(vOpp.position - circumcenter);
		return squaredDistance >= squaredRadius;
	}
	//Expects bound to be a CW list of edges surrounding the hole to be triangulated.
	//Triangulate handles both closed and open edge contours
	//Open contours occur when triangulating the first side of an edge pair.
	static Index triangulate(Mesh & mesh, Oryol::Array<Index> & bound, bool real) {

	}
 */

};
//TODO: Investigate a better approach for initialising everything.
void Delaunay::Mesh::Setup(double width, double height)
{
    enum vIndices : Index {
        vInfinite, vBottomLeft, vBottomRight, vTopRight, vTopLeft
    };
    enum eIndices : Index {
        eBR_TL = 1, eTL_BL, eBL_BR,
        eTL_BR = 5, eBR_TR, eTR_TL,
        eInf_TL = 9, eTL_TR, eTR_Inf,
        eInf_BL = 13, eBL_TL, eTL_Inf,
        eInf_BR = 17, eBR_BL, eBL_Inf,
        eInf_TR = 21, eTR_BR, eBR_Inf
    };
	//Clear everything
	vertices.Clear();
	faces.Clear();
    constraints.Clear();

    vertices.Add({ {width * 0.5f, height * 0.5}, eTR_Inf,0,0 });
    vertices.Add({ {0,0}, eTL_BL, 2,2 });
    vertices.Add({ {width, 0}, eBL_BR, 2, 2 });
    vertices.Add({ {width, height}, eBR_TR, 2, 2 });
    vertices.Add({ {0, height}, eTR_TL, 2, 2 });
    
    faces.Add({ 0b1100, 0, {{vTopLeft,eTL_BR}, {vBottomLeft,eBL_TL}, {vBottomRight,eBR_BL}}}); //fTL_BL_BR = eBR_TL & eTL_BL & eBL_BR [0]
    faces.Add({ 0b1100, 0, {{vBottomRight,eBR_TL},{vTopRight,eTR_BR},{ vTopLeft,eTL_TR}}}); //fBR_TR_TL = eTL_BR & eBR_TR & eTR_TL [1]
    faces.Add({ 0b0101, 0, {{vTopLeft, eTL_Inf},{vTopRight,eTR_TL}, {vInfinite, eInf_TR}}}); //fTL_TR_vInf = eInf_TL & eTL_TR & eTR_Inf [2]
    faces.Add({ 0b0101, 0, {{vBottomLeft, eBL_Inf},{vTopLeft, eTL_BL},{vInfinite, eInf_TL}}}); //fBL_TL_vInf = eInf_BL & eBL_TL & eTL_Inf [3]
    faces.Add({ 0b0101, 0, {{vBottomRight, eBR_Inf},{vBottomLeft,eBL_BR},{vInfinite, eInf_BL}}}); //fBR_BL_vInf = eInf_BR & eBR_BL & eBL_Inf [4]
    faces.Add({ 0b0101, 0, {{vTopRight, eTR_Inf},{vBottomRight, eBR_TR},{vInfinite, eInf_BR}}}); //fTR_BR_vInf = eInf_TR & eTR_BR & eBR_Inf [5]

    //Add edge constraints
    //BottomLeft To BottomRight
    {
        auto & cBL_BR = constraints.Add();
        cBL_BR.startVertex = vBottomLeft;
        cBL_BR.endVertex = vBottomRight;
        cBL_BR.edgePairs.Add(0);
        halfEdgeToEdgePairMapping.Add(eBL_BR,0);
        halfEdgeToEdgePairMapping.Add(eBR_BL,0);
        edgePairToConstraintMapping.Add().Add(0);
    }
    //BottomRight To TopRight
    {
        auto & cBR_TR = constraints.Add();
        cBR_TR.startVertex = vBottomRight;
        cBR_TR.endVertex = vTopRight;
        cBR_TR.edgePairs.Add(1);
        halfEdgeToEdgePairMapping.Add(eBR_TR,1);
        halfEdgeToEdgePairMapping.Add(eTR_BR,1);
        edgePairToConstraintMapping.Add().Add(1);
    }
    //TopRight to TopLeft
    {
        auto & cTR_TL = constraints.Add();
        cTR_TL.startVertex = vTopRight;
        cTR_TL.endVertex = vTopLeft;
        cTR_TL.edgePairs.Add(2);
        halfEdgeToEdgePairMapping.Add(eTR_TL,2);
        halfEdgeToEdgePairMapping.Add(eTL_TR,2);
        edgePairToConstraintMapping.Add().Add(2);
    }
    //TopLeft to BottomLeft
    {
        auto & cTL_BL = constraints.Add();
        cTL_BL.startVertex = vTopLeft;
        cTL_BL.endVertex = vBottomLeft;
        cTL_BL.edgePairs.Add(3);
        halfEdgeToEdgePairMapping.Add(eTL_BL,3);
        halfEdgeToEdgePairMapping.Add(eBL_TL,3);
        edgePairToConstraintMapping.Add().Add(3);
    }
}

/*
Index Delaunay::Mesh::InsertVertex(double x, double y)
{
	Index centerVertex;
	Oryol::Array<Index> edgesToCheck;
	vertices.Reserve(1);
    HalfEdge::Index vertex = HalfEdge::Invalid;
	//Make sure the vertex is inside the bounding box the mesh was initialised with.
	//Locate the primitive the vertex falls on
	LocateResult result = this->Locate(x, y);
	switch (result.type) {
	case LocateResult::Vertex:
		vertex = result.object;
		break;
	case LocateResult::Edge:
		vertex = Impl::SplitEdge(*this,result.object, x, y, centerVertex, edgesToCheck);
		break;
	case LocateResult::Face:
		vertex = Impl::SplitFace(*this,result.object, x, y, edgesToCheck);
		break;
	}
	//Restore delaunay condition
	//TODO: In rare cases this code can fail and cause a crash because edge references are not stable and may occasionally hit an invalid half-edge index
	while (!edgesToCheck.Empty()) {
		Index h = edgesToCheck.PopFront();
		
		if (!isHalfEdgeConstrained(h) && !Impl::IsDelaunay(*this,h)) {
			h = Impl::FlipEdge(*this,h);
			HalfEdge & current = edgeAt(h);
			HalfEdge & opposite = edgeAt(current.oppositeHalfEdge);
			//If the newly flipped edge is now pointing towards the center, it was not before
			if (opposite.destinationVertex == centerVertex || current.destinationVertex == centerVertex) { 
				edgesToCheck.Add(Face::prevHalfEdge(indexFor(opposite)));
				edgesToCheck.Add(Face::nextHalfEdge(indexFor(current)));
			}
			else {
				edgesToCheck.Add(Face::prevHalfEdge(h));
				edgesToCheck.Add(Face::nextHalfEdge(h));
			}
			
		}
	}
	

	return vertex;
}
*/
/*
bool Delaunay::Mesh::DeleteVertex(Index index)
{
	// Few typical cases
	//If vertex.constraintCount == 0 it is completely free of constraint --> Remove immediately
	//If vertex.constraintCount > 0 && vertex.endPointCount != 0 it is the end point of a constraint --> Cannot be removed until the constraint is broken first
	//If vertex.constraintCount == 2 && vertex.endPointCount == 0 it is the midpoint of a constraint and is no longer needed
	Vertex & vertex = this->vertices[index];
	Oryol::Array<Index> bound; //Outer bound is required by untriangulate
	Oryol::Array<Index> taggedFaces; //Contains the face indices that will be removed from the mesh
	if (vertex.constraintCount == 0 && vertex.endPointCount == 0) {
		for (Index h : vertex.OutgoingEdges(*this)) {
			bound.Add(edgeAt(Face::nextHalfEdge(h)).oppositeHalfEdge);
			taggedFaces.Add(h / 4);
		}
	}
	else if (vertex.constraintCount == 2 && vertex.endPointCount == 0) {

	}

	return false;
}

*/

Delaunay::Mesh::LocateResult Delaunay::Mesh::Locate(double x, double y)
{
	LocateResult result{ HalfEdge::InvalidIndex, LocateResult::None };
	Index bestVertex = HalfEdge::InvalidIndex;
	{
		//Seed the random generator
		srand(x * 10 + y * 4);
		//Find closest vertex by randomly sampling the vertices (excluding the infinite vertex)
		int vertexSampleCount = std::pow(this->vertices.Size(), 1 / 3.);
		double minDistanceSquared = std::numeric_limits<double>::infinity();
		for (int i = 0; i < vertexSampleCount; i++) {
			Index index = rand_range(1, this->vertices.Size() - 1);
			Vertex & vertex = this->vertices[index];
			double distanceSquared = Geo2D::DistanceSquared(glm::dvec2(x,y)-vertex.position);
			if (distanceSquared < minDistanceSquared) {
				minDistanceSquared = distanceSquared;
				bestVertex = index;
			}
		}
	}
	//Start jump-and-walk search with the first face associated with the best vertex
    Index currentFace = -1;
    for(Index h : this->vertices[bestVertex].OutgoingEdges(*this)){
        if(isHalfEdgeReal(h)){
            currentFace = h / 4;
            break;
        }
    }
	Oryol::Set<Index> visitedFaces;
	int iterations = 0;
	while (visitedFaces.Contains(currentFace) || !(result = Impl::IsInFace(*this,x, y, this->faces[currentFace]))) {
		visitedFaces.Add(currentFace);
		iterations++;
		if (iterations == 50) {
			//Log this as it is taking longer than expected
		}
		else if (iterations > 1000) {
			//Bail out if too many iterations have elapsed
			result.type = LocateResult::None;
			break;
		}
		Index nextFace = HalfEdge::InvalidIndex;
		//Find the best direction to look in.
		for (HalfEdge & h : this->faces[currentFace].edges) {
			//Determine if the position falls to the right of the current half edge (thus outside of the current face)
			Vertex & originVertex = this->vertices[edgeAt(h.oppositeHalfEdge).destinationVertex];
			Vertex & destinationVertex = this->vertices[h.destinationVertex];
			if (Geo2D::Sign(originVertex.position, destinationVertex.position, { x,y }) < 0) {
				nextFace = h.oppositeHalfEdge / 4;
				break;
			}
		}
		if (nextFace != HalfEdge::InvalidIndex) {
			currentFace = nextFace;
		}
		else {
			result.type = LocateResult::None;
			break; //Something has gone wrong so log it and bail
		}
	}
	
	return result;
}

void Delaunay::Mesh::SetDebugDraw(DebugDraw * debug)
{
	this->debugDraw = debug;
}

void Delaunay::Mesh::DrawDebugData()
{
	if (debugDraw == nullptr) return;
	Oryol::Set<Index> visitedFaces;
    for (int i = 1; i < vertices.Size(); i++) {
		Vertex & vertex = this->vertices[i];
		debugDraw->DrawVertex(vertex.position);
		const Index first = vertex.edge;
		Index h = first;
		do {
			Face & face = this->faces[h / 4];
			HalfEdge & current = edgeAt(h);
			//Render Edges
			if ((h < current.oppositeHalfEdge) && (current.destinationVertex != 0) && (edgeAt(current.oppositeHalfEdge).destinationVertex != 0)) {
				Vertex & origin = this->vertices[current.destinationVertex];
				Vertex & destination = this->vertices[edgeAt(current.oppositeHalfEdge).destinationVertex];
				debugDraw->DrawEdge(origin.position, destination.position, isHalfEdgeConstrained(h));
			}
			//TODO: Render Faces
			if(!visitedFaces.Contains(h/4)){
				visitedFaces.Add(h / 4);
			}
			h = edgeAt(Face::nextHalfEdge(h)).oppositeHalfEdge;
		} while (h != first);
		
	}
}




inline Delaunay::Mesh::HalfEdge & Delaunay::Mesh::edgeAt(Index index) {
	return *(reinterpret_cast<HalfEdge*>(&faces[0])+ index);
}
inline const Delaunay::Mesh::HalfEdge & Delaunay::Mesh::edgeAt(Index index) const {
	return *(reinterpret_cast<const HalfEdge*>(&faces[0]) + index);
}

inline Index Delaunay::Mesh::Face::nextHalfEdge(Index h) {
	++h;
	return (h & 3) != 0 ? h : h - 3;
}

inline Index Delaunay::Mesh::Face::prevHalfEdge(Index h) {
	--h;
	return (h & 3) != 0 ? h : h + 3;
}


Index Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::operator*()
{
	return current;
}

void Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::operator++()
{	
	current = mesh.edgeAt(Mesh::Face::nextHalfEdge(current)).oppositeHalfEdge;
	if (current == first)
		current = Mesh::HalfEdge::InvalidIndex;
}

bool Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::operator!=(const IncomingHalfEdgeIterator & rhs)
{
	return current != Mesh::HalfEdge::InvalidIndex;
}

Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator & Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::begin()
{
	return *this;
}
Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator & Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::end()
{
	return *this;
}

Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::IncomingHalfEdgeIterator(Mesh & mesh, Index first)
	: mesh(mesh), current(first), first(first)
{
	/*
	current = mesh.vertices[vertex].edge;
	if (mesh.edgeAt(current).destinationVertex != vertex)
	current = mesh.edgeAt(current).oppositeHalfEdge;
	first = current;
	*/
	
}

Index Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::operator*()
{
	return current;
}

Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::OutgoingHalfEdgeIterator(Mesh & mesh, Index first)
	: mesh(mesh), current(first), first(first)
{
	/*
	current = mesh.vertices[vertex].edge;
	if (mesh.edgeAt(current).destinationVertex == vertex)
	current = mesh.edgeAt(current).oppositeHalfEdge;
	first = current;
	*/
}
void Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::operator++() {
	current = Mesh::Face::nextHalfEdge(mesh.edgeAt(current).oppositeHalfEdge);
	if (current == first)
		current = Mesh::HalfEdge::InvalidIndex;
}
bool Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::operator!=(const OutgoingHalfEdgeIterator & rhs)
{
	return current != Mesh::HalfEdge::InvalidIndex;
}
Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator & Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::begin()
{
	return *this;
}
Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator & Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::end()
{
	return *this;
}

Mesh::HalfEdge::IncomingHalfEdgeIterator Mesh::Vertex::IncomingEdges(Mesh & mesh) {
	return HalfEdge::IncomingHalfEdgeIterator(mesh, edge);
}
Mesh::HalfEdge::OutgoingHalfEdgeIterator Mesh::Vertex::OutgoingEdges(Mesh & mesh) {
	return HalfEdge::OutgoingHalfEdgeIterator(mesh, edge);
}
