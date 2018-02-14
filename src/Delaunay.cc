//------------------------------------------------------------------------------
//  Triangle.cc
//------------------------------------------------------------------------------
#include "Pre.h"
#include "Core/Main.h"
#include "Gfx/Gfx.h"
#include "Input/Input.h"
#include "DebugBatch.h"
#include "glm/gtc/matrix_transform.hpp"
#include "Mesh.h"
#include "Geo2D.h"
#include "Path.h"
using namespace Oryol;
using namespace Delaunay;

class MeshDraw : public DebugBatch {
public:
	void DrawVertex(glm::vec2 position) {
		this->Point(position.x, position.y, 5, { 1,1,1 });
	}
    void DrawEdge(const Mesh::Vertex & origin, const Mesh::Vertex & destination, bool constrained) {
		if (constrained)
			this->Line(origin.position.x, origin.position.y, destination.position.x, destination.position.y, { 1,0,0,0.8f });
		else
			this->Line(origin.position.x, origin.position.y, destination.position.x, destination.position.y, { 1,1,1,0.8f });
	}
    void DrawFace(const Mesh & mesh,uint32_t fIndex,Color color){
        const Mesh::Face & face = mesh.FaceAt(fIndex);
        const Mesh::Vertex & vA = mesh.VertexAt(face.edges[0].destinationVertex);
        const Mesh::Vertex & vB = mesh.VertexAt(face.edges[1].destinationVertex);
        const Mesh::Vertex & vC = mesh.VertexAt(face.edges[2].destinationVertex);
        this->Triangle(vA.position.x, vA.position.y, vB.position.x, vB.position.y, vC.position.x, vC.position.y, color);
    }
};

class DelaunayApp : public App {
public:
    AppState::Code OnRunning();
    AppState::Code OnInit();
    AppState::Code OnCleanup();
	glm::mat4 projectionMatrix;
	MeshDraw debug;
	Delaunay::Mesh mesh;
    Array<uint32_t> pathFaces;
    Array<uint32_t> pathEdges;
};
OryolMain(DelaunayApp);

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnInit() {
    // setup rendering system
    Gfx::Setup(GfxSetup::Window(800, 600, "Oryol Delaunay Sample"));
	Input::Setup();
	debug.Setup(GfxSetup());

	mesh.Setup(600, 600);
    //mesh.InsertConstraintSegment({0,100}, {175,100});
    //mesh.InsertConstraintSegment({400,100}, {225,100});
    //mesh.InsertConstraintSegment({100,300}, {300,300});
    //mesh.InsertConstraintSegment({0,200}, {300,200});
    //mesh.InsertConstraintSegment({300,200}, {500,400});
	mesh.InsertConstraintSegment({450,500}, {500,350});
	mesh.InsertConstraintSegment({ 150,500 }, { 100,350 });
	mesh.InsertConstraintSegment({150,500}, {450,500});
	mesh.InsertConstraintSegment({ 100,350 }, { 300,200 });
	auto vertex = mesh.InsertVertex({ 300,300 });
	mesh.RemoveVertex(vertex);
    //mesh.InsertVertex({400,150});
    /*
    mesh.InsertConstraintSegment({ 50,100 }, { 350,100 });
    //auto segment = mesh.InsertConstraintSegment({100,300}, {200,150});
    mesh.InsertConstraintSegment({50,300}, {350,300});
    mesh.InsertConstraintSegment( { 0,200 }, { 400,200 } );
    mesh.InsertConstraintSegment({25,250}, {375,250});
    mesh.InsertConstraintSegment({350,100}, {200,300});
    mesh.InsertConstraintSegment({200,300}, {50,100});
    
    mesh.InsertConstraintSegment({50,50}, {350,50});
    //mesh.RemoveConstraintSegment(segment);
    */
    
    //Path::FindPath(mesh, {85,65}, {100,225}, 0, pathFaces, pathEdges);
    
	projectionMatrix = glm::ortho<float>(0, 800, 600, 0, -10, 10);
    return App::OnInit();
}
static const char * names[5] = { "None","Vertex","Edge","Face" };
//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnRunning() {
    
	if (Input::MouseButtonDown(MouseButton::Left)) {
        auto pos = Input::MousePosition();
        pos.x = std::round(pos.x);
        pos.y = std::round(pos.y);
        if(Geo2D::IsInRange<double>(0, pos.x, 600) && Geo2D::IsInRange<double>(0, pos.y, 600)){
            //auto result = mesh.InsertVertex(pos);
            //Log::Info("Added vertex: %lu\n", result);
            auto result = mesh.Locate({pos.x,pos.y});
            Log::Info("Got object: %s at (%.2f,%.2f) = (id: %u)\n",names[result.type],pos.x,pos.y, result.object);
        } else {
            Log::Info("Outside bounding box at (%d,%d)\n",(int)pos.x,(int)pos.y);
        }

	}
    
    Gfx::BeginPass();
	
    for(uint32_t vIndex : mesh.ActiveVertexIndices()){
        if(vIndex != 0){
            const Delaunay::Mesh::Vertex & vertex = mesh.VertexAt(vIndex);
            debug.DrawVertex(vertex.position);
            const uint32_t first = mesh.GetOutgoingEdgeFor(vIndex);
            uint32_t h = first;
            do {
                const Delaunay::Mesh::HalfEdge & edge = mesh.EdgeAt(h);
                const Delaunay::Mesh::HalfEdge & opposite = mesh.EdgeAt(edge.oppositeHalfEdge);
                if((h < edge.oppositeHalfEdge) && edge.destinationVertex != 0 && opposite.destinationVertex != 0){
                    const Mesh::Vertex & destination = mesh.VertexAt(edge.destinationVertex);
                    debug.DrawEdge(vertex,destination,edge.constrained);
                }
            } while((h = mesh.GetNextOutgoingEdge(h)) != first);
        }
    }
    for(uint32_t fIndex : pathFaces){
        debug.DrawFace(mesh,fIndex,{0,1,0,0.3f});
    }
    
    debug.Draw(projectionMatrix);

    Gfx::EndPass();
    Gfx::CommitFrame();
    
    // continue running or quit?
    return Gfx::QuitRequested() ? AppState::Cleanup : AppState::Running;
}

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnCleanup() {
	Input::Discard();
    Gfx::Discard();
	
    return App::OnCleanup();
}
