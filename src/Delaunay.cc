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
#include "ObjectPool.h"
using namespace Oryol;

class MeshDraw : public Delaunay::DebugDraw, public DebugBatch {
	virtual void DrawVertex(glm::vec2 position) override {
		this->Point(position.x, position.y, 5, { 1,1,1 });
	}
	virtual void DrawEdge(glm::vec2 origin, glm::vec2 destination, bool constrained) override {
		if (constrained)
			this->Line(origin.x, origin.y, destination.x, destination.y, { 1,0,0,0.8f });
		else
			this->Line(origin.x, origin.y, destination.x, destination.y, { 1,1,1,0.8f });
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
};
OryolMain(DelaunayApp);

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnInit() {
    // setup rendering system
    Gfx::Setup(GfxSetup::Window(600, 600, "Oryol Delaunay Sample"));
	Input::Setup();
	debug.Setup(GfxSetup());

	mesh.Setup(400, 400);
	mesh.SetDebugDraw(&debug);
    
    mesh.InsertConstraintSegment({ 50,100 }, { 350,100 });
    //auto segment = mesh.InsertConstraintSegment({100,300}, {200,150});
    mesh.InsertConstraintSegment({50,300}, {350,300});
    mesh.InsertConstraintSegment( { 0,200 }, { 400,200 } );
    mesh.InsertConstraintSegment({25,250}, {375,250});
    mesh.InsertConstraintSegment({350,100}, {200,300});
    mesh.InsertConstraintSegment({200,300}, {50,100});
    
    mesh.InsertConstraintSegment({50,50}, {350,50});
    
    
	projectionMatrix = glm::ortho<float>(-100, 500, 500, -100, -10, 10);
    return App::OnInit();
}
static const char * names[5] = { "None","Vertex","Edge","Face","Segment" };
//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnRunning() {
    
	if (Input::MouseButtonDown(MouseButton::Left)) {
        auto pos = Input::MousePosition() - glm::vec2{100,100};
        if(Geo2D::IsInRange<double>(0, pos.x, 401) && Geo2D::IsInRange<double>(0, pos.y, 401)){
            //auto result = mesh.InsertVertex(pos);
            //Log::Info("Added vertex: %lu\n", result);
            auto result = mesh.Locate({(int)pos.x,(int)pos.y});
            Log::Info("Got object: %s at (%d,%d) = (id: %u, generation: %u)\n",names[result.type],(int)pos.x,(int)pos.y, result.index, result.generation);
        } else {
            Log::Info("Outside bounding box at (%d,%d)\n",(int)pos.x,(int)pos.y);
        }

	}
    
    Gfx::BeginPass();
	
	mesh.DrawDebugData();
	
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
