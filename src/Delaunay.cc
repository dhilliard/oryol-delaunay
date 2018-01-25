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

struct Entity {
	double x, y;
	double width, height;
};
class DelaunayApp : public App {
public:
    AppState::Code OnRunning();
    AppState::Code OnInit();
    AppState::Code OnCleanup();
	glm::mat4 projectionMatrix;
	MeshDraw debug;
	Delaunay::Mesh mesh;
	ObjectPool<Entity> entities;
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
    

	mesh.InsertConstraintSegment({50,300}, {350,300});
	mesh.InsertConstraintSegment({ 200,350 }, { 100,300 });
	//mesh.SplitFace(l.object , 300, 300);
	projectionMatrix = glm::ortho<float>(-100, 500, -100, 500, -10, 10);
    return App::OnInit();
}
static const char * names[4] = { "None","Vertex","Edge","Face" };
//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnRunning() {
	if (Input::MouseButtonDown(MouseButton::Left)) {
		auto pos = Input::MousePosition() - glm::vec2{100, 100};
		auto result = mesh.Locate(pos);
		Log::Info("Got Result: %s (x: %.2f, y: %.2f, Index: %d, Generation: %d)\n", names[result.type], pos.x, pos.y, result.object, result.generation);
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
