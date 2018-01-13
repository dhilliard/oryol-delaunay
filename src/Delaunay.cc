//------------------------------------------------------------------------------
//  Triangle.cc
//------------------------------------------------------------------------------
#include "Pre.h"
#include "Core/Main.h"
#include "Gfx/Gfx.h"
#include "DebugDraw.h"
#include "glm/gtc/matrix_transform.hpp"
#include "Mesh.h"

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
    Gfx::Setup(GfxSetup::Window(400, 400, "Oryol Delaunay Sample"));
    
	debug.Setup(GfxSetup());
	mesh.Setup(400, 400);
	
	projectionMatrix = glm::ortho<double>(0, 400, 400, 0, -10, 10);
    return App::OnInit();
}

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnRunning() {
    
    Gfx::BeginPass();
	
	mesh.SetDebugDraw(&debug);
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
    Gfx::Discard();
    return App::OnCleanup();
}
