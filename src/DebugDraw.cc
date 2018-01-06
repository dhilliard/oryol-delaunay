#include "DebugDraw.h"
#include "shaders.h"

using namespace Oryol;

uint32_t CompactColor(const Color & color) {
	union {
		struct { uint8_t r, g, b, a; };
		uint32_t value;
	} conv = { (uint8_t)(color.r * 255), (uint8_t)(color.g * 255), (uint8_t)(color.b * 255), (uint8_t)(color.a * 255) };
	return conv.value;
}
Oryol::Id BuildPipeline(const Oryol::GfxSetup & setup,Oryol::Id shd,const Oryol::VertexLayout & layout,Oryol::PrimitiveType::Code prim) {
	auto pipSetup = PipelineSetup::FromLayoutAndShader(layout, shd);
	pipSetup.RasterizerState.SampleCount = setup.SampleCount;
	pipSetup.BlendState.ColorFormat = setup.ColorFormat;
	pipSetup.BlendState.DepthFormat = setup.DepthFormat;
	pipSetup.PrimType = prim;
	pipSetup.BlendState.BlendEnabled = true;
	pipSetup.BlendState.SrcFactorRGB = BlendFactor::SrcAlpha;
	pipSetup.BlendState.DstFactorRGB = BlendFactor::OneMinusSrcAlpha;
	return Gfx::CreateResource(pipSetup);
}
void DebugDraw::Setup(const Oryol::GfxSetup & setup)
{
	Gfx::PushResourceLabel();
	{
		auto meshSetup = MeshSetup::Empty(MaxNumTriangleVertices, Usage::Stream);
		meshSetup.Layout = {
			{ VertexAttr::Position, VertexFormat::Float2 },
			{ VertexAttr::Color0, VertexFormat::UByte4N }
		};
		this->triangleDrawState.Mesh[0] = Gfx::CreateResource(meshSetup);
		//Setup up the shader and pipeline for rendering
		Id shd = Gfx::CreateResource(DebugGeometryShader::Setup());
		this->triangleDrawState.Pipeline = BuildPipeline(setup,shd,meshSetup.Layout,PrimitiveType::Triangles);
	}
	{
		auto meshSetup = MeshSetup::Empty(MaxNumLineVertices, Usage::Stream);
		meshSetup.Layout = {
			{ VertexAttr::Position, VertexFormat::Float2 },
			{ VertexAttr::Color0, VertexFormat::UByte4N }
		};
		this->lineDrawState.Mesh[0] = Gfx::CreateResource(meshSetup);
		//Setup up the shader and pipeline for rendering
		Id shd = Gfx::CreateResource(DebugGeometryShader::Setup());

		this->lineDrawState.Pipeline = BuildPipeline(setup, shd, meshSetup.Layout, PrimitiveType::Lines);
	}
	{
		auto meshSetup = MeshSetup::Empty(MaxNumPointVertices, Usage::Stream);
		meshSetup.Layout = {
			{ VertexAttr::Position, VertexFormat::Float3 },
			{ VertexAttr::Color0, VertexFormat::UByte4N }
		};
		this->pointDrawState.Mesh[0] = Gfx::CreateResource(meshSetup);
		//Setup up the shader and pipeline for rendering
		Id shd = Gfx::CreateResource(DebugPointShader::Setup());

		this->pointDrawState.Pipeline = BuildPipeline(setup, shd, meshSetup.Layout, PrimitiveType::Points);
	}
	this->resourceLabel = Gfx::PopResourceLabel();
}
void DebugDraw::Discard() {
	Gfx::DestroyResources(this->resourceLabel);
	this->resourceLabel.Invalidate();
}

void DebugDraw::Triangle(float x1, float y1, float x2, float y2, float x3, float y3, const Color & color)
{
	if (this->triangles.Size() + 3 < MaxNumTriangleVertices) {
		uint32_t uColor = CompactColor(color);
		this->triangles.Add({ x1,y1,uColor });
		this->triangles.Add({ x2,y2,uColor });
		this->triangles.Add({ x3,y3,uColor });
	}
}

void DebugDraw::Line(float x1, float y1, float x2, float y2, const Color & color)
{
	if (this->lines.Size() + 2 < MaxNumLineVertices) {
		uint32_t uColor = CompactColor(color);
		this->lines.Add({ x1,y1,uColor });
		this->lines.Add({ x2,y2,uColor });
	}
}

void DebugDraw::Point(float x, float y, float size, const Color & color)
{
	if (this->points.Size() + 1 < MaxNumPointVertices) {
		this->points.Add({ x,y,size,CompactColor(color) });
	}
}

void DebugDraw::Draw(glm::mat4x4 projectionMatrix)
{
	DebugGeometryShader::vsParams params{ projectionMatrix };
	if(!triangles.Empty()){
		Gfx::UpdateVertices(this->triangleDrawState.Mesh[0], this->triangles.begin(), this->triangles.Size() * sizeof(vertex_t));
		Gfx::ApplyDrawState(this->triangleDrawState);
		Gfx::ApplyUniformBlock(params);
		Gfx::Draw({ 0,this->triangles.Size() });
		this->triangles.Clear();
	}
	if(!lines.Empty()) {
		Gfx::UpdateVertices(this->lineDrawState.Mesh[0], this->lines.begin(), this->lines.Size() * sizeof(vertex_t));
		Gfx::ApplyDrawState(this->lineDrawState);
		Gfx::ApplyUniformBlock(params);
		Gfx::Draw({ 0,this->lines.Size() });
		this->lines.Clear();
	}
	if(!points.Empty()) {
		Gfx::UpdateVertices(this->pointDrawState.Mesh[0], this->points.begin(), this->points.Size() * sizeof(point_t));
		Gfx::ApplyDrawState(this->pointDrawState);
		Gfx::ApplyUniformBlock(params);
		Gfx::Draw({ 0,this->points.Size() });
		this->points.Clear();
	}
}
