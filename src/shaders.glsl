
@vs debugGeometryVS
uniform vsParams {
	mat4 mvp;
};
in vec2 position;
in vec4 color0;
out vec4 color;
void main() {
    gl_Position = mvp * vec4(position,0,1);
    color = color0;
}
@end

@fs debugGeometryFS
in vec4 color;
out vec4 fragColor;
void main() {
    fragColor = color;
}
@end

@vs debugPointVS
uniform vsParams {
	mat4 mvp;
};
in vec4 color0;
in vec3 position;
out vec4 color;
void main() {
    gl_Position = mvp * vec4(position.xy,0,1);
	gl_PointSize = position.z;
	color = color0;
}
@end

@program DebugPointShader debugPointVS debugGeometryFS
@program DebugGeometryShader debugGeometryVS debugGeometryFS
