#version 150 

//these are the interpolated values out of the rasterizer, so you can't know
//their specific values without knowing the vertices that contributed to them


in vec3 fs_normal;
in vec3 fs_light_vector;
in vec3 fs_color;
in vec4 vs_position;
in vec4 vs_normal;

uniform mat4 u_modelMatrix;
uniform vec3 u_lightColor;


out vec4 out_Color;

void main() {
    //base colors for materials
    vec4 diffuseColor = vec4(fs_color, 1.0);
    vec4 specularColor=vec4(u_lightColor,1.0);
    vec4 ambientColor=vec4(0,0,0.8,1.0);

    vec3 v=vec3((u_modelMatrix * vs_position).xyz);
    vec3 L = fs_light_vector;   
    vec3 E = normalize(-v); // we are in Eye Coordinates, so EyePos is (0,0,0) 
    vec3 R = normalize(-reflect(L,fs_normal));  

    //calculate diffuse term and clamp to the range [0, 1]
    float diffuseTerm = clamp(dot(normalize(fs_normal), normalize(fs_light_vector)), 0.0, 1.0);
	float specularTerm=clamp(pow(dot(R, E),4), 0.0, 1.0);

    out_Color = 0.5*diffuseColor * diffuseTerm+0.1*ambientColor+specularColor*specularTerm;
 //	out_Color = 0.9*diffuseColor * diffuseTerm+0.1*ambientColor;





//http://www.opengl.org/sdk/docs/tutorials/ClockworkCoders/lighting.php   OPENGL PHONG SHADER
}