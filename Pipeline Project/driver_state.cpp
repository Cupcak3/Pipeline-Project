#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
	delete [] image_color;
	delete [] image_depth;
	image_color = NULL;
	image_depth = NULL;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
	state.image_width=width;
	state.image_height=height;
	state.image_color=0;
	state.image_depth=0;
	// Allocate memory for image_color. Initialize all pixels to black. Ignore depth_image (Used in z-buffer later)
	state.image_color = new pixel[width*height];
	state.image_depth = new float[width*height];
	
	for (int i = 0; i < state.image_height * state.image_width; ++i)
	{
		//0-255 with 0,0,0 being black and 255,255,255 being white
		state.image_color[i] = make_pixel(0, 0, 0); // Initialize all pixels to black
		state.image_depth[i] = 2; // Each pixel's default depth is 2 > 1 so each will always be set when turning on pixels
	}
}


static data_geometry perspective_divide(const data_geometry &g)
{
	float w = g.gl_Position[3];
	data_geometry g_return = g;
	
	g_return.gl_Position[0] = g.gl_Position[0] / w;
	g_return.gl_Position[1] = g.gl_Position[1] / w;
	g_return.gl_Position[2] = g.gl_Position[2] / w;
	//g_return.gl_Position[3] =g.gl_Position[3] / w;
	
	
	
	return g_return;
}

static void translate_to_pixel_space(data_geometry &g, driver_state &state)
{
	//Calculate pixel coordinates
	// 		(x,y) in 2D NDC -------> (i,j) in pixel space
	// 		i = (w/2)x + (w/2 - .5)
	// 		j = (h/2)y + (h/2 - .5)
	g.gl_Position[0] = abs(((state.image_width/2)  * g.gl_Position[0]) + ((state.image_width/2)  - .5));
	g.gl_Position[1] = abs(((state.image_height/2) * g.gl_Position[1]) + ((state.image_height/2) - .5));
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
	switch (type)
	{
		case render_type::triangle:
		{
			//std::cout<<"Type = triangle"<<std::endl;
			//Read every 3 verticies into a data_geometry array, call rasterize_triangle
			const data_geometry* g_array[3];
			//Go over all vertices
			for (int i = 0; i < state.num_vertices*state.floats_per_vertex; i+=3*state.floats_per_vertex)
			{
				//Process 3 verticies at a time
				data_geometry g1, g2, g3;
				data_vertex v1, v2, v3;
				
				v1.data = &state.vertex_data[i];
				v2.data = &state.vertex_data[i+state.floats_per_vertex];
				v3.data = &state.vertex_data[i+2*state.floats_per_vertex];
				
				g1.data = v1.data;
				g2.data = v2.data;
				g3.data = v3.data;
				
				state.vertex_shader(v1,g1,state.uniform_data);
				state.vertex_shader(v2,g2,state.uniform_data);
				state.vertex_shader(v3,g3,state.uniform_data);
				
				
				g_array[0] = &g1;
				g_array[1] = &g2;
				g_array[2] = &g3;
				
				//Render each triangle
				//clip_triangle(state, g_array, 6);
				clip_triangle(state, g_array, 0);
			}
			g_array[0] = 0;
			g_array[1] = 0;
			g_array[2] = 0;
			
			break;
		}
		case render_type::indexed:
		{
			const data_geometry* g_array[3];
			for (int i = 0; i < 3*state.num_triangles; i+=3)
			{
				//Process 3 verticies at a time
				data_geometry g1, g2, g3;
				data_vertex v1, v2, v3;
				
				v1.data = &state.vertex_data[state.index_data[i]  *state.floats_per_vertex];
				v2.data = &state.vertex_data[state.index_data[i+1]*state.floats_per_vertex];
				v3.data = &state.vertex_data[state.index_data[i+2]*state.floats_per_vertex];
				
				g1.data = v1.data;
				g2.data = v2.data;
				g3.data = v3.data;
				
				state.vertex_shader(v1,g1,state.uniform_data);
				state.vertex_shader(v2,g2,state.uniform_data);
				state.vertex_shader(v3,g3,state.uniform_data);
				
				
				g_array[0] = &g1;
				g_array[1] = &g2;
				g_array[2] = &g3;
				
				//Render each triangle
				//clip_triangle(state, g_array, 6);
				clip_triangle(state, g_array, 0);
			}
			
			g_array[0] = 0;
			g_array[1] = 0;
			g_array[2] = 0;
			
			break;
		}
		case render_type::fan:
		{
			const data_geometry* g_array[3];
			data_geometry base;
			data_vertex base_data;
			
			base_data.data = &state.vertex_data[0];
			base.data = base_data.data;
			state.vertex_shader(base_data, base, state.uniform_data);
			
			g_array[0] = &base;
			
			//translate_to_pixel_space(base, state);
			
			for (int i = state.floats_per_vertex; i < state.num_vertices*state.floats_per_vertex-state.floats_per_vertex; i+=state.floats_per_vertex)
			{
				//Process 3 verticies at a time
				data_geometry g2, g3;
				data_vertex   v2, v3;
				
				v2.data = &state.vertex_data[i];
				v3.data = &state.vertex_data[i+state.floats_per_vertex];
				
				g2.data = v2.data;
				g3.data = v3.data;
				
				state.vertex_shader(v2,g2,state.uniform_data);
				state.vertex_shader(v3,g3,state.uniform_data);
				
				
				g_array[1] = &g2;
				g_array[2] = &g3;
				
				//Render each triangle
				//clip_triangle(state, g_array, 6);
				clip_triangle(state, g_array, 0);
			}
			g_array[0] = 0;
			g_array[1] = 0;
			g_array[2] = 0;
			
			break;
		}
		case render_type::strip:
		{
			const data_geometry* g_array[3];
			int flip = 0;
			for (int i = 0; i < state.floats_per_vertex*state.num_vertices-2*state.floats_per_vertex; i+=state.floats_per_vertex)
			{
				data_geometry g1, g2, g3;
				data_vertex v1, v2, v3;
				
				v1.data = &state.vertex_data[i];
				v2.data = &state.vertex_data[i+state.floats_per_vertex];
				v3.data = &state.vertex_data[i+2*state.floats_per_vertex];
				
				if (flip % 2 == 1)
				{
					data_vertex temp;
					temp = v1;
					v1 = v2;
					v2 = temp;
				}
				
				g1.data = v1.data;
				g2.data = v2.data;
				g3.data = v3.data;
				
				state.vertex_shader(v1,g1,state.uniform_data);
				state.vertex_shader(v2,g2,state.uniform_data);
				state.vertex_shader(v3,g3,state.uniform_data);
				
				
				g_array[0] = &g1;
				g_array[1] = &g2;
				g_array[2] = &g3;
				
				//Render each triangle
				//clip_triangle(state, g_array, 6);
				clip_triangle(state, g_array, 0);
				++flip;
			}
			
			g_array[0] = 0;
			g_array[1] = 0;
			g_array[2] = 0;
			
			break;
		}
		default:
			break;
	}
	//std::cout<<"TODO: implement rendering."<<std::endl;
}


static void interpolation_helper(vec4 &point_1, float alpha_1, float alpha_2, vec4 &point_2, vec4 &point_3, data_geometry *data_array, const data_geometry **in, driver_state &state, int tri_case) {
	for (int i = 0; i < 3; ++i)
	{
		data_array[i].data = new float[MAX_FLOATS_PER_VERTEX];
		for (int j = 0; j < MAX_FLOATS_PER_VERTEX; ++j)
		{
			data_array[i].data[j] = 0; //Zero everything out
		}
	}
	for(int a = 0; a < 3; ++a)
	{
		for (int i = 0; i < state.floats_per_vertex; ++i)
		{
			switch (state.interp_rules[a*state.floats_per_vertex+i])
			{
				case interp_type::smooth:
				{
					switch (tri_case)
					{
						case 1: //A AB AC
						{
							float k_1 = 1.0 / (alpha_1 * point_1[3] + (1-alpha_1) * point_2[3]);
							float k_2 = 1.0 / (alpha_2 * point_1[3] + (1-alpha_2) * point_3[3]);
							float alpha_1_corrected = alpha_1 * point_1[3] * k_1;
							float alpha_2_corrected = alpha_2 * point_1[3] * k_2;
							
							data_array[0].data = in[0]->data;
							data_array[1].data[i] = alpha_1_corrected * in[0]->data[i] + (1-alpha_1_corrected) * in[1]->data[i];
							data_array[2].data[i] = alpha_2_corrected * in[0]->data[i] + (1-alpha_2_corrected) * in[2]->data[i];
							continue;
						}
						case 2:
						{
							float k_1 = 1.0 / (alpha_1 * point_2[3] + (1-alpha_1) * point_3[3]);
							float k_2 = 1.0 / (alpha_2 * point_1[3] + (1-alpha_2) * point_2[3]);
							float alpha_1_corrected = alpha_1 * point_2[3] * k_1;
							float alpha_2_corrected = alpha_2 * point_1[3] * k_2;
							
							data_array[0].data = in[1]->data;
							data_array[1].data[i] = alpha_1_corrected * in[1]->data[i] + (1-alpha_1_corrected) * in[2]->data[i];
							data_array[2].data[i] = alpha_2_corrected * in[0]->data[i] + (1-alpha_2_corrected) * in[1]->data[i];
							continue;
						}
						case 3:
						{
							float k_1 = 1.0 / (alpha_1 * point_3[3] + (1-alpha_1) * point_1[3]);
							float k_2 = 1.0 / (alpha_2 * point_3[3] + (1-alpha_2) * point_2[3]);
							float alpha_1_corrected = alpha_1 * point_3[3] * k_1;
							float alpha_2_corrected = alpha_2 * point_3[3] * k_2;
							
							data_array[0].data = in[0]->data;
							data_array[1].data[i] = alpha_1_corrected * in[2]->data[i] + (1-alpha_1_corrected) * in[0]->data[i];
							data_array[2].data[i] = alpha_2_corrected * in[2]->data[i] + (1-alpha_2_corrected) * in[1]->data[i];
							continue;
						}
						case 4:
						{
							float k_2 = 1.0 / (alpha_2 * point_2[3] + (1-alpha_2) * point_3[3]);
							float alpha_2_corrected = alpha_2 * point_2[3] * k_2;
							
							data_array[0].data = in[0]->data;
							data_array[1].data = in[1]->data;
							data_array[2].data[i] = alpha_2_corrected * in[1]->data[i] + (1-alpha_2_corrected) * in[2]->data[i];
							continue;
						}
						case 5:
						{
							float k_1 = 1.0 / (alpha_1 * point_2[3] + (1-alpha_1) * point_3[3]);
							float k_2 = 1.0 / (alpha_2 * point_1[3] + (1-alpha_2) * point_3[3]);
							float alpha_1_corrected = alpha_1 * point_2[3] * k_1;
							float alpha_2_corrected = alpha_2 * point_1[3] * k_2;
							
							data_array[0].data = in[0]->data;
							data_array[1].data[i] = alpha_1_corrected * in[1]->data[i] + (1-alpha_1_corrected) * in[2]->data[i];
							data_array[2].data[i] = alpha_2_corrected * in[0]->data[i] + (1-alpha_2_corrected) * in[2]->data[i];
							continue;
						}
						case 6:
						{
							float k_2 = 1.0 / (alpha_2 * point_1[3] + (1-alpha_2) * point_2[3]); // FIX POINTS
							float alpha_2_corrected = alpha_2 * point_1[3] * k_2; // FIX POINTS
							
							data_array[0].data = in[2]->data;
							data_array[1].data = in[0]->data;
							data_array[2].data[i] = alpha_2_corrected * in[0]->data[i] + (1-alpha_2_corrected) * in[1]->data[i];
							continue;
						}
						case 7:
						{
							float k_1 = 1.0 / (alpha_1 * point_1[3] + (1-alpha_1) * point_2[3]);
							float k_2 = 1.0 / (alpha_2 * point_2[3] + (1-alpha_2) * point_3[3]);
							float alpha_1_corrected = alpha_1 * point_1[3] * k_1;
							float alpha_2_corrected = alpha_2 * point_2[3] * k_2;
							
							data_array[0].data = in[0]->data;
							data_array[1].data[i] = alpha_1_corrected * in[0]->data[i] + (1-alpha_1_corrected) * in[1]->data[i];
							data_array[2].data[i] = alpha_2_corrected * in[1]->data[i] + (1-alpha_2_corrected) * in[2]->data[i];
							continue;
						}
						case 8:
						{
							float k_2 = 1.0 / (alpha_2 * point_1[3] + (1-alpha_2) * point_3[3]);
							float alpha_2_corrected = alpha_2 * point_1[3] * k_2;
							
							data_array[0].data = in[1]->data;
							data_array[1].data = in[2]->data;
							data_array[2].data[i] = alpha_2_corrected * in[0]->data[i] + (1-alpha_2_corrected) * in[2]->data[i];
							continue;
						}
						case 9:
						{
							float k_1 = 1.0 / (alpha_1 * point_1[3] + (1-alpha_1) * point_3[3]);
							float k_2 = 1.0 / (alpha_2 * point_1[3] + (1-alpha_2) * point_2[3]);
							float alpha_1_corrected = alpha_1 * point_1[3] * k_1;
							float alpha_2_corrected = alpha_2 * point_1[3] * k_2;
							
							data_array[0].data = in[1]->data;
							data_array[1].data[i] = alpha_1_corrected * in[0]->data[i] + (1-alpha_1_corrected) * in[2]->data[i];
							data_array[2].data[i] = alpha_2_corrected * in[0]->data[i] + (1-alpha_2_corrected) * in[1]->data[i];
							continue;
						}
						default:
							continue;
					}
				}
				case interp_type::noperspective:
				{
					switch (tri_case)
					{
						case 1: //A AB AC
						{
							data_array[0].data = in[0]->data;
							data_array[1].data[i] = alpha_1 * in[0]->data[i] + (1-alpha_1) * in[1]->data[i];
							data_array[2].data[i] = alpha_2 * in[0]->data[i] + (1-alpha_2) * in[2]->data[i];
							continue;
						}
						case 2: //B BC AB
						{
							data_array[0].data = in[1]->data;
							data_array[1].data[i] = alpha_1 * in[1]->data[i] + (1-alpha_1) * in[2]->data[i];
							data_array[2].data[i] = alpha_2 * in[0]->data[i] + (1-alpha_2) * in[1]->data[i];
							continue;
						}
						case 3: //C CA CB
						{
							data_array[0].data = in[2]->data;
							data_array[1].data[i] = alpha_1 * in[2]->data[i] + (1-alpha_1) * in[0]->data[i];
							data_array[2].data[i] = alpha_2 * in[2]->data[i] + (1-alpha_2) * in[1]->data[i];
							continue;
						}
						case 4: //A B BC
						{
							data_array[0].data = in[0]->data;
							data_array[1].data = in[1]->data;
							data_array[2].data[i] = alpha_2 * in[1]->data[i] + (1-alpha_2) * in[2]->data[i];
							continue;
						}
						case 5:
						{
							data_array[0].data = in[0]->data;
							data_array[1].data[i] = alpha_1 * in[1]->data[i] + (1-alpha_1) * in[2]->data[i];
							data_array[2].data[i] = alpha_2 * in[0]->data[i] + (1-alpha_2) * in[2]->data[i];
							continue;
						}
						case 6:
						{
							data_array[0].data = in[2]->data;
							data_array[1].data = in[0]->data;
							data_array[2].data[i] = alpha_2 * in[0]->data[i] + (1-alpha_2) * in[1]->data[i];
							continue;
						}
						case 7:
						{
							data_array[0].data = in[2]->data;
							data_array[1].data[i] = alpha_1 * in[0]->data[i] + (1-alpha_1) * in[1]->data[i];
							data_array[2].data[i] = alpha_2 * in[1]->data[i] + (1-alpha_2) * in[2]->data[i];
							continue;
						}
						case 8:
						{
							data_array[0].data = in[1]->data;
							data_array[1].data = in[2]->data;
							data_array[2].data[i] = alpha_2 * in[0]->data[i] + (1-alpha_2) * in[2]->data[i];
							continue;
						}
						case 9:
						{
							data_array[0].data = in[1]->data;
							data_array[1].data[i] = alpha_1 * in[0]->data[i] + (1-alpha_1) * in[2]->data[i];
							data_array[2].data[i] = alpha_2 * in[0]->data[i] + (1-alpha_2) * in[1]->data[i];
							continue;
						}
						default:
							continue;
					}
				}
				case interp_type::flat:
				{
					for (int j = 0; j < 3; ++j)
					{
						data_array[j].data = in[0]->data;
					}
					continue;
				}
				default:
					continue;
			}
		}
	}
}

static float midpoint_weight(vec4 &A, vec4 &B, int axis, int sign) {
	return (sign*B[3] - B[axis]) / (A[axis] - sign*A[3] + sign*B[3] - B[axis]);
}

static vec<float, 4> midpoint(const vec4 &first, float midpoint_weight, const vec4 &last) {
	return midpoint_weight * first + (1-midpoint_weight) * last;
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3], int face)
{
	if (face > 6) return;
	if (face == 6)
	{
		const data_geometry *tri[3];
		
		data_geometry g_a = perspective_divide(*in[0]);
		data_geometry g_b = perspective_divide(*in[1]);
		data_geometry g_c = perspective_divide(*in[2]);
		
		tri[0] = &g_a;
		tri[1] = &g_b;
		tri[2] = &g_c;
		
		translate_to_pixel_space(g_a, state);
		translate_to_pixel_space(g_b, state);
		translate_to_pixel_space(g_c, state);
		
		rasterize_triangle(state, tri);
		return;
	}
	vec4 A = in[0]->gl_Position;
	vec4 B = in[1]->gl_Position;
	vec4 C = in[2]->gl_Position;
	
	//std::cout<<"New clip face: "<<face<<std::endl;
	
	const data_geometry * in1[3];
	const data_geometry * in2[3];
	
	int inside_points[3] = {0,0,0};  // Track A,B,C inside axis. 1 true, 0 false
	
	int axis = face % 3; //X Y Z
	int sign = (face%2==0) ? 1 : -1; // Positive or negative sign
	
	if (sign == 1)
	{
		if (A[axis] < A[3])
			inside_points[0] = 1;
		if (B[axis] < B[3])
			inside_points[1] = 1;
		if (C[axis] < C[3])
			inside_points[2] = 1;
	}
	else
	{
		if (A[axis] > -A[3])
			inside_points[0] = 1;
		if (B[axis] > -B[3])
			inside_points[1] = 1;
		if (C[axis] > -C[3])
			inside_points[2] = 1;
	}
	
	for (int i = 0; i < 3; ++i)
	{
		//std::cout<< inside_points[i] << " ";
		//std::cout<<in[i]->gl_Position<<std::endl;
	}//std::cout<<std::endl;
	
	
	
	//std::cout<<"Clipping face "<<face<<std::endl;
	if (!inside_points[0] && !inside_points[1] && !inside_points[2])   return;  //No points inside
	else if (inside_points[0] && inside_points[1] && inside_points[2]) clip_triangle(state, in, face+1); // All points inside
	
	//ONE POINT INSIDE
	else if (inside_points[0] && !inside_points[1] && !inside_points[2])
	{
		//A INSIDE   BC OUTSIDE   A AB AC
		float AB_alpha, AC_alpha;
		
		AB_alpha = midpoint_weight(A, B, axis, sign); // alpha 1
		AC_alpha = midpoint_weight(A, C, axis, sign); // alpha 2
		
		vec4 AB = midpoint(A, AB_alpha, B);           // Midpoint with alpha 1
		vec4 AC = midpoint(A, AC_alpha, C);           // Midpoint with alpha 2
		
		data_geometry * data_array = new data_geometry[3];
		data_array[0].gl_Position = A;
		data_array[1].gl_Position = AB;
		data_array[2].gl_Position = AC;
		
		interpolation_helper(A, AB_alpha, AC_alpha, B, C, data_array, in, state, 1);
		
		for (int j = 0; j < 3; ++j)
			in1[j] = &data_array[j];
		clip_triangle(state, in1, face+1);
	}
	else if (!inside_points[0] && inside_points[1] && !inside_points[2])
	{
		//std::cout<<"CASE A BC AB"<<std::endl;
		//B INSIDE   AC OUTSIDE   B BC AB
		float BC_t, AB_t;
		
		BC_t = midpoint_weight(B, C, axis, sign);
		AB_t = midpoint_weight(A, B, axis, sign);
		
		vec4 BC = midpoint(B, BC_t, C);
		vec4 AB = midpoint(A, AB_t, B);
		
		data_geometry * data_array = new data_geometry[3];
		data_array[0].gl_Position = B;
		data_array[1].gl_Position = BC;
		data_array[2].gl_Position = AB;
		
		interpolation_helper(A, BC_t, AB_t, B, C, data_array, in, state, 2);
		
		for (int j = 0; j < 3; ++j)
			in1[j] = &data_array[j];
		clip_triangle(state, in1, face+1);
	}
	else if (!inside_points[0] && !inside_points[1] && inside_points[2])
	{
		//std::cout<<"CASE C CA CB"<<std::endl;
		//C INSIDE   AB OUTSIDE   C CA CB <--
		
		float CA_t, CB_t;
		
		CA_t = midpoint_weight(C, A, axis, sign); // alpha 1
		CB_t = midpoint_weight(C, B, axis, sign); // alpha 2
		
		vec4 CA = midpoint(C, CA_t, A);           // Midpoint with alpha 1
		vec4 CB = midpoint(C, CB_t, B);           // Midpoint with alpha 2
		
		data_geometry * data_array = new data_geometry[3];
		data_array[0].gl_Position = C;
		data_array[1].gl_Position = CA;
		data_array[2].gl_Position = CB;
		
		interpolation_helper(A, CA_t, CB_t, B, C, data_array, in, state, 3);
		
		for (int j = 0; j < 3; ++j)
			in1[j] = &data_array[j];
		clip_triangle(state, in1, face+1);
		
	}
	
	
	//TWO POINTS INSIDE
	else if (inside_points[0] && inside_points[1] && !inside_points[2])
	{
		//TRIANGLE         A B BC
		float BC_t;
		
		BC_t = midpoint_weight(B, C, axis, sign); // alpha 2
		
		vec4 BC = midpoint(B, BC_t, C);           // Midpoint with alpha 2
		
		data_geometry * data_array = new data_geometry[3];
		data_array[0].gl_Position = A;
		data_array[1].gl_Position = B;
		data_array[2].gl_Position = BC;
		
		interpolation_helper(A, 1, BC_t, B, C, data_array, in, state, 4);
		
		for (int j = 0; j < 3; ++j)
			in1[j] = &data_array[j];
		//clip
		clip_triangle(state, in1, face+1);
		
		//std::cout<<"CASE A BC AC (2)"<<std::endl;
		//TRIANGLE               A BC AC
		//A INSIDE   BC OUTSIDE  A AB AC
		float AC_t;
		
		AC_t = midpoint_weight(A, C, axis, sign); // alpha 2
		
		vec4 AC = midpoint(A, AC_t, C);           // Midpoint with alpha 2
		
		data_array[0].gl_Position = A;
		data_array[1].gl_Position = BC;
		data_array[2].gl_Position = AC;
		
		interpolation_helper(A, BC_t, AC_t, B, C, data_array, in, state, 5);
		
		
		for (int j = 0; j < 3; ++j)
			in2[j] = &data_array[j];
		//clip
		clip_triangle(state, in2, face+1);
	}
	else if (inside_points[0] && inside_points[2] && !inside_points[1])
	{
		//std::cout<<"CASE C A AB (2)"<<std::endl;
		//AC INSIDE   B OUTSIDE
		
		//TRIANGLE    C A AB
		float AB_t;
		
		AB_t = midpoint_weight(A, B, axis, sign);
		
		vec4 AB = midpoint(A, AB_t, B);
		
		data_geometry * data_array = new data_geometry[3];
		data_array[0].gl_Position = C;
		data_array[1].gl_Position = A;
		data_array[2].gl_Position = AB;
		
		interpolation_helper(A, 1, AB_t, B, C, data_array, in, state, 6);
		for (int j = 0; j < 3; ++j)
			in1[j] = &data_array[j];
		//clip
		clip_triangle(state, in1, face+1);
		
		//TRIANGLE    C AB BC
		float BC_t;
		
		BC_t = midpoint_weight(B, C, axis, sign);
		
		vec4 BC = midpoint(A, AB_t, B);
		
		data_array[0].gl_Position = C;
		data_array[1].gl_Position = AB;
		data_array[2].gl_Position = BC;
		
		interpolation_helper(A, AB_t, BC_t, B, C, data_array, in, state, 7);
		
		for (int j = 0; j < 3; ++j)
			in2[j] = &data_array[j];
		//clip
		clip_triangle(state, in2, face+1);
	}
	else if (inside_points[1] && inside_points[2] && !inside_points[0])
	{
		//BC INSIDE   A OUTSIDE
		
		//TRIANGLE    B C AC
		float AC_t;
		
		AC_t = midpoint_weight(A, C, axis, sign);
		
		vec4 AC = midpoint(A, AC_t, C);
		
		data_geometry * data_array = new data_geometry[3];
		data_array[0].gl_Position = B;
		data_array[1].gl_Position = C;
		data_array[2].gl_Position = AC;
		
		interpolation_helper(A, 1, AC_t, B, C, data_array, in, state, 8);
		for (int j = 0; j < 3; ++j)
			in1[j] = &data_array[j];
		//clip
		clip_triangle(state, in1, face+1);
		
		//std::cout<<"CASE B AC AB (2)"<<std::endl;
		//Triangle    B AC AB
		float AB_t;
		
		AB_t = midpoint_weight(A, B, axis, sign);
		
		vec4 AB = midpoint(A, AB_t, B);
		
		data_array[0].gl_Position = B;
		data_array[1].gl_Position = AC;
		data_array[2].gl_Position = AB;
		
		interpolation_helper(A, AC_t, AB_t, B, C, data_array, in, state, 9);
		
		for (int j = 0; j < 3; ++j)
			in2[j] = &data_array[j];
		//clip
		clip_triangle(state, in2, face+1);
	}
	
	return;
	
}

static double calc_area(double Ax, double Ay, double Bx, double By, double Cx, double Cy)
{
	return (0.5 * ((Bx*Cy - Cx*By) - (Ax*Cy - Cx*Ay) + (Ax*By - Bx*Ay)));
}


static float calc_depth(float alpha, float beta, float gamma, const data_geometry **in)
{
	//Z-buffer vec4 p' = alpha'*a' + beta'*b' + gamma'*c'
	vec4 point_on_frag = alpha * in[0]->gl_Position + beta * in[1]->gl_Position + gamma * in[2]->gl_Position;
	return point_on_frag[2];
}

static void shade_pixel(int i, int j, const data_geometry **in, driver_state &state, float alpha, float beta, float gamma)
{
	data_output output;
	float * frag_data = new float[MAX_FLOATS_PER_VERTEX];
	
	for (int a = 0; a < 3; ++a)
	{
		for (int b = 0; b < state.floats_per_vertex; ++b)
		{
			// Pixel Fragment coloring rules
			switch (state.interp_rules[a*state.floats_per_vertex+b])
			{
				case interp_type::noperspective:
				{
					frag_data[b] = alpha * in[0]->data[b] + beta * in[1]->data[b] + gamma * in[2]->data[b];
					continue;
				}
				case interp_type::flat:
				{
					frag_data[b] = in[a]->data[b];
					continue;
				}
				case interp_type::smooth: //Perspective correct
				{
					//alpha = alpha' / (Wa*k);
					//beta  = beta'  / (Wb*k);
					//gamma = gamma' / (Wc*k);
					float k = (alpha / in[0]->gl_Position[3]) + (beta / in[1]->gl_Position[3]) + (gamma / in[2]->gl_Position[3]);
					float alpha_corrected = (alpha / (in[0]->gl_Position[3] * k));
					float beta_corrected  = (beta  / (in[1]->gl_Position[3] * k));
					float gamma_corrected = (gamma / (in[2]->gl_Position[3] * k));
					frag_data[b] = alpha_corrected * in[0]->data[b] + beta_corrected * in[1]->data[b] + gamma_corrected * in[2]->data[b];
					continue;
				}
				default:
					continue;
			}
		}
	}
	const data_fragment fragment {frag_data};
	
	state.fragment_shader(fragment, output, state.uniform_data);
	output.output_color *= 255;
	
	
	float depth = calc_depth(alpha, beta, gamma, in);
	//std::cout<<i+j*state.image_width<<std::endl;
	//std::cout<<state.image_depth[i+j*state.image_width]<<std::endl;
	if(i+j*state.image_width < state.image_width*state.image_height)
		if(depth <= state.image_depth[i+j*state.image_width])
		{
			state.image_color[i+j*state.image_width] = make_pixel(output.output_color[0], output.output_color[1], output.output_color[2]);
			state.image_depth[i+j*state.image_width] = depth;
		}
	
	delete [] frag_data;
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
	/*
	 for (int k = 0; k < 3; ++k)
	 {
	 int i = in[k]->gl_Position[0];
	 int j = in[k]->gl_Position[1];
	 
	 std::cout<<"("<<i<<", "<<j<<")"<<std::endl;
	 }
	 */
	
	//Rasterize triangle
	//	Iterate over all pixels
	//		At pixel (i,j) use barycentric coordinates of pixel to determine if pixel is inside triangle or not.
	//		If inside set to white
	
	float Ax = in[0]->gl_Position[0];
	float Ay = in[0]->gl_Position[1];
	
	float Bx = in[1]->gl_Position[0];
	float By = in[1]->gl_Position[1];
	
	float Cx = in[2]->gl_Position[0];
	float Cy = in[2]->gl_Position[1];
	
	float A_Triangle_Total = (calc_area(Ax, Ay, Bx, By, Cx, Cy));
	
	float maxX = std::max(std::max(Ax, Bx), Cx);
	float minX = std::min(std::min(Ax, Bx), Cx);
	
	float maxY = std::max(std::max(Ay, By), Cy);
	float minY = std::min(std::min(Ay, By), Cy);
	
	float alpha = -1, beta = -1, gamma = -1;
	
	for (int j = minY; j <= maxY; ++j)
	{
		for (int i = minX; i <= maxX; ++i)
		{
			float tri_A = calc_area(i, j, Bx, By, Cx, Cy);
			float tri_B = calc_area(Ax, Ay, i, j, Cx, Cy);
			float tri_C = calc_area(Ax, Ay, Bx, By, i, j);
			
			alpha = tri_A / A_Triangle_Total;
			beta  = tri_B / A_Triangle_Total;
			gamma = tri_C / A_Triangle_Total;
			
			if ((0 <= alpha) && (0 <= beta) && (0 <= gamma) && (alpha <= 1) && (beta <= 1) && (gamma <= 1))
			{
				shade_pixel(i, j, in, state, alpha, beta, gamma);
			}
		}
	}
	
	//std::cout<<"TODO: implement rasterization"<<std::endl;
}
