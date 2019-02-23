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


static void divide_coordinates(data_geometry &g)
{
	g.gl_Position[0] /= g.gl_Position[3];
	g.gl_Position[1] /= g.gl_Position[3];
}

static void translate_to_pixel_space(data_geometry &g, driver_state &state)
{
	//Calculate pixel coordinates
	// 		(x,y) in 2D NDC -------> (i,j) in pixel space
	// 		i = (w/2)x + (w/2 - .5)
	// 		j = (h/2)y + (h/2 - .5)
	g.gl_Position[0] = ((state.image_width/2)  * g.gl_Position[0]) + ((state.image_width/2)  - .5);
	g.gl_Position[1] = ((state.image_height/2) * g.gl_Position[1]) + ((state.image_height/2) - .5);
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
				
				divide_coordinates(g1);
				divide_coordinates(g2);
				divide_coordinates(g3);
			
				translate_to_pixel_space(g1, state);
				translate_to_pixel_space(g2, state);
				translate_to_pixel_space(g3, state);
				
				g_array[0] = &g1;
				g_array[1] = &g2;
				g_array[2] = &g3;
				
				//Render each triangle
				rasterize_triangle(state, g_array);
			}
			g_array[0] = 0;
			g_array[1] = 0;
			g_array[2] = 0;
			
			break;
		}
		case render_type::indexed:
		{
			break;
		}
		case render_type::fan:
		{
			break;
		}
		case render_type::strip:
		{
			break;
		}
		default:
			break;
	}
    //std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

static double calc_area(int Ax, int Ay, int Bx, int By, int Cx, int Cy)
{
	return (0.5 * ((Bx*Cy - Cx*By) - (Ax*Cy - Cx*Ay) + (Ax*By - Bx*Ay)));
}


static void shade_pixel(int i, int j, const data_geometry **in, driver_state &state, float alpha, float beta, float gamma) {
	data_output output;
	float * frag_data = new float[MAX_FLOATS_PER_VERTEX];
	
	for (int m = 0; m < 3; ++m)
	{
		for (int n = 0; n < state.floats_per_vertex; ++n)
		{
			switch (state.interp_rules[m*state.floats_per_vertex+n])
			{
				case interp_type::noperspective:
				{
					frag_data[n] = alpha * in[0]->data[n] + beta * in[1]->data[n] + gamma * in[2]->data[n];
					continue;
				}
				case interp_type::flat:
				{
					frag_data[m] = in[m]->data[n];
					continue;
				}
				case interp_type::smooth:
				{
					
				}
				default:
					continue;
			}
		}
	}
	const data_fragment fragment {frag_data};
	
	state.fragment_shader(fragment, output, state.uniform_data);
	output.output_color *= 255;
	
	
	state.image_color[i+j*state.image_width] = make_pixel(output.output_color[0], output.output_color[1], output.output_color[2]);
	
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
	
	int Ax = in[0]->gl_Position[0];
	int Ay = in[0]->gl_Position[1];
	
	int Bx = in[1]->gl_Position[0];
	int By = in[1]->gl_Position[1];
	
	int Cx = in[2]->gl_Position[0];
	int Cy = in[2]->gl_Position[1];
	
	float A_Triangle_Total = (calc_area(Ax, Ay, Bx, By, Cx, Cy));
	
	int maxX = std::max(std::max(Ax, Bx), Cx);
	int minX = std::min(std::min(Ax, Bx), Cx);
	
	int maxY = std::max(std::max(Ay, By), Cy);
	int minY = std::min(std::min(Ay, By), Cy);
	
	float alpha = -1, beta = -1, gamma = -1;
	
	// Attempted optimizations
	/*
	float a_k0 = calc_area(0, 0, Bx, By, Cx, Cy);
	float a_k1 = calc_area(1, 0, Bx, By, Cx, Cy) - a_k0;
	float a_k2 = calc_area(0, 1, Bx, By, Cx, Cy) - a_k0;
	
	float b_k0 = calc_area(Ax, Ay, 0, 0, Cx, Cy);
	float b_k1 = calc_area(Ax, Ay, 1, 0, Cx, Cy) - b_k0;
	float b_k2 = calc_area(Ax, Ay, 0, 1, Cx, Cy) - b_k0;
	
	float g_k0 = calc_area(Ax, Ay, Bx, By, 0, 0);
	float g_k1 = calc_area(Ax, Ay, Bx, By, 1, 0) - g_k0;
	float g_k2 = calc_area(Ax, Ay, Bx, By, 0, 1) - g_k0;
	
	alpha = a_k0 + a_k1*minX + a_k2*minY;
	beta  = b_k0 + b_k1*minX + b_k2*minY;
	gamma = g_k0 + g_k1*minX + g_k2*minY;
	*/

	
	for (int j = minY; j <= maxY; ++j)
	{
		for (int i = minX; i <= maxX; ++i)
		{
			// Unoptimized
			float tri_A = calc_area(i, j, Bx, By, Cx, Cy);
			float tri_B = calc_area(Ax, Ay, i, j, Cx, Cy);
			float tri_C = calc_area(Ax, Ay, Bx, By, i, j);
			
			alpha = tri_A / A_Triangle_Total;
			beta  = tri_B / A_Triangle_Total;
			gamma = tri_C / A_Triangle_Total;
			
			if ((0 <= alpha) && (0 <= beta) && (0 <= gamma) && (alpha <= 1) && (beta <= 1) && (gamma <= 1))
			{
				//color = alpha * color_v0 + beta * color_v1 + gamma * color_v3
				//Attempted optimization
				shade_pixel(i, j, in, state, alpha, beta, gamma);
				//state.image_color[(i+j*state.image_width)] = make_pixel(255, 255, 255);
			}
		}
	}
	
		//Extras
		//	Implement color interpolation by checking interp_rules in driver_state before sending color to the fragment_shader.
		// 		Only one interp_rule for each float in vertex.data. If the rule type is noperspective, interpolate float from
		//		3 vertices using barycentric coordinates.
	
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

