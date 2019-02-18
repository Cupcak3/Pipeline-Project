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

	for (int i = 0; i < state.image_height * state.image_width; ++i)
	{
		//0-255 with 0,0,0 being black and 255,255,255 being white
		state.image_color[i] = make_pixel(0, 0, 0); // Initialize all pixels to black
	}
	
	//std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
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
				
				state.vertex_shader(v1,g1,state.uniform_data);
				state.vertex_shader(v2,g2,state.uniform_data);
				state.vertex_shader(v3,g3,state.uniform_data);
			
				g1.gl_Position[0] *= state.image_width/2;
				g1.gl_Position[0] += ((state.image_width/2)-.5);
				g1.gl_Position[1] *= state.image_height/2;
				g1.gl_Position[1] += ((state.image_height/2)-.5);
				
				g2.gl_Position[0] *= state.image_width/2;
				g2.gl_Position[0] += ((state.image_width/2)-.5);
				g2.gl_Position[1] *= state.image_height/2;
				g2.gl_Position[1] += ((state.image_height/2)-.5);
				
				g3.gl_Position[0] *= state.image_width/2;
				g3.gl_Position[0] += ((state.image_width/2)-.5);
				g3.gl_Position[1] *= state.image_height/2;
				g3.gl_Position[1] += ((state.image_height/2)-.5);
				
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

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
	/*
	for (int k = 0; k < 3; ++k)
	{
		//Calculate pixel coordinates
		//	X and Y positions in NDC (each from -1 to 1)
		//		X from 0 to width, Y from 0 to height
		//  	NDC(-1,-1) is bottom left corner but not center of bottom left pixel
		// 		(x,y) in 2D NDC -------> (i,j) in pixel space
		// 		i = (w/2)x + (w/2 - .5)
		// 		j = (h/2)y + (h/2 - .5)
		
		int i = (state.image_width/2) * in[k]->gl_Position[0] + ((state.image_width/2) - .5);
		int j = (state.image_height/2) * in[k]->gl_Position[1] + ((state.image_height/2) - .5);
		
		// Draw verticies in image (using image_color in driver_state)
		// 	Make sure they fall on vertices of 00.png
		//	Already have (i,j) of pixel position. Determine specific pixel of color_image to set using (i,j)
		// (i,j) --> 0 1 2 3 4 5 ... ij ... width*height
	
		//int pixel_index = i + j * state.image_width;
		//state.image_color[pixel_index] = make_pixel(255, 255, 255);
		//std::cout<<"("<<i<<", "<<j<<")"<<std::endl;
	}*/
	

	//return;
	//Rasterize triangle
	//	Iterate over all pixels
	//		At pixel (i,j) use barycentric coordinates of pixel to determine if pixel is inside triangle or not.
	//		If inside set to white
	/*
	int Ax = (state.image_width/2) * in[0]->gl_Position[0] + ((state.image_width/2) - .5);
	int Ay = (state.image_height/2) * in[0]->gl_Position[1] + ((state.image_height/2) - .5);
	
	int Bx = (state.image_width/2) * in[1]->gl_Position[0] + ((state.image_width/2) - .5);
	int By = (state.image_height/2) * in[1]->gl_Position[1] + ((state.image_height/2) - .5);
	
	int Cx = (state.image_width/2) * in[2]->gl_Position[0] + ((state.image_width/2) - .5);
	int Cy = (state.image_height/2) * in[2]->gl_Position[1] + ((state.image_height/2) - .5);
	*/
	
	int Ax = in[0]->gl_Position[0];
	int Ay = in[0]->gl_Position[1];
	
	int Bx = in[1]->gl_Position[0];
	int By = in[1]->gl_Position[1];
	
	int Cx = in[2]->gl_Position[0];
	int Cy = in[2]->gl_Position[1];
	
	
	float A_Triangle_Total = 0.5 * ((Bx*Cy - Cx*By) - (Ax*Cy - Cx*Ay) - (Ax*By - Bx*Ax));
	
	int maxX = std::max({Ax, Bx, Cx});
	int minX = std::min({Ax, Bx, Cx});
	
	int maxY = std::max({Ay, By, Cy});
	int minY = std::min({Ay, By, Cy});
	
	float alpha = -1, beta = -1, gamma = -1;
	/* Attempted optimizations
	float k0 = (0.5 * (Bx*Cy - Cx*By)) / A_Triangle_Total;
	float k1 = ((0.5 * ((Bx*Cy - Cx*By) - (Cy) - (By - Bx))) / A_Triangle_Total) - k0;
	float k2 = ((0.5 * ((Bx*Cy - Cx*By) - (Cx))) / A_Triangle_Total) - k0;
	*/
	
	for (int j = minY; j < maxY; ++j)
	{
		
		for (int i = minX+1; i < maxX; ++i)
		{
			// Unoptimized
			float tri_A = 0.5 * ((Bx*Cy - Cx*By) - (i*Cy - Cx*j) - (i*By - Bx*i));
			float tri_B = 0.5 * ((i*Cy - Cx*j) - (Ax*Cy - Cx*Ay) - (Ax*j - j*Ax));
			float tri_C = 0.5 * ((Bx*j - j*By) - (Ax*j - i*Ay) - (Ax*By - Bx*Ax));
			
			alpha = tri_A / A_Triangle_Total;
			beta  = tri_B / A_Triangle_Total;
			gamma = tri_C / A_Triangle_Total;
			
			
			if (alpha >= 0 && beta >= 0 && gamma >= 0)
			{
				//Attempted optimization
				data_output output;
				const data_fragment fragment{};
				float * frag_data = new float[MAX_FLOATS_PER_VERTEX];
				state.fragment_shader(fragment, output, state.uniform_data);
				
				output.output_color *= 255;
				state.image_color[i+j*state.image_width] = make_pixel(output.output_color[0], output.output_color[1], output.output_color[2]);
				delete [] frag_data;
			}
		}
	}
		
		//Extras
		//	Use Fragment shader to calculate pixel color rather than setting to white explicitly
		//		Use data_output in common.h and fragment_shader function in driver_state.h
		//	Implement color interpolation by checking interp_rules in driver_state before sending color to the fragment_shader.
		// 		Only one interp_rule for each float in vertex.data. If the rule type is noperspective, interpolate float from
		//		3 vertices using barycentric coordinates.
		
    std::cout<<"TODO: implement rasterization"<<std::endl;
}

