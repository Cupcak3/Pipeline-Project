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
			std::cout<<"Type = triangle"<<std::endl;
			/*
			 *  v# are float arrays
			 *  v1 = x1 y1 z1
			 *  v2 = x2 y2 z2
			 *  v3 = x3 y3 z3
			 *  g_array = {v1, v2, v3} = (data) = {x1 y1 z1, x2 y2 z2, x3 y3 z3}
			 */
			
			//Read every 3 verticies into a data_geometry array, call rasterize_triangle
			data_geometry* g_array = new data_geometry[3];
			int k = 0;
			for (int i = 0; i < state.num_vertices*state.floats_per_vertex; i+=state.floats_per_vertex)
			{
				//g_array[k]'s data now points to the first vertex coordinate data address
				g_array[k].data = &state.vertex_data[i];
				++k;
			}
			/* Debugging for correctness
			std::cout<<"array data:\n{";
			for (int i = 0; i < 3; ++i)
			{
				std::cout<<"( ";
				for (int j = 0; j < state.floats_per_vertex; ++j)
				{
					std::cout<<g_array[i].data[j]<<" ";
				}
				std::cout<<")";
			}
			std::cout<<"}"<<std::endl;
			 */
			
			const data_geometry* g = g_array;
			//Call rasterize_triangle
			rasterize_triangle(state, &g);
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
    std::cout<<"TODO: implement rendering."<<std::endl;
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
	for (int i = 0; i < 3; ++i)
	{
		data_vertex * v_data = new data_vertex;
		data_geometry vertex = *in[i];
		
		v_data->data = in[i]->data;
		
		void (*shader) (const data_vertex&, data_geometry&, const float*) = state.vertex_shader;
		shader(*v_data, vertex, state.uniform_data); //vertex now has position filled in (?) as (x y z w)
		vertex.gl_Position /= vertex.gl_Position[3]; // Divide position by w (?)
		
		//Calculate pixel coordinates
		//	X and Y positions in NDC (each from -1 to 1)
		//		X from 0 to width, Y from 0 to height
		//  	NDC(-1,-1) is bottom left corner but not center of bottom left pixel
		// 		(x,y) in 2D NDC -------> (i,j) in pixel space
		// 		i =
		// 		j =
		
		//Draw verticies in image (useing image_color in driver_state)
		//	Make sure they fall on vertices of 00.png
		//	Already have (i,j) of pixel position. Determine specific pixel of color_image to set using (i,j)
		
		//Rasterize triangle
		//	Iterate over all pixels
		//		At pixel (i,j) use barycentric coordinates of pixel to determine if pixel is inside triangle or not.
		//		If inside set to white
		
		//Extras
		//	Determine square containing triangle. Only scan that square
		//		min/max x/y coordinates (x_min, x_max, y_min, y_max)
		//	Use Fragment shader to calculate pixel color rather than setting to white explicitly
		//		Use data_output in common.h and fragment_shader function in driver_state.h
		//	Implement color interpolation by checking interp_rules in driver_state before sending color to the fragment_shader.
		// 		Only one interp_rule for each float in vertex.data. If the rule type is noperspective, interpolate float from
		//		3 vertices using barycentric coordinates.
		
	}
    std::cout<<"TODO: implement rasterization"<<std::endl;
}

