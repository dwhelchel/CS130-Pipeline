#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
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

    // Set image_color and image_depth to size width*height pixels
    unsigned long total_pixels = width * height;
    state.image_color = new pixel[total_pixels];
    state.image_depth = new float[total_pixels];

    // Set all pixels to black
    for (unsigned int i = 0; i < total_pixels; ++i) {
        state.image_color[i] = make_pixel(0, 0, 0);
    }

    // std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
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

    // Determine the type of rendering
    switch (type) {
        case render_type::invalid:
        break;

        /* Triangle */
        case render_type::triangle:
        {
            // Vertex array
            const data_geometry * geo[3];

            // Data geometry and data vertex variables
            data_geometry dg1;
            data_geometry dg2;
            data_geometry dg3;
            data_vertex dv1;
            data_vertex dv2;
            data_vertex dv3;

            // Allocate memory of floats_per_vertex size per data array
            dg1.data = new float[state.floats_per_vertex];
            dg2.data = new float[state.floats_per_vertex];
            dg3.data = new float[state.floats_per_vertex];
            dv1.data = new float[state.floats_per_vertex];
            dv2.data = new float[state.floats_per_vertex];
            dv3.data = new float[state.floats_per_vertex];

            for (int v = 0; v < state.num_vertices * state.floats_per_vertex; v += state.floats_per_vertex * 3) {

                // Assign appropriate data values from vertex_data into data_vertex
                for (int i = 0; i < state.floats_per_vertex; ++i) {
                    dv1.data[i] = state.vertex_data[v*state.floats_per_vertex+i];
                    dg1.data[i] = state.vertex_data[v*state.floats_per_vertex+i];
                }
                for (int i = 0; i < state.floats_per_vertex; ++i) {
                    dv2.data[i] = state.vertex_data[v+1*state.floats_per_vertex+i];
                    dg2.data[i] = state.vertex_data[v+1*state.floats_per_vertex+i];
                }
                for (int i = 0; i < state.floats_per_vertex; ++i) {
                    dv3.data[i] = state.vertex_data[v+2*state.floats_per_vertex+i];
                    dg3.data[i] = state.vertex_data[v+2*state.floats_per_vertex+i];
                }

                // Divide position by w
                dg1.gl_Position = dg1.gl_Position / dg1.gl_Position[3];
                dg2.gl_Position = dg2.gl_Position / dg2.gl_Position[3];
                dg3.gl_Position = dg3.gl_Position / dg3.gl_Position[3];

                // Set the pointers in geo to the proper locations
                geo[0] = &dg1;
                geo[1] = &dg2;
                geo[2] = &dg3;

                // Call vertex shader with data_vertex, data_geometry, and uniform_data
                state.vertex_shader(dv1, dg1, state.uniform_data);
                state.vertex_shader(dv2, dg2, state.uniform_data);
                state.vertex_shader(dv3, dg3, state.uniform_data);

                // Rasterize the triangle with state and new vertex array
                rasterize_triangle(state, geo);

            }

            // Deallocate memory used
            delete [] dg1.data;
            delete [] dg2.data;
            delete [] dg3.data;
            delete [] dv1.data;
            delete [] dv2.data;
            delete [] dv3.data;
        }
        // Exit
        break;

        /* Indexed */
        case render_type::indexed:
        break;

        /* Fan */
        case render_type::fan:
        break;

        /* Strip */
        case render_type::strip:
        break;
    }

    // std::cout<<"TODO: implement rendering."<<std::endl;
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

    // Vertex A
    double Ax = in[0]->gl_Position[0] * (state.image_width / 2) + ((state.image_width / 2) - 0.5);
    double Ay = in[0]->gl_Position[1] * (state.image_height / 2) + ((state.image_height / 2) - 0.5);

    // Vertex B
    double Bx = in[1]->gl_Position[0] * (state.image_width / 2) + ((state.image_width / 2) - 0.5);
    double By = in[1]->gl_Position[1] * (state.image_height / 2) + ((state.image_height / 2) - 0.5);

    // Vertex C
    double Cx = in[2]->gl_Position[0] * (state.image_width / 2) + ((state.image_width / 2) - 0.5);
    double Cy = in[2]->gl_Position[1] * (state.image_height / 2) + ((state.image_height / 2) - 0.5);

    // barycentric areas
    double totalArea = 0.5 * ((Bx*Cy - Cx*By) - (Ax*Cy - Cx*Ay) + (Ax*By - Bx*Ay));
    double alphaA = 0;
    double betaA = 0;
    double gammaA = 0;

    // loop over all pixels and do barycentric calculations
    for (int i = 0; i < state.image_width; ++i){
        for (int j = 0; j < state.image_height; ++j) {

            // Calculating alpha
            alphaA = 0.5 * ((Bx*Cy - Cx*By) + (By - Cy)*i + (Cx-Bx)*j);
            double alpha = alphaA / totalArea;

            // Calulcating beta
            betaA = 0.5 * ((Cx*Ay - Ax*Cy) + (Cy - Ay)*i + (Ax - Cx)*j);
            double beta = betaA / totalArea;

            // Calulating gamma
            gammaA = 0.5 * ((Ax*By - Bx*Ay) + (Ay - By)*i + (Bx - Ax)*j);
            double gamma = gammaA / totalArea;

            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                unsigned int index = i+j*state.image_width;
                state.image_color[index] = make_pixel(255, 255, 255);
            }
        }
    }

    // std::cout<<"TODO: implement rasterization"<<std::endl;
}
