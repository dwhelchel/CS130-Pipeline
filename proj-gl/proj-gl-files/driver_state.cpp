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
        state.image_depth[i] = 1.0f;
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

            for (int v = 0; v < state.num_vertices * state.floats_per_vertex; v += 3 * state.floats_per_vertex) {

                // Assign appropriate data values from vertex_data into data_vertex
                for (int i = 0; i < state.floats_per_vertex; ++i) {
                    dv1.data[i] = state.vertex_data[i + 0*state.floats_per_vertex+v];
                    dg1.data[i] = state.vertex_data[i + 0*state.floats_per_vertex+v];
                }
                for (int i = 0; i < state.floats_per_vertex; ++i) {
                    dv2.data[i] = state.vertex_data[i + 1*state.floats_per_vertex+v];
                    dg2.data[i] = state.vertex_data[i + 1*state.floats_per_vertex+v];
                }
                for (int i = 0; i < state.floats_per_vertex; ++i) {
                    dv3.data[i] = state.vertex_data[i + 2*state.floats_per_vertex+v];
                    dg3.data[i] = state.vertex_data[i + 2*state.floats_per_vertex+v];
                }

                // Set the pointers in geo to the proper locations
                geo[0] = &dg1;
                geo[1] = &dg2;
                geo[2] = &dg3;

                // Call vertex shader with data_vertex, data_geometry, and uniform_data
                state.vertex_shader(dv1, dg1, state.uniform_data);
                state.vertex_shader(dv2, dg2, state.uniform_data);
                state.vertex_shader(dv3, dg3, state.uniform_data);

                // Rasterize the triangle with state and new vertex array
                // rasterize_triangle(state, geo);
                clip_triangle(state, geo, 0);

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
        {

            const data_geometry * geo[3];
            data_geometry dg1;
            data_geometry dg2;
            data_geometry dg3;
            data_vertex dv1;
            data_vertex dv2;
            data_vertex dv3;

            for (int i = 0; i < state.num_triangles * 3; i += 3) {

                dv1.data = &state.vertex_data[state.index_data[i+0] * state.floats_per_vertex];
                dg1.data = dv1.data;

                dv2.data = &state.vertex_data[state.index_data[i+1] * state.floats_per_vertex];
                dg2.data = dv2.data;

                dv3.data = &state.vertex_data[state.index_data[i+2] * state.floats_per_vertex];
                dg3.data = dv3.data;

                state.vertex_shader(dv1, dg1, state.uniform_data);
                state.vertex_shader(dv2, dg2, state.uniform_data);
                state.vertex_shader(dv3, dg3, state.uniform_data);

                geo[0] = &dg1;
                geo[1] = &dg2;
                geo[2] = &dg3;

                clip_triangle(state, geo, 0);

            }
        }
        break;

        /* Fan */
        case render_type::fan:
        {

            const data_geometry * geo[3];
            data_geometry dg[3];
            data_vertex dv[3];

            int count = 0;

            for (int i = 0; i < state.num_vertices; ++i) {
                for (int j = 0; j < 3; ++j) {
                    count = i+j;
                    if (j == 0) count = 0;

                    dv[j].data = &state.vertex_data[count * state.floats_per_vertex];

                    dg[j].data = dv[j].data;
                    state.vertex_shader(dv[j], dg[j], state.uniform_data);
                    geo[j] = &dg[j];
                }
                clip_triangle(state, geo, 0);
            }

        }
        break;

        /* Strip */
        case render_type::strip:
        {

            const data_geometry * geo[3];
            data_geometry dg1;
            data_geometry dg2;
            data_geometry dg3;
            data_vertex dv1;
            data_vertex dv2;
            data_vertex dv3;

            for (int i = 0; i < state.num_vertices - 2; ++i) {

                dv1.data = &state.vertex_data[i * state.floats_per_vertex];
                dv2.data = &state.vertex_data[(i+1) * state.floats_per_vertex];
                dv3.data = &state.vertex_data[(i+2) * state.floats_per_vertex];
                dg1.data = dv1.data;
                dg2.data = dv2.data;
                dg3.data = dv3.data;

                state.vertex_shader(dv1, dg1, state.uniform_data);
                state.vertex_shader(dv2, dg2, state.uniform_data);
                state.vertex_shader(dv3, dg3, state.uniform_data);

                geo[0] = &dg1;
                geo[1] = &dg2;
                geo[2] = &dg3;

                clip_triangle(state, geo, 0);

            }

        }
        break;
    }

    // std::cout<<"TODO: implement rendering."<<std::endl;
}

float generate_alpha(driver_state& state, bool sign, int position, const data_geometry* in[3], int a, int b) {

    float alpha = 0;

    if (sign) {
        alpha = (in[b]->gl_Position[3] - in[b]->gl_Position[position])
                / (in[a]->gl_Position[position] - in[a]->gl_Position[3] + in[b]->gl_Position[3] - in[b]->gl_Position[position]);
    } else {
        alpha = (-1 * in[b]->gl_Position[3] - in[b]->gl_Position[position])
                / (in[a]->gl_Position[position] + in[a]->gl_Position[3] - in[b]->gl_Position[3] - in[b]->gl_Position[position]);
    }

    return alpha;

}

// Checks to see if vertices of the triangle are on opposite sides
void check_vertices(driver_state& state, bool sign, int position, const data_geometry* in[3], int face) {

    // True = inside, false = outside
    bool vertexA = false;
    bool vertexB = false;
    bool vertexC = false;

    float alpha_0;
    float alpha_1;

    // New data geometries
    const data_geometry * geo[3];
    const data_geometry * geo2[3];
    data_geometry dg1;
    data_geometry dg2;
    data_geometry dg3;
    data_geometry dg_1;

    dg1.data = new float[state.floats_per_vertex];
    dg2.data = new float[state.floats_per_vertex];
    dg3.data = new float[state.floats_per_vertex];
    dg_1.data = new float[state.floats_per_vertex];

    if (sign) {
        if (in[0]->gl_Position[position] <= in[0]->gl_Position[3]) vertexA = true;
        if (in[1]->gl_Position[position] <= in[1]->gl_Position[3]) vertexB = true;
        if (in[2]->gl_Position[position] <= in[2]->gl_Position[3]) vertexC = true;
    } else {
        if (in[0]->gl_Position[position] >= -1 * in[0]->gl_Position[3]) vertexA = true;
        if (in[1]->gl_Position[position] >= -1 * in[1]->gl_Position[3]) vertexB = true;
        if (in[2]->gl_Position[position] >= -1 * in[2]->gl_Position[3]) vertexC = true;
    }

    // If all inside nothing to clip against this plane
    if (vertexA && vertexB && vertexC) {
        clip_triangle(state, in, face+1);
    }

    // For middle cases
    // A, B outside C inside
    if (!vertexA && !vertexB && vertexC) { // bc, ca ; call once
        alpha_0 = generate_alpha(state, sign, position, in, 2, 0);
        alpha_1 = generate_alpha(state, sign, position, in, 2, 1);

        vec4 position_0 = alpha_0 * in[2]->gl_Position + (1 - alpha_0) * in[0]->gl_Position; // ca
        vec4 position_1 = alpha_1 * in[2]->gl_Position + (1 - alpha_1) * in[1]->gl_Position; // cb

        dg1.gl_Position = in[2]->gl_Position;
        dg2.gl_Position = position_0;
        dg3.gl_Position = position_1;

        for (int i = 0; i < state.floats_per_vertex; ++i) {
            if (state.interp_rules[i] == interp_type::flat) {
                dg1.data[i] = in[2]->data[i];
                dg2.data[i] = in[0]->data[i];
                dg3.data[i] = in[1]->data[i];
            } else if (state.interp_rules[i] == interp_type::smooth) {
                dg1.data[i] = in[2]->data[i];
                dg2.data[i] = alpha_0 * in[2]->data[i] * (1 - alpha_0) * in[0]->data[i]; // ca
                dg3.data[i] = alpha_1 * in[2]->data[i] * (1 - alpha_1) * in[1]->data[i]; // cb
            } else if (state.interp_rules[i] == interp_type::noperspective) {
                float ca_w = 1.0 / (alpha_0 * in[2]->gl_Position[3] + (1 - alpha_0) * in[0]->gl_Position[3]);
                float cb_w = 1.0 / (alpha_1 * in[2]->gl_Position[3] + (1 - alpha_1) * in[1]->gl_Position[3]);
                float ca_nonalpha = alpha_0 * in[2]->gl_Position[3] * ca_w;
                float cb_nonalpha = alpha_1 * in[2]->gl_Position[3] * cb_w;
                dg1.data[i] = in[2]->data[i];
                dg2.data[i] = ca_nonalpha * in[2]->data[i] + (1 - ca_nonalpha) * in[0]->data[i];
                dg3.data[i] = cb_nonalpha * in[2]->data[i] + (1 - cb_nonalpha) * in[1]->data[i];
            }
        }

        geo[0] = &dg1;
        geo[1] = &dg2;
        geo[2] = &dg3;

        clip_triangle(state, geo, face+1);
    }
    // A, C outside B inside
    if (!vertexA && vertexB && !vertexC) { // ab, bc ; call once
        alpha_0 = generate_alpha(state, sign, position, in, 1, 2); // bc
        alpha_1 = generate_alpha(state, sign, position, in, 1, 0); // ba

        vec4 position_0 = alpha_0 * in[1]->gl_Position + (1 - alpha_0) * in[2]->gl_Position; //bc
        vec4 position_1 = alpha_1 * in[1]->gl_Position + (1 - alpha_1) * in[0]->gl_Position; //ba

        dg1.gl_Position = in[1]->gl_Position; // b position
        dg2.gl_Position = position_0; // bc position
        dg3.gl_Position = position_1; // ba position

        for (int i = 0; i < state.floats_per_vertex; ++i) {
            if (state.interp_rules[i] == interp_type::flat) {
                dg1.data[i] = in[1]->data[i];
                dg2.data[i] = in[0]->data[i];
                dg3.data[i] = in[2]->data[i];
            } else if (state.interp_rules[i] == interp_type::smooth) {
                dg1.data[i] = in[1]->data[i];
                dg2.data[i] = alpha_0 * in[1]->data[i] + (1 - alpha_0) * in[2]->data[i];
                dg3.data[i] = alpha_1 * in[1]->data[i] + (1 - alpha_1) * in[0]->data[i];
            } else if (state.interp_rules[i] == interp_type::noperspective) {
                float bc_w = 1.0 / (alpha_0 * in[1]->gl_Position[3] + (1 - alpha_0) * in[2]->gl_Position[3]);
                float ba_w = 1.0 / (alpha_1 * in[1]->gl_Position[3] + (1 - alpha_1) * in[0]->gl_Position[3]);
                float bc_nonalpha = alpha_0 * in[1]->gl_Position[3] * bc_w;
                float ba_nonalpha = alpha_1 * in[1]->gl_Position[3] * ba_w;
                dg1.data[i] = in[1]->data[i];
                dg2.data[i] = bc_nonalpha * in[1]->data[i] + (1 - bc_nonalpha) * in[2]->data[i];
                dg3.data[i] = ba_nonalpha * in[1]->data[i] + (1 - ba_nonalpha) * in[0]->data[i];
            }
        }

        geo[0] = &dg1;
        geo[1] = &dg2;
        geo[2] = &dg3;

        clip_triangle(state, geo, face+1);
    }
    // A outside B, C inside
    if (!vertexA && vertexB && vertexC) { // ab, ca ; call twice
        alpha_0 = generate_alpha(state, sign, position, in, 2, 0); //ca
        alpha_1 = generate_alpha(state, sign, position, in, 1, 0); //ba

        vec4 position_0 = alpha_0 * in[2]->gl_Position + (1 - alpha_0) * in[0]->gl_Position; // ca
        vec4 position_1 = alpha_1 * in[1]->gl_Position + (1 - alpha_1) * in[0]->gl_Position; // ba

        dg1.gl_Position = in[1]->gl_Position;
        dg2.gl_Position = position_0;
        dg3.gl_Position = position_1;
        dg_1.gl_Position = in[2]->gl_Position;

        for (int i = 0; i < state.floats_per_vertex; ++i) {
            if (state.interp_rules[i] == interp_type::flat) {
                dg1.data[i] = in[1]->data[i];
                dg2.data[i] = in[0]->data[i];
                dg3.data[i] = in[0]->data[i];
                dg_1.data[i] = in[2]->data[i];
            } else if (state.interp_rules[i] == interp_type::smooth) {
                dg1.data[i] = in[1]->data[i]; // b data
                dg_1.data[i] = in[2]->data[i]; // c data
                dg2.data[i] = alpha_0 * in[2]->data[i] + (1 - alpha_0) * in[0]->data[i];
                dg3.data[i] = alpha_1 * in[1]->data[i] + (1 - alpha_1) * in[0]->data[i];
            } else if (state.interp_rules[i] == interp_type::noperspective) {
                float ca_w = 1.0 / (alpha_0 * in[2]->gl_Position[3] + (1 - alpha_0) * in[0]->gl_Position[3]);
                float ba_w = 1.0 / (alpha_1 * in[1]->gl_Position[3] + (1 - alpha_1) * in[0]->gl_Position[3]);
                float ca_noalpha = alpha_0 * in[0]->gl_Position[3] * ca_w;
                float ba_nonalpha = alpha_1 * in[0]->gl_Position[3] * ba_w;
                dg1.data[i] = in[1]->data[i];
                dg_1.data[i] = in[2]->data[i];
                dg2.data[i] = ca_noalpha * in[2]->data[i] + (1 - ca_noalpha) * in[0]->data[i];
                dg3.data[i] = ba_nonalpha * in[1]->data[i] + (1 - ba_nonalpha) * in[0]->data[i];
            }
        }

        geo[0] = &dg1; // b
        geo[1] = &dg2; // ca
        geo[2] = &dg3; // ba

        geo2[0] = &dg1; // b
        geo2[1] = &dg_1; // c
        geo2[2] = &dg2; // ca

        clip_triangle(state, geo, face+1);
        clip_triangle(state, geo2, face+1);
    }
    // A inside B, C outside
    if (vertexA && !vertexB && !vertexC) { // ab, ca ; call once
        alpha_0 = generate_alpha(state, sign, position, in, 0, 1);
        alpha_1 = generate_alpha(state, sign, position, in, 0, 2);

        vec4 position_0 = alpha_0 * in[0]->gl_Position + (1 - alpha_0) * in[1]->gl_Position; // ab
        vec4 position_1 = alpha_1 * in[0]->gl_Position + (1 - alpha_1) * in[2]->gl_Position; // ac

        dg1.gl_Position = in[0]->gl_Position;
        dg2.gl_Position = position_0;
        dg3.gl_Position = position_1;

        for (int i = 0; i < state.floats_per_vertex; ++i) {
            if (state.interp_rules[i] == interp_type::flat) {
                dg1.data[i] = in[0]->data[i];
                dg2.data[i] = in[1]->data[i];
                dg3.data[i] = in[2]->data[i];
            } else if (state.interp_rules[i] == interp_type::smooth) {
                dg1.data[i] = in[0]->data[i]; // a data
                dg2.data[i] = alpha_0 * in[0]->data[i] + (1 - alpha_0) * in[1]->data[i];
                dg3.data[i] = alpha_1 * in[0]->data[i] + (1 - alpha_1) * in[2]->data[i];
            } else if (state.interp_rules[i] == interp_type::noperspective) {
                float ab_w = 1.0 / (alpha_0 * in[0]->gl_Position[3] + (1 - alpha_0) * in[1]->gl_Position[3]);
                float ac_w = 1.0 / (alpha_1 * in[0]->gl_Position[3] + (1 - alpha_1) * in[2]->gl_Position[3]);
                float ab_nonalpha = alpha_0 * in[0]->gl_Position[3] * ab_w;
                float ac_nonalpha = alpha_1 * in[0]->gl_Position[3] * ac_w;
                dg1.data[i] = in[0]->data[i];
                dg2.data[i] = ab_nonalpha * in[0]->data[i] + (1 - ab_nonalpha) * in[1]->data[i];
                dg3.data[i] = ac_nonalpha * in[0]->data[i] + (1 - ac_nonalpha) * in[2]->data[i];
            }
        }

        geo[0] = &dg1;
        geo[1] = &dg2;
        geo[2] = &dg3;

        clip_triangle(state, geo, face+1);
    }
    // A, C inside B outside
    if (vertexA && !vertexB && vertexC) { // bc, ab ; call twice
        alpha_0 = generate_alpha(state, sign, position, in, 0, 1); //ab
        alpha_1 = generate_alpha(state, sign, position, in, 2, 1); //cb

        vec4 position_0 = alpha_0 * in[0]->gl_Position + (1 - alpha_0) * in[1]->gl_Position; // ab
        vec4 position_1 = alpha_1 * in[2]->gl_Position + (1 - alpha_1) * in[1]->gl_Position; // cb

        dg1.gl_Position = in[2]->gl_Position; // c position
        dg2.gl_Position = position_0; // ab
        dg3.gl_Position = position_1; // cb
        dg_1.gl_Position = in[0]->gl_Position; // a

        for (int i = 0; i < state.floats_per_vertex; ++i) {
            if (state.interp_rules[i] == interp_type::flat) {
                dg1.data[i] = in[2]->data[i];
                dg2.data[i] = in[1]->data[i];
                dg3.data[i] = in[1]->data[i];
                dg_1.data[i] = in[0]->data[i];
            } else if (state.interp_rules[i] == interp_type::smooth) {
                dg1.data[i] = in[2]->data[i]; // c data
                dg_1.data[i] = in[0]->data[i]; // a data
                dg2.data[i] = alpha_0 * in[0]->data[i] + (1 - alpha_0) * in[1]->data[i]; //ab
                dg3.data[i] = alpha_1 * in[2]->data[i] + (1 - alpha_1) * in[1]->data[i]; //cb
            } else if (state.interp_rules[i] == interp_type::noperspective) {
                float ab_w = 1.0 / (alpha_0 * in[0]->gl_Position[3] + (1 - alpha_0) * in[1]->gl_Position[3]);
                float cb_w = 1.0 / (alpha_1 * in[2]->gl_Position[3] + (1 - alpha_1) * in[1]->gl_Position[3]);
                float ab_noalpha = alpha_0 * in[1]->gl_Position[3] * ab_w;
                float cb_nonalpha = alpha_1 * in[1]->gl_Position[3] * cb_w;
                dg1.data[i] = in[2]->data[i];
                dg_1.data[i] = in[0]->data[i];
                dg2.data[i] = ab_noalpha * in[0]->data[i] + (1 - ab_noalpha) * in[1]->data[i];
                dg3.data[i] = cb_nonalpha * in[2]->data[i] + (1 - cb_nonalpha) * in[1]->data[i];
            }
        }

        geo[0] = &dg1; // c
        geo[1] = &dg2; // ab
        geo[2] = &dg3; // cb

        geo2[0] = &dg1; // c
        geo2[1] = &dg_1; // a
        geo2[2] = &dg2; // ab

        clip_triangle(state, geo, face+1);
        clip_triangle(state, geo2, face+1);
    }
    // A, B inside C outside
    if (vertexA && vertexB && !vertexC) { // ca, bc ; call twice
        alpha_0 = generate_alpha(state, sign, position, in, 1, 2); // bc
        alpha_1 = generate_alpha(state, sign, position, in, 0, 2); // ac

        vec4 position_0 = alpha_0 * in[1]->gl_Position + (1 - alpha_0) * in[2]->gl_Position; // bc
        vec4 position_1 = alpha_1 * in[0]->gl_Position + (1 - alpha_1) * in[2]->gl_Position; // ac

        dg1.gl_Position = in[0]->gl_Position; // a position
        dg2.gl_Position = position_0; // bc
        dg3.gl_Position = position_1; // ac
        dg_1.gl_Position = in[1]->gl_Position; // b

        for (int i = 0; i < state.floats_per_vertex; ++i) {
            if (state.interp_rules[i] == interp_type::flat) {
                dg1.data[i] = in[0]->data[i];
                dg2.data[i] = in[2]->data[i];
                dg3.data[i] = in[2]->data[i];
                dg_1.data[i] = in[1]->data[i];
            } else if (state.interp_rules[i] == interp_type::smooth) {
                dg1.data[i] = in[0]->data[i]; // a data
                dg_1.data[i] = in[1]->data[i]; // b data
                dg2.data[i] = alpha_0 * in[1]->data[i] + (1 - alpha_0) * in[2]->data[i]; // bc
                dg3.data[i] = alpha_1 * in[0]->data[i] + (1 - alpha_1) * in[2]->data[i]; // ac
            } else if (state.interp_rules[i] == interp_type::noperspective) {
                float bc_w = 1.0 / (alpha_0 * in[1]->gl_Position[3] + (1 - alpha_0) * in[2]->gl_Position[3]);
                float ac_w = 1.0 / (alpha_1 * in[0]->gl_Position[3] + (1 - alpha_1) * in[2]->gl_Position[3]);
                float bc_noalpha = alpha_0 * in[2]->gl_Position[3] * bc_w;
                float ac_nonalpha = alpha_1 * in[2]->gl_Position[3] * ac_w;
                dg1.data[i] = in[0]->data[i];
                dg_1.data[i] = in[1]->data[i];
                dg2.data[i] = bc_noalpha * in[1]->data[i] + (1 - bc_noalpha) * in[2]->data[i];
                dg3.data[i] = ac_nonalpha * in[0]->data[i] + (1 - ac_nonalpha) * in[2]->data[i];
            }
        }

        geo[0] = &dg1; // a
        geo[1] = &dg2; // bc
        geo[2] = &dg3; // ac

        geo2[0] = &dg1; // a
        geo2[1] = &dg_1; // b
        geo2[2] = &dg2; // bc

        clip_triangle(state, geo, face+1);
        clip_triangle(state, geo2, face+1);
    }

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

    // True = inside, false = outside
    bool sign = false;
    int position = 0;

    // clip right plane (+1)
    if (face == 0) {
        sign = true;
        position = 0;
        check_vertices(state, sign, position, in, face);
    }
    // clip left plane (-1)
    else if (face == 1) {
        sign = false;
        position = 0;
        check_vertices(state, sign, position, in, face);
    }
    // clip top plane (+1)
    else if (face == 2) {
        sign = true;
        position = 1;
        check_vertices(state, sign, position, in, face);
    }
    // clip bottom plane (-1)
    else if (face == 3) {
        sign = false;
        position = 1;
        check_vertices(state, sign, position, in, face);
    }
    // clip far plane (+1)
    else if (face == 4) {
        sign = true;
        position = 2;
        check_vertices(state, sign, position, in, face);
    }
    // clip near plane (-1)
    else if (face == 5){
        sign = false;
        position = 2;
        check_vertices(state, sign, position, in, face);
    }
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{

    // Store w values
    const float a_w = in[0]->gl_Position[3];
    const float b_w = in[1]->gl_Position[3];
    const float c_w = in[2]->gl_Position[3];

    // Divide position by w
    vec4 a_position = in[0]->gl_Position / in[0]->gl_Position[3];
    vec4 b_position = in[1]->gl_Position / in[1]->gl_Position[3];
    vec4 c_position = in[2]->gl_Position / in[2]->gl_Position[3];

    // Vertex A
    double Ax = a_position[0] * (state.image_width / 2) + ((state.image_width / 2) - 0.5);
    double Ay = a_position[1] * (state.image_height / 2) + ((state.image_height / 2) - 0.5);

    // Vertex B
    double Bx = b_position[0] * (state.image_width / 2) + ((state.image_width / 2) - 0.5);
    double By = b_position[1] * (state.image_height / 2) + ((state.image_height / 2) - 0.5);

    // Vertex C
    double Cx = c_position[0] * (state.image_width / 2) + ((state.image_width / 2) - 0.5);
    double Cy = c_position[1] * (state.image_height / 2) + ((state.image_height / 2) - 0.5);

    // barycentric areas
    double totalArea = 0.5 * ((Bx*Cy - Cx*By) - (Ax*Cy - Cx*Ay) + (Ax*By - Bx*Ay));
    double alphaA = 0;
    double betaA = 0;
    double gammaA = 0;

    // Data fragment and Data output
    data_fragment df;
    data_output dout;

    // Allocate memory temp array
    df.data = new float[MAX_FLOATS_PER_VERTEX];

    // New barycentric coordinates and k value
    double newAlpha, newBeta, newGamma, k;

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

            // std::cout << alphaA << " " << betaA << " " << gammaA << std::endl;

            // Check if point lies within triangle
            if (alpha >= 0 && beta >= 0 && gamma >= 0) {

                unsigned int index = i+j*state.image_width;

                float depth = (alpha * a_position[2]) +
                               (beta * b_position[2]) +
                               (gamma * c_position[2]);

                if (depth < state.image_depth[index]) {
                    for (int i = 0; i < state.floats_per_vertex; ++i) {
                        if (state.interp_rules[i] == interp_type::flat) {
                            df.data[i] = in[0]->data[i];
                        }
                        else if (state.interp_rules[i] == interp_type::smooth) {

                            k = (alpha / a_w) + (beta / b_w) + (gamma / c_w);

                            newAlpha = alpha / (a_w * k);
                            newBeta = beta / (b_w * k);
                            newGamma = gamma / (c_w * k);

                            df.data[i] = newAlpha * in[0]->data[i] +
                                           newBeta * in[1]->data[i] +
                                           newGamma * in[2]->data[i];
                        }
                        else if (state.interp_rules[i] == interp_type::noperspective) {
                            df.data[i] = alpha * in[0]->data[i] +
                                           beta * in[1]->data[i] +
                                           gamma * in[2]->data[i];
                        }
                    }

                    state.image_depth[index] = depth;
                    state.fragment_shader(df, dout, state.uniform_data);
                    state.image_color[index] = make_pixel(dout.output_color[0] * 255,
                                                          dout.output_color[1] * 255,
                                                          dout.output_color[2] * 255);
                }

            }
        }
    }

    delete [] df.data;

}
