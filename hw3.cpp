#define COW_CRASH_ON_FLOATING_POINT_EXCEPTIONS
#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_stretchy_buffer.cpp" // bare-bones templatized stretchy buffer

////////////////////////////////////////////////////////////////////////////////
// documentation                                                              //
////////////////////////////////////////////////////////////////////////////////

// mouse wheel to zoom

// HUGE is a real big number

// vec3 cross(vec3 a, vec3 b); // a x b
// double squaredNorm(vec3 v); // |v|^2
// double norm(vec3 v);        // |v|
// vec3 normalized(vec3 v);    // v / |v|

// int3 is three contiguous ints
//      and works a lot like vec3
#if 0
int3 triangle = {};   // (0, 0, 0)
triangle.i = 4;       // (4, 0, 0)
triangle[1] = 5;      // (4, 5, 0)
triangle.data[2] = 6; // (4, 5, 6)
#endif

// // soup mesh
//
// struct BasicTriangleMesh3D {
//     int num_vertices;
//     vec3 *vertex_positions;
// };

// // indexed mesh
//
// struct FancyTriangleMesh3D {
//     int num_vertices;
//     int num_triangles;
//     vec3 *vertex_positions;
//     int3 *triangle_indices;
//     vec3 *vertex_normals;
// };

// https://cplusplus.com/reference/cstdio/sscanf/
// https://cplusplus.com/reference/cstring/memcpy/

////////////////////////////////////////////////////////////////////////////////
// hw                                                                         //
////////////////////////////////////////////////////////////////////////////////

// begin please ignore these lines
void mesh_transform_vertex_positions_to_double_unit_box(int num_vertices, vec3 *vertex_positions);
void fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals(FancyTriangleMesh3D *fancy_mesh);
void fancy_mesh_merge_duplicated_vertices(FancyTriangleMesh3D *fancy_mesh);
// end please ignore these lines

// this function is already complete; you are free to use it as reference for load_basic_mesh(...)
FancyTriangleMesh3D load_fancy_mesh(char *filename, bool transform_vertex_positions_to_double_unit_box, bool compute_normals, bool unify_duplicated_vertices)
{
    FancyTriangleMesh3D fancy_mesh = {};
    {
        StretchyBuffer<vec3> vertex_positions = {};
        StretchyBuffer<int3> triangle_indices = {};
        {
            FILE *fp = fopen(filename, "r");
            ASSERT(fp);
            char buffer[4096];
            while (fgets(buffer, NELEMS(buffer), fp) != NULL)
            {
                char prefix[16] = {};
                sscanf(buffer, "%s", prefix);
                if (strcmp(prefix, "v") == 0)
                {
                    double x, y, z;
                    ASSERT(sscanf(buffer, "%s %lf %lf %lf", prefix, &x, &y, &z) == 4);
                    sbuff_push_back(&vertex_positions, {x, y, z});
                }
                else if (strcmp(prefix, "f") == 0)
                {
                    int i, j, k;
                    ASSERT(sscanf(buffer, "%s %d %d %d", prefix, &i, &j, &k) == 4);
                    sbuff_push_back(&triangle_indices, {i - 1, j - 1, k - 1});
                }
            }
            fclose(fp);
        }
        // note: don't free the data pointers! (we're stealing them)
        fancy_mesh.num_triangles = triangle_indices.length;
        fancy_mesh.triangle_indices = triangle_indices.data;
        fancy_mesh.num_vertices = vertex_positions.length;
        fancy_mesh.vertex_positions = vertex_positions.data;
    }
    if (transform_vertex_positions_to_double_unit_box)
    {
        mesh_transform_vertex_positions_to_double_unit_box(fancy_mesh.num_vertices, fancy_mesh.vertex_positions);
    }
    if (unify_duplicated_vertices)
    {
        fancy_mesh_merge_duplicated_vertices(&fancy_mesh);
    }
    if (compute_normals)
    {
        fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals(&fancy_mesh);
    }
    return fancy_mesh;
}

// begin submission

BasicTriangleMesh3D load_basic_mesh(char *filename, bool transform_vertex_positions_to_double_unit_box)
{
    BasicTriangleMesh3D basic_mesh = {};
    StretchyBuffer<vec3> vertex_positions = {};
    {
        FILE *fp = fopen(filename, "r");
        ASSERT(fp);
        char buffer[4096];
        while (fgets(buffer, NELEMS(buffer), fp) != NULL)
        {
            double x, y, z;
            ASSERT(sscanf(buffer, "%lf %lf %lf", &x, &y, &z) == 3);
            sbuff_push_back(&vertex_positions, {x, y, z});
        }
        fclose(fp);
    }
    basic_mesh.num_vertices = vertex_positions.length;
    basic_mesh.vertex_positions = vertex_positions.data; // NOTE stealing data pointer

    if (transform_vertex_positions_to_double_unit_box)
    {
        mesh_transform_vertex_positions_to_double_unit_box(basic_mesh.num_vertices, basic_mesh.vertex_positions);
    }
    return basic_mesh;
}

BasicTriangleMesh3D fancy2basic(FancyTriangleMesh3D fancy_mesh)
{
    BasicTriangleMesh3D basic_mesh = {};
    {
        basic_mesh.num_vertices = 3 * fancy_mesh.num_triangles;
        basic_mesh.vertex_positions = (vec3 *)malloc(basic_mesh.num_vertices * sizeof(vec3));
        for (int i = 0; i < fancy_mesh.num_triangles; i++)
        {
            int3 triangle = fancy_mesh.triangle_indices[i];
            for (int k = 0; k < 3; k++) // !!!loop through vertex positions!!!
            {
                basic_mesh.vertex_positions[k + 3 * i] = fancy_mesh.vertex_positions[triangle[k]];
            }
        }
    }
    return basic_mesh;
}

void mesh_transform_vertex_positions_to_double_unit_box(int num_vertices, vec3 *vertex_positions)
{
    { // () mesh_transform_vertex_positions_to_double_unit_box
        // TODO overwrite entries of vertex_positions
        vec3 left = V3(HUGE, HUGE, HUGE);
        vec3 right = V3(-HUGE, -HUGE, -HUGE);
        for (int i = 0; i < num_vertices; i++)
        {
            for (int d = 0; d < 3; ++d)
            {
                left[d] = MIN(left[d], vertex_positions[i][d]);
                right[d] = MAX(right[d], vertex_positions[i][d]);
            }
        }
        vec3 center = (left + right) / 2;
        vec3 length = right - left;

        double scaling_factor = 0;

        for (int l = 0; l < 2; l++)
        {
            scaling_factor = length[l];

            if (length[l] < length[l + 1])
            {
                scaling_factor = length[l + 1];
            }
        }

        for (int n = 0; n < num_vertices; n++)
        {
            for (int m = 0; m < 3; m++)
            {
                vertex_positions[n][m] = vertex_positions[n][m] - center[m];
                vertex_positions[n][m] = vertex_positions[n][m] * (2 / scaling_factor);
            }
        }
    }
}

vec3 calc(int3 tri, vec3 *vertex_positions)
{
    vec3 point_x = vertex_positions[tri[0]];
    vec3 point_y = vertex_positions[tri[1]];
    vec3 point_z = vertex_positions[tri[2]];

    vec3 cross_1 = point_y - point_x;
    vec3 cross_2 = point_z - point_x;

    double area = norm(cross(cross_1, cross_2) / 2);
    vec3 normal = cross(cross_1, cross_2) / norm(cross(cross_1, cross_2));

    return area * normal;
}

void fancy_mesh_alloc_compute_and_store_area_weighted_vertex_normals(FancyTriangleMesh3D *fancy_mesh)
{
    ASSERT(fancy_mesh->vertex_normals == NULL);
    { /*Pseudocode for calculating area-weighted vertex normals for an indexed triangle mesh:
# start with a zero vector for every vertex
# e.g. use calloc(num_vertices, sizeof(vec3))


normals = [ V3(0, 0, 0) for _ in range(num_vertices) ]
for tri in triangle_indices:
    area   = calculate_triangle_area(...)
    normal = calculate_triangle_unit_normal(...)
    for vertex_index in tri
        normals[vertex_index] += area * normal
for normal in normals:
    normal = normalized(normal) */
        // TODO allocate fancy_mesh->vertex_normals
        // TODO write entries of fancy_mesh->vertex_normals
        fancy_mesh->vertex_normals = (vec3 *)calloc(fancy_mesh->num_vertices, sizeof(vec3));

        for (int i = 0; i < fancy_mesh->num_triangles; i++)
        {
            vec3 normal_calc = calc(fancy_mesh->triangle_indices[i], fancy_mesh->vertex_positions);

            for (int d = 0; d < 3; d++)
            {
                fancy_mesh->vertex_normals[fancy_mesh->triangle_indices[i][d]] += normal_calc;
            }
        }

        for (int j = 0; j < fancy_mesh->num_vertices; j++)
        {
            fancy_mesh->vertex_normals[j] = normalized(fancy_mesh->vertex_normals[j]);
        }
    }
}

void fancy_mesh_merge_duplicated_vertices(FancyTriangleMesh3D *fancy_mesh)
{
    int new_num_vertices = 0;
    vec3 *new_vertex_positions = (vec3 *)calloc(fancy_mesh->num_vertices, sizeof(vec3)); // (more space than we'll need)
    {                                                                                    // [] fancy_mesh_merge_duplicated_vertices
                                                                                         // TODO set new_num_vertices
                                                                                         // TODO wrie entries of new_vertex_positions
                                                                                         // TODO overwrite entries of fancy_mesh->triangle_indices with new triangle indices
                                                                                         // NOTE it is OK if your implementation is slow (mine takes ~5 seconds to fix up the teapot in debug mode)
                                                                                         // NOTE please don't worry about space efficiency at all
    }
    if (new_num_vertices)
    {
        fancy_mesh->num_vertices = new_num_vertices;
        free(fancy_mesh->vertex_positions);
        fancy_mesh->vertex_positions = new_vertex_positions;
    }
}

void hw3a()
{
    init();

    // no workarounds allowed :)

    BasicTriangleMesh3D basic_box = load_basic_mesh("data_basic_box", true);
    FancyTriangleMesh3D fancy_bunny = load_fancy_mesh("data_fancy_bunny", true, true, false);
    BasicTriangleMesh3D basic_bunny = fancy2basic(fancy_bunny);
    FancyTriangleMesh3D fancy_teapot_with_seams = load_fancy_mesh("data_fancy_teapot_with_seams", true, true, false);
    FancyTriangleMesh3D fancy_teapot_no_seams = load_fancy_mesh("data_fancy_teapot_with_seams", true, true, true);

    Camera3D camera = {5, RAD(45), RAD(0), RAD(-10), 0, .1};
    double t = 0;
    bool paused = false;
    int part = 0;
    while (begin_frame())
    {
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 M = RotationY(.7 * t);
        mat4 PVM = P * V * M;

        int num_parts = 0;
        if (part == num_parts++)
        {
            basic_draw(TRIANGLE_MESH, PVM, basic_box, monokai.green, 3, false, AVG(monokai.green, monokai.white));
        }
        if (part == num_parts++)
        {
            fancy_draw(P, V, M, fancy_bunny, monokai.blue);
        }
        if (part == num_parts++)
        {
            basic_draw(TRIANGLE_MESH, PVM, basic_bunny, monokai.blue, 3, false, AVG(monokai.blue, monokai.white));
        }
        if (part == num_parts++)
        {
            fancy_draw(P, V, M, fancy_teapot_with_seams, monokai.red);
        }
        if (part == num_parts++)
        {
            fancy_draw(P, V, M, fancy_teapot_no_seams, monokai.red);
        }
        imgui_slider("part", &part, 0, num_parts - 1, 'j', 'k', true);

        { // bounding [-1, 1]^3 box
            double tmp[] = {
                -1,
                -1,
                -1,
                -1,
                1,
                -1,
                -1,
                1,
                1,
                -1,
                -1,
                1,
                -1,
                -1,
                -1,
                1,
                -1,
                -1,
                1,
                1,
                -1,
                1,
                1,
                1,
                1,
                -1,
                1,
                1,
                -1,
                -1,
                1,
                1,
                -1,
                -1,
                1,
                -1,
                1,
                1,
                -1,
                1,
                1,
                1,
                -1,
                1,
                1,
                1,
                1,
                1,
                1,
                -1,
                1,
                -1,
                -1,
                1,
                1,
                -1,
                1,
            };
            basic_draw(LINE_STRIP, PVM.data, XYZ, RGB, NELEMS(tmp) / 3, tmp, NULL, 1, 1, 1, .1, 2);
        }

        imgui_checkbox("paused", &paused, 'p');
        if (!paused)
        {
            t += .0167;
        }
        if (imgui_button("reset", 'r'))
        {
            t = 0;
        }
    }
}
FancyTriangleMesh3D hw1()
{
    FancyTriangleMesh3D fancy_mesh = {};
    {

        StretchyBuffer<vec3> vertex_positions = {};
        StretchyBuffer<int3> triange_indices = {};

        int N = 64;
        int n = 16;
        for (int i = 0; i < N; i++)
        {
            double theta = (double(i) / (N - 1)) * (2 * PI);
            for (int j = 0; j < n; ++j)
            {
                double phi = (double(j) / (n - 1)) * (2 * PI);
                vec3 p = V3(3, 0, 0) + V3(cos(phi), sin(phi), 0);
                p = transformPoint(RotationY(theta), p);
                sbuff_push_back(&vertex_positions, p);
                int3 tri_a = {
                    i * n + j,
                    ((i + 1) % N) * n + j,
                    ((i + 1) % N) * n + ((j + 1) % n)};

                int3 tri_b = {
                    i * n + j,
                    ((i + 1) % N) * n + ((j + 1) % n),
                    i * n + ((j + 1) % n)};
                sbuff_push_back(&triange_indices, tri_a);
                sbuff_push_back(&triange_indices, tri_b);
            }
        }
    }
    return fancy_mesh;
}

void hw3b()
{
    init();

    mat4 M[10] = {};
    {
        for (int i = 0; i < NELEMS(M); ++i)
        {
            M[i] = Identity4x4;
        }
    }

    vec2 L_quads[] = {{0, 0}, {2, 0}, {2, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 3}, {0, 3}};

    Camera2D camera = {20};
    while (begin_frame())
    {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        // // 2D transform API
        // mat4 Translation(double t_x, double t_y); // translation by (t_x, t_y)^T
        // mat4 Scaling(double s_x, double s_y) // scaling by s_x in x and x_y in y
        // mat4 Rotation(double theta); // counter-clockwise rotation about origin by theta radians
        // - NOTE RAD(deg) converts from degrees to radians

#if 0
        M[0] = Translation(5, 5);
#endif

        // () yellow
        M[0] = Translation(2, -1);

        // () purple
        M[1] = Rotation(3.14159 / 2);

        // () orange
        M[2] = Scaling(1, -1);

        // () lightblue
        M[3] = Rotation(3.14159 / 2) * Scaling(-1, 1);

        // () red
        M[4] = Rotation(3.14159 * 3 / 2) * Translation(-3, 1);

        // () buff (tan)
        // M[5] = Translation(7, -3) * Scaling(-3, 2);
        M[5] = Scaling(3, 2);

        // () gray
        M[6] = Rotation(3.14159 * 3 / 2) * Translation(-3, -4) * Scaling(2, 1);

        // [] green
        M[7] = Rotation(3.14159) * Translation(-4, -2) * Scaling(0.5, 0.5);

        // [] purplishpink
        // NOTE only a perfect solution will score full credit
        //      i.e., no hard-coded constants with 10 digits after the decimal place :)
        M[8] = Translation(-4, -1) * Rotation(3.14159 * 3 / 4) * inverse(Translation(-4, -1)) * Rotation(3.14159 * 3 / 2) * Translation(-3, -4) * Scaling(2, 1);

        // <> blue
        // NOTE only a perfect solution will score full credit
        //      i.e., no hard-coded constants with 10 digits after the decimal place :)
        M[9] = {};

        { // draw L-block's
            basic_draw(QUADS, PV, NELEMS(L_quads), L_quads, monokai.white);
            for (int i = 0; i < NELEMS(M); ++i)
            {
                basic_draw(QUADS, PV * M[i], NELEMS(L_quads), L_quads, color_get_kelly(i));
            }
        }
        { // bespoke widget
            vec2 mouse_position = input_get_mouse_position_in_world_coordinates(PV);
            int x = (int)roundf((float)mouse_position.x);
            int y = (int)roundf((float)mouse_position.y);
            imgui_readout("x", &x);
            imgui_readout("y", &y);
            gl_PV(PV);
            gl_color(0, 1, 0, 1);
            gl_begin(GL_POINTS);
            gl_vertex(x, y);
            gl_end();
            gl_begin(GL_LINE_LOOP);
            gl_vertex(x, y);
            gl_color(0, 1, 0, .5);
            gl_vertex(0, y);
            gl_vertex(0, 0);
            gl_vertex(x, 0);
            gl_end();
        }

        // NOTE if you want to draw other stuff to help you debug, do it down here
    }
}

// // 3D transform API
// mat4 Translation(double t_x, double t_y, double t_z);
// mat4 Scaling(double s_x, double s_y, double s_z);
// mat4 RotationX(double theta);
// mat4 RotationY(double theta);
// mat4 RotationZ(double theta);
// mat4 Rotation(vec3 axis, dbouel theta); // probably not so useful?

void hw3c()
{
    init();
    Camera3D camera = {10, RAD(45)};
    double t = 0;
    bool playing = false;
    while (begin_frame())
    {
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);

        mat4 S = Scaling(1.25 - .25 * cos(5 * t), .722 + .278 * cos(5 * t), 1.25 - .25 * cos(5 * t));
        basic_draw(P * V, meshlib.basic_axes);

        fancy_draw(P, V, Translation(4.5, -1, -1) * S * Scaling(0.5, 0.5, 0.5), meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(3.5, -1.5, -0.5) * S * Scaling(0.9, 0.9, 0.9), meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(3.5, -2, 0) * S, meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(2.5, -1, 0) * S * Scaling(0.8, 0.8, 0.8), meshlib.fancy_sphere, monokai.white);

        fancy_draw(P, V, Translation(-4.5, -1, -1) * S * Scaling(0.5, 0.5, 0.5), meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(-3.5, -1.5, -0.5) * S * Scaling(0.9, 0.9, 0.9), meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(-3.5, -2, 0) * S, meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(-2.5, -1, 0) * S * Scaling(0.8, 0.8, 0.8), meshlib.fancy_sphere, monokai.white);

        fancy_draw(P, V, Translation(-3.5, 1, -1) * S * Scaling(1, 0.5, 1), meshlib.fancy_cylinder, monokai.red);
        fancy_draw(P, V, Translation(3.5, 1, -1) * S * Scaling(1, 0.5, 1), meshlib.fancy_cylinder, monokai.red);
        fancy_draw(P, V, Translation(-3.5, 0, -1) * S * Scaling(1, 0.5, 1), meshlib.fancy_cylinder, monokai.orange);
        fancy_draw(P, V, Translation(3.5, 0, -1) * S * Scaling(1, 0.5, 1), meshlib.fancy_cylinder, monokai.orange);
        fancy_draw(P, V, Translation(-3.5, -1, -1) * S * Scaling(0.8, 0.8, 0.8), meshlib.fancy_box, monokai.gray);
        fancy_draw(P, V, Translation(3.5, -1, -1) * S * Scaling(0.8, 0.8, 0.8), meshlib.fancy_box, monokai.gray);
        fancy_draw(P, V, Translation(-3.5, 1, -1) * S * Scaling(0.8, 0.8, 0.8), meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(-3.5, -0.5, -1) * S * Scaling(0.5, 0.5, 0.5), meshlib.fancy_cylinder, monokai.white);
        fancy_draw(P, V, Translation(3.5, -0.5, -1) * S * Scaling(0.5, 0.5, 0.5), meshlib.fancy_cylinder, monokai.white);
        fancy_draw(P, V, Translation(3.5, 1, -1) * S * Scaling(0.8, 0.8, 0.8), meshlib.fancy_sphere, monokai.white);

        fancy_draw(P, V, Translation(-3.5, 1.5, -1) * S, meshlib.fancy_cone, monokai.blue);
        fancy_draw(P, V, Translation(3.5, 1.5, -1) * S, meshlib.fancy_cone, monokai.blue);

        fancy_draw(P, V, Translation(0, -1, -1) * S * Scaling(2.5, 2.5, 2.5), meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(0, -1, -1) * S * Scaling(1.8, 1.8, 1.8), meshlib.fancy_box, monokai.gray);
        fancy_draw(P, V, Translation(0, 1, -1) * S * Scaling(1.2, 0.8, 1.2), meshlib.fancy_cylinder, monokai.orange);

        fancy_draw(P, V, Translation(0, 2.7, -1) * S * Scaling(1.2, 1.2, 1.2), meshlib.fancy_sphere, monokai.white);
        fancy_draw(P, V, Translation(0, 0, -1) * S * Scaling(1, 1.5, 1), meshlib.fancy_cylinder, monokai.white);
        // fancy_draw(P, V, Translation(0, 2, -1) * S * Scaling(1.5, 0.2, 1.5), meshlib.fancy_cylinder, monokai.purple);

        fancy_draw(P, V, Translation(0, 2.2, -1) * S * Scaling(1.5, 0.5, 1.5), meshlib.fancy_cylinder, monokai.red);
        fancy_draw(P, V, Translation(0, 3, -1) * S * Scaling(1.5, 1.5, 1.5), meshlib.fancy_cone, monokai.blue);

        // fancy_draw(P, V, Translation(-1.5, 0, 0) * S, meshlib.fancy_cone, monokai.yellow);
        // fancy_draw(P, V, Translation(1.5, 0, 0) * S, meshlib.fancy_cylinder, monokai.blue);
        // fancy_draw(P, V, Translation(4.5, 0, 0) * S, meshlib.fancy_sphere, monokai.purple);

        { // floor
            gl_PV(P * V);
            gl_begin(QUADS);
            gl_color(1, 1, 1, .5);
            gl_vertex(10, 0, 10);
            gl_vertex(-10, 0, 10);
            gl_vertex(-10, 0, -10);
            gl_vertex(10, 0, -10);
            gl_end();
        }

        imgui_checkbox("playing", &playing, 'p');
        if (playing)
        {
            t += .0167;
        }
    }
}

void hw()
{
    init(true, "", 540, 0, 100);
    hw1();
    hw3a();
    hw3b();
    hw3c();
}

// end submission

int main()
{
    hw();
    return 0;
}
