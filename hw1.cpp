// !! This is the hw.cpp for hw1.
//
// If you cloned today expecting to start working on hw0, please...
//
// - open the file explorer / finder
// - navigate into folder old_hw\hw0
// - copy hw.cpp
// - navigate back to the root folder (CSCI-371)
// - paste hw.cpp inside, replacing the file that's already there (this file)
//
// You should then be good to go :)
//
// Please ask me ASAP if you have any questions.

// if you find a bug in the hw, please report it on Github Issues or email
//                                        the class will be rewarded with snacks

#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include <vector>
#include <iostream>

//
// documentation
//

#if 0
// usage code
#endif

// press the \ key to display an fps counter
// (macbook users see README.txt if you don't have a 60 fps cap)
// press the / key to uncap the framerate

// here is the actual basic_draw function signature (at least for this week)
// if vertex_colors is NULL, then it draws all vertices with the same color,
// specifically { r_fallback, g_fallback, b_fallback, a_fallback }
// if overlay is true, then we draw on top of _everything_
//
// void basic_draw(
//         int primitive,
//         double *transform,
//         int dimension_of_positions,
//         int dimension_of_colors,
//         int num_vertices,
//         double *vertex_positions,
//         double *vertex_colors = NULL,
//         double r_fallback = 1,
//         double g_fallback = 1,
//         double b_fallback = 1,
//         double a_fallback = 1,
//         double size_in_pixels = 0,
//         bool overlay = false,
//         double r_wireframe = 1,
//         double g_wireframe = 1,
//         double b_wireframe = 1,
//         double a_wireframe = 1);
//
// snail gives us a couple convenient wrappers, which we can think of as
//
//   // single color
//   void basic_draw(
//           int primtive,
//           mat4 transform,
//           int num_vertices,
//           vecX *vertex_positions,
//           vec3 color,
//           double size_in_pixels = 0,
//           bool overlay = false);
//
//   // per-vertex color
//   void basic_draw(
//           int primtive,
//           mat4 transform,
//           int num_vertices,
//           vecX *vertex_positions,
//           vec3 *vertex_colors,
//           double size_in_pixels = 0,
//           bool overlay = false);
//
// here vecX means you can use vec2 or vec3 (cow.cpp uses templates for this)

#if 0
static Camera2D camera = { 5 };
mat4 PV = camera_get_PV(&camera);
vec2 foo[3] = { { 0, 0 }, { 1, 0}, { 0, 1} };
basic_draw(TRIANGLES, PV, 3, foo, monokai.white);
#endif

// this function lets you drag 2D vertices
// the arguments with default arguments specify the size and color of the dot
// that pops up when you hover over a point
//
// bool widget_drag(
//         double *PV,
//         int num_vertices,
//         double *vertex_positions,
//         double size_in_pixels = 0,
//         double r = 1,
//         double g = 1,
//         double b = 1,
//         double a = 1);
//
// if we include snail we can call this wrapper if we prefer
//
// bool widget_drag(
//         mat4 PV,
//         int num_vertices,
//         vec2 *vertex_positions,
//         double size_in_pixels = 0,
//         vec3 color = monokai.white);

#if 0
static Camera2D camera = { 5 };
mat4 PV = camera_get_PV(&camera);
static vec2 foo[3] = { { 0, 0 }, { 1, 0 }, { 0, 1 } };
basic_draw(LINE_LOOP, PV, 3, foo, monokai.white);
widget_drag(PV, 3, foo);
#endif

// imgui_slider lets us scrub an int t between bounds a and b.
// pressing the key j will decrement t, pressing the key k will increment it
// if loop is true, e.g. going above b will take us back to a
//
// void imgui_slider(
//         char *name,
//         int *t,
//         int a,
//         int b,
//         char j = 0,
//         char k = 0,
//         bool loop = false);
//
// there is also a version for doubles
//
// void imgui_slider(char *name, double *t, double a, double b);

#if 0
static int foo = 0;
imgui_slider("foo", &foo, 0, 100, 'j', 'k', true);
#endif

// here is how we can access user input
//
// input.key_pressed[...]
// input.key_released[...]
// input.key_held[...]
// input.key_toggle[...]
// input.mouse_left_pressed
// input.mouse_left_held
// input.mouse_left_released
// input.mouse_right_pressed
// input.mouse_right_held
// input.mouse_right_released
//
// void input_get_mouse_position_and_change_in_position_in_world_coordinates(
//         double *PV,
//         double *mouse_x_world,
//         double *mouse_y_world,
//         double *mouse_dx_world = NULL,
//         double *mouse_dy_world = NULL);
//
// which has these snail wrappers
//
// vec2 input_get_mouse_position_in_world_coordinates(mat4 PV);
// vec2 input_get_mouse_change_in_position_in_world_coordinates(mat4 PV);

#if 0
if (input.key_pressed['a']) { printf("you pressed the A key!\n"); }

static Camera2D camera = { 5 };
camera_move(&camera);
mat4 PV = camera_get_PV(&camera);
if (input.mouse_left_pressed) {
    vec2 s_mouse = input_get_mouse_position_in_world_coordinates(PV);
    printf("you lefted clicked at (%lf, %lf) in world coordinates\n", s_mouse.x, s_mouse.y);
}
#endif

// NELEMS(fixed_size_array) gives the number of elements in a fixed-size array
// use at your own risk; careful not to call it on a pointer

//
// implementation
//

// begin submission

void hw1a()
{

    while (begin_frame())
    {
        static Camera2D camera = {5};
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        // draggable red square outline
        // static array contains the
        static vec2 square[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        basic_draw(LINE_LOOP, PV, 4, square, monokai.red);
        widget_drag(PV, 4, square);

        // green circle outline with slider to change number of vertices
        gl_PV(PV);
        gl_begin(LINE_LOOP);
        gl_color(166. / 255, 226. / 255, 46. / 255);

        // slider
        static int N = 64;
        imgui_slider("N", &N, 0, 64, 'j', 'k', true);
        // i is the number of times the loop runs and should go up to N (number of vertex)
        for (double i = 0; i <= N; ++i)
        {
            // gl_vertex((1.5 * sqrt(2) * cos((i / N) * 6.28)), 1.5 * sqrt(2) * sin((i / N) * 6.28));

            double t = double(i) / (N - 1);
            double theta = double(i) / N * 2 * PI;

            gl_vertex((LERP(t, 2, 1) * cos(theta)), (LERP(t, 2, 1) * sin(theta)));
        }

        gl_end();
    }
}

double random_number()
{
    return 2 * (rand() / (double)RAND_MAX) - 1;
}

double random_color()
{
    return rand() % 256;
}

double random_velocity()
{
    return (rand() / (double)RAND_MAX) / 50;
}

void hw1b()
{

    vec3 *vertex_positions = (vec3 *)malloc(100000 * sizeof(vec3));
    for (int p = 0; p <= 100000; p++)
    {
        vertex_positions[p] = {random_number(), random_number(), random_number()};
    }

    vec3 *vertex_colors = (vec3 *)malloc(100000 * sizeof(vec3));
    for (int c = 0; c <= 100000; c++)
    {
        vertex_colors[c] = {random_color() / 255, random_color() / 255, random_color() / 255};
    }

    vec3 *vertex_velocities = (vec3 *)malloc(100000 * sizeof(vec3));
    for (int v = 0; v <= 100000; v++)
    {
        vertex_velocities[v] = {random_velocity(), random_velocity(), random_velocity()};
    }

    while (begin_frame())
    {
        static Camera3D camera = {5, RAD(45)};
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        gl_PV(PV);

        bool playing = false;
        double t = 0;
        double a = 0;
        double b = 2;

        static vec3 cube[8] = {{-1, 1, 1}, {-1, 1, -1}, {-1, -1, -1}, {-1, -1, 1}, {1, -1, 1}, {1, -1, -1}, {1, 1, -1}, {1, 1, 1}};
        basic_draw(LINE_LOOP, PV, 8, cube, monokai.white);
        static vec3 cubefill_horizontal[4] = {{-1, 1, -1}, {1, 1, -1}, {-1, -1, -1}, {1, -1, -1}};
        basic_draw(LINES, PV, 4, cubefill_horizontal, monokai.white);
        static vec3 cubefill_vertical[4] = {{1, 1, 1}, {1, -1, 1}, {-1, 1, 1}, {-1, -1, 1}};
        basic_draw(LINES, PV, 4, cubefill_vertical, monokai.white);

        static int num = 12;
        imgui_slider("num_vertices_to_draw", &num, 4, 100000);

        std::vector<int> primitive_types = {POINTS, LINES, LINE_STRIP, LINE_LOOP, TRIANGLES, TRIANGLE_STRIP, TRIANGLE_FAN, QUADS, TRIANGLE_MESH, QUAD_MESH};
        static int primitive_num = 0;
        imgui_slider("i_primitive", &primitive_num, 0, 9, 'j', 'k', true);

        static double size = 0.000000;
        imgui_slider("size", &size, 0.0, 50.0);

        for (int n = 0; n <= 100000; n++)
        {
            if ((vertex_positions[n].data[0] >= 1) || (vertex_positions[n].data[0] <= -1))
            {
                vertex_velocities[n].data[0] *= -1;
            }
            if ((vertex_positions[n].data[1] >= 1) || (vertex_positions[n].data[1] <= -1))
            {
                vertex_velocities[n].data[1] *= -1;
            }
            if ((vertex_positions[n].data[2] >= 1) || (vertex_positions[n].data[2] <= -1))
            {
                vertex_velocities[n].data[2] *= -1;
            }
            vertex_positions[n].data[0] += vertex_velocities[n].data[0];
            vertex_positions[n].data[1] += vertex_velocities[n].data[1];
            vertex_positions[n].data[2] += vertex_velocities[n].data[2];
        }
        basic_draw(primitive_types[primitive_num], PV, num, vertex_positions, vertex_colors, size);
    }
}

struct StretchyBuffer
{
    int length;   // the amount of stuff the buffer currently holds
    int capacity; // the amount of stuff the buffer currently CAN hold
    vec2 *data;   // stretchybuffer's elements are of vec2 type
};

void sbuff_push_back(StretchyBuffer *buffer, vec2 point)
{
    if (buffer->length == 0)
    {
        buffer->data = (vec2 *)malloc(16 * sizeof(vec2));
        buffer->capacity = 16;
    }

    if (buffer->length == buffer->capacity)
    {
        buffer->capacity = 2 * buffer->capacity;
        vec2 *new_buffer = (vec2 *)realloc(buffer->data, buffer->capacity * sizeof(vec2));
        buffer->data = new_buffer;
    }
    buffer->data[buffer->length] = point;
    buffer->length++;
}

void sbuff_free(StretchyBuffer *buffer)
{
    free(buffer->data);
    buffer->data = NULL;
    buffer->length = 0;
    buffer->capacity = 0;
}

void hw1c()
{

#if 1
    StretchyBuffer buffer = {};

    ASSERT(buffer.length == 0);
    ASSERT(buffer.capacity == 0);

    sbuff_push_back(&buffer, {});

    ASSERT(buffer.length == 1);
    ASSERT(buffer.capacity == 16);

    int N = 2048;
    for (int i = 0; i < N; ++i)
    {
        double theta = double(i) / 64 * 2 * PI;
        double r = double(i) / N;
        sbuff_push_back(&buffer, r * V2(cos(theta), sin(theta)));
    }

    while (begin_frame())
    {
        static Camera2D camera = {.8};
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        static double time = 0;
        time += .0167;

        basic_draw(LINE_STRIP, PV * RotationZ(-10 * time), buffer.length, buffer.data, color_rainbow_swirl(.1 * time));
    }

    // note: we don't technically have to free here if the program is going
    // to exit anyway; i just needed to call free somewhere to test it :)
    sbuff_free(&buffer);
    ASSERT(buffer.length == 0);
    ASSERT(buffer.capacity == 0);
    ASSERT(buffer.data == NULL);
#endif
}

void hw1d()
{
    std::vector<std::vector<vec2>> strokes;

    while (begin_frame())
    {
        static Camera2D camera = {5};
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        gl_PV(PV);

        static double time = 0;
        time += .1;

        if (input.mouse_left_pressed)
        {
            strokes.push_back(std::vector<vec2>());
        }
        if (input.mouse_left_held)
        {
            strokes.back().push_back(input_get_mouse_position_in_world_coordinates(PV));
        }

        for (int i = 0; i < int(strokes.size()); i++)
        {
            basic_draw(LINE_STRIP, PV, strokes[i].size(), strokes[i].data(), color_rainbow_swirl(.1 * time));
        }

        if (input.key_pressed['x'] == true)
        {
            strokes.clear();
        }
    }
}

void red_spiral()
{
    init();
    Camera2D camera = {5};
    mat4 PV = camera_get_PV(&camera);
    gl_PV(PV);
    while (begin_frame(1, 1, 1, 1))
    {
        gl_begin(LINES);
        gl_color(0, 0, 0);
        gl_vertex(1, 0);
        gl_vertex(2, 0);
        gl_end();
    }

    gl_begin(LINE_LOOP);
    gl_color(1, 0, 0);
    int N = 64;
    for (int i = 0; i < N; ++i)
    {
        double t = double(i) / (N - 1);
        double theta = double(i) / 64 * 2 * PI;
        gl_vertex((LERP(t, 1, 2) * sqrt(2) * cos((i / N) * 6.28)), (LERP(t, 1, 2) * sqrt(2) * sin((i / N) * 6.28)));
    }
    gl_end();
}

void hw0()
{
    Camera2D camera = {5};
    while (begin_frame())
    {
        camera_move(&camera);
        gl_PV(camera_get_PV(&camera));
        gl_begin(LINE_LOOP);
        // red square
        gl_color(249. / 255, 38. / 255, 114. / 255); // cow.cpp 382
        gl_vertex(-1.5, -1.5);
        gl_vertex(-1.5, 1.5);
        gl_vertex(1.5, 1.5);
        gl_vertex(1.5, -1.5);
        gl_end();
        // green circle
        gl_begin(LINE_LOOP);
        gl_color(166. / 255, 226. / 255, 46. / 255);
        // approximating each points along the circle using trig
        for (double angle = 0; angle < 360; ++angle)
        {
            gl_vertex(1.5 * sqrt(2) * cos(RAD(angle)), 1.5 * sqrt(2) * sin(RAD(angle)));
        }
        gl_end();
    }
}

void extracredit_v()
{
    static Camera3D camera = {5, RAD(45)};
    camera_move(&camera);
    mat4 PV = camera_get_PV(&camera);
    gl_PV(PV);

    while (begin_frame())
    {
        static vec3 big_grid[8] = {{-1, 1, -10}, {-1, 1, -10}, {-1, -1, -10}, {-1, -1, 10}, {1, -1, 10}, {1, -1, -10}, {1, 1, -10}, {1, 1, 10}};
        basic_draw(LINES, PV, 8, big_grid, monokai.white);
    }
}

void hw()
{
    hw1a();
    hw1b();
    hw1c();
    hw1d();
    red_spiral();
    hw0();
    extracredit_v();

    // you can add additional apps to demo extra credit here
}

// end submission

int main()
{
    init(false, "hw1 :D");
    hw();
    return 0;
}
