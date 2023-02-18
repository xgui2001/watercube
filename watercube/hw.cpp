#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_cow_supplement.cpp"
#include <vector>
#include <iostream>

double random_position()
{
    return rand() % 100;
}

void hw6c_draw_textured_square(mat4 P, mat4 V, mat4 M, char *texture_filename)
{
    // please ignore; this function is hack-a-saurus rex
    static FancyTriangleMesh3D square;
    if (!square.num_vertices)
    {
        square = meshlib.fancy_square;
        square.vertex_normals = NULL;
    }
    fancy_draw(P, V, M,
               square.num_triangles, square.triangle_indices, square.num_vertices, square.vertex_positions,
               NULL, NULL, {},
               square.vertex_texCoords, texture_filename);
};

void world()
{
    init();

    // vec2 point[i] = {0, 1};
    // vec2 point2 = {0, 2};

    const int NUM_POINTS = 300;
    vec2 point[NUM_POINTS] = {};
    for (int i = 0; i < NUM_POINTS; i++)
    {
        point[i] = {util_random_double(-2, 2), util_random_double(-2, 2)};
    }

    double gravity_constant = -9.81;
    double h = 1. / 90.;
    double t = 0;

    // vec2 floor[2] = {{-5, -5}, {5, -5}};
    // vec2 wall_left[2] = {{-5, 0}, {-5, 10}};
    // vec2 wall_right[2] = {{5, 0}, {5, 10}};
    static vec2 speed[NUM_POINTS] = {};

    vec2 gravity = gravity_constant * V2(cos(-RAD(90)), sin(-RAD(90)));

    // texture

    int side_length_in_pixels = 16;
    int nrChannels = 4;                                                                                        // r0 g0 b0 a0 r1 g1 b1 a1 ...
    unsigned char *data = (unsigned char *)malloc(side_length_in_pixels * side_length_in_pixels * nrChannels); // the amount of data it takes to represent everything
    memset(data, 255, side_length_in_pixels * side_length_in_pixels * nrChannels);
    fancy_texture_create("my_custom_texture", side_length_in_pixels, side_length_in_pixels, nrChannels, data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); // HACK so we can *see* the pixels
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); // ""

    int3 triangle_indices[] = {
        {0, 1, 2},
        {0, 2, 3}};
    vec3 vertex_positions[] = {
        {-10, -10, -1},
        {10, -10, -1},
        {10, 10, -1},
        {-10, 10, -1},
    };
    vec2 vertex_texCoords[] = {
        {0, 0},
        {1, 0},
        {1, 1},
        {0, 1}};

    StretchyBuffer<vec2> drawings = {};

    while (begin_frame())
    {
        t += 1. / 60;
        static Camera3D camera = {5, RAD(75)};
        // camera_move(&camera);

        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);

        mat4 PV = camera_get_PV(&camera);
        gl_PV(PV);
        mat4 M_new = Translation(0, 0) * Scaling(1, 1);
        {
            if (input.key_held['a'] == true)
            {
                M_new *= Translation(-1, 0);
            }

            if (input.key_held['w'] == true)
            {
                M_new *= Translation(0, 1);
            }

            if (input.key_held['d'] == true)
            {
                M_new *= Translation(1, 0);
            }

            if (input.key_held['s'] == true)
            {
                M_new *= Translation(0, -1);
            }
        }
        fancy_draw(P, V, Identity4x4,
                   2, triangle_indices, 4, vertex_positions,
                   NULL, NULL, {},
                   vertex_texCoords, "underwater.jpeg");

        { // update texture

            int pixelIndex = 0;
            for (int i = 0; i < side_length_in_pixels; ++i)
            {
                for (int j = 0; j < side_length_in_pixels; ++j)
                {
                    vec2 UV_cords = V2(double(i) / (side_length_in_pixels - 1), double(j) / (side_length_in_pixels - 1));
                    data[4 * pixelIndex + 0] = (0.5 + 0.5 * cos(t + UV_cords.x + 0)) * 255 / 2; // r
                    data[4 * pixelIndex + 1] = (0.5 + 0.5 * cos(t + UV_cords.y + 2)) * 255 / 2; // g
                    data[4 * pixelIndex + 2] = (0.5 + 0.5 * cos(t + UV_cords.x + 4)) * 255;     // b
                    data[4 * pixelIndex + 3] = 255 / 8;                                         // brightness
                    pixelIndex++;
                }
            }
            fancy_texture_update("my_custom_texture", side_length_in_pixels, side_length_in_pixels, nrChannels, data);
        }
        fancy_draw(P, V, Identity4x4,
                   2, triangle_indices,
                   4, vertex_positions,
                   NULL, NULL, {},
                   vertex_texCoords, "my_custom_texture");
        hw6c_draw_textured_square(PV, Identity4x4, M_new, "pixil-frame-0-removebg-preview.png");

        vec2 s_mouse = input_get_mouse_position_in_world_coordinates(PV);

        std::vector<int> primitive_types = {POINTS, LINES, LINE_STRIP, LINE_LOOP, TRIANGLES, TRIANGLE_STRIP, TRIANGLE_FAN, QUADS, TRIANGLE_MESH, QUAD_MESH};
        static int primitive_num = 0;
        imgui_slider("primitive_type", &primitive_num, 0, 9, 'j', 'k', true);

        static double size = 10.0;
        imgui_slider("size", &size, 1.0, 30.0);

        static double radius = 2.5;
        imgui_slider("radius", &radius, 2.0, 4.0);

        static int num = 200;
        imgui_slider("points_to_draw", &num, 1, 300);

        double pos = 0.0 + radius;
        double neg = 0.0 - radius;

        vec2 box[4] = {{pos, pos}, {pos, neg}, {neg, neg}, {neg, pos}};

        static vec2 particle_normalized[NUM_POINTS][NUM_POINTS] = {};
        double particle_weight[NUM_POINTS] = {};
        double particle_pressure[NUM_POINTS] = {};

        // drawings CURRENTLY ONLY WORK WHILE CAMERA_MOVE IS ENABLED
        /*if (input.mouse_left_pressed)
        {
            sbuff_push_back(&drawings, {});
        }

        if (input.mouse_left_held)
        {
            if ((-s_mouse.x <= radius) && (-s_mouse.x >= -radius) && (-s_mouse.y <= radius) && (-s_mouse.y >= -radius))
            {
                sbuff_push_back(&drawings, s_mouse);
            }
        }

        if (input.key_pressed['x'] == true)
        {
            sbuff_free(&drawings);
        }
        */

        // update all points
        for (int i = 0; i < NUM_POINTS; i++)
        {

            // determine speed
            speed[i] -= h * gravity / 100.;
            for (int j = 0; j < NUM_POINTS; j++)
            {
                if (i != j) // cannot be while (i != j)
                {
                    double D = norm(point[i] - point[j]);
                    if (D < 1 && !IS_ZERO(D))
                    {
                        // 1 and 2 are colliding particles
                        // n[0][1] is normalized relative position
                        particle_normalized[i][j] = normalized(point[i] - point[j]);
                        particle_normalized[j][i] = normalized(point[j] - point[i]);

                        // vec2 particle_normalized0 = normalized(point1 - point2);
                        // vec2 particle_normalized1 = normalized(point2 - point1);

                        // w[0][1] is weight
                        particle_weight[i] = (1 - (norm(point[i] - point[j]) / 1.0));
                        particle_weight[j] = (1 - (norm(point[j] - point[i]) / 1.0));

                        // double particle_weight0 = (1 - (norm(point1 - point2) / 1.0));
                        // double particle_weight1 = (1 - (norm(point2 - point1]) / 1.0));

                        // particle weight sum is the sum of particle weights around that particle
                        // particle_weight_sum[0] = particle_weight[0][1];
                        // particle_weight_sum[1] = particle_weight[1][0];

                        // double particle_weight_avg = (particle_weight_sum[0] + particle_weight_sum[1]) / 2.0;
                        particle_pressure[i] = fmax(0, 1 * (particle_weight[i] - 0.4));
                        particle_pressure[j] = fmax(0, 1 * (particle_weight[j] - 0.4));

                        // double particle_pressure0 = fmax(0, 2 * (particle_weight0 - 0.5));
                        // double particle_pressure1 = fmax(0, 2 * (particle_weight1 - 0.5));

                        speed[i] += h * 0.2 * (particle_pressure[i] + particle_pressure[j]) * particle_weight[i] * particle_normalized[i][j];
                        speed[j] += h * 0.2 * (particle_pressure[j] + particle_pressure[i]) * particle_weight[j] * particle_normalized[j][i];
                        //   speed[0] += h * 2 * (particle_pressure0 + particle_pressure1) * particle_weight0 * particle_normalized0;
                        //   speed[1] += h * 2 * (particle_pressure1 + particle_pressure0) * particle_weight1 * particle_normalized1;

                        // viscosity
                        speed[i] += h * 0.0001 * (speed[i] - speed[j]);
                        speed[j] += h * 0.0001 * (speed[j] - speed[i]);

                        // tensile WIP
                    }
                }
            }

            /*
            { // drawing
                if (drawings.length > 1)
                {
                    for (int l = 0; l < drawings.length; l++)
                    {
                        if (norm(drawings.data[l] - point[i]) < 0.1 && norm(drawings.data[l] - point[i]) != 0)
                        {
                            speed[i].x *= 0;
                            speed[i].y *= 0;
                        }
                    }
                }
            }
            */

            // move point
            point[i] += speed[i];

            // enforce bounds
            {
                if (point[i].y <= -radius)
                {
                    speed[i].y *= -.8;
                    point[i].y = -radius + TINY;
                }

                if (point[i].y >= radius)
                {
                    speed[i].y *= -.8;
                    point[i].y = radius - TINY;
                }

                if (point[i].x <= -radius)
                {
                    speed[i].x *= -.8;
                    point[i].x = -radius + TINY;
                }

                if (point[i].x >= radius)
                {
                    speed[i].x *= -.8;
                    point[i].x = radius - TINY;
                }

                // duck
                vec2 center = {0, 0};

                if (input.key_held['a'] == true)
                {
                    center.x = center.x - 1;
                }

                if (input.key_held['w'] == true)
                {
                    center.y = center.y + 1;
                }

                if (input.key_held['d'] == true)
                {
                    center.x = center.x + 1;
                }

                if (input.key_held['s'] == true)
                {
                    center.y = center.y - 1;
                }

                if (norm(point[i] - center) < 0.3 && norm(point[i] - center) != 0)
                {
                    speed[i].x *= -0.8;
                    speed[i].y *= -0.8;
                    if (point[i].x < center.x)
                    {
                        point[i].x -= 0.03;
                    }
                    if (point[i].x > center.x)
                    {
                        point[i].x += 0.03;
                    }
                    if (point[i].y > center.y)
                    {
                        point[i].y += 0.03;
                    }
                    if (point[i].y < center.y)
                    {
                        point[i].y -= 0.03;
                    }
                }
            }
        }

        // draw everything
        basic_draw(LINE_LOOP, PV, 4, box, monokai.blue);
        basic_draw(primitive_types[primitive_num], PV, num, point, COS_LERP(t, monokai.blue, monokai.blue2), size);
        // basic_draw(LINE_STRIP, PV, drawings.length, drawings.data, color_rainbow_swirl(.1 * t));
    }
}

int main()
{
    world();
    return 0;
}