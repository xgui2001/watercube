#define COW_CRASH_ON_FLOATING_POINT_EXCEPTIONS

#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"

// // documentation ////////////////////////////////////////////////////////////
//
// // random double in interval [a, b]
// double util_random_double(double a = 0, double b = 1);
//
// double norm(vecX v); // length of a vector v
//
// matX inverse(matX M); // inverse of matrix M
//
// // useful macros (see hw2 writeup for explanation)
// #define LERP(t, a, b)          ((1 - (t)) * (a) + (t) * (b))
// #define INVERSE_LERP(c, a, b)  (((c) - (a)) / double((b) - (a)))
// #define CLAMP(t, a, b)         MIN(MAX(t, a), b)
// #define CLAMPED_LERP(t, a, b)  LERP(CLAMP(t, 0, 1), a, b)
// #define COS_LERP(t, a, b)      LERP(.5 - .5 * cos((t)*PI), a, b)

// begin submission

void hw2a()
{
    init();
    mat3 M_2a = {3, 0, -1, 1, -1, 0, 0, -1, 1};
    vec3 V_2a = {0, -1, 1};
    vec3 p_answer = inverse(M_2a) * V_2a;
    pprint(p_answer);

    vec3 q = {-0.282550, 3.282550, 5.565100};
    double q_answer = norm(q - p_answer);
    printf("%lf\n", q_answer);
}

void hw2b()
{
    init();
    Camera2D camera = {10};
    vec2 circle_center = V2(0, 0);
    double circle_radius = 2;
    vec2 test_point = V2(0, 0);
    while (begin_frame())
    {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        imgui_slider("r", &circle_radius, 0, 5);
        imgui_slider("x", &circle_center.x, -5, 5);
        imgui_slider("y", &circle_center.y, -5, 5);

        bool inside = false;

        if (norm(test_point - circle_center) < circle_radius)
        {
            inside = true;
        }

        basic_draw(POINTS, PV, 1, &test_point, (inside) ? monokai.green : monokai.red);
        widget_drag(PV, 1, &test_point, 0, V4(1, 1, 1, .5));

        gl_PV(PV);
        gl_begin(LINE_LOOP);
        gl_color(monokai.yellow);
        for (double theta = 0; theta < 2 * PI; theta += .1)
        {
            // gl_vertex(circle_center.x + circle_radius * cos(theta), circle_center.y + circle_radius * sin(theta));
            // gl_vertex(circle_center + circle_radius * V2(cos(theta), sin(theta)));
            gl_vertex(circle_center + circle_radius * e_theta(theta));
        }
        gl_end();
    }
}

void hw2c()
{
    init();
    Camera2D camera = {10};
    vec2 circle_center = V2(0, 0);
    double circle_radius = 2;
    vec2 *test_points = (vec2 *)malloc(4096 * sizeof(vec2));
    vec3 *test_points_colors = (vec3 *)malloc(4096 * sizeof(vec3));
    for (int i = 0; i < 4096; i++)
    {
        test_points[i] = V2(util_random_double(5, -5), util_random_double(5, -5));
    }

    while (begin_frame())
    {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        imgui_slider("r", &circle_radius, 0, 5);
        imgui_slider("x", &circle_center.x, -5, 5);
        imgui_slider("y", &circle_center.y, -5, 5);

        for (int i = 0; i < 4096; i++)
        {
            if (norm(test_points[i] - circle_center) < circle_radius)
            {
                test_points_colors[i] = monokai.green;
            }
            else
            {
                test_points_colors[i] = monokai.red;
            }
        }

        basic_draw(POINTS, PV, 4096, test_points, test_points_colors);
        widget_drag(PV, 4096, test_points, 0, V4(1, 1, 1, .5));

        gl_PV(PV);
        gl_begin(LINE_LOOP);
        gl_color(monokai.yellow);
        for (double theta = 0; theta < 2 * PI; theta += .1)
        {
            // gl_vertex(circle_center.x + circle_radius * cos(theta), circle_center.y + circle_radius * sin(theta));
            // gl_vertex(circle_center + circle_radius * V2(cos(theta), sin(theta)));
            gl_vertex(circle_center + circle_radius * e_theta(theta));
        }
        gl_end();
    }
}

void hw2d()
{
    init();
    int the_first_frame_that_k_ends_in_47 = 9; // not 39

    // begin don't modify this code
    if (!the_first_frame_that_k_ends_in_47)
    {
        unsigned int k = 0;
        int frame = 0;
        while (begin_frame())
        {
            unsigned int m_w = 1 + k;
            unsigned int m_z = 2 + k;
            m_z = 36969 * (m_z & 65535) + (m_z >> 16);
            m_w = 18000 * (m_w & 65535) + (m_w >> 16);
            k += (m_z << 17) + m_w;
            ++frame;
        }
    }
    else
    {
        while (begin_frame())
        {
            char r = char(the_first_frame_that_k_ends_in_47 * 233.111111111);
            char g = char(the_first_frame_that_k_ends_in_47 * 2833.33333333);
            g = char(g + 1);
            char b = char(the_first_frame_that_k_ends_in_47 * 210.666666667);
            clear_draw_buffer(r / 255., g / 255., b / 255., 1);
        }
    }
    // end don't modify this code
}

void hw2e()
{
    init();
    Camera2D camera = {5};
    bool playing = false;
    double t = 0;
    double a = 0;
    double b = 2;

    while (begin_frame())
    {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        imgui_checkbox("playing", &playing, 'p');
        if (playing)
        {
            t += 1. / 60;
        }
        if (imgui_button("reset", 'r'))
        {
            t = 0;
        }
        imgui_slider("t", &t, -2, 2);
        imgui_slider("a", &a, -4, 4);
        imgui_slider("b", &b, -4, 4);

        vec2 green_dot = V2(LERP(t, a, b), 1);
        vec2 red_dot = V2(CLAMPED_LERP(t, a, b), 0);
        vec2 purple_dot = V2(COS_LERP(t, a, b), -1);

        gl_PV(PV);
        gl_begin(LINES);
        gl_color(monokai.blue);
        gl_vertex(a, 10);
        gl_vertex(a, -10);
        gl_color(monokai.orange);
        gl_vertex(b, 10);
        gl_vertex(b, -10);
        gl_end();

        basic_draw(POINTS, PV, 1, &green_dot, monokai.green);
        basic_draw(POINTS, PV, 1, &red_dot, monokai.red);
        basic_draw(POINTS, PV, 1, &purple_dot, monokai.purple);
    }
}

void hw2f()
{
    init();
    Camera2D camera = {3};
    vec2 a = V2(-1, -1);
    vec2 b = V2(1, -1);
    vec2 c = V2(0, .5 * sqrt(3));
    vec2 p = V2(0, 0);
    vec3 alpha_beta_gamma = V3(1, 1, 1); // (alpha, beta, gamma)
    while (begin_frame())
    {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        mat3 M = {a.x, b.x, c.x, a.y, b.y, c.y, 1, 1, 1};
        alpha_beta_gamma = inverse(M) * V3(p.x, p.y, 1);

        vec3 color = alpha_beta_gamma.x * V3(1, 0, 0) + alpha_beta_gamma.y * V3(0, 1, 0) + alpha_beta_gamma.z * V3(0, 0, 1);
        basic_draw(POINTS, PV, 1, &p, color);
        widget_drag(PV, 1, &p, 15, color);

        imgui_readout("alpha", &alpha_beta_gamma.x);
        imgui_readout("beta", &alpha_beta_gamma.y);
        imgui_readout("gamma", &alpha_beta_gamma.z);

        gl_PV(PV);
        gl_begin(LINE_LOOP);
        gl_color(V3(1, 0, 0));
        gl_vertex(a);
        gl_color(V3(0, 1, 0));
        gl_vertex(b);
        gl_color(V3(0, 0, 1));
        gl_vertex(c);
        gl_end();
        widget_drag(PV, 1, &a);
        widget_drag(PV, 1, &b);
        widget_drag(PV, 1, &c);
    }
}

double random_number()
{
    return 2 * (rand() / (double)RAND_MAX) - 2;
}

double random_velocity()
{
    return (rand() / (double)RAND_MAX) / 50 + 0.03;
}

void hw2g()
{
    init();
    Camera2D camera = {5};
    bool playing = false;
    double t = 0;

    vec2 *vertex_positions = (vec2 *)malloc(100 * sizeof(vec2));
    for (int p = 0; p <= 100; p++)
    {
        vertex_positions[p] = {random_number(), 2 + random_number()};
    }

    vec2 *vertex_velocities = (vec2 *)malloc(100 * sizeof(vec2));
    for (int v = 0; v <= 100; v++)
    {
        vertex_velocities[v] = {random_velocity(), random_velocity()};
    }

    struct
    {
        vec3 brown = {58.8 / 255, 29.4 / 255, 0 / 255};
        vec3 light_gray = {225.0 / 255, 225.0 / 255, 225.0 / 255};
        vec3 lighter_brown = {90.0 / 255, 50.0 / 255, 0 / 255};
        vec3 dark_green0 = {0.0 / 255, 60.0 / 255, 0.0 / 255};
        vec3 dark_green1 = {0.0 / 255, 80.0 / 255, 0.0 / 255};
        vec3 dark_green2 = {0.0 / 255, 100.0 / 255, 0.0 / 255};
        vec3 sky_blue = {4.0 / 255, 99.0 / 255, 202.0 / 255};
        vec3 deep_sky_blue = {12.0 / 255, 20.0 / 255, 69.0 / 255};
        vec3 lighter_red = {249. / 255, 80. / 255, 150. / 255};
    } new_colors;

    while (begin_frame())
    {
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);

        vec3 sky_color = {COS_LERP(t, new_colors.sky_blue, new_colors.deep_sky_blue)};
        vec3 snow_color = {COS_LERP(t, monokai.white, new_colors.light_gray)};
        vec3 planet_color = {COS_LERP(t, monokai.orange, new_colors.light_gray)};
        vec3 cloud_color = {COS_LERP(t, monokai.white, monokai.gray)};

        imgui_checkbox("playing", &playing, 'p');
        if (playing)
        {
            t += 1. / 100;
            for (int n = 0; n < 100; n++)
            {
                if (vertex_positions[n].data[1] <= -2)
                {
                    vertex_velocities[n].data[0] = 0;
                    vertex_velocities[n].data[1] = 0;
                }
                vertex_positions[n].data[1] -= vertex_velocities[n].data[1];
            }
        }

        vec2 square[4] = {{-2, -2}, {0, -2}, {0, 2}, {-2, 2}};
        basic_draw(QUADS, PV, 4, square, sky_color);

        vec2 planet[4] = {{-0.6, 1.3}, {-0.1, 1.3}, {-0.1, 1.8}, {-0.6, 1.8}};
        vec2 planet_new[4] = {{-1.6, 1.3}, {-1.1, 1.3}, {-1.1, 1.8}, {-1.6, 1.8}};

        vec2 *planet_position = (vec2 *)malloc(4 * sizeof(vec2));
        for (int p = 0; p < 4; p++)
        {
            planet_position[p] = {COS_LERP(t, planet[p].x, planet_new[p].x),
                                  COS_LERP(t, planet[p].y, planet_new[p].y)};
        }

        basic_draw(QUADS, PV, 4, planet_position, planet_color);

        vec2 cloud1[4] = {{-1.9, 1.2}, {-1.1, 1.2}, {-1.1, 1.5}, {-1.9, 1.5}};
        basic_draw(QUADS, PV, 4, cloud1, cloud_color);

        vec2 cloud2[4] = {{-2, 1.0}, {-0.7, 1.0}, {-0.7, 1.2}, {-2, 1.2}};
        basic_draw(QUADS, PV, 4, cloud2, cloud_color);

        vec2 cloud3[4] = {{-0.4, 0.8}, {0, 0.8}, {0, 1.0}, {-0.4, 1.0}};
        basic_draw(QUADS, PV, 4, cloud3, cloud_color);

        vec2 tree_trunk1[4] = {{-2, -2}, {-1.7, -2}, {-1.7, 0.5}, {-2, 0.5}};
        basic_draw(QUADS, PV, 4, tree_trunk1, new_colors.lighter_brown);

        vec2 tree_trunk2[4] = {{-2, -2}, {-1.9, -2}, {-1.9, 0.5}, {-2, 0.5}};
        basic_draw(QUADS, PV, 4, tree_trunk2, new_colors.brown);

        vec2 tree_body1[4] = {{-2, -1.1}, {-1.1, -1.1}, {-1.1, -0.1}, {-2, -0.1}};
        basic_draw(QUADS, PV, 4, tree_body1, new_colors.dark_green0);

        vec2 tree_body2[4] = {{-1.8, -0.8}, {-1.1, -0.8}, {-1.1, -0.1}, {-1.8, -0.1}};
        basic_draw(QUADS, PV, 4, tree_body2, new_colors.dark_green1);

        vec2 tree_top1[4] = {{-2, 0.1}, {-1.5, 0.1}, {-1.5, 0.7}, {-2, 0.7}};
        basic_draw(QUADS, PV, 4, tree_top1, new_colors.dark_green1);

        vec2 tree_top2[4] = {{-1.9, 0.2}, {-1.5, 0.2}, {-1.5, 0.7}, {-1.9, 0.7}};
        basic_draw(QUADS, PV, 4, tree_top2, new_colors.dark_green2);

        vec2 snowman_body[4] = {{-0.8, -2}, {-0.1, -2}, {-0.1, -1.3}, {-0.8, -1.3}};
        basic_draw(QUADS, PV, 4, snowman_body, snow_color);

        vec2 snowman_head[4] = {{-0.6, -1.3}, {-0.3, -1.3}, {-0.3, -1.0}, {-0.6, -1.0}};
        basic_draw(QUADS, PV, 4, snowman_head, snow_color);

        vec2 snowman_hat1[4] = {{-0.15, -1.0}, {-0.75, -1.0}, {-0.75, -1.1}, {-0.15, -1.1}};
        basic_draw(QUADS, PV, 4, snowman_hat1, monokai.red);

        vec2 snowman_hat2[4] = {{-0.6, -1.1}, {-0.3, -1.1}, {-0.3, -0.6}, {-0.6, -0.6}};
        basic_draw(QUADS, PV, 4, snowman_hat2, monokai.red);

        vec2 snowman_hat3[4] = {{-0.5, -1.0}, {-0.3, -1.0}, {-0.3, -0.6}, {-0.5, -0.6}};
        basic_draw(QUADS, PV, 4, snowman_hat3, new_colors.lighter_red);

        vec2 snowman_nose[4] = {{-0.4, -1.25}, {-0.5, -1.25}, {-0.5, -1.2}, {-0.4, -1.2}};
        basic_draw(QUADS, PV, 4, snowman_nose, monokai.orange);

        vec2 snowman_dot1[1] = {{-0.45, -1.8}};
        basic_draw(POINTS, PV, 1, snowman_dot1, monokai.black);

        vec2 snowman_dot2[1] = {{-0.45, -1.6}};
        basic_draw(POINTS, PV, 1, snowman_dot2, monokai.black);

        vec2 snowman_hand1[2] = {{-0.8, -1.6}, {-1.2, -1.4}};
        basic_draw(LINE_STRIP, PV, 2, snowman_hand1, new_colors.brown);

        basic_draw(POINTS, PV, 100, vertex_positions, snow_color);
    }
}

void hw()
{
    hw2a();
    hw2b();
    hw2c();
    hw2d();
    hw2e();
    hw2f();
    hw2g();
}

// end submission

int main()
{
    // please change this call to init(...)
    // so you like where the popup window spawns :)
    //
    // void init(
    //     bool transparent_framebuffer,
    //     char *window_title,
    //     int screen_height_in_pixels,
    //     int window_top_left_init_x_in_pixels,
    //     int window_top_left_init_y_in_pixels
    // );
    init(true, "hw2", 540, 0, 100);
    hw();
    return 0;
}
