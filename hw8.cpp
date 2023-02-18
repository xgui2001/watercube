#define _CRT_SECURE_NO_WARNINGS
#include "snail.cpp"
#include "cow.cpp"
#include "_cow_supplement.cpp"
#include <vector>
#include <iostream>

// BEGIN fine to ignore
void hw8a_draw_textured_square(mat4 P, mat4 V, mat4 M, char *texture_filename)
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

struct
{
#define Q(deg) (2 * V3(cos(RAD(deg)), sin(RAD(deg)), -2))
    vec3 trivial_vertex_positions[3] = {Q(0), Q(120), Q(240)};
#undef Q
    vec3 trivial_vertex_colors[3] = {V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1)};
    BasicTriangleMesh3D trivial = {NELEMS(trivial_vertex_positions), trivial_vertex_positions, trivial_vertex_colors};

#define Q(deg, is_point) (2 * V3(cos(RAD(deg)), sin(RAD(deg)), (is_point) ? .5 : -.5))
    vec3 cycle_vertex_positions[12] = {Q(0, false), Q(20, false), Q(200, true), Q(120, false), Q(140, false), Q(320, true), Q(240, false), Q(260, false), Q(100, true), V3(5 * cos(RAD(240)), -1.5, 5 * sin(RAD(240))), V3(5 * cos(RAD(120)), -1.5, 5 * sin(RAD(120))), V3(5 * cos(RAD(0)), -1.5, 5 * sin(RAD(0)))};
#undef Q
    vec3 cycle_vertex_colors[12] = {monokai.brown, monokai.brown, monokai.brown, monokai.yellow, monokai.yellow, monokai.yellow, monokai.purple, monokai.purple, monokai.purple, monokai.red, monokai.red, monokai.red};
    BasicTriangleMesh3D cycle = {NELEMS(cycle_vertex_positions), cycle_vertex_positions, cycle_vertex_colors};

    vec3 tilt_vertex_positions[6] = {
        V3(-1, -.2, 1),
        V3(1, -.2, 1),
        V3(1, .2, -1),
        V3(-1, -.2, 1),
        V3(1, .2, -1),
        V3(-1, .2, -1),
    };
    vec3 tilt_vertex_colors[6] = {V3(1, 1, 1), V3(1, 1, 1), V3(0, 0, 1), V3(1, 1, 1), V3(0, 0, 1), V3(0, 0, 1)};
    BasicTriangleMesh3D tilt = {NELEMS(tilt_vertex_positions), tilt_vertex_positions, tilt_vertex_colors};

    vec3 clip2_vertex_positions[6] = {V3(-1, -1, 5), V3(1, -1, 5), V3(0, -1, 0)};
    vec3 clip2_vertex_colors[6] = {V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1)};
    BasicTriangleMesh3D clip2 = {NELEMS(clip2_vertex_positions), clip2_vertex_positions, clip2_vertex_colors};

    vec3 clip1_vertex_positions[6] = {V3(-1, -1, 0), V3(1, -1, 0), V3(0, -1, 5)};
    vec3 clip1_vertex_colors[6] = {V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1)};
    BasicTriangleMesh3D clip1 = {NELEMS(clip1_vertex_positions), clip1_vertex_positions, clip1_vertex_colors};
} hw8a_meshes;

struct
{
    bool draw_rays;
    bool draw_scene_3D = true;
    bool fully_transparent_film_plane_bg;
    bool draw_film_plane = true;
    bool draw_cube_at_observer;
    bool bunny_stress_test;
    double renderer_distance_to_film_plane = 1; // for visualization only (doesn't impact rendering)
} hw8a_tweaks;
// END fine to ignore

// mesh    -- the current scene
// light_p -- the position of the point light in world coordinates
BasicTriangleMesh3D *mesh;
vec3 light_p = V3(0, 2.5, 0);

// S                  -- the film plane side length in pixels
// color_buffer.data  -- an unsigned char array in r_0 g_0 b_0 a_0 r_1 ... format
//                       where, e.g., g_k is the green component of the k-th pixel
// hw8a_set_pixel(...) -- write color to pixel (i, j); !! does clamping for you
int S = 64;
Texture color_buffer;
void hw8a_set_pixel(int i, int j, vec3 color)
{
    color_buffer.data[4 * (j * S + i) + 0] = (unsigned char)(255 * CLAMP(color.r, 0, 1));
    color_buffer.data[4 * (j * S + i) + 1] = (unsigned char)(255 * CLAMP(color.g, 0, 1));
    color_buffer.data[4 * (j * S + i) + 2] = (unsigned char)(255 * CLAMP(color.b, 0, 1));
    color_buffer.data[4 * (j * S + i) + 3] = 255;
}

struct CastRayResult
{
    bool hit_at_least_one_triangle;
    vec3 base_color;
    double min_t;
    vec3 n_hit;
};
CastRayResult cast_ray(vec3 dir, vec3 o_renderer)
{
    bool hit_at_least_one_triangle = false;
    double min_t = INFINITY;
    vec3 base_color = {};
    vec3 n_hit;

    int num_triangles = mesh->num_vertices / 3;
    for (int triangle_i = 0; triangle_i < num_triangles; ++triangle_i)
    {
        // (a, b, c) -- triangle vertex positions in world coordinates
        // color_*   -- vertex colors
        // n         -- unit vector normal to triangle
        vec3 a, b, c;
        vec3 color_a, color_b, color_c;
        vec3 n;
        {
            a = mesh->vertex_positions[3 * triangle_i + 0];
            b = mesh->vertex_positions[3 * triangle_i + 1];
            c = mesh->vertex_positions[3 * triangle_i + 2];
            vec3 e1 = b - a;
            vec3 e2 = c - a;
            n = normalized(cross(e1, e2));
            if (mesh->vertex_colors != NULL)
            {
                color_a = mesh->vertex_colors[3 * triangle_i + 0];
                color_b = mesh->vertex_colors[3 * triangle_i + 1];
                color_c = mesh->vertex_colors[3 * triangle_i + 2];
            }
            else
            {
                vec3 fallback_color = V3(.5, .5, .5) + .5 * n;
                color_a = fallback_color;
                color_b = fallback_color;
                color_c = fallback_color;
            }

            vec4 answer_cord = inverse(M4(a.x, b.x, c.x, -dir.x, a.y, b.y, c.y, -dir.y, a.z, b.z, c.z, -dir.z, 1, 1, 1, 0)) * V4(o_renderer.x, o_renderer.y, o_renderer.z, 1);

            // vec3 p = o_renderer + answer_cord.w * dir; //point we hit on the triangle
            bool hit = false;

            if (answer_cord.x > 0 && answer_cord.y > 0 && answer_cord.z > 0 && answer_cord.w > TINY)
            {
                hit = true;
            }

            if (hit)
            {
                hit_at_least_one_triangle = true;
                if (answer_cord.w < min_t)
                {
                    min_t = answer_cord.w;

                    // the absolute of -t tells you if triangle is in the front
                    // smallest is closest to the camera
                    // if hit value less than min_t draw; only want to draw if absolute value less than min t
                    base_color = answer_cord.x * color_a + answer_cord.y * color_b + answer_cord.z * color_c;
                    n_hit = n;
                }
            }
        }
    }

    return {hit_at_least_one_triangle, base_color, min_t, n_hit};
};
void hw8a()
{
    init();

    // renderer -- the camera doing the actual rendering (raycasting)
    // observer -- the camera that's observing _everything_ (scene and ray-casting camera)
    Camera3D renderer = {4, RAD(60)};
    Camera3D observer = {6.5, RAD(60), RAD(30), RAD(-15), -2};
    { // (fine to ignore)
        color_buffer = {"color_buffer", S, S, 4, (unsigned char *)malloc(S * S * 4)};
        fancy_texture_create(&color_buffer);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }

    while (begin_frame())
    {
        { // tweaks (fine to ignore)
            imgui_checkbox("draw_rays", &hw8a_tweaks.draw_rays, 'z');
            imgui_checkbox("draw_film_plane", &hw8a_tweaks.draw_film_plane, 'x');
            imgui_checkbox("draw_scene_3D", &hw8a_tweaks.draw_scene_3D, 's');
            imgui_checkbox("draw_cube_at_observer", &hw8a_tweaks.draw_cube_at_observer, 'c');
            imgui_checkbox("fully_transparent_film_plane_bg", &hw8a_tweaks.fully_transparent_film_plane_bg, 't');
            imgui_checkbox("bunny_stress_test", &hw8a_tweaks.bunny_stress_test);

            imgui_slider("hw8a_tweaks.renderer_distance_to_film_plane", &hw8a_tweaks.renderer_distance_to_film_plane, 0, 5);
            imgui_slider("renderer.distance_to_origin", &renderer.distance_to_origin, 1, 10);
            imgui_slider("renderer.theta", &renderer.theta, RAD(-30), RAD(30), true);
            imgui_slider("renderer.phi", &renderer.phi, RAD(-30), RAD(30), true);
            imgui_slider("renderer.angle_of_view", &renderer.angle_of_view, RAD(1), RAD(178), true);
        }

        // *_renderer                   -- matrices for the renderer
        // *_observere                  -- matrices for the observer
        // film_plane_side_length_world -- the length of the film plane in world coordinates
        mat4 C_renderer, P_renderer, V_renderer;
        mat4 P_observer, V_observer, PV_observer;
        { // (fine to ignore)
         {// mesh
          static BasicTriangleMesh3D * examples[] = {&hw8a_meshes.trivial, &hw8a_meshes.cycle, &hw8a_meshes.tilt, &hw8a_meshes.clip2, &hw8a_meshes.clip1};
        static char *example_names[] = {"trivial", "cycle", "tilt", "clip2", "clip1"};
        static BasicTriangleMesh3D bunny = load_basic_mesh("data_basic_bunny", true);
        static int mesh_index = 0;
        imgui_slider("mesh_index", &mesh_index, 0, NELEMS(examples) - 1, 'j', 'k', true);
        mesh = examples[mesh_index];
        if (!hw8a_tweaks.bunny_stress_test)
        {
            _imgui_printf("example: %s", example_names[mesh_index]);
        }
        else
        {
            mesh = &bunny;
            _imgui_printf("i'm starting to worry our implementation wasn't very efficient");
        }
    }

    {     // C_*, P_*, PV_*
        { // reset
            static Camera3D _renderer_0 = renderer;
            if (imgui_button("reset C_renderer", 'r'))
            {
                renderer = _renderer_0;
            }
        }
        camera_move(&observer);
        C_renderer = camera_get_C(&renderer);
        P_renderer = tform_get_P_perspective(renderer.angle_of_view); // aspect <- 1
        V_renderer = inverse(C_renderer);
        { // C_observer <- C_renderer
            static bool clicked;
            bool clicked_this_frame;
            bool selected;
            bool released_this_frame = false;
            {
                char *name = "hold to toggle C_observer <- C_renderer";
                clicked_this_frame = imgui_button(name, KEY_TAB);
                if (clicked_this_frame)
                {
                    clicked = true;
                }
                selected = (imgui.selected_widget_ID == (void *)name);
                if (clicked && !selected)
                {
                    clicked = false;
                    released_this_frame = true;
                }
            }

            { // memcpy's
                static Camera3D safe;
                static bool draw_cube_push;
                if (clicked_this_frame)
                {
                    memcpy(&safe, &observer, sizeof(Camera3D));
                    draw_cube_push = hw8a_tweaks.draw_cube_at_observer;
                    hw8a_tweaks.draw_cube_at_observer = false;
                }
                if (selected)
                {
                    memcpy(&observer, &renderer, sizeof(Camera3D));
                }
                if (released_this_frame)
                {
                    memcpy(&observer, &safe, sizeof(Camera3D));
                    hw8a_tweaks.draw_cube_at_observer = draw_cube_push;
                }
            }
        }
        P_observer = camera_get_P(&observer);
        V_observer = camera_get_V(&observer);
        PV_observer = P_observer * V_observer;
    }
}

{  // render (your work here!)
 { // prep (fine to ignore)
  {// clear color_buffer
   if (hw8a_tweaks.fully_transparent_film_plane_bg){
       memset(color_buffer.data, 0, S *S * 4);
}
else
{
    memset(color_buffer.data, 255, S * S * 4);
    int pixelIndex = 0;
    for (int j = 0; j < S; ++j)
    {
        for (int i = 0; i < S; ++i)
        {
            color_buffer.data[4 * pixelIndex++ + 3] = (((j + i) % 2) == 0) ? 200 : 220;
        }
    }
}
}

{ // light_p
    jank_widget_translate3D(PV_observer, 1, &light_p);
    basic_draw(POINTS, PV_observer, 1, &light_p, monokai.white);
}
}

{ // render (your work here!)
    // renderer.angle_of_view -- _full_ (not half) angle of view of the renderer
    // *_renderer             -- the axes and origin of the renderer
    // NOTE: o_renderer is where rays originate from
    // NOTE: -z_renderer points from o_renderer to the center of the film plane
    vec3 x_renderer, y_renderer, z_renderer, o_renderer;
    {
        camera_get_coordinate_system(&renderer, NULL, x_renderer.data, y_renderer.data, z_renderer.data, o_renderer.data);
    }

    { // write to color_buffer.data (your work here!)

        // gl_* will be useful for debugging direction (dir)
        gl_begin(LINES);
        gl_PV(PV_observer);
        {

            for (int i = 0; i < S; ++i)
            {
                for (int j = 0; j < S; ++j)
                {

                    double theta = renderer.angle_of_view;
                    vec3 tmp = V3(i - double(S) / 2, j - double(S) / 2, -(double(S) / 2) / tan(theta / 2)); // student answer from board
                    vec3 dir = tmp.x * x_renderer + tmp.y * y_renderer + tmp.z * z_renderer;
                    // if (hw8a_tweaks.draw_rays)
                    //{
                    // gl_color(monokai.red);
                    // gl_vertex(o_renderer);
                    // gl_vertex(o_renderer + dir);
                    //}
                    CastRayResult result = cast_ray(dir, o_renderer);
                    if (result.hit_at_least_one_triangle)
                    {
                        vec3 p_hit = o_renderer + result.min_t * dir;
                        vec3 light_dir = light_p - p_hit;
                        vec3 normal = result.n_hit;

                        vec3 norm = normalized(normal);
                        vec3 light_dir_norm = normalized(light_dir);
                        float diff = fmax(dot(norm, light_dir_norm), 0.0);
                        vec3 diffuse = diff * result.base_color;

                        float specular_strength = 0.1;
                        vec3 view_dir = normalized(o_renderer - p_hit);
                        vec3 reflect_dir = light_dir_norm - 2 * dot(light_dir_norm, normal) * normal;
                        float spec = pow(fmax(dot(view_dir, reflect_dir), 0.0), 32);
                        vec3 specular = specular_strength * spec * result.base_color;

                        bool triangle_is_facing_the_light = dot(normal, light_dir) > 0;
                        CastRayResult shadow = cast_ray(light_dir, p_hit);
                        bool lit = !shadow.hit_at_least_one_triangle;
                        // if in the light:specular and diffuse;
                        if (triangle_is_facing_the_light && lit)
                        {

                            hw8a_set_pixel(i, j, diffuse + specular + result.base_color);
                        }
                        else
                        {
                            hw8a_set_pixel(i, j, 0.5 * result.base_color);
                        }
                        if (hw8a_tweaks.draw_rays)
                        {
                            gl_color(result.base_color);
                            gl_vertex(o_renderer);
                            gl_vertex(p_hit);
                        }
                    }
                }

                // TODO intersect with the scene
                // for (...)
                //     if (...) {
                //         hw8a_set_pixel(i, j, color);
                //     }
                // }
            }
        }
    }
    gl_end();

    { // send updated texture to the GPU (fine to ignore)
        fancy_texture_update(&color_buffer);
    }
}
}

{ // observe (fine to ignore)
    if (hw8a_tweaks.draw_scene_3D)
    {
        basic_draw(TRIANGLES, PV_observer, *mesh, V3(1, 0, 1));
    }
    { // bespoke widget
        basic_draw(PV_observer * C_renderer, meshlib.basic_axes);

        if (hw8a_tweaks.draw_film_plane)
        {
            double film_plane_side_length_world = 2 * hw8a_tweaks.renderer_distance_to_film_plane * tan(renderer.angle_of_view / 2);
            mat4 M = C_renderer * Translation(0, 0, -hw8a_tweaks.renderer_distance_to_film_plane) * Scaling(film_plane_side_length_world / 2);
            hw8a_draw_textured_square(P_observer, V_observer, M, color_buffer.filename);
            { // outline
                vec3 tmp[] = {
                    {-1, -1, 0},
                    {-1, 1, 0},
                    {1, 1, 0},
                    {1, -1, 0},
                };
                basic_draw(LINE_LOOP, P_observer * V_observer * M, 4, tmp);
            }
        }
        if (hw8a_tweaks.draw_cube_at_observer)
        {
            basic_draw(PV_observer * C_renderer * Scaling(.2), meshlib.basic_box, .5 * monokai.gray);
        }
    }
}
}
}

char *hw8b_frag = R""""(
    #version 330 core

    vec3 plasma(float t) {
        const vec3 c0 = vec3(0.05873234392399702, 0.02333670892565664, 0.5433401826748754);
        const vec3 c1 = vec3(2.176514634195958, 0.2383834171260182, 0.7539604599784036);
        const vec3 c2 = vec3(-2.689460476458034, -7.455851135738909, 3.110799939717086);
        const vec3 c3 = vec3(6.130348345893603, 42.3461881477227, -28.51885465332158);
        const vec3 c4 = vec3(-11.10743619062271, -82.66631109428045, 60.13984767418263);
        const vec3 c5 = vec3(10.02306557647065, 71.41361770095349, -54.07218655560067);
        const vec3 c6 = vec3(-3.658713842777788, -22.93153465461149, 18.19190778539828);
        return c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6)))));
    }

    vec3 inferno(float t) {

    const vec3 c0 = vec3(0.0002189403691192265, 0.001651004631001012, -0.01948089843709184);
    const vec3 c1 = vec3(0.1065134194856116, 0.5639564367884091, 3.932712388889277);
    const vec3 c2 = vec3(11.60249308247187, -3.972853965665698, -15.9423941062914);
    const vec3 c3 = vec3(-41.70399613139459, 17.43639888205313, 44.35414519872813);
    const vec3 c4 = vec3(77.162935699427, -33.40235894210092, -81.80730925738993);
    const vec3 c5 = vec3(-71.31942824499214, 32.62606426397723, 73.20951985803202);
    const vec3 c6 = vec3(25.13112622477341, -12.24266895238567, -23.07032500287172);

    return c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6)))));

}


    // https://iquilezles.org/articles/distfunctions/
    float dot2(vec2 v) { return dot(v,v); }
    float dot2(vec3 v) { return dot(v,v); }
    float ndot(vec2 a, vec2 b) { return a.x*b.x - a.y*b.y; }

    // for computing ray directions
    uniform vec3 x_renderer;
    uniform vec3 y_renderer;
    uniform vec3 z_renderer;
    uniform vec3 o_renderer;
    uniform float renderer_angle_of_view;
    uniform vec2 iResolution;

    uniform float time; // for time-varying distance fields

    out vec4 fragColor;

    // begin https://iquilezles.org/articles/distfunctions/
    float sdSphere(vec3 p, float r) {
        return length(p) - r;
    }
    float sdTorus(vec3 p, vec2 t) {
        vec2 q = vec2(length(p.xz)-t.x,p.y);
        return length(q)-t.y;
    }
    float sdCone(vec3 p, vec2 c, float h ){
  // c is the sin/cos of the angle, h is height
  // Alternatively pass q instead of (c,h),
  // which is the point at the base in 2D
  vec2 q = h*vec2(c.x/c.y,-1.0);
    
  vec2 w = vec2( length(p.xz), p.y );
  vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
  vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
  float k = sign( q.y );
  float d = min(dot( a, a ),dot(b, b));
  float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
  return sqrt(d)*sign(s);
}
    float sdSphere2(vec3 p, float r) {
        return length(p) - r;
    }

    float sdCutSphere( vec3 p, float r, float h )
{
  // sampling independent computations (only depend on shape)
  float w = sqrt(r*r-h*h);

  // sampling dependant computations
  vec2 q = vec2( length(p.xz), p.y );
  float s = max( (h-r)*q.x*q.x+w*w*(h+r-2.0*q.y), h*q.x-w*q.y );
  return (s<0.0) ? length(q)-r :
         (q.x<w) ? h - q.y     :
                   length(q-vec2(w,h));
}

float sdRoundedCylinder( vec3 p, float ra, float rb, float h )
{
  vec2 d = vec2( length(p.xz)-2.0*ra+rb, abs(p.y) - h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0)) - rb;
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
  vec3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

float sdCappedCone( vec3 p, float h, float r1, float r2 )
{
  vec2 q = vec2( length(p.xz), p.y );
  vec2 k1 = vec2(r2,h);
  vec2 k2 = vec2(r2-r1,2.0*h);
  vec2 ca = vec2(q.x-min(q.x,(q.y<0.0)?r1:r2), abs(q.y)-h);
  vec2 cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot2(k2), 0.0, 1.0 );
  float s = (cb.x<0.0 && ca.y<0.0) ? -1.0 : 1.0;
  return s*sqrt( min(dot2(ca),dot2(cb)) );
}

float sdVerticalCapsule( vec3 p, float h, float r )
{
  p.y -= clamp( p.y, 0.0, h );
  return length( p ) - r;
}

    // end

    vec3 march(vec3 o, vec3 dir) {
        // https://michaelwalczyk.com/blog-ray-marching.html

        const int MAX_STEPS = 64;
        const float HIT_TOLERANCE = 0.001;
        const float MAX_MARCH_DISTANCE = 100.0;

        // p   -- current position along ray
        // o   -- camera origin             
        // t   -- distance marched along ray
        // dir -- camera direction          
        // f   -- distance to surface       

        float t = 0.0;
        int step = 0;

        while (step++ < MAX_STEPS && t < MAX_MARCH_DISTANCE) {
            // get current position of ray's head
            vec3 p = o + t * dir;
                
                vec3 a = vec3(8, 8, 8);
                vec3 q_new = mod(p+0.5*a,a)-0.5*a;

                const float k2 = 10.0; // or some other amount
                float c2 = cos(k2*p.y);
                float s2 = sin(k2*p.y);
                mat2  m2 = mat2(c2,-s2,s2,c2);
                vec3  q2 = vec3(m2*p.xz,p.y);

            // compute distance to implicit surface
            float f = MAX_MARCH_DISTANCE; {
                {
                    //vec3 sphere_position = vec3(0.0, 0.0, 0.0);
                   // float sphere_radius = 1.0;
                    //float distance_to_sphere = sdSphere(-(q_new - sphere_position), sphere_radius);
                   // f = min(f, distance_to_sphere);
                }
                {
                    // f = min(f, sdTorus(p - vec3(0.0, 1.0 * sin(time), 0.0), vec2(1.0, 0.25)));
                    vec3 torus_position = vec3(0.0, sin(time), 0.0);
                    float torus_major_radius = 1.5;
                    float torus_minor_radius = 0.1;
                    vec2 torus_radii = vec2(torus_major_radius, torus_minor_radius);
                    float distance_to_torus = sdTorus(-(q_new - torus_position), torus_radii);
                    f = min(f, distance_to_torus);
                }
                {
                    vec3 half_position = vec3(0.0, 0.6, 0.0);
                    float half_radius = 1.03;
                    float half_height = 0.3;

                    float distance_to_half = sdCutSphere((q_new - half_position), half_radius, half_height);
                    f = min(f, distance_to_half);
                }
                {
                    vec3 cylinder_position = vec3(0.0, 0.0, 0.0);
                    float ra_cylinder = 0.5;
                    float rb_cylinder = 0.1;
                    float height_cylinder = 0.8;

                    float distance_to_cylinder = sdRoundedCylinder((q_new - cylinder_position), ra_cylinder, rb_cylinder, height_cylinder);
                    f = min(f, distance_to_cylinder);

                }
                {
                    vec3 cylinder_position1 = vec3(0.5, -1.0, 0.0);
                    float ra_cylinder1 = 0.21;
                    float rb_cylinder1 = 0.1;
                    float height_cylinder1 = 0.6;

                    float distance_to_cylinder1 = sdRoundedCylinder((p - cylinder_position1), ra_cylinder1, rb_cylinder1, height_cylinder1);
                    f = min(f, distance_to_cylinder1);
                }
                {
                    vec3 cylinder_position2 = vec3(-0.5, -1.0, 0.0);
                    float ra_cylinder2 = 0.21;
                    float rb_cylinder2 = 0.1;
                    float height_cylinder2 = 0.6;

                    float distance_to_cylinder2 = sdRoundedCylinder((q_new - cylinder_position2), ra_cylinder2, rb_cylinder2, height_cylinder2);
                    f = min(f, distance_to_cylinder2);
                }
                {
                    vec3 window_position = vec3(0.0, 0.45, 0.35);
                    vec3 window_a = vec3(0.4, 0.45, 0.35);
                    vec3 window_b = vec3(-0.4, 0.45, 0.35);
                    float window_radius = 0.35;

                    float distance_to_window = sdCapsule((q_new -window_position), window_a,  window_b, window_radius);
                    f = min(f, distance_to_window);
                }
                {
                    vec3 bag_position = vec3(0.0, 0.3, -0.8);
                    float bag_height = 0.4;
                    float bag_radius = 0.6;
                    
                    float distance_to_bag = sdVerticalCapsule(q_new -bag_position, bag_height, bag_radius);
                    f = min(f, distance_to_bag);
                }
                {
                    vec3 bay_position = vec3(0.0, -1.8, 0.0);
                    float bay_height = 0.15;
                    float bay_r1 = 1.5;
                    float bay_r2 = 1.2;

                    float distance_to_bay = sdCappedCone(q_new -bay_position, bay_height, bay_r1, bay_r2);
                    f = min(f, distance_to_bay);
                }
            }

            if (f < HIT_TOLERANCE) { // hit!
                return inferno(.5 + .5 * cos(t));
            }

            // NOTE if you're getting weird "overstepping" artifacts
            // (weird missing slices in the geometry)               
            // a (hacky) solution is to replace t += f; with e.g.   
            //t += min(f, .1);
            t += f;
        }
        return vec3(0.0);
    }

    void main() {
        vec3 o = o_renderer;
        vec3 dir; {
            // NOTE assume unit distance to film plane
            vec2 ds; { // [-R, R]
                float theta = renderer_angle_of_view / 2;
                float _R = tan(theta);
                ds = gl_FragCoord.xy;
                ds -= vec2(iResolution.x / 2, iResolution.y / 2);
                ds *= _R * 2. / iResolution.y;
            }
            // vec3 p_world = o_renderer - z_renderer + dx * x_renderer + dy * y_renderer;
            // dir = p_world - o_renderer;
            dir = -z_renderer + ds.x * x_renderer + ds.y * y_renderer;
        }
        vec3 col = march(o, dir);
        fragColor = vec4(col, 1);
    }
)"""";

/*
void hw8b()
{
    init();

    // mesh
    int num_vertices = 4;
    int num_triangles = 2;
    vec3 vertex_positions[] = {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
    int3 triangle_indices[] = {{0, 1, 2}, {0, 2, 3}};

    // shaders
    char *vert = R""""(
        #version 330 core
        layout (location = 0) in vec3 _p_model;
        void main() {
            gl_Position = vec4(_p_model, 1);
        }
    )"""";

    int shader_program = shader_build_program(vert, hw8b_frag);
    ASSERT(shader_program);

    // misc opengl
    GLuint VAO, VBO, EBO;
    {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);
    }

    // misc cow
    bool playing = false;
    double time = 0;
    Camera3D renderer = {8, RAD(45)};

    while (begin_frame())
    {
        camera_move(&renderer);

        {     // imgui
            { // reset
                static Camera3D _camera_0 = renderer;
                if (imgui_button("reset", 'r'))
                {
                    renderer = _camera_0;
                    time = 0;
                }
            }
            imgui_checkbox("playing", &playing, 'p');
        }

        { // draw
            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, num_vertices * 3 * sizeof(double), vertex_positions, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);
            glEnableVertexAttribArray(0);

            glUseProgram(shader_program);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * num_triangles * sizeof(int), triangle_indices, GL_DYNAMIC_DRAW);

            shader_set_uniform_vec2(shader_program, "iResolution", window_get_dimensions_in_pixels());
            {
                shader_set_uniform_double(shader_program, "time", time);

                shader_set_uniform_double(shader_program, "renderer_angle_of_view", renderer.angle_of_view);
                vec3 x_renderer, y_renderer, z_renderer, o_renderer;
                {
                    camera_get_coordinate_system(&renderer, NULL, x_renderer.data, y_renderer.data, z_renderer.data, o_renderer.data);
                }
                shader_set_uniform_vec3(shader_program, "x_renderer", x_renderer);
                shader_set_uniform_vec3(shader_program, "y_renderer", y_renderer);
                shader_set_uniform_vec3(shader_program, "z_renderer", z_renderer);
                shader_set_uniform_vec3(shader_program, "o_renderer", o_renderer);
            }

            glDrawElements(GL_TRIANGLES, 3 * num_triangles, GL_UNSIGNED_INT, NULL);
        }

        if (playing || input.key_pressed['.'])
        {
            time += .0167;
        }
    }
}
*/

double random_position()
{
    return rand() % 100;
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
void hwfinal()
{
    init();

    bool playing = false;

    double h = 1. / 60.;
    double g = -9.81; // gravitational acceleration
    vec2 gravity = h * g * V2(cos(-RAD(90)), sin(-RAD(90)));

    double diameter = 100;
    vec2 *box = (vec2 *)calloc(4, sizeof(vec2));

    vec2 *particle_positions = (vec2 *)calloc(100, sizeof(vec2));
    for (int p = 0; p < 100; p++)
    {
        particle_positions[p] = {random_position(), random_position()};
    }

    vec2 *vertex_velocities = (vec2 *)calloc(100, sizeof(vec2));

    while (begin_frame())
    {
        static Camera2D camera = {200, 100, 50};
        camera_move(&camera);
        mat4 PV = camera_get_PV(&camera);
        gl_PV(PV);

        vec2 s_mouse = input_get_mouse_position_in_world_coordinates(PV);

        vec2 particle_norm[500][500] = {};
        double particle_weight[500][500] = {};
        double particle_weight_sum[500] = {};
        double particle_pressure[500] = {};

        static double size = 10.0;
        static double time = 0;
        static int new_particle_count = 0;

        imgui_checkbox("playing", &playing, 'p');
        imgui_slider("diameter", &diameter, 100.0, 200.0);

        box[0] = {diameter, diameter};
        box[1] = {diameter, 0};
        box[2] = {0, 0};
        box[3] = {0, diameter};

        double h = 1. / 60.;

        basic_draw(LINE_LOOP, PV, 4, box, monokai.white);

        if (playing)
        {
            time += .01;

            for (int g = 0; g < 100; g++)
            {
                // gravity
                vertex_velocities[g].y -= h * gravity.y * 1 / 100;
            }

            // collision detection inbetween particles

            for (int n = 0; n < 100; n++)
            {
                particle_pressure[n] = 0.;

                for (int m = 0; m < 100; m++)
                {
                    particle_norm[n][m] = {0.0, 0.0};
                    particle_weight[n][m] = 0.;

                    if (m != n)
                    {
                        // if particle close to each other
                        if (norm(particle_positions[n] - particle_positions[m]) < 1 && !IS_ZERO(norm(particle_positions[n] - particle_positions[m])))
                        {

                            // keep track of the normalized position for all particles in regards to the specific particle
                            particle_norm[n][m] = normalized(particle_positions[n] - particle_positions[m]);
                            // calculate the weight of the contacts of each particle
                            particle_weight[n][m] += (1.0 - norm(particle_positions[n] - particle_positions[m]) / 1.0);
                            // calculate pressure of each particle (constants: coefficient of weight, avg weight <- assumed to be 0.5)
                            particle_pressure[n] += fmax(0, 0.2 * (particle_weight[n][m] - 0.3));
                        }
                    }
                    // apply repulsive force based on pressure
                    vertex_velocities[n] += (h * 0.5 * (particle_pressure[n] + particle_pressure[m]) * particle_weight[n][m] * particle_norm[n][m]) / 2;
                }
            }

            // if the particle is outside of box (??)
            // dragging particle inside box?
            /*
            if ((particle_positions[n].x >= 100))
            {
                particle_positions[n].x = 100 - TINY;
            }
            if ((particle_positions[n].y >= 100))
            {
                particle_positions[n].y = 100 - TINY;
            }
            if ((particle_positions[n].x <= 0))
            {
                particle_positions[n].x = 0 + TINY;
            }
            if ((particle_positions[n].y <= 0))
            {
                particle_positions[n].y = 0 + TINY;
            }
            */

            // collision detection with the box: if particle will be going outside of the box, set velocity = 0

            for (int n = 0; n < 100; n++)
            {
                if ((particle_positions[n].x >= diameter) && (vertex_velocities[n].x >= 0))
                {
                    particle_positions[n].x = diameter - TINY;
                    vertex_velocities[n].x *= 0;
                }

                if ((particle_positions[n].x <= 0) && (vertex_velocities[n].x <= 0))
                {
                    particle_positions[n].x = 0 + TINY;
                    vertex_velocities[n].x *= -0;
                }

                if ((particle_positions[n].y >= diameter) && (vertex_velocities[n].y >= 0))
                {
                    particle_positions[n].y = diameter - TINY;
                    vertex_velocities[n].y *= -0;
                }

                if ((particle_positions[n].y <= 0) && (vertex_velocities[n].y <= 0))
                {
                    particle_positions[n].y = 0 + TINY;
                    vertex_velocities[n].y *= -0;
                }

                particle_positions[n] += vertex_velocities[n];
            }
        }

        basic_draw(POINTS, PV, 100, particle_positions, color_rainbow_swirl(.5 * time), size);
    }
    free(particle_positions);
    // free(particle_weight);
    // free(particle_pressure);
    free(vertex_velocities);
}

int main()
{
    // hw8a();
    //  hw8b();
    // hwfinal();
    return 0;
}
