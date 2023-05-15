extern "C"
{
#include "raylib.h"
#include "raymath.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
}
#if defined(PLATFORM_WEB)
#include <emscripten/emscripten.h>
#endif

#include "common.h"
#include "vec.h"
#include "quat.h"
#include "spring.h"
#include "array.h"
#include "character.h"
#include "database.h"

#include <initializer_list>
#include <functional>

//--------------------------------------

static inline Vector3 to_Vector3(vec3 v)
{
    return (Vector3){ v.x, v.y, v.z };
}

//--------------------------------------

// Perform linear blend skinning and copy
// result into mesh data. Update and upload
// deformed vertex positions and normals to GPU
void deform_character_mesh(
  Mesh& mesh,
  const character& c,
  const slice1d<vec3> bone_anim_positions,
  const slice1d<quat> bone_anim_rotations,
  const slice1d<int> bone_parents)
{
    linear_blend_skinning_positions(
        slice1d<vec3>(mesh.vertexCount, (vec3*)mesh.vertices),
        c.positions,
        c.bone_weights,
        c.bone_indices,
        c.bone_rest_positions,
        c.bone_rest_rotations,
        bone_anim_positions,
        bone_anim_rotations);

    linear_blend_skinning_normals(
        slice1d<vec3>(mesh.vertexCount, (vec3*)mesh.normals),
        c.normals,
        c.bone_weights,
        c.bone_indices,
        c.bone_rest_rotations,
        bone_anim_rotations);

    UpdateMeshBuffer(mesh, 0, mesh.vertices, mesh.vertexCount * 3 * sizeof(float), 0);
    UpdateMeshBuffer(mesh, 2, mesh.normals, mesh.vertexCount * 3 * sizeof(float), 0);
}

Mesh make_character_mesh(character& c)
{
    Mesh mesh = { 0 };

    mesh.vertexCount = c.positions.size;
    mesh.triangleCount = c.triangles.size / 3;
    mesh.vertices = (float*)MemAlloc(c.positions.size * 3 * sizeof(float));
    mesh.texcoords = (float*)MemAlloc(c.texcoords.size * 2 * sizeof(float));
    mesh.normals = (float*)MemAlloc(c.normals.size * 3 * sizeof(float));
    mesh.indices = (unsigned short*)MemAlloc(c.triangles.size * sizeof(unsigned short));

    memcpy(mesh.vertices, c.positions.data, c.positions.size * 3 * sizeof(float));
    memcpy(mesh.texcoords, c.texcoords.data, c.texcoords.size * 2 * sizeof(float));
    memcpy(mesh.normals, c.normals.data, c.normals.size * 3 * sizeof(float));
    memcpy(mesh.indices, c.triangles.data, c.triangles.size * sizeof(unsigned short));

    UploadMesh(&mesh, true);

    return mesh;
}

//--------------------------------------

// Basic functionality to get gamepad input including deadzone and
// squaring of the stick location to increase sensitivity. To make
// all the other code that uses this easier, we assume stick is
// oriented on floor (i.e. y-axis is zero)

enum
{
    GAMEPAD_PLAYER = 0,
};

enum
{
    GAMEPAD_STICK_LEFT,
    GAMEPAD_STICK_RIGHT,
};

vec3 gamepad_get_stick(int stick, const float deadzone = 0.2f)
{

    float gamepadx = int(IsKeyDown(KEY_D)) - int(IsKeyDown(KEY_A));//GetGamepadAxisMovement(GAMEPAD_PLAYER, stick == GAMEPAD_STICK_LEFT ? GAMEPAD_AXIS_LEFT_X : GAMEPAD_AXIS_RIGHT_X);
    float gamepady = int(IsKeyDown(KEY_S)) - int(IsKeyDown(KEY_W));//GetGamepadAxisMovement(GAMEPAD_PLAYER, stick == GAMEPAD_STICK_LEFT ? GAMEPAD_AXIS_LEFT_Y : GAMEPAD_AXIS_RIGHT_Y);



    float gamepadmag = sqrtf(gamepadx*gamepadx + gamepady*gamepady);

    if (gamepadmag > deadzone)
    {
        float gamepaddirx = gamepadx / gamepadmag;
        float gamepaddiry = gamepady / gamepadmag;
        float gamepadclippedmag = gamepadmag > 1.0f ? 1.0f : gamepadmag*gamepadmag;
        gamepadx = gamepaddirx * gamepadclippedmag;
        gamepady = gamepaddiry * gamepadclippedmag;
    }
    else
    {
        gamepadx = 0.0f;
        gamepady = 0.0f;
    }

    return vec3(gamepadx, 0.0f, gamepady);
}

//--------------------------------------

float orbit_camera_update_azimuth(
    const float azimuth,
    const vec3 gamepadstick_right,
    const bool desired_strafe,
    const float dt)
{
    vec3 gamepadaxis = desired_strafe ? vec3() : gamepadstick_right;
    return azimuth + 2.0f * dt * -gamepadaxis.x;
}

float orbit_camera_update_altitude(
    const float altitude,
    const vec3 gamepadstick_right,
    const bool desired_strafe,
    const float dt)
{
    vec3 gamepadaxis = desired_strafe ? vec3() : gamepadstick_right;
    return clampf(altitude + 2.0f * dt * gamepadaxis.z, 0.0, 0.4f * PIf);
}

float orbit_camera_update_distance(
    const float distance,
    const float dt)
{
    float gamepadzoom =
        IsGamepadButtonDown(GAMEPAD_PLAYER, GAMEPAD_BUTTON_LEFT_TRIGGER_1)  ? +1.0f :
        IsGamepadButtonDown(GAMEPAD_PLAYER, GAMEPAD_BUTTON_RIGHT_TRIGGER_1) ? -1.0f : 0.0f;

    return clampf(distance +  10.0f * dt * gamepadzoom, 0.1f, 100.0f);
}

// Updates the camera using the orbit cam controls
void orbit_camera_update(
    Camera3D& cam,
    float& camera_azimuth,
    float& camera_altitude,
    float& camera_distance,
    const vec3 target,
    const vec3 gamepadstick_right,
    const bool desired_strafe,
    const float dt)
{

    camera_azimuth = orbit_camera_update_azimuth(camera_azimuth, gamepadstick_right, desired_strafe, dt);
    camera_altitude = orbit_camera_update_altitude(camera_altitude, gamepadstick_right, desired_strafe, dt);
    camera_distance = orbit_camera_update_distance(camera_distance, dt);

    quat rotation_azimuth = quat_from_angle_axis(camera_azimuth, vec3(0, 1, 0));
    vec3 position = quat_mul_vec3(rotation_azimuth, vec3(0, 0, camera_distance));
    vec3 axis = normalize(cross(position, vec3(0, 1, 0)));

    quat rotation_altitude = quat_from_angle_axis(camera_altitude, axis);

    vec3 eye = target + quat_mul_vec3(rotation_altitude, position);

    cam.target = (Vector3){ target.x, target.y, target.z };
    cam.position = (Vector3){ eye.x, eye.y, eye.z };

    UpdateCamera(&cam, CAMERA_ORBITAL);
}

//--------------------------------------

bool desired_strafe_update()
{
    return IsGamepadButtonDown(GAMEPAD_PLAYER, GAMEPAD_BUTTON_LEFT_TRIGGER_2) > 0.5f;
}

void desired_gait_update(
    float& desired_gait,
    float& desired_gait_velocity,
    const float dt,
    const float gait_change_halflife = 0.1f)
{

    simple_spring_damper_exact(
        desired_gait,
        desired_gait_velocity,
        IsGamepadButtonDown(GAMEPAD_PLAYER, GAMEPAD_BUTTON_RIGHT_FACE_DOWN) ? 1.0f : 0.0f,
        gait_change_halflife,
        dt);
}

vec3 desired_velocity_update(
    const vec3 gamepadstick_left,
    const float camera_azimuth,
    const quat simulation_rotation,
    const float fwrd_speed,
    const float side_speed,
    const float back_speed)
{

    // Find stick position in world space by rotating using camera azimuth
    vec3 global_stick_direction = quat_mul_vec3(
        quat_from_angle_axis(camera_azimuth, vec3(0, 1, 0)), gamepadstick_left);

    // Find stick position local to current facing direction
    vec3 local_stick_direction = quat_inv_mul_vec3(
        simulation_rotation, global_stick_direction);

    // Scale stick by forward, sideways and backwards speeds
    vec3 local_desired_velocity = local_stick_direction.z > 0.0 ?
        vec3(side_speed, 0.0f, fwrd_speed) * local_stick_direction :
        vec3(side_speed, 0.0f, back_speed) * local_stick_direction;

    // Re-orientate into the world space
    return quat_mul_vec3(simulation_rotation, local_desired_velocity);
}

quat desired_rotation_update(
    const quat desired_rotation,
    const vec3 gamepadstick_left,
    const vec3 gamepadstick_right,
    const float camera_azimuth,
    const bool desired_strafe,
    const vec3 desired_velocity)
{

    quat desired_rotation_curr = desired_rotation;

    // If strafe is active then desired direction is coming from right
    // stick as long as that stick is being used, otherwise we assume
    // forward facing
    if (desired_strafe)
    {
        vec3 desired_direction = quat_mul_vec3(quat_from_angle_axis(camera_azimuth, vec3(0, 1, 0)), vec3(0, 0, -1));

        if (length(gamepadstick_right) > 0.01f)
        {
            desired_direction = quat_mul_vec3(quat_from_angle_axis(camera_azimuth, vec3(0, 1, 0)), normalize(gamepadstick_right));
        }

        return quat_from_angle_axis(atan2f(desired_direction.x, desired_direction.z), vec3(0, 1, 0));
    }

    // If strafe is not active the desired direction comes from the left
    // stick as long as that stick is being used
    else if (length(gamepadstick_left) > 0.01f)
    {

        vec3 desired_direction = normalize(desired_velocity);
        return quat_from_angle_axis(atan2f(desired_direction.x, desired_direction.z), vec3(0, 1, 0));
    }

    // Otherwise desired direction remains the same
    else
    {
        return desired_rotation_curr;
    }
}

//--------------------------------------

// Moving the root is a little bit difficult when we have the
// inertializer set up in the way we do. Essentially we need
// to also make sure to adjust all of the locations where
// we are transforming the data to and from as well as the
// offsets being blended out
void inertialize_root_adjust(
    vec3& offset_position,
    vec3& transition_src_position,
    quat& transition_src_rotation,
    vec3& transition_dst_position,
    quat& transition_dst_rotation,
    vec3& position,
    quat& rotation,
    const vec3 input_position,
    const quat input_rotation)
{

    // Find the position difference and add it to the state and transition location
    vec3 position_difference = input_position - position;
    position = position_difference + position;
    transition_dst_position = position_difference + transition_dst_position;

    // Find the point at which we want to now transition from in the src data
    transition_src_position = transition_src_position + quat_mul_vec3(transition_src_rotation,
        quat_inv_mul_vec3(transition_dst_rotation, position - offset_position - transition_dst_position));
    transition_dst_position = position;
    offset_position = vec3();

    // Find the rotation difference. We need to normalize here or some error can accumulate
    // over time during adjustment.
    quat rotation_difference = quat_normalize(quat_mul_inv(input_rotation, rotation));

    // Apply the rotation difference to the current rotation and transition location
    rotation = quat_mul(rotation_difference, rotation);
    transition_dst_rotation = quat_mul(rotation_difference, transition_dst_rotation);
}

void inertialize_pose_reset(
    slice1d<vec3> bone_offset_positions,
    slice1d<vec3> bone_offset_velocities,
    slice1d<quat> bone_offset_rotations,
    slice1d<vec3> bone_offset_angular_velocities,
    vec3& transition_src_position,
    quat& transition_src_rotation,
    vec3& transition_dst_position,
    quat& transition_dst_rotation,
    const vec3 root_position,
    const quat root_rotation)
{

    bone_offset_positions.zero();
    bone_offset_velocities.zero();
    bone_offset_rotations.set(quat());
    bone_offset_angular_velocities.zero();

    transition_src_position = root_position;
    transition_src_rotation = root_rotation;
    transition_dst_position = vec3();
    transition_dst_rotation = quat();
}

// This function transitions the inertializer for
// the full character. It takes as input the current
// offsets, as well as the root transition locations,
// current root state, and the full pose information
// for the pose being transitioned from (src) as well
// as the pose being transitioned to (dst) in their
// own animation spaces.
void inertialize_pose_transition(
    slice1d<vec3> bone_offset_positions,
    slice1d<vec3> bone_offset_velocities,
    slice1d<quat> bone_offset_rotations,
    slice1d<vec3> bone_offset_angular_velocities,
    vec3& transition_src_position,
    quat& transition_src_rotation,
    vec3& transition_dst_position,
    quat& transition_dst_rotation,
    const vec3 root_position,
    const vec3 root_velocity,
    const quat root_rotation,
    const vec3 root_angular_velocity,
    const slice1d<vec3> bone_src_positions,
    const slice1d<vec3> bone_src_velocities,
    const slice1d<quat> bone_src_rotations,
    const slice1d<vec3> bone_src_angular_velocities,
    const slice1d<vec3> bone_dst_positions,
    const slice1d<vec3> bone_dst_velocities,
    const slice1d<quat> bone_dst_rotations,
    const slice1d<vec3> bone_dst_angular_velocities)
{

    // First we record the root position and rotation
    // in the animation data for the source and destination
    // animation
    transition_dst_position = root_position;
    transition_dst_rotation = root_rotation;
    transition_src_position = bone_dst_positions(0);
    transition_src_rotation = bone_dst_rotations(0);

    // We then find the velocities so we can transition the
    // root inertiaizers
    vec3 world_space_dst_velocity = quat_mul_vec3(transition_dst_rotation,
        quat_inv_mul_vec3(transition_src_rotation, bone_dst_velocities(0)));

    vec3 world_space_dst_angular_velocity = quat_mul_vec3(transition_dst_rotation,
        quat_inv_mul_vec3(transition_src_rotation, bone_dst_angular_velocities(0)));

    // Transition inertializers recording the offsets for
    // the root joint
    inertialize_transition(
        bone_offset_positions(0),
        bone_offset_velocities(0),
        root_position,
        root_velocity,
        root_position,
        world_space_dst_velocity);

    inertialize_transition(
        bone_offset_rotations(0),
        bone_offset_angular_velocities(0),
        root_rotation,
        root_angular_velocity,
        root_rotation,
        world_space_dst_angular_velocity);

    // Transition all the inertializers for each other bone
    for (int i = 1; i < bone_offset_positions.size; i++)
    {
        inertialize_transition(
            bone_offset_positions(i),
            bone_offset_velocities(i),
            bone_src_positions(i),
            bone_src_velocities(i),
            bone_dst_positions(i),
            bone_dst_velocities(i));

        inertialize_transition(
            bone_offset_rotations(i),
            bone_offset_angular_velocities(i),
            bone_src_rotations(i),
            bone_src_angular_velocities(i),
            bone_dst_rotations(i),
            bone_dst_angular_velocities(i));
    }
}

// This function updates the inertializer states. Here
// it outputs the smoothed animation (input plus offset)
// as well as updating the offsets themselves. It takes
// as input the current playing animation as well as the
// root transition locations, a halflife, and a dt
void inertialize_pose_update(
    slice1d<vec3> bone_positions,
    slice1d<vec3> bone_velocities,
    slice1d<quat> bone_rotations,
    slice1d<vec3> bone_angular_velocities,
    slice1d<vec3> bone_offset_positions,
    slice1d<vec3> bone_offset_velocities,
    slice1d<quat> bone_offset_rotations,
    slice1d<vec3> bone_offset_angular_velocities,
    const slice1d<vec3> bone_input_positions,
    const slice1d<vec3> bone_input_velocities,
    const slice1d<quat> bone_input_rotations,
    const slice1d<vec3> bone_input_angular_velocities,
    const vec3 transition_src_position,
    const quat transition_src_rotation,
    const vec3 transition_dst_position,
    const quat transition_dst_rotation,
    const float halflife,
    const float dt)
{

    // First we find the next root position, velocity, rotation
    // and rotational velocity in the world space by transforming
    // the input animation from it's animation space into the
    // space of the currently playing animation.
    vec3 world_space_position = quat_mul_vec3(transition_dst_rotation,
        quat_inv_mul_vec3(transition_src_rotation,
            bone_input_positions(0) - transition_src_position)) + transition_dst_position;

    vec3 world_space_velocity = quat_mul_vec3(transition_dst_rotation,
        quat_inv_mul_vec3(transition_src_rotation, bone_input_velocities(0)));

    // Normalize here because quat inv mul can sometimes produce
    // unstable returns when the two rotations are very close.
    quat world_space_rotation = quat_normalize(quat_mul(transition_dst_rotation,
        quat_inv_mul(transition_src_rotation, bone_input_rotations(0))));

    vec3 world_space_angular_velocity = quat_mul_vec3(transition_dst_rotation,
        quat_inv_mul_vec3(transition_src_rotation, bone_input_angular_velocities(0)));

    // Then we update these two inertializers with these new world space inputs
    inertialize_update(
        bone_positions(0),
        bone_velocities(0),
        bone_offset_positions(0),
        bone_offset_velocities(0),
        world_space_position,
        world_space_velocity,
        halflife,
        dt);

    inertialize_update(
        bone_rotations(0),
        bone_angular_velocities(0),
        bone_offset_rotations(0),
        bone_offset_angular_velocities(0),
        world_space_rotation,
        world_space_angular_velocity,
        halflife,
        dt);

    // Then we update the inertializers for the rest of the bones
    for (int i = 1; i < bone_positions.size; i++)
    {
        inertialize_update(
            bone_positions(i),
            bone_velocities(i),
            bone_offset_positions(i),
            bone_offset_velocities(i),
            bone_input_positions(i),
            bone_input_velocities(i),
            halflife,
            dt);

        inertialize_update(
            bone_rotations(i),
            bone_angular_velocities(i),
            bone_offset_rotations(i),
            bone_offset_angular_velocities(i),
            bone_input_rotations(i),
            bone_input_angular_velocities(i),
            halflife,
            dt);
    }
}

//--------------------------------------

// Copy a part of a feature vector from the
// matching database into the query feature vector
void query_copy_denormalized_feature(
    slice1d<float> query,
    int& offset,
    const int size,
    const slice1d<float> features,
    const slice1d<float> features_offset,
    const slice1d<float> features_scale)
{

    for (int i = 0; i < size; i++)
    {
        query(offset + i) = features(offset + i) * features_scale(offset + i) + features_offset(offset + i);
    }

    offset += size;
}

// Compute the query feature vector for the current
// trajectory controlled by the gamepad.
void query_compute_trajectory_position_feature(
    slice1d<float> query,
    int& offset,
    const vec3 root_position,
    const quat root_rotation,
    const slice1d<vec3> trajectory_positions)
{

    vec3 traj0 = quat_inv_mul_vec3(root_rotation, trajectory_positions(1) - root_position);
    vec3 traj1 = quat_inv_mul_vec3(root_rotation, trajectory_positions(2) - root_position);
    vec3 traj2 = quat_inv_mul_vec3(root_rotation, trajectory_positions(3) - root_position);

    query(offset + 0) = traj0.x;
    query(offset + 1) = traj0.z;
    query(offset + 2) = traj1.x;
    query(offset + 3) = traj1.z;
    query(offset + 4) = traj2.x;
    query(offset + 5) = traj2.z;

    offset += 6;
}

// Same but for the trajectory direction
void query_compute_trajectory_direction_feature(
    slice1d<float> query,
    int& offset,
    const quat root_rotation,
    const slice1d<quat> trajectory_rotations)
{

    vec3 traj0 = quat_inv_mul_vec3(root_rotation, quat_mul_vec3(trajectory_rotations(1), vec3(0, 0, 1)));
    vec3 traj1 = quat_inv_mul_vec3(root_rotation, quat_mul_vec3(trajectory_rotations(2), vec3(0, 0, 1)));
    vec3 traj2 = quat_inv_mul_vec3(root_rotation, quat_mul_vec3(trajectory_rotations(3), vec3(0, 0, 1)));

    query(offset + 0) = traj0.x;
    query(offset + 1) = traj0.z;
    query(offset + 2) = traj1.x;
    query(offset + 3) = traj1.z;
    query(offset + 4) = traj2.x;
    query(offset + 5) = traj2.z;

    offset += 6;
}

//--------------------------------------

// Collide against the obscales which are
// essentially bounding boxes of a given size
vec3 simulation_collide_obstacles(
    const vec3 prev_pos,
    const vec3 next_pos,
    const slice1d<vec3> obstacles_positions,
    const slice1d<vec3> obstacles_scales,
    const float radius = 0.6f)
{

    vec3 dx = next_pos - prev_pos;
    vec3 proj_pos = prev_pos;

    // Substep because I'm too lazy to implement CCD
    int substeps = 1 + (int)(length(dx) * 5.0f);

    for (int j = 0; j < substeps; j++)
    {
        proj_pos = proj_pos + dx / substeps;

        for (int i = 0; i < obstacles_positions.size; i++)
        {
            // Find nearest point inside obscale and push out
            vec3 nearest = clamp(proj_pos,
              obstacles_positions(i) - 0.5f * obstacles_scales(i),
              obstacles_positions(i) + 0.5f * obstacles_scales(i));

            if (length(nearest - proj_pos) < radius)
            {
                proj_pos = radius * normalize(proj_pos - nearest) + nearest;
            }
        }
    }

    return proj_pos;
}

// Taken from https://theorangeduck.com/page/spring-roll-call#controllers
void simulation_positions_update(
    vec3& position,
    vec3& velocity,
    vec3& acceleration,
    const vec3 desired_velocity,
    const float halflife,
    const float dt,
    const slice1d<vec3> obstacles_positions,
    const slice1d<vec3> obstacles_scales)
{

    float y = halflife_to_damping(halflife) / 2.0f;
    vec3 j0 = velocity - desired_velocity;
    vec3 j1 = acceleration + j0*y;
    float eydt = fast_negexpf(y*dt);

    vec3 position_prev = position;

    position = eydt*(((-j1)/(y*y)) + ((-j0 - j1*dt)/y)) +
        (j1/(y*y)) + j0/y + desired_velocity * dt + position_prev;
    velocity = eydt*(j0 + j1*dt) + desired_velocity;
    acceleration = eydt*(acceleration - j1*y*dt);

    position = simulation_collide_obstacles(
        position_prev,
        position,
        obstacles_positions,
        obstacles_scales);
}

void simulation_rotations_update(
    quat& rotation,
    vec3& angular_velocity,
    const quat desired_rotation,
    const float halflife,
    const float dt)
{
    simple_spring_damper_exact(
        rotation,
        angular_velocity,
        desired_rotation,
        halflife, dt);
}

// Predict what the desired velocity will be in the
// future. Here we need to use the future trajectory
// rotation as well as predicted future camera
// position to find an accurate desired velocity in
// the world space
void trajectory_desired_velocities_predict(
  slice1d<vec3> desired_velocities,
  const slice1d<quat> trajectory_rotations,
  const vec3 desired_velocity,
  const float camera_azimuth,
  const vec3 gamepadstick_left,
  const vec3 gamepadstick_right,
  const bool desired_strafe,
  const float fwrd_speed,
  const float side_speed,
  const float back_speed,
  const float dt)
{

    desired_velocities(0) = desired_velocity;

    for (int i = 1; i < desired_velocities.size; i++)
    {
        desired_velocities(i) = desired_velocity_update(
            gamepadstick_left,
            orbit_camera_update_azimuth(
                camera_azimuth, gamepadstick_right, desired_strafe, i * dt),
            trajectory_rotations(i),
            fwrd_speed,
            side_speed,
            back_speed);
    }
}

void trajectory_positions_predict(
    slice1d<vec3> positions,
    slice1d<vec3> velocities,
    slice1d<vec3> accelerations,
    const vec3 position,
    const vec3 velocity,
    const vec3 acceleration,
    const slice1d<vec3> desired_velocities,
    const float halflife,
    const float dt,
    const slice1d<vec3> obstacles_positions,
    const slice1d<vec3> obstacles_scales)
{

    positions(0) = position;
    velocities(0) = velocity;
    accelerations(0) = acceleration;

    for (int i = 1; i < positions.size; i++)
    {
        positions(i) = positions(i-1);
        velocities(i) = velocities(i-1);
        accelerations(i) = accelerations(i-1);

        simulation_positions_update(
            positions(i),
            velocities(i),
            accelerations(i),
            desired_velocities(i),
            halflife,
            dt,
            obstacles_positions,
            obstacles_scales);
    }
}

// Predict desired rotations given the estimated future
// camera rotation and other parameters
void trajectory_desired_rotations_predict(
  slice1d<quat> desired_rotations,
  const slice1d<vec3> desired_velocities,
  const quat desired_rotation,
  const float camera_azimuth,
  const vec3 gamepadstick_left,
  const vec3 gamepadstick_right,
  const bool desired_strafe,
  const float dt)
{

    desired_rotations(0) = desired_rotation;

    for (int i = 1; i < desired_rotations.size; i++)
    {
        desired_rotations(i) = desired_rotation_update(
            desired_rotations(i-1),
            gamepadstick_left,
            gamepadstick_right,
            orbit_camera_update_azimuth(
                camera_azimuth, gamepadstick_right, desired_strafe, i * dt),
            desired_strafe,
            desired_velocities(i));
    }
}

void trajectory_rotations_predict(
    slice1d<quat> rotations,
    slice1d<vec3> angular_velocities,
    const quat rotation,
    const vec3 angular_velocity,
    const slice1d<quat> desired_rotations,
    const float halflife,
    const float dt)
{

    rotations.set(rotation);
    angular_velocities.set(angular_velocity);

    for (int i = 1; i < rotations.size; i++)
    {
        simulation_rotations_update(
            rotations(i),
            angular_velocities(i),
            desired_rotations(i),
            halflife,
            i * dt);
    }
}



//--------------------------------------

void draw_axis(const vec3 pos, const quat rot, const float scale = 1.0f)
{

    vec3 axis0 = pos + quat_mul_vec3(rot, scale * vec3(1.0f, 0.0f, 0.0f));
    vec3 axis1 = pos + quat_mul_vec3(rot, scale * vec3(0.0f, 1.0f, 0.0f));
    vec3 axis2 = pos + quat_mul_vec3(rot, scale * vec3(0.0f, 0.0f, 1.0f));

    DrawLine3D(to_Vector3(pos), to_Vector3(axis0), RED);
    DrawLine3D(to_Vector3(pos), to_Vector3(axis1), GREEN);
    DrawLine3D(to_Vector3(pos), to_Vector3(axis2), BLUE);
}

void draw_features(const slice1d<float> features, const vec3 pos, const quat rot, const Color color)
{

    vec3 lfoot_pos = quat_mul_vec3(rot, vec3(features( 0), features( 1), features( 2))) + pos;
    vec3 rfoot_pos = quat_mul_vec3(rot, vec3(features( 3), features( 4), features( 5))) + pos;
    vec3 lfoot_vel = quat_mul_vec3(rot, vec3(features( 6), features( 7), features( 8)));
    vec3 rfoot_vel = quat_mul_vec3(rot, vec3(features( 9), features(10), features(11)));
    //vec3 hip_vel   = quat_mul_vec3(rot, vec3(features(12), features(13), features(14)));
    vec3 traj0_pos = quat_mul_vec3(rot, vec3(features(15),         0.0f, features(16))) + pos;
    vec3 traj1_pos = quat_mul_vec3(rot, vec3(features(17),         0.0f, features(18))) + pos;
    vec3 traj2_pos = quat_mul_vec3(rot, vec3(features(19),         0.0f, features(20))) + pos;
    vec3 traj0_dir = quat_mul_vec3(rot, vec3(features(21),         0.0f, features(22)));
    vec3 traj1_dir = quat_mul_vec3(rot, vec3(features(23),         0.0f, features(24)));
    vec3 traj2_dir = quat_mul_vec3(rot, vec3(features(25),         0.0f, features(26)));

    DrawSphereWires(to_Vector3(lfoot_pos), 0.05f, 4, 10, color);
    DrawSphereWires(to_Vector3(rfoot_pos), 0.05f, 4, 10, color);
    DrawSphereWires(to_Vector3(traj0_pos), 0.05f, 4, 10, color);
    DrawSphereWires(to_Vector3(traj1_pos), 0.05f, 4, 10, color);
    DrawSphereWires(to_Vector3(traj2_pos), 0.05f, 4, 10, color);

    DrawLine3D(to_Vector3(lfoot_pos), to_Vector3(lfoot_pos + 0.1f * lfoot_vel), color);
    DrawLine3D(to_Vector3(rfoot_pos), to_Vector3(rfoot_pos + 0.1f * rfoot_vel), color);

    DrawLine3D(to_Vector3(traj0_pos), to_Vector3(traj0_pos + 0.25f * traj0_dir), color);
    DrawLine3D(to_Vector3(traj1_pos), to_Vector3(traj1_pos + 0.25f * traj1_dir), color);
    DrawLine3D(to_Vector3(traj2_pos), to_Vector3(traj2_pos + 0.25f * traj2_dir), color);
}

void draw_trajectory(
    const slice1d<vec3> trajectory_positions,
    const slice1d<quat> trajectory_rotations,
    const Color color)
{

    for (int i = 1; i < trajectory_positions.size; i++)
    {
        DrawSphereWires(to_Vector3(trajectory_positions(i)), 0.05f, 4, 10, color);
        DrawLine3D(to_Vector3(trajectory_positions(i)), to_Vector3(
            trajectory_positions(i) + 0.6f * quat_mul_vec3(trajectory_rotations(i), vec3(0, 0, 1.0f))), color);
        DrawLine3D(to_Vector3(trajectory_positions(i-1)), to_Vector3(trajectory_positions(i)), color);
    }
}

void draw_obstacles(
    const slice1d<vec3> obstacles_positions,
    const slice1d<vec3> obstacles_scales)
{

    for (int i = 0; i < obstacles_positions.size; i++)
    {
        vec3 position = vec3(
            obstacles_positions(i).x,
            obstacles_positions(i).y + 0.5f * obstacles_scales(i).y + 0.01f,
            obstacles_positions(i).z);

        DrawCube(
            to_Vector3(position),
            obstacles_scales(i).x,
            obstacles_scales(i).y,
            obstacles_scales(i).z,
            LIGHTGRAY);

        DrawCubeWires(
            to_Vector3(position),
            obstacles_scales(i).x,
            obstacles_scales(i).y,
            obstacles_scales(i).z,
            GRAY);
    }
}


//--------------------------------------

vec3 clamp_character_position(
    const vec3 character_position,
    const vec3 simulation_position,
    const float max_distance)
{

    // If the character deviates too far from the simulation
    // position we need to clamp it to within the max distance
    if (length(character_position - simulation_position) > max_distance)
    {
        return max_distance *
            normalize(character_position - simulation_position) +
            simulation_position;
    }
    else
    {
        return character_position;
    }
}

quat clamp_character_rotation(
    const quat character_rotation,
    const quat simulation_rotation,
    const float max_angle)
{

    // If the angle between the character rotation and simulation
    // rotation exceeds the threshold we need to clamp it back
    if (quat_angle_between(character_rotation, simulation_rotation) > max_angle)
    {
        // First, find the rotational difference between the two
        quat diff = quat_abs(quat_mul_inv(
            character_rotation, simulation_rotation));

        // We can then decompose it into angle and axis
        float diff_angle; vec3 diff_axis;
        quat_to_angle_axis(diff, diff_angle, diff_axis);

        // We then clamp the angle to within our bounds
        diff_angle = clampf(diff_angle, -max_angle, max_angle);

        // And apply back the clamped rotation
        return quat_mul(
          quat_from_angle_axis(diff_angle, diff_axis), simulation_rotation);
    }
    else
    {
        return character_rotation;
    }
}

//--------------------------------------

void update_callback(void* args)
{
    ((std::function<void()>*)args)->operator()();
}

int main(void)
{
    // Init Window

    const int screen_width = 1280;
    const int screen_height = 720;
    int TheBest = -1, TheCur = -1, TheFound = -1;

    SetConfigFlags(FLAG_VSYNC_HINT);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(screen_width, screen_height, "raylib [data vs code driven displacement]");
    SetTargetFPS(60);

    // Camera

    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 10.0f, 10.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    float camera_azimuth = 0.0f;
    float camera_altitude = 0.4f;
    float camera_distance = 4.0f;

    // Scene Obstacles

    array1d<vec3> obstacles_positions(3);
    array1d<vec3> obstacles_scales(3);

    obstacles_positions(0) = vec3(5.0f, 0.0f, 6.0f);
    obstacles_positions(1) = vec3(-3.0f, 0.0f, -5.0f);
    obstacles_positions(2) = vec3(-8.0f, 0.0f, 3.0f);

    obstacles_scales(0) = vec3(2.0f, 1.0f, 5.0f);
    obstacles_scales(1) = vec3(4.0f, 1.0f, 4.0f);
    obstacles_scales(2) = vec3(2.0f, 1.0f, 2.0f);

    // Ground Plane

    Shader ground_plane_shader = LoadShader("./resources/checkerboard.vs", "./resources/checkerboard.fs");
    Mesh ground_plane_mesh = GenMeshPlane(20.0f, 20.0f, 10, 10);
    Model ground_plane_model = LoadModelFromMesh(ground_plane_mesh);
    ground_plane_model.materials[0].shader = ground_plane_shader;

    // Character

    character character_data;
    character_load(character_data, "./resources/character.bin");

    Shader character_shader = LoadShader("./resources/character.vs", "./resources/character.fs");
    Mesh character_mesh = make_character_mesh(character_data);
    Model character_model = LoadModelFromMesh(character_mesh);
    character_model.materials[0].shader = character_shader;

    // Load Animation Data and build Matching Database

    database db;
    database_load(db, "./resources/database.bin");

    float feature_weight_foot_position = 0.75f;
    float feature_weight_foot_velocity = 1.0f;
    float feature_weight_hip_velocity = 1.0f;
    float feature_weight_trajectory_positions = 1.0f;
    float feature_weight_trajectory_directions = 1.5f;

    database_build_matching_features(
        db,
        feature_weight_foot_position,
        feature_weight_foot_velocity,
        feature_weight_hip_velocity,
        feature_weight_trajectory_positions,
        feature_weight_trajectory_directions);

    database_save_matching_features(db, "./resources/features.bin");

    // Pose & Inertializer Data

    int frame_index = db.range_starts(0);
    float inertialize_blending_halflife = 0.1f;

    array1d<vec3> curr_bone_positions = db.bone_positions(frame_index);
    array1d<vec3> curr_bone_velocities = db.bone_velocities(frame_index);
    array1d<quat> curr_bone_rotations = db.bone_rotations(frame_index);
    array1d<vec3> curr_bone_angular_velocities = db.bone_angular_velocities(frame_index);
    array1d<bool> curr_bone_contacts = db.contact_states(frame_index);

    array1d<vec3> trns_bone_positions = db.bone_positions(frame_index);
    array1d<vec3> trns_bone_velocities = db.bone_velocities(frame_index);
    array1d<quat> trns_bone_rotations = db.bone_rotations(frame_index);
    array1d<vec3> trns_bone_angular_velocities = db.bone_angular_velocities(frame_index);
    array1d<bool> trns_bone_contacts = db.contact_states(frame_index);

    array1d<vec3> bone_positions = db.bone_positions(frame_index);
    array1d<vec3> bone_velocities = db.bone_velocities(frame_index);
    array1d<quat> bone_rotations = db.bone_rotations(frame_index);
    array1d<vec3> bone_angular_velocities = db.bone_angular_velocities(frame_index);

    array1d<vec3> bone_offset_positions(db.nbones());
    array1d<vec3> bone_offset_velocities(db.nbones());
    array1d<quat> bone_offset_rotations(db.nbones());
    array1d<vec3> bone_offset_angular_velocities(db.nbones());

    array1d<vec3> global_bone_positions(db.nbones());
    array1d<vec3> global_bone_velocities(db.nbones());
    array1d<quat> global_bone_rotations(db.nbones());
    array1d<vec3> global_bone_angular_velocities(db.nbones());
    array1d<bool> global_bone_computed(db.nbones());

    vec3 transition_src_position;
    quat transition_src_rotation;
    vec3 transition_dst_position;
    quat transition_dst_rotation;

    inertialize_pose_reset(
        bone_offset_positions,
        bone_offset_velocities,
        bone_offset_rotations,
        bone_offset_angular_velocities,
        transition_src_position,
        transition_src_rotation,
        transition_dst_position,
        transition_dst_rotation,
        bone_positions(0),
        bone_rotations(0));

    inertialize_pose_update(
        bone_positions,
        bone_velocities,
        bone_rotations,
        bone_angular_velocities,
        bone_offset_positions,
        bone_offset_velocities,
        bone_offset_rotations,
        bone_offset_angular_velocities,
        db.bone_positions(frame_index),
        db.bone_velocities(frame_index),
        db.bone_rotations(frame_index),
        db.bone_angular_velocities(frame_index),
        transition_src_position,
        transition_src_rotation,
        transition_dst_position,
        transition_dst_rotation,
        inertialize_blending_halflife,
        0.0f);

    // Trajectory & Gameplay Data

    float search_time = 0.1f;
    float search_timer = search_time;
    float force_search_timer = search_time;

    vec3 desired_velocity;
    vec3 desired_velocity_change_curr;
    vec3 desired_velocity_change_prev;
    float desired_velocity_change_threshold = 50.0;

    quat desired_rotation;
    vec3 desired_rotation_change_curr;
    vec3 desired_rotation_change_prev;
    float desired_rotation_change_threshold = 50.0;

    float desired_gait = 0.0f;
    float desired_gait_velocity = 0.0f;

    vec3 simulation_position;
    vec3 simulation_velocity;
    vec3 simulation_acceleration;
    quat simulation_rotation;
    vec3 simulation_angular_velocity;

    float simulation_velocity_halflife = 0.27f;
    float simulation_rotation_halflife = 0.27f;

    // All speeds in m/s
    float simulation_run_fwrd_speed = 4.0f;
    float simulation_run_side_speed = 3.0f;
    float simulation_run_back_speed = 2.5f;

    float simulation_walk_fwrd_speed = 1.75f;
    float simulation_walk_side_speed = 1.5f;
    float simulation_walk_back_speed = 1.25f;

    array1d<vec3> trajectory_desired_velocities(4);
    array1d<quat> trajectory_desired_rotations(4);
    array1d<vec3> trajectory_positions(4);
    array1d<vec3> trajectory_velocities(4);
    array1d<vec3> trajectory_accelerations(4);
    array1d<quat> trajectory_rotations(4);
    array1d<vec3> trajectory_angular_velocities(4);



    // Clamping

    bool clamping_enabled = true;
    float clamping_max_distance = 0.15f;
    float clamping_max_angle = 0.5f * PIf;




    array1d<vec3> adjusted_bone_positions = bone_positions;
    array1d<quat> adjusted_bone_rotations = bone_rotations;

    // Learned Motion Matching

    bool brute_force = false;


    array1d<float> features_proj = db.features(frame_index);
    array1d<float> features_curr = db.features(frame_index);
    array1d<float> latent_proj(32); latent_proj.zero();
    array1d<float> latent_curr(32); latent_curr.zero();

    // Go

    float dt = 1.0f / 60.0f;

    auto update_func = [&]()
    {
        // Get gamepad stick states
        vec3 gamepadstick_left = gamepad_get_stick(GAMEPAD_STICK_LEFT);
        vec3 gamepadstick_right = gamepad_get_stick(GAMEPAD_STICK_RIGHT);

        // Get if strafe is desired
        bool desired_strafe = desired_strafe_update();

        // Get the desired gait (walk / run)
        desired_gait_update(
            desired_gait,
            desired_gait_velocity,
            dt);

        // Get the desired simulation speeds based on the gait
        float simulation_fwrd_speed = lerpf(simulation_run_fwrd_speed, simulation_walk_fwrd_speed, desired_gait);
        float simulation_side_speed = lerpf(simulation_run_side_speed, simulation_walk_side_speed, desired_gait);
        float simulation_back_speed = lerpf(simulation_run_back_speed, simulation_walk_back_speed, desired_gait);

        // Get the desired velocity
        vec3 desired_velocity_curr = desired_velocity_update(
            gamepadstick_left,
            camera_azimuth,
            simulation_rotation,
            simulation_fwrd_speed,
            simulation_side_speed,
            simulation_back_speed);

        // Get the desired rotation/direction
        quat desired_rotation_curr = desired_rotation_update(
            desired_rotation,
            gamepadstick_left,
            gamepadstick_right,
            camera_azimuth,
            desired_strafe,
            desired_velocity_curr);

        // Check if we should force a search because input changed quickly
        desired_velocity_change_prev = desired_velocity_change_curr;
        desired_velocity_change_curr =  (desired_velocity_curr - desired_velocity) / dt;
        desired_velocity = desired_velocity_curr;

        desired_rotation_change_prev = desired_rotation_change_curr;
        desired_rotation_change_curr = quat_to_scaled_angle_axis(quat_abs(quat_mul_inv(desired_rotation_curr, desired_rotation))) / dt;
        desired_rotation =  desired_rotation_curr;

        bool force_search = false;

        if (force_search_timer <= 0.0f && (
            (length(desired_velocity_change_prev) >= desired_velocity_change_threshold &&
             length(desired_velocity_change_curr)  < desired_velocity_change_threshold)
        ||  (length(desired_rotation_change_prev) >= desired_rotation_change_threshold &&
             length(desired_rotation_change_curr)  < desired_rotation_change_threshold)))
        {
            force_search = true;
            force_search_timer = search_time;
        }
        else if (force_search_timer > 0)
        {
            force_search_timer -= dt;
        }

        // Predict Future Trajectory

        trajectory_desired_rotations_predict(
          trajectory_desired_rotations,
          trajectory_desired_velocities,
          desired_rotation,
          camera_azimuth,
          gamepadstick_left,
          gamepadstick_right,
          desired_strafe,
          20.0f * dt);

        trajectory_rotations_predict(
            trajectory_rotations,
            trajectory_angular_velocities,
            simulation_rotation,
            simulation_angular_velocity,
            trajectory_desired_rotations,
            simulation_rotation_halflife,
            20.0f * dt);

        trajectory_desired_velocities_predict(
          trajectory_desired_velocities,
          trajectory_rotations,
          desired_velocity,
          camera_azimuth,
          gamepadstick_left,
          gamepadstick_right,
          desired_strafe,
          simulation_fwrd_speed,
          simulation_side_speed,
          simulation_back_speed,
          20.0f * dt);

        trajectory_positions_predict(
            trajectory_positions,
            trajectory_velocities,
            trajectory_accelerations,
            simulation_position,
            simulation_velocity,
            simulation_acceleration,
            trajectory_desired_velocities,
            simulation_velocity_halflife,
            20.0f * dt,
            obstacles_positions,
            obstacles_scales);

        // Make query vector for search.
        // In theory this only needs to be done when a search is
        // actually required however for visualization purposes it
        // can be nice to do it every frame
        array1d<float> query(db.nfeatures());

        // Compute the features of the query vector

        slice1d<float> query_features = db.features(frame_index);

        int offset = 0;
        query_copy_denormalized_feature(query, offset, 3, query_features, db.features_offset, db.features_scale); // Left Foot Position
        query_copy_denormalized_feature(query, offset, 3, query_features, db.features_offset, db.features_scale); // Right Foot Position
        query_copy_denormalized_feature(query, offset, 3, query_features, db.features_offset, db.features_scale); // Left Foot Velocity
        query_copy_denormalized_feature(query, offset, 3, query_features, db.features_offset, db.features_scale); // Right Foot Velocity
        query_copy_denormalized_feature(query, offset, 3, query_features, db.features_offset, db.features_scale); // Hip Velocity
        query_compute_trajectory_position_feature(query, offset, bone_positions(0), bone_rotations(0), trajectory_positions);
        query_compute_trajectory_direction_feature(query, offset, bone_rotations(0), trajectory_rotations);

        assert(offset == db.nfeatures());

        // Check if we reached the end of the current anim
        bool end_of_anim = database_trajectory_index_clamp(db, frame_index, 1) == frame_index;

        // Do we need to search?
        if (force_search || search_timer <= 0.0f || end_of_anim)
        {

                // Search

                int best_index = end_of_anim ? -1 : frame_index;
                float best_cost = FLT_MAX;

                database_search(
                    best_index,
                    best_cost,
                    db,
                    query,
                    brute_force);
                TheFound = best_index;
                TheCur = frame_index;

                // Transition if better frame found

                if (best_index != frame_index)
                {
                    TheBest = best_index;

                    inertialize_pose_transition(
                        bone_offset_positions,
                        bone_offset_velocities,
                        bone_offset_rotations,
                        bone_offset_angular_velocities,
                        transition_src_position,
                        transition_src_rotation,
                        transition_dst_position,
                        transition_dst_rotation,
                        bone_positions(0),
                        bone_velocities(0),
                        bone_rotations(0),
                        bone_angular_velocities(0),
                        curr_bone_positions,
                        curr_bone_velocities,
                        curr_bone_rotations,
                        curr_bone_angular_velocities,
                        db.bone_positions(best_index),
                        db.bone_velocities(best_index),
                        db.bone_rotations(best_index),
                        db.bone_angular_velocities(best_index));

                    frame_index = best_index;
                }

            // Reset search timer
            search_timer = search_time;
        }

        // Tick down search timer
        search_timer -= dt;

        // Tick frame
        frame_index++; // Assumes dt is fixed to 60fps

        // Look-up Next Pose
        curr_bone_positions = db.bone_positions(frame_index);
        curr_bone_velocities = db.bone_velocities(frame_index);
        curr_bone_rotations = db.bone_rotations(frame_index);
        curr_bone_angular_velocities = db.bone_angular_velocities(frame_index);
        curr_bone_contacts = db.contact_states(frame_index);

        // Update inertializer

        inertialize_pose_update(
            bone_positions,
            bone_velocities,
            bone_rotations,
            bone_angular_velocities,
            bone_offset_positions,
            bone_offset_velocities,
            bone_offset_rotations,
            bone_offset_angular_velocities,
            curr_bone_positions,
            curr_bone_velocities,
            curr_bone_rotations,
            curr_bone_angular_velocities,
            transition_src_position,
            transition_src_rotation,
            transition_dst_position,
            transition_dst_rotation,
            inertialize_blending_halflife,
            dt);

        // Update Simulation

        vec3 simulation_position_prev = simulation_position;

        simulation_positions_update(
            simulation_position,
            simulation_velocity,
            simulation_acceleration,
            desired_velocity,
            simulation_velocity_halflife,
            dt,
            obstacles_positions,
            obstacles_scales);

        simulation_rotations_update(
            simulation_rotation,
            simulation_angular_velocity,
            desired_rotation,
            simulation_rotation_halflife,
            dt);

        // Clamping

        if (clamping_enabled)
        {
            vec3 adjusted_position = bone_positions(0);
            quat adjusted_rotation = bone_rotations(0);

            adjusted_position = clamp_character_position(
                adjusted_position,
                simulation_position,
                clamping_max_distance);

            adjusted_rotation = clamp_character_rotation(
                adjusted_rotation,
                simulation_rotation,
                clamping_max_angle);

            inertialize_root_adjust(
                bone_offset_positions(0),
                transition_src_position,
                transition_src_rotation,
                transition_dst_position,
                transition_dst_rotation,
                bone_positions(0),
                bone_rotations(0),
                adjusted_position,
                adjusted_rotation);
        }

        // Contact fixup with foot locking and IK

        adjusted_bone_positions = bone_positions;
        adjusted_bone_rotations = bone_rotations;


        // Full pass of forward kinematics to compute
        // all bone positions and rotations in the world
        // space ready for rendering

        forward_kinematics_full(
            global_bone_positions,
            global_bone_rotations,
            adjusted_bone_positions,
            adjusted_bone_rotations,
            db.bone_parents);

        // Update camera

        orbit_camera_update(
            camera,
            camera_azimuth,
            camera_altitude,
            camera_distance,
            bone_positions(0) + vec3(0, 1, 0),
            // simulation_position + vec3(0, 1, 0),
            gamepadstick_right,
            desired_strafe,
            dt);

        // Render

        BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);

        // Draw Simulation Object

        DrawCylinderWires(to_Vector3(simulation_position), 0.6f, 0.6f, 0.001f, 17, ORANGE);
        DrawSphereWires(to_Vector3(simulation_position), 0.05f, 4, 10, ORANGE);
        DrawLine3D(to_Vector3(simulation_position), to_Vector3(
            simulation_position + 0.6f * quat_mul_vec3(simulation_rotation, vec3(0.0f, 0.0f, 1.0f))), ORANGE);

        // Draw Clamping Radius/Angles

        if (clamping_enabled)
        {
            DrawCylinderWires(
                to_Vector3(simulation_position),
                clamping_max_distance,
                clamping_max_distance,
                0.001f, 17, SKYBLUE);

            quat rotation_clamp_0 = quat_mul(quat_from_angle_axis(+clamping_max_angle, vec3(0.0f, 1.0f, 0.0f)), simulation_rotation);
            quat rotation_clamp_1 = quat_mul(quat_from_angle_axis(-clamping_max_angle, vec3(0.0f, 1.0f, 0.0f)), simulation_rotation);

            vec3 rotation_clamp_0_dir = simulation_position + 0.6f * quat_mul_vec3(rotation_clamp_0, vec3(0.0f, 0.0f, 1.0f));
            vec3 rotation_clamp_1_dir = simulation_position + 0.6f * quat_mul_vec3(rotation_clamp_1, vec3(0.0f, 0.0f, 1.0f));

            DrawLine3D(to_Vector3(simulation_position), to_Vector3(rotation_clamp_0_dir), SKYBLUE);
            DrawLine3D(to_Vector3(simulation_position), to_Vector3(rotation_clamp_1_dir), SKYBLUE);
        }



        draw_trajectory(
            trajectory_positions,
            trajectory_rotations,
            ORANGE);

        draw_obstacles(
            obstacles_positions,
            obstacles_scales);

        deform_character_mesh(
            character_mesh,
            character_data,
            global_bone_positions,
            global_bone_rotations,
            db.bone_parents);

        DrawModel(character_model, (Vector3){0.0f, 0.0f, 0.0f}, 1.0f, RAYWHITE);

        // Draw matched features

        array1d<float> current_features = db.features(frame_index);
        denormalize_features(current_features, db.features_offset, db.features_scale);
        draw_features(current_features, bone_positions(0), bone_rotations(0), MAROON);

        // Draw Simuation Bone

        DrawSphereWires(to_Vector3(bone_positions(0)), 0.05f, 4, 10, MAROON);
        DrawLine3D(to_Vector3(bone_positions(0)), to_Vector3(
            bone_positions(0) + 0.6f * quat_mul_vec3(bone_rotations(0), vec3(0.0f, 0.0f, 1.0f))), MAROON);

        // Draw Ground Plane

        DrawModel(ground_plane_model, (Vector3){0.0f, -0.01f, 0.0f}, 1.0f, WHITE);
        DrawGrid(20, 1.0f);
        draw_axis(vec3(), quat());

        EndMode3D();

        // UI

        //---------

        float ui_sim_hei = 20;

        GuiGroupBox((Rectangle){ 970, ui_sim_hei, 290, 250 }, "simulation object");

        simulation_velocity_halflife = GuiSliderBar(
            (Rectangle){ 1100, ui_sim_hei + 10, 120, 20 },
            "velocity halflife",
            TextFormat("%5.3f", simulation_velocity_halflife),
            simulation_velocity_halflife, 0.0f, 0.5f);

        simulation_rotation_halflife = GuiSliderBar(
            (Rectangle){ 1100, ui_sim_hei + 40, 120, 20 },
            "rotation halflife",
            TextFormat("%5.3f", simulation_rotation_halflife),
            simulation_rotation_halflife, 0.0f, 0.5f);

        simulation_run_fwrd_speed = GuiSliderBar(
            (Rectangle){ 1100, ui_sim_hei + 70, 120, 20 },
            "run forward speed",
            TextFormat("%5.3f", simulation_run_fwrd_speed),
            simulation_run_fwrd_speed, 0.0f, 10.0f);

        simulation_run_side_speed = GuiSliderBar(
            (Rectangle){ 1100, ui_sim_hei + 100, 120, 20 },
            "run sideways speed",
            TextFormat("%5.3f", simulation_run_side_speed),
            simulation_run_side_speed, 0.0f, 10.0f);

        simulation_run_back_speed = GuiSliderBar(
            (Rectangle){ 1100, ui_sim_hei + 130, 120, 20 },
            "run backwards speed",
            TextFormat("%5.3f", simulation_run_back_speed),
            simulation_run_back_speed, 0.0f, 10.0f);

        simulation_walk_fwrd_speed = GuiSliderBar(
            (Rectangle){ 1100, ui_sim_hei + 160, 120, 20 },
            "walk forward speed",
            TextFormat("%5.3f", simulation_walk_fwrd_speed),
            simulation_walk_fwrd_speed, 0.0f, 5.0f);

        simulation_walk_side_speed = GuiSliderBar(
            (Rectangle){ 1100, ui_sim_hei + 190, 120, 20 },
            "walk sideways speed",
            TextFormat("%5.3f", simulation_walk_side_speed),
            simulation_walk_side_speed, 0.0f, 5.0f);

        simulation_walk_back_speed = GuiSliderBar(
            (Rectangle){ 1100, ui_sim_hei + 220, 120, 20 },
            "walk backwards speed",
            TextFormat("%5.3f", simulation_walk_back_speed),
            simulation_walk_back_speed, 0.0f, 5.0f);

        //---------

        float ui_inert_hei = 280;

        GuiGroupBox((Rectangle){ 970, ui_inert_hei, 290, 40 }, "inertiaization blending");

        inertialize_blending_halflife = GuiSliderBar(
            (Rectangle){ 1100, ui_inert_hei + 10, 120, 20 },
            "halflife",
            TextFormat("%5.3f", inertialize_blending_halflife),
            inertialize_blending_halflife, 0.0f, 0.3f);
        //---------

        float ui_ctrl_hei = 380;

        GuiGroupBox((Rectangle){ 1010, ui_ctrl_hei, 250, 140 }, "controls");

        GuiLabel((Rectangle){ 1030, ui_ctrl_hei +  10, 200, 20 }, "WASD - Move");

        //---------

        GuiGroupBox((Rectangle){ 20, 20, 290, 190 }, "feature weights");

        feature_weight_foot_position = GuiSliderBar(
            (Rectangle){ 150, 30, 120, 20 },
            "foot position",
            TextFormat("%5.3f", feature_weight_foot_position),
            feature_weight_foot_position, 0.001f, 3.0f);

        feature_weight_foot_velocity = GuiSliderBar(
            (Rectangle){ 150, 60, 120, 20 },
            "foot velocity",
            TextFormat("%5.3f", feature_weight_foot_velocity),
            feature_weight_foot_velocity, 0.001f, 3.0f);

        feature_weight_hip_velocity = GuiSliderBar(
            (Rectangle){ 150, 90, 120, 20 },
            "hip velocity",
            TextFormat("%5.3f", feature_weight_hip_velocity),
            feature_weight_hip_velocity, 0.001f, 3.0f);

        feature_weight_trajectory_positions = GuiSliderBar(
            (Rectangle){ 150, 120, 120, 20 },
            "trajectory positions",
            TextFormat("%5.3f", feature_weight_trajectory_positions),
            feature_weight_trajectory_positions, 0.001f, 3.0f);

        feature_weight_trajectory_directions = GuiSliderBar(
            (Rectangle){ 150, 150, 120, 20 },
            "trajectory directions",
            TextFormat("%5.3f", feature_weight_trajectory_directions),
            feature_weight_trajectory_directions, 0.001f, 3.0f);

        if (GuiButton((Rectangle){ 150, 180, 120, 20 }, "rebuild database"))
        {
            database_build_matching_features(
                db,
                feature_weight_foot_position,
                feature_weight_foot_velocity,
                feature_weight_hip_velocity,
                feature_weight_trajectory_positions,
                feature_weight_trajectory_directions);
        }


        //---------

        float ui_clamp_hei = 440;

        GuiGroupBox((Rectangle){ 20, ui_clamp_hei, 290, 100 }, "clamping");

        clamping_enabled = GuiCheckBox(
            (Rectangle){ 50, ui_clamp_hei + 10, 20, 20 },
            "enabled",
            clamping_enabled);

        clamping_max_distance = GuiSliderBar(
            (Rectangle){ 150, ui_clamp_hei + 40, 120, 20 },
            "distance",
            TextFormat("%5.3f", clamping_max_distance),
            clamping_max_distance, 0.0f, 0.5f);

        clamping_max_angle = GuiSliderBar(
            (Rectangle){ 150, ui_clamp_hei + 70, 120, 20 },
            "angle",
            TextFormat("%5.3f", clamping_max_angle),
            clamping_max_angle, 0.0f, PIf);

        //---------

        float ui_ik_hei = 550;

        //---------

        EndDrawing();

    };

#if defined(PLATFORM_WEB)
    std::function<void()> u{update_func};
    emscripten_set_main_loop_arg(update_callback, &u, 0, 1);
#else
    while (!WindowShouldClose())
    {
        update_func();
    }
#endif

    // Unload stuff and finish
    UnloadModel(character_model);
    UnloadModel(ground_plane_model);
    UnloadShader(character_shader);
    UnloadShader(ground_plane_shader);

    CloseWindow();

    return 0;
}
