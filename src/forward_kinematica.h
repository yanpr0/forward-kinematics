#pragma once

#include <functional>

#include "common.h"
#include "vec.h"
#include "quat.h"
#include "array.h"

//--------------------------------------

// Compute forward kinematics for all joints
void forward_kinematics_full(
    slice1d<vec3> global_bone_positions,
    slice1d<quat> global_bone_rotations,
    const slice1d<vec3> local_bone_positions,
    const slice1d<quat> local_bone_rotations,
    const slice1d<int> bone_parents)
{
    array1d<bool> is_set(global_bone_positions.size);
    is_set.zero();
    global_bone_positions(0) = local_bone_positions(0);
    global_bone_rotations(0) = quat_normalize(local_bone_rotations(0));
    is_set(0) = true;

    std::function<void(int)> set = [&](int bone)
    {
        int parent_bone = bone_parents(bone);
        if (!is_set(parent_bone))
        {
            set(parent_bone);
        }
        global_bone_positions(bone) = global_bone_positions(parent_bone) + quat_mul_vec3(global_bone_rotations(parent_bone), local_bone_positions(bone));
        //global_bone_positions(bone) = global_bone_positions(parent_bone) + local_bone_positions(bone); // uncomment for flex
        global_bone_rotations(bone) = quat_mul(global_bone_rotations(parent_bone), local_bone_rotations(bone));
        //global_bone_rotations(bone) = quat_mul(local_bone_rotations(bone), global_bone_rotations(parent_bone)); // uncomment for horror
        global_bone_rotations(bone) = quat_normalize(global_bone_rotations(bone));
        is_set(bone) = true;
    };

    for (int i = 1; i < is_set.size; ++i)
    {
        if (!is_set(i))
        {
            set(i);
        }
    }
}

// Here we using a simple recursive version of forward kinematics for calculation only one bone
void forward_kinematics(
    vec3& bone_position,
    quat& bone_rotation,
    const slice1d<vec3> bone_positions,
    const slice1d<quat> bone_rotations,
    const slice1d<int> bone_parents,
    const int bone)
{
    if (bone == 0)
    {
        bone_position = bone_positions(0);
        bone_rotation = bone_rotations(0);
        return;
    }

    vec3 parent_pos;
    quat parent_rot;
    forward_kinematics(parent_pos, parent_rot, bone_positions, bone_rotations, bone_parents, bone_parents(bone));

    bone_position = parent_pos + quat_mul_vec3(parent_rot, bone_positions(bone));
    bone_rotation = quat_mul(parent_rot, bone_rotations(bone));
    bone_rotation = quat_normalize(bone_rotation);
}

// Forward kinematics but also compute the velocities
void forward_kinematics_velocity(
    vec3& bone_position,
    vec3& bone_velocity,
    quat& bone_rotation,
    vec3& bone_angular_velocity,
    const slice1d<vec3> bone_positions,
    const slice1d<vec3> bone_velocities,
    const slice1d<quat> bone_rotations,
    const slice1d<vec3> bone_angular_velocities,
    const slice1d<int> bone_parents,
    const int bone)
{
    if (bone == 0)
    {
        bone_position = bone_positions(0);
        bone_velocity = bone_velocities(0);
        bone_rotation = bone_rotations(0);
        bone_angular_velocity = bone_angular_velocities(0);
        return;
    }

    vec3 parent_pos;
    vec3 parent_vel;
    quat parent_rot;
    vec3 parent_ang_vel;
    forward_kinematics_velocity(
            parent_pos, parent_vel, parent_rot, parent_ang_vel,
            bone_positions, bone_velocities, bone_rotations, bone_angular_velocities,
            bone_parents, bone_parents(bone));

    bone_position = parent_pos + quat_mul_vec3(parent_rot, bone_positions(bone));
    bone_velocity = parent_vel + quat_mul_vec3(parent_rot, bone_velocities(bone)) + cross(parent_ang_vel, quat_mul_vec3(parent_rot, bone_positions(bone)));
    bone_rotation = quat_mul(parent_rot, bone_rotations(bone));
    bone_rotation = quat_normalize(bone_rotation);
    bone_angular_velocity = parent_ang_vel + quat_mul_vec3(parent_rot, bone_angular_velocities(bone));
}

