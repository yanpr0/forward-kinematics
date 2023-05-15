#pragma once

#include "common.h"
#include "vec.h"
#include "quat.h"
#include "array.h"
#include "forward_kinematica.h"

#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <math.h>
//--------------------------------------

enum
{
    BOUND_SM_SIZE = 16,
    BOUND_LR_SIZE = 64,
};

struct database
{
    array2d<vec3> bone_positions;
    array2d<vec3> bone_velocities;
    array2d<quat> bone_rotations;
    array2d<vec3> bone_angular_velocities;
    array1d<int> bone_parents;

    array1d<int> range_starts;
    array1d<int> range_stops;

    array2d<float> features;
    array1d<float> features_offset;
    array1d<float> features_scale;

    array2d<bool> contact_states;

    array2d<float> bound_sm_min;
    array2d<float> bound_sm_max;
    array2d<float> bound_lr_min;
    array2d<float> bound_lr_max;

    int nframes() const { return bone_positions.rows; }
    int nbones() const { return bone_positions.cols; }
    int nranges() const { return range_starts.size; }
    int nfeatures() const { return features.cols; }
    int ncontacts() const { return contact_states.cols; }
};

void database_load(database& db, const char* filename)
{
    FILE* f = fopen(filename, "rb");
    assert(f != NULL);

    array2d_read(db.bone_positions, f);
    array2d_read(db.bone_velocities, f);
    array2d_read(db.bone_rotations, f);
    array2d_read(db.bone_angular_velocities, f);
    array1d_read(db.bone_parents, f);

    array1d_read(db.range_starts, f);
    array1d_read(db.range_stops, f);

    array2d_read(db.contact_states, f);

    fclose(f);
    printf("Frames count %d\n", db.nframes());//53500
}

void database_save_matching_features(const database& db, const char* filename)
{
    FILE* f = fopen(filename, "wb");
    assert(f != NULL);

    array2d_write(db.features, f);
    array1d_write(db.features_offset, f);
    array1d_write(db.features_scale, f);

    fclose(f);
}

// When we add an offset to a frame in the database there is a chance
// it will go out of the relevant range so here we can clamp it to
// the last frame of that range.
int database_trajectory_index_clamp(database& db, int frame, int offset)
{
    for (int i = 0; i < db.nranges(); i++)
    {
        if (frame >= db.range_starts(i) && frame < db.range_stops(i))
        {
            return clamp(frame + offset, db.range_starts(i), db.range_stops(i) - 1);
        }
    }

    assert(false);
    return -1;
}

//--------------------------------------

void normalize_feature(
    slice2d<float> features,
    slice1d<float> features_offset,
    slice1d<float> features_scale,
    const int offset,
    const int size,
    const float weight = 1.0f)
{
    // First compute what is essentially the mean
    // value for each feature dimension
    for (int j = 0; j < size; j++)
    {
        features_offset(offset + j) = 0.0f;
    }

    for (int i = 0; i < features.rows; i++)
    {
        for (int j = 0; j < size; j++)
        {
            features_offset(offset + j) += features(i, offset + j) / features.rows;
        }
    }

    // Now compute the variance of each feature dimension
    array1d<float> vars(size);
    vars.zero();

    for (int i = 0; i < features.rows; i++)
    {
        for (int j = 0; j < size; j++)
        {
            vars(j) += squaref(features(i, offset + j) - features_offset(offset + j)) / features.rows;
        }
    }

    // We compute the overall std of the feature as the average
    // std across all dimensions
    float std = 0.0f;
    for (int j = 0; j < size; j++)
    {
        std += sqrtf(vars(j)) / size;
    }

    // Features with no variation can have zero std which is
    // almost always a bug.
    assert(std > 0.0);

    // The scale of a feature is just the std divided by the weight
    for (int j = 0; j < size; j++)
    {
        features_scale(offset + j) = std / weight;
    }

    // Using the offset and scale we can then normalize the features
    for (int i = 0; i < features.rows; i++)
    {
        for (int j = 0; j < size; j++)
        {
            features(i, offset + j) = (features(i, offset + j) - features_offset(offset + j)) / features_scale(offset + j);
        }
    }
}

void denormalize_features(
    slice1d<float> features,
    const slice1d<float> features_offset,
    const slice1d<float> features_scale)
{
    for (int i = 0; i < features.size; i++)
    {
        features(i) = (features(i) * features_scale(i)) + features_offset(i);
    }
}



//--------------------------------------

// Compute a feature for the position of a bone relative to the simulation/root bone
void compute_bone_position_feature(database& db, int& offset, int bone, float weight = 1.0f)
{
    for (int i = 0; i < db.nframes(); i++)
    {
        vec3 bone_position;
        quat bone_rotation;

        forward_kinematics(
            bone_position,
            bone_rotation,
            db.bone_positions(i),
            db.bone_rotations(i),
            db.bone_parents,
            bone);

        bone_position = quat_mul_vec3(quat_inv(db.bone_rotations(i, 0)), bone_position - db.bone_positions(i, 0));

        db.features(i, offset + 0) = bone_position.x;
        db.features(i, offset + 1) = bone_position.y;
        db.features(i, offset + 2) = bone_position.z;
    }

    normalize_feature(db.features, db.features_offset, db.features_scale, offset, 3, weight);

    offset += 3;
}

// Similar but for a bone's velocity
void compute_bone_velocity_feature(database& db, int& offset, int bone, float weight = 1.0f)
{
    for (int i = 0; i < db.nframes(); i++)
    {
        vec3 bone_position;
        vec3 bone_velocity;
        quat bone_rotation;
        vec3 bone_angular_velocity;

        forward_kinematics_velocity(
            bone_position,
            bone_velocity,
            bone_rotation,
            bone_angular_velocity,
            db.bone_positions(i),
            db.bone_velocities(i),
            db.bone_rotations(i),
            db.bone_angular_velocities(i),
            db.bone_parents,
            bone);

        bone_velocity = quat_mul_vec3(quat_inv(db.bone_rotations(i, 0)), bone_velocity);

        db.features(i, offset + 0) = bone_velocity.x;
        db.features(i, offset + 1) = bone_velocity.y;
        db.features(i, offset + 2) = bone_velocity.z;
    }

    normalize_feature(db.features, db.features_offset, db.features_scale, offset, 3, weight);

    offset += 3;
}

// Compute the trajectory at 20, 40, and 60 frames in the future
void compute_trajectory_position_feature(database& db, int& offset, float weight = 1.0f)
{
    for (int i = 0; i < db.nframes(); i++)
    {
        int t0 = database_trajectory_index_clamp(db, i, 20);
        int t1 = database_trajectory_index_clamp(db, i, 40);
        int t2 = database_trajectory_index_clamp(db, i, 60);

        vec3 trajectory_pos0 = quat_mul_vec3(quat_inv(db.bone_rotations(i, 0)), db.bone_positions(t0, 0) - db.bone_positions(i, 0));
        vec3 trajectory_pos1 = quat_mul_vec3(quat_inv(db.bone_rotations(i, 0)), db.bone_positions(t1, 0) - db.bone_positions(i, 0));
        vec3 trajectory_pos2 = quat_mul_vec3(quat_inv(db.bone_rotations(i, 0)), db.bone_positions(t2, 0) - db.bone_positions(i, 0));

        db.features(i, offset + 0) = trajectory_pos0.x;
        db.features(i, offset + 1) = trajectory_pos0.z;
        db.features(i, offset + 2) = trajectory_pos1.x;
        db.features(i, offset + 3) = trajectory_pos1.z;
        db.features(i, offset + 4) = trajectory_pos2.x;
        db.features(i, offset + 5) = trajectory_pos2.z;
    }

    normalize_feature(db.features, db.features_offset, db.features_scale, offset, 6, weight);

    offset += 6;
}

// Same for direction
void compute_trajectory_direction_feature(database& db, int& offset, float weight = 1.0f)
{
    for (int i = 0; i < db.nframes(); i++)
    {
        int t0 = database_trajectory_index_clamp(db, i, 20);
        int t1 = database_trajectory_index_clamp(db, i, 40);
        int t2 = database_trajectory_index_clamp(db, i, 60);

        vec3 trajectory_dir0 = quat_mul_vec3(quat_inv(db.bone_rotations(i, 0)), quat_mul_vec3(db.bone_rotations(t0, 0), vec3(0, 0, 1)));
        vec3 trajectory_dir1 = quat_mul_vec3(quat_inv(db.bone_rotations(i, 0)), quat_mul_vec3(db.bone_rotations(t1, 0), vec3(0, 0, 1)));
        vec3 trajectory_dir2 = quat_mul_vec3(quat_inv(db.bone_rotations(i, 0)), quat_mul_vec3(db.bone_rotations(t2, 0), vec3(0, 0, 1)));

        db.features(i, offset + 0) = trajectory_dir0.x;
        db.features(i, offset + 1) = trajectory_dir0.z;
        db.features(i, offset + 2) = trajectory_dir1.x;
        db.features(i, offset + 3) = trajectory_dir1.z;
        db.features(i, offset + 4) = trajectory_dir2.x;
        db.features(i, offset + 5) = trajectory_dir2.z;
    }

    normalize_feature(db.features, db.features_offset, db.features_scale, offset, 6, weight);

    offset += 6;
}

// Build the Motion Matching search acceleration structure. Here we
// just use axis aligned bounding boxes regularly spaced at BOUND_SM_SIZE
// and BOUND_LR_SIZE frames
void database_build_bounds(database& db)
{
    int nbound_sm = ((db.nframes() + BOUND_SM_SIZE - 1) / BOUND_SM_SIZE);
    int nbound_lr = ((db.nframes() + BOUND_LR_SIZE - 1) / BOUND_LR_SIZE);

    db.bound_sm_min.resize(nbound_sm, db.nfeatures());
    db.bound_sm_max.resize(nbound_sm, db.nfeatures());
    db.bound_lr_min.resize(nbound_lr, db.nfeatures());
    db.bound_lr_max.resize(nbound_lr, db.nfeatures());

    db.bound_sm_min.set(FLT_MAX);
    db.bound_sm_max.set(FLT_MIN);
    db.bound_lr_min.set(FLT_MAX);
    db.bound_lr_max.set(FLT_MIN);

    for (int i = 0; i < db.nframes(); i++)
    {
        int i_sm = i / BOUND_SM_SIZE;
        int i_lr = i / BOUND_LR_SIZE;

        for (int j = 0; j < db.nfeatures(); j++)
        {
            db.bound_sm_min(i_sm, j) = minf(db.bound_sm_min(i_sm, j), db.features(i, j));
            db.bound_sm_max(i_sm, j) = maxf(db.bound_sm_max(i_sm, j), db.features(i, j));
            db.bound_lr_min(i_lr, j) = minf(db.bound_lr_min(i_lr, j), db.features(i, j));
            db.bound_lr_max(i_lr, j) = maxf(db.bound_lr_max(i_lr, j), db.features(i, j));
        }
    }
}

// Build all motion matching features and acceleration structure
void database_build_matching_features(
    database& db,
    const float feature_weight_foot_position,
    const float feature_weight_foot_velocity,
    const float feature_weight_hip_velocity,
    const float feature_weight_trajectory_positions,
    const float feature_weight_trajectory_directions)
{
    int nfeatures =
        3 + // Left Foot Position
        3 + // Right Foot Position
        3 + // Left Foot Velocity
        3 + // Right Foot Velocity
        3 + // Hip Velocity
        6 + // Trajectory Positions 2D
        6 ; // Trajectory Directions 2D

    db.features.resize(db.nframes(), nfeatures);
    db.features_offset.resize(nfeatures);
    db.features_scale.resize(nfeatures);

    int offset = 0;
    compute_bone_position_feature(db, offset, Bone_LeftFoot, feature_weight_foot_position);
    compute_bone_position_feature(db, offset, Bone_RightFoot, feature_weight_foot_position);
    compute_bone_velocity_feature(db, offset, Bone_LeftFoot, feature_weight_foot_velocity);
    compute_bone_velocity_feature(db, offset, Bone_RightFoot, feature_weight_foot_velocity);
    compute_bone_velocity_feature(db, offset, Bone_Hips, feature_weight_hip_velocity);
    compute_trajectory_position_feature(db, offset, feature_weight_trajectory_positions);
    compute_trajectory_direction_feature(db, offset, feature_weight_trajectory_directions);

    assert(offset == nfeatures);

    database_build_bounds(db);
}

// Motion Matching search function essentially consists
// of comparing every feature vector in the database,
// against the query feature vector, first checking the
// query distance to the axis aligned bounding boxes used
// for the acceleration structure.
void motion_matching_search(
    int& __restrict best_index,
    float& __restrict best_cost,
    const slice1d<int> range_starts,
    const slice1d<int> range_stops,
    const slice2d<float> features,
    const slice1d<float> features_offset,
    const slice1d<float> features_scale,
    const slice2d<float> bound_sm_min,
    const slice2d<float> bound_sm_max,
    const slice2d<float> bound_lr_min,
    const slice2d<float> bound_lr_max,
    const slice1d<float> query_normalized,
    const bool brute_force,
    const float transition_cost,
    const int ignore_range_end,
    const int ignore_surrounding)
{

    int nfeatures = query_normalized.size;
    int nranges = range_starts.size;

    int curr_index = best_index;

    // Find cost for current frame
    if (best_index != -1)
    {
        best_cost = 0.0;
        for (int i = 0; i < nfeatures; i++)
        {
            best_cost += squaref(query_normalized(i) - features(best_index, i));
        }
    }

    float curr_cost = 0.0f;

    // Search rest of database
    for (int r = 0; r < nranges; r++)
    {
        // Exclude end of ranges from search
        int i = range_starts(r);
        int range_end = range_stops(r) - ignore_range_end;

        if (brute_force)
        {
            while (i < range_end)
            {
                // Skip surrounding frames
                if (curr_index != - 1 && abs(i - curr_index) < ignore_surrounding)
                {
                    i++;
                    continue;
                }

                // Check against each frame inside small box
                curr_cost = transition_cost;
                for (int j = 0; j < nfeatures; j++)
                {
                    curr_cost += squaref(query_normalized(j) - features(i, j));
                    if (curr_cost >= best_cost)
                    {
                        break;
                    }
                }

                // If cost is lower than current best then update best
                if (curr_cost < best_cost)
                {
                    best_index = i;
                    best_cost = curr_cost;
                }

                i++;
            }
            continue;
        }
        while (i < range_end)
        {
            // Find index of current and next large box
            int i_lr = i / BOUND_LR_SIZE;
            int i_lr_next = (i_lr + 1) * BOUND_LR_SIZE;

            // Find distance to box
            curr_cost = transition_cost;
            for (int j = 0; j < nfeatures; j++)
            {
                curr_cost += squaref(query_normalized(j) - clampf(query_normalized(j),
                    bound_lr_min(i_lr, j), bound_lr_max(i_lr, j)));

                if (curr_cost >= best_cost)
                {
                    break;
                }
            }

            // If distance is greater than current best jump to next box
            if (curr_cost >= best_cost)
            {
                i = i_lr_next;
                continue;
            }

            // Check against small box
            while (i < i_lr_next && i < range_end)
            {
                // Find index of current and next small box
                int i_sm = i / BOUND_SM_SIZE;
                int i_sm_next = (i_sm + 1) * BOUND_SM_SIZE;

                // Find distance to box
                curr_cost = transition_cost;
                for (int j = 0; j < nfeatures; j++)
                {
                    curr_cost += squaref(query_normalized(j) - clampf(query_normalized(j),
                        bound_sm_min(i_sm, j), bound_sm_max(i_sm, j)));

                    if (curr_cost >= best_cost)
                    {
                        break;
                    }
                }

                // If distance is greater than current best jump to next box
                if (curr_cost >= best_cost)
                {
                    i = i_sm_next;
                    continue;
                }

                // Search inside small box
                while (i < i_sm_next && i < range_end)
                {
                    // Skip surrounding frames
                    if (curr_index != - 1 && abs(i - curr_index) < ignore_surrounding)
                    {
                        i++;
                        continue;
                    }

                    // Check against each frame inside small box
                    curr_cost = transition_cost;
                    for (int j = 0; j < nfeatures; j++)
                    {
                        curr_cost += squaref(query_normalized(j) - features(i, j));
                        if (curr_cost >= best_cost)
                        {
                            break;
                        }
                    }

                    // If cost is lower than current best then update best
                    if (curr_cost < best_cost)
                    {
                        best_index = i;
                        best_cost = curr_cost;
                    }

                    i++;
                }
            }
        }
    }
}

// Search database
void database_search(
    int& best_index,
    float& best_cost,
    const database& db,
    const slice1d<float> query,
    const bool brute_force,
    const float transition_cost = 0.0f,
    const int ignore_range_end = 20,
    const int ignore_surrounding = 20)
{
    // Normalize Query
    array1d<float> query_normalized(db.nfeatures());
    for (int i = 0; i < db.nfeatures(); i++)
    {
        query_normalized(i) = (query(i) - db.features_offset(i)) / db.features_scale(i);
    }

    // Search
    motion_matching_search(
        best_index,
        best_cost,
        db.range_starts,
        db.range_stops,
        db.features,
        db.features_offset,
        db.features_scale,
        db.bound_sm_min,
        db.bound_sm_max,
        db.bound_lr_min,
        db.bound_lr_max,
        query_normalized,
        brute_force,
        transition_cost,
        ignore_range_end,
        ignore_surrounding);
}
