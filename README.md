# Parallel-Game-Engine
modified this single-threaded game engine to make it run as fast  as possible.
Instructions
 *
 * Your job is to modify this single-threaded game engine to make it run as fast
 * as possible. You can use a combination of threading and code optimization to
 * improve the runtime. Speed is measured in terms of the number of frames
 * processed per second. I have inserted some timing code into the main game
 * loop (calls to clock_gettime()). You might have to move the timing code
 * around, but you cannot cheat (for example by starting the timer after your
 * computations start running).
 *
 * The 5 fastest submissions will get extra credit.
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define NPOLYS 1000000
#define NROTATIONS 10

// Create timespec structs, which hold the start and end times of our computation
struct timespec tstart={0,0}, tend={0,0};

struct point {
  double x, y, z;
};

struct polygon {
  struct point pt[3];
};


/*
 * matrix_vec_mult
 *
 * Multiplies square matrix mat (of size dim_size x dim_size) by vector vec.
 * Output written to new_vec. mat is stored as an array of doubles in row-major
 * order (ie rows are contiguous in memory).
 */
void matrix_vec_mult(double* mat, uint32_t dim_size, double* vec, double* new_vec) {
    for (int i = 0; i < dim_size; i++) {
        new_vec[i] = 0.;
        for (int j = 0; j < dim_size; j++) {
            new_vec[i] += mat[i*dim_size + j] * vec[j];
        }
    }
}

/*
 * build_Rx
 *
 * Computes the elements of 3x3 rotation matrix Rx, which rotates a point by
 * angle theta_x around the x-axis. Rx is a 9-element array of doubles
 * allocated by the caller.
 *
 * See the link below for a more detailed description:
 *
 * https://en.wikipedia.org/wiki/Rotation_matrix
 */
void build_Rx(double *Rx, double theta_x) {
    Rx[0] = 1.;
    Rx[1] = 0.;
    Rx[2] = 0.;
    Rx[3] = 0.;
    Rx[4] = cos(theta_x);
    Rx[5] = -sin(theta_x);
    Rx[6] = 0.;
    Rx[7] = sin(theta_x);
    Rx[8] = cos(theta_x);
}

void build_Ry(double *Ry, double theta_y) {
    Ry[0] = cos(theta_y);
    Ry[1] = 0.;
    Ry[2] = sin(theta_y);
    Ry[3] = 0.;
    Ry[4] = 1.;
    Ry[5] = 0.;
    Ry[6] = -sin(theta_y);
    Ry[7] = 0.;
    Ry[8] = cos(theta_y);
}

void build_Rz(double *Rz, double theta_z) {
    Rz[0] = cos(theta_z);
    Rz[1] = -sin(theta_z);
    Rz[2] = 0.;
    Rz[3] = sin(theta_z);
    Rz[4] = cos(theta_z);
    Rz[5] = 0.;
    Rz[6] = 0.;
    Rz[7] = 0.;
    Rz[8] = 1.;
}

int main() {
    // Create an array of polygons to work with
    struct polygon *polys = malloc(NPOLYS *sizeof(struct polygon));
    double Rx[9], Ry[9], Rz[9]; // Rotation matrices
    double theta_x = 0., theta_y = 0., theta_z = 0.; // Angles to rotate each point by

    // Get time from the operating system. This is reported as the number of nanoseconds
    // since midnight Jan 1, 1970
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    // MAIN GAME LOOP BELOW
    for(int frame = 0; frame < NROTATIONS; frame++) {
        printf("frame %d\n", frame);
        // 1. Build rotation matrices
        build_Rx(Rx, theta_x);
        build_Ry(Ry, theta_y);
        build_Rz(Rz, theta_z);

        // 2. Rotate every point in every polygon
        for(int p = 0; p < NPOLYS; p++) {
            for(int q = 0; q < 3; q++) { // 3 points per polygon
                struct point rot_x, rot_y, rot_z;

                // Rotate current point around the x-axis, stashing the result in rot_x
                matrix_vec_mult(Rx, 3, (double*) (&polys[p].pt[q]), (double*)&rot_x);

                // Rotate current point around y-axis, stashing result in rot_y
                matrix_vec_mult(Ry, 3, (double*) (&rot_x), (double*)&rot_y);

                // Rotate current point around z-axis, stashing result in rot_z.
                // This is the rotated point that will get projected onto the screen.
                matrix_vec_mult(Rz, 3, (double*) (&rot_y), (double*)&rot_z);
            }
        }

        // Update the rotation angle for each point in the scene. This is
        // usually done from user input, but we'll just hardcode the updates.
        theta_x += 0.1;
        theta_y += 0.1;
        theta_z += 0.1;
    }
    // Get time from the operating system. This is reported as the number of
    // nanoseconds since midnight Jan 1, 1970. tstart and tend are global
    // variables declared at the top of the file.
    clock_gettime(CLOCK_MONOTONIC, &tend);

    // compute time difference between tstart and tend
    double elapsed_time = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - 
                          ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
    printf("drew %f frames per second\n", NROTATIONS/elapsed_time);
}
