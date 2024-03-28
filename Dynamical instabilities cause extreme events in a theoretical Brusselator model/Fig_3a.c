/*
Dynamical instabilities cause extreme events in a theoretical Brusselator model.
Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)
*/

#include<stdio.h>
#include<math.h>

// Define constants
#define M 150

// Declare global variables
double omega, A, B, f, y[10];

// Function prototypes
void RK4(int, double, double, double[], void(*DGL)(double, double[], double[]));
void DGL(double, double[], double[]);

// Main function
int main() {
    int i, nn = 3;
    double t, h, pi, total_time, step, skip_time, skip;

    FILE *fp;

    // Open output file
    fp = fopen("bo_poin.dat", "w");
    if (fp == NULL) {
        printf("File - Open Error\n");
        return 0;
    }

    // Get input from the user
    printf("Can you enter the B value???\n");
    scanf("%lf", &B);

    // Set initial values
    A = 0.2;
    f = 0.06;
    omega = 0.7;
    pi = 4.0 * atan(1.0);
    h = (2.0 * pi) / (omega * 1000);
    total_time = 2000;
    step = total_time / h;
    skip_time = 5000;
    skip = skip_time / h;
    y[1] = 0.01;
    y[2] = 0.02;

    // Perform integration using RK4 method
    for (i = 0; i < step; i++) {
        t = h * (double)(i);
        RK4(nn, h, t, y, DGL);
        if (i > skip) {
            fprintf(fp, "%20.12lf%20.12lf%20.12lf\n", t, y[1], y[2]);
        }
    }

    // Close output file
    fclose(fp);
    return 0;
}

// RK4 function definition
void RK4(int nn, double h, double t, double y[10], void(*DGL)(double, double[], double[])) {
    int i;
    double k1[10], k2[10], k3[10], k4[10];
    double yaux[10];

    // Calculate k1
    DGL(t, y, k1);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k1[i] / 2.0;
    }

    // Calculate k2
    DGL(t + h / 2.0, yaux, k2);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k2[i] / 2.0;
    }

    // Calculate k3
    DGL(t + h / 2.0, yaux, k3);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k3[i];
    }

    // Calculate k4 and update y using RK4 formula
    DGL(t + h, yaux, k4);
    for (i = 1; i <= nn; i++) {
        y[i] = y[i] + h * (k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0);
    }
}

// DGL function definition
void DGL(double t, double y[10], double F[10]) {
    F[1] = A - B * y[1] - y[1] + y[1] * y[1] * y[2] + f * sin(y[3]);
    F[2] = B * y[1] - y[1] * y[1] * y[2];
    F[3] = omega;
}
