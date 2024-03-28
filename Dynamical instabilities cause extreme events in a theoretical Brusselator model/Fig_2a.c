/*
Dynamical instabilities cause extreme events in a theoretical Brusselator model.
Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)
*/

#include <stdio.h>
#include <math.h>

// Define a constant M with value 150
#define M 150

// Declare global variables
double omega, A, B, f, y[10];

// Function prototypes
void RK4(int, double, double, double[], void (*DGL)(double, double[], double[]));
void DGL(double, double[], double[]);

// Main function
int main() {
    int i, j, nn = 3;
    double t, h, pi;

    // File pointer for writing data to a file
    FILE *fp;

    // Open the file "bif_b.dat" for writing
    fp = fopen("bif_b.dat", "w");
    if (fp == NULL) {
        printf("File - Open Error\n");
        return 0;
    }

    // Value of pi
    pi = 4.0 * atan(1.0);

    // Initialize B
    B = 0.9;

    // Loop for different values of B
    do {
        // Initialize A, f, omega, and h
        A = 0.219;
        f = 0.06;
        omega = 0.7;
        h = (2.0 * pi) / (omega * 1000);

        // Initialize initial values for y
        y[1] = 0.01;
        y[2] = 0.02;

        // Iterate for M steps
        for (j = 0; j < M; j++) {
            for (i = 0; i < 1000; i++) {
                t = h * (double)(i);
                RK4(nn, h, t, y, DGL);
            }
            // Write data to the file for the last 60 steps
            if (j > (M - 60))
                fprintf(fp, "%20.12lf%20.12lf%20.12lf\n", B, y[1], y[2]);
            printf("%20.12lf\n", B);
        }

        // Increment B
        B += 0.00001;
    } while (B < 1.2);

    // Close the file
    fclose(fp);

    return 0;
}

// Function to perform Runge-Kutta 4th order integration
void RK4(int nn, double h, double t, double y[10], void (*DGL)(double, double[], double[])) {
    int i;
    double k1[10], k2[10], k3[10], k4[10];
    double yaux[10];

    // Compute k1
    DGL(t, y, k1);

    // Compute k2
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k1[i] / 2.0;
    }
    DGL(t + h / 2.0, yaux, k2);

    // Compute k3
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k2[i] / 2.0;
    }
    DGL(t + h / 2.0, yaux, k3);

    // Compute k4
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k3[i];
    }
    DGL(t + h, yaux, k4);

    // Update y using the RK4 formula
    for (i = 1; i <= nn; i++) {
        y[i] = y[i] + h * (k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0);
    }
}

// Function defining the system of differential equations
void DGL(double t, double y[10], double F[10]) {
    F[1] = A - B * y[1] - y[1] + y[1] * y[1] * y[2] + f * sin(y[3]);
    F[2] = B * y[1] - y[1] * y[1] * y[2];
    F[3] = omega;
}






