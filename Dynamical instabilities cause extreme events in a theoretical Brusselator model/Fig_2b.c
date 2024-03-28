/*
Dynamical instabilities cause extreme events in a theoretical Brusselator model.
Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)
*/

#include <stdio.h>
#include <math.h>

// Define global variables
int n = 3; // Dimension
double A, B, t, omega, f, y[25]; // Parameters and state variables
double znorm[10], gsc[10], cum[10]; // Arrays for calculations

// Function prototypes
void RK4(int, double, double, double[], void (*DGL)(double, double[], double[]));
void DGL(double, double[], double[]);
void RK4_1(int, double, double, double[], void (*DGL_1)(double, double[], double[]));
void DGL_1(double, double[], double[]);

// Main function
int main() {
    int i, j, k, l, nn = 12, tf = 1000; // Loop variables and parameters
    double h, pi, lexp, t1; // Step size, pi, and auxiliary variables

    FILE *fp; // File pointer for writing data to a file

    // Open the file "bo_lyp_revised2.dat" for writing
    fp = fopen("bo_lyp_revised2.dat", "w");
    if (fp == NULL) {
        printf("File - Open Error\n");
        return 0;
    }

    // Initial value of B
    B = 1.095;

    // Loop for different values of B
    while (B < 1.099) {
        // Value of pi
        pi = 4.0 * atan(1.0);
        t = 0.0;

        // Initialize parameters
        A = 0.2;
        f = 0.06;
        omega = 0.7;

        // Calculate step size
        h = (2.0 * pi) / (omega * 1098);

        // Initial conditions for the nonlinear system
        y[1] = 0.01;
        y[2] = 0.02;

        // Initial conditions for a linear system (orthonormal frame)
        for (i = n + 1; i <= nn; i++)
            y[i] = 0.0;
        for (i = 1; i <= n; i++) {
            y[(n + 1) * i] = 1.0;
            cum[i] = 0.0;
        }

        // Integration loop
        for (i = 0; i <= tf; i++) {
            for (j = 1; j <= 800; j++) {
                t1 = t + h * (double)(j);
                RK4_1(nn, h, t1, y, DGL_1);
            }

            // Construction of a new orthonormal basis by Gram-Schmidt process
            // Normalization of the first vector
            znorm[1] = 0.0;
            for (j = 1; j <= n; j++)
                znorm[1] += y[n * j + 1] * y[n * j + 1];
            znorm[1] = sqrt(znorm[1]);
            for (j = 1; j <= n; j++)
                y[n * j + 1] /= znorm[1];

            // Generation of a new set of orthonormal vectors
            for (j = 2; j <= n; j++) {
                for (k = 1; k <= (j - 1); k++) {
                    gsc[k] = 0.0;
                    for (l = 1; l <= n; l++)
                        gsc[k] += y[n * l + j] * y[n * l + k];
                }
                for (k = 1; k <= n; k++) {
                    for (l = 1; l <= (j - 1); l++)
                        y[n * k + j] -= gsc[l] * y[n * k + l];
                }
                znorm[j] = 0.0;
                for (k = 1; k <= n; k++)
                    znorm[j] += y[n * k + j] * y[n * k + j];
                znorm[j] = sqrt(znorm[j]);
                for (k = 1; k <= n; k++)
                    y[n * k + j] /= znorm[j];
            }

            // Compute cumulative average
            for (k = 1; k <= n; k++)
                cum[k] += log(znorm[k]) / log(2.0);
            t = t1;
        }

        // Compute final cumulative average
        for (k = 1; k <= n; k++)
            cum[k] /= t;

        // Write data to the file
        fprintf(fp, "%20.12lf%20.12lf%20.12lf%20.12lf\n", B, cum[1], cum[2], cum[3]);
        printf("%20.12lf\n", B);

        // Increment B
        B += 0.00001;
    }

    // Close the file
    fclose(fp);

    return 0;
}

// Function to perform Runge-Kutta 4th order integration
void RK4(int nn, double h, double t, double y[25], void (*DGL)(double, double[], double[])) {
    int i;
    double k1[25], k2[25], k3[25], k4[25];
    double yaux[25];

    DGL(t, y, k1);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k1[i] / 2.0;
    }
    DGL(t + h / 2.0, yaux, k2);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k2[i] / 2.0;
    }
    DGL(t + h / 2.0, yaux, k3);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k3[i];
    }
    DGL(t + h, yaux, k4);
    for (i = 1; i <= nn; i++) {
        y[i] = y[i] + h * (k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0);
    }
}

// Function defining the system of differential equations
void DGL(double t, double y[25], double F[25]) {
    F[1] = A - B * y[1] - y[1] + y[1] * y[1] * y[2] + f * sin(y[3]);
    F[2] = B * y[1] - y[1] * y[1] * y[2];
    F[3] = omega;
}

// Function to perform Runge-Kutta 4th order integration for the linear system
void RK4_1(int nn, double h, double t, double y[25], void (*DGL_1)(double, double[], double[])) {
    int i;
    double k1[25], k2[25], k3[25], k4[25];
    double yaux[25];

    DGL_1(t, y, k1);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k1[i] / 2.0;
    }
    DGL_1(t + h / 2.0, yaux, k2);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k2[i] / 2.0;
    }
    DGL_1(t + h / 2.0, yaux, k3);
    for (i = 1; i <= nn; i++) {
        yaux[i] = y[i] + h * k3[i];
    }
    DGL_1(t + h, yaux, k4);
    for (i = 1; i <= nn; i++) {
        y[i] = y[i] + h * (k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0);
    }
}

// Function defining the system of differential equations for the linear system
void DGL_1(double t, double y[25], double F[25]) {
    int i;
    double yy1;

    F[1] = A - B * y[1] - y[1] + y[1] * y[1] * y[2] + f * sin(y[3]);
    F[2] = B * y[1] - y[1] * y[1] * y[2];
    F[3] = omega;

    for (i = 0; i < n; i++) {
        F[4 + i] = -(B * y[4 + i]) - (y[4 + i]) + (2.0 * y[1] * y[2] * y[4 + i]) + (y[1] * y[1] * y[7 + i]) + (f * cos(y[3]) * y[10 + i]);
        F[7 + i] = (B * y[4 + i]) - (2.0 * y[1] * y[2] * y[4 + i]) - (y[1] * y[1] * y[7 + i]);
        F[10 + i] = 0.0;
    }
}
