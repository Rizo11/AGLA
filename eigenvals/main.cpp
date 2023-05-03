#include <iostream>
#include <math.h>
#include <iomanip>
#include <bits/stdc++.h>
#include <cstdio>

using namespace std;

#define GNUPLOT_NAME "gnuplot -persist"

void prnt(double d) {
    cout << fixed << setprecision(2) << d << " ";
}

int main() {
    double v0, k0, a1, b1, a2, b2, T, N;

    cin >> v0 >> k0 >> a1 >> b1 >> a2 >> b2 >> T >> N;
    double step = T/N;
    double eigenvalue = sqrt(a1*a2);

    vector<double> time = vector<double>();
    vector<double> victimPopulation = vector<double>();
    vector<double> killerPopulation = vector<double>();


    double v_t, k_t, p1, p2;
    double rad;
    double sinn, coss;
    v0 -= a2/b2;
    k0 -= a1/b1;
    for (double t = 0; t <= T; t += step) {
        time.push_back(t);
        rad = eigenvalue * t;

        sinn = sin(rad);
        coss = cos(rad);

        p1 = v0 * coss;
        p2 = k0 * sqrt(a2) * b1 * sinn;
        p2 /= b2* sqrt(a1);

        v_t = p1 - p2 + a2/b2;

        p1 = k0 * coss;
        p2 = v0 * sqrt(a1) * b2 * sinn;
        p2 /= b1 * sqrt(a2);
        k_t = p2 + p1 + a1/b1;

        victimPopulation.push_back(v_t);
        killerPopulation.push_back(k_t);

    }

    FILE* pipe = popen(GNUPLOT_NAME, "w");

    if (pipe != NULL) {
        fprintf(pipe, "set xlabel 'Fox population'\n");
        fprintf(pipe, "set ylabel 'Rabbit population'\n");
        ::fprintf(pipe, "%s\n", "plot '-' using 1:2 title 'Rabbit' with lines");

        for (int i = 0; i < time.size(); ++i) {
            fprintf(pipe, "%f\t%f\n", killerPopulation[i], victimPopulation[i]);
        }

        fflush(pipe);


        pclose(pipe);
    } else {
        cout << "Could not open file" << endl;
    }

    cout << "t:" << endl;
    for (int i = 0; i < time.size(); ++i) {
        prnt(time[i]);
    }
    cout << endl;

    cout << "v:" << endl;
    for (int i = 0; i < time.size(); ++i) {
        prnt(victimPopulation[i]);
    }
    cout << endl;

    cout << "k:" << endl;
    for (int i = 0; i < time.size(); ++i) {
        prnt(killerPopulation[i]);
    }
    cout << endl;

    return 0;
}
