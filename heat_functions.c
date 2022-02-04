//
// Created by zyuhao on 2022-01-18.
//

// functions used in analytical heat transfer computation

#include "CA_heat.h"

void import_data(double* X, double* Y, double* L, double* D, double* T, int nrows) {
    // read_csv
    char filename[] = "DAQ_T50us_MS_Cube_fused_P04_layer0101.csv";
    FILE* the_file;
    int err;
    err = fopen(&the_file, filename, "r");
    if (err != 0) {
        perror("Unable to open the file.");
    }
    char line[nrows];
    int row = 0;
    while ((fgets(line, sizeof(line), the_file)) && (row < nrows)) {
        char* token;
        token = strtok(line, ",");
        X[row] = atof(token);
        token = strtok(NULL, ",");
        Y[row] = atof(token);
        token = strtok(NULL, ",");
        L[row] = atof(token);
        token = strtok(NULL, ",");
        D[row] = atof(token);
        token = strtok(NULL, ",");
        T[row] = atof(token);
        token = strtok(NULL, ",");

        printf("%lf  ", X[row]);
        printf("%lf  ", Y[row]);
        printf("%lf  ", L[row]);
        printf("%lf  ", D[row]);
        printf("%lf  ", T[row]);
        printf("\n");
        row++;
    }
    fclose(the_file);
}

void Update_Temperature_Field(double ***temp, double ***temp_old, int xsize, int ysize, int zsize, double *yc, double timeca, double Gy, double Vy){
    int i, j, k;
    double SF;

    SF = 0.0 + Vy*timeca;
    for (i = 0; i < xsize; i++) {
        for (j = 0; j < ysize; j++) {
            for (k = 0; k < zsize; k++) {

                temp_old[i][j][k] = temp[i][j][k];
                temp[i][j][k] = metalts + (yc[j]-SF)*Gy;

            }
        }
    }

}


