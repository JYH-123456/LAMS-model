#include "snow_model.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>

int initialize_atmos_array(atmos_array *atmos, int num) {
    atmos->size = num;
    atmos->air_temp = (double *)calloc(num, sizeof(double));
    atmos->snow = (double *)calloc(num, sizeof(double));
    atmos->rain = (double *)calloc(num, sizeof(double));
    atmos->in_short = (double *)calloc(num, sizeof(double));
    atmos->in_long = (double *)calloc(num, sizeof(double));
    atmos->wind = (double *)calloc(num, sizeof(double));
    atmos->air_pressure = (double *)calloc(num, sizeof(double));
    atmos->vapor_pressure = (double *)calloc(num, sizeof(double));
    atmos->humidity = (double *)calloc(num, sizeof(double));
    atmos->Tg1 = (double *)calloc(num, sizeof(double));
    atmos->Tg2 = (double *)calloc(num, sizeof(double));
    atmos->Tg3 = (double *)calloc(num, sizeof(double));

    if (!atmos->air_temp || !atmos->rain || !atmos->snow || !atmos->in_short || 
        !atmos->in_long || !atmos->wind || !atmos->air_pressure || !atmos->vapor_pressure || 
        !atmos->humidity || !atmos->Tg1 || !atmos->Tg1 || !atmos->Tg3) {
        printf("Memory allocation failed.\n");
        free_atmos_array(atmos);
        return 1;
    }
    return 0;
}

void free_atmos_array(atmos_array *atmos) {
    free(atmos->air_temp);
    free(atmos->rain);
    free(atmos->snow);
    free(atmos->in_short);
    free(atmos->in_long);
    free(atmos->wind);
    free(atmos->air_pressure);
    free(atmos->vapor_pressure);
    free(atmos->humidity);
    free(atmos->Tg1);
    free(atmos->Tg2);
    free(atmos->Tg3);
}


int read_atmos_data(atmos_array *atmos, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        return -1;
    }

    for (int i = 0; i<atmos->size; i++) {
        int result = fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                            &atmos->rain[i], &atmos->snow[i], &atmos->air_temp[i], &atmos->in_short[i],
                            &atmos->in_long[i], &atmos->air_pressure[i], &atmos->vapor_pressure[i],
                            &atmos->wind[i], &atmos->humidity[i], &atmos->Tg1[i], atmos->Tg2, atmos->Tg3);
        if (result != 12) {
            fprintf(stderr, "Error reading data at line %d: expected 12 values, got %d\n", i + 1, result);
            fclose(file);
            return -1;
        }        
    }
    fclose(file);
    return 0;
}

int read_param_file(const char *filename, ParamData *params) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "无法打开文件: %s\n", filename);
        return 0;
    }
    int count = fscanf(file, "%lf\t%lf\t%lf\t%lf", &params->Tpmax, &params->Tpmin, &params->C_h, &params->Tf);
    fclose(file);
    if (count != param_num) {
        fprintf(stderr, "文件格式错误: %s\n", filename);
        return 0;
    }
    return 1;
}

// int allocate_soils(soillayers *soils, int num){
//     soils->num_grids = num;
//     soils->soil_temp = (double *)calloc(num, sizeof(double));
//     soils->flow_mass = (double *)calloc(num, sizeof(double));
//     soils->flow_heat = (double *)calloc(num, sizeof(double));
//     soils->thermal_cond =(double *)calloc(num, sizeof(double));
//     soils->heat_cap = (double *)calloc(num, sizeof(double));
//     soils->water_content = (double *)calloc(num, sizeof(double));
//     soils->porosity = (double *)calloc(num, sizeof(double));
//     soils->hydro_cond = (double *)calloc(num, sizeof(double));
//     soils->water_pot = (double *)calloc(num, sizeof(double));
//     soils->pram_b = (double *)calloc(num, sizeof(double));
//     soils->D = (double *)calloc(num, sizeof(double));
//     soils->water_cont_dt = (double *)calloc(num, sizeof(double));
//     soils->satur_k = (double *)calloc(num, sizeof(double));
//     soils->satur_pot = (double *)calloc(num, sizeof(double));
//     soils->ice_cont = (double *)calloc(num, sizeof(double));
//     soils->volumn_fra =(double *)calloc(num, sizeof(double));
//     soils->frozen = (double *)calloc(num, sizeof(double));
//     soils->melt =  (double *)calloc(num, sizeof(double));
//     return 0;
// }

// void initialize_soils(soillayers *soils, char *filename){
//     FILE *file = fopen(filename, "r");
//     if (!file) {
//         perror("Failed to open file");
//         return -1;
//     }
//     int num = soils->num_grids;
//     for (int i = 0; i<num; i++) {
//         int result = fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
//                             &soils->soil_temp[i], &soils->thermal_cond[i], &soils->heat_cap[i],
//                             &soils->water_content[i], &soils->porosity[i], &soils->hydro_cond[i],
//                             &soils->water_pot[i], &soils->pram_b[i], &soils->satur_k);
//         if (result != 9) {
//             fprintf(stderr, "Error reading data at line %d: expected 9 values, got %d\n", i + 1, result);
//             fclose(file);
//             return -1;
//         }
//     }
//     soils->soil_depth = ground_depth/(num+1);
//     soils->density = soil_density;
//     for (int i=0; i<num; i++){
//         soils->hydro_cond[i] = soils->satur_k[i]*pow(soils->water_content[i]/soils->porosity[i], 2*soils->pram_b[i]+3);
//         soils->water_pot[i] = soils->satur_pot[i]*pow(soils->water_content[i]/soils->porosity[i], -soils->pram_b[i]);
//         soils-> D[i] = -soils->hydro_cond[i]*soils->porosity[i]/soils->pram_b[i]/soils->satur_pot[i]*
//                 pow(soils->water_pot[i]/soils->satur_pot[i], -1/soils->pram_b[i]-1);
//     }
//     soils->water_flow_bottom =soils->hydro_cond[num-1]*(1+(soils->satur_pot[num-1]-soils->water_pot[num-1])/(0.5*soils->soil_depth));
// }

// void free_soils(soillayers *soils) {
//     free(soils->soil_temp);
//     free(soils->flow_mass);
//     free(soils->flow_heat);
//     free(soils->thermal_cond);
//     free(soils->heat_cap);
//     free(soils->water_content);
//     free(soils->porosity);
//     free(soils->hydro_cond);
//     free(soils->water_pot);
//     free(soils->pram_b);
//     free(soils->hydro_cond);
//     free(soils->D);
//     free(soils->water_cont_dt);
//     free(soils->satur_k);
//     free(soils->satur_pot);
//     free(soils->ice_cont);
//     free(soils->volumn_fra);
//     free(soils->frozen);
//     free(soils->melt);
// }
