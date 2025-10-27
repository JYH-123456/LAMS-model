#include "snow_model.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
static FILE *snow_file = NULL;
static char snow_filename[100] = "";
static char output_filename[100] = "";


int main(){
        
    // initialize date time list
    atmos_array atmos;
    DateList datelist;
    ParamData params;
    int total_steps;
    double delta_time; // Time step in seconds

    // initialize date list
    char station[10] = "cdp";
    int start_year = 2006;
    int start_month = 10;
    int start_day = 1;
    int start_hour = 0;
    int end_year = 2007;
    int end_month = 10;
    int end_day = 1;
    int end_hour = 0;
    delta_time = delta_t*3600;
    datelist = generate_datetime_list(start_year,start_month, start_day, start_hour,
                                        end_year, end_month, end_day, end_hour,
                                        delta_t);
    total_steps = datelist.count;
    // initialize atmos array
    int array_size = datelist.count;
    char atmos_file[50];
    snprintf(atmos_file, sizeof(atmos_file), "data/%s_%d_%02d-%d_%02d_input.txt", 
            station, start_year, start_month, end_year, end_month);
    
    if (initialize_atmos_array(&atmos, array_size) != 0) {
        printf("Failed to initialize atmospheric data.\n");
        return 1;
    }
    if (read_atmos_data(&atmos, atmos_file)==0){
        printf("Atmospheric data read successfully.\n");
    }

    // set corresponding parameters
    char param_path[256];
    params.Tpmin = T_min_rain;
    params.Tpmax = T_max_snow;
    params.C_h = 0.003;
    params.albedo_max = 0.90;
    params.albedo_min = 0.35;
    params.Tf = 273.15;
    
    // initialize soil layers
    // int num_grids = 10;
    // char soil_file[50] = "soils.txt";
    // soillayers soils;
    // allocate_soils(&soils, num_grids);
    // if (initialize_soils(&soils, soil_file)!= 0){
    //     printf("Failed to initialize soils data.\n");
    //     return 1;
    // }

    // initialize snow model
    SnowModel model;
    init_model(&model);
    snprintf(snow_filename, sizeof(snow_filename), "output3/snow_%s_%d-%d.txt", 
            station, start_year, end_year);
    snprintf(output_filename, sizeof(output_filename), "output3/output_%s_%d-%d.txt", 
            station, start_year, end_year);
    initialize_output_file(output_filename);

    // start simulation loop
    int time_step = 0;
    while (time_step < total_steps) { // Simulate for 24 hours
        // printf("time step:%d\n", time_step);
        update_prec(&model, &atmos, &params, time_step, delta_time);
        update_snow_layer(&model, &atmos, &params, datelist, time_step, delta_time);
        // 记录单独的降雪积雪事件
        snow_file = fopen(snow_filename, "w");
        if (snow_file == NULL) {
            printf("无法打开文件 %s 进行写入\n", snow_filename);
        } else {
            // 写入文件头
            snow_file = fopen(snow_filename, "a");
            fprintf(snow_file, "Time Step: %d\n", time_step);
            print_layers_to_file(&model, snow_file);
        }   
        // 输出output文件
        fprintf(outputFile, "%d\t%.3f\t%.3f\t%.3f\t%.7f\t%.7f\t%.7f\t%.2f\t%.2f\t%.7f\n",
            model.num_layers, model.total_thick*100, model.swe*100, model.Tsfc,
            model.melt, model.refreeze, model.outflow, model.density, model.Tgrnd-273.15, model.temp);
        //}
        time_step += 1;
    }
    close_output_file();
    fclose(snow_file);
    return 0;
}
