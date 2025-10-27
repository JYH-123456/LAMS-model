#include <stdio.h>
#include <stdlib.h>
#include "snow_model.h"
FILE *outputFile;

void init_model(SnowModel *model) {
    model->head = NULL;
    model->tail = NULL;
    model->num_layers = 0;
    model->total_thick = 0;
    model->ave_thick = 0;
    model->snow_cover = 0;
    model->Tgrnd = -1;
    model->Tsfc = -1;
    model->T1 = -1;
    model->outflow = 0;
    model->fi = 0;
    model->fl = 0;
    model->density = 0;
    model-> temp = -1;
    model->refreeze = 0;
    model->melt = 0;
    model->swe = 0;
    model->thick_add = 0;
}

void add_node(SnowModel *model, 
            atmos_array *atmos,
            ParamData *params,
            double snowfall, 
            double rainfall, 
            double h, 
            int current_step) {
    double new_density;
    double ice_ratio;
    double airtemp = atmos->air_temp[current_step]+273.15;
    double wind = atmos->wind[current_step];
    SnowLayerNode *new_node = (SnowLayerNode *)calloc(1, sizeof(SnowLayerNode));
    if (!new_node) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    new_density = new_snow_density(airtemp, wind, params);
    new_node->thick = snowfall*density_water/new_density;
    model->thick_add = new_node->thick;
    new_node->ice_ratio = snowfall*density_water/density_ice/new_node->thick; 
    new_node->ice_content = snowfall*density_water; 
    new_node->water_content = rainfall*density_water;
    new_node->water_ratio = rainfall*density_water/density_water/new_node->thick;
    new_node->density = (new_node->ice_content+new_node->water_content)/new_node->thick;
    new_node->refreeze = 0;    
    new_node->rain_in = 0; // rainfall/h
    new_node->density_dt = 0;
    new_node->water_out = 0;
    new_node->delta_temp = 0;    
    new_node->ice_in = 0;
    new_node->temp = atmos->air_temp[current_step]+273.15;
    new_node->snow_time = current_step;
    new_node->Cs = (new_node->water_ratio*density_water*C_l+new_node->ice_ratio*density_ice*C_i)/new_node->density;
    ice_ratio = new_node->ice_ratio;
    new_node->f_l = A1*ice_ratio*(1-ice_ratio);
    new_node -> short_in = atmos->in_short[current_step];
    new_node -> long_in = atmos->in_long[current_step];
    if (model->num_layers == 0) {
        new_node->layer_id = 1;
        new_node->upper_layer = NULL;
        new_node->lower_layer = NULL;
        model->head = new_node;
        model->tail = new_node;
    } else {
        new_node->layer_id = model->head->layer_id + 1;
        new_node->lower_layer = model->head;
        new_node->upper_layer = NULL;
        model->head->upper_layer = new_node;
        model->head = new_node;
    }
    model->num_layers++;
    model->total_thick += new_node->thick;
    // printf("add a new node with depth: %.4fcm\n", new_node->thick*100);
}

void remove_node(SnowModel *model, SnowLayerNode *node) {
    if (node == NULL) return;

    if (node->upper_layer != NULL) {
        node->upper_layer->lower_layer = node->lower_layer;
    } else {
        model->head = node->lower_layer;
    }

    if (node->lower_layer != NULL) {
        // node->lower_layer->ice_content += node->ice_content;
        // node->lower_layer->water_content += node->water_content;
        node->lower_layer->upper_layer = node->upper_layer;
    } else {
        model->tail = node->upper_layer;
        if (node->upper_layer!=NULL){
            model->outflow += node->upper_layer->water_out;
        }
    }
    model->outflow += node->density*node->thick;
    model->total_thick -= node->thick;
    free(node);
    model->num_layers--;
}

void print_layers_to_file(const SnowModel *model, FILE *file) {
    SnowLayerNode *current = model->head;
    fprintf(file, "ID\tThickness(cm)\tDensity(kg/m^3)\tTemp(K)\tRefreeze(mm)\tSnow Time\twater ratio\tice ratio\tCs\n");
    while (current != NULL) {
        fprintf(file, "%d\t%.2f\t%.2f\t%.2f\t%.4f\t%d\t%.3f\t%.3f\t%.2f\n",
               current->layer_id,
               current->thick*100,
               current->density,
               current->temp,
               current->refreeze*3600,
               current->snow_time,
            //    current->upper_layer ? current->upper_layer->layer_id : -1,
            //    current->lower_layer ? current->lower_layer->layer_id : -1,
               current->water_ratio,
               current->ice_ratio,
               current->Cs);
        current = current->lower_layer;
    }
    fprintf(file, "Ground Temperature: %.2fC, Total Thickness: %.2fcm, Number of Layers: %d\n", 
            model->Tgrnd-273.15, model->total_thick*100, model->num_layers);
    fprintf(file, "------------------------------------------------------------\n");
}

void initialize_output_file(const char *filename) {
    // 打开文件以追加模式写入
    outputFile = fopen(filename, "w");
    if (outputFile == NULL) {
        fprintf(stderr, "无法打开文件 %s\n", filename);
        exit(EXIT_FAILURE);
    }
}

void close_output_file() {
    // 关闭文件
    if (outputFile != NULL) {
        fclose(outputFile);
        outputFile = NULL;
    }
}



