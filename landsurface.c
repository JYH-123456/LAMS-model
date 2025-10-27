#include "snow_model.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double rain_snow_seperate(ParamData *params, double prec_amount, double air_temp, double humidity, double air_pressure){
    double es;
    double Tw;
    double delta_es;
    double psnow;
    double psleet;
    double prain;
    if (air_temp> params->Tpmax){
        return 0;
    }else if(air_temp< params->Tpmin){
        return 1;
    }
    else{
        double ratio = (params->Tpmax-air_temp)/(params->Tpmax-params->Tpmin);
        return ratio;
    }
    // es = 6.1078*exp(17.27*air_temp/(air_temp+237.3));
    // delta_es = es*(17.27*237.3)/pow(air_temp+237.3, 2);
    // Tw = air_temp-es*(1-humidity)/(0.000643*air_pressure*10+delta_es);
    // psnow = 1/(1+exp((Tw+4.035)/0.9127));
    // psleet = 1/(1+exp((Tw+1.633)/0.9127)) - psnow;
    // prain = 1-psnow-psleet;
    // if (psnow>= psleet && psnow>=prain){
    //     return 1;
    // }
    // if(prain>= psleet && prain>=psnow){
    //     return 0;
    // }
    // if(psleet>= prain && psleet>=psnow){
    //     return psnow/(psnow+prain);
    // }
}


void update_prec(SnowModel *model, atmos_array *atmos, ParamData *params, int t_count, double h){
    double snowfall;
    double rainfall;
    double air_temp;
    double humidity;
    double wind;
    double thickness;
    double new_density;
    SnowLayerNode *top_node;
    SnowLayerNode *new_node;
    air_temp = atmos->air_temp[t_count]+273.15;
    model->Tgrnd = atmos->Tg1[t_count]+273.15;  // convert C to K   
    humidity = atmos->humidity[t_count];
    wind = atmos->wind[t_count];
    snowfall = atmos->snow[t_count];
    rainfall = atmos->rain[t_count];
    snowfall = snowfall/1000;  // convert mm to m
    rainfall = rainfall/1000;  // convert mm to m
    new_density = new_snow_density(air_temp, wind, params);
    thickness = snowfall*density_water/new_density;
    model->outflow = 0;
    if (model->num_layers>0){  //snowfall+rainfall>0 &&
        top_node = model->head;
        if((top_node->thick+thickness<SNOW_THRESHOLD && t_count-top_node->snow_time<=TIME_THRESHOLD) || thickness<thick_min){
            top_node->rain_in = rainfall/h;
            top_node->ice_in = snowfall/h;
            // double water_content = top_node->water_ratio*density_water*top_node->thick+rainfall*density_water;
            // double ice_content = top_node->ice_ratio*density_ice*top_node->thick+snowfall*density_water;
            // top_node->ice_content = ice_content;
            // top_node->water_content = water_content;
            top_node->thick += thickness;
            // top_node->ice_ratio = ice_content/top_node->thick/density_ice;
            // top_node->water_ratio = water_content/density_water/(top_node->thick);
            // top_node->density = (ice_content+water_content)/(top_node->thick);
            model->total_thick += thickness;
            model->thick_add = thickness;
        }else{
            if (snowfall>0) {
                add_node(model, atmos, params, snowfall, rainfall, h, t_count);
            }else{
                top_node->rain_in = rainfall/h;
                model->thick_add = 0;
            }   
        }
    }else{
        if (snowfall>0 && thickness>thick_min){
            add_node(model, atmos, params, snowfall, rainfall, h, t_count);
        }else{
            model->outflow = (snowfall+rainfall)*density_water;
            model->thick_add = 0;
        }   
    }
}
