#include "snow_model.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// double snow_albedo(SnowLayerNode *top_node, 
//                 atmos_array *atmos,
//                 ParamData *params,
//                 double h,
//                 int t_count){
//     double albedo;
//     double value;
//     double Tc = params->Tf-2;
//     value = (top_node->temp-Tc)/(params->Tf-Tc);
//     if (value>0){
//         albedo = ALBEDO_MAX+(ALBEDO_MIN-ALBEDO_MAX)*value;
//     }else{
//         albedo =  ALBEDO_MAX;
//     }
//     return albedo;
// }

double snow_albedo(SnowModel *model,
                SnowLayerNode *top_node,
                atmos_array *atmos,
                ParamData *params,
                double h,
                int t_count){
    double refreeze;
    double old_albedo=model->albedo;
    double new_albedo;
    double temp = top_node->temp;
    double albedo_dt=0;
    double Sf;
    refreeze = top_node->refreeze;

    if(top_node->ice_in>0){
        Sf = top_node->ice_in*h*density_water;
        // printf("Sf:%.2f\n", Sf);
        albedo_dt = (0.95-old_albedo)*Sf/10;
        old_albedo = old_albedo + albedo_dt;
    }
    if (model->head->refreeze<0){
        new_albedo = 0.45+(old_albedo-0.45)*exp(-1/360000*h);
    }else{
        new_albedo = old_albedo-0.0000002*h;
    }
    if (new_albedo < 0.2){
        new_albedo = 0.2;
    }

    return new_albedo;
}


// double new_snow_density(double air_temp){
//     double density_new;
//     double airtemp;
//     airtemp = C_TO_F(air_temp);
//     if (airtemp > 0) {
//         density_new = SNOW_NEW_SNOW_DENSITY + 1000 *
//                         (airtemp / SNOW_NEW_BRAS_DENOM) *
//                         (airtemp / SNOW_NEW_BRAS_DENOM);
//     }
//     else {
//         density_new = SNOW_NEW_SNOW_DENSITY;
//     }
// }

double new_snow_density(double air_temp, double wind, ParamData *params){
    double new_density;
    double Tm = params->Tf;
    new_density = 109 + 6*(air_temp-Tm)+ 26*sqrt(wind);
    if (new_density < Density_min){
        new_density = Density_min;
    }
    return new_density;
}

void snow_compaction(SnowModel *model,
                    ParamData *params,
                    double h){
    double vis;
    double d_density;
    double stress_f;
    double value;
    double new_density;
    double total_thick=0;
    double air_temp;
    double ice_ratio;
    SnowLayerNode *top_node;
    SnowLayerNode *next_node;
    SnowLayerNode *node;
    if (model->num_layers>0){
        top_node = model->head;
        top_node->gravity_stress = top_node->density*top_node->thick*gravity;
        vis = 37000000*exp(0.081*(params->Tf-top_node->temp)+0.018*top_node->density);
        value = top_node->density-100;
        if (value<0){
            value = 0;
        }
        d_density = top_node->gravity_stress/vis+0.0000028*exp(-0.042*(params->Tf-top_node->temp)-0.046*value);
        d_density = d_density * top_node->density;
        new_density = top_node->density + d_density*h;
        // if (isnan(d_density)){
        //     printf("gravity:%.3f, ice ratio:%.3f, water ratio:%.3f, vis:%.3f, temp:%.3f, old density:%.3f, vsfc:%.4f\n", 
        //         top_node->gravity_stress, top_node->ice_ratio, top_node->water_ratio, vis, top_node->temp, top_node->density, top_node->vi_sfc);
        // }
        top_node->thick = top_node->thick*top_node->density/new_density;
        top_node->density_dt = (new_density-top_node->density)/h;
        top_node->density = new_density;
        
        if (top_node->thick>THICKNESS_THRESHOLD){
            top_node->water_ratio = top_node->water_content/density_water/top_node->thick;
            if (top_node->water_ratio < 0){
                top_node->water_ratio = 0;
            }else if(top_node->water_ratio>1){
                top_node->water_ratio = 1;
            }
            top_node->ice_ratio = top_node->ice_content/density_ice/top_node->thick;
            
            if (top_node->ice_ratio < 0){
                top_node->ice_ratio = 0;
            }else if(top_node->ice_ratio>1){
                top_node->ice_ratio = 1;
            }
            top_node->f_l = A1*top_node->ice_ratio*(1-top_node->ice_ratio);  //A1*ice_ratio*(1-ice_ratio)
            top_node->Cs = (top_node->ice_ratio*density_ice*C_i+top_node->water_ratio*density_water*C_l)/top_node->density;
        }
        if (top_node->Cs > C_l){
            top_node->Cs = C_l;
        }else if(top_node->Cs<C_i){
            top_node->Cs = C_i;
        }
        next_node = top_node->lower_layer;
        while(next_node != NULL){
            next_node->gravity_stress = next_node->upper_layer->gravity_stress+next_node->density*next_node->thick*gravity;
            vis = 37000000*exp(0.081*(params->Tf-next_node->temp)+0.018*next_node->density);
            value = next_node->density - 100;
            if (value<0){
                value = 0;
            }
            d_density = next_node->gravity_stress/vis+0.0000028*exp(-0.042*(params->Tf-next_node->temp)-0.046*value);
            d_density = d_density * next_node->density;
            new_density = next_node->density + d_density*h;
            next_node->thick = next_node->thick*next_node->density/new_density;
            next_node->density_dt = (new_density-next_node->density)/h;
            next_node->density = new_density;
            if (next_node->thick>THICKNESS_THRESHOLD){
                next_node->ice_ratio = next_node->ice_content/density_ice/next_node->thick;
                if (next_node->ice_ratio < 0){
                    next_node->ice_ratio = 0;
                }else if(next_node->ice_ratio>1){
                    next_node->ice_ratio = 1;
                }
                next_node->water_ratio = next_node->water_content/density_water/next_node->thick;
                if (next_node->water_ratio < 0){
                    next_node->water_ratio = 0;
                }else if(next_node->water_ratio>1){
                    next_node->water_ratio = 1;
                }
                ice_ratio = next_node->ice_ratio;
                next_node->f_l=  A1*ice_ratio*(1-ice_ratio);  // A1*ice_ratio*(1-ice_ratio)
                next_node->Cs = (next_node->ice_ratio*density_ice*C_i+next_node->water_ratio*density_water*C_l)/next_node->density;
            }
            if (next_node->Cs >C_l){
                next_node->Cs = C_l;
            }else if(next_node->Cs<C_i){
                next_node->Cs = C_i;
            }
            next_node= next_node->lower_layer;
        }    
    }
    node = model->head;
    while(node != NULL){
        next_node = node->lower_layer;
        if (node->thick<=THICKNESS_THRESHOLD || node->ice_ratio<= 0.01){
            remove_node(model, node);
        }else{
            total_thick += node->thick;
        }
        node = next_node;
    }
    model->total_thick = total_thick;

}

double snow_cover_fraction(double snow_depth){
    double snow_cover;
    snow_cover = snow_depth/(snow_depth+0.25);
    // snow_cover = tanh(snow_depth/0.1);
    // snow_cover = snow_depth/0.1;
    // if (snow_cover>1){
    //     snow_cover = 1;
    // }
    return snow_cover;
}

double sensible_heat_efficient(
            double air_temp,
            double surf_temp,
            double wind,
            double snow_depth){
    double RiB;
    double Chn;
    double Fh;
    double C_h;
    if (wind>0){
        RiB = gravity*wind_height*(air_temp-surf_temp)/(air_temp*wind*wind);
        Chn = 0.4*0.4/(log(wind_height/0.01))/(log(air_height/0.1));
        if (RiB>0){
            Fh = 1/(1+2*5*RiB/sqrt(1+RiB));
        }else{
            Fh = 1 - 3*5*RiB/(1+3*5*5*Chn*sqrt(-RiB*wind_height/0.01));
        }
        C_h = Fh*Chn;
        // printf("surf_temp:%.2f\tRiB:%.4f\tChn:%.4f\tFh:%.4f\tC_h:%.4f\n", surf_temp, RiB, Chn, Fh, C_h);
    }else{
        C_h = 0;
    }
    return C_h;
}

double snow_conductivity(double density){
    if (density < 100) return 0.2; // 低密度雪
    if (density > 500) return 0.3; // 高密度雪
    // 线性近似
    return 0.2 + (density - 100) * (0.1 / 400.0);
}
