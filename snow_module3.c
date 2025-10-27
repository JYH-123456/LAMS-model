#include "snow_model.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// mass transfer related functions
void mass_balance(
    SnowModel *model, 
    double h,
    int pr){
    SnowLayerNode *top_node;
    SnowLayerNode *current_node;
    SnowLayerNode *node;
    SnowLayerNode *next_node;
    double total_thick;
    double refreeze;
    double vi;
    double vi_sfc;
    double vl_sfc;
    double thick;
    double fi;
    double fl;
    // calculate mass balance for top layer/node
    if (model->num_layers>0){
        top_node = model->head;
        thick = top_node->thick;
        refreeze = top_node->refreeze; 
        vi = top_node-> vi;
        vl_sfc = top_node->vl_sfc;
        vi_sfc = top_node->vi_sfc;
        fi = (top_node->ice_in*density_water+refreeze-vi-vi_sfc)*h + top_node->ice_content;
        top_node->ice_content = fi;
        if (fi<0){
            top_node->ice_content = 0;
        }
        fl = (top_node->rain_in*density_water-top_node->water_out-refreeze-vl_sfc)*h + top_node->water_content;  //top_node->rain_in*density_water
        top_node->water_content = fl;
        // printf("top node vi_sfc: %.4f, vl_sfc:%.4f\n", vi_sfc*1000000, vl_sfc*1000000);
        if (fl<0) {
            top_node->water_content = 0;
        }
        // if (top_node->ice_content>0){
        //     top_node->f_l = A1*top_node->ice_ratio*(1-top_node->ice_ratio);
        // }else{
        //     top_node->f_l = 0;
        // }

        // calculate mass balance for inner layer/node
        current_node = top_node->lower_layer;
        while (current_node != NULL){
            refreeze = current_node->refreeze; 
            fi = (current_node->refreeze-current_node->vi)*h+current_node->ice_content;
            current_node->ice_content = fi;
            if (fi<0){
                current_node->ice_content = 0;
            }
            fl = (current_node->water_in-current_node->water_out-current_node->refreeze)*h+current_node->water_content;
            current_node->water_content = fl;
            if (fl<0){
                current_node->water_content =0;
            }
            // if (current_node->ice_ratio>0){
            //     current_node->f_l = A1*current_node->ice_ratio*(1-current_node->ice_ratio);
            // }else{
            //     current_node->f_l = 0;
            // }
            current_node = current_node->lower_layer;
        }
    }  
}

// void mass_balance_step(
//         SnowLayerNode *node, 
//         double dt,
//         int top){
//     double ice_ratio;
//     double water_ratio;
//     double refreeze = node->refreeze;
//     double thick = node->thick;
//     double density = node ->density;
//     double vi = node-> vi;
//     if (top==1){
//         double vl_sfc = node->vl_sfc;
//         double vi_sfc = node->vi_sfc;
//         ice_ratio = (refreeze-vi-vi_sfc)/thick/density_ice*dt+node->ice_ratio;
//         water_ratio = (node->rain_in-node->water_out-refreeze-vl_sfc)/density_water/thick*dt+node->water_ratio;    
//     }else{
//         ice_ratio = (node->refreeze-node->vi)/node->thick/density_ice*dt+node->ice_ratio;
//         water_ratio = (node->water_in-node->water_out-node->refreeze)/node->thick/density_water*dt+node->water_ratio;
//     }
//     if (ice_ratio<0){
//             node->ice_ratio = 0;
//         }else{
//             node->ice_ratio = ice_ratio;
//         }
//     if (water_ratio < 0){
//         node->water_ratio = 0;
//     }else{
//         node->water_ratio = water_ratio;
//     }
//     node->f_l = A1*node->ice_ratio*(1-node->ice_ratio);
//     node->air_ratio = 1-node->ice_ratio-node->water_ratio;
//     node->density = node->ice_ratio*density_ice+node->air_ratio*density_air+node->water_ratio*density_water;
//     node->Cs=(node->ice_ratio*density_ice*C_i+node->water_ratio*density_water*C_l)/node->density;
// }


double calculate_stress_threshold(double snow_density){
    double stress_f;
    if (snow_density<=500){
        stress_f = stress_cr1*exp(snow_density/rho_snow_cr1);
    }else{
        stress_f = stress_cr2*exp(snow_density/rho_snow_cr2);
    }
    return stress_f;
}

// snow heat trasnfer related functions
double f_top(SnowLayerNode *node, 
        atmos_array *atmos, 
        int t_count, 
        double layer_temp, 
        double Tgrnd,
        double albedo,
        double C_h,
        int print) {
    // atmospheric variables
    double temp_prec = atmos->air_temp[t_count]+273.15;  // in unit K
    double in_short = atmos->in_short[t_count];  // in unit W/m2
    double in_long = atmos->in_long[t_count]; // in unit W/m2
    double wind = atmos->wind[t_count]; // wind speed at 10 m, m/s
    double air_temp = atmos->air_temp[t_count]+273.15; // in unit K
    double air_pressure  = atmos->air_pressure[t_count]*1000; // in Pa
    double vapor_pressure = atmos->vapor_pressure[t_count]*1000;  //in Pa
    // node layer related varaibles
    double temp_grad;
    double lower_layer_temp;
    double temp_grad_upper;
    double thick = node->thick;
    double heat_net;
    double vi;
    double keff;
    double V_sfc;
    double V_sfc_heat;
    double fl;  //water ratio
    double snow_density = node->density;
    double Cs = node-> Cs;
    double out_value; 
    int snow_mode;
    double right_term;
    double heat_grad=0;
    double snow_cond = snow_conductivity(node->density);
    double density_dt = node->density_dt;
    if (node->lower_layer==NULL){
        lower_layer_temp = Tgrnd;
        heat_grad = snow_cond * (Tgrnd+air_temp-2*layer_temp) / thick;
    } else{
        lower_layer_temp = node->lower_layer->temp;
        // keff = calculate_keff(snow_density);
        heat_grad = snow_cond*(lower_layer_temp+air_temp-2*layer_temp)/thick;
    }
    temp_grad = (layer_temp-lower_layer_temp)/thick;
    
    vi = calculate_vi(density_dt, thick, snow_density);
    node -> vi = vi;
    fl = node->water_ratio;
    if (fl>0){
        snow_mode = 0;
        V_sfc = calculate_sfc_v(node, layer_temp, wind, air_pressure, vapor_pressure, C_h, snow_mode);
        node->vl_sfc = V_sfc;
        node ->vi_sfc = 0;
        V_sfc_heat = V_sfc*(layer_temp*C_l-Lv);
    }
    else{
        snow_mode = 1;
        V_sfc = calculate_sfc_v(node, layer_temp, wind, air_pressure, vapor_pressure, C_h, snow_mode);
        node->vi_sfc = V_sfc;
        node->vl_sfc = 0;
        V_sfc_heat = V_sfc*(layer_temp*C_i-Ld);
    }
    double heat_prec = (C_l*density_water*node->rain_in+C_i*density_water*node->ice_in)*(temp_prec-layer_temp);
    double heat_vi = -vi*(Ld-layer_temp*C_i);
    double heat_solar = (in_short-albedo*in_short)*(1-exp(-thick/LAMDA));
    double heat_h = density_air*C_a*C_h*wind*(air_temp-layer_temp);
    double heat_lw_out = 0.976*0.0000000567*pow(layer_temp, 4);
    double heat_lw_in = in_long;
    double heat_vsfc = V_sfc_heat;
    double heat_melt = node->refreeze*(Lf+layer_temp*(C_l-C_i));
    right_term = heat_prec+heat_solar+heat_h+heat_vsfc+heat_vi+heat_grad+heat_melt+heat_lw_in-heat_lw_out;  
    out_value = right_term/snow_density/Cs/thick;
    // if (print==1){
    //     fprintf(outputFile, "%.4f\t%.6f\t%.7f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.7f\n",
    //     layer_temp, Tgrnd, node->refreeze, heat_melt, heat_prec, V_sfc, heat_grad, heat_h, heat_lw_in-heat_lw_out, heat_vsfc, heat_vi, heat_solar, out_value);
    //     // fprintf(outputFile, "top layer temp:%.4f\tground temp:%.6f\trefreeze:%.7f\tsnow density:%.3f\tsnow depth:%.3fcm\n",
    //     // layer_temp, Tgrnd, node->refreeze, node->density,node->thick*100);
    //     // fflush(outputFile); 
    // }
    return out_value;
}

double f_inner(SnowLayerNode *node, 
        double layer_temp, 
        double Tgrnd,
        int print){
    if (node == NULL) {
        printf("Error: Node is NULL\n");
    }
    // 然后检查 node->upper_layer 是否为 NULL
    if (node->upper_layer == NULL) {
        printf("Error: upper_layer is NULL for layer_id %d\n", node->layer_id);
        // 处理上层节点为空的情况，比如设置默认值或返回错误码。
    }
    double lower_layer_temp;
    double right_term;
    double water_in = node->water_in;
    double temp_in = node->upper_layer->temp;
    double vi;
    double keff;
    double short_in = node->short_in;
    double thick = node->thick;
    double out_value;
    double snow_density = node->density;
    double Cs = node-> Cs; 
    double heat_grad;
    double heat_grnd;
    double density_dt = node->density_dt;
    double heat_melt = node->refreeze*(Lf+layer_temp*(C_l-C_i));
    if(node->lower_layer==NULL){
        // keff = calculate_keff(snow_density);
        heat_grad = SNOW_CONDUCT*(temp_in+Tgrnd-2*layer_temp)/thick;
    }else{
        lower_layer_temp = node->lower_layer->temp;
        // keff = calculate_keff(snow_density);
        heat_grad = SNOW_CONDUCT*(temp_in+lower_layer_temp-2*layer_temp)/thick;
    }
    
    vi = calculate_vi(density_dt, thick, snow_density);
    node->vi = vi;
    double heat_vi = -vi*(Ld-layer_temp*C_i);
    double heat_water = C_l*water_in*(temp_in-layer_temp);
    
    double heat_solar = fabs(short_in)*(1-exp(-thick/LAMDA));
    right_term = heat_water+heat_solar+heat_vi+heat_grad+heat_melt;
    out_value = right_term/snow_density/Cs/thick;
    // if (print==1){
    //     fprintf(outputFile, "layerid: %d\ttemp:%.3f\tground temp:%.3f\theat_vi:%.7f\theat_water:%.7f\theat_grad:%.7f\theat_solar:%.7f\theat melt:%.7f\tout_value:%.7f\n", 
    //     node->layer_id, layer_temp, Tgrnd, heat_vi, heat_water, heat_grad, heat_solar, heat_melt, out_value);
    //     fflush(outputFile); 
    // }
    return out_value;
}


double min(double a, double b) {
    return (a < b) ? a : b;
}

double calculate_keff(double snow_density){
    double keff;
    keff = 0.021+2.5/1000000*pow(snow_density, 2);
}

double calculate_vi(double density_dt,
        double thick,
        double snow_density){
    double vi;
    double density_ratio;
    density_ratio = 1-snow_density/density_ice;
    double De = diffusion_water*10*pow((1-density_ratio),0.51);
    vi = -De*density_dt*thick; //
    if (isnan(vi)){
        printf("De:%.4f, density_dt:%.4f, porosity:%.4f, thick:%.4f\n", De, density_dt, density_ratio, thick);
    }
    return vi;
}

double calculate_sfc_v(SnowLayerNode *node, 
        double layer_temp,
        double wind, 
        double air_pressure, 
        double vapor_pressure,
        double C_h, 
        int snow_mode){
    double V_sfc;
    double V_max;
    double Tsfc = layer_temp;
    double qair; //the specific humidity of the air at height 2 m
    double q; //the specific humidity of snow surface layer
    qair = calculate_qair(vapor_pressure, air_pressure);
    q = calculate_q(Tsfc, snow_mode, air_pressure);
    V_sfc = density_air*C_h*wind*(q-qair);
    // if (fabs(qair-q)<=0.0001){
    //     printf("Tsfc:%.4f, qair:%.4f, q:%.4f\n", layer_temp, qair, q);
    // }
    if(snow_mode==0){
        V_max = density_water*node->thick*node->water_ratio/delta_t/3600;
    }else{
        V_max = density_ice*node->thick*node->ice_ratio/delta_t/3600;
    }
    // if (V_sfc<0){
    //     V_sfc=0;
    // }
    if(V_sfc>V_max){
        V_sfc = V_max;
    }
    // printf("qair:%.7f, q:%.7f\n", q, qair);
    return V_sfc;  // unit: kg/m2/s
}

double calculate_q(double Tsfc, 
        int snow_mode, 
        double air_pressure){
    double es;
    double q;
    es = calculate_es(Tsfc, snow_mode);
    q = 0.622*es/(air_pressure-(1-0.622)*es);
    return q;
}

double calculate_es(double temp, // in K
        int snow_mode){
    double es;
    double log_es;
    double temp_c = temp-273.15;
    if (snow_mode) {
        // Goff-Gratch equation
        if(temp<273.15){
            log_es = -9.09718*(273.16/temp-1)-3.56654*log10(273.16/temp)+0.876793*(1-temp/273.16)+log10(6.1071);
        }else{
            log_es = -7.90298*(373.16/temp-1)+5.02808*log10(373.16/temp)+
            0.00000013816*(pow(10, 11.344*(1-temp/373.16))-1)+
            0.0081328*(pow(10, -3.49149*(373.16/temp-1))-1) +
            log10(1013.246);
        }
        es = pow(10, log_es);
    }
    else{
        // Magnus-Tetens equation
        es = 6.1121*exp(17.67*temp_c/(temp_c+243.5));
    }
    return es*100;  // in Pa
}

double calculate_qair(double e, 
        double air_pressure){
    double qair;
    qair = 0.622*e/(air_pressure-(1-0.622)*e);
    return qair;
}

// Fourth Order Runge-Kutta method
void update_snow_layer(SnowModel *model, 
        atmos_array *atmos,
        ParamData *params, 
        DateList datelist,
        int t_count, 
        double h) {
    SnowLayerNode *current_node;
    SnowLayerNode *top_node;
    SnowLayerNode *node;
    SnowLayerNode *lower_node;
    int month = datelist.months[t_count];
    double layer_temp;
    double refreeze;
    double k1;
    double k2;
    double k3;
    double k4;
    double short_down;
    double albedo;
    double ml_max;
    double fi_mix=0;
    double fl_mix=0;
    double density_mix=0;
    double temp_mix=-1;
    double ml;
    double new_temp;
    int t=0;
    int pr = 0;
    int count = 0;
    double ddt = 1;
    double total_refreeze= 0;
    double total_melt = 0;
    double res_tol=999;
    double melt_max;
    double last_temp;
    double refreeze_max;
    int num_layer = model->num_layers;
    double store_temp[num_layer];
    double max_temp;
    double min_temp;
    double C_h;
    double Tgrnd = atmos->Tg1[t_count]+273.15;  // convert C to K
    model->Tgrnd = Tgrnd;
    // store current temp
    node = model->head;
    for (int i = 0; i < num_layer && node != NULL; i++) {
        store_temp[i] = node->temp;
        node = node->lower_layer;
    }

    // assign atmospheric data
    double air_temp = atmos->air_temp[t_count]+273.15;
    double humidity = atmos->humidity[t_count];
    double air_pressure = atmos->air_pressure[t_count];
    double wind = atmos->wind[t_count];
    double rainfall = atmos->rain[t_count]/1000;
    double snowfall = atmos->snow[t_count]/1000;
    double last_snow_depth = model->total_thick;
    double if_melt=0;
    if (last_snow_depth>0){
        if_melt = -model->total_thick*model->density/h;
    }

    if (num_layer>0 && model->total_thick>thick_min){
        if (air_temp<273.15 && model->total_thick>0.1){ // air_temp<273.15 && 
            albedo = params->albedo_max;
        }else{
            albedo = params->albedo_min;
        }
        // if (wind>3){
        //     C_h = 0.0005;
        // }else{
        //     C_h = 0.0005;
        // }
        C_h = params->C_h;
        // albedo = snow_albedo(model, model->head, atmos, params, h, t_count);
        // C_h = sensible_heat_efficient(air_temp, model->Tsfc, wind, model->total_thick);
        model-> albedo = albedo; 
        printf("t_count:%d, albedo:%.3f, C_H:%.5f\n", t_count, albedo, C_h);  
        // calculate temperature for top layer
        top_node = model->head;                                     
        top_node -> long_in = atmos->in_long[t_count];
        top_node -> short_in = atmos->in_short[t_count];
        short_down = (1-albedo)*atmos->in_short[t_count] * exp(-top_node->thick/LAMDA);
        top_node -> short_down = short_down;
        Tgrnd = model->Tgrnd;
        while (t < h){
            layer_temp = top_node->temp;
            // update top layer temperature by Fourth-order Runge-Kutta
            k1 = ddt*f_top(top_node, atmos, t_count, layer_temp, Tgrnd, albedo, C_h, 0);
            k2 = ddt*f_top(top_node, atmos, t_count, layer_temp + 0.5 * k1, Tgrnd, albedo, C_h, 0);
            k3 = ddt*f_top(top_node, atmos, t_count, layer_temp + 0.5 * k2, Tgrnd, albedo, C_h, 0);
            k4 = ddt*f_top(top_node, atmos, t_count, layer_temp + k3, Tgrnd, albedo, C_h, 0);
            new_temp = layer_temp + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
            max_temp = melt_max*Lf*h/top_node->density/top_node->Cs/top_node->thick+params->Tf;
            if (isnan(new_temp) || new_temp>280 || new_temp<223){  // || new_temp<223
                new_temp = layer_temp;
            }
            top_node->temp = new_temp;
            // calculate top layer refreeze flux
            refreeze_max = density_water*top_node->thick*top_node->water_ratio/h;
            melt_max = density_ice*top_node->thick*top_node->ice_ratio/h;
            refreeze = top_node->density*top_node->Cs*top_node->thick*(params->Tf-new_temp)/Lf/h; //273.15
            top_node->refreeze = refreeze;
            if (refreeze<0 && fabs(refreeze)>melt_max){
                top_node->refreeze = -melt_max;
            }
            if (refreeze>0 && refreeze>refreeze_max){
                top_node->refreeze = refreeze_max;
            }
            // calculate down flow 
            ml_max = density_water*top_node->thick*top_node->f_l/h;
            ml = density_water*top_node->thick*top_node->water_ratio/h-top_node->refreeze;
            if (ml>ml_max){
                top_node->water_out = ml - ml_max;
            }else{
                top_node->water_out = 0;
            }    
            // total_refreeze_flux += refreeze*(Lf-layer_temp*C_i+layer_temp*C_l);
            current_node = top_node->lower_layer;
            while(current_node != NULL) {
                layer_temp = current_node->temp;
                current_node->short_in = current_node->upper_layer->short_down;
                current_node->water_in = current_node->upper_layer->water_out;
                current_node->short_down = current_node->short_in*exp(-current_node->thick/LAMDA);
                // update temperature
                k1 = ddt * f_inner(current_node, layer_temp, Tgrnd, 0);
                k2 = ddt * f_inner(current_node, layer_temp + 0.5 * k1, Tgrnd, 0);
                k3 = ddt * f_inner(current_node, layer_temp + 0.5 * k2, Tgrnd, 0);
                k4 = ddt * f_inner(current_node, layer_temp + k3, Tgrnd, 0);
                new_temp = layer_temp+(k1 + 2*k2 + 2*k3 + k4) / 6.0;
                
                if (isnan(new_temp) || new_temp>280|| new_temp<223){
                    new_temp = layer_temp;
                }
                current_node->temp = new_temp;
                // calculate refreeze mass flux
                refreeze_max = density_water*current_node->thick*current_node->water_ratio/h;
                melt_max = density_ice*current_node->thick*current_node->ice_ratio/h;
                refreeze = current_node->density*current_node->Cs*current_node->thick*(params->Tf-new_temp)/Lf/h;
                current_node->refreeze = refreeze;
                if (refreeze<0 && fabs(refreeze)>melt_max){
                        current_node->refreeze = -melt_max;
                }
                if (refreeze>0 && refreeze>refreeze_max){
                        current_node->refreeze = refreeze_max;
                }
                // calculate down flow
                ml_max = density_water*current_node->thick*current_node->f_l/h;
                ml = density_water*current_node->thick*current_node->water_ratio/h-current_node->refreeze;  //+current_node->water_in
                if (ml>ml_max){
                    current_node->water_out = ml - ml_max;
                }
                else{
                    current_node->water_out = 0;
                }
                current_node = current_node->lower_layer;            
            }
            //calculate temperature for ground layer
            t += ddt;
            //printf("initial temp:%.5f\tresidual:%.5f\trefreeze:%.5f\tground temp:%.5f\n", atmos->temp_prec[t_count]+273.15,residual, refreeze, Tgrnd);
        }
       
        // update snow layer temp and calculate total refreeze and melt flux
        current_node = model->head;
        int k = 0;
        while(current_node != NULL){
            if(isnan(current_node->temp)){
                current_node->temp= store_temp[k];
            }
            if(current_node->temp>params->Tf){
                current_node->temp = params->Tf;
            }
            if(current_node->refreeze>0){
                total_refreeze += current_node->refreeze;
            }else{
                total_melt += current_node->refreeze;
            }
            current_node = current_node->lower_layer;
            k++;
        }
        model->refreeze = total_refreeze;
        model->melt = total_melt;
        model->outflow += model->tail->water_out*h;

    }else{
        model->refreeze = 0;
        model->melt = 0;
        model->outflow += 0;        
    }
    // printf("time step:%d\tTgrnd:%.3f\n", t_count, model->Tgrnd);
    mass_balance(model, h, 0);
    snow_compaction(model, params, h);

    // update snow characteristic as a whole
    current_node = model->head;
    if (current_node == NULL){
        model->fl = 0;
        model->fi = 0;
        model->density = 0;
        model->temp = -1;
        model->swe = 0;
        model->Tsfc = -1;
        
    }else{
        model->Tsfc = model->head->temp;
        double pack_temp=0;
        double swe = 0;
        while (current_node != NULL){
            fi_mix += current_node->thick*current_node->ice_ratio/model->total_thick;
            fl_mix += current_node->thick*current_node->water_ratio/model->total_thick;
            swe += current_node->thick*current_node->density/density_water;
            current_node = current_node -> lower_layer;
        }
        model->swe = swe;
        density_mix = fi_mix*density_ice+fl_mix*density_water;
        model->fi = fi_mix;
        model->fl = fl_mix;
        model->density = density_mix;
        double C_mix = fi_mix*density_ice*C_i+fl_mix*density_water*C_l;  //+fa_mix*density_air*C_a

        current_node = model->head;
        while (current_node != NULL){
            pack_temp += current_node->thick*current_node->Cs*current_node->temp/model->total_thick/C_mix;
            current_node= current_node->lower_layer;
        }
        model->temp = pack_temp;
    }
}
