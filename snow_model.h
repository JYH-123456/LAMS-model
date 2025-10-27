#ifndef SNOW_MODEL_H
#define SNOW_MODEL_H
#include <time.h>
#include <stdbool.h>
#include <gsl/gsl_multifit_nlin.h>
#include <math.h>

// Define constant
#define C_l 4217.6  // specific heat of liquid water, J/kg/K
#define C_i 2106   // specific heat of ice, J/kg/K
#define C_a 1007.6  // specific heat of air, J/kg/K
#define C_g  1200 // predefined specific heat of ground soil, J/kg/K
#define Ld 2834000 // the latent heat of deposition, J/kg
#define Lf 334000  // the latent heat of fusion, J/kg
#define Lv 2500000  //the latent heat of vaporization, J/kg
#define Rv 461  // the specific gas constant for water vapor, Pa/K*m3/kg
#define boltzmann 0.0000000567   // Stefan-Boltzmann constant
#define gravity 9.8  // m/s2
#define unit_g 0.0098  // kpa(kg/m2)-1
#define wind_height  10 // the near-surface reference height where wind is measured, m
#define C_w  4180000//volumetric heat capacity (J·m−3·K−1) of liquid water.

#define param_num 4 // number of grid specific parameters
#define snow_emissity 0.976  // Snow emissivity, dimensionless
#define diffusion_water 0.000022  //the water vapor diffusion coefficient in air, m2/s
#define density_ice 916.68  // the density of solid ice at 0℃, kg/m3
#define density_water 999.84  //the density of liquid water at 0℃, kg/m3
#define density_air 1.14 // density of air, kg/m3
#define density_grnd  2800 // predefined ground soil layer density, kg/m3
#define soil_density 2800 //soil layers density, kg/m3
// #define Ce 0.004  // Bulk transfer coefficient for water vapor, Dimensionless
#define snow_rough   0.005// snow surface roughness, m
#define SNOW_CONDUCT 0.05  // thermal conductivity of snow, W·m−1·K−1;
#define SOIL_CONDUCT 2 // soil heat conductivity, W·m−1·K−1;
#define air_cond 0.2  // thermal conductivity of air, W·m−1·K−1;
#define air_height 1.5

// key adjustable parameters
#define SNOW_THRESHOLD 0.1 // thickness threshold beyond which a new node is added, m
#define TIME_THRESHOLD 24 // in number of time steps, 24
#define THICKNESS_THRESHOLD 0.001 // thickness threshold when a node is removed, m
#define thick_min 0.001 // the minumum thickness for snow model
#define delta_t 1  // time step in hours
#define T_max_snow 0  // maximum temperature for snowfall, C
#define T_min_rain -2   //minimum temperature for rainfall, C
#define ALBEDO_K 0.02 //albedo coefficient for emprical function
#define LAMDA 0.05  // solar absorption length in snow, m
#define stress_cr1 0.2  // kpa
#define stress_cr2 56.234  // kpa
#define rho_snow_cr1 58.7048  // kg/m3
#define rho_snow_cr2 173.7178  // kg/m3
#define A1 0.3 // empirical parameter a1
#define wind_hight 10 // wind speed hight
#define SNOW_NEW_SNOW_DENSITY 100 // new snow density, kg/m3
#define SNOW_NEW_BRAS_DENOM 100 // parameter used to calculate new snow density when air temp>0, kg/m3
#define Density_min 70
#define ATTEN_k 0.01  //attenuation coefficient for snow albedo
#define SURF_DEP  1 // soil surface depth, m
#define CONST_DEP 5  // m
#define C_TO_F(t) (t * 9. / 5.) + 32.
#define Tcri 273.15  // maximum temperature for snow layer

#define water_pot_bottom  0.1// underwater saturation water potential
#define param_c 0.4 //soil texture related parameters
#define param_d -0.38 //soil texture related parameters
#define CONST_T2  270  // constant soil temperature at CONST_DEP, K

// #define underwater_prosity 0.2  //underwater soil porosity == water content
// #define underwater_water_pot 0.1 // water potential for subsurface water
// #define underwater_k 0.7 // subsurface water soil layer hydraulic conductivity, m/s
// #define underwater_b 0.5 // parameter b for groundwater soil layer

typedef struct {
    int *years;
    int *months;
    int *days;
    int *hours;
    int count;
} DateList;

// typedef struct veg_params
// {
//     double LAI[12];
//     double canopy;
//     double aero_resist;
//     double roughness;
//     double max_intercept;
//     double Tcanopy;
//     double inter_rain;
//     double inter_snow;
// }veg_params;

// typedef struct soillayers{
//     int num_grids;  // invariable
//     double density;  // time-invariant, kg/m3
//     double soil_depth;  // invariable
//     double *soil_temp;  // in K
//     double *flow_heat;  // heat transfer in W/m2 
//     double *flow_mass;  // in m/s
//     double *thermal_cond;  //thermal conductivity (W·m−1·K−1) 
//     double *heat_cap;  //volumetric heat capacity, J·m·−3·K−1
//     double *water_content; // volumetric liquid water content (m3/m3) 
//     double *porosity;  // soil porosity (m3/m3)
//     double *hydro_cond; //hydraulic conductivity, m/s
//     double *water_pot;  //soil water potential (m) 
//     double *pram_b;  // soil related paramter for saturated hydraulic conductivity water potential calculation
//     double *D;  // intermediate function
//     double *water_cont_dt;  // water content change during dt
//     double *satur_k;  // saturated hydraulic conductivity, m/s
//     double *satur_pot;  //satuated water potential (m)
//     double *ice_cont;  // soil ice content, (m3/m3)
//     double *volumn_fra; // volumetric portion of liquid water existing in frozen soil
//     double *frozen;  // frozen water fraction in soil, (m3/m3)
//     double *melt;  // melt water fraction in soil, (m3/m3)
//     double water_ratio_top; // interface water content
//     double water_flow_bottom; // groundwater water flow, in m/s
//     double T_0;  // interface temperature, if snow not exist, Tair, K
//     double Q_0;  // interface heat transfer, in W/m2
// } soillayers;

typedef struct atmos_array {
    int size;
    double *rain; // snow amount in mm
    double *air_temp; // in ℃
    double *snow;   // snow amount in mm
    double *in_short; //incoming short wave radiation
    double *in_long; //incoming short wave radiation
    double *wind; // wind speed at 10m, m/s
    double *air_pressure; //air pressure, in kPa
    double *vapor_pressure;  // vapor pressure, in kPa
    double *humidity;
    double *Tg1;  // soil temperature in 0.1m depth, C
    double *Tg2;  // soil temperature in 0.2m depth, C
    double *Tg3;  // soil temperature in 0.5m depth, C 
} atmos_array;


typedef struct {
    double Tpmin;
    double Tpmax;
    double C_h;
    double albedo_max;
    double albedo_min;
    double Tf;
} ParamData;


typedef struct SnowLayerNode {
    int layer_id;
    double thick; // in cm
    double density;   // in kg/m^3
    double density_dt; // in kg/m^3/t
    double ice_ratio;  // volume fraction
    double water_ratio;  // volume fraction
    double water_content;  // mass fraction,kg/m2/s
    double ice_content;  // mass fraction,kg/m2/s
    double gravity_stress;
    double temp;  // in K
    int snow_time; // time when last snow happened, in time step
    double Cs;  // specific heat of snow
    struct SnowLayerNode *upper_layer;
    struct SnowLayerNode *lower_layer;
    double short_in;   // short radiation transported from upper layer
    double short_down;  // short radiation emitted from current layer
    double long_in; // incoming longwave radiation;
    double refreeze;      // liquid water refreeze in current layer
    double water_out;  // liquid water flow to the next layer, in kg/m^2/s;
    double water_in; // liquid water come from upper layer, for inner_layers;
    double rain_in;  // rainfall amount come from prec, only for top node, in m/s;
    double ice_in;  //snowfall amount come from prec, only for top node, in m/s;
    double f_l;  // water holding capacity of current layer;
    double vi;
    double vi_sfc;
    double vl_sfc;
    double delta_temp; // temperature change by a timestep
} SnowLayerNode;


typedef struct SnowModel{
    SnowLayerNode *head;
    SnowLayerNode *tail;
    int num_layers;
    double total_thick;
    double ave_thick;
    double snow_cover;
    double Tgrnd;
    double Tsfc;  // snow surface temperature
    double outflow;  // infiltrated liquid water, in kg/m2/s
    double T1; // snow node tempe, K
    double fi;  // mass content of ice, kg/m2
    double fl;  // mass content of water, kg/m2
    double density;
    double temp;
    double refreeze;
    double melt;
    double swe; //m
    double albedo; 
    double thick_add;
} SnowModel;

// Function declarations
void init_model(SnowModel *model);
void add_node(SnowModel *model, atmos_array *atmos, ParamData *params, double snowfall, double rainfall, double h, int current_step);
void remove_node(SnowModel *model, SnowLayerNode *node);
void print_layers_to_file(const SnowModel *model, FILE *file);

double f_top(SnowLayerNode *node, atmos_array *atmos, int t_count, double layer_temp, double Tgrnd, double albedo, double C_h, int print);
double f_inner(SnowLayerNode *node, double layer_temp, double Tgrnd, int print);
double calculate_Tgrnd(SnowModel *model, ParamData *params, double shortwave, double snowfall, double rainfall, double airtemp, double wind, 
                        double et, double dt, int month, int print);
double calculate_keff(double snow_density);
double calculate_vi(double density_dt, double thick, double snow_density);
double calculate_sfc_v(SnowLayerNode *node, double layer_temp, double wind, double air_pressure, double vapor_pressure, double C_h, int snow_mode);
double calculate_q(double Tsfc, int snow_mode, double air_pressure);
double calculate_es(double temp, int snow_mode);
double calculate_qair(double e, double air_pressure);
void update_snow_layer(SnowModel *model, atmos_array *atmos, ParamData *params, DateList datelist, int t_count, double h);

double rain_snow_seperate(ParamData *params, double prec_amount, double air_temp, double humidity, double air_pressure);
void update_prec(SnowModel *model, atmos_array *atmos,  ParamData *params, int t_count, double h);
void mass_balance(SnowModel *model, double h, int pr);
void mass_balance_step(SnowLayerNode *node, double dt, int top);
void calculate_thick(SnowModel *model, double dt);

double snow_albedo(SnowModel *model,
                SnowLayerNode *top_node,
                atmos_array *atmos,
                ParamData *params,
                double h,
                int t_count);
void snow_compaction(SnowModel *model,
                    ParamData *params,
                    double h);
double snow_cover_fraction(double snow_depth);
double sensible_heat_efficient(
            double air_temp,
            double surf_temp,
            double wind,
            double snow_depth);
double snow_conductivity(double density);
double new_snow_density(double air_temp, double wind, ParamData *params);
// double new_snow_density(double air_temp);

int initialize_atmos_array(atmos_array *atmos, int num);
int read_atmos_data(atmos_array *atmos, const char *filename);
void free_atmos_array(atmos_array *atmos);
int read_param_file(const char *filename, ParamData *params);

int calculate_time_steps(struct tm start, struct tm end, int hour_interval);
DateList generate_datetime_list(int start_year, int start_month, int start_day, int start_hour,
                                   int end_year, int end_month, int end_day, int end_hour,
                                   int hour_interval);
void free_datetime_list(DateList* datetime_list);

// void snow_run(SnowModel *model, atmos_array *atmos, int t_count, double h, int dt);
// void snow_run_step(SnowModel *model, atmos_array *atmos, int t_count, double h, int dt);
// int expb_f(const gsl_vector *x, void *data, gsl_vector *f);
// int fit_paramter(double *T, double *Z, size_t n, double *params);
extern FILE *outputFile;
void initialize_output_file(const char *filename);
void close_output_file();
double snow_density(double air_temp);

#endif // SNOW_MODEL_H
