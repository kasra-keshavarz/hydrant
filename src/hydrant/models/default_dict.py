"""
Default dictionaries
"""

# header
header= {
    'section_head' :"""!!	=======================================================================================================																
!! Parameter file for: Case Name																	
!!	------------------------																
!!																	
!!	=======================================================================================================																"""
}

# general met param
General_meteo_param = {
    'section_head' :"""!! METEOROLOGICAL INPUT - general parameters related to temperature and precipitation corrections																	
!!	-----																
!! All of these will be kept as 0 because we are not correcting the temperature or the precpitation																""",
    'tcobselev'    :{'value': 0, 'comment': "!! parameter for temperature correction due to observation elevation deviation from subbasin elevation (deg C)"},
    'tcalt'        :{'value': 0, 'comment': "!! parameter for temperature’s elevation dependence"},
    'tcelevadd'    :{'value': 0, 'comment': "!! parameter for temperature’s elevation dependence"},
    'pcaddg'       :{'value': 0, 'comment': "!! correction parameter for precipitation"},
    'pcusnow'      :{'value': 0, 'comment': "!! undercatch correction for snowfall"},
    'pcurain'      :{'value': 0, 'comment': "!! undercatch correction for rainfall"}
}

# snow param
Snow_param = {
    'section_head' :"""!!	=======================================================================================================																
!!	"SNOWMELT, DISTRIBUTION, DENSITY, HEAT; sublimation is sorted under Evapotranspiration"																
!!	-----																
!!	"General snow accumulation and melt related parameters (baseline values from SHYPE, unless noted otherwise)"																
!! Snow distribution submodel: 0, snow melt submodel: 2, snow density submodel: 0, snow heat submodel: 1 (including snkika below)""",
    'ttpi'         :{'value': 1.7083, 'comment': "!! half of temperature interval with mixed snow and rainfall"},
    'sdnsnew'      :{'value': 0.13,   'comment': "!! density of new-fallen snow (kg/dm3)"},
    'snowdensdt'   :{'value': 0.0016, 'comment': "!! increase of snow density per day"},
    'fsceff'       :{'value': 0.99,   'comment': "!! efficiency of snow cover to influence snow melt and snow evaporation"},
    'cmrefr'       :{'value': 0.05,   'comment': "!! refreeze efficiency compared to the degree-day snow melt factor Used for second snow melt model"},
    'whcsnow'      :{'value': 0.08,   'comment': "!! water holding capacity of snow"}
}

# Regional or deep Groundwater outflow
Deep_groundwater = {
    'section_head' :"""!!	-----																
!!	Regional or deep groundwater outflow																
!! for deepground submodel: 0""",
    'rcgrw'        :{'value': 0.0,
                     'comment': "!! recession coefficient for regional groundwater outflow from soil layers (set to zero because we are not considering deepground)"}
}


# river routing
River_routing ={
    'section_head' :"""!!																	
!!	=======================================================================================================																
!!	RIVER ROUTING																
!!	-----			
!! for riverflow submodel:1?""",
    'damp'         :{'value': 0.1614,
                     'comment': "!! fraction of delay in the watercourse which also causes damping"},
    'rivvel'       :{'value': 9.9267,
                     'comment': "!! celerity of flood in watercourse"},
    'qmean'        :{'value': 200,
                     'comment': "!! initial value for calculation of mean flow "}
}


# Glacier
Glacier  = {
    'section_head' :"""!!	=======================================================================================================																
!!	"GLACIER - parameters for volume-area scaling, accumulation and melt (sublimation sorted under Evapotranspiration)"																
!!	-----																
!!	Glacier volume-area scaling	
!! the parameters used calculate the area of the glaciers""",
    'glacvexp'     :{'value'  : 1.395    , 'comment': "!! exponent of glacier area-volume relationship for glacier of type zero"},
    'glacvcoef'    :{'value'  : 0.17157  , 'comment': "!! coefficient of glacier area-volume relationship for glacier of type zero"},
    'glacvexp1'    :{'value'  : 1.25     , 'comment': "!! exponent of glacier area-volume relationship for glacier of type one)"},
    'glacvcoef1'   :{'value'  : 2.88364  , 'comment': "!! coefficient of glacier area-volume relationship for glacier of type one"},
    'glac2arlim'   :{'value'  : 25000000 , 'comment': "!! area limit for determine glacier type which is used only if glacier type is given in GlacierData.txt"},
    'glacannmb'    :{'value'  : 0        , 'comment': "!! annual mass balance for correction of initial glacier volume"}
}

# Glacier melt
Glacier_melt = {
    'section_head' :"""!!	-----																
!!	Glacier melt parameters 																
!!	----																
!! considered with snowevaporation submodel: 1, snowmelt submodel 2""",
    'glacttmp'     :{'value'  : 0          , 'comment': "!! threshold temperature for glacier melt"},
    'glaccmlt'     :{'value'  : 1.58595482 , 'comment': "!! melting parameter for glacier"},
    'glaccmrad'    :{'value'  : 0.19090136 , 'comment': "!! coefficient for radiation glacier melt parameter for second snowmelt model"},
    'glaccmrefr'   :{'value'  : 0.06259448 , 'comment': "!! refreeze efficiency compared to the degree-day glacier melt factor parameter for second snow meltmodel"},
    'glacalb'      :{'value'  : 0.35       , 'comment': "!! albedo for glacier ice"},
    'fepotglac'    :{'value'  : 0          , 'comment': "!! fraction of snow-free potential evapotranspiration for first snowevaporation model"}
} 

# evaporation
Evap = {
    'section_head' :"""!!	=======================================================================================================																
!!	EVAPOTRANSPIRATION PARAMETERS																
!!	-----																
!!	General evapotranspiration parameters																
!! used for petmodel""",
    'lp'           :{'value'  : 0.546 , 'comment': "!! Threshold for water content reduction of transpiration as fraction of field capacity"},
    'epotdist'     :{'value'  : 0.546 , 'comment': "!! Coefficient in exponential function for potential evapotranspiration's depth dependency"},
    'krs'          :{'value'  : 0.546 , 'comment': "!! parameter for estimating shortwave radiation used in the third petmodel"},
    'jhtadd'       :{'value'  : 0.546 , 'comment': "!! parameter for second petmodel"},
    'jhtscale'     :{'value'  : 0.546 , 'comment': "!! parameter for second petmodel"}
}

# for 19 land covers
Snow_land_submodel_1_cec = {
    'section_head' :"""!!	-----																
!!	SNOWMELT Landuse dependent parameters													
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19
!! snowmelt submodel:2, snow heat submodel: 1""",
    'ttmp'         :{'value'  : [-9.7740,-2.4419,2.5478,-3.8724,2.9143,-7.2759,-6.1012,-6.5266,\
                                 -1.8872,-1.2143,-9.9603,-5.4364,-9.774,-9.774,-9.7740,-9.7740,\
                                 -9.7740,-9.7740,-9.7740],
                     'comment': "!! threshold temperature for snow melt snow density and evapotranspiration"},
    'cmlt'         :{'value'  : [9.7021,6.0035,1.1786,9.3525,1.7176,5.8523,4.1957,8.6383,\
                                 8.0090,5.4865,1.1010,5.5150,9.7021,9.7021,9.7021,9.7021,\
                                 9.7021,9.7021,9.7021],
                     'comment': "!! melting parameter for snow"},
    'cmrad'        :{'value'  : [0.249065876,0.249065876,1.5,0.176534861,0.685361445,0.174564317,\
                                 0.174564317,0.174564317,0.174564317,0.685361445,0.501842737,\
                                 0.011482887,0.249065876,0.249065876,0.249065876,0.249065876,\
                                 0.249065876,0.249065876,0.249065876],
                     'comment': "!! coefficient for radiation snow melt, parameter for second snowmelt model"},
    'snalbmin'     :{'value'  : [0.524781764,0.524781764,0.45,0.250044137,0.243925437,0.251664609,\
                                0.251664609,0.251664609,0.251664609,0.243925437,0.409460604,0.22856541,\
                                0.524781764,0.524781764,0.524781764,0.524781764,0.524781764,\
                                0.524781764,0.524781764],
                     'comment': "!! parameter for second snowmelt model"},
    'snalbmax'     :{'value'  : [0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,\
                                0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,\
                                0.85,0.85,0.85],
                     'comment': "!! parameter for second snowmelt model"},
    'snalbkexp'    :{'value'  : [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,\
                                0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1],
                     'comment': "!! parameter for second snowmelt model"},
    'snkika'       :{'value'  : [50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50],
                     'comment': "!! snow heat model, relation between snow thermal conductivity and surface heat exchange coefficient"}
}

# for 19 land covers
Snow_land_submodel_2_cec = {
    'section_head' :"""!!	-----																
!!	SNOWCOVER parameters (general and landuse) - baseline from Rossby RCA model (Samuelsson et al 2006;Lindström et al)																
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19
!! used in SNOWMELT submodel:2 """,
    'fscdistmax'   :{'value'  : [0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
                                 0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8],
                     'comment': "!! maximum snow distribution factor"},
    'fscdist0'     :{'value'  : [0.571998656,0.571998656,0.6,0.672227979,0.718591213,\
                                 0.672161579,0.672161579,0.672161579,0.672161579,\
                                 0.718591213,0.302164137,0.663832068,0.663832068,\
                                 0.663832068,0.663832068,0.663832068,0.663832068,\
                                 0.663832068,0.663832068],
                     'comment': "!! minimum snow distribution factor"},
    'fscdist1'     :{'value'  : [0.001,0.001,0.001,0.001,0.001,0.001,0.001,\
                                 0.001,0.001,0.001,0.001,0.001,0.001,0.001,\
                                 0.001,0.001,0.001,0.001,0.001],
                     'comment': "!! std coefficient for snow distribution factor parameter for second snowmelt model"},
    'fscmax'       :{'value'  : 1.00     , 'comment': "!! maximum fractional snow cover area"},
    'fscmin'       :{'value'  : 0.01     , 'comment': "!! minimum fractional snow cover area"},
    'fsclim'       :{'value'  : 0.001    , 'comment': "!! limit of fractional snow cover area for onset of snowmax"},
    'fsck1'        :{'value'  : 0.2      , 'comment': "!! Snowmass threshold to start decreasing the internal snowmax variable towards the end of the melt season"},
    'fsckexp'      :{'value'  : 0.000001 , 'comment': "!! Coefficient in exponential decrease of the internal Snowmax variable"}
}


# evapoation land use
Evap_land_cec = {
    'section_head' :"""!!	-----																
!!																	
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19""",
    'kc3'          :{'value'  :[1.017511845,1.017511845,1.201224208,1.334493399,1.265059352,\
                                1.020708799,1.020708799,1.020708799,1.020708799,1.265059352,\
                                1.342448354,1.024959087,1.024959087,1.024959087,1.024959087,\
                                1.024959087,1.024959087,1.024959087,1.024959087],
                     'comment': "!! crop coefficient for third petmodel"},
    'alb'          :{'value'  :[0.476534277,0.476534277,0.7,0.45542863,0.669192433,0.799822092,\
                                0.799822092,0.799822092,0.799822092,0.669192433,0.400103867,\
                                0.479658425,0.479658425,0.479658425,0.479658425,0.479658425,\
                                0.479658425,0.479658425,0.479658425],
                     'comment': "!! albedo for petmodels"},
    'ttrig'        :{'value'  :[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     'comment': "!! temperature threshold for soil temperature control on soil evapotranspiration"},
    'treda'        :{'value'  :[0.84,0.84,0.84,0.84,0.95,0.95,0.95,\
                                0.95,0.95,0.7,0.9,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8],
                     'comment': "!! soil temperature control on soil evapotranspiration"},
    'tredb'        :{'value'  :[0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
                                0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4],
                     'comment': "!! soil temperature control on soil evapotranspiration"},
    'cevp'         :{'value'  :[0.22,0.22,1.6,1.9,0.17,0.17,0.17,0.17,\
                                0.17,0.1,0.21,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07],
                     'comment': "!! evapotranspiration parameter"},
    'fepotsnow'    :{'value'  :[0.912879467,0.912879467,0.18,0.533387661,\
                                0.460848987,0.12002416,0.12002416,0.12002416,\
                                0.12002416,0.460848987,0.206956849,0.197802201,\
                                0.197802201,0.197802201,0.197802201,0.197802201,\
                                0.197802201,0.197802201,0.197802201],
                     'comment': "!! fraction of snow-free potential evapotranspiration, used for calculation of snow evaporation"}
}

# land use recession coeffcient
Land_recession_cec = {
    'section_head' :"""!! ------------
!! Frozen soil infiltration parameters															
!!Land use parameters																	
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19""",
    'srrcs'        :{'value': [0.1259,0.0701,0.187,0.1977,0.0951,0.1208,0.1594,0.0694,0.1136,0.0575,1,0.1213,0.1213,0.1213,0.1213,0.1213,0.1213,0.1213,0.1213],
                     'comment': "!! (landuse) recession coefficient for surface runoff should be set to one for lake and riverclasses with floodplains"}
}

# Frozen soil infiltration for land use
Forzen_soil_infil_LU_cec = {
    'section_head' :"""!! ------------
!! Frozen soil infiltration parameters															
!!Land use parameters																	
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19""",
    'surfmem'      :{'value': [17.8,17.8,17.8,17.8,5.15,5.15,5.15,\
                               5.15,5.15,5.15,5.15,5.15,5.15,5.15,\
                               5.15,5.15,5.15,5.15,5.15],
                     'comment': "!! upper soil layer soil temperature memory"},
    'depthrel'     :{'value': [1.1152,1.1152,1.1152,1.1152,2.47,2.47,\
                               2.47,2.47,2.47,2.47,2.47,2.47,2.47,2.47,\
                               2.47,2.47,2.47,2.47,2.47],
                     'comment': "!! depth relation for soil temperature memory"},
    'frost'        :{'value': [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
                     'comment': "!! frost depth parameter"}
}


# Frozen soil infiltration for soil
Forzen_soil_infil_soil_USDA = {
    'section_head' :"""!!======================================================																	
!! Frozen soil infiltration parameters																	
!! General and specific frozen soil 
!! SOIL:	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10	S11	S12					
!! for frozen soil submodel: 2""",
    'deepmem'      :{'value': 1000,  'comment': "!! deep soil temperature memory"},
    'bfroznsoil'   :{'value': [2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1],
                     'comment': "!! ??"},
    'logsatmp'     :{'value': [1.15,1.15,1.88,1.59,1.88,2.17,2.46,2.46,2.46,2.46,2.46,2.46],
                     'comment': "!! coefficient in unfrozen soil water content function"},
    'bcosby'       :{'value': [4.74,4.74,5.33,7.22,5.33,3.44,1.55,1.55,1.55,1.55,1.55,1.55],
                     'comment': "!! coefficient in unfrozen soil water content function"}
}


# soil class parameters
Soil_class_param_USDA = {
    'section_head' :"""!!	=======================================================================================================																
!!	"SOIL/LAND HYDRAULIC RESPONSE PARAMETERS - recession coef., water retention, infiltration, macropore, surface runoff; etc."																
!!	-----																
!!	Soil-class parameters																
!!	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10	S11	S12					
!! surfacerunoff submodel: 0, soilleakage submodel: 0
!! recession coefficient for surface runoff should be set to one for lake and riverclasses with floodplains""",
    'rrcs1'        :{'value': [0.3201,0.2982,0.2663,0.451,0.1637,0.1637,0.1637,0.1637,0.1637,0.1637,0.1637,0.1637],
                     'comment': "!! recession coefficient for uppermost soil layer"},
    'rrcs2'        :{'value': [0.1612,0.0858,0.1422,0.0112,0.1914,0.1914,0.1914,0.1914,0.1914,0.1914,0.1914,0.1914],
                     'comment': "!! recession coefficient for lowest soil layer"},
    'trrcs'        :{'value': [0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15],
                     'comment': "!! recession coefficient for tile drains"},
    'mperc1'       :{'value': [63.9842,119.5863,93.9854,111.8318,20.1177,20.1177,20.1177,20.1177,20.1177,20.1177,20.1177,20.1177],
                     'comment': "!! maximum percolation capacity from soil layer one to soil layer two"},
    'mperc2'       :{'value': [97.6492,12.5429,20.0276,79.481,12.0754,12.0754,12.0754,12.0754,12.0754,12.0754,12.0754,12.0754],
                     'comment': "!! maximum percolation capacity from soil layer two to soil layer three"},
    'sfrost'       :{'value': [1,1,1,1,1,1,1,1,1,1,1,1],
                     'comment': "!! frost depth parameter"},
    'srrate'       :{'value': [0.4975,0.4489,0.1874,0.4499,0.4956,0.4956,0.4956,0.4956,0.4956,0.4956,0.4956,0.4956],
                     'comment': "!! fraction for surface runoff"},
    'wcwp1'        :{'value': [0.2732,0.214,0.1479,0.4233,0.4941,0.4941,0.4941,0.4941,0.4941,0.4941,0.4941,0.4941],
                     'comment': "!! wilting point as a fraction for uppermost soil layer"},
    'wcwp2'        :{'value': [0.2293,0.256,0.0984,0.4303,0.1308,0.1308,0.1308,0.1308,0.1308,0.1308,0.1308,0.1308],
                     'comment': "!! wilting point as a fraction for second soil layer"},
    'wcwp3'        :{'value': [0.3257,0.058,0.3639,0.0433,0.3632,0.3632,0.3632,0.3632,0.3632,0.3632,0.3632,0.3632],
                     'comment': "!! wilting point as a fraction for lowest soil layer"},
    'wcfc1'        :{'value': [0.4344,0.1758,0.3526,0.4812,0.4894,0.4894,0.4894,0.4894,0.4894,0.4894,0.4894,0.4894],
                     'comment': "!! fraction of soil available for evapotranspiration but not for runoff for uppermost soil layer"},
    'wcfc2'        :{'value': [0.1392,0.1966,0.3818,0.1163,0.2385,0.2385,0.2385,0.2385,0.2385,0.2385,0.2385,0.2385],
                     'comment': "!! fraction of soil available for evapotranspiration but not for runoff for second soil layer"},
    'wcfc3'        :{'value': [0.2307,0.2075,0.4055,0.077,0.343,0.343,0.343,0.343,0.343,0.343,0.343,0.343],
                     'comment': "!! fraction of soil available for evapotranspiration but not for runoff for lowest soil layer"},
    'wcep1'        :{'value': [0.8729,0.4168,0.8743,0.5142,0.0117,0.0117,0.0117,0.0117,0.0117,0.0117,0.0117,0.0117],
                     'comment': "!! effective porosity as a fraction for uppermost soil layer"},
    'wcep2'        :{'value': [0.1177,0.2773,0.0329,0.8547,0.087,0.087,0.087,0.087,0.087,0.087,0.087,0.087],
                     'comment': "!! effective porosity as a fraction for second soil layer"},
    'wcep3'        :{'value': [0.3064,0.8004,0.5832,0.474,0.8299,0.8299,0.8299,0.8299,0.8299,0.8299,0.8299,0.8299],
                     'comment': "!! effective porosity as a fraction for lowest soil layer"},
    'rrcs3'        :{'value': 0.1612,
                     'comment': "!! recession coefficient for slope dependence (upper layer) ???"}
}




# general info that is shared between the soil and land cover
General = {'header':header,
           'General_meteo_param': General_meteo_param,
           'Snow_param': Snow_param,
           'Deep_groundwater': Deep_groundwater,
           'River_routing': River_routing,
           'Glacier': Glacier,
           'Glacier_melt': Glacier_melt, 
           'Evap': Evap}

# USDA soil type from triangle
# USDA 12 soil types
USDA = {'Forzen_soil_infil_soil': Forzen_soil_infil_soil_USDA,
        'Soil_class_param': Soil_class_param_USDA}

# land use related parameter
# CEC 19 land covers
CEC = {'Snow_land_submodel_1': Snow_land_submodel_1_cec,
       'Snow_land_submodel_2': Snow_land_submodel_2_cec,
       'Evap_land': Evap_land_cec,
       'Land_recession': Land_recession_cec,
       'Forzen_soil_infil_LU': Forzen_soil_infil_LU_cec}