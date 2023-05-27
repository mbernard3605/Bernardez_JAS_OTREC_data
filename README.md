# Bernardez_JAS_OTREC_data

MAT files:
anglecolormap5.mat - colormap for top-heaviness angle
ERA5_EOFS.mat - vertical motion EOFs calculated from ERA5 climatology
OTREC_B1.mat - More advanced thermodynamic variables from all B1 flights
OTREC_B1_zbp_start_ent.mat - simple thermodynamic variables from all B1 flights
OTREC_B2.mat - Same for B2 flights
OTREC_B2_zbp_start_ent.mat - Same for B2 flights
pres.mat - Mean pressure profile from all 3dvar data

.m files
Buoyancy_analysis_b1.m  - Creates OTREC_B1.mat
Buoyancy_analysis_b2.m  - Creates OTREC_B2.mat
calc_ent.m - Calculates entropy from temperature and pressure
calc_grad.m - Calculates gradient of 2D field
calc_ii.m - Calculates the instability index from profile of entropy and height
calc_qvstar.m - Calculates the saturation specific humidity 
calc_sf.m - Calculates the saturation fraction
calc_sq.m - Calculates the MDC
figure_2.m - Creates the base of figure 2
figure_3.m - Creates the base of figure 3
figure_4.m - Creates the base of figure 4
figure_5.m - Creates the base of figure 5
figure_6.m - Creates the base of figure 6
figure_8.m - Creates the base of figure 7
figure_10.m - Creates the base of figure 10
redblue.m - creates a red-blue colormap
Zero_Buoyancy_Plume_analysis_b1.m - Creates OTREC_B1_zbp_start_ent.mat
Zero_Buoyancy_Plume_analysis_b2.m - Creates OTREC_B2_zbp_start_ent.mat
zero_plume_buoyancy_multi_plume.m

