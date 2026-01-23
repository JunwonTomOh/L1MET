#include "puppimet.h"
#include <hls_math.h>
// #include <hls_stream.h>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

poly2_t cos2_par0[NPOINT] = {-1.00007,-0.924181,-0.707596,-0.382902,-0.000618262,0.382137,0.707056,0.923708,1.00007,0.924181,0.707594,0.383285,0.000188727,-0.382139,-0.706719,-0.923708};
poly2_t cos2_par1[NPOINT] = {9.164680268990924e-06, 0.0017064607695524156, 0.0031441321076514446, 0.004079929656016374, 0.004437063290882583, 0.004095969231842202, 0.0031107221424451436, 0.001689531075808071, -9.161756842493832e-06, -0.001706456406229286, -0.003143961938049376, -0.004103015998697129, -0.004411145151490469, -0.0040958165155326525, -0.0031310072316764474, -0.001689531075808071};
poly2_t cos2_par2[NPOINT] = {9.319674765430664e-06, 7.871694899063284e-06, 5.222989318251642e-06, 2.0256106486379287e-06, -1.9299417402361656e-06, -5.35167113952279e-06, -7.740062096537953e-06, -9.348822844786505e-06, -9.319674765430664e-06, -7.871694899063284e-06, -5.225331064666252e-06, -1.780776301343235e-06, 1.6556927733433181e-06, 5.3495197789955455e-06, 7.954684107366423e-06, 9.348822844786505e-06};

poly2_t sin2_par0[NPOINT] = {0.000524872,-0.382229,-0.706791,-0.923959,-1.00008,-0.924156,-0.707264,-0.383199,-0.000525527,0.382228,0.706792,0.923752,1.00013,0.924155,0.707535,0.3832};
poly2_t sin2_par1[NPOINT] = {-0.004431478237276202, -0.00409041472149773, -0.0031267268116859314, -0.00167440343451641, 9.741773386162849e-06, 0.0017049641497188307, 0.00312406082125351, 0.0040978672774037465, 0.004431478237276202, 0.00409041472149773, 0.0031266351819002015, 0.0016868781753450394, -1.249302315254411e-05, -0.001704846339994321, -0.003140405829698437, -0.0040978672774037465};
poly2_t sin2_par2[NPOINT] = {1.870674613498914e-06, 5.292404012785538e-06, 7.909829192302831e-06, 9.188746390688592e-06, 9.313525301268721e-06, 7.887020962996302e-06, 5.435897856093815e-06, 1.8358587462761668e-06, -1.870668901922293e-06, -5.292404012785538e-06, -7.908420336736317e-06, -9.320836119343602e-06, -9.284396260501616e-06, -7.88869635880513e-06, -5.262894200243701e-06, -1.835864457852788e-06};

phi_t phi2_edges[NPOINT+1] = {-720, -630, -540, -450, -360, -270, -180, -90, 0, 90, 180, 270, 360, 450, 540, 630, 720};


void Get_xy(Particle_T in_particles, Particle_xy &proj_xy) {
    // This function calculates x, y projection values using 2nd Polynomial interpolation
    
    #pragma HLS pipeline

    int phibin = 0;
    if      (in_particles.hwPhi < phi2_edges[1]) phibin = 0;
    else if (in_particles.hwPhi < phi2_edges[2]) phibin = 1;
    else if (in_particles.hwPhi < phi2_edges[3]) phibin = 2;
    else if (in_particles.hwPhi < phi2_edges[4]) phibin = 3;
    else if (in_particles.hwPhi < phi2_edges[5]) phibin = 4;
    else if (in_particles.hwPhi < phi2_edges[6]) phibin = 5;
    else if (in_particles.hwPhi < phi2_edges[7]) phibin = 6;
    else if (in_particles.hwPhi < phi2_edges[8]) phibin = 7;
    else if (in_particles.hwPhi < phi2_edges[9]) phibin = 8;
    else if (in_particles.hwPhi < phi2_edges[10]) phibin = 9;
    else if (in_particles.hwPhi < phi2_edges[11]) phibin = 10;
    else if (in_particles.hwPhi < phi2_edges[12]) phibin = 11;
    else if (in_particles.hwPhi < phi2_edges[13]) phibin = 12;
    else if (in_particles.hwPhi < phi2_edges[14]) phibin = 13;
    else if (in_particles.hwPhi < phi2_edges[15]) phibin = 14;
    else if (in_particles.hwPhi < phi2_edges[16]) phibin = 15;
    
    poly_t cos_var = cos2_par0[phibin] + cos2_par1[phibin] * (in_particles.hwPhi - phi2_edges[phibin]) + cos2_par2[phibin] * (in_particles.hwPhi - phi2_edges[phibin]) * (in_particles.hwPhi - phi2_edges[phibin]);
    poly_t sin_var = sin2_par0[phibin] + sin2_par1[phibin] * (in_particles.hwPhi - phi2_edges[phibin]) + sin2_par2[phibin] * (in_particles.hwPhi - phi2_edges[phibin]) * (in_particles.hwPhi - phi2_edges[phibin]);
    proj_xy.hwPx = in_particles.hwPt * cos_var;
    proj_xy.hwPy = in_particles.hwPt * sin_var;
}


void Sum_Particles(Particle_xy proj_xy[NPARTICLES], Particle_xy &met_xy) {
  // Sum of maximum 128 vectors in parallel

  #pragma HLS pipeline

  proj_t proj_x[NPARTICLES];
  proj_t proj_y[NPARTICLES];
  proj_t met_x = 0;
  proj_t met_y = 0;

  met_xy.hwPx = 0;
  met_xy.hwPy = 0;

  #pragma HLS ARRAY_PARTITION variable=proj_xy complete
  #pragma HLS ARRAY_PARTITION variable=proj_x complete
  #pragma HLS ARRAY_PARTITION variable=proj_y complete

  for(int i=0; i<NPARTICLES; ++i) {
    #pragma HLS unroll
    proj_x[i] = proj_xy[i].hwPx;
    proj_y[i] = proj_xy[i].hwPy;
  }

  CALC_LOOP:
  for(int i=0; i<NPARTICLES; ++i) {
    #pragma HLS unroll
    // Note: If Proj value has more float bits than metx&y, it makes lots of latency
    met_xy.hwPx -= proj_x[i];
    met_xy.hwPy -= proj_y[i];
  }
}


void puppimet_xy(Particle_T in_particles[N_INPUT_LINKS], Particle_xy &met_xy, METCtrlToken token_d, METCtrlToken& token_q) {
    #pragma HLS PIPELINE

    Particle_xy proj_xy[N_INPUT_LINKS];

    #pragma HLS ARRAY_PARTITION variable=in_particles complete
    #pragma HLS ARRAY_PARTITION variable=proj_xy complete

    #pragma HLS aggregate variable=token_d compact=bit
    #pragma HLS aggregate variable=token_q compact=bit
    #pragma HLS aggregate variable=in_particles compact=bit
    #pragma HLS aggregate variable=met_xy compact=bit
    
    #pragma HLS INTERFACE ap_none port=in_particles
    #pragma HLS INTERFACE ap_none port=met_xy
    #pragma HLS INTERFACE ap_none port=token_d
    #pragma HLS INTERFACE ap_none port=token_q

    static proj_t    global_acc_px;
    static proj_t    global_acc_py;
    METCtrlToken token_i;
    
    if (token_d.start == 1){
        // Reset
        global_acc_px      = 0;
        global_acc_py      = 0;
        
        token_i.start_of_orbit = 0;
        token_i.start = 1;
        token_i.last = 0;
        token_i.dataValid = false;
        token_i.frameValid = false;
    } else if (token_d.last) {
        token_i.start_of_orbit = token_d.start_of_orbit;
        token_i.start = 0;
        token_i.last = 1;
        token_i.dataValid = true;
        token_i.frameValid = true;
    } else {
        token_i.start_of_orbit = 0;
        token_i.start = 0;
        token_i.last = 0;
        token_i.dataValid = false;
        token_i.frameValid = false;
    }

    PROJ_LOOP:
    for (int i = 0; i < N_INPUT_LINKS; ++i) {
        #pragma HLS UNROLL
        Get_xy(in_particles[i], proj_xy[i]);
    }

    Particle_xy slice_sum;

    Sum_Particles(proj_xy, slice_sum);
    global_acc_px += slice_sum.hwPx;
    global_acc_py += slice_sum.hwPy;
    met_xy.hwPx = global_acc_px;
    met_xy.hwPy = global_acc_py;

    token_q = token_i;
  return;
}


void pxpy_to_ptphi(const Particle_xy met_xy, Sum &out_met, METCtrlToken token_d, METCtrlToken& token_q) {
  // convert x, y coordinate to pt, phi coordinate using HLS math library
  
  #pragma HLS pipeline
  
  #pragma HLS aggregate variable=met_xy compact=bit
  #pragma HLS aggregate variable=out_met compact=bit
  #pragma HLS aggregate variable=token_d compact=bit
  #pragma HLS aggregate variable=token_q compact=bit

  #pragma HLS interface ap_none port=met_xy
  #pragma HLS interface ap_none port=out_met
  #pragma HLS interface ap_none port=token_d
  #pragma HLS interface ap_none port=token_q

  out_met.clear();
  out_met.hwPt = hls::hypot(met_xy.hwPx, met_xy.hwPy);
  // Reduce Latency by not-using division.
  out_met.hwPhi = phi_t(ap_fixed<26, 11>(hls::atan2(met_xy.hwPy, met_xy.hwPx)) * ap_fixed<26, 11>(229.29936));
  // out_met.hwPhi = l1ct::Scales::makeGlbPhi(hls::atan2(met_xy.hwPy, met_xy.hwPx));
  

  token_q = token_d;

  return;
}