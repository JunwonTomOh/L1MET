#include "puppimet.h"
#include "lut.h"
#include <hls_math.h>
// #include <hls_stream.h>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#include <ap_int.h>
#include <ap_fixed.h>
#endif


void Get_LUT(phi_t in_phi, LUT_tri_T &out_sin, LUT_tri_T &out_cos){

    int sin_index = 0;
    int cos_index = 0;

    if (in_phi <= -360)      sin_index = in_phi + 720;
    else if (in_phi < 0)     sin_index = -in_phi;
    else if (in_phi <= 360)  sin_index = in_phi;
    else                  sin_index = 720 - in_phi;
    
    if (in_phi <= -360)      cos_index = in_phi + 720;
    else if (in_phi < 0)     cos_index = -in_phi;
    else if (in_phi <= 360)  cos_index = in_phi;
    else                  cos_index = 720 - in_phi;
    

    LUT_tri_T sin_mag = sin_LUT[sin_index];
    LUT_tri_T cos_mag = cos_LUT[cos_index];

    
    LUT_tri_T sin_cal = sin_mag;
    if (in_phi < 0) sin_cal = -sin_cal;

    LUT_tri_T cos_cal = cos_mag;
    if (in_phi < -360 || in_phi > 360) cos_cal = -cos_cal;

    out_sin = sin_cal;
    out_cos = cos_cal;
}

void Get_xy(Particle_T in_particles, Particle_xy &proj_xy) {
    // This function calculates x, y projection values using LUT
    #pragma HLS pipeline
    
    int sin_index = 0;
    int cos_index = 0;

    const int phi = (int)in_particles.hwPhi;


    if (phi <= -360)      sin_index = phi + 720;
    else if (phi < 0)     sin_index = -phi;
    else if (phi <= 360)  sin_index = phi;
    else                  sin_index = 720 - phi;
    
    if (phi <= -360)      cos_index = phi + 720;
    else if (phi < 0)     cos_index = -phi;
    else if (phi <= 360)  cos_index = phi;
    else                  cos_index = 720 - phi;

    LUT_tri_T sin_mag = sin_LUT[sin_index];
    LUT_tri_T cos_mag = cos_LUT[cos_index];

    
    LUT_tri_T sin_cal = sin_mag;
    if (phi < 0) sin_cal = -sin_cal;

    LUT_tri_T cos_cal = cos_mag;
    if (phi < -360 || phi > 360) cos_cal = -cos_cal;

    
    proj_xy.hwPx = in_particles.hwPt * cos_cal;
    proj_xy.hwPy = in_particles.hwPt * sin_cal;
}


void Sum_Particles(Particle_xy proj_xy[N_INPUT_LINKS], Particle_xy &met_xy) {
  // Sum of maximum 128 vectors in parallel

  #pragma HLS pipeline

  proj_t proj_x[N_INPUT_LINKS];
  proj_t proj_y[N_INPUT_LINKS];
  proj_t met_x = 0;
  proj_t met_y = 0;

  met_xy.hwPx = 0;
  met_xy.hwPy = 0;

  #pragma HLS ARRAY_PARTITION variable=proj_xy complete
  #pragma HLS ARRAY_PARTITION variable=proj_x complete
  #pragma HLS ARRAY_PARTITION variable=proj_y complete

  for(int i=0; i<N_INPUT_LINKS; ++i) {
    #pragma HLS unroll
    proj_x[i] = proj_xy[i].hwPx;
    proj_y[i] = proj_xy[i].hwPy;
  }

  CALC_LOOP:
  for(int i=0; i<N_INPUT_LINKS; ++i) {
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