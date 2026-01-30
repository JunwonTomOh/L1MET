#ifndef PUPPIMET_H
#define PUPPIMET_H

#include "data.h"
#define NPOINT 16
#define POLY_MAP 0.004363323129985824 // pi / 720

void puppimet_xy(Particle_T in_particles[N_INPUT_LINKS], Particle_xy &met_xy, METCtrlToken token_d, METCtrlToken& token_q);
void pxpy_to_ptphi(Particle_xy met_xy, Sum &out_met, METCtrlToken token_d, METCtrlToken& token_q);
void Get_LUT(phi_t in_phi, LUT_tri_T &out_sin, LUT_tri_T &out_cos);


#endif
