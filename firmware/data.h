#ifndef PUPPIMET_DATA_H
#define PUPPIMET_DATA_H

#include "ap_int.h"
#include "ap_fixed.h"
#include <hls_math.h>

#include "../utils/bit_encoding.h"

#define NPARTICLES 128
#define N_INPUT_LINKS 36
#define N_Frames 54

typedef ap_ufixed<14, 12, AP_TRN, AP_SAT> pt_t;
typedef ap_ufixed<28, 24, AP_TRN, AP_SAT> pt2_t;
typedef ap_fixed<16, 14, AP_TRN, AP_SAT> pxy_t;
typedef ap_int<12> eta_t;
typedef ap_int<11> phi_t;   


typedef ap_fixed<18, 14> acc_t; // for hls math
typedef ap_fixed<22, 12> proj_t; // for Projection values
typedef ap_fixed<32, 22> proj2_t; // for Projection values
typedef ap_fixed<16, 5> eta_var; // for eta values, set more bits for floats (float eta, eta edge, etc.)
typedef ap_fixed<8, 1> tri_t; // for trigonometric values, -1 ~ 1
// typedef ap_fixed<21, 3> poly_t;
typedef ap_fixed<32, 2> poly_t;
typedef ap_fixed<32, 2> poly2_t;

#define LUT_tri_T_Float 10
typedef ap_fixed<2+LUT_tri_T_Float, 2> LUT_tri_T;

typedef ap_uint<3> ParticleID;

constexpr int INTPHI_PI = 720;
constexpr float ETAPHI_LSB = M_PI / INTPHI_PI;
inline float floatPhi(phi_t phi) { return phi.to_float() * ETAPHI_LSB; }


struct Sum {
    pt_t hwPt;
    phi_t hwPhi;
    pt_t hwSumPt;

    inline bool operator==(const Sum &other) const {
      return hwPt == other.hwPt && hwPhi == other.hwPhi && hwSumPt == other.hwSumPt;
    }

    inline void clear() {
      hwPt = 0;
      hwPhi = 0;
      hwSumPt = 0;
    }
    
    static const int BITWIDTH = pt_t::width + phi_t::width + pt_t::width;
    inline ap_uint<BITWIDTH> pack() const {
      ap_uint<BITWIDTH> ret;
      unsigned int start = 0;
      pack_into_bits(ret, start, hwPt);
      pack_into_bits(ret, start, hwPhi);
      pack_into_bits(ret, start, hwSumPt);
      return ret;
    }
    inline static Sum unpack(const ap_uint<BITWIDTH> &src) {
      Sum ret;
      unsigned int start = 0;
      unpack_from_bits(src, start, ret.hwPt);
      unpack_from_bits(src, start, ret.hwPhi);
      unpack_from_bits(src, start, ret.hwSumPt);
      return ret;
    }
};

  struct PuppiObj {
    pt_t hwPt;
    eta_t hwEta;  // wider range to support global coordinates
    phi_t hwPhi;
    ParticleID hwId;

    static const int DATA_BITS_TOTAL = 21;
    ap_uint<DATA_BITS_TOTAL> hwData;

    inline bool operator==(const PuppiObj &other) const {
      return hwPt == other.hwPt && hwEta == other.hwEta && hwPhi == other.hwPhi && hwId == other.hwId &&
             hwData == other.hwData;
    }

    inline bool operator>(const PuppiObj &other) const { return hwPt > other.hwPt; }
    inline bool operator<(const PuppiObj &other) const { return hwPt < other.hwPt; }

    inline void clear() {
      hwPt = 0;
      hwEta = 0;
      hwPhi = 0;
      hwId = 0;
      hwData = 0;
    }


    static const int BITWIDTH = pt_t::width + eta_t::width + phi_t::width + 3 + DATA_BITS_TOTAL;
    inline ap_uint<BITWIDTH> pack() const {
      ap_uint<BITWIDTH> ret;
      unsigned int start = 0;
      pack_into_bits(ret, start, hwPt);
      pack_into_bits(ret, start, hwEta);
      pack_into_bits(ret, start, hwPhi);
      pack_into_bits(ret, start, hwId);
      pack_into_bits(ret, start, hwData);
      return ret;
    }
    inline void initFromBits(const ap_uint<BITWIDTH> &src) {
      unsigned int start = 0;
      unpack_from_bits(src, start, hwPt);
      unpack_from_bits(src, start, hwEta);
      unpack_from_bits(src, start, hwPhi);
      unpack_from_bits(src, start, hwId);
      unpack_from_bits(src, start, hwData);
    }
    inline static PuppiObj unpack(const ap_uint<BITWIDTH> &src) {
      PuppiObj ret;
      ret.initFromBits(src);
      return ret;
    }
  };


class Particle_T {
public:
// 25bits + 12bits
    pt_t hwPt; // 14bits
    eta_t hwEta; // 12bits
    phi_t hwPhi; // 11bits
    // bool operator >= (const Particle_T &b){
    //     return hwPt >= b.hwPt;
    // }
};

struct Particle_xy{
// 44bits
    proj_t hwPx; // 22bits
    proj_t hwPy; // 22bits
};

struct METCtrlToken{
    ap_uint<1> start_of_orbit;
    ap_uint<1> start;
    ap_uint<1> last;
    bool dataValid;
    bool frameValid;
};

#endif
