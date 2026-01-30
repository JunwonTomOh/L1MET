#include <cassert>
#include <cmath>
#include <cstdio>

#include "utils/rufl_io.h"
#include "firmware/puppimet.h"

#include "hls_math.h"

int main() {
    
    int matched_sin = 0;
    int unmatched_sin = 0;
    int matched_cos = 0;
    int unmatched_cos = 0;
    float total_diff_sin = 0;
    float total_diff_cos = 0;

    for(int i=-720; i<721; i++) {
        phi_t in_phi = phi_t(i);
        
        LUT_tri_T hls_sin, hls_cos;
        LUT_tri_T pred_sin, pred_cos;
        
        Get_LUT(in_phi, hls_sin, hls_cos);
        
        float pred_phi = floatPhi(phi_t(i));
        pred_cos = LUT_tri_T(cos(pred_phi));
        pred_sin = LUT_tri_T(sin(pred_phi));

        typedef ap_fixed<2+LUT_tri_T_Float*2, 2+LUT_tri_T_Float> scaled_tri_t;
        bool is_sin_Same, is_cos_Same;
        scaled_tri_t sin_diff, cos_diff;

        sin_diff = scaled_tri_t(hls_sin) - scaled_tri_t(pred_sin);
        cos_diff = scaled_tri_t(hls_cos) - scaled_tri_t(pred_cos);

        is_sin_Same = (sin_diff<<LUT_tri_T_Float) < 1;
        // is_sin_Same = hls_sin == pred_sin;
        is_cos_Same = (cos_diff<<LUT_tri_T_Float) < 1;
        // is_cos_Same = hls_cos == pred_cos;


        if (!is_sin_Same){
            unmatched_sin += 1;
            total_diff_sin += float(sin_diff);
            std::cout << "Sin value Unmatched!!\n"
                      << "Input Phi Value: " << i
                      << " / " << floatPhi(in_phi) << "\n"
                      << "HLS LUT Sin: " << hls_sin.to_string()
                      << " / " << hls_sin << "\n"
                      << "Predicted Sin: " << pred_sin.to_string()
                      << " / " << pred_sin
                      << std::endl;
        } else matched_sin += 1;

        if (!is_cos_Same){
            unmatched_cos += 1;
            total_diff_cos += float(cos_diff);
            std::cout << "Cos value Unmatched!!\n"
                      << "Input Phi Value: " << i
                      << " / " << floatPhi(in_phi) << "\n"
                      << "HLS LUT Cos: " << hls_cos.to_string()
                      << " / " << hls_cos << "\n"
                      << "Predicted Cos: " << pred_cos.to_string()
                      << " / " << pred_cos
                      << std::endl;
        } else matched_cos += 1;

    }
    std::cout << "Sin Matched: " << matched_sin
              << ", Unmatched: " << unmatched_sin << "\n"
              << "Mean Diff: " << total_diff_sin / unmatched_sin
              << std::endl;
    std::cout << "Cos Matched: " << matched_cos
              << ", Unmatched: " << unmatched_cos << "\n"
              << "Mean Diff: " << total_diff_cos / unmatched_cos
              << std::endl;
    return 0;
}