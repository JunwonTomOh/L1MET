from math import sin, cos, pi


with open("./firmware/lut.h", "w") as f:
    f.write('#include "puppimet.h"\n')
    f.write('#include <ap_int.h>\n')
    f.write('#include <ap_fixed.h>\n\n')

    f.write('static const LUT_tri_T sin_LUT[361] = {\n')
    for i in range(361):
        f.write(f'{sin(i*pi/720)}, \n')
    f.write('};\n')

    f.write('static const LUT_tri_T cos_LUT[361] = {\n')
    for i in range(361):
        f.write(f'{cos(i*pi/720):.15f}, \n')
    f.write('};\n')