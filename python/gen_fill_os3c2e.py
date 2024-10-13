#!/usr/bin/env python3


import argparse
import sys


def UniqueTriplet(quartet):
    if quartet[0] < quartet[1]:
        return False
    return True


# some helpers
amchar = "spdfghiklmnoqrtuvwxyzabceSPDFGHIKLMNOQRTUVWXYZABCE0123456789"

parser = argparse.ArgumentParser()
parser.add_argument("-l", type=int, required=True, help="Maximum AM")
parser.add_argument("-d", type=int, default=0, help="Derivative level")
parser.add_argument("outfile", type=str, help="Output file")

args = parser.parse_args()


maxam = args.l
der = args.d

print("-------------------------------")
print("Generating Array Filling")
print(f"Maximum AM: {maxam}")
print("-------------------------------")

valid = {}
for L in range(0, maxam+1):
    valid[L] = set()


for i in range(0, maxam + 1):
    for j in range(0, maxam + 1):
        for k in range(0, maxam + 1):
            q = (i, j, k)
            maxL = max(q)
            valid[maxL].add(q)

# Start the source file
sourcefile = args.outfile
with open(sourcefile, 'w', encoding='ascii') as f:
    f.write("/*\n")
    f.write(" Generated with:\n")
    f.write("   " + " ".join(sys.argv[:]))
    f.write("\n")
    f.write("*/\n\n")
    f.write("\n\n")
    f.write("#include \"simint/os3c2e/os3c2e.h\"\n")
    f.write("#include \"simint/os3c2e/os3c2e_init.h\"\n")
    f.write("#include \"simint/os3c2e/os3c2e_config.h\"\n")
    if der > 0:
        f.write(
            f"#include \"simint/os3c2e/gen/os3c2e_deriv{der}_generated.h\"\n")
    else:
        f.write("#include \"simint/os3c2e/gen/os3c2e_generated.h\"\n")
    f.write("\n\n")

    f.write("// Stores pointers to the os3c2e functions\n")
    f.write("#define AMSIZE   SIMINT_OS3C2E_MAXAM+1\n")
    f.write("#define DERSIZE  SIMINT_OS3C2E_DERIV+1\n")
    f.write(
        "extern simint_os3c2efunc simint_os3c2efunc_array[DERSIZE][AMSIZE][AMSIZE][AMSIZE][AMSIZE];\n")
    f.write("\n\n\n")

    # write the finalize functions
    if der > 0:
        f.write(f"void simint_os3c2e_deriv{der}_finalize(void)\n")
    else:
        f.write("void simint_os3c2e_finalize(void)\n")
    f.write("{\n")
    f.write("    // nothing to do\n")
    f.write("}\n\n\n")

    # Now write the filling function
    if der > 0:
        f.write(f"void simint_os3c2e_deriv{der}_init(void)\n")
    else:
        f.write("void simint_os3c2e_init(void)\n")
    f.write("{\n")

    for L, qset in valid.items():
        q1, q2, q3 = q
        l1, l2, l3 = [amchar[z] for z in q]
        if der > 0:
            f.write(f"    #if SIMINT_OS3C2E_DERIV{der}_MAXAM >= {L}\n")
        else:
            f.write(f"    #if SIMINT_OS3C2E_MAXAM >= {L}\n")

        for q in sorted(list(qset)):
            if der > 0:
                fname = f"os3c2e_deriv{der}_{l1}_{l2}_{l3}"
            else:
                fname = f"os3c2e_{l1}_{l2}_{l3}".format(
                    amchar[q[0]], amchar[q[1]], amchar[q[2]])

            f.write(f"    simint_os3c2efunc_array[{der}][{q1}][{q2}][{q3}] = {fname};\n")

        f.write("    #endif\n\n")

    f.write("}\n\n")
