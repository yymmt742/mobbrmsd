import sys
import numpy as np

"""
gzmat templates are given by
$ echo 'CC' | obabel -i smi -o gzmat --gen3D > ethane.gzmat
$ echo 'CC(C)C' | obabel -i smi -o gzmat --gen3D > isobutane.gzmat
$ echo 'CC(C)(C)C' | obabel -i smi -o gzmat --gen3D > neopentane.gzmat
$ echo 'CC(C)(C)C(C)(C)C' | obabel -i smi -o gzmat --gen3D > tetramethylbutane.gzmat
$ echo 'CC(C)(C)CCC(C)(C)C' | obabel -i smi -o gzmat --gen3D > tetramethylhexane.gzmat
$ echo 'CC(C)(C)C(C)(C)C(C)(C)C' | obabel -i smi -o gzmat --gen3D > hexamethylpentane.gzmat
$ echo 'CC(C)(C)C(C)(C)C(C)(C)C(C)(C)C' | obabel -i smi -o gzmat --gen3D > octamethylhexane.gzmat
"""

sin = lambda x: np.sin(x * np.pi / 360.0)  # sin(x/2) for degree
cos = lambda x: np.cos(x * np.pi / 360.0)  # cos(x/2) for degree

rCH = 1.09
rCC = 1.54
rBB = 1.62  # bond length for tertiary carbon to tertiary carbon
aCCC = 111.66
aCCH = 111.09
aTTH = 180 * 2 * np.arctan(np.sqrt(2.0)) / np.pi

aWID = 135.0  # wide CCC angle
aNAR = 100.0  # narrow CCC angle

aPSI = 180.0 * np.arccos(-cos(aWID) * cos(aNAR)) / np.pi  # CCC angle to aWID
dPHI = (
    180
    * np.arccos(
        -sin(aWID)
        * cos(aNAR)
        / np.sqrt(
            (sin(aWID) * cos(aNAR)) ** 2
            + (cos(aWID) * sin(aNAR)) ** 2
            + (sin(aWID) * sin(aNAR)) ** 2
        )
    )
    / np.pi
)  # CCCC dihedral angle to aWID


def print_bonds(t, r, base, n):
    for i in range(base, base + n):
        print(f"{t}{i}={r:7.2f}")


def print_triplet_dihedrals(base, n):
    for i, ri in enumerate(360 * np.random.rand(n)):
        j = 3 * i + base
        print(f"d{j + 0}= {ri:6.2f}")
        print(f"d{j + 1}= 120.00")
        print(f"d{j + 2}= 240.00")


def mod360(v):
    if v < 0.0:
        return str(v + 360.0)
    elif v > 360.0:
        return str(v - 360.0)
    else:
        return str(v)


def print_ethane():
    print(
        """

#

 

0  1
C
C  1  r2
H  1  r3  2  a3
H  1  r4  2  a4  3  d4
H  1  r5  2  a5  3  d5
H  2  r6  1  a6  3  d6
H  2  r7  1  a7  6  d7
H  2  r8  1  a8  6  d8
Variables:"""
    )
    print_bonds("r", rCC, 2, 1)
    print_bonds("r", rCH, 3, 6)
    print_bonds("a", aCCH, 3, 6)
    print_triplet_dihedrals(3, 2)
    print()


def print_propane():
    print(
        """

#

 

0  1
C
C  1  r2
C  1  r3   2  a3
H  1  r4   2  a4   3  d4
H  1  r5   3  a5   2  d5
H  2  r6   1  a6   3  d6
H  2  r7   1  a7   6  d7
H  2  r8   1  a8   6  d8
H  3  r9   1  a9   2  d9
H  3  r10  1  a10  9  d10
H  3  r11  1  a11  9  d11
Variables:"""
    )
    print_bonds("r", rCC, 2, 2)
    print_bonds("r", rCH, 4, 9)

    print_bonds("a", aCCC, 3, 1)
    print_bonds("a", aCCH, 4, 9)

    print_bonds("d", 120.0, 4, 2)
    print_triplet_dihedrals(6, 2)


def print_isobutane():
    print(
        """

#

 

0  1
H
C  1  r2
C  2  r3   1  a3
C  2  r4   1  a4   3  d4
C  2  r5   1  a5   3  d5
H  3  r6   2  a6   1  d6
H  3  r7   2  a7   6  d7
H  3  r8   2  a8   6  d8
H  4  r9   2  a9   1  d9
H  4  r10  2  a10  9  d10
H  4  r11  2  a11  9  d11
H  5  r12  2  a12  1  d12
H  5  r13  2  a13 12  d13
H  5  r14  2  a14 12  d14
Variables:"""
    )
    print_bonds("r", rCC, 2, 4)
    print_bonds("r", rCH, 6, 9)

    print_bonds("a", aCCC, 3, 3)
    print_bonds("a", aCCH, 6, 9)

    print_bonds("d", 120.0, 4, 1)
    print_bonds("d", 240.0, 5, 1)
    print_triplet_dihedrals(6, 3)
    print()


def print_neopentane():
    print(
        """

#

 

0  1
C
C  1  r2
C  1  r3   2  a3
C  1  r4   2  a4   3  d4
C  1  r5   2  a5   3  d5
H  2  r6   1  a6   3  d6
H  2  r7   1  a7   6  d7
H  2  r8   1  a8   6  d8
H  3  r9   1  a9   2  d9
H  3  r10  1  a10  9  d10
H  3  r11  1  a11  9  d11
H  4  r12  1  a12  2  d12
H  4  r13  1  a13 12  d13
H  4  r14  1  a14 12  d14
H  5  r15  1  a15  2  d15
H  5  r16  1  a16 15  d16
H  5  r17  1  a17 15  d17
Variables:"""
    )
    print_bonds("r", rCC, 2, 5)
    print_bonds("r", rCH, 6, 12)

    print_bonds("a", aTTH, 3, 4)
    print_bonds("a", aCCH, 6, 12)

    print_bonds("d", 120.0, 4, 1)
    print_bonds("d", 240.0, 5, 1)
    print_triplet_dihedrals(6, 4)


def print_tetramethylbutane():
    print(
        """
#

 

0  1
C
C  1  r2
C  1  r3   2  a3
C  1  r4   2  a4   3  d4
C  1  r5   2  a5   3  d5
C  2  r6   1  a6   3  d6
C  2  r7   1  a7   6  d7
C  2  r8   1  a8   6  d8
H  3  r9   1  a9   2  d9
H  3  r10  1  a10  9  d10
H  3  r11  1  a11  9  d11
H  4  r12  1  a12  2  d12
H  4  r13  1  a13 12  d13
H  4  r14  1  a14 12  d14
H  5  r15  1  a15  2  d15
H  5  r16  1  a16 15  d16
H  5  r17  1  a17 15  d17
H  6  r18  2  a18  1  d18
H  6  r19  2  a19 18  d19
H  6  r20  2  a20 18  d20
H  7  r21  2  a21  1  d21
H  7  r22  2  a22 21  d22
H  7  r23  2  a23 21  d23
H  8  r24  2  a24  1  d24
H  8  r25  2  a25 24  d25
H  8  r26  2  a26 24  d26
Variables:"""
    )

    print_bonds("r", rBB, 2, 3)
    print_bonds("r", rCC, 5, 4)
    print_bonds("r", rCH, 9, 18)

    print_bonds("a", aCCC, 3, 6)
    print_bonds("a", aCCH, 9, 26)

    print_triplet_dihedrals(3, 8)
    print()


def print_tetramethylpentane():
    print(
        """
#

 

0  1
C
C  1   r2
C  1   r3  2  a3
C  2   r4  1  a4   3  d4
C  2   r5  1  a5   4  d5
C  2   r6  1  a6   4  d6
C  3   r7  1  a7   2  d7
C  3   r8  1  a8   7  d8
C  3   r9  1  a9   7  d9
H  1  r10  2  a10  3  d10
H  1  r11  3  a11  2  d11
H  4  r12  2  a12  1  d12
H  4  r13  2  a13 12  d13
H  4  r14  2  a14 12  d14
H  5  r15  2  a15  1  d15
H  5  r16  2  a16 15  d16
H  5  r17  2  a17 15  d17
H  6  r18  2  a18  1  d18
H  6  r19  2  a19 18  d19
H  6  r20  2  a20 18  d20
H  7  r21  3  a21  1  d21
H  7  r22  3  a22 21  d22
H  7  r23  3  a23 21  d23
H  8  r24  3  a24  1  d24
H  8  r25  3  a25 24  d25
H  8  r26  3  a26 24  d26
H  9  r27  3  a27  1  d27
H  9  r28  3  a28 27  d28
H  9  r29  3  a29 27  d29
Variables:"""
    )

    print_bonds("r", rBB, 2, 2)
    print_bonds("r", rCC, 4, 6)
    print_bonds("r", rCH, 10, 20)

    print_bonds("a", aWID, 3, 1)
    print_bonds("a", aCCC, 4, 6)
    print_bonds("a", aCCH, 10, 20)

    print_triplet_dihedrals(4, 2)
    print_bonds("d", 120.0, 10, 2)
    print_triplet_dihedrals(12, 6)
    print()


def print_tetramethylhexane():
    print(
        """
#

 

0  1
C
C  1   r2
C  1   r3   2  a3
C  2   r4   1  a4   3  d4
C  3   r5   1  a5   2  d5
C  3   r6   1  a6   5  d6
C  3   r7   1  a7   5  d7
C  4   r8   2  a8   1  d8
C  4   r9   2  a9   8  d9
C  4   r10  2  a10  8  d10
H  1   r11  2  a11  3  d11
H  1   r12  3  a12  2  d12
H  2   r13  1  a13  4  d13
H  2   r14  4  a14  1  d14
H  5   r15  3  a15  1  d15
H  5   r16  3  a16 15  d16
H  5   r17  3  a17 15  d17
H  6   r18  3  a18  1  d18
H  6   r19  3  a19 18  d19
H  6   r20  3  a20 18  d20
H  7   r21  3  a21  1  d21
H  7   r22  3  a22 21  d22
H  7   r23  3  a23 21  d23
H  8   r24  4  a24  2  d24
H  8   r25  4  a25 24  d25
H  8   r26  4  a26 24  d26
H  9   r27  4  a27  2  d27
H  9   r28  4  a28 27  d28
H  9   r29  4  a29 27  d29
H  10  r30  4  a30  2  d30
H  10  r31  4  a31 30  d31
H  10  r32  4  a32 30  d32
Variables:"""
    )

    print_bonds("r", rCC, 2, 9)
    print_bonds("r", rCH, 11, 22)

    print_bonds("a", aWID, 3, 2)
    print_bonds("a", aCCC, 5, 6)
    print_bonds("a", aCCH, 11, 22)

    print_bonds("d", 360 * np.random.rand(1)[0], 4, 1)
    print_triplet_dihedrals(5, 6)
    print_bonds("d", 120.0, 11, 4)
    print_triplet_dihedrals(15, 6)

    print()


def print_hexamethylpentane():
    print(
        """
#

 

0  1
C
C  1   r2
C  1   r3   2  a3
C  1   r4   2  a4   3  d4
C  1   r5   3  a5   2  d5
C  2   r6   1  a6   3  d6
C  2   r7   1  a7   6  d7
C  2   r8   1  a8   6  d8
C  3   r9   1  a9   2  d9
C  3   r10  1  a10  9  d10
C  3   r11  1  a11  9  d11
H  4   r12  1  a12  5  d12
H  4   r13  1  a13 12  d13
H  4   r14  1  a14 12  d14
H  5   r15  1  a15  4  d15
H  5   r16  1  a16 15  d16
H  5   r17  1  a17 15  d17
H  6   r18  2  a18  1  d18
H  6   r19  2  a19 18  d19
H  6   r20  2  a20 18  d20
H  7   r21  2  a21  1  d21
H  7   r22  2  a22 21  d22
H  7   r23  2  a23 21  d23
H  8   r24  2  a24  1  d24
H  8   r25  2  a25 24  d25
H  8   r26  2  a26 24  d26
H  9   r27  3  a27  1  d27
H  9   r28  3  a28 27  d28
H  9   r29  3  a29 27  d29
H  10  r30  3  a30  1  d30
H  10  r31  3  a31 30  d31
H  10  r32  3  a32 30  d32
H  11  r33  3  a33  1  d33
H  11  r34  3  a34 33  d34
H  11  r35  3  a35 33  d35
Variables:"""
    )

    print_bonds("r", rBB, 2, 2)
    print_bonds("r", rCC, 4, 8)
    print_bonds("r", rCH, 12, 24)

    print_bonds("a", aWID, 3, 1)
    print_bonds("a", aPSI, 4, 2)
    print_bonds("a", aCCC, 6, 6)
    print_bonds("a", aCCH, 12, 24)

    print_bonds("d", dPHI, 4, 2)
    print_triplet_dihedrals(6, 10)

    print()


def print_octamethylhexane():
    print(
        """
#

 

0  1
C
C  1   r2
C  1   r3   2   a3
C  2   r4   1   a4   3  d4
C  1   r5   2   a5   3  d5
C  1   r6   3   a6   2  d6
C  2   r7   1   a7   4  d7
C  2   r8   4   a8   1  d8
C  3   r9   1   a9   2  d9
C  3   r10  1   a10  9  d10
C  3   r11  1   a11  9  d11
C  4   r12  2   a12  1  d12
C  4   r13  2   a13 12  d13
C  4   r14  2   a14 12  d14
H  5   r15  1   a15  2  d15
H  5   r16  1   a16 15  d16
H  5   r17  1   a17 15  d17
H  6   r18  1   a18  2  d18
H  6   r19  1   a19 18  d19
H  6   r20  1   a20 18  d20
H  7   r21  2   a21  1  d21
H  7   r22  2   a22 21  d22
H  7   r23  2   a23 21  d23
H  8   r24  2   a24  1  d24
H  8   r25  2   a25 24  d25
H  8   r26  2   a26 24  d26
H  9   r27  3   a27  1  d27
H  9   r28  3   a28 27  d28
H  9   r29  3   a29 27  d29
H  10  r30  3   a30  1  d30
H  10  r31  3   a31 30  d31
H  10  r32  3   a32 30  d32
H  11  r33  3   a33  1  d33
H  11  r34  3   a34 33  d34
H  11  r35  3   a35 33  d35
H  12  r36  4   a36  2  d36
H  12  r37  4   a37 36  d37
H  12  r38  4   a38 36  d38
H  13  r39  4   a39  2  d39
H  13  r40  4   a40 39  d40
H  13  r41  4   a41 39  d41
H  14  r42  4   a42  2  d42
H  14  r43  4   a43 42  d43
H  14  r44  4   a44 42  d44
Variables:"""
    )

    print_bonds("r", rBB, 2, 3)
    print_bonds("r", rCC, 5, 10)
    print_bonds("r", rCH, 15, 30)

    print_bonds("a", aWID, 3, 2)
    print_bonds("a", aPSI, 5, 4)
    print_bonds("a", aCCC, 9, 6)
    print_bonds("a", aCCH, 15, 30)

    print_bonds("d", 360 * np.random.rand(1)[0], 4, 1)
    print_bonds("d", dPHI, 5, 4)
    print_triplet_dihedrals(9, 12)

    print()


if sys.argv[1] == "ethane":
    gen = print_ethane
elif sys.argv[1] == "propane":
    gen = print_propane
elif sys.argv[1] == "isobutane":
    gen = print_isobutane
elif sys.argv[1] == "neopentane":
    gen = print_neopentane
elif sys.argv[1] == "tetramethylbutane":
    gen = print_tetramethylbutane
elif sys.argv[1] == "tetramethylpentane":
    gen = print_tetramethylpentane
elif sys.argv[1] == "tetramethylhexane":
    gen = print_tetramethylhexane
elif sys.argv[1] == "hexamethylpentane":
    gen = print_hexamethylpentane
elif sys.argv[1] == "octamethylhexane":
    gen = print_octamethylhexane
else:
    exit()


for n in range(int(sys.argv[2]) if len(sys.argv) > 2 else 1):
    gen()
