import sys

macros = (
    ('"""DEFAULT_ISTATE"""', "0"),
    ('"""DEFAULT_RSTATE"""', "1.0,0.0,0.0,0.0,0.0,0.0"),
    ('"""RECIPROCAL_OF_N"""', "0"),
    ('"""INDEX_TO_AUTOCORR"""', "1"),
    ('"""INDEX_TO_UPPERBOUND"""', "2"),
    ('"""INDEX_TO_LOWERBOUND"""', "3"),
    ('"""INDEX_TO_N_EVAL"""', "4"),
    ('"""INDEX_TO_LOG_RATIO"""', "5"),
    ('"""INDEX_TO_ROTMAT"""', "6"),
    ('"""NOT_YET_RUN_FLAG"""', "-1"),
    ('"""IS_FINISHED_FLAG"""', "0"),
    ('"""ISNOT_FINISHED_FLAG"""', "1"),
)


def replace_line(l):
    ret = l
    for f, t in macros:
        ret = ret.replace(f, t)
    return ret


with open(sys.argv[1], "r") as f:
    for l in f.readlines():
        print(replace_line(l), end="")
