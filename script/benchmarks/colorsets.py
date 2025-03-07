import numpy as np
import matplotlib.colors as colors


def hex_to_RGB(hex_str):
    """#FFFFFF -> [255,255,255]"""
    # Pass 16 to the integer function for change of base
    return [int(hex_str[i : i + 2], 16) for i in range(1, 6, 2)]


def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    assert n > 1
    c1_rgb = np.array(hex_to_RGB(c1)) / 255
    c2_rgb = np.array(hex_to_RGB(c2)) / 255
    mix_pcts = [x / (n - 1) for x in range(n)]
    rgb_colors = [((1 - mix) * c1_rgb + (mix * c2_rgb)) for mix in mix_pcts]
    return [
        "#" + "".join([format(int(round(val * 255)), "02x") for val in item])
        for item in rgb_colors
    ]


uwht = "#FFFFFF"
ugrb = "#C8C8C8"
ugrd = "#7F878F"
ublk = "#000000"

ured = "#FF2800"
ublu = "#0041FF"
ugrn = "#35A16B"
uylw = "#FAF500"
usbu = "#66CCFF"
upnk = "#FF99A0"
uorg = "#FF9900"
uppl = "#9A0079"
ubrw = "#663300"

sred = "#FFD1D1"
sblu = "#B4EBFA"
sgrn = "#87E7B0"
sygr = "#CBF266"
sylw = "#FFFF99"
sorg = "#EDC58F"
sppl = "#C7B2DE"

mycmap1 = colors.LinearSegmentedColormap.from_list(
    "mycmap1", [ublk, ublu, sblu, sylw, uwht]
)
