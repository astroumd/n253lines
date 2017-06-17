
#
# See also Table 2 in  http://adsabs.harvard.edu/abs/2015ApJ...801...63M
#

def pos253(i):
    """ return key positions in N253 (1..10) from Meier's Table 2:
            0 = map center, the reference position of N253
        1..10 = reference positions in the table
    """
    pos = [ ['00h47m33.100s',   '-25d17m17.50s' ],  # map reference 
            ['00h47m33.041s',	'-25d17m26.61s'	],  # pos 1
            ['00h47m32.290s',	'-25d17m19.10s'	],  #     2
            ['00h47m31.936s',	'-25d17m29.10s'	],  #     3
            ['00h47m32.792s',	'-25d17m21.10s'	],  #     4
            ['00h47m32.969s',	'-25d17m19.50s'	],  #     5
            ['00h47m33.159s',	'-25d17m17.41s'	],  #     6
            ['00h47m33.323s',	'-25d17m15.50s'	],  #     7
            ['00h47m33.647s',	'-25d17m13.10s'	],  #     8
            ['00h47m33.942s',	'-25d17m11.10s'	],  #     9
            ['00h47m34.148s',	'-25d17m12.30s'	],  # pos 10
          ]
    if i<0: return pos       # special case
    return pos[i]
