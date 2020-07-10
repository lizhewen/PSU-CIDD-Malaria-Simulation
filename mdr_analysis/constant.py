HEADLINE = "time_elapsed\tsystem_time\tyear\tmonth\tday\tseasonal_factor\ttreatment_coverage_0_1\ttreatment_coverage_0_10\tpopulation\tsep\teir\tsep\tbsp_2_10\tbsp_0_5\tblood_slide_prevalence\tsep\tmonthly_new_infection\tsep\tmonthly_new_treatment\tsep\tmonthly_clinical_episode\tsep\tmonthly_ntf_raw\tsep\tKNY--C1x\tKNY--C1X\tKNY--C2x\tKNY--C2X\tKNY--Y1x\tKNY--Y1X\tKNY--Y2x\tKNY--Y2X\tKYY--C1x\tKYY--C1X\tKYY--C2x\tKYY--C2X\tKYY--Y1x\tKYY--Y1X\tKYY--Y2x\tKYY--Y2X\tKNF--C1x\tKNF--C1X\tKNF--C2x\tKNF--C2X\tKNF--Y1x\tKNF--Y1X\tKNF--Y2x\tKNF--Y2X\tKYF--C1x\tKYF--C1X\tKYF--C2x\tKYF--C2X\tKYF--Y1x\tKYF--Y1X\tKYF--Y2x\tKYF--Y2X\tKNYNYC1x\tKNYNYC1X\tKNYNYC2x\tKNYNYC2X\tKNYNYY1x\tKNYNYY1X\tKNYNYY2x\tKNYNYY2X\tKYYYYC1x\tKYYYYC1X\tKYYYYC2x\tKYYYYC2X\tKYYYYY1x\tKYYYYY1X\tKYYYYY2x\tKYYYYY2X\tKNFNFC1x\tKNFNFC1X\tKNFNFC2x\tKNFNFC2X\tKNFNFY1x\tKNFNFY1X\tKNFNFY2x\tKNFNFY2X\tKYFYFC1x\tKYFYFC1X\tKYFYFC2x\tKYFYFC2X\tKYFYFY1x\tKYFYFY1X\tKYFYFY2x\tKYFYFY2X\tTNY--C1x\tTNY--C1X\tTNY--C2x\tTNY--C2X\tTNY--Y1x\tTNY--Y1X\tTNY--Y2x\tTNY--Y2X\tTYY--C1x\tTYY--C1X\tTYY--C2x\tTYY--C2X\tTYY--Y1x\tTYY--Y1X\tTYY--Y2x\tTYY--Y2X\tTNF--C1x\tTNF--C1X\tTNF--C2x\tTNF--C2X\tTNF--Y1x\tTNF--Y1X\tTNF--Y2x\tTNF--Y2X\tTYF--C1x\tTYF--C1X\tTYF--C2x\tTYF--C2X\tTYF--Y1x\tTYF--Y1X\tTYF--Y2x\tTYF--Y2X\tTNYNYC1x\tTNYNYC1X\tTNYNYC2x\tTNYNYC2X\tTNYNYY1x\tTNYNYY1X\tTNYNYY2x\tTNYNYY2X\tTYYYYC1x\tTYYYYC1X\tTYYYYC2x\tTYYYYC2X\tTYYYYY1x\tTYYYYY1X\tTYYYYY2x\tTYYYYY2X\tTNFNFC1x\tTNFNFC1X\tTNFNFC2x\tTNFNFC2X\tTNFNFY1x\tTNFNFY1X\tTNFNFY2x\tTNFNFY2X\tTYFYFC1x\tTYFYFC1X\tTYFYFC2x\tTYFYFC2X\tTYFYFY1x\tTYFYFY1X\tTYFYFY2x\tTYFYFY2X\t"

ENCODINGDB = [
    'KNY--C1x', 'KNY--C1X', 'KNY--C2x', 'KNY--C2X', 'KNY--Y1x', 'KNY--Y1X',
    'KNY--Y2x', 'KNY--Y2X', 'KYY--C1x', 'KYY--C1X', 'KYY--C2x', 'KYY--C2X',
    'KYY--Y1x', 'KYY--Y1X', 'KYY--Y2x', 'KYY--Y2X', 'KNF--C1x', 'KNF--C1X',
    'KNF--C2x', 'KNF--C2X', 'KNF--Y1x', 'KNF--Y1X', 'KNF--Y2x', 'KNF--Y2X',
    'KYF--C1x', 'KYF--C1X', 'KYF--C2x', 'KYF--C2X', 'KYF--Y1x', 'KYF--Y1X',
    'KYF--Y2x', 'KYF--Y2X', 'KNYNYC1x', 'KNYNYC1X', 'KNYNYC2x', 'KNYNYC2X',
    'KNYNYY1x', 'KNYNYY1X', 'KNYNYY2x', 'KNYNYY2X', 'KYYYYC1x', 'KYYYYC1X',
    'KYYYYC2x', 'KYYYYC2X', 'KYYYYY1x', 'KYYYYY1X', 'KYYYYY2x', 'KYYYYY2X',
    'KNFNFC1x', 'KNFNFC1X', 'KNFNFC2x', 'KNFNFC2X', 'KNFNFY1x', 'KNFNFY1X',
    'KNFNFY2x', 'KNFNFY2X', 'KYFYFC1x', 'KYFYFC1X', 'KYFYFC2x', 'KYFYFC2X',
    'KYFYFY1x', 'KYFYFY1X', 'KYFYFY2x', 'KYFYFY2X', 'TNY--C1x', 'TNY--C1X',
    'TNY--C2x', 'TNY--C2X', 'TNY--Y1x', 'TNY--Y1X', 'TNY--Y2x', 'TNY--Y2X',
    'TYY--C1x', 'TYY--C1X', 'TYY--C2x', 'TYY--C2X', 'TYY--Y1x', 'TYY--Y1X',
    'TYY--Y2x', 'TYY--Y2X', 'TNF--C1x', 'TNF--C1X', 'TNF--C2x', 'TNF--C2X',
    'TNF--Y1x', 'TNF--Y1X', 'TNF--Y2x', 'TNF--Y2X', 'TYF--C1x', 'TYF--C1X',
    'TYF--C2x', 'TYF--C2X', 'TYF--Y1x', 'TYF--Y1X', 'TYF--Y2x', 'TYF--Y2X',
    'TNYNYC1x', 'TNYNYC1X', 'TNYNYC2x', 'TNYNYC2X', 'TNYNYY1x', 'TNYNYY1X',
    'TNYNYY2x', 'TNYNYY2X', 'TYYYYC1x', 'TYYYYC1X', 'TYYYYC2x', 'TYYYYC2X',
    'TYYYYY1x', 'TYYYYY1X', 'TYYYYY2x', 'TYYYYY2X', 'TNFNFC1x', 'TNFNFC1X',
    'TNFNFC2x', 'TNFNFC2X', 'TNFNFY1x', 'TNFNFY1X', 'TNFNFY2x', 'TNFNFY2X',
    'TYFYFC1x', 'TYFYFC1X', 'TYFYFC2x', 'TYFYFC2X', 'TYFYFY1x', 'TYFYFY1X',
    'TYFYFY2x', 'TYFYFY2X'
]

REPORTDAYS = [
    0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366, 397, 425, 456,
    486, 517, 547, 578, 609, 639, 670, 700, 731, 762, 790, 821, 851, 882, 912,
    943, 974, 1004, 1035, 1065, 1096, 1127, 1155, 1186, 1216, 1247, 1277, 1308,
    1339, 1369, 1400, 1430, 1461, 1492, 1521, 1552, 1582, 1613, 1643, 1674,
    1705, 1735, 1766, 1796, 1827, 1858, 1886, 1917, 1947, 1978, 2008, 2039,
    2070, 2100, 2131, 2161, 2192, 2223, 2251, 2282, 2312, 2343, 2373, 2404,
    2435, 2465, 2496, 2526, 2557, 2588, 2616, 2647, 2677, 2708, 2738, 2769,
    2800, 2830, 2861, 2891, 2922, 2953, 2982, 3013, 3043, 3074, 3104, 3135,
    3166, 3196, 3227, 3257, 3288, 3319, 3347, 3378, 3408, 3439, 3469, 3500,
    3531, 3561, 3592, 3622, 3653, 3684, 3712, 3743, 3773, 3804, 3834, 3865,
    3896, 3926, 3957, 3987, 4018, 4049, 4077, 4108, 4138, 4169, 4199, 4230,
    4261, 4291, 4322, 4352, 4383, 4414, 4443, 4474, 4504, 4535, 4565, 4596,
    4627, 4657, 4688, 4718, 4749, 4780, 4808, 4839, 4869, 4900, 4930, 4961,
    4992, 5022, 5053, 5083, 5114, 5145, 5173, 5204, 5234, 5265, 5295, 5326,
    5357, 5387, 5418, 5448, 5479, 5510, 5538, 5569, 5599, 5630, 5660, 5691,
    5722, 5752, 5783, 5813, 5844, 5875, 5904, 5935, 5965, 5996, 6026, 6057,
    6088, 6118, 6149, 6179, 6210, 6241, 6269, 6300, 6330, 6361, 6391, 6422,
    6453, 6483, 6514, 6544, 6575, 6606, 6634, 6665, 6695, 6726, 6756, 6787,
    6818, 6848, 6879, 6909, 6940, 6971, 6999, 7030, 7060, 7091, 7121, 7152,
    7183, 7213, 7244, 7274, 7305, 7336, 7365, 7396, 7426, 7457, 7487, 7518,
    7549, 7579, 7610, 7640, 7671, 7702, 7730, 7761, 7791, 7822, 7852, 7883,
    7914, 7944, 7975, 8005, 8036, 8067, 8095, 8126, 8156, 8187, 8217, 8248,
    8279, 8309, 8340, 8370, 8401, 8432, 8460, 8491, 8521, 8552, 8582, 8613,
    8644, 8674, 8705, 8735, 8766, 8797, 8826, 8857, 8887, 8918, 8948, 8979,
    9010, 9040, 9071, 9101, 9132, 9163, 9191, 9222, 9252, 9283, 9313, 9344,
    9375, 9405, 9436, 9466, 9497, 9528, 9556, 9587, 9617, 9648, 9678, 9709,
    9740, 9770, 9801, 9831, 9862, 9893, 9921, 9952, 9982, 10013, 10043, 10074,
    10105, 10135, 10166, 10196, 10227, 10258, 10287, 10318, 10348, 10379,
    10409, 10440, 10471, 10501, 10532, 10562, 10593, 10624, 10652, 10683,
    10713, 10744, 10774, 10805, 10836, 10866, 10897, 10927, 10958
]

FIVE_YR_CYCLING_SWITCH_DAYS = [
    0, 3653, # Burn-in
    5478,    # First Five-Year Period
    7303,    # Second Period
    9128,    # Third
    10953,   # Fourth / Last
    10958    # Residue due to 365-and-calendar difference
]

