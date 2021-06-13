def find_matches(text, patterns):
	matches = []
	for i in range(len(text)):
		for pattern in patterns:
			if i+len(pattern) <= len(text) and text[i:i+len(pattern)] == pattern:
				matches.append(i)
	return matches

def 

if __name__ == "__main__":
	with open("dataset_317416_4 (2).txt", "r") as f:
		content = f.read().strip("\n").split("\n")
		text = content[0]
		patterns = content[1:]
	# print (text)
	# print (patterns)
	# print (find_matches(text, patterns))

	text = [9, 24, 33, 39, 51, 55, 59, 60, 61, 62, 63, 64, 65, 66, 79, 85, 93, 94, 95, 96, 97, 98, 107, 115, 116, 126, 141, 142, 156, 157, 158, 159, 165, 179, 180, 181, 182, 183, 184, 185, 186, 192, 195, 196, 197, 198, 202, 209, 214, 215, 216, 217, 218, 219, 220, 221, 242, 243, 244, 246, 247, 248, 249, 256, 257, 258, 259, 260, 261, 262, 263, 264, 283, 284, 285, 286, 296, 297, 298, 299, 300, 301, 302, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 326, 327, 331, 344, 348, 349, 350, 356, 357, 359, 372, 382, 383, 384, 385, 386, 387, 388, 389, 398, 413, 423, 425, 426, 427, 428, 429, 443, 444, 459, 468, 469, 470, 471, 472, 473, 502, 506, 507, 514, 515, 516, 517, 530, 531, 541, 548, 549, 555, 556, 557, 558, 559, 560, 561, 562, 576, 577, 586, 590, 591, 592, 593, 594, 595, 597, 598, 606, 607, 608, 609, 617, 618, 619, 620, 621, 622, 623, 624, 635, 636, 637, 638, 639, 640, 662, 663, 664, 665, 666, 667, 668, 669, 685, 686, 695, 696, 697, 698, 699, 702, 712, 713, 714, 715, 716, 717, 718, 719, 737, 744, 745, 746, 755, 768, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 787, 788, 797, 798, 799, 809, 810, 816, 822, 832, 848, 849, 850, 851, 852, 853, 854, 855, 856, 872, 884, 885, 886, 894, 904, 905, 927, 928, 929, 930, 931, 932, 933, 934, 959, 972, 973, 981, 982, 999, 1000, 1001, 1002, 1003, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1041, 1042, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1079, 1080, 1087, 1088, 1089, 1090, 1091, 1099, 1107, 1117, 1119, 1120, 1121, 1122, 1123, 1124, 1138, 1142, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1168, 1169, 1184, 1185, 1186, 1195, 1196, 1197, 1198, 1199, 1200, 1201, 1202, 1213, 1221, 1222, 1230, 1231, 1256, 1259, 1278, 1280, 1283, 1289, 1290, 1291, 1292, 1293, 1294, 1305, 1306, 1308, 1309, 1310, 1311, 1312, 1313, 1336, 1340, 1341, 1342, 1343, 1344, 1345, 1346, 1347, 1366, 1367, 1368, 1369, 1370, 1371, 1372, 1373, 1381, 1382, 1383, 1384, 1385, 1386, 1388, 1404, 1405, 1408, 1409, 1410, 1411, 1412, 1425, 1426, 1428, 1429, 1440, 1441, 1442, 1443, 1444, 1445, 1446, 1447, 1456, 1457, 1458, 1459, 1460, 1461, 1462, 1463, 1464, 1465, 1468, 1486, 1487, 1488, 1490, 1491, 1492, 1493, 1494, 1514, 1518, 1519, 1524, 1528, 1545, 1546, 1547, 1548, 1549, 1550, 1551, 1552, 1564, 1565, 1566, 1567, 1568, 1569, 1570, 1571, 1572, 1600, 1601, 1602, 1603, 1604, 1608, 1617, 1618, 1619, 1628, 1641, 1650, 1651, 1652, 1661, 1662, 1668, 1678, 1682, 1686, 1687, 1688, 1689, 1690, 1705, 1706, 1707, 1716, 1719, 1720, 1721, 1722, 1723, 1727, 1731, 1732, 1736, 1741, 1742, 1743, 1744, 1745, 1746, 1747, 1748, 1749, 1750, 1760, 1761, 1762, 1763, 1764, 1765, 1774, 1775, 1776, 1777, 1778, 1779, 1780, 1781, 1796, 1804, 1805, 1806, 1807, 1808, 1809, 1810, 1811, 1830, 1840, 1841, 1842, 1843, 1844, 1857, 1858, 1883, 1884, 1885, 1886, 1887, 1888, 1889, 1890, 1906, 1907, 1908, 1909, 1910, 1911, 1912, 1913, 1940, 1947, 1948, 1949, 1950, 1958, 1959, 1960, 1961, 1962, 1968, 1981, 1982, 1993, 2005, 2006, 2007, 2008, 2010, 2018, 2019, 2020, 2021, 2030, 2031, 2032, 2033, 2034, 2035, 2036, 2037, 2046, 2047, 2048, 2049, 2050, 2051, 2053, 2057, 2060, 2061, 2062, 2063, 2064, 2065, 2066, 2067, 2074, 2075, 2086, 2102, 2105, 2111, 2112, 2113, 2114, 2115, 2116, 2117, 2118, 2119, 2120, 2134, 2135, 2136, 2137, 2138, 2139, 2140, 2141, 2142, 2146, 2157, 2158, 2159, 2160, 2161, 2162, 2171, 2180, 2181, 2182, 2183, 2184, 2185, 2186, 2187, 2190, 2194, 2199, 2202, 2222, 2223, 2234, 2235, 2236, 2237, 2238, 2239, 2240, 2241, 2242, 2251, 2252, 2253, 2254, 2255, 2256, 2268, 2269, 2270, 2271, 2272, 2273, 2274, 2275, 2291, 2293, 2294, 2295, 2297, 2298, 2301, 2307, 2319, 2329, 2330, 2331, 2332, 2333, 2334, 2335, 2336, 2337, 2346, 2360, 2375, 2376, 2377, 2378, 2379, 2380, 2381, 2382, 2385, 2386, 2396, 2397, 2428, 2429, 2430, 2431, 2432, 2433, 2434, 2435, 2436, 2439, 2442, 2443, 2446, 2447, 2448, 2449, 2450, 2451, 2452, 2453, 2462, 2463, 2472, 2473, 2474, 2476, 2477, 2478, 2479, 2491, 2492, 2493, 2494, 2495, 2496, 2497, 2498, 2503, 2504, 2509, 2531, 2532, 2533, 2536, 2540, 2555, 2556, 2557, 2558, 2559, 2560, 2561, 2562, 2572, 2573, 2574, 2575, 2576, 2577, 2578, 2579, 2596, 2597, 2598, 2599, 2604, 2605, 2606, 2607, 2608, 2609, 2610, 2611, 2622, 2623, 2624, 2636, 2653, 2657, 2658, 2662, 2663, 2664, 2665, 2674, 2676, 2677, 2678, 2679, 2680, 2681, 2699, 2705, 2712, 2713, 2719, 2720, 2721, 2722, 2723, 2731, 2732, 2733, 2734, 2735, 2736, 2737, 2738, 2739, 2740, 2748, 2762, 2764, 2773, 2782, 2783, 2784, 2785, 2786, 2787, 2788, 2804, 2805, 2806, 2818, 2824, 2838, 2850, 2851, 2852, 2853, 2854, 2855, 2857, 2858, 2861, 2862, 2865, 2868, 2869, 2872, 2873, 2882, 2883, 2884, 2885, 2886, 2887, 2888, 2889, 2898, 2907, 2908, 2909, 2910, 2911, 2913, 2914, 2920, 2921, 2925, 2926, 2927, 2928, 2929, 2930, 2931, 2932, 2933, 2934, 2947, 2961, 2964, 2965, 2966, 2967, 2968, 2969, 2970, 2971, 2972, 2979, 2980, 2981, 2982, 2987, 2991, 2992, 2993, 2994, 2995, 2996, 2997, 2998, 3017, 3018, 3019, 3020, 3021, 3022, 3023, 3024, 3030, 3031, 3034, 3051, 3052, 3053, 3054, 3055, 3056, 3057, 3058, 3059, 3060, 3068, 3072, 3083, 3087, 3091, 3103, 3104, 3105, 3106, 3107, 3108, 3109, 3114, 3118, 3129, 3130, 3131, 3132, 3133, 3134, 3135, 3136, 3137, 3138, 3160, 3161, 3162, 3163, 3164, 3165, 3166, 3167, 3168, 3169, 3187, 3188, 3189, 3190, 3191, 3192, 3193, 3194, 3213, 3226, 3227, 3254, 3255, 3256, 3257, 3258, 3267, 3268, 3269, 3270, 3271, 3272, 3273, 3274, 3282, 3288, 3289, 3290, 3291, 3292, 3293, 3294, 3295, 3296, 3309, 3310, 3311, 3312, 3313, 3314, 3315, 3316, 3325, 3326, 3329, 3338, 3340, 3341, 3342, 3343, 3344, 3345, 3347, 3360, 3361, 3362, 3363, 3364, 3365, 3366, 3374, 3375, 3376, 3377, 3378, 3379, 3380, 3381, 3390, 3391, 3392, 3393, 3394, 3395, 3396, 3397, 3405, 3406, 3407, 3408, 3409, 3423, 3437, 3445, 3447, 3448, 3454, 3464, 3473, 3474, 3475, 3476, 3477, 3478, 3479, 3480, 3481, 3482, 3483, 3490, 3493, 3500, 3501, 3502, 3503, 3504, 3513, 3527, 3534, 3535, 3547, 3548, 3549, 3550, 3551, 3552, 3553, 3554, 3566, 3598, 3599, 3600, 3602, 3603, 3604, 3605, 3606, 3607, 3612, 3613, 3614, 3623, 3624, 3634, 3635, 3636, 3647, 3655, 3658, 3659, 3667, 3668, 3669, 3670, 3671, 3672, 3673, 3674, 3687, 3701, 3706, 3707, 3712, 3719, 3720, 3723, 3729, 3730, 3731, 3732, 3733, 3734, 3743, 3751, 3762, 3763, 3764, 3765, 3766, 3767, 3768, 3769, 3795, 3796, 3797, 3798, 3799, 3800, 3801, 3802, 3814, 3815, 3816, 3826, 3827, 3828, 3829, 3830, 3831, 3832, 3833, 3852, 3853, 3861, 3870, 3888, 3889, 3891, 3892, 3893, 3894, 3895, 3896, 3904, 3905, 3906, 3907, 3908, 3909, 3910, 3911, 3912, 3931, 3934, 3935, 3936, 3937, 3938, 3940, 3941, 3952, 3953, 3968, 3981, 3982, 3983, 3984, 3985, 3986, 3988, 3995, 3996, 3999, 4010, 4011, 4031, 4032, 4044, 4057, 4067, 4068, 4069, 4070, 4071, 4082, 4091, 4092, 4093, 4094, 4095, 4096, 4097, 4109, 4110, 4111, 4112, 4113, 4114, 4115, 4116, 4119, 4127, 4130, 4131, 4132, 4133, 4134, 4135, 4136, 4137, 4159, 4172, 4173, 4174, 4175, 4176, 4177, 4178, 4179, 4180, 4181, 4182, 4185, 4188, 4189, 4195, 4209, 4210, 4218, 4219, 4220, 4221, 4222, 4231, 4234, 4243, 4248, 4251, 4252, 4261, 4262, 4263, 4264, 4270, 4271, 4272, 4273, 4275, 4276, 4290, 4305, 4306, 4307, 4324, 4337, 4351, 4360, 4361, 4362, 4363, 4375, 4388, 4389, 4390, 4391, 4392, 4407, 4408, 4409, 4410, 4411, 4412, 4413, 4414, 4419, 4428, 4436, 4437, 4450, 4465, 4466, 4468, 4469, 4470, 4471, 4472, 4495, 4496, 4497, 4509, 4510, 4511, 4512, 4513, 4514, 4516, 4519, 4520, 4525, 4526, 4534, 4535, 4536, 4537, 4538, 4539, 4540, 4541, 4558, 4559, 4560, 4561, 4562, 4563, 4564, 4565, 4574, 4575, 4576, 4577, 4578, 4579, 4580, 4581, 4594, 4602, 4605, 4606, 4616, 4617, 4618, 4619, 4620, 4621, 4651, 4652, 4661, 4662, 4663, 4664, 4665, 4676, 4677, 4700, 4701, 4702, 4703, 4704, 4705, 4707, 4708, 4709, 4710, 4719, 4720, 4721, 4722, 4723, 4724, 4726, 4727, 4744, 4751, 4752, 4757, 4758, 4759, 4760, 4761, 4762, 4763, 4764, 4773, 4774, 4780, 4787, 4788, 4814, 4815, 4824, 4829, 4830, 4831, 4832, 4833, 4834, 4835, 4836, 4837, 4858, 4881, 4882, 4883, 4884, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4896, 4897, 4898, 4899, 4900, 4901, 4917, 4918, 4919, 4920, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4942, 4954, 4965, 4966, 4967, 4968, 4969, 4970, 4971, 4972, 4977, 4981, 4982, 4983, 4984, 5009, 5027, 5030, 5031, 5032, 5033, 5041, 5042, 5047, 5056, 5057, 5069, 5070, 5071, 5072, 5073, 5074, 5075, 5076, 5077, 5088, 5089, 5090, 5091, 5092, 5093, 5094, 5095, 5096, 5097, 5110, 5111, 5112, 5113, 5114, 5115, 5116, 5117, 5118, 5130, 5131, 5132, 5133, 5134, 5135, 5136, 5137, 5138, 5147, 5155, 5171, 5181, 5182, 5188, 5193, 5194, 5195, 5196, 5204, 5205, 5207, 5208, 5209, 5210, 5219, 5220, 5222, 5223, 5224, 5225, 5226, 5236, 5237, 5243, 5248, 5251, 5261, 5262, 5281, 5291, 5303, 5320, 5321, 5322, 5323, 5324, 5326, 5327, 5342, 5343, 5344, 5345, 5346, 5347, 5348, 5349, 5350, 5358, 5363, 5370, 5371, 5373, 5374, 5375, 5376, 5377, 5387, 5392, 5405, 5419, 5428, 5429, 5430, 5431, 5432, 5440, 5441, 5442, 5443, 5444, 5445, 5449, 5450, 5454, 5455, 5456, 5465, 5475, 5476, 5477, 5478, 5479, 5480, 5481, 5482, 5490, 5491, 5492, 5493, 5501, 5502, 5503, 5513, 5514, 5515, 5529, 5535, 5536, 5537, 5538, 5539, 5540, 5541, 5542, 5548, 5551, 5559, 5560, 5561, 5562, 5563, 5564, 5565, 5566, 5567, 5587, 5588, 5589, 5611, 5626, 5637, 5638, 5639, 5640, 5641, 5642, 5643, 5644, 5654, 5655, 5659, 5684, 5693, 5694, 5695, 5696, 5697, 5698, 5699, 5711, 5712, 5713, 5714, 5715, 5716, 5717, 5718, 5739, 5760, 5769, 5770, 5778, 5779, 5782, 5798, 5809, 5810, 5811, 5812, 5813, 5814, 5815, 5816, 5817, 5829, 5830, 5831, 5832, 5836, 5843, 5855, 5856, 5857, 5867, 5868, 5881, 5882, 5883, 5891, 5892, 5893, 5894, 5895, 5896, 5897, 5898, 5899, 5900, 5907, 5908, 5909, 5910, 5911, 5912, 5920, 5921, 5922, 5924, 5925, 5926, 5927, 5954, 5957, 5966, 5968, 5969, 5970, 5991, 6000, 6001, 6002, 6003, 6004, 6005, 6008, 6017, 6018, 6022, 6025, 6026, 6027, 6028, 6029, 6030, 6031, 6032, 6033, 6054, 6055, 6064, 6065, 6066, 6067, 6068, 6069, 6071, 6080, 6083, 6087, 6088, 6089, 6094, 6095, 6102, 6103, 6104, 6105, 6114, 6117, 6118, 6119, 6120, 6121, 6122, 6123, 6124, 6125, 6132, 6133, 6134, 6135, 6151, 6152, 6153, 6154, 6155, 6156, 6157, 6158, 6163, 6164, 6167, 6168, 6170, 6171, 6172, 6173, 6186, 6187, 6188, 6189, 6198, 6206, 6210, 6215, 6216, 6217, 6218, 6219, 6220, 6229, 6230, 6231, 6234, 6240, 6241, 6242, 6243, 6244, 6245, 6246, 6247, 6255, 6266, 6285, 6291, 6296, 6304, 6308, 6309, 6310, 6311, 6312, 6313, 6314, 6315, 6316, 6329, 6330, 6333, 6336, 6337, 6343, 6352, 6353, 6354, 6364, 6365, 6366, 6367, 6368, 6369, 6370, 6371, 6372, 6373, 6389, 6390, 6391, 6392, 6393, 6394, 6395, 6403, 6418, 6419, 6427, 6435, 6436, 6437, 6438, 6439, 6440, 6445, 6449, 6450, 6468, 6469, 6470, 6478, 6487, 6488, 6499, 6500, 6501, 6502, 6510, 6511, 6512, 6513, 6514, 6515, 6516, 6517, 6518, 6519, 6522, 6527, 6528, 6535, 6536, 6537, 6538, 6539, 6559, 6560, 6578, 6579, 6592, 6603, 6604, 6612, 6613, 6614, 6620, 6631, 6632, 6633, 6635, 6636, 6637, 6638, 6639, 6645, 6646, 6651, 6655, 6663, 6674, 6675, 6676, 6677, 6685, 6697, 6698, 6714, 6715, 6716, 6717, 6718, 6719, 6720, 6721, 6730, 6731, 6732, 6733, 6755, 6759, 6763, 6768, 6774, 6775, 6778, 6781, 6796, 6801, 6802, 6803, 6804, 6805, 6815, 6821, 6822, 6841, 6842, 6843, 6844, 6845, 6846, 6847, 6848, 6857, 6863, 6864, 6865, 6866, 6867, 6868, 6869, 6870, 6871, 6890, 6897, 6898, 6899, 6900, 6901, 6902, 6903, 6904, 6913, 6914, 6921, 6929, 6935, 6939, 6940, 6941, 6942, 6951, 6955, 6960, 6970, 6971, 6972, 6994, 6999, 7000, 7013, 7014, 7015, 7016, 7017, 7018, 7019, 7020, 7027, 7032, 7045, 7046, 7047, 7048, 7049, 7050, 7051, 7052, 7072, 7073, 7074, 7075, 7076, 7077, 7078, 7079, 7085, 7089, 7090, 7091, 7092, 7093, 7094, 7095, 7096, 7097, 7098, 7099, 7100, 7113, 7114, 7115, 7116, 7117, 7118, 7119, 7120, 7127, 7134, 7135, 7136, 7137, 7138, 7139, 7140, 7141, 7161, 7165, 7170, 7177, 7183, 7184, 7185, 7186, 7187, 7188, 7190, 7199, 7205, 7206, 7207, 7208, 7216, 7217, 7235, 7236, 7237, 7238, 7239, 7240, 7242, 7269, 7282, 7283, 7295, 7309, 7322, 7333, 7336, 7337, 7338, 7340, 7341, 7342, 7343, 7344, 7367, 7368, 7369, 7370, 7371, 7372, 7373, 7374, 7384, 7401, 7402, 7403, 7404, 7405, 7406, 7407, 7408, 7415, 7420, 7421, 7422, 7423, 7424, 7425, 7426, 7427, 7428, 7429, 7455, 7456, 7457, 7458, 7459, 7460, 7461, 7462, 7463, 7471, 7472, 7473, 7474, 7475, 7476, 7477, 7485, 7494, 7495, 7496, 7497, 7510, 7519, 7520, 7521, 7530, 7536, 7537, 7538, 7539, 7540, 7541, 7544, 7548, 7558, 7565, 7578, 7579, 7580, 7581, 7589, 7591, 7592, 7593, 7594, 7595, 7596, 7606, 7608, 7609, 7610, 7611, 7612, 7613, 7622, 7625, 7634, 7642, 7650, 7651, 7665, 7666, 7669, 7680, 7681, 7682, 7690, 7695, 7699, 7700, 7708, 7709, 7716, 7732, 7740, 7741, 7742, 7743, 7744, 7745, 7746, 7747, 7748, 7758, 7764, 7765, 7775, 7776, 7785, 7786, 7787, 7788, 7797, 7798, 7799, 7800, 7801, 7802, 7811, 7812, 7813, 7814, 7815, 7816, 7817, 7818, 7827, 7853, 7860, 7867, 7878, 7879, 7886, 7887, 7888, 7889, 7893, 7896, 7897, 7898, 7910, 7913, 7914, 7915, 7923, 7924, 7926, 7933, 7934, 7935, 7936, 7937, 7938, 7939, 7940, 7948, 7953, 7954, 7958, 7973, 7978, 7979, 7986, 7987, 7988, 7990, 7991, 7992, 7993, 8011, 8016, 8025, 8026, 8027, 8028, 8029, 8030, 8041, 8054, 8055, 8056, 8057, 8058, 8059, 8060, 8061, 8062, 8079, 8080, 8081, 8082, 8083, 8084, 8085, 8096, 8105, 8121, 8134, 8135, 8136, 8137, 8138, 8139, 8140, 8141, 8154, 8156, 8157, 8163, 8164, 8165, 8166, 8167, 8168, 8175, 8180, 8184, 8185, 8186, 8187, 8188, 8189, 8190, 8191, 8202, 8203, 8210, 8211, 8212, 8213, 8214, 8215, 8216, 8217, 8223, 8226, 8227, 8245, 8246, 8265, 8266, 8267, 8268, 8269, 8270, 8271, 8272, 8281, 8282, 8283, 8292, 8304, 8307, 8311, 8320, 8339, 8343, 8348, 8349, 8350, 8375, 8385, 8392, 8393, 8394, 8395, 8396, 8397, 8398, 8399, 8405, 8406, 8410, 8411, 8412, 8416, 8417, 8420, 8421, 8422, 8423, 8424, 8442, 8443, 8444, 8451, 8452, 8453, 8454, 8455, 8456, 8457, 8458, 8459, 8460, 8472, 8473, 8474, 8475, 8476, 8477, 8478, 8479, 8480, 8492, 8493, 8494, 8502, 8503, 8504, 8505, 8506, 8507, 8508, 8509, 8518, 8519, 8520, 8521, 8522, 8523, 8531, 8532, 8538, 8548, 8549, 8550, 8566, 8567, 8568, 8569, 8570, 8571, 8573, 8574, 8575, 8584, 8585, 8586, 8588, 8589, 8593, 8598, 8599, 8600, 8601, 8602, 8603, 8604, 8605, 8606, 8621, 8623, 8624, 8630, 8631, 8632, 8633, 8634, 8635]
	for i in range(len(text)):
		text[i] = str(text[i])
	print (" ".join(text))