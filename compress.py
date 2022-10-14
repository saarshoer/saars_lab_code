import os
import glob
import pandas as pd
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers

input_dir = '/net/mraid08/export/mb/MBPipeline/Analyses/MBSNP/Gut/CountsS'
output_dir = '/net/mraid08/export/jasmine/CountsS_compressed_new'
file_name = 'mb_snp_counts_Rep_*.h5'

not_updated = ['Rep_2697', 'Rep_417', 'Rep_330', 'Rep_1316', 'Rep_1317', 'Rep_880', 'Rep_2696', 'Rep_1287', 'Rep_1335', 'Rep_1305', 'Rep_347', 'Rep_1300', 'Rep_1813', 'Rep_1787', 'Rep_1292', 'Rep_1315', 'Rep_376', 'Rep_883', 'Rep_1306', 'Rep_868', 'Rep_866', 'Rep_864', 'Rep_858', 'Rep_10', 'Rep_867', 'Rep_897', 'Rep_1802', 'Rep_2860', 'Rep_296', 'Rep_1339', 'Rep_11', 'Rep_1302', 'Rep_1295', 'Rep_855', 'Rep_1289', 'Rep_1245', 'Rep_817', 'Rep_2274', 'Rep_795', 'Rep_2871', 'Rep_1067', 'Rep_1803', 'Rep_1056', 'Rep_1079', 'Rep_794', 'Rep_2468', 'Rep_317', 'Rep_1155', 'Rep_2008', 'Rep_649', 'Rep_1046', 'Rep_882', 'Rep_1303', 'Rep_343', 'Rep_1877', 'Rep_1057', 'Rep_2858', 'Rep_814', 'Rep_1725', 'Rep_2893', 'Rep_2471', 'Rep_2193', 'Rep_2890', 'Rep_1183', 'Rep_2244', 'Rep_2797', 'Rep_2866', 'Rep_1703', 'Rep_2800', 'Rep_793', 'Rep_2845', 'Rep_3519', 'Rep_1016', 'Rep_1707', 'Rep_2985', 'Rep_1041', 'Rep_2966', 'Rep_329', 'Rep_2186', 'Rep_1088', 'Rep_1973', 'Rep_297', 'Rep_1267', 'Rep_2867', 'Rep_1967', 'Rep_166', 'Rep_1719', 'Rep_2905', 'Rep_1809', 'Rep_892', 'Rep_2785', 'Rep_861', 'Rep_2913', 'Rep_1452', 'Rep_415', 'Rep_828', 'Rep_1615', 'Rep_2505', 'Rep_1185', 'Rep_2097', 'Rep_3506', 'Rep_2501', 'Rep_57', 'Rep_203', 'Rep_3453', 'Rep_1294', 'Rep_2096', 'Rep_573', 'Rep_2894', 'Rep_3531', 'Rep_170', 'Rep_826', 'Rep_2457', 'Rep_1126', 'Rep_2973', 'Rep_2001', 'Rep_574', 'Rep_3476', 'Rep_798', 'Rep_1354', 'Rep_545', 'Rep_2803', 'Rep_268', 'Rep_2815', 'Rep_2811', 'Rep_2488', 'Rep_1282', 'Rep_73', 'Rep_1940', 'Rep_2901', 'Rep_2910', 'Rep_3518', 'Rep_1408', 'Rep_3091', 'Rep_2569', 'Rep_447', 'Rep_1701', 'Rep_2210', 'Rep_3143', 'Rep_3488', 'Rep_2891', 'Rep_390', 'Rep_1007', 'Rep_1054', 'Rep_3445', 'Rep_803', 'Rep_3537', 'Rep_17', 'Rep_1169', 'Rep_2819', 'Rep_2859', 'Rep_254', 'Rep_1201', 'Rep_2825', 'Rep_1952', 'Rep_2880', 'Rep_16', 'Rep_2752', 'Rep_373', 'Rep_1949', 'Rep_2877', 'Rep_768', 'Rep_1796', 'Rep_2887', 'Rep_2884', 'Rep_2888', 'Rep_2927', 'Rep_2916', 'Rep_77', 'Rep_2907', 'Rep_167', 'Rep_2915', 'Rep_1218', 'Rep_1076', 'Rep_2904', 'Rep_2892', 'Rep_914', 'Rep_2812', 'Rep_2801', 'Rep_13', 'Rep_2428', 'Rep_1595', 'Rep_3379', 'Rep_1332', 'Rep_3204', 'Rep_2224', 'Rep_1702', 'Rep_26', 'Rep_2799', 'Rep_223', 'Rep_697', 'Rep_2805', 'Rep_1665', 'Rep_257', 'Rep_1105', 'Rep_325', 'Rep_2375', 'Rep_2900', 'Rep_1194', 'Rep_265', 'Rep_2154', 'Rep_1055', 'Rep_1106', 'Rep_878', 'Rep_145', 'Rep_1290', 'Rep_1724', 'Rep_245', 'Rep_833', 'Rep_3462', 'Rep_1664', 'Rep_1099', 'Rep_3159', 'Rep_3538', 'Rep_1504', 'Rep_2016', 'Rep_1027', 'Rep_406', 'Rep_24', 'Rep_2827', 'Rep_218', 'Rep_292', 'Rep_3482', 'Rep_1883', 'Rep_1090', 'Rep_2430', 'Rep_1814', 'Rep_270', 'Rep_2872', 'Rep_361', 'Rep_3373', 'Rep_2490', 'Rep_1307', 'Rep_165', 'Rep_2412', 'Rep_3543', 'Rep_1050', 'Rep_1274', 'Rep_2971', 'Rep_15', 'Rep_1100', 'Rep_569', 'Rep_3565', 'Rep_1668', 'Rep_1177', 'Rep_1286', 'Rep_3487', 'Rep_2885', 'Rep_2699', 'Rep_3564', 'Rep_1165', 'Rep_1301', 'Rep_1309', 'Rep_3553', 'Rep_159', 'Rep_168', 'Rep_1505', 'Rep_3536', 'Rep_1529', 'Rep_2127', 'Rep_3124', 'Rep_1863', 'Rep_2153', 'Rep_1972', 'Rep_1254', 'Rep_2199', 'Rep_1971', 'Rep_69', 'Rep_2146', 'Rep_3503', 'Rep_2145', 'Rep_354', 'Rep_1981', 'Rep_2198', 'Rep_896', 'Rep_1795', 'Rep_1784', 'Rep_1690', 'Rep_3594', 'Rep_1742', 'Rep_2842', 'Rep_114', 'Rep_1333', 'Rep_1336', 'Rep_2817', 'Rep_1605', 'Rep_2802', 'Rep_147', 'Rep_2798', 'Rep_2419', 'Rep_1606', 'Rep_3470', 'Rep_2794', 'Rep_2863', 'Rep_3426', 'Rep_176', 'Rep_479', 'Rep_907', 'Rep_575', 'Rep_570', 'Rep_568', 'Rep_1304', 'Rep_865', 'Rep_1608', 'Rep_1104', 'Rep_1188', 'Rep_872', 'Rep_944', 'Rep_1153', 'Rep_1375', 'Rep_1310', 'Rep_1371', 'Rep_2830', 'Rep_2816', 'Rep_2754', 'Rep_1901', 'Rep_1320', 'Rep_1314', 'Rep_1323', 'Rep_1319', 'Rep_1187', 'Rep_911', 'Rep_1900', 'Rep_1326', 'Rep_1313', 'Rep_1324', 'Rep_1329', 'Rep_1331', 'Rep_1102', 'Rep_1330', 'Rep_2965', 'Rep_910', 'Rep_894', 'Rep_1192', 'Rep_808', 'Rep_3361', 'Rep_198', 'Rep_1318', 'Rep_857', 'Rep_856', 'Rep_812', 'Rep_1036', 'Rep_1806', 'Rep_1346', 'Rep_2862', 'Rep_1311', 'Rep_1350', 'Rep_1321', 'Rep_1066', 'Rep_1808', 'Rep_2059', 'Rep_1327', 'Rep_1807', 'Rep_3360', 'Rep_3407', 'Rep_2695', 'Rep_2060', 'Rep_154', 'Rep_1045', 'Rep_1322', 'Rep_152', 'Rep_811', 'Rep_2824', 'Rep_271', 'Rep_1969', 'Rep_2809', 'Rep_1547', 'Rep_2829', 'Rep_1569', 'Rep_906', 'Rep_1366', 'Rep_1308', 'Rep_2831', 'Rep_1383', 'Rep_1291', 'Rep_936', 'Rep_1006', 'Rep_860', 'Rep_1299', 'Rep_1348', 'Rep_2492', 'Rep_1035', 'Rep_1285', 'Rep_1152', 'Rep_862', 'Rep_863', 'Rep_920', 'Rep_3563', 'Rep_810', 'Rep_155', 'Rep_809', 'Rep_1288', 'Rep_149', 'Rep_1887', 'Rep_3319', 'Rep_2349', 'Rep_68', 'Rep_1811', 'Rep_1097', 'Rep_1334', 'Rep_1667', 'Rep_2010', 'Rep_824', 'Rep_1328', 'Rep_362', 'Rep_1781', 'Rep_2045', 'Rep_1968', 'Rep_1970', 'Rep_481', 'Rep_1098', 'Rep_1149', 'Rep_1812', 'Rep_1794', 'Rep_1951', 'Rep_815', 'Rep_85', 'Rep_157']
skip_species = not_updated


def compress(in_file, out_file):
    with pd.HDFStore(in_file, 'r') as hdf:
        for key in hdf.keys():
            hdf[key].to_hdf(out_file, key=key, complevel=9)


os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)
sethandlers()

with qp(jobname='comp') as q:

    q.startpermanentrun()
    tkttores = {}

    for input_file in glob.glob(os.path.join(input_dir, file_name)):
        output_file = os.path.join(output_dir, os.path.basename(input_file))
        species = os.path.basename(input_file).replace('mb_snp_counts_', '').replace('.h5', '')
        if species not in skip_species:
            tkttores[input_file] = q.method(compress, (input_file, output_file))

    for k, v in tkttores.items():
        q.waitforresult(v)
