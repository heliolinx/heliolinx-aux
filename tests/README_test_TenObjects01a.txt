INVOCATION FOR make_tracklets:

time make_tracklets -dets test_TenObjects01a.csv -outimgs imgs_test_TenObjects01a.txt -pairdets pairdets_test_TenObjects01a.csv -tracklets tracklets_test_TenObjects01a.csv -trk2det trk2det_test_TenObjects01a.csv -colformat colformat_LSST_02.txt -imrad 2.0 -maxtime 1.5 -maxGCR 0.5 -maxvel 1.5 -earth Earth1day2020s_02a.csv -obscode ObsCodesNew.txt

EXPECTED FINAL OUTPUT AND TIMING:

Output image catalog imgs_test_TenObjects01a.txt, with 123 lines, has been written
Writing paired detection file pairdets_test_TenObjects01a.csv with 109 lines
Writing tracklet file tracklets_test_TenObjects01a.csv with 42 lines
Writing trk2det file trk2det_test_TenObjects01a.csv with 109 lines

real	0m0.063s
user	0m0.020s
sys	0m0.003s


INVOCATION FOR heliolinc:

time heliolinc -imgs imgs_test_TenObjects01a.txt -pairdets pairdets_test_TenObjects01a.csv -tracklets tracklets_test_TenObjects01a.csv -trk2det trk2det_test_TenObjects01a.csv -mjd 61109.09 -obspos Earth1day2020s_02a.csv -heliodist heliohyp_rmb00a.txt -clustrad 200000.0 -outsum sum_test_TenObjects01a.csv -clust2det clust2det_test_TenObjects01a.csv

EXPECTED FINAL OUTPUT AND TIMING:

De-duplicating output set of 53 candidate linkages
Final de-duplicated set contains 26 linkages
Writing 26 lines to output cluster-summary file sum_test_TenObjects01a.csv
Writing 269 lines to output clust2det file clust2det_test_TenObjects01a.csv

real	0m0.058s
user	0m0.010s
sys	0m0.007s


SETUP AND INVOCATION FOR link_planarity:

printf "sum_test_TenObjects01a.csv clust2det_test_TenObjects01a.csv\n" > clusterlist_test_TenObjects01a

time link_planarity -imgs imgs_test_TenObjects01a.txt -pairdet pairdets_test_TenObjects01a.csv -lflist clusterlist_test_TenObjects01a -mjd 61109.09 -simptype 1 -max_astrom_rms 0.5 -oop 10000.0 -ptpow 3 -maxrms 400000.0 -outsum LPLsum_test_TenObjects01a.csv -clust2det LPLclust2det_test_TenObjects01a.csv

EXPECTED FINAL OUTPUT AND TIMING:

Writing 10 lines to output cluster-summary file LPLsum_test_TenObjects01a.csv
Writing 105 lines to output clust2det file LPLclust2det_test_TenObjects01a.csv

real	0m0.041s
user	0m0.023s
sys	0m0.005s
