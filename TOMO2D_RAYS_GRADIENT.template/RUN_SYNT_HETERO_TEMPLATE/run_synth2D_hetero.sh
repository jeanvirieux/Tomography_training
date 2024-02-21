##########################################################
#
# simple example of computing synthetic "observed" travel times ...
#
# a way to compute observed data
#
##########################################################
# create of a synthetic model: constant velocity + exponential anomaly
# model_hetero_ref.dat (for tomo) and model_hetero_ref.bin (for plotting)
##########################################################
if [ $# -eq 2 ]
then

#########################################################
# create gnu file for plotting
#########################################################
echo "set term jpeg" >dessin.gnu
echo "set size ratio 2" >> dessin.gnu
echo "set palette functions (sin(gray*2*pi)+1)/2, (sin(gray*2*pi + 2*pi/3)+1)/2,(sin(gray*2*pi + 4*pi/3)+1)/2" >> dessin.gnu
echo "unset key" >> dessin.gnu
echo "set output 'MAP.jpg" >> dessin.gnu
echo "set cbrange [1980.:2020.]" >> dessin.gnu
echo "set xrange [0:$1]" >> dessin.gnu
echo "set yrange [0:$2]" >> dessin.gnu
echo "set pm3d map" >> dessin.gnu
echo "set rmargin at screen 0.85" >> dessin.gnu
echo "set lmargin at screen 0.15" >> dessin.gnu
echo "splot 'model_hetero_ref.bin' binary array=$1x$2"  >> dessin.gnu

nsrc=190
nrec=190

../BIN/model_hetero <<EOD
$1 $2       ! nx,nz
1. 1.         ! dx,dz
2000.         ! constant velocity
50.          ! amplitude of exp
50. 100.      ! xpos,zpos of anomaly
20. 20.       ! xleng,zleng
EOD
gnuplot dessin.gnu
#
#
#
../BIN/line_array <<EOD
$nsrc       ! nber of sources 
2. 2.    ! first position (x,z)
1.       ! dz vertical step
0.       ! dx horizontal step
EOD
mv fort.7 sources.dat
#
#
#
../BIN/line_array <<EOD
$nrec       ! nber of receivers
98. 2.   ! first position  (x,z)
1.       ! dz vertical step
0.       ! dx horizontal step
EOD
mv fort.7 receivers.dat

#rm toto
#echo "1" >toto
#ndof=`echo $nsrc $nrec | awk '{print $1*$2+1}' `
#awk -v n=$ndof '{for(i=1;i<=n;i++) print 1 1 1}' toto  >data_init.dat
#rm toto

#
# build a fake data file
#
../BIN/temps_homo <<EOD
0. 0. 1. 1.
$1 $2
receivers.dat
sources.dat
model_hetero_ref.dat
EOD
mv fcal.dat data_init.dat

#########################################################
#  construction of synthetic data using format of data.dat
#########################################################
echo "0                      ! option inversion 1 forward 0"  > inputs
echo "0. 0. 1. 1.            ! xorg,zorg, dx,dz            " >> inputs
echo "$1 $2                  ! nx,nz                       " >> inputs 
echo "model_hetero_ref.dat                                 " >> inputs
echo "receivers.dat                                        " >> inputs
echo "sources.dat                                          " >> inputs
echo "data_init.dat                                        " >> inputs
echo "10 0.5                 ! iter_max, damp              " >> inputs

../BIN/tomo2D_gradient
mv fcal.dat data_fwd_hetero.dat

else

##########################
# warning messages here 
##########################

if [ $# -eq 1 ]
then
if [ $1 = -h ]
then
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo please enter number of points along x and number of points along z
echo for this example please specify 101 201 
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
exit 0
fi
fi

#########################
# missing arguments
#########################

echo we need two integers ... nx and nz: type -h for more hints

fi
