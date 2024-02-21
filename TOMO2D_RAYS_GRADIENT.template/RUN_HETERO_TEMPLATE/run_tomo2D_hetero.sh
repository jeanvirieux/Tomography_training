##########################################################
#
# simple example of computing travel times ...
#
##########################################################
# create of a synthetic model: constant velocity + exponential anomaly
# model_hetero_ref.dat (for tomo) and model_hetero_ref.bin (for plotting)
##########################################################
if [ $# -eq 2 ]
then

echo "set term jpeg" >target.gnu
echo "set size ratio 2" >> target.gnu
echo "set palette functions (sin(gray*2*pi)+1)/2, (sin(gray*2*pi + 2*pi/3)+1)/2,(sin(gray*2*pi + 4*pi/3)+1)/2" >> target.gnu
echo "unset key" >> target.gnu
echo "set output 'MAP_TARGET.jpg" >> target.gnu
echo "set cbrange [1980.:2020.]" >> target.gnu
echo "set xrange [0:$1]" >> target.gnu
echo " set yrange [0:$2]" >> target.gnu
echo "set pm3d map" >> target.gnu
echo "splot 'model_hetero_ref.bin' binary array=$1x$2"  >> target.gnu

echo "set term jpeg" >final.gnu
echo "set size ratio 2" >> final.gnu
echo "set palette functions (sin(gray*2*pi)+1)/2, (sin(gray*2*pi + 2*pi/3)+1)/2,(sin(gray*2*pi + 4*pi/3)+1)/2" >> final.gnu
echo "unset key" >> final.gnu
echo "set output 'MAP_FIN.jpg" >> final.gnu
echo "set cbrange [1980.:2020.]" >> final.gnu
echo "set xrange [0:$1]" >> final.gnu
echo " set yrange [0:$2]" >> final.gnu
echo "set pm3d map" >> final.gnu
echo "splot 'model.bin' binary array=$1x$2"  >> final.gnu

echo "set term jpeg" >init.gnu
echo "set size ratio 2" >> init.gnu
echo "set palette functions (sin(gray*2*pi)+1)/2, (sin(gray*2*pi + 2*pi/3)+1)/2,(sin(gray*2*pi + 4*pi/3)+1)/2" >> init.gnu
echo "unset key" >> init.gnu
echo "set output 'MAP_INI.jpg" >> init.gnu
echo "set cbrange [1980.:2020.]" >> init.gnu
echo "set xrange [0:$1]" >> init.gnu
echo " set yrange [0:$2]" >> init.gnu
echo "set pm3d map" >> init.gnu
echo "splot 'model_homo.bin' binary array=$1x$2"  >> init.gnu


nsrc=190
nrec=190
#######################################
#   compute the model to be recovered
#######################################
../BIN/model_hetero <<EOD
$1 $2       ! nx,nz
1. 1.         ! dx,dz
2000.         ! constant velocity
50.          ! amplitude of exp
50. 100.      ! xpos,zpos of anomaly
20. 20.       ! xleng,zleng
EOD
############################# plot the target model
gnuplot target.gnu
#######################################
#   compute the initial model
#######################################
../BIN/model_homo <<EOD
$1 $2         ! nx,nz
2000.         ! constant velocity
EOD
mv model_homo.dat model_ini.dat
############################# plot the initial model
gnuplot init.gnu
#######################################
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

#########################################################
# copy the synthetic data from RUN_SYNT_HETERO
#########################################################
cp ../RUN_SYNT_HETERO/data_fwd_hetero.dat data.dat
#########################################################
#  construction of synthetic data using format of data.dat
#########################################################
echo "1                      ! option inversion 1 forward 0"  > inputs
echo "0. 0. 1. 1.            ! xorg,zorg, dx,dz            " >> inputs
echo "$1 $2                  ! nx,nz                       " >> inputs 
echo "model_ini.dat                                        " >> inputs
echo "receivers.dat                                        " >> inputs
echo "sources.dat                                          " >> inputs
echo "data.dat                                             " >> inputs
echo "10 0.5                 ! iter_max, damp              " >> inputs

../BIN/tomo2D_gradient
mv fcal.dat data_final.dat
gnuplot final.gnu

else

##########################
# warning messages here 
##########################

if [ $# -eq 1 ]
then
if [ $1 = h ]
then
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo please enter number of points along x and number of points along z
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
exit 0
fi
fi

#########################
# missing arguments
#########################

echo we need two integers ... nx and nz

fi
