#
#
#
cp ../HamiltonFastMarching/bin/FileHFM_Isotropic2.exe .
#
#  python HFM_homo.py
python HFM_homo.py

#  python HFM.py 'file name' nz nx dcarre zorg xorg zsrc xsrc    (first vertical second horizontal)
rm input_* output_*
python HFM.py nankai_double.bin 67 365 468.75 -468.75 -468.75 1000. 100000.

#  receivers are described in HFM_all.py
rm input_* output_*
python HFM_all.py nankai_double.bin 67 365 468.75 -468.75 -468.75 1000. 100000.

# marmousi
rm input_* output_*
python HFM.py marmousi_double.bin 131 395 23.4375 -23.4375 -23.4375 1000. 1000.
#manipulate HFM_all.py for doing the two-points ray computation for marmousi
