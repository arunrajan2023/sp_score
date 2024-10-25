

XXX="0"
for YYY in  {1..5..1} {7..14..1}
do
	echo $XXX $YYY
	#create py file
	cd tests
		sed  "s/XXX/$XXX/g" tmp_test_pscore.py > tmp                     
		sed  "s/YYY/$YYY/g" tmp > "$XXX_$YYY_test_pscore.py"
	cd ..
	#create slurm file
	sed "s/XXX/$XXX/g" tmp_script_compute.sh > tmp
	sed "s/YYY/$YYY/g" tmp > "$XXX_$YYY_script_compute.sh"
done


XXX="1"
for YYY in  {1..16..1}
do
	echo $XXX $YYY
	#create py file
	cd tests
		sed  "s/XXX/$XXX/g" tmp_test_pscore.py > tmp                     
		sed  "s/YYY/$YYY/g" tmp > "$XXX_$YYY_test_pscore.py"
	cd ..
	#create slurm file
	sed  "s/XXX/$XXX/g" tmp_script_compute.sh > tmp
	sed  "s/YYY/$YYY/g" tmp > "$XXX_$YYY_script_compute.sh"
done
