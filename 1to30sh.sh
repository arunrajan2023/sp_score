for i in {1..30..1}; do sed "s/PPP/$i/g" PPP_ether_script_compute.sh >"$i"_ether_script_compute.sh
done
