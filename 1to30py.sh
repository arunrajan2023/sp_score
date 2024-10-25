for i in {1..30..1}; do sed "s/PPP/$i/g" PPP__ether_test_pscore.py >"$i"_ether_test_pscore.py; done
