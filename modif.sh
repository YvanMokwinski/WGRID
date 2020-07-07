who=`ls src/*.cpp`
for i in $who;do
    grep "wmeshspacedg.hpp" $i 2>&1 >/dev/null;
    a=$?
    if [ $a -ne 1 ];then
	echo "mod "$i;
	cat $i | sed 's/wmeshspacedg.hpp/wmeshspacedg_t.hpp/g' > tmp;
	mv tmp $i;
    fi
done
	 
