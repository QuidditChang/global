age1=$1
age2=$2
incr=$3
cap=0
while [ $cap -le 11 ]
do
    if [ $cap -le 9 ]; then
	str=0$cap
    else
	str=$cap
    fi
    echo $str
    k=$age1
    while [ $k -le $age2 ]
    do
        awk '{f1=1.0; for(i=1;i<25;i++) f1=f1-$i; print f1,$0}' global.opt$str.$k > temp00.$k
        awk '{max=$1; col=1; for(i=2;i<=25;i++) {if(max<$i) {max=$i; col=i};} print col-1}' temp00.$k | sed '1 s/.*//' > tracer00.$k
        paste global.cap$str.$k tracer00.$k > temp01.$k
	tail -n +2 temp01.$k > global.comp$str.$k
        mv global.comp$str.$k global.$str.$k.txt
    
        rm temp00.$k temp01.$k tracer00.$k
    
        #lat=1.57
        #echo $lat | xargs -I {} awk '{if($1=={}) print $2/3.1415926*180.0,$3,$9}' sub.comp00.$k > comp00
    
        #gmtset PAPER_MEDIA letter
        #gmtset PAGE_ORIENTATION portrait
    
        #LONMIN=0
        #LONMAX=60
        #RADMIN=0.75
        #RADMAX=1.0
        #REGION=$LONMIN/$LONMAX/$RADMIN/$RADMAX
    
        #pivot=`echo $LONMIN $LONMAX | awk '{print ($1+$2)/2}'`
        #Proj=Pa15/$pivot
    
        #makecpt -Cpolar -T0/7/1 -Z -I > tmp.cpt
        #surface comp00 -R$REGION -I0.1/0.002 -Gcomp.grd  -Ll-3000.0 -Lu3000.0 -T0.8
        #grdimage comp.grd -J$Proj -R$REGION -Ctmp.cpt -Ba10f5/a10f5 -X3 -Y5.0 -K -P > composition.$((age)).eps
        #psscale -Ctmp.cpt -D7/-1/5.0/0.50h -Ba100f50:"Composition":/:%%: -O -K -P >> composition.$((age)).eps
    
        let "k+=incr"
    
    done
    let "cap=cap+1"
done
