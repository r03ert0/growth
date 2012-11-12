for name in $(cat log2/name.txt)
do
	name=$name
	cat $name | awk 'Begin{i=0;flagEnd=0}{i++;if($1=="Collision"||$1=="Saving")flagEnd=1; if(i>11 && !flagEnd)printf("%f ",$6)}END{printf("\n")}'
done
