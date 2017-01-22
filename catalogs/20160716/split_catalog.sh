master=clusterUDGlenses.list
clean=clusterUDGlenses_cleaned.list

# read header
head=`awk '{if ($1=="#") print}' $master`
echo $head

# create master cleaned file
awk -v x=$clean '{if ($1!="#" && $11>=5) print > x}' $master
# insert header
sed -i "1i $head" $clean

# read clusters, sort and take unique entries
clusters=`awk '{if ($1!="#") print $1}' $master | sort | uniq`
echo $clusters

# print content to each file
awk '{if ($1!="#") print > "all/"$1"_udgs.cat"}' $master
# also print cleaned files
awk '{if ($1!="#") print > "good/"$1"_good_udgs.cat"}' $clean
# this one will not print the first column
#awk '{for (i=2; i<=NF; i++) if($1!="#") print $i " " > "all/"$1"_udgs.cat"}' $master

# insert header
for cluster in $clusters
do
    # must use double quotes for the variable to be interpreted properly
    sed -i "1i $head" "all/${cluster}_udgs.cat"
    sed -i "1i $head" "good/${cluster}_good_udgs.cat"
done
