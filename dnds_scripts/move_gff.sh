in_dir=../results/dnds_out/agat_fix

for clique_path in $in_dir/*
do
        clique=$(echo $clique_path|sed "s|$in_dir/||g")

        mkdir $out_basedir/$clique

        # Bring gffs to main clique directory
        for i in $clique_path/*
        do
                cp $i/*.gff $clique_path
       	done
done
