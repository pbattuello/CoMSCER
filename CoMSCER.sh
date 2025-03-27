#!/usr/bin/bash


echo -e "\n\n\t=====> CoMSCER Analysis <=====\n\n"

if [ $# -ne 6 ]; then
	echo "==============================================================="
	echo -e "\n\n!!! Needed 6 arguments !!!\n--> 1: Mutational Matrix1\n--> 2: Mutational Matrix2\n--> 3: Output Directory\n--> 4: Cancer Specific SBS mutational signatures from Cosmic2\n--> 5: Cancer Specific SBS mutational signatures from Cosmic3.2\n--> 6: coding non-coding [0-1]\n\n"
	echo "==============================================================="
	exit 1
fi


mutational_file1=$1
mutational_file2=$2
working_dir=$3
sbs_to_evaluate_C2=$4
sbs_to_evaluate_C3_2=$5
non_coding=$6
script_path=$(dirname $(realpath $0))



sbs_cosmic2=$(echo $sbs_to_evaluate_C2 | awk '{split($1, a, ","); for (i in a) sub(/SBS/, "Signature_"); print}')
sbs_cosmic3_2=$(echo $sbs_to_evaluate_C3_2) 

echo -e "--> Start analysis..."

mkdir -p $working_dir/TMP
mkdir -p $working_dir/data
mkdir -p $working_dir/Figures
mkdir -p $working_dir/data/logs

mkdir -p $working_dir/Cosine_Similarity
mkdir -p $working_dir/Cosine_Similarity/matrix1
mkdir -p $working_dir/Cosine_Similarity/matrix2

mkdir -p $working_dir/SBS_signature_contributions
mkdir -p $working_dir/SBS_signature_contributions/matrix1
mkdir -p $working_dir/SBS_signature_contributions/matrix2

mkdir -p $working_dir/TMP/matrix1
mkdir -p $working_dir/TMP/matrix2

if [ $non_coding -eq 1 ]; then
	mkdir -p $working_dir/Cosine_Similarity/matrix1_coding
	mkdir -p $working_dir/Cosine_Similarity/matrix2_coding
	mkdir -p $working_dir/Cosine_Similarity/matrix1_non_coding
	mkdir -p $working_dir/Cosine_Similarity/matrix2_non_coding

	mkdir -p $working_dir/SBS_signature_contributions/matrix1_coding
	mkdir -p $working_dir/SBS_signature_contributions/matrix2_coding
	mkdir -p $working_dir/SBS_signature_contributions/matrix1_non_coding
	mkdir -p $working_dir/SBS_signature_contributions/matrix2_non_coding

	mkdir -p $working_dir/TMP/matrix1_coding
	mkdir -p $working_dir/TMP/matrix2_coding
	mkdir -p $working_dir/TMP/matrix1_non_coding
	mkdir -p $working_dir/TMP/matrix2_non_coding
fi

echo -e "--> Creating mutational matrix..."

awk 'BEGIN{FS=OFS="\t"}{print "signature", $1, ".", $4, "SNP", $2, $3, $3, $5, $6, "SOMATIC"}' $mutational_file1 > $working_dir/TMP/matrix1/mut_matrix1.tmp.txt
awk 'BEGIN{FS=OFS="\t"}{print "signature", $1, ".", $4, "SNP", $2, $3, $3, $5, $6, "SOMATIC"}' $mutational_file2 > $working_dir/TMP/matrix2/mut_matrix2.tmp.txt

genome=$(cut -f 4 $working_dir/TMP/matrix1/mut_matrix1.tmp.txt | sort | uniq)

if [ $non_coding -eq 1 ]; then
	
	#subset coding and non coding mutations
	awk -v wd=$working_dir 'BEGIN{FS=OFS="\t"}{\
		if($7=="non_coding"){print "signature", $1, ".", $4, "SNP", $2, $3, $3, $5, $6, "SOMATIC" > wd"/TMP/matrix1_non_coding/mut_matrix1_non_coding.tmp.txt"}\
		if($7=="coding"){print "signature", $1, ".", $4, "SNP", $2, $3, $3, $5, $6, "SOMATIC" > wd"/TMP/matrix1_coding/mut_matrix1_coding.tmp.txt"}\
		}' $mutational_file1
	
	awk -v wd=$working_dir 'BEGIN{FS=OFS="\t"}{\
		if($7=="non_coding"){print "signature", $1, ".", $4, "SNP", $2, $3, $3, $5, $6, "SOMATIC" > wd"/TMP/matrix2_non_coding/mut_matrix2_non_coding.tmp.txt"}\
		if($7=="coding"){print "signature", $1, ".", $4, "SNP", $2, $3, $3, $5, $6, "SOMATIC" > wd"/TMP/matrix2_coding/mut_matrix2_coding.tmp.txt"}
		}' $mutational_file2

fi




$script_path/bin/matrix_generator.py "matrix1" $working_dir/TMP/matrix1 $genome > $working_dir/data/logs/MatrixGenerator.matrix1.log 2> $working_dir/data/logs/MatrixGenerator.matrix1.err
$script_path/bin/matrix_generator.py "matrix2" $working_dir/TMP/matrix2 $genome > $working_dir/data/logs/MatrixGenerator.matrix2.log 2> $working_dir/data/logs/MatrixGenerator.matrix2.err

cp $working_dir/TMP/matrix1/output/SBS/matrix1.SBS96.all $working_dir/data/.
cp $working_dir/TMP/matrix2/output/SBS/matrix2.SBS96.all $working_dir/data/.

mutational_matrix1=$working_dir/data/matrix1.SBS96.all
mutational_matrix2=$working_dir/data/matrix2.SBS96.all


if [ $non_coding -eq 1 ]; then

	#creating mutational matrix for coding and non-coding mutations
	$script_path/bin/matrix_generator.py "matrix1_coding" $working_dir/TMP/matrix1_coding $genome > $working_dir/data/logs/MatrixGenerator.matrix1_coding.log 2> $working_dir/data/logs/MatrixGenerator.matrix1_coding.err
	$script_path/bin/matrix_generator.py "matrix2_coding" $working_dir/TMP/matrix2_coding $genome > $working_dir/data/logs/MatrixGenerator.matrix2_coding.log 2> $working_dir/data/logs/MatrixGenerator.matrix2_coding.err

	$script_path/bin/matrix_generator.py "matrix1_non_coding" $working_dir/TMP/matrix1_non_coding $genome > $working_dir/data/logs/MatrixGenerator.matrix1_non_coding.log 2> $working_dir/data/logs/MatrixGenerator.matrix1_non_coding.err
	$script_path/bin/matrix_generator.py "matrix2_non_coding" $working_dir/TMP/matrix2_non_coding $genome > $working_dir/data/logs/MatrixGenerator.matrix2_non_coding.log 2> $working_dir/data/logs/MatrixGenerator.matrix2_non_coding.err


	cp $working_dir/TMP/matrix1_coding/output/SBS/matrix1_coding.SBS96.all $working_dir/data/.
	cp $working_dir/TMP/matrix2_coding/output/SBS/matrix2_coding.SBS96.all $working_dir/data/.

	cp $working_dir/TMP/matrix1_non_coding/output/SBS/matrix1_non_coding.SBS96.all $working_dir/data/.
	cp $working_dir/TMP/matrix2_non_coding/output/SBS/matrix2_non_coding.SBS96.all $working_dir/data/.


	mutational_matrix1_coding=$working_dir/data/matrix1_coding.SBS96.all
	mutational_matrix2_coding=$working_dir/data/matrix2_coding.SBS96.all

	mutational_matrix1_non_coding=$working_dir/data/matrix1_non_coding.SBS96.all
	mutational_matrix2_non_coding=$working_dir/data/matrix2_non_coding.SBS96.all

fi

echo -e "--> Done"
echo -e "--> Running mutational signature analysis with multiple tools..."
echo -e "--> Whole Mutational Matrix..."

Rscript $script_path/bin/mutationalpatterns.R $mutational_matrix1 $script_path $working_dir "matrix1" $genome > $working_dir/data/logs/MutationalPatterns.matrix1.log 2> $working_dir/data/logs/MutationalPatterns.matrix1.err
Rscript $script_path/bin/mutationalpatterns.R $mutational_matrix2 $script_path $working_dir "matrix2" $genome > $working_dir/data/logs/MutationalPatterns.matrix2.log 2> $working_dir/data/logs/MutationalPatterns.matrix2.err

Rscript $script_path/bin/deconstructsig.R $mutational_matrix1 $script_path $working_dir "matrix1" $genome > $working_dir/data/logs/deconstructSigs.matrix1.log 2> $working_dir/data/logs/deconstructSigs.matrix1.err
Rscript $script_path/bin/deconstructsig.R $mutational_matrix2 $script_path $working_dir "matrix2" $genome > $working_dir/data/logs/deconstructSigs.matrix2.log 2> $working_dir/data/logs/deconstructSigs.matrix2.err

Rscript $script_path/bin/signature.tools.lib.R $mutational_matrix1 $script_path $working_dir "matrix1" $genome > $working_dir/data/logs/signature.tools.lib.matrix1.log 2> $working_dir/data/logs/signature.tools.lib.matrix1.err
Rscript $script_path/bin/signature.tools.lib.R $mutational_matrix2 $script_path $working_dir "matrix2" $genome > $working_dir/data/logs/signature.tools.lib.matrix2.log 2> $working_dir/data/logs/signature.tools.lib.matrix2.err

python $script_path/bin/sigprofilerassignment.py $mutational_matrix1 $script_path $working_dir "matrix1" $genome > $working_dir/data/logs/sigprofilerassignment.matrix1.log 2> $working_dir/data/logs/sigprofilerassignment.matrix1.err
python $script_path/bin/sigprofilerassignment.py $mutational_matrix2 $script_path $working_dir "matrix2" $genome > $working_dir/data/logs/sigprofilerassignment.matrix2.log 2> $working_dir/data/logs/sigprofilerassignment.matrix2.err

mv $working_dir/TMP/matrix1/C2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix1/sbs_sig_C2_sigprofilerassignment.txt
mv $working_dir/TMP/matrix1/C3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix1/sbs_sig_C3.2_sigprofilerassignment.txt

mv $working_dir/TMP/matrix2/C2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix2/sbs_sig_C2_sigprofilerassignment.txt
mv $working_dir/TMP/matrix2/C3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix2/sbs_sig_C3.2_sigprofilerassignment.txt

Rscript $script_path/bin/cosine_sim_SPA.R $mutational_matrix1 $script_path $working_dir "matrix1"
Rscript $script_path/bin/cosine_sim_SPA.R $mutational_matrix2 $script_path $working_dir "matrix2"


###$script_path/bin/siganalyzer.sh $mutational_matrix $script_path $working_dir

###$script_path/bin/siganalyzer.py $mutational_matrix $script_path $working_dir

echo -e "--> Done"

if [ $non_coding -eq 1 ]; then

echo -e "--> Coding Genome Mutational Matrix..."

	Rscript $script_path/bin/mutationalpatterns.R $mutational_matrix1_coding $script_path $working_dir "matrix1_coding" $genome > $working_dir/data/logs/MutationalPatterns.matrix1_coding.log 2> $working_dir/data/logs/MutationalPatterns.matrix1_coding.err
	Rscript $script_path/bin/mutationalpatterns.R $mutational_matrix2_coding $script_path $working_dir "matrix2_coding" $genome > $working_dir/data/logs/MutationalPatterns.matrix2_coding.log 2> $working_dir/data/logs/MutationalPatterns.matrix2_coding.err

	Rscript $script_path/bin/deconstructsig.R $mutational_matrix1_coding $script_path $working_dir "matrix1_coding" $genome > $working_dir/data/logs/deconstructSigs.matrix1_coding.log 2> $working_dir/data/logs/deconstructSigs.matrix1_coding.err
	Rscript $script_path/bin/deconstructsig.R $mutational_matrix2_coding $script_path $working_dir "matrix2_coding" $genome > $working_dir/data/logs/deconstructSigs.matrix2_coding.log 2> $working_dir/data/logs/deconstructSigs.matrix2_coding.err

	Rscript $script_path/bin/signature.tools.lib.R $mutational_matrix1_coding $script_path $working_dir "matrix1_coding" $genome > $working_dir/data/logs/signature.tools.lib.matrix1_coding.log 2> $working_dir/data/logs/signature.tools.lib.matrix1_coding.err
	Rscript $script_path/bin/signature.tools.lib.R $mutational_matrix2_coding $script_path $working_dir "matrix2_coding" $genome > $working_dir/data/logs/signature.tools.lib.matrix2_coding.log 2> $working_dir/data/logs/signature.tools.lib.matrix2_coding.err

	python $script_path/bin/sigprofilerassignment.py $mutational_matrix1_coding $script_path $working_dir "matrix1_coding" $genome > $working_dir/data/logs/sigprofilerassignment.matrix1_coding.log 2> $working_dir/data/logs/sigprofilerassignment.matrix1_coding.err
	python $script_path/bin/sigprofilerassignment.py $mutational_matrix2_coding $script_path $working_dir "matrix2_coding" $genome > $working_dir/data/logs/sigprofilerassignment.matrix2_coding.log 2> $working_dir/data/logs/sigprofilerassignment.matrix2_coding.err

	mv $working_dir/TMP/matrix1_coding/C2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix1_coding/sbs_sig_C2_sigprofilerassignment.txt 
	mv $working_dir/TMP/matrix1_coding/C3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix1_coding/sbs_sig_C3.2_sigprofilerassignment.txt

	mv $working_dir/TMP/matrix2_coding/C2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix2_coding/sbs_sig_C2_sigprofilerassignment.txt
	mv $working_dir/TMP/matrix2_coding/C3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix2_coding/sbs_sig_C3.2_sigprofilerassignment.txt

	Rscript $script_path/bin/cosine_sim_SPA.R $mutational_matrix1_coding $script_path $working_dir "matrix1_coding"
	Rscript $script_path/bin/cosine_sim_SPA.R $mutational_matrix2_coding $script_path $working_dir "matrix2_coding"

	echo -e "--> Done"
	echo -e "--> Non - Coding Genome Mutational Matrix..."

	Rscript $script_path/bin/mutationalpatterns.R $mutational_matrix1_non_coding $script_path $working_dir "matrix1_non_coding" $genome > $working_dir/data/logs/MutationalPatterns.matrix1_non_coding.log 2> $working_dir/data/logs/MutationalPatterns.matrix1_non_coding.err
	Rscript $script_path/bin/mutationalpatterns.R $mutational_matrix2_non_coding $script_path $working_dir "matrix2_non_coding" $genome > $working_dir/data/logs/MutationalPatterns.matrix2_non_coding.log 2> $working_dir/data/logs/MutationalPatterns.matrix2_non_coding.err

	Rscript $script_path/bin/deconstructsig.R $mutational_matrix1_non_coding $script_path $working_dir "matrix1_non_coding" $genome > $working_dir/data/logs/deconstructSigs.matrix1_non_coding.log 2> $working_dir/data/logs/deconstructSigs.matrix1_non_coding.err
	Rscript $script_path/bin/deconstructsig.R $mutational_matrix2_non_coding $script_path $working_dir "matrix2_non_coding" $genome > $working_dir/data/logs/deconstructSigs.matrix2_non_coding.log 2> $working_dir/data/logs/deconstructSigs.matrix2_non_coding.err

	Rscript $script_path/bin/signature.tools.lib.R $mutational_matrix1_non_coding $script_path $working_dir "matrix1_non_coding" $genome > $working_dir/data/logs/signature.tools.lib.matrix1_non_coding.log 2> $working_dir/data/logs/signature.tools.lib.matrix1_non_coding.err
	Rscript $script_path/bin/signature.tools.lib.R $mutational_matrix2_non_coding $script_path $working_dir "matrix2_non_coding" $genome > $working_dir/data/logs/signature.tools.lib.matrix2_non_coding.log 2> $working_dir/data/logs/signature.tools.lib.matrix2_non_coding.err

	python $script_path/bin/sigprofilerassignment.py $mutational_matrix1_non_coding $script_path $working_dir "matrix1_non_coding" $genome > $working_dir/data/logs/sigprofilerassignment.matrix1_non_coding.log 2> $working_dir/data/logs/sigprofilerassignment.matrix1_non_coding.err
	python $script_path/bin/sigprofilerassignment.py $mutational_matrix2_non_coding $script_path $working_dir "matrix2_non_coding" $genome > $working_dir/data/logs/sigprofilerassignment.matrix2_non_coding.log 2> $working_dir/data/logs/sigprofilerassignment.matrix2_non_coding.err

	mv $working_dir/TMP/matrix1_non_coding/C2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix1_non_coding/sbs_sig_C2_sigprofilerassignment.txt
	mv $working_dir/TMP/matrix1_non_coding/C3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix1_non_coding/sbs_sig_C3.2_sigprofilerassignment.txt

	mv $working_dir/TMP/matrix2_non_coding/C2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix2_non_coding/sbs_sig_C2_sigprofilerassignment.txt
	mv $working_dir/TMP/matrix2_non_coding/C3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt $working_dir/TMP/matrix2_non_coding/sbs_sig_C3.2_sigprofilerassignment.txt

	Rscript $script_path/bin/cosine_sim_SPA.R $mutational_matrix1_non_coding $script_path $working_dir "matrix1_non_coding"
	Rscript $script_path/bin/cosine_sim_SPA.R $mutational_matrix2_non_coding $script_path $working_dir "matrix2_non_coding"

	echo -e "--> Done"

fi

echo -e "--> Integrating mutational signature results..."

if [ $non_coding -eq 1 ]; then

	Rscript $script_path/bin/integrative_signature_analysis_coding.R $working_dir $sbs_cosmic2 $sbs_cosmic3_2 > $working_dir/data/logs/integrative_signature_analysis_coding.log 2> $working_dir/data/logs/integrative_signature_analysis_coding.err

fi


if [ $non_coding -eq 0 ]; then
	
	Rscript $script_path/bin/integrative_signature_analysis.R $working_dir $sbs_cosmic2 $sbs_cosmic3_2 > $working_dir/data/logs/integrative_signature_analysis.log 2> $working_dir/data/logs/integrative_signature_analysis.err

fi


rm Rplots.pdf
rm -r $working_dir/TMP

echo -e "\n\n\t==> Analysis Completed <==\n\n"
