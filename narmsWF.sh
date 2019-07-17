#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 16-18
#$ -cwd -V
#$ -o narmsWF.log
#$ -j y
#$ -N JessFinder
#$ -q all.q

### Jessica Corron Chen - Developed 05/15/2019 ###
### Last update - 05/15/2019 ###
### Usage: sh narmsWF.sh organismName ###
### Supported bugs: salmonella, escherichia, campylobacter, vibrio ###
### Update paths and modules before usage to local system prior to usage ###

# Do stuff

#create time stamp
time_stamp=$(date +%Y_%m_%d_%H_%M_%S)

if [ "$1"  == "salmonella" ]; then

	#configure environment for CG pipeline
	module purge
	export PATH=/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/bin/CG-Pipeline/scripts:$PATH
	module load perl/5.16.1-MT

	#compute assembly metrics
	run_assembly_readMetrics.pl *.fastq.gz --fast --numcpus 16 -e 5000000 | sort -k3,3n > readMetrics.tsv
	
	#sum across R1 and R2 read decks; divide by ten to get cov-cutoff
	cat readMetrics.tsv | awk 'BEGIN { FS = "\t"} {if ($1 != "File"){print $1, $9/10}}' | sed -e 's/_[1,2]\.fastq\.gz//' \
	| awk '{a[$1] += $2} END{for (i in a) print i, a[i]}' | awk '{print $1 ,"\t" , $2}' > readMetrics2.tsv

	#configure environment to run shovill
	module purge
	module load SPAdes/3.13.0 Skesa/2.3.0 megahit/1.1.2 velvet/1.2.10 lighter/1.1.1 flash/1.2.11 samtools/1.8 bwa/0.7.17  Mash/2.0  seqtk/1.2 pilon/1.22 trimmomatic/0.35
	export PATH=/scicomp/home/lly3/bin/pigz-2.4:$PATH
	export PATH=/scicomp/home/lly3/bin/samclip:$PATH
	export PATH=/scicomp/home/lly3/bin/shovill/bin:$PATH
	unset PERL5LIB

	#run shovill - trims, assembles, runs read correction; formats names, shovill downsamples to 100X, so we set the cov cutoff to 10 unless it's calculated to be less
	cat readMetrics2.tsv | awk '{if($2>=10) print $0}' | xargs -P 1 -n 2 sh -c 'shovill --outdir $0 --R1 $0_1.fastq.gz --R2 $0_2.fastq.gz --gsize 5000000 \
	--mincov 10 --trim --namefmt $0_contig%05d; mv $0/contigs.fa $0/$0.contigs.fa; mv $0/contigs.gfa $0/$0.contigs.gfa; \
	mv $0/shovill.corrections $0/$0.shovill.corrections; mv $0/shovill.log $0/$0.shovill.log; mv $0/spades.fasta $0/$0.spades.fasta'

	cat readMetrics2.tsv | awk '{if($2<10 && $2>=3) print $0}' | xargs -P 1 -n 2 sh -c 'shovill --outdir $0 --R1 $0_1.fastq.gz --R2 $0_2.fastq.gz --gsize 5000000 \
	--mincov $1 --trim --namefmt $0_contig%05d; mv $0/contigs.fa $0/$0.contigs.fa; mv $0/contigs.gfa $0/$0.contigs.gfa; \
	 mv $0/shovill.corrections $0/$0.shovill.corrections; mv $0/shovill.log $0/$0.shovill.log; mv $0/spades.fasta $0/$0.spades.fasta'

	cat readMetrics2.tsv | awk '{if($2<3) print $0}' > lowCoverageIsolates.salm."$time_stamp".tsv

	#move files; cleanup
	mkdir graphs; mv */*.contigs.gfa ./graphs
	mkdir shovill-logs; mv */*.shovill.log ./shovill-logs; mv */*.shovill.corrections ./shovill-logs
	mkdir assemblies; mv */*.contigs.fa ./assemblies; 
	mkdir assemblies-raw; mv */*.spades.fasta ./assemblies-raw
	cut -f 1 readMetrics2.tsv | xargs -P 1 -n 1 sh -c 'rmdir $0'
	cd assemblies
	
	#configure environment to run staramr
	module purge
	export PATH=/scicomp/home/lly3/miniconda3/bin:$PATH
	source activate staramr

	#run staramr 
	staramr search --pointfinder-organism salmonella --exclude-genes-file /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/NARST_research/scripts/databases/salm_genes_to_exclude.tsv --pid-threshold 90 --percent-length-overlap-resfinder 50  -o out.staramr *.fa
	conda deactivate

	#clean up output files
	cd out.staramr
	cat resfinder.tsv | cut -f 1-2 |  sed -e 's/\.contigs//' | sed -e 's/\ //' | sed '1d' > import.resfinder.tsv
	cat summary.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($2 == "None") print $1,$2 }' | sed -e 's/\.contigs//' >> import.resfinder.tsv 
	cat pointfinder.tsv | cut -f 1-2 | sed -e 's/\.contigs//' | sed -e 's/\ //' | sed '1d'| sed -e 's/([A-Z]/(/' |sed -e 's/[A-Z])/)/' > import.pointfinder.tsv
	cd ..
	
	#run plasmid finder with abricate, formats import sheet
	module purge
   	unset PERL5LIB
	export PATH=/scicomp/home/lly3/miniconda3/bin:$PATH
	mkdir out.plasmids
	ls *.fa | sed -e 's/\.fa//' | xargs -I one sh -c 'abricate --db KTPF --threads 12 --minid 90 --mincov 60 --noheader one.fa  > out.plasmids/one.pf.tsv'
	cd out.plasmids
	abricate --summary *.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | awk 'BEGIN {FS = "\t"} {OFS="\t"}{if ($2 == "0") {print $1,"None"}}' > import.plasmidfinder.tsv 
	cat *.contigs.pf.tsv | sed -e 's/\.contigs\.fa//' | tee plasmidfinder.tsv | awk 'BEGIN {FS = "\t"}{OFS="\t"}{print $1,$5}' >> import.plasmidfinder.tsv
	ls *.contigs.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | xargs -I one sh -c 'cat one.contigs.pf.tsv | sed -e 's/\.contigs\.fa//' | echo "$(cut -f 1 | head -n 1)"\
	 "$(cut -f 5 one.contigs.pf.tsv | paste -d, -s)">> plasmidsummary.tsv'
	abricate --summary *.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | awk 'BEGIN {FS = "\t"}{OFS="\t"}{if ($2 == "0") {print $1,"None"}}' >> plasmidsummary.tsv
	sed -i '/^ *$/d' plasmidsummary.tsv; 
	mkdir hits; mv *.contigs.pf.tsv ./hits; cd ..
	cut -f 1-2 out.staramr/summary.tsv | sed '1d' | sed -e 's/\.contigs//' | sed -e 's/\ //g' > out.staramr/summary2.tsv
	sort -k 1 -o out.staramr/summary2.tsv out.staramr/summary2.tsv; sort -k 1 -o out.plasmids/plasmidsummary.tsv out.plasmids/plasmidsummary.tsv; 
	join out.staramr/summary2.tsv out.plasmids/plasmidsummary.tsv | sed -e 's/\ /\t/g' > final.summary.tsv
	mv final.summary.tsv ../final.summary.salm."$time_stamp".tsv; mv out.plasmids/import.plasmidfinder.tsv ../import.plasmidfinder.salm."$time_stamp".tsv; mv \
	out.staramr/import.resfinder.tsv ../import.resfinder.salm."$time_stamp".tsv; mv out.staramr/import.pointfinder.tsv ../import.pointfinder.salm."$time_stamp".tsv
	
	#set up environment for seqsero2
	
	module purge
	export PATH='/bin:/sbin':"$PATH"
	export PATH=$PATH:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/out/Salmonella/SeqSero2_new/SeqSero2-master
	export PATH=$PATH:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/out/Salmonella/SeqSero2_new/SeqSero2-master/SalmID
	module load SeqSero/1.0 BEDTools/2.17.0 SPAdes/3.8.0 Python/3.5
	
	#run seqsero2
	cd .. ; ls *_1.fastq.gz | sed -e 's/_1\.fastq\.gz//' | xargs -I one sh -c 'SeqSero2_package.py -t 2 -p 16 -m a -i one_1.fastq.gz one_2.fastq.gz';

	for d in SeqSero_result*

	do 

		cd $d;
		cat Seqsero_result.txt | awk 'BEGIN { FS = "\t"} {if ($1=="Input files:" || $1 == "Predicted serotype:" || $1 == "Predicted antigenic profile:") {print $2}}' > Seqsero_result.1.tsv; 
		cd ..

	done
	paste -s  SeqSero_result*/Seqsero_result.1.tsv > allseqsero."$time_stamp".tsv
	
	
elif [ "$1" == "escherichia" ]; then
	
	#configure environment for CG pipeline
	module purge
	export PATH=/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/bin/CG-Pipeline/scripts:$PATH
	module load perl/5.16.1-MT
	
	#compute assembly metrics
	run_assembly_readMetrics.pl *.fastq.gz --fast --numcpus 16 -e 5000000 | sort -k3,3n > readMetrics.tsv
	
	#sum across R1 and R2 read decks; divide by two to get cov-cutoff
	cat readMetrics.tsv | awk 'BEGIN { FS = "\t"} {if ($1 != "File"){print $1, $9/10}}' | sed -e 's/_[1,2]\.fastq\.gz//' \
	| awk '{a[$1] += $2} END{for (i in a) print i, a[i]}' | awk '{print $1 ,"\t" , $2}' > readMetrics2.tsv
	
	#configure environment to run shovill
	module purge
	module load SPAdes/3.13.0 Skesa/2.3.0 megahit/1.1.2 velvet/1.2.10 lighter/1.1.1 flash/1.2.11 samtools/1.8 bwa/0.7.17  Mash/2.0  seqtk/1.2 pilon/1.22 trimmomatic/0.35
	export PATH=/scicomp/home/lly3/bin/pigz-2.4:$PATH
	export PATH=/scicomp/home/lly3/bin/samclip:$PATH
	export PATH=/scicomp/home/lly3/bin/shovill/bin:$PATH
	unset PERL5LIB

	#run shovill - trims, assembles, runs read correction; formats names
	cat readMetrics2.tsv | awk '{if($2>=10) print $0}' | xargs -P 1 -n 2 sh -c 'shovill --outdir $0 --R1 $0_1.fastq.gz --R2 $0_2.fastq.gz --gsize 5000000 \
	--mincov 10 --trim --namefmt $0_contig%05d; mv $0/contigs.fa $0/$0.contigs.fa; mv $0/contigs.gfa $0/$0.contigs.gfa; \
	mv $0/shovill.corrections $0/$0.shovill.corrections; mv $0/shovill.log $0/$0.shovill.log; mv $0/spades.fasta $0/$0.spades.fasta'

	cat readMetrics2.tsv | awk '{if($2<10 && $2>=4) print $0}' | xargs -P 1 -n 2 sh -c 'shovill --outdir $0 --R1 $0_1.fastq.gz --R2 $0_2.fastq.gz --gsize 5000000 \
	--mincov $1 --trim --namefmt $0_contig%05d; mv $0/contigs.fa $0/$0.contigs.fa; mv $0/contigs.gfa $0/$0.contigs.gfa; \
	mv $0/shovill.corrections $0/$0.shovill.corrections; mv $0/shovill.log $0/$0.shovill.log; mv $0/spades.fasta $0/$0.spades.fasta'
	
	cat readMetrics2.tsv | awk '{if($2<4) print $0}' > lowCoverageIsolates.ec."$time_stamp".tsv

	#move files; cleanup
	mkdir graphs; mv */*.contigs.gfa ./graphs
	mkdir shovill-logs; mv */*.shovill.log ./shovill-logs; mv */*.shovill.corrections ./shovill-logs
	mkdir assemblies; mv */*.contigs.fa ./assemblies; 
	mkdir assemblies-raw; mv */*.spades.fasta ./assemblies-raw
	cut -f 1 readMetrics2.tsv | xargs -P 1 -n 1 sh -c 'rmdir $0'
	cd assemblies
	
	#configure environment to run staramr
	module purge
	export PATH=/scicomp/home/lly3/miniconda3/bin:$PATH
	source activate staramr

	#run staramr 
	staramr search --exclude-genes-file /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/NARST_research/scripts/databases/ec_genes_to_exclude.tsv --pid-threshold 90 --percent-length-overlap-resfinder 50 -o out.staramr *.fa

	conda deactivate

	#clean up output files
	cd out.staramr
	cat resfinder.tsv | cut -f 1-2 |  sed -e 's/\.contigs//' | sed -e 's/\ //' | sed '1d' > import.resfinder.tsv
	cat summary.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($2 == "None") print $1,$2 }' | sed -e 's/\.contigs//' >> import.resfinder.tsv 
	#cat pointfinder.tsv | cut -f 1-2 | sed -e 's/\.contigs//' | sed -e 's/\ //' | sed '1d' > import.pointfinder.tsv
	cd ..
	
	#run plasmid finder with abricate, formats import sheet
	module purge
   	unset PERL5LIB
	export PATH=/scicomp/home/lly3/miniconda3/bin:$PATH
	mkdir out.plasmids
	ls *.fa | sed -e 's/\.fa//' | xargs -I one sh -c 'abricate --db KTPF --threads 12 --minid 90 --mincov 60 --noheader one.fa  > out.plasmids/one.pf.tsv'
	cd out.plasmids
	abricate --summary *.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | awk 'BEGIN {FS = "\t"}{OFS="\t"}{if ($2 == "0") {print $1,"None"}}' > import.plasmidfinder.tsv 
	cat *.contigs.pf.tsv | sed -e 's/\.contigs\.fa//' | tee plasmidfinder.tsv | awk 'BEGIN {FS = "\t"}{OFS="\t"}{print $1,$5}' >> import.plasmidfinder.tsv
	ls *.contigs.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | xargs -I one sh -c 'cat one.contigs.pf.tsv | sed -e 's/\.contigs\.fa//' | echo "$(cut -f 1 | head -n 1)" \
	"$(cut -f 5 one.contigs.pf.tsv | paste -d, -s)">> plasmidsummary.tsv'
	abricate --summary *.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | awk 'BEGIN {FS = "\t"} {OFS="\t"} {if ($2 == "0") {print $1,"None"}}' >> plasmidsummary.tsv
	sed -i '/^ *$/d' plasmidsummary.tsv; 
	mkdir hits; mv *.contigs.pf.tsv ./hits; cd ..
	cut -f 1-2 out.staramr/summary.tsv | sed '1d' | sed -e 's/\.contigs//' | sed -e 's/\ //g' > out.staramr/summary2.tsv
	sort -k 1 -o out.staramr/summary2.tsv out.staramr/summary2.tsv; sort -k 1 -o out.plasmids/plasmidsummary.tsv out.plasmids/plasmidsummary.tsv; 
	
	#run ariba to detect pointfinder mutations
	module purge
	module load ariba/2.12.0
	cd .. 
	mkdir mutational
	ls *_1.fastq.gz | sed -e 's/_1\.fastq\.gz//' | xargs -I one sh -c 'ariba run  /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/NARST_research/scripts/databases/ecpoints/ecpoints one_1.fastq.gz one_2.fastq.gz mutational/one.mut; mv mutational/one.mut/report.tsv mutational/one.mut/one.report.tsv' 
	cd mutational
	ls | sed -e 's/\.mut//' | while read id 
	do 
		cat $id.mut/$id.report.tsv | awk -v id="$id" 'BEGIN{OFS="\t"} {if ($18 == "1") {print id,$1"("$19")"}}' | sed -e 's/(./(/' | sed -e 's/.)/)/' >> import.mutational.tsv
		cat $id.mut/$id.report.tsv | awk -v id="$id" 'BEGIN{OFS="\t"} {if ($18 == "1") {print id,$1"("$19")"}}' >> $id.summary.mutational.tsv
		
	done 
	
	for g in *.summary.mutational.tsv

        do
                (if [ -s "$g" ]
                then
                        g1="${g%.}";
                        cat "$g" | echo -e  "${g1%.summary.mutational.tsv}\t""$(cut -f 2 | paste -d, -s)" >> all.summary.mutational.tsv
                else
                        g1="${g%.}";
                        cd ${g1%.summary.mutational.tsv}.mut; [ -f log.clusters.gz ] && { cd .. ; echo -e "${g1%.summary.mutational.tsv}\tNone" >> \
			all.summary.mutational.tsv; } || { cd ..; echo -e "${g1%.summary.mutational.tsv}\tError" >> all.summary.mutational.tsv; 
			echo -e "${g1%.summary.mutational.tsv}\tError" >> import.mutational.tsv; }
                fi)
        done

	sort -k 1 -o all.summary.mutational.tsv all.summary.mutational.tsv
	
	cd ../assemblies; 
	join out.staramr/summary2.tsv ../mutational/all.summary.mutational.tsv > determinants.tsv
	join determinants.tsv out.plasmids/plasmidsummary.tsv | sed -e 's/\ /\t/g' > final.summary.tsv
	
	mv final.summary.tsv ../final.summary.ec."$time_stamp".tsv; mv out.plasmids/import.plasmidfinder.tsv ../import.plasmidfinder.ec."$time_stamp".tsv; 
	mv out.staramr/import.resfinder.tsv ../import.resfinder.ec."$time_stamp".tsv ; mv ../mutational/import.mutational.tsv ../import.pointfinder.ec."$time_stamp".tsv

elif [ "$1" == "campylobacter" ]; then

	#configure environment for CG pipeline
	module purge
	export PATH=/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/bin/CG-Pipeline/scripts:$PATH
	module load perl/5.16.1-MT
	
	#compute assembly metrics
	run_assembly_readMetrics.pl *.fastq.gz --fast --numcpus 16 -e 1800000 | sort -k3,3n > readMetrics.tsv
	
	#sum across R1 and R2 read decks; divide by two to get cov-cutoff
	cat readMetrics.tsv | awk 'BEGIN { FS = "\t"} {if ($1 != "File"){print $1, $9/10}}' | sed -e 's/_[1,2]\.fastq\.gz//' \
	| awk '{a[$1] += $2} END{for (i in a) print i, a[i]}' | awk '{print $1 ,"\t" , $2}' > readMetrics2.tsv

	#configure environment to run shovill
	module purge
	module load SPAdes/3.13.0 Skesa/2.3.0 megahit/1.1.2 velvet/1.2.10 lighter/1.1.1 flash/1.2.11 samtools/1.8 bwa/0.7.17  Mash/2.0  seqtk/1.2 pilon/1.22 trimmomatic/0.35
	export PATH=/scicomp/home/lly3/bin/pigz-2.4:$PATH
	export PATH=/scicomp/home/lly3/bin/samclip:$PATH
	export PATH=/scicomp/home/lly3/bin/shovill/bin:$PATH
	unset PERL5LIB

	#run shovill - trims, assembles, runs read correction; formats names
	cat readMetrics2.tsv | awk '{if($2>=10) print $0}' | xargs -P 1 -n 2 sh -c 'shovill --outdir $0 --R1 $0_1.fastq.gz --R2 $0_2.fastq.gz --gsize 1800000 \
	--mincov 10 --trim --namefmt $0_contig%05d; mv $0/contigs.fa $0/$0.contigs.fa; mv $0/contigs.gfa $0/$0.contigs.gfa; mv \
	$0/shovill.corrections $0/$0.shovill.corrections; mv $0/shovill.log $0/$0.shovill.log; mv $0/spades.fasta $0/$0.spades.fasta'
	
	cat readMetrics2.tsv | awk '{if($2<10 && $2>=2) print $0}' | xargs -P 1 -n 2 sh -c 'shovill --outdir $0 --R1 $0_1.fastq.gz --R2 $0_2.fastq.gz --gsize 1800000 \
	--mincov $1 --trim --namefmt $0_contig%05d; mv $0/contigs.fa $0/$0.contigs.fa; mv $0/contigs.gfa $0/$0.contigs.gfa; mv \
	$0/shovill.corrections $0/$0.shovill.corrections; mv $0/shovill.log $0/$0.shovill.log; mv $0/spades.fasta $0/$0.spades.fasta'

	cat readMetrics2.tsv | awk '{if($2<2) print $0}' > lowCoverageIsolates.campy."$time_stamp".tsv

	#move files; cleanup
	mkdir graphs; mv */*.contigs.gfa ./graphs
	mkdir shovill-logs; mv */*.shovill.log ./shovill-logs; mv */*.shovill.corrections ./shovill-logs
	mkdir assemblies; mv */*.contigs.fa ./assemblies; 
	mkdir assemblies-raw; mv */*.spades.fasta ./assemblies-raw
	cut -f 1 readMetrics2.tsv | xargs -P 1 -n 1 sh -c 'rmdir $0'
	cd assemblies
	
	#configure environment to run staramr
	module purge
	export PATH=/scicomp/home/lly3/miniconda3/bin:$PATH
	source activate staramr

	#run staramr 
	staramr search --pointfinder-organism campylobacter --pid-threshold 90 --percent-length-overlap-resfinder 50 --exclude-genes-file /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/NARST_research/scripts/databases/campy_genes_to_exclude.tsv -o out.staramr *.fa
	
	conda deactivate

	#clean up output files
	cd out.staramr
	cat resfinder.tsv | cut -f 1-2 |  sed -e 's/\.contigs//' | sed -e 's/\ //' | sed '1d'> import.resfinder.tsv
	cat summary.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($2 == "None") print $1,$2 }' | sed -e 's/\.contigs//' >> import.resfinder.tsv 
	cat pointfinder.tsv | cut -f 1-2 | sed -e 's/\.contigs//' | sed -e 's/\ //' | sed '1d'| sed -e 's/([A-Z]/(/' |sed -e 's/[A-Z])/)/' > staramr.pointfinder.tsv
	cd ..
	
	#run plasmid finder with abricate, formats import sheet
	module purge
   	unset PERL5LIB
	export PATH=/scicomp/home/lly3/miniconda3/bin:$PATH
	mkdir out.plasmids
	ls *.fa | sed -e 's/\.fa//' | xargs -I one sh -c 'abricate --db KTPF --threads 12 --minid 90 --mincov 60 --noheader one.fa  > out.plasmids/one.pf.tsv'
	cd out.plasmids
	abricate --summary *.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | awk 'BEGIN {FS = "\t"} {OFS="\t"}{if ($2 == "0") {print $1,"None"}}' > import.plasmidfinder.tsv 
	cat *.contigs.pf.tsv | sed -e 's/\.contigs\.fa//' | tee plasmidfinder.tsv | awk 'BEGIN {FS = "\t"} {OFS="\t"}{print $1,$5}' >> import.plasmidfinder.tsv
	ls *.contigs.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | xargs -I one sh -c 'cat one.contigs.pf.tsv | sed -e 's/\.contigs\.fa//' | echo "$(cut -f 1 | head -n 1)" \
	"$(cut -f 5 one.contigs.pf.tsv | paste -d, -s)">> plasmidsummary.tsv'
	abricate --summary *.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | awk 'BEGIN {FS = "\t"}{OFS="\t"} {if ($2 == "0") {print $1,"None"}}' >> plasmidsummary.tsv
	sed -i '/^ *$/d' plasmidsummary.tsv; 
	mkdir hits; mv *.contigs.pf.tsv ./hits; cd ..
	cut -f 1-2 out.staramr/summary.tsv | sed '1d' | sed -e 's/\.contigs//' | sed -e 's/\ //g' > out.staramr/summary2.tsv
	sort -k 1 -o out.staramr/summary2.tsv out.staramr/summary2.tsv; sort -k 1 -o out.plasmids/plasmidsummary.tsv out.plasmids/plasmidsummary.tsv;  
	
	#run ariba to detect 23S mutations
	module purge
	module load ariba/2.12.0
	cd .. ; mkdir mutational
	ls *_1.fastq.gz | sed -e 's/_1\.fastq\.gz//' | xargs -I one sh -c 'ariba run /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/NARST_research/scripts/databases/23S/out.prepareref.23s one_1.fastq.gz one_2.fastq.gz mutational/one.mut; mv mutational/one.mut/report.tsv mutational/one.mut/one.report.tsv' 
	cd mutational
	ls | sed -e 's/\.mut//' | while read id 
	do 
		cat $id.mut/$id.report.tsv | awk -v id="$id" 'BEGIN {OFS="\t"}{if ($18 == "1") {print id,$1"("$21")"}}' >> import.23Sdetection.tsv
		cat $id.mut/$id.report.tsv | awk -v id="$id" 'BEGIN {OFS="\t"}{if ($18 == "1") {print id,$1"("$19")"}}' >> $id.summary.23Sdetection.tsv
	
	
	done 
	sort -k 1 -o import.23Sdetection.tsv import.23Sdetection.tsv
	
	for g in *.summary.23Sdetection.tsv
	do
		(if [ -s "$g" ]
		then
			g1="${g%.}";
			cat "$g" | echo -e "${g1%.summary.23Sdetection.tsv}\t""$(cut -f 2 | paste -d, -s)" >> all.summary.23Sdetection.tsv
		else 
			g1="${g%.}";
			cd ${g1%.summary.23Sdetection.tsv}.mut; [ -f log.clusters.gz ] && { cd .. ; echo -e "${g1%.summary.23Sdetection.tsv}\tNone" >> all.summary.23Sdetection.tsv; } || { cd ..; echo -e "${g1%.summary.23Sdetection.tsv}\tError" >> all.summary.23Sdetection.tsv;  echo -e "${g1%.summary.23Sdetection.tsv}\tError" >> import.23Sdetection.tsv;}
		fi)
	done
	
	sort -k 1 -o all.summary.23Sdetection.tsv all.summary.23Sdetection.tsv
	
	cd ../assemblies; cat out.staramr/staramr.pointfinder.tsv ../mutational/import.23Sdetection.tsv > out.staramr/import.pointfinder.tsv
	join out.staramr/summary2.tsv ../mutational/all.summary.23Sdetection.tsv > determinants.tsv
	join determinants.tsv out.plasmids/plasmidsummary.tsv | sed -e 's/\ /\t/g' > final.summary.tsv
	
	mv final.summary.tsv ../final.summary.campy."$time_stamp".tsv; mv out.plasmids/import.plasmidfinder.tsv ../import.plasmidfinder.campy."$time_stamp".tsv;
	mv out.staramr/import.resfinder.tsv ../import.resfinder.campy."$time_stamp".tsv; mv out.staramr/import.pointfinder.tsv ../import.pointfinder.campy."$time_stamp".tsv	

elif [ "$1" == "vibrio" ]; then

	#configure environment for CG pipeline
	module purge
	export PATH=/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/bin/CG-Pipeline/scripts:$PATH
	module load perl/5.16.1-MT
	
	#compute assembly metrics
	run_assembly_readMetrics.pl *.fastq.gz --fast --numcpus 16 -e 5000000 | sort -k3,3n > readMetrics.tsv
	
	#sum across R1 and R2 read decks; divide by two to get cov-cutoff
	cat readMetrics.tsv | awk 'BEGIN { FS = "\t"} {if ($1 != "File"){print $1, $9/10}}' | sed -e 's/_[1,2]\.fastq\.gz//' \
	| awk '{a[$1] += $2} END{for (i in a) print i, a[i]}' | awk '{print $1 ,"\t" , $2}' > readMetrics2.tsv

	#configure environment to run shovill
	module purge
	module load SPAdes/3.13.0 Skesa/2.3.0 megahit/1.1.2 velvet/1.2.10 lighter/1.1.1 flash/1.2.11 samtools/1.8 bwa/0.7.17  Mash/2.0  seqtk/1.2 pilon/1.22 trimmomatic/0.35
	export PATH=/scicomp/home/lly3/bin/pigz-2.4:$PATH
	export PATH=/scicomp/home/lly3/bin/samclip:$PATH
	export PATH=/scicomp/home/lly3/bin/shovill/bin:$PATH
	unset PERL5LIB

	#run shovill - trims, assembles, runs read correction; formats names
	cat readMetrics2.tsv | awk '{if($2>=10) print $0}' | xargs -P 1 -n 2 sh -c 'shovill --outdir $0 --R1 $0_1.fastq.gz --R2 $0_2.fastq.gz --gsize 5000000 \
	--mincov 10 --trim --namefmt $0_contig%05d; mv $0/contigs.fa $0/$0.contigs.fa; mv $0/contigs.gfa $0/$0.contigs.gfa; \
	mv $0/shovill.corrections $0/$0.shovill.corrections; mv $0/shovill.log $0/$0.shovill.log; mv $0/spades.fasta $0/$0.spades.fasta'
	
	cat readMetrics2.tsv | awk '{if($2<10 && $2>=4) print $0}' | xargs -P 1 -n 2 sh -c 'shovill --outdir $0 --R1 $0_1.fastq.gz --R2 $0_2.fastq.gz --gsize 5000000 \
	--mincov $1 --trim --namefmt $0_contig%05d; mv $0/contigs.fa $0/$0.contigs.fa; mv $0/contigs.gfa $0/$0.contigs.gfa; \
	mv $0/shovill.corrections $0/$0.shovill.corrections; mv $0/shovill.log $0/$0.shovill.log; mv $0/spades.fasta $0/$0.spades.fasta'

	cat readMetrics2.tsv | awk '{if($2<4) print $0}' > lowCoverageIsolates.vibrio."$time_stamp".tsv

	#move files; cleanup
	mkdir graphs; mv */*.contigs.gfa ./graphs
	mkdir shovill-logs; mv */*.shovill.log ./shovill-logs; mv */*.shovill.corrections ./shovill-logs
	mkdir assemblies; mv */*.contigs.fa ./assemblies; 
	mkdir assemblies-raw; mv */*.spades.fasta ./assemblies-raw
	cut -f 1 readMetrics2.tsv | xargs -P 1 -n 1 sh -c 'rmdir $0'
	cd assemblies
	
	#configure environment to run staramr
	module purge
	export PATH=/scicomp/home/lly3/miniconda3/bin:$PATH
	source activate staramr

	#run staramr 
	staramr search --pid-threshold 90 --percent-length-overlap-resfinder 50 --exclude-genes-file /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/NARST_research/scripts/databases/campy_genes_to_exclude.tsv -o out.staramr *.fa
	
	conda deactivate

	#clean up output files
	cd out.staramr
	cat resfinder.tsv | cut -f 1-2 |  sed -e 's/\.contigs//' | sed -e 's/\ //' | sed '1d' > import.resfinder.tsv
	cat summary.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($2 == "None") print $1,$2 }' | sed -e 's/\.contigs//'  >> import.resfinder.tsv 
	cd ..
	
	#run plasmid finder with abricate, formats import sheet
	module purge
   	unset PERL5LIB
	export PATH=/scicomp/home/lly3/miniconda3/bin:$PATH
	mkdir out.plasmids
	ls *.fa | sed -e 's/\.fa//' | xargs -I one sh -c 'abricate --db KTPF --threads 12 --minid 90 --mincov 60 --noheader one.fa  > out.plasmids/one.pf.tsv'
	cd out.plasmids
	abricate --summary *.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | awk 'BEGIN {FS = "\t"} {OFS="\t"}{if ($2 == "0") {print $1,"None"}}' > import.plasmidfinder.tsv 
	cat *.contigs.pf.tsv | sed -e 's/\.contigs\.fa//' | tee plasmidfinder.tsv | awk 'BEGIN {FS = "\t"}{OFS="\t"}{print $1,$5}' >> import.plasmidfinder.tsv
	ls *.contigs.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | xargs -I one sh -c 'cat one.contigs.pf.tsv | sed -e 's/\.contigs\.fa//' | echo "$(cut -f 1 | head -n 1)" \
	"$(cut -f 5 one.contigs.pf.tsv | paste -d, -s)">> plasmidsummary.tsv'
	abricate --summary *.pf.tsv | sed -e 's/\.contigs\.pf\.tsv//' | awk 'BEGIN {FS = "\t"}{OFS="\t"} {if ($2 == "0") {print $1,"None"}}' >> plasmidsummary.tsv
	sed -i '/^ *$/d' plasmidsummary.tsv; 
	mkdir hits; mv *.contigs.pf.tsv ./hits; cd ..
	cut -f 1-2 out.staramr/summary.tsv | sed '1d' | sed -e 's/\.contigs//' | sed -e 's/\ //g' > out.staramr/summary2.tsv
	sort -k 1 -o out.staramr/summary2.tsv out.staramr/summary2.tsv; sort -k 1 -o out.plasmids/plasmidsummary.tsv out.plasmids/plasmidsummary.tsv; 

	
	#run ariba to detect qrdr 
	module purge
	module load ariba/2.12.0
	cd .. ; mkdir mutational
	ls *_1.fastq.gz | sed -e 's/_1\.fastq\.gz//' | xargs -I one sh -c 'ariba run /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/NARST_research/scripts/databases/vcpoints/vcpoints one_1.fastq.gz one_2.fastq.gz mutational/one.mut; mv mutational/one.mut/report.tsv mutational/one.mut/one.report.tsv' 
	cd mutational
	ls | sed -e 's/\.mut//' | while read id 
	do 
		cat $id.mut/$id.report.tsv | awk -v id="$id" 'BEGIN {OFS="\t"}{if ($18 == "1") {print id,$1"("$19")"}}' | sed -e 's/(./(/' | sed -e 's/.)/)/' >> import.mutational.tsv
		cat $id.mut/$id.report.tsv | awk -v id="$id" 'BEGIN {OFS="\t"}{if ($18 == "1") {print id,$1"("$19")"}}' >> $id.summary.mutational.tsv
		
	done 
	
	for g in *.summary.mutational.tsv

        do
                (if [ -s "$g" ]
                then
                        g1="${g%.}";
                        cat "$g" | echo -e "${g1%.summary.mutational.tsv}\t""$(cut -f 2 | paste -d, -s)" >> all.summary.mutational.tsv
                else
                        g1="${g%.}";
                        cd ${g1%.summary.mutational.tsv}.mut; [ -f log.clusters.gz ] && { cd .. ; echo -e "${g1%.summary.mutational.tsv}\tNone" >> \
			all.summary.mutational.tsv; } || { cd ..; echo -e "${g1%.summary.mutational.tsv}\tError" >> all.summary.mutational.tsv;
			echo -e "${g1%.summary.mutational.tsv}\tError" >> import.mutational.tsv; }
                fi)
        done
	
	cd ../assemblies; 
	join out.staramr/summary2.tsv ../mutational/all.summary.mutational.tsv > determinants.tsv
	join determinants.tsv out.plasmids/plasmidsummary.tsv | sed -e 's/\ /\t/g' > final.summary.tsv
	
	mv final.summary.tsv ../final.summary.vibrio."$time_stamp".tsv; mv out.plasmids/import.plasmidfinder.tsv ../import.plasmidfinder.vibrio."$time_stamp".tsv; 
	mv out.staramr/import.resfinder.tsv ../import.resfinder.vibrio."$time_stamp".tsv; 
	mv ../mutational/import.mutational.tsv ../import.pointfinder.vibrio."$time_stamp".tsv	
	
else 

	echo "organism not recognized"

fi

