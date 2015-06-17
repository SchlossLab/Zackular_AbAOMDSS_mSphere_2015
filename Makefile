################################################################################
#
#	Part 1: Get the reference files
#
#	Here we give instructions on how to get the necessary reference files that
#	are used throughout the rest of the analysis. These are used to calculate
#	error rates, generate alignments, and classify sequences. We will pull down
#	the mock community reference (HMP_MOCK.fasta), the silva reference alignment
#	(silva.bacteria.align), and the RDP training set data (trainset9_032012).
#	Finally, we use the HMP_MOCK.align to get the alignment coordinates for the
#	V3-V4, V4, and V4-V5 data. These data will be stored in the data/references/
#	folder.
#
#	The targets in this part are all used as dependencies in other rules
#
################################################################################

#Location of reference files
REFS = data/references/

#get the silva reference alignment
$(REFS)silva.bacteria.align :
	wget -N -P $(REFS) http://www.mothur.org/w/images/2/27/Silva.nr_v119.tgz; \
	tar xvzf $(REFS)Silva.nr_v119.tgz -C $(REFS);
	mothur "#get.lineage(fasta=$(REFS)silva.nr_v119.align, taxonomy=$(REFS)silva.nr_v119.tax, taxon=Bacteria)";
	mv $(REFS)silva.nr_v119.pick.align $(REFS)silva.bacteria.align; \
	rm $(REFS)README.html; \
	rm $(REFS)README.Rmd; \
	rm $(REFS)silva.nr_v119.*

#get the v4 region of the alignment
$(REFS)silva.v4.align : $(REFS)silva.bacteria.align
	mothur "#pcr.seqs(fasta=$(REFS)silva.bacteria.align, start=11894, end=25319, keepdots=F, processors=8);\
			unique.seqs(fasta=current);"; \
	mv $(REFS)silva.bacteria.pcr.unique.align $(REFS)silva.v4.align; \
	rm $(REFS)silva.bacteria.pcr.*

#get the rdp training set data
$(REFS)trainset10_082014.pds.tax $(REFS)trainset10_082014.pds.fasta :
	wget -N -P $(REFS) http://www.mothur.org/w/images/2/24/Trainset10_082014.pds.tgz; \
	tar xvzf $(REFS)Trainset10_082014.pds.tgz -C $(REFS);\
	mv $(REFS)trainset10_082014.pds/trainset10_082014.* $(REFS);\
	rm -rf $(REFS)trainset10_082014.pds

#get the V4 region of the RDP training set
$(REFS)trainset10_082014.v4.tax $(REFS)trainset10_082014.v4.fasta : \
						$(REFS)trainset10_082014.pds.tax \
						$(REFS)trainset10_082014.pds.fasta \
						$(REFS)silva.v4.align
	mothur "#align.seqs(fasta=$(REFS)trainset10_082014.pds.fasta, reference=$(REFS)silva.v4.align, processors=8);\
		screen.seqs(fasta=current, taxonomy=$(REFS)trainset10_082014.pds.tax, start=1968, end=11550);\
		degap.seqs(fasta=current)"; \
	mv $(REFS)trainset10_082014.pds.good.ng.fasta $(REFS)trainset10_082014.v4.fasta; \
	mv $(REFS)trainset10_082014.pds.good.tax $(REFS)trainset10_082014.v4.tax;\
	rm data/references/trainset10_082014.pds.align*;\
	rm data/references/trainset10_082014.pds.bad.accnos;\
	rm data/references/trainset10_082014.pds.flip.accnos;

$(REFS)HMP_MOCK.fasta :
	wget --no-check-certificate -N -P $(REFS) https://raw.githubusercontent.com/SchlossLab/Kozich_MiSeqSOP_AEM_2013/master/data/references/HMP_MOCK.fasta

#align the mock community reference sequeces
$(REFS)HMP_MOCK.v4.fasta : $(REFS)HMP_MOCK.fasta $(REFS)silva.v4.align
	mothur "#align.seqs(fasta=$(REFS)HMP_MOCK.fasta, reference=$(REFS)silva.v4.align);\
			degap.seqs()";\
	mv $(REFS)HMP_MOCK.ng.fasta $(REFS)HMP_MOCK.v4.fasta;\
	rm $(REFS)HMP_MOCK.align;\
	rm $(REFS)HMP_MOCK.align.report;\
	rm $(REFS)HMP_MOCK.flip.accnos


################################################################################
#
#	Part 2: Run data through mothur
#
#
################################################################################

data/raw/get_data : code/get_fastqs.sh data/raw/ab_aomdss.files
	bash code/get_fastqs.sh data/raw/ab_aomdss.files;\
	touch data/raw/get_data

BASIC_STEM = data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster



# here we go from the raw fastq files and the files file to generate a fasta,
# taxonomy, and count_table file that has had the chimeras removed as well as
# any non bacterial sequences
$(BASIC_STEM).denovo.uchime.pick.pick.count_table $(BASIC_STEM).pick.pick.fasta $(BASIC_STEM).pick.v4.wang.pick.taxonomy : code/get_good_seqs.batch\
										data/raw/get_data\
										data/references/silva.v4.align\
										data/references/trainset10_082014.v4.fasta\
										data/references/trainset10_082014.v4.tax
	mothur code/get_good_seqs.batch;\
	rm data/process/*.map



# here we go from the good sequences and generate a shared file and a
# cons.taxonomy file based on OTU data
$(BASIC_STEM).pick.pick.pick.an.unique_list.shared $(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.cons.taxonomy : code/get_shared_otus.batch\
										$(BASIC_STEM).denovo.uchime.pick.pick.count_table\
										$(BASIC_STEM).pick.pick.fasta\
										$(BASIC_STEM).pick.v4.wang.pick.taxonomy
	mothur code/get_shared_otus.batch;\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table;\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta;\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.taxonomy;



# here we go from the good sequences and generate a shared file and a
# cons.taxonomy file based on phylum-level data
$(BASIC_STEM).pick.v4.wang.pick.pick.tx.5.cons.taxonomy $(BASIC_STEM).pick.v4.wang.pick.pick.tx.shared : code/get_shared_phyla.batch\
										$(BASIC_STEM).denovo.uchime.pick.pick.count_table\
										$(BASIC_STEM).pick.pick.fasta\
										$(BASIC_STEM).pick.v4.wang.pick.taxonomy
	mothur code/get_shared_phyla.batch;\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table;\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta;\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.taxonomy;\
	rm data/process/*.tx.*rabund;


# now we want to get the sequencing error as seen in the mock community samples
$(BASIC_STEM).pick.pick.pick.error.summary : code/get_error.batch\
										$(BASIC_STEM).denovo.uchime.pick.pick.count_table\
										$(BASIC_STEM).pick.pick.fasta\
										$(REFS)HMP_MOCK.v4.fasta
	mothur code/get_error.batch


# rarefy the number of reads to 2500 sequences per library for the alpha and beta diversity analyses
$(BASIC_STEM).pick.pick.pick.an.unique_list.groups.ave-std.summary $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist $(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared : $(BASIC_STEM).pick.pick.pick.an.unique_list.shared
	mothur "#set.dir(seed=19760620);dist.shared(shared=$^, calc=thetayc, subsample=2500, iters=100); summary.single(shared=$^, subsample=2500, calc=nseqs-sobs-shannon-invsimpson, iters=100); sub.sample(shared=$^, size=2500)";\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.summary;\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.dist;\
	rm data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.std.dist


# pull out the day 0 submatrix from above and the design file
$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.dist $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.design : code/get_day0_matrix_design.R $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist
	R -e "source('code/get_day0_matrix_design.R')"


# build NMDS plot for day 0 submatrix
$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.nmds.stress $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.nmds.axes : $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.dist
	mothur "#set.dir(seed=19760620);nmds(phylip=$<, maxdim=2)";\
	rm data/process/*day0.nmds.iters


# run AMOVA on subsetted data
$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.amova : $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.dist $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.design
	mothur "#set.dir(seed=19760620);amova(phylip=data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.dist, design=data/process/ab_aomdss.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.design, iters=100000)"


################################################################################
#
#	Part 3: Write the paper
#
#
################################################################################

$(FIGURES) = results/figures
$(FIGURES)/figure_1.pdf : code/build_figure_1.R\
							code/rf_baseline_analysis.R\
 							code/aomdss_timeline.R\
							$(BASIC_STEM).pick.v4.wang.pick.pick.tx.shared\
							$(BASIC_STEM).pick.v4.wang.pick.pick.tx.5.cons.taxonomy\ $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.nmds.axes\
							results/pictures/tumor_images.png\
							data/process/tumor_counts.tsv
	source("code/build_figure_1.R")


$(FIGURES)/figure_S1.pdf $(FIGURES)/figure_2.pdf $(FIGURES)/figure_3.pdf data/process/baseline_model.Rdata : code/rf_baseline_analysis.R\
							$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared\
							$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.cons.taxonomy\
							data/process/tumor_counts.tsv
	source("code/run_baseline_model.R")


$(FIGURES)/figure_S2.pdf $(FIGURES)/figure_S3.pdf $(FIGURES)/figure_S4.pdf data/process/final_model.Rdata : code/rf_final_analysis.R\
								$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared\
								$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.cons.taxonomy\
								data/process/tumor_counts.tsv
	source("code/run_final_model.R")


$(FIGURES)/figure_4.pdf data/process/distance_tumor_correlation.Rdata : code/build_figure_4.R\
				code/make_timeline_plot.R\
				code/rf_baseline_analysis.R\
				$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist\
				$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared
	source("code/build_figure_4.R")

$(FIGURES)/figure_5.pdf : code/build_figure_5.R data/process/tumor_counts.tsv\
							data/process/intervention_tumor_counts.tsv\
							results/pictures/tumor_images_intervention.png
	source("code/build_figure_5.R")



write.paper :   data/process/baseline_model.Rdata\
				data/process/final_model.Rdata\
				data/process/distance_tumor_correlation.Rdata\
				$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.amova\
				$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.nmds.stress\
				$(FIGURES)/figure_1.pdf\
				$(FIGURES)/figure_2.pdf\
				$(FIGURES)/figure_3.pdf\
				$(FIGURES)/figure_4.pdf\
				$(FIGURES)/figure_5.pdf\
				$(FIGURES)/figure_S1.pdf\
				$(FIGURES)/figure_S2.pdf\
				$(FIGURES)/figure_S3.pdf\
				$(FIGURES)/figure_S4.pdf


				$(BASIC_STEM).pick.pick.pick.error.summary\
				                $(BASIC_STEM).pick.pick.pick.an.unique_list.groups.ave-std.summary\
				                $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist\
				                $(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared\
				                $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.amova\

	R -e "library(knitr);knit2html('Zackular_AbAOMDSS_SciReports_2015.Rmd', 'Zackular_AbAOMDSS_SciReports_2015.html')"; \
	pandoc -f markdown -t docx Zackular_AbAOMDSS_SciReports_2015.md -o Zackular_AbAOMDSS_SciReports_2015.docx






################################################################################
#
#	Part 4: Notebook entries
#
#
################################################################################

doc/notebook/2015_02_02-random_forest.html : $(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared\
											$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.cons.taxonomy\
											code/rf_baseline_analysis.R
	R -e "library(knitr);knit2html('doc/notebook/2015_02_02-random_forest.Rmd')";\
	mv 2015_02_02-random_forest.* doc/notebook/

notebooks : doc/notebook/2015_02_02-random_forest.html



################################################################################
#
#	Part 5: Everything
#
#
################################################################################

all : notebooks\
		$(BASIC_STEM).pick.pick.pick.error.summary\
		$(BASIC_STEM).pick.pick.pick.an.unique_list.groups.ave-std.summary\
		$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist\
		$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared\
		$(BASIC_STEM).pick.v4.wang.pick.pick.tx.5.cons.taxonomy\
		$(BASIC_STEM).pick.v4.wang.pick.pick.tx.shared\
		$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.amova\
		$(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.day0.nmds.stress
