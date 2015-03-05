Manipulation of the Gut Microbiome Reveals Role for Microbial Community Structure in Colon Tumorigenesis
=======

There is growing evidence that individuals with colonic adenomas and carcinomas
harbor a distinct microbiota. Alterations in the gut microbiota may allow the
outgrowth of bacterial populations that induce genomic mutations or exacerbate
tumor-promoting inflammation. In addition, it is likely that the loss of key
bacterial populations may result in the loss of protective functions that are
normally provided by a healthy microbiota. Using an inflammation-based murine
model of colorectal cancer, we explored the role of the gut microbiota in colon
tumorigenesis. We observed that the number of tumors that developed in the model
could be altered by perturbing the microbiota with different combinations of
antibiotics. One particular antibiotic combination, specifically, metronidazole
and streptomycin, was superior to other antibiotic combinations tested and was
sufficient to completely suppress tumor development that was not due to changes
in bacterial load. Using the random forest machine learning algorithm we were
able to create a regression model that predicted the number of tumors that
developed over the course of 73 days based on the composition of the microbiota
on the first day. Finally, the timing of antibiotic treatment was an important
determinant of tumor outcome as colon tumorigenesis was arrested with the use of
antibiotics during the inflammation period of the murine model. These results
suggest that it is possible to predict colon tumorigenesis based on knowledge of
the microbiota and that altering the structure of the gut microbiota can alter
the course of colorectal cancer.


Overview
--------

    project
    |- README          # the top level description of content
    |
    |- doc/            # documentation for the study
    |  |- notebook/    # preliminary analyses (dead branches of analysis)
    |  +- paper/       # manuscript(s), whether generated or not
    |
    |- data            # raw and primary data, are not changed once created
    |  |- references/  # reference files to be used in analysis
    |  |- raw/         # raw data, will not be altered
    |  +- process/     # cleaned data, will not be altered once created
    |
    |- code/           # any programmatic code
    |- results         # all output from workflows and analyses
    |  |- tables/      # text version of tables to be rendered with kable in R
    |  |- figures/     # graphs, likely designated for manuscript figures
    |  +- pictures/    # diagrams, images, and other non-graph graphics
    |
    |- scratch/        # temporary files that can be safely deleted or lost
    |
    |- Zackular_AbAOMDSS_GutMicrobes_2015.Rmd  # executable Rmarkdown for this study, if applicable
    |- Zackular_AbAOMDSS_GutMicrobes_2015.md   # Markdown (GitHub) version of the *Rmd file
    |- Zackular_AbAOMDSS_GutMicrobes_2015.html # HTML version of *.Rmd file
    |
    +- Makefile        # executable Makefile for this study
