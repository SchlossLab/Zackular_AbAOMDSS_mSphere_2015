Manipulation of the Gut Microbiome Reveals Role for Microbial Community Structure in Colon Tumorigenesis
=======

There is growing evidence that individuals with colonic adenomas and carcinomas harbor a distinct microbiota. Alterations to the gut microbiota may allow the outgrowth of bacterial populations that induce genomic mutations or exacerbate tumor-promoting inflammation. In addition, it is likely that the loss of key bacterial populations may result in the loss of protective functions that are normally provided by the microbiota. We explored the role of the gut microbiota in colon tumorigenesis using an inflammation-based murine model. We observed that perturbing the microbiota with different combinations of antibiotics did not change the bacterial load but reduced the number of tumors at the end of the model. Using the random forest machine learning algorithm we successfully modeled the number of tumors that developed over the course of the model based on the composition of the microbiota at the beginning. The timing of antibiotic treatment was an important determinant of tumor outcome as colon tumorigenesis was arrested with the use of antibiotics during the inflammation period of the murine model. Together, these results indicate that it is possible to predict colon tumorigenesis based on the composition of the microbiota and that altering the gut microbiota can alter the course of tumorigenesis.


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
    |- Zackular_AbAOMDSS_mSphere_2015.Rmd  # executable Rmarkdown for this study, if applicable
    |- Zackular_AbAOMDSS_mSphere_2015.md   # Markdown (GitHub) version of the *Rmd file
    |- Zackular_AbAOMDSS_mSphere_2015.html # HTML version of *.Rmd file
    |- Zackular_AbAOMDSS_mSphere_2015.docx # Word version of *.Rmd file
    |
    +- Makefile        # executable Makefile for this study
