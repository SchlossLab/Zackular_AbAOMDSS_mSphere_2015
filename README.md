Manipulation of the Gut Microbiome Reveals Role for Microbial Community Structure in Colon Tumorigenesis
=======

There is growing evidence that individuals with colonic adenomas and carcinomas
harbor a distinct microbiota. These alterations may allow the outgrowth of
populations that induce mutations or exacerbate inflammation. In addition, it is
likely that the loss of key populations may result in the loss of protective
functions that are provided for by a healthy microbiota. Using an
inflammation-based murine model of colorectal cancer we explored the
host-microbiota relationship to better understand the role of various
populations through the process of tumorigenesis. By perturbing the microbiota
with mixtures of antibiotics that targeted distinct groups of bacteria we
observed that it was possible to predict the number of tumors that the animals
would harbor by the end of the model. It was apparent that distinct microbiota
could lead to similar numbers of tumors and that variation in the composition of
the microbiota could also lead to wide variation in the number of colonic tumors
that formed. Furthermore, without altering the number of bacteria in the colon,
we were able to fully suppress tumor formation using a combination of
metronidazole and streptomycin. Finally, by altering when the antibiotics were
given to the model we showed that the role of the microbiota in tumorigenesis is
most pronounced during the period of inflammation rather than in the processing
of the mutagen. These results suggest that altering the structure and function
of the gut microbiota can arrest the course of colorectal cancer.


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
