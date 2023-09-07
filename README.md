# MLseq_evaluation

This project contains the code to my BSc thesis.

## Summary
RNA-Sequencing (RNA-Seq) is a method to acquire information on the transcriptome of a cell or tissue
at a specific time point or condition. Due to its high resolution and breadth, it could be a promising tool
for revealing the genetic mechanisms behind many diseases. Moreover, it has the potential to discover
biomarkers that could be used to create new cures. However, since RNA-Seq data typically consist of
hundreds or thousands of transcripts, the analysis remains a key challenge, and approved methods often
fail to extract the biomarkers most important for the condition. Recently a tool, MLSeq [7], has been
proposed to overcome this challenge. It offers statistical learning models that can be trained to classify the
data. Some models, so-called sparse classifiers, rely on algorithms selecting a subset of genes for improved
training. I hypothesized that extracting this subset may be a way of discovering potential biomarkers, as
proposed by the authors.
I tested this hypothesis by comparing the approved method of DESeq2 to the novel method of MLSeq
on three datasets, two of the datasets used in the original MLSeq evaluation and one in-house dataset.
Thus, I aimed to evaluate which package is more suitable for biomarker discovery in terms of results,
comprehensiveness, and reproducibility.
I found that while DESeq2 performed well throughout the datasets, MLSeq showed a different performance
for each. While it was able to predict biomarkers well in the first and third datasets, it failed to select a
sufficiently small number of genes in the second. Since MLSeq also shows further limitations, I concluded
that MLSeq does not outperform DESeq2 in the analysis of RNA-Seq data. Instead, the results indicate
the best selection of biomarkers can be achieved by using both methods, ideally determining differentially
expressed genes with DESeq2 first and then performing additional gene selection with MLSeq.
