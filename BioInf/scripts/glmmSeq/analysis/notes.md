Notes on Process

GOseq -- when changing the gene2cat variable in the goseq() function from a list of GO terms in the frontloaded data to a list of GO terms in all 9011 genes, the data makes much more sense to me. 



Smaller list --> Whole List ()


Normal p-values

  For Frontloading (2631/9011 genes): We go from 6493 (out of 12,434 (unique out of duplicated 217,506)) GO terms with p value < 0.05 to 621 (out of 17,741 (unique out of duplicated 705,483)).

  For Origin DEGs (840/9011 genes): We go from 3732 (out of 6,767 (unique out of duplicated 46,011)) GO terms with p value < 0.05 to 219 (out of 17,741 (unique out of duplicated 705,483)).



If we use adjusted p-values, we have 0 enriched GO terms for the frontloaded data set and for the DEGs

  For Frontloading (2631/9011 genes): We go from 5618 (out of 12,434 (unique out of duplicated 217,506)) GO terms with adjusted p value < 0.05 to 0 (out of 17,741 (unique out of duplicated 705,483)).

  For Origin DEGs (840/9011 genes): We go from 3732 (out of 6,767 (unique out of duplicated 46,011)) GO terms with adjusted p value < 0.05 to 0 (out of 17,741 (unique out of duplicated 705,483)).

maybe we should remove the filter(numInCat>10) filter??

