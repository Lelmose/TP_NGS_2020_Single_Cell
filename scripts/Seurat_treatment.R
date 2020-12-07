library(SeuratWrappers)
library(velocyto.R)

# Read the loom file, the resulting object contains 4 tabs, spliced, unspliced, ambiguous and spanning
tooth <- ReadVelocity('/ifb/data/mydatalocal/loom_files/onefilepercell_SRR11201752_and_others_DFMGU.loom')


# Normalization step
norm_spliced <- NormalizeData(tooth$spliced)
norm_unspliced <- NormalizeData(tooth$unspliced)
norm_ambiguous <- NormalizeData(tooth$ambiguous)
