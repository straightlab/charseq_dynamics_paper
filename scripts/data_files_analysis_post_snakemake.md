This notebook logs which files are produced or used by analysis scripts (post snakemake preprocessing) for reproducibility. This notebook can be used to convert file paths in analysis notebooks into appropriate paths for a user specific system.

# Data files used or produced by `scripts/prepare_charseq_data_for_python.ipynb`
 
## Support files
```bash
rsync -avh -P /oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/ENSG_1_X.csv /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/ENSG_1_X.csv

rsync -avh -P /oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/resources/chrNameLength_1-22X.txt /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/chrNameLength_1-22X.txt

rsync -avh -P /home/groups/astraigh/differentiation_charseq_paper/notebooks/charles/13_intergenic/data/charOnly/strngtie.1to22X.NOdup.csv /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/strngtie.1to22X.NOdup.csv

```

## ChARseq data in stored in chartable objects
These data files are all what is necessary to produce ChAR-seq maps as in Fig. 1

```bash
#pickled contact maps for annotated genes
parallel -j8 'rsync -avh -P /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/{1}/data/{2}/pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/dna.Q255Q40.{4}.chartable.compress.10.pickle /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/charseq_data_as_chartables/{3}_dna.Q255Q40.{4}.chartable.compress.10.pickle' ::: novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun :::+ TFONM2_ES TGNQ5_DE TMF2_Rep2B_ES TMF3_Rep2A_DE :::+ ES_rep1 DE_rep1 ES_rep2 DE_rep2 ::: exons introns all

#same for UTLs
parallel -j8 'rsync -avh -P /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/{1}/data/{2}/pairs/strngtieCharOnly_hg38/exons/filtered/DNAq15-blacklisted_RNAunq/dna.Q255Q40.all.chartable.compress.10.pickle /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/charseq_data_as_chartables/{3}_dna.Q255Q40.UTLs.chartable.compress.10.pickle' ::: novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun :::+ TFONM2_ES TGNQ5_DE TMF2_Rep2B_ES TMF3_Rep2A_DE :::+ ES_rep1 DE_rep1 ES_rep2 DE_rep2
```

## Location of RNA-side 3' "source" locus for RNAs binding in CIS
```bash
#pickled chartables containing RNA coordinates of the RNA-DNA contact point for each RNA binding to CIS chrosomome.
parallel -j8 'rsync -avh -P /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/{1}/data/{2}/pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/rna.Q255Q40.cis.{4}.chartable.1.pickle /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/charseq_data_as_chartables/{3}_rna.Q255Q40.cis.{4}.chartable.1.pickle' ::: novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun :::+ TFONM2_ES TGNQ5_DE TMF2_Rep2B_ES TMF3_Rep2A_DE :::+ ES_rep1 DE_rep1 ES_rep2 DE_rep2 ::: exons introns all

#same for UTLs
parallel -j8 'rsync -avh -P /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/{1}/data/{2}/pairs/strngtieCharOnly_hg38/exons/filtered/DNAq15-blacklisted_RNAunq/rna.Q255Q40.cis.all.chartable.1.pickle /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/charseq_data_as_chartables/{3}_rna.Q255Q40.cis.UTLs.chartable.1.pickle' ::: novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun :::+ TFONM2_ES TGNQ5_DE TMF2_Rep2B_ES TMF3_Rep2A_DE :::+ ES_rep1 DE_rep1 ES_rep2 DE_rep2
```


## Datafiles for modeling background using TRANS mRNAs predictive models 
```bash
#delocalization score files used for masking trans-delocalized RNAs from models
rsync -avh -P /home/groups/astraigh/differentiation_charseq_paper/notebooks/charles/17_tscores/data/tscoresAllSamples.exons.tab /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/tscoresAllSamples.exons.tab

#pickled chartables containing charseqability background models
parallel --dry-run -j8 'rsync -avh -P /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/{1}/data/{2}/pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/bkg.Q255Q40.chartable.compress.10.pickle /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/{3}_bkg.Q255Q40.chartable.compress.10.pickle' ::: novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun :::+ TFONM2_ES TGNQ5_DE TMF2_Rep2B_ES TMF3_Rep2A_DE :::+ ES_rep1 DE_rep1 ES_rep2 DE_rep2

```

# Data files used or produced by `scripts/compute_cis_background_model.ipynb`
```bash

rsync -avh -P  /home/groups/astraigh/differentiation_charseq_paper/notebooks/charles/17_tscores/data/tscoresAllSamples.introns.tab /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/tscoresAllSamples.introns.tab

rsync -avh -P /home/groups/astraigh/differentiation_charseq_paper/notebooks/charles/02_broadBinders/data/tscores_cep002bb.csv /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/tscores_cep002bb.csv



#pickled chartables containing RNA travel properties by RNA
parallel --dry-run -j8 'rsync -avh -P /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/{1}/data/{2}/pairs/gencondeV29_hg38/{4}/analysis/01_flight/fly.stranded.pickle /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/{3}_{4}_fly.stranded.pickled' ::: novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun :::+ TFONM2_ES TGNQ5_DE TMF2_Rep2B_ES TMF3_Rep2A_DE :::+ ES_rep1 DE_rep1 ES_rep2 DE_rep2 ::: exons introns


#pickled CIS model files 
parallel -j8 'rsync -avh -P /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/{1}/data/{2}/pairs/gencondeV29_hg38/{4}/analysis/01_flight/cisprofiles.pickle /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/{3}_{4}_cisprofiles.pickle' ::: novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun :::+ TFONM2_ES TGNQ5_DE TMF2_Rep2B_ES TMF3_Rep2A_DE :::+ ES_rep1 DE_rep1 ES_rep2 DE_rep2 ::: exons introns

parallel -j8 'rsync -avh -P /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/{1}/data/{2}/pairs/gencondeV29_hg38/{4}/analysis/01_flight/cismodel.pickle /oak/stanford/groups/astraigh/differentiation_paper_data/for_zenodo/predictive_models_data_files/{3}_{4}_cismodel.pickle.pickle' ::: novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun novchar2/CL_2020-05-06_rerun :::+ TFONM2_ES TGNQ5_DE TMF2_Rep2B_ES TMF3_Rep2A_DE :::+ ES_rep1 DE_rep1 ES_rep2 DE_rep2 ::: exons introns

```

