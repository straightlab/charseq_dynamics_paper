Visualization : chromatin associated RNAs basic statistics
================

Load info

``` r
allgenes <- read_parquet(here('../rdana/genes/data-output/allgenes_final.parquet'))
exprdata <- read_parquet(here('../rdana/carnas/data-output/expression.rnacharIndep.exonsScaling.Q255Q40.with_unannotated.parquet'))
```

# RNA census

## By annotation type (exons, introns, intergenic)

### exons vs introns

``` r
type_tallies <- exprdata %>%
  
dplyr::select(GeneID, annotation_type, starts_with("FPM")) %>%
  dplyr::mutate(annotation_type = factor(annotation_type, levels = c("introns","exons","intergenic"))) %>%
pivot_longer(cols=-c(GeneID, annotation_type), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM") %>%
group_by(cell, sequencing, annotation_type) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup()
```

    ## `summarise()` has grouped output by 'cell', 'sequencing'. You can override
    ## using the `.groups` argument.

``` r
p<- type_tallies %>%
  group_by(cell, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  dplyr::select(-cell) %>%
  group_by(sequencing, annotation_type) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() %>%
  ggplot(aes(x=sequencing, y=FPM, fill=annotation_type))+
  geom_col(stats="identity") +
  scale_fill_manual(values=c('#8856a7','#43a2ca','#bdbdbd'))+
  #ggsci::scale_fill_nejm()+
  #1f78b4
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.avg.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.1 * nrow(p$data), panel.height=1.5)
```

    ## [1] "fig.width=2.5, fig.height=2.1"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
p<- type_tallies %>%
  mutate(cell= factor(cell, levels = c("ES","DE"))) %>%
  group_by(cell, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  ggplot(aes(x=cell, y=FPM, fill=annotation_type))+
  geom_col(stats="identity") +
  scale_fill_manual(values=c('#8856a7','#43a2ca','#bdbdbd'))+
  #ggsci::scale_fill_nejm()+
  #1f78b4
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~sequencing)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.bycell.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.1/2 * nrow(p$data), panel.height=1.5)
```

    ## [1] "fig.width=3.2, fig.height=2.3"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
type_tallies %>% group_by(cell, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  dplyr::select(-cell) %>%
  group_by(sequencing, annotation_type) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() 
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `sequencing`

    ## # A tibble: 6 × 3
    ##   sequencing annotation_type   FPM
    ##   <chr>      <fct>           <dbl>
    ## 1 char       introns         72.1 
    ## 2 char       exons           15.2 
    ## 3 char       intergenic      12.7 
    ## 4 rna        introns         18.4 
    ## 5 rna        exons           74.2 
    ## 6 rna        intergenic       7.42

## By rna type

### Exons, introns

``` r
types_to_plot <- c("mRNA","lncRNA","ncRNA")
type_tallies <- exprdata %>%
        dplyr::filter(annotation_type %in% c('exons','introns'), rna_type %in% types_to_plot) %>%

dplyr::select(GeneID, annotation_type, rna_type, starts_with("FPM")) %>%
  dplyr::mutate(rna_type = factor(rna_type, levels = c("mRNA","lncRNA","ncRNA"))) %>%
pivot_longer(cols=-c(GeneID, annotation_type, rna_type), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM") %>%
group_by(cell, sequencing, annotation_type, rna_type) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup()
```

    ## `summarise()` has grouped output by 'cell', 'sequencing', 'annotation_type'.
    ## You can override using the `.groups` argument.

``` r
p<- type_tallies %>%
    mutate(cell= factor(cell, levels = c("ES","DE"))) %>%
  group_by(cell, annotation_type, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  dplyr::select(-cell) %>%
  group_by(sequencing, annotation_type, rna_type) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() %>%
  ggplot(aes(x=sequencing, y=FPM, fill=rna_type))+
  geom_col(stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_manual(values=c(mRNA='#cccccc', lncRNA='#1f78b4', ncRNA='#b2df8a'))+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~annotation_type)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `annotation_type`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.avg.simple.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.7, panel.height=1.5)
```

    ## [1] "fig.width=3.3, fig.height=2.4"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
p<- type_tallies %>%
   mutate(cell= factor(cell, levels = c("ES","DE"))) %>%
  group_by(cell, annotation_type, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  ggplot(aes(x=cell, y=FPM, fill=rna_type))+
  geom_col(stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_manual(values=c(mRNA='#cccccc', lncRNA='#1f78b4', ncRNA='#b2df8a'))+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~annotation_type+sequencing)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `annotation_type`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.bycell.simple.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.7, panel.height=1.5)
```

    ## [1] "fig.width=3.3, fig.height=4.9"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### UNAs

``` r
type_tallies <- exprdata %>%
        dplyr::filter(annotation_type=='intergenic') %>%
        dplyr::select(GeneID, starts_with("FPM")) %>%
#   dplyr::mutate(rna_type = if_else(rna_type=='protein_coding', 'mRNA', rna_type)) %>%
#   dplyr::mutate(rna_type = factor(rna_type, levels = c("mRNA","lncRNA","ncRNA", "repeat","tRNAderived", "sRNAderived", "cre","geneproximal","antisense","intergenic"))) %>%
pivot_longer(cols=-c(GeneID), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM") %>%
  inner_join(allgenes %>%
          dplyr::select(GeneID, rna_type, rna_subtype), by=c('GeneID')) %>%
  
group_by(cell, sequencing, rna_type) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup() %>%
  dplyr::mutate(rna_type = factor(rna_type, levels = c("repeat","tRNAderived", "snRNAderived", "cre","readthrough","antisense","intergenic")))
```

    ## `summarise()` has grouped output by 'cell', 'sequencing'. You can override
    ## using the `.groups` argument.

``` r
p<- type_tallies %>%
  mutate(cell = factor(cell, levels = c("ES","DE"))) %>%
  group_by(cell, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  # dplyr::select(-cell) %>%
  # group_by(sequencing, rna_type) %>%
  # summarise_if(is.numeric, mean) %>%
  # ungroup() %>%
  ggplot(aes(x=sequencing, y=FPM, fill=rna_type))+
  geom_col(stats="identity") +
  scale_fill_brewer(palette = 'Paired')+
  #1f78b4
  #scale_fill_manual(values=c(mRNA='#cccccc', lncRNA='#1f78b4', ncRNA='#b2df8a'))+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~cell)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.bysequencing.simple.INTERGENIC.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.7, panel.height=1.5)
```

    ## [1] "fig.width=3.8, fig.height=2.4"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
p<- type_tallies %>%
  group_by(cell, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  mutate(cell=factor(cell, levels = c("ES","DE"))) %>%
  #dplyr::select(-cell) %>%
  #group_by(sequencing, rna_type) %>%
  #summarise_if(is.numeric, mean) %>%
  #ungroup() %>%
  ggplot(aes(x=cell, y=FPM, fill=rna_type))+
  geom_col(stats="identity") +
  scale_fill_brewer(palette = 'Paired')+
  #1f78b4
  #scale_fill_manual(values=c(mRNA='#cccccc', lncRNA='#1f78b4', ncRNA='#b2df8a'))+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~sequencing)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.bycell.simple.INTERGENIC.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.7, panel.height=1.5)
```

    ## [1] "fig.width=3.8, fig.height=2.3"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## By rna type

### Exons, introns

``` r
type_tallies <- exprdata %>%
        dplyr::filter(annotation_type %in% c('exons','introns'), rna_type=='ncRNA') %>%

dplyr::select(GeneID, annotation_type, rna_subtype, starts_with("FPM")) %>%
pivot_longer(cols=-c(GeneID, annotation_type, rna_subtype), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM") %>%
group_by(cell, sequencing, annotation_type, rna_subtype) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup()
```

    ## `summarise()` has grouped output by 'cell', 'sequencing', 'annotation_type'.
    ## You can override using the `.groups` argument.

``` r
p<- type_tallies %>%
  group_by(cell, annotation_type, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  dplyr::select(-cell) %>%
  group_by(sequencing, annotation_type, rna_subtype) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() %>%
  ggplot(aes(x=sequencing, y=FPM, fill=rna_subtype))+
  geom_col(stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_brewer(palette='Paired')+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~annotation_type)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `annotation_type`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.avg.ncRNAsDetails.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.7, panel.height=1.5)
```

    ## [1] "fig.width=3.6, fig.height=2.4"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
type_tallies <- exprdata %>%
        dplyr::filter(annotation_type %in% c('exons','introns')) %>%
          dplyr::filter(rna_type=='ncRNA') %>%
dplyr::select(GeneID, annotation_type, rna_subtype, starts_with("FPM")) %>%
pivot_longer(cols=-c(GeneID, annotation_type, rna_subtype), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM") %>%
group_by(cell, sequencing, annotation_type, rna_subtype) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup()
```

    ## `summarise()` has grouped output by 'cell', 'sequencing', 'annotation_type'.
    ## You can override using the `.groups` argument.

``` r
p<- type_tallies %>%
   mutate(cell= factor(cell, levels = c("ES","DE"))) %>%
  group_by(cell, annotation_type, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  ggplot(aes(x=cell, y=FPM, fill=rna_subtype))+
  geom_col(stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_brewer(palette='Paired')+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~annotation_type + sequencing)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `annotation_type`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.bycell.ncRNAsDetails.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.7, panel.height=1.5)
```

    ## [1] "fig.width=3.6, fig.height=4.9"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Exons, introns

``` r
type_tallies <- exprdata %>%
        dplyr::filter(annotation_type %in% c('exons','introns')) %>%
        left_join(allgenes %>% dplyr::select(GeneID, rna_type_gencode)) %>%
          dplyr::filter(rna_subtype=='ncRNA', rna_type=='ncRNA') %>%
dplyr::select(GeneID, annotation_type, rna_type_gencode , starts_with("FPM")) %>%
pivot_longer(cols=-c(GeneID, annotation_type, rna_type_gencode), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM") %>%
group_by(cell, sequencing, annotation_type, rna_type_gencode) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup()
```

    ## Joining, by = "GeneID"
    ## `summarise()` has grouped output by 'cell', 'sequencing', 'annotation_type'.
    ## You can override using the `.groups` argument.

``` r
p<- type_tallies %>%
  group_by(cell, annotation_type, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  dplyr::select(-cell) %>%
  group_by(sequencing, annotation_type, rna_type_gencode) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() %>%
  ggplot(aes(x=sequencing, y=FPM, fill=rna_type_gencode))+
  geom_col(stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_brewer(palette='Paired')+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~annotation_type)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `annotation_type`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.avg.ncRNAsDetails_more.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.7, panel.height=1.5)
```

    ## [1] "fig.width=2.6, fig.height=2.4"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### CREs

``` r
type_tallies <- exprdata %>%
        dplyr::filter(annotation_type=='intergenic') %>%
        dplyr::select(GeneID, starts_with("FPM")) %>%
#   dplyr::mutate(rna_type = if_else(rna_type=='protein_coding', 'mRNA', rna_type)) %>%
#   dplyr::mutate(rna_type = factor(rna_type, levels = c("mRNA","lncRNA","ncRNA", "repeat","tRNAderived", "sRNAderived", "cre","geneproximal","antisense","intergenic"))) %>%
pivot_longer(cols=-c(GeneID), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM") %>%
  inner_join(allgenes %>%
          dplyr::select(GeneID, rna_type, rna_subtype), by=c('GeneID')) %>%
  
group_by(cell, sequencing, rna_type, rna_subtype) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup() %>%
  #dplyr::mutate(rna_type = factor(rna_type, levels = c("repeat","tRNAderived", "snRNAderived", "cre","readthrough","antisense","intergenic"))) %>%
  dplyr::filter(rna_type=='cre') %>%
   mutate(rna_subtype = factor(rna_subtype, levels = c("dELS","pELS","PLS","DNase-H3K4me3","CTCF-only"))) 
```

    ## `summarise()` has grouped output by 'cell', 'sequencing', 'rna_type'. You can
    ## override using the `.groups` argument.

``` r
p<- type_tallies %>%
  group_by(cell, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  #dplyr::select(-cell) %>%
  group_by(sequencing, rna_subtype, cell) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() %>%
  mutate(cell=factor(cell, levels= c("ES","DE"))) %>%
  ggplot(aes(x=rna_subtype, y=FPM, fill=cell))+
  geom_col(stats="identity", position = position_dodge()) +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_brewer(palette='Paired')+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~sequencing)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.avg.INTERGENICdetails.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 1.2, panel.height=1.5)
```

    ## [1] "fig.width=3.8, fig.height=3.1"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
\### Repeats

``` r
type_tallies <- exprdata %>%
        dplyr::filter(annotation_type=='intergenic') %>%
        dplyr::select(GeneID, rna_type, rna_subtype, starts_with("FPM")) %>%
pivot_longer(cols=-c(GeneID, rna_type, rna_subtype), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM") %>%
group_by(cell, sequencing, rna_type, rna_subtype) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup() %>%
  dplyr::filter(rna_type=='repeat') %>%
group_by(cell, sequencing, rna_subtype) %>%
summarise(across(starts_with("FPM"), function(x) sum(x, na.rm=T))) %>%
ungroup() %>%
            dplyr::mutate(rna_subtype = forcats::fct_reorder(rna_subtype, FPM, .fun=mean, .desc=F)) %>%
  mutate(cell =factor(cell, levels = c("ES","DE")))
```

    ## `summarise()` has grouped output by 'cell', 'sequencing', 'rna_type'. You can
    ## override using the `.groups` argument.
    ## `summarise()` has grouped output by 'cell', 'sequencing'. You can override
    ## using the `.groups` argument.

``` r
p<- type_tallies %>%
  group_by(cell, sequencing) %>%
  mutate_if(is.numeric, ~./sum(., na.rm = T)*100) %>%
  ungroup() %>%
  ggplot(aes(x=cell, y=FPM, fill=rna_subtype))+
  geom_col(stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_brewer(palette='Paired')+
scale_y_continuous(expand = c(0,0)) + theme_publish()+

  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% FPM", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~sequencing)
```

    ## `mutate_if()` ignored the following grouping variables:
    ## • Columns `cell`, `sequencing`

    ## Warning: Ignoring unknown parameters: stats

``` r
fname='annotation_types.barplot.reptypes.pdf'
p_fixed<- prettysave(p, here('figures/carnas', fname), panel.width= 0.8, panel.height=1.5)
```

    ## [1] "fig.width=4, fig.height=2.3"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

# Diversity

## Prepare data

``` r
expr_genes <- exprdata %>%
        #dplyr::select(-type) %>%
        dplyr::filter(annotation_type %in% c('exons','introns')) %>%
dplyr::select(GeneID, annotation_type, rna_type, starts_with("FPM")) %>%
pivot_longer(cols=-c(GeneID, annotation_type, rna_type), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM")


expr_intergenes <- exprdata %>%
        dplyr::filter(annotation_type %in% c('intergenic')) %>%
        inner_join(allgenes %>%
          dplyr::select(GeneID, annotation_type, top_repeat.name) %>%
            dplyr::mutate(top_repeat.name = case_when(stringr::str_detect(top_repeat.name, "\\|") ~ "ambiguous",top_repeat.name =="-1" ~ "ambiguous",top_repeat.name =="Unknown" ~ "ambiguous", T ~ top_repeat.name)) , by=c('GeneID','annotation_type')) %>%
dplyr::select(GeneID, rna_type, top_repeat.name, starts_with("FPM")) %>%
  #mutate(rna_type = factor(rna_type, levels = c("dELS","pELS","PLS","H3K4me3","CTCF"))) %>%
pivot_longer(cols=-c(GeneID, rna_type, top_repeat.name), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM")
```

## Percent of RNAs at threshold

``` r
g0 <- allgenes %>% dplyr::select(GeneID, annotation_type, rna_type, rna_subtype, rna_type_gencode) %>%
          dplyr::filter(annotation_type=='exons', rna_type %in% c('ncRNA','mRNA','lncRNA')) %>%
  dplyr::select(GeneID, rna_type) %>%
  dplyr::count(rna_type) %>%
  dplyr::rename(Ntotal=n) %>%
  mutate(rna_type = factor(rna_type, levels = c("mRNA","lncRNA","ncRNA")))

g0introns <- allgenes %>% 
  dplyr::filter(annotation_type=='introns', ilen>0) %>%
  dplyr::select(GeneID, rna_type, rna_subtype) %>%
          dplyr::filter(rna_type %in% c('ncRNA','mRNA','lncRNA')) %>%
  dplyr::select(GeneID, rna_type) %>%
  dplyr::count(rna_type) %>%
  dplyr::rename(Ntotal=n) %>%
  mutate(rna_type = factor(rna_type, levels = c("mRNA","lncRNA","ncRNA")))

g0intergenes <- allgenes %>%  
   dplyr::filter(annotation_type=='intergenic') %>%
  dplyr::select(GeneID, rna_type) %>%
  dplyr::count(rna_type) %>%
  dplyr::rename(Ntotal=n) %>%
   mutate(rna_type = factor(rna_type, levels = c( "repeat","tRNAderived", "snRNAderived", "cre","readthrough","antisense","intergenic")))


plot_census <- function(exprdata, g0, cmap, fname=NULL){

exprdata %>%
  group_by(GeneID, sequencing) %>%
  summarise(rna_type = rna_type[1], FPM=max(FPM)) %>%
  ungroup() %>%
  dplyr::filter(FPM>0.1)%>%
  mutate(FPM_above=case_when(FPM>10 ~ "10", FPM>1 ~ "1", T ~ "0.1")) %>%
  dplyr::count(rna_type, sequencing, FPM_above) %>%
  inner_join(g0, by=c('rna_type')) %>%
    mutate(rna_type = factor(rna_type, levels = levels(g0$rna_type))) %>%
  mutate(nper = n/Ntotal*100) %>%
  ggplot(aes(x=rna_type, y=nper, fill=FPM_above))+
  geom_col(stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_brewer(palette=cmap)+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% RNAs", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~sequencing) +
    labs(fill="FPM cutoff")->p
#print(p$data)
  
  fname2=NULL
  if (! is.null(fname)){
    fname2=here('figures/carnas', fname)
  }
  p_fixed<- prettysave(p,fname2 , panel.width= 0.2/6*nrow(p$data), panel.height=1.5)
plot_grid(p_fixed)
}
```

### by rna type

``` r
plot_census(expr_intergenes, g0intergenes, "Reds", fname='rnadiversity.barplot.percent.intergenic.pdf')
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=4.6, fig.height=2.9"

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
plot_census(expr_genes %>% dplyr::filter(annotation_type=='introns'), g0introns, "Purples", fname="rnadiversity.barplot.percent.introns.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=3, fig.height=2.6"

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
plot_census(expr_genes %>% dplyr::filter(annotation_type=='exons'), g0, "Blues", fname="rnadiversity.barplot.percent.exons.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=3, fig.height=2.6"

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->
\### subtype CREs

``` r
expr_cre <-
   exprdata %>%
        dplyr::filter(annotation_type %in% c('intergenic'), rna_type=='cre') %>%
  dplyr::select(-rna_type) %>%
  dplyr::rename(rna_type = rna_subtype) %>%
dplyr::select(GeneID, rna_type, starts_with("FPM")) %>%
  #mutate(rna_type = factor(rna_type, levels = c("dELS","pELS","PLS","H3K4me3","CTCF"))) %>%
pivot_longer(cols=-c(GeneID, rna_type), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM")

g0cre <- allgenes %>%  
  dplyr::filter(rna_type=="cre") %>%
  dplyr::select(GeneID, rna_subtype) %>%
  dplyr::rename(rna_type = rna_subtype) %>%
  dplyr::count(rna_type) %>%
  dplyr::rename(Ntotal=n) %>%
    mutate(rna_type = factor(rna_type, levels = c("dELS","pELS","PLS","DNase-H3K4me3","CTCF-only"))) 


plot_census(expr_cre, g0cre, "Reds", fname="rnadiversity.barplot.percent.CREs.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=3.8, fig.height=3.1"

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->
\### subtype repeats

``` r
expr_rep <-
   exprdata %>%
        dplyr::filter(annotation_type %in% c('intergenic')) %>%
        inner_join(allgenes %>%
                     dplyr::filter(annotation_type=='intergenic') %>%
          dplyr::select(GeneID, top_repeat.name)) %>%
            dplyr::filter(rna_type=="repeat") %>%
            dplyr::mutate(rna_type = case_when(stringr::str_detect(top_repeat.name, "\\|") ~ "ambiguous",top_repeat.name =="-1" ~ "ambiguous",top_repeat.name =="Unknown" ~ "ambiguous", T ~ top_repeat.name))  %>%
dplyr::select(GeneID, rna_type, starts_with("FPM")) %>%
  #mutate(rna_type = factor(rna_type, levels = c("dELS","pELS","PLS","H3K4me3","CTCF"))) %>%
pivot_longer(cols=-c(GeneID, rna_type), names_to = c("cell", "sequencing"),
names_pattern = "FPM\\.(.*)\\.(.*)",
values_to = "FPM")
```

    ## Joining, by = "GeneID"

``` r
g0rep <- allgenes %>%  
  dplyr::filter(rna_type=="repeat") %>%
  dplyr::select(GeneID, top_repeat.name) %>%
   dplyr::mutate(rna_type = case_when(stringr::str_detect(top_repeat.name, "\\|") ~ "ambiguous",top_repeat.name =="-1" ~ "ambiguous",top_repeat.name =="Unknown" ~ "ambiguous", T ~ top_repeat.name)) %>%
  dplyr::count(rna_type) %>%
  dplyr::rename(Ntotal=n) %>%
    mutate(rna_type = factor(rna_type, levels = c("LTR","LINE","SINE","DNA","Satellite","Simple_repeat","rRNA","Low_complexity","Retroposon","ambiguous"))) 


plot_census(expr_rep, g0rep, "Reds", fname="rnadiversity.barplot.percent.repeats.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=4.4, fig.height=2.8"

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->
\## Absolute number of RNAs at threshold

``` r
plot_N <- function(exprdata, g0, cmap, fname=NULL){

s <- exprdata %>%
  group_by(GeneID, sequencing) %>%
  summarise(rna_type = rna_type[1], FPM=max(FPM)) %>%
  ungroup() %>%
  dplyr::filter(FPM>0.1)%>%
  dplyr::count(rna_type, sequencing) %>%
    dplyr::rename(Nabove = n) %>%
  inner_join(g0, by=c('rna_type')) %>%
    mutate(rna_type = factor(rna_type, levels = levels(g0$rna_type)))


   

s %>%
   dplyr::mutate(Ntotal = Ntotal-Nabove) %>%
  pivot_longer(starts_with("N"), names_prefix = "N", values_to = "N", names_to="measured") %>%
    dplyr::mutate(measured = factor(case_when(measured=="above"~"Expressed",T~"In transcriptome"), levels = c("In transcriptome", "Expressed"))) %>%
  ggplot(aes(x=rna_type, y=N, fill=measured))+
  geom_col(color = 'black', stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_manual(values = c(Expressed = RColorBrewer::brewer.pal(3,cmap)[[1]], `In transcriptome` = "white"))+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="# RNAs", fill=NULL)+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~sequencing) ->p
#print(p$data)

 fname2=NULL
  if (! is.null(fname)){
    fname2=here('figures/carnas', fname)
  }
  p_fixed<- prettysave(p,fname2 , panel.width= 0.2/4*nrow(p$data), panel.height=1.5)
#plot_grid(p_fixed)

return(list(data=s, plt=p_fixed))
}


plot_N2 <- function(exprdata, g0, cmap){

exprdata %>%
  group_by(GeneID, sequencing) %>%
  summarise(rna_type = rna_type[1], FPM=max(FPM)) %>%
  ungroup() %>%
  dplyr::filter(FPM>0.1)%>%
  mutate(FPM_above=case_when(FPM>10 ~ "10", FPM>1 ~ "1", T ~ "0.1")) %>%
  dplyr::count(rna_type, sequencing, FPM_above) %>%
  inner_join(g0, by=c('rna_type')) %>%
    mutate(rna_type = factor(rna_type, levels = levels(g0$rna_type))) %>%
  mutate(nper = n/Ntotal*100) %>%
  ggplot(aes(x=rna_type, y=nper, fill=FPM_above))+
  geom_col(stats="identity") +
  #ggsci::scale_fill_nejm()+
  #1f78b4
  scale_fill_brewer(palette=cmap)+
scale_y_continuous(expand = c(0,0)) + theme_publish()+
  
  theme(legend.position = "right", legend.direction = "vertical")+
  labs(x=NULL, y="% RNAs", fill=NULL)+
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~sequencing) +
    labs(fill="FPM cutoff")->p
#print(p$data)
p_fixed <- egg::set_panel_size(p=p, margin = unit(0, "in"), width = unit(0.2/6*nrow(p$data), "in"), height = unit(1.5, "in"))
plot_grid(p_fixed)
}
```

``` r
expr_genes %>% dplyr::filter(annotation_type=='exons') %>%
  group_by(GeneID, sequencing) %>%
  summarise(rna_type = rna_type[1], FPM=max(FPM)) %>%
  ungroup() %>%
  dplyr::filter(FPM>0.1)%>%
  mutate(FPM_above=case_when(FPM>10 ~ "10", FPM>1 ~ "1", T ~ "0.1")) %>%
  dplyr::count(rna_type, sequencing, FPM_above) %>%
  left_join(g0, by=c('rna_type')) %>%
    mutate(rna_type = factor(rna_type, levels = levels(g0$rna_type))) %>%
  mutate(nper = n/Ntotal*100) 
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 24 × 6
    ##    rna_type sequencing FPM_above     n Ntotal  nper
    ##    <fct>    <chr>      <chr>     <int>  <int> <dbl>
    ##  1 lncRNA   char       0.1        3474  17142 20.3 
    ##  2 lncRNA   char       1          1143  17142  6.67
    ##  3 lncRNA   char       10          173  17142  1.01
    ##  4 lncRNA   rna        0.1        2580  17142 15.1 
    ##  5 lncRNA   rna        1           893  17142  5.21
    ##  6 lncRNA   rna        10          173  17142  1.01
    ##  7 mRNA     char       0.1        3292  19864 16.6 
    ##  8 mRNA     char       1          7579  19864 38.2 
    ##  9 mRNA     char       10         4143  19864 20.9 
    ## 10 mRNA     rna        0.1        2292  19864 11.5 
    ## # … with 14 more rows

### by rna type

``` r
out = plot_N(expr_genes %>% dplyr::filter(annotation_type=='exons'), g0, "Blues","rnadiversity.barplot.N.exons.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=3.9, fig.height=2.6"

``` r
print(out$data)
```

    ## # A tibble: 6 × 4
    ##   rna_type sequencing Nabove Ntotal
    ##   <fct>    <chr>       <int>  <int>
    ## 1 lncRNA   char         4790  17142
    ## 2 lncRNA   rna          3646  17142
    ## 3 mRNA     char        15014  19864
    ## 4 mRNA     rna         15199  19864
    ## 5 ncRNA    char         2667  20469
    ## 6 ncRNA    rna          1999  20469

``` r
plot_grid(out$plt)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
out = plot_N(expr_genes %>% dplyr::filter(annotation_type=='introns'), g0introns, "Purples", fname="rnadiversity.barplot.N.introns.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=3.9, fig.height=2.6"

``` r
print(out$data)
```

    ## # A tibble: 6 × 4
    ##   rna_type sequencing Nabove Ntotal
    ##   <fct>    <chr>       <int>  <int>
    ## 1 lncRNA   char         6931  12493
    ## 2 lncRNA   rna          4834  12493
    ## 3 mRNA     char        15896  18772
    ## 4 mRNA     rna         14625  18772
    ## 5 ncRNA    char          993   3774
    ## 6 ncRNA    rna          1365   3774

``` r
plot_grid(out$plt)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
out = plot_N(expr_intergenes, g0intergenes, "Reds", fname="rnadiversity.barplot.N.intergenic.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=5.5, fig.height=2.9"

``` r
print(out$data)
```

    ## # A tibble: 14 × 4
    ##    rna_type     sequencing Nabove Ntotal
    ##    <fct>        <chr>       <int>  <int>
    ##  1 antisense    char         3444  16571
    ##  2 antisense    rna          1398  16571
    ##  3 cre          char         5114   6974
    ##  4 cre          rna          4204   6974
    ##  5 intergenic   char         5454  24487
    ##  6 intergenic   rna          2338  24487
    ##  7 readthrough  char         7671   8075
    ##  8 readthrough  rna          7140   8075
    ##  9 repeat       char         8360  21752
    ## 10 repeat       rna          2101  21752
    ## 11 snRNAderived char          134    184
    ## 12 snRNAderived rna            51    184
    ## 13 tRNAderived  char          265    271
    ## 14 tRNAderived  rna           217    271

``` r
plot_grid(out$plt)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### cre subtypes

``` r
out = plot_N(expr_cre, g0cre, "Reds", fname="rnadiversity.barplot.N.cre.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=4.6, fig.height=3.1"

``` r
print(out$data)
```

    ## # A tibble: 10 × 4
    ##    rna_type      sequencing Nabove Ntotal
    ##    <fct>         <chr>       <int>  <int>
    ##  1 CTCF-only     char           87    192
    ##  2 CTCF-only     rna            34    192
    ##  3 dELS          char         1930   2801
    ##  4 dELS          rna          1453   2801
    ##  5 DNase-H3K4me3 char           56    110
    ##  6 DNase-H3K4me3 rna            29    110
    ##  7 pELS          char         1539   2085
    ##  8 pELS          rna          1273   2085
    ##  9 PLS           char         1502   1786
    ## 10 PLS           rna          1415   1786

``` r
plot_grid(out$plt)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->
\### repeats subtypes

``` r
out = plot_N(expr_rep, g0rep, "Reds", fname="rnadiversity.barplot.N.repeats.pdf")
```

    ## `summarise()` has grouped output by 'GeneID'. You can override using the
    ## `.groups` argument.

    ## Warning: Ignoring unknown parameters: stats

    ## [1] "fig.width=5.4, fig.height=2.8"

``` r
print(out$data)
```

    ## # A tibble: 14 × 4
    ##    rna_type   sequencing Nabove Ntotal
    ##    <fct>      <chr>       <int>  <int>
    ##  1 ambiguous  char         1960   5125
    ##  2 ambiguous  rna           322   5125
    ##  3 DNA        char          186    483
    ##  4 DNA        rna            33    483
    ##  5 LINE       char         1059   3571
    ##  6 LINE       rna           423   3571
    ##  7 LTR        char          830   1968
    ##  8 LTR        rna           331   1968
    ##  9 Retroposon char            7     97
    ## 10 Retroposon rna             4     97
    ## 11 Satellite  char           53    330
    ## 12 Satellite  rna            50    330
    ## 13 SINE       char         4174   9807
    ## 14 SINE       rna           892   9807

``` r
plot_grid(out$plt)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

## Expression distribution

### exons, introns, UNAs

``` r
expr_genes %>%
dplyr::select(GeneID, cell, sequencing, annotation_type, FPM) %>%
  bind_rows(expr_intergenes %>% dplyr::select(GeneID, cell, sequencing, FPM) %>% mutate(annotation_type="UNA")) %>%
  ungroup() %>%
  mutate(annotation_type=factor(annotation_type, levels = c("exons","introns","UNA")), sequencing= factor(sequencing, levels = c("char","rna")), cell = factor(cell, levels = c("ES","DE"))) %>%
 dplyr::filter(FPM>0.1) %>%
 ggplot(aes(y = annotation_type , x = log10(FPM), color=sequencing)) +
  stat_slab(aes(slab_color=sequencing), fill=NA, normalize="panels") +
 stat_pointinterval(point_interval = median_qi, .width = c(.75), interval_size=1, position = position_dodge(width = .15, preserve = "single"))+
    facet_wrap(~cell)+
   #scale_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  #scale_color_manual(aesthetics = "slab_color", values=c(ES="cornflowerblue", DE="gold3"))+
 # scale_fill_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  #scale_slab_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  theme_publish() + 
  scale_x_continuous(breaks = seq(-2,3,1),expand=c(0,0), labels = scales::math_format(10^.x))+ #scales::math_format
  annotation_logticks(sides = "b")+
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = "dotted"))+
  labs(x="FPM", y=NULL)->p

p_fixed<- prettysave(p,here("figures/carnas/annotation_types.exprdist.eix_compareSequencing.pdf") , panel.width= 2.1, panel.height=0.5*3)
```

    ## [1] "fig.width=5, fig.height=2.8"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
expr_genes %>%
dplyr::select(GeneID, cell, sequencing, annotation_type, FPM) %>%
  bind_rows(expr_intergenes %>% dplyr::select(GeneID, cell, sequencing, FPM) %>% mutate(annotation_type="UNA")) %>%
  ungroup() %>%
  mutate(annotation_type=factor(annotation_type, levels = c("exons","introns","UNA")), sequencing= factor(sequencing, levels = c("char","rna")), cell = factor(cell, levels = c("ES","DE"))) %>%
dplyr::filter(FPM>0.1) %>%
 ggplot(aes(y = annotation_type , x = log10(FPM), color=cell)) +
  stat_slab(aes(slab_color=cell), fill=NA, normalize="panels") +
 stat_pointinterval(point_interval = median_qi, .width = c(.75), interval_size=1, position = position_dodge(width = .15, preserve = "single"))+
    facet_wrap(~sequencing)+
  scale_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  scale_color_manual(aesthetics = "slab_color", values=c(ES="cornflowerblue", DE="gold3"))+
  theme_publish() + 
  scale_x_continuous(breaks = seq(-2,3,1),expand=c(0,0), labels = scales::math_format(10^.x))+ #scales::math_format
  annotation_logticks(sides = "b")+
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = "dotted"))+
  labs(x="FPM", y=NULL)->p

p_fixed<- prettysave(p,here("figures/carnas/annotation_types.exprdist.eix_compareCcell.pdf") , panel.width= 2.1, panel.height=0.5*3)
```

    ## [1] "fig.width=5, fig.height=2.8"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

### by RNA type (all UNAs combined)

``` r
expr_genes %>%
  dplyr::filter(annotation_type=="exons") %>%
  #mutate(rna_type = paste0(rna_type,"_",annotation_type)) %>%
  dplyr::select(-annotation_type) %>%
  bind_rows(expr_intergenes %>% dplyr::select(GeneID, cell, sequencing, FPM) %>% mutate(rna_type="UNA")) %>%
  ungroup() %>%
  mutate(rna_type = factor(rna_type, levels = c("mRNA", "lncRNA","ncRNA","UNA")), sequencing= factor(sequencing, levels = c("char","rna")), cell = factor(cell, levels = c("ES","DE"))) %>%
dplyr::filter(FPM>0.1, !is.na(rna_type)) %>%
 ggplot(aes(y = rna_type , x = log10(FPM), color=cell)) +
  stat_slab(aes(slab_color=cell), fill=NA, normalize="panels") +
  stat_pointinterval(point_interval = median_qi, .width = c(.75), interval_size=1, position = position_dodge(width = .15, preserve = "single"))+
    facet_wrap(~sequencing)+
   scale_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  scale_color_manual(aesthetics = "slab_color", values=c(ES="cornflowerblue", DE="gold3"))+
 # scale_fill_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  #scale_slab_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  theme_publish() + 
  scale_x_continuous(breaks = seq(-2,3,1),expand=c(0,0), labels = scales::math_format(10^.x))+ #scales::math_format
  annotation_logticks(sides = "b")+
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = "dotted"))+
  labs(x="FPM", y=NULL)->p

p_fixed<- prettysave(p,here("figures/carnas/annotation_types.exprdist.exons_compareCells.pdf") , panel.width= 1.8, panel.height=0.5*4)
```

    ## [1] "fig.width=4.5, fig.height=3.3"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
expr_genes %>%
  dplyr::filter(annotation_type=="exons") %>%
  #mutate(rna_type = paste0(rna_type,"_",annotation_type)) %>%
  dplyr::select(-annotation_type) %>%
  bind_rows(expr_intergenes %>% dplyr::select(GeneID, cell, sequencing, FPM) %>% mutate(rna_type="UNA")) %>%
  ungroup() %>%
  mutate(rna_type = factor(rna_type, levels = c("mRNA", "lncRNA","ncRNA","UNA")), sequencing= factor(sequencing, levels = c("char","rna")), cell = factor(cell, levels = c("ES","DE"))) %>%
dplyr::filter(FPM>0.1, !is.na(rna_type)) %>%
 ggplot(aes(y = rna_type , x = log10(FPM), color=sequencing)) +
  stat_slab(aes(slab_color=sequencing), fill=NA, normalize="panels") +
  stat_pointinterval(point_interval = median_qi, .width = c(.75), interval_size=1, position = position_dodge(width = .15, preserve = "single"))+
    facet_wrap(~cell)+
   #scale_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  #scale_color_manual(aesthetics = "slab_color", values=c(ES="cornflowerblue", DE="gold3"))+
 # scale_fill_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  #scale_slab_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  theme_publish()+ 
  scale_x_continuous(breaks = seq(-2,3,1),expand=c(0,0), labels = scales::math_format(10^.x))+ #scales::math_format
  annotation_logticks(sides = "b")+
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = "dotted"))+
  labs(x="FPM", y=NULL)->p

p_fixed<- prettysave(p,here("figures/carnas/annotation_types.exprdist.exons_compareSequencing.pdf") , panel.width= 1.8, panel.height=0.5*4)
```

    ## [1] "fig.width=4.5, fig.height=3.3"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

### UNAs details

``` r
expr_intergenes  %>%
  mutate(sequencing= factor(sequencing, levels = c("char","rna")), cell = factor(cell, levels = c("ES","DE"))) %>%
  
dplyr::mutate(rna_type = factor(rna_type, levels = c("repeat","tRNAderived", "snRNAderived", "cre","readthrough","antisense","intergenic")))  %>%
  dplyr::filter(FPM>0.1, !is.na(rna_type)) %>%
 ggplot(aes(y = rna_type , x = log10(FPM), color=cell)) +
  stat_slab(aes(slab_color=cell), fill=NA, normalize="panels") +
  stat_pointinterval(point_interval = median_qi, .width = c(.75), interval_size=1, position = position_dodge(width = .15, preserve = "single"))+
    facet_wrap(~sequencing)+
   scale_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  scale_color_manual(aesthetics = "slab_color", values=c(ES="cornflowerblue", DE="gold3"))+
 # scale_fill_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  #scale_slab_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  theme_publish()+ 
  scale_x_continuous(breaks = seq(-2,3,1),expand=c(0,0), labels = scales::math_format(10^.x))+ #scales::math_format
  annotation_logticks(sides = "b")+
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = "dotted"))+
  labs(x="FPM", y=NULL)->p

p_fixed<- prettysave(p,here("figures/carnas/annotation_types.exprdist.intergenic.pdf") , panel.width= 1.8, panel.height=0.5*7)
```

    ## [1] "fig.width=5, fig.height=4.8"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

### repeat details

``` r
expr_intergenes  %>%
  mutate(sequencing= factor(sequencing, levels = c("char","rna")), cell = factor(cell, levels = c("ES","DE"))) %>%
  dplyr::filter(rna_type=="repeat") %>%
   mutate(top_repeat.name = factor(top_repeat.name, levels = c("LTR","LINE","SINE","DNA","Satellite","Simple_repeat","rRNA","Low_complexity","Retroposon","ambiguous"))) %>%
  dplyr::filter(FPM>0.1, !is.na(top_repeat.name)) %>%
#dplyr::mutate(rna_type = factor(rna_type, levels = c("repeat","tRNAderived", "sRNAderived", "cre","geneproximal","antisense","intergenic")))  %>%
 ggplot(aes(y = top_repeat.name , x = log10(FPM), color=cell)) +
  stat_slab(aes(slab_color=cell), fill=NA, normalize="panels") +
  stat_pointinterval(point_interval = median_qi, .width = c(.75), interval_size=1, position = position_dodge(width = .15, preserve = "single"))+
    facet_wrap(~sequencing)+
   scale_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  scale_color_manual(aesthetics = "slab_color", values=c(ES="cornflowerblue", DE="gold3"))+
 # scale_fill_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  #scale_slab_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  theme_publish()+ 
  scale_x_continuous(breaks = seq(-2,3,1),expand=c(0,0), labels = scales::math_format(10^.x))+ #scales::math_format
  annotation_logticks(sides = "b")+
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = "dotted"))+
  labs(x="FPM", y=NULL)->p

p_fixed<- prettysave(p,here("figures/carnas/annotation_types.exprdist.intergenicREPEATS.pdf") , panel.width= 1.8, panel.height=0.5*7)
```

    ## [1] "fig.width=4.8, fig.height=4.8"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

### CRE details

``` r
expr_intergenes  %>%
  inner_join(allgenes %>% dplyr::select(GeneID, rna_subtype)) %>%
  mutate(sequencing= factor(sequencing, levels = c("char","rna")), cell = factor(cell, levels = c("ES","DE"))) %>%
  dplyr::filter(rna_type=="cre") %>%
  dplyr::filter(FPM>0.1) %>%
#dplyr::mutate(rna_type = factor(rna_type, levels = c("repeat","tRNAderived", "sRNAderived", "cre","geneproximal","antisense","intergenic")))  %>%
 ggplot(aes(y = rna_subtype , x = log10(FPM), color=cell)) +
  stat_slab(aes(slab_color=cell), fill=NA, normalize="panels") +
  stat_pointinterval(point_interval = median_qi, .width = c(.75), interval_size=1, position = position_dodge(width = .15, preserve = "single"))+
    facet_wrap(~sequencing)+
   scale_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  scale_color_manual(aesthetics = "slab_color", values=c(ES="cornflowerblue", DE="gold3"))+
 # scale_fill_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  #scale_slab_color_manual(values=c(ES="cornflowerblue", DE="gold3"))+
  theme_publish()+ 
  scale_x_continuous(breaks = seq(-2,3,1),expand=c(0,0), labels = scales::math_format(10^.x))+ #scales::math_format
  annotation_logticks(sides = "b")+
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = "dotted"))+
  labs(x="FPM", y=NULL)->p
```

    ## Joining, by = "GeneID"

``` r
p_fixed<- prettysave(p,here("figures/carnas/annotation_types.exprdist.intergenicCRE.pdf") , panel.width= 1.8, panel.height=0.5*5)
```

    ## [1] "fig.width=5.2, fig.height=3.8"

``` r
plot_grid(p_fixed)
```

![](visu_carnas_basicstats_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->
