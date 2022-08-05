import os
import pickle
from tqdm import tqdm


def load_RD_data(libs, genes,  smpls=[], loadrna=False, loaddna=True, loadbkg=True, loadmodels=False, sfx = "Q255ok", sfx_rna="Q255"):
    cur_msg='GENES'
    nstep=1+int(loadmodels)
    

    prfx="compress"
    
    paths_dna={'exons':'pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/dna.%s.%s.chartable.%s.%g.pickle'%(sfx, 'exons',prfx,10), 
          'introns':'pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/dna.%s.%s.chartable.%s.%g.pickle'%(sfx, 'introns',prfx,10),
              'all':'pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/dna.%s.%s.chartable.%s.%g.pickle'%(sfx, 'all',prfx,10),
              'intergenic': 'pairs/strngtieCharOnly_hg38/exons/filtered/DNAq15-blacklisted_RNAunq/dna.%s.%s.chartable.%s.%g.pickle'%(sfx, 'all',prfx,10)}
    paths_rna={'exons':'pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/rna.%s.cis.%s.chartable.%g.pickle'%(sfx_rna,'exons',1), 
                   'introns':'pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/rna.%s.cis.%s.chartable.%g.pickle'%(sfx_rna,'introns',1),
                       'all':'pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/rna.%s.cis.%s.chartable.%g.pickle'%(sfx_rna,'all',1),
                       'intergenic':'pairs/strngtieCharOnly_hg38/exons/filtered/DNAq15-blacklisted_RNAunq/rna.%s.cis.%s.chartable.%g.pickle'%(sfx_rna,'all',1)}
    path_bk='pairs/gencondeV29_hg38/all/filtered/DNAq15-blacklisted_RNAunq/bkg.%s.chartable.%s.%g.pickle'%(sfx, prfx,10)

    if (len(smpls)==0):
        smpls = list(libs.keys())

    dna_to_load = {}
    if ((type(loaddna)==bool) and loaddna):
        dna_to_load=paths_dna
    else:
        if type(loaddna)!=bool:
            if type(loaddna)==str:
                loaddna = [loaddna]
            dna_to_load = {k: paths_dna[k] for k in loaddna}

    rna_to_load = {}
    if ((type(loadrna)==bool) and loadrna):
        rna_to_load=paths_rna
    else:
        if type(loadrna)!=bool:
            if type(loadrna)==str:
                loadrna = [loadrna]
            rna_to_load = {k: paths_rna[k] for k in loadrna}
    
    for s in smpls:
        nstep+=int(loadbkg)
        nstep+=len(rna_to_load)
        nstep+=len(dna_to_load)
    
    
    with tqdm(total=nstep) as pbar:
        pbar.set_postfix(step=cur_msg, refresh=True)
        
        
    # Load the RNAcis and DNA chartables

        dna={}
        if len(dna_to_load)>0:
            dna={k:{} for k in smpls}
            pbar.update(1)
            for sname in smpls:
                for k, v in dna_to_load.items():
                    cur_msg='DNA:'+sname+":"+k
                    pbar.set_postfix(step=cur_msg, refresh=True)
                    dna[sname][k]=pickle.load(open(os.path.join(libs[sname], v), "rb")) #libs[s]
                    pbar.update(1)

        rna={}
        if len(rna_to_load)>0:
            rna={k:{} for k in smpls}
            for sname in smpls:
                for k, v in rna_to_load.items():
                    cur_msg='RNA:'+sname+":"+k
                    pbar.set_postfix(step=cur_msg, refresh=True)
                    rna[sname][k]=pickle.load(open(os.path.join(libs[sname], v), "rb"))
                    pbar.update(1)

    #Load bkg
        bkg={}
        if loadbkg:
            bkg={}
            B={}


            for sname in smpls:
                cur_msg='bkg:'+sname
                pbar.set_postfix(step=cur_msg, refresh=True)
                bkg[sname]=pickle.load(open(os.path.join(libs[sname], path_bk), "rb"))
                B[sname]=bkg[sname]['all'].normalize_row_bytrans().counts.tocsr()
                pbar.update(1)

    #Load fly data
        cismodels={}
        if loadmodels:
            cur_msg='Models'
            pbar.set_postfix(step=cur_msg, refresh=True)
            cismodels={k:{} for k in smpls}

            for sname in smpls:
                for t in ['exons','introns']:
                    cismodels[sname][t]=pickle.load(open(os.path.join(libs[sname],'pairs/gencondeV29_hg38/%s/analysis/01_flight/cismodel.pickle'%t),'rb'))

        genetypes = set(genes.keys())
        gene_strand_dict = {gt : {k: (True if v=="-" else False) for k, v in genes[gt].strand.to_dict().items()} for gt in genetypes} 
        rd_data={'genes': genes, 'gene_strand_dict':gene_strand_dict, 'dna':dna,'rna':rna,'bkg':bkg, 'cismodel':cismodels, 'libs':libs}
        pbar.update(1)
        pbar.set_postfix(refresh=True)
    
    return(rd_data)