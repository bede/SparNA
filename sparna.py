#!/usr/bin/env python3

# Author: Bede Constantinides - b|at|bede|dot|im

# TODO
# | Decide on BWA vs Bowtie2
# | GZIP support
# | Report % read alignment in mapping to ref and contig
# | Interleaved reads (ONE TRUE FORMAT)
# | add minimum similarity threshold for reference selection
# | report on trimming, %remapped
# | Bootstrap/dogfoood assemblies with --trusted-contigs etc ?

# DEPENDENCIES
# | python packages:
# |    argh, biopython, khmer, plotly
# | others, expected inside $PATH:
# |    bwa, bowtie2, samtools, vcftools, bcftools, bedtools, seqtk, spades, quast
# | others, bundled inside res/ directory:
# |    trimmomatic

import os
import io
import sys
import argh
import time
import json
import pandas
import logging
import requests
import subprocess
import multiprocessing
import concurrent.futures

import pprint

from collections import OrderedDict

from Bio import SeqIO
from Bio import SeqUtils

import plotly.offline as py
import plotly.graph_objs as go

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run(cmd):
    return subprocess.run(cmd,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)


def name_sample(fwd_fq):
    fwd_fq_prefix = os.path.splitext(os.path.split(fwd_fq)[1])[0]
    return fwd_fq_prefix


def import_reads(fwd_fq, rev_fq, params):
    print('Importing reads...')
    cmd = (
    'cat {fwd_fq} > {out}/raw/{name}.f.fastq '
    '&& cat {rev_fq} > {out}/raw/{name}.r.fastq '
    '&& interleave-reads.py {out}/raw/{name}.f.fastq {out}/raw/{name}.r.fastq '
    '> {out}/raw/{name}.fr.fastq'
    .format(**params, fwd_fq=fwd_fq, rev_fq=rev_fq))
    logger.info(cmd)
    cmd_run = run(cmd)
    logger.info(cmd_run.stdout, cmd_run.stderr)
    print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_IMPORT')


def trim(norm_k_list, params):
    print('Trimming...')
    # Fetch smallest norm_k for trimming with Trimmomatic min_len - Screed bug workaround
    params['min_len'] = max(map(int, norm_k_list.split(',')))
    cmd = (
    'java -jar {pipe}/res/trimmomatic-0.33.jar PE'
    ' {out}/raw/{name}.f.fastq'
    ' {out}/raw/{name}.r.fastq'
    ' {out}/trim/{name}.f_pe.fastq'
    ' {out}/trim/{name}.f_se.fastq'
    ' {out}/trim/{name}.r_pe.fastq'
    ' {out}/trim/{name}.r_se.fastq'
    ' ILLUMINACLIP:{pipe}/res/illumina_adapters.fa:2:30:10'.format(**params))
    if params['qual_trim']:
        print('\tQuality trimming...')
        logger.info('Quality trimming...')
        cmd += ' SLIDINGWINDOW:4:20'.format(**params)
    cmd += (
    ' MINLEN:{min_len}'
    ' && cat {out}/trim/{name}.f_se.fastq {out}/trim/{name}.r_se.fastq'
    ' > {out}/trim/{name}.se.fastq'
    ' && interleave-reads.py {out}/trim/{name}.f_pe.fastq {out}/trim/{name}.r_pe.fastq'
    ' > {out}/trim/{name}.fr.fastq'.format(**params))
    logger.info(cmd)
    cmd_run = run(cmd)
    logger.info(cmd_run.stderr)
    print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_TRIM')


def normalise(norm_perms, params):
    print('Normalising...')
    cmds = []
    for norm_perm in norm_perms:
        cmd_vars = dict(**params,
                        k=str(norm_perm['k']),
                        c=str(norm_perm['c']))
        cmd = (
        'normalize-by-median.py -C {c} -k {k} -N 4 -x 1e8 -p'
        ' {out}/trim/{name}.fr.fastq'
        ' -o {out}/norm/{name}.norm_k{k}c{c}.fr.fastq'
        ' && normalize-by-median.py -C {c} -k {k} -N 4 -x 1e9'
        ' {out}/trim/{name}.se.fastq'
        ' -o {out}/norm/{name}.norm_k{k}c{c}.se.fastq'
        ' && split-paired-reads.py'
        ' -1 {out}/norm/{name}.norm_k{k}c{c}.f_pe.fastq'
        ' -2 {out}/norm/{name}.norm_k{k}c{c}.r_pe.fastq'
        ' {out}/norm/{name}.norm_k{k}c{c}.fr.fastq'
        ' && cat {out}/norm/{name}.norm_k{k}c{c}.fr.fastq'
        ' {out}/norm/{name}.norm_k{k}c{c}.se.fastq >'
        ' {out}/norm/{name}.norm_k{k}c{c}.pe_and_se.fastq'.format(**cmd_vars))
        cmds.append(cmd)
        print('\tNormalising norm_c={c}, norm_k={k}'.format(**cmd_vars))
        logger.info('Normalising norm_c={c}, norm_k={k}'.format(**cmd_vars))
    with multiprocessing.Pool(params['threads']) as pool:
        results = pool.map(run, cmds)
    logger.info([result.stdout + result.stdout for result in results])
    print('\tAll done') if not max([r.returncode for r in results]) else sys.exit('ERR_NORM')
    return norm_perms


def assemble(asm_perms, params):
    '''
    Performs multiple assemblies and returns and OrderedDict of assembly names and paths
    '''
    print('Assembling...')

    if params['asm_k']:
        asm_k_fmt = 'k' + 'k'.join(params['asm_k'])
    else:
        asm_k_fmt = 'k'

    cmds_asm = []
    cmd_vars = dict(**params,
                    asm_k_fmt=asm_k_fmt)
    
    if params['no_norm']:
        cmd_asm = (
        'spades.py -m 8 -t 12'
        ' --12 {out}/trim/{name}.fr.fastq'
        ' -s {out}/trim/{name}.se.fastq'
        ' -o {out}/asm/{name}.no_norm.asm_{asm_k_fmt} --careful'.format(**cmd_vars))
        print('\tAssembling without prior normalisation'.format(**cmd_vars))
        cmds_asm.append(cmd_asm)

    for asm_perm in asm_perms:
        cmd_vars['k'] = str(asm_perm['k'])
        cmd_vars['c'] = str(asm_perm['c'])
        cmd_asm = (
        'spades.py -m 8 -t {threads}'
        ' --pe1-1 {out}/norm/{name}.norm_k{k}c{c}.f_pe.fastq'
        ' --pe1-2 {out}/norm/{name}.norm_k{k}c{c}.r_pe.fastq'.format(**cmd_vars))
        if params['asm_k']:
            cmd_asm += ' -k {asm_k}'.format(**cmd_vars)
        cmd_asm += (
        ' --s1 {out}/norm/{name}.norm_k{k}c{c}.se.fastq'
        ' -o {out}/asm/{name}.norm_k{k}c{c}.asm_{asm_k_fmt} --careful'.format(**cmd_vars))
        cmds_asm.append(cmd_asm)
        print('\tAssembling norm_c={c}, norm_k={k}, asm_k={asm_k}'.format(**cmd_vars))

    with multiprocessing.Pool(params['threads']) as pool:
        results = pool.map(run, cmds_asm)
    logger.info([result.stdout for result in results])
    print('\tAll done') if not max([r.returncode for r in results]) else sys.exit('ERR_ASM')
    asms = os.listdir(params['out'] + '/asm')
    asm_paths = [params['out'] + '/asm/' + asm + '/contigs.fasta' for asm in asms]
    return OrderedDict(zip(asms, asm_paths))


def prune_assemblies(asms_paths, min_len, params):
    asms_paths_pruned = OrderedDict()
    for asm, path in asms_paths.items():
        asms_paths_pruned[asm] = path.replace('/asm/', '/asm_prune/')
        records = (r for r in SeqIO.parse(path, 'fasta') if len(r.seq) >= min_len)
        os.makedirs(asms_paths_pruned[asm].replace('/contigs.fasta', ''))
        SeqIO.write(records, asms_paths_pruned[asm], 'fasta')
    return asms_paths_pruned 


def gc_content(asms_paths):
    asms_gc = {}
    for asm, path in asms_paths.items():
        asm_gc = []
        for record in SeqIO.parse(path, 'fasta'):
            asm_gc.append(SeqUtils.GC(record.seq)/100)
        asms_gc[asm] = asm_gc
    return asms_gc


def map_to_assemblies(asms_paths, params):
    '''
    Map original reads to each assembly with Bowtie2
    Record mapping statistics
    Screen uniquely mapped reads and quantify reads mapped per contig
    Returns dict of tuples containing contig_len and n_reads_mapped for each contig
    '''
    print('Aligning to assemblies... (Bowtie2)')
    asms_coverages = {}
    for asm, asm_path in asms_paths.items():
        i = 0
        print(asm, i)
        cmd_vars = {**params,
                    'asm':asm,
                    'asm_path':asm_path}
        cmds = [
        'bowtie2-build -q {asm_path} {out}/remap/{asm}',
        'bowtie2 -x {out}/remap/{asm} --no-unal --very-sensitive-local --threads {threads}'
        ' -1 {out}/raw/{name}.f.fastq'
        ' -2 {out}/raw/{name}.r.fastq'
        ' -S {out}/remap/{asm}.sam'
        ' 2> {out}/remap/{asm}.bt2.stats',
        'grep -v XS:i: {out}/remap/{asm}.sam > {out}/remap/{asm}.uniq.sam',
        'samtools view -bS {out}/remap/{asm}.uniq.sam'
        ' | samtools sort - -o {out}/remap/{asm}.uniq.bam',
        'samtools index {out}/remap/{asm}.uniq.bam',
        'samtools idxstats {out}/remap/{asm}.uniq.bam'
        ' > {out}/remap/{asm}.uniq.bam.stats',
        'rm {out}/remap/{asm}.sam {out}/remap/{asm}.uniq.sam']
        cmds = [cmd.format(**cmd_vars) for cmd in cmds]
        for cmd in cmds:
            logger.info(cmd)
            cmd_run = run(cmd)
            logger.info(cmd_run.stdout)
            cmd_prefix = cmd.split(' ')[0]
            print('\tDone (' +cmd_prefix+ ')') if cmd_run.returncode == 0 else sys.exit('ERR_REMAP')
        
        with open('{out}/remap/{asm}.bt2.stats'.format(**cmd_vars), 'r') as bt2_stats:
            map_prop = float(bt2_stats.read().partition('% overall')[0].split('\n')[-1].strip())/100
        
        asm_coverages = []
        with open('{out}/remap/{asm}.uniq.bam.stats'.format(**cmd_vars), 'r') as bam_stats:
            for line in bam_stats:
                if not line.startswith('*'):
                    reads_mapped = int(line.strip().split('\t')[2])
                    asm_coverages.append(int(reads_mapped))
        
        asms_coverages[asm] = asm_coverages
    return asms_coverages


def onecodex_lca(seq, onecodex_api_key):
    '''
    Returns dict of OneCodex real-time API k-mer hits for a given sequence
    e.g. {'elapsed_secs':'0.0003','k': 31,'n_hits': 97,'n_lookups': 128,'tax_id': 9606}
    '''
    url = 'https://app.onecodex.com/api/v0/search'
    payload = {'sequence':str(seq)}
    auth = requests.auth.HTTPBasicAuth(onecodex_api_key, '')
    response = requests.post(url, payload, auth=auth, timeout=5)
    result = json.loads(response.text)
    result['prop_hits'] = round(int(result['n_hits'])/int(result['n_lookups']), 3)
    return result

def ebi_taxid_to_lineage(tax_id):
    '''
    Returns scientific name and lineage for a given taxid using EBI's taxonomy API
    e.g.('Retroviridae', ['Viruses', 'Retro-transcribing viruses'])
    '''
    url = 'http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{}'
    if tax_id == 0 or tax_id == 1:
        return None, None
    response = requests.get(url.format(tax_id), timeout=5)
    result = json.loads(response.text)
    sciname = result['scientificName']
    taxonomy = [x for x in result['lineage'].split('; ') if x]
    return sciname, taxonomy

def onecodex_lca_taxa(seqrecord, onecodex_api_key):
    '''
    Returns a scientific name and lineage for a SeqRecord using OneCodex and EBI APIs
    e.g. ('NODE_3_length_4481_cov_46.6129_ID_7947',
     ('Hepatitis C virus genotype 3',
      ['Viruses',
       'ssRNA viruses',
       'ssRNA positive-strand viruses, no DNA stage',
       'Flaviviridae',
       'Hepacivirus']))
    '''
    hits = onecodex_lca(str(seqrecord.seq), onecodex_api_key)
    sciname, taxonomy = ebi_taxid_to_lineage(hits['tax_id'])
    result = (sciname, taxonomy, hits)
    return result

def fasta_onecodex_lca_taxa(fasta_path, onecodex_api_key):
    '''
    Executes onecodex_lca_taxa() in parallel for a multifasta file
    '''
    seqrecords = SeqIO.parse(fasta_path, 'fasta')
    taxa = {}
    with concurrent.futures.ThreadPoolExecutor(50) as executor:
        futures = {executor.submit(onecodex_lca_taxa, seqrecord, onecodex_api_key): seqrecord for seqrecord in seqrecords}
        for future in concurrent.futures.as_completed(futures):
            seqrecord = futures[future]
            try:
                data = future.result()
            except Exception as exception:
                taxa[seqrecord.id] = (None, None, None)
                logger.info('Skipping {}'.format(seqrecord.id))
            else:
                taxa[seqrecord.id] = future.result()
    return taxa

def onecodex_assemblies(asms_paths, onecodex_api_key):
    '''
    Returns OneCodex hits for a dict of assembly names and corresponding paths 
    '''
    print('Fetching LCA taxonomic assignments... (requires network access)')
    sample_results = OrderedDict()
    for asm_name, asm_path in asms_paths.items():
        print('\tAssembly {}'.format(asm_name))
        sample_results[asm_name] = fasta_onecodex_lca_taxa(asm_path, onecodex_api_key)
    return sample_results

def seqrecords(fasta_path):
    '''
    Accepts path to multifasta, returns list of Biopython SeqRecords
    '''
    return SeqIO.parse(fasta_path, 'fasta')

def lengths(asms_paths):
    '''
    Accepts path to multifasta, returns OrderedDict of sequence lengths 
    '''
    lengths = {}
    for asm_name, asm_path in asms_paths.items():
        lengths[asm_name] = []
        records = seqrecords(asm_path)
        for record in records:
            lengths[asm_name].append(len(record.seq))
    return lengths

def gc_contents(asms_paths):
    '''
    Accepts path to multifasta, returns OrderedDict of sequence GC content
    '''
    gc_contents = {}
    for asm_name, asm_path in asms_paths.items():
        records = seqrecords(asm_path)
        gc_contents[asm_name] = []
        for record in seqrecords:
            gc_contents[record.id].append(SeqUtils.GC(record.seq)/100)
    return gc_contents

def marker_metadata(asms_paths, lengths, gc_contents, taxa):
    '''
    Returns summary metadata for each sequence
    '''
    metadata = {}
    for asm_name, asm_path in asms_paths.items():
        metadata[asm_name] = []
        records = seqrecords(asm_path)
        for i, record in enumerate(records):
            lineage = ';'.join(taxa[asm_name][record.id][1]) if taxa[asm_name][record.id][1] else ''
            lineage_fmt = (lineage[:40] + '..') if len(lineage) > 50 else lineage
            text = (
                '{}<br>'
                'lca: {} (taxid: {})<br>'
                'lineage: {}<br>'
                'length: {}<br>'
                'gc_content: {}<br>'
                ''.format(record.id,
                          taxa[asm_name][record.id][0],
                          0,
                          #taxa[asm_name][record.id][2]['tax_id'],
                          lineage_fmt,
                          lengths[asm_name][i],
                          round(float(gc_contents[asm_name][i]), 3)))
            metadata[asm_name].append(text)
    return metadata


def build_ebi_blast_query(title, sequence, database):
    '''
    Returns dict of REST params for the EBI BLAST API
    '''
    logger.info('building query')
    return { 'email': 'bede.constantinides@manchester.ac.uk',
             'program': 'blastn',
             'stype': 'dna',
             'database': database,
             'align': 6,
             'match_scores': '1,-3',
             'gapopen': 5, 
             'gapextend': 2,
             'exp': '1e-10',
             'filter': 'T',
             'dropoff': 0,
             'scores': 5,
             'alignments': 5,
             'title': title,
             'sequence': str(sequence) }

def parse_hits(title, raw_hits):
    '''
    Returns list of tuples of BLAST hits
    [(blast, tab, output, fields), (blast, tab, output, fields)]
    '''
    hits = []
    for line in io.StringIO(raw_hits):
        if ':' in line:
            fields = [field.strip() for field in line.split('\t')]
            hit = (title, ) + tuple(fields[1].split(':')) + tuple(fields[2:])
            hits.append(hit)
    return hits

def fetch_annotation(database, accession):
    '''
    Return SeqRecord of annotation for given EMBL accession number
    '''
    query = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/{}/{}'.format(database, accession)
    request = requests.get(query)
    annotation = SeqIO.read(io.StringIO(request.text), 'embl')
    return annotation

def ebi_blast(query):
    '''
    Returns BLAST hits as a tuple containing a list of tuples for each hit
    ('seq', [(blast, tab, output, fields),
             (blast, tab, output, fields)])
    '''
    run_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/'
    status_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/'
    results_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/'

    call = requests.post(run_url, data=query)
    start_time = time.time()
    logger.info('dispatched blast jobid: ' + call.text)
    while True:
        status = requests.get(status_url + call.text)
        if status.text == 'FINISHED':
            hits_r = requests.get(results_url + call.text + '/out')
            hits = parse_hits(query['title'], hits_r.text)
            logger.info(status.text + ' ' + call.text)
            logger.info('Job completed in ' + str(time.time() - start_time))
            break
        elif time.time() - start_time > 180:
            print('blast timeout')
            logger.error('blast timeout')
            break
        elif status.text == 'RUNNING':
            time.sleep(2)
        else:
            print('status: ' + status.text)
            logger.error('status: ' + status.text)
            break
    return (query['title'], hits)

def ebi_annotated_blast(query):
    '''
    Returns BLAST hits as a tuple containing a list of tuples of hit tuples and SeqRecord annotations
    ('seq', [((blast, tab, output, fields), SeqRecord),
             ((blast, tab, output, fields), SeqRecord)])
    '''
    run_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/'
    status_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/'
    results_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/'

    call = requests.post(run_url, data=query)
    start_time = time.time()
    logger.info('dispatched blast jobid: ' + call.text)
    while True:
        status = requests.get(status_url + call.text)
        if status.text == 'FINISHED':
            hits_r = requests.get(results_url + call.text + '/out')
            hits = parse_hits(query['title'], hits_r.text)
            annotations_items = [hit[1] + hit[2] for hit in hits] # all there
            annotations = list(fetch_annotation(hit[1], hit[2]) for hit in hits)
            hits_annotations = list(zip(hits, annotations))
            # hits_annotations = list(zip(hits)) TESTING WITHOUT SEQRECORD
            print('\t\tQuery ' + query['title'])
            logger.info(status.text + ' ' + call.text)
            # print(time.time() - start_time)
            break
        elif time.time() - start_time > 180:
            logger.error('blast timeout')
            hits_annotations = None
            break
        elif status.text == 'RUNNING':
            time.sleep(5)
        else:
            logger.error('status: ' + status.text)
            hits_annotations = None
            break
    return (query['title'], hits_annotations)

def fasta_blaster(fasta, database, max_seqs, min_len):
    '''
    NEEDS UPDATING FOR NESTED ORDEREDDICTS
    MIN_LEN NEEDS IMPLEMENTING
    Returns BLAST results as an OrderedDict of ebi_blast() or ebi_annotated_blast() output
    ebi_blast():
    OrderedDict([('seq_1', [(blast, tab, output, fields),
                            (blast, tab, output, fields)],
                 ('seq_2', [(blast, tab, output, fields),
                            (blast, tab, output, fields)])])
    ebi_annotated_blast():
    OrderedDict([('seq_1', [((blast, tab, output, fields), SeqRecord),
                            ((blast, tab, output, fields), SeqRecord)],
                 ('seq_2', [((blast, tab, output, fields), SeqRecord),
                            ((blast, tab, output, fields), SeqRecord)])])
    '''
    records = OrderedDict()
    with open(fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if len(record.seq) >= min_len:
                records[record.id] = record.seq

    queries = [build_ebi_blast_query(title, seq, database) for title, seq in records.items()]
    with multiprocessing.Pool(30) as pool:
        results = pool.map(ebi_annotated_blast, queries[0:max_seqs+1])
        if len(queries) > max_seqs:
            results += zip([q['title'] for q in queries[max_seqs+1:]], [None]*len(queries[max_seqs+1:]))
    return OrderedDict(results)

def blast_assemblies(asms_paths, database, max_seqs, min_len):
    '''
    Returns BLAST hit information for a dict of assembly names and corresponding paths 
    '''
    print('BLASTing assemblies...')
    sample_results = OrderedDict()
    for asm_name, asm_path in asms_paths.items():
        print('\tAssembly {}'.format(asm_name))
        sample_results[asm_name] = fasta_blaster(asm_path, database, max_seqs, min_len)
    return sample_results

def blast_superkingdoms(blast_results):
    asms_superkingdoms = {}
    for asm, contigs in blast_results.items():
        asm_superkingdoms = []
        for contig, hits in contigs.items():
            if hits: # Hits found
                top_hit_superkingdom = hits[0][1].annotations['taxonomy'][0]
                asm_superkingdoms.append(top_hit_superkingdom)
            elif type(hits) is list: # Zero hits
                asm_superkingdoms.append(False)
            elif hits is None: # Not searched
                asm_superkingdoms.append(None)
            assert hits or type(hits) is list or hits is None
        asms_superkingdoms[asm] = asm_superkingdoms
    return asms_superkingdoms

def blast_summary(blast_results, asms_covs):
    asms_summaries = {}
    for asm, contigs in blast_results.items():
        asm_summaries = []
        for i, (contig, hits) in enumerate(contigs.items()):
            if hits: # Hits found
                description = hits[0][1].description[:40] + (hits[0][1].description[40:] and '…')
                top_hit_summary = (''
                '{2}<br>Coverage: {0} reads<br>{1}<br>Accession: {3}:{4}'
                '<br>Identity: {5}%<br>Alignment length: {6}<br>Mismatches: {7}<br>'
                'E-value: {13}'.format(asms_covs[asm][i], description, *hits[0][0]))
                asm_summaries.append(top_hit_summary)
            elif type(hits) is list: # Zero hits
                asm_summaries.append(False)
            elif hits is None: # Not searched
                asm_summaries.append(None)
            assert hits or type(hits) is list or hits is None
        asms_summaries[asm] = asm_summaries
    return asms_summaries


def plotly(asms_names, asms_stats, lca, blast, params):
    cov_max = max(sum([i for i in asms_stats['covs'].values()], []))
    cov_scale_factor = round(cov_max/5000, 1) # For bubble scaling

    traces = []
    for asm_name in asms_names:
        if lca:
            traces.append(
                go.Scatter(
                    x=asms_stats['lens'][asm_name],
                    y=asms_stats['gc'][asm_name],
                    mode='markers',
                    name=asm_name,
                    text=asms_stats['legend'][asm_name],
                    line=dict(shape='spline'),
                    marker=dict(
                        opacity=0.5,
                        symbol='circle',
                        sizemode='area',
                        sizeref=cov_scale_factor,
                        size=asms_stats['covs'][asm_name],
                        line=dict(width=1))))
        elif blast:
            traces.append(
                go.Scatter(
                    x=asms_stats['lens'][asm_name],
                    y=asms_stats['gc'][asm_name],
                    mode='markers',
                    name=asm_name,
                    text=asms_stats['blast_summary'][asm_name],
                    line=dict(shape='spline'),
                    marker=dict(
                        opacity=0.5,
                        symbol='circle',
                        sizemode='area',
                        sizeref=cov_scale_factor,
                        size=asms_stats['covs'][asm_name],
                        line=dict(width=1))))
        else:
            traces.append(
                go.Scatter(
                    x=asms_stats['lens'][asm_name],
                    y=asms_stats['gc'][asm_name],
                    mode='markers',
                    name=asm_name,
                    text=asms_stats['names'][asm_name],
                    line=dict(shape='spline'),
                    marker=dict(
                        opacity=0.5,
                        symbol='circle',
                        sizemode='area',
                        sizeref=cov_scale_factor,
                        size=asms_stats['covs'][asm_name],
                        line=dict(width=1))))

    layout = go.Layout(
        title='Contig length vs. GC content vs. coverage',
        xaxis=dict(
            title='Contig length',
            gridcolor='rgb(255, 255, 255)',
            zerolinewidth=1,
            type='log',
            gridwidth=2),
        yaxis=dict(
            title='GC content',
            gridcolor='rgb(255, 255, 255)',
            zerolinewidth=1,
            gridwidth=2),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)')

    fig = go.Figure(data=traces, layout=layout)
    return py.plot(fig, filename=params['out'] + '/plot.html')


def report(chart_url, start_time, end_time, params):
    elapsed_time = end_time - start_time
    report_content = 'wall_time\t{}\nchart_url\t{}'.format(elapsed_time, chart_url)
    with open(params['out'] + '/summary.txt', 'w') as report:
        report.write(report_content)
        print(report_content)


def main(
    debug=False,
    fwd_fq=None, rev_fq=None,
    qual_trim=False,
    lca=False, blast=False,
    norm_c_list=None, norm_k_list=None,
    asm_k_list=None, no_norm=False,
    onecodex_api_key='a1d32ce32583468192101cc1d0cf27ec',
    blast_db='em_rel', blast_max_seqs=5, min_len=500,
    out_prefix='sparna', threads=4):

    start_time = int(time.time())
    
    params = dict(name=name_sample(fwd_fq),
                  out=out_prefix + '_' + name_sample(fwd_fq),
                  pipe=os.path.dirname(os.path.realpath(__file__)),
                  qual_trim=qual_trim,
                  no_norm=no_norm,
                  norm_c=norm_c_list.split(',') if norm_c_list else None,
                  norm_k=norm_k_list.split(',') if norm_k_list else None,
                  asm_k=asm_k_list if asm_k_list else 0,
                  threads=threads)

    if norm_k_list and norm_c_list:
        norm_perms = [{'k':k, 'c':c} for k in params['norm_k'] for c in params['norm_c']]
        asm_perms = [{'k':p['k'],'c':p['c']} for p in norm_perms]
    else:
        norm_perms = [{'k':'0', 'c':'0'}]
        asm_perms = [{'k':'0', 'c':'0'}]

    for dir in ['raw', 'trim', 'norm', 'asm', 'asm_prune', 'remap', 'eval']:
        if not os.path.exists(params['out'] + '/' + dir):
            os.makedirs(params['out'] + '/' + dir)

    import_reads(fwd_fq, rev_fq, params)
    trim(norm_k_list, params)
    if norm_k_list and norm_c_list:
        norm_perms = normalise(norm_perms, params)
    else:
        norm_perms = None
    asms_paths_full = assemble(asm_perms, params)
    asms_paths = prune_assemblies(asms_paths_full, min_len, params)
    
    asms_names = {a: [r.id for r in SeqIO.parse(p, 'fasta')] for a, p in asms_paths.items()}
    asms_lens = lengths(asms_paths)
    asms_covs = map_to_assemblies(asms_paths, params)
    asms_gc = gc_content(asms_paths)
    
    if lca:
        lca_taxa = onecodex_assemblies(asms_paths, onecodex_api_key)
        metadata_summaries = marker_metadata(asms_paths, asms_lens, asms_gc, lca_taxa)
        asms_stats = dict(names=asms_names,
                          lens=asms_lens,
                          covs=asms_covs,
                          legend=metadata_summaries,
                          gc=asms_gc,
                          cpg=None)

        # print('\npaths', asms_paths)
        # print('\nasm_names', asms_names)
        # print(len(asms_names['060-660_r1_Cap1_F.norm_k21c1.asm_k']))
        # print(len(asms_names['060-660_r1_Cap1_F.norm_k21c10.asm_k']))
        # print('\nlengths', asms_lens)
        # print(len(asms_lens['060-660_r1_Cap1_F.norm_k21c1.asm_k']))
        # print(len(asms_lens['060-660_r1_Cap1_F.norm_k21c10.asm_k']))
        # print('\ncovs', asms_covs)
        # print(len(asms_covs['060-660_r1_Cap1_F.norm_k21c1.asm_k']))
        # print(len(asms_covs['060-660_r1_Cap1_F.norm_k21c10.asm_k']))
        # print('\ngc', asms_gc)
        # print(len(asms_gc['060-660_r1_Cap1_F.norm_k21c1.asm_k']))
        # print(len(asms_gc['060-660_r1_Cap1_F.norm_k21c10.asm_k']))
        # pprint.pprint(lca_taxa)
        # pprint.pprint(metadata_summaries)

    elif blast:
        blast_results = blast_assemblies(asms_paths, blast_db, blast_max_seqs, min_len)
        asms_stats = dict(names=asms_names,
                          lens=asms_lens,
                          covs=asms_covs,
                          blast_summary=blast_summary(blast_results, asms_covs), 
                          blast_superkingdoms=blast_superkingdoms(blast_results),
                          gc=gc_content(asms_paths),
                          cpg=None)

    else:
        asms_stats = dict(names=asms_names,
                          lens=asms_lens,
                          covs=asms_covs,
                          gc=asms_gc,
                          cpg=None)

    chart_url = plotly(asms_names, asms_stats, lca, blast, params)
    report(chart_url, start_time, time.time(), params)


argh.dispatch_command(main)