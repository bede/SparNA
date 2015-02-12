#!/usr/bin/env python

# Author: Bede Constantinides
# Python (2.7+) pipeline for paired-end HepC assembly developed during a placement at PHE
# Please seek permission prior to use or distribution

# TODO
# | include reference_path inside paths?
# | if no reads are aligned in the blast step, try blasting every single read
# | itop appeasing Insanely Bad Format - keep reads interleaved except as required
# | use khmer for deinterleaving (split-paired-reads.py)
# | sort out assess_coverage() (currently disused)
# | add minimum similarity threshold for reference selection
# | mauve/nucmer integration for when a contiguous assembly is unavailable
# | parellelise normalisation and assembly for speed
# | report on trimming, %remapped
# | increase khmer table size
# | improve command line interface
# | check exit code of every call to os.system?
# | pep8

# DEPENDENCIES
# | python packages:
# |    argh, biopython, envoy, khmer
# | others, expected inside $PATH:
# |    blast, samtools, seqtk, spades, quast
# | others, bundled inside res/ directory:
# |    trimmomatic, fastq_deinterleave
# | others, bundled and requiring compilation: segemehl

# USAGE: ./pipeline.py --threads 12 --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-c-list 5 --asm-k-list 33 --multiple-samples --in-dir /path/to/fastqs --out-dir /path/to/output
# Input fastq filenames should have an extension and a signature to allow identification of forward and reverse reads

# min_cov
# segemehl_sensitivity_pc
# min_depth

from __future__ import division, print_function
import os
import sys
import time
import argh
import envoy
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

def list_fastqs(fwd_reads_sig, rev_reads_sig, paths):
    print('Identifying input... ')
    fastqs = {'f':[], 'r':[]}
    for fastq in os.listdir(paths['in']):
        if fastq.endswith('.fastq') or fastq.endswith('.fastq'):
            if fwd_reads_sig in fastq:
                fastqs['f'].append(paths['in'] + '/' + fastq)
            elif rev_reads_sig in fastq:
                fastqs['r'].append(paths['in'] + '/' + fastq)
    fastq_pairs = zip(fastqs['f'], fastqs['r'])
    fastq_pairs = {os.path.splitext(p[0].replace(fwd_reads_sig,''))[0]: p for p in fastq_pairs}
    print('\tDone') if fastq_pairs else sys.exit('ERR_READS')
    return fastqs, fastq_pairs
    
def import_reads(multiple_samples, fastqs, fastq_pair, paths, i=1):
    print('Importing reads... ')
    cmd_vars = {
     'i':str(i),
     'fq_pair_f':fastq_pair[0],
     'fq_pair_r':fastq_pair[1],
     'path_out':paths['out'],
     'fastqs_f':' '.join(fastqs['f']),
     'fastqs_r':' '.join(fastqs['r'])}
    if multiple_samples:
        cmd_import = (
         'cp {fq_pair_f} {path_out}/merge/{i}.raw.r1.fastq && '
         'cp {fq_pair_r} {path_out}/merge/{i}.raw.r2.fastq && '
         .format(**cmd_vars))
    else:
        cmd_import = (
         'cat {fastqs_f} > {path_out}/merge/{i}.raw.r1.fastq && '
         'cat {fastqs_r} > {path_out}/merge/{i}.raw.r2.fastq && '
         'interleave-reads.py {path_out}/merge/{i}.raw.r1.fastq '
         .format(**cmd_vars))
    cmd_import += (
     'interleave-reads.py {path_out}/merge/{i}.raw.r1.fastq '
     '{path_out}/merge/{i}.raw.r2.fastq > {path_out}/merge/{i}.raw.r12.fastq'
     .format(**cmd_vars))
    cmd_import = os.system(cmd_import)
    print('\tDone') if cmd_import == 0 else sys.exit('ERR_IMPORT')


def count_reads(paths, i=1):
    print('Counting reads...')
    cmd_count = ('wc -l {path_out}/merge/{i}.raw.r12.fastq'
    .format(i=str(i),
            path_out=paths['out']))
    cmd_count = envoy.run(cmd_count)
    n_reads = int(cmd_count.std_out.replace(' ','').split('/')[0])/4
    print('\tDone') if cmd_count.status_code == 0 else sys.exit('ERR_COUNT')
    return n_reads

def sample_reads(n_reads, paths, i=1):
    print('Sampling reads...')
    n_reads_sample = n_reads/100
    cmd_sample = (
     'cat {path_out}/merge/{i}.raw.r12.fastq | seqtk sample - '
     '{n_reads_sample} | seqtk seq -a - > {path_out}/sample/{i}.sample.fasta'
     .format(i=str(i),
             path_out=paths['out'],
             n_reads_sample=n_reads_sample))
    cmd_sample = os.system(cmd_sample)
    print('\tDone') if cmd_sample == 0 else sys.exit('ERR_SAMPLE')
    return n_reads_sample

def blast_references(paths, threads, i=1):
    print('BLASTing reference sequences...')
    if not os.path.exists(paths['pipe'] + '/res/hcv_db/db.nhr'):
        cmd_blastn_index = (
         'makeblastdb -dbtype nucl -input_type fasta '
         '-in {path_pipe}/res/hcv_db/db.fasta -title db'
         .format(path_pipe=paths['pipe']))
        cmd_blastn_index = os.system(cmd_blastn_index)
    cmd_blastn = NcbiblastnCommandline(
        query = paths['out'] + '/sample/' + str(i) + '.sample.fasta', num_alignments = 1,
        db = paths['pipe'] + '/res/hcv_db/db.fasta', evalue = 1e-4, outfmt = 7,
        out = paths['out'] + '/blast/' + str(i) + '.blast.tsv', 
        num_threads = threads)
    cmd_blastn()
    print('\tDone')

def choose_reference(paths, i=1):
    print('Choosing reference sequence...')
    accession_freqs = {}
    with open(paths['out'] + '/blast/' + str(i) + '.blast.tsv', 'r') as blast_out:
        for line in blast_out:
            if not line.startswith('#'):
                accession = line.split('\t')[1]
                if accession in accession_freqs.keys():
                    accession_freqs[accession] += 1
                else: accession_freqs[accession] = 1
    if accession_freqs:
        reference_found = True
        top_accession = max(accession_freqs, key=accession_freqs.get)
    else:
        reference_found = False
        top_accession = None
        print('\tWARNING: failed to identify a similar reference sequence')
    print('\tDone')
    return reference_found, top_accession

def extract_reference(top_accession, paths, i=1):
    print('\tExtracting ' + top_accession + '...')
    reference = ''
    with open(paths['pipe'] + '/res/hcv_db/db.fasta', 'r') as references_fa:
        inside_best_reference = False
        for line in references_fa:
            if line.startswith('>'):
                if top_accession in line:
                    inside_best_reference = True
                else: inside_best_reference = False
            elif inside_best_reference:
                reference += line.strip()
    reference_len = len(reference)
    reference_path = paths['out'] + '/ref/' + str(i) + '.ref.fasta'
    with open(reference_path, 'w') as reference_fa:
        reference_fa.write('>' + top_accession + '\n' + reference)
    print('\tDone')
    return reference_path, reference_len

def genotype(n_reads, paths, i=1):
    print('Genotyping...')
    genotype_freqs = {}
    with open(paths['out'] + '/blast/' + str(i) + '.blast.tsv', 'r') as blast_out:
        for line in blast_out:
            if not line.startswith('#'):
                genotype = line.split('\t')[1].split('_')[1].split('.')[0]
                if genotype in genotype_freqs.keys():
                    genotype_freqs[genotype] += 1
                else: genotype_freqs[genotype] = 1
    top_genotype = max(genotype_freqs, key=genotype_freqs.get) if genotype_freqs else None
    genotype_props = {k: genotype_freqs[k]/n_reads*100 for k in genotype_freqs.keys()}
    genotype_props_pc = {k: round(genotype_freqs[k]/n_reads*1e4, 2) for k in genotype_freqs.keys()}
    genotype_props_pc_sorted = []
    for genotype, proportion in reversed(sorted(genotype_props_pc.items(), key=lambda(k,v):(v,k))):
        record = (genotype + ': ' + str(proportion) + '% (' + str(genotype_freqs[genotype]) + ')')
        genotype_props_pc_sorted.append(record)
        print('\t' + record)
    with open(paths['out'] + '/blast/' + str(i) + '.genotypes.txt', 'w') as genotypes_file:
        for item in genotype_props_pc_sorted:
             genotypes_file.write(item + '\n')
    print('\tDone')
    return top_genotype

def map_reads(reference_path, paths, threads, i=1):
    print('Aligning with Segemehl... ')
    cmd_vars = {
     'i':str(i),
     'path_pipe':paths['pipe'],
     'path_out':paths['out'],
     'reference_path':reference_path,
     'threads':threads}
    cmds_map = [
     '{path_pipe}/res/segemehl/segemehl.x -d {reference_path} -x {reference_path}.idx',
     '{path_pipe}/res/segemehl/segemehl.x -d {reference_path} -x {reference_path}.idx -q '
     '{path_out}/merge/{i}.raw.r12.fastq --threads {threads} -A 60 > '
     '{path_out}/map/{i}.segemehl_mapped.sam',
     'samtools view -bS {path_out}/map/{i}.segemehl_mapped.sam | samtools sort - '
     '{path_out}/map/{i}.segemehl_mapped',
     'samtools index {path_out}/map/{i}.segemehl_mapped.bam',
     'samtools mpileup -d 1000 -f {reference_path} {path_out}/map/{i}.segemehl_mapped.bam > '
     '{path_out}/map/{i}.segemehl_mapped.pile',
     'samtools mpileup -ud 1000 -f {reference_path} /map/{i}.segemehl_mapped.bam | '
     'bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -a - | fasta_formatter -o '
     '{path_out}/map/{i}.consensus.fasta']
    for i, cmd_map in enumerate(cmds_map):
        cmd_map = os.system(cmd_map.format(**cmd_vars))
        print('\tDone (step ' + str(i) + ')') if cmd_map == 0 else sys.exit('ERR_MAP')

def assess_coverage(reference_len, paths, i=1):
    print('Identifying low coverage regions... ')
    min_depth = 1
    min_coverage = 0.9
    depths = {}
    uncovered_sites = []
    with open(paths['out'] + '/map/' + str(i) + '.segemehl_mapped.pile', 'r') as pileup:
        bases_covered = 0
        for line in pileup:
            site = int(line.split('\t')[1])
            depths[site] = int(line.split('\t')[3])
            if depths[site] < min_depth:
                uncovered_sites.append(site)
    uncovered_region = 0
    uncovered_regions = []
    last_uncovered_site = 0
    largest_uncovered_region = 0
    for uncovered_site in uncovered_sites:
        if uncovered_site == last_uncovered_site + 1:
            uncovered_region += 1
            if uncovered_region > largest_uncovered_region:
                largest_uncovered_region = uncovered_region
        else:
            if uncovered_region > 0:
                uncovered_regions.append(uncovered_region)
            uncovered_region = 1
        last_uncovered_site = uncovered_site

    print('\tUncovered sites: ' + str(len(uncovered_sites)))
    print('\tUncovered regions: ' + str(len(uncovered_regions)))
    if not uncovered_sites:
        print('\tAll reference bases covered!')
    elif len(uncovered_sites) < (1-min_coverage)*reference_len:
        print('\tReference coverage safely above threshold')
    else:
        print('\tReference coverage below threshold')
    print('\tDone')

def trim(paths, i=1):
    print('Trimming... ')
    cmd_trim = (
     'java -jar {path_pipe}/res/trimmomatic-0.32.jar PE '
     '{path_out}/merge/{i}.raw.r1.fastq {path_out}/merge/{i}.raw.r2.fastq '
     '{path_out}/trim/{i}.trim.r1_pe.fastq {path_out}/trim/{i}.trim.r1_se.fastq '
     '{path_out}/trim/{i}.trim.r2_pe.fastq {path_out}/trim/{i}.trim.r2_se.fastq '
     'ILLUMINACLIP:{path_pipe}/res/illumina_adapters.fa:2:30:10 MINLEN:25'
     .format(i=str(i),
             path_pipe=paths['pipe'],
             path_out=paths['out']))
    cmd_trim_pp = (
     'cat {path_out}/trim/{i}.trim.r1_se.fastq {path_out}/trim/{i}.trim.r2_se.fastq > '
     '{path_out}/trim/{i}.trim.se.fastq '
     '&& interleave-reads.py {path_out}/trim/{i}.trim.r1_pe.fastq '
     '{path_out}/trim/{i}.trim.r2_pe.fastq > {path_out}/trim/{i}.trim.r12_pe.fastq'
     .format(i=str(i), 
             path_pipe=paths['pipe'],
             path_out=paths['out']))
    cmd_trim = envoy.run(cmd_trim)
    print(cmd_trim.std_err)
    cmd_trim_stats = ''.join(cmd_trim.std_err).split('\n')[25]
    cmd_trim_pp = os.system(cmd_trim_pp)
    print('\tDone') if cmd_trim.status_code == 0 else sys.exit('ERR_TRIM')

def normalise(norm_k_list, norm_c_list, paths, i=1):
    print('Normalising... ')
    ks = norm_k_list.split(',')
    cs = norm_c_list.split(',')
    norm_perms = [{'k':k, 'c':c} for k in ks for c in cs]
    for norm_perm in norm_perms:
        cmd_norm = (
         'normalize-by-median.py -C {c} -k {k} -N 1 -x 1e9 -p '
         '{path_out}/trim/{i}.trim.r12_pe.fastq '
         '-o {path_out}/norm/{i}.norm_k{k}c{c}.r12_pe.fastq '
         '&& normalize-by-median.py -C {c} -k {k} -N 1 -x 1e9 '
         '{path_out}/trim/{i}.trim.se.fastq '
         '-o {path_out}/norm/{i}.norm_k{k}c{c}.se.fastq '
         '&& {path_pipe}/res/fastq_deinterleave '
         '{path_out}/norm/{i}.norm_k{k}c{c}.r12_pe.fastq '
         '{path_out}/norm/{i}.norm_k{k}c{c}.r1_pe.fastq '
         '{path_out}/norm/{i}.norm_k{k}c{c}.r2_pe.fastq '
         '&& cat {path_out}/norm/{i}.norm_k{k}c{c}.r12_pe.fastq '
         '{path_out}/norm/{i}.norm_k{k}c{c}.se.fastq > '
         '{path_out}/norm/{i}.norm_k{k}c{c}.pe_and_se.fastq'
         .format(i=str(i),
                 k=str(norm_perm['k']),
                 c=str(norm_perm['c']),
                 path_pipe=paths['pipe'],
                 path_out=paths['out']))
        cmd_norm = os.system(cmd_norm)
        print('\tDone (k=' + k + ', c=' + c + ')') if cmd_norm == 0 else sys.exit('ERR_NORM')
    return norm_perms

def assemble(norm_perms, asm_k_list, asm_untrusted_contigs, reference_found, paths, threads, i=1):
    print('Assembling... ')
    if reference_found and asm_untrusted_contigs:
        asm_perms = [{'k':p['k'],'c':p['c'],'uc':uc} for p in norm_perms for uc in [1, 0]]
    else:
        asm_perms = [{'k':p['k'],'c':p['c'],'uc':uc} for p in norm_perms for uc in [0]]
    for asm_perm in asm_perms:
        cmd_vars = {
         'i':str(i),
         'k':str(asm_perm['k']),
         'c':str(asm_perm['c']),
         'uc':str(asm_perm['uc']),
         'asm_k_list':asm_k_list,
         'path_out':paths['out'],
         'threads':threads}
        cmd_asm = (
         'spades.py -m 8 -t {threads} -k {asm_k_list} '
         '--pe1-1 {path_out}/norm/{i}.norm_k{k}c{c}.r1_pe.fastq '
         '--pe1-2 {path_out}/norm/{i}.norm_k{k}c{c}.r2_pe.fastq '
         '--s1 {path_out}/norm/{i}.norm_k{k}c{c}.se.fastq '
         '-o {path_out}/asm/{i}.norm_k{k}c{c}.asm_k{asm_k_list}.uc{uc} --careful'
         .format(**cmd_vars))
        if asm_perm['uc']:
            cmd_asm += ' --untrusted-contigs ' + paths['out'] + '/ref/' + str(i) + '.ref.fasta'
        cmd_asm = envoy.run(cmd_asm)
        print(cmd_asm.std_out, cmd_asm.std_err)
        print('\tDone (k=' + asm_k_list + ')') if cmd_asm.status_code == 0 else sys.exit('ERR_ASM')

def evaluate_assemblies(reference_found, paths, threads, i=1):
    print('Comparing assemblies... ')
    asm_dirs = (
     [paths['out'] + '/asm/' + dir + '/contigs.fasta' for dir in os.listdir(paths['out'] + '/asm')
     if not dir.startswith('.')])
    cmd_vars = {
     'i':str(i),
     'asm_dirs':' '.join(asm_dirs),
     'path_out':paths['out'],
     'threads':threads}
    eval_cmd = ('quast.py {asm_dirs} -o {path_out}/eval/{i} --threads {threads}'.format(**cmd_vars))
    if reference_found:
        eval_cmd += ' -R {path_out}/ref/{i}.ref.fasta'.format(**cmd_vars)
    eval_cmd += ' &> /dev/null'
    os.system(eval_cmd)

def remap_reads():
    pass

def report(paths, i):
    os.makedirs(paths['out'] + '/eval/summary/')
    cmd_vars = {
     'i':str(i),
     'path_out':paths['out']}
    cmd_report = (
     'cp -R {path_out}/eval/{i}/report.html {path_out}/eval/{i}/transposed_report.tsv '
     '{path_out}/eval/{i}/report_html_aux {path_out}/eval/sumary/'.format(**cmd_vars))
    cmd_report = os.system(cmd_report)
    print('\tQUAST report: ' + paths['out'] + '/eval/summary/')

def main(in_dir=None, out_dir=None, fwd_reads_sig=None, rev_reads_sig=None, norm_k_list=None,
    norm_c_list=None, asm_k_list=None, asm_untrusted_contigs=False, multiple_samples=False,
    threads=1):
    paths = {
     'in':in_dir,
     'pipe':os.path.dirname(os.path.realpath(__file__)),
     'out':out_dir + '/run_' + str(int(time.time()))}
    job_dirs = ['merge', 'sample', 'blast', 'ref', 'map', 'trim', 'norm', 'asm', 'remap', 'eval']
    for dir in job_dirs:
        os.makedirs(paths['out'] + '/' + dir)

    fastqs, fastq_pairs = list_fastqs(fwd_reads_sig, rev_reads_sig, paths)
    for i, fastq_pair in enumerate(fastq_pairs, start=1):
        import_reads(multiple_samples, fastqs, fastq_pairs[fastq_pair], paths, i)
        n_reads = count_reads(paths, i)
        sample_reads(n_reads, paths, i)
        blast_references(paths, threads, i)
        reference_found, top_accession = choose_reference(paths, i)
        if reference_found:
            reference_path, reference_len = extract_reference(top_accession, paths, i)
            top_genotype = genotype(n_reads, paths, i)
            map_reads(reference_path, paths, i)
            assess_coverage(reference_len, paths, i)
        trim(paths, i)
        assemble(normalise(norm_k_list, norm_c_list, paths, i), asm_k_list, asm_untrusted_contigs,
                 reference_found, paths, threads, i)
        evaluate_assemblies(reference_found, paths, threads, i)
        remap_reads()
    report(paths, i)

argh.dispatch_command(main)