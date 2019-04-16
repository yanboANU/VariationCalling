# This file contains wrapper functions to run various haplotype assembly tools

import time
import sys
import os
import fileIO
import subprocess
import shutil

# Runs (cmd) as a process, kills it after (timeout) seconds
# kills process if memory hits 8 GB / 8M ??
#def run_process(cmd, timeout=None, memlimit=7700000):
def run_process(cmd, timeout=None, memlimit=1000000000):
    cmd += ' 2>&1'
    print(cmd)
    # credit to roland smith for method of retrieving stdout and stderr: http://stackoverflow.com/questions/14059558/why-is-python-no-longer-waiting-for-os-system-to-finish
    t1 = time.time()
    if timeout == None:
        cmd = "ulimit -v {}; {}".format(memlimit, cmd)
        prog = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    else:
        timeout_cmd = 'ulimit -v {}; timeout {} {}'.format(memlimit, timeout, cmd)
        prog = subprocess.Popen(timeout_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = prog.communicate()
    t2 = time.time()
    runtime = t2 -t1
    print ("running time:", runtime)
    o = out.decode("ISO-8859-1")
    print(o)
    return runtime

#yanbo 
def run_whatsHap():
    return

#path_to_bin = /home/yulin/liyanbo/Tools/freebayes/bin/freebayes
def run_freebyes(path_to_bin, ref_file, bam_file, out_file):
    cmd = "{} --fasta-reference {} {} > {}".format(path_to_bin, ref_file, bam_file, out_file)
    return run_process(cmd,timeout)


def run_GATK():
     
    return run_process(cmd,timeout)

#bwa index my.fasta
#bwa aln [opts] my.fasta my.fastq > my.sai
#bwa samse my.fasta my.sai my.fastq > my.sam
#samtools view -S -b my.sam > my.bam
#samtools sort my.bam my-sorted
#samtools faidx my.fasta
#samtools mpileup -g -f my.fasta my-sorted1.bam my-sorted-2.bam my-sorted-n.bam > myraw.bcf
#bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf


def run_samtools():
    return run_process(cmd,timeout)





def run_hapcut2(path_to_bin, frag_file, vcf_file, output_file, timeout=None):
    cmd = "{} --fragments {} --vcf {} --output {}".format(path_to_bin, frag_file, vcf_file, output_file)
    return run_process(cmd,timeout)


def run_evaluator(path_to_bin, valid_master, out_file, chrID, stat_file, timeout=None):
    cmd = "{} --master {} --alg-phase {} --chr {} > {}".format(path_to_bin, valid_master, out_file, chrID, stat_file)
    #print cmd
    return run_process(cmd,timeout)

# a wrapper function for calling the HapCUT program
def run_hapcut(path_to_bin, frag_file, vcf_file, output_file, maxiter, maxcutiter, longreads=0, timeout=None):
    cmd = "{} --fragments {} --VCF {} --output {} --maxiter {} --maxcutiter {} --longreads {} --maxmem 32000".format(path_to_bin, frag_file, vcf_file, output_file, maxiter, maxcutiter, longreads)
    return run_process(cmd,timeout)

def run_refhap(path_to_jar, frag_file, output_file,timeout=None):
    cmd = "java -cp {} mpg.molgen.sih.main.SIH {} {}".format(path_to_jar, frag_file, output_file)
    return run_process(cmd,timeout)

#not to java, path to mpg.molgen.sih.main.SIH
def run_dgs(path_to_jar, frag_file, output_file,timeout=None):
    cmd = "java -cp {} mpg.molgen.sih.main.SIH -a DGS {} {}".format(path_to_jar, frag_file, output_file)
    return run_process(cmd,timeout)

def run_fasthare(path_to_jar, frag_file, output_file,timeout=None):
    cmd = "java -cp {} mpg.molgen.sih.main.SIH -a FastHare {} {}".format(path_to_jar, frag_file, output_file)
    return run_process(cmd,timeout)

# a wrapper function for calling the MixSIH program
def run_mixsih(path_to_bin, frag_file, vcf_file, output_file, mixsih_frag_file,timeout=None):

    num_snps  = fileIO.count_SNPs(vcf_file)

    # mixsih fragment file has no qual scores in the last column
    # it also has the num SNPs at top
    # we will simultaneously estimate average miscall rate and remove the qual scores
    total_qual = 0
    num_qual   = 0
    mixsih_frag_file_unsorted = mixsih_frag_file + ".unsorted"
    with open (mixsih_frag_file_unsorted, 'w') as mff:
        print(str(num_snps), file=mff)
        with open(frag_file, 'r') as ff:
            for line in ff:
                el = line.strip().split()
                if len(el) > 2:
                    # tally up the present qual scores
                    for q in el[-1]:
                        num_qual += 1
                        total_qual += 10**((ord(q)-33)/-10)

                    # remove qual string
                    el = el[:-1]
                    print(' '.join(el), file=mff)
                else:
                    print(line.strip(), file=mff)

    # mixsih requires we input a miscall rate, so we take the average
    miscall_rate = total_qual/num_qual
    # need to be sorted by value of third column...
    os.system("sort -n -k 3 {} > {}".format(mixsih_frag_file_unsorted, mixsih_frag_file))

    cmd = "{0} -a {1} {2} {3}.profile {3}".format(path_to_bin, miscall_rate, mixsih_frag_file, output_file)
    return run_process(cmd,timeout)

# a wrapper function for calling the Haptree program
def run_haptree(path_to_bin, frag_file, vcf_file, output_file, haptree_vcf_file,timeout=None):

    # haptree vcf requires no header
    with open (haptree_vcf_file, 'w') as ht_vcf:
        with open(vcf_file, 'r') as vcf:
            for line in vcf:
                if line[:1] != '#':
                    print(line.strip(), file=ht_vcf)

    haptree_output_dir = output_file + ".dir"

    cmd = "{} {} {} {}".format(path_to_bin, frag_file, haptree_vcf_file, haptree_output_dir)
    haptree_solution = os.path.join(haptree_output_dir,"HapTreeSolution")
    if os.path.isfile(haptree_solution):
        shutil.move(haptree_solution, output_file)
        os.rmdir(haptree_output_dir)

    return run_process(cmd,timeout)


