#!/usr/bin/env python3
#import pysam
import sys
import re
import pyfaidx
from collections import defaultdict
import datetime
from whatshap.args import HelpfulArgumentParser as ArgumentParser

#def zipped_pileup(filenames, chromosome):
	#samfiles = [pysam.AlignmentFile(f, "rb") for f in filenames]
	#iterators = [samfile.pileup(chromosome) for samfile in samfiles]
	#nexts = [None] * len(samfiles)
	#def fetch_next(i):
		#try:
			#nexts[i] = iterators[i].__next__()
		#except StopIteration:
			#nexts[i] = None
	#def all_none():
		#return sum(n is None for n in nexts) == len(nexts)
	#for i in range(len(nexts)):
		#fetch_next(i)
	#while not all_none():
		#y = [None] * len(samfiles)
		#position = min(col.pos for col in nexts if col is not None)
		#for i, col in enumerate(nexts):
			#if col is None:
				#continue
			#if col.pos != position:
				#assert col.pos > position
				#continue
			#y[i] = col
			#fetch_next(i)
		#yield position, y
	#for samfile in samfiles:
		#samfile.close()

#for position, columns in zipped_pileup(['pacbio/NA19238.chr1.bam','pacbio/NA19239.chr1.bam', 'pacbio/NA19240.chr1.bam'], 'chr1'):
	#if position %100 == 0:
		#print('Position:', position)
	#for i, column in enumerate(columns):
		#if column is None:
			##print('  column {}, None'.format(i))
			#pass
		#else:
			#bases = []
			#for pileupread in column.pileups:
				#if not pileupread.is_del and not pileupread.is_refskip:
					#bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
			##print('  column {}, coverage: {}, bases: {}'.format(i,column.n,''.join(bases)))
	#if position == 38606:
		#break
	
def main():
	parser = ArgumentParser(prog='simple-snp-caller.py', description=__doc__)
	parser.add_argument('--minabs', metavar='MIN_ABS', default=3, type=int,
		help='Minimum absolute ALT depth to call a SNP (default: %(default)s).')
	parser.add_argument('--minrel', metavar='MIN_REL', default=0.25, type=float,
		help='Minimum relative ALT depth to call a SNP (default: %(default)s).')
	parser.add_argument('--multi-allelics', default=False, action='store_true',
		help='Also output multi-allelic sites, if not given only the best ALT allele is reported (if unique).')
	parser.add_argument('--sample', metavar='SAMPLE', default=None, 
		help='Put this sample column into VCF (default: output sites-only VCF).')
	parser.add_argument('ref', metavar='REF', help='FASTA with reference genome')
	args = parser.parse_args()
	
	fasta = pyfaidx.Fasta(args.ref, as_raw=True)

	print('##fileformat=VCFv4.2')
	print('##fileDate={}'.format(datetime.datetime.now().strftime('%Y%m%d')))
	print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
	print('##FILTER=<ID=PASS,Description="All filters passed">')
	header_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
	if args.sample is not None:
		header_columns += ['FORMAT', args.sample]
	print(*header_columns, sep='\t')

	re_nucleotide = re.compile('[ACGTNacgtn]')
	re_indel = re.compile('[-\\+]([0-9]+)')
	re_ref = re.compile('[,\\.]')
	re_ignore = re.compile('([\\$\\*]|\\^.)')
	for line in sys.stdin:
		fields = line.split()
		if len(fields) == 4:
			continue
		assert len(fields) == 6, fields
		chromosome = fields[0]
		position = int(fields[1])
		pileup = fields[4]
		i = 0
		bases = defaultdict(int)
		#print('digesting pileup', pileup)
		n = 0
		ref = fasta[chromosome][position-1].upper()
		while i < len(pileup):
			m = re_nucleotide.match(pileup[i:])
			if m is not None:
				bases[pileup[i].upper()] += 1
				n += 1
				#print('  found base:', pileup[i])
				i += 1
				continue
			m = re_indel.match(pileup[i:])
			if m is not None:
				l = int(m.group(1))
				skip = (m.end()-m.start()) + l
				#print('  found indel:', pileup[i:i+skip])
				i += skip
				continue
			m = re_ref.match(pileup[i:])
			if m is not None:
				bases[ref] += 1
				n += 1
				#print('  found REF:', pileup[i:i+skip])
				i += 1
				continue
			m = re_ignore.match(pileup[i:])
			if m is not None:
				skip = (m.end()-m.start())
				#print('  found other things to ignore:', pileup[i:i+skip])
				i += skip
				continue
			assert False
		#print(bases)
		ref_count = bases[ref]
		alts = []
		for base, count in bases.items():
			if base == ref:
				continue
			if (count >= args.minabs) and (count / (count+ref_count) >= args.minrel):
				alts.append( (count, base) )
		alts.sort(reverse=True)
		if len(alts) > 0:
			columns = [chromosome, position, '.',  ref, '.', '.', 'PASS', '.']
			if args.sample is not None:
				columns += ['GT', '.']
			if args.multi_allelics:
				columns[4] = ','.join(base for count, base in alts)
			else:
				# Do we have two equally supported ALT alleles
				if len(alts) > 1 and (alts[0][0] == alts[1][0]):
					columns[4] = 'N'
				else:
					columns[4] = alts[0][1]
			print(*columns, sep='\t')

if __name__ == '__main__':
	main()
