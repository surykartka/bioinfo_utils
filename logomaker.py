from Bio import SeqIO, Seq
import os, logomaker, pandas, numpy, sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['arial']
matplotlib.rcParams['font.size'] = 14

def get_ntDist(seq):
	ntDist = {}
	for nt in seq:
		if nt not in '+X':
			ntDist[nt] = ntDist.get(nt,0) + 1
	total = sum(ntDist.values())
	ntDist['G'] = (ntDist['G'] + ntDist['C']) / 2 / total
	ntDist['C'] = ntDist['G']
	ntDist['A'] = (ntDist['A'] + ntDist['T']) / 2 / total
	ntDist['T'] = ntDist['A']
	print(ntDist)
	return ntDist

def print_weblogo(file, start=None):
	seqs = {}
	for rec in SeqIO.parse(open(file),'fasta'):
	    r,seq = rec.description,str(rec.seq)
	    seqs[r] = seq
	ntDist = get_ntDist(''.join(seqs.values()))
	
	cmd = 'weblogo --first-index -%d --color-scheme classic -n 100 ' % (start)
	cmd += '--scale-width NO -F pdf --errorbars NO --yaxis 0.5 '
	cmd += '--composition "%s" ' % ntDist
	cmd += '< %s > weblogo.pdf' % file
	os.system(cmd)

def print_EDlogo(file, start=None):
	#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2489-3
	# EDLogo plot
	# start denotes the first position on x axis

	seqs = {}
	for rec in SeqIO.parse(open(file),'fasta'):
	    r,seq = rec.description,str(rec.seq)
	    seqs[r] = seq

	background_seq = ''.join(seqs.values())
	q_i = get_ntDist(background_seq)
	n_seqs = len(seqs)
	acgt = sorted(q_i.keys())
	mat = []
	for j in range(len(seq)):
		nts = [x[j] for x in seqs.values()]
		r_i = {i: numpy.log2(nts.count(i)/n_seqs/q_i[i]) for i in acgt}
		r_i_median = numpy.median(list(r_i.values()))
		r_i_corr = {i: r_i[i] - r_i_median for i in acgt} 
		## without corrections for median
		#r_i_corr = r_i
		vector = [r_i_corr[i] for i in acgt]
		mat.append(vector)
	
	if not start:
		start = 1
	end = len(seq)+start

	logo_df = pandas.DataFrame(mat, columns=acgt, index=range(start, end))
	blue='#234379';green='#4D8738';yellow='#DF902F';red='#AC1A24'

	logomaker.Logo(logo_df,
				   color_scheme='colorblind_safe', #{'A':green,'C':blue,'G':yellow,'T':red}
				   figsize=(7,2.5))
				   #center_values=True)
	plt.xlabel('Distance [nt]')
	plt.ylabel('Enrichment [bits]')
	plt.title('n = %d' % (n_seqs))
	#plt.xticks(range(-upLen, downLen, 5))
	plt.tight_layout()
	plt.savefig('TTS_EDlogo.pdf')


if __name__ == '__main__':
	file = sys.argv[1]
	print_EDlogo(file)
	#print_weblogo(file)