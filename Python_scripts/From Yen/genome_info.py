"""
This file contain all genome information from yeast to Drosophila

"""
#
# constant for yeast genome 
#
tss_total = 4792

nucleotides = ["A", "C", "G", "T"]

#
# amino acid codon
#
dict_codon2aa = {

	"ATT" : ("I", "Isoleucine"),
	"ATA" : ("I", "Isoleucine"),
	"ATC" : ("I", "Isoleucine"),
	
	"CTT" : ("L", "Leucine"),
	"CTC" : ("L", "Leucine"),
	"CTA" : ("L", "Leucine"),
	"CTG" : ("L", "Leucine"),
	"TTA" : ("L", "Leucine"),
	"TTG" : ("L", "Leucine"),
	
	"GTT" : ("V", "Valine"),
	"GTC" : ("V", "Valine"),
	"GTA" : ("V", "Valine"),
	"GTG" : ("V", "Valine"),
	
	"TTT" : ("F", "Phenylalanine"),
	"TTC" : ("F", "Phenylalanine"),
	
	"ATG" : ("M", "Methionine"),
	
	"TGT" : ("C", "Cysteine"),
	"TGC" : ("C", "Cysteine"),
	
	"GCT" : ("A", "Alanine"),
	"GCC" : ("A", "Alanine"),
	"GCA" : ("A", "Alanine"),
	"GCG" : ("A", "Alanine"),
	
	"GGT" : ("G", "Glycine"),
	"GGC" : ("G", "Glycine"),
	"GGA" : ("G", "Glycine"),
	"GGG" : ("G", "Glycine"),
	
	"CCT" : ("P", "Proline"),
	"CCC" : ("P", "Proline"),
	"CCA" : ("P", "Proline"),
	"CCG" : ("P", "Proline"),
	
	"ACT" : ("T", "Threonine"),
	"ACC" : ("T", "Threonine"),
	"ACA" : ("T", "Threonine"),
	"ACG" : ("T", "Threonine"),
	
	"TCT" : ("S", "Serine"),
	"TCC" : ("S", "Serine"),
	"TCA" : ("S", "Serine"),
	"TCG" : ("S", "Serine"),
	"AGT" : ("S", "Serine"),
	"AGC" : ("S", "Serine"),
	
	"TAT" : ("Y", "Tyrosine"),
	"TAC" : ("Y", "Tyrosine"),
	
	"TGG" : ("W", "Tryptophan"),
	
	"CAA" : ("Q", "Glutamine"),
	"CAG" : ("Q", "Glutamine"),
	
	"AAT" : ("N", "Asparagine"), 
	"AAC" : ("N", "Asparagine"),
	
	"CAT" : ("H", "Histidine"), 
	"CAC" : ("H", "Histidine"),
	
	"GAA" : ("E", "Glutamic acid"),
	"GAG" : ("E", "Glutamic acid"),
	
	"GAT" : ("D", "Aspartic acid"),
	"GAC" : ("D", "Aspartic acid"),
	
	"AAA" : ("K", "Lysine"),
	"AAG" : ("K", "Lysine"),
	
	"CGT" : ("R", "Arginine"),
	"CGC" : ("R", "Arginine"),
	"CGA" : ("R", "Arginine"),
	"CGG" : ("R", "Arginine"),
	"AGA" : ("R", "Arginine"),
	"AGG" : ("R", "Arginine"),
	
	"TAA" : ("STOP", "Stop codons"),
	"TAG" : ("STOP", "Stop codons"),
	"TGA" : ("STOP", "Stop codons")

	}

dict_aa2codon = {

	"I" : ("ATT", "ATC", "ATA"),
	"L" : ("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"),
	"V" : ("GTT", "GTC", "GTA", "GTG"),
	"F" : ("TTT", "TTC"),
	"M" : ("ATG"),
	"C" : ("TGT", "TGC"),
	"A" : ("GCT", "GCC", "GCA", "GCG"),
	"G" : ("GGT", "GGC", "GGA", "GGG"),
	"P" : ("CCT", "CCC", "CCA", "CCG"),
	"T" : ("ACT", "ACC", "ACG", "ACA"),
	"S" : ("TCT", "TCC", "TCG", "TCA", "AGT", "AGC"),
	"Y" : ("TAT", "TAC"),
	"W" : ("TGG"),
	"Q" : ("CAA", "CAG"),
	"N" : ("AAT", "AAC"),
	"H" : ("CAT", "CAC"),
	"E" : ("GAA", "GAG"),
	"D" : ("GAT", "GAC"),
	"K" : ("AAA", "AAG"),
	"R" : ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
	"STOP" : ("TAA", "TAG", "TGA")

	}
dict_aa2name = {

	"I" : "Isoleucine",
	"L" : "Leucine",
	"V" : "Valine",
	"F" : "Phenylalanine",
	"M" : "Methionine",
	"C" : "Cysteine",
	"A" : "Alanine",
	"G" : "Glycine",
	"P" : "Proline",
	"T" : "Threonine",
	"S" : "Serine",
	"Y" : "Tyrosine",
	"W" : "Tryptophan",
	"Q" : "Glutamine",
	"N" : "Asparagine",
	"H" : "Histidine",
	"E" : "Glutamic acid",
	"D" : "Aspartic acid",
	"K" : "Lysine",
	"R" : "Arginine",
	"STOP" : "Stop condons"
}
# the second letter of yeast gene id corresponds to chromosome numbers
dict_id_chrom = {

	"A": "chr01", "B": "chr02", "C": "chr03", "D": "chr04", 
	
	"E": "chr05", "F": "chr06", "G": "chr07", "H": "chr08", 
	
	"I": "chr09", "J": "chr10",  "K": "chr11", "L": "chr12", 
							
	"M": "chr13", "N": "chr14",  "O": "chr15", "P": "chr16"							
	
	}

# degenerated nucleotide IUB-Code:
IUB_CODE = {
	"A" : ["A"], 
	"C" : ["C"], 
	"G" : ["G"], 
	"T" : ["T"], 
	"R" : ["A", "G"], 
	"Y" : ["C", "T"], 
	"K" : ["G", "T"], 
	"M" : ["A", "C"], 
	"S" : ["G", "C"], 
	"W" : ["A", "T"], 
	"B" : ["C", "G", "T"], 
	"D" : ["A", "G", "T"], 
	"H" : ["A", "C", "T"], 
	"V" : ["A", "C", "G"], 
	"N" : ["A", "C", "G", "T"]
	}

# roman chromosome name to normal name
dict_roman2chrom = {

	"I": "chr1", 
	"II": "chr2", 
	"III": "chr3",
	"IV": "chr4", 
	"V": "chr5", 
	"VI": "chr6", 
	"VII": "chr7", 
	"VIII": "chr8", 
	"IX": "chr9", 
	"X": "chr10", 
	"XI": "chr11", 
	"XII": "chr12", 
	"XIII": "chr13", 
	"XIV": "chr14", 
	"XV": "chr15", 
	"XVI": "chr16", 
	"mitochondrion": "chrMito", 
	"plasmid": "2-micron"

	}

#
# dictionary related to nucleotide to color-coded convertion
#	

## nucleotide to color-code dictionary
dict_barcode = {			

		"AA": "0", "AC": "1", "AG": "2", "AT": "3", 

		"CA": "1", "CC": "0", "CG": "3", "CT": "2", 

		"GA": "2", "GC": "3", "GG": "0", "GT": "1", 

		"TA": "3", "TC": "2", "TG": "1", "TT": "0"

		}

dict_nuc2color = {

		"A0": "A", "A1": "C", "A2": "G", "A3": "T", 

		"C1": "A", "C0": "C", "C3": "G", "C2": "T", 

		"G2": "A", "G3": "C", "G0": "G", "G1": "T", 

		"T3": "A", "T2": "C", "T1": "G", "T0": "T"

		}



#E. coli genome dinucleotide frequency
di_freq_ecoli = {"AA": 5.5, "AC": 5.53, "GT": 5.51, "AG": 5.13, "CC": 4.99, "CA": 7.01, "CG": 7.47, "TT": 5.52, "GG": 4.96, "GC": 8.27, "AT": 6.68, "GA": 5.76, "TG": 6.95, "TA": 4.57, "TC": 5.76, "CT": 5.09}

#yeast genome dinucleotide frequency
di_freq_yeast = {"AA": 8.0, "AC": 5.26, "GT": 5.24, "AG": 5.83, "CC": 3.32, "CA": 6.47, "CG": 2.92, "TT": 7.95, "GG": 3.3, "GC": 3.73, "AT": 9.02, "GA": 6.22, "TG": 6.44, "TA": 7.43, "TC": 6.2, "CT": 5.8}

#drosophila mel dinucleotide frequency
di_freq_dmel_het = {"AA": 6.82, "AC": 4.71, "GT": 4.76, "AG": 4.91, "CC": 2.97, "CA": 5.51, "CG": 3.23, "TT": 6.9, "GG": 2.96, "GC": 3.96, "AT": 7.5, "GA": 5.01, "TG": 5.59, "TA": 6.6, "TC": 5.05, "CT": 4.98}

#yeast trinucleotide frequency
tri_freq_yeast = {"AAA" :2.71,"AAC" :1.88,"AAG" :2.26,"AAT" :3.08,"ACA" :1.71,"ACC" :1.21,"ACG" :0.91,"ACT" :1.57,"AGA" :1.98,"AGC" :1.19,"AGG" :1.19,"AGT" :1.58,"ATA" :2.35,"ATC" :1.84,"ATG" :1.9,"ATT" :3.07,"CAA" :2.41,"CAC" :1.06,"CAG" :1.31,"CAT" :1.91,"CCA" :1.55,"CCC" :0.59,"CCG" :0.61,"CCT" :1.19,"CGA" :0.95,"CGC" :0.57,"CGG" :0.6,"CGT" :0.91,"CTA" :1.34,"CTC" :1.1,"CTG" :1.3,"CTT" :2.24,"GAA" :2.48,"GAC" :1.01,"GAG" :1.1,"GAT" :1.83,"GCA" :1.29,"GCC" :0.82,"GCG" :0.57,"GCT" :1.19,"GGA" :1.32,"GGC" :0.82,"GGG" :0.57,"GGT" :1.2,"GTA" :1.48,"GTC" :1.01,"GTG" :1.05,"GTT" :1.87,"TAA" :2.33,"TAC" :1.48,"TAG" :1.34,"TAT" :2.35,"TCA" :2.1,"TCC" :1.32,"TCG" :0.95,"TCT" :1.96,"TGA" :2.1,"TGC" :1.29,"TGG" :1.54,"TGT" :1.69,"TTA" :2.33,"TTC" :2.46,"TTG" :2.39,"TTT" :2.7}

#drosophila trinucleotide including mitochondira and U
tri_freq_dmel = {"AAA" :2.67,"AAC" :1.79,"AAG" :1.82,"AAT" :2.95,"ACA" :1.78,"ACC" :1.08,"ACG" :0.97,"ACT" :1.52,"AGA" :1.43,"AGC" :1.54,"AGG" :1.12,"AGT" :1.52,"ATA" :2.01,"ATC" :1.51,"ATG" :1.8,"ATT" :2.95,"CAA" :2.31,"CAC" :1.34,"CAG" :1.57,"CAT" :1.8,"CCA" :1.73,"CCC" :0.81,"CCG" :0.98,"CCT" :1.12,"CGA" :1.29,"CGC" :1.05,"CGG" :0.98,"CGT" :0.98,"CTA" :1.04,"CTC" :1.19,"CTG" :1.57,"CTT" :1.82,"GAA" :1.99,"GAC" :1.06,"GAG" :1.18,"GAT" :1.51,"GCA" :1.84,"GCC" :1.37,"GCG" :1.05,"GCT" :1.53,"GGA" :1.37,"GGC" :1.37,"GGG" :0.81,"GGT" :1.08,"GTA" :1.18,"GTC" :1.06,"GTG" :1.34,"GTT" :1.79,"TAA" :2.26,"TAC" :1.19,"TAG" :1.04,"TAT" :2.01,"TCA" :1.64,"TCC" :1.38,"TCG" :1.29,"TCT" :1.43,"TGA" :1.64,"TGC" :1.84,"TGG" :1.73,"TGT" :1.78,"TTA" :2.26,"TTC" :1.99,"TTG" :2.31,"TTT" :2.67}

#drosophila trinucleotide frequency in het regions
tri_freq_dmel_het = {"AAA" :2.76,"AAC" :1.87,"AAG" :1.96,"AAT" :3.12,"ACA" :1.81,"ACC" :1.12,"ACG" :0.99,"ACT" :1.58,"AGA" :1.69,"AGC" :1.36,"AGG" :1.12,"AGT" :1.57,"ATA" :2.33,"ATC" :1.56,"ATG" :1.69,"ATT" :3.12,"CAA" :2.27,"CAC" :1.16,"CAG" :1.37,"CAT" :1.68,"CCA" :1.45,"CCC" :0.67,"CCG" :0.81,"CCT" :1.14,"CGA" :1.18,"CGC" :0.84,"CGG" :0.8,"CGT" :1.01,"CTA" :1.32,"CTC" :1.13,"CTG" :1.4,"CTT" :2.01,"GAA" :2.07,"GAC" :1.14,"GAG" :1.12,"GAT" :1.56,"GCA" :1.46,"GCC" :1.01,"GCG" :0.84,"GCT" :1.38,"GGA" :1.26,"GGC" :1.01,"GGG" :0.66,"GGT" :1.12,"GTA" :1.37,"GTC" :1.14,"GTG" :1.17,"GTT" :1.9,"TAA" :2.61,"TAC" :1.36,"TAG" :1.32,"TAT" :2.34,"TCA" :1.71,"TCC" :1.27,"TCG" :1.19,"TCT" :1.73,"TGA" :1.74,"TGC" :1.48,"TGG" :1.47,"TGT" :1.85,"TTA" :2.61,"TTC" :2.1,"TTG" :2.32,"TTT" :2.78}

#chromsome length dictionary: yeast
#chr_len_yeast = { "chr13": 924429 }

chr_len_yeast = {
	"chr1" :230208, 
	"chr2" :813178, 
	"chr3" :316617, 
	"chr4" :1531918, 
	"chr5" :576869, 
	"chr6" :270148, 
	"chr7" :1090946, 
	"chr8" :562643, 
	"chr9" :439885,
	"chr10" :745745, 
	"chr11" :666454, 
	"chr12" :1078175, 
	"chr13" :924429, 
	"chr14" :784333, 
	"chr15" :1091290, 
	"chr16" :948062
} #, "chrMito" :85779, "2-micron" :6318}

"""
#chr_len_yeast = { "chr13": 924429 }
chr_len_yeast = {
	"chr01" :230208, 
	"chr02" :813178, 
	"chr03" :316617, 
	"chr04" :1531918, 
	"chr05" :576869, 
	"chr06" :270148, 
	"chr07" :1090946, 
	"chr08" :562643, 
	"chr09" :439885, 
	"chr10" :745745, 
	"chr11" :666454, 
	"chr12" :1078175, 
	"chr13" :924429, 
	"chr14" :784333, 
	"chr15" :1091290, 
	"chr16" :948062
} #, "chrMito" :85779, "2-micron" :6318}
"""

#chromsome length dictionary: pombe: 
chr_len_sp11 = {
	"chr01":	5579133,
	"chr02":	4539804,	
	"chr03":	2452883
} 

# 2009 new genome length
chr_len_ye_09 = {
	"chr1":	230208,
	"chr2":	813178,	
	"chr3":	316617,	
	"chr4":	1531919,	
	"chr5":	576869,	
	"chr6":	270148,	
	"chr7":	1090947,	
	"chr8":	562643,	
	"chr9":	439885,	
	"chr10":	745741,	
	"chr11":	666454,	
	"chr12":	1078175,	
	"chr13":	924429,	
	"chr14":	784333,	
	"chr15":	1091289,	
	"chr16":	948062,	
	"chrm":	85779
}

# chromsome length dictionary: Kluyveromyces lactis
chr_len_klactis = {
	"chrA":	1062590,
	"chrB":	1320834,	
	"chrC":	1753957,	
	"chrD":	1715506,	
	"chrE":	2234072,	
	"chrF":	2602197	
}



#D. mel's chromsome len
chr_len_dmel_eu = {"2R": 21146708,"3R": 27905053, "4": 1351857,"3L": 24543557,"2L": 23011544,"X": 22422827}
chr_len_dmel_het = {"2RHet": 3288761,"3RHet": 2517507,"3LHet": 2555491,"YHet": 347038,"2LHet": 368872,"XHet": 204112}
chr_len_dmel = {"2R": 21146708,"3R": 27905053, "4": 1351857,"3L": 24543557,"2L": 23011544,"X": 22422827, "2RHet": 3288761,"3RHet": 2517507,"3LHet": 2555491,"YHet": 347038,"2LHet": 368872,"XHet": 204112}

## fly chromsome len; assembly dm3 from ucsc browser;
chr_len_dm3 = {
	"chr2L"    :23011544,
	"chr2LHet" :368872,
	"chr2R"    :21146708,
	"chr2RHet" :3288761,
	"chr3L"    :24543557,
	"chr3LHet" :2555491,
	"chr3R"    :27905053,
	"chr3RHet" :2517507,
	"chr4"     :1351857,
	"chrU"     :10049037,
	"chrUextra":29004656,
	"chrX"     :22422872,
	"chrXHet"  :204112,
	"chrYHet"  :347038,
	"chrM"     :19517
}

GC_dmel = {"2RHet": 39.2029805202, "3RHet": 39.4713435375, "3LHet": 39.4357010778, "YHet": 38.4467434907, "2LHet": 39.0781946608, "XHet": 40.0907344987, "2R": 43.2730724474, "3R": 42.9149444726, "dmel_mitochondrion_genome": 17.840856689, "U": 39.9763128435, "4": 35.523692498, "3L": 41.9463077267, "2L": 41.8353617242, "X": 42.4998120471}

# human chromsome len; assembly hg18 from ucsc browser; chromosome length were modified based on last tag coor from Zhao KJ's paper
# chr3, chr4, 19, 
chr_len_hg18 = {

	"chr1" :249250621,
	"chr10" :135534747,
	"chr11" :135006516,
	"chr12" :133851895,
	"chr13" :115169878,
	"chr14" :107349540,
	"chr15" :102531392,
	"chr16" :90354753,
	"chr17" :81195210,
	"chr18" :78077248,
	"chr19" :63806188,
	"chr2" :243199373,
	"chr20" :63025520,
	"chr21" :48129895,
	"chr22" :51304566,
	"chr3" :199431662,
	"chr4" :191263063,
	"chr5" :180915260,
	"chr6" :171115067,
	"chr7" :159138663,
	"chr8" :146364022,
	"chr9" :141213431,
	"chrX" :155270560,
	"chrY" :59373566

}


#
#chromsome length dictionary: mouse; got it from Celine Han
#

chr_len_mm9 = {
   "chr1" : 197195432,
   "chr2" : 181748087,
   "chr3" : 159599783,
   "chr4" : 155630120,
   "chr5" : 152537259,
   "chr6" : 149517037,
   "chr7" : 152524553,
   "chr8" : 131738871,
   "chr9" : 124076172,
   "chr10" : 129993255,
   "chr11" : 121843856,
   "chr12" : 121257530,
   "chr13" : 120284312,
   "chr14" : 125194864,
   "chr15" : 103494974,
   "chr16" : 98319150,
   "chr17" : 95272651,
   "chr18" : 90772031,
   "chr19" : 61342430,
   "chrX"  : 166650296,
   "chrY"  : 58682461
}


tri_nuc_bendability = {

	"AAT" : - 0.280, "ATT" : - 0.280,
	"AAA" : - 0.274, "TTT" : - 0.274, 
	"CCA" : - 0.246, "TGG" : - 0.246,
	"AAC" : - 0.205, "GTT" : - 0.205,
	"ACT" : - 0.183, "AGT" : - 0.183,
	"CCG" : - 0.136, "CGG" : - 0.136,
	"ATC" : - 0.110, "GAT" : - 0.110, 
	"AAG" : - 0.081, "CTT" : - 0.081,
	"CGC" : - 0.077, "GCG" : - 0.077,
	"AGG" : - 0.057, "CCT" : - 0.057, 
	"GAA" : - 0.037, "TTC" : - 0.037, 
	"ACG" : - 0.033, "CGT" : - 0.033, 
	"ACC" : - 0.032, "GGT" : - 0.032,
	"GAC" : - 0.013, "GTC" : - 0.013, 
	"CCC" : - 0.012, "GGG" : - 0.012, 
	"ACA" : - 0.006, "TGT" : - 0.006, 
	"CGA" : - 0.003, "TCG" : - 0.003, 
	
	"GGA" : 0.013, "TCC" : 0.013, 
	"CAA" : 0.015, "TTG" : 0.015, 
	"AGC" : 0.017, "GCT" : 0.017, 
	"GTA" : 0.025, "TAC" : 0.025, 
	"AGA" : 0.027, "TCT" : 0.027, 
	"CTC" : 0.031, "GAG" : 0.031, 
	"CAC" : 0.040, "GTG" : 0.040, 
	"TAA" : 0.068, "TTA" : 0.068, 
	"GCA" : 0.076, "TGC" : 0.076, 
	"CTA" : 0.090, "TAG" : 0.090, 
	"GCC" : 0.107, "GGC" : 0.107,
	"ATG" : 0.134, "CAT" : 0.134, 
	"CAG" : 0.175, "CTG" : 0.175, 
	"ATA" : 0.182, "TAT" : 0.182,
	"TCA" : 0.194, "TGA" : 0.194
}

pared_trinucs = [
	("AAT", "ATT"),
	("AAA", "TTT"), 
	("CCA", "TGG"),
	("AAC", "GTT"),
	("ACT", "AGT"),
	("CCG", "CGG"),
	("ATC", "GAT"), 
	("AAG", "CTT"),
	("CGC", "GCG"),
	("AGG", "CCT"), 
	("GAA", "TTC"), 
	("ACG", "CGT"), 
	("ACC", "GGT"),
	("GAC", "GTC"), 
	("CCC", "GGG"), 
	("ACA", "TGT"), 
	("CGA", "TCG"), 
	("GGA", "TCC"), 
	("CAA", "TTG"), 
	("AGC", "GCT"), 
	("GTA", "TAC"), 
	("AGA", "TCT"), 
	("CTC", "GAG"), 
	("CAC", "GTG"), 
	("TAA", "TTA"), 
	("GCA", "TGC"), 
	("CTA", "TAG"), 
	("GCC", "GGC"),
	("ATG", "CAT"), 
	("CAG", "CTG"), 
	("ATA", "TAT"),
	("TCA", "TGA")
]

set_dinuc = set([
	"AA", "AC", "AG", "AT", 
	"CA", "CC", "CG", "CT",
	"GA", "GC", "GG", "GT",
	"TA", "TC", "TG", "TT"
])


dict_dinuc_pair = {
	"AA|TT" : ("AA", "TT"),
	"AC|GT" : ("AC", "GT"),
	"AG|CT" : ("AG", "CT"),
	# AT self reverse conmplementary
	
	"CA|TG" : ("CA", "TG"),
	"CC|GG" : ("CC", "GG"),
	# CG self reverse complementary
	"CT|AG" : ("CT", "AG"),
	
	"GA|TC" : ("GA", "TC")
}

def testGenomeInfo():
	print "genome_info is working fine!!"

