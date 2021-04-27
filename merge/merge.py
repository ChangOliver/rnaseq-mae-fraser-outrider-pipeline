import os, glob
import pandas as pd
import numpy as np
import sqlite3
from optparse import OptionParser

# argument parsing 
def parse_arg():
	parser = OptionParser()
	parser.add_option("-m", "--mae", type="string", default=None,
					help="input MAE.result.csv", metavar="FILE")
	parser.add_option("-f", "--fraser", type="string", default=None, 
					help="input FRASER.result.csv", metavar="FILE")
	parser.add_option("-u", "--outrider", type="string", default=None, 
					help="input OUTRIDER.result.csv", metavar="FILE")
	parser.add_option("-d", "--dataset", type="int", default=None, 
					help="37 for GRCh37 and 38 for GRCh38", metavar="INT")
	parser.add_option("-o", "--output", type="string", default=None, 
					help="output directory to store results", metavar="FILE")

	(opt, args) = parser.parse_args()

	if opt.mae is None or opt.fraser is None or opt.outrider is None or opt.output is None:
		parser.error("Input & output must be provided.")

	if opt.dataset is None:
		parser.error("Choose a dataset to use.")

	if opt.dataset != 37 and opt.dataset != 38:
		parser.error("Choose 37 for GRCh37 and 38 for GRCh38.")

	if not os.path.isdir(opt.output):
		os.mkdir(opt.output)

	opt.output = opt.output if opt.output[-1] == '/' else opt.output + '/' 
	return opt

def read_input(opt):
	# input mart
	mart = pd.read_csv("mart37_export.csv" if opt.dataset == 37 else "mart38_export.csv", dtype={'chr': "string"})
	mart.rename(columns={"GeneName": "OUTRIDER_geneID", "Chromosome": "OUTRIDER_chr", "GeneStart": "OUTRIDER_start", "GeneEnd":"OUTRIDER_end"}, inplace=True)

	# input MAE
	mae = pd.read_csv(opt.mae, dtype={'contig': "string"})
	mae.drop('sampleID', inplace=True, axis=1)
	mae = mae.add_prefix("MAE_")
	mae.sort_values(by=["MAE_contig", "MAE_position"], inplace=True)
	
	# input FRASER
	fraser = pd.read_csv(opt.fraser, dtype={'seqnames': "string"})
	fraser.drop('sampleID', inplace=True, axis=1)
	fraser = fraser.add_prefix("FRASER_")
	fraser.sort_values(by=["FRASER_seqnames", "FRASER_start", "FRASER_end"], inplace=True)

	# input OUTRIDER
	outrider = pd.read_csv(opt.outrider)
	outrider.drop('sampleID', inplace=True, axis=1)
	outrider = outrider.add_prefix("OUTRIDER_")
	outrider = outrider.merge(mart, how="left", on="OUTRIDER_geneID")
	cols = outrider.columns.tolist()
	cols = cols[0:1] + cols[-3:] + cols[1:-3]
	outrider = outrider[cols]
	outrider.sort_values(by=["OUTRIDER_chr", "OUTRIDER_start", "OUTRIDER_end"], inplace=True)

	return {"mae": mae, "fraser": fraser, "outrider":outrider}

def merge(mae, fraser, outrider):

	# use sqlite3
	conn = sqlite3.connect(':memory:')
	mae.to_sql('mae', conn, index=False)
	fraser.to_sql('fraser', conn, index=False)
	outrider.to_sql('outrider', conn, index=False)

	# build query (full outer join)
	query = '''
		select  
			mae.*, fraser.*, outrider.*
		from
			mae left join fraser on MAE_position between FRASER_start and FRASER_end and MAE_contig = FRASER_seqnames
				left join outrider on MAE_position between OUTRIDER_start and OUTRIDER_end and MAE_contig = OUTRIDER_chr
		union
		select
			mae.*, fraser.*, outrider.*
		from
			fraser left join mae on MAE_position between FRASER_start and FRASER_end and FRASER_seqnames = MAE_contig
					left join outrider on FRASER_start >= OUTRIDER_start and FRASER_end <= OUTRIDER_end and FRASER_seqnames = OUTRIDER_chr
		union
		select
			mae.*, fraser.*, outrider.*
		from
			outrider left join mae on MAE_position between OUTRIDER_start and OUTRIDER_end and MAE_contig = OUTRIDER_chr
					left join fraser on FRASER_start >= OUTRIDER_start and FRASER_end <= OUTRIDER_end and FRASER_seqnames = OUTRIDER_chr
		'''

	res = pd.read_sql_query(query, conn)

	# typecast
	res = res.astype({ 	'MAE_position': 'Int64', 'MAE_refCount': 'Int64', 'MAE_altCount': 'Int64',
						'MAE_totalCount': 'Int64', 'MAE_signif': 'boolean', 'MAE_signif_ALT': 'boolean',
						'MAE_rare': 'boolean', 'FRASER_start': 'Int64', 'FRASER_end': 'Int64',
						'FRASER_width': 'Int64', 'FRASER_counts': 'Int64', 'FRASER_totalCounts': 'Int64',
						'OUTRIDER_start': 'Int64', 'OUTRIDER_end': 'Int64', 'OUTRIDER_rawcounts': 'Int64',
						'OUTRIDER_aberrant': 'boolean', 'OUTRIDER_AberrantBySample': 'Int64', 'OUTRIDER_AberrantByGene': 'Int64'
					})

	res[["FRASER_hgncSymbol", "OUTRIDER_geneID", "MAE_contig", "FRASER_seqnames", "OUTRIDER_chr"]] = res[["FRASER_hgncSymbol", "OUTRIDER_geneID", "MAE_contig", "FRASER_seqnames", "OUTRIDER_chr"]].fillna(".")

	res["sampleID"] = sample
	res["seqnames"] = np.select([ res.MAE_contig != ".", res.FRASER_seqnames != ".", res.OUTRIDER_chr != "."], [res.MAE_contig, res.FRASER_seqnames, res.OUTRIDER_chr], default=None)
	res["geneID"] = np.select([ res.FRASER_hgncSymbol != ".", res.OUTRIDER_geneID != "."], [res.FRASER_hgncSymbol, res.OUTRIDER_geneID], default=None)
	
	res.drop(['MAE_contig', 'FRASER_seqnames', 'OUTRIDER_chr', 'OUTRIDER_geneID'], inplace=True, axis=1)
	
	# sort & arrange
	res.sort_values(by=["seqnames", "geneID", "MAE_position", "FRASER_start", "OUTRIDER_start"], inplace=True, na_position="last")
	cols = res.columns.tolist()
	cols = cols[-3:] + cols[:-3]
	res = res[cols]

	return res

if __name__ == "__main__":
	# parse arguments
	opt = parse_arg()

	# extract sampleID
	sample = os.path.basename(opt.fraser)[:-18]

	print("Merging " + sample)

	data = read_input(opt)
	res = merge(data["mae"], data["fraser"], data["outrider"])
	res.to_csv(opt.output + sample + ".rnaseq.merge.tsv", index=False, na_rep=".")
