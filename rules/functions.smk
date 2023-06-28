from snakemake.utils import min_version

import sys
import os
import platform
import re
import csv
import gzip
import warnings

def build_dict_from_key_value_pairs(list_of_dicts):
	'''
	Each dict in the list MUST contain the keys 'id' and 'path'.
	path will be converted to absolute paths.
	Takes the list of dictionaries and makes it into a new dictionary where the keys are the id values from each dictionary and the values are each dictionary
	e.g. [{"id":"a", "value": 1}, {"id":"b","value":2}] ->
	{"a":{"id":"a", "value": 1}, "b":{"id":"b","value":2}}
	'''
	out_dict = {}
	for d in list_of_dicts:
		d['path'] = os.path.abspath(d['path'])
		out_dict[d['id']] = d
	return(out_dict)

def build_gwas_type_dict_str(list_of_dicts):
	gwas_dict_list = []
	for d in list_of_dicts:
		gwas_dict_list.append("\""+d['id']+"\""+":"+"\""+d['type']+"\"")
	gwas_dict_str = "'{"+",".join(gwas_dict_list)+"}'"
	return(gwas_dict_str)



##############################
GWAS_SUMSTATS = build_dict_from_key_value_pairs(config['GWAS_SUMSTATS'])
GWAS_TYPE_DICT_STR = build_gwas_type_dict_str(config['GWAS_SUMSTATS'])
eQTL_SUMSTATS = build_dict_from_key_value_pairs(config['eQTL_INPUT'])
