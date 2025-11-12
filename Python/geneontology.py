#!/bin/env python

import re
import gm

def load_OBO(filename='go-basic.obo'):
	"""
	parse the OBO file and returns the graph
	obsolete terms are discarded
	only is_a and part_of relationships are loaded

	Extract of a file to be parsed:
	[Term]
	id: GO:0000028
	name: ribosomal small subunit assembly
	namespace: biological_process
	def: "The aggregation, arrangement and bonding together of constituent RNAs and proteins to form the small ribosomal subunit." [GOC:jl]
	subset: gosubset_prok
	synonym: "30S ribosomal subunit assembly" NARROW [GOC:mah]
	synonym: "40S ribosomal subunit assembly" NARROW [GOC:mah]
	is_a: GO:0022618 ! ribonucleoprotein complex assembly
	relationship: part_of GO:0042255 ! ribosome assembly
	relationship: part_of GO:0042274 ! ribosomal small subunit biogenesis
	"""
	def parseTerm(lines):
		# search for obsolete
		for l in lines:
			if l.startswith('is_obsolete: true'):
				return
		# otherwise create node
		go_id = re_go_id.match(lines.pop(0)).group(1)
		go_attr = add_node(go_graph, go_id) # add node to graph and get the node attribute dict
		go_attr['type'] = 'GOTerm'
		for line in lines:
			if re_go_name.match(line): go_attr['name'] = re_go_name.match(line).group(1)
			elif re_go_namespace.match(line): go_attr['namespace'] = re_go_namespace.match(line).group(1)
			elif re_go_def.match(line): go_attr['def'] = re_go_def.match(line).group(1)
			elif re_go_alt_id.match(line): go_graph['alt_id'][ re_go_alt_id.match(line).group(1) ] = go_id  # alt_id → go_id
			elif re_go_is_a.match(line):
				parent_id = re_go_is_a.match(line).group(1)
				add_edge(go_graph, go_id, parent_id, { 'relationship': 'is a' })
			elif re_go_part_of.match(line):
				parent_id = re_go_part_of.match(line).group(1)
				add_edge(go_graph, go_id, parent_id, { 'relationship': 'part of' })
	# method main
	go_graph          = create_graph(directed=True, weighted=False)
	go_graph['alt_id'] = {} # alternate GO ids
	# regexp to parse term lines
	re_go_id          = re.compile(r'^id:\s+(GO:\d+)\s*$')
	re_go_name        = re.compile(r'^name:\s+(.+)\s*$')
	re_go_namespace   = re.compile(r'^namespace:\s+(.+)\s*$')
	re_go_def         = re.compile(r'^def:\s+(.+)\s*$')
	re_go_alt_id      = re.compile(r'^alt_id:\s+(GO:\d+)\s*$')
	re_go_is_a        = re.compile(r'^is_a:\s+(GO:\d+)\s')
	re_go_xref        = re.compile(r'^xref:\s+(\S+)\s*$')
	re_go_part_of      = re.compile(r'^relationship:\s+part_of\s+(GO:\d+)\s')
	# buffer each term lines, then parse lines to create GOTerm node
	with open(filename) as f:
		line = f.readline().rstrip()
		# skip header until first [Term] is reached
		while not line.startswith('[Term]'):
			line = f.readline().rstrip()
		buff = []
		line = f.readline()
		stop = False
		while line and not stop:
			line = line.rstrip()
			# new Term
			if line.startswith('[Term]'):
				parseTerm(buff)
				buff=[]
			# last Term
			elif line.startswith('[Typedef]'):
				parseTerm(buff)
				stop=True
			# or append to buffer
			else:
				buff.append(line)
			line = f.readline()
	return go_graph

def load_GOA(go, filename, warnings=True):
	"""
	parse GOA file and add annotated gene products to previsouly loaded graph go

	Extract of a file to be parsed:
	gaf-version: 2.1
	!GO-version: http://purl.obolibrary.org/obo/go/releases/2020-11-28/extensions/go-plus.owl
	UniProtKB       O05154  tagX            GO:0008360      GO_REF:0000043  IEA     UniProtKB-KW:KW-0133    P       Putative glycosyltransferase TagX       tagX|SAOUHSC_00644      protein 93061   20201128        UniProt

	UniProtKB       O05154  tagX            GO:0016740      GO_REF:0000043  IEA     UniProtKB-KW:KW-0808    F       Putative glycosyltransferase TagX       tagX|SAOUHSC_00644      protein 93061   20201128        UniProt

	UniProtKB       O05204  ahpF            GO:0000302      GO_REF:0000002  IEA     InterPro:IPR012081      P       Alkyl hydroperoxide reductase subunit F ahpF|SAOUHSC_00364      protein 93061   20201128        InterPro

		0        1       2   3       4             5          6        7      8             9                              10
				id    name        go_id               evidence-codes                     desc                           aliases

	GAF spec: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
	Column 	Content 						Required? 	Cardinality 	Example
	1 		DB 								required 	1 				UniProtKB
	2 		DB Object ID 					required 	1 				P12345
	3 		DB Object Symbol 				required 	1 				PHO3
	4 		Qualifier 						optional 	0 or greater 	NOT
	5 		GO ID 							required 	1 				GO:0003993
	6 		DB:Reference (|DB:Reference) 	required 	1 or greater 	PMID:2676709
	7 		Evidence Code 					required 	1 				IMP
	8 		With (or) From 					optional 	0 or greater 	GO:0000346
	9 		Aspect 							required 	1 				F
	10 		DB Object Name 					optional 	0 or 1 			Toll-like receptor 4
	11 		DB Object Synonym (|Synonym) 	optional 	0 or greater 	hToll 	Tollbooth
	12 		DB Object Type 					required 	1 				protein
	13 		Taxon(|taxon) 					required 	1 or 2 			taxon:9606
	14 		Date 							required 	1 				20090118
	15 		Assigned By 					required 	1 				SGD
	16 		Annotation Extension 			optional 	0 or greater 	part_of(CL:0000576)
	17 		Gene Product Form ID 			optional 	0 or 1 			UniProtKB:P12345-2
	"""
	with open(filename) as f:
		line = f.readline()
		while line:
			if not line.startswith('!'): # skip comments
				cols = line.rstrip().split('\t')
				gp_id = cols[1]
				gt_id = cols[4]
				if gt_id not in go['nodes']: # GOTerm not found search alternate ids
					while gt_id not in go['nodes'] and gt_id in go['alt_id']:
						gt_id = go['alt_id'][gt_id] # replace term by alternate
				if gt_id not in go['nodes']: # failure: warn user
					if warnings:
						print(f'Warning: could not attach a gene product ({gp_id}) to a non existing GO Term ({gt_id})')
				else: # success: GOTerm to attach to was found
					# create node for gene product if not already present
					if gp_id not in go['nodes']:
						gp_attr = add_node(go, gp_id, { 'id': gp_id, 'type': 'GeneProduct'})
					# create or update gene product attributes
					gp_attr = go['nodes'][gp_id]
					gp_attr['name'] = cols[2]
					gp_attr['desc'] = cols[9]
					gp_attr['aliases'] = cols[10].split('|')
					# attach gene product to GOTerm
					gt_attr = go['nodes'][gt_id]
					e_attr = add_edge(go, gp_id, gt_id)
					e_attr['relationship'] = 'annotation'
					if 'evidence-codes' not in e_attr:
						e_attr['evidence-codes'] = []
					e_attr['evidence-codes'].append( cols[6] )
			line = f.readline()

def GOTerms(go, gp_id, recursive=False):
	if not recursive:
		if gp_id in go['nodes']:
			return neighbors(go, gp_id)
	else:
		# TO DO: GOTerms directly linked (successors) AND their descendants should be returned
		return None
	return None

def GeneProducts(go, go_id, recursive=False):
	# TO DO: a list of GeneProducts directly linked (predecessors of type GeneProduct) should be returned if all is False
	# if all is True, a list GeneProducts connected by a path (ancestors of type GeneProduct) should be returned
	return None

def max_depth(go):
	# TO DO
	return None

##### main → tests #####
if __name__ == "__main__":
    print("# Gene Ontology module tests")
