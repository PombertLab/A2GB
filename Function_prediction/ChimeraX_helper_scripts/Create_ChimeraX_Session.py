#!/usr/bin/python
from chimerax.core.commands import run
from sys import argv
import argparse
import re
import os

name = 'Create_ChimeraX_Session.py'
version = '0.2'
updated = '2021-04-26'

usage = f'''
NAME		{name}
VERSION		{version}
UPDATED		{updated}

SYNOPSIS	This script is used to align a reference .pdb to a predicted .pdb,
		changes the predicted .pdb color, hides all atoms, shows only matched
		chains, and saves the result as a ChimeraX session (.cxs)

COMMAND		{name} \\
			-p ...preference \\
			-r ...reference

OPTIONS

-p_f (--pred)		Predicted pdb file
-p_m (--pfam)		PFAM match file
-r_m (--rcsb)		RCSD match file
-o_d (--outdir)		Output directory for .cxs files [Default: ./CXS]
'''

if len(argv) < 2:
	print(usage)
	exit()

outdir = './CXS'

parser = argparse.ArgumentParser(usage=usage)
parser.add_argument('-p_f','--pred',type=str,required=True)
parser.add_argument('-p_m','--pfam',type=str)
parser.add_argument('-r_m','--rcsb',type=str)
parser.add_argument('-o_d','--outdir',type=str)
parser.add_argument('--nogui')

args = parser.parse_args()
pred = args.pred
pfam = args.pfam
rcsb = args.rcsb
if(args.outdir):
	outdir = args.outdir

locus_tag = os.path.basename(pred)
locus_tag = re.findall("\w+",locus_tag)

## Load pdb files
model_pfam = False
model_pfam_name = False
model_rcsb = False
model_rcsb_name = False

model_pred = run(session,f"open {pred}")[0][0]
model_pred_name = (model_pred.id_string)

if(pfam):
	model_pfam = run(session,f"open {pfam}")[0][0]
	model_pfam_name = (model_pfam.id_string)
if(rcsb):
	model_rcsb = run(session,f"open {rcsb}")[0][0]
	model_rcsb_name = (model_rcsb.id_string)


## Prepare file for display by hiding everything
run(session,"hide atoms")
run(session,"hide ribbons")

## Superimpose match structure against prediction
chain_pred_pfam = False
chain_pred_rcsb = False
chain_pfam = False
chain_rcsb = False

if(model_pfam and model_rcsb):
	match = run(session,f"match #{model_pred_name} to #{model_pfam_name}")
	chain_pred_pfam = (match[0][0].unique_chain_ids)[0]
	chain_pfam = (match[0][1].unique_chain_ids)[0]
	match = run(session,f"match #{model_pred_name} to #{model_rcsb_name}")
	chain_pred_rcsb = (match[0][0].unique_chain_ids)[0]
	chain_rcsb = (match[0][1].unique_chain_ids)[0]
elif(model_pred and model_pfam):
	match = run(session,f"match #{model_pred_name} to #{model_pfam_name}")
	chain_pred_pfam = (match[0][0].unique_chain_ids)[0]
	chain_pfam = (match[0][1].unique_chain_ids)[0]
elif(model_pred and model_rcsb):
	match = run(session,f"match #{model_pred_name} to #{model_rcsb_name}")
	chain_pred_rcsb = (match[0][0].unique_chain_ids)[0]
	chain_rcsb = (match[0][1].unique_chain_ids)[0]

## Color reference structure a diferrent color
if(chain_pfam):
	run(session,f"color #{model_pfam_name}/{chain_pfam} #FF00FF ribbons")
if(chain_rcsb):
	run(session,f"color #{model_rcsb_name}/{chain_rcsb} #00FFFF ribbons")

## Show only matching chains
run(session,f"show #{model_pred_name} ribbons")
if(chain_pfam):
	run(session,f"show #{model_pfam_name}/{chain_pfam} ribbons")
if(chain_rcsb):
	run(session,f"show #{model_rcsb_name}/{chain_rcsb} ribbons")

## Orient the chain to view
run(session,"view")

## Save match as a new file
run(session,f"save {outdir}/{locus_tag[0]}.cxs format session")

quit()