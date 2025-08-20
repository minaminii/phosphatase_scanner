import os, sys, tempfile, string, random
import re, getopt
import csv
from subprocess import call
import subprocess
from config import *

# load biopython
#import sys
#sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/')
from Bio.Blast import NCBIXML
from Bio import SeqIO

class PhosphataseScanner:
    def __init__(self, fasta_file):
        self.prefix = fasta_file
        self.fasta_file = fasta_file
    def wrapper(self, job='', fasta=''):
        self.hmmscan_general = Hmmscan(self.fasta_file, HMM_GENERAL)
        self.hmmscan_general.wrapper()
        self.hmmscan_specific = Hmmscan(self.fasta_file, HMM_PD)
        self.hmmscan_specific.wrapper()
        self.summary()
        self.write_summary()
        #self.hmmscan_protein = Hmmscan(self.fasta_file, HMM_PROTEIN)
        #self.hmmscan_protein.wrapper()
    def summary(self):
        q2h_general = self.hmmscan_general.tab_query2hits
        q2h_specific = self.hmmscan_specific.tab_query2hits
        query = set(q2h_general.keys())
        query.update(set(q2h_specific.keys()))
        rows = []
        for q in query:
            row = {'Query': q}
            row.update(self._best_hit(q2h_general.get(q, []), 'General HMM'))
            row.update(self._best_hit(q2h_specific.get(q, []), 'Specific HMM'))
            rows.append(row)
        self.summary_rows = rows
        self.summary_fieldnames = ['Query', 'General HMM Name', 'General HMM E-value', 'Specific HMM Name', 'Specific HMM E-value']
    def write_summary(self):
        self.summary_file = self.prefix + '.summary.tab'
        sys.stderr.write('#LOG: run write_summary. ' + self.summary_file + '\n')
        with open(self.summary_file,'wb') as fou:
            dw = csv.DictWriter(fou, delimiter='\t', fieldnames=self.summary_fieldnames)
            headers = {}
            for n in dw.fieldnames:
                headers[n] = n
            dw.writerow(headers)
            for row in self.summary_rows:
                dw.writerow(row)
    def _best_hit(self, hits, name):
        name_key = name + ' Name'
        evalue_key = name + ' E-value'
        bh = {name_key: '', evalue_key: ''}
        if len(hits) > 0:
            h0 = hits[0]
            if h0.has_key('target_name'): bh[name_key] = h0.get('target_name')
            if h0.has_key('evalue_1st_domain'): bh[evalue_key] = h0.get('evalue_1st_domain')
        return bh

class Hmmscan:
    def __init__(self, fasta_file, hmmdb0):
        self.hmmscan_bin = HMMSCAN_BIN
        self.fasta_file = fasta_file
        self.hmmdb = hmmdb0['file']
        self.hmmdb_id = hmmdb0['id']
    def wrapper(self):
        self.run()
        self.list_hit_sequence_ids()
        self.get_hit_sequence_fasta()
        self.parse()
    def run(self):
        self.out_file = self.fasta_file + '.' + self.hmmdb_id + '.hmmscan.out'
        self.tab_file = self.fasta_file + '.' + self.hmmdb_id + '.hmmscan.tab'
        self.dom_tab_file = self.fasta_file + '.' + self.hmmdb_id + '.hmmscan.dom.tab'
        if os.path.exists(self.out_file):
            sys.stderr.write("#LOG: hmm result file exists. skip hmmscan database " + self.hmmdb_id + ".\n")
            return 0
        if not os.path.exists(self.hmmdb):
            sys.stderr.write('#ERR: cannot find hmmdb ' + self.hmmdb + '\n')
        # Example:
        # hmmscan -o all_protein_sequneces_vs_PD.out --tblout all_protein_sequneces_vs_PD.tab --domtblout all_protein_sequneces_vs_PD.dom.tab /Users/bioinfo/workspace/phosphatase/web/kp1/colt/static/colt/data/HMM/specific/PD.hmm3 all_protein_sequneces.fasta &> all_protein_sequneces_vs_PD.log &
        cline = ' '.join((self.hmmscan_bin, '-o', self.out_file, '--tblout', self.tab_file, '--domtblout', self.dom_tab_file, self.hmmdb, self.fasta_file))
        sys.stderr.write('#LOG: run hmmscan. ' + cline + '\n')
        p = subprocess.Popen(cline, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        out, err = p.communicate()
    def list_hit_sequence_ids(self):
        if hasattr(self, 'hit_sequence_ids'): return self.hit_sequence_ids
        ids = set()
        with open(self.tab_file, 'r') as f:
            for row in f:
                if re.search(r'^#', row): continue
                fields = row.split()
                ids.add(fields[2])
        self.hit_sequence_ids = list(ids)
        return list(ids)
    def get_hit_sequence_fasta(self):
        ids = self.list_hit_sequence_ids()
        records = []
        for seq_record in SeqIO.parse(self.fasta_file, "fasta"):
            if seq_record.id in ids:
                records.append(seq_record)
        self.hit_sequence_fasta_file = self.fasta_file + '.' + self.hmmdb_id + '.hit_sequence.fasta'
        SeqIO.write(records, self.hit_sequence_fasta_file, "fasta")
    def parse(self):
        self.tab_to_dict()
    def tab_to_string(self):
        sys.stderr.write('#LOG: read hmmscan table out ' + self.tab_file + '\n')
        with open(self.tab_file, 'r') as f:
            self.tab_string = f.read()
    def dom_tab_to_string(self):
        sys.stderr.write('#LOG: read hmmscan domain table out ' + self.dom_tab_file + '\n')
        with open(self.dom_tab_file, 'r') as f:
            self.dom_tab_string = f.read()
    def tab_to_dict(self):
        query2hits = {}
        with open(self.tab_file, 'r') as f:
            for row in f:
                if re.search(r'^#', row): continue
                d = self._tab_row_to_dict(row)
                if query2hits.has_key(d['query_name']):
                    query2hits[d['query_name']].append(d)
                else:
                    query2hits[d['query_name']] = [d]
        self.tab_query2hits = query2hits
        return query2hits
    def _tab_row_to_dict(self, row):
        f = row.split()
        target_name = f[0]
        query_name = f[2]
        evalue_fullseq = f[4]
        score_fullseq = f[5]
        evalue_1st_domain = f[7]
        score_1st_domain = f[8]
        p = re.compile('^(\w+)\.(\S+)\.(\w+)$')
        m = p.match(target_name)
        d = {'target_name': target_name,
             'query_name': query_name,
             'evalue_fullseq': evalue_fullseq,
             'score_fullseq': score_fullseq,
             'evalue_1st_domain': evalue_1st_domain,
             'score_1st_domain': score_1st_domain,}
        return d

# parameters
def get_argument(argv, usage="usage"):
    try:
        opts, args = getopt.getopt(argv, "hi:", ['infile='])
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    if args or not opts:
        print usage
        sys.exit(3)
    for opt, arg in opts:
        if opt == '-h':
            print usage
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
            return infile
        else:
            print usage


def main(argv):
    usage = "python " + sys.argv[0] + " -i infile"
    infile = get_argument(argv, usage=usage)
    ps = PhosphataseScanner(infile)
    ps.wrapper()

if __name__ == '__main__':
    main(sys.argv[1:])
