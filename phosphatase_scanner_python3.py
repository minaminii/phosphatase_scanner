import os
import sys
import logging
import getopt
import csv
import subprocess

from Bio import SeqIO

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.NOTSET)

# Set default directory to script's directory
SCANNER_DIR = os.path.realpath(os.path.dirname(__file__))

KLASS_LEVELS = ['Group', 'Family', 'Subfamily']
HMM_GENERAL = {
    'id': 'general',
    'file': f'{SCANNER_DIR}/data/general.hmm3'
}
HMM_PD = {
    'id': 'PD',
    'file': f'{SCANNER_DIR}/data/PD.hmm3'
}
HMM_PROTEIN = {
    'id': 'protein',
    'file': f'{SCANNER_DIR}/data/protein.hmm3'
}
HMMSCAN_BIN = 'hmmscan'

class PhosphataseScanner:
    def __init__(self, fasta_file):
        self.prefix = fasta_file
        self.fasta_file = fasta_file

        self.hmmscan_general = Hmmscan(self.fasta_file, HMM_GENERAL)
        self.hmmscan_specific = Hmmscan(self.hmmscan_general.hit_sequence_fasta_file, HMM_PD)

        self.summary_rows = None
        self.summary_file = self.prefix + '.summary.tab'

    def wrapper(self):
        self.hmmscan_general.wrapper()
        self.hmmscan_specific.wrapper()
        self.summary()
        self._write_summary()
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

    def _write_summary(self):
        logger.log(logging.INFO, 'Run write_summary %s', self.summary_file)
        summary_fieldnames = [
            'Query',
            'General HMM Name',
            'General HMM E-value',
            'Specific HMM Name',
            'Specific HMM E-value'
        ]
        with open(self.summary_file,'w', encoding='utf-8') as fou:
            dw = csv.DictWriter(fou, delimiter='\t', fieldnames=summary_fieldnames)
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
            bh[name_key] = h0.get('target_name', None)
            bh[evalue_key] = h0.get('evalue_1st_domain', None)
        return bh

class Hmmscan:
    def __init__(self, fasta_file, hmmdb0):
        self.hmmscan_bin = HMMSCAN_BIN
        self.fasta_file = fasta_file
        self.hmmdb = hmmdb0['file']
        self.hmmdb_id = hmmdb0['id']

        self.hit_sequence_ids = None
        self.hit_sequence_fasta_file = self.fasta_file + '.' + self.hmmdb_id + '.hit_sequence.fasta'

        self._out_file = self.fasta_file + '.' + self.hmmdb_id + '.hmmscan.out'
        self._tab_file = self.fasta_file + '.' + self.hmmdb_id + '.hmmscan.tab'
        self._dom_tab_file = self.fasta_file + '.' + self.hmmdb_id + '.hmmscan.dom.tab'

        self.tab_query2hits = None

    def wrapper(self):
        self.run()
        self.list_hit_sequence_ids()
        self.get_hit_sequence_fasta()
        self._tab_to_dict()

    def run(self):
        if os.path.exists(self._out_file):
            logging.info('hmm result file exists. skip hmmscan database %s', self.hmmdb_id)
            return 0
        if not os.path.exists(self.hmmdb):
            logging.error('cannot find hmmdb %s', self.hmmdb)

        # Example:
        # hmmscan \
        # -o all_protein_sequneces_vs_PD.out \
        # --tblout all_protein_sequneces_vs_PD.tab \
        # --domtblout all_protein_sequneces_vs_PD.dom.tab \
        # /Users/bioinfo/workspace/phosphatase/web/kp1/colt/static/colt/data/HMM/specific/PD.hmm3 \
        # all_protein_sequneces.fasta &> all_protein_sequneces_vs_PD.log &

        hmm_command = ' '.join((
            self.hmmscan_bin,
            '-o', self._out_file,
            '--tblout', self._tab_file,
            '--domtblout', self._dom_tab_file,
            self.hmmdb, self.fasta_file
        ))
        logger.info('Run hmmscan: %s', hmm_command)
        with subprocess.Popen(
            hmm_command,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            shell=True
        ) as proc:
            out, err = proc.communicate()
            if out != b'':
                logger.debug(out.decode())
            if err != b'':
                logger.error(err.decode())

    def list_hit_sequence_ids(self):
        '''
        Read the tab output from `.hmmscan.tab`
        '''
        ids = set()
        # Open xxx.hmmscan.tab
        with open(self._tab_file, 'r', encoding='utf-8') as f:
            # iterate each row
            for row in f:
                # Ignore rows start with #
                if row.startswith('#'):
                    continue
                # Split by idk what but I think it tried to obtain the 3rd column
                fields = row.split()
                # Get 3rd column from the record - which is `query name``
                ids.add(fields[2])
        self.hit_sequence_ids = list(ids)
        # Return a list of ids, in this case is likely to be a query name that contains a target
        return list(ids)

    def get_hit_sequence_fasta(self):
        '''
        Extract sequences with target hit from the input fasta.
        '''
        ids = self.list_hit_sequence_ids()
        records = []
        for seq_record in SeqIO.parse(self.fasta_file, "fasta"):
            if seq_record.id in ids:
                records.append(seq_record)
        SeqIO.write(records, self.hit_sequence_fasta_file, "fasta")

    def _tab_to_dict(self):
        query2hits = {}
        with open(self._tab_file, 'r', encoding='utf-8') as f:
            for row in f:
                if row.startswith('#'):
                    continue
                d = self._tab_row_to_dict(row)
                if query2hits.get(d['query_name'], None):
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
        return {
            'target_name': target_name,
            'query_name': query_name,
            'evalue_fullseq': evalue_fullseq,
            'score_fullseq': score_fullseq,
            'evalue_1st_domain': evalue_1st_domain,
            'score_1st_domain': score_1st_domain
        }

def get_argument(argv, usage="usage"):
    '''
    Process arguments
    '''
    usage = "python " + sys.argv[0] + " -i infile"
    try:
        opts, args = getopt.getopt(argv, "hi:", ['infile='])
    except getopt.GetoptError:
        print(usage)
        sys.exit(1)
    if args or not opts:
        print(usage)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit(0)
        elif opt in ("-i", "--infile"):
            infile = arg
            return infile
        else:
            print(usage)


def main(argv):
    '''
    Main entrypoint
    '''
    infile = get_argument(argv)
    ps = PhosphataseScanner(infile)
    ps.wrapper()

if __name__ == '__main__':
    main(sys.argv[1:])
