KLASS_LEVELS = ['Group', 'Family', 'Subfamily']

HMM_GENERAL = {'id': 'general',
               'file': 'data/general.hmm3',}
HMM_PD = {'id': 'PD', 
          'file': 'data/PD.hmm3',}
HMM_PROTEIN = {'id': 'protein',
               'file': 'data/protein.hmm3',}
HMMSCAN_BIN = '/usr/local/bin/hmmscan'

BLASTALL = '/opt/local/bin/blastall'
BLASTDB_PD = {'id': 'pd',
              'name': 'Phosphatase Domain',
              'file': 'colt/static/colt/data/blastdb/PD.fasta',}
BLASTDB_PROTEIN = {'id': 'protein',
                   'name': 'Protein',
                   'file': 'colt/static/colt/data/blastdb/protein.fasta',}
