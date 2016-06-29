from __future__ import print_function
import argparse
import csv
import sys
from collections import defaultdict

def chrom_key(chrom):
  chrom = chrom.lower()
  if chrom == 'x':
    chrom = 23
  elif chrom == 'y':
    chrom = 24
  elif chrom.isdigit():
    chrom = int(chrom)
  else:
    chrom = 999
  return chrom

def chrom_copy_number(chrom):
  chrom = chrom.lower()
  if chrom == 'x':
    cn = 1
  elif chrom == 'y':
    cn = 2
  elif chrom.isdigit():
    cn = 2
  else:
    cn = 2
  return cn
  
class CopyNumberWriter(object):
  def __init__(self, cn_output_fn):
    self._cn_output_fn = cn_output_fn
    # cellular_prevalence represents fraction of *all* cells that are affected
    # by CNAs, *not* just tumor cells.
    self._keys = ('chromosome', 'start', 'end', 'copy_number', 'minor_cn', 'major_cn', 'cellular_prevalence')

  def _write_header(self):
    self._cn_output.write('' + '\t'.join(self._keys) + '\n')

  def _write_cn_record(self, region):
    vals = [str(region[k]) for k in self._keys]
    self._cn_output.write('\t'.join(vals) + '\n')

  def write_cn(self, cn_regions):
    self._cn_output = open(self._cn_output_fn, 'w')
    self._write_header()

    chroms = sorted(cn_regions.keys(), key = chrom_key)
    for chrom in chroms:
      chrom_regions = cn_regions[chrom]
      chrom_regions.sort(key = lambda r: r['start'])
      for region in chrom_regions:
        # Insert chromosome into record, as including it originally would have
        # duplicated the dictionary key corresponding to per-chromosome CNAs.
        region['chromosome'] = chrom_key(chrom)
        region['copy_number'] = region['major_cn'] + region['minor_cn']
        self._write_cn_record(region)

    self._cn_output.close()

class MeanTotalCopyNumberWriter(object):
  def __init__(self, cn_output_fn):
    self._cn_output_fn = cn_output_fn
    # cellular_prevalence represents fraction of *all* cells that are affected
    # by CNAs, *not* just tumor cells.
    self._keys = ('chromosome', 'start', 'n_loci', 'end', 'mean_tcn')

  def _write_header(self):
    self._cn_output.write('#' + '\t'.join(self._keys) + '\n')

  def _write_cn_record(self, region):
    vals = [str(region[k]) for k in self._keys]
    self._cn_output.write('\t'.join(vals) + '\n')
    
  def write_mean_tcn(self, cn_regions):
    self._cn_output = open(self._cn_output_fn, 'w')
    self._write_header()

    chroms = sorted(cn_regions.keys(), key = chrom_key)
    for chrom in chroms:
      chrom_regions = cn_regions[chrom]
      chrom_regions.sort(key = lambda r: r['start'])
      for region in chrom_regions:
        # Insert chromosome into record, as including it originally would have
        # duplicated the dictionary key corresponding to per-chromosome CNAs.
        region['chromosome'] = chrom_key(chrom)
        self._write_cn_record(region)

    self._cn_output.close()
    
class AvailableCopyNumberWriter(object):
  def __init__(self, cn_output_fn):
    self._cn_output_fn = cn_output_fn
    # cellular_prevalence represents fraction of *all* cells that are affected
    # by CNAs, *not* just tumor cells.
    self._keys = ('chromosome', 'start', 'n_loci', 'end', 'avail_cn baf')

  def _write_header(self):
    self._cn_output.write('#' + '\t'.join(self._keys) + '\n')

  def _write_cn_record(self, region):
    vals = [str(region[k]) for k in self._keys]
    self._cn_output.write('\t'.join(vals) + '\n')
    
  def write_avail_cn(self, cn_regions):
    self._cn_output = open(self._cn_output_fn, 'w')
    self._write_header()

    chroms = sorted(cn_regions.keys(), key = chrom_key)
    for chrom in chroms:
      chrom_regions = cn_regions[chrom]
      chrom_regions.sort(key = lambda r: r['start'])
      for region in chrom_regions:
        # Insert chromosome into record, as including it originally would have
        # duplicated the dictionary key corresponding to per-chromosome CNAs.
        region['chromosome'] = chrom_key(chrom)            
        region['avail_cn baf'] = '\t'.join(map(str, region['avail_cn']))
        
        self._write_cn_record(region)

    self._cn_output.close()

class CnaParser(object):
  def parse(self):
    raise Exception('Not implemented')

class BattenbergParser(CnaParser):
  def __init__(self, bb_filename, cellularity):
    self._bb_filename = bb_filename
    self._cellularity = cellularity
    # Used by SMC-Het parser, which has fields shifted by 1.
    self._field_offset = 0
  
  def _compute_total_cn(self, cna1, cna2):
    cn1 = (cna1['major_cn'] + cna1['minor_cn']) * cna1['cellular_prevalence']
    if cna2:
      cn2 = (cna2['major_cn'] + cna2['minor_cn']) * cna2['cellular_prevalence']
    else:
      cn2 = 0
    total_cn = cn1 + cn2
    # print(cn1,cn2)
    return total_cn

  def _compute_avail_cn(self, cna1, cna2):
    max_cn=8
    avail_cn=[]
    for i in range(max_cn+1):
      if i <= cna1['major_cn']:
        avail_cn.append(1.0)
    else:
        avail_cn.append(0.0)
    return avail_cn

  def parse(self):
    cn_regions = defaultdict(list)
    pval_threshold = 0.05

    with open(self._bb_filename) as bbf:
      header = bbf.next()
      for line in bbf:
        fields = line.strip().split()
        chrom = fields[1 + self._field_offset].lower()
        start = int(fields[2 + self._field_offset])
        end = int(fields[3 + self._field_offset])
        pval = float(fields[5 + self._field_offset])
        logr = float(fields[6 + self._field_offset])

        cna1 = {}
        cna1['start'] = start
        cna1['end'] = end
        cna1['n_loci'] = 1
        cna1['logr'] = logr
        cna1['major_cn'] = int(fields[8 + self._field_offset])
        cna1['minor_cn'] = int(fields[9 + self._field_offset])
        cna1['cellular_prevalence'] = float(fields[10 + self._field_offset]) * self._cellularity

        cna2 = None
        # Stefan's comment on p values: The p-values correspond "to whether a
        # segment should be clonal or subclonal copynumber. We first fit a
        # clonal copynumber profile for the whole sample and then perform a
        # simple two-sided t-test twhere the null hypothesis is: A particular
        # segment is clonal. And the alternative: It is subclonal."
        #
        # Thus: if t-test falls below significance threshold, we push cna1 to
        # clonal frequency.
        if pval <= pval_threshold:
          cna2 = {}
          cna2['start'] = start
          cna2['end'] = end
          cna2['n_loci'] = 1
          cna2['logr'] = logr
          cna2['major_cn'] = int(fields[11 + self._field_offset])
          cna2['minor_cn'] = int(fields[12 + self._field_offset])
          cna2['cellular_prevalence'] = float(fields[13 + self._field_offset]) * self._cellularity
        else:
          cna1['cellular_prevalence'] = self._cellularity

        cna1['total_cn'] = self._compute_total_cn(cna1, cna2)
        cna1['mean_tcn'] = cna1['total_cn'] + (1 - self._cellularity) * chrom_copy_number(chrom)
        cna1['avail_cn'] = self._compute_avail_cn(cna1, cna2)
         
        cn_regions[chrom].append(cna1)
        # if cna2 is not None:
        #   cn_regions[chrom].append(cna2)
    return cn_regions
    
class BattenbergSmchetParser(BattenbergParser):
  def __init__(self, bb_filename, cellularity):
    super(BattenbergSmchetParser, self).__init__(bb_filename, cellularity)
    # SMC-Het Battenberg files lack the initial index column.
    self._field_offset = -1

def restricted_float(x):
  x = float(x)
  if x < 0.0 or x > 1.0:
    raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
  return x
  
def parse_cellularity(args):
    
  try:
    i = float(args.cellularity)
    cellularity = args.cellularity
    if args.cellularity > 1.0:
      print('Cellularity for %s is %s. Setting to 1.0.' % (args.cna_file, args.cellularity), file=sys.stderr)
      cellularity = 1.0
    else:
      pass
  except (ValueError, TypeError):
    with open(args.cellularity) as f:
      header = f.readline()
      for line in f.readlines():
        cellularity, ploidy = map(float, line.split('\t'))
    pass # can ignore the  process value or raise a error
  
  return cellularity

def main():
  parser = argparse.ArgumentParser(
    description='Create CNA input file for parser from Battenberg data',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('-f', '--cna-format', dest='cna_type', required=True, choices=('battenberg-smchet',),
    help='Type of CNA input')
  parser.add_argument('-c', '--cellularity', dest='cellularity', required=True,
    help='File containing newline-separated cellularity and ploidy')
  parser.add_argument('--cn', dest='cn_filename', default='cn.txt',
    help='Output destination for copy number')
  parser.add_argument('--mean-tcn', dest='mean_tcn_filename', default='mean_tcn.txt',
    help='Output destination for mean total copy number')
  parser.add_argument('--avail-cn', dest='avail_cn_filename', default='avail_cn.txt',
    help='Output destination for available copy number')
  parser.add_argument('cna_file')
  args = parser.parse_args()
  
  cellularity = parse_cellularity(args)

  if args.cna_type == 'battenberg-smchet':
    parser = BattenbergSmchetParser(args.cna_file, cellularity)
  else:
    raise Exception('Unknown input type')

  regions = parser.parse()
  cn = CopyNumberWriter(args.cn_filename)
  cn.write_cn(regions)
  mean_tcn = MeanTotalCopyNumberWriter(args.mean_tcn_filename)
  mean_tcn.write_mean_tcn(regions)
  avail_cn = AvailableCopyNumberWriter(args.avail_cn_filename)
  avail_cn.write_avail_cn(regions)

if __name__ == '__main__':
  main()
