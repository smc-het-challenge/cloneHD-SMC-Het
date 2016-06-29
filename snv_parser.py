#!/usr/bin/env python2
from __future__ import print_function

# Requires PyVCF. To install: pip2 install pyvcf
import vcf
import argparse
import csv
from collections import defaultdict
import random
import sys

class VariantParser(object):
  def __init__(self):
    # Child classes must give the following variables sensible values in
    # constructor so that list_variants() works subsequently.
    self._cnas = None
    self._vcf_filename = None

  def list_variants(self):
    variants = self._filter(self._vcf_filename)
    variants_and_reads = []
    for variant in variants:
      ref_reads, alt_reads, total_reads = self._calc_read_counts(variant)
      variants_and_reads.append((variant, ref_reads, alt_reads, total_reads))
    return variants_and_reads

  def _calc_read_counts(self, variant):
    raise Exception('Not implemented -- use child class')

  def _parse_vcf(self, vcf_filename):
    vcfr = vcf.Reader(filename=vcf_filename)
    records = []
    for variant in vcfr:
      variant.CHROM = variant.CHROM.lower()
      # Some VCF dialects prepend "chr", some don't. Remove the prefix to
      # standardize.
      if variant.CHROM.startswith('chr'):
        variant.CHROM = variant.CHROM[3:]
      records.append(variant)
    return records

  def _is_good_chrom(self, chrom):
    # Ignore the following:
    #   * Variants unmapped ('chrUn') or mapped to fragmented chromosome ('_random')
    #   * Weird chromosomes from Mutect (e.g., "chr17_ctg5_hap1").
    #   * Mitochondrial ("mt" or "m"), which are weird
    #   * Sex chromosomes difficult to deal with, as expected frequency depends on
    #     whether patient is male or female, so ignore them for now. TODO: fix this.
    if chrom in [str(i) for i in range(1, 23)]:
      return True
    elif chrom in ['x', 'y']:
      return True
    else:
      return False

  def _does_variant_pass_filters(self, variant):
    if variant.FILTER is None:
      return True
    if len(variant.FILTER) > 0:
      # Variant failed one or more filters.
      return False
    return True

  def _filter(self, vcf_filename):
    variants = []

    all_variants = self._parse_vcf(vcf_filename)

    for variant in all_variants:
      if not self._is_good_chrom(variant.CHROM):
        continue
      if not self._does_variant_pass_filters(variant):
        continue
      variants.append(variant)
    return variants

  def _get_tumor_index(self, variant, tumor_sample=None):
    """Find the index of the tumor sample.

    Currently hardcodes tumour sample as the last column if name not specified.
    Might not always be true
    """
    if self._tumor_sample:
      tumor_is = [i for i, s in enumerate(variant.samples) if s.sample == tumor_sample]
      assert len(tumor_is) == 1, "Did not find tumor name %s in samples" % tumor_sample
      return tumor_is[0]
    else:
      # Don't make this -1, as some code assumes it will be >= 0.
      return len(variant.samples) - 1

class MutectSmchetParser(VariantParser):
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample

  def _calc_read_counts(self, variant):
    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    ref_reads = int(variant.samples[tumor_i]['AD'][0])
    alt_reads = int(variant.samples[tumor_i]['AD'][1])
    total_reads = ref_reads + alt_reads

    return (ref_reads, alt_reads, total_reads)

class CnaFormatter(object):
  def __init__(self, cna_confidence, cellularity, read_depth, read_length):
    self._cna_confidence = cna_confidence
    self._cellularity = cellularity
    self._read_depth = read_depth
    self._read_length = read_length

  def _max_reads(self):
    return 1e6 * self._read_depth

  def _find_overlapping_variants(self, chrom, cna, variants):
    overlapping = []

    start = cna['start']
    end = cna['end']
    for variant in variants:
      if chrom.lower() == variant['chrom'].lower():
        if start <= variant['pos'] <= end:
          overlapping.append(variant['snv_id'])
    return overlapping

  def _calc_ref_reads(self, cellular_prev, total_reads):
    vaf = cellular_prev / 2
    ref_reads = int((1 - vaf) * total_reads)
    return ref_reads
    
  def _calc_alt_reads(self, cellular_prev, total_reads):
    vaf = cellular_prev / 2
    alt_reads = int(vaf * total_reads)
    return alt_reads

  def _calc_total_reads(self, cellular_prev, locus_start, locus_end, new_cn):
    # Proportion of all cells carrying CNA.
    P = cellular_prev
    if new_cn == 2:
      # If no net change in copy number -- e.g., because (major, minor) went
      # from (1, 1) to (2, 0) -- force the delta_cn to be 1.
      delta_cn = 1.
      no_net_change = True
    else:
      delta_cn = float(new_cn - 2)
      no_net_change = False

    region_length = locus_end - locus_start + 1
    fn = (self._read_depth * region_length) / self._read_length

    # This is a hack to prevent division by zero (when delta_cn = -2). Its
    # effect will be to make d large.
    if P == 1.0:
      P = 0.999

    d = (delta_cn**2 / 4) * (fn * P * (2 - P)) / (1 + (delta_cn  * P) / 2)

    if no_net_change:
      # If no net change in CN occurred, the estimate was just based on BAFs,
      # meaning we have lower confidence in it. Indicate this lack of
      # confidence via d by multiplying it by (read length / distance between
      # common SNPs), with the "distance between common SNPs" taken to be 1000 bp.
      d *= (self._read_length / 1000.)

    # Cap at 1e6 * read_depth.
    return int(round(min(d, self._max_reads())))

  def _format_overlapping_variants(self, variants, maj_cn, min_cn):
      variants = [(snv_id, str(min_cn), str(maj_cn)) for snv_id in variants]
      return variants

  def _format_cnas(self, cnas, variants):
    log('Estimated read depth: %s' % self._read_depth)

    for chrom, chrom_cnas in cnas.items():
      for cna in chrom_cnas:
        overlapping_variants = self._find_overlapping_variants(chrom, cna, variants)
        total_reads = self._calc_total_reads(
          cna['cellular_prevalence'],
          cna['start'],
          cna['end'],
          cna['major_cn'] + cna['minor_cn'],
        )
        yield {
          'chrom': chrom,
          'start': cna['start'],
          'end': cna['end'],
          'major_cn': cna['major_cn'],
          'minor_cn': cna['minor_cn'],
          'cellular_prevalence': cna['cellular_prevalence'],
          'ref_reads': self._calc_ref_reads(cna['cellular_prevalence'], total_reads),
          'alt_reads': self._calc_alt_reads(cna['cellular_prevalence'], total_reads),
          'total_reads': total_reads,
          'overlapping_variants': self._format_overlapping_variants(overlapping_variants, cna['major_cn'], cna['minor_cn']),
        }

  def _merge_variants(self, cna1, cna2):
    cna1_variant_names = set([v[0] for v in cna1['overlapping_variants']])
    for variant in cna2['overlapping_variants']:
      variant_name = variant[0]
      if variant_name not in cna1_variant_names:
        cna1['overlapping_variants'].append(variant)
    else:
        # If variant already in cna1's list, ignore it. This should only occur
        # if two subclonal CNAs have close to 0.5 frequency each. In this case,
        # we lose information about major/minor status of the cna2 relative to
        # its SNVs.
        log('%s already in %s' % (variant, cna1['cna_id']))

  # CNAs with similar a/d values should not be free to move around the
  # phylogeny independently, and so we merge them into a single entity. We may
  # do the same with SNVs bearing similar frequencies later on.
  def format_and_merge_cnas(self, cnas, variants):
    formatted = list(self._format_cnas(cnas, variants))
    formatted.sort(key = lambda f: f['cellular_prevalence'])
    if len(formatted) == 0:
      return []

    merged, formatted = formatted[:1], formatted[1:]
    merged[0]['cna_id'] = 'c0'
    counter = 1

    cellularity = find_cellularity(cnas)

    for current in formatted:
      last = merged[-1]

      # Only merge CNAs if they're clonal. If they're subclonal, leave them
      # free to move around the tree.
      if current['cellular_prevalence'] == last['cellular_prevalence'] == cellularity:
        # Merge the CNAs.
        log('Merging %s_%s and %s_%s' % (current['chrom'], current['start'], last['chrom'], last['start']))
        last['total_reads'] = current['total_reads'] + last['total_reads']
        last['ref_reads'] = self._calc_ref_reads(last['cellular_prevalence'], last['total_reads'])
        last['alt_reads'] = self._calc_alt_reads(last['cellular_prevalence'], last['total_reads'])
        self._merge_variants(last, current)
      else:
        # Do not merge the CNAs.
        current['cna_id'] = 'c%s' % counter
        merged.append(current)
        counter += 1

    for cna in merged:
      cna['ref_reads'] = int(round(cna['ref_reads'] * self._cna_confidence))
      cna['alt_reads'] = int(round(cna['alt_reads'] * self._cna_confidence))
      cna['total_reads'] = int(round(cna['total_reads'] * self._cna_confidence))

    return merged

class VariantFormatter(object):
  def __init__(self):
    self._counter = 0

  def _split_types(self, genotype):
    types = [int(e) for e in genotype.split('/')]
    if len(types) != 2:
      raise Exception('Not diploid: %s' % types)
    return types

  def _calc_ref_freq(self, ref_genotype, error_rate):
    types = self._split_types(ref_genotype)
    num_ref = len([t for t in types if t == 0])
    freq = (num_ref / 2) - error_rate
    if freq < 0:
      freq = 0.0
    if freq > 1:
      raise Exception('Nonsensical frequency: %s' % freq)
    return freq

  def format_variants(self, variant_list, error_rate):
    variant_list.sort(key = lambda v: variant_key(v[0]))
    
    for variant, alt_reads, total_reads in variant_list:
      snv_id = 's%s' % self._counter
      chrom, pos = variant_key(variant)
      variant_name = '%s_%s' % (chrom, pos)

      # TODO: switch back to using calc_ref_freq() when we no longer want mu_r
      # and mu_v fixed.
      # This is mu_r in PhyloWGS.
      expected_ref_freq = 1 - error_rate
      if chrom in ('x', 'y', 'm'):
        # Haploid, so should only see non-variants when sequencing error
        # occurred. Note that chrY and chrM are always haploid; chrX is haploid
        # only in men, so script must know sex of patient to choose correct
        # value. Currently, I just assume that all data comes from men.
        #
        # This is mu_v in PhyloWGS.
        expected_alt_freq = error_rate
      else:
        # Diploid, so should see variants in (0.5 - error_rate) proportion of
        # reads.
        #
        # This is mu_v in PhyloWGS.
        expected_alt_freq = 0.5 - error_rate

      yield {
        'snv_id': snv_id,
        'chrom': chrom,
        'pos': pos,
        'variant_name': variant_name,
        # 'ref_reads': ref_reads,
        'alt_reads': alt_reads,
        'total_reads': total_reads,
        'expected_ref_freq': expected_ref_freq,
        'expected_alt_freq': expected_alt_freq,
      }
      self._counter += 1

def restricted_float(x):
  x = float(x)
  if x < 0.0 or x > 1.0:
    raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
  return x

def variant_key(var):
  chrom = var.CHROM
  if chrom == 'x':
    chrom = 23
  elif chrom == 'y':
    chrom = 24
  elif chrom.isdigit():
    chrom = int(chrom)
  else:
    chrom = 999
  return (chrom, var.POS)

def find_cellularity(cnas):
  max_cellular_prev = 0
  for chrom, chrom_regions in cnas.items():
    for cnr in chrom_regions:
      if cnr['cellular_prevalence'] > max_cellular_prev:
        max_cellular_prev = cnr['cellular_prevalence']
  return max_cellular_prev

class VariantAndCnvGroup(object):
  def __init__(self):
    self._cn_regions = None
    self._cellularity = None

  def add_variants(self, variants_and_reads):
    self._variants_and_reads = variants_and_reads
    # Estimate read depth before any filtering of variants is performed, in
    # case no SNVs remain afterward.
    self._estimated_read_depth = self._estimate_read_depth()

  def add_cnas(self, cn_regions):
    self._cn_regions = cn_regions
    self._cellularity = find_cellularity(self._cn_regions)

  def has_cnas(self):
    return self._cn_regions is not None

  def _filter_variants_outside_regions(self, variants_and_reads, regions):
    filtered = []

    for variant, ref_reads, alt_reads, total_reads in self._variants_and_reads:
      for region in regions[variant.CHROM]:
        if region['start'] <= variant.POS <= region['end']:
          filtered.append((variant, ref_reads, alt_reads, total_reads))
          break

    return filtered

  def _is_region_normal_cn(self, region):
    return region['major_cn'] == region['minor_cn'] == 1

  def _print_variant_differences(self, before, after, before_label, after_label):
    before = set([var for (var, _, _, _) in before])
    after = set([var for (var, _, _, _) in after])
    log('%s=%s %s=%s delta=%s' % (before_label, len(before), after_label, len(after), len(before) - len(after)))


    assert after.issubset(before)
    removed = list(before - after)
    removed.sort(key = variant_key)

    for var in removed:
      var_name = '%s_%s' % (var.CHROM, var.POS)
      for region in self._cn_regions[var.CHROM]:
        if region['start'] <= var.POS <= region['end']:
          region_type = (self._is_region_normal_cn(region) and 'normal') or 'abnormal'
          log('%s\t[in %s-CN region chr%s(%s, %s)]' % (var_name, region_type, var.CHROM, region['start'], region['end']))
          break
      else:
        log('%s\t[outside all regions]' % var_name)

  def retain_only_variants_in_normal_cn_regions(self):
    if not self.has_cnas():
      raise Exception('CN regions not yet provided')

    normal_cn = defaultdict(list)

    for chrom, regions in self._cn_regions.items():
      for region in regions:
        if self._is_region_normal_cn(region) and region['cellular_prevalence'] == self._cellularity:
          normal_cn[chrom].append(region)

    filtered = self._filter_variants_outside_regions(self._variants_and_reads, normal_cn)
    self._print_variant_differences(self._variants_and_reads, filtered, 'all_variants', 'only_normal_cn')
    self._variants_and_reads = filtered

  def _filter_multiple_abnormal_cn_regions(self, regions):
    good_regions = defaultdict(list)
    for chrom, reg in regions.items():
      idx = 0
      while idx < len(reg):
        region = reg[idx]

        # Accept clonal regions unconditonally, whether normal or abnormal CN.
        if region['cellular_prevalence'] == self._cellularity:
          good_regions[chrom].append(region)
          idx += 1

        else:
          regions_at_same_coords = [region]

          i = idx + 1
          while i < len(reg) and reg[i]['start'] == region['start']:
            # Battenerg either has entire regions at same coords, or they have
            # no overlap whatsoever. Thus, this assertion maintains sanity for
            # Battenberg, but may fail for other CN callers.
            assert reg[i]['end'] == region['end']
            regions_at_same_coords.append(reg[i])
            i += 1

          abnormal_regions = [r for r in regions_at_same_coords if not self._is_region_normal_cn(r)]
          # In Battenberg, either one region is normal and the other abnormal,
          # or both are abnormal.
          # In TITAN, only one abnormal region will be listed, without a
          # corresponding normal region.
          # Ignore normal region(s) and add only one abnormal one. We do
          # this so PWGS can recalculate the frequencies based on the
          # major/minor CN of the region, according to the a & d values we will
          # assign to the region.
          if len(abnormal_regions) == 1:
            good_regions[chrom].append(abnormal_regions[0])
          else:
            # Ignore CNA regions with multiple abnormal CN states, as we don't
            # know what order the CN events occurred in.
            log('Multiple abnormal regions: chrom=%s %s' % (chrom, abnormal_regions))
          idx += len(regions_at_same_coords)

    return good_regions

  def exclude_variants_in_subclonal_cnas(self):
    # Battenberg:
    #   Five possible placements for variant in Battenberg according to CN records:
    #   1 record:
    #     That record has normal CN: include
    #     That record has abnormal CN: include
    #   2 records:
    #     One record is normal CN, one record is abnormal CN: include
    #     Both records are abnormal CN: exclude (as we don't know what order the CN events occurred in)
    # TITAN:
    #   In output seen to date, TITAN will only list one record per region. If
    #   the CN state is abnormal and clonal_frac < 1, this implies the
    #   remainder of the region will be normal CN. Multiple abnormal records
    #   for the same region are likely possible, but I haven't yet seen any.
    #   Regardless, when they occur, they should be properly handled by the
    #   code.
    if not self.has_cnas():
      raise Exception('CN regions not yet provided')

    good_regions = self._filter_multiple_abnormal_cn_regions(self._cn_regions)
    # If variant isn't listed in *any* region: exclude (as we suspect CNA
    # caller didn't know what to do with the region).
    filtered = self._filter_variants_outside_regions(self._variants_and_reads, good_regions)
    self._print_variant_differences(self._variants_and_reads, filtered, 'all_variants', 'outside_subclonal_cn')
    self._variants_and_reads = filtered

  def format_variants(self, sample_size, error_rate, priority_snvs):
    if sample_size is None:
      sample_size = len(self._variants_and_reads)
    random.shuffle(self._variants_and_reads)

    subsampled, remaining = [], []
    
    for variant, ref_reads, alt_reads, total_reads in self._variants_and_reads:
      if len(subsampled) < sample_size and (variant.CHROM, variant.POS) in priority_snvs:
        subsampled.append((variant, alt_reads, total_reads))
    else:
        remaining.append((variant, alt_reads, total_reads))

    assert len(subsampled) <= sample_size
    needed = sample_size - len(subsampled)
    subsampled = subsampled + remaining[:needed]
    nonsubsampled = remaining[needed:]

    formatter = VariantFormatter()
    subsampled = list(formatter.format_variants(subsampled, error_rate))
    nonsubsampled = list(formatter.format_variants(nonsubsampled, error_rate))

    return (subsampled, nonsubsampled)

  def write_variants(self, variants, outfn):
    with open(outfn, 'w') as outf:
      for variant in variants:
        vals = (
          'chrom',
          'pos',
          'alt_reads',
          'total_reads',
        )
        vals = [variant[k] for k in vals]
        print('\t'.join([str(v) for v in vals]), file=outf)

  def _estimate_read_depth(self):
    read_sum = 0
    if len(self._variants_and_reads) == 0:
      default_read_depth = 50
      log('No variants available, so fixing read depth at %s.' % default_read_depth)
      return default_read_depth
    for variant, ref_reads, alt_reads, total_reads in self._variants_and_reads:
      read_sum += total_reads
    return float(read_sum) / len(self._variants_and_reads)

  def write_cnas(self, variants, outfn, cna_confidence, read_length):
    abnormal_regions = {}
    filtered_regions = self._filter_multiple_abnormal_cn_regions(self._cn_regions)
    for chrom, regions in filtered_regions.items():
      abnormal_regions[chrom] = [r for r in regions if not self._is_region_normal_cn(r)]

    with open(outfn, 'w') as outf:
      # print('\t'.join(('cna', 'a', 'd', 'snvs')), file=outf)
      formatter = CnvFormatter(cna_confidence, self._cellularity, self._estimated_read_depth, read_length)
      for cna in formatter.format_and_merge_cnas(abnormal_regions, variants):
        overlapping = [','.join(o) for o in cna['overlapping_variants']]
        vals = (
          cna['cna_id'],
          str(cna['alt_reads']),
          str(cna['total_reads']),
          ';'.join(overlapping),
        )
        print('\t'.join(vals), file=outf)

def log(msg):
  if log.verbose:
    print(msg, file=sys.stderr)
log.verbose = False

class CnaParser(object):
  def __init__(self, cn_filename):
    self._cn_filename = cn_filename

  def parse(self):
    cn_regions = defaultdict(list)

    with open(self._cn_filename) as cnf:
      reader = csv.DictReader(cnf, delimiter='\t')
      for record in reader:

        chrom = record['chromosome']
        del record['chromosome']
        for key in ('start', 'end', 'major_cn', 'minor_cn'):
          record[key] = int(record[key])
        record['cellular_prevalence'] = float(record['cellular_prevalence'])
        cn_regions[chrom].append(record)

    # Ensure CN regions are properly sorted, which we later rely on when
    # filtering out regions with multiple abnormal CN states.
    for chrom, regions in cn_regions.items():
      cn_regions[chrom] = sorted(regions, key = lambda r: r['start'])

    return cn_regions

def parse_priority_snvs(priority_snv_filename):
  if priority_snv_filename is None:
    return set()
  with open(priority_snv_filename) as priof:
    lines = priof.readlines()

  priority_snvs = set()
  for line in lines:
    chrom, pos = line.strip().split('_', 1)
    priority_snvs.add((chrom.lower(), int(pos)))
  return set(priority_snvs)

def main():
  parser = argparse.ArgumentParser(
    description='Create snv_data.txt and cna_data.txt input files for cloneHD from VCF and CNA data.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('-e', '--error-rate', dest='error_rate', type=restricted_float, default=0.001,
    help='Expected error rate of sequencing platform')
  parser.add_argument('-s', '--sample-size', dest='sample_size', type=int,
    help='Subsample SNVs to reduce cloneHD runtime')
  parser.add_argument('-P', '--priority-snvs', dest='priority_snv_filename',
    help='File containing newline-separated list of SNVs in "<chr>_<locus>" format to prioritize for inclusion')
  parser.add_argument('--cnas', dest='cna_file',
    help='Path to CNA list created with parse_cnas.py')
  parser.add_argument('--only-normal-cn', dest='only_normal_cn', action='store_true', default=False,
      help='Only output variants lying in normal CN regions. Do not output CNA data directly.')
  parser.add_argument('--output-cnas', dest='output_cnas', default='cna_data.txt',
    help='Output destination for CNAs')
  parser.add_argument('--output-snvs', dest='output_snvs', default='snv_data.txt',
    help='Output destination for variants')
  parser.add_argument('-v', '--variant-type', dest='input_type', required=True, choices=('mutect-smchet',),
    help='Type of VCF file')
  parser.add_argument('--tumor-sample', dest='tumor_sample',
    help='Name of the tumor sample in the input VCF file. Defaults to last sample if not specified.')
  parser.add_argument('--cna-confidence', dest='cna_confidence', type=restricted_float, default=0.5,
    help='Confidence in CNAs. Set to < 1 to scale "d" values used in CNA output file')
  parser.add_argument('--read-length', dest='read_length', type=int, default=100,
    help='Approximate length of reads. Used to calculate confidence in CNA frequencies')
  # parser.add_argument('--nonsubsampled-variants', dest='output_nonsubsampled_variants',
  #   help='If subsampling, write nonsubsampled variants to separate file, in addition to subsampled variants')
  # parser.add_argument('--nonsubsampled-variants-cnas', dest='output_nonsubsampled_variants_cnas',
  #   help='If subsampling, write CNAs for nonsubsampled variants to separate file')
  parser.add_argument('--verbose', dest='verbose', action='store_true')
  parser.add_argument('vcf_file')
  args = parser.parse_args()

  log.verbose = args.verbose

  # Fix random seed to ensure same set of SNVs chosen when subsampling on each
  # invocation.
  random.seed(1)

  grouper = VariantAndCnvGroup()
  if args.input_type == 'mutect-smchet':
    variant_parser = MutectSmchetParser(args.vcf_file, args.tumor_sample)
  variants_and_reads = variant_parser.list_variants()
  grouper.add_variants(variants_and_reads)

  if args.cna_file:
    cna_parser = CnaParser(args.cna_file)
    cn_regions = cna_parser.parse()
    grouper.add_cnas(cn_regions)

  if args.only_normal_cn:
    grouper.retain_only_variants_in_normal_cn_regions()
  elif grouper.has_cnas():
    grouper.exclude_variants_in_subclonal_cnas()

  priority_snvs = parse_priority_snvs(args.priority_snv_filename)
  subsampled_vars, nonsubsampled_vars = grouper.format_variants(args.sample_size, args.error_rate, priority_snvs)
  grouper.write_variants(subsampled_vars, args.output_snvs)
  # if args.output_nonsubsampled_variants:
  #   grouper.write_variants(nonsubsampled_vars, args.output_nonsubsampled_variants)

  if not args.only_normal_cn and grouper.has_cnas():
    grouper.write_cnas(subsampled_vars, args.output_cnas, args.cna_confidence, args.read_length)
    # if args.output_nonsubsampled_variants and args.output_nonsubsampled_variants_cnas:
    #   grouper.write_cnas(nonsubsampled_vars, args.output_nonsubsampled_variants_cnas, args.cna_confidence, args.read_length)
  else:
    # Write empty CNA file.
    with open(args.output_cnas, 'w'):
      pass

if __name__ == '__main__':
  main()
