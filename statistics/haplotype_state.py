import vcf 
import argparse
from record import Record, PhaseSet, ChromosomoHaplotype
from stats import PhaseSetStats, HapStats
import traceback


def get_phase_set_stats(template_phase_set:PhaseSet, phase_set:PhaseSet):
    prev_record: Record
    record: Record
    t_record: Record
    t_prev_record: Record 
    record_count = 0
    switch_error_count = 0
    mismatch_error_count = 0
    total_record = len(phase_set.records_idx)
    prev_switch_error = False
    last_record_pos = 0
    last_record_idx = 0
    first_record_idx = 0
    bnd_count = 0
    bnd_switched_error_count = 0
    bnd_missmatch_error_count = 0
    switchs = []
    mismatchs = []
    for record_pos in phase_set.records.keys():
        if record_pos not in template_phase_set.records.keys():
            continue
        record = phase_set.records[record_pos]
        if record.bnd:
            bnd_count = bnd_count + 1
            # continue
        record_count += 1
        # print(template_phase_set.records)

        t_record = template_phase_set.records[record_pos]
        if record_count == total_record:
            last_record_idx = record.idx
            last_record_pos = record.pos
        if record_count == 1:
            prev_record = record
            first_record_idx = record.idx
            t_prev_record = t_record
        else:
            switched = record.switched(prev_record)
            t_switched = t_record.switched(t_prev_record)
            if switched != t_switched:  # switch error
                if record_count > 2 and record_count < total_record:
                    switch_error_count += 1
                    switchs.append(record)
                    # if record.bnd:
                    #     bnd_switched_error_count += 1
                if prev_switch_error:   # mismatch error
                    mismatch_error_count += 1
                    switch_error_count -= 2
                    prev_switch_error = False
                    if len(switchs) > 0:
                        switchs.pop(-1)
                    if len(switchs) > 0:
                        switchs.pop(-1)
                    mismatchs.append(prev_record)
                        # bnd_switched_error_count -= 2
                        # bnd_missmatch_error_count += 1
                else:
                    prev_switch_error = True
            else:                       #no switch error for ajunct pos, reset
                prev_switch_error = False
            prev_record = record
            t_prev_record = t_record
    sws = []
    mis = []
    for i in switchs:
        if i.bnd == 1:
            bnd_switched_error_count +=1
            sws.append(i.pos)
    for i in mismatchs:
        if i.bnd == 1:
            bnd_missmatch_error_count +=1
            mis.append(i.pos)
    S50 = total_record
    N50 = last_record_pos - phase_set.starting_pos
    spaned_record = last_record_idx - first_record_idx + 1
    AN50 = N50/spaned_record * S50
    if bnd_switched_error_count < 0:
        bnd_switched_error_count = 0
    return AN50, S50, N50, switch_error_count, mismatch_error_count, spaned_record, bnd_count, bnd_switched_error_count, bnd_missmatch_error_count,sws, mis


def get_haplotype_stats_chromo(template_chromo:ChromosomoHaplotype, in_chromo:ChromosomoHaplotype, out, contig):
    template_phase_set:PhaseSet
    phase_set : PhaseSet
    template_phase_set = list(template_chromo.chromo_phase_set.values() )[0]
    chromo_snp_count = len(template_phase_set.records_idx)
    chromo_span = max(template_phase_set.records_idx) - min(template_phase_set.records_idx)
    hap_stats = HapStats(chromo_snp_count, chromo_span)
    index = 0
    for phase_set in in_chromo.chromo_phase_set.values():
        AN50, S50, N50, switch_error_count, mismatch_error_count, spanned_snp, bnd_count, bnd_switched_error_count, bnd_missmatch_error_count,switchs,mismatchs = get_phase_set_stats(template_phase_set, phase_set)
        phase_set_stats = PhaseSetStats(switch_error_count, mismatch_error_count, S50, N50, AN50, spanned_snp, bnd_count, bnd_switched_error_count, bnd_missmatch_error_count)
        if S50 < 2 or S50 - bnd_count == 0:
            continue
        hap_stats.insert_phase_set_stats(0, phase_set_stats)
        index += 1
        out.write("%s\t%d\t%d\t%d\t%d\t%.8f\t%.8f\t%d\t%d\t%d\t%s\t%s\n" % (contig, phase_set_stats.get_AN50(), phase_set_stats.get_N50(), phase_set_stats.get_phased_snp(), spanned_snp, phase_set_stats.get_switch_error(), phase_set_stats.get_mismatch_error(),
                                                                    phase_set_stats.bnd_count, phase_set_stats.bnd_switched_error_count, phase_set_stats.bnd_missmatch_error_count,switchs,mismatchs))
    out.write("%s\t%d\t%d\t%d\t%d\t%.8f\t%.8f\t%d\t%d\t%d\n" % (contig + "_total", hap_stats.get_AN50(), hap_stats.get_N50(), hap_stats.get_total_phased(), hap_stats.get_total_spanned(), hap_stats.get_switch_error(), hap_stats.get_mismatch_error(), hap_stats.get_bnd_count(), hap_stats.get_bnd_switch(), hap_stats.get_bnd_mismatch()))
    return hap_stats


def get_haplotype_stats(template_vcf:vcf.Reader, in_vcf:vcf.Reader, out):
    contigs = in_vcf.contigs.keys()
    hap_stats = HapStats()
    for contig in contigs:
        try: 
            template_vcf.fetch(contig)
            template_chromo = ChromosomoHaplotype(template_vcf, contig)
            in_chromo = ChromosomoHaplotype(in_vcf, contig)
            print(template_chromo)
            chromo_hap_stats = get_haplotype_stats_chromo(template_chromo, in_chromo, out, contig)
            hap_stats.insert_hap_stats(chromo_hap_stats)
        except Exception as e:
            traceback.print_exc()
        break
    out.write("%s\t%d\t%d\t%d\t%d\t%.8f\t%.8f\t%d\t%d\t%d\n" % ("total", hap_stats.get_AN50(), hap_stats.get_N50(), hap_stats.get_total_phased(), hap_stats.get_total_spanned(),hap_stats.get_switch_error(), hap_stats.get_mismatch_error(), hap_stats.get_bnd_count(), hap_stats.get_bnd_switch(), hap_stats.get_bnd_mismatch()))

def main():
    parser = argparse.ArgumentParser('phaseset_to_vcf.py')
    parser.add_argument('-t', '--template',  help='template vcf, indexed', required=True)
    parser.add_argument('-v', '--vcf', help='input vcf, indexed', required=True)
    parser.add_argument('-o', '--out', help='output stats', required=True)
    options = parser.parse_args()

    in_vcf = vcf.Reader(filename=options.vcf)
    template_vcf = vcf.Reader(filename=options.template)
    outf = open(options.out, 'w')
    outf.write("Chromosome\tAN50\tN50\tphased_snp\ttotal_snp\tswitch_error_rate\tmismatch_error_rate\tbnd_count\tbnd_mismatch\tbnd_switch\n")
    get_haplotype_stats(template_vcf, in_vcf, outf)    
        
    outf.close()

    return

    
if __name__ == '__main__':
    main()
