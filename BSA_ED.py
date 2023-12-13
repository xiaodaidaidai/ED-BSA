import argparse
from collections import defaultdict
import math
from statistics import mean

def parse_vcf_header(vcf_file):
    col1 = 0
    col2 = 0

    for line in vcf_file:
        if line.startswith("#CHROM"):
            header = line.strip().split("\t")
            col1 = header.index(args.bulk1)
            col2 = header.index(args.bulk2)
            return col1, col2
    return col1, col2

def calculate_edpower(ad1, ad2, power):
    dp1 = sum(ad1)
    dp2 = sum(ad2)

    if dp1 == 0 or dp2 == 0 or dp1 < args.minDP or dp1 > args.maxDP or dp2 < args.minDP or dp2 > args.maxDP:
        return None

    freq1 = [ad / dp1 for ad in ad1]
    freq2 = [ad / dp2 for ad in ad2]

    sqr_diff = [(freq1[i] - freq2[i]) ** 2 for i in range(len(ad1))]
    ed = math.sqrt(sum(sqr_diff))

    ed_power = ed ** power
    return ed_power

def main():
    with open(args.vcf, 'r') as vcf_file, open(args.outpre + f".snp_EDpower{args.power}.tsv", 'w') as snp_output:
        col1, col2 = parse_vcf_header(vcf_file)

        snp_output.write("CHROM\tPOS\t{0}_AD\t{1}_AD\tED{2}\n".format(args.bulk1, args.bulk2, args.power))
        print(args.bulk1 + ": " + str(col1) + "  " + args.bulk2 + ": " + str(col2))
        ed_dict = defaultdict(dict)

        for line in vcf_file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")

            ad_field = fields[8].split(":").index("AD")
            ad1 = [int(ad) for ad in fields[col1].split(":")[ad_field].split(",")]
            ad2 = [int(ad) for ad in fields[col2].split(":")[ad_field].split(",")]

            ed_power = calculate_edpower(ad1, ad2, args.power)

            if ad1.count(0) == 1 and ad2.count(0) == 1 and ed_power == 0:
                continue

            if ed_power is not None:
                ed_dict[fields[0]][fields[1]] = ed_power
                snp_output.write(f"{fields[0]}\t{fields[1]}\t{','.join(map(str, ad1))}\t{','.join(map(str, ad2))}\t{ed_power}\n")

    with open(args.outpre + f".sliding_EDpower{args.power}.tsv", 'w') as sliding_output:
        sliding_output.write("CHROM\tPOS\tSTART\tEND\tnSNP\tmean_ED{0}\n".format(args.power))

        for chrom, pos_dict in ed_dict.items():
            pos_list = sorted(pos_dict.keys())

            max_pos = max(int(pos) for pos in pos_list)

            window = args.window
            step = args.step

            window = window * 1000
            step = step * 1000

            if int(max_pos) <= window:
                mean_ed = mean(ed_dict[chrom][pos] for pos in pos_list)
                sliding_output.write(f"{chrom}\t{int(max_pos) / 2}\t1\t{max_pos}\t{len(pos_list)}\t{mean_ed}\n")
            else:
                win_edges = list(range(window, int(max_pos) + 1, step))

                for edge in win_edges:
                    start = edge - window + 1
                    sub_pos = [pos for pos in pos_list if start <= int(pos) < edge]

                    if sub_pos:
                        mean_ed = mean(ed_dict[chrom][pos] for pos in sub_pos)
                    else:
                        mean_ed = "nan"
                    mean_pos = edge - window / 2
                    sliding_output.write(f"{chrom}\t{int(mean_pos):d}\t{start}\t{edge}\t{len(sub_pos)}\t{mean_ed}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate ED power from VCF file")
    parser.add_argument('-vcf', required=True, help='Input VCF file')
    parser.add_argument('-bulk1', required=True, help='Sample name for bulk1')
    parser.add_argument('-bulk2', required=True, help='Sample name for bulk2')
    parser.add_argument('-power', type=int, default=5, help='ED power (default is 5)')
    parser.add_argument('-minDP', type=int, default=4, help='Minimum depth of variation (per sample) (default is 4)')
    parser.add_argument('-maxDP', type=int, default=250, help='Maximum depth of variation (per sample) (default is 250)')
    parser.add_argument('-window', type=int, default=2000, help='Window size (kilobases) (default is 2000)')
    parser.add_argument('-step', type=int, default=100, help='Step size (kilobases) (default is 100)')
    parser.add_argument('-outpre', required=True, help='Output file prefix')

    args = parser.parse_args()
    main()
