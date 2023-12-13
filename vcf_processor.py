import sys

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} <input_file_vcf> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

genotype_mapping = {}
sample_names = []
output_data = []

with open(input_file, 'r') as infile:
    for line in infile:
        if line.startswith('#CHROM'):
            parts = line.strip().split('\t')
            sample_names = parts[9:]

            header = ['#CHROM', 'POS', 'REF', 'ALT']
            for sample_name in sample_names:
                header.extend([sample_name, 'Depth', 'AlleDp'])
            output_data.append(header)

        if not line.startswith('#'):
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]

            sample_data = parts[9:]
            new_sample_data = []

            for data in sample_data:
                sample_parts = data.split(':')

                genotype = sample_parts[0]
                ad = sample_parts[1].split(',')
                depth = sample_parts[2]
                alle_dp = ','.join(ad)

                genotype_mapping[genotype] = f'{ref},{alt}'

                if genotype in ('0/0', '0|0'):
                    base_type = ref
                elif genotype in ('1/1', '1|1'):
                    base_type = alt
                elif genotype in ('./.', '.|.'):
                    base_type = 'NA'
                    depth = 'NA'
                    alle_dp = 'NA'
                else:
                    base_type = genotype_mapping.get(genotype, '')

                new_sample_info = [base_type, depth, alle_dp]
                new_sample_data.extend(new_sample_info)

            data_row = [chrom, pos, ref, alt] + new_sample_data
            output_data.append(data_row)

with open(output_file, 'w') as outfile:
    for data_row in output_data:
        outfile.write('\t'.join(data_row) + '\n')
