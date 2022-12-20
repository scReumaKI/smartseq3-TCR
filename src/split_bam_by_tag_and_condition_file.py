#!/usr/bin/env python3
import pandas, argparse, pysam, tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_in')
    parser.add_argument('bam_tag_flag',type=str)
    parser.add_argument('condition_csv')
    parser.add_argument('condition_tag_col',type=str)
    parser.add_argument('condition_name_col',type=str)
    parser.add_argument('bam_out_prefix')
    parser.add_argument('--name_part_filter')
    o = parser.parse_args()

    CT = pandas.read_csv(o.condition_csv)
    if o.name_part_filter is not None:
        CT = CT.loc[CT[o.condition_name_col].str.contains(o.name_part_filter),:]
    tag_to_name = dict(zip(CT.loc[:,o.condition_tag_col], CT.loc[:,o.condition_name_col].astype(str)))

    bam_in = pysam.AlignmentFile(o.bam_in, "rb",check_sq=False)
    name_to_bam = {name:pysam.AlignmentFile(o.bam_out_prefix + name + '.bam', "wb", template=bam_in) for name in set(tag_to_name.values())}

    for aln in tqdm.tqdm(bam_in.fetch(until_eof=True)):
        barcode = aln.get_tag(o.bam_tag_flag)
        if barcode in tag_to_name.keys():
            name_to_bam[tag_to_name[barcode]].write(aln)
