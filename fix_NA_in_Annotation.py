#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description='fix problem with N/A in Diamond annotation: remove the hypothetical proteins and recover taxonomy (taxid, superkingdom, kingdom and phylum) for the rest')
    parser.add_argument('-i', type=argparse.FileType(), help='input csv from diamond annotation')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output file without N/A')
    parser.add_argument('--tax_nodes', type=argparse.FileType(), help='nodes.nmd')
    parser.add_argument('--tax_names', type=argparse.FileType(), help='names.nmd, only scientific names')
    args = parser.parse_args()
    
    outp = args.output
    
    name_id_db = {}
    id_name_db = {}
    for name in args.tax_names:
        name = name.split('\t|\t')
        name_id_db[name[1]] = name[0]
        id_name_db[int(name[0])] = name[1] 
    
    kid_parent_rank = {}
    for ids_pair in args.tax_nodes:
        ids_pair = ids_pair.split('\t|\t')
        kid_parent_rank[int(ids_pair[0])] = [int(ids_pair[1]), ids_pair[2]] # ids_pair[2] - the rank of ids_pair[0]

    for record in args.i:
        record = record.strip().split('\t')
        if record[13] == 'N/A':
            if 'hypothetical protein' not in record[17]:
                short_lineage = {'sskingdoms': 0, 'skingdoms': 0, 'sphylums': 0}
                sci_name = record[17][record[17].rfind('[')+1 : record[17].rfind(']')] # scientific name between brakets []
                short_lineage['sscinames'] = sci_name
                if sci_name in name_id_db:
                    taxid = name_id_db[sci_name] # get taxonomic id of scientific names from []
                    short_lineage['staxids'] = taxid
                    taxid = int(taxid)

                    if taxid in kid_parent_rank:
                        while kid_parent_rank[taxid][0] != 1:
                            if kid_parent_rank[taxid][1] == 'superkingdom':
                                short_lineage['sskingdoms'] = id_name_db[taxid]
                            elif kid_parent_rank[taxid][1] == 'kingdom':
                                short_lineage['skingdoms'] = id_name_db[taxid]
                            elif kid_parent_rank[taxid][1] == 'phylum':
                                short_lineage['sphylums'] = id_name_db[taxid]

                            taxid = kid_parent_rank[taxid][0]

                    outp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(record[:12]), short_lineage['staxids'], \
                        short_lineage['sscinames'], short_lineage['sskingdoms'], short_lineage['skingdoms'], short_lineage['sphylums'], '\t'.join(record[17:])))
                else:
                    print('Could not find %s in names.nmd db' % (sci_name))
        else:
            outp.write('%s\n' % ('\t'.join(record)))
            # args.unmatched.write('%s\n' % (taxid))

if __name__ == '__main__':
    main()
