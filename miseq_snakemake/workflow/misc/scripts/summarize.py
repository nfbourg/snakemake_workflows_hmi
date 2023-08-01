#!/usr/bin/env python
import pandas as pd

def remove_ends(row):
    indices = row.index[:-2]

    # remove beginning
    first_ind = indices[0]
    if row[first_ind] == 'N':
        for i in range(1,len(indices)):
            previous_ind = indices[i-1]
            row[previous_ind] = '-'
            next_nt = indices[i]
            if row[next_nt] != 'N':
                break

    # remove end
    last_ind = indices[-1]
    if row[last_ind] == 'N':
        for i in range(1,len(indices)):
            previous_ind = indices[(i*-1)]
            row[previous_ind] = '-'
            next_i = (i*-1)-1
            next_nt = indices[next_i]
            if row[next_nt] != 'N':
                break
    return row

def main(allele_txt_list,output,lowq):
    
    df = pd.DataFrame()

    for allele_loc in allele_txt_list:
        try:
            data = pd.read_csv(allele_loc, sep=' ', header = None)
            data.drop(data.columns[[0, 2, 4, 6, 8, 10, 11, 12, 13]], axis=1, inplace = True)
            data.columns = ["Read_ID", "Chromosome", "Coordinate", "Strand", "Base"]
            data.insert(0, "Sample", allele_loc.split('.')[0][8:])
            data.insert(0, "Read", allele_loc.split('.')[1])
            df = pd.concat([df, data])
        except pd.errors.EmptyDataError as e:
            if not lowq:
                assert(False)
            print(e)
            print(f"{allele_loc} has no data.")
            continue
    df2 = df[["Read_ID", "Coordinate", "Base"]]
    df2.drop_duplicates(subset=['Read_ID', 'Coordinate'], keep = 'first', inplace = True)
    df2.reset_index(inplace = True)
    data = df2.pivot(index='Read_ID', columns='Coordinate', values = "Base")
    df3 = df.loc[:, ['Read_ID','Sample', 'Read']]
    df3.drop_duplicates(inplace = True)
    sample_dict = pd.Series(df3.Sample.values,index=df3.Read_ID).to_dict()
    read_dict = pd.Series(df3.Read.values,index=df3.Read_ID).to_dict()

    data['Sample'] = data.index.map(sample_dict)
    data['Read'] = data.index.map(read_dict)
    
    if lowq:
        data = data.fillna('N')
        data.apply(lambda row: remove_ends(row), axis=1)
    else:
        data.dropna(inplace=True,axis=0)
    data.to_csv(output, sep="\t", index=False)

if __name__ == "__main__":
    allele_txt_list = snakemake.input.allele
    output = snakemake.output.res
    lowq = snakemake.params.lowq
    main(allele_txt_list,output,lowq)



