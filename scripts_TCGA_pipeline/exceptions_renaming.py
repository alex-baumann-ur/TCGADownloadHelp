def file_rename_exceptions(df_input, tcga_study):
    df = df_input.copy()

    # 36 characters file ID
    df['File Suffix'] = ''

    # if file name has the 36 characters file ID as prefix and after that the suffix with a "." separated
    df.loc[(df['File Name'].str[36]=='.')&(df['File Name'].str[8]=='-')&(df['File Name'].str[13]=='-'), 'File Suffix'] = df['File Name'].str[36:]

    # if file name has a 36 characters file ID as prefix and after that the suffix with a "_" separated
    df.loc[(df['File Name'].str[36]=='_')&(df['File Name'].str[8]=='-'), 'File Suffix'] = '.' + df['File Name'].str[37:]

    # if file name starts with study abbreviation (such as TCGA_LUAD)
    for study in tcga_study:
        s1 = study.replace('_','-')
        s2 = study.replace('-','_')
        df.loc[((df['File Name'].str.startswith(s1))|(df['File Name'].str.startswith(s2)))&
               (df['File Name'].str[(len(study)+9)]=='-'), 'File Suffix'] = df['File Name'].str[(len(study)+37):]
    
    # if file name contains "WholeGenome" (GATK4 CNV)
    df.loc[(df['File Name'].str.contains('.WholeGenome.')), 'File Suffix'] = df['File Name'].str[28:]
    
    # for every other file that is not included in the previous rules
    df.loc[df['File Suffix']=='', 'File Suffix'] = '.' + df['File Name']

    return df
