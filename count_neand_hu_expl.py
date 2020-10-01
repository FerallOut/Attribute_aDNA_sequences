# IMPORT NECESSARY MODULES
#******************************************************************************
# import 'glob' for Unix style pathname pattern expansion. 
# 'glob.glob('*.csv')' function will give you all the .csv files in a directory as a list.
import glob 

# will use 'pandas' library for manipulating data
import pandas as pd

# to retrieve the path to the current working directory & play sound when the script ends
import os

# to be able to track how long the script takes
from datetime import datetime 

# to use regex 
import re



start_time = datetime.now()    # initiate the time stamp

#FUNCTIONS
#******************************************************************************

# print-check:
  # how many df you have, how many rows and columns they have, 
  # how they look like, what kind of columns they have
def verif_print(possible_list_of_dataframes):
    if isinstance(possible_list_of_dataframes, list) ==True:
        for count in range(len(possible_list_of_dataframes)):
            return print(f'There are {len(possible_list_of_dataframes)} files', 
                             f'\n{possible_list_of_dataframes[count].head(5)}\n', 
                             f'DataFrame{count} dimensions/axes: {possible_list_of_dataframes[count].ndim}', 
                             f'DataFrame{count} rows :{possible_list_of_dataframes[count].shape[0]}', 
                             f'DataFrame{count} columns :{possible_list_of_dataframes[count].shape[1]}', 
                             sep='\n', end='\n*******************\n')
    elif isinstance(possible_list_of_dataframes, list) ==False:
        return print(f'{possible_list_of_dataframes.head(5)}\n', 
                     f'Dimensions/axes of the reference file: {possible_list_of_dataframes.ndim}', 
                     f'Number of genomic positions (=rows): {possible_list_of_dataframes.shape[0]}', 
                     f'Number of columns: {possible_list_of_dataframes.shape[1]}', 
                     sep='\n', end='\n*******************\n')

# define the path to the files you use
def open_path(file_name):
    path_ref=os.path.join(os.getcwd(), file_name)
    return path_ref

# the 'dtype=object' so that py already knwos what kind of data is looking at so it will not waste
#time and resources deciding on what data type it is working with
#don't add it if you have columns made of only one data type
# if you have mixed, this prevents py from exiting with an error

# the na_filter was necessary to keep the NA as empty space
# header prevents the erronious read of files (loss of first line)

# open reference file and name columns
def input_file(file_name):
    with open(open_path(file_name)) as my_fileRef:
        my_ref = pd.read_csv(my_fileRef, sep='\t', header=None, dtype='object', na_filter= False)
        my_ref.columns = ['chr', 'position', 'base_orig', 'base_deriv', 'group', 'empty']  #giving names to the columns to make it easier to $
        return my_ref
    
# evaluate the results you get
def evaluate(df_N, df_hu):
    if df_N.shape[0] < df_hu.shape[0]:
        return print(f'Neand:{df_N.shape[0]} - Hu:{df_hu.shape[0]} \nThe individual is most likely to be related to modern humans', end='\n*******************\n')
    elif df_N.shape[0] > df_hu.shape[0]:
        return print(f'\nNeand:{df_N.shape[0]} - Hu:{df_hu.shape[0]} \nThe individual is most likely to be related to Neanderthals', end='\n*******************\n')
    else:
        return print(f'\nNeand:{df_N.shape[0]} - Hu:{df_hu.shape[0]} \nEqual nr of sequences', end='\n*******************\n')

# function to determine which base to choose if there are multiple reads with different results
def most_common_base(my_list):
    counts,values = pd.Series(my_list).value_counts().values, pd.Series(my_list).value_counts().index
    df_results = pd.DataFrame(list(zip(values,counts)),columns=["value","count"])
    #print(f'the df for counts: {df_results}', end='\n******************\n')
    find_max_val = df_results['count'].max(axis=0, skipna = True)               # axis=0 takes the column!
    most_common_baseS = df_results.loc[(df_results['count'] == find_max_val), 'value']      # df.loc[mask, column]
    #print(f'Max value is \n {find_max_val}', end='\n******************\n')
    #print(f'The most common base is(are) \n {most_common_baseS}', end='\n******************\n')
    #print(f'df_results is \n {df_results}', end='\n******************\n')
    if len(most_common_baseS) == 1:
        return most_common_baseS[0]
    elif len(most_common_baseS) > 1:
        random_sampled = most_common_baseS.sample(n=1, axis=0) # Number of items from axis to randomly return
        return random_sampled.iloc[0]    #'iloc' function is indexing based on integer-location for selection by position 
    
# clean the 'base_read' column and find out the base that was read
def base_stripper(row):
    matchesList = re.findall('[ACGT]{1}', row['base_read'], flags=re.IGNORECASE)    #strips each row of the col of non-letters
    #print(matchesList)
    matchesList = [item.upper() for item in matchesList]    #make all uppercase
    if matchesList and len(matchesList) == 1:           #one base     #means empty list; For sequences, (strings, lists, tuples), empty sequences are false.
        #row['base_read'] = matchesList[0] 
        return matchesList[0]
    elif matchesList and len(matchesList) > 1:          #2 different bases, or same but different case
        base=most_common_base(matchesList)
        return base
    elif not matchesList:                           #(if list 'a' is empty, 'not a' is True)
        return row['base_reference']


    
# CODE
#******************************************************************************
# READ IN the reference file with genomic locations that differ between Neanderthals and modern humans

# Ref paper: Phylogenetic analysis in 'Nuclear DNA sequences from the Middle Pleistocene Sima de los Huesos hominins' - Meyer, Nature, 2016
# DL link: https://bioinf.eva.mpg.de/sima/

# CHANGE path if your reference file is named a different way
# mine is in folder 'new'
file_ref='new/duplicatedNOT_rows_tab_5k'
my_ref=input_file(file_ref)
#print(f'The reference file:')
#my_ref_verif=verif_print(my_ref)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# OPEN the files containing the sequenced data, after going through 'mpileup'

# CHANGE REGEX if your sequence files are named a different way
# mine are in folder 'new' and all are named 'jar00X.XXXXXX'
# read in all necessary '.mpileup' files as a list of dataframes
df = [pd.read_csv(files, sep='\t', dtype='object') for files in glob.glob(open_path('new/jar*'))]

# extract the name of each file to use later
file_name = [os.path.basename(f) for f in glob.glob(os.path.join(os.getcwd(), 'new/jar*'))]    
#print(file_name)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# loop through the sequence files and get needed data
for a in range(len(df)):

    # name each column of the sequence file
    df[a].columns = ["chr", "position", "base_reference", " nr_reads", "base_read", "quality"]  #giving names to the columns to make it easier to work with
    print(f'\n The DataFrame from sequence file "{file_name[a]}":')
    #df_verif=verif_print(df[a])
    
    # merge the reference file and one of the sequence file to get all needed data in one dataframe
    # this reduces the length of the file because it searches for rows with common ('chr','position')    
    common_df=pd.merge(left=my_ref[['chr','position','base_orig', 'base_deriv', 'group']], right=df[a][['chr','position', "base_reference","base_read"]], how='inner', on=['chr','position'])
    print(f'"{file_name[a]}": The dataframe with the common values (reference file and sequence):', f'{common_df.head(5)}', f'{common_df.head(5)}', sep='\n', end='\n******************\n')
    
    # drop na rows and those that only have one reading of '*'
    # replace missing data with "*" (not needed if data was checked for missing data and NA)
    common_df=common_df.replace('^\s*$', '*', regex=True)          
    print(f'"{file_name[a]}": The common dataframe after replacing missing data with "*":', f'{common_df.shape}', f'{common_df.head(5)}', sep='\n', end='\n***************\n') 
   
    indexNames = common_df[common_df['base_read'] == r'*'].index    #drop all rows in column 'base_read' that strictly have the value '*'
    common_df.drop(indexNames, inplace=True)
    print(f'"{file_name[a]}": The common dataframe after drop of "*" rows:', f'{common_df.shape}', f'{common_df.head(5)}', sep='\n', end='\n*************\n')
    
    #strip the 'base_read' column and assign the final base for each position
    common_df['base_read'] = common_df.apply(lambda row: base_stripper(row), axis=1)
    print(f'"{file_name[a]}": common dataframe after final filtering', f'{common_df.shape}', f'{common_df.head(5)}', sep='\n', end='\n***********\n') 
    
    #get all the rows for which ('base_deriv','base_read') are identical 
    mask=common_df[['base_deriv','base_read']].nunique(axis=1)==1 #on the df made out of the columns ('base_deriv','base_read') taken together,
                                                                  #I need to work only on the columns where these 2 columns have the same value (the sequenced nucleotide and the one in reference are the same)
                                                                  #so find the rows where it is true (==1) that the column values (axis =1) are equal

    # count only the rows you are interested in: check the next 2 blocks of code:
    df_N = common_df[mask & (common_df['group']=='Neandertal')] #now only save the rows where the 'group' values are only 'Nearderthal'
    print(f'"{file_name[a]}": The dataframe with the results only for Neanderthal', f'{df_N.head(5)}', f'{df_N.shape}', sep='\n', end ='\n********************\n')

    df_hu = common_df[mask & (common_df['group']=='human')]
    print(f'"{file_name[a]}": The dataframe with the results only for human', f'{df_hu.head(5)}', f'{df_hu.shape}', sep='\n', end ='\n*************\n')
   
    print(f'The results, wihout Denisovan contribution: ')
    # always put the dataframe for Neanderthals first!
    e=evaluate(df_N, df_hu)

    # delete this part if NOT interested in counting the Denisova contribution to the final result
    # same as before
    df_ND = common_df[mask & ((common_df['group']=='Neandertal-Denisova') | (common_df['group']=='Neandertal'))]    
    df_huD = common_df[mask & ((common_df['group']=='Human-Denisova') | (common_df['group']=='human'))]
    print(f'The results, when considering Denisovan contribution:')
    # always put the dataframe for Neanderthals first!
    e_D=evaluate(df_ND, df_huD)

end_time=datetime.now()
Duration=end_time-start_time   #gives you a timedelta object!
print ('Duration: {}'.format(end_time - start_time))

# !!!plays a sound when the script is done
# must have 'sox' installed (sudo apt install sox)

duration = 1   # seconds
freq = 440     # Hz
os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the unique names in the 'group' columns are: 
#Human-Denisova
#Neandertal-Human
#Denisova
#Neandertal
#Neandertal-Denisova
#human
#all
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~